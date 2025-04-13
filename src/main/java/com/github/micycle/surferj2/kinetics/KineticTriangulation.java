package com.github.micycle.surferj2.kinetics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeSubdivision;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeTriangle;

/**
 * The static structure of the triangulation at time t=0. The algorithm now
 * needs to determine how this structure will evolve over time. This evolution
 * is driven by events, which signify topological changes in the wavefront.
 * <P>
 * Uses a kinetic data structure to witness events: Triangulate the
 * not-yet-swept plane; triangles witness events
 */
public class KineticTriangulation {

	/*-
	 * Relevant C++ Classes:
		KineticTriangulation: Overall manager and initialization logic.
		KineticTriangle: The core node storing connectivity and vertex/edge references.
		WavefrontVertex: Represents the kinetic vertices.
		WavefrontEdge: Represents the kinetic edges corresponding to original constraints.
		BasicTriangulation: Performs the initial Constrained Delaunay Triangulation.
		BasicInput: Holds the initial polygon data.
	 */

	private final List<KineticTriangle> triangles = new ArrayList<>();
	private final List<WavefrontVertex> vertices = new ArrayList<>(); // Includes INFINITE_VERTEX potentially
	private final List<WavefrontEdge> wavefrontEdges = new ArrayList<>();

	// Lookup maps built during initialization
	private final Map<Coordinate, WavefrontVertex> vertexMap = new HashMap<>();
	private final Map<CanonicalSegment, WavefrontEdge> edgeMap = new HashMap<>();
	// form initial wavefront
	private final Set<CanonicalSegment> constraintSegments = new HashSet<>();

	// Add lists to hold newly created elements during an event step if needed
	private final List<WavefrontVertex> newVerticesThisStep = new ArrayList<>();
	private final List<WavefrontEdge> newEdgesThisStep = new ArrayList<>();

	public List<KineticTriangle> getTriangles() {
		return Collections.unmodifiableList(triangles);
	}

	public List<WavefrontVertex> getVertices() {
		return Collections.unmodifiableList(vertices);
	}

	public List<WavefrontEdge> getWavefrontEdges() {
		return Collections.unmodifiableList(wavefrontEdges);
	}

	/**
	 * Initializes the KineticTriangulation from a JTS Polygon.
	 *
	 * @param inputPolygon    The input polygon (may have holes).
	 * @param geometryFactory The geometry factory to use.
	 */
	public KineticTriangulation(Polygon inputPolygon, GeometryFactory geometryFactory) {
		if (inputPolygon == null || inputPolygon.isEmpty()) {
			// Or throw exception? For now, allow empty init.
			System.err.println("Warning: Initializing KineticTriangulation with empty polygon.");
			return;
		}
		if (!inputPolygon.isValid()) {
			throw new IllegalArgumentException("Input polygon is not valid according to JTS rules.");
		}

		// 1. Store Constraint Segments
		populateConstraintSegments(inputPolygon);

		// 2. Perform Constrained Delaunay Triangulation
		// Use QuadEdge implementation for better triangle access
		var dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(inputPolygon);
		QuadEdgeSubdivision subdivision = dtb.getSubdivision();

		// 3. Build Initial Data Structures (Triangles, Vertices, Edges)
		Map<QuadEdgeTriangle, KineticTriangle> qetToKt = new HashMap<>();
		List<QuadEdgeTriangle> jtsTriangles = QuadEdgeTriangle.createOn(subdivision);

		// First Pass: Create KineticTriangle objects and WavefrontVertex placeholders
		for (QuadEdgeTriangle jtsTriangle : jtsTriangles) {
			if (!isTriangleInside(jtsTriangle, inputPolygon)) {
				continue; // Skip triangles outside the polygon area (e.g., in holes or outside convex
							// hull)
			}

			KineticTriangle kt = new KineticTriangle();
			triangles.add(kt);
			qetToKt.put(jtsTriangle, kt); // map jts geom triangle to kinetic triangle.

			for (int i = 0; i < 3; i++) {
				Coordinate coord = jtsTriangle.getCoordinate(i);
				// Create vertex placeholder if not seen before
				WavefrontVertex wv = vertexMap.computeIfAbsent(coord, k -> {
					WavefrontVertex newWv = new WavefrontVertex(k);
					this.vertices.add(newWv);
					return newWv;
				});
				kt.setVertex(i, wv);
			}
		}

		// Second Pass: Link Neighbors and Create/Link WavefrontEdges
		// --- Pass 2: Link Neighbors and Create/Identify WavefrontEdges ---
		for (Map.Entry<QuadEdgeTriangle, KineticTriangle> entry : qetToKt.entrySet()) {
			QuadEdgeTriangle jtsTriangle = entry.getKey();
			KineticTriangle currentKt = entry.getValue();

			for (int i = 0; i < 3; i++) {
				Coordinate p0 = jtsTriangle.getCoordinate(i);
				Coordinate p1 = jtsTriangle.getCoordinate((i + 1) % 3);
				if (p0 == null || p1 == null || p0.equals2D(p1)) {
					continue; // Skip null or degenerate edges
				}

				CanonicalSegment edgeSegment = new CanonicalSegment(p0, p1);
				int edgeIndexInTriangle = (i + 2) % 3; // Edge i is opposite vertex (i+2)%3

				boolean isConstraint = constraintSegments.contains(edgeSegment);
				QuadEdgeTriangle neighborQet = findNeighbor(jtsTriangle, i, jtsTriangles);
				KineticTriangle neighborKt = (neighborQet != null) ? qetToKt.get(neighborQet) : null;

				if (isConstraint) {
					// This edge is on the input boundary. Create/get the WavefrontEdge.
					WavefrontEdge wfEdge = createOrGetWavefrontEdge(edgeSegment);
					currentKt.setWavefront(edgeIndexInTriangle, wfEdge);
					// Check if neighbor also exists and borders this constraint
					if (neighborKt != null) {
						// This might happen if the constraint is very thin or has collinear points
						// The neighbor should also see this as a constraint
						System.err.println("Warning: Constraint edge " + edgeSegment + " has an internal neighbor triangle " + neighborKt.id);
					}
				} else {
					// This edge is internal to the polygon's triangulation area.
					if (neighborKt != null) {
						// Link neighbors
						currentKt.setNeighbor(edgeIndexInTriangle, neighborKt);
					} else {
						// Internal edge with no valid neighbor inside the polygon area.
						// This means the edge must lie on the boundary (e.g., hole boundary
						// that wasn't explicitly marked as constraint or outer boundary near convex
						// hull).
						// Treat it as a constraint.
						if (!constraintSegments.contains(edgeSegment)) {
							System.err
									.println("Warning: Internal edge " + edgeSegment + " borders outside, treating as constraint for triangle " + currentKt.id);
							// Add to constraints if this happens, maybe triangulation issue?
							// constraintSegments.add(edgeSegment); // Use cautiously
						}
						WavefrontEdge wfEdge = createOrGetWavefrontEdge(edgeSegment);
						currentKt.setWavefront(edgeIndexInTriangle, wfEdge);
					}
				}
			}
		}

		// Third Pass: Link WavefrontEdges to their WavefrontVertices
		for (WavefrontEdge wfEdge : this.wavefrontEdges) {
			Coordinate c0 = wfEdge.canonicalSegment.getP0();
			Coordinate c1 = wfEdge.canonicalSegment.getP1();
			WavefrontVertex v0 = vertexMap.get(c0);
			WavefrontVertex v1 = vertexMap.get(c1);

			if (v0 == null || v1 == null) {
				throw new IllegalStateException("Vertex not found for WavefrontEdge segment: " + wfEdge.canonicalSegment);
			}

			// Determine orientation relative to the incident triangle
			KineticTriangle tri = wfEdge.getIncidentTriangle();
			if (tri == null) {
				// This might happen if an edge is only bordered by outside triangles
				// Find the *other* triangle potentially if needed, or handle based on algorithm
				// reqs.
				// For now, search all triangles (less efficient)
				tri = findTriangleWithWavefrontEdge(wfEdge);
				if (tri == null) {
					System.err.println("Warning: Could not find incident triangle for edge " + wfEdge);
					continue; // Skip linking if no triangle is found
				}
				wfEdge.setIncidentTriangle(tri); // Link back if found this way
			}

			int edgeIndexInTri = -1;
			for (int i = 0; i < 3; ++i) {
				if (tri.getWavefront(i) == wfEdge) {
					edgeIndexInTri = i;
					break;
				}
			}
			if (edgeIndexInTri == -1) {
				throw new IllegalStateException("WavefrontEdge " + wfEdge + " not found in its incident triangle " + tri);
			}

			// Edge edgeIndexInTri is opposite vertex edgeIndexInTri
			// Vertices for this edge are (edgeIndexInTri+1)%3 and (edgeIndexInTri+2)%3
			WavefrontVertex triV1 = tri.getVertex((edgeIndexInTri + 1) % 3);
			WavefrontVertex triV2 = tri.getVertex((edgeIndexInTri + 2) % 3);

			// Assign vertices to edge based on triangle's vertex order
			// Vertex 0 of edge should correspond to triV1 (CCW order)
			// Vertex 1 of edge should correspond to triV2 (CW order)
			if (triV1 == v0 && triV2 == v1) {
				wfEdge.setVerticesAndUpdateAdj(v0, v1);
			} else if (triV1 == v1 && triV2 == v0) {
				wfEdge.setVerticesAndUpdateAdj(v1, v0); // Assign in order matching triangle
			} else {
				throw new IllegalStateException("Vertices mismatch between edge " + wfEdge + " and triangle " + tri);
			}

			// Link vertices back to this edge
			// Vertex 0 (wfEdge.getVertex(0)) should have this as its edge1 (CW)
			// Vertex 1 (wfEdge.getVertex(1)) should have this as its edge0 (CCW)
			// NOTE handled by AndUpdateAdj?
//			if (wfEdge.getVertex(0) != null) {
//				wfEdge.getVertex(0).setIncidentEdge(1, wfEdge); // Set as right/CW edge
//			}
//			if (wfEdge.getVertex(1) != null) {
//				wfEdge.getVertex(1).setIncidentEdge(0, wfEdge); // Set as left/CCW edge
//			}

		}

		// Final check for vertex incident edge consistency (optional, debug)
		// for (WavefrontVertex wv : this.vertices) {
		// if (!wv.isInfinite && (wv.getIncidentEdge(0) == null || wv.getIncidentEdge(1)
		// == null)) {
		// System.err.println("Warning: Vertex " + wv + " has incomplete incident
		// edges.");
		// }
		// }

	}

	/** Creates a new WavefrontVertex, adds it to the main list, and returns it. */
	public WavefrontVertex createAndAddVertex(Coordinate posZero, Coordinate posStart, double timeStart, WavefrontEdge edgeA, WavefrontEdge edgeB) {
		// Use the factory method logic from WavefrontVertex.make_vertex C++ equivalent
		// Need compute_intersection and potentially compute_velocity logic here or
		// accessible
		// For now, simplify: assume posZero is correctly calculated externally or is
		// posStart
		WavefrontVertex newVertex = new WavefrontVertex(posStart); // Simplified creation
		newVertex.setIncidentEdges(edgeA, edgeB); // Link to edges

		// Link edges back to vertex - IMPORTANT
		// NOTE setVertexRawAndUpdateAdj instead?
		if (edgeA != null) {
			edgeA.setVertexRaw(1, newVertex); // edgeA is edge0 (CCW), its vertex1 is the new vertex
		}
		if (edgeB != null) {
			edgeB.setVertexRaw(0, newVertex); // edgeB is edge1 (CW), its vertex0 is the new vertex
		}

		this.vertices.add(newVertex);
		this.newVerticesThisStep.add(newVertex); // Track for potential later processing if needed
		System.out.println("Created Vertex: " + newVertex);
		return newVertex;
	}

	/**
	 * Splits a wavefront edge due to a split event. Marks the original edge dead,
	 * creates two new edges, adds them to the main list, and returns them. The
	 * caller is responsible for setting the new vertices on the returned edges.
	 *
	 * @param edgeToSplit The original WavefrontEdge.
	 * @return A Pair containing the two new WavefrontEdges (edge corresponding to
	 *         original vertex0, edge corresponding to original vertex1).
	 */
	public Pair<WavefrontEdge, WavefrontEdge> splitWavefrontEdge(WavefrontEdge edgeToSplit) {
		if (edgeToSplit == null || edgeToSplit.isDead()) {
			throw new IllegalArgumentException("Cannot split null or dead edge.");
		}
		System.out.println("Splitting Edge: " + edgeToSplit);

		edgeToSplit.markDead();

		// Create new edges with same underlying segment and weight

		WavefrontEdge newEdgeA = new WavefrontEdge(new WavefrontSupportingLine(edgeToSplit.getSegment(), edgeToSplit.weight));
		WavefrontEdge newEdgeB = new WavefrontEdge(new WavefrontSupportingLine(edgeToSplit.getSegment(), edgeToSplit.weight));

		// Add to main list
		this.wavefrontEdges.add(newEdgeA);
		this.wavefrontEdges.add(newEdgeB);
		this.newEdgesThisStep.add(newEdgeA); // Track if needed
		this.newEdgesThisStep.add(newEdgeB);

		// Set initial vertex links (caller will set the split vertex)
		newEdgeA.setVerticesAndUpdateAdj(edgeToSplit.getVertex(0), null); // Keeps original vertex 0
		newEdgeB.setVerticesAndUpdateAdj(null, edgeToSplit.getVertex(1)); // Keeps original vertex 1

		// Update vertex linkage for the retained vertices
//		if (edgeToSplit.getVertex(0) != null && !edgeToSplit.getVertex(0).isInfinite) {
//			edgeToSplit.getVertex(0).setIncidentEdge(1, newEdgeA); // Original v0's edge1 (CW) is now newEdgeA
//		}
//		if (edgeToSplit.getVertex(1) != null && !edgeToSplit.getVertex(1).isInfinite) {
//			edgeToSplit.getVertex(1).setIncidentEdge(0, newEdgeB); // Original v1's edge0 (CCW) is now newEdgeB
//		}

		return Pair.of(newEdgeA, newEdgeB);
	}

	/**
	 * Performs an edge flip between two adjacent triangles. Assumes the edge is
	 * internal (not a constraint) and validity checks passed.
	 *
	 * @param triangle  The first triangle.
	 * @param edgeIndex The index (0, 1, 2) of the internal edge in the first
	 *                  triangle to flip.
	 */
	public void flipEdge(KineticTriangle triangle, int edgeIndex) {
		int i = edgeIndex % 3;
		if (triangle == null || triangle.isDead()) {
			throw new IllegalArgumentException("Cannot flip edge in null or dead triangle");
		}
		if (triangle.isConstrained(i)) {
			throw new IllegalArgumentException("Cannot flip a constrained edge (index " + i + " in triangle " + triangle.id + ")");
		}

		KineticTriangle neighbor = triangle.getNeighbor(i);
		if (neighbor == null || neighbor.isDead()) {
			throw new IllegalArgumentException("Cannot flip edge with null or dead neighbor (triangle " + triangle.id + " edge " + i + ")");
		}

		int neighborEdgeIndex = neighbor.indexOfNeighbor(triangle); // Index of edge in neighbor opposite triangle
		if (neighbor.isConstrained(neighborEdgeIndex)) {
			throw new IllegalStateException("Neighbor " + neighbor.id + " has constraint where triangle " + triangle.id + " expects internal edge."); // Should
																																						// not
																																						// happen
																																						// if
																																						// structure
																																						// is
																																						// consistent
		}

		System.out.println("Flipping Edge between KT" + triangle.id + "[" + i + "] and KT" + neighbor.id + "[" + neighborEdgeIndex + "]");

		// --- Port C++ do_raw_flip_inner logic ---

		// 1. Identify vertices
		WavefrontVertex v_tri = triangle.getVertex(i); // Vertex in tri opposite the edge
		WavefrontVertex v_neigh = neighbor.getVertex(neighborEdgeIndex); // Vertex in neighbor opposite the edge
		WavefrontVertex v_shared1 = triangle.getVertex((i + 1) % 3); // Shared vertex 1 (CCW in tri)
		WavefrontVertex v_shared2 = triangle.getVertex((i + 2) % 3); // Shared vertex 2 (CW in tri)

		// Check consistency
		if (v_shared1 != neighbor.getVertex((neighborEdgeIndex + 2) % 3) || // Should be CW in neighbor
				v_shared2 != neighbor.getVertex((neighborEdgeIndex + 1) % 3)) { // Should be CCW in neighbor
			throw new IllegalStateException("Shared vertex mismatch during flip preparation.");
		}

		// 2. Update vertex assignments in triangles
		triangle.setVertex((i + 2) % 3, v_neigh); // Replace v_shared2 with v_neigh in triangle
		neighbor.setVertex((neighborEdgeIndex + 2) % 3, v_tri); // Replace v_shared1 with v_tri in neighbor

		// 3. Update neighbor/wavefront links - Careful reassignment!
		// Let edgeN1 = neighbor's edge previously adjacent to v_shared1 (index
		// (neighborEdgeIndex + 1) % 3)
		// Let edgeN2 = neighbor's edge previously adjacent to v_shared2 (index
		// (neighborEdgeIndex + 2) % 3)
		// Let edgeT1 = triangle's edge previously adjacent to v_shared1 (index (i + 2)
		// % 3)
		// Let edgeT2 = triangle's edge previously adjacent to v_shared2 (index (i + 1)
		// % 3)

		KineticTriangle neighborOfN2 = neighbor.getNeighbor((neighborEdgeIndex + 2) % 3);
		WavefrontEdge wavefrontOfN2 = neighbor.getWavefront((neighborEdgeIndex + 2) % 3);
		KineticTriangle neighborOfT1 = triangle.getNeighbor((i + 1) % 3);
		WavefrontEdge wavefrontOfT1 = triangle.getWavefront((i + 1) % 3);

		// Update edge i in triangle (now borders neighborOfN2 or has wavefrontOfN2)
		triangle.setNeighbor(i, neighborOfN2);
		triangle.setWavefront(i, wavefrontOfN2);
		if (neighborOfN2 != null) {
			neighborOfN2.setNeighbor(neighborOfN2.indexOfNeighbor(neighbor), triangle); // Update backlink
		}
		if (wavefrontOfN2 != null) {
			wavefrontOfN2.setIncidentTriangle(triangle);
		}

		// Update edge neighborEdgeIndex in neighbor (now borders neighborOfT1 or has
		// wavefrontOfT1)
		neighbor.setNeighbor(neighborEdgeIndex, neighborOfT1);
		neighbor.setWavefront(neighborEdgeIndex, wavefrontOfT1);
		if (neighborOfT1 != null) {
			neighborOfT1.setNeighbor(neighborOfT1.indexOfNeighbor(triangle), neighbor); // Update backlink
		}
		if (wavefrontOfT1 != null) {
			wavefrontOfT1.setIncidentTriangle(neighbor);
		}

		// Update the edge between triangle and neighbor (indices (i+1)%3 in tri,
		// (neighborEdgeIndex+1)%3 in neighbor)
		triangle.setNeighbor((i + 1) % 3, neighbor);
		triangle.setWavefront((i + 1) % 3, null);
		neighbor.setNeighbor((neighborEdgeIndex + 1) % 3, triangle);
		neighbor.setWavefront((neighborEdgeIndex + 1) % 3, null);

		// 4. Update WavefrontVertex edge incidence if needed (though setVertex should
		// handle this?)
		// Re-linking vertex-to-edge is tricky here, might need dedicated logic
		// after the flip if setVertex doesn't fully capture it.
		// For now, assume setVertex updates are sufficient, but verify.
		// Crucially, v_tri and v_neigh now have different adjacent edges.

		System.out.println(
				" Flipped: KT" + triangle.id + " now vertices [" + triangle.getVertex(0) + "," + triangle.getVertex(1) + "," + triangle.getVertex(2) + "]");
		System.out.println(
				" Flipped: KT" + neighbor.id + " now vertices [" + neighbor.getVertex(0) + "," + neighbor.getVertex(1) + "," + neighbor.getVertex(2) + "]");

	}

	/** Helper to find incident triangles for a vertex (needed by handlers). */
	public List<KineticTriangle> findIncidentTriangles(WavefrontVertex vertex) {
		List<KineticTriangle> incident = new ArrayList<>();
		for (KineticTriangle kt : triangles) {
			if (kt.isDead()) {
				continue;
			}
			for (int i = 0; i < 3; i++) {
				if (kt.getVertex(i) == vertex) { // Object identity
					incident.add(kt);
					break;
				}
			}
		}
		return incident;
	}

	/** Clear lists tracking elements created during the current step. */
	public void clearNewElementsList() {
		newVerticesThisStep.clear();
		newEdgesThisStep.clear();
	}

	/** Helper to create or retrieve a WavefrontEdge, ensuring uniqueness. */
	private WavefrontEdge createOrGetWavefrontEdge(CanonicalSegment segment) {
		return edgeMap.computeIfAbsent(segment, k -> {
			// TODO: Get weight from actual input source if available
			double weight = 1.0;
			// Use the constructor that creates the WavefrontSupportingLine
			WavefrontEdge newEdge = new WavefrontEdge(k.getSegment(), weight);
			this.wavefrontEdges.add(newEdge);
			return newEdge;
		});
	}

	/** Final linking step after all triangles/edges are processed. */
	private void linkWavefrontEdgesAndVertices() {
		for (WavefrontEdge wfEdge : this.wavefrontEdges) {
			LineSegment seg = wfEdge.getSegment();
			Coordinate c0 = seg.p0; // Use original segment coords for lookup
			Coordinate c1 = seg.p1;
			WavefrontVertex mapV0 = vertexMap.get(new CanonicalSegment(c0, c0).getP0()); // Lookup canonical coord
			WavefrontVertex mapV1 = vertexMap.get(new CanonicalSegment(c1, c1).getP0()); // Lookup canonical coord

			if (mapV0 == null || mapV1 == null) {
				throw new IllegalStateException("Vertex not found in vertexMap for WavefrontEdge segment: " + seg + " c0Key:"
						+ new CanonicalSegment(c0, c0).getP0() + " c1Key:" + new CanonicalSegment(c1, c1).getP0());
			}

			// Determine orientation relative to the incident triangle(s)
			KineticTriangle tri = wfEdge.getIncidentTriangle(); // Get the one set during Pass 2

			// If triangle is null, try finding the other side (less common)
			if (tri == null) {
				tri = findOtherTriangleWithWavefrontEdge(wfEdge);
				if (tri != null) {
					wfEdge.setIncidentTriangle(tri); // Set if found
				} else {
					System.err.println("Warning: Could not find any incident triangle for edge " + wfEdge.id + " seg=" + seg);
					continue; // Skip linking if no triangle found
				}
			}

			int edgeIndexInTri = tri.indexOfWavefront(wfEdge);
			if (edgeIndexInTri == -1) {
				// This might happen if findOtherTriangleWithWavefrontEdge was called and set
				// the wrong tri?
				// Or if the initial setIncidentTriangle failed.
				System.err.println("Error: WavefrontEdge " + wfEdge.id + " not found in its supposedly incident triangle " + tri.id);
				// Attempt to find the correct triangle again?
				boolean found = false;
				for (KineticTriangle ktCheck : triangles) {
					for (int i = 0; i < 3; ++i) {
						if (ktCheck.getWavefront(i) == wfEdge) {
							tri = ktCheck;
							edgeIndexInTri = i;
							wfEdge.setIncidentTriangle(tri);
							found = true;
							System.err.println("Found edge in triangle " + tri.id + " edge index " + i);
							break;
						}
					}
					if (found) {
						break;
					}
				}
				if (!found) {
					System.err.println("FATAL: Could not resolve incident triangle for edge " + wfEdge.id);
					continue; // Skip this edge if unresolvable
				}
			}

			// Edge edgeIndexInTri is opposite vertex edgeIndexInTri
			// Vertices forming this edge in the triangle are (edgeIndexInTri+1)%3 and
			// (edgeIndexInTri+2)%3
			// These indices give us the *expected* WavefrontVertex objects from the
			// triangle's perspective
			WavefrontVertex triV_CCW = tri.getVertex((edgeIndexInTri + 1) % 3); // This vertex is CCW of edge
			WavefrontVertex triV_CW = tri.getVertex((edgeIndexInTri + 2) % 3); // This vertex is CW of edge

			// Assign vertices to edge in the order (CCW_Vertex, CW_Vertex) which
			// corresponds to (vertex0, vertex1)
			// This matches the C++ convention where edge->vertex(0) is the "left" vertex
			// seen from the triangle.
			if ((mapV0 == triV_CCW && mapV1 == triV_CW)) {
				wfEdge.setVerticesAndUpdateAdj(mapV0, mapV1);
			} else if ((mapV1 == triV_CCW && mapV0 == triV_CW)) {
				wfEdge.setVerticesAndUpdateAdj(mapV1, mapV0); // Assign in (CCW, CW) order
			} else {
				// This can happen if the vertexMap lookup failed or triangle vertices are wrong
				System.err.println("TRIANGLE: " + tri);
				System.err.println("EDGE: " + wfEdge);
				System.err.println("EDGE SEG: " + seg);
				System.err.println("MAP V0: " + mapV0 + " (" + mapV0.initialPosition + ")");
				System.err.println("MAP V1: " + mapV1 + " (" + mapV1.initialPosition + ")");
				System.err.println("TRI V_CCW (idx " + ((edgeIndexInTri + 1) % 3) + "): " + triV_CCW + " ("
						+ (triV_CCW != null ? triV_CCW.initialPosition : "null") + ")");
				System.err.println(
						"TRI V_CW  (idx " + ((edgeIndexInTri + 2) % 3) + "): " + triV_CW + " (" + (triV_CW != null ? triV_CW.initialPosition : "null") + ")");

				throw new IllegalStateException("Vertex mismatch linking edge " + wfEdge.id + " to triangle " + tri.id);
			}

			// Link vertices back to this edge
			// Vertex 0 (CCW) should have this as its edge1 (CW edge *relative to the
			// vertex*)
			// Vertex 1 (CW) should have this as its edge0 (CCW edge *relative to the
			// vertex*)
//			WavefrontVertex v0 = wfEdge.getVertex(0);
//			WavefrontVertex v1 = wfEdge.getVertex(1);
//
//			if (v0 != null && !v0.isInfinite) {
//				if (v0.getIncidentEdge(1) != null && v0.getIncidentEdge(1) != wfEdge) {
//					System.err
//							.println("Warning: Overwriting incident edge 1 for vertex " + v0.id + ". Old: " + v0.getIncidentEdge(1).id + " New: " + wfEdge.id);
//				}
//				v0.setIncidentEdge(1, wfEdge);
//			}
//			if (v1 != null && !v1.isInfinite) {
//				if (v1.getIncidentEdge(0) != null && v1.getIncidentEdge(0) != wfEdge) {
//					System.err
//							.println("Warning: Overwriting incident edge 0 for vertex " + v1.id + ". Old: " + v1.getIncidentEdge(0).id + " New: " + wfEdge.id);
//				}
//				v1.setIncidentEdge(0, wfEdge);
//			}
		}

		// Final check: Ensure all non-infinite vertices have both incident edges set
		for (WavefrontVertex v : this.vertices) {
			if (!v.isInfinite) {
				if (v.getIncidentEdge(0) == null) {
					System.err.println("Warning: Vertex " + v.id + " (" + v.initialPosition + ") missing incident edge 0 (CCW) after linking.");
				}
				if (v.getIncidentEdge(1) == null) {
					System.err.println("Warning: Vertex " + v.id + " (" + v.initialPosition + ") missing incident edge 1 (CW) after linking.");
				}
			}
		}
	}

	// Helper needed for linking pass 3
	public KineticTriangle findOtherTriangleWithWavefrontEdge(WavefrontEdge edge) {
		KineticTriangle first = edge.getIncidentTriangle();
		for (KineticTriangle kt : triangles) {
			if (kt == first) {
				continue; // Skip the one we already know
			}
			for (int i = 0; i < 3; i++) {
				if (kt.getWavefront(i) == edge) {
					return kt;
				}
			}
		}
		return null; // Only one triangle borders this edge (or none were found)
	}

	private KineticTriangle findTriangleWithWavefrontEdge(WavefrontEdge edge) {
		for (KineticTriangle kt : triangles) {
			for (int i = 0; i < 3; i++) {
				if (kt.getWavefront(i) == edge) {
					return kt;
				}
			}
		}
		return null;
	}

	private void populateConstraintSegments(Polygon polygon) {
		extractSegments(polygon.getExteriorRing(), constraintSegments);
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			extractSegments(polygon.getInteriorRingN(i), constraintSegments);
		}
	}

	private void extractSegments(LineString ring, Set<CanonicalSegment> segmentSet) {
		Coordinate[] coords = ring.getCoordinates();
		for (int i = 0; i < coords.length - 1; i++) {
			// Skip zero-length segments
			if (!coords[i].equals2D(coords[i + 1])) {
				segmentSet.add(new CanonicalSegment(coords[i], coords[i + 1]));
			}
		}
	}

	// Helper to check if a triangle is roughly inside the polygon
	// Uses centroid test, which is usually sufficient for CDT results
	private boolean isTriangleInside(QuadEdgeTriangle triangle, Polygon polygon) {
		// Calculate centroid
		Coordinate p0 = triangle.getCoordinate(0);
		Coordinate p1 = triangle.getCoordinate(1);
		Coordinate p2 = triangle.getCoordinate(2);
		

		// Check for degenerate triangles (collinear vertices) before centroid
		// calculation

		if (p0 == null || p1 == null || p2 == null || (Orientation.index(p0, p1, p2) == Orientation.COLLINEAR)) {
			return false; // Or handle as needed, CDT usually avoids this
		}

		Coordinate centroid = new Coordinate((p0.x + p1.x + p2.x) / 3.0, (p0.y + p1.y + p2.y) / 3.0);
		Point centroidPoint = polygon.getFactory().createPoint(centroid);
		return polygon.contains(centroidPoint); // Use 'contains' which is robust
	}

	// Helper to find neighbor using QuadEdgeTriangle structure (simpler than raw
	// QuadEdge)
	private QuadEdgeTriangle findNeighbor(QuadEdgeTriangle triangle, int edgeIndex, List<QuadEdgeTriangle> allTriangles) {
		Coordinate p0 = triangle.getCoordinate(edgeIndex);
		Coordinate p1 = triangle.getCoordinate((edgeIndex + 1) % 3);
		LineSegment edge = new LineSegment(p0, p1);

		for (QuadEdgeTriangle potentialNeighbor : allTriangles) {
			if (potentialNeighbor == triangle) {
				continue; // Skip self
			}

			for (int j = 0; j < 3; j++) {
				Coordinate n0 = potentialNeighbor.getCoordinate(j);
				Coordinate n1 = potentialNeighbor.getCoordinate((j + 1) % 3);
				LineSegment neighborEdge = new LineSegment(n0, n1);

				// Check if edges are the same (ignoring direction)
				if ((edge.p0.equals2D(neighborEdge.p0) && edge.p1.equals2D(neighborEdge.p1))
						|| (edge.p0.equals2D(neighborEdge.p1) && edge.p1.equals2D(neighborEdge.p0))) {
					return potentialNeighbor;
				}
			}
		}
		return null; // No neighbor found
	}

	public Map<CanonicalSegment, WavefrontEdge> getEdgeMap() {
		return edgeMap;
	}

}