package com.github.micycle.surferj2.kinetics;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.triangulate.DelaunayTriangulationBuilder;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeSubdivision;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeTriangle;

import java.util.*;

/**
 * The static structure of the triangulation at time t=0. The algorithm now
 * needs to determine how this structure will evolve over time. This evolution
 * is driven by events, which signify topological changes in the wavefront.
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
			qetToKt.put(jtsTriangle, kt);

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
		for (Map.Entry<QuadEdgeTriangle, KineticTriangle> entry : qetToKt.entrySet()) {
			QuadEdgeTriangle jtsTriangle = entry.getKey();
			KineticTriangle currentKt = entry.getValue();

			for (int i = 0; i < 3; i++) {
				Coordinate p0 = jtsTriangle.getCoordinate(i);
				Coordinate p1 = jtsTriangle.getCoordinate((i + 1) % 3);
				CanonicalSegment edgeSegment = new CanonicalSegment(p0, p1);

				boolean isConstraint = constraintSegments.contains(edgeSegment);

				if (isConstraint) {
					// Find or create the WavefrontEdge
					WavefrontEdge wfEdge = edgeMap.computeIfAbsent(edgeSegment, k -> {
						// TODO: Get weight from input if available, default 1.0
						WavefrontEdge newWfEdge = new WavefrontEdge(k, 1.0);
						this.wavefrontEdges.add(newWfEdge);
						return newWfEdge;
					});
					currentKt.setWavefront((i + 2) % 3, wfEdge); // Edge i is opposite vertex (i+2)%3
				} else {
					// Find neighbor triangle
					QuadEdgeTriangle neighborQet = findNeighbor(jtsTriangle, i, jtsTriangles);
					if (neighborQet != null) {
						KineticTriangle neighborKt = qetToKt.get(neighborQet);
						if (neighborKt != null) { // Check if neighbor is inside polygon
							currentKt.setNeighbor((i + 2) % 3, neighborKt);
						} else {
							// Neighbor exists in triangulation but is outside polygon area
							// Treat this edge like a constraint from the perspective of currentKt
							// This can happen along hole boundaries or the outer boundary
							// Ensure this segment IS actually a constraint segment in these cases
							if (!constraintSegments.contains(edgeSegment)) {
								// This case might indicate an issue or a very specific geometry
								// For simplicity now, we might ignore it or treat as constraint
								System.err.println("Warning: Internal edge " + edgeSegment + " borders an external triangle for " + currentKt);
								// Optionally create a WF edge here if needed by algorithm logic
							} else {
								// It's a constraint segment bordering outside, set WF edge
								WavefrontEdge wfEdge = edgeMap.computeIfAbsent(edgeSegment, k -> {
									WavefrontEdge newWfEdge = new WavefrontEdge(k, 1.0);
									this.wavefrontEdges.add(newWfEdge);
									return newWfEdge;
								});
								currentKt.setWavefront((i + 2) % 3, wfEdge);
							}
						}
					} else {
						// No neighbor found via QuadEdgeTriangle logic, potentially boundary
						if (!constraintSegments.contains(edgeSegment)) {
							System.err.println("Warning: Non-constraint edge " + edgeSegment + " has no neighbor for " + currentKt);
						} else {
							// Should have been caught by isConstraint check, but safety check
							WavefrontEdge wfEdge = edgeMap.computeIfAbsent(edgeSegment, k -> {
								WavefrontEdge newWfEdge = new WavefrontEdge(k, 1.0);
								this.wavefrontEdges.add(newWfEdge);
								return newWfEdge;
							});
							currentKt.setWavefront((i + 2) % 3, wfEdge);
						}
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
				wfEdge.setVertices(v0, v1);
			} else if (triV1 == v1 && triV2 == v0) {
				wfEdge.setVertices(v1, v0); // Assign in order matching triangle
			} else {
				throw new IllegalStateException("Vertices mismatch between edge " + wfEdge + " and triangle " + tri);
			}

			// Link vertices back to this edge
			// Vertex 0 (wfEdge.getVertex(0)) should have this as its edge1 (CW)
			// Vertex 1 (wfEdge.getVertex(1)) should have this as its edge0 (CCW)
			if (wfEdge.getVertex(0) != null)
				wfEdge.getVertex(0).setIncidentEdge(1, wfEdge); // Set as right/CW edge
			if (wfEdge.getVertex(1) != null)
				wfEdge.getVertex(1).setIncidentEdge(0, wfEdge); // Set as left/CCW edge

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
		if (p0 == null || p1 == null || p2 == null)
			return false; // Should not happen for valid triangles

		// Check for degenerate triangles (collinear vertices) before centroid
		// calculation

		if (Orientation.index(p0, p1, p2) == Orientation.COLLINEAR) {
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
			if (potentialNeighbor == triangle)
				continue; // Skip self

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

}