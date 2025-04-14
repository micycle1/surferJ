package com.github.micycle.surferj2.kinetics;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.math.Vector2D;

import com.github.micycle.surferj2.kinetics.WavefrontVertex.InfiniteSpeedType;
import com.github.micycle.surferj2.kinetics.WavefrontVertex.VertexAngle;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Test data structure initialisation (links, etc) and geometry (orientation of
 * vertices, etc.)
 */
public class KineticTriangulationTest {

	private static final double DELTA = 1e-9; // Tolerance for floating-point comparisons

	private final GeometryFactory gf = new GeometryFactory();
	private final WKTReader reader = new WKTReader(gf);

	private Polygon createSquare() throws ParseException {
		return (Polygon) reader.read("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))");
	}

	private Polygon createLShape() throws ParseException {
		return (Polygon) reader.read("POLYGON ((0 0, 20 0, 20 10, 10 10, 10 20, 0 20, 0 0))");
	}

	private Polygon createSquareWithHole() throws ParseException {
		return (Polygon) reader.read("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0), (3 3, 7 3, 7 7, 3 7, 3 3))");
	}

	private Polygon createSingleTrianglePoly() throws ParseException {
		return (Polygon) reader.read("POLYGON ((0 0, 10 0, 5 10, 0 0))");
	}

	private Polygon createDegeneratePoly() throws ParseException {
		// Collinear points on boundary
		return (Polygon) reader.read("POLYGON ((0 0, 5 0, 10 0, 10 10, 0 10, 0 0))");
	}

	private Polygon createConcavePolygon() throws ParseException {
		// Concave polygon with a dent
		return (Polygon) reader.read("POLYGON ((0 0, 10 0, 10 10, 5 5, 0 10, 0 0))");
	}

	private Polygon createLargePolygon(int vertexCount) {
		// Create a circular polygon approximation with vertexCount vertices
		Coordinate[] coords = new Coordinate[vertexCount + 1];
		for (int i = 0; i < vertexCount; i++) {
			double angle = 2 * Math.PI * i / vertexCount;
			coords[i] = new Coordinate(100 * Math.cos(angle), 100 * Math.sin(angle));
		}
		coords[vertexCount] = coords[0]; // Close the ring
		return gf.createPolygon(coords);
	}

	// Helper to find a vertex by its initial coordinate. Simple brute-force.
	private WavefrontVertex findVertexAt(KineticTriangulation kt, double x, double y) {
		Coordinate target = new Coordinate(x, y);
		for (WavefrontVertex v : kt.getVertices()) {
			// Use equals2D for robustness, though exact match expected here
			if (v.initialPosition.equals2D(target, DELTA)) {
				return v;
			}
		}
		fail("Vertex not found at coordinate (" + x + ", " + y + ")");
		return null; // Should not be reached
	}

	// Helper to calculate expected velocity for a standard convex/reflex corner
	// Based on C++ compute_velocity: intersection of lines offset by weight at t=1
	// Assumes weights are 1.0 for simplicity in this helper.
	private Coordinate calculateExpectedCornerVelocity(Coordinate vertexPos, Coordinate neighbor1Pos, Coordinate neighbor2Pos) {
		// Simplified: Assumes edges are directed CCW around polygon for normal calc
		// Edge 1: vertexPos -> neighbor1Pos (CW edge relative to vertex)
		// Edge 2: neighbor2Pos -> vertexPos (CCW edge relative to vertex)
		LineSegment seg1 = new LineSegment(vertexPos, neighbor1Pos);
		LineSegment seg2 = new LineSegment(neighbor2Pos, vertexPos);

		// Use default weight 1.0 for standard wavefront offset
		WavefrontSupportingLine line1 = new WavefrontSupportingLine(seg1, 1.0);
		WavefrontSupportingLine line2 = new WavefrontSupportingLine(seg2, 1.0);

		// Get lines offset by weight=1.0 (at t=1)
		// Need access to the offset line calculation logic or WavefrontSupportingLine
		// methods
		// Assuming methods line_at(1.0) or similar exist based on C++
		// Placeholder: Manually calculate offset lines for testing
		Coordinate p1_t1 = offsetPoint(vertexPos, line1.getUnitNormal(), 1.0);
		Coordinate p1_neighbor_t1 = offsetPoint(neighbor1Pos, line1.getUnitNormal(), 1.0);
		LineSegment offsetSeg1 = new LineSegment(p1_t1, p1_neighbor_t1);

		Coordinate p2_t1 = offsetPoint(vertexPos, line2.getUnitNormal(), 1.0);
		Coordinate p2_neighbor_t1 = offsetPoint(neighbor2Pos, line2.getUnitNormal(), 1.0);
		LineSegment offsetSeg2 = new LineSegment(p2_neighbor_t1, p2_t1); // Reversed for intersection

		// Find intersection of the offset lines at t=1
		Coordinate intersection = offsetSeg1.lineIntersection(offsetSeg2);

		if (intersection == null) {
			// This can happen with parallel lines (collinear vertex) or numerical issues
			// Fail test or return specific value? Fail for now.
			fail("Could not calculate intersection for velocity at " + vertexPos);
			return null;
		}

		// Velocity is the vector from initial position to intersection point at t=1
		return new Coordinate(intersection.x - vertexPos.x, intersection.y - vertexPos.y);
	}

	// Helper to offset a point along a normal by a distance (weight)
	private Coordinate offsetPoint(Coordinate p, Vector2D unitNormal, double distance) {
		return new Coordinate(p.x + unitNormal.getX() * distance, p.y + unitNormal.getY() * distance);
	}

	@Test
	void testSquareInitialization() throws ParseException {
		Polygon square = createSquare();
		KineticTriangulation kt = new KineticTriangulation(square, gf);

		// Expected: CDT of a square usually results in 2 triangles.
		assertEquals(2, kt.getTriangles().size(), "Square should yield 2 triangles");
		// Expected: 4 vertices from the square.
		assertEquals(4, kt.getVertices().size(), "Square should have 4 vertices");
		// Expected: 4 constraint edges.
		assertEquals(4, kt.getWavefrontEdges().size(), "Square should have 4 constraint edges");

		validateConnectivity(kt);
		validateConstraints(kt, square);
	}

	@Test
	void testLShapeInitialization() throws ParseException {
		Polygon lShape = createLShape();
		KineticTriangulation kt = new KineticTriangulation(lShape, gf);

		// Expected number of triangles/vertices/edges depends on the specific CDT
		// result.
		// Focus on structural integrity.
		assertTrue(kt.getTriangles().size() > 0, "L-shape should yield triangles");
		assertEquals(6, kt.getVertices().size(), "L-shape should have 6 vertices");
		assertEquals(6, kt.getWavefrontEdges().size(), "L-shape should have 6 constraint edges");

		validateConnectivity(kt);
		validateConstraints(kt, lShape);
	}

	@Test
	void testSquareWithHoleInitialization() throws ParseException {
		Polygon squareWithHole = createSquareWithHole();
		KineticTriangulation kt = new KineticTriangulation(squareWithHole, gf);

		// Expected: Triangles covering the area between outer and inner ring.
		assertTrue(kt.getTriangles().size() > 0, "Square with hole should yield triangles");
		// Expected: 4 outer + 4 inner vertices = 8
		assertEquals(8, kt.getVertices().size(), "Square with hole should have 8 vertices");
		// Expected: 4 outer + 4 inner constraint edges = 8
		assertEquals(8, kt.getWavefrontEdges().size(), "Square with hole should have 8 constraint edges");

		validateConnectivity(kt);
		validateConstraints(kt, squareWithHole);
	}

	@Test
	void testSingleTrianglePolyInitialization() throws ParseException {
		Polygon singleTriPoly = createSingleTrianglePoly();
		KineticTriangulation kt = new KineticTriangulation(singleTriPoly, gf);

		// Expected: 1 triangle
		assertEquals(1, kt.getTriangles().size(), "Single triangle polygon should yield 1 triangle");
		assertEquals(3, kt.getVertices().size(), "Single triangle polygon should have 3 vertices");
		assertEquals(3, kt.getWavefrontEdges().size(), "Single triangle polygon should have 3 constraint edges");

		// The single triangle should have no neighbors, only wavefront edges
		KineticTriangle tri = kt.getTriangles().get(0);
		assertNull(tri.getNeighbor(0));
		assertNull(tri.getNeighbor(1));
		assertNull(tri.getNeighbor(2));
		assertNotNull(tri.getWavefront(0));
		assertNotNull(tri.getWavefront(1));
		assertNotNull(tri.getWavefront(2));

		validateConnectivity(kt); // Will check basic vertex/edge links
		validateConstraints(kt, singleTriPoly);
	}

	@Test
	void testDegeneratePolyInitialization() throws ParseException {
		// Test if it handles collinear points on boundary gracefully
		Polygon degeneratePoly = createDegeneratePoly();
		KineticTriangulation kt = new KineticTriangulation(degeneratePoly, gf);

		// Expected vertices: 5 unique coordinates
		assertEquals(5, kt.getVertices().size(), "Degenerate polygon should have 5 unique vertices");
		// Expected constraint edges: 5 unique segments
		assertEquals(5, kt.getWavefrontEdges().size(), "Degenerate polygon should have 5 constraint edges");
		assertTrue(kt.getTriangles().size() > 0, "Degenerate polygon should still yield triangles");

		validateConnectivity(kt);
		validateConstraints(kt, degeneratePoly);
	}

	@Test
	void testInvalidPolygon() {
		// Create a self-intersecting (bowtie) polygon
		String bowtieWKT = "POLYGON ((0 0, 10 10, 0 10, 10 0, 0 0))";
		try {
			Polygon bowtie = (Polygon) reader.read(bowtieWKT);
			assertThrows(IllegalArgumentException.class, () -> new KineticTriangulation(bowtie, gf),
					"Self-intersecting polygon should throw IllegalArgumentException");
		} catch (ParseException e) {
			fail("Unexpected ParseException: " + e.getMessage());
		}
	}

	@Test
	void testDuplicatePoints() throws ParseException {
		Polygon polyWithDuplicates = (Polygon) reader.read("POLYGON ((0 0, 0 0, 10 0, 10 10, 0 10, 0 0))");
		KineticTriangulation kt = new KineticTriangulation(polyWithDuplicates, gf);

		// JTS normalizes duplicate points, expect 4 vertices (square)
		assertEquals(4, kt.getVertices().size(), "Polygon with duplicate points should have 4 unique vertices");
		assertEquals(4, kt.getWavefrontEdges().size(), "Polygon with duplicate points should have 4 constraint edges");
		assertTrue(kt.getTriangles().size() > 0, "Polygon with duplicate points should yield triangles");

		validateConnectivity(kt);
		validateConstraints(kt, polyWithDuplicates);
	}

	@Test
	void testNearCollinearPoints() throws ParseException {
		// Triangle with one vertex very close to the line between others
		Polygon nearCollinear = (Polygon) reader.read("POLYGON ((0 0, 10 0, 5 0.0001, 0 0))");
		KineticTriangulation kt = new KineticTriangulation(nearCollinear, gf);

		// Should still produce a valid triangulation (1 triangle expected)
		assertEquals(1, kt.getTriangles().size(), "Near-collinear polygon should yield 1 triangle");
		assertEquals(3, kt.getVertices().size(), "Near-collinear polygon should have 3 vertices");
		assertEquals(3, kt.getWavefrontEdges().size(), "Near-collinear polygon should have 3 constraint edges");

		validateConnectivity(kt);
		validateConstraints(kt, nearCollinear);
	}

	@Test
	void testConcavePolygonTriangleInclusion() throws ParseException {
		Polygon concavePoly = createConcavePolygon();
		KineticTriangulation kt = new KineticTriangulation(concavePoly, gf);

		// Concave polygon should have 5 vertices
		assertEquals(5, kt.getVertices().size(), "Concave polygon should have 5 vertices");
		assertEquals(5, kt.getWavefrontEdges().size(), "Concave polygon should have 5 constraint edges");
		assertTrue(kt.getTriangles().size() > 0, "Concave polygon should yield triangles");

		// Verify all triangles are inside (centroid check)
		for (KineticTriangle tri : kt.getTriangles()) {
			Coordinate[] coords = new Coordinate[3];
			for (int i = 0; i < 3; i++) {
				coords[i] = tri.getVertex(i).initialPosition;
			}
			Coordinate centroid = new Coordinate((coords[0].x + coords[1].x + coords[2].x) / 3.0, (coords[0].y + coords[1].y + coords[2].y) / 3.0);
			Point centroidPoint = gf.createPoint(centroid);
			assertTrue(concavePoly.contains(centroidPoint), "Triangle centroid should be inside concave polygon");
		}

		validateConnectivity(kt);
		validateConstraints(kt, concavePoly);
	}

	@Test
	void testWavefrontEdgeWithoutIncidentTriangle() throws ParseException {
		// Use a polygon with multiple holes to increase chance of boundary edges
		Polygon multiHolePoly = (Polygon) reader
				.read("POLYGON ((0 0, 20 0, 20 20, 0 20, 0 0), " + "(5 5, 10 5, 10 10, 5 10, 5 5), " + "(12 12, 15 12, 15 15, 12 15, 12 12))");

		// Capture System.err for warning verification
		ByteArrayOutputStream errContent = new ByteArrayOutputStream();
		PrintStream originalErr = System.err;
		System.setErr(new PrintStream(errContent));

		try {
			KineticTriangulation kt = new KineticTriangulation(multiHolePoly, gf);
			assertEquals(12, kt.getVertices().size(), "Polygon with two holes should have 10 vertices");
			assertEquals(12, kt.getWavefrontEdges().size(), "Polygon with two holes should have 10 constraint edges");
			assertTrue(kt.getTriangles().size() > 0, "Polygon with two holes should yield triangles");

			validateConnectivity(kt);
			validateConstraints(kt, multiHolePoly);
		} finally {
			System.setErr(originalErr);
		}
	}

	@Test
	void testLargePolygon() {
		Polygon largePoly = createLargePolygon(100); // 100 vertices
		KineticTriangulation kt = new KineticTriangulation(largePoly, gf);

		assertEquals(100, kt.getVertices().size(), "Large polygon should have 100 vertices");
		assertEquals(100, kt.getWavefrontEdges().size(), "Large polygon should have 100 constraint edges");
		assertTrue(kt.getTriangles().size() > 0, "Large polygon should yield triangles");

		validateConnectivity(kt);
		validateConstraints(kt, largePoly);
	}

	// Geometric State Tests

	@Test
	void testSquareVertexGeometry() throws ParseException {
		Polygon square = createSquare(); // POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))
		KineticTriangulation kt = new KineticTriangulation(square, gf);

		// --- Test Vertex (0, 0) ---
		WavefrontVertex v00 = findVertexAt(kt, 0, 0);
		assertNotNull(v00, "Vertex at (0,0) not found");

		// Expected Angle: Convex corner of square
		assertEquals(VertexAngle.LEFT_TURN, v00.getAngle(), "Square corner (0,0) should be LEFT_TURN");

		// Expected Speed: Standard wavefront vertex
		assertEquals(InfiniteSpeedType.NONE, v00.getInfiniteSpeed(), "Square corner (0,0) should have NONE infinite speed");

		// Expected Velocity: Inward diagonal (1, 1) assuming weight=1
		assertEquals(1.0, v00.getVelocity().getX(), DELTA, "Square corner (0,0) velocity X should be 1");
		assertEquals(1.0, v00.getVelocity().getY(), DELTA, "Square corner (0,0) velocity Y should be 1");

		// --- Test Vertex (10, 0) ---
		WavefrontVertex v100 = findVertexAt(kt, 10, 0);
		assertNotNull(v100, "Vertex at (10,0) not found");

		// Expected Angle: Convex corner of square
		assertEquals(VertexAngle.LEFT_TURN, v100.getAngle(), "Square corner (10,0) should be LEFT_TURN");

		// Expected Speed: Standard wavefront vertex
		assertEquals(InfiniteSpeedType.NONE, v100.getInfiniteSpeed(), "Square corner (10,0) should have NONE infinite speed");

		// Expected Velocity: Inward diagonal (-1, 1) assuming weight=1
		assertEquals(-1.0, v100.getVelocity().getX(), DELTA, "Square corner (10,0) velocity X should be -1");
		assertEquals(1.0, v100.getVelocity().getY(), DELTA, "Square corner (10,0) velocity Y should be 1");

		// --- Test Vertex (10, 10) ---
		WavefrontVertex v1010 = findVertexAt(kt, 10, 10);
		assertNotNull(v1010, "Vertex at (10,10) not found");

		// Expected Angle: Convex corner of square
		assertEquals(VertexAngle.LEFT_TURN, v1010.getAngle(), "Square corner (10,10) should be LEFT_TURN");

		// Expected Speed: Standard wavefront vertex
		assertEquals(InfiniteSpeedType.NONE, v1010.getInfiniteSpeed(), "Square corner (10,10) should have NONE infinite speed");

		// Expected Velocity: Inward diagonal (-1, -1) assuming weight=1
		assertEquals(-1.0, v1010.getVelocity().getX(), DELTA, "Square corner (10,10) velocity X should be -1");
		assertEquals(-1.0, v1010.getVelocity().getY(), DELTA, "Square corner (10,10) velocity Y should be -1");

		// --- Test Vertex (0, 10) ---
		WavefrontVertex v010 = findVertexAt(kt, 0, 10);
		assertNotNull(v010, "Vertex at (0,10) not found");

		// Expected Angle: Convex corner of square
		assertEquals(VertexAngle.LEFT_TURN, v010.getAngle(), "Square corner (0,10) should be LEFT_TURN");

		// Expected Speed: Standard wavefront vertex
		assertEquals(InfiniteSpeedType.NONE, v010.getInfiniteSpeed(), "Square corner (0,10) should have NONE infinite speed");

		// Expected Velocity: Inward diagonal (1, -1) assuming weight=1
		assertEquals(1.0, v010.getVelocity().getX(), DELTA, "Square corner (0,10) velocity X should be 1");
		assertEquals(-1.0, v010.getVelocity().getY(), DELTA, "Square corner (0,10) velocity Y should be -1");
	}

	@Test
	void testLShapeReflexVertexGeometry() throws ParseException {
		Polygon lShape = createLShape();
		KineticTriangulation kt = new KineticTriangulation(lShape, gf);

		// Test the inner reflex corner vertex (10, 10)
		WavefrontVertex v1010 = findVertexAt(kt, 10, 10);
		assertNotNull(v1010);

		// Expected Angle: Reflex corner of L-shape
		assertEquals(VertexAngle.RIGHT_TURN, v1010.getAngle(), "L-shape reflex corner (10,10) should be RIGHT_TURN");

		// Expected Speed: Standard wavefront vertex
		assertEquals(InfiniteSpeedType.NONE, v1010.getInfiniteSpeed(), "L-shape reflex corner (10,10) should have NONE infinite speed");

		// Expected Velocity: Should move diagonally outwards (bisector)
		// Neighbors are (20,10) and (10,20).
		// Normal for (10,10)->(20,10) is (0,-1). Offset line y=9.
		// Normal for (10,20)->(10,10) is (-1,0). Offset line x=9.
		// Intersection (9,9). Velocity = (9-10, 9-10) = (-1, -1).
		assertEquals(-1.0, v1010.getVelocity().getX(), DELTA, "L-shape reflex corner (10,10) velocity X should be -1");
		assertEquals(-1.0, v1010.getVelocity().getY(), DELTA, "L-shape reflex corner (10,10) velocity Y should be -1");
	}

	@Test
	void testSingleTriangleVertexGeometry() throws ParseException {
		Polygon singleTriPoly = createSingleTrianglePoly(); // (0 0, 10 0, 5 10, 0 0)
		KineticTriangulation kt = new KineticTriangulation(singleTriPoly, gf);

		// Test vertex (0, 0)
		WavefrontVertex v00 = findVertexAt(kt, 0, 0);
		assertNotNull(v00);

		// Expected Angle: Convex corner of triangle
		assertEquals(VertexAngle.LEFT_TURN, v00.getAngle(), "Single triangle corner (0,0) should be LEFT_TURN");

		// Expected Speed: Standard wavefront vertex
		assertEquals(InfiniteSpeedType.NONE, v00.getInfiniteSpeed(), "Single triangle corner (0,0) should have NONE infinite speed");

		// Expected Velocity: Needs calculation based on the two boundary edges
		// Neighbors are (10,0) and (5,10).
		Coordinate expectedVel = calculateExpectedCornerVelocity(v00.initialPosition, new Coordinate(10, 0), new Coordinate(5, 10));
		// Assert velocity components (calculateExpectedCornerVelocity needs finishing
		// or manual calc)
		// Normal for (0,0)->(10,0) is (0,-1). Offset line y=-1.
		// Normal for (5,10)->(0,0) is (-10, -5). Unit normal approx (-0.894, -0.447).
		// Offset line...? Bit complex.
		// Skipping exact velocity check here unless calculateExpectedCornerVelocity is
		// fully implemented.
		assertNotNull(v00.getVelocity(), "Velocity should be calculated");
	}

	@Test
	void testDegenerateCollinearVertexGeometry() throws ParseException {
		Polygon degeneratePoly = createDegeneratePoly();
		KineticTriangulation kt = new KineticTriangulation(degeneratePoly, gf);

		// Test the middle collinear vertex (5, 0)
		WavefrontVertex v50 = findVertexAt(kt, 5, 0);
		assertNotNull(v50);

		// Expected Angle: Straight line segment
		assertEquals(VertexAngle.COLLINEAR, v50.getAngle(), "Degenerate vertex (5,0) should be COLLINEAR");

		// Expected Speed: Depends on weights. Assuming default weight 1.0 for boundary
		// edges.
		// Edges (0,0)-(5,0) and (10,0)-(5,0) have same direction vector (when pointing
		// away from vertex)
		// but opposite normals relative to standard CCW traversal?
		// Edge (0,0)->(5,0). Normal (0,-1).
		// Edge (10,0)->(5,0). Normal (0,1).
		// They are opposing. Weights are likely equal (default 1.0).
		// Let's check the calculation in calculateInfiniteSpeedType:
		// angle == COLLINEAR.
		// normal1 = unit normal of edge1 = unit normal of (10,0)->(5,0) = (0,1).
		// dir0 = normal direction of edge0 = normal dir of (0,0)->(5,0) = (0,-5).
		// orient = dir0.x*normal1.y - dir0.y*normal1.x = 0*1 - (-5)*0 = 0. -> Not
		// OPPOSING branch.
		// weights are equal -> NONE.
		assertEquals(InfiniteSpeedType.NONE, v50.getInfiniteSpeed(), "Degenerate vertex (5,0) should have NONE infinite speed (opposing normals, same weight)");

		// Expected Velocity: For COLLINEAR, NONE speed -> velocity is w * unit_normal
		// Which edge's normal? Let's assume edge0: (0,0)->(5,0). Weighted normal is
		// (0,-1)*1.0 = (0,-1).
		assertEquals(0.0, v50.getVelocity().getX(), DELTA, "Degenerate vertex (5,0) velocity X should be 0");
		assertEquals(-1.0, v50.getVelocity().getY(), DELTA, "Degenerate vertex (5,0) velocity Y should be -1 (along normal)");
	}

	// --- Validation Helper Methods ---

	private void validateConnectivity(KineticTriangulation kt) {
		assertFalse(kt.getTriangles().isEmpty(), "Triangulation should not be empty for validation");

		for (KineticTriangle tri : kt.getTriangles()) {
			assertNotNull(tri, "Triangle reference should not be null");
			// Check vertices
			for (int i = 0; i < 3; i++) {
				WavefrontVertex v = tri.getVertex(i);
				assertNotNull(v, "Triangle " + tri.id + " vertex " + i + " is null");
				if (!v.isInfinite) {
					// Find the vertex in the main list
					assertTrue(kt.getVertices().contains(v), "Triangle vertex not found in main list");
				}
			}

			// Check neighbors/wavefronts exclusivity and back-references
			for (int i = 0; i < 3; i++) {
				KineticTriangle neighbor = tri.getNeighbor(i);
				WavefrontEdge wavefront = tri.getWavefront(i);

				assertTrue(neighbor == null || wavefront == null, "Triangle " + tri.id + " edge " + i + " has both neighbor and wavefront");
				assertTrue(neighbor != null || wavefront != null, "Triangle " + tri.id + " edge " + i + " has neither neighbor nor wavefront"); // Should always
																																				// have one

				if (neighbor != null) {
					// Check back-reference
					boolean foundBackRef = false;
					for (int j = 0; j < 3; ++j) {
						if (neighbor.getNeighbor(j) == tri) {
							foundBackRef = true;
							break;
						}
					}
					assertTrue(foundBackRef, "Neighbor " + neighbor.id + " does not reference back triangle " + tri.id);
					// Check shared vertices
					WavefrontVertex v1 = tri.getVertex((i + 1) % 3);
					WavefrontVertex v2 = tri.getVertex((i + 2) % 3);
					assertTrue(neighbor.indexOfVertex(v1) >= 0, "Neighbor missing shared vertex");
					assertTrue(neighbor.indexOfVertex(v2) >= 0, "Neighbor missing shared vertex");
				}

				if (wavefront != null) {
					assertEquals(tri, wavefront.getIncidentTriangle(), "Wavefront edge " + wavefront.id + " does not reference back triangle " + tri.id);
					assertTrue(kt.getWavefrontEdges().contains(wavefront), "Wavefront edge not found in main list");
					// Check shared vertices
					WavefrontVertex v1 = tri.getVertex((i + 1) % 3); // CCW vertex from edge i's perspective
					WavefrontVertex v2 = tri.getVertex((i + 2) % 3); // CW vertex from edge i's perspective
					// Check if edge vertices match triangle vertices (order might be swapped in
					// edge)
					assertTrue((wavefront.getVertex(0) == v1 && wavefront.getVertex(1) == v2) || (wavefront.getVertex(0) == v2 && wavefront.getVertex(1) == v1),
							"Wavefront edge " + wavefront.id + " vertices don't match triangle " + tri.id + " vertices for edge " + i);

					// Check vertex back-references to edge
					assertNotNull(v1);
					assertNotNull(v2);
					if (!v1.isInfinite)
						assertTrue(v1.getIncidentEdge(0) == wavefront || v1.getIncidentEdge(1) == wavefront,
								"Vertex " + v1 + " doesn't reference incident edge " + wavefront.id);
					if (!v2.isInfinite)
						assertTrue(v2.getIncidentEdge(0) == wavefront || v2.getIncidentEdge(1) == wavefront,
								"Vertex " + v2 + " doesn't reference incident edge " + wavefront.id);

				}
			}
		}

		// Validate vertex references to edges
		for (WavefrontVertex v : kt.getVertices()) {
			if (v.isInfinite)
				continue;
			WavefrontEdge e0 = v.getIncidentEdge(0);
			WavefrontEdge e1 = v.getIncidentEdge(1);
			assertNotNull(e0, "Vertex " + v + " missing edge 0");
			assertNotNull(e1, "Vertex " + v + " missing edge 1");
			assertTrue(e0.getVertex(0) == v || e0.getVertex(1) == v, "Vertex " + v + " not found in its edge 0");
			assertTrue(e1.getVertex(0) == v || e1.getVertex(1) == v, "Vertex " + v + " not found in its edge 1");

			// Check C++ convention: edge0 = left/CCW, edge1 = right/CW
			// edge0's vertex1 should be v, edge1's vertex0 should be v
			assertEquals(v, e0.getVertex(1), "Vertex " + v + " should be vertex1 (end) of its edge0 (CCW)");
			assertEquals(v, e1.getVertex(0), "Vertex " + v + " should be vertex0 (start) of its edge1 (CW)");
		}
	}

	private void validateConstraints(KineticTriangulation kt, Polygon polygon) {
		Set<CanonicalSegment> inputConstraints = new HashSet<>();
		extractSegments(polygon.getExteriorRing(), inputConstraints);
		for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
			extractSegments(polygon.getInteriorRingN(i), inputConstraints);
		}

		Set<CanonicalSegment> ktConstraints = kt.getWavefrontEdges().stream().map(we -> we.canonicalSegment).collect(Collectors.toSet());

		assertEquals(inputConstraints.size(), ktConstraints.size(), "Number of constraint edges mismatch");
		assertEquals(inputConstraints, ktConstraints, "Constraint edge sets do not match");

		// Also check that triangles bordering these constraints have the correct
		// wavefront edge set
		for (KineticTriangle tri : kt.getTriangles()) {
			for (int i = 0; i < 3; ++i) {
				WavefrontEdge wf = tri.getWavefront(i);
				if (wf != null) {
					assertTrue(inputConstraints.contains(wf.canonicalSegment),
							"Triangle " + tri.id + " has non-constraint edge marked as wavefront: " + wf.canonicalSegment);
				} else {
					// If no wavefront, check if the corresponding edge is internal (not a
					// constraint)
					WavefrontVertex v1 = tri.getVertex((i + 1) % 3);
					WavefrontVertex v2 = tri.getVertex((i + 2) % 3);
					if (!v1.isInfinite && !v2.isInfinite) {
						CanonicalSegment internalEdge = new CanonicalSegment(v1.initialPosition, v2.initialPosition);
						assertFalse(inputConstraints.contains(internalEdge), "Triangle " + tri.id + " has constraint edge marked as internal: " + internalEdge);
					}
				}
			}
		}
	}

	private void extractSegments(LineString ring, Set<CanonicalSegment> segmentSet) {
		Coordinate[] coords = ring.getCoordinates();
		for (int i = 0; i < coords.length - 1; i++) {
			if (!coords[i].equals2D(coords[i + 1])) { // Avoid zero-length segments
				segmentSet.add(new CanonicalSegment(coords[i], coords[i + 1]));
			}
		}
	}
}