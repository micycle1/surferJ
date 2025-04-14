package com.github.micycle.surferj2.kinetics;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;

import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class SupportLineTest {

	private final GeometryFactory gf = new GeometryFactory();
	private final WKTReader reader = new WKTReader(gf);

	// --- Test Polygons (reuse from previous test class if desired) ---
	private Polygon createSquare() throws ParseException {
		return (Polygon) reader.read("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))");
	}

	private Polygon createSquareWithHole() throws ParseException {
		return (Polygon) reader.read("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0), (3 3, 7 3, 7 7, 3 7, 3 3))");
	}

	private Polygon createSingleTrianglePoly() throws ParseException {
		return (Polygon) reader.read("POLYGON ((0 0, 10 0, 5 10, 0 0))");
	}

	// --- Tests ---

	@Test
	void testSquareWavefrontEdgeCreation() throws ParseException {
		Polygon square = createSquare();
		KineticTriangulation kt = new KineticTriangulation(square, gf);

		assertEquals(4, kt.getWavefrontEdges().size(), "Should create 4 wavefront edges for square boundary");
		assertEquals(4, kt.getEdgeMap().size(), "Edge map should contain 4 unique constraint segments");

		// Verify segments in the edge map
		Set<CanonicalSegment> expectedSegments = Set.of(new CanonicalSegment(new Coordinate(0, 0), new Coordinate(10, 0)),
				new CanonicalSegment(new Coordinate(10, 0), new Coordinate(10, 10)), new CanonicalSegment(new Coordinate(10, 10), new Coordinate(0, 10)),
				new CanonicalSegment(new Coordinate(0, 10), new Coordinate(0, 0)));
		assertEquals(expectedSegments, kt.getEdgeMap().keySet());

		// Verify edges in the list correspond to map entries
		for (WavefrontEdge we : kt.getWavefrontEdges()) {
			assertNotNull(we.supportingLine, "WavefrontEdge must have a supporting line");
			assertEquals(1.0, we.getWeight(), "Default weight should be 1.0"); // Assuming default
			CanonicalSegment cs = we.canonicalSegment;
			assertTrue(kt.getEdgeMap().containsKey(cs), "Edge segment not found in map");
			assertSame(we, kt.getEdgeMap().get(cs), "Edge in list differs from edge in map");
		}
	}

	@Test
	void testSquareTriangleWavefrontLinks() throws ParseException {
		Polygon square = createSquare();
		KineticTriangulation kt = new KineticTriangulation(square, gf);

		assertEquals(2, kt.getTriangles().size()); // Expect 2 triangles

		for (KineticTriangle tri : kt.getTriangles()) {
			int constraintCount = 0;
			int neighborCount = 0;
			for (int i = 0; i < 3; i++) {
				WavefrontEdge wf = tri.getWavefront(i);
				KineticTriangle neighbor = tri.getNeighbor(i);

				if (wf != null) {
					constraintCount++;
					assertNull(neighbor, "Triangle " + tri.id + " edge " + i + " has both wavefront and neighbor");
					// Check back reference
					assertTrue(wf.getIncidentTriangle() == tri || kt.findOtherTriangleWithWavefrontEdge(wf) == tri,
							"Wavefront edge " + wf.id + " does not reference triangle " + tri.id);
					// Check if the segment is actually a boundary segment
					assertTrue(kt.getEdgeMap().containsKey(wf.canonicalSegment), "Triangle " + tri.id + " has non-boundary edge marked as wavefront");
				} else {
					neighborCount++;
					assertNotNull(neighbor, "Triangle " + tri.id + " edge " + i + " has neither wavefront nor neighbor");
					// Check if the edge is internal
					WavefrontVertex v1 = tri.getVertex((i + 1) % 3);
					WavefrontVertex v2 = tri.getVertex((i + 2) % 3);
					assertFalse(kt.getEdgeMap().containsKey(new CanonicalSegment(v1.initialPosition, v2.initialPosition)),
							"Triangle " + tri.id + " has boundary edge marked as internal neighbor edge");
				}
			}
			// Each triangle in a simple square CDT should have 2 constraint edges and 1
			// neighbor edge
			assertEquals(2, constraintCount, "Triangle in square CDT should have 2 constraint edges");
			assertEquals(1, neighborCount, "Triangle in square CDT should have 1 neighbor edge");
		}
	}

	@Test
	void testSquareEdgeVertexLinks() throws ParseException {
		Polygon square = createSquare();
		KineticTriangulation kt = new KineticTriangulation(square, gf);

		assertEquals(4, kt.getWavefrontEdges().size());
		assertEquals(4, kt.getVertices().size());

		Map<Coordinate, WavefrontVertex> coordToVertex = kt.getVertices().stream().collect(Collectors.toMap(v -> v.initialPosition, v -> v));

		for (WavefrontEdge we : kt.getWavefrontEdges()) {
			WavefrontVertex v0 = we.getVertex(0); // CCW vertex from triangle's perspective
			WavefrontVertex v1 = we.getVertex(1); // CW vertex from triangle's perspective

			assertNotNull(v0, "WavefrontEdge " + we.id + " vertex 0 is null");
			assertNotNull(v1, "WavefrontEdge " + we.id + " vertex 1 is null");
			assertNotEquals(v0, v1, "WavefrontEdge vertices are the same");

			// Check consistency with segment coordinates
			LineSegment seg = we.getSegment();
			assertTrue((v0.initialPosition.equals2D(seg.p0) && v1.initialPosition.equals2D(seg.p1))
					|| (v0.initialPosition.equals2D(seg.p1) && v1.initialPosition.equals2D(seg.p0)), "Edge vertices don't match segment endpoints");

			// Check vertex back-references
			// v0 (CCW from tri) is END of edge 0 (CCW from vert), START of edge 1 (CW from
			// vert)
			// v1 (CW from tri) is START of edge 0 (CCW from vert), END of edge 1 (CW from
			// vert)
			assertSame(we, v0.getIncidentEdge(1), "Vertex " + v0.id + " edge 1 (CW) should be " + we.id);
			assertSame(we, v1.getIncidentEdge(0), "Vertex " + v1.id + " edge 0 (CCW) should be " + we.id);
		}

		// Check all vertices have exactly two edges
		for (WavefrontVertex v : kt.getVertices()) {
			assertNotNull(v.getIncidentEdge(0), "Vertex " + v.id + " missing edge 0");
			assertNotNull(v.getIncidentEdge(1), "Vertex " + v.id + " missing edge 1");
			assertNotEquals(v.getIncidentEdge(0), v.getIncidentEdge(1), "Vertex " + v.id + " incident edges are the same");
		}
	}

	@Test
	void testHoleWavefrontEdgeCreation() throws ParseException {
		Polygon squareWithHole = createSquareWithHole();
		KineticTriangulation kt = new KineticTriangulation(squareWithHole, gf);

		// 4 outer + 4 inner edges
		assertEquals(8, kt.getWavefrontEdges().size(), "Should create 8 wavefront edges for square with hole");
		assertEquals(8, kt.getEdgeMap().size(), "Edge map should contain 8 unique constraint segments");

		// Verify some segments (example)
		assertTrue(kt.getEdgeMap().containsKey(new CanonicalSegment(new Coordinate(0, 0), new Coordinate(10, 0))), "Outer edge missing");
		assertTrue(kt.getEdgeMap().containsKey(new CanonicalSegment(new Coordinate(3, 3), new Coordinate(7, 3))), "Inner (hole) edge missing");

		// Check links (similar structure to square test, but more complex)
		validateWavefrontLinks(kt);
		validateEdgeVertexLinks(kt);
	}

	@Test
	void testSingleTriangleWavefrontEdgeCreation() throws ParseException {
		Polygon singleTriPoly = createSingleTrianglePoly();
		KineticTriangulation kt = new KineticTriangulation(singleTriPoly, gf);

		assertEquals(3, kt.getWavefrontEdges().size());
		assertEquals(3, kt.getEdgeMap().size());
		assertEquals(1, kt.getTriangles().size());

		KineticTriangle tri = kt.getTriangles().get(0);
		assertNotNull(tri.getWavefront(0));
		assertNotNull(tri.getWavefront(1));
		assertNotNull(tri.getWavefront(2));
		assertNull(tri.getNeighbor(0));
		assertNull(tri.getNeighbor(1));
		assertNull(tri.getNeighbor(2));

		validateWavefrontLinks(kt);
		validateEdgeVertexLinks(kt);
	}

	// --- Helper Validation Methods (can be shared or adapted) ---

	private void validateWavefrontLinks(KineticTriangulation kt) {
		for (KineticTriangle tri : kt.getTriangles()) {
			for (int i = 0; i < 3; i++) {
				WavefrontEdge wf = tri.getWavefront(i);
				if (wf != null) {
					assertTrue(wf.getIncidentTriangle() == tri || kt.findOtherTriangleWithWavefrontEdge(wf) == tri,
							"Wavefront edge " + wf.id + " does not reference triangle " + tri.id);
					assertTrue(kt.getEdgeMap().containsKey(wf.canonicalSegment), "Triangle " + tri.id + " has non-boundary edge marked as wavefront");
				}
			}
		}
	}

	private void validateEdgeVertexLinks(KineticTriangulation kt) {
		for (WavefrontEdge we : kt.getWavefrontEdges()) {
			WavefrontVertex v0 = we.getVertex(0);
			WavefrontVertex v1 = we.getVertex(1);
			assertNotNull(v0);
			assertNotNull(v1);
			assertSame(we, v0.getIncidentEdge(1), "Vertex " + v0.id + " edge 1 (CW) link wrong for " + we.id);
			assertSame(we, v1.getIncidentEdge(0), "Vertex " + v1.id + " edge 0 (CCW) link wrong for " + we.id);
		}
		for (WavefrontVertex v : kt.getVertices()) {
			if (!v.isInfinite) {
				assertNotNull(v.getIncidentEdge(0));
				assertNotNull(v.getIncidentEdge(1));
			}
		}
	}

}