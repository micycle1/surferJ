package com.github.micycle1.surferj.collapse;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;

import com.github.micycle1.surferj.kinetics.WavefrontEdge;
import com.github.micycle1.surferj.kinetics.WavefrontVertex;
import com.github.micycle1.surferj.kinetics.WavefrontVertex.InfiniteSpeedType;
import com.github.micycle1.surferj.kinetics.WavefrontVertex.VertexAngle;

/**
 * Unit tests for the WavefrontVertex class.
 */
public class WavefrontVertexTest2 {

	private static final double DELTA = 1e-9;

	// Counter for unique vertex IDs if WavefrontVertex doesn't manage it itself for
	// tests
	private static long testIdCounter = 1;

	// Helper to create vertex using the test constructor
	private WavefrontVertex createTestVertex(double x, double y, double vx, double vy) {
		// If your WavefrontVertex constructor doesn't take ID, manage it externally or
		// use reflection if necessary
		WavefrontVertex v = new WavefrontVertex(x, y, vx, vy);
		// Reset id via reflection if needed for predictability, otherwise rely on
		// internal counter
		// setIdViaReflection(v, testIdCounter++); // Example if needed
		return v;
	}

	// Helper to create vertex using the test constructor
	private WavefrontVertex createTestVertex(Coordinate pos) {
		WavefrontVertex v = new WavefrontVertex(pos);
		return v;
	}

	// Helper to create a basic edge (doesn't handle incident triangle etc.)
	// Assumes WavefrontEdge constructor takes LineSegment and weight
	private WavefrontEdge createTestEdge(WavefrontVertex v0, WavefrontVertex v1, double weight) {
		if (v0 == null || v1 == null) {
			throw new NullPointerException("Null vertex for edge");
		}
		LineSegment seg = new LineSegment(v0.initialPosition, v1.initialPosition);
		WavefrontEdge edge = new WavefrontEdge(seg, weight);
		// Crucially link vertices back to edge FOR THE TEST SETUP
		// This mimics what KineticTriangulation Pass 3 does before geometry calc.
		edge.setVerticesRaw(v0, v1); // Assumes a method that ONLY sets vertex refs
		return edge;
	}

	// Helper to fully link two edges to a vertex, triggering recalculateGeometry
	private void setupVertexGeometry(WavefrontVertex vertex, WavefrontEdge edge0, WavefrontEdge edge1) {
		// Ensure edges are linked to the vertex
		// Order matters: edge0 is incoming CCW, edge1 is outgoing CCW
		vertex.setIncidentEdge(0, edge0); // This should trigger recalculateGeometry
		vertex.setIncidentEdge(1, edge1); // This should trigger recalculateGeometry again
	}

	// Reset counter before each test if necessary
	// @BeforeEach
	// void resetCounter() {
	// testIdCounter = 1;
	// }

	@Nested
	@DisplayName("Constructor and Factory Tests")
	class ConstructorTests {

		@Test
		void testSimpleConstructorWithVelocity() {
			WavefrontVertex v = createTestVertex(1.0, 2.0, 3.0, 4.0);
			assertEquals(1.0, v.initialPosition.getX(), "Initial X mismatch");
			assertEquals(2.0, v.initialPosition.getY(), "Initial Y mismatch");
			assertEquals(3.0, v.getVelocity().getX(), DELTA, "Velocity X mismatch");
			assertEquals(4.0, v.getVelocity().getY(), DELTA, "Velocity Y mismatch");
			assertFalse(v.isInfinite(), "Should be finite");
			assertNull(v.getIncidentEdge(0), "Edge 0 should be null initially");
			assertNull(v.getIncidentEdge(1), "Edge 1 should be null initially");
			// Check default calculated values (before recalculateGeometry runs)
			assertEquals(VertexAngle.COLLINEAR, v.getAngle(), "Default angle should be COLLINEAR");
			assertEquals(InfiniteSpeedType.NONE, v.getInfiniteSpeed(), "Default infinite speed should be NONE");
			assertFalse(v.hasStopped(), "Should not be stopped initially");
		}

		@Test
		void testSimpleConstructorWithPosition() {
			WavefrontVertex v = createTestVertex(new Coordinate(5, 6));
			assertEquals(5.0, v.initialPosition.getX(), "Initial X mismatch");
			assertEquals(6.0, v.initialPosition.getY(), "Initial Y mismatch");
			assertEquals(0.0, v.getVelocity().getX(), DELTA, "Default Velocity X should be 0");
			assertEquals(0.0, v.getVelocity().getY(), DELTA, "Default Velocity Y should be 0");
			assertFalse(v.isInfinite(), "Should be finite");
			assertNull(v.getIncidentEdge(0), "Edge 0 should be null initially");
			assertNull(v.getIncidentEdge(1), "Edge 1 should be null initially");
			assertEquals(VertexAngle.COLLINEAR, v.getAngle(), "Default angle should be COLLINEAR");
			assertEquals(InfiniteSpeedType.NONE, v.getInfiniteSpeed(), "Default infinite speed should be NONE");
		}

		@Test
		void testMakeInfinite_None() {
			// Assuming a static factory method like in C++
			// If not present, test the constructor/method that sets isInfinite=true
			WavefrontVertex v = WavefrontVertex.makeInfinite(InfiniteSpeedType.NONE); // Adjust if factory name differs
			assertTrue(v.isInfinite(), "Should be infinite");
			assertEquals(InfiniteSpeedType.NONE, v.getInfiniteSpeed(), "Infinite speed should be NONE");
			assertEquals(0.0, v.getVelocity().getX(), DELTA, "Infinite vertex velocity X should be 0");
			assertEquals(0.0, v.getVelocity().getY(), DELTA, "Infinite vertex velocity Y should be 0");
			// Angle might be default/undefined for infinite, check expectation
			// assertEquals(VertexAngle.COLLINEAR, v.getAngle(), "Infinite vertex default
			// angle");
		}

		@Test
		void testMakeInfinite_Opposing() {
			WavefrontVertex v = WavefrontVertex.makeInfinite(InfiniteSpeedType.OPPOSING);
			assertTrue(v.isInfinite(), "Should be infinite");
			assertEquals(InfiniteSpeedType.OPPOSING, v.getInfiniteSpeed(), "Infinite speed should be OPPOSING");
			assertEquals(0.0, v.getVelocity().getX(), DELTA, "Infinite vertex velocity X should be 0");
			assertEquals(0.0, v.getVelocity().getY(), DELTA, "Infinite vertex velocity Y should be 0");
		}

		@Test
		void testMakeInfinite_Weighted() {
			WavefrontVertex v = WavefrontVertex.makeInfinite(InfiniteSpeedType.WEIGHTED);
			assertTrue(v.isInfinite(), "Should be infinite");
			assertEquals(InfiniteSpeedType.WEIGHTED, v.getInfiniteSpeed(), "Infinite speed should be WEIGHTED");
			assertEquals(0.0, v.getVelocity().getX(), DELTA, "Infinite vertex velocity X should be 0");
			assertEquals(0.0, v.getVelocity().getY(), DELTA, "Infinite vertex velocity Y should be 0");
		}
	}

	@Nested
	@DisplayName("Geometric Property Calculation Tests")
	class GeometryCalculationTests {

		// Setup vertices for reuse
		WavefrontVertex v00 = createTestVertex(0, 0, 0, 0);
		WavefrontVertex v10 = createTestVertex(1, 0, 0, 0);
		WavefrontVertex v_10 = createTestVertex(-1, 0, 0, 0);
		WavefrontVertex v01 = createTestVertex(0, 1, 0, 0);
		WavefrontVertex v11 = createTestVertex(1, 1, 0, 0);
		WavefrontVertex v_11 = createTestVertex(-1, 1, 0, 0);
		WavefrontVertex v20 = createTestVertex(2, 0, 0, 0);

		@Test
		void testConvexCorner_AngleVelocitySpeed() {
			// Square corner at (0,0), edges have weight 1.0
			WavefrontEdge edge_v01_v00 = createTestEdge(v01, v00, 1.0); // Arriving CCW (edge0)
			WavefrontEdge edge_v00_v10 = createTestEdge(v00, v10, 1.0); // Leaving CCW (edge1)

			setupVertexGeometry(v00, edge_v01_v00, edge_v00_v10);

			assertEquals(VertexAngle.LEFT_TURN, v00.getAngle(), "Convex angle should be LEFT_TURN");
			assertEquals(InfiniteSpeedType.NONE, v00.getInfiniteSpeed(), "Convex angle should have NONE speed");
			assertEquals(1.0, v00.getVelocity().getX(), DELTA, "Convex angle velocity X mismatch");
			assertEquals(1.0, v00.getVelocity().getY(), DELTA, "Convex angle velocity Y mismatch");
		}

		@Test
		void testReflexCorner_AngleVelocitySpeed() {
			// Inner corner of L-shape at (0,0), edges have weight 1.0
			// Polygon points: (-1, 0), (0, 0), (0, -1) ...
			WavefrontVertex v0_1 = createTestVertex(0, -1, 0, 0); // Need this point

			WavefrontEdge edge_v_10_v00 = createTestEdge(v_10, v00, 1.0); // Arriving CCW (edge0)
			WavefrontEdge edge_v00_v0_1 = createTestEdge(v00, v0_1, 1.0); // Leaving CCW (edge1)

			setupVertexGeometry(v00, edge_v_10_v00, edge_v00_v0_1);

			assertEquals(VertexAngle.RIGHT_TURN, v00.getAngle(), "Reflex angle should be RIGHT_TURN");
			assertEquals(InfiniteSpeedType.NONE, v00.getInfiniteSpeed(), "Reflex angle should have NONE speed");
			// Velocity: offset x=0 by (-1,0)*1 -> x=-1; offset y=0 by (0,-1)*1 -> y=-1.
			// Intersect (-1,-1). Vel = (-1,-1)-(0,0)=(-1,-1)
			assertEquals(-1.0, v00.getVelocity().getX(), DELTA, "Reflex angle velocity X mismatch");
			assertEquals(-1.0, v00.getVelocity().getY(), DELTA, "Reflex angle velocity Y mismatch");
		}

		@Test
		void testCollinear_NoneSpeed_AngleVelocitySpeed() {
			// Collinear vertex at (1,0) between (0,0) and (2,0), weight 1.0
			WavefrontEdge edge_v00_v10 = createTestEdge(v00, v10, 1.0); // Arriving CCW (edge0)
			WavefrontEdge edge_v10_v20 = createTestEdge(v10, v20, 1.0); // Leaving CCW (edge1)

			setupVertexGeometry(v10, edge_v00_v10, edge_v10_v20);

			assertEquals(VertexAngle.COLLINEAR, v10.getAngle(), "Collinear angle should be COLLINEAR");
			// Check infinite speed calculation: angle==COLLINEAR, edge directions same,
			// weights same -> NONE
			assertEquals(InfiniteSpeedType.NONE, v10.getInfiniteSpeed(), "Collinear (same dir, same weight) should have NONE speed");
			// Check velocity calculation: COLLINEAR+NONE -> velocity = weighted outward
			// normal of edge0
			// edge0 segment (0,0)->(1,0). Dir(1,0). Outward Normal (0,-1). Weighted Normal
			// (0,-1)*1=(0,-1)
			assertEquals(0.0, v10.getVelocity().getX(), DELTA, "Collinear NONE velocity X mismatch");
			assertEquals(-1.0, v10.getVelocity().getY(), DELTA, "Collinear NONE velocity Y mismatch");
		}

		@Test
		void testCollinear_OpposingSpeed_AngleVelocitySpeed() {
			// Collinear vertex at (0,0), edges point opposite ways, weight 1.0
			WavefrontEdge edge_v_10_v00 = createTestEdge(v_10, v00, 1.0); // Arriving CCW (edge0) segment (-1,0)->(0,0)
			WavefrontEdge edge_v10_v00 = createTestEdge(v10, v00, 1.0); // This is the second edge incident on v00. Should be leaving?
																		// Let's use edge v00 -> v10 instead for standard setup.
			WavefrontEdge edge_v00_v10 = createTestEdge(v00, v10, 1.0); // Leaving CCW (edge1) segment (0,0)->(1,0)

			// We need edges arriving/leaving v00 that are collinear but opposing directions
			// e.g., edge (-1,0)->(0,0) and edge (1,0)->(0,0)
			// This setup requires modifying the edge direction assumptions slightly or the
			// test setup.
			// Let's simulate by setting up the edges directly for the vertex
			WavefrontEdge edge_in = createTestEdge(v_10, v00, 1.0); // Normal edge
			// Create an edge pointing inwards but conceptually linked as the "other" edge
			WavefrontEdge edge_in_opp = createTestEdge(v10, v00, 1.0); // Segment (1,0) -> (0,0)

			// Setup: edge_in is edge0, edge_in_opp is edge1 for calculation purposes
			// Need to bypass normal setupVertexGeometry if it assumes standard vertex
			// ordering
			v00.setIncidentEdge(0, edge_in);
			v00.setIncidentEdge(1, edge_in_opp); // This should trigger recalc
			v00.recalculateGeometry(); // Force recalc if setIncidentEdge doesn't guarantee it

			assertEquals(VertexAngle.COLLINEAR, v00.getAngle(), "Collinear opposing angle should be COLLINEAR");
			// Check infinite speed calc: angle==COLLINEAR, edge directions opposite ->
			// OPPOSING
			assertEquals(InfiniteSpeedType.OPPOSING, v00.getInfiniteSpeed(), "Collinear (opposing dir) should have OPPOSING speed");
			assertEquals(0.0, v00.getVelocity().getX(), DELTA, "Collinear OPPOSING velocity X should be 0");
			assertEquals(0.0, v00.getVelocity().getY(), DELTA, "Collinear OPPOSING velocity Y should be 0");
		}

		@Test
		void testCollinear_WeightedSpeed_AngleVelocitySpeed() {
			// Collinear vertex at (1,0) between (0,0) and (2,0), different weights
			WavefrontEdge edge_v00_v10 = createTestEdge(v00, v10, 1.0); // Arriving CCW (edge0), weight 1
			WavefrontEdge edge_v10_v20 = createTestEdge(v10, v20, 2.0); // Leaving CCW (edge1), weight 2

			setupVertexGeometry(v10, edge_v00_v10, edge_v10_v20);

			assertEquals(VertexAngle.COLLINEAR, v10.getAngle(), "Collinear weighted angle should be COLLINEAR");
			// Check infinite speed calculation: angle==COLLINEAR, edge directions same,
			// weights DIFFERENT -> WEIGHTED
			assertEquals(InfiniteSpeedType.WEIGHTED, v10.getInfiniteSpeed(), "Collinear (same dir, diff weight) should have WEIGHTED speed");
			assertEquals(0.0, v10.getVelocity().getX(), DELTA, "Collinear WEIGHTED velocity X should be 0");
			assertEquals(0.0, v10.getVelocity().getY(), DELTA, "Collinear WEIGHTED velocity Y should be 0");
		}

		@Test
		void testInfiniteVertex_IgnoresFiniteEdgeRecalculation() {
			WavefrontVertex vInf = WavefrontVertex.makeInfinite(InfiniteSpeedType.OPPOSING);
			WavefrontVertex v00 = createTestVertex(0, 0, 0, 0);
			WavefrontVertex v10 = createTestVertex(1, 0, 0, 0);

			WavefrontEdge edge0 = createTestEdge(v00, vInf, 1.0); // Attach finite edges
			WavefrontEdge edge1 = createTestEdge(vInf, v10, 1.0);

			// These calls should store refs but NOT trigger recalculateGeometry on vInf
			vInf.setIncidentEdge(0, edge0);
			vInf.setIncidentEdge(1, edge1);

			assertTrue(vInf.isInfinite());
			assertEquals(InfiniteSpeedType.OPPOSING, vInf.getInfiniteSpeed(), "Infinite speed should NOT change");
			assertEquals(0.0, vInf.getVelocity().getX(), DELTA, "Infinite velocity X should NOT change");
			assertEquals(0.0, vInf.getVelocity().getY(), DELTA, "Infinite velocity Y should NOT change");
			// Angle may or may not be defined, depends on makeInfinite implementation
		}
	}

	@Nested
	@DisplayName("State and Position Tests")
	class StateTests {
		WavefrontVertex v = createTestVertex(1, 2, 3, 4); // Starts at (1,2), vel (3,4)

		@Test
		void testGetPositionAt() {
			Coordinate p0 = v.getPositionAt(0.0);
			assertEquals(1.0, p0.getX(), DELTA);
			assertEquals(2.0, p0.getY(), DELTA);

			Coordinate p1 = v.getPositionAt(1.0);
			assertEquals(1.0 + 3.0 * 1.0, p1.getX(), DELTA); // 4.0
			assertEquals(2.0 + 4.0 * 1.0, p1.getY(), DELTA); // 6.0

			Coordinate p2 = v.getPositionAt(2.5);
			assertEquals(1.0 + 3.0 * 2.5, p2.getX(), DELTA); // 1.0 + 7.5 = 8.5
			assertEquals(2.0 + 4.0 * 2.5, p2.getY(), DELTA); // 2.0 + 10.0 = 12.0

			Coordinate p_neg = v.getPositionAt(-1.0); // Test negative time
			assertEquals(1.0 + 3.0 * (-1.0), p_neg.getX(), DELTA); // -2.0
			assertEquals(2.0 + 4.0 * (-1.0), p_neg.getY(), DELTA); // -2.0
		}

		@Nested
		@DisplayName("Predicate Tests")
		class PredicateTests {
			// Reuse vertices from geometry tests if needed, or setup new ones

			@Test
			void testIsReflexOrStraight() {
				WavefrontVertex convex = createTestVertex(0, 0, 0, 0);
				WavefrontVertex reflex = createTestVertex(0, 0, 0, 0);
				WavefrontVertex collinear = createTestVertex(0, 0, 0, 0);

				// Setup convex (LEFT_TURN)
				WavefrontEdge e0_c = createTestEdge(createTestVertex(0, 1, 0, 0), convex, 1.0);
				WavefrontEdge e1_c = createTestEdge(convex, createTestVertex(1, 0, 0, 0), 1.0);
				setupVertexGeometry(convex, e0_c, e1_c);
				assertEquals(VertexAngle.LEFT_TURN, convex.getAngle());
				assertFalse(convex.isReflexOrStraight(), "Convex vertex should NOT be reflex or straight");

				// Setup reflex (RIGHT_TURN)
				WavefrontEdge e0_r = createTestEdge(createTestVertex(-1, 0, 0, 0), reflex, 1.0);
				WavefrontEdge e1_r = createTestEdge(reflex, createTestVertex(0, -1, 0, 0), 1.0);
				setupVertexGeometry(reflex, e0_r, e1_r);
				assertEquals(VertexAngle.RIGHT_TURN, reflex.getAngle());
				assertTrue(reflex.isReflexOrStraight(), "Reflex vertex SHOULD be reflex or straight");

				// Setup collinear (COLLINEAR)
				WavefrontEdge e0_s = createTestEdge(createTestVertex(-1, 0, 0, 0), collinear, 1.0);
				WavefrontEdge e1_s = createTestEdge(collinear, createTestVertex(1, 0, 0, 0), 1.0);
				setupVertexGeometry(collinear, e0_s, e1_s);
				assertEquals(VertexAngle.COLLINEAR, collinear.getAngle());
				assertTrue(collinear.isReflexOrStraight(), "Collinear vertex SHOULD be reflex or straight");
			}

			// Add testIsConvexOrStraight() similar to above if needed
		}

		// --- Helper Methods ---

	}
}