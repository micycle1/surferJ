package com.github.micycle1.surferj.kinetics;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;

public class WavefrontVertexTest {

	private static final double DELTA = 1e-9;
	private WavefrontVertex centerVertex;
	private WavefrontVertex vertexP0; // Represents vertex before centerVertex on edge0
	private WavefrontVertex vertexP2; // Represents vertex after centerVertex on edge1

	// Helper to create vertices (position only, velocity not relevant for these
	// tests)
	private WavefrontVertex createVertex(double x, double y) {
		WavefrontVertex v = new WavefrontVertex(new Coordinate(x, y));
		return v;
	}

	// Helper to create edges and link incident pointers implicitly via
	// setVerticesAndUpdateAdj
	private WavefrontEdge createEdge(WavefrontVertex vA, WavefrontVertex vB, double weight) {
		if (vA == null || vB == null) {
			throw new NullPointerException("Cannot create edge with null vertex");
		}
		LineSegment seg = new LineSegment(vA.getInitialPosition(), vB.getInitialPosition());
		WavefrontEdge edge = new WavefrontEdge(seg, weight);
		edge.setVerticesAndUpdateAdj(vA, vB); // Links edge to vertices and vice-versa
		return edge;
	}

	@BeforeEach
	public void setUp() {
		centerVertex = createVertex(0, 0); // Central vertex for testing
	}

	/**
	 * Test convex angle (LEFT_TURN). Setup: edge0 arrives from left (-1,0) ->
	 * (0,0). edge1 leaves upwards (0,0) -> (0,1). Points for Orientation.index:
	 * p0=(-1,0), p1=(0,0), p2=(0,1). Turn from edge0 vector (1,0) to edge1 vector
	 * (0,1) is CCW -> LEFT_TURN.
	 */
	@Test
	public void testConvexAngle_LeftTurn() {
		vertexP0 = createVertex(-1, 0); // << Swapped geometry from original test
		vertexP2 = createVertex(0, 1);

		WavefrontEdge edge0 = createEdge(vertexP0, centerVertex, 1.0);
		WavefrontEdge edge1 = createEdge(centerVertex, vertexP2, 1.0);

		centerVertex.setIncidentEdges(edge0, edge1);

		assertEquals(WavefrontVertex.VertexAngle.LEFT_TURN, centerVertex.getAngle(), "Angle should be LEFT_TURN"); // << Correct expected angle
		assertEquals(WavefrontVertex.InfiniteSpeedType.NONE, centerVertex.getInfiniteSpeed(), "Infinite speed should be NONE");

		// Velocity calculation (previously in Reflex test):
		// line0: dir(1,0), n(0,1), wn(0,1). Shifted point (0,1). Line y=1.
		// line1: dir(0,1), n(-1,0), wn(-1,0). Shifted point (-1,0). Line x=-1.
		// Intersection = (-1,1). Velocity = intersect - initial = (-1,1).
		Coordinate vel = centerVertex.getVelocity();
		assertEquals(-1.0, vel.getX(), DELTA, "Velocity X should be -1");
		assertEquals(1.0, vel.getY(), DELTA, "Velocity Y should be 1");
	}

	/**
	 * Test reflex angle (RIGHT_TURN). Setup: edge0 arrives from right (1,0) ->
	 * (0,0). edge1 leaves upwards (0,0) -> (0,1). Points for Orientation.index:
	 * p0=(1,0), p1=(0,0), p2=(0,1). Turn from edge0 vector (-1,0) to edge1 vector
	 * (0,1) is CW -> RIGHT_TURN.
	 */
	@Test
	public void testReflexAngle_RightTurn() {
		vertexP0 = createVertex(1, 0); // << Swapped geometry from original test
		vertexP2 = createVertex(0, 1);

		WavefrontEdge edge0 = createEdge(vertexP0, centerVertex, 1.0);
		WavefrontEdge edge1 = createEdge(centerVertex, vertexP2, 1.0);

		centerVertex.setIncidentEdges(edge0, edge1);

		assertEquals(WavefrontVertex.VertexAngle.RIGHT_TURN, centerVertex.getAngle(), "Angle should be RIGHT_TURN"); // << Correct expected angle
		assertEquals(WavefrontVertex.InfiniteSpeedType.NONE, centerVertex.getInfiniteSpeed(), "Infinite speed should be NONE");

		// Velocity calculation (previously in Convex test):
		// line0: dir(-1,0), n(0,-1), wn(0,-1). Shifted point (0,-1). Line y=-1.
		// line1: dir(0,1), n(-1,0), wn(-1,0). Shifted point (-1,0). Line x=-1.
		// Intersection = (-1,-1). Velocity = intersect - initial = (-1,-1).
		Coordinate vel = centerVertex.getVelocity();
		assertEquals(-1.0, vel.getX(), DELTA, "Velocity X should be -1");
		assertEquals(-1.0, vel.getY(), DELTA, "Velocity Y should be -1");
	}

	// ... Keep the other collinear and position tests as they were in the last
	// revision ...
	/**
	 * Test collinear edges where normals should be opposing but code calculates
	 * NONE. edge0 arrives from right (1,0) -> (0,0). Dir(-1,0), Normal(0,-1). edge1
	 * leaves towards right (-1,0) <- (0,0). Dir(-1,0), Normal(0,-1). The actual
	 * segments used for calculation are (1,0)->(0,0) and (0,0)->(-1,0)
	 */
	@Test
	public void testCollinearOpposing_CalculatedAsNone() {
		vertexP0 = createVertex(1, 0);
		vertexP2 = createVertex(-1, 0); // Vertex defining end of edge1

		WavefrontEdge edge0 = createEdge(vertexP0, centerVertex, 1.0); // weight 1
		WavefrontEdge edge1 = createEdge(centerVertex, vertexP2, 1.0); // weight 1

		centerVertex.setIncidentEdges(edge0, edge1);

		assertEquals(WavefrontVertex.VertexAngle.COLLINEAR, centerVertex.getAngle(), "Angle should be COLLINEAR");
		// wn0=(0,-1), un1=(0,-1), orient=0, weights equal -> NONE
		assertEquals(WavefrontVertex.InfiniteSpeedType.NONE, centerVertex.getInfiniteSpeed(), "Infinite speed should be NONE (based on current code logic)");

		// Velocity for COLLINEAR/NONE is weighted normal of edge0
		Coordinate vel = centerVertex.getVelocity();
		assertEquals(0.0, vel.getX(), DELTA, "Velocity X should be 0");
		assertEquals(-1.0, vel.getY(), DELTA, "Velocity Y should be -1 (wn0)");
	}

	/** Test collinear edges, same direction, different weights. */
	@Test
	public void testCollinearWeighted() {
		vertexP0 = createVertex(-1, 0);
		vertexP2 = createVertex(1, 0); // Vertex defining end of edge1

		WavefrontEdge edge0 = createEdge(vertexP0, centerVertex, 1.0); // weight 1
		WavefrontEdge edge1 = createEdge(centerVertex, vertexP2, 2.0); // weight 2

		centerVertex.setIncidentEdges(edge0, edge1);

		assertEquals(WavefrontVertex.VertexAngle.COLLINEAR, centerVertex.getAngle(), "Angle should be COLLINEAR");
		// wn0=(0,1), un1=(0,1), orient=0. Weights differ (1!=2) -> WEIGHTED.
		assertEquals(WavefrontVertex.InfiniteSpeedType.WEIGHTED, centerVertex.getInfiniteSpeed(), "Infinite speed should be WEIGHTED");

		Coordinate vel = centerVertex.getVelocity();
		assertEquals(0.0, vel.getX(), DELTA, "Velocity X should be 0");
		assertEquals(0.0, vel.getY(), DELTA, "Velocity Y should be 0");
	}

	/** Test collinear edges, same direction, same weight. */
	@Test
	public void testCollinearSame() {
		vertexP0 = createVertex(-1, 0);
		vertexP2 = createVertex(1, 0);

		WavefrontEdge edge0 = createEdge(vertexP0, centerVertex, 1.0); // weight 1
		WavefrontEdge edge1 = createEdge(centerVertex, vertexP2, 1.0); // weight 1

		centerVertex.setIncidentEdges(edge0, edge1);

		assertEquals(WavefrontVertex.VertexAngle.COLLINEAR, centerVertex.getAngle(), "Angle should be COLLINEAR");
		// wn0=(0,1), un1=(0,1), orient=0. Weights equal (1==1) -> NONE.
		assertEquals(WavefrontVertex.InfiniteSpeedType.NONE, centerVertex.getInfiniteSpeed(), "Infinite speed should be NONE");

		// Velocity for COLLINEAR/NONE is weighted normal of edge0
		Coordinate vel = centerVertex.getVelocity();
		assertEquals(0.0, vel.getX(), DELTA, "Velocity X should be 0");
		assertEquals(1.0, vel.getY(), DELTA, "Velocity Y should be 1 (wn0)");
	}

	/** Test position changes over time based on calculated velocity. */
	@Test
	public void testPositionOverTime() {
		// --- Moving Vertex (Use LEFT_TURN setup) ---
		centerVertex = createVertex(0, 0); // Reset center vertex
		vertexP0 = createVertex(-1, 0); // Geometry for LEFT_TURN
		vertexP2 = createVertex(0, 1);
		WavefrontEdge edge0_move = createEdge(vertexP0, centerVertex, 1.0);
		WavefrontEdge edge1_move = createEdge(centerVertex, vertexP2, 1.0);
		centerVertex.setIncidentEdges(edge0_move, edge1_move);

		// Expected velocity (-1, 1) for LEFT_TURN setup
		Coordinate expectedVelMove = new Coordinate(-1, 1);
		assertEquals(expectedVelMove.getX(), centerVertex.getVelocity().getX(), DELTA);
		assertEquals(expectedVelMove.getY(), centerVertex.getVelocity().getY(), DELTA);

		Coordinate posAt1 = centerVertex.getPositionAt(1.0);
		Coordinate expectedPosAt1 = new Coordinate(centerVertex.getInitialPosition().getX() + expectedVelMove.getX() * 1.0,
				centerVertex.getInitialPosition().getY() + expectedVelMove.getY() * 1.0);
		assertEquals(expectedPosAt1.getX(), posAt1.getX(), DELTA, "Moving vertex pos X at t=1");
		assertEquals(expectedPosAt1.getY(), posAt1.getY(), DELTA, "Moving vertex pos Y at t=1");
		assertEquals(-1.0, posAt1.getX(), DELTA); // Check specific value
		assertEquals(1.0, posAt1.getY(), DELTA); // Check specific value

		// --- Stationary Vertex (Use WEIGHTED setup) ---
		centerVertex = createVertex(0, 0); // Reset center vertex
		vertexP0 = createVertex(-1, 0);
		vertexP2 = createVertex(1, 0);
		WavefrontEdge edge0_stat = createEdge(vertexP0, centerVertex, 1.0);
		WavefrontEdge edge1_stat = createEdge(centerVertex, vertexP2, 2.0); // Different weight
		centerVertex.setIncidentEdges(edge0_stat, edge1_stat);

		// Expected velocity (0, 0)
		Coordinate expectedVelStat = new Coordinate(0, 0);
		assertEquals(expectedVelStat.getX(), centerVertex.getVelocity().getX(), DELTA);
		assertEquals(expectedVelStat.getY(), centerVertex.getVelocity().getY(), DELTA);
		assertEquals(WavefrontVertex.InfiniteSpeedType.WEIGHTED, centerVertex.getInfiniteSpeed());

		Coordinate posAt100 = centerVertex.getPositionAt(100.0);
		Coordinate expectedPosAt100 = new Coordinate(centerVertex.getInitialPosition().getX() + expectedVelStat.getX() * 100.0,
				centerVertex.getInitialPosition().getY() + expectedVelStat.getY() * 100.0);
		assertEquals(expectedPosAt100.getX(), posAt100.getX(), DELTA, "Stationary vertex pos X at t=100");
		assertEquals(expectedPosAt100.getY(), posAt100.getY(), DELTA, "Stationary vertex pos Y at t=100");
		assertEquals(0.0, posAt100.getX(), DELTA); // Check specific value
		assertEquals(0.0, posAt100.getY(), DELTA); // Check specific value
	}
}