package com.github.micycle1.surferj;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.math.Vector2D;

import com.github.micycle1.surferj.kinetics.WavefrontEdge;
import com.github.micycle1.surferj.kinetics.WavefrontSupportingLine;
import com.github.micycle1.surferj.kinetics.WavefrontVertex;
import com.github.micycle1.surferj.kinetics.WavefrontVertex.InfiniteSpeedType;
import com.github.micycle1.surferj.kinetics.WavefrontVertex.VertexAngle;

import org.locationtech.jts.geom.Coordinate;

public class WavefrontVertexTest {
    private WavefrontVertex vertex;
    private WavefrontEdge edge0, edge1;
    private WavefrontSupportingLine line0, line1;

    // Test convex angle (left turn)
    @Test
    public void testConvexAngle() {
        line0 = new WavefrontSupportingLine(1, 0, 0, 0); // Horizontal, rightward
        line1 = new WavefrontSupportingLine(0, 1, 0, 0); // Vertical, upward
        edge0 = new WavefrontEdge(line0);
        edge1 = new WavefrontEdge(line1);
        vertex = new WavefrontVertex(new Coordinate(0, 0));
        vertex.setIncidentEdges(edge0, edge1);

        assertEquals(VertexAngle.LEFT_TURN, vertex.getAngle());
        assertEquals(InfiniteSpeedType.NONE, vertex.getInfiniteSpeed());
        Coordinate vel = vertex.getVelocity();
        assertTrue(vel.getX() > 0, "Velocity X should be positive");
        assertTrue(vel.getY() > 0, "Velocity Y should be positive");
    }

    // Test reflex angle (right turn)
    @Test
    public void testReflexAngle() {
        line0 = new WavefrontSupportingLine(1, 0, 0, 1); // Horizontal, rightward
        line1 = new WavefrontSupportingLine(0, -1, 0, 1); // Vertical, downward
        edge0 = new WavefrontEdge(line0);
        edge1 = new WavefrontEdge(line1);
        vertex = new WavefrontVertex(new Coordinate(0, 0));
        vertex.setIncidentEdges(edge0, edge1);

        assertEquals(VertexAngle.RIGHT_TURN, vertex.getAngle());
        assertEquals(InfiniteSpeedType.NONE, vertex.getInfiniteSpeed());
    }

    // Test collinear edges moving towards each other (infinite speed)
    @Test
    public void testCollinearOpposing() {
        line0 = new WavefrontSupportingLine(1, 0, 0, 1); // Rightward
        line1 = new WavefrontSupportingLine(-1, 0, 0, 1); // Leftward
        edge0 = new WavefrontEdge(line0);
        edge1 = new WavefrontEdge(line1);
        vertex = new WavefrontVertex(new Coordinate(0, 0));
        vertex.setIncidentEdges(edge0, edge1);

        assertEquals(VertexAngle.COLLINEAR, vertex.getAngle());
        assertEquals(InfiniteSpeedType.OPPOSING, vertex.getInfiniteSpeed());
        Coordinate vel = vertex.getVelocity();
        assertEquals( 0, vel.getX(), 1e-12, "Velocity X should be zero");
        assertEquals( 0, vel.getY(), 1e-12, "Velocity Y should be zero");
    }

    // Test collinear edges with different weights (infinite speed)
    @Test
    public void testCollinearWeighted() {
        line0 = new WavefrontSupportingLine(1, 0, 0, 1); // Rightward, weight 1
        line1 = new WavefrontSupportingLine(1, 0, 0, 2); // Rightward, weight 2
        edge0 = new WavefrontEdge(line0);
        edge1 = new WavefrontEdge(line1);
        vertex = new WavefrontVertex(new Coordinate(0, 0));
        vertex.setIncidentEdges(edge0, edge1);

        assertEquals(VertexAngle.COLLINEAR, vertex.getAngle());
        assertEquals(InfiniteSpeedType.WEIGHTED, vertex.getInfiniteSpeed());
        Coordinate vel = vertex.getVelocity();
        assertEquals( 0, vel.getX(), 1e-12, "Velocity X should be zero");
        assertEquals( 0, vel.getY(), 1e-12, "Velocity Y should be zero");
    }

    // Test collinear edges with same weight and direction
    @Test
    public void testCollinearSame() {
        line0 = new WavefrontSupportingLine(0, 0, 1, 0); // Rightward, weight 1
        line1 = new WavefrontSupportingLine(1, 0, 0, 1); // Rightward, weight 1
        edge0 = new WavefrontEdge(line0);
        edge1 = new WavefrontEdge(line1);
        vertex = new WavefrontVertex(new Coordinate(0, 0));
        vertex.setIncidentEdges(edge0, edge1);

        assertEquals(VertexAngle.COLLINEAR, vertex.getAngle());
        assertEquals(InfiniteSpeedType.NONE, vertex.getInfiniteSpeed());
        Coordinate vel = vertex.getVelocity();
        assertEquals(0, vel.getX(), 1e-12, "Velocity X should be zero");
        assertTrue(vel.getY() != 0, "Velocity Y should be non-zero (normal)");
    }

    // Test position changes over time
    @Test
    public void testPositionOverTime() {
        // Moving vertex
        line0 = new WavefrontSupportingLine(1, 0, 0, 1); // Horizontal, rightward
        line1 = new WavefrontSupportingLine(0, 1, 0, 1); // Vertical, upward
        edge0 = new WavefrontEdge(line0);
        edge1 = new WavefrontEdge(line1);
        vertex = new WavefrontVertex(new Coordinate(0, 0));
        vertex.setIncidentEdges(edge0, edge1);
        Coordinate posAt1 = vertex.getPositionAt(1);
        assertNotEquals(0, posAt1.getX(), 1e-12);
        assertNotEquals(0, posAt1.getY(), 1e-12);

        // Stationary vertex (infinite speed)
        line0 = new WavefrontSupportingLine(1, 0, 0, 1); // Rightward
        line1 = new WavefrontSupportingLine(-1, 0, 0, 1); // Leftward
        edge0 = new WavefrontEdge(line0);
        edge1 = new WavefrontEdge(line1);
        vertex = new WavefrontVertex(new Coordinate(0, 0));
        vertex.setIncidentEdges(edge0, edge1);
        Coordinate posAt100 = vertex.getPositionAt(100);
        assertEquals(0, posAt100.getX(), 1e-12);
        assertEquals(0, posAt100.getY(), 1e-12);
    }
}