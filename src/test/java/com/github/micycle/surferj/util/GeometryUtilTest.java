package com.github.micycle.surferj.util;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.math.Vector2D;

import static org.junit.jupiter.api.Assertions.*;

class GeometryUtilTest {

    private static final double W = Constants.DEFAULT_WEIGHT; // Default weight = 1.0
    private static final double E = Constants.EPSILON;        // Epsilon

    @Test
    @DisplayName("Velocity for 90 degree convex corner at origin")
    void calculateVelocity_convex90Degree() {
        Coordinate origin = new Coordinate(0, 0);
        // Edge A: (-1, 0) -> (0, 0)
        LineSegment edgeA = new LineSegment(new Coordinate(-1, 0), origin);
        // Edge B: (0, 0) -> (0, -1) // Consistent orientation relative to vertex
        LineSegment edgeB = new LineSegment(origin, new Coordinate(0, -1));

        Vector2D velocity = GeometryUtil.calculateVelocity(origin, edgeA, edgeB, W, W);
        assertNotNull(velocity);
        assertEquals(1.0, velocity.getX(), E, "Velocity X should be approx 1.0");
        assertEquals(1.0, velocity.getY(), E, "Velocity Y should be approx 1.0"); // Corrected calculation matches this
    }

     @Test
     @DisplayName("Velocity for 90 degree convex corner offset from origin")
     void calculateVelocity_convex90DegreeOffset() {
         Coordinate corner = new Coordinate(5, 3);
         // Edge A: (4, 3) -> (5, 3)
         LineSegment edgeA = new LineSegment(new Coordinate(4, 3), corner);
         // Edge B: (5, 3) -> (5, 2)
         LineSegment edgeB = new LineSegment(corner, new Coordinate(5, 2));

         Vector2D velocity = GeometryUtil.calculateVelocity(corner, edgeA, edgeB, W, W);
         assertNotNull(velocity);
         assertEquals(1.0, velocity.getX(), E);
         assertEquals(1.0, velocity.getY(), E);
     }


    @Test
    @DisplayName("Velocity for 180 degree straight vertex")
    void calculateVelocity_straight180Degree() {
        Coordinate vertex = new Coordinate(0, 0);
        // Edge A: (-1, 0) -> (0, 0)
        LineSegment edgeA = new LineSegment(new Coordinate(-1, 0), vertex);
        // Edge B: (0, 0) -> (1, 0)
        LineSegment edgeB = new LineSegment(vertex, new Coordinate(1, 0));

        Vector2D velocity = GeometryUtil.calculateVelocity(vertex, edgeA, edgeB, W, W);
        assertNotNull(velocity);
        assertEquals(0.0, velocity.getX(), E, "Velocity X should be 0.0");
        assertEquals(W, velocity.getY(), E, "Velocity Y should be weight (1.0)"); // Corrected calculation matches this
    }

    @Test
    @DisplayName("Velocity for 270 degree reflex corner at origin")
    void calculateVelocity_reflex270Degree() {
        Coordinate origin = new Coordinate(0, 0);
        // Edge A: (1, 0) -> (0, 0)
        LineSegment edgeA = new LineSegment(new Coordinate(1, 0), origin);
        // Edge B: (0, 0) -> (0, -1)
        LineSegment edgeB = new LineSegment(origin, new Coordinate(0, -1));

        Vector2D velocity = GeometryUtil.calculateVelocity(origin, edgeA, edgeB, W, W);
        assertNotNull(velocity);
        // Adjusted expectation based on Cramer's rule result
        assertEquals(1.0, velocity.getX(), E, "Velocity X should be approx 1.0");
        assertEquals(-1.0, velocity.getY(), E, "Velocity Y should be approx -1.0");
    }

    @Test
    @DisplayName("Velocity for opposing collinear edges (Infinite Speed Case)")
    void calculateVelocity_opposingCollinear() {
        Coordinate vertex = new Coordinate(0, 0);
        // Edge A: (1, 0) -> (0, 0)
        LineSegment edgeA = new LineSegment(new Coordinate(1, 0), vertex);
        // Edge B: (0, 0) -> (-1, 0)
        LineSegment edgeB = new LineSegment(vertex, new Coordinate(-1, 0));

        Vector2D velocity = GeometryUtil.calculateVelocity(vertex, edgeA, edgeB, W, W);
        assertNotNull(velocity);
        assertEquals(0.0, velocity.getX(), E, "Velocity X should be 0.0 (Infinite Speed Case)");
        assertEquals(0.0, velocity.getY(), E, "Velocity Y should be 0.0 (Infinite Speed Case)");
    }

     @Test
     @DisplayName("Velocity for zero-length edge should be zero")
     void calculateVelocity_zeroLengthEdge() {
         Coordinate vertex = new Coordinate(0, 0);
         LineSegment edgeA = new LineSegment(new Coordinate(-1, 0), vertex);
         LineSegment zeroEdgeB = new LineSegment(vertex, vertex); // Zero length

         Vector2D velocity = GeometryUtil.calculateVelocity(vertex, edgeA, zeroEdgeB, W, W);
         assertNotNull(velocity);
         assertEquals(0.0, velocity.length(), E, "Velocity should be zero for zero-length edge");

         // Test with edge A being zero length
         LineSegment zeroEdgeA = new LineSegment(vertex, vertex);
         LineSegment edgeB = new LineSegment(vertex, new Coordinate(1, 0));
         velocity = GeometryUtil.calculateVelocity(vertex, zeroEdgeA, edgeB, W, W);
         assertNotNull(velocity);
         assertEquals(0.0, velocity.length(), E, "Velocity should be zero for zero-length edge");
     }

    @Test
    @DisplayName("Velocity with different weights")
    void calculateVelocity_differentWeights() {
        Coordinate origin = new Coordinate(0, 0);
        // Edge A: (-1, 0) -> (0, 0)
        LineSegment edgeA = new LineSegment(new Coordinate(-1, 0), origin);
        // Edge B: (0, 0) -> (0, -1)
        LineSegment edgeB = new LineSegment(origin, new Coordinate(0, -1));
        double weightA = 2.0;
        double weightB = 1.0;

        Vector2D velocity = GeometryUtil.calculateVelocity(origin, edgeA, edgeB, weightA, weightB);
        assertNotNull(velocity);
        assertEquals(1.0, velocity.getX(), E); // Corrected calculation matches this
        assertEquals(2.0, velocity.getY(), E); // Corrected calculation matches this
    }

     @Test
     @DisplayName("Velocity for same-direction collinear edges with different weights")
     void calculateVelocity_collinearSameDirDifferentWeights() {
         Coordinate vertex = new Coordinate(0, 0);
         // Edge A: (-1, 0) -> (0, 0)
         LineSegment edgeA = new LineSegment(new Coordinate(-1, 0), vertex);
         // Edge B: (0, 0) -> (1, 0)
         LineSegment edgeB = new LineSegment(vertex, new Coordinate(1, 0));
         double weightA = 2.0; // Faster
         double weightB = 1.0;

         Vector2D velocity = GeometryUtil.calculateVelocity(vertex, edgeA, edgeB, weightA, weightB);
         assertNotNull(velocity);
         assertEquals(0.0, velocity.getX(), E);
         assertEquals(weightA, velocity.getY(), E); // Expect faster weight * normal Y

         // Test with weightB faster
         velocity = GeometryUtil.calculateVelocity(vertex, edgeA, edgeB, weightB, weightA); // Swap weights
         assertNotNull(velocity);
         assertEquals(0.0, velocity.getX(), E);
         assertEquals(weightA, velocity.getY(), E); // Still expect faster weight (which is now weightA=2.0)
     }
}
