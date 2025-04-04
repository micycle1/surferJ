package com.github.micycle.surferj.util;

import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Polygon;

import java.util.ArrayList;
import java.util.List;

public class TinFourUtil {

    public static IIncrementalTin triangulatePolygon(Polygon polygon) {
        IIncrementalTin tin = new IncrementalTin(1.0); // Adjust nominal spacing if needed

        // Add polygon vertices
        List<Vertex> vertices = new ArrayList<>();
        for (Coordinate coord : polygon.getExteriorRing().getCoordinates()) {
            // Avoid duplicate closing coordinate if present
            if (vertices.isEmpty() || !vertices.get(vertices.size() - 1).getCoordinate().equals(coord)) {
                vertices.add(new Vertex(coord.x, coord.y, 0)); // Z is not used here
            }
        }
        tin.add(vertices, null); // No constraints added yet

        // Add polygon edges as constraints
        List<IQuadEdge> constraintEdges = new ArrayList<>();
        for (int i = 0; i < vertices.size(); i++) {
            Vertex v1 = vertices.get(i);
            Vertex v2 = vertices.get((i + 1) % vertices.size()); // Wrap around
            // Find the edge in the TIN corresponding to this segment
             IQuadEdge edge = tin.getEdge(v1, v2);
             if (edge != null) {
                 // Check if it's already constrained (TinFour might add constraints implicitly sometimes)
                 if (!edge.isConstrainedRegionMember()) {
                    tin.addConstraintEdge(edge);
                    constraintEdges.add(edge);
                 }
             } else {
                 System.err.println("Warning: Could not find TIN edge for constraint: " + v1 + " -> " + v2);
                  // This might happen if the vertex list wasn't clean or TIN failed
             }
        }
         // TinFour might require constraints to form a polygon
         // tin.addConstraints(constraintEdges, true); // Mark as polygon constraints

        // TODO: Handle interior holes if needed

        tin.isBootstrapped(); // Ensure triangulation is built
        return tin;
    }
}