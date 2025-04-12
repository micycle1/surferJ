package com.github.micycle.surferj.util;

import org.tinfour.common.IConstraint;
import org.tinfour.common.PolygonConstraint;
import org.tinfour.common.Vertex;
import org.tinfour.standard.IncrementalTin;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class TinFourUtil {

//	public static IIncrementalTin triangulatePolygonOld(Polygon polygon) {
//		IIncrementalTin tin = new IncrementalTin(1.0); // Adjust nominal spacing if needed
//
//		// Add polygon vertices
//		List<Vertex> vertices = new ArrayList<>();
//		for (Coordinate coord : polygon.getExteriorRing().getCoordinates()) {
//			// Avoid duplicate closing coordinate if present
//			if (vertices.isEmpty() || !vertices.get(vertices.size() - 1).getCoordinate().equals(coord)) {
//				vertices.add(new Vertex(coord.x, coord.y, 0)); // Z is not used here
//			}
//		}
//		tin.add(vertices, null); // No constraints added yet
//
//		// Add polygon edges as constraints
//		List<IQuadEdge> constraintEdges = new ArrayList<>();
//		for (int i = 0; i < vertices.size(); i++) {
//			Vertex v1 = vertices.get(i);
//			Vertex v2 = vertices.get((i + 1) % vertices.size()); // Wrap around
//			// Find the edge in the TIN corresponding to this segment
//			IQuadEdge edge = tin.getEdge(v1, v2);
//			if (edge != null) {
//				// Check if it's already constrained (TinFour might add constraints implicitly
//				// sometimes)
//				if (!edge.isConstrainedRegionMember()) {
//					tin.addConstraintEdge(edge);
//					constraintEdges.add(edge);
//				}
//			} else {
//				System.err.println("Warning: Could not find TIN edge for constraint: " + v1 + " -> " + v2);
//				// This might happen if the vertex list wasn't clean or TIN failed
//			}
//		}
//		// TinFour might require constraints to form a polygon
//		// tin.addConstraints(constraintEdges, true); // Mark as polygon constraints
//
//		// TODO: Handle interior holes if needed
//
//		tin.isBootstrapped(); // Ensure triangulation is built
//		return tin;
//	}

	/**
	 * Creates an incremental tin from a JTS Polygon, adding vertices and
	 * constraints for both exterior and interior rings with proper orientation.
	 * 
	 * @param polygon The JTS Polygon to convert to a tin
	 * @param pretty  Whether to apply pretty constraint options
	 * @return A fully constructed IncrementalTin with vertices and constraints
	 */
	public static IncrementalTin triangulatePolygon(Polygon polygon) {
		final IncrementalTin tin = new IncrementalTin(1);

		// Extract all coordinates from the geometry to add as vertices
		final Coordinate[] coords = polygon.getCoordinates();
		final List<Vertex> vertices = new ArrayList<>();
		for (int i = 0; i < coords.length; i++) {
			vertices.add(new Vertex(coords[i].x, coords[i].y, Double.NaN, i));
		}

		// Sort vertices to prevent degenerate insertion
//        HilbertSort hs = new HilbertSort();
//        hs.sort(vertices);

		// Add all vertices to the tin
		tin.add(vertices, null);

		// Create constraints from the linear components
		List<IConstraint> constraints = new ArrayList<>();

		for (int n = 0; n < polygon.getNumGeometries(); n++) {
			Polygon geom = (Polygon) polygon.getGeometryN(n);
			if (geom instanceof Polygon) { // in case multipolygon?
				Polygon poly = (Polygon) geom;

				// Handle exterior ring
				LineString exteriorRing = poly.getExteriorRing();
				if (exteriorRing != null && !exteriorRing.isEmpty()) {
					final List<Vertex> exteriorPoints = new ArrayList<>();
					Coordinate[] exteriorCoords = exteriorRing.getCoordinates();

					for (Coordinate coord : exteriorCoords) {
						exteriorPoints.add(new Vertex(coord.x, coord.y, Double.NaN));
					}

					// Ensure exterior ring is CCW
					if (!Orientation.isCCWArea(exteriorCoords)) {
						Collections.reverse(exteriorPoints);
					}

					constraints.add(new PolygonConstraint(exteriorPoints));
				}

				// Handle interior rings (holes)
				for (int i = 0; i < poly.getNumInteriorRing(); i++) {
					LineString interiorRing = poly.getInteriorRingN(i);
					if (interiorRing != null && !interiorRing.isEmpty()) {
						final List<Vertex> interiorPoints = new ArrayList<>();
						Coordinate[] interiorCoords = interiorRing.getCoordinates();

						for (Coordinate coord : interiorCoords) {
							interiorPoints.add(new Vertex(coord.x, coord.y, Double.NaN));
						}

						// Ensure interior ring is CW
						if (Orientation.isCCWArea(interiorCoords)) {
							Collections.reverse(interiorPoints);
						}

						constraints.add(new PolygonConstraint(interiorPoints));
					}
				}
			}
		}

		// Add constraints to the tin
		if (!constraints.isEmpty()) {
			tin.addConstraints(constraints, false); // NOTE false
		}

		return tin;
	}

}