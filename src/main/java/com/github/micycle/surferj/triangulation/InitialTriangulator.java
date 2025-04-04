package com.github.micycle.surferj.triangulation;

import org.locationtech.jts.geom.Polygon;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.Vertex;

import com.github.micycle.surferj.util.Constants;
import com.github.micycle.surferj.util.TinFourUtil;

import org.tinfour.common.SimpleTriangle;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Triangulate a JTS polygon using TinFour.
 */
public class InitialTriangulator {

	private final Polygon inputPolygon;
	private IIncrementalTin tin;
	private Map<IQuadEdge, Double> constraintWeights; // Store weights for original polygon edges
	private List<SimpleTriangle> triangles;

	public InitialTriangulator(Polygon polygon) {
		this.inputPolygon = Objects.requireNonNull(polygon);
		this.constraintWeights = new HashMap<>();
	}

	public void triangulate() {
		this.tin = TinFourUtil.triangulatePolygon(inputPolygon);
//		tin.triangles().spliterator()
		this.triangles = tin.getTriangles();
		storeConstraintWeights();
	}

	private void storeConstraintWeights() {
		// In this simplified version, all weights are 1.0
		// A real implementation would need to map TIN constraint edges
		// back to the original Polygon segments and store their associated weights.
		for (IQuadEdge edge : tin.getEdges()) {
			if (edge.isConstraint()) {
				constraintWeights.put(edge, Constants.DEFAULT_WEIGHT);
				// Also store for the twin, TinFour edges might not be consistently oriented
				constraintWeights.put(edge.getTwin(), Constants.DEFAULT_WEIGHT);
			}
		}
	}

	public IIncrementalTin getTin() {
		return tin;
	}

	public List<SimpleTriangle> getTriangles() {
		return triangles;
	}

	public double getConstraintWeight(IQuadEdge edge) {
		// Default to 1.0 if not found (should only happen for non-constraint edges)
		return constraintWeights.getOrDefault(edge, Constants.DEFAULT_WEIGHT);
	}

	public Map<IQuadEdge, Double> getConstraintWeightsMap() {
		return constraintWeights;
	}

	// Helper to get the Tinfour edge corresponding to two kinetic vertices
	// (Requires KineticVertex to store its original Tinfour Vertex)
	public IQuadEdge findTinEdge(Vertex v1, Vertex v2) {
		if (tin == null || v1 == null || v2 == null)
			return null;
		return tin.getEdge(v1, v2);
	}
}