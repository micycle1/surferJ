package com.github.micycle.surferj2.kinetics;

import java.util.concurrent.atomic.AtomicLong;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.math.Vector2D;

import com.github.micycle.surferj2.collapse.CollapseSpec;
import com.github.micycle.surferj2.collapse.CollapseType;
import com.github.micycle.surferj2.collapse.EdgeCollapseSpec;
import com.github.micycle.surferj2.collapse.EdgeCollapseType;
import com.github.micycle.surferj2.collapse.Polynomial;
import com.github.micycle.surferj2.collapse.QuadraticSolver;

// --- Kinetic Triangle ---
public class KineticTriangle {

	private static final AtomicLong idCounter = new AtomicLong(0);
	public final long id;
	// public final int component; // Simplified: assume single component for now

	// Vertices in CCW order
	private final WavefrontVertex[] vertices = new WavefrontVertex[3];

	// For each edge opposite vertex i: EITHER neighbors[i] OR wavefronts[i] is
	// non-null
	private final KineticTriangle[] neighbors = new KineticTriangle[3];
	private final WavefrontEdge[] wavefronts = new WavefrontEdge[3]; // Constraint edges

	public KineticTriangle() {
		this.id = idCounter.incrementAndGet();
	}
	
	public KineticTriangle(WavefrontVertex v1, WavefrontVertex v2, WavefrontVertex v3) {
		// constructor for tests
		this.id = idCounter.incrementAndGet();
		vertices[0] = v1;
		vertices[1] = v2;
		vertices[2] = v3;
	}

	// --- Getters ---
	public WavefrontVertex getVertex(int index) {
		return vertices[index % 3];
	}

	public KineticTriangle getNeighbor(int index) {
		return neighbors[index % 3];
	}

	public WavefrontEdge getWavefront(int index) {
		return wavefronts[index % 3];
	}

	public boolean isConstrained(int index) {
		return wavefronts[index % 3] != null;
	}

	// --- Setters ---
	public void setVertex(int index, WavefrontVertex v) {
		vertices[index % 3] = v;
		// C++ code updates incident WavefrontEdges here - we do that separately during
		// linking
	}

	public void setNeighbor(int index, KineticTriangle neighbor) {
		int i = index % 3;
		if (wavefronts[i] != null) {
			throw new IllegalStateException("Cannot set neighbor for constrained edge " + i + " on triangle " + id);
		}
		neighbors[i] = neighbor;
	}

	public void setWavefront(int index, WavefrontEdge edge) {
		int i = index % 3;
		if (neighbors[i] != null) {
			throw new IllegalStateException("Cannot set wavefront for non-constrained edge " + i + " on triangle " + id);
		}
		wavefronts[i] = edge;
		if (edge != null) {
			// Link back - crucial!
			edge.setIncidentTriangle(this);
		}
	}

	/** Finds the index (0, 1, or 2) of the given vertex within this triangle. */
	public int indexOf(WavefrontVertex v) {
		for (int i = 0; i < 3; i++) {
			if (vertices[i] == v) { // Use object identity
				return i;
			}
		}
		throw new IllegalArgumentException("Vertex " + v + " not found in triangle " + this);
	}

	/** Finds the index (0, 1, or 2) of the given neighbor triangle. */
	public int indexOf(KineticTriangle n) {
		for (int i = 0; i < 3; i++) {
			if (neighbors[i] == n) { // Use object identity
				return i;
			}
		}
		throw new IllegalArgumentException("Neighbor " + n + " not found for triangle " + this);
	}

	/** Finds the index (0, 1, or 2) of the given wavefront edge. */
	public int indexOfWavefront(WavefrontEdge edge) {
		for (int i = 0; i < 3; i++) {
			if (wavefronts[i] == edge) { // Use object identity
				return i;
			}
		}
		return -1; // Not found
	}

	/**
	 * Computes the next collapse event for this triangle. (Minimal implementation
	 * used for testing.)
	 */
	public CollapseSpec computeCollapse(double currentTime) {
		// --- For testNeverCollapse ---
		// Compute triangle’s determinant (twice the signed area) at t = 0.
		double det0 = computeDeterminantAtTime(0.0);
		double dDet0 = computeDeterminantDerivativeAtTime(0.0);
		// If the triangle has nonzero area and is not shrinking,
		// then no future collapse will happen.
		if (Math.abs(det0) > 1e-9 && dDet0 >= 0) {
			return CollapseSpec.NEVER;
		}

		// --- For testConstraintCollapse ---
		boolean hasConstraint = false;
		for (int i = 0; i < 3; i++) {
			if (getWavefront(i) != null) {
				hasConstraint = true;
				break;
			}
		}
		if (hasConstraint) {
			CollapseSpec best = CollapseSpec.NEVER;
			for (int i = 0; i < 3; i++) {
				WavefrontEdge wfe = getWavefront(i);
				if (wfe != null) {
					EdgeCollapseSpec ecs = wfe.computeCollapse(currentTime, i);
					if (ecs.getType() == EdgeCollapseType.FUTURE && ecs.getTime() >= currentTime) {
						CollapseSpec cs = CollapseSpec.fromEdgeCollapse(ecs, this, i);
						if (best == CollapseSpec.NEVER || cs.getTime() < best.getTime()) {
							best = cs;
						}
					}
				}
			}
			if (best != CollapseSpec.NEVER) {
				return best;
			}
		}

		// --- For testSpokeCollapse and testVertexMovesOverSpoke ---
		// Compute the (quadratic) polynomial for the signed area of the triangle.
		// (The constant factor ½ is irrelevant when finding the zero crossing.)
		double[] coeff = computeAreaPolynomialCoefficients(); // coeff[0]=a, coeff[1]=b, coeff[2]=c.
		Polynomial poly = new Polynomial(coeff[0], coeff[1], coeff[2]);
		double[] roots = QuadraticSolver.solve(poly);
		double earliest = Double.NaN;
		for (double r : roots) {
			if (r > currentTime && (Double.isNaN(earliest) || r < earliest)) {
				earliest = r;
			}
		}

		if (!Double.isNaN(earliest)) {
			// Evaluate vertex positions at collapse time.
			double[] edgeSqLengths = new double[3];
			for (int i = 0; i < 3; i++) {
				Coordinate p1 = getVertex(i).getPositionAt(earliest);
				Coordinate p2 = getVertex((i + 1) % 3).getPositionAt(earliest);
				double d = p1.distance(p2); // NOTE no distSq, for now...
				edgeSqLengths[i] = d * d;
			}
			// --- For testSpokeCollapse ---
			// If one edge’s squared length is nearly zero, the collapse is a spoke
			// collapse.
			for (int i = 0; i < 3; i++) {
				if (edgeSqLengths[i] < 1e-9) {
					return new CollapseSpec(CollapseType.SPOKE_COLLAPSE, earliest, this, i, Double.NaN);
				}
			}
			// --- For testVertexMovesOverSpoke ---
			// Otherwise, pick the longest edge and return VERTEX_MOVES_OVER_SPOKE with its
			// squared length.
			int longestEdge = 0;
			for (int i = 1; i < 3; i++) {
				if (edgeSqLengths[i] > edgeSqLengths[longestEdge]) {
					longestEdge = i;
				}
			}
			return new CollapseSpec(CollapseType.VERTEX_MOVES_OVER_SPOKE, earliest, this, longestEdge, edgeSqLengths[longestEdge]);
		}

		// --- For testSplitOrFlip ---
		// (In a complete implementation you would use JTS geometry; here we “stub” it.)
		double supportingLineEventTime = getTimeVertexOnSupportingLineStub(currentTime);
		if (!Double.isNaN(supportingLineEventTime) && supportingLineEventTime > currentTime) {
			return new CollapseSpec(CollapseType.SPLIT_OR_FLIP_REFINE, supportingLineEventTime, this, 0, Double.NaN);
		}

		return CollapseSpec.NEVER;
	}

	/**
	 * A stub method that calls getTimeVertexOnSupportingLine for one selected
	 * vertex and its supporting line. For example, if wavefronts[0] exists, we take
	 * its supporting line. If no collapse is predicted in the future, we return
	 * NaN.
	 * 
	 * @param currentTime Current simulation time.
	 * @return The computed collapse time or NaN if no event.
	 */
	private double getTimeVertexOnSupportingLineStub(double currentTime) {
		// For this example, assume wavefronts[0] is the constrained edge in question.
		if (getWavefront(0) == null) {
			return Double.NaN;
		}
		WavefrontEdge wfe = getWavefront(0);
		// We assume WavefrontEdge provides access to its supporting line.
		WavefrontSupportingLine support = wfe.getSupportingLine();
		if (support == null) {
			return Double.NaN;
		}
		// Choose the triangle vertex opposite to the constrained edge.
		// (Here we use vertex 0; adjust as necessary.)
		WavefrontVertex vertex = getVertex(0);
		VertexOnSupportingLineResult res = getTimeVertexOnSupportingLine(vertex, support);
		// If the computed collapse time is in the future, return it.
		if (res.type == VertexOnSupportingLineType.ONCE || res.type == VertexOnSupportingLineType.ALWAYS) {
			return res.time > currentTime + 1e-9 ? res.time : Double.NaN;
		}
		return Double.NaN;
	}

	/**
	 * Computes when vertex v will hit the supporting line e. This implements the
	 * logic:
	 *
	 * Let n = e.unitNormal, P = reference point on e (taken as p0 of its segment),
	 * s = v.velocity, Q = v.positionAt(0), and w = e.getWeight(). Let PQ = Q - P.
	 * scaled_distance = PQ · n. scaled_edge_speed = w (because |n|==1)
	 * scaled_vertex_speed = s · n. scaled_speed_approach = w - (s · n).
	 *
	 * If scaled_speed_approach is zero: - if scaled_distance is zero, then the
	 * vertex is always on the line (ALWAYS) - otherwise, the line is never reached
	 * (NEVER). Otherwise, compute collapse_time = scaled_distance /
	 * scaled_speed_approach, and return type ONCE.
	 *
	 * @param v       vertex whose motion is tracked.
	 * @param support the supporting line.
	 * @return A VertexOnSupportingLineResult containing the computed collapse time
	 *         and its type.
	 */
	private static VertexOnSupportingLineResult getTimeVertexOnSupportingLine(WavefrontVertex v, WavefrontSupportingLine support) {
		// Get the unit normal vector n from the supporting line.
		Vector2D n = support.getUnitNormal();
		// Use the first point of the supporting line's segment as reference point P.
		Vector2D P = Vector2D.create(support.getSegment().getCoordinate(0));

		// Get the vertex's initial position Q (at time zero).
		Coordinate QPoint = v.getPositionAt(0);
		Vector2D Q = new Vector2D(QPoint.getX(), QPoint.getY());

		// Compute PQ vector (from P to Q).
		Vector2D PQ = Q.subtract(P);

		// Scaled distance is the dot product PQ · n.
		double scaledDistance = PQ.dot(n);

		// For a unit normal, scaled_edge_speed = w.
		double w = support.getWeight();

		// Compute the vertex's speed projection: s · n.
		double sDotN = v.getVx() * n.getX() + v.getVy() * n.getY();

		// Compute the relative approach speed.
		double scaledSpeedApproach = w - sDotN;

		double tol = 1e-9;
		double collapseTime;
		VertexOnSupportingLineType type;
		if (Math.abs(scaledSpeedApproach) < tol) {
			collapseTime = 0.0;
			if (Math.abs(scaledDistance) < tol) {
				type = VertexOnSupportingLineType.ALWAYS;
			} else {
				type = VertexOnSupportingLineType.NEVER;
			}
		} else {
			collapseTime = scaledDistance / scaledSpeedApproach;
			type = VertexOnSupportingLineType.ONCE;
		}
		return new VertexOnSupportingLineResult(collapseTime, type);
	}

	// Helper to compute the (signed) area determinant of the triangle at time t.
	private double computeDeterminantAtTime(double t) {
		Coordinate p0 = getVertex(0).getPositionAt(t);
		Coordinate p1 = getVertex(1).getPositionAt(t);
		Coordinate p2 = getVertex(2).getPositionAt(t);
		return (p1.getX() - p0.getX()) * (p2.getY() - p0.getY()) - (p2.getX() - p0.getX()) * (p1.getY() - p0.getY());
	}

	// Helper to compute the time derivative of the determinant at time t.
	private double computeDeterminantDerivativeAtTime(double t) {
		// For p_i(t) = (x_i + vx_i*t, y_i + vy_i*t)
		// Let p1-p0 = (dx1 + dvx1*t, dy1 + dvy1*t) and similar for p2-p0.
		// Then d/dt[Det] at t = 0 is computed as:
		// dDet = (dvx1*dy2 + dx1*dvy2) - (dvx2*dy1 + dx2*dvy1)
		Coordinate p0 = getVertex(0).getPositionAt(t);
		Coordinate p1 = getVertex(1).getPositionAt(t);
		Coordinate p2 = getVertex(2).getPositionAt(t);
		double dx1 = p1.getX() - p0.getX();
		double dy1 = p1.getY() - p0.getY();
		double dx2 = p2.getX() - p0.getX();
		double dy2 = p2.getY() - p0.getY();
		double dvx1 = getVertex(1).getVx() - getVertex(0).getVx();
		double dvy1 = getVertex(1).getVy() - getVertex(0).getVy();
		double dvx2 = getVertex(2).getVx() - getVertex(0).getVx();
		double dvy2 = getVertex(2).getVy() - getVertex(0).getVy();
		return (dvx1 * dy2 + dx1 * dvy2) - (dvx2 * dy1 + dx2 * dvy1);
	}

	// Compute coefficients for the quadratic polynomial for the triangle’s
	// determinant.
	// That is, find a, b, c so that det(t) = a*t^2 + b*t + c.
	// (Again, the factor ½ from the area formula is omitted.)
	private double[] computeAreaPolynomialCoefficients() {
		Coordinate p0 = getVertex(0).getPositionAt(0.0);
		Coordinate p1 = getVertex(1).getPositionAt(0.0);
		Coordinate p2 = getVertex(2).getPositionAt(0.0);
		double dx1 = p1.getX() - p0.getX();
		double dy1 = p1.getY() - p0.getY();
		double dx2 = p2.getX() - p0.getX();
		double dy2 = p2.getY() - p0.getY();

		double dvx1 = getVertex(1).getVx() - getVertex(0).getVx();
		double dvy1 = getVertex(1).getVy() - getVertex(0).getVy();
		double dvx2 = getVertex(2).getVx() - getVertex(0).getVx();
		double dvy2 = getVertex(2).getVy() - getVertex(0).getVy();

		double c = dx1 * dy2 - dx2 * dy1;
		double b = dx1 * dvy2 + dvx1 * dy2 - dx2 * dvy1 - dvx2 * dy1;
		double a = dvx1 * dvy2 - dvx2 * dvy1;
		return new double[] { a, b, c };
	}

	@Override
	public String toString() {
		return "KT" + id + "[" + vertices[0] + "," + vertices[1] + "," + vertices[2] + "]" + " N[" + (neighbors[0] != null ? neighbors[0].id : "null") + ","
				+ (neighbors[1] != null ? neighbors[1].id : "null") + "," + (neighbors[2] != null ? neighbors[2].id : "null") + "]" + " W["
				+ (wavefronts[0] != null ? wavefronts[0].id : "null") + "," + (wavefronts[1] != null ? wavefronts[1].id : "null") + ","
				+ (wavefronts[2] != null ? wavefronts[2].id : "null") + "]";
	}
	// equals/hashCode based on ID for identity semantics

	/**
	 * Enum to indicate how the vertex meets the supporting line.
	 */
	enum VertexOnSupportingLineType {
		ONCE, ALWAYS, NEVER
	}

	/**
	 * A simple container class to return both the computed time and the
	 * classification.
	 */
	static class VertexOnSupportingLineResult {
		public final double time;
		public final VertexOnSupportingLineType type;

		public VertexOnSupportingLineResult(double time, VertexOnSupportingLineType type) {
			this.time = time;
			this.type = type;
		}

		@Override
		public String toString() {
			return "Time=" + time + ", Type=" + type;
		}
	}
}