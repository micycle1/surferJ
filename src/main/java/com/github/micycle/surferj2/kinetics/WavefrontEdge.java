package com.github.micycle.surferj2.kinetics;

import java.util.concurrent.atomic.AtomicLong;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;

public class WavefrontEdge {
	
	private static final AtomicLong idCounter = new AtomicLong(0);
	public final long id;
	public CanonicalSegment canonicalSegment; // Original segment
	public final double weight; // From input, default 1.0
    public final WavefrontSupportingLine supportingLine;


	// Endpoints: vertex0 is start, vertex1 is end relative to KineticTriangle's
	// perspective
	private WavefrontVertex vertex0; // Corresponds to C++ vertices[0] (left)
	private WavefrontVertex vertex1; // Corresponds to C++ vertices[1] (right)

	// The triangle this edge forms a boundary for
	private KineticTriangle incidentTriangle;
	
    public WavefrontEdge(WavefrontSupportingLine supportingLine) {
        this.id = idCounter.incrementAndGet();
        this.weight = supportingLine.getWeight();
        this.supportingLine = supportingLine;
        // Validation: Check if supportingLine is valid? (e.g., non-zero length)
        if (supportingLine == null) {
             throw new IllegalArgumentException("WavefrontSupportingLine cannot be null");
        }
        this.canonicalSegment = new CanonicalSegment(supportingLine.getSegment());
    }

//	public WavefrontEdge(CanonicalSegment segment, double weight) {
//		this.id = idCounter.incrementAndGet();
//		this.canonicalSegment = segment;
//		this.weight = weight;
//	}
	
	 public WavefrontEdge(LineSegment segment, double weight) {
		 this(new WavefrontSupportingLine(segment, weight));
	    }
	 
		public WavefrontSupportingLine getSupportingLine() {
			return supportingLine;
		}

		public LineSegment getSegment() {
			return supportingLine.getSegment();
		}

		public double getWeight() {
			return supportingLine.getWeight();
		}

		public void setVertices(WavefrontVertex v0, WavefrontVertex v1) {
			// Ensure the vertices match the segment coordinates conceptually
			if (v0 != null && v1 != null && !v0.isInfinite && !v1.isInfinite) {
				Coordinate c0 = getSegment().p0;
				Coordinate c1 = getSegment().p1;
				// Check if vertices align with segment ends, allowing for swapped order
				boolean match1 = (v0.initialPosition.equals2D(c0) && v1.initialPosition.equals2D(c1));
				boolean match2 = (v0.initialPosition.equals2D(c1) && v1.initialPosition.equals2D(c0));
				if (!match1 && !match2) {
					throw new IllegalArgumentException("Vertices " + v0 + ", " + v1 + " do not match edge segment ends " + getSegment());
				}
			} else if (v0 == null || v1 == null) {
				throw new IllegalArgumentException("Cannot set null vertices for WavefrontEdge");
			}
			this.vertex0 = v0;
			this.vertex1 = v1;
		}

		public void setVertex(int index, WavefrontVertex v) {
			if (index == 0)
				this.vertex0 = v;
			else if (index == 1)
				this.vertex1 = v;
			else
				throw new IndexOutOfBoundsException("Index must be 0 or 1");
		}

	public WavefrontVertex getVertex(int index) {
		if (index == 0)
			return vertex0;
		if (index == 1)
			return vertex1;
		throw new IndexOutOfBoundsException("Index must be 0 or 1");
	}

	public void setIncidentTriangle(KineticTriangle triangle) {
		this.incidentTriangle = triangle;
	}

	public KineticTriangle getIncidentTriangle() {
		return incidentTriangle;
	}
	
    /**
     * Computes when this edge collapses to zero length.
     * Stub implementation for TDD.
     *
     * @param currentTime The current time in the simulation.
     * @return The EdgeCollapseSpec describing the edge collapse.
     */
//	public EdgeCollapseSpec computeCollapse(double currentTime) {
//        double x0 = vertex0.getX0(), y0 = vertex0.getY0();
//        double vx0 = vertex0.getVx(), vy0 = vertex0.getVy();
//        double x1 = vertex1.getX0(), y1 = vertex1.getY0();
//        double vx1 = vertex1.getVx(), vy1 = vertex1.getVy();
//
//        double dx = x0 - x1, dy = y0 - y1;
//        double dvx = vx0 - vx1, dvy = vy0 - vy1;
//
//        if (Math.abs(dvx) < 1e-12 && Math.abs(dvy) < 1e-12) {
//            return (Math.abs(dx) < 1e-12 && Math.abs(dy) < 1e-12) ?
//                new EdgeCollapseSpec(EdgeCollapseType.FUTURE, currentTime) :
//                new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
//        }
//
//        double t_x = (Math.abs(dvx) > 1e-12) ? -dx / dvx : Double.POSITIVE_INFINITY;
//        double t_y = (Math.abs(dvy) > 1e-12) ? -dy / dvy : Double.POSITIVE_INFINITY;
//
//        if (Math.abs(t_x - t_y) < 1e-9 || (Double.isInfinite(t_y) && !Double.isInfinite(t_x))) {
//            double t = t_x;
//            if (t > currentTime) return new EdgeCollapseSpec(EdgeCollapseType.FUTURE, t);
//        } else if (Double.isInfinite(t_x) && !Double.isInfinite(t_y)) {
//            double t = t_y;
//            if (t > currentTime) return new EdgeCollapseSpec(EdgeCollapseType.FUTURE, t);
//        }
//        return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
//    }

	@Override
	public String toString() {
		return "WE" + id + "[" + vertex0 + "->" + vertex1 + "]";
	}
	// equals/hashCode based on ID for identity semantics
}