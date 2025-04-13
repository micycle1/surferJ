package com.github.micycle.surferj2.kinetics;

import java.util.concurrent.atomic.AtomicLong;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;

import com.github.micycle.surferj2.collapse.EdgeCollapseSpec;
import com.github.micycle.surferj2.collapse.EdgeCollapseType;

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

	/**
	 * Compute collapse time for this constraint edge. For now we determine the
	 * current length and a “relative speed.” If the edge is shrinking (i.e. length
	 * > 0 and relative speed > 0) we return a FUTURE collapse event.
	 */
	public EdgeCollapseSpec computeCollapse(double currentTime, int edgeIndex) {
		Coordinate p0 = vertex0.getPositionAt(currentTime);
		Coordinate p1 = vertex1.getPositionAt(currentTime);
		double length = p0.distance(p1);
		// approximate relative speed along the edge:
		double dvx = vertex1.getVx() - vertex0.getVx();
		double dvy = vertex1.getVy() - vertex0.getVy();
		double relSpeed = Math.hypot(dvx, dvy); // NOTE hypot (slower)
		if (length < 1e-9) {
			return new EdgeCollapseSpec(EdgeCollapseType.ALWAYS, currentTime);
		}
		if (relSpeed < 1e-9) {
			return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
		}
		double dt = length / relSpeed; // simplistic estimate
		double collapseTime = currentTime + dt;
		return new EdgeCollapseSpec(EdgeCollapseType.FUTURE, collapseTime);
	}


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

	@Override
	public String toString() {
		return "WE" + id + "[" + vertex0 + "->" + vertex1 + "]";
	}
	// equals/hashCode based on ID for identity semantics
}