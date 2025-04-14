package com.github.micycle1.surferj.kinetics;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.math.Vector2D;

import com.github.micycle1.surferj.SurfConstants;

/**
 * When are WavefrontSupportingLine created in C++?
 * 
 * During KineticTriangulation::create_supporting_lines: This is the primary
 * creation point. When the code iterates through the faces of the initial
 * BasicTriangulation (the CDT result) and identifies a constrained edge (i.e.,
 * an edge corresponding to an original input polygon segment), it creates a
 * WavefrontEdge. Inside the WavefrontEdge constructor (the one taking Point_2
 * u, Point_2 v, NT weight), a std::make_shared<const
 * WavefrontSupportingLine>(u, v, weight) is created. This
 * WavefrontSupportingLine object captures the infinite line containing the
 * segment (u, v) and its associated weight.
 * 
 * During KineticTriangulation::create_bevels_at_vertex (Potentially): When
 * handling degree-1 vertices (and potentially other reflex vertices requiring
 * beveling, though that part is complex and we're ignoring it for now), the
 * code might create new, temporary WavefrontEdges that don't correspond
 * directly to input segments. These also get a WavefrontSupportingLine created
 * in their respective WavefrontEdge constructor (the one taking a pre-existing
 * shared_ptr<const WavefrontSupportingLine>). For the degree-1 case mentioned
 * in create_bevels_at_vertex, it explicitly creates a new supporting line
 * perpendicular to the incident edge. During WavefrontEdge::split: When a
 * WavefrontEdge is split during the simulation (due to a split event), the two
 * new WavefrontEdge objects created reuse the same shared_ptr<const
 * WavefrontSupportingLine> as the original edge they came from. The underlying
 * supporting line (and its motion rule) doesn't change, only the endpoints
 * tracked by the WavefrontEdge and its WavefrontVertex objects change.
 * 
 * 
 * Relationship with WavefrontEdge:
 * 
 * Yes, the other LLM is correct. In the C++ code, WavefrontEdge contains a
 * shared_ptr<const WavefrontSupportingLine> (accessed via the l() method). This
 * is the standard pattern.
 * 
 * Your Java WavefrontEdge class absolutely needs to contain a
 * WavefrontSupportingLine member. It should be created within the WavefrontEdge
 * constructor using the initial segment data and weight.
 */
public class WavefrontSupportingLine {

	private final LineSegment segment; // Initial segment at t=0
	private final double weight;
	private final Vector2D normalDirection; /* arbitrary length, perpendicular to line direction */
	private final Vector2D normalUnit; /* unit length */
	private final Vector2D normal; // WEIGHTED unit-normal

	public WavefrontSupportingLine(double p1x, double p1y, double p2x, double p2y) {
		this(new LineSegment(p1x, p1y, p2x, p2y), 1.0); // NOTE default weight of 1.0
	}

	public WavefrontSupportingLine(LineSegment segment, double weight) {
		this.segment = segment;
		this.weight = weight;

		// Compute direction vector of the line (from p0 to p1)
		Vector2D directionVector = new Vector2D(segment.p0, segment.p1);
		if (directionVector.lengthSquared() < SurfConstants.ZERO_DIST_SQ) {
			throw new IllegalArgumentException("Cannot create supporting line for zero-length segment: " + segment);
		}

		// Calculate normal direction (rotate direction vector 90 degrees
		// counterclockwise)
		double dx = directionVector.getX();
		double dy = directionVector.getY();
		this.normalDirection = new Vector2D(-dy, dx);

		// Calculate unit normal
		this.normalUnit = this.normalDirection.normalize();

		// Calculate weighted normal
		this.normal = this.normalUnit.multiply(this.weight);
	}

	public Coordinate lineIntersection(WavefrontSupportingLine o) {
		return this.segment.lineIntersection(o.segment);
	}

	public LineSegment getSegment() {
		return segment;
	}

	public double getWeight() {
		return weight;
	}

	public Vector2D getNormalDirection() {
		return normalDirection;
	}

	/**
	 * Unit-length normal.
	 */
	public Vector2D getUnitNormal() {
		return normalUnit;
	}

	/**
	 * Weighted.
	 */
	public Vector2D getNormal() {
		return normal;
	}

	/**
	 * The wavefront line after it has moved for time t=1.
	 */
	public LineSegment lineAtOne() {
		return getOffsetSegment(1);
	}

	/**
	 * Gets the position of the supporting line offset by time t. The offset is
	 * applied along the weighted normal direction.
	 * 
	 * @param t time
	 * @return The offset line segment (endpoints translated by normal * t)
	 */
	private LineSegment getOffsetSegment(double t) {
		Vector2D offset = normal.multiply(t);
		Coordinate p0Offset = Vector2D.create(segment.p0).add(offset).toCoordinate();
		Coordinate p1Offset = Vector2D.create(segment.p1).add(offset).toCoordinate();
		return new LineSegment(p0Offset, p1Offset);
	}

	@Override
	public String toString() {
		return "WSL{" + segment + ", w=" + weight + ", n=" + normalUnit + '}';
	}

	public boolean isVertical() {
		return segment.isVertical();
	}

	public boolean isVertical(double epsilon) {
		return Math.abs(segment.p0.x - segment.p1.x) < epsilon;
	}
}