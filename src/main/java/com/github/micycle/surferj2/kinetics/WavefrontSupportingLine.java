package com.github.micycle.surferj2.kinetics;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.math.Vector2D;

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
	// Calculated properties
	private final Vector2D directionVector; // Vector along the segment
	private final Vector2D unitNormal; // Unit normal vector (consistent side, e.g., "left" relative to direction)
	private final Vector2D weightedNormal; // unitNormal * weight
	
	public WavefrontSupportingLine(double p1x, double p1y, double p2x, double p2y) {
		this(new LineSegment(p1x, p1y, p2x, p2y), 1);
	}

	public WavefrontSupportingLine(LineSegment segment, double weight) {
		this.segment = segment; // Store the original segment
		this.weight = weight;

		// Calculate direction (ensure consistent orientation if needed, e.g., p0 to p1)
		this.directionVector = new Vector2D(segment.p0, segment.p1);
		if (this.directionVector.length() < 1e-12) {
			throw new IllegalArgumentException("Cannot create supporting line for zero-length segment: " + segment);
		}

		// Calculate unit normal (e.g., rotate direction vector -90 degrees for "left"
		// normal)
		// JTS doesn't have a direct Vector2D normal, calculate manually:
		// Normal to (dx, dy) is (-dy, dx) or (dy, -dx)
		double dx = directionVector.getX();
		double dy = directionVector.getY();
		Vector2D normal = new Vector2D(-dy, dx); // Left-hand normal relative to segment direction p0->p1
		this.unitNormal = normal.normalize(); // Make it unit length

		// Calculate weighted normal
		this.weightedNormal = this.unitNormal.multiply(this.weight);
		
	}
	
	public Coordinate lineIntersection(WavefrontSupportingLine o) {
		return this.segment.lineIntersection(o.segment); // robust -- CGAlgorithmsDD
	}

	// Getters
	public LineSegment getSegment() {
		return segment;
	}

	public double getWeight() {
		return weight;
	}
	
	public Vector2D getDirection() {
		return directionVector;
	}

	public Vector2D getUnitNormal() {
		return unitNormal;
	}

	public Vector2D getWeightedNormal() {
		return weightedNormal;
	}

	/**
	 * Gets the position of the supporting line offset by time t. Note: Calculating
	 * the infinite line representation might be more useful for intersections than
	 * just offsetting the segment endpoints.
	 * 
	 * @param t time
	 * @return The offset line segment (for visualization or basic checks)
	 */
	public LineSegment getOffsetSegment(double t) {
		Vector2D offset = weightedNormal.multiply(t);
		Coordinate p0Offset = Vector2D.create(segment.p0).add(offset).toCoordinate();
		Coordinate p1Offset = Vector2D.create(segment.p1).add(offset).toCoordinate();
		return new LineSegment(p0Offset, p1Offset);
	}

	// TODO: Potentially add methods to get an infinite Line representation (e.g.,
	// using parameters a,b,c for ax+by+c=0)
	// for more robust intersection calculations, especially for offset lines.

	@Override
	public String toString() {
		return "WSL{" + segment + ", w=" + weight + ", n=" + unitNormal + '}';
	}
}