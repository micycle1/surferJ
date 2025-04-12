package com.github.micycle.surferj2.kinetics;

import java.util.concurrent.atomic.AtomicLong;

import org.locationtech.jts.geom.Coordinate;

public class WavefrontEdge {
	
	private static final AtomicLong idCounter = new AtomicLong(0);
	public final long id;
	public final CanonicalSegment canonicalSegment; // Original segment
	public final double weight; // From input, default 1.0

	// Endpoints: vertex0 is start, vertex1 is end relative to KineticTriangle's
	// perspective
	private WavefrontVertex vertex0; // Corresponds to C++ vertices[0] (left)
	private WavefrontVertex vertex1; // Corresponds to C++ vertices[1] (right)

	// The triangle this edge forms a boundary for
	private KineticTriangle incidentTriangle;

	public WavefrontEdge(CanonicalSegment segment, double weight) {
		this.id = idCounter.incrementAndGet();
		this.canonicalSegment = segment;
		this.weight = weight;
	}

	public void setVertices(WavefrontVertex v0, WavefrontVertex v1) {
		// Ensure the vertices match the segment coordinates conceptually
		if (v0 != null && v1 != null) {
			Coordinate c0 = canonicalSegment.getP0();
			Coordinate c1 = canonicalSegment.getP1();
			// Check if vertices align with segment ends, allowing for swapped order
			boolean match1 = (v0.initialPosition.equals(c0) && v1.initialPosition.equals(c1));
			boolean match2 = (v0.initialPosition.equals(c1) && v1.initialPosition.equals(c0));
			if (!match1 && !match2) {
				throw new IllegalArgumentException("Vertices " + v0 + ", " + v1 + " do not match segment ends " + canonicalSegment);
			}
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