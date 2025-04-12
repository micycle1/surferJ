package com.github.micycle.surferj2.kinetics;

import java.util.concurrent.atomic.AtomicLong;

import org.locationtech.jts.geom.Coordinate;

// --- Wavefront Vertex ---
public class WavefrontVertex {
	private static final AtomicLong idCounter = new AtomicLong(0);
	public final long id;
	public final Coordinate initialPosition; // Corresponds to pos_zero/pos_start
	// Velocity, time_start, time_stop etc. added later for kinetic part
	public final boolean isInfinite;

	// References to the two incident wavefront edges forming this vertex
	// Order matters: edge0 is CCW, edge1 is CW when looking from inside the
	// skeleton towards the vertex
	private WavefrontEdge edge0; // Left/CCW in C++ incident_wavefront_edges[0] (a)
	private WavefrontEdge edge1; // Right/CW in C++ incident_wavefront_edges[1] (b)

	// Placeholder for the infinite vertex
	public static final WavefrontVertex INFINITE_VERTEX = new WavefrontVertex();

	private WavefrontVertex() {
		this.id = -1; // Special ID for infinite
		this.initialPosition = null;
		this.isInfinite = true;
	}

	public WavefrontVertex(Coordinate initialPosition) {
		this.id = idCounter.incrementAndGet();
		this.initialPosition = initialPosition;
		this.isInfinite = false;
	}

	public void setIncidentEdges(WavefrontEdge edge0, WavefrontEdge edge1) {
		if (this.isInfinite)
			throw new IllegalStateException("Cannot set edges for infinite vertex");
		// Basic check: the edges should share this vertex's coordinate
		Coordinate sharedCoord = null;
		if (edge0 != null && edge1 != null) {
			if (edge0.getVertex(0) == this)
				sharedCoord = edge0.getVertex(0).initialPosition;
			else if (edge0.getVertex(1) == this)
				sharedCoord = edge0.getVertex(1).initialPosition;

			if (sharedCoord == null || !sharedCoord.equals(this.initialPosition)) {
				throw new IllegalArgumentException(
						"Edges do not consistently share this vertex coordinate " + this.initialPosition + " edge0: " + edge0 + " edge1: " + edge1);
			}
			if ((edge1.getVertex(0) != this && edge1.getVertex(1) != this)) {
				throw new IllegalArgumentException("Edge1 does not contain this vertex coordinate " + this.initialPosition + " edge1: " + edge1);
			}
		}
		// We might set them one by one during initialization
		this.edge0 = edge0;
		this.edge1 = edge1;
	}

	public WavefrontEdge getIncidentEdge(int index) {
		if (index == 0)
			return edge0;
		if (index == 1)
			return edge1;
		throw new IndexOutOfBoundsException("Index must be 0 or 1");
	}

	public void setIncidentEdge(int index, WavefrontEdge edge) {
		if (isInfinite)
			throw new IllegalStateException("Cannot set edges for infinite vertex");
		if (index == 0)
			this.edge0 = edge;
		else if (index == 1)
			this.edge1 = edge;
		else
			throw new IndexOutOfBoundsException("Index must be 0 or 1");
	}

	@Override
	public String toString() {
		if (isInfinite)
			return "WV(Inf)";
		return "WV" + id + "(" + initialPosition.getX() + "," + initialPosition.getY() + ")";
	}

	// equals/hashCode based on ID for identity semantics if needed,
	// but using object identity is often sufficient.
	// Coordinate equality is based on value.
}