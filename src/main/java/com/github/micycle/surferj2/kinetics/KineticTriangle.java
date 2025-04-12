package com.github.micycle.surferj2.kinetics;

import java.util.concurrent.atomic.AtomicLong;

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
		// this.component = component;
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

	@Override
	public String toString() {
		return "KT" + id + "[" + vertices[0] + "," + vertices[1] + "," + vertices[2] + "]" + " N[" + (neighbors[0] != null ? neighbors[0].id : "null") + ","
				+ (neighbors[1] != null ? neighbors[1].id : "null") + "," + (neighbors[2] != null ? neighbors[2].id : "null") + "]" + " W["
				+ (wavefronts[0] != null ? wavefronts[0].id : "null") + "," + (wavefronts[1] != null ? wavefronts[1].id : "null") + ","
				+ (wavefronts[2] != null ? wavefronts[2].id : "null") + "]";
	}
	// equals/hashCode based on ID for identity semantics
}