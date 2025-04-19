package com.github.micycle1.surferj.kinetics;

public class AroundVertexIterator {
	
	// TODO swap for TriangulationUtils approach?
	private static final int[] CW = { 2, 0, 1 };
	private static final int[] CCW = { 1, 2, 0 };

	private KineticTriangle t;
	private int vInTIdx;

	/** Create an iterator at triangle t, pointing into corner vInTIdx. */
	public AroundVertexIterator(KineticTriangle t, int vInTIdx) {
		this.t = t;
		this.vInTIdx = vInTIdx;
	}

	/** The triangle we are currently on (or null if end). */
	public KineticTriangle t() {
		return t;
	}

	/** The corner index in t (0,1,2) that we are iterating around. */
	public int vInTIdx() {
		return vInTIdx;
	}

	/** True if we have walked off the mesh (end of fan). */
	public boolean isEnd() {
		return t == null;
	}

	/** Return the neighbor across edge = direction[vInTIdx], or null. */
	KineticTriangle nextTriangle(int[] direction) {
		if (t == null)
			return null;
		return t.getNeighbor(direction[vInTIdx]);
	}

	/** Move one step CW (or CCW) around the vertex. */
	private AroundVertexIterator walkDir(int[] direction) {
		if (t == null) {
			// already at end
			return this;
		}
		KineticTriangle next = nextTriangle(direction);
		int newVInNext = 0;
		if (next != null) {
			// find index of current t in next.neighbors
			int idxInNext = next.indexOfNeighbor(t);
			// then apply direction[] to compute the new corner
			newVInNext = direction[idxInNext];
		}
		this.t = next;
		this.vInTIdx = newVInNext;
		return this;
	}

	/** Walk one step clockwise around the shared vertex. */
	public AroundVertexIterator walkCw() {
		return walkDir(CW);
	}

	/** Walk one step counter‐clockwise around the shared vertex. */
	public AroundVertexIterator walkCcw() {
		return walkDir(CCW);
	}

	/** Return an iterator advanced CW until the last non‐null triangle. */
	public AroundVertexIterator mostCw() {
		AroundVertexIterator res = new AroundVertexIterator(t, vInTIdx);
		while (res.nextTriangle(CW) != null) {
			res.walkDir(CW);
		}
		return res;
	}

	/** Return an iterator advanced CCW until the last non‐null triangle. */
	public AroundVertexIterator mostCcw() {
		AroundVertexIterator res = new AroundVertexIterator(t, vInTIdx);
		while (res.nextTriangle(CCW) != null) {
			res.walkDir(CCW);
		}
		return res;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!(obj instanceof AroundVertexIterator))
			return false;
		AroundVertexIterator o = (AroundVertexIterator) obj;
		return this.t == o.t && this.vInTIdx == o.vInTIdx;
	}

	@Override
	public int hashCode() {
		int result = 17;
		result = 31 * result + System.identityHashCode(t);
		result = 31 * result + vInTIdx;
		return result;
	}

	@Override
	public String toString() {
		if (t == null) {
			return "AVI(end)";
		}
		return "AVI(T" + t.getId() + ", corner=" + vInTIdx + ")";
	}

	public AroundVertexIterator incidentFacesIterator(KineticTriangle t, int vInT) {
		return new AroundVertexIterator(t, vInT);
	}

	public AroundVertexIterator incidentFacesEnd() {
		return new AroundVertexIterator(null, -1);
	}
}