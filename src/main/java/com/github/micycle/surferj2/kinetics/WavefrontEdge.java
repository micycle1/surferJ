package com.github.micycle.surferj2.kinetics;

import java.util.concurrent.atomic.AtomicLong;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;

import com.github.micycle.surferj2.collapse.CollapseSpec;
import com.github.micycle.surferj2.collapse.EdgeCollapseSpec;
import com.github.micycle.surferj2.collapse.EdgeCollapseType;

public class WavefrontEdge {

	private static final AtomicLong idCounter = new AtomicLong(0);
	public final long id;
	public CanonicalSegment canonicalSegment; // Original segment representation
	public final double weight; // Speed, from input, default 1.0
	public final WavefrontSupportingLine supportingLine; // Geometric and kinetic line info

	// Endpoints: vertex0 is CCW, vertex1 is CW relative to incidentTriangle edge
	private WavefrontVertex vertex0; // Corresponds to C++ vertices[0]
	private WavefrontVertex vertex1; // Corresponds to C++ vertices[1]

	// The triangle this edge forms a boundary for
	private KineticTriangle incidentTriangle;

	private boolean isDead = false;

	// Fields related to initial state and skeleton tracking
	public final boolean isInitial; // True for edges from input polygon
	public final boolean isBeveling; // True for beveling edges (zero length initial)
	private final WavefrontVertex[] initialVertices = new WavefrontVertex[2]; // For isInitial edges, store original vertices
	public final SkeletonDCELFace skeletonFace; // Straight skeleton face traced by this edge (if not beveling)

	// Collapse caching
	private EdgeCollapseSpec cachedCollapseSpec; // Cached collapse computation
	private boolean collapseSpecValid = false; // Cache validity flag

	// --- Constructors ---

	// Constructor for initial edges from input coordinates
	public WavefrontEdge(Coordinate u, Coordinate v, double weight, KineticTriangle incidentTriangle, SkeletonDCELFace skeletonFace) {
		this.id = idCounter.incrementAndGet();
		// Ensure segment direction matches typical CCW winding if possible, although
		// supporting line handles it
		this.supportingLine = new WavefrontSupportingLine(new LineSegment(u, v), weight);
		this.canonicalSegment = new CanonicalSegment(this.supportingLine.getSegment()); // Use canonical segment
		this.weight = weight; // Store weight directly too? Redundant but maybe convenient
		this.incidentTriangle = incidentTriangle;
		this.isInitial = true;
		this.isBeveling = u.equals2D(v); // Beveling if start/end points are the same
		this.skeletonFace = skeletonFace; // Skeleton face associated with this edge
		if (skeletonFace != null && isBeveling) {
			throw new IllegalStateException("skeletonFace must be null for beveling edges");
		}
		// Initial vertices are set later after vertex objects are created and assigned
	}

	// Constructor from a supporting line (used internally?)
	public WavefrontEdge(WavefrontSupportingLine supportingLine) {
		this.id = idCounter.incrementAndGet();
		if (supportingLine == null) {
			throw new IllegalArgumentException("WavefrontSupportingLine cannot be null");
		}
		this.supportingLine = supportingLine;
		this.weight = supportingLine.getWeight();
		this.canonicalSegment = new CanonicalSegment(supportingLine.getSegment());
		LineSegment segment = supportingLine.getSegment();
		this.isInitial = false; // Typically not initial if created this way? Or depends on context? Assume
								// false.
		this.isBeveling = segment.p0.equals2D(segment.p1);
		this.skeletonFace = null; // No skeleton face assigned by default here
		if (isBeveling) {
			// Beveling edges should ideally have a skeleton face assigned later if needed
		}
	}

	// Constructor from basic segment and weight (maybe less common)
	public WavefrontEdge(LineSegment segment, double weight) {
		this(new WavefrontSupportingLine(segment, weight));
	}

	// Private constructor for internal operations like split
	private WavefrontEdge(WavefrontVertex v0, WavefrontVertex v1, WavefrontSupportingLine supportingLine, KineticTriangle incidentTriangle,
			SkeletonDCELFace skeletonFace, boolean isBeveling, boolean isInitial) {
		this.id = idCounter.incrementAndGet();
		this.vertex0 = v0;
		this.vertex1 = v1;
		this.supportingLine = supportingLine;
		this.canonicalSegment = new CanonicalSegment(supportingLine.getSegment());
		this.weight = supportingLine.getWeight();
		this.incidentTriangle = incidentTriangle;
		this.isInitial = isInitial; // Propagate initial status if needed (e.g. split of initial edge?)
		this.isBeveling = isBeveling;
		this.skeletonFace = skeletonFace;
		if (skeletonFace != null && isBeveling) {
			throw new IllegalStateException("skeletonFace must be null for beveling edges");
		}
		// Update vertex back pointers if vertices are provided
		updateVertexIncidentEdges();
	}

	// --- Getters ---

	public WavefrontSupportingLine getSupportingLine() {
		return supportingLine;
	}

	public LineSegment getSegment() {
		// Return the segment from the supporting line, which might be updated
		return supportingLine.getSegment();
	}

	public double getWeight() {
		return weight; // Use the stored weight
	}

	public WavefrontVertex getVertex(int index) {
		if (index == 0)
			return vertex0;
		if (index == 1)
			return vertex1;
		throw new IndexOutOfBoundsException("WavefrontEdge vertex index must be 0 or 1, got: " + index);
	}

	public KineticTriangle getIncidentTriangle() {
		return incidentTriangle;
	}

	public boolean isDead() {
		return isDead;
	}

	public WavefrontVertex getInitialVertex(int index) {
		if (!isInitial) {
			throw new IllegalStateException("Initial vertices only exist for initial edges (Edge " + id + ")");
		}
		if (index < 0 || index > 1) {
			throw new IndexOutOfBoundsException("Initial vertex index must be 0 or 1, got: " + index);
		}
		if (initialVertices[index] == null) {
			// This might happen if setInitialVertices wasn't called after assigning final
			// vertices
			System.err.println("Warning: Accessing unset initialVertex[" + index + "] for initial edge " + id);
			// Fallback to current vertex? Or throw? Let's throw for now.
			throw new IllegalStateException("Initial vertex " + index + " not set for initial edge " + id);
		}
		return initialVertices[index];
	}

	// --- Setters and State Changers ---

	/**
	 * Sets the incident triangle for this edge. Also invalidates the collapse spec.
	 * 
	 * @param triangle The new incident triangle.
	 */
	public void setIncidentTriangle(KineticTriangle triangle) {
		if (this.incidentTriangle != triangle) { // Avoid invalidation if no change
			this.incidentTriangle = triangle;
			invalidateCollapseSpec();
		}
	}

	/**
	 * Sets one of the edge's vertices (0 or 1).
	 * <p>
	 * IMPORTANT: This method only updates the edge's internal pointer. It does NOT
	 * update the vertex's back-pointer (`wavefronts` array). Use
	 * {@link #setVerticesAndUpdateAdj(WavefrontVertex, WavefrontVertex)} or
	 * {@link #setVertexAndUpdateAdj(int, WavefrontVertex)} for full updates.
	 * Invalidates the collapse spec.
	 *
	 * @param index The index (0 or 1) of the vertex to set.
	 * @param v     The vertex to set. Must not be null.
	 */
	public void setVertexRaw(int index, WavefrontVertex v) {
		if (v == null)
			throw new NullPointerException("Cannot set null vertex for edge " + id);
		boolean changed = false;
		if (index == 0) {
			if (this.vertex0 != v) {
				this.vertex0 = v;
				changed = true;
			}
		} else if (index == 1) {
			if (this.vertex1 != v) {
				this.vertex1 = v;
				changed = true;
			}
		} else {
			throw new IndexOutOfBoundsException("WavefrontEdge vertex index must be 0 or 1, got: " + index);
		}
		if (changed) {
			invalidateCollapseSpec();
		}
	}

	/**
	 * Sets one of the edge's vertices (0 or 1) AND updates that vertex's incident
	 * edge pointer (`wavefronts` array) back to this edge. Invalidates the collapse
	 * spec.
	 *
	 * @param index The index (0 or 1) of the vertex to set.
	 * @param v     The vertex to set. Must not be null.
	 */
	public void setVertexAndUpdateAdj(int index, WavefrontVertex v) {
		if (v == null)
			throw new NullPointerException("Cannot set null vertex for edge " + id);

		// Determine the corresponding index in the vertex's wavefront array
		// vertex0 (edge index 0) corresponds to vertex.wavefronts[1] (CW edge)
		// vertex1 (edge index 1) corresponds to vertex.wavefronts[0] (CCW edge)
		int vertexWavefrontIndex = (index == 0) ? 1 : 0;

		// Get the old vertex at this index to potentially clear its pointer
		WavefrontVertex oldVertex = (index == 0) ? this.vertex0 : this.vertex1;

		// Check if vertex is actually changing
		if (oldVertex == v) {
			// Vertex isn't changing, but ensure back pointer is correct
			if (!v.isInfinite()) {
				v.setIncidentEdge(vertexWavefrontIndex, this);
			}
			return; // No further action needed
		}

		// Clear old vertex's back pointer if it's not null and not infinite
		if (oldVertex != null && !oldVertex.isInfinite()) {
			// Only clear if the pointer actually points to this edge
			if (oldVertex.getWavefront(vertexWavefrontIndex) == this) {
				oldVertex.setIncidentEdge(vertexWavefrontIndex, null);
			} else {
				// This might indicate an inconsistency elsewhere
				System.err.println("Warning: Old vertex V" + oldVertex.id + " at index " + index + " of edge " + id
						+ " did not have correct back-pointer (expected edge " + id + ", found edge "
						+ (oldVertex.getWavefront(vertexWavefrontIndex) != null ? oldVertex.getWavefront(vertexWavefrontIndex).id : "null") + ")");
			}
		}

		// Update edge's internal pointer
		if (index == 0) {
			this.vertex0 = v;
		} else {
			this.vertex1 = v;
		}

		// Update new vertex's back pointer if it's not infinite
		if (!v.isInfinite()) {
			v.setIncidentEdge(vertexWavefrontIndex, this);
		}

		// Invalidate cache since vertex changed
		invalidateCollapseSpec();
	}

	/**
	 * Sets both vertices of the edge AND updates their incident edge pointers
	 * (`wavefronts` array) back to this edge. This is the recommended method for
	 * ensuring consistency after operations like flips or constraint moves.
	 * Invalidates the collapse spec.
	 *
	 * @param v0 The vertex for index 0 (CCW end relative to triangle). Must not be
	 *           null.
	 * @param v1 The vertex for index 1 (CW end relative to triangle). Must not be
	 *           null.
	 */
	public void setVerticesAndUpdateAdj(WavefrontVertex v0, WavefrontVertex v1) {
		if (v0 == null || v1 == null) {
			throw new NullPointerException("Cannot set null vertices for edge " + id);
		}
		if (v0 == v1) {
			throw new IllegalArgumentException("Cannot set the same vertex for both ends of edge " + id);
		}

		// Check if vertices actually changed to avoid unnecessary updates/invalidations
		boolean changed = (this.vertex0 != v0 || this.vertex1 != v1);

		// --- Update Vertex 0 ---
		// Clear old vertex0's back pointer (if needed)
		if (this.vertex0 != null && this.vertex0 != v0 && !this.vertex0.isInfinite()) {
			if (this.vertex0.getWavefront(1) == this) { // vertex0 corresponds to wavefronts[1]
				this.vertex0.setIncidentEdge(1, null);
			}
		}
		// Set new vertex0
		this.vertex0 = v0;
		// Set new vertex0's back pointer
		if (!v0.isInfinite()) {
			v0.setIncidentEdge(1, this); // vertex0 corresponds to wavefronts[1]
		}

		// --- Update Vertex 1 ---
		// Clear old vertex1's back pointer (if needed)
		if (this.vertex1 != null && this.vertex1 != v1 && !this.vertex1.isInfinite()) {
			if (this.vertex1.getWavefront(0) == this) { // vertex1 corresponds to wavefronts[0]
				this.vertex1.setIncidentEdge(0, null);
			}
		}
		// Set new vertex1
		this.vertex1 = v1;
		// Set new vertex1's back pointer
		if (!v1.isInfinite()) {
			v1.setIncidentEdge(0, this); // vertex1 corresponds to wavefronts[0]
		}

		// Invalidate collapse spec if vertices changed
		if (changed) {
			invalidateCollapseSpec();
		}
	}

	/**
	 * Updates the incident edge pointers (`wavefronts` array) for the current
	 * `vertex0` and `vertex1` to point back to this edge. Call this when the vertex
	 * objects themselves haven't changed but their relationship to this edge might
	 * have been established externally (e.g., in a constructor).
	 */
	private void updateVertexIncidentEdges() {
		if (vertex0 != null && !vertex0.isInfinite()) {
			vertex0.setIncidentEdge(1, this); // vertex0 -> wavefronts[1]
		}
		if (vertex1 != null && !vertex1.isInfinite()) {
			vertex1.setIncidentEdge(0, this); // vertex1 -> wavefronts[0]
		}
	}

	/**
	 * Marks the edge as dead.
	 */
	public void markDead() {
		if (!this.isDead) {
			this.isDead = true;
			invalidateCollapseSpec(); // Cannot collapse dead edge
			// Optional: Clear references to break cycles and help GC
			// clearReferences();
		}
	}

	/**
	 * Stores the current vertices as the initial vertices. Should only be called
	 * once for initial edges after vertices are finalized.
	 */
	public void setInitialVertices() {
		if (!isInitial) {
			throw new IllegalStateException("Can only set initial vertices for initial edges (Edge " + id + ")");
		}
		if (vertex0 == null || vertex1 == null) {
			throw new IllegalStateException("Vertices must be set before setting initial vertices for edge " + id);
		}
		// Check if already set?
		if (initialVertices[0] != null || initialVertices[1] != null) {
			System.err.println("Warning: Overwriting initial vertices for edge " + id);
		}
		initialVertices[0] = vertex0;
		initialVertices[1] = vertex1;
	}

	// --- Collapse Logic ---

	/** Invalidates the cached collapse spec. */
	private void invalidateCollapseSpec() {
		// Check flag first to avoid redundant work if already invalid
		if (collapseSpecValid) {
			collapseSpecValid = false;
			cachedCollapseSpec = null;
		}
	}

	/** Gets the cached collapse spec, computing it if necessary. */
	private EdgeCollapseSpec getCachedCollapseSpec() {
		if (isDead) {
			// Return a spec indicating no collapse for dead edges
			return EdgeCollapseSpec.NEVER;
		}
		if (!collapseSpecValid) {
			// Compute and cache the spec
			// Pass time=0 initially? Or rely on caller passing correct time in
			// getEdgeCollapse?
			// Let's assume computeCollapse needs a time. If cache is invalid, we don't know
			// 'currentTime'.
			// This suggests caching should happen inside getEdgeCollapse(timeNow).
			throw new IllegalStateException("Attempted to access invalid collapse spec cache for edge " + id);
		}
		if (cachedCollapseSpec == null) {
			// This implies collapseSpecValid was true, but cache is null - internal error
			throw new IllegalStateException("Internal Error: Collapse spec cache null despite being marked valid for edge " + id);
		}
		return cachedCollapseSpec;
	}

	/**
	 * Checks if the edge endpoints are moving parallel (resulting in ALWAYS or
	 * NEVER collapse).
	 * 
	 * @param timeNow The current simulation time.
	 * @return true if the edge collapse type is ALWAYS or NEVER.
	 */
	public boolean parallelEndpoints(double timeNow) {
		// Compute or get cached spec first
		EdgeCollapseSpec spec = getEdgeCollapse(timeNow);
		return spec.getType() == EdgeCollapseType.ALWAYS || spec.getType() == EdgeCollapseType.NEVER;
	}

	/**
	 * Gets the EdgeCollapseSpec, using cache if valid for the given time, otherwise
	 * computes and caches it.
	 * 
	 * @param timeNow The current simulation time.
	 * @return The EdgeCollapseSpec.
	 */
	public EdgeCollapseSpec getEdgeCollapse(double timeNow) {
		if (isDead) {
			return EdgeCollapseSpec.NEVER;
		}
		// Simple caching: invalidate always if time advances? No, time is input.
		// Caching is only useful if called multiple times *at the same logical step*
		// or if computation is expensive. Let's assume basic caching is desired.
		if (!collapseSpecValid) {
			// log.debug("Edge {}: Cache miss, computing collapse at time {}", id, timeNow);
			// // Use logger
			cachedCollapseSpec = computeCollapse(timeNow);
			collapseSpecValid = true; // Mark cache valid after computation
		} else {
			// log.debug("Edge {}: Cache hit.", id); // Use logger
		}
		// Should we check if cached spec time matches timeNow? Depends on how cache is
		// used.
		// If invalidated externally on state changes, this should be fine.
		if (cachedCollapseSpec == null) {
			throw new IllegalStateException("Internal Error: Collapse spec cache null after check/compute for edge " + id);
		}
		return cachedCollapseSpec;
	}

	/**
	 * Creates a CollapseSpec suitable for the EventQueue based on this edge's
	 * potential collapse.
	 * 
	 * @param component           The component ID.
	 * @param timeNow             Current simulation time.
	 * @param collapsingEdgeIndex The index (0, 1, or 2) this edge corresponds to
	 *                            within its incident triangle.
	 * @return The CollapseSpec.
	 */
	public CollapseSpec getCollapse(int component, double timeNow, int collapsingEdgeIndex) {
		assertEdgeSane(collapsingEdgeIndex);
		// Compute or get cached edge collapse spec
		EdgeCollapseSpec edgeSpec = getEdgeCollapse(timeNow);
		// Create the triangle-level collapse spec from the edge spec
		return CollapseSpec.fromEdgeCollapse(edgeSpec, incidentTriangle, collapsingEdgeIndex);
	}

	// --- Utility and Validation ---

	/** Basic sanity check before using the edge in collapse calculations. */
	private void assertEdgeSane(int collapsingEdgeIndex) {
		if (isDead) {
			throw new IllegalStateException("Operation on dead edge " + id);
		}
		if (collapsingEdgeIndex < 0 || collapsingEdgeIndex >= 3) {
			throw new IllegalArgumentException("Collapsing edge index must be 0, 1, or 2, got " + collapsingEdgeIndex);
		}
		if (incidentTriangle == null) {
			throw new IllegalStateException("Incident triangle must be set for edge " + id);
		}
		// Check if the triangle still points to this edge at the expected index
		if (incidentTriangle.getWavefront(collapsingEdgeIndex) != this) {
			throw new IllegalStateException("Incident triangle " + incidentTriangle.getName() + " does not reference this edge " + id + " at index "
					+ collapsingEdgeIndex + " (found: "
					+ (incidentTriangle.getWavefront(collapsingEdgeIndex) != null ? incidentTriangle.getWavefront(collapsingEdgeIndex).id : "null") + ")");
		}
		if (vertex0 == null || vertex1 == null) {
			throw new IllegalStateException("Vertices must be set for edge " + id);
		}
		// Additional check: ensure vertices match triangle's vertices for that edge
		// index
		WavefrontVertex triVCCW = incidentTriangle.getVertex(KineticTriangle.ccw(collapsingEdgeIndex));
		WavefrontVertex triVCW = incidentTriangle.getVertex(KineticTriangle.cw(collapsingEdgeIndex));
		if (vertex0 != triVCCW || vertex1 != triVCW) {
			throw new IllegalStateException("Edge " + id + " vertices (V" + vertex0.id + ", V" + vertex1.id + ")" + " do not match incident triangle "
					+ incidentTriangle.getName() + " vertices for edge index " + collapsingEdgeIndex + " (expected V" + (triVCCW != null ? triVCCW.id : "null")
					+ ", V" + (triVCW != null ? triVCW.id : "null") + ")");
		}
	}

	/** Optional: Clears references to potentially help GC. Use with caution. */
	private void clearReferences() {
		// Only clear if dead?
		if (isDead) {
			// Clear vertex back pointers if they point here
			if (vertex0 != null && !vertex0.isInfinite() && vertex0.getWavefront(1) == this) {
				vertex0.setIncidentEdge(1, null);
			}
			if (vertex1 != null && !vertex1.isInfinite() && vertex1.getWavefront(0) == this) {
				vertex1.setIncidentEdge(0, null);
			}
			// Clear internal pointers
			this.vertex0 = null;
			this.vertex1 = null;
			this.incidentTriangle = null;
			this.cachedCollapseSpec = null;
			// Keep initialVertices? Probably useful for history/debugging even if dead.
		}
	}

	@Override
	public String toString() {
		return "WE" + id + "[" + (vertex0 != null ? "V" + vertex0.id : "null") + "->" + (vertex1 != null ? "V" + vertex1.id : "null") + "]"
				+ (isInitial ? " (Initial)" : "") + (isBeveling ? " (Bevel)" : "") + (isDead ? " (DEAD)" : "");
	}

	@Override
	public int hashCode() {
		return Long.hashCode(id);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		WavefrontEdge other = (WavefrontEdge) obj;
		return id == other.id;
	}

	// Enhanced computeCollapse to match C++ logic
	public EdgeCollapseSpec computeCollapse(double currentTime) {
		WavefrontVertex wfv0 = vertex0;
		WavefrontVertex wfv1 = vertex1;
		if (wfv0 == null || wfv1 == null) {
			throw new IllegalStateException("Vertices must be set to compute collapse");
		}

		Coordinate v0 = wfv0.getVelocity();
		Coordinate v1 = wfv1.getVelocity();
		Coordinate p0 = wfv0.getInitialPosition();
		Coordinate p1 = wfv1.getInitialPosition();

		// Use JTS Orientation.index() to determine orientation
		int orient = Orientation.index(p0, p1, new Coordinate(p0.x + v0.x, p0.y + v0.y));
		int orientation;
		switch (orient) {
			case org.locationtech.jts.algorithm.Orientation.LEFT :
				orientation = Orientation.LEFT;
				break;
			case org.locationtech.jts.algorithm.Orientation.RIGHT :
				orientation = Orientation.RIGHT;
				break;
			case org.locationtech.jts.algorithm.Orientation.COLLINEAR :
			default :
				orientation = Orientation.COLLINEAR;
				break;
		}

		if (orientation != Orientation.LEFT) {
			if (orientation == Orientation.RIGHT) {
				return new EdgeCollapseSpec(EdgeCollapseType.PAST, Double.NaN);
			} else { // COLLINEAR
				double sqDist = p0.distance(p1);
				if (sqDist < 1e-9) {
					return new EdgeCollapseSpec(EdgeCollapseType.ALWAYS, currentTime);
				} else {
					return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
				}
			}
		} else {
			// Project onto edge direction
			double edgeDelta, wfvsDelta;
			if (!supportingLine.isVertical()) {
				edgeDelta = p1.x - p0.x;
				wfvsDelta = v0.x - v1.x;
			} else {
				edgeDelta = p1.y - p0.y;
				wfvsDelta = v0.y - v1.y;
			}

			if (Math.abs(edgeDelta) < 1e-9 || Math.abs(wfvsDelta) < 1e-9) {
				throw new IllegalStateException("Invalid edge or velocity delta");
			}

			double time = edgeDelta / wfvsDelta;
			if (time < currentTime - 1e-9) {
				throw new IllegalStateException("Computed collapse time is in the past");
			}
			return new EdgeCollapseSpec(EdgeCollapseType.FUTURE, time);
		}
	}

	class SkeletonDCELFace {
		// TODO Define as needed
	}
}