package com.github.micycle1.surferj.kinetics;

import java.util.List;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicLong;

import org.apache.commons.lang3.tuple.Pair;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.math.Vector2D;

import com.github.micycle1.surferj.SurfConstants;
import com.github.micycle1.surferj.TriangulationUtils;
import com.github.micycle1.surferj.collapse.CollapseSpec;
import com.github.micycle1.surferj.collapse.EdgeCollapseSpec;
import com.github.micycle1.surferj.collapse.EdgeCollapseType;

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

	// --- NEW Private Constructor (for split) ---
	/**
	 * Private constructor used internally, primarily by the split() method, to
	 * create new edge halves. Allows setting initial vertex references and
	 * inheriting properties from the edge being split.
	 *
	 * @param v0               The first vertex (can be null initially).
	 * @param v1               The second vertex (can be null initially).
	 * @param supportingLine   The supporting line (shared with the original edge).
	 * @param incidentTriangle The initial incident triangle (may change).
	 * @param skeletonFace     The associated skeleton face (shared with the
	 *                         original edge).
	 * @param isBeveling       Whether this edge represents a bevel.
	 */
	private WavefrontEdge(WavefrontVertex v0, WavefrontVertex v1, WavefrontSupportingLine supportingLine, KineticTriangle incidentTriangle,
			SkeletonDCELFace skeletonFace, boolean isBeveling) {
		this.id = idCounter.incrementAndGet();
		this.vertex0 = v0;
		this.vertex1 = v1;
		this.supportingLine = Objects.requireNonNull(supportingLine, "SupportingLine cannot be null in private constructor");
		this.canonicalSegment = new CanonicalSegment(supportingLine.getSegment()); // Re-derive canonical segment
		this.weight = supportingLine.getWeight(); // Get weight from supporting line
		this.incidentTriangle = incidentTriangle;
		this.skeletonFace = skeletonFace;
		this.isBeveling = isBeveling;
		this.isInitial = false; // Split edges are never 'initial' input edges

		// Important: Link edge back to vertices if needed by vertex logic.
		// Be careful not to double-link if the public constructors also do this.
		// This linking might be better done when the vertices are finalized.
		// if (v0 != null) v0.addIncidentEdge(this); // Example
		// if (v1 != null) v1.addIncidentEdge(this); // Example
		updateVertexIncidentEdges(); // NOTE correct???
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
		if (index == 0) {
			return vertex0;
		}
		if (index == 1) {
			return vertex1;
		}
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
	    if (index < 0 || index > 1) {
	        throw new IndexOutOfBoundsException("WavefrontEdge vertex index must be 0 or 1, got: " + index);
	    }
	    // Allow null vertex? Let's assume yes for now, caller beware.
	    // if (v == null) {
	    //     LOGGER.warn("Setting null vertex at index {} for edge {}", index, this.id);
	    // }

	    WavefrontVertex oldVertex = (index == 0) ? this.vertex0 : this.vertex1;
	    boolean changed = (oldVertex != v);

	    // 1. Update edge's internal pointer
	    if (index == 0) {
	        this.vertex0 = v;
	    } else { // index == 1
	        this.vertex1 = v;
	    }

	    // 2. Update NEW vertex's back pointer (if vertex is not null)
	    if (v != null && !v.isInfinite()) {
	         // Index 0 of edge corresponds to vertex's side 1 (outgoing CW)
	         // Index 1 of edge corresponds to vertex's side 0 (incoming CCW)
	         // Make sure WavefrontVertex.setIncidentEdge handles nulls if needed
	         v.setIncidentEdge(1 - index, this); // <<< THIS IS THE MISSING LINK
	    }

	    // 3. CRITICAL: DO NOT MODIFY oldVertex's pointers TO edges here.

	    // 4. Invalidate spec only if a vertex actually changed
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
		if (v == null) {
			throw new NullPointerException("Cannot set null vertex for edge " + id);
		}

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

	    boolean changed = (this.vertex0 != v0 || this.vertex1 != v1);

	    // Only update internal state and NEW vertices' back pointers
	    this.vertex0 = v0;
	    if (v0 != null && !v0.isInfinite()) {
	        // Associate this edge with the correct side of the new vertex v0
	        // Assuming vertex0 maps to incoming/side 1 for v0
	        v0.setIncidentEdge(1, this);
	    }

	    this.vertex1 = v1;
	    if (v1 != null && !v1.isInfinite()) {
	        // Associate this edge with the correct side of the new vertex v1
	        // Assuming vertex1 maps to outgoing/side 0 for v1
	        v1.setIncidentEdge(0, this);
	    }

	    // *** DO NOT MODIFY OLD VERTICES ***
	    // The logic that cleared the old vertices' pointers is removed.

	    if (changed) {
	        invalidateCollapseSpec();
	    }
	}

	public void setVerticesRaw(WavefrontVertex v0, WavefrontVertex v1) {
		// Ensure the vertices match the segment coordinates conceptually
		if (v0 != null && v1 != null && !v0.isInfinite && !v1.isInfinite) {
			Coordinate c0 = getSegment().p0;
			Coordinate c1 = getSegment().p1;
			// Check if vertices align with segment ends, allowing for swapped order
			boolean match1 = (v0.posStart.equals2D(c0) && v1.posStart.equals2D(c1));
			boolean match2 = (v0.posStart.equals2D(c1) && v1.posStart.equals2D(c0));
			if (!match1 && !match2) {
				throw new IllegalArgumentException("Vertices " + v0 + ", " + v1 + " do not match edge segment ends " + getSegment());
			}
		} else if (v0 == null || v1 == null) {
			throw new IllegalArgumentException("Cannot set null vertices for WavefrontEdge");
		}
		this.vertex0 = v0;
		this.vertex1 = v1;
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

	/**
	 * Splits this wavefront edge into two new edges at an implicit point. Marks the
	 * current edge as dead and adds the two new edge halves to the provided list.
	 * The new vertices connecting the halves must be set by the caller.
	 *
	 * @param wavefrontEdges The list to which the new edge halves will be added.
	 * @param splitPosition  The geometric location of the split (used by caller,
	 *                       not directly here).
	 * @return A Pair containing the two new WavefrontEdge instances (first = part
	 *         connected to original vertex 0, second = part connected to original
	 *         vertex 1).
	 */
	public Pair<WavefrontEdge, WavefrontEdge> split(List<WavefrontEdge> wavefrontEdges, Coordinate splitPosition) {
		// NOTE: splitPosition isn't used in the C++ code provided for *creating* the
		// edges,
		// but it's essential context for the caller
		// (KineticEventHandler.handleSplitEvent)
		// to create the new WavefrontVertices *at* that position.

		if (isDead) {
			// Handle appropriately - maybe return null or throw exception?
			// C++ doesn't show a check, but it's good practice.
			System.err.println("Warning: Attempting to split an already dead WavefrontEdge.");
			// Returning null might be problematic downstream. Maybe return existing dead
			// parts if tracked?
			// For now, let's proceed but log heavily. If split is called, it implies an
			// event occurred.
		}

		markDead(); // Mark this edge instance as dead

		// Create the first new edge half: Original vertex[0] ----> NULL (new vertex to
		// be set by caller)
		// It uses the same supporting line, incident triangle, skeleton face, and
		// beveling status.
		WavefrontEdge edgeA = new WavefrontEdge(vertex0, // Start vertex is the same as original start
				null, // End vertex will be the new split vertex (set by caller)
				this.supportingLine, this.incidentTriangle, // Initially assumes same incident triangle
				this.skeletonFace, this.isBeveling);

		// Create the second new edge half: NULL (new vertex to be set by caller) ---->
		// Original vertex[1]
		WavefrontEdge edgeB = new WavefrontEdge(null, // Start vertex will be the new split vertex (set by caller)
				vertex1, // End vertex is the same as original end
				this.supportingLine, this.incidentTriangle, // Initially assumes same incident triangle
				this.skeletonFace, // Assumes same skeleton face
				this.isBeveling);

		// Add the new edges to the master list
		wavefrontEdges.add(edgeA);
		wavefrontEdges.add(edgeB);

		// Return the pair of new edges
		return Pair.of(edgeA, edgeB);
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
		WavefrontVertex triVCCW = incidentTriangle.getVertex(TriangulationUtils.ccw(collapsingEdgeIndex));
		WavefrontVertex triVCW = incidentTriangle.getVertex(TriangulationUtils.cw(collapsingEdgeIndex));
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
		if (this == obj) {
			return true;
		}
		if (obj == null || getClass() != obj.getClass()) {
			return false;
		}
		WavefrontEdge other = (WavefrontEdge) obj;
		return id == other.id;
	}

	/**
	 * Computes the collapse specification for this edge based on the relative
	 * movement of its endpoint vertices. Mirrors the logic of the C++
	 * WavefrontEdge::compute_collapse. Calculates the physical collapse time, which
	 * might be in the past relative to currentTime. The caller is responsible for
	 * interpreting this.
	 * <p>
	 * The computeCollapse method determines when or if an edge, defined by two
	 * vertices (vertex0 and vertex1), will collapse in a 2D simulation, likely
	 * related to computational geometry or computer graphics. It takes a
	 * currentTime parameter (the current simulation time) and returns an
	 * EdgeCollapseSpec object, which specifies:
	 *
	 * Type: NEVER (won’t collapse), ALWAYS (currently collapsed), or FUTURE (will
	 * collapse at a specific time). Time: The time of collapse (or Double.NaN if
	 * not applicable). The method uses the initial positions (p0, p1) and
	 * velocities (v0, v1) of the vertices to predict collapse behavior, considering
	 * their relative orientation and motion.
	 *
	 * @param currentTime The current simulation time (used for context, e.g., in
	 *                    ALWAYS case, but NOT for filtering).
	 * @return An EdgeCollapseSpec indicating the type and time of collapse.
	 * @throws IllegalStateException if vertices are not set.
	 */
	public EdgeCollapseSpec computeCollapse(double currentTime) {
		WavefrontVertex wfv0 = vertex0;
		WavefrontVertex wfv1 = vertex1;
		if (wfv0 == null || wfv1 == null) {
			throw new IllegalStateException("Vertices must be set to compute collapse for Edge " + id);
		}

		// Use initial positions and velocities for the core calculation,
		// matching the C++ logic.
		Coordinate p0 = wfv0.getInitialPosition();
		Coordinate p1 = wfv1.getInitialPosition();
		Vector2D v0 = Vector2D.create(wfv0.getVelocity()); // Assuming getVelocity() returns (0,0) for static/infinite
		Vector2D v1 = Vector2D.create(wfv1.getVelocity());

		// Check relative orientation of velocities (v0 relative to v1, or vice-versa?)
		// C++ CGAL::orientation(v0, v1) checks orientation of vector v1 relative to v0
		// (origin implicit).
		// This checks if v1 is left/right/collinear with v0.
		// Let's use cross product: v0.x * v1.y - v0.y * v1.x
		double crossProduct = v0.getX() * v1.getY() - v0.getY() * v1.getX();
		int orientation;

		if (Math.abs(crossProduct) < SurfConstants.ZERO_AREA) { // Use an appropriate tolerance for cross product (area)
			orientation = Orientation.COLLINEAR;
		} else if (crossProduct > 0) {
			orientation = Orientation.LEFT; // v1 is left of v0
		} else {
			orientation = Orientation.RIGHT; // v1 is right of v0
		}

		if (orientation != Orientation.LEFT) {
			// Velocities are collinear or v1 is to the right of v0 (implying divergence or
			// parallel movement away?)
			// C++ logic: If not LEFT_TURN -> PAST or (NEVER/ALWAYS if COLLINEAR)
			if (orientation == Orientation.RIGHT) {
				// C++ returns PAST. This means they will never converge *in the future*
				// according to this velocity orientation check. Let's map to NEVER
				// as PAST isn't very informative without a time.
				return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
				// Or maybe keep PAST if needed upstream? Let's use NEVER for simplicity.
				// return new EdgeCollapseSpec(EdgeCollapseType.PAST, Double.NaN);
			} else { // COLLINEAR
				// Check initial distance. If initially collapsed, it's ALWAYS. Otherwise NEVER.
				// Note: C++ checks sqdist == CORE_ZERO based on initial positions.
				double initialSqDist = p0.distanceSq(p1);
				if (initialSqDist < SurfConstants.ZERO_DIST_SQ) { // Use squared distance tolerance
					// Need to verify they aren't moving apart despite being collinear
					Vector2D relPos = Vector2D.create(p1).subtract(Vector2D.create(p0));
					Vector2D relVel = v1.subtract(v0);
					if (relPos.dot(relVel) >= -SurfConstants.ZERO_SPEED_SQ) { // Check if moving apart (dot >= 0) or static
						// If they start collapsed and are static or moving apart (relative vel along
						// axis), they stay collapsed
						return new EdgeCollapseSpec(EdgeCollapseType.ALWAYS, currentTime); // Return current time as reference
					} else {
						// Started collapsed but moving towards each other (impossible if truly
						// collapsed?)
						// This case is tricky. Assume ALWAYS if initially zero dist.
						System.err.println("Warning: Edge " + id + " initially collapsed but velocities are convergent?");
						return new EdgeCollapseSpec(EdgeCollapseType.ALWAYS, currentTime);
					}
				} else {
					// Check if moving towards each other
					Vector2D relPos = Vector2D.create(p1).subtract(Vector2D.create(p0));
					Vector2D relVel = v1.subtract(v0);
					// If velocities are zero, never collapse. If non-zero but parallel, check if
					// converging.
					if (relVel.lengthSquared() < SurfConstants.ZERO_SPEED_SQ) { // Velocities are effectively the same
						return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
					}
					if (relPos.dot(relVel) < -SurfConstants.ZERO_SPEED_SQ) { // Dot product negative -> Converging
						// Calculate time to collapse based on 1D motion along the line
						// Time = - (relPos . relVel) / (relVel . relVel) ? No, that's projection time.
						// Time = distance / relative_speed = |relPos| / |relVel| (if moving directly
						// towards)
						// Let's use the C++ projection method for consistency if needed, otherwise
						// assume NEVER for now.
						// The C++ code relies on the projection logic below even for collinear cases?
						// If v0 and v1 are collinear, wfvs_delta might be zero.
						// Revisit this: Let's assume for now, if collinear and not initially collapsed,
						// NEVER.
						// This might be wrong if they are collinear and moving towards each other.
						// Let's recalculate using the projection logic from the C++ `else` block.
						// Fall through to the projection logic.
					} else {
						// Collinear, not collapsed, not converging -> NEVER
						return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
					}
				}
			}
		}

		// --- Case: Orientation == LEFT (Potentially Converging) ---
		// Calculate collapse time using projection method from C++.
		// This finds time 't' where Pi(t) projected onto some axis are equal.
		// t = (P1_proj(0) - P0_proj(0)) / (V0_proj - V1_proj)
		// t = edge_delta / wfvs_delta

		double edge_delta; // Difference in projected initial positions
		double wfvs_delta; // Difference in projected velocities (V0_proj - V1_proj)

		// Use projection axis based on supporting line (prefer non-vertical)
		boolean isVertical = supportingLine.isVertical(SurfConstants.VERTICAL_SLOPE_TOL); // Use tolerance

		if (!isVertical) {
			// Project onto X-axis
			edge_delta = p1.x - p0.x;
			wfvs_delta = v0.getX() - v1.getX(); // Note: C++ uses v0.x() - v1.x()
		} else {
			// Project onto Y-axis
			edge_delta = p1.y - p0.y;
			wfvs_delta = v0.getY() - v1.getY(); // Note: C++ uses v0.y() - v1.y()
		}

		// Check for zero relative velocity along projection axis
		if (Math.abs(wfvs_delta) < SurfConstants.ZERO_SPEED) {
			// If projected velocities are the same, they maintain their projected distance.
			// Check if they were initially collapsed along this projection.
			if (Math.abs(edge_delta) < SurfConstants.ZERO_DIST) {
				// They start at same projected position and move at same projected speed.
				// This doesn't guarantee collapse, only that they are relatively static along
				// this axis.
				// Need to check the *other* axis or actual distance.
				// If orientation was LEFT, they can't be fully collapsed initially unless
				// edge_delta=0 AND crossProduct=0.
				// This case implies parallel movement. Return NEVER.
				return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
			} else {
				// Maintain separation along this axis -> NEVER collapse by this projection
				// method.
				return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
			}
		}

		// Calculate potential collapse time
		double time = edge_delta / wfvs_delta;

		// *** CHANGE HERE: Do not compare with currentTime or throw exception ***
		// The calculated 'time' is the physical time the projected coordinates
		// coincide.

		// Determine midpoint at collapse time 'time'
		Coordinate collapseP0 = wfv0.getPositionAt(time);
		Coordinate collapseP1 = wfv1.getPositionAt(time);

		// Sanity check: the distance at this time should be very small IF the
		// orientation was LEFT
		// and wfvs_delta was non-zero. This projection method finds when they align on
		// one axis.
		double collapseDistSq = collapseP0.distanceSq(collapseP1);
		if (collapseDistSq > SurfConstants.ZERO_DIST_SQ * 100) { // Use larger tolerance?
			System.err.println(
					"Warning: Edge " + id + " collapse time " + time + " calculated via projection results in non-zero distance squared: " + collapseDistSq);
			// This might happen if velocities were nearly parallel or initial orientation
			// check was borderline.
			// Treat as NEVER? Or proceed cautiously? Let's return NEVER.
			return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
		}

		Coordinate collapsePoint = new Coordinate((collapseP0.x + collapseP1.x) / 2.0, (collapseP0.y + collapseP1.y) / 2.0);
		boolean sameInitial = wfv0.getInitialPosition().equals2D(wfv1.getInitialPosition()); // Needed for spec

		// Return the calculated time, regardless of whether it's past/present/future.
		return new EdgeCollapseSpec(EdgeCollapseType.FUTURE, time); // Use FUTURE type as placeholder
	}

	/**
	 * The goal is to maintain strict equivalence in outcomes (e.g., orientation
	 * handling, collinear cases, and time computation) while allowing minor
	 * improvements like better error handling or precision tolerance, as long as
	 * they don’t deviate from the C++ logic.
	 *
	 * @param currentTime
	 * @return
	 */
	public EdgeCollapseSpec computeCollapseSimple(double currentTime) {
		WavefrontVertex wfv0 = vertex0;
		WavefrontVertex wfv1 = vertex1;
		if (wfv0 == null || wfv1 == null) {
			throw new IllegalStateException("Vertices must be set to compute collapse for Edge " + id);
		}

		Coordinate p0 = wfv0.getInitialPosition();
		Coordinate p1 = wfv1.getInitialPosition();
		Vector2D v0 = Vector2D.create(wfv0.getVelocity());
		Vector2D v1 = Vector2D.create(wfv1.getVelocity());

		double crossProduct = v0.getX() * v1.getY() - v0.getY() * v1.getX();
		int orientation;

		if (Math.abs(crossProduct) < SurfConstants.ZERO_AREA) {
			orientation = Orientation.COLLINEAR;
		} else if (crossProduct > 0) {
			orientation = Orientation.LEFT;
		} else {
			orientation = Orientation.RIGHT;
		}

		if (orientation != Orientation.LEFT) {
			if (orientation == Orientation.RIGHT) {
				return new EdgeCollapseSpec(EdgeCollapseType.PAST, Double.NaN);
			} else { // COLLINEAR
				double initialSqDist = p0.distanceSq(p1);
				if (initialSqDist < SurfConstants.ZERO_DIST_SQ) {
					return new EdgeCollapseSpec(EdgeCollapseType.ALWAYS, currentTime);
				} else {
					return new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
				}
			}
		} else {
			double edge_delta;
			double wfvs_delta;
			boolean isVertical = supportingLine.isVertical(SurfConstants.VERTICAL_SLOPE_TOL);

			if (!isVertical) {
				edge_delta = p1.x - p0.x;
				wfvs_delta = v0.getX() - v1.getX();
			} else {
				edge_delta = p1.y - p0.y;
				wfvs_delta = v0.getY() - v1.getY();
			}

			if (Math.abs(wfvs_delta) < SurfConstants.ZERO_SPEED) {
				throw new IllegalStateException("Invalid edge or velocity delta");
			}

			double time = edge_delta / wfvs_delta;
			if (time < currentTime - SurfConstants.TIME_TOL) {
				return new EdgeCollapseSpec(EdgeCollapseType.PAST, Double.NaN);
			}
			return new EdgeCollapseSpec(EdgeCollapseType.FUTURE, time);
		}
	}

	class SkeletonDCELFace {
		// TODO Define as needed
	}
}