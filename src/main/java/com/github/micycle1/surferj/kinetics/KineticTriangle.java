package com.github.micycle1.surferj.kinetics;

import static com.github.micycle1.surferj.TriangulationUtils.ccw;
import static com.github.micycle1.surferj.TriangulationUtils.cw;

import java.util.Arrays;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicLong;

import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Pair;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.math.Vector2D;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

// Import logging framework if desired, e.g. import org.slf4j.Logger; import org.slf4j.LoggerFactory;
import com.github.micycle1.surferj.SurfConstants;
import com.github.micycle1.surferj.collapse.CollapseSpec;
import com.github.micycle1.surferj.collapse.CollapseType;
import com.github.micycle1.surferj.collapse.EdgeCollapseSpec;
import com.github.micycle1.surferj.collapse.EdgeCollapseType;
import com.github.micycle1.surferj.collapse.Polynomial;
import com.github.micycle1.surferj.collapse.QuadraticSolver;
import com.github.micycle1.surferj.kinetics.WavefrontVertex.InfiniteSpeedType;

/**
 * A triangle, part of the kinetic triangulation.
 * <p>
 * Maintains a kinetic triangulation of the area not yet swept by the wavefront.
 * <p>
 * - A triangle always has three vertices (one of which may be "infinite"). -
 * Triangle edges may be part of the wavefront (constrained), pointing to the
 * original input edge. - Each side either has a neighbor (KineticTriangle) or a
 * wavefront (WavefrontEdge constraint).
 *
 * Mirrors the C++ KineticTriangle class.
 */
public class KineticTriangle {

	private static final Logger LOGGER = LoggerFactory.getLogger(KineticTriangle.class);

	// NOTE quite confident this is thoroughly implemented!
	// NOTE good starting point to make other code

	// Static variables for stuck loop detection in get_generic_collapse_time() (as
	// in C++)
	private static double lastTime = Double.NEGATIVE_INFINITY;
	private static int loopDetectionCount = 0;
	private static final int LOOP_DETECTION_LIMIT = 1000;

	private static final AtomicLong ID_COUNTER = new AtomicLong(0);

	public final long id;
	public int component; // Represents the polygon component the triangle belongs to. Simplified in some
							// logic.

	// Vertices in CCW order
	private final WavefrontVertex[] vertices = new WavefrontVertex[3];

	// For each edge opposite vertex i: EITHER neighbors[i] OR wavefronts[i] is
	// non-null
	private final KineticTriangle[] neighbors = new KineticTriangle[3];
	private final WavefrontEdge[] wavefronts = new WavefrontEdge[3]; // Constraint edges (final removed to allow modification in flip/move)

	private boolean isDead = false;
	private boolean isDying = false; // Marked when scheduled for death

	// --- Collapse Spec Caching ---
	private CollapseSpec cachedCollapseSpec;
	private boolean isCollapseSpecValid = false;
	// --- Debugging aid (optional, like C++) ---
	// private final WavefrontVertex[] collapseSpecComputedWithVertices = new
	// WavefrontVertex[3]; // Make final if used

	public KineticTriangle() {
		// NOTE this constructor, not in original.
		// required now because KineticTriangulation initalises vertices lazily.
		this.id = ID_COUNTER.getAndIncrement();
	}

	/**
	 * Creates a new KineticTriangle with a unique ID.
	 *
	 * @param component The component ID this triangle belongs to.
	 */
	public KineticTriangle(int component) {
		this.id = ID_COUNTER.getAndIncrement();
		this.component = component;
	}

	/**
	 * Constructor primarily for testing purposes.
	 *
	 * @param component The component ID.
	 * @param v0        Vertex 0.
	 * @param v1        Vertex 1.
	 * @param v2        Vertex 2.
	 */
	public KineticTriangle(int component, WavefrontVertex v0, WavefrontVertex v1, WavefrontVertex v2) {
		this(component);
		// Use the proper setter to ensure adjacent edges are updated if constrained
		// later
		// However, for simple test setup, direct assignment might be okay if
		// constraints aren't involved yet.
		// Let's use direct assignment here assuming it's for initial setup before
		// linking.
		this.vertices[0] = v0;
		this.vertices[1] = v1;
		this.vertices[2] = v2;
	}
	// --- Getters ---

	public long getId() {
		return id;
	}

	public int getComponent() {
		return component;
	}

	public WavefrontVertex getVertex(int index) {
		if (index < 0 || index >= 3) {
			throw new IllegalArgumentException("Invalid vertex index: " + index);
		}
		return vertices[index];
	}

	public KineticTriangle getNeighbor(int index) {
		if (index < 0 || index >= 3) {
			throw new IllegalArgumentException("Invalid neighbor index: " + index);
		}
		return neighbors[index];
	}

	public WavefrontEdge getWavefront(int index) {
		if (index < 0 || index >= 3) {
			throw new IllegalArgumentException("Invalid wavefront index: " + index);
		}
		return wavefronts[index];
	}

	public boolean isConstrained(int index) {
		if (index < 0 || index >= 3) {
			throw new IllegalArgumentException("Invalid constraint check index: " + index);
		}
		return wavefronts[index] != null;
	}

	public boolean isDead() {
		return isDead;
	}

	public boolean isDying() {
		return isDying;
	}

	public boolean isCollapseSpecValid() {
		return isCollapseSpecValid;
	}

	public boolean isUnbounded() {
		// Check if vertices array is initialized before accessing
		return (vertices[0] != null && vertices[0].isInfinite()) || (vertices[1] != null && vertices[1].isInfinite())
				|| (vertices[2] != null && vertices[2].isInfinite());
	}

	public int getInfiniteVertexIndex() {
		int infiniteCount = 0;
		int index = -1;
		for (int i = 0; i < 3; i++) {
			if (vertices[i] != null && vertices[i].isInfinite()) {
				infiniteCount++;
				index = i;
			}
		}
		if (infiniteCount != 1) {
			throw new IllegalStateException("Expected exactly one infinite vertex in unbounded triangle " + getName() + ", found " + infiniteCount);
		}
		return index;
	}

	public InfiniteSpeedType hasVertexInfiniteSpeed() {
		boolean hasOpposing = false;
		boolean hasWeighted = false;
		for (int i = 0; i < 3; i++) {
			// Ensure vertex exists before checking speed type
			if (vertices[i] != null) {
				InfiniteSpeedType speedType = vertices[i].getInfiniteSpeed();
				if (speedType == InfiniteSpeedType.OPPOSING) {
					hasOpposing = true;
					break; // Opposing takes precedence
				} else if (speedType == InfiniteSpeedType.WEIGHTED) {
					hasWeighted = true;
				}
			}
		}
		if (hasOpposing) {
			return InfiniteSpeedType.OPPOSING;
		} else if (hasWeighted) {
			return InfiniteSpeedType.WEIGHTED;
		} else {
			return InfiniteSpeedType.NONE;
		}
	}

	public int getInfiniteSpeedOpposingVertexIndex() {
		for (int i = 0; i < 3; i++) {
			if (vertices[i] != null && vertices[i].getInfiniteSpeed() == InfiniteSpeedType.OPPOSING) {
				return i;
			}
		}
		throw new IllegalStateException("No OPPOSING infinite speed vertex found in " + getName());
	}

	public String getName() {
		return "KT" + id;
	}

	// --- Setters and State Changers ---

	public void setVertex(int index, WavefrontVertex v) {
		if (index < 0 || index >= 3) {
			throw new IllegalArgumentException("Invalid vertex index: " + index);
		}
		if (v == null) {
			throw new NullPointerException("Cannot set null vertex for " + getName());
		}
		vertices[index] = v;

		// Update adjacent wavefront edges if they exist.
		int edgeIdxCW = cw(index); // Index of edge whose CCW vertex is 'index'
		if (isConstrained(edgeIdxCW)) {
			WavefrontEdge edge = wavefronts[edgeIdxCW];
			// Check edge exists before calling method
			if (edge != null) {
				edge.setVertexAndUpdateAdj(0, v); // Vertex 'index' is vertex 0 (CCW) for edge 'edgeIdxCW'
			} else {
				// This case implies isConstrained was true but wavefronts[edgeIdxCW] became
				// null concurrently? Unlikely but possible.
				LOGGER.warn("Constraint flag mismatch during setVertex for " + getName() + " edge " + edgeIdxCW);
			}
		}

		int edgeIdxCCW = ccw(index); // Index of edge whose CW vertex is 'index'
		if (isConstrained(edgeIdxCCW)) {
			WavefrontEdge edge = wavefronts[edgeIdxCCW];
			if (edge != null) {
				edge.setVertexAndUpdateAdj(1, v); // Vertex 'index' is vertex 1 (CW) for edge 'edgeIdxCCW'
			} else {
				LOGGER.warn("Constraint flag mismatch during setVertex for " + getName() + " edge " + edgeIdxCCW);
			}
		}

		invalidateCollapseSpec();
	}

	public void setNeighborRaw(int index, KineticTriangle neighbor) {
		if (index < 0 || index >= 3) {
			throw new IllegalArgumentException("Invalid neighbor index: " + index);
		}
		// Allow setting null neighbor even if constrained (e.g., during cleanup)
		// But prevent setting non-null neighbor if constrained
		if (neighbor != null && wavefronts[index] != null) {
			throw new IllegalStateException("Cannot set neighbor " + neighbor.getName() + " for constrained edge " + index + " on triangle " + getName());
		}
		// Check component compatibility if setting a non-null neighbor
		if (neighbor != null && neighbor.component != this.component) {
			LOGGER.warn("Setting neighbor " + neighbor.getName() + " (comp " + neighbor.component + ") for triangle " + getName() + " (comp " + this.component
					+ ") with different component.");
			// Depending on strictness, could throw an exception here.
		}

		neighbors[index] = neighbor;
		if (neighbor != null) {
			wavefronts[index] = null; // Ensure only one is set if neighbor is non-null
		}
		invalidateCollapseSpec();
	}

	public void setWavefront(int index, WavefrontEdge edge) {
		if (index < 0 || index >= 3) {
			throw new IllegalArgumentException("Invalid wavefront index: " + index);
		}
		// Allow setting null wavefront even if neighbor exists (e.g., during cleanup)
		// But prevent setting non-null wavefront if neighbor exists
		if (edge != null && neighbors[index] != null) {
			throw new IllegalStateException("Cannot set wavefront edge " + edge.id + " for edge " + index + " which has neighbor " + neighbors[index].getName()
					+ " on triangle " + getName());
		}
		wavefronts[index] = edge;
		if (edge != null) {
			neighbors[index] = null; // Ensure only one is set if edge is non-null
			edge.setIncidentTriangle(this);
		}
		invalidateCollapseSpec();
	}

	public void setNeighbors(KineticTriangle n0, KineticTriangle n1, KineticTriangle n2) {
		setNeighborRaw(0, n0);
		setNeighborRaw(1, n1);
		setNeighborRaw(2, n2);
	}

	public void setWavefronts(WavefrontEdge w0, WavefrontEdge w1, WavefrontEdge w2) {
		setWavefront(0, w0);
		setWavefront(1, w1);
		setWavefront(2, w2);
	}

	public void invalidateCollapseSpec() {
		if (isCollapseSpecValid) { // Avoid redundant invalidation if already invalid
			this.isCollapseSpecValid = false;
			this.cachedCollapseSpec = null;
			// Optional: Invalidate debug state
			// Arrays.fill(collapseSpecComputedWithVertices, null);
		}
	}

	public void markDying() {
		if (isDead || isDying) {
			return; // Already dead or dying
		}
		this.isDying = true;
		// Do NOT invalidate collapse spec here, event queue relies on it being valid
		// until the event is processed or discarded.
	}

	public void markDead() {
		if (isDead) {
			// log.warn("markDead() called multiple times for {}", getName()); // Use logger
			// if available
			LOGGER.warn("markDead() called multiple times for " + getName());
			return;
		}
		this.isDead = true;
		this.isDying = true; // Dead implies Dying
		invalidateCollapseSpec(); // Cannot compute collapse for dead triangle

		// Null out references only if strict memory management is needed and
		// ownership logic guarantees no external references remain. Might hide bugs.
		// Generally safer not to null out unless profiling shows benefit.
		// Arrays.fill(vertices, null);
		// Arrays.fill(neighbors, null);
		// Arrays.fill(wavefronts, null);
		// cachedCollapseSpec = null;
	}

	// --- Index Finding ---

	public int indexOfVertex(WavefrontVertex v) {
		if (v == null) {
			throw new NullPointerException("Cannot find index of null vertex in " + getName());
		}
		for (int i = 0; i < 3; i++) {
			if (vertices[i] == v) { // Object identity
				return i;
			}
		}
		throw new IllegalArgumentException("Vertex " + v.id + " not found in triangle " + getName());
	}

	public int indexOfNeighbor(KineticTriangle n) {
		if (n == null) {
			throw new NullPointerException("Cannot find index of null neighbor for " + getName());
		}
		for (int i = 0; i < 3; i++) {
			if (neighbors[i] == n) { // Object identity
				return i;
			}
		}
		throw new IllegalArgumentException("Neighbor " + n.getName() + " not found for triangle " + getName());
	}

	public int indexOfWavefront(WavefrontEdge edge) {
		if (edge == null) {
			throw new NullPointerException("Cannot find index of null wavefront edge for " + getName());
		}
		for (int i = 0; i < 3; i++) {
			if (wavefronts[i] == edge) { // Object identity
				return i;
			}
		}
		throw new IllegalArgumentException("Wavefront edge " + edge.id + " not found on triangle " + getName());
	}

	public boolean hasNeighbor(KineticTriangle n) {
		if (n == null) {
			return false;
		}
		return neighbors[0] == n || neighbors[1] == n || neighbors[2] == n;
	}

	public boolean hasVertex(WavefrontVertex v) {
		if (v == null) {
			return false;
		}
		return vertices[0] == v || vertices[1] == v || vertices[2] == v;
	}

	public boolean hasWavefront(WavefrontEdge e) {
		if (e == null) {
			return false;
		}
		return wavefronts[0] == e || wavefronts[1] == e || wavefronts[2] == e;
	}

	// --- Collapse Calculation ---

	public CollapseSpec getCollapseSpec(double currentTime) {
		if (isDead()) {
			// Explicitly handle dead case, return NEVER
			return CollapseSpec.NEVER;
		}
		if (!isCollapseSpecValid) {
			// --- Debug checks (optional) ---
			// assertVerticesMatchCache(); // Implement if using cache debug check

			cachedCollapseSpec = calculateCollapseSpec(currentTime);
			isCollapseSpecValid = true;

			// --- Debug state update (optional) ---
			// cacheCurrentVertices(); // Implement if using cache debug check
		}
		if (isDying()) {
			// If dying, the event might still need to be processed.
			// If cache is valid, return it. If not, compute it.
			// The event queue handler should ultimately decide based on isDying.
			// log.debug("getCollapseSpec called on dying triangle {}", getName()); // Use
			// logger
		}
		// --- Debug checks (optional) ---
		// assertVerticesMatchCache();

		if (cachedCollapseSpec == null) {
			// This should not happen if logic is correct, calculateCollapseSpec should
			// always return something (even NEVER)
			LOGGER.error("Cached collapse spec is null after calculation/cache hit for " + getName() + " at time " + currentTime);
			// Force recalculation or return NEVER as fallback?
			cachedCollapseSpec = calculateCollapseSpec(currentTime); // Try again
			if (cachedCollapseSpec == null) { // Still null? Major issue.
				return CollapseSpec.NEVER;
			}
		}
		return cachedCollapseSpec;
	}

	public CollapseSpec refineCollapseSpec(CollapseSpec refinedSpec) {
		if (isDead || isDying) {
			throw new IllegalStateException("Cannot refine collapse spec for dead/dying triangle " + getName());
		}
		if (!isCollapseSpecValid) {
			throw new IllegalStateException("Cannot refine collapse spec for " + getName() + ": current spec is not valid.");
		}
		if (cachedCollapseSpec == null) {
			// Should not happen if isCollapseSpecValid is true
			throw new IllegalStateException("Cannot refine collapse spec for " + getName() + ": cached spec is null despite being valid.");
		}
		if (refinedSpec == null) {
			throw new NullPointerException("Cannot refine collapse spec with null spec for " + getName());
		}

		// --- Debug checks (optional) ---
		// assertVerticesMatchCache();

		if (!cachedCollapseSpec.allowsRefinementTo(refinedSpec)) {
			throw new IllegalStateException("Refinement from " + cachedCollapseSpec + " to " + refinedSpec + " is not allowed for " + getName());
		}

		// log.debug("Refining collapse spec for {} from {} to {}", getName(),
		// cachedCollapseSpec.getType(), refinedSpec.getType()); // Use logger
		LOGGER.debug("Refining collapse spec for " + getName() + " from " + cachedCollapseSpec.getType() + " to " + refinedSpec.getType());
		this.cachedCollapseSpec = refinedSpec;
		// isCollapseSpecValid remains true
		return this.cachedCollapseSpec;
	}

	CollapseSpec calculateCollapseSpec(double currentTime) {
		// Added null check for vertices early on, might be needed if construction is
		// delayed
		if (vertices[0] == null || vertices[1] == null || vertices[2] == null) {
			LOGGER.warn("Calculating collapse for incomplete triangle " + getName());
			return CollapseSpec.NEVER; // Cannot calculate
		}

		final InfiniteSpeedType infiniteSpeedType = hasVertexInfiniteSpeed();
		if (infiniteSpeedType != InfiniteSpeedType.NONE) {
			// log.debug("KT{}: Calculating collapse - Has infinitely fast vertex ({})", id,
			// infiniteSpeedType); // Use logger
			LOGGER.debug("  KT" + id + ": Calculating collapse - Has infinitely fast vertex (" + infiniteSpeedType + ")");
			if (infiniteSpeedType == InfiniteSpeedType.OPPOSING) {
				return new CollapseSpec(CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING, currentTime, this);
			} else { // WEIGHTED
				Pair<Integer, Double> edgeInfo = findFastestWeightedEdgeInfo();
				int relevantEdge = edgeInfo.getLeft();
				if (relevantEdge != -1) {
					double relevantSpeedFactor = edgeInfo.getRight(); // Use the factor as secondary key
					return new CollapseSpec(CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED, currentTime, this, relevantEdge, relevantSpeedFactor);
				} else {
					LOGGER.error("KT" + id + " reported WEIGHTED infinite speed, but no relevant edge found.");
					return new CollapseSpec(CollapseType.INVALID_EVENT, currentTime, this);
				}
			}
		}

		if (isUnbounded()) {
			// log.debug("KT{}: Calculating collapse for UNBOUNDED triangle.", id); // Use
			// logger
			LOGGER.debug("  KT" + id + ": Calculating collapse for UNBOUNDED triangle.");
			return calculateCollapseUnbounded(currentTime);
		} else {
			// log.debug("KT{}: Calculating collapse for BOUNDED triangle.", id); // Use
			// logger
			LOGGER.debug("  KT" + id + ": Calculating collapse for BOUNDED triangle.");
			return calculateCollapseBounded(currentTime);
		}
	}

	private Pair<Integer, Double> findFastestWeightedEdgeInfo() {
		int relevantEdge = -1;
		double relevantEdgeSpeedFactor = Double.NEGATIVE_INFINITY;
		final double factor = SurfConstants.FASTER_EDGE_WINS_IN_COLLINEAR_CASES;

		for (int i = 0; i < 3; i++) {
			WavefrontEdge wfe = wavefronts[i];
			if (wfe == null) {
				continue;
			}

			WavefrontVertex v_ccw = getVertex(ccw(i));
			WavefrontVertex v_cw = getVertex(cw(i));

			// Ensure vertices exist
			if (v_ccw == null || v_cw == null) {
				continue;
			}

			boolean ccwIsWeighted = v_ccw.getInfiniteSpeed() == InfiniteSpeedType.WEIGHTED;
			boolean cwIsWeighted = v_cw.getInfiniteSpeed() == InfiniteSpeedType.WEIGHTED;

			if (!ccwIsWeighted && !cwIsWeighted) {
				continue;
			}

			double edgeSpeedFactor = wfe.getWeight() * factor;

			if (relevantEdge < 0 || edgeSpeedFactor > relevantEdgeSpeedFactor) {
				relevantEdge = i;
				relevantEdgeSpeedFactor = edgeSpeedFactor;
			}
		}
		return Pair.of(relevantEdge, relevantEdgeSpeedFactor);
	}

	// --- Bounded Triangle Calculation ---

	private CollapseSpec calculateCollapseBounded(double currentTime) {
		int numWavefronts = countWavefronts();
		// log.debug(" KT{} Bounded: Num Wavefronts: {}", id, numWavefronts); // Use
		// logger
		LOGGER.debug("    KT" + id + " Bounded: Num Wavefronts: " + numWavefronts);

		switch (numWavefronts) {
			case 3 :
				return calculateCollapseBoundedConstrained3(currentTime);
			case 2 :
				return calculateCollapseBoundedConstrained2(currentTime);
			case 1 :
				return calculateCollapseBoundedConstrained1(currentTime);
			case 0 :
				return calculateCollapseBoundedConstrained0(currentTime);
			default :
				throw new IllegalStateException("Invalid number of wavefronts: " + numWavefronts + " for KT" + id);
		}
	}

	private CollapseSpec calculateCollapseBoundedConstrained3(double currentTime) {
		// log.debug(" KT{} Calculating bounded constrained 3", id); // Use logger
		LOGGER.debug("      KT" + id + " Calculating bounded constrained 3");
		// Fetch collapse spec for one edge (assumes others match)
		if (wavefronts[0] == null) {
			throw new IllegalStateException("Constrained 3 triangle KT" + id + " missing wavefront 0."); // Should not happen
		}
		CollapseSpec cs0 = wavefronts[0].getCollapse(this.component, currentTime, 0);

		// Optional: Add assertions to verify other edges give the same result (within
		// tolerance)
		// CollapseSpec cs1 = wavefronts[1].getCollapse(this.component, currentTime, 1);
		// CollapseSpec cs2 = wavefronts[2].getCollapse(this.component, currentTime, 2);
		// assert cs0.isEquivalent(cs1, SurfConstants.ZERO_NT) && cs0.isEquivalent(cs2,
		// SurfConstants.ZERO_NT) : "Constrained 3 edges don't collapse together for KT"
		// + id;

		if (cs0.getType() == CollapseType.CONSTRAINT_COLLAPSE) {
			return new CollapseSpec(CollapseType.TRIANGLE_COLLAPSE, cs0.getTime(), this);
		} else if (cs0.getType() == CollapseType.NEVER) {
			return CollapseSpec.NEVER;
		} else {
			LOGGER.warn("Unexpected collapse type " + cs0.getType() + " from edge 0 in constrained 3 KT" + id);
			return CollapseSpec.NEVER;
		}
	}

	private CollapseSpec calculateCollapseBoundedConstrained2(double currentTime) {
		// log.debug(" KT{} Calculating bounded constrained 2", id); // Use logger
		LOGGER.debug("      KT" + id + " Calculating bounded constrained 2");
		int c1_idx = -1, c2_idx = -1;
		for (int i = 0; i < 3; i++) {
			if (wavefronts[i] != null) {
				if (c1_idx == -1) {
					c1_idx = i;
				} else {
					c2_idx = i;
				}
			}
		}
		if (c1_idx == -1 || c2_idx == -1) {
			throw new IllegalStateException("Constrained 2 triangle KT" + id + " has incorrect wavefront count.");
		}

		CollapseSpec cs1 = wavefronts[c1_idx].getCollapse(this.component, currentTime, c1_idx);
		CollapseSpec cs2 = wavefronts[c2_idx].getCollapse(this.component, currentTime, c2_idx);

		// log.debug(" KT{} Edge {} collapse: {}", id, c1_idx, cs1); // Use logger
		// log.debug(" KT{} Edge {} collapse: {}", id, c2_idx, cs2); // Use logger
		LOGGER.debug("        KT" + id + " Edge " + c1_idx + " collapse: " + cs1);
		LOGGER.debug("        KT" + id + " Edge " + c2_idx + " collapse: " + cs2);

		final boolean cs1ValidFuture = cs1.getType() == CollapseType.CONSTRAINT_COLLAPSE && cs1.getTime() >= currentTime - SurfConstants.ZERO_NT;
		final boolean cs2ValidFuture = cs2.getType() == CollapseType.CONSTRAINT_COLLAPSE && cs2.getTime() >= currentTime - SurfConstants.ZERO_NT;

		if (cs1ValidFuture && cs2ValidFuture) {
			if (Math.abs(cs1.getTime() - cs2.getTime()) < SurfConstants.ZERO_NT) {
				// log.debug(" KT{} Simultaneous collapse -> TRIANGLE_COLLAPSE @ {}", id,
				// cs1.getTime()); // Use logger
				LOGGER.debug("        KT" + id + " Simultaneous collapse -> TRIANGLE_COLLAPSE @ " + cs1.getTime());
				return new CollapseSpec(CollapseType.TRIANGLE_COLLAPSE, cs1.getTime(), this);
			} else {

				CollapseSpec earlier = ObjectUtils.min(cs1, cs2);
				// log.debug(" KT{} Returning earlier edge collapse: {}", id, earlier); // Use
				// logger
				LOGGER.debug("        KT" + id + " Returning earlier edge collapse: " + earlier);
				return earlier;
			}
		} else if (cs1ValidFuture) {
			// log.debug(" KT{} Returning edge 1 collapse: {}", id, cs1); // Use logger
			LOGGER.debug("        KT" + id + " Returning edge 1 collapse: " + cs1);
			return cs1;
		} else if (cs2ValidFuture) {
			// log.debug(" KT{} Returning edge 2 collapse: {}", id, cs2); // Use logger
			LOGGER.debug("        KT" + id + " Returning edge 2 collapse: " + cs2);
			return cs2;
		} else {
			// log.debug(" KT{} No future edge collapses -> NEVER", id); // Use logger
			LOGGER.debug("        KT" + id + " No future edge collapses -> NEVER");
			return CollapseSpec.NEVER;
		}
	}

	private CollapseSpec calculateCollapseBoundedConstrained1(double currentTime) {
		// log.debug(" KT{} Calculating bounded constrained 1", id); // Use logger
		LOGGER.debug("      KT" + id + " Calculating bounded constrained 1");
		int c_idx = -1;
		WavefrontEdge wf = null;
		for (int i = 0; i < 3; i++) {
			if (wavefronts[i] != null) {
				c_idx = i;
				wf = wavefronts[i];
				break;
			}
		}
		if (wf == null) {
			throw new IllegalStateException("Constrained 1 triangle KT" + id + " has no wavefront edge.");
		}

		// log.debug(" KT{} Constraint edge index: {} (Edge {})", id, c_idx, wf.id); //
		// Use logger
		LOGGER.debug("        KT" + id + " Constraint edge index: " + c_idx + " (Edge " + wf.id + ")");

		final Polynomial determinant = computeDeterminantPolynomial();
		// log.debug(" KT{} Determinant: {}", id, determinant); // Use logger
		LOGGER.debug("        KT" + id + " Determinant: " + determinant);
		// Check determinant sign at current time (should be >= 0 within tolerance)
		double detNow = determinant.evaluate(currentTime);
		if (detNow < -SurfConstants.ZERO_AREA_SQ) {
			LOGGER.warn("Triangle " + getName() + " has negative area " + detNow + " at time " + currentTime + ". (triangle's orientation fipped earlier?)");
			// Depending on strictness, could return INVALID or proceed cautiously.
			// can be ok, if event occurred in past
		}

		if (wf.parallelEndpoints(currentTime)) {
			// log.debug(" KT{} Endpoints are parallel.", id); // Use logger
			LOGGER.debug("        KT" + id + " Endpoints are parallel.");
			EdgeCollapseSpec edgeCollapse = wf.getEdgeCollapse(currentTime);
			if (edgeCollapse.getType() == EdgeCollapseType.ALWAYS) {
				// log.debug(" KT{} Edge collapses ALWAYS (now).", id); // Use logger
				LOGGER.debug("          KT" + id + " Edge collapses ALWAYS (now).");
				return CollapseSpec.fromEdgeCollapse(edgeCollapse, this, c_idx);
			} else { // NEVER
				// log.debug(" KT{} Edge collapses NEVER.", id); // Use logger
				LOGGER.debug("          KT" + id + " Edge collapses NEVER.");
				// Check if determinant suggests shrinkage -> split/flip
				if (determinant.getDegree() == 1 && determinant.b < 0) { // Check if linear and decreasing
					LOGGER.debug("          KT" + id + " Determinant linear & decreasing, checking split/flip.");
					return calculateSplitOrFlipEventBoundedConstrained1(currentTime, c_idx, determinant);
				} else {
					// log.debug(" KT{} Not linear & decreasing, returning NEVER.", id); // Use
					// logger
					LOGGER.debug("          KT" + id + " Not linear & decreasing, returning NEVER.");
					return CollapseSpec.NEVER;
				}
			}
		} else { // Endpoints not parallel
			// log.debug(" KT{} Endpoints are not parallel.", id); // Use logger
			LOGGER.debug("        KT" + id + " Endpoints are not parallel.");
			CollapseSpec candidate = wf.getCollapse(this.component, currentTime, c_idx);
			// log.debug(" KT{} Edge collapse candidate: {}", id, candidate); // Use logger
			LOGGER.debug("          KT" + id + " Edge collapse candidate: " + candidate);

			boolean candidateIsFuture = candidate.getType() == CollapseType.CONSTRAINT_COLLAPSE && candidate.getTime() >= currentTime - SurfConstants.ZERO_NT;

			if (determinant.getDegree() == 2) {
				// log.debug(" KT{} Determinant is quadratic.", id); // Use logger
				LOGGER.debug("          KT" + id + " Determinant is quadratic.");
				boolean acceptCandidate = false;
				if (candidateIsFuture) {
					acceptCandidate = acceptCollapseBoundedConstrained1(candidate.getTime(), determinant, true);
					// log.debug(" KT{} Candidate acceptable based on determinant roots? {}", id,
					// acceptCandidate); // Use logger
					LOGGER.debug("          KT" + id + " Candidate acceptable based on determinant roots? " + acceptCandidate);
				}

				if (acceptCandidate) {
					return candidate;
				} else {
					// log.debug(" KT{} Edge collapse rejected or invalid, checking split/flip.",
					// id); // Use logger
					LOGGER.debug("          KT" + id + " Edge collapse rejected or invalid, checking split/flip.");
					return calculateSplitOrFlipEventBoundedConstrained1(currentTime, c_idx, determinant);
				}
			} else { // Determinant degree <= 1
				// log.debug(" KT{} Determinant degree <= 1.", id); // Use logger
				LOGGER.debug("          KT" + id + " Determinant degree <= 1.");
				if (candidateIsFuture) {
					// log.debug(" KT{} Using edge collapse candidate.", id); // Use logger
					LOGGER.debug("          KT" + id + " Using edge collapse candidate.");
					return candidate;
				} else {
					// log.debug(" KT{} Edge collapse is past/never, checking split/flip.", id); //
					// Use logger
					LOGGER.debug("          KT" + id + " Edge collapse is past/never, checking split/flip.");
					return calculateSplitOrFlipEventBoundedConstrained1(currentTime, c_idx, determinant);
				}
			}
		}
	}

	private CollapseSpec calculateCollapseBoundedConstrained0(double currentTime) {
		// log.debug(" KT{} Calculating bounded constrained 0 (Flip/Spoke/Meet)", id);
		// // Use logger
		LOGGER.debug("      KT" + id + " Calculating bounded constrained 0 (Flip/Spoke/Meet)");
		final Polynomial determinant = computeDeterminantPolynomial();
		// log.debug(" KT{} Determinant: {}", id, determinant); // Use logger
		LOGGER.debug("        KT" + id + " Determinant: " + determinant);
		// Check determinant sign at current time
		double detNow = determinant.evaluate(currentTime);
		if (detNow < -SurfConstants.ZERO_AREA_SQ) {
			LOGGER.warn("Triangle " + getName() + " has negative area " + detNow + " at time " + currentTime + ". (triangle's orientation fipped earlier?)");
		}

		return calculateFlipEvent(currentTime, determinant);
	}

	// --- Unbounded Triangle Calculation ---

	private CollapseSpec calculateCollapseUnbounded(double currentTime) {
		// log.debug(" KT{} Calculating unbounded", id); // Use logger
		LOGGER.debug("      KT" + id + " Calculating unbounded");
		if (!isUnbounded()) {
			throw new IllegalStateException("calculateCollapseUnbounded called on bounded triangle " + getName());
		}
		if (hasVertexInfiniteSpeed() != InfiniteSpeedType.NONE) {
			throw new IllegalStateException("calculateCollapseUnbounded called on triangle " + getName() + " with infinite speed vertex");
		}

		final int infIdx = getInfiniteVertexIndex();
		CollapseSpec edgeCollapse = CollapseSpec.NEVER;

		// 1. Check for collapse of the triangle's own finite edge (if constrained)
		if (isConstrained(infIdx)) {
			WavefrontEdge boundedEdge = wavefronts[infIdx];
			if (boundedEdge == null) {
				throw new IllegalStateException("Unbounded KT" + id + " constraint flag mismatch edge " + infIdx);
			}
			edgeCollapse = boundedEdge.getCollapse(this.component, currentTime, infIdx);
			// log.debug(" KT{} Bounded edge {} (Edge {}) collapse: {}", id, infIdx,
			// boundedEdge.id, edgeCollapse); // Use logger
			LOGGER.debug("        KT" + id + " Bounded edge " + infIdx + " (Edge " + (boundedEdge != null ? boundedEdge.id : "null") + ") collapse: "
					+ edgeCollapse);
		} else {
			// log.debug(" KT{} Bounded edge {} is not constrained (spoke).", id, infIdx);
			// // Use logger
			LOGGER.debug("        KT" + id + " Bounded edge " + infIdx + " is not constrained (spoke).");
		}

		// 2. Check if CCW vertex from neighbor leaves the convex hull boundary
		final int neighborEdgeIdx = ccw(infIdx);
		final KineticTriangle neighbor = neighbors[neighborEdgeIdx];
		if (neighbor == null) {
			throw new IllegalStateException("Unbounded triangle " + getName() + " missing neighbor on edge " + neighborEdgeIdx);
		}
		if (!neighbor.isUnbounded()) {
			throw new IllegalStateException("Neighbor " + neighbor.getName() + " of unbounded " + getName() + " must be unbounded.");
		}

		final int neighborInfIdx = neighbor.getInfiniteVertexIndex();
		// Get the common finite vertex 'v' (vertex opposite the edge towards the
		// neighbor in this triangle)
		final WavefrontVertex v = getVertex(cw(infIdx));
		if (v == null || v.isInfinite()) {
			throw new IllegalStateException("Invalid common finite vertex v for unbounded KT" + id);
		}
		// Verification: Check if 'v' is the vertex opposite the edge towards this
		// triangle in the neighbor
		if (v != neighbor.getVertex(ccw(neighborInfIdx))) {
			throw new IllegalStateException("Finite vertex mismatch between unbounded neighbors " + getName() + " and " + neighbor.getName());
		}

		// Get the vertex 'w' from the neighbor (the one that might leave the CH)
		final WavefrontVertex w = neighbor.getVertex(cw(neighborInfIdx)); // Vertex 'w' in neighbor (ccw from vInf->v edge)
		if (w == null || w.isInfinite()) {
			throw new IllegalStateException("Invalid vertex w for unbounded neighbor " + neighbor.getName());
		}

		CollapseSpec vertexLeavesCH = CollapseSpec.NEVER;

		// Check if the boundary edge (either in this triangle or the neighbor) is
		// constrained
		if (isConstrained(infIdx) || neighbor.isConstrained(neighborInfIdx)) {
			// --- Constrained Boundary Case ---
			final WavefrontEdge edgeE;
			final WavefrontVertex definingVertexForLog; // Vertex used to determine which edge is 'edgeE'

			if (isConstrained(infIdx)) {
				edgeE = wavefronts[infIdx];
				if (edgeE == null) {
					throw new IllegalStateException("Unbounded KT" + id + " constraint flag mismatch edge " + infIdx);
				}
				definingVertexForLog = neighbor.getVertex(ccw(neighborInfIdx)); // This is 'v'
				LOGGER.debug(" KT" + id + " Using constrained edge from this triangle (Edge " + edgeE.id + ")");

			} else { // neighbor is constrained
				edgeE = neighbor.wavefronts[neighborInfIdx];
				if (edgeE == null) {
					throw new IllegalStateException("Unbounded neighbor KT" + neighbor.id + " constraint flag mismatch edge " + neighborInfIdx);
				}
				definingVertexForLog = getVertex(cw(infIdx)); // This is 'v'
				LOGGER.debug(" KT" + id + " Using constrained edge from neighbor (Edge " + edgeE.id + ")");
			}
			LOGGER.debug("      (Constrained check uses common vertex V" + v.id + " against Edge " + edgeE.id + ", following C++ logic)");

			final WavefrontSupportingLine lineE = edgeE.getSupportingLine();

			// *** KEY CHANGE: Use common vertex 'v' for checks, following C++ logic ***
			LOGGER.debug(" KT" + id + " Checking relative speed of common vertex V" + v.id + " and constrained edge (Edge " + edgeE.id + ")");
			final Sign edgeFaster = edgeIsFasterThanVertex(v, lineE); // Use common vertex 'v'

			if (edgeFaster == Sign.POSITIVE) {
				// C++ logic: If edge pulls away from common vertex 'v' faster, assume 'w' won't
				// cross inwards.
				LOGGER.debug(" Edge " + edgeE.id + " is faster than common vertex V" + v.id + " (along normal). Treating as NEVER leaves CH.");
				vertexLeavesCH = CollapseSpec.NEVER; // Match C++ behavior
			} else if (edgeFaster == Sign.ZERO) {
				// Parallel motion. Check if already collinear (distance is zero).
				LOGGER.debug(" Edge " + edgeE.id + " and common vertex V" + v.id + " have same speed along normal (parallel).");
				VertexOnSupportingLineResult hitResult = getTimeVertexOnSupportingLine(v, lineE); // Use 'v'
				if (hitResult.type == VertexOnSupportingLineType.ALWAYS) {
					LOGGER.debug(" Common vertex V" + v.id + " is always on the line (collinear). Treating as NEVER leaves CH.");
				} else { // NEVER type
					LOGGER.debug(" Common vertex V" + v.id + " is never on the line (parallel, offset). Treating as NEVER leaves CH.");
				}
				vertexLeavesCH = CollapseSpec.NEVER; // Match C++ behavior for ZERO sign
			} else { // edgeFaster == Sign.NEGATIVE
				// C++ logic: Edge is NOT pulling away faster than common vertex 'v'.
				// Calculate time when 'v' would hit the line 'lineE'.
				// This time is attributed to the CCW_VERTEX_LEAVES_CH event.
				LOGGER.debug(" Common vertex V" + v.id + " is faster than or equal speed to edge " + edgeE.id + " (along normal).");
				VertexOnSupportingLineResult hitResult = getTimeVertexOnSupportingLine(v, lineE); // Use 'v'

				if (hitResult.type == VertexOnSupportingLineType.ONCE && hitResult.time >= currentTime - SurfConstants.ZERO_NT) {
					// Check if the event V hits Line E occurs in the future
					LOGGER.debug(" Common vertex V" + v.id + " hits line of edge " + edgeE.id + " at future time: " + hitResult.time
							+ ". Creating CCW_VERTEX_LEAVES_CH event.");
					// Assign this time to the 'w' leaving event, using the original triangle's
					// infIdx
					vertexLeavesCH = new CollapseSpec(CollapseType.CCW_VERTEX_LEAVES_CH, hitResult.time, this, infIdx);
				} else {
					LOGGER.debug(" Common vertex V" + v.id + " hits line of edge " + edgeE.id + " in past or never/always. Treating as NEVER leaves CH.");
					// vertexLeavesCH remains NEVER
				}
			}
		} else {
			// --- Unconstrained Boundary (Spoke) Case ---
			// This logic remains the same, using the determinant of u, v, w
			final WavefrontVertex u = getVertex(ccw(infIdx)); // Other finite vertex of this triangle
			if (u == null || u.isInfinite()) {
				throw new IllegalStateException("Invalid vertex u for unbounded KT" + id);
			}
			LOGGER.debug(" KT" + id + " Checking if vertex w (V" + w.id + ") leaves CH across spoke (u=V" + u.id + ", v=V" + v.id + ")");
			boolean sameVel = u.getVelocity().equals(v.getVelocity()) && v.getVelocity().equals(w.getVelocity());
			if (sameVel) {
				// Match C++ logic: Check for coincident vertices. If not coincident, result is
				// NEVER.
				// Note: C++ checks pos_zero. Using getInitialPosition() assuming it's
				// equivalent.
				boolean coincident = u.getInitialPosition().equals(v.getInitialPosition()) || v.getInitialPosition().equals(w.getInitialPosition());
				// C++ doesn't check u == w, but might be implied? Added for completeness.
				// || u.getInitialPosition().equals(w.getInitialPosition());

				if (coincident) {
					// C++ aborts here. We should probably throw or handle distinct from NEVER.
					LOGGER.warn("Vertices u,v,w have same velocity AND are coincident. C++ would abort.");
					// Decide on appropriate error handling or specific CollapseSpec type if needed.
					// For now, mirroring the 'NEVER' from the non-coincident case, but this is
					// ambiguous.
					vertexLeavesCH = CollapseSpec.NEVER; // Or throw new IllegalStateException("...");
				} else {
					LOGGER.debug(" Vertices u,v,w have same velocity, not coincident -> NEVER leaves CH.");
					vertexLeavesCH = CollapseSpec.NEVER; // This matches the C++ 'else' branch for non-coincident.
				}
			} else {
				// Velocities differ, calculate time of collinearity using determinant
				final Polynomial determinantUVw = computeDeterminantFromVertices(u, v, w);
				Optional<Double> collapseTimeOpt = getGenericCollapseTime(currentTime, determinantUVw);
				if (collapseTimeOpt.isPresent()) {
					double collapseTime = collapseTimeOpt.get();
					LOGGER.debug(" Vertices u,v,w become collinear at time: " + collapseTime);
					vertexLeavesCH = new CollapseSpec(CollapseType.CCW_VERTEX_LEAVES_CH, collapseTime, this, infIdx);
				} else {
					LOGGER.debug(" Vertices u,v,w never become collinear in the future.");
					// vertexLeavesCH remains NEVER
				}
			}
		}

		// log.debug(" KT{} Edge collapse: {}", id, edgeCollapse); // Use logger
		// log.debug(" KT{} Vertex leaves CH: {}", id, vertexLeavesCH); // Use logger
		LOGGER.debug("        KT" + id + " Final Edge collapse: " + edgeCollapse);
		LOGGER.debug("        KT" + id + " Final Vertex leaves CH: " + vertexLeavesCH);

		// Return the earliest of the two possible events
		return ObjectUtils.min(edgeCollapse, vertexLeavesCH);
	}

	public double cross(Vector2D v1, Vector2D v2) {
		return v1.getX() * v2.getY() - v1.getY() * v2.getX();
	}

	// --- Helper Methods for Collapse Calculation ---

	private int countWavefronts() {
		int count = 0;
		if (wavefronts[0] != null) {
			count++;
		}
		if (wavefronts[1] != null) {
			count++;
		}
		if (wavefronts[2] != null) {
			count++;
		}
		return count;
	}

	/**
	 * Computes the polynomial representing twice the signed area of the triangle as
	 * a function of time. det(t) = a*t^2 + b*t + c. Derived from the determinant
	 * formula: det = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0) where xi(t) = xi(0) + vxi *
	 * t, etc.
	 *
	 * @return The Polynomial for 2 * signed area.
	 * @throws IllegalStateException if any vertex is null.
	 */
	Polynomial computeDeterminantPolynomial() {
		if (vertices[0] == null || vertices[1] == null || vertices[2] == null) {
			throw new IllegalStateException("Cannot compute determinant for incomplete triangle " + getName());
		}
		// Note: This assumes finite vertices or that WavefrontVertex handles infinite
		// cases gracefully (e.g., zero velocity).
		// If infinite vertices have special non-zero "velocities", this needs
		// adjustment.
		// Based on the provided getPositionAt, infinite vertices have zero velocity.

		// Initial positions
		final Coordinate p0 = vertices[0].getInitialPosition();
		final Coordinate p1 = vertices[1].getInitialPosition();
		final Coordinate p2 = vertices[2].getInitialPosition();

		// Velocities
		final Coordinate v0 = vertices[0].getVelocity();
		final Coordinate v1 = vertices[1].getVelocity();
		final Coordinate v2 = vertices[2].getVelocity();

		// Calculate differences in initial positions
		final double dx10 = p1.x - p0.x;
		final double dy10 = p1.y - p0.y;
		final double dx20 = p2.x - p0.x;
		final double dy20 = p2.y - p0.y;

		// Calculate differences in velocities
		final double dvx10 = v1.getX() - v0.getX();
		final double dvy10 = v1.getY() - v0.getY();
		final double dvx20 = v2.getX() - v0.getX();
		final double dvy20 = v2.getY() - v0.getY();

		// Calculate coefficients a, b, c for det(t) = a*t^2 + b*t + c
		// a = (dvx10 * dvy20) - (dvx20 * dvy10)
		final double a = (dvx10 * dvy20) - (dvx20 * dvy10);

		// b = (dx10*dvy20 + dvx10*dy20) - (dx20*dvy10 + dvx20*dy10)
		final double b = (dx10 * dvy20 + dvx10 * dy20) - (dx20 * dvy10 + dvx20 * dy10);

		// c = (dx10 * dy20) - (dx20 * dy10) (This is 2 * signed area at t=0)
		final double c = (dx10 * dy20) - (dx20 * dy10);

		return new Polynomial(a, b, c);
	}

	/** Static version to compute determinant polynomial from three vertices. */
	public static Polynomial computeDeterminantFromVertices(WavefrontVertex vert0, WavefrontVertex vert1, WavefrontVertex vert2) {
		if (vert0 == null || vert1 == null || vert2 == null) {
			throw new NullPointerException("Cannot compute determinant from null vertices.");
		}
		// Initial positions
		final Coordinate p0 = vert0.getInitialPosition();
		final Coordinate p1 = vert1.getInitialPosition();
		final Coordinate p2 = vert2.getInitialPosition();

		// Velocities
		final Coordinate v0 = vert0.getVelocity();
		final Coordinate v1 = vert1.getVelocity();
		final Coordinate v2 = vert2.getVelocity();

		final double dx10 = p1.x - p0.x;
		final double dy10 = p1.y - p0.y;
		final double dx20 = p2.x - p0.x;
		final double dy20 = p2.y - p0.y;

		final double dvx10 = v1.getX() - v0.getX();
		final double dvy10 = v1.getY() - v0.getY();
		final double dvx20 = v2.getX() - v0.getX();
		final double dvy20 = v2.getY() - v0.getY();

		final double a = (dvx10 * dvy20) - (dvx20 * dvy10);
		final double b = (dx10 * dvy20 + dvx10 * dy20) - (dx20 * dvy10 + dvx20 * dy10);
		final double c = (dx10 * dy20) - (dx20 * dy10);

		return new Polynomial(a, b, c);
	}

	/**
	 * Calculates the time of the next generic collapse event (determinant becoming
	 * zero). Ported logic from C++ implementation.
	 *
	 * @param currentTime The current time.
	 * @param det         The polynomial representing the determinant (area * 2)
	 *                    over time.
	 * @return An Optional containing the time of the *first* future collapse event,
	 *         or Optional.empty() if no collapse occurs at or after currentTime.
	 */
	private Optional<Double> getGenericCollapseTime(double currentTime, Polynomial det) {

		final int degree = det.getDegree();
		final double detNow = det.evaluate(currentTime); // Evaluate determinant at current time

		// --- Degree 0 ---
		if (degree == 0) {
			LOGGER.warn("Have a polynomial of degree zero. Can we catch this sooner? det={}, time={}", det.c, currentTime);
			if (Math.abs(det.c) < SurfConstants.ZERO_AREA_SQ) { // Check if constant is zero
				LOGGER.warn("Degree 0 polynomial is zero -> collapses now (and always). time={}", currentTime);

				// C++ loop detection logic
				if (currentTime == lastTime) {
					loopDetectionCount++;
					if (loopDetectionCount > LOOP_DETECTION_LIMIT) {
						LOGGER.error("Potential infinite loop detected at time {} (Degree 0 Zero). Aborting check.", currentTime);
						// Consider throwing an exception in critical systems instead of just returning
						return Optional.empty(); // Avoid returning a time that might perpetuate loop
					}
				} else {
					loopDetectionCount = 0;
					lastTime = currentTime;
				}
				// End C++ loop detection logic

				return Optional.of(currentTime); // Collapse happens now (or always)
			} else {
				LOGGER.debug("        Degree 0, non-zero constant -> NEVER");
				return Optional.empty(); // Non-zero constant, never collapses
			}
		}

		// --- Degree 1 ---
		else if (degree == 1) {
			final double b = det.b;
			final double c = det.c;

			if (Math.abs(b) < SurfConstants.ZERO_POLYNOMIAL_COEFF) {
				// Should not happen if degree is truly 1, but handle degeneracies
				LOGGER.warn("Degree 1 polynomial has near-zero slope b={}, c={}. Treating as degree 0.", b, c);
				// Re-evaluate as degree 0
				if (Math.abs(c) < SurfConstants.ZERO_AREA_SQ) {
					LOGGER.warn("Effective Degree 0 is zero -> collapses now (and always). time={}", currentTime);
					// Add loop detection here too if desired, similar to degree 0 case
					return Optional.of(currentTime);
				} else {
					return Optional.empty();
				}
			}

			final double root = -c / b; // Calculate the single root
			LOGGER.debug("        GenericCollapseTime: Degree 1, root={}", root);

			// Compare root with currentTime using tolerance ZERO_NT
			if (Math.abs(root - currentTime) < SurfConstants.ZERO_NT) { // Root is effectively now
				// C++ logic: collapse now only if sign is negative (slope b < 0)
				if (b < -SurfConstants.ZERO_POLYNOMIAL_COEFF) { // Check sign of b (slope)
					LOGGER.warn("Polynomial (degree 1) has a zero right now and decreasing. Collapse NOW. time={}", currentTime);

					// C++ loop detection logic
					if (currentTime == lastTime) {
						loopDetectionCount++;
						if (loopDetectionCount > LOOP_DETECTION_LIMIT) {
							LOGGER.error("Potential infinite loop detected at time {} (Degree 1 Zero Decreasing). Aborting check.", currentTime);
							return Optional.empty();
						}
					} else {
						loopDetectionCount = 0;
						lastTime = currentTime;
					}
					// End C++ loop detection logic

					return Optional.of(currentTime); // Collapse now because it's decreasing
				} else {
					LOGGER.debug("          Root is now, but slope is non-negative (b={}) -> No future event", b);
					return Optional.empty(); // Not collapsing, or passing through zero while increasing/flat
				}
			} else if (root > currentTime /*- SurfConstants.ZERO_NT removed, compare directly */) { // Root is strictly in the future
				LOGGER.debug("          Root is future -> Event at {}", root);
				// Optional sanity check (sign mismatch) - useful for debugging state issues
				if (detNow * b >= SurfConstants.ZERO_AREA_SQ * SurfConstants.ZERO_POLYNOMIAL_COEFF && Math.abs(detNow) > SurfConstants.ZERO_AREA_SQ) {
					LOGGER.warn("Sign mismatch? Future linear root={}, but detNow={} and slope={} have same sign.", root, detNow, b);
				}
				return Optional.of(root); // Future collapse
			} else { // Root is in the past
				LOGGER.debug("          Root is past -> No future event");
				return Optional.empty();
			}
		}

		// --- Degree 2 ---
		else { // degree == 2
			LOGGER.debug("        GenericCollapseTime: Degree 2: {}", det);

			// Use the revised solver that returns ALL real roots
			double[] roots = QuadraticSolver.solve(det); // Returns 0, 1, or 2 roots, sorted.

			LOGGER.debug("          Mathematical Roots: {}", Arrays.toString(roots));

			if (roots.length == 0) {
				LOGGER.debug("          No real roots.");
				// Optional sanity check: if no roots, determinant should not cross zero.
				// If a>0, det should always be positive (if starting positive).
				// If a<0, det should always be negative (if starting negative).
				if (det.a > SurfConstants.ZERO_POLYNOMIAL_COEFF && detNow < -SurfConstants.ZERO_AREA_SQ) {
					LOGGER.warn("Quadratic opens up (a>0), but detNow is negative ({}) with no real roots!", detNow);
				} else if (det.a < -SurfConstants.ZERO_POLYNOMIAL_COEFF && detNow > SurfConstants.ZERO_AREA_SQ) {
					LOGGER.warn("Quadratic opens down (a<0), but detNow is positive ({}) with no real roots!", detNow);
				}
				return Optional.empty(); // No roots means no collapse
			} else {
				// We have one or two roots. Let's call them r0 and r1 (r0 <= r1).
				// If only one root, r0 = r1.
				double r0 = roots[0];
				double r1 = (roots.length == 1) ? r0 : roots[1];

				double a = det.a;
				Optional<Double> collapseTime = Optional.empty();

				// Apply C++ logic based on sign of 'a'
				if (a < -SurfConstants.ZERO_POLYNOMIAL_COEFF) { // Parabola opens downward (a < 0)
					// C++ logic: collapse time is the *larger* root (x1).
					// We only report it if it's a *future* event relative to currentTime.
					LOGGER.debug("          Parabola opens down (a < 0). Considering larger root r1 = {}", r1);
					if (r1 >= currentTime - SurfConstants.ZERO_NT) {
						// Clamp to currentTime if it's very close to avoid past times due to precision
						collapseTime = Optional.of(Math.max(currentTime, r1));
						LOGGER.debug("          Selected future collapse time (a<0): {}", collapseTime.get());
					} else {
						LOGGER.debug("          Larger root r1 is in the past. No future collapse.");
					}
				} else if (a > SurfConstants.ZERO_POLYNOMIAL_COEFF) { // Parabola opens upward (a > 0)
					// C++ logic: collapse time is the *smaller* root (x0), but ONLY if it's >=
					// time_now.
					LOGGER.debug("          Parabola opens up (a > 0). Considering smaller root r0 = {}", r0);
					if (r0 >= currentTime - SurfConstants.ZERO_NT) {
						// Clamp to currentTime if very close
						collapseTime = Optional.of(Math.max(currentTime, r0));
						LOGGER.debug("          Selected future collapse time (a>0): {}", collapseTime.get());
					} else {
						// r0 is in the past. The next root r1 is when the determinant becomes positive
						// again,
						// so it's not considered a "collapse" in the C++ logic's context here.
						LOGGER.debug("          Smaller root r0 is in the past. No future collapse.");
						// We also need to check if r1 is a future time, because if r0 is past and r1 is
						// future,
						// we are currently *inside* the interval where det < 0 (if detNow < 0).
						// However, the C++ code explicitly says `result = false` if x0 < time_now when
						// a > 0.
						// We will adhere to that direct port. If the system requires handling the exit
						// from a
						// negative determinant state when a>0, the logic needs further adaptation
						// beyond the C++ source.
					}
				} else {
					// a is effectively zero. Should have been caught by degree check.
					LOGGER.error("Polynomial has degree 2 but 'a' coefficient is near zero: {}", det);
					// Fallback or return empty ? Returning empty is safer.
					return Optional.empty();
				}

				// Final sanity check (optional but good):
				// If we found a collapse time, detNow should ideally have the "correct" sign.
				// If a < 0, we expect detNow >= 0 for a future collapse at r1.
				// If a > 0, we expect detNow >= 0 for a future collapse at r0.
				if (collapseTime.isPresent()) {
					if (a < -SurfConstants.ZERO_POLYNOMIAL_COEFF && detNow < -SurfConstants.ZERO_AREA_SQ) {
						LOGGER.warn("Logic Warning (a<0): Found future collapse time {} but detNow={} is negative.", collapseTime.get(), detNow);
					} else if (a > SurfConstants.ZERO_POLYNOMIAL_COEFF && detNow < -SurfConstants.ZERO_AREA_SQ) {
						LOGGER.warn("Logic Warning (a>0): Found future collapse time {} but detNow={} is negative.", collapseTime.get(), detNow);
					}
				} else {
					// If no future collapse time was found, log why
					if (roots.length == 1)
						LOGGER.debug("          Single root r0={} is in the past.", r0);
					// else handled above (r1 past for a<0, r0 past for a>0)
				}

				return collapseTime;
			}
		}
	} // end getGenericCollapseTime

	private boolean acceptCollapseBoundedConstrained1(double collapseTime, Polynomial determinant, boolean isEdgeCollapseCandidate) {
		if (determinant.getDegree() != 2) {
			LOGGER.warn("acceptCollapseBoundedConstrained1 called with non-quadratic determinant for " + getName());
			// Fallback for non-quadratic: If it's a future time, accept it? Or always
			// reject?
			// Let's return true if it's a future edge collapse, mirroring linear case.
			return isEdgeCollapseCandidate; // Accept any future edge collapse if determinant isn't quadratic
		}

		final double a = determinant.a;
		int signLead = (int) Math.signum(a);
		if (Math.abs(a) < SurfConstants.ZERO_POLYNOMIAL_COEFF) { // close to 0 -> 0
			signLead = 0;
		} // Handle degenerate quadratic

		if (signLead < 0) { // Opens down
			// log.debug(" AcceptCollapse: a < 0 -> Accepting"); // Use logger
			LOGGER.debug("          AcceptCollapse: a < 0 -> Accepting");
			return true;
		} else if (signLead == 0) { // Degenerate (linear/constant)
			LOGGER.warn("acceptCollapse logic encountered degenerate quadratic for " + getName());
			return isEdgeCollapseCandidate; // Accept if edge collapse, reject otherwise?
		} else { // signLead > 0, Opens up.
			// log.debug(" AcceptCollapse: a > 0 -> Checking derivative"); // Use logger
			LOGGER.debug("          AcceptCollapse: a > 0 -> Checking derivative");
			final Polynomial derivative = determinant.differentiate();
			final double derivAtCollapse = derivative.evaluate(collapseTime);
			// log.debug(" Derivative at collapse time {}: {}", collapseTime,
			// derivAtCollapse); // Use logger
			LOGGER.debug("            Derivative at collapse time " + collapseTime + ": " + derivAtCollapse);

			if (Math.abs(derivAtCollapse) < SurfConstants.ZERO_POLYNOMIAL_VALUE) { // Derivative is zero
				// log.debug(" Derivative is zero. Accepting only if edge collapse: {}",
				// isEdgeCollapseCandidate); // Use logger
				LOGGER.debug("            Derivative is zero. Accepting only if edge collapse: " + isEdgeCollapseCandidate);
				return isEdgeCollapseCandidate;
			} else if (derivAtCollapse < 0) { // Derivative is negative (first root)
				// log.debug(" Derivative is negative -> Accepting (first root)"); // Use logger
				LOGGER.debug("            Derivative is negative -> Accepting (first root)");
				return true;
			} else { // Derivative is positive (second root)
				// log.debug(" Derivative is positive -> Rejecting (second root)"); // Use
				// logger
				LOGGER.debug("            Derivative is positive -> Rejecting (second root)");
				return false;
			}
		}
	}

	private CollapseSpec calculateSplitOrFlipEventBoundedConstrained1(double currentTime, int c_idx, Polynomial determinant) {
		// log.debug(" KT{} Calculating Split/Flip for Constrained 1 edge {}", id,
		// c_idx); // Use logger
		LOGGER.debug("          KT" + id + " Calculating Split/Flip for Constrained 1 edge " + c_idx);

		boolean anyReflex = false;
		for (int i = 0; i < 3; i++) {
			if (vertices[i] != null && vertices[i].isReflexOrStraight()) {
				anyReflex = true;
				break;
			}
		}
		if (!anyReflex) {
			// log.debug(" All vertices convex, returning NEVER.", id); // Use logger
			LOGGER.debug("            All vertices convex, returning NEVER.");
			return eventThatWillNotHappen(currentTime, determinant);
		}

		final WavefrontVertex vOpposite = vertices[c_idx];
		final WavefrontEdge edge = wavefronts[c_idx];
		if (vOpposite == null || edge == null) {
			throw new IllegalStateException("Missing vertex/edge for split/flip calc in KT" + id);
		}
		final WavefrontSupportingLine line = edge.getSupportingLine();

		final VertexOnSupportingLineResult hitResult = getTimeVertexOnSupportingLine(vOpposite, line);
		// log.debug(" Vertex V{} hits supporting line of Edge {}: {}", id,
		// vOpposite.id, edge.id, hitResult); // Use logger
		LOGGER.debug("            Vertex V" + vOpposite.id + " hits supporting line of Edge " + edge.id + ": " + hitResult);

		switch (hitResult.type) {
			case ONCE :
				final double hitTime = hitResult.time; // aka collapse_time
				if (hitTime >= currentTime - SurfConstants.ZERO_NT) { // Hit is now or future
					if (Math.abs(hitTime - currentTime) < SurfConstants.ZERO_NT) { // Hit is now
						// log.debug(" Hit time is now.", id); // Use logger
						LOGGER.debug("            Hit time is now.");
						// Check if shrinking (using acceptCollapse logic without edge flag)
						boolean shrinkingOrStagnant = false;
						if (determinant.getDegree() == 2) {
							shrinkingOrStagnant = acceptCollapseBoundedConstrained1(currentTime, determinant, false); // false -> check derivative only
						} else if (determinant.getDegree() == 1) {
							shrinkingOrStagnant = determinant.b <= SurfConstants.ZERO_POLYNOMIAL_COEFF; // negative -- Shrinking or flat
						} else { // Constant degree 0
							shrinkingOrStagnant = Math.abs(determinant.c) < SurfConstants.ZERO_AREA_SQ; // Only if already zero
						}

						if (shrinkingOrStagnant) {
							// log.debug(" Triangle shrinking/stagnant now -> SPLIT_OR_FLIP_REFINE @ {}",
							// id, currentTime); // Use logger
							LOGGER.debug("            Triangle shrinking/stagnant now -> SPLIT_OR_FLIP_REFINE @ " + currentTime);
							return new CollapseSpec(CollapseType.SPLIT_OR_FLIP_REFINE, currentTime, this, c_idx);
						} else {
							// log.debug(" Triangle expanding now -> NEVER", id); // Use logger
							LOGGER.debug("            Triangle expanding now -> NEVER");
							return eventThatWillNotHappen(currentTime, determinant);
						}
					} else { // Hit is strictly future
						// log.debug(" Hit time is in the future -> SPLIT_OR_FLIP_REFINE @ {}", id,
						// hitTime); // Use logger
						LOGGER.debug("            Hit time is in the future -> SPLIT_OR_FLIP_REFINE @ " + hitTime);
						return new CollapseSpec(CollapseType.SPLIT_OR_FLIP_REFINE, hitTime, this, c_idx);
					}
				} else { // hitTime < currentTime (past)
					// log.debug(" Hit time is in the past -> NEVER", id); // Use logger
					LOGGER.debug("            Hit time is in the past -> NEVER");
					return eventThatWillNotHappen(currentTime, determinant);
				}
			case ALWAYS :
				// log.debug(" Vertex always on line -> SPLIT_OR_FLIP_REFINE @ {}", id,
				// currentTime); // Use logger
				LOGGER.debug("            Vertex always on line -> SPLIT_OR_FLIP_REFINE @ " + currentTime);
				return new CollapseSpec(CollapseType.SPLIT_OR_FLIP_REFINE, currentTime, this, c_idx);
			case NEVER :
				// log.debug(" Vertex never hits line -> NEVER", id); // Use logger
				LOGGER.debug("            Vertex never hits line -> NEVER");
				return new CollapseSpec(CollapseType.NEVER, currentTime, this, c_idx);
//				return eventThatWillNotHappen(currentTime, determinant); // NOTE
			default :
				throw new AssertionError("Unknown VertexOnSupportingLineType: " + hitResult.type);
		}
	}

	private CollapseSpec calculateFlipEvent(double currentTime, Polynomial determinant) {
		// log.debug(" KT{} Calculating Flip/Spoke Event", id); // Use logger
		LOGGER.debug("        KT" + id + " Calculating Flip/Spoke Event");

		boolean potentiallyCouldFlip = false;
		for (int i = 0; i < 3; i++) {
			if (vertices[i] != null && vertices[i].isReflexOrStraight()) {
				potentiallyCouldFlip = true;
				break;
			}
		}
		if (!potentiallyCouldFlip) {
			// log.debug(" No reflex/straight vertices, cannot flip. Returning NEVER.", id);
			// // Use logger
			LOGGER.debug("          No reflex/straight vertices, cannot flip. Returning NEVER.");
			return eventThatWillNotHappen(currentTime, determinant);
		}

		CollapseSpec genericCollapse = getGenericCollapse(currentTime, determinant);
		// log.debug(" Generic collapse result: {}", id, genericCollapse); // Use logger
		LOGGER.debug("          Generic collapse result: " + genericCollapse);

		switch (genericCollapse.getType()) {
			case NEVER :
				return CollapseSpec.NEVER; // Already handled by eventThatWillNotHappen if needed
			case TRIANGLE_COLLAPSE :
			case SPOKE_COLLAPSE :
				// log.info(" Flip calculation resulted in Triangle/Spoke collapse: {}",
				// genericCollapse); // Use info level?
				LOGGER.debug("          Flip calculation resulted in Triangle/Spoke collapse: " + genericCollapse);
				return genericCollapse;

			case VERTEX_MOVES_OVER_SPOKE :
				final int movingVertexIndex = genericCollapse.getRelevantEdge();
				final WavefrontVertex movingVertex = vertices[movingVertexIndex];
				if (movingVertex == null) {
					throw new IllegalStateException("Null moving vertex " + movingVertexIndex + " during flip check for " + getName());
				}

				if (movingVertex.isReflexOrStraight()) {
					// log.debug(" Reflex vertex V{} moves over spoke edge opposite index {}. Valid
					// flip/spoke event.", id, movingVertex.id, movingVertexIndex); // Use logger
					LOGGER.debug("          Reflex vertex V" + movingVertex.id + " moves over spoke edge opposite index " + movingVertexIndex
							+ ". Valid flip/spoke event.");
					return genericCollapse;
				} else {
					// log.debug(" Convex vertex V{} would move over spoke edge opposite index {}.
					// Invalidating event -> NEVER.", id, movingVertex.id, movingVertexIndex); //
					// Use logger
					LOGGER.debug("          Convex vertex V" + movingVertex.id + " would move over spoke edge opposite index " + movingVertexIndex
							+ ". Invalidating event -> NEVER.");
					return eventThatWillNotHappen(currentTime, determinant);
				}

			default :
				LOGGER.error("Unexpected type " + genericCollapse.getType() + " from getGenericCollapse in calculateFlipEvent for " + getName());
				System.exit(1); // abort!
				return null;
		}
	}

	private CollapseSpec getGenericCollapse(double currentTime, Polynomial determinant) {
		// log.debug(" KT{} Calculating Generic Collapse (for Flip/Spoke)", id); // Use
		// logger
		LOGGER.debug("          KT" + id + " Calculating Generic Collapse (for Flip/Spoke)");
		Optional<Double> collapseTimeOpt = getGenericCollapseTime(currentTime, determinant);

		if (!collapseTimeOpt.isPresent()) {
			// log.debug(" No generic collapse time found -> NEVER.", id); // Use logger
			LOGGER.debug("            No generic collapse time found -> NEVER.");
			return CollapseSpec.NEVER;
		}

		final double collapseTime = collapseTimeOpt.get();
		// log.debug(" Generic collapse time: {}", id, collapseTime); // Use logger
		LOGGER.debug("            Generic collapse time: {}", collapseTime);

		// Calculate squared lengths of edges at collapse time
		final Coordinate p0 = vertices[0].getPositionAt(collapseTime);
		final Coordinate p1 = vertices[1].getPositionAt(collapseTime);
		final Coordinate p2 = vertices[2].getPositionAt(collapseTime);

		final double sqLen0 = p1.distanceSq(p2); // Opposite vertex 0
		final double sqLen1 = p2.distanceSq(p0); // Opposite vertex 1
		final double sqLen2 = p0.distanceSq(p1); // Opposite vertex 2

		final double[] sqLengths = { sqLen0, sqLen1, sqLen2 };
		// log.debug(" SqLengths at collapse: {}", id, Arrays.toString(sqLengths)); //
		// Use logger
		LOGGER.debug("            SqLengths at collapse: [" + sqLen0 + ", " + sqLen1 + ", " + sqLen2 + "]");

		int zeroCount = 0;
		int firstZeroIndex = -1; // Index of the *edge* opposite the vertex index
		boolean[] isZero = new boolean[3];
		for (int i = 0; i < 3; i++) {
			if (sqLengths[i] < SurfConstants.ZERO_AREA_SQ) {
				isZero[i] = true;
				zeroCount++;
				if (firstZeroIndex == -1) {
					firstZeroIndex = i;
				}
			}
		}

		if (zeroCount == 2) {
			int thirdIndex = -1; // Find the index not in {firstZeroIndex, (firstZeroIndex+1)%3 or
									// (firstZeroIndex+2)%3}
			for (int i = 0; i < 3; ++i) {
				if (!isZero[i]) {
					thirdIndex = i;
				}
			}

			if (thirdIndex != -1 && sqLengths[thirdIndex] < SurfConstants.ZERO_AREA_SQ) {
				isZero[thirdIndex] = true;
				zeroCount = 3;
				// log.debug(" Zero count corrected from 2 to 3.", id); // Use logger
				LOGGER.debug("            Zero count corrected from 2 to 3.");
			} else {
				LOGGER.warn("Exactly 2 edges have zero length at generic collapse time " + collapseTime + " for KT" + id + ". Treating as full collapse.");
				zeroCount = 3;
			}
		}

		if (zeroCount == 3) {
			// log.debug(" All edges zero -> TRIANGLE_COLLAPSE", id); // Use logger
			LOGGER.debug("            All edges zero -> TRIANGLE_COLLAPSE");
			return new CollapseSpec(CollapseType.TRIANGLE_COLLAPSE, collapseTime, this);
		} else if (zeroCount == 1) {
			// log.debug(" One edge ({}) zero -> SPOKE_COLLAPSE", id, firstZeroIndex); //
			// Use logger
			LOGGER.debug("            One edge (" + firstZeroIndex + ") zero -> SPOKE_COLLAPSE");
			// Optional: Check other two lengths are equal within tolerance
			// double lenA = sqLengths[(firstZeroIndex + 1) % 3];
			// double lenB = sqLengths[(firstZeroIndex + 2) % 3];
			// if (Math.abs(lenA - lenB) > SurfConstants.ZERO_AREA_SQ * 10) { // Use larger
			// tolerance?
			// LOGGER.warn("Spoke collapse non-zero edges differ
			// significantly for KT" + id);
			// }
			return new CollapseSpec(CollapseType.SPOKE_COLLAPSE, collapseTime, this, firstZeroIndex);
		} else { // zeroCount == 0
			// log.debug(" No edges zero. Vertices collinear -> VERTEX_MOVES_OVER_SPOKE",
			// id); // Use logger
			LOGGER.debug("            No edges zero. Vertices collinear -> VERTEX_MOVES_OVER_SPOKE");

			ImmutableTriple<Integer, Integer, Integer> sortedIndices = indirectSort3(sqLengths);
			int longestEdgeIndex = sortedIndices.getRight(); // Index of the longest edge
			double longestSqLength = sqLengths[longestEdgeIndex];
			// log.debug(" Longest edge index: {} (SqLen: {})", id, longestEdgeIndex,
			// longestSqLength); // Use logger
			LOGGER.debug("            Longest edge index: " + longestEdgeIndex + " (SqLen: " + longestSqLength + ")");

			int relevantVertexIndex = longestEdgeIndex; // Vertex opposite longest edge moves over spoke
			return new CollapseSpec(CollapseType.VERTEX_MOVES_OVER_SPOKE, collapseTime, this, relevantVertexIndex, longestSqLength);
		}
	}

	static VertexOnSupportingLineResult getTimeVertexOnSupportingLine(WavefrontVertex v, WavefrontSupportingLine line) {
		if (v == null || line == null) {
			throw new NullPointerException("Null vertex or line in getTimeVertexOnSupportingLine");
		}

		final Vector2D n = line.getNormalDirection();
		final double nSqLen = n.lengthSquared();
		if (nSqLen < SurfConstants.ZERO_NORM_SQ) {
			LOGGER.warn("Supporting line has near-zero normal vector: " + line + " in vertex-line time calc.");
			return new VertexOnSupportingLineResult(Double.NaN, VertexOnSupportingLineType.NEVER);
		}
		final double nLen = Math.sqrt(nSqLen);

		final Coordinate Q0 = v.getInitialPosition(); // Position at t=0
		final Coordinate P0 = line.getSegment().getCoordinate(0); // Reference point on line at t=0
		final Vector2D s = Vector2D.create(v.getVelocity());
		final double w = line.getWeight();

		final Vector2D PQ = Vector2D.create(Q0).subtract(Vector2D.create(P0));

		// NOTE dot, since CGAL::Vector_2<Kernel> overloads the multiplication operator
		// *
		final double scaledDistance = n.dot(PQ);
		final double scaledVertexSpeedNormal = n.dot(s);
		final double scaledEdgeSpeedNormal = w * nLen;

		final double scaledSpeedApproach = scaledEdgeSpeedNormal - scaledVertexSpeedNormal;

		// log.debug(" TimeVertexOnLine: V{} vs Line ({})", v.id, line.getSegment()); //
		// Use logger
		// log.debug(" n={}, P0={}, Q0={}, s={}, w={}", n, P0, Q0, s, w); // Use logger
		// log.debug(" PQ={}, ScaledDist (num)={}, ScaledSpeedApproach (den)={}", PQ,
		// scaledDistance, scaledSpeedApproach); // Use logger
		LOGGER.debug("          TimeVertexOnLine: V" + v.id + " vs Line (" + line.getSegment().p0 + "->" + line.getSegment().p1 + ")");
		LOGGER.debug("            n=" + n + ", P0=" + P0 + ", Q0=" + Q0 + ", s=" + s + ", w=" + w);
		LOGGER.debug("            PQ=" + PQ);
		LOGGER.debug("            ScaledDist (num) = " + scaledDistance);
		LOGGER.debug("            ScaledSpeedApproach (den) = " + scaledSpeedApproach);

		double collapseTime;
		VertexOnSupportingLineType type;

		// Use a relative tolerance for speed approach? Or absolute? Using absolute for
		// now.
		if (Math.abs(scaledSpeedApproach) < SurfConstants.ZERO_SPEED) {
			// Check distance using appropriate tolerance (related to coordinates, not
			// speed)
			if (Math.abs(scaledDistance) < SurfConstants.ZERO_DIST * nLen) { // Scale tolerance by normal length?
				// log.debug(" Speeds equal, distance zero -> ALWAYS"); // Use logger
				LOGGER.debug("            Speeds equal, distance zero -> ALWAYS");
				type = VertexOnSupportingLineType.ALWAYS;
				collapseTime = 0.0;
			} else {
				// log.debug(" Speeds equal, distance non-zero -> NEVER"); // Use logger
				LOGGER.debug("            Speeds equal, distance non-zero -> NEVER");
				type = VertexOnSupportingLineType.NEVER;
				collapseTime = Double.NaN;
			}
		} else {
			collapseTime = scaledDistance / scaledSpeedApproach;
			// log.debug(" Speeds differ -> ONCE @ {}", collapseTime); // Use logger
			LOGGER.debug("            Speeds differ -> ONCE @ " + collapseTime);
			type = VertexOnSupportingLineType.ONCE;
		}
		return new VertexOnSupportingLineResult(collapseTime, type);
	}

	static Sign edgeIsFasterThanVertex(WavefrontVertex v, WavefrontSupportingLine line) {
		if (v == null || line == null) {
			throw new NullPointerException("Null vertex or line in edgeIsFasterThanVertex");
		}

		final Vector2D n = line.getNormalDirection();
		final double nSqLen = n.lengthSquared();
		if (nSqLen < SurfConstants.ZERO_NORM_SQ) {
			LOGGER.warn("Supporting line has near-zero normal vector in speed comparison: " + line);
			return Sign.ZERO;
		}
		final double nLen = Math.sqrt(nSqLen);

		final Vector2D s = Vector2D.create(v.getVelocity());
		final double w = line.getWeight();

		final double scaledEdgeSpeed = w * nLen;
		final double scaledVertexSpeed = n.dot(s);
		final double diff = scaledEdgeSpeed - scaledVertexSpeed;

		// log.debug(" EdgeFasterThanVertex: V{} vs Line ({})", v.id,
		// line.getSegment()); // Use logger
		// log.debug(" ScaledEdgeSpeed={}, ScaledVertexSpeed={}, Diff={}",
		// scaledEdgeSpeed, scaledVertexSpeed, diff); // Use logger
		LOGGER.debug("          EdgeFasterThanVertex: V" + v.id + " vs Line (" + line.getSegment().p0 + "->" + line.getSegment().p1 + ")");
		LOGGER.debug("            ScaledEdgeSpeed=" + scaledEdgeSpeed + ", ScaledVertexSpeed=" + scaledVertexSpeed + ", Diff=" + diff);

		// Use tolerance for comparison
		if (Math.abs(diff) < SurfConstants.ZERO_SPEED) {
			return Sign.ZERO;
		} else if (diff > 0) {
			return Sign.POSITIVE;
		} else {
			return Sign.NEGATIVE;
		}
	}

	// --- Triangulation Manipulation ---

	/**
	 * Moves a constraint (WavefrontEdge) from a dying source triangle's edge onto
	 * the edge of this triangle that bordered the source. Updates the edge's
	 * incident triangle and vertices. Does NOT modify this triangle's neighbor
	 * pointers (that happens later in the collapse handler).
	 *
	 * @param targetEdgeIndex Index in this triangle (0,1,2) corresponding to the
	 *                        edge that bordered the sourceTriangle.
	 * @param sourceTriangle  The dying triangle from which the constraint
	 *                        originates.
	 * @param sourceEdgeIndex Index in the sourceTriangle (0,1,2) where the
	 *                        constraint edge resides.
	 */
	public void moveConstraintFrom(int targetEdgeIndex, KineticTriangle sourceTriangle, int sourceEdgeIndex) {
		if (targetEdgeIndex < 0 || targetEdgeIndex >= 3) {
			throw new IllegalArgumentException("Invalid target index: " + targetEdgeIndex);
		}
		if (sourceEdgeIndex < 0 || sourceEdgeIndex >= 3) {
			throw new IllegalArgumentException("Invalid source index: " + sourceEdgeIndex);
		}
		if (sourceTriangle == null) {
			throw new NullPointerException("Source triangle is null");
		}
		// Allow source to be dying, that's the typical case
		// if (!sourceTriangle.isDying()) {
		// throw new IllegalStateException("Source triangle " + sourceTriangle.getName()
		// + " is not dying");
		// }
		// Check consistency: this triangle should have sourceTriangle as neighbor at
		// targetEdgeIndex *before* this call
		if (neighbors[targetEdgeIndex] != sourceTriangle) {
			LOGGER.warn("Consistency warning in moveConstraintFrom: Target edge " + targetEdgeIndex + " of " + getName() + " does not border source "
					+ sourceTriangle.getName() + " (Neighbor is: " + neighbors[targetEdgeIndex] + ")");
			// Don't throw exception here, as topology might be in flux during multi-event
			// steps.
		}
		if (isConstrained(targetEdgeIndex)) {
			throw new IllegalStateException("Target edge " + targetEdgeIndex + " of " + getName() + " is already constrained.");
		}
		if (!sourceTriangle.isConstrained(sourceEdgeIndex)) {
			throw new IllegalStateException("Source triangle " + sourceTriangle.getName() + " edge " + sourceEdgeIndex + " is not constrained.");
		}

		final WavefrontEdge edgeToMove = sourceTriangle.getWavefront(sourceEdgeIndex);
		if (edgeToMove == null) {
			// This might happen if the source triangle was already processed differently?
			LOGGER.error("Constraint edge to move is null from source T{} edge {}", sourceTriangle.getId(), sourceEdgeIndex);
			// Cannot proceed without the edge
			return;
			// Or throw new NullPointerException("Constraint edge to move is null from
			// source " + sourceTriangle.getName());
		}
		// Optional: Check if edge's incident triangle is indeed the source
		if (edgeToMove.getIncidentTriangle() != sourceTriangle) {
			LOGGER.warn("Edge {} incident triangle mismatch (expected T{}, found T{}) during moveConstraintFrom T{} edge {} -> T{} edge {}", edgeToMove.id,
					sourceTriangle.getId(), (edgeToMove.getIncidentTriangle() != null ? edgeToMove.getIncidentTriangle().getId() : "null"),
					sourceTriangle.getId(), sourceEdgeIndex, getId(), targetEdgeIndex);
			// Proceed cautiously, but this indicates a potential state issue.
		}

		LOGGER.debug("Moving constraint Edge " + edgeToMove.id + " from dying " + sourceTriangle.getName() + "(edge " + sourceEdgeIndex + ") to " + getName()
				+ "(edge " + targetEdgeIndex + ")");

		// 1. Update this triangle's links: Set wavefront, clear neighbor
		// ***** FIX: ONLY SET THE WAVEFRONT, DO NOT CLEAR NEIGHBOR YET *****
		// this.neighbors[targetEdgeIndex] = null; // <<< REMOVED THIS LINE
		this.wavefronts[targetEdgeIndex] = edgeToMove;

		// 2. Update the edge's incident triangle
		edgeToMove.setIncidentTriangle(this);

		// 3. Update the edge's vertices to match this triangle's perspective
		final WavefrontVertex vCCW = this.getVertex(ccw(targetEdgeIndex));
		final WavefrontVertex vCW = this.getVertex(cw(targetEdgeIndex));
		if (vCCW == null || vCW == null) {
			throw new IllegalStateException("Null vertex found during moveConstraintFrom for " + getName());
		}
		edgeToMove.setVerticesAndUpdateAdj(vCCW, vCW);

		// 4. Null out the constraint on the source triangle
		sourceTriangle.wavefronts[sourceEdgeIndex] = null;

		// 5. Invalidate collapse specs
		this.invalidateCollapseSpec();
		// Source triangle spec doesn't matter (it's dying)
		// Edge spec was invalidated by setVerticesAndUpdateAdj
	}

	public void doRawFlip(int edgeIndex) {
		// --- Input Validation ---
		if (edgeIndex < 0 || edgeIndex >= 3) {
			throw new IllegalArgumentException("Invalid flip index: " + edgeIndex);
		}
		if (isConstrained(edgeIndex)) {
			throw new IllegalStateException("Cannot flip constrained edge " + edgeIndex + " of " + getName());
		}
		final KineticTriangle n = neighbors[edgeIndex]; // The neighbor triangle
		if (n == null) {
			throw new IllegalStateException("Cannot flip edge " + edgeIndex + " of " + getName() + ": no neighbor");
		}
		if (n.isDead() || n.isDying()) {
			throw new IllegalStateException("Cannot flip edge " + edgeIndex + " of " + getName() + ": neighbor " + n.getName() + " is dead/dying");
		}
		final int nEdgeIndex = n.indexOfNeighbor(this); // Edge index in neighbor
		if (nEdgeIndex < 0) {
			// Should not happen if topology is consistent before flip
			throw new IllegalStateException("Neighbor " + n.getName() + " does not point back to " + getName());
		}
		if (n.isConstrained(nEdgeIndex)) {
			throw new IllegalStateException("Cannot flip edge " + edgeIndex + ": neighbor " + n.getName() + " edge " + nEdgeIndex + " is constrained");
		}

		LOGGER.debug("Flipping edge between {} (idx {}) and {} (idx {})", getName(), edgeIndex, n.getName(), nEdgeIndex);

		// --- Identify Vertices ---
		final WavefrontVertex v_opposite_this = this.getVertex(edgeIndex); // Vertex 'v' in C++
		final WavefrontVertex v_opposite_n = n.getVertex(nEdgeIndex); // Vertex 'o' in C++
		// Shared vertices (relative to 'this' triangle's edgeIndex)
		final WavefrontVertex v_shared_ccw = this.getVertex(ccw(edgeIndex)); // Vertex 'v1' in C++
		final WavefrontVertex v_shared_cw = this.getVertex(cw(edgeIndex)); // Vertex 'v2' in C++

		// --- Store Original Neighbor/Wavefront Info (needed for backlink updates) ---
		KineticTriangle n_neighbor_across_n_cw = n.getNeighbor(cw(nEdgeIndex));
		WavefrontEdge w_wavefront_across_n_cw = n.getWavefront(cw(nEdgeIndex));
		KineticTriangle n_neighbor_across_this_cw = this.getNeighbor(cw(edgeIndex));
		WavefrontEdge w_wavefront_across_this_cw = this.getWavefront(cw(edgeIndex));

		// --- Step 1: Update Vertices ---
		// Update using setVertex to handle potential edge endpoint updates if needed
		this.setVertex(ccw(edgeIndex), v_opposite_n); // this.vertex[v1_idx] = o
		n.setVertex(ccw(nEdgeIndex), v_opposite_this); // n.vertex[v_shared_cw_idx_in_n] = v

		// --- Step 2: Update Neighbor/Wavefront pointers (Following C++) ---
		// Edge 'edgeIndex' in 'this' gets properties from edge 'cw(nEdgeIndex)' in 'n'
		this.neighbors[edgeIndex] = n_neighbor_across_n_cw;
		this.wavefronts[edgeIndex] = w_wavefront_across_n_cw;

		// Edge 'nEdgeIndex' in 'n' gets properties from edge 'cw(edgeIndex)' in 'this'
		n.neighbors[nEdgeIndex] = n_neighbor_across_this_cw;
		n.wavefronts[nEdgeIndex] = w_wavefront_across_this_cw;

		// Edge 'cw(nEdgeIndex)' in 'n' now points to 'this' (new shared edge)
		n.neighbors[cw(nEdgeIndex)] = this;
		n.wavefronts[cw(nEdgeIndex)] = null;

		// Edge 'cw(edgeIndex)' in 'this' now points to 'n' (new shared edge)
		this.neighbors[cw(edgeIndex)] = n;
		this.wavefronts[cw(edgeIndex)] = null;

		// --- Step 3: Update Backlinks ---
		// Update neighbor that is now across edge 'this[edgeIndex]'
		if (this.wavefronts[edgeIndex] != null) { // If it's now constrained
			this.wavefronts[edgeIndex].setIncidentTriangle(this);
		} else if (this.neighbors[edgeIndex] != null) { // If it's an internal edge
			KineticTriangle newNeighbor = this.neighbors[edgeIndex];
			// Find where newNeighbor pointed to 'n' and update it to point to 'this'
			int idxInNewNeighbor = newNeighbor.indexOfNeighbor(n);
			if (idxInNewNeighbor != -1) {
				newNeighbor.setNeighborRaw(idxInNewNeighbor, this);
			} else {
				// This could happen if topology was inconsistent *before* the flip
				LOGGER.warn("Could not find backlink from new neighbor {} to old neighbor {} during flip of T{}/T{}", newNeighbor.getName(), n.getName(),
						getName(), n.getName());
			}
		}

		// Update neighbor that is now across edge 'n[nEdgeIndex]'
		if (n.wavefronts[nEdgeIndex] != null) { // If it's now constrained
			n.wavefronts[nEdgeIndex].setIncidentTriangle(n);
		} else if (n.neighbors[nEdgeIndex] != null) { // If it's an internal edge
			KineticTriangle newNeighborN = n.neighbors[nEdgeIndex];
			// Find where newNeighborN pointed to 'this' and update it to point to 'n'
			int idxInNewNeighborN = newNeighborN.indexOfNeighbor(this);
			if (idxInNewNeighborN != -1) {
				newNeighborN.setNeighborRaw(idxInNewNeighborN, n);
			} else {
				LOGGER.warn("Could not find backlink from new neighbor {} to old neighbor {} during flip of T{}/T{}", newNeighborN.getName(), getName(),
						getName(), n.getName());
			}
		}

		// --- Step 4: Invalidate Collapse Specs ---
		// Invalidation should be handled by setVertex calls and neighbor updates.
		// Explicitly invalidate just in case.
		this.invalidateCollapseSpec();
		n.invalidateCollapseSpec();
		// Also invalidate neighbors whose pointers changed
		if (n_neighbor_across_n_cw != null)
			n_neighbor_across_n_cw.invalidateCollapseSpec();
		if (n_neighbor_across_this_cw != null)
			n_neighbor_across_this_cw.invalidateCollapseSpec();

		LOGGER.trace("Flip complete.");
	}

	// --- Validation ---

	/**
	 * Checks internal consistency. Comment out for production builds.
	 *
	 * @throws AssertionError if inconsistent.
	 */
	/**
	 * Checks that this triangle’s local pointers (vertices, neighbors,
	 * wavefront‐edges and their back‐links) form a consistent fan. Throws
	 * IllegalStateException on any invariant violation.
	 */
	public void assertValid() {
		if (isDead()) {
			throw new IllegalStateException("assertValid(): triangle is dead");
		}

		// For each corner i = 0,1,2
		for (int i = 0; i < 3; i++) {
			WavefrontVertex vi = vertices[i];
			KineticTriangle nbr = neighbors[i];
			WavefrontEdge wf = wavefronts[i];

			// must have exactly one of (neighbor) XOR (wavefront)
			if ((nbr != null) == (wf != null)) {
				throw new IllegalStateException(
						String.format("Vertex/neigh vs. wavefront mismatch at corner %d of T%d: hasNbr=%b, hasWf=%b", i, getId(), nbr != null, wf != null));
			}

			// vertex must exist
			if (vi == null) {
				throw new IllegalStateException(String.format("Missing vertex %d in T%d", i, getId()));
			}

			if (nbr != null) {
				// interior edge: neighbor must point back
				if (!nbr.hasNeighbor(this)) {
					throw new IllegalStateException(String.format("Neighborhood inconsistency: T%d not found in neighbor[%d]=T%d", getId(), i, nbr.getId()));
				}
				int idxInNbr = nbr.indexOfNeighbor(this);

				// the two “other” vertices must match across the shared edge
				int cwI = cw(i);
				int ccwI = ccw(i);
				int nbrCcw = ccw(idxInNbr);
				int nbrCw = cw(idxInNbr);

				WavefrontVertex myLeft = vertices[cwI];
				WavefrontVertex nbrLeft = nbr.getVertex(nbrCcw);
				if (myLeft != nbrLeft) {
					throw new IllegalStateException(String.format("Edge‐vertex mismatch (left) at T%d[%d] vs nbr T%d[%d]: %s vs %s", getId(), cwI, nbr.getId(),
							nbrCcw, myLeft, nbrLeft));
				}

				WavefrontVertex myRight = vertices[ccwI];
				WavefrontVertex nbrRight = nbr.getVertex(nbrCw);
				if (myRight != nbrRight) {
					throw new IllegalStateException(String.format("Edge‐vertex mismatch (right) at T%d[%d] vs nbr T%d[%d]: %s vs %s", getId(), ccwI,
							nbr.getId(), nbrCw, myRight, nbrRight));
				}

			} else {
				// boundary edge: must have a wavefront
				if (wf == null) {
					throw new IllegalStateException(String.format("Missing wavefront on boundary edge %d of T%d", i, getId()));
				}
				// the wavefront must reference this triangle
				if (wf.getIncidentTriangle() != this) {
					throw new IllegalStateException(
							String.format("Wavefront‐triangle mismatch: wf.incident=%s vs this=T%d", wf.getIncidentTriangle(), getId()));
				}
				if (wf.isDead()) {
					throw new IllegalStateException(String.format("Wavefront on boundary edge %d of T%d is dead", i, getId()));
				}

				// the wavefront’s endpoints must be the two opposite vertices:
				// wf.vertex(0) == vertices[ccw(i)]
				// wf.vertex(1) == vertices[cw(i)]
				int cwI = cw(i);
				int ccwI = ccw(i);
				WavefrontVertex wf0 = wf.getVertex(0);
				WavefrontVertex wf1 = wf.getVertex(1);

				if (wf0 != vertices[ccwI] || wf1 != vertices[cwI]) {
					throw new IllegalStateException(String.format("Wavefront endpoint mismatch in T%d[%d]: got (%s,%s) expected (%s,%s)", getId(), i, wf0, wf1,
							vertices[ccwI], vertices[cwI]));
				}

				// each vertex must have this as its incident wavefront on the correct side
				// (side 1 for ccw‐vertex, side 0 for cw‐vertex)
				WavefrontEdge v0Side1 = vertices[ccwI].getIncidentEdge(1);
				WavefrontEdge v1Side0 = vertices[cwI].getIncidentEdge(0);
				if (v0Side1 != wf || v1Side0 != wf) {
					throw new IllegalStateException(String.format("Vertex<‐>wavefront back‐links mismatch in T%d[%d]: v[%d].e1=%s, v[%d].e0=%s, wf=%s", getId(),
							i, ccwI, v0Side1, cwI, v1Side0, wf));
				}
			}
		}
	}

	// --- toString, equals, hashCode ---

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getName());
		sb.append("[");
		sb.append(vertices[0] != null ? "V" + vertices[0].id : "null").append(",");
		sb.append(vertices[1] != null ? "V" + vertices[1].id : "null").append(",");
		sb.append(vertices[2] != null ? "V" + vertices[2].id : "null");
		sb.append("] N[");
		sb.append(neighbors[0] != null ? neighbors[0].id : "null").append(",");
		sb.append(neighbors[1] != null ? neighbors[1].id : "null").append(",");
		sb.append(neighbors[2] != null ? neighbors[2].id : "null");
		sb.append("] W[");
		sb.append(wavefronts[0] != null ? "E" + wavefronts[0].id : "null").append(",");
		sb.append(wavefronts[1] != null ? "E" + wavefronts[1].id : "null").append(",");
		sb.append(wavefronts[2] != null ? "E" + wavefronts[2].id : "null");
		sb.append("]");
		if (isDying()) {
			sb.append(" (dying)");
		}
		if (isDead()) {
			sb.append(" (DEAD)");
		}
		// sb.append(" C").append(component);
		return sb.toString();
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
		KineticTriangle other = (KineticTriangle) obj;
		return id == other.id; // Identity based on ID
	}

	// --- Helper Classes and Enums ---

	enum VertexOnSupportingLineType {
		ONCE, ALWAYS, NEVER
	}

	static class VertexOnSupportingLineResult {
		public final double time;
		public final VertexOnSupportingLineType type;

		public VertexOnSupportingLineResult(double time, VertexOnSupportingLineType type) {
			this.time = time;
			this.type = type;
		}

		@Override
		public String toString() {
			return "Time=" + time + ", Type=" + type;
		}
	}

	public enum Sign {
		NEGATIVE, ZERO, POSITIVE
	}

	private static ImmutableTriple<Integer, Integer, Integer> indirectSort3(double[] values) {
		if (values.length != 3) {
			throw new IllegalArgumentException("Requires array of length 3");
		}
		int i0 = 0, i1 = 1, i2 = 2;
		// Simple bubble sort for 3 elements
		if (values[i0] > values[i1]) {
			i0 = 1;
			i1 = 0;
		} // Swap indices if needed
		if (values[i1] > values[i2]) {
			int tmp = i1;
			i1 = i2;
			i2 = tmp;
		}
		if (values[i0] > values[i1]) {
			int tmp = i0;
			i0 = i1;
			i1 = tmp;
		}
		return ImmutableTriple.of(i0, i1, i2);
	}

	private CollapseSpec eventThatWillNotHappen(double currentTime, Polynomial determinant) {
		// C++ Debug build: return INVALID_EVENT at determinant's zero time if it exists
		// future
		// C++ Release build: return NEVER
		/*
		 * // Mimic Debug build: Optional<Double> timeOpt =
		 * getGenericCollapseTime(currentTime, determinant); if (timeOpt.isPresent()) {
		 * return new CollapseSpec(CollapseType.INVALID_EVENT, timeOpt.get(), this); }
		 * else { return CollapseSpec.NEVER; }
		 */
		// Mimic Release build:
		return CollapseSpec.NEVER;
	}

}