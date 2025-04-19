package com.github.micycle1.surferj.collapse;

import java.util.Objects;

import com.github.micycle1.surferj.SurfConstants;
import com.github.micycle1.surferj.kinetics.KineticTriangle;

/**
 * Represents the specification for a potential future event (usually a
 * collapse) for a kinetic triangle. Includes the event type, time, and
 * potentially relevant edge or secondary key for tie-breaking. Implements
 * Comparable for use in priority queues.
 */
public class CollapseSpec implements Comparable<CollapseSpec> {

	protected CollapseType type;
	protected final double time; // Event time, NaN if NEVER
	// NOTE: Kept triangle field in Java CollapseSpec, unlike C++ where it's added
	// in Event
	protected final KineticTriangle triangle;
	protected final int relevantEdge; // Index (0,1,2) in 'triangle', -1 if not applicable
	protected final double secondaryKey; // For tie-breaking, NaN if not applicable
	protected final int component; // Component ID, used for priority comparison

	protected boolean isValid = true; // Flag for EventQueue strategy (Java specific)

	// --- Static NEVER Instance ---
	// Assign a default component (e.g., -1)
	public static final CollapseSpec NEVER = new CollapseSpec(CollapseType.NEVER, Double.NaN, null, -1, Double.NaN, -1);

	// --- Base Constructor (Takes all fields including component) ---
	// This is the constructor all others should ideally delegate to.
	public CollapseSpec(CollapseType type, double time, KineticTriangle triangle, int relevantEdge, double secondaryKey, int component) {
		this.type = type;
		this.triangle = triangle; // Can be null for NEVER type
		this.relevantEdge = relevantEdge;
		this.secondaryKey = secondaryKey;
		this.component = component;

		// Time validation and clamping logic (same as before)
		if (type != CollapseType.NEVER && type != CollapseType.UNDEFINED) {
			double clampedTime = (time < 0.0 && time >= -SurfConstants.TIME_TOL) ? 0.0 : time; // Clamp near-zero negatives
			if (Double.isNaN(clampedTime) || clampedTime < -SurfConstants.TIME_TOL) {
				System.err.println("Warning: CollapseSpec time " + time + " for type " + type + " is invalid. Setting to NEVER.");
				this.time = Double.NaN;
				this.type = CollapseType.NEVER; // Force NEVER if time is truly invalid
			} else {
				this.time = clampedTime;
			}
		} else {
			this.time = Double.NaN; // Ensure NEVER/UNDEFINED always have NaN time
		}

		// Basic validation based on type requirements (logging warnings)
		if (requiresRelevantEdge(type) && relevantEdge < 0) {
			System.err.println("Warning: CollapseType " + type + " requires relevantEdge >= 0, but got " + relevantEdge);
		}
		if (requiresSecondaryKey(type) && Double.isNaN(secondaryKey)) {
			System.err.println("Warning: CollapseType " + type + " requires secondaryKey, but it's NaN");
		}
	}

	// --- NEW Compatibility Constructors (Extract component from triangle) ---

	/**
	 * Compatibility constructor. Infers component ID from the triangle.
	 *
	 * @param type     Event type.
	 * @param time     Event time.
	 * @param triangle Associated triangle (must provide component ID via
	 *                 getComponent()). Can be null only if type is NEVER/UNDEFINED.
	 */
	public CollapseSpec(CollapseType type, double time, KineticTriangle triangle) {
		this(type, time, triangle, -1, Double.NaN, // Provide defaults for edge/key
				(triangle != null) ? triangle.getComponent() : -1); // Extract component
		// Add validation check if needed (e.g., ensuring triangle is not null for types
		// requiring it)
		if (triangle == null && !(type == CollapseType.NEVER || type == CollapseType.UNDEFINED)) {
			System.err.println("Warning: Non-NEVER/UNDEFINED CollapseSpec created with null triangle.");
			// Consider throwing an exception based on type requirements
		}
	}

	/**
	 * Compatibility constructor. Infers component ID from the triangle.
	 *
	 * @param type         Event type.
	 * @param time         Event time.
	 * @param triangle     Associated triangle (must provide component ID via
	 *                     getComponent()).
	 * @param relevantEdge Relevant edge index.
	 */
	public CollapseSpec(CollapseType type, double time, KineticTriangle triangle, int relevantEdge) {
		this(type, time, triangle, relevantEdge, Double.NaN, // Provide default for key
				(triangle != null) ? triangle.getComponent() : -1); // Extract component
		if (triangle == null) { // Add null check
			System.err.println("Warning: CollapseSpec created with null triangle.");
		}
	}

	/**
	 * Compatibility constructor. Infers component ID from the triangle.
	 *
	 * @param type         Event type.
	 * @param time         Event time.
	 * @param triangle     Associated triangle (must provide component ID via
	 *                     getComponent()).
	 * @param relevantEdge Relevant edge index.
	 * @param secondaryKey Secondary key for tie-breaking.
	 */
	public CollapseSpec(CollapseType type, double time, KineticTriangle triangle, int relevantEdge, double secondaryKey) {
		this(type, time, triangle, relevantEdge, secondaryKey, // Pass through values
				(triangle != null) ? triangle.getComponent() : -1); // Extract component
		if (triangle == null) { // Add null check
			System.err.println("Warning: CollapseSpec created with null triangle.");
		}
	}

	public CollapseSpec(CollapseType type, double time, KineticTriangle triangle, int relevantEdge, int component) {
		this(type, time, triangle, relevantEdge, Double.NaN, component);
	}

	// --- Constructor from EdgeCollapseSpec (No longer needs explicit component)
	// ---
	/**
	 * Constructor mapping from an EdgeCollapseSpec. Infers component from triangle.
	 * Mimics C++ constructor CollapseSpec(int, const EdgeCollapseSpec&, int).
	 *
	 * @param edgeCollapseSpec The edge collapse info.
	 * @param triangle         The triangle this spec pertains to (used to get
	 *                         component).
	 * @param collapsingEdge   The index (0,1,2) of the edge within the triangle.
	 */
	public CollapseSpec(EdgeCollapseSpec edgeCollapseSpec, KineticTriangle triangle, int collapsingEdge) {
		this( // Delegate to base constructor
				mapEdgeCollapseType(edgeCollapseSpec.getType()), // Map type
				(edgeCollapseSpec.getType() == EdgeCollapseType.FUTURE || edgeCollapseSpec.getType() == EdgeCollapseType.ALWAYS) ? edgeCollapseSpec.getTime()
						: Double.NaN, // Set time or NaN
				triangle, (mapEdgeCollapseType(edgeCollapseSpec.getType()) == CollapseType.CONSTRAINT_COLLAPSE) ? collapsingEdge : -1, // Set edge only if
																																		// relevant type
				Double.NaN, // No secondary key from EdgeCollapseSpec
				(triangle != null) ? triangle.getComponent() : -1 // Extract component
		);
		// C++ assertion check: edge_collapse.type() must be FUTURE, ALWAYS, NEVER, or
		// PAST. Handled by mapEdgeCollapseType.
	}

	// Inside CollapseSpec class...

	/**
	 * Internal constructor used for delegation, typically from subclasses like
	 * Event. It copies properties from a source CollapseSpec but explicitly sets
	 * the triangle reference for this new instance.
	 *
	 * @param sourceSpec    The pre-calculated CollapseSpec whose properties to
	 *                      copy.
	 * @param eventTriangle The specific KineticTriangle this new CollapseSpec
	 *                      instance pertains to.
	 */
	protected CollapseSpec(CollapseSpec sourceSpec, KineticTriangle eventTriangle) {
		// Initialize fields by copying from the source spec
		this( // Delegate to the main constructor
				sourceSpec.getType(), sourceSpec.getTime(), eventTriangle, // Use the explicitly provided triangle for this instance
				sourceSpec.getRelevantEdge(), sourceSpec.getSecondaryKey(), sourceSpec.getComponent());
		// Ensure the source spec wasn't null
		Objects.requireNonNull(sourceSpec, "Source CollapseSpec cannot be null in delegating constructor");
		// Optional: Add an assertion to check if eventTriangle matches
		// sourceSpec.getTriangle()
		// assert eventTriangle == sourceSpec.getTriangle() : "Triangle mismatch during
		// delegation";
	}

	// Helper to map EdgeCollapseType (Unchanged)
	private static CollapseType mapEdgeCollapseType(EdgeCollapseType edgeType) {
		switch (edgeType) {
			case FUTURE :
			case ALWAYS :
				return CollapseType.CONSTRAINT_COLLAPSE;
			case PAST :
			case NEVER :
				return CollapseType.NEVER;
			case UNDEFINED :
			default :
				return CollapseType.UNDEFINED;
		}
	}

	// --- Factory Method (Now extracts component internally) ---
	// This is now consistent with the constructor `CollapseSpec(EdgeCollapseSpec,
	// KineticTriangle, int)`
	public static CollapseSpec fromEdgeCollapse(EdgeCollapseSpec edgeSpec, KineticTriangle triangle, int relevantEdgeInTriangle) {
		// Extracts component internally via the constructor it calls
		return new CollapseSpec(edgeSpec, triangle, relevantEdgeInTriangle);
	}

	// --- Getters ---
	public CollapseType getType() {
		return type;
	}

	public double getTime() {
		return time;
	}

	public KineticTriangle getTriangle() {
		return triangle;
	}

	public int getRelevantEdge() {
		return relevantEdge;
	}

	public double getSecondaryKey() {
		return secondaryKey;
	}

	public int getComponent() {
		return component;
	}

	// --- Validity Methods (Java Specific) ---
	public void invalidate() {
		this.isValid = false;
	}

	public boolean isValid() {
		return this.isValid;
	}

	public boolean allowsRefinementTo(CollapseSpec o) {
		if (type == CollapseType.SPLIT_OR_FLIP_REFINE) {
			if (o.type == CollapseType.VERTEX_MOVES_OVER_SPOKE || o.type == CollapseType.SPOKE_COLLAPSE) {
				return relevantEdge != o.relevantEdge;
			}
		}
		return false;
	}

	// --- Static Helper Methods (Requirement Checks) ---
	// Visibility changed to package-private or public if needed outside
	static boolean requiresRelevantEdge(CollapseType type) {
		// Matches C++ list exactly
		return type == CollapseType.CONSTRAINT_COLLAPSE || type == CollapseType.SPOKE_COLLAPSE || type == CollapseType.SPLIT_OR_FLIP_REFINE
				|| type == CollapseType.VERTEX_MOVES_OVER_SPOKE || type == CollapseType.CCW_VERTEX_LEAVES_CH
				|| type == CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED;
		// C++ includes || false which is redundant
	}

	static boolean requiresSecondaryKey(CollapseType type) {
		// Matches C++ list exactly
		return type == CollapseType.VERTEX_MOVES_OVER_SPOKE || type == CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED;
		// C++ includes || false which is redundant
	}

	// --- Comparison Logic (Crucial for Event Queue) ---
	@Override
	public int compareTo(CollapseSpec other) {
		// Handle NEVER cases first (always highest/latest)
		boolean thisIsNever = (this.type == CollapseType.NEVER);
		boolean otherIsNever = (other.type == CollapseType.NEVER);

		if (thisIsNever && otherIsNever) {
			return 0;
		}
		if (thisIsNever) {
			return 1; // this is later
		}
		if (otherIsNever) {
			return -1; // other is later
		}

		// Compare component ID (lower component ID has higher priority/comes earlier)
		int compCompare = Integer.compare(this.component, other.component);
		if (compCompare != 0) {
			return compCompare;
		}

		// Compare time with tolerance
		// Using Double.compare is generally safer than subtraction for
		// NaN/Infinity/-0.0
		int timeCompare = Double.compare(this.time, other.time);
		// Apply tolerance ONLY if Double.compare results in 0
		if (timeCompare == 0 || Math.abs(this.time - other.time) <= SurfConstants.TIME_TOL) {
			// Times are effectively equal, compare by type priority (lower enum ordinal =
			// higher priority)
			int typeCompare = Integer.compare(this.type.ordinal(), other.type.ordinal());
			if (typeCompare != 0) {
				return typeCompare;
			}

			// If times and types are equal, use secondary key for tie-breaking if
			// applicable
			// Higher secondary key means HIGHER priority (comes earlier) -> reverse
			// comparison
			if (requiresSecondaryKey(this.type)) {
				// Check for NaN before comparing secondary keys
				boolean thisHasKey = !Double.isNaN(this.secondaryKey);
				boolean otherHasKey = !Double.isNaN(other.secondaryKey);
				if (thisHasKey && otherHasKey) {
					// Reverse compare: higher key comes first
					return Double.compare(other.secondaryKey, this.secondaryKey);
				} else if (thisHasKey) {
					return -1; // This has a key, other doesn't -> this comes first
				} else if (otherHasKey) {
					return 1; // Other has a key, this doesn't -> other comes first
				}
				// If neither has a valid secondary key, they remain equal at this stage.
			}
			// If still equal, consider them equal for priority purposes
			return 0;

		} else {
			// Times are different beyond tolerance
			return timeCompare;
		}
	}

	// --- toString, equals, hashCode ---
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("CollapseSpec{");
		sb.append(type);
		if (triangle != null) {
			sb.append(", tri=").append(triangle.getId()); // Use ID for brevity
		}
		if (component >= 0) { // Only print component if non-negative? Or always? C++ always prints.
			sb.append(", comp=").append(component);
		}
		// Use isNaN check for time, as NEVER type might not guarantee NaN if
		// constructed improperly
		if (!Double.isNaN(time)) {
			sb.append(", time=").append(String.format("%.6f", time));
		}
		if (requiresRelevantEdge(type) && relevantEdge != -1) { // Only print edge if required and valid
			sb.append(", edge=").append(relevantEdge);
		}
		if (requiresSecondaryKey(type) && !Double.isNaN(secondaryKey)) { // Only print key if required and valid
			sb.append(", key=").append(String.format("%.4f", secondaryKey));
		}
		sb.append('}');
		return sb.toString();
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		// Use instanceof to allow comparison with subclasses like Event
		if (!(o instanceof CollapseSpec)) {
			return false;
		}
		CollapseSpec that = (CollapseSpec) o;

		// Equality based on the comparison logic (compareTo == 0)
		// AND potentially stricter checks if needed (e.g., exact triangle instance or
		// ID).
		// C++ operator== just uses compare(o) == CGAL::EQUAL.
		// Let's stick to the C++ definition: equality means they have the same
		// priority.
		// Adding triangle ID check makes it stricter than C++ priority equality.
		// Reverting to priority-based equality:
		return this.compareTo(that) == 0;
		// If stricter equality needed (e.g., for HashMap keys if used):
		// return this.compareTo(that) == 0 &&
		// this.component == that.component && // Redundant if compareTo checks
		// component
		// Objects.equals(this.triangle == null ? null : this.triangle.getId(),
		// that.triangle == null ? null : that.triangle.getId());
	}

	@Override
	public int hashCode() {
		// Hash based on fields determining priority in compareTo.
		// Must be consistent with the fields used when compareTo returns 0.
		// Primary fields: component, time (discretized), type, secondaryKey (if
		// applicable).
		// Include triangle ID for better distribution if used in HashMaps/Sets for
		// non-priority reasons.

		// Discretize time based on tolerance for hashing
		long timeBits = Double.doubleToLongBits(time);
		long timeHash = Double.isNaN(time) ? 0 : (long) (time / SurfConstants.TIME_TOL); // Discretize

		long secondaryKeyBits = Double.doubleToLongBits(secondaryKey);
		long secondaryKeyHash = (requiresSecondaryKey(type) && !Double.isNaN(secondaryKey)) ? (long) (secondaryKey / SurfConstants.TIME_TOL) // Discretize if
																																				// used
				: 0;

		// Combine hashes. Use fields consistent with compareTo returning 0.
		return Objects.hash(component, timeHash, type, secondaryKeyHash // Only relevant if type requires it
		// Optional: Add triangle ID if stricter equality/hashing is desired beyond
		// priority
		// , (triangle != null ? triangle.getId() : 0)
		);
	}
}