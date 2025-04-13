package com.github.micycle.surferj2.collapse;

import com.github.micycle.surferj2.kinetics.KineticTriangle;

import java.util.Objects;

public class CollapseSpec implements Comparable<CollapseSpec> {
	private final CollapseType type;
	private final double time; // Event time, NaN if NEVER
	private final KineticTriangle triangle; // The triangle this event pertains to
	private final int relevantEdge; // Index (0,1,2) in 'triangle', -1 if not applicable
	private final double secondaryKey; // For tie-breaking (e.g., edge length, speed), NaN if not applicable

	// Flag for the EventQueue strategy
	// Use AtomicBoolean if there's any potential for concurrent modification,
	// otherwise, a simple volatile boolean is fine. Sticking with simple boolean
	// for now.
	private boolean isValid = true; // Default to true

	public static final CollapseSpec NEVER = new CollapseSpec(CollapseType.NEVER, Double.NaN, null, -1, Double.NaN);

	// Constructor variations based on C++
	public CollapseSpec(CollapseType type, double time, KineticTriangle triangle, int relevantEdge, double secondaryKey) {
		this.type = type;
		this.triangle = triangle;
		this.relevantEdge = relevantEdge; // Assume validation happens elsewhere or is implied by type
		this.secondaryKey = secondaryKey;

		// Ensure time is non-negative for relevant types
		if (type != CollapseType.NEVER && type != CollapseType.UNDEFINED && (Double.isNaN(time) || time < -1e-9)) { // Allow small negatives near zero
			throw new IllegalArgumentException("Time must be valid and non-negative for most collapse types, got: " + time + " for type " + type);
		}
		// Clamp small negative times resulting from floating point errors
		this.time = (type != CollapseType.NEVER && type != CollapseType.UNDEFINED && time < 0.0) ? 0.0 : time;

		if (type == CollapseType.NEVER || type == CollapseType.UNDEFINED) {
			if (!Double.isNaN(this.time)) {
				// System.err.println("Warning: Time should be NaN for NEVER/UNDEFINED
				// CollapseSpec, but was " + this.time);
			}
		}
	}

	// Simplified constructors
	public CollapseSpec(CollapseType type, double time, KineticTriangle triangle) {
		this(type, time, triangle, -1, Double.NaN);
		if (requiresRelevantEdge(type)) {
			throw new IllegalArgumentException("CollapseType " + type + " requires a relevant edge index.");
		}
	}

	public CollapseSpec(CollapseType type, double time, KineticTriangle triangle, int relevantEdge) {
		this(type, time, triangle, relevantEdge, Double.NaN);
		if (!requiresRelevantEdge(type)) {
			// Allow edge for types that *might* use it later, but warn if strictly not
			// needed?
			// System.err.println("Warning: Relevant edge provided for CollapseType " + type
			// + " which might not use it.");
		}
		if (requiresSecondaryKey(type)) {
			throw new IllegalArgumentException("CollapseType " + type + " requires a secondary key.");
		}
	}

	// New constructor to match C++ pattern
	public CollapseSpec(int component, EdgeCollapseSpec edgeCollapseSpec, int collapsingEdge, KineticTriangle triangle) {
		this.component = component;
		this.triangle = triangle;
		this.relevantEdge = collapsingEdge;

		// Map EdgeCollapseType to CollapseType
		switch (edgeCollapseSpec.getType()) {
			case FUTURE :
			case ALWAYS :
				this.type = CollapseType.CONSTRAINT_COLLAPSE;
				this.time = edgeCollapseSpec.getTime();
				break;
			case PAST :
			case NEVER :
				this.type = CollapseType.NEVER;
				this.time = Double.NaN;
				break;
			case UNDEFINED :
				this.type = CollapseType.UNDEFINED;
				this.time = Double.NaN;
				break;
			default :
				throw new IllegalArgumentException("Unknown EdgeCollapseType: " + edgeCollapseSpec.getType());
		}
		this.secondaryKey = Double.NaN; // Default; adjust if needed
	}

	// Factory from EdgeCollapseSpec (simulates C++ constructor)
	public static CollapseSpec fromEdgeCollapse(EdgeCollapseSpec edgeSpec, KineticTriangle triangle, int relevantEdgeInTriangle) {
		if (edgeSpec.getType() == EdgeCollapseType.FUTURE || edgeSpec.getType() == EdgeCollapseType.ALWAYS) {
			return new CollapseSpec(CollapseType.CONSTRAINT_COLLAPSE, edgeSpec.getTime(), triangle, relevantEdgeInTriangle);
		} else {
			// Map NEVER, PAST, UNDEFINED EdgeCollapseTypes to NEVER CollapseSpec
			return new CollapseSpec(CollapseType.NEVER, Double.NaN, triangle, -1, Double.NaN);
		}
	}

	private int component; // Added to store the component from C++ NOTE should be final
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

	/** Marks this event instance as invalid (it will be skipped by the queue). */
	public void invalidate() {
		this.isValid = false;
	}

	/** Checks if this event instance is still considered valid. */
	public boolean isValid() {
		return this.isValid;
	}

	public boolean allowsRefinementTo(CollapseSpec o) {
		if (type == CollapseType.SPLIT_OR_FLIP_REFINE) {
			if (o.type == CollapseType.VERTEX_MOVES_OVER_SPOKE || o.type == CollapseType.SPOKE_COLLAPSE) {
				if (relevantEdge != o.relevantEdge) {
					return true;
				}
			}
		}
		return false;
	}

	// --- Helper Methods ---
	private static boolean requiresRelevantEdge(CollapseType type) {
		return type == CollapseType.CONSTRAINT_COLLAPSE || type == CollapseType.SPOKE_COLLAPSE || type == CollapseType.SPLIT_OR_FLIP_REFINE
				|| type == CollapseType.VERTEX_MOVES_OVER_SPOKE || type == CollapseType.CCW_VERTEX_LEAVES_CH
				|| type == CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED;
	}

	private static boolean requiresSecondaryKey(CollapseType type) {
		return type == CollapseType.VERTEX_MOVES_OVER_SPOKE || type == CollapseType.FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED;
	}

	// --- Comparison Logic (Crucial for Event Queue) ---
	@Override
	public int compareTo(CollapseSpec other) {
		// Handle NEVER cases first (always highest time)
		if (this.type == CollapseType.NEVER && other.type == CollapseType.NEVER)
			return 0;
		if (this.type == CollapseType.NEVER)
			return 1; // this is later
		if (other.type == CollapseType.NEVER)
			return -1; // other is later

		// Compare time with tolerance
		double timeDiff = this.time - other.time;
		if (Math.abs(timeDiff) > 1e-9) { // Use tolerance
			return Double.compare(this.time, other.time);
		}

		// If times are effectively equal, compare by type priority (lower enum ordinal
		// = higher priority)
		int typeCompare = Integer.compare(this.type.ordinal(), other.type.ordinal());
		if (typeCompare != 0) {
			return typeCompare;
		}

		// If times and types are equal, use secondary key for tie-breaking if
		// applicable
		// Higher secondary key means HIGHER priority (comes earlier) -> reverse
		// comparison
		if (requiresSecondaryKey(this.type) && !Double.isNaN(this.secondaryKey) && !Double.isNaN(other.secondaryKey)) {
			return Double.compare(other.secondaryKey, this.secondaryKey); // Reverse compare
		}

		// If still equal, they are considered equal for priority purposes
		return 0;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("CollapseSpec{");
		sb.append(type);
		if (triangle != null)
			sb.append(", tri=").append(triangle.id);
		if (!Double.isNaN(time))
			sb.append(", time=").append(String.format("%.6f", time));
		if (relevantEdge != -1)
			sb.append(", edge=").append(relevantEdge);
		if (!Double.isNaN(secondaryKey))
			sb.append(", key=").append(String.format("%.4f", secondaryKey));
		sb.append('}');
		return sb.toString();
	}

	// Basic equals for testing (uses compareTo == 0)
	@Override
	public boolean equals(Object o) {
		if (this == o)
			return true;
		if (o == null || getClass() != o.getClass())
			return false;
		CollapseSpec that = (CollapseSpec) o;
		// Equality based on the comparison logic defined in compareTo
		return this.compareTo(that) == 0 && Objects.equals(triangle, that.triangle); // Include triangle ID check for stricter test equality
	}

	@Override
	public int hashCode() {
		// Hash based on fields relevant to priority, excluding secondary key due to
		// tolerance issues
		// Triangle ID is included for uniqueness in tests
		return Objects.hash(type, (long) (time / 1e-9), triangle != null ? triangle.id : 0, relevantEdge);
	}
}