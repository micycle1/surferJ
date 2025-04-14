package com.github.micycle.surferj2.collapse;

import java.util.Objects;

public class EdgeCollapseSpec {

	private final EdgeCollapseType type;
	private final double time; // Meaningful only for FUTURE and ALWAYS

	public static final EdgeCollapseSpec NEVER = new EdgeCollapseSpec(EdgeCollapseType.NEVER, Double.NaN);
	public static final EdgeCollapseSpec PAST = new EdgeCollapseSpec(EdgeCollapseType.PAST, Double.NaN);
	public static final EdgeCollapseSpec UNDEFINED = new EdgeCollapseSpec(EdgeCollapseType.UNDEFINED, Double.NaN);

	public EdgeCollapseSpec(EdgeCollapseType type, double time) {
		this.type = type;
		// Ensure time is non-negative for relevant types
		if ((type == EdgeCollapseType.FUTURE || type == EdgeCollapseType.ALWAYS) && (Double.isNaN(time) || time < 0.0)) {
			// Allow a small tolerance for floating point comparisons near zero
			if (Math.abs(time) > 1e-9) {
				throw new IllegalArgumentException("Time must be non-negative for FUTURE/ALWAYS collapse types, got: " + time);
			}
			this.time = 0.0; // Clamp small negative values to zero
		} else {
			this.time = time;
		}

		if (type == EdgeCollapseType.NEVER || type == EdgeCollapseType.PAST || type == EdgeCollapseType.UNDEFINED) {
			if (!Double.isNaN(this.time)) {
				// System.err.println("Warning: Time should be NaN for NEVER/PAST/UNDEFINED
				// EdgeCollapseType, but was " + this.time);
				// Allow time for debugging, but it's not strictly used for these types
			}
		}
	}

	public EdgeCollapseType getType() {
		return type;
	}

	public double getTime() {
		if (type != EdgeCollapseType.FUTURE && type != EdgeCollapseType.ALWAYS) {
			// Or throw exception? Depends on usage.
			// System.err.println("Warning: Accessing time for EdgeCollapseType " + type);
			return Double.NaN;
		}
		return time;
	}

	@Override
	public String toString() {
		if (Double.isNaN(time))
			return "EdgeCollapseSpec{" + type + '}';
		return "EdgeCollapseSpec{" + type + ", time=" + String.format("%.6f", time) + '}';
	}

	// Basic equals/hashCode for testing
	@Override
	public boolean equals(Object o) {
		if (this == o)
			return true;
		if (o == null || getClass() != o.getClass())
			return false;
		EdgeCollapseSpec that = (EdgeCollapseSpec) o;
		// Use tolerance for time comparison if type requires time
		boolean timeEquals;
		if (type == EdgeCollapseType.FUTURE || type == EdgeCollapseType.ALWAYS) {
			timeEquals = Math.abs(time - that.time) < 1e-9; // Tolerance
		} else {
			timeEquals = true; // Time doesn't matter for other types for equality check
		}
		return type == that.type && timeEquals;
	}

	@Override
	public int hashCode() {
		// Time hash code isn't reliable with tolerance, focus on type for basic hash
		return Objects.hash(type);
	}
}