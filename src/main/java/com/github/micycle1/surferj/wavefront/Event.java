package com.github.micycle1.surferj.wavefront;

import java.util.Objects;

import com.github.micycle1.surferj.collapse.CollapseSpec;
import com.github.micycle1.surferj.kinetics.KineticTriangle;

/**
 * Represents a specific event instance occurring for a particular kinetic
 * triangle. It extends CollapseSpec, inheriting its properties (type, time,
 * relevant edge, secondary key, component, triangle) and comparison logic.
 * Primarily used within the EventQueueItem.
 */
public class Event extends CollapseSpec {

	/**
	 * Creates a new Event based on the collapse specification calculated for the
	 * given triangle at the specified time. The necessary fields (type, time, edge,
	 * key, component, triangle) are initialized by the superclass constructor using
	 * the calculated spec.
	 *
	 * @param triangle The kinetic triangle this event pertains to. Must not be
	 *                 null.
	 * @param now      The current time, used to calculate the initial collapse
	 *                 spec.
	 */
	public Event(KineticTriangle triangle, double now) {
		super(calculateInitialSpec(triangle, now), triangle);

	}

	private static CollapseSpec calculateInitialSpec(KineticTriangle triangle, double now) {
		Objects.requireNonNull(triangle, "Event triangle cannot be null during initial spec calculation");
		CollapseSpec spec = triangle.getCollapseSpec(now);
		Objects.requireNonNull(spec, "KineticTriangle.getCollapseSpec returned null for triangle " + triangle.getId());
		return spec;
	}

	@Override
	public String toString() {
		// Refine toString slightly to indicate it's an Event wrapping a CollapseSpec's
		// info
		return "Event" + super.toString().substring("CollapseSpec".length()); // Removes "CollapseSpec" prefix
		// Example output: Event{TRIANGLE_COLLAPSE, tri=12, comp=0, time=1.234567}
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		// Check if 'o' is an Event instance specifically
		if (!(o instanceof Event)) {
			return false;
		}
		// Rely on CollapseSpec's equals definition (which uses compareTo == 0 for
		// priority equality)
		return super.equals(o);
		// If stricter equality (same triangle instance/ID needed beyond priority):
		// return super.equals(o) && Objects.equals(this.getTriangle(), ((Event)
		// o).getTriangle());
	}

	@Override
	public int hashCode() {
		// Rely on CollapseSpec's hashCode for priority-based hashing.
		// If stricter hashing needed, combine with triangle ID.
		return super.hashCode();
		// Stricter version:
		// KineticTriangle tri = getTriangle();
		// return Objects.hash(super.hashCode(), (tri != null ? tri.getId() : 0));
	}

}