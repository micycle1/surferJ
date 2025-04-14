package com.github.micycle1.surferj.wavefront;

import java.util.Objects;

import com.github.micycle1.surferj.kinetics.KineticTriangle;

/**
 * An item placed in the EventQueue. It wraps the current Event specification
 * for a triangle and allows its priority to be updated by recalculating the
 * underlying Event. Implements Comparable based on the wrapped Event.
 */
public class EventQueueItem implements Comparable<EventQueueItem> {

	// Make the Event reference non-final so it can be replaced
	private Event event;

	/**
	 * Creates a new EventQueueItem for a given triangle at a specific time.
	 *
	 * @param t   The kinetic triangle. Must not be null.
	 * @param now The current time.
	 */
	public EventQueueItem(KineticTriangle t, double now) {
		// Initial event calculation
		this.event = new Event(Objects.requireNonNull(t, "Triangle cannot be null for EventQueueItem"), now);
	}

	/**
	 * Gets the underlying Event object reflecting the current specification.
	 *
	 * @return The current Event associated with this queue item.
	 */
	public Event getEvent() {
		return event;
	}

	/**
	 * Updates the priority of this item by creating a NEW underlying Event based on
	 * the triangle's state at the given time. The internal reference is replaced.
	 *
	 * @param now The current time.
	 */
	public void updatePriority(double now) {
		// Get the triangle from the *current* event object
		KineticTriangle triangle = this.event.getTriangle();
		if (triangle == null) {
			// Should not happen if constructor enforces non-null triangle
			throw new IllegalStateException("Cannot update priority for EventQueueItem with a null triangle reference in its event.");
		}
		// Create a completely new Event instance with the updated calculation
		// This new Event will have the correct time, type, edge, key, component based
		// on 'now'.
		this.event = new Event(triangle, now);
		// The EventQueue will handle removing this EventQueueItem and re-adding it
		// to reflect the changed priority derived from the new internal 'event'.
	}

	/**
	 * Compares this item to another based on their *current* underlying Event's
	 * comparison logic.
	 */
	@Override
	public int compareTo(EventQueueItem other) {
		// Comparison uses the currently held 'event' reference
		return this.event.compareTo(other.getEvent());
	}

	@Override
	public String toString() {
		// Show the currently held event
		return "QueueItem{" + event + '}';
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}
		EventQueueItem that = (EventQueueItem) o;
		// Equality depends on the currently held event instance's equality
		return Objects.equals(event, that.event);
	}

	@Override
	public int hashCode() {
		// Hash code depends on the currently held event instance
		return Objects.hash(event);
	}
}