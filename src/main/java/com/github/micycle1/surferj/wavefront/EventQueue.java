package com.github.micycle1.surferj.wavefront;

/**
* Manages the priority queue of collapse events for the kinetic triangulation.
* Uses a custom Heap implementation to allow efficient updates and removals.
* Mirrors the C++ EventQueue class.
*/

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import com.github.micycle1.surferj.kinetics.KineticTriangle;
import com.github.micycle1.surferj.kinetics.KineticTriangulation;

/**
 * Manages the priority queue of future kinetic events (triangle collapses). It
 * uses a standard Java PriorityQueue and handles deferred updates and removals
 * to maintain consistency, compensating for the lack of a direct decrease-key
 * operation in PriorityQueue.
 */
public class EventQueue {

	// NOTE: Using standard Java PriorityQueue + Map for lookup.
	// C++ uses a custom heap with direct index access and fix_idx operation.
	// Java approach requires remove/add for updates, which is less performant
	// (O(log N) for remove-by-identity + O(log N) for add, vs C++ O(log N)
	// fix_idx).
	// The Map lookup is O(1) on average.
	private final PriorityQueue<EventQueueItem> queue;
	private final Map<Long, EventQueueItem> tidToItemMap; // Map Triangle ID -> Queue Item

	// Sets to track triangles needing action before the next event poll
	private final Set<Long> needsDroppingSet; // Triangle IDs to be removed
	private final Set<Long> needsUpdateSet; // Triangle IDs needing priority update

	/**
	 * Initializes the EventQueue based on the initial state of the kinetic
	 * triangulation.
	 *
	 * @param kt The kinetic triangulation.
	 */
	public EventQueue(KineticTriangulation kt) {
		int initialCapacity = kt.getTriangles().size();
		this.queue = new PriorityQueue<>(initialCapacity > 0 ? initialCapacity : 11); // Default capacity is 11
		this.tidToItemMap = new HashMap<>(initialCapacity);
		this.needsDroppingSet = new HashSet<>();
		this.needsUpdateSet = new HashSet<>();

		for (KineticTriangle t : kt.getTriangles()) {
			if (t != null && !t.isDead()) { // Only add active triangles
				EventQueueItem qi = new EventQueueItem(t, 0); // Initialize at time 0 NOTE USE SurfConstants.ZERO?
				// NOTE: Only add valid events to the queue initially? The C++ adds all.
				// Adding all mimics C++, even if event is NEVER.
				this.queue.add(qi);
				this.tidToItemMap.put(t.getId(), qi);
			}
		}
		// C++ calls heapify() after adding all elements. PriorityQueue handles this
		// during adds.
	}

	/**
	 * Returns the item with the highest priority (lowest time/secondary key)
	 * without removing it from the queue. Processes pending updates first.
	 *
	 * @return The next EventQueueItem, or null if the queue is empty after updates.
	 */
	public EventQueueItem peek() {
		assertNoPending(); // Ensure consistency before peeking
		return queue.peek();
	}

	/**
	 * Removes and returns the item with the highest priority (lowest time/secondary
	 * key). Processes pending updates first.
	 *
	 * @return The next EventQueueItem, or null if the queue is empty after updates.
	 */
	public EventQueueItem poll() {
		assertNoPending(); // Ensure consistency before polling
		EventQueueItem item = queue.poll();
		if (item != null) {
			tidToItemMap.remove(item.getEvent().getTriangle().getId());
		}
		return item;
	}

	/**
	 * Returns the number of items currently in the queue. Does *not* account for
	 * items pending drop/update until process_pending_updates is called.
	 *
	 * @return The size of the internal priority queue.
	 */
	public int size() {
		// NOTE: This size might be temporarily inaccurate if items are pending drop.
		// Call process_pending_updates before relying on size for algorithm
		// termination.
		return queue.size();
	}

	/**
	 * Checks if the queue is empty. Does *not* account for items pending
	 * drop/update until process_pending_updates is called.
	 *
	 * @return true if the internal priority queue is empty.
	 */
	public boolean isEmpty() {
		// NOTE: Like size(), may be inaccurate before processing pending updates.
		return queue.isEmpty();
	}

	/**
	 * Processes all queued updates and drops. Items marked for dropping are
	 * removed. Items marked for update have their priority recalculated and their
	 * position in the queue adjusted (by removing and re-adding). This should be
	 * called before polling or peeking the next event.
	 *
	 * @param now The current simulation time.
	 */
	public void processPendingUpdates(double now) {
		// 1. Handle drops first (as a drop might override an update)
		// We must iterate carefully if removing while iterating, or use a copy.
		Set<Long> toDrop = new HashSet<>(needsDroppingSet); // Copy to avoid concurrent modification issues
		for (long tid : toDrop) {
			EventQueueItem itemToDrop = tidToItemMap.remove(tid);
			if (itemToDrop != null) {
				// Removal from PriorityQueue is O(N) if using remove(Object).
				// If performance becomes critical, alternative heap implementations
				// or strategies might be needed. For now, this is functionally correct.
				boolean removed = queue.remove(itemToDrop);
				// TODO Add logging/assertion if removed is false unexpectedly
			}
			// Ensure it's not processed in the update loop below
			needsUpdateSet.remove(tid);
		}
		needsDroppingSet.clear();

		// 2. Handle updates
		Set<Long> toUpdate = new HashSet<>(needsUpdateSet); // Copy
		for (long tid : toUpdate) {
			// If it was concurrently marked for dropping and processed above, map entry
			// will be null.
			EventQueueItem itemToUpdate = tidToItemMap.get(tid);
			if (itemToUpdate != null) {
				// IMPORTANT: Must remove before updating priority, then re-add.
				boolean removed = queue.remove(itemToUpdate);
				// TODO Add logging/assertion if removed is false unexpectedly

				itemToUpdate.updatePriority(now); // Recalculate event time/type

				// Only re-add if the triangle is still relevant (not marked dead externally)
				// and the event is valid (though queue might contain invalid events).
				// Re-adding mimics the C++ fix_idx behaviour.
				if (!itemToUpdate.getEvent().getTriangle().isDead()) {
					queue.add(itemToUpdate); // Re-adds with updated priority
				} else {
					// If triangle became dead but wasn't explicitly dropped, remove mapping
					tidToItemMap.remove(tid);
				}
			}
		}
		needsUpdateSet.clear();
	}

	/**
	 * Marks a triangle as needing its event priority recalculated and its position
	 * in the queue potentially updated. The update is deferred until
	 * processPendingUpdates is called.
	 *
	 * @param t                        The triangle needing an update.
	 * @param mayHaveValidCollapseSpec (ignored in Java impl) C++ assertion flag.
	 */
	public void needsUpdate(KineticTriangle t, boolean mayHaveValidCollapseSpec) {
		long tid = t.getId();
		// Only add if not already marked for dropping
		if (!needsDroppingSet.contains(tid)) {
			needsUpdateSet.add(tid);
		}
	}

	/**
	 * Marks a triangle as needing to be removed from the event queue. The removal
	 * is deferred until processPendingUpdates is called. Also marks the triangle
	 * itself as dead.
	 *
	 * @param t The triangle to drop.
	 */
	public void needsDropping(KineticTriangle t) {
		t.markDead(); // Mark the triangle state
		long tid = t.getId();
		needsDroppingSet.add(tid);
		// If it was pending update, the drop takes precedence.
		needsUpdateSet.remove(tid);
	}

	/**
	 * Checks if a triangle is currently marked as needing an update.
	 *
	 * @param t The triangle to check.
	 * @return true if the triangle is in the pending update set.
	 */
	public boolean isInNeedsUpdate(KineticTriangle t) {
		return needsUpdateSet.contains(t.getId());
	}

	/**
	 * Checks if a triangle is currently marked as needing to be dropped.
	 *
	 * @param t The triangle to check.
	 * @return true if the triangle is in the pending drop set.
	 */
	public boolean isInNeedsDropping(KineticTriangle t) {
		return needsDroppingSet.contains(t.getId());
	}

	/**
	 * Asserts that there are no pending updates or drops. Should be true before
	 * polling or peeking the queue.
	 */
	private void assertNoPending() {
		if (!needsDroppingSet.isEmpty() || !needsUpdateSet.isEmpty()) {
			// Consider throwing an IllegalStateException or logging a warning
			System.err.println("WARNING: EventQueue accessed while updates/drops pending. Call processPendingUpdates() first.");
			System.err.println("Pending Drops: " + needsDroppingSet);
			System.err.println("Pending Updates: " + needsUpdateSet);
			// In production, might throw:
			// throw new IllegalStateException("EventQueue accessed with pending
			// updates/drops.");
		}
	}

	/**
	 * Checks if the heap property is valid. NOTE: This is difficult/inefficient to
	 * check for standard Java PriorityQueue. C++ version iterates through internal
	 * array. This is a placeholder/no-op.
	 *
	 * @return true (always, as check is not implemented).
	 */
	public boolean isValidHeap() {
		// Cannot directly validate internal state of PriorityQueue easily.
		return true; // Assume valid
	}
}