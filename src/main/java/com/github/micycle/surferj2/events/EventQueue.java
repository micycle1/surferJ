package com.github.micycle.surferj2.events;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.PriorityQueue;

import com.github.micycle.surferj2.collapse.CollapseSpec;
import com.github.micycle.surferj2.collapse.CollapseType;
import com.github.micycle.surferj2.kinetics.KineticTriangle;

/**
 * Manages future collapse events for the Kinetic Triangulation simulation.
 * Uses a PriorityQueue ordered by event time and type.
 * Implements a "Mark as Invalid + Re-add" strategy for handling updates,
 * requiring an external map to track the currently active event for each triangle.
 */
public class EventQueue {

    private final PriorityQueue<CollapseSpec> queue;
    // Tracks the single, currently active CollapseSpec associated with each triangle in the queue.
    private final Map<KineticTriangle, CollapseSpec> activeEvents; // NOTE check equality and hashCode are implemented properly?

    /**
     * Initializes the Event Queue by calculating the initial collapse event
     * for every triangle in the provided list.
     *
     * @param initialTriangles The list of KineticTriangles at the start (t=0).
     * @param initialTime The time at which to calculate initial events (usually 0.0).
     */
    public EventQueue(List<KineticTriangle> initialTriangles, double initialTime) {
        this.queue = new PriorityQueue<>();
        this.activeEvents = new HashMap<>();

        if (initialTriangles == null) {
            return; // Handle null input gracefully
        }

        for (KineticTriangle triangle : initialTriangles) {
            if (triangle == null) continue; // Skip null triangles if any
            // Calculate initial event
            CollapseSpec initialSpec = triangle.computeCollapse(initialTime);
            // Add if it's a valid future event (not NEVER)
            addOrUpdate(initialSpec);
        }
    }

    /**
     * Adds a new potential event or updates the existing event for the
     * associated triangle. If an older event for the same triangle exists
     * in the queue, it is marked as invalid. Only non-NEVER events are
     * actually added.
     *
     * @param newSpec The new CollapseSpec to consider. Must not be null.
     */
    public void addOrUpdate(CollapseSpec newSpec) {
        Objects.requireNonNull(newSpec, "Cannot add a null CollapseSpec");
        KineticTriangle triangle = newSpec.getTriangle();

        if (triangle == null || newSpec.getType() == CollapseType.NEVER) {
             // If the event is NEVER or has no associated triangle, ensure any
             // previous event for this (potentially null) triangle is invalidated.
             if (triangle != null) {
                 invalidateEventFor(triangle); // Remove from map, mark old spec invalid
             }
            return; // Don't add NEVER events to the queue
        }

        // 1. Invalidate any existing active event for this triangle
        invalidateEventFor(triangle); // Marks old spec invalid, removes from map

        // 2. Add the new spec to the queue and the tracking map
        if (newSpec.getType() != CollapseType.NEVER) { // Double-check, though handled above
             if (!newSpec.isValid()) {
                // This shouldn't happen if invalidate logic is correct, but safety check
                System.err.println("Warning: Attempting to add an already invalid spec: " + newSpec);
                newSpec = triangle.computeCollapse(newSpec.getTime()); // Recompute maybe? Or just ignore? For now ignore.
                if (newSpec.getType() == CollapseType.NEVER || !newSpec.isValid()) return;
            }
            queue.add(newSpec);
            activeEvents.put(triangle, newSpec);
        }
    }

    /**
     * Marks the currently active event associated with the given triangle
     * as invalid and removes it from the active tracking map. The invalid
     * event object remains in the priority queue but will be skipped during polling.
     *
     * @param triangle The triangle whose event should be invalidated.
     */
    public void invalidateEventFor(KineticTriangle triangle) {
        if (triangle == null) return;

        CollapseSpec currentlyActiveSpec = activeEvents.remove(triangle);
        if (currentlyActiveSpec != null) {
            currentlyActiveSpec.invalidate();
            // Note: The invalidated spec stays in the PriorityQueue until polled.
        }
    }

    /**
     * Retrieves and removes the next valid event (earliest time, highest priority)
     * from the queue. Invalidated events are skipped.
     *
     * @return The next valid CollapseSpec, or null if the queue is effectively empty.
     */
    public CollapseSpec pollNextValidEvent() {
        while (!queue.isEmpty()) {
            CollapseSpec polledSpec = queue.poll();

            // 1. Check if the polled instance itself was marked invalid
            if (!polledSpec.isValid()) {
                // This spec was invalidated before being polled, skip it.
                continue;
            }

            KineticTriangle triangle = polledSpec.getTriangle();
            if (triangle == null) {
                 // Should not happen for valid specs added, but safety check
                 System.err.println("Warning: Polled valid spec with null triangle: " + polledSpec);
                 continue;
            }

            // 2. Check if this polled spec is still the *active* one for its triangle
            CollapseSpec activeSpecForTriangle = activeEvents.get(triangle);
            if (polledSpec == activeSpecForTriangle) {
                // It's valid and it's the one we are tracking! Process it.
                activeEvents.remove(triangle); // Remove from tracking map as it's being processed
                return polledSpec;
            } else {
                // This polled spec is valid, but it's NOT the one currently tracked
                // in the activeEvents map. This means a *newer* spec for the same
                // triangle was added via addOrUpdate after this one, but before
                // this older one could be marked invalid and polled.
                // Discard this older (but technically still valid) instance.
                // The newer one (in the map) is still in the queue and will be polled later.
                continue;
            }
        }
        return null; // Queue is empty or only contains invalid/superseded events
    }


    /**
     * Peeks at the next valid event without removing it.
     * WARNING: This can be less efficient than polling as it might need to
     * iterate past several invalid events at the top of the queue.
     *
     * @return The next valid CollapseSpec, or null if the queue is effectively empty.
     */
    public CollapseSpec peekNextValidEvent() {
        // NOTE This requires iterating without removing, which PriorityQueue doesn't support efficiently while checking validity.
        // A less efficient simulation:
        PriorityQueue<CollapseSpec> tempQueue = new PriorityQueue<>(queue); // Copy O(n)
        while (!tempQueue.isEmpty()) {
            CollapseSpec peekedSpec = tempQueue.poll(); // Removes from temp queue
             if (!peekedSpec.isValid()) {
                continue;
            }
            KineticTriangle triangle = peekedSpec.getTriangle();
             if (triangle == null) {
                 continue;
            }
             CollapseSpec activeSpecForTriangle = activeEvents.get(triangle);
            if (peekedSpec == activeSpecForTriangle) {
                 // Found the next valid event
                return peekedSpec;
            }
            // Otherwise, it's superseded, continue checking the temp queue
        }
        return null; // No valid event found
    }

    /**
     * Checks if there are likely any valid events remaining.
     * Based on the tracking map size.
     *
     * @return true if the active event map is not empty.
     */
    public boolean hasMoreEvents() {
        return !activeEvents.isEmpty();
    }

    /**
     * Returns the number of triangles that currently have an active,
     * pending event in the queue.
     *
     * @return The count of active events.
     */
    public int size() {
        return activeEvents.size();
    }

     /**
     * Clears the event queue and the tracking map.
     */
    public void clear() {
        queue.clear();
        activeEvents.clear();
    }
}