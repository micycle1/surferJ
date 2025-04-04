package com.github.micycle.surferj.event;

import java.util.HashMap;
import java.util.Map;
import java.util.PriorityQueue;

import com.github.micycle.surferj.kinetics.KineticTriangle;

public class EventQueue {

    private PriorityQueue<CollapseEvent> queue;
    // Map to quickly find the current event for a triangle if it's in the queue
    private Map<Integer, CollapseEvent> triangleEventMap;

    public EventQueue() {
        this.queue = new PriorityQueue<>();
        this.triangleEventMap = new HashMap<>();
    }

    public void add(CollapseEvent event) {
        if (event.getType() != EventType.NONE) {
            // Remove existing event for this triangle before adding new one
             removeEventForTriangle(event.getTriangle());
            if (queue.add(event)) { // add returns boolean in PriorityQueue
                triangleEventMap.put(event.getTriangle().id, event);
            }
        }
    }

    public CollapseEvent peek() {
        // Need to ensure the top event is still valid (triangle not stopped)
        while (queue.peek() != null && queue.peek().getTriangle().hasStopped()) {
            removeInternal(queue.poll()); // Remove stopped triangle's event
        }
        return queue.peek();
    }

    public CollapseEvent poll() {
         // Ensure top event is valid before polling
         while (queue.peek() != null && queue.peek().getTriangle().hasStopped()) {
             removeInternal(queue.poll());
         }
        CollapseEvent event = queue.poll();
        removeInternal(event);
        return event;
    }

     // Update event for a given triangle
     public void update(KineticTriangle triangle, double timeNow) {
         if (triangle.hasStopped()) {
              removeEventForTriangle(triangle); // Remove if stopped
             return;
         }
         // Calculate the new event
         CollapseEvent newEvent = triangle.computeCollapseEvent(timeNow);
         // Add the new event (add method handles removing the old one)
         add(newEvent);
     }

      // Remove event associated with a specific triangle
      private void removeEventForTriangle(KineticTriangle triangle) {
           CollapseEvent existingEvent = triangleEventMap.get(triangle.id);
           if (existingEvent != null) {
               removeInternal(existingEvent);
               // Standard PriorityQueue doesn't support efficient removal by element.
               // This requires iterating or using a different queue implementation
               // (like pairing heap) or marking events as invalid.
               // HACK: Mark as invalid - peek/poll will skip them.
               // OR: Rebuild queue (inefficient) or use remove (O(n)).
               // Using remove() for simplicity here, but be aware of performance.
               queue.remove(existingEvent);
           }
      }

      // Internal removal helper
      private void removeInternal(CollapseEvent event) {
           if (event != null) {
               triangleEventMap.remove(event.getTriangle().id);
           }
      }


    public boolean isEmpty() {
        // Need to account for potentially invalid (stopped) events at the top
        return peek() == null;
    }

    public int size() {
        // This size includes potentially invalid events until peek/poll cleans them
        return queue.size();
    }
}