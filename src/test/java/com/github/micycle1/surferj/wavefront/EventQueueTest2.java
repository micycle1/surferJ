package com.github.micycle1.surferj.wavefront;

import static org.junit.jupiter.api.Assertions.*;

import java.util.List;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import com.github.micycle1.surferj.collapse.CollapseSpec;
import com.github.micycle1.surferj.collapse.CollapseType;
import com.github.micycle1.surferj.kinetics.KineticTriangle;
import com.github.micycle1.surferj.kinetics.KineticTriangulation;
import com.github.micycle1.surferj.wavefront.EventQueue;
import com.github.micycle1.surferj.wavefront.EventQueueItem;

// o4-mini
public class EventQueueTest2 {

    /**
     * A minimal stub of KineticTriangle that:
     *  - carries an ID
     *  - can be marked dead
     *  - returns a CollapseSpec with a fixed time
     */
    private static class FakeTriangle extends KineticTriangle {

        private final long id;
        private boolean dead = false;
        private double collapseTime;

        FakeTriangle(long id, double initialTime) {
            super(); // if your real KineticTriangle has arguments, adjust here
            this.id = id;
            this.collapseTime = initialTime;
        }

        void setCollapseTime(double t) {
            this.collapseTime = t;
        }

        public void markDead() {
            this.dead = true;
        }

        @Override
        public long getId() {
            return id;
        }

        @Override
        public boolean isDead() {
            return dead;
        }

        /**
         * Return a real CollapseSpec built from your ported constructor.
         * We pick a benign CollapseType (e.g. EDGE if you prefer VERTEX).
         */
        @Override
        public CollapseSpec getCollapseSpec(double now) {
            // We ignore 'now' in this fake: the spec time is fixed.
            // Use the constructor CollapseSpec(CollapseType, double, KineticTriangle).
            return new CollapseSpec(CollapseType.VERTEX_MOVES_OVER_SPOKE, this.collapseTime, this);
        }
    }

    private FakeTriangle t1, t2, t3;
    private EventQueue queue;

    @BeforeEach
    void before() {
        // Create three fake triangles with distinct collapse times
        t1 = new FakeTriangle(1, 3.0);
        t2 = new FakeTriangle(2, 1.0);
        t3 = new FakeTriangle(3, 2.0);

        // Fake KineticTriangulation that just returns our list
        KineticTriangulation kt = new KineticTriangulation() {
            @Override
            public List<KineticTriangle> getTriangles() {
                return List.of(t1, t2, t3);
            }
            // other methods not used by EventQueue can remain unimplemented
        };

        queue = new EventQueue(kt);
    }

    @Test
    void testInitialPeekAndPollOrder() {
        // We expect times: t2@1.0, t3@2.0, t1@3.0
        assertEquals(3, queue.size());
        assertFalse(queue.isEmpty());

        EventQueueItem head = queue.peek();
        assertNotNull(head);
        assertEquals(1.0, head.getEvent().getTime());
        assertEquals(2, head.getEvent().getTriangle().getId());

        // poll off the first
        EventQueueItem polled = queue.poll();
        assertEquals(1.0, polled.getEvent().getTime());
        assertEquals(2, polled.getEvent().getTriangle().getId());

        // remaining size
        assertEquals(2, queue.size());
        // next peek is t3@2.0
        EventQueueItem next = queue.peek();
        assertEquals(2.0, next.getEvent().getTime());
        assertEquals(3, next.getEvent().getTriangle().getId());
    }

    @Test
    void testDeferredDropping() {
        // Mark t2 for drop
        queue.needsDropping(t2);
        assertTrue(queue.isInNeedsDropping(t2));

        // Until we flush, size remains 3
        assertEquals(3, queue.size());

        queue.processPendingUpdates(0.0);
        // t2 should vanish
        assertFalse(queue.isInNeedsDropping(t2));
        assertEquals(2, queue.size());

        // Now peek() is t3
        EventQueueItem head = queue.peek();
        assertEquals(3, head.getEvent().getTriangle().getId());
    }

    @Test
    void testDeferredUpdateReordersQueue() {
        // Initial collapse times: t1@3.0, t2@1.0, t3@2.0
        // Remove the earliest (t2), leaving t3@2.0, t1@3.0
        EventQueueItem removed = queue.poll();
        assertEquals(2L, removed.getEvent().getTriangle().getId());
        assertEquals(2, queue.size());

        // Before any update, peek() must be t3 (ID=3, time=2.0)
        EventQueueItem beforeUpdate = queue.peek();
        assertEquals(3L, beforeUpdate.getEvent().getTriangle().getId());
        assertEquals(2.0, beforeUpdate.getEvent().getTime(), 1e-9);

        // Now lower t1’s collapse time to 1.5 and mark for update
        t1.setCollapseTime(1.5);
        queue.needsUpdate(t1, true);
        assertTrue(queue.isInNeedsUpdate(t1));

        // STILL, until we call processPendingUpdates, peek() remains t3
        EventQueueItem stillBefore = queue.peek();
        assertEquals(3L, stillBefore.getEvent().getTriangle().getId());

        // Flush the update
        queue.processPendingUpdates(0.0);

        // After flush, t1@1.5 < t3@2.0, so head must be t1 (ID=1)
        EventQueueItem head = queue.peek();
        assertEquals(1L, head.getEvent().getTriangle().getId());
        assertEquals(1.5, head.getEvent().getTime(), 1e-9);
    }

    @Test
    void testPollEmpties() {
        // Poll three times
        assertNotNull(queue.poll());
        assertNotNull(queue.poll());
        assertNotNull(queue.poll());
        assertTrue(queue.isEmpty());
        assertEquals(0, queue.size());
        // Further poll/peek yield null
        assertNull(queue.peek());
        assertNull(queue.poll());
    }

    @Test
    void testMixedDropAndUpdate() {
        // Setup: drop t2, update t3 to 0.5
        t3.setCollapseTime(0.5);
        queue.needsUpdate(t3, true);
        queue.needsDropping(t2);

        queue.processPendingUpdates(0.0);
        // After flush we have t3@0.5 and t1@3.0
        assertEquals(2, queue.size());
        EventQueueItem head = queue.peek();
        assertEquals(0.5, head.getEvent().getTime());
        assertEquals(3, head.getEvent().getTriangle().getId());
    }
}