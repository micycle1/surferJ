//package com.github.micycle.surferj2.events;
//
//import org.junit.jupiter.api.BeforeEach;
//import org.junit.jupiter.api.Test;
//import org.locationtech.jts.geom.Coordinate;
//
//import com.github.micycle.surferj2.collapse.CollapseSpec;
//import com.github.micycle.surferj2.collapse.CollapseType;
//import com.github.micycle.surferj2.kinetics.KineticTriangle;
//
//import java.util.ArrayList;
//import java.util.List;
//
//import static org.junit.jupiter.api.Assertions.*;
//
//public class EventQueueTest {
//
//    private EventQueue eventQueue;
//    private List<KineticTriangle> triangles;
//    private KineticTriangle t1, t2, t3;
//
//    @BeforeEach
//    void setUp() {
//        // Create some dummy triangles for testing
//        t1 = new KineticTriangle(); // ID will be 1 (assuming fresh counter)
//        t2 = new KineticTriangle(); // ID will be 2
//        t3 = new KineticTriangle(); // ID will be 3
//        triangles = new ArrayList<>(List.of(t1, t2, t3));
//
//        // Initialize EventQueue with an empty list, we'll add events manually
//        eventQueue = new EventQueue(new ArrayList<>(), 0.0);
//    }
//
//    // Helper to create a basic spec
//    private CollapseSpec createSpec(CollapseType type, double time, KineticTriangle tri, int edge) {
//        return new CollapseSpec(type, time, tri, edge, Double.NaN); // Secondary key NaN unless specified
//    }
//     private CollapseSpec createSpecWithKey(CollapseType type, double time, KineticTriangle tri, int edge, double key) {
//        return new CollapseSpec(type, time, tri, edge, key);
//    }
//
//    @Test
//    void testEmptyQueueOnInit() {
//        EventQueue eq = new EventQueue(new ArrayList<>(), 0.0);
//        assertNull(eq.pollNextValidEvent());
//        assertFalse(eq.hasMoreEvents());
//        assertEquals(0, eq.size());
//    }
//
//    @Test
//    void testAddAndPollSimpleOrder() {
//        CollapseSpec specA = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t1, 0);
//        CollapseSpec specB = createSpec(CollapseType.SPOKE_COLLAPSE, 2.0, t2, 1);
//        CollapseSpec specC = createSpec(CollapseType.TRIANGLE_COLLAPSE, 8.0, t3, -1); // No edge needed
//
//        eventQueue.addOrUpdate(specA);
//        eventQueue.addOrUpdate(specB);
//        eventQueue.addOrUpdate(specC);
//
//        assertEquals(3, eventQueue.size());
//        assertTrue(eventQueue.hasMoreEvents());
//
//        assertEquals(specB, eventQueue.pollNextValidEvent());
//        assertEquals(2, eventQueue.size());
//        assertEquals(specA, eventQueue.pollNextValidEvent());
//        assertEquals(1, eventQueue.size());
//        assertEquals(specC, eventQueue.pollNextValidEvent());
//        assertEquals(0, eventQueue.size());
//        assertFalse(eventQueue.hasMoreEvents());
//        assertNull(eventQueue.pollNextValidEvent());
//    }
//
//    @Test
//    void testTimeTieBreakingByType() {
//        // CONSTRAINT_COLLAPSE has lower ordinal (higher priority) than SPOKE_COLLAPSE
//        CollapseSpec specA = createSpec(CollapseType.SPOKE_COLLAPSE, 5.0, t1, 1);
//        CollapseSpec specB = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t2, 0);
//
//        eventQueue.addOrUpdate(specA);
//        eventQueue.addOrUpdate(specB);
//
//        assertEquals(specB, eventQueue.pollNextValidEvent(), "CONSTRAINT_COLLAPSE should have higher priority");
//        assertEquals(specA, eventQueue.pollNextValidEvent());
//        assertNull(eventQueue.pollNextValidEvent());
//    }
//
//     @Test
//    void testTimeAndTypeTieBreakingBySecondaryKey() {
//        // VERTEX_MOVES_OVER_SPOKE uses secondary key (higher key = higher priority)
//        CollapseSpec specA = createSpecWithKey(CollapseType.VERTEX_MOVES_OVER_SPOKE, 5.0, t1, 2, 10.0); // Lower key
//        CollapseSpec specB = createSpecWithKey(CollapseType.VERTEX_MOVES_OVER_SPOKE, 5.0, t2, 1, 20.0); // Higher key
//
//        eventQueue.addOrUpdate(specA);
//        eventQueue.addOrUpdate(specB);
//
//        assertEquals(specB, eventQueue.pollNextValidEvent(), "Higher secondary key should have higher priority");
//        assertEquals(specA, eventQueue.pollNextValidEvent());
//        assertNull(eventQueue.pollNextValidEvent());
//    }
//
//     @Test
//    void testAddNeverEvent() {
//        CollapseSpec specA = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t1, 0);
//        CollapseSpec specNever = new CollapseSpec(CollapseType.NEVER, Double.NaN, t2, -1, Double.NaN);
//
//        eventQueue.addOrUpdate(specA);
//        eventQueue.addOrUpdate(specNever); // Should not be added
//
//        assertEquals(1, eventQueue.size());
//        assertEquals(specA, eventQueue.pollNextValidEvent());
//        assertNull(eventQueue.pollNextValidEvent());
//        assertEquals(0, eventQueue.size());
//    }
//
//     @Test
//    void testInvalidateEvent() {
//        CollapseSpec specA = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t1, 0);
//        CollapseSpec specB = createSpec(CollapseType.SPOKE_COLLAPSE, 2.0, t2, 1);
//
//        eventQueue.addOrUpdate(specA);
//        eventQueue.addOrUpdate(specB);
//        assertEquals(2, eventQueue.size());
//
//        // Invalidate specB (the earlier one)
//        eventQueue.invalidateEventFor(t2); // This marks specB instance as invalid
//
//        assertEquals(1, eventQueue.size(), "Size should reflect only active events after invalidation"); // t1 is still active
//        assertFalse(specB.isValid(), "SpecB instance should be marked invalid");
//        assertTrue(specA.isValid(), "SpecA instance should still be valid");
//
//        // Polling should skip the invalidated specB and return specA
//        assertEquals(specA, eventQueue.pollNextValidEvent());
//        assertEquals(0, eventQueue.size(), "Queue should be empty after polling the only valid event");
//        assertNull(eventQueue.pollNextValidEvent());
//    }
//
//    @Test
//    void testUpdateEvent_ReplaceWithLaterTime() {
//        CollapseSpec specA_T1_time2 = createSpec(CollapseType.SPOKE_COLLAPSE, 2.0, t1, 1);
//        CollapseSpec specB_T2_time5 = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t2, 0);
//        CollapseSpec specC_T1_time8 = createSpec(CollapseType.TRIANGLE_COLLAPSE, 8.0, t1, -1); // New event for t1
//
//        eventQueue.addOrUpdate(specA_T1_time2);
//        eventQueue.addOrUpdate(specB_T2_time5);
//        assertEquals(2, eventQueue.size());
//        assertEquals(specA_T1_time2, eventQueue.peekNextValidEvent()); // T1 at t=2 is next
//
//        // Update event for T1
//        eventQueue.addOrUpdate(specC_T1_time8); // This invalidates specA, adds specC
//
//        assertEquals(2, eventQueue.size(), "Size should remain 2 (B and C active)");
//        assertFalse(specA_T1_time2.isValid(), "Original specA should be invalid after update");
//        assertTrue(specC_T1_time8.isValid(), "New specC should be valid");
//
//        // Now specB should be the first event
//        assertEquals(specB_T2_time5, eventQueue.pollNextValidEvent());
//        assertEquals(1, eventQueue.size());
//        // Then the updated specC for T1
//        assertEquals(specC_T1_time8, eventQueue.pollNextValidEvent());
//        assertEquals(0, eventQueue.size());
//        assertNull(eventQueue.pollNextValidEvent());
//    }
//
//     @Test
//    void testUpdateEvent_ReplaceWithEarlierTime() {
//        CollapseSpec specA_T1_time8 = createSpec(CollapseType.TRIANGLE_COLLAPSE, 8.0, t1, -1);
//        CollapseSpec specB_T2_time5 = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t2, 0);
//        CollapseSpec specC_T1_time2 = createSpec(CollapseType.SPOKE_COLLAPSE, 2.0, t1, 1); // New, earlier event for t1
//
//        eventQueue.addOrUpdate(specA_T1_time8);
//        eventQueue.addOrUpdate(specB_T2_time5);
//        assertEquals(2, eventQueue.size());
//        assertEquals(specB_T2_time5, eventQueue.peekNextValidEvent()); // T2 at t=5 is next
//
//        // Update event for T1
//        eventQueue.addOrUpdate(specC_T1_time2); // This invalidates specA, adds specC
//
//        assertEquals(2, eventQueue.size(), "Size should remain 2 (B and C active)");
//        assertFalse(specA_T1_time8.isValid(), "Original specA should be invalid after update");
//        assertTrue(specC_T1_time2.isValid(), "New specC should be valid");
//
//        // Now specC (new T1 event) should be the first event
//        assertEquals(specC_T1_time2, eventQueue.pollNextValidEvent());
//        assertEquals(1, eventQueue.size());
//        // Then the original specB for T2
//        assertEquals(specB_T2_time5, eventQueue.pollNextValidEvent());
//        assertEquals(0, eventQueue.size());
//        assertNull(eventQueue.pollNextValidEvent());
//    }
//
//     @Test
//    void testUpdateWithNeverInvalidatesPrevious() {
//        CollapseSpec specA_T1_time5 = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t1, 0);
//        CollapseSpec specNever_T1 = new CollapseSpec(CollapseType.NEVER, Double.NaN, t1, -1, Double.NaN);
//
//        eventQueue.addOrUpdate(specA_T1_time5);
//        assertEquals(1, eventQueue.size());
//        assertTrue(specA_T1_time5.isValid());
//
//        eventQueue.addOrUpdate(specNever_T1); // Add NEVER event for T1
//
//        assertEquals(0, eventQueue.size(), "Adding NEVER should remove active event");
//        assertFalse(specA_T1_time5.isValid(), "Original spec should be invalidated by NEVER");
//        assertNull(eventQueue.pollNextValidEvent(), "Queue should be empty");
//    }
//
//     @Test
//    void testPollSkipsInvalidatedAndSuperseded() {
//        // Sequence: Add A(t1, t=2), Add B(t2, t=5), Add C(t3, t=8)
//        // Then: Invalidate A, Update C with D(t3, t=1)
//        CollapseSpec specA_T1_time2 = createSpec(CollapseType.SPOKE_COLLAPSE, 2.0, t1, 1);
//        CollapseSpec specB_T2_time5 = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t2, 0);
//        CollapseSpec specC_T3_time8 = createSpec(CollapseType.TRIANGLE_COLLAPSE, 8.0, t3, -1);
//        CollapseSpec specD_T3_time1 = createSpec(CollapseType.SPOKE_COLLAPSE, 1.0, t3, 0); // New, earlier for T3
//
//        eventQueue.addOrUpdate(specA_T1_time2); // Queue: [A(2)] Map: {t1:A}
//        eventQueue.addOrUpdate(specB_T2_time5); // Queue: [A(2), B(5)] Map: {t1:A, t2:B}
//        eventQueue.addOrUpdate(specC_T3_time8); // Queue: [A(2), B(5), C(8)] Map: {t1:A, t2:B, t3:C}
//
//        eventQueue.invalidateEventFor(t1);      // Queue: [A(2, invalid), B(5), C(8)] Map: {t2:B, t3:C}
//        assertFalse(specA_T1_time2.isValid());
//        assertEquals(2, eventQueue.size());
//
//        eventQueue.addOrUpdate(specD_T3_time1); // Queue: [D(1), A(2, invalid), B(5), C(8, invalid)] Map: {t2:B, t3:D}
//        assertFalse(specC_T3_time8.isValid());
//        assertTrue(specD_T3_time1.isValid());
//        assertEquals(2, eventQueue.size()); // Still only t2 and t3 have active events
//
//        // Expected polling order: D(1), B(5)
//        assertEquals(specD_T3_time1, eventQueue.pollNextValidEvent()); // Polls D, skips nothing
//        assertEquals(1, eventQueue.size());
//        assertEquals(specB_T2_time5, eventQueue.pollNextValidEvent()); // Polls B, skips A(invalid) and C(invalid) implicitly
//        assertEquals(0, eventQueue.size());
//        assertNull(eventQueue.pollNextValidEvent()); // Queue is now effectively empty
//    }
//
//     @Test
//    void testPeekNextValidEvent() {
//        CollapseSpec specA_T1_time2_inv = createSpec(CollapseType.SPOKE_COLLAPSE, 2.0, t1, 1); // Will invalidate
//        CollapseSpec specB_T2_time5 = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t2, 0);
//        CollapseSpec specC_T1_time8 = createSpec(CollapseType.TRIANGLE_COLLAPSE, 8.0, t1, -1); // Supersedes A
//
//        eventQueue.addOrUpdate(specA_T1_time2_inv);
//        eventQueue.addOrUpdate(specB_T2_time5);
//        eventQueue.addOrUpdate(specC_T1_time8); // Invalidates A, adds C
//
//        assertEquals(2, eventQueue.size());
//
//        // Peek should skip invalidated A and return B
//        assertEquals(specB_T2_time5, eventQueue.peekNextValidEvent());
//        assertEquals(2, eventQueue.size(), "Peek should not change size"); // Ensure peek didn't remove
//
//        // Poll B
//        assertEquals(specB_T2_time5, eventQueue.pollNextValidEvent());
//        assertEquals(1, eventQueue.size());
//
//        // Peek should now return C
//        assertEquals(specC_T1_time8, eventQueue.peekNextValidEvent());
//        assertEquals(1, eventQueue.size());
//
//        // Poll C
//        assertEquals(specC_T1_time8, eventQueue.pollNextValidEvent());
//        assertEquals(0, eventQueue.size());
//
//        // Peek should return null
//        assertNull(eventQueue.peekNextValidEvent());
//    }
//
//    @Test
//    void testClearQueue() {
//        CollapseSpec specA = createSpec(CollapseType.CONSTRAINT_COLLAPSE, 5.0, t1, 0);
//        CollapseSpec specB = createSpec(CollapseType.SPOKE_COLLAPSE, 2.0, t2, 1);
//        eventQueue.addOrUpdate(specA);
//        eventQueue.addOrUpdate(specB);
//        assertEquals(2, eventQueue.size());
//
//        eventQueue.clear();
//
//        assertEquals(0, eventQueue.size());
//        assertFalse(eventQueue.hasMoreEvents());
//        assertNull(eventQueue.pollNextValidEvent());
//        // Check internal state if possible (e.g., if activeEvents was accessible)
//    }
//}