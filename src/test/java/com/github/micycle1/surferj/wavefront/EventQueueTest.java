package com.github.micycle1.surferj.wavefront;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.anyDouble;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.lenient;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.never;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.Mock;
import org.mockito.junit.jupiter.MockitoExtension;

import com.github.micycle1.surferj.collapse.CollapseSpec;
import com.github.micycle1.surferj.collapse.CollapseType;
import com.github.micycle1.surferj.kinetics.KineticTriangle;
import com.github.micycle1.surferj.kinetics.KineticTriangulation;

// Using MockitoExtension to handle mock lifecycle
@ExtendWith(MockitoExtension.class)
class EventQueueTest {

	private static final int T_0 = 0; // time zero

	@Mock
	private KineticTriangulation mockTriangulation;

	// Mocks for individual triangles will be created per test
	private EventQueue eventQueue;

	// Helper: Create mock triangle with lenient defaults
	private KineticTriangle createMockTriangle(long id, int component) {
		KineticTriangle mockTriangle = mock(KineticTriangle.class);
		// Use lenient() for default setups that might be overridden or not always used
		lenient().when(mockTriangle.getId()).thenReturn(id);
		lenient().when(mockTriangle.getComponent()).thenReturn(component);
		lenient().when(mockTriangle.isDead()).thenReturn(false);
		// **Crucial:** Provide a default fallback for getCollapseSpec to avoid NPEs
		// Return NEVER spec by default if no specific time is matched.
		lenient().when(mockTriangle.getCollapseSpec(anyDouble()))
				.thenReturn(new CollapseSpec(CollapseType.NEVER, Double.NaN, mockTriangle, -1, Double.NaN, component));
		return mockTriangle;
	}

	// Helper: Setup specific getCollapseSpec mocking (overrides the default
	// anyDouble)
	private void mockGetCollapseSpecTime(KineticTriangle mockTriangle, double time, CollapseSpec specToReturn) {
		// No need for lenient() here if we *expect* this specific call
		// However, using lenient() makes it more robust if test logic changes slightly
		lenient().when(mockTriangle.getCollapseSpec(eq(time))).thenReturn(specToReturn);
	}

	// Helper: Create CollapseSpec (simplified)
	private CollapseSpec createCollapseSpec(CollapseType type, double time, KineticTriangle triangle, int component) {
		return new CollapseSpec(type, time, triangle, component); // Use constructor inferring component from triangle
	}

	private CollapseSpec createCollapseSpec(CollapseType type, double time, KineticTriangle triangle, int relevantEdge, double secondaryKey) {
		// Use constructor that infers component
		return new CollapseSpec(type, time, triangle, relevantEdge, secondaryKey);
	}

	@BeforeEach
	void setUp() {
		// Test-specific mock setups happen within each @Test method
	}

	@Test
	void testInitialization_EmptyTriangulation() {
		when(mockTriangulation.getTriangles()).thenReturn(Collections.emptyList());
		eventQueue = new EventQueue(mockTriangulation);
		assertTrue(eventQueue.isEmpty());
	}

	@Test
	void testInitialization_PopulatesQueueCorrectly() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		KineticTriangle t2 = createMockTriangle(2L, 0);
		CollapseSpec spec1 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 10.0, t1, 0);
		CollapseSpec spec2 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 5.0, t2, 0);

		// *** Explicitly mock for time = CORE_ZERO for init ***
		mockGetCollapseSpecTime(t1, T_0, spec1);
		mockGetCollapseSpecTime(t2, T_0, spec2);

		when(mockTriangulation.getTriangles()).thenReturn(Arrays.asList(t1, t2));
		eventQueue = new EventQueue(mockTriangulation); // This calls getCollapseSpec(CORE_ZERO)

		assertFalse(eventQueue.isEmpty());
		assertEquals(2, eventQueue.size());
		assertNotNull(eventQueue.peek());
		assertEquals(2L, eventQueue.peek().getEvent().getTriangle().getId()); // t2 is earlier
	}

	// --- Ordering Tests --- (Need CORE_ZERO mocking)
	@Test
	void testOrdering_ByTime() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		KineticTriangle t2 = createMockTriangle(2L, 0);
		CollapseSpec spec1 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 10.0, t1, 0);
		CollapseSpec spec2 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 5.0, t2, 0);
		mockGetCollapseSpecTime(t1, T_0, spec1);
		mockGetCollapseSpecTime(t2, T_0, spec2);
		when(mockTriangulation.getTriangles()).thenReturn(Arrays.asList(t1, t2));
		eventQueue = new EventQueue(mockTriangulation);

		assertEquals(2L, eventQueue.poll().getEvent().getTriangle().getId());
		assertEquals(1L, eventQueue.poll().getEvent().getTriangle().getId());
	}

	@Test
	void testOrdering_ByComponent() {
		KineticTriangle t1 = createMockTriangle(1L, 1); // Higher component
		KineticTriangle t2 = createMockTriangle(2L, 0); // Lower component
		CollapseSpec spec1 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 5.0, t1, 1);
		CollapseSpec spec2 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 5.0, t2, 0);
		mockGetCollapseSpecTime(t1, T_0, spec1);
		mockGetCollapseSpecTime(t2, T_0, spec2);
		when(mockTriangulation.getTriangles()).thenReturn(Arrays.asList(t1, t2));
		eventQueue = new EventQueue(mockTriangulation);
		assertEquals(2L, eventQueue.poll().getEvent().getTriangle().getId()); // t2 first (lower component)
	}

	@Test
	void testOrdering_ByType() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		KineticTriangle t2 = createMockTriangle(2L, 0);
		CollapseSpec spec1 = createCollapseSpec(CollapseType.SPOKE_COLLAPSE, 5.0, t1, 0); // Higher ordinal
		CollapseSpec spec2 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 5.0, t2, 0); // Lower ordinal
		mockGetCollapseSpecTime(t1, T_0, spec1);
		mockGetCollapseSpecTime(t2, T_0, spec2);
		when(mockTriangulation.getTriangles()).thenReturn(Arrays.asList(t1, t2));
		eventQueue = new EventQueue(mockTriangulation);
		assertEquals(2L, eventQueue.poll().getEvent().getTriangle().getId()); // t2 first (lower type ordinal)
	}

	@Test
	void testOrdering_BySecondaryKey() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		KineticTriangle t2 = createMockTriangle(2L, 0);
		CollapseSpec spec1 = createCollapseSpec(CollapseType.VERTEX_MOVES_OVER_SPOKE, 5.0, t1, 0, 10.0);
		CollapseSpec spec2 = createCollapseSpec(CollapseType.VERTEX_MOVES_OVER_SPOKE, 5.0, t2, 0, 20.0); // Higher key -> higher priority
		mockGetCollapseSpecTime(t1, T_0, spec1);
		mockGetCollapseSpecTime(t2, T_0, spec2);
		when(mockTriangulation.getTriangles()).thenReturn(Arrays.asList(t1, t2));
		eventQueue = new EventQueue(mockTriangulation);
		assertEquals(2L, eventQueue.poll().getEvent().getTriangle().getId()); // t2 first
	}

	// --- Update/Drop Tests ---

	@Test
	void testNeedsDropping_RemovesItem() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		KineticTriangle t2 = createMockTriangle(2L, 0);
		CollapseSpec spec1 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 10.0, t1, 0);
		CollapseSpec spec2 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 5.0, t2, 0);
		mockGetCollapseSpecTime(t1, T_0, spec1);
		mockGetCollapseSpecTime(t2, T_0, spec2);
		when(mockTriangulation.getTriangles()).thenReturn(Arrays.asList(t1, t2));
		eventQueue = new EventQueue(mockTriangulation);

		eventQueue.needsDropping(t2); // Request drop for t2

		// No need to mock getCollapseSpec for the processing time here,
		// as t2 should be removed without its spec being recalculated.
		// The default anyDouble() mock for t1 ensures no NPE if EQ checks it.

		eventQueue.processPendingUpdates(T_0);

		assertEquals(1, eventQueue.size());
		assertEquals(1L, eventQueue.peek().getEvent().getTriangle().getId());
		verify(t2).markDead();
	}

	@Test
	void testNeedsUpdate_UpdatesPriority() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		KineticTriangle t2 = createMockTriangle(2L, 0);
		CollapseSpec spec1_initial = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 10.0, t1, 0);
		CollapseSpec spec2_initial = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 5.0, t2, 0); // t2 initially first
		mockGetCollapseSpecTime(t1, T_0, spec1_initial);
		mockGetCollapseSpecTime(t2, T_0, spec2_initial);
		when(mockTriangulation.getTriangles()).thenReturn(Arrays.asList(t1, t2));
		eventQueue = new EventQueue(mockTriangulation);

		assertEquals(2L, eventQueue.peek().getEvent().getTriangle().getId()); // Verify initial state

		// Prepare updated spec for t2 making it later
		CollapseSpec spec2_updated = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 15.0, t2, 0);
		double updateTime = 6.0;

		// *** Mock specific time for the update ***
		mockGetCollapseSpecTime(t2, updateTime, spec2_updated);
		// The default anyDouble() mock handles t1 if it were checked at updateTime

		eventQueue.needsUpdate(t2, false);
		eventQueue.processPendingUpdates(updateTime); // This triggers getCollapseSpec(updateTime) for t2

		assertEquals(2, eventQueue.size());
		assertEquals(1L, eventQueue.peek().getEvent().getTriangle().getId()); // t1 should now be first
	}

	@Test
	void testUpdateAndDrop_DropWins() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		CollapseSpec spec1 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 10.0, t1, 0);
		mockGetCollapseSpecTime(t1, T_0, spec1);
		when(mockTriangulation.getTriangles()).thenReturn(List.of(t1));
		eventQueue = new EventQueue(mockTriangulation);

		eventQueue.needsUpdate(t1, false);
		eventQueue.needsDropping(t1); // Drop requested after update

		eventQueue.processPendingUpdates(T_0);
		assertTrue(eventQueue.isEmpty());
		verify(t1).markDead();
	}

	@Test
	void testDropAndUpdate_DropWins() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		CollapseSpec spec1 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 10.0, t1, 0);
		mockGetCollapseSpecTime(t1, T_0, spec1);
		when(mockTriangulation.getTriangles()).thenReturn(List.of(t1));
		eventQueue = new EventQueue(mockTriangulation);

		eventQueue.needsDropping(t1);
		eventQueue.needsUpdate(t1, false); // Update requested after drop

		eventQueue.processPendingUpdates(T_0);
		assertTrue(eventQueue.isEmpty());
		verify(t1).markDead();
	}

	@Test
	void testProcessPendingUpdates_Idempotent() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		CollapseSpec spec1 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 10.0, t1, 0);
		mockGetCollapseSpecTime(t1, T_0, spec1);
		when(mockTriangulation.getTriangles()).thenReturn(List.of(t1));
		eventQueue = new EventQueue(mockTriangulation);

		eventQueue.processPendingUpdates(T_0);
		assertEquals(1, eventQueue.size());
		eventQueue.processPendingUpdates(T_0); // Call again
		assertEquals(1, eventQueue.size());
	}

	@Test
	void testNeverEvent_HandledCorrectly() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		KineticTriangle t2 = createMockTriangle(2L, 0);
		CollapseSpec spec1 = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 10.0, t1, 0);
		// Create a NEVER spec specifically for t2
		CollapseSpec spec2_never = new CollapseSpec(CollapseType.NEVER, Double.NaN, t2, -1, Double.NaN, 0); // Ensure it has t2

		mockGetCollapseSpecTime(t1, T_0, spec1);
		mockGetCollapseSpecTime(t2, T_0, spec2_never); // Mock t2 returning its NEVER spec

		when(mockTriangulation.getTriangles()).thenReturn(Arrays.asList(t1, t2));
		eventQueue = new EventQueue(mockTriangulation);

		assertEquals(2, eventQueue.size());
		assertEquals(1L, eventQueue.poll().getEvent().getTriangle().getId()); // t1 first

		EventQueueItem item2 = eventQueue.poll(); // Should be t2's item
		assertEquals(2L, item2.getEvent().getTriangle().getId());
		assertEquals(CollapseType.NEVER, item2.getEvent().getType());
		assertTrue(Double.isNaN(item2.getEvent().getTime()));

		assertTrue(eventQueue.isEmpty());
	}

	@Test
	void testUpdate_TriangleBecomesDead() {
		KineticTriangle t1 = createMockTriangle(1L, 0);
		CollapseSpec spec1_initial = createCollapseSpec(CollapseType.TRIANGLE_COLLAPSE, 10.0, t1, 0);
		mockGetCollapseSpecTime(t1, T_0, spec1_initial); // Mock initial state
		when(mockTriangulation.getTriangles()).thenReturn(List.of(t1));
		eventQueue = new EventQueue(mockTriangulation);

		eventQueue.needsUpdate(t1, false); // Request update

		// Simulate t1 becoming dead *before* processing updates
		double updateTime = 5.0;
		lenient().when(t1.isDead()).thenReturn(true); // t1 is now dead
		// Mock getCollapseSpec for the update time (might still be called before isDead
		// check)
		// Return a NEVER spec to be safe if it's called
		mockGetCollapseSpecTime(t1, updateTime, new CollapseSpec(CollapseType.NEVER, Double.NaN, t1, 0));

		eventQueue.processPendingUpdates(updateTime);

		// The item should be removed because the triangle was dead
		assertTrue(eventQueue.isEmpty(), "Queue should be empty as updated triangle was dead");
		verify(t1, never()).markDead(); // EQ should not call markDead, it just observes isDead
	}
}