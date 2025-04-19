package com.github.micycle1.surferj.wavefront;

import static org.junit.jupiter.api.Assertions.*;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import com.github.micycle1.surferj.SurfConstants;
import com.github.micycle1.surferj.collapse.CollapseSpec;
import com.github.micycle1.surferj.collapse.CollapseType;
import com.github.micycle1.surferj.kinetics.KineticEventHandler;
import com.github.micycle1.surferj.kinetics.KineticTriangle;
import com.github.micycle1.surferj.kinetics.KineticTriangulation;

/**
 * Exercise WavefrontPropagator’s stepping logic, injecting a stubbed
 * KineticEventHandler so we can capture the times at which it is invoked.
 */
class WavefrontPropagatorTest {

	/**
	 * Stub triangle whose collapse‐time we can set arbitrarily.
	 */
	private static class FakeTriangle extends KineticTriangle {
		private final long id;
		private double collapseTime;
		private boolean dead = false;

		FakeTriangle(long id, double collapseTime) {
			super(); // adapt if your real KineticTriangle has args
			this.id = id;
			this.collapseTime = collapseTime;
		}

		@Override
		public long getId() {
			return id;
		}

		@Override
		public boolean isDead() {
			return dead;
		}

		public void markDead() {
			dead = true;
		}

		@Override
		public CollapseSpec getCollapseSpec(double now) {
			// ignore 'now'; just return a real CollapseSpec with our fixed time
			return new CollapseSpec(CollapseType.SPLIT_OR_FLIP_REFINE, collapseTime, this);
		}
	}

	/**
	 * A stub KineticTriangulation that simply holds our FakeTriangles. We override
	 * only what the propagator uses.
	 */
	private static class FakeKt extends KineticTriangulation {

		final List<FakeTriangle> tris;

		FakeKt(List<FakeTriangle> tris) {
			super();
			this.tris = tris;
		}

		@Override
		public List<KineticTriangle> getTriangles() {
			// unchecked cast; feeding in our FakeTriangles
			@SuppressWarnings("unchecked")
			List<KineticTriangle> result = (List) tris;
			return result;
		}

		@Override
		public boolean isRestrictComponent() {
			return false;
		}

		@Override
		public void createRemainingSkeletonDcel() {
			// no-op for tests
		}
	}

	/**
	 * A stub event handler which simply records the sequence of times it was asked
	 * to handle.
	 */
	private static class StubEventHandler extends KineticEventHandler {
		final List<Double> handledTimes = new ArrayList<>();

		StubEventHandler(KineticTriangulation kt) {
			super(kt);
		}

		@Override
		public void handleEvent(Event e) {
			// kinetics.Event extends CollapseSpec, so getTime() is available
			handledTimes.add(e.getTime());
		}
	}

	private FakeTriangle t1, t2, t3;
	private FakeKt kt;
	private WavefrontPropagator wp;
	private StubEventHandler stubHandler;

	@BeforeEach
	void setUp() throws Exception {
		// 3 triangles with collapse times 1.0, 2.0, 3.0
		t1 = new FakeTriangle(1, 1.0);
		t2 = new FakeTriangle(2, 2.0);
		t3 = new FakeTriangle(3, 3.0);

		kt = new FakeKt(List.of(t1, t2, t3));

		// Create a dummy SkeletonStructure just to satisfy the constructor:
		SkeletonStructure dummySk = new SkeletonStructure() {
			@Override
			public KineticTriangulation getKt() {
				return kt;
			}
		};

		wp = new WavefrontPropagator(dummySk);
		wp.setupQueue(kt);

		// Now inject our stub handler via reflection
		stubHandler = new StubEventHandler(kt);
		Field fh = WavefrontPropagator.class.getDeclaredField("eventHandler");
		fh.setAccessible(true);
		fh.set(wp, stubHandler);
	}

	@Test
	void testAdvanceStepProcessesOneEventAtATime() {
		assertEquals(0.0, wp.getTime(), 1e-9);
		assertTrue(stubHandler.handledTimes.isEmpty());

		wp.advanceStep();
		assertEquals(1.0, wp.getTime(), 1e-9);
		assertEquals(List.of(1.0), stubHandler.handledTimes);

		wp.advanceStep();
		assertEquals(2.0, wp.getTime(), 1e-9);
		assertEquals(List.of(1.0, 2.0), stubHandler.handledTimes);

		wp.advanceStep();
		assertEquals(3.0, wp.getTime(), 1e-9);
		assertEquals(List.of(1.0, 2.0, 3.0), stubHandler.handledTimes);
	}

	@Test
	void testAdvanceTimeNextJumpsClockButDoesNotFire() {
		wp.advanceTimeNext();
		assertEquals(1.0, wp.getTime(), 1e-9);
		assertTrue(stubHandler.handledTimes.isEmpty());
	}

	@Test
	void testAdvanceTimeIgnoreEventWithTarget() {
		wp.advanceTimeIgnoreEvent(2.1);
		assertEquals(2.1, wp.getTime(), 1e-9);
		// events @1.0 and @2.0 should have fired
		assertEquals(List.of(1.0, 2.0), stubHandler.handledTimes);

		// one more step should process the last event
		wp.advanceStep();
		assertEquals(3.0, wp.getTime(), 1e-9);
		assertEquals(List.of(1.0, 2.0, 3.0), stubHandler.handledTimes);
	}

	@Test
	void testAdvanceTimeIgnoreEventNoArgsJustBumpsTime() {
		final double before = wp.getTime();
		final double inc = wp.increment;
		wp.advanceTimeIgnoreEvent();
		assertEquals(before + inc, wp.getTime(), 1e-9);
		// no events should have fired
		assertTrue(stubHandler.handledTimes.isEmpty());
	}

	@Test
	void testAdvanceTimeDrainsEventsWhenPresent() {
		wp.setIncrement(2.5);
		wp.advanceTime(); // from 0 → 2.5
		assertEquals(2.5, wp.getTime(), 1e-9);
		assertEquals(List.of(1.0, 2.0), stubHandler.handledTimes);
	}

	@Test
	void testAdvanceTimeSingleStepIfNoMoreEvents() {
		// drain manually
		wp.advanceStep();
		wp.advanceStep();
		wp.advanceStep();
		assertEquals(List.of(1.0, 2.0, 3.0), stubHandler.handledTimes);

		// now queue empty: advanceTime should bump & do exactly one advanceStep()
		wp.setIncrement(0.7);
		double before = wp.getTime();
		wp.advanceTime();
		assertEquals(before + 0.7, wp.getTime(), 1e-9);
		// no new events
		assertEquals(3, stubHandler.handledTimes.size());
	}

	@Test
	void testAdvanceToEndFinalizes() {
		assertFalse(wp.isPropagationComplete());
		wp.advanceToEnd();
		assertTrue(wp.isPropagationComplete());

		// all three must have fired
		assertEquals(List.of(1.0, 2.0, 3.0), stubHandler.handledTimes);

		// the KineticTriangulation should have built the DCEL
//    assertTrue(kt.createdDcel);
	}

	@Test
	void testDoInitialSkipsSkipAll() {
		wp.doInitialSkips(true, 0, SurfConstants.ZERO);
		assertTrue(wp.isPropagationComplete());
		assertEquals(List.of(1.0, 2.0, 3.0), stubHandler.handledTimes);
	}

	@Test
	void testDoInitialSkipsByCount() {
		// skip_to = 2 → C++ processes exactly 1 event
		wp.doInitialSkips(false, 2, SurfConstants.ZERO);

		assertFalse(wp.isPropagationComplete());
		// only event#1 at t=1.0 should have fired
		assertEquals(List.of(1.0), stubHandler.handledTimes);

		// next event in the queue is event#2 at t=2.0
		assertEquals(2.0, peekNextTime(), 1e-9);
	}

	@Test
	void testDoInitialSkipsByTime() {
		wp.doInitialSkips(false, 0, 2.1);
		assertFalse(wp.isPropagationComplete());
		assertEquals(List.of(1.0, 2.0), stubHandler.handledTimes);
		assertEquals(2.1, wp.getTime(), 1e-9);
		assertEquals(3.0, peekNextTime(), 1e-9);
	}

	@Test
	void testResetTimeToLastEvent() {
		wp.advanceStep(); // 1.0
		wp.advanceStep(); // 2.0
		assertEquals(2.0, wp.getTime(), 1e-9);

		// roll back to that last event
		wp.resetTimeToLastEvent();
		assertEquals(2.0, wp.getTime(), 1e-9);

		// then step again → 3.0
		wp.advanceStep();
		assertEquals(3.0, wp.getTime(), 1e-9);
	}

	/** Helper to inspect the next‐event time in the queue. */
	private double peekNextTime() {
		WavefrontPropagator w = wp;
		EventQueueItem next = w.peekNextEvent();
		return next == null ? Double.NaN : next.getEvent().getTime();
	}
}