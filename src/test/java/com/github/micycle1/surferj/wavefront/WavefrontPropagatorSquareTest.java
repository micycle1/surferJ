package com.github.micycle1.surferj.wavefront;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;

import com.github.micycle1.surferj.collapse.CollapseType;
import com.github.micycle1.surferj.kinetics.KineticTriangle;
import com.github.micycle1.surferj.kinetics.KineticTriangulation;
// Import other necessary classes (Vec2d, Coordinate, etc.)

// Assume SkeletonDCEL is sufficiently stubbed for initialization
// Assume BasicTriangulation and initialization work correctly

@Disabled
public class WavefrontPropagatorSquareTest {

	private SkeletonStructure sk;
	private WavefrontPropagator wp;
	private KineticTriangulation kt;
	// private EventQueue eq; // Access queue via wp if needed

	// Define test tolerance
	private static final double TIME_DELTA = 1e-7;
	private static final double POS_DELTA = 1e-7;

	// JTS GeometryFactory
	private static final GeometryFactory geometryFactory = new GeometryFactory();

	// Helper to create a square JTS Polygon centered at origin, size 2x2
	private Polygon createSquareJtsPolygon() {
		Coordinate[] coords = new Coordinate[] { new Coordinate(-1, -1), new Coordinate(1, -1), new Coordinate(1, 1), new Coordinate(-1, 1),
				new Coordinate(-1, -1) // Close the ring
		};
		LinearRing shell = geometryFactory.createLinearRing(coords);
		return geometryFactory.createPolygon(shell);
	}

	@BeforeEach
	void setUp() {
		Polygon squarePolygon = createSquareJtsPolygon();

		// *** ADAPTATION POINT ***
		// Assume SkeletonStructure constructor now takes Polygon,
		// or that there's a clear way to initialize it from a Polygon.
		// Option 1: SkeletonStructure(Polygon) constructor exists
		sk = new SkeletonStructure(squarePolygon);

		// Option 2: Initialize method takes Polygon (less likely based on C++)
		// sk = new SkeletonStructure(...); // May need default constructor or different
		// args
		// sk.initializeFromPolygon(squarePolygon, -1); // Example method

		// Option 3: Convert Polygon to intermediate format if SkeletonStructure MUST
		// take BasicInput
		// BasicInput squareInput = convertPolygonToBasicInput(squarePolygon); // Need
		// this conversion
		// sk = new SkeletonStructure(squareInput);

		// The rest of the initialization happens inside initialize()
		sk.initialise(-1); // Process all components

		wp = sk.getWp();
		kt = sk.getKt();
		// eq = wp.getEqForTesting(); // Get queue reference if needed
	}

	// Helper method to get EQ (add to WavefrontPropagator or SkeletonStructure for
	// testing)
	// Example: add this to WavefrontPropagator
	// public EventQueue getEqForTesting() { return eq; }

	@Test
	void testSquareInitializationState() {
		assertNotNull(wp, "WavefrontPropagator should be initialized");
		assertNotNull(kt, "KineticTriangulation should be initialized");

		assertEquals(0.0, wp.getTime(), "Initial time should be 0.0");
		assertFalse(wp.isPropagationComplete(), "Simulation should not be complete initially");

		// --- Compare with C++ Log (using equivalent JTS Polygon input) ---
		EventQueueItem initialTopItem = wp.peekNextEvent();
		assertNotNull(initialTopItem, "Initial event queue should not be empty");

		Event initialEvent = initialTopItem.getEvent();
		// *** UPDATE EXPECTED VALUES BASED ON C++ LOG WITH JTS-EQUIVALENT INPUT ***
		assertEquals(CollapseType.TRIANGLE_COLLAPSE, initialEvent.getType(), "Initial event type mismatch");
		assertEquals(1.0, initialEvent.getTime(), TIME_DELTA, "Initial event time mismatch"); // Example value
		// Check triangle ID (might change based on triangulation from JTS)
		assertTrue(initialEvent.getTriangle().getId() >= 0, "Initial event triangle ID should be valid");
		// assertEquals(EXPECTED_ID_FROM_LOG, initialEvent.getTriangle().getId(),
		// "Initial event triangle ID mismatch");

		// Check number of initial triangles (should still be 2 for simple square)
		// long activeCount = kt.getTriangles().stream().filter(t -> t != null &&
		// !t.isDead()).count();
		// assertEquals(2, activeCount, "Should have 2 initial active triangles");
	}

	@Test
	void testSquareSingleStep() {
		// --- Get Initial State (Update expected values from C++ log) ---
		EventQueueItem initialTopItem = wp.peekNextEvent();
		assertNotNull(initialTopItem);
		long expectedInitialTriangleId = initialTopItem.getEvent().getTriangle().getId(); // Use actual ID
		double expectedEventTime = 1.0; // Example value, update from log
		assertEquals(expectedEventTime, initialTopItem.getEvent().getTime(), TIME_DELTA);

		// --- Execute Step 1 ---
		wp.advanceStep();

		// --- Verify State After Step 1 (Update expected values from C++ log) ---
		assertFalse(wp.isPropagationComplete(), "Simulation should not be complete after 1 step");
		assertEquals(expectedEventTime, wp.getTime(), TIME_DELTA, "Time should advance to event time");
		assertEquals(1, wp.getEventCtrForTesting(), "Event counter should be 1"); // Need getter

		// Check triangulation state changes
		KineticTriangle processedTriangle = kt.getTriangles().get((int) expectedInitialTriangleId); // Assuming ID is index
		assertTrue(processedTriangle.isDying() || processedTriangle.isDead(), "Processed triangle should be marked dying/dead");

		// Check queue state after updates
		EventQueueItem nextTopItem = wp.peekNextEvent();
		assertNotNull(nextTopItem, "Queue should have next event after step 1");
		Event nextEvent = nextTopItem.getEvent();
		// Find the *other* initial triangle ID (assuming 2 triangles)
		long expectedNextTriangleId = kt.getTriangles().stream().filter(t -> t != null && t.getId() != expectedInitialTriangleId).map(KineticTriangle::getId)
				.findFirst().orElse(-1L);
		assertTrue(expectedNextTriangleId != -1, "Could not find the other initial triangle");

		assertEquals(CollapseType.TRIANGLE_COLLAPSE, nextEvent.getType(), "Next event type mismatch");
		assertEquals(expectedEventTime, nextEvent.getTime(), TIME_DELTA, "Next event time mismatch");
		assertEquals(expectedNextTriangleId, nextEvent.getTriangle().getId(), "Next event triangle ID mismatch");
	}

	@Test
	void testSquareRunToEnd() {
		// --- Execute to End ---
		wp.advanceToEnd();

		// --- Verify Final State (Update expected values from C++ log) ---
		assertTrue(wp.isPropagationComplete(), "Simulation should be complete");
		assertEquals(1.0, wp.getTime(), TIME_DELTA, "Final time should be the last event time"); // Example value
		assertEquals(2, wp.getEventCtrForTesting(), "Final event count mismatch"); // Example value

		EventQueueItem finalTopItem = wp.peekNextEvent();
		assertNull(finalTopItem, "Queue should be empty at the end");
	}

	// Add more tests...

	// Helper method to get event counter (add to WavefrontPropagator for testing)
	// Example: add this to WavefrontPropagator
	// public long getEventCtrForTesting() { return eventCtr; }

	// Getter for EventQueue (if needed)
	// Example: add this to WavefrontPropagator
	// public EventQueue getEqForTesting() { return eq; }
}