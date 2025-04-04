package com.github.micycle.surferj;

import java.util.List;

import org.locationtech.jts.geom.Polygon;

import com.github.micycle.surferj.event.CollapseEvent;
import com.github.micycle.surferj.event.EventQueue;
import com.github.micycle.surferj.kinetics.KineticTriangulation;
import com.github.micycle.surferj.triangulation.InitialTriangulator;
import com.github.micycle.surferj.util.Constants;

public class StraightSkeletonGenerator {

	private final Polygon inputPolygon;
	private InitialTriangulator initialTriangulator;
	private KineticTriangulation kineticTriangulation;
	private EventQueue eventQueue;
	private double currentTime = 0.0;
	private int eventCounter = 0;
	private final int maxEvents; // Safety break

	public StraightSkeletonGenerator(Polygon polygon, int maxEvents) {
		this.inputPolygon = polygon;
		this.maxEvents = maxEvents;
	}

	public SkeletonOutput generate() {
		System.out.println("1. Initial Triangulation...");
		initialTriangulator = new InitialTriangulator(inputPolygon);
		initialTriangulator.triangulate();
		if (initialTriangulator.getTin() == null || initialTriangulator.getTriangles().isEmpty()) {
			System.err.println("Initial triangulation failed.");
			return new SkeletonOutput(List.of());
		}
		System.out.println("   Triangulation complete. Triangles: " + initialTriangulator.getTriangles().size());

		System.out.println("2. Initializing Kinetic Structure...");
		kineticTriangulation = new KineticTriangulation(initialTriangulator);
		eventQueue = new EventQueue();
		kineticTriangulation.initialize(eventQueue);
		System.out.println("   Kinetic structure initialized.");

		System.out.println("3. Running Simulation...");
		while (!eventQueue.isEmpty() && eventCounter < maxEvents) {
			CollapseEvent nextEvent = eventQueue.peek(); // Don't remove yet

			if (nextEvent.getTime() < currentTime - Constants.EPSILON) {
				System.err.println("Error: Event time regression! Current: " + currentTime + ", Event: " + nextEvent);
				// This indicates a logic error or precision issue. Stop processing.
				break;
			}

			// Check for very small time steps - potential precision issue / infinite loop
			if (Math.abs(nextEvent.getTime() - currentTime) < Constants.EPSILON && eventCounter > 0) {
				System.out.println("Warning: Very small time step or duplicate event time at " + currentTime);
				// Might need more sophisticated tie-breaking or tolerance handling.
			}

			currentTime = nextEvent.getTime();
			eventCounter++;

			// Actually remove the event now
			nextEvent = eventQueue.poll();
			if (nextEvent == null)
				break; // Should not happen if isEmpty check passed

			System.out.println("\n--- Event " + eventCounter + " at time " + String.format("%.5f", currentTime) + " ---");
			kineticTriangulation.handleEvent(nextEvent);

			// Basic validation after event (optional, can be slow)
			// validateStructure();
		}

		if (eventCounter >= maxEvents) {
			System.err.println("Warning: Reached maximum event limit (" + maxEvents + "). Skeleton may be incomplete.");
		}
		System.out.println("   Simulation finished after " + eventCounter + " events. Time: " + String.format("%.5f", currentTime));

		System.out.println("4. Generating Output Skeleton...");
		SkeletonOutput output = kineticTriangulation.generateSkeletonOutput();
		System.out.println("   Output generated.");

		return output;
	}

	// Basic validation (placeholder)
	private void validateStructure() {
		// Check for overlapping triangles, consistent neighbor pointers etc.
		// This would be quite involved.
	}

}
