package com.github.micycle1.surferj.wavefront;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.github.micycle1.surferj.SurfConstants;
import com.github.micycle1.surferj.collapse.CollapseType;
import com.github.micycle1.surferj.kinetics.KineticEventHandler;
import com.github.micycle1.surferj.kinetics.KineticTriangulation;

/**
 * Manages the propagation of the wavefront by processing events from an
 * EventQueue. It advances simulation time, handles events by delegating to
 * KineticTriangulation, and manages the overall simulation lifecycle.
 */

public class WavefrontPropagator {

	private static final Logger LOGGER = LoggerFactory.getLogger(WavefrontPropagator.class);

	private final SkeletonStructure sk; // Reference back to the main structure
	private EventQueue eq; // The event queue, set during initialization
	private long eventCtr = 0; // Event counter
	private boolean finalized = false; // Has the simulation finished and finalized?

	private double time = SurfConstants.ZERO; // Current simulation time
	private double lastEventTime = SurfConstants.ZERO; // Time of the last processed event
	private double increment = 0.0005; // Default time step for manual advance

	// Component tracking (optional, for visualization/debugging)
	private int currentComponent = -1;
	private int lastEventComponent = -1;

	private KineticEventHandler eventHandler; // NOTE new in Java version

	/**
	 * Creates a WavefrontPropagator associated with a SkeletonStructure.
	 *
	 * @param sk The parent SkeletonStructure.
	 */
	public WavefrontPropagator(SkeletonStructure sk) {
		this.sk = sk;
		eventHandler = new KineticEventHandler(sk.getKt());
	}

	/**
	 * Sets up the event queue using the provided kinetic triangulation. This is
	 * typically called after the triangulation is initialized.
	 *
	 * @param kt The initialized kinetic triangulation.
	 */
	void setupQueue(KineticTriangulation kt) {
		this.eq = new EventQueue(kt);
		kt.setQueue(this.eq); // Provide the queue back to the triangulation
	}

	/**
	 * Checks if the simulation has been finalized.
	 *
	 * @return true if finalize() has been called.
	 */
	public boolean isPropagationComplete() {
		return finalized;
	}

	/**
	 * Gets the current simulation time.
	 *
	 * @return The current time.
	 */
	public double getTime() {
		return time;
	}

	/**
	 * Sets the time increment used by advanceTime().
	 *
	 * @param i The new time increment.
	 */
	public void setIncrement(double i) {
		this.increment = i;
	}

	/**
	 * Gets the current component being processed (for visualization/debugging). -1
	 * means all components or component tracking is inactive.
	 *
	 * @return The current component index.
	 */
	public int getCurrentComponent() {
		return currentComponent;
	}

	/**
	 * Checks if there are no more processable events in the queue. An empty queue
	 * or a queue where the next event is 'NEVER' means no more events. Ensures
	 * pending updates are processed before checking.
	 *
	 * @return true if no more events are left.
	 */
	private boolean noMoreEvents() {
		if (eq == null) {
			return true; // Not initialized yet
		}
		eq.processPendingUpdates(this.time); // Process updates before checking state
		EventQueueItem next = eq.peek();
		return next == null || next.getEvent().getType() == CollapseType.NEVER;
	}

	/**
	 * Advances simulation time by the current increment. If events occur during
	 * this time step, they are processed sequentially until the target time is
	 * reached.
	 */
	public void advanceTime() {
		double targetTime = this.time + this.increment;
		advanceTimeIgnoreEvent(targetTime); // Process events up to targetTime
	}

	/**
	 * Advances simulation time to the specified target time. Processes all events
	 * that occur before or at this target time. Use this instead of setting time
	 * directly if events need processing.
	 *
	 * @param targetTime The time to advance to.
	 */
	public void advanceTimeIgnoreEvent(double targetTime) {
		if (finalized || eq == null) {
			this.time = targetTime; // Just set time if finalized or not setup
			return;
		}

		while (!isPropagationComplete()) {
			eq.processPendingUpdates(this.time); // Ensure queue is up-to-date
			EventQueueItem nextItem = eq.peek();

			if (nextItem == null || nextItem.getEvent().getTime() > targetTime) {
				// No more events before target time
				break;
			}
			// Next event is at or before target time, process it
			advanceStep();
		}
		// After processing all events <= targetTime, set time to the target
		this.time = targetTime;
	}

	/**
	 * Moves simulation time forward to the time of the very next event, without
	 * processing the event itself.
	 */
	public void advanceTimeNext() {
		if (finalized || eq == null) {
			return;
		}
		eq.processPendingUpdates(this.time); // Ensure queue is up-to-date
		EventQueueItem next = eq.peek();
		if (next != null && next.getEvent().getType() != CollapseType.NEVER) {
			this.time = next.getEvent().getTime();
			if (sk.getKt().isRestrictComponent()) {
				// Component tracking might be based on the event's triangle
				this.currentComponent = next.getEvent().getTriangle().getComponent();
			}
		}
		// If next is null or NEVER, time doesn't advance
	}

	/**
	 * Processes the next single event from the queue. Advances simulation time to
	 * the event time and delegates event handling to the KineticTriangulation.
	 */
	public void advanceStep() {
		// C++ DBG_INDENT_LEVEL_STORE / CHECK omitted

		if (finalized || eq == null) {
			return;
		}

		// CRITICAL: Process pending updates *before* polling the next event
		eq.processPendingUpdates(this.time);

		if (noMoreEvents()) {
			LOGGER.info("No more events to process.");
			finalizeSim(); // Automatically finalize if queue becomes empty
			return;
		}

		EventQueueItem nextItem = eq.poll(); // Poll retrieves and removes the head
		if (nextItem == null) { // Should not happen if noMoreEvents was false, but safety check
			LOGGER.warn("Polled null item unexpectedly.");
			finalizeSim();
			return;
		}

		Event event = nextItem.getEvent();

		// Advance time to the event time
		this.time = event.getTime();
		if (sk.getKt().isRestrictComponent()) {
			this.currentComponent = event.getTriangle().getComponent();
		}

		this.eventCtr++;
		LOGGER.info("Processing Event #{} @ t={}: {}", eventCtr, String.format("%.6f", time), event);

		eventHandler.handleEvent(event);

		// CRITICAL: Process updates generated *by handling this event*.
		// This ensures the queue is correct before the *next* advanceStep or check.
		eq.processPendingUpdates(this.time);

		LOGGER.info("Finished Event #{}. Queue size: {}", eventCtr, eq.size());

		// Update last event info
		this.lastEventTime = this.time;
		this.lastEventComponent = this.currentComponent;

		// Check again if simulation should end
		if (noMoreEvents()) {
			LOGGER.info("All events processed after event #{}", eventCtr);
			finalizeSim();
		} else {
			// Log info about the *next* event
			EventQueueItem nextEvent = eq.peek(); // Peek after updates
			if (nextEvent != null) {
				LOGGER.debug("    Next event: {} @ t={}", nextEvent.getEvent(), String.format("%.6f", nextEvent.getEvent().getTime()));
			}
		}
	}

	/**
	 * Processes all remaining events in the queue until the simulation is complete.
	 */
	public void advanceToEnd() {
		LOGGER.info("Advancing simulation to the end...");
		while (!isPropagationComplete()) {
			if (noMoreEvents()) { // Check required in case queue empties mid-process
				finalizeSim();
				break;
			}
			advanceStep();
		}
		LOGGER.info("Simulation finished. Total events processed: {}", eventCtr);
	}

	/**
	 * Skips initial events based on count or time.
	 *
	 * @param skipAll       If true, advances to the end immediately.
	 * @param skipToEvent   Process events until the event counter reaches this
	 *                      number (1-based).
	 * @param skipUntilTime Process events until the simulation time reaches this
	 *                      value.
	 */
	public void doInitialSkips(boolean skipAll, long skipToEvent, double skipUntilTime) {
		if (skipAll) {
			advanceToEnd();
			return;
		}

		// Skip by event count (process N events)
		// Note: eventCtr is 0-based internally but increments before processing event
		// #1
		while (!isPropagationComplete() && skipToEvent > eventCtr + 1) {
			if (noMoreEvents()) {
				break;
			}
			advanceStep();
		}

		// Skip by time
		if (skipUntilTime > SurfConstants.ZERO && !isPropagationComplete()) {
			// Process events strictly *before* the skip time
			while (!isPropagationComplete()) {
				eq.processPendingUpdates(this.time); // Update before peeking
				EventQueueItem next = eq.peek();
				if (next == null || next.getEvent().getTime() >= skipUntilTime) {
					break; // Next event is at or after skip time
				}
				advanceStep();
			}
			// Now, ensure current time is at least the skip time
			if (skipUntilTime > getTime()) {
				// Use advanceTimeIgnoreEvent to set time correctly without processing
				// events *at* skipUntilTime if they weren't already processed.
				advanceTimeIgnoreEvent(skipUntilTime);
			}
		}
	}

	/**
	 * Finalizes the simulation. This typically involves creating the remaining
	 * skeleton structures in the DCEL based on the final state of the
	 * triangulation. Prevents further event processing.
	 */
	public void finalizeSim() {
		if (!finalized) {
			LOGGER.info("Finalizing simulation...");
			// Ask KineticTriangulation to build the final DCEL parts
			sk.getKt().createRemainingSkeletonDcel();
			finalized = true;
			currentComponent = -1; // Reset component tracking

			// Update stats if that system is ported
			// sk.getKt().updateEventTimingStats(-1);
			LOGGER.info("Simulation finalized.");
		}
	}

	/**
	 * Resets the simulation time and current component back to the state right
	 * after the last processed event. Useful for stepping backwards in a debugger
	 * or visualization.
	 */
	public void resetTimeToLastEvent() {
		this.time = this.lastEventTime;
		this.currentComponent = this.lastEventComponent;
		// NOTE: This does NOT undo the state changes caused by the last event.
		// It only resets the time/component tracker. True rollback is complex.
	}

	/**
	 * Returns the next event item in the queue without removing it. Processes
	 * pending updates first.
	 *
	 * @return The next EventQueueItem or null if none.
	 */
	public EventQueueItem peekNextItem() {
		if (eq == null) {
			return null;
		}
		eq.processPendingUpdates(this.time);
		return eq.peek();
	}

	long getEventCtrForTesting() {
		return eventCtr;
	}
}