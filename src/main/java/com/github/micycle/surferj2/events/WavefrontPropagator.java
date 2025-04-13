package com.github.micycle.surferj2.events;

import org.apache.commons.lang3.tuple.Pair;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;

import com.github.micycle.surferj2.SurfConstants;
import com.github.micycle.surferj2.collapse.CollapseSpec;
import com.github.micycle.surferj2.collapse.CollapseType;
import com.github.micycle.surferj2.kinetics.KineticTriangle;
import com.github.micycle.surferj2.kinetics.KineticTriangulation;
import com.github.micycle.surferj2.kinetics.WavefrontEdge;
import com.github.micycle.surferj2.kinetics.WavefrontVertex;

import java.util.*;

/**
 * Orchestrates the straight skeleton wavefront propagation algorithm. It
 * processes events from an EventQueue in chronological order, updating the
 * KineticTriangulation at each step.
 *
 * Corresponds mostly to C++ WavefrontPropagator and event handling logic
 * previously within KineticTriangulation.cpp.
 */
public class WavefrontPropagator {

	private static final double TIME_COMPARISON_TOLERANCE = 1e-9; // Tolerance for time comparisons
	private static final int MAX_EVENT_COUNT = 500000; // Safety break

	private final KineticTriangulation kineticTriangulation;
	private final EventQueue eventQueue;
	private double currentTime;
	private double lastEventTime; // For statistics/debugging C++ porting
	private int currentComponent = -1; // For statistics/debugging C++ porting
	private int lastEventComponent = -1; // For statistics/debugging C++ porting
	private long eventCounter; // For logging/debugging
	private boolean finalized; // To prevent running multiple times

	// Statistics mirroring C++
	private final long[] eventTypeCounter = new long[CollapseType.values().length]; // Use Enum ordinals
	private int maxTrianglesPerEdgeEvent = 0;
	private long avgTrianglesPerEdgeEventSum = 0;
	private long avgTrianglesPerEdgeEventCtr = 0;
	private int maxTrianglesPerSplitEvent = 0;
	private long avgTrianglesPerSplitEventSum = 0;
	private long avgTrianglesPerSplitEventCtr = 0;
	private int eventsPerCurrentEventTime = 0;
	private int maxEventsPerTime = 0;
	private long avgEventsPerTimeSum = 0;
	private long avgEventsPerTimeCtr = 0;

	public WavefrontPropagator(KineticTriangulation kineticTriangulation) {
		this.kineticTriangulation = Objects.requireNonNull(kineticTriangulation, "KineticTriangulation cannot be null");
		// Initialize EventQueue potentially requires the initial set of triangles
		this.eventQueue = new EventQueue(kineticTriangulation.getTriangles(), 0.0);
		this.kineticTriangulation.setEventQueue(this.eventQueue); // Link back if needed by KT methods
		this.currentTime = 0.0;
		this.lastEventTime = 0.0;
		this.eventCounter = 0;
		this.finalized = false;
	}

	/**
	 * Runs the simulation until no more events are left in the queue.
	 */
	public void propagateToEnd() {
		if (finalized) {
			System.err.println("Warning: Propagation already finalized. Cannot run again.");
			return;
		}
		System.out.println("Starting wavefront propagation...");

		while (eventQueue.hasMoreEvents()) {
			kineticTriangulation.clearNewElementsList(); // Clear tracking lists at start of step

			// Get the absolute next event based on time, handling potential invalid ones
			CollapseSpec nextEvent = eventQueue.pollNextValidEvent();

			if (nextEvent == null) {
				System.out.println("No more valid events polled. Propagation finished.");
				break; // Simulation ends
			}

			// Validate time progression
			if (nextEvent.getTime() < currentTime - TIME_COMPARISON_TOLERANCE) {
				System.err.printf("FATAL ERROR: Event time %.12f is before current time %.12f. Aborting.%nEvent: %s%n", nextEvent.getTime(), currentTime,
						nextEvent);
				throw new IllegalStateException("Event time regression detected.");
			}

			// Advance time
			currentTime = nextEvent.getTime();
			updateEventTimingStats(currentTime, nextEvent.getTriangle().getComponent()); // C++ stats port
			eventCounter++;

			System.out.printf("--- Processing Event #%d @ time %.9f (Comp %d): %s ---%n", eventCounter, currentTime, nextEvent.getTriangle().getComponent(),
					nextEvent);

			// Process the event - This modifies the triangulation AND updates the queue
			handleEvent(nextEvent);

			// Update events for newly created elements (mirrors C++ post-event update
			// logic)
			updateEventsForNewElements();

			if (eventCounter > MAX_EVENT_COUNT) { // Safety break
				System.err.println("Error: Exceeded maximum event count (" + MAX_EVENT_COUNT + "). Stopping simulation.");
				throw new RuntimeException("Exceeded maximum event count");
				// break;
			}
		}

		finalizePropagation();
		System.out.printf("Propagation finalized at time: %.9f after %d events.%n", currentTime, eventCounter);
		printStatistics();
	}

	/** Dispatcher mirroring C++ KineticTriangulation::handle_event */
	private void handleEvent(CollapseSpec event) {
		KineticTriangle triangle = event.getTriangle();
		if (triangle == null) { // Should not happen if pollNextValidEvent works
			System.err.println("Critical Error: Event has null triangle: " + event);
			return;
		}
		if (triangle.isDead()) { // Should not happen if pollNextValidEvent works
			System.err.println("Warning: Skipping event for already dead triangle: " + event);
			return;
		}
		// Check if event is still valid (matches current computed event for triangle)
		// This handles cases where a previous event invalidated this one
		CollapseSpec currentSpec = triangle.getCollapseSpec(currentTime); // Recompute or get cached
		if (!event.equals(currentSpec)) {
			System.out.println("  Skipping stale event: " + event + " (Current is: " + currentSpec + ")");
			// No need to invalidate, pollNextValidEvent handles the map cleanup
			return;
		}

		// Event is valid and for a live triangle, proceed.
		// Event queue invalidation is handled by pollNextValidEvent/addOrUpdate

		// C++ Stats port
		eventTypeCounter[CollapseType.UNDEFINED.ordinal()]++; // Overall counter
		eventTypeCounter[event.getType().ordinal()]++;

		switch (event.getType()) {
			case TRIANGLE_COLLAPSE :
				handleTriangleCollapse(event);
				break;
			case CONSTRAINT_COLLAPSE :
				handleConstraintCollapse(event);
				break;
			case SPOKE_COLLAPSE :
				handleSpokeCollapse(event);
				break;
			case SPLIT_OR_FLIP_REFINE :
				handleSplitOrFlipRefine(event);
				break;
			case VERTEX_MOVES_OVER_SPOKE :
				handleVertexMovesOverSpoke(event);
				break;
			case CCW_VERTEX_LEAVES_CH :
				handleCcwVertexLeavesCH(event);
				break;
			case FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING :
				handleFaceWithInfinitelyFastOpposingVertex(event);
				break;
			case FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED :
				handleFaceWithInfinitelyFastWeightedVertex(event);
				break;
			case INVALID_EVENT :
				System.err.println("FATAL ERROR: Handling INVALID_EVENT. Event: " + event);
				throw new IllegalStateException("Invalid event encountered: " + event);
			case NEVER :
			case UNDEFINED :
			default :
				System.err.println("FATAL ERROR: Unexpected event type in handleEvent: " + event.getType() + ". Event: " + event);
				throw new IllegalStateException("Unexpected event type: " + event.getType());
		}
	}

	// --- Event Handler Implementations ---

	/** Mirrors C++ KineticTriangulation::handle_triangle_collapse_event */
	private void handleTriangleCollapse(CollapseSpec event) {
		System.out.println("  Handling TRIANGLE_COLLAPSE for KT" + event.getTriangle().id);
		KineticTriangle tri = event.getTriangle();
		double time = event.getTime();

		// 1. Stop Vertices
		WavefrontVertex v0 = tri.getVertex(0);
		WavefrontVertex v1 = tri.getVertex(1);
		WavefrontVertex v2 = tri.getVertex(2);

		// Need a reliable stop position. Calculate from one vertex.
		Coordinate stopPos = v0.getPositionAt(time); // Assume getPositionAt handles stopped vertices correctly
		if (v0.hasStopped())
			stopPos = v0.getPosStop();

		System.out.printf("    Stopping vertices %d, %d, %d at (%.6f, %.6f)%n", v0.id, v1.id, v2.id, stopPos.x, stopPos.y);
		v0.stop(time, stopPos);
		v1.stop(time, stopPos);
		v2.stop(time, stopPos);

		// 2. Invalidate neighbors' events & update their structure
		List<KineticTriangle> neighborsToUpdate = new ArrayList<>(3);
		for (int i = 0; i < 3; i++) {
			WavefrontEdge wfe = tri.getWavefront(i);
			KineticTriangle neighbor = tri.getNeighbor(i);

			if (wfe != null) { // Constraint Edge
				System.out.println("    Edge " + i + ": Constraint WFE" + wfe.id);
				if (!wfe.isDead()) {
					wfe.markDead();
					// DCEL Link: Head-to-Head across the dead edge
					WavefrontVertex edgeV_ccw = tri.getVertex((i + 1) % 3); // CCW vertex is edge's v0 perspective in Tri
					WavefrontVertex edgeV_cw = tri.getVertex((i + 2) % 3); // CW vertex is edge's v1 perspective in Tri
					// Link edgeV_cw side 0 (the side facing the dead edge) to edgeV_ccw side 1 (the
					// side facing the dead edge)
					System.out.printf("      DCEL Link H2H: V%d (side 0) <-> V%d (side 1)%n", edgeV_cw.id, edgeV_ccw.id);
					edgeV_cw.linkHeadToHead(0, edgeV_ccw, 1);
				} else {
					System.out.println("      Constraint WFE" + wfe.id + " already dead.");
				}
			} else if (neighbor != null) { // Internal Edge (Spoke)
				System.out.println("    Edge " + i + ": Spoke to KT" + neighbor.id);
				if (!neighbor.isDead()) {
					int neighborEdgeIndex = neighbor.indexOf(tri);
					if (neighborEdgeIndex == -1) {
						throw new IllegalStateException("Neighbor " + neighbor.id + " doesn't know triangle " + tri.id);
					}
					System.out.println("      Removing neighbor link from KT" + neighbor.id + " edge " + neighborEdgeIndex);
					neighbor.setNeighbor(neighborEdgeIndex, null); // Remove link to dead triangle
					// This neighbor now effectively experiences a spoke collapse at its boundary.
					// Call the helper, which might merge vertices or kill the neighbor too.
					System.out.println("      Triggering spoke collapse part 2 for neighbor KT" + neighbor.id + " edge " + neighborEdgeIndex);
					handleSpokeCollapsePart2(neighbor, neighborEdgeIndex, time);
					// Don't add neighbor to update list here; handleSpokeCollapsePart2 will trigger
					// updates if needed.
				} else {
					System.out.println("      Neighbor KT" + neighbor.id + " already dead.");
				}
			} else {
				System.err.println("Warning: Triangle " + tri.id + " edge " + i + " has neither wavefront nor neighbor!");
				// This might happen if the neighbor died in the same step?
			}
		}

		// 3. Mark triangle dead
		System.out.println("    Marking KT" + tri.id + " dead.");
		tri.markDead();
		// Event queue update: Invalidation implicitly done by polling.
		// Recalculations triggered by handleSpokeCollapsePart2 for neighbors.
	}

	/** Mirrors C++ KineticTriangulation::handle_constraint_event */
	private void handleConstraintCollapse(CollapseSpec event) {
		System.out.println("  Handling CONSTRAINT_COLLAPSE for KT" + event.getTriangle().id + " edge " + event.getRelevantEdge());
		KineticTriangle tri = event.getTriangle();
		int edgeIndex = event.getRelevantEdge();
		double time = event.getTime();

		WavefrontVertex va = tri.getVertex((edgeIndex + 1) % 3); // CCW vertex of edge (Edge's v0)
		WavefrontVertex vb = tri.getVertex((edgeIndex + 2) % 3); // CW vertex of edge (Edge's v1)
		WavefrontEdge edge = tri.getWavefront(edgeIndex);

		if (edge == null || edge.isDead()) {
			System.err.println("Warning: Constraint edge already null/dead for event " + event + ". Skipping.");
			return;
		}
		if (va.hasStopped() && vb.hasStopped()) {
			System.err.println("Warning: Vertices already stopped for constraint collapse event " + event + ". Skipping.");
			return;
		}
		// Additional check: Are the vertices *supposed* to collide here?
		Coordinate posA = va.getPositionAt(time);
		Coordinate posB = vb.getPositionAt(time);
		if (posA.distanceSq(posB) > SurfConstants.ZERO_NT_SQ) {
			System.err.printf(
					"Warning: Constraint collapse event %s scheduled but vertices V%d (%.6f, %.6f) and V%d (%.6f, %.6f) are not coincident at t=%.9f. DistSq=%.12e. Skipping.%n",
					event, va.id, posA.x, posA.y, vb.id, posB.x, posB.y, time, posA.distanceSq(posB));
			// This suggests an issue with event time calculation vs actual vertex movement.
			// Recompute event for triangle t and update queue instead of processing?
			updateEventsFor(Collections.singletonList(tri));
			return;
		}

		// 1. Stop vertices
		Coordinate stopPos = posA; // Use the calculated coincident position
		System.out.printf("    Stopping vertices V%d, V%d at (%.6f, %.6f)%n", va.id, vb.id, stopPos.x, stopPos.y);
		va.stop(time, stopPos);
		vb.stop(time, stopPos);

		// 2. DCEL Link: Head-to-Head between colliding vertices across the edge
		// Link va's side 1 (CW, facing collapsing edge) to vb's side 0 (CCW, facing
		// collapsing edge)
		System.out.printf("    DCEL Link H2H: V%d (side 1) <-> V%d (side 0)%n", va.id, vb.id);
		va.linkHeadToHead(1, vb, 0);

		// 3. Handle the triangle removal and vertex merging
		System.out.println("    Triggering constraint collapse part 2 for KT" + tri.id + " edge " + edgeIndex);
		handleConstraintCollapsePart2(tri, edgeIndex, time);

		// 4. Mark edge dead (done inside the helper)
		// Event queue update handled within handleConstraintCollapsePart2
	}

	/** Mirrors C++ KineticTriangulation::handle_spoke_collapse_event */
	private void handleSpokeCollapse(CollapseSpec event) {
		System.out.println("  Handling SPOKE_COLLAPSE for KT" + event.getTriangle().id + " edge " + event.getRelevantEdge());
		KineticTriangle tri = event.getTriangle();
		int edgeIndex = event.getRelevantEdge();
		double time = event.getTime();

		if (tri.isConstrained(edgeIndex)) {
			System.err.println("Error: Spoke collapse called on constrained edge " + edgeIndex + " for KT" + tri.id);
			// Recompute event?
			updateEventsFor(Collections.singletonList(tri));
			return;
		}
		KineticTriangle neighbor = tri.getNeighbor(edgeIndex);
		if (neighbor == null || neighbor.isDead()) {
			System.err.println("Warning: Spoke collapse neighbor null or dead for event " + event + ". Skipping.");
			return;
		}

		WavefrontVertex va = tri.getVertex((edgeIndex + 1) % 3); // CCW vertex of edge
		WavefrontVertex vb = tri.getVertex((edgeIndex + 2) % 3); // CW vertex of edge
		if (va.hasStopped() && vb.hasStopped()) {
			System.err.println("Warning: Vertices already stopped for spoke collapse event " + event + ". Skipping.");
			return;
		}
		// Check for coincidence
		Coordinate posA = va.getPositionAt(time);
		Coordinate posB = vb.getPositionAt(time);
		if (posA.distanceSq(posB) > SurfConstants.ZERO_NT_SQ) {
			System.err.printf(
					"Warning: Spoke collapse event %s scheduled but vertices V%d (%.6f, %.6f) and V%d (%.6f, %.6f) are not coincident at t=%.9f. DistSq=%.12e. Skipping.%n",
					event, va.id, posA.x, posA.y, vb.id, posB.x, posB.y, time, posA.distanceSq(posB));
			updateEventsFor(Arrays.asList(tri, neighbor));
			return;
		}

		// 1. Stop vertices
		Coordinate stopPos = posA; // Use the calculated coincident position
		System.out.printf("    Stopping vertices V%d, V%d at (%.6f, %.6f)%n", va.id, vb.id, stopPos.x, stopPos.y);
		va.stop(time, stopPos);
		vb.stop(time, stopPos);

		// 2. Remove neighbor links
		int neighborEdgeIndex = neighbor.indexOf(tri);
		if (neighborEdgeIndex == -1) {
			throw new IllegalStateException("Spoke collapse neighbor " + neighbor.id + " doesn't know triangle " + tri.id);
		}
		System.out.println("    Removing neighbor link between KT" + tri.id + " edge " + edgeIndex + " and KT" + neighbor.id + " edge " + neighborEdgeIndex);
		tri.setNeighbor(edgeIndex, null);
		neighbor.setNeighbor(neighborEdgeIndex, null);

		// 3. Handle the two triangles involved using the helper function
		System.out.println("    Triggering spoke collapse part 2 for KT" + tri.id + " edge " + edgeIndex);
		handleSpokeCollapsePart2(tri, edgeIndex, time);
		System.out.println("    Triggering spoke collapse part 2 for KT" + neighbor.id + " edge " + neighborEdgeIndex);
		handleSpokeCollapsePart2(neighbor, neighborEdgeIndex, time);

		// DCEL Link: Done implicitly inside
		// handleSpokeCollapsePart2/handleConstraintCollapsePart2
		// when the new vertex is created and linked.

		// Event queue updates handled within handleSpokeCollapsePart2 calls
	}

	/** Mirrors C++ KineticTriangulation::handle_split_or_flip_refine_event */
	private void handleSplitOrFlipRefine(CollapseSpec event) {
		System.out.println("  Handling SPLIT_OR_FLIP_REFINE for KT" + event.getTriangle().id + " edge " + event.getRelevantEdge());
		KineticTriangle tri = event.getTriangle();
		int edgeIndex = event.getRelevantEdge(); // The *constraint* edge index
		double time = event.getTime();
		WavefrontVertex vertexV = tri.getVertex(edgeIndex); // The vertex *opposite* the constraint

		if (!tri.isConstrained(edgeIndex)) {
			System.err.println("Error: SPLIT_OR_FLIP called on non-constrained edge " + edgeIndex + " for KT" + tri.id);
			updateEventsFor(Collections.singletonList(tri));
			return;
		}
		if (vertexV.hasStopped()) {
			System.err.println("Warning: Vertex V" + vertexV.id + " already stopped for SPLIT_OR_FLIP event " + event + ". Skipping.");
			return;
		}

		WavefrontVertex vertexA = tri.getVertex((edgeIndex + 1) % 3); // CCW end of constraint
		WavefrontVertex vertexB = tri.getVertex((edgeIndex + 2) % 3); // CW end of constraint

		Coordinate posV = vertexV.getPositionAt(time);
		Coordinate posA = vertexA.getPositionAt(time);
		Coordinate posB = vertexB.getPositionAt(time);
		LineSegment constraintSeg = new LineSegment(posA, posB);

		// Refine based on position of V relative to segment AB
		double distSqV_A = posV.distanceSq(posA);
		double distSqV_B = posV.distanceSq(posB);
		// Use projection factor or check collinearity and betweenness
		double projFactor = GeoUtils.projectionFactor(constraintSeg, posV);
		boolean onLine = Math.abs(GeoUtils.ptLineDist(posA, posB, posV)) < SurfConstants.ZERO_NT; // Check collinearity

		CollapseSpec refinedSpec = null;

		if (!onLine) {
			// This shouldn't happen if the event time calculation is correct.
			// The vertex should be *on* the supporting line at the event time.
			System.err.printf(
					"ERROR: SPLIT_OR_FLIP event %s scheduled, but vertex V%d (%.6f, %.6f) is not on supporting line of V%d-V%d at t=%.9f. Dist=%.12e. Recomputing event.%n",
					event, vertexV.id, posV.x, posV.y, vertexA.id, vertexB.id, time, GeoUtils.ptLineDist(posA, posB, posV));
			updateEventsFor(Collections.singletonList(tri));
			return;
		}

		// Vertex V is on the supporting line. Check its position relative to the
		// segment.
		if (distSqV_A < SurfConstants.ZERO_NT_SQ) { // V hits endpoint A (within tolerance)
			System.out.println("    Refining SPLIT_OR_FLIP to SPOKE_COLLAPSE (V hits A)");
			int spokeEdgeIndex = (edgeIndex + 2) % 3; // Edge opposite vertex B (edge VA)
			refinedSpec = tri.refineCollapseSpec(new CollapseSpec(CollapseType.SPOKE_COLLAPSE, time, tri, spokeEdgeIndex));
		} else if (distSqV_B < SurfConstants.ZERO_NT_SQ) { // V hits endpoint B (within tolerance)
			System.out.println("    Refining SPLIT_OR_FLIP to SPOKE_COLLAPSE (V hits B)");
			int spokeEdgeIndex = (edgeIndex + 1) % 3; // Edge opposite vertex A (edge VB)
			refinedSpec = tri.refineCollapseSpec(new CollapseSpec(CollapseType.SPOKE_COLLAPSE, time, tri, spokeEdgeIndex));
		} else if (projFactor > SurfConstants.ZERO_NT && projFactor < 1.0 - SurfConstants.ZERO_NT) { // V hits interior (robust check)
			System.out.println("    Refining SPLIT_OR_FLIP to actual SPLIT_EVENT");
			// Handle the split directly
			handleSplitEvent(event); // Pass original event spec
			return; // Split handled, no need to re-queue
		} else { // V hits supporting line outside segment AB -> Flip
			System.out.println("    Refining SPLIT_OR_FLIP to VERTEX_MOVES_OVER_SPOKE (Flip)");
			int flipEdgeIndex;
			double longestSpokeSqLen;
			// Determine longest spoke: edge opposite the endpoint V is *further* from
			if (distSqV_A > distSqV_B) { // V is further from A, closer to B. Flip edge VA (opposite B).
				flipEdgeIndex = (edgeIndex + 2) % 3; // Flip edge opposite B
				longestSpokeSqLen = distSqV_A;
				System.out.println("      V is further from A. Flipping edge " + flipEdgeIndex + " (opp B). Longest spoke sq: " + longestSpokeSqLen);
			} else { // V is further from B, closer to A. Flip edge VB (opposite A).
				flipEdgeIndex = (edgeIndex + 1) % 3; // Flip edge opposite A
				longestSpokeSqLen = distSqV_B;
				System.out.println("      V is further from B. Flipping edge " + flipEdgeIndex + " (opp A). Longest spoke sq: " + longestSpokeSqLen);
			}

			// Check if the edge to flip is internal (a spoke)
			if (tri.isConstrained(flipEdgeIndex)) {
				System.err.println("Error: Cannot refine SPLIT to FLIP because the required flip edge " + flipEdgeIndex + " is constrained! Event: " + event);
				// This indicates a potential issue in event calculation or geometry.
				refinedSpec = tri.refineCollapseSpec(new CollapseSpec(CollapseType.INVALID_EVENT, time, tri));
			} else {
				refinedSpec = tri.refineCollapseSpec(new CollapseSpec(CollapseType.VERTEX_MOVES_OVER_SPOKE, time, tri, flipEdgeIndex, longestSpokeSqLen));
			}
		}

		// If refinement occurred and didn't handle directly (i.e., not split), re-queue
		// the refined event.
		if (refinedSpec != null) {
			System.out.println("    Re-queuing refined event: " + refinedSpec);
			eventQueue.addOrUpdate(refinedSpec); // addOrUpdate handles replacing old event
		}
	}

	/**
	 * Handler for actual Split Events, after refinement. Mirrors C++
	 * KineticTriangulation::handle_split_event
	 */
	private void handleSplitEvent(CollapseSpec event) {
		System.out.println("    Handling actual SPLIT_EVENT for KT" + event.getTriangle().id + " edge " + event.getRelevantEdge());
		KineticTriangle tri = event.getTriangle();
		int constraintEdgeIndex = event.getRelevantEdge();
		double time = event.getTime();
		WavefrontVertex vertexV = tri.getVertex(constraintEdgeIndex); // Vertex opposite constraint
		WavefrontEdge edgeE = tri.getWavefront(constraintEdgeIndex);

		if (edgeE == null || edgeE.isDead()) {
			System.err.println("Warning: Split edge null or dead for event " + event + ". Skipping.");
			return;
		}
		if (vertexV.hasStopped()) {
			System.err.println("Warning: Split vertex V" + vertexV.id + " already stopped for event " + event + ". Skipping.");
			return;
		}

		// 1. Stop vertex V
		Coordinate stopPos = vertexV.getPositionAt(time);
		System.out.printf("      Stopping vertex V%d at (%.6f, %.6f)%n", vertexV.id, stopPos.x, stopPos.y);
		vertexV.stop(time, stopPos);

		// 2. Split the edge E
		System.out.println("      Splitting edge WFE" + edgeE.id);
		Pair<WavefrontEdge, WavefrontEdge> newEdges = kineticTriangulation.splitWavefrontEdge(edgeE);
		WavefrontEdge newEdgeA = newEdges.getLeft(); // Corresponds to original vertex A side (CCW side of constraint)
		WavefrontEdge newEdgeB = newEdges.getRight(); // Corresponds to original vertex B side (CW side of constraint)
		System.out.println("      New edges: WFE" + newEdgeA.id + " (A side), WFE" + newEdgeB.id + " (B side)");

		// 3. Create two new vertices at the split point
		WavefrontEdge edgeIncidentA = vertexV.getIncidentEdge(1); // Edge CW of V (towards original A side)
		WavefrontEdge edgeIncidentB = vertexV.getIncidentEdge(0); // Edge CCW of V (towards original B side)

		System.out.println("      Creating new vertices at split point.");
		WavefrontVertex newVertexA = kineticTriangulation.createAndAddVertex(stopPos, // pos_zero (needs back-calculation?)
				stopPos, time, newEdgeA, edgeIncidentA, true // Mark as split vertex
		);
		WavefrontVertex newVertexB = kineticTriangulation.createAndAddVertex(stopPos, // pos_zero (needs back-calculation?)
				stopPos, time, edgeIncidentB, newEdgeB, true // Mark as split vertex
		);
		System.out.println("      New vertices: V" + newVertexA.id + " (A side), V" + newVertexB.id + " (B side)");

		// 4. Link new vertices on the split edge (Update edge endpoints)
		System.out.println("      Linking new vertices to split edge ends.");
		newEdgeA.setVertex(1, newVertexA); // newEdgeA ends at newVertexA
		newEdgeB.setVertex(0, newVertexB); // newEdgeB starts at newVertexB
		// Link original endpoints to the new edges
		if (newEdgeA.getVertex(0) != null)
			newEdgeA.getVertex(0).setIncidentWavefrontEdge(1, newEdgeA);
		if (newEdgeB.getVertex(1) != null)
			newEdgeB.getVertex(1).setIncidentWavefrontEdge(0, newEdgeB);

		// 5. Update incident triangles around vertex V
		System.out.println("      Updating triangles incident to V" + vertexV.id);
		List<KineticTriangle> trianglesToUpdate = new ArrayList<>();
		// This requires iterating the triangle fan around V. Assume KT provides this.
		// The logic is: walk from edgeIncidentB (CCW) towards edgeIncidentA (CW).
		// Triangles encountered before reaching the split edge E belong to newVertexB.
		// Triangles encountered after reaching the split edge E belong to newVertexA.
		// The original triangle 'tri' is skipped as it's dying.
		Collection<KineticTriangle> incidentTriangles = kineticTriangulation.findIncidentTriangles(vertexV);
		boolean passedSplit = false;
		WavefrontEdge currentEdge = edgeIncidentB; // Start CCW
		// Hypothetical loop structure:
		for (KineticTriangle incidentT : /* Ordered list from KT.getOrderedIncidentTriangles(vertexV, edgeIncidentB) */ incidentTriangles) {
			if (incidentT == tri || incidentT.isDead())
				continue;

			int vIndex = incidentT.indexOf(vertexV);
			if (vIndex == -1)
				continue; // Should not happen

			// Check if the edge CW of V in incidentT is edgeE (or its successors
			// newEdgeA/B?)
			// This logic is complex and needs careful geometry/neighbor checks.
			// Simplified: Assume we know based on which neighbor (na/nb) the triangle is.
			KineticTriangle neighborA = tri.getNeighbor((constraintEdgeIndex + 2) % 3); // Neighbor across spoke VA
			KineticTriangle neighborB = tri.getNeighbor((constraintEdgeIndex + 1) % 3); // Neighbor across spoke VB

			if (/* incidentT is between edgeIncidentB and edgeE */ kineticTriangulation.isTriangleBetweenEdges(incidentT, vertexV, edgeIncidentB, edgeE)) { // Hypothetical
																																							// check
				System.out.println("        Updating KT" + incidentT.id + " with V" + newVertexB.id);
				incidentT.setVertex(vIndex, newVertexB);
				trianglesToUpdate.add(incidentT);
			} else { // incidentT is between edgeE and edgeIncidentA
				System.out.println("        Updating KT" + incidentT.id + " with V" + newVertexA.id);
				incidentT.setVertex(vIndex, newVertexA);
				trianglesToUpdate.add(incidentT);
			}
		}
		System.err.println("TODO: Implement robust triangle fan update in handleSplitEvent!");

		// 6. Update the neighbors that bordered the original triangle 'tri' along
		// spokes
		KineticTriangle neighborA = tri.getNeighbor((constraintEdgeIndex + 2) % 3); // Neighbor across spoke VA
		KineticTriangle neighborB = tri.getNeighbor((constraintEdgeIndex + 1) % 3); // Neighbor across spoke VB

		if (neighborA != null && !neighborA.isDead()) {
			int idx = neighborA.indexOf(tri);
			System.out.println("      Updating neighbor KT" + neighborA.id + " edge " + idx + " to border WFE" + newEdgeA.id);
			neighborA.setWavefront(idx, newEdgeA); // Was spoke, now constraint newEdgeA
			if (!trianglesToUpdate.contains(neighborA))
				trianglesToUpdate.add(neighborA);
		}
		if (neighborB != null && !neighborB.isDead()) {
			int idx = neighborB.indexOf(tri);
			System.out.println("      Updating neighbor KT" + neighborB.id + " edge " + idx + " to border WFE" + newEdgeB.id);
			neighborB.setWavefront(idx, newEdgeB); // Was spoke, now constraint newEdgeB
			if (!trianglesToUpdate.contains(neighborB))
				trianglesToUpdate.add(neighborB);
		}

		// 7. Mark original triangle dead
		System.out.println("      Marking original triangle KT" + tri.id + " dead.");
		tri.markDead();

		// 8. DCEL Links
		System.out.printf("      DCEL Link H2T: V%d (side 0) -> V%d%n", vertexV.id, newVertexB.id);
		vertexV.linkHeadToTail(0, newVertexB); // V's CCW side links to newVertexB's tail
		System.out.printf("      DCEL Link H2T: V%d (side 1) -> V%d%n", vertexV.id, newVertexA.id);
		vertexV.linkHeadToTail(1, newVertexA); // V's CW side links to newVertexA's tail
		System.out.printf("      DCEL Link T2T: V%d <-> V%d%n", newVertexA.id, newVertexB.id);
		newVertexA.linkTailToTail(newVertexB); // Link the new vertices together tail-to-tail

		// 9. Update event queue
		updateEventsFor(trianglesToUpdate);

		// Increment C++ stats counters
		avgTrianglesPerSplitEventSum += trianglesToUpdate.size(); // Approx count
		avgTrianglesPerSplitEventCtr++;
		maxTrianglesPerSplitEvent = Math.max(maxTrianglesPerSplitEvent, trianglesToUpdate.size());
	}

	/** Mirrors C++ KineticTriangulation::handle_vertex_moves_over_spoke_event */
	private void handleVertexMovesOverSpoke(CollapseSpec event) {
		System.out.println("  Handling VERTEX_MOVES_OVER_SPOKE for KT" + event.getTriangle().id + " edge " + event.getRelevantEdge());
		KineticTriangle tri = event.getTriangle();
		int edgeIndex = event.getRelevantEdge(); // The spoke edge index being flipped
		double time = event.getTime();

		if (tri.isConstrained(edgeIndex)) {
			System.err.println("Error: Attempting to flip constrained edge " + edgeIndex + " in VERTEX_MOVES_OVER_SPOKE for KT" + tri.id + ". Event " + event);
			updateEventsFor(Collections.singletonList(tri));
			return;
		}
		KineticTriangle neighbor = tri.getNeighbor(edgeIndex);
		if (neighbor == null || neighbor.isDead()) {
			System.err.println("Warning: Flip neighbor null or dead for event " + event + ". Skipping.");
			return;
		}

		// Add geometric check: Is the flip valid? (Check if vertex opposite spoke is
		// inside circumcircle of the other 3)
		// Or simpler check: Is the resulting diagonal valid (orientation checks)?
		// The C++ code seems to rely on the event calculation being correct.

		// 1. Perform the flip
		System.out.println("    Performing flip on KT" + tri.id + " edge " + edgeIndex);
		kineticTriangulation.flipEdge(tri, edgeIndex);

		// 2. Update events for the two affected triangles
		List<KineticTriangle> trianglesToUpdate = Arrays.asList(tri, neighbor);
		updateEventsFor(trianglesToUpdate);
		System.out.println("    Completed VERTEX_MOVES_OVER_SPOKE (Flip).");
	}

	/** Mirrors C++ KineticTriangulation::handle_ccw_vertex_leaves_ch_event */
	private void handleCcwVertexLeavesCH(CollapseSpec event) {
		System.out.println("  Handling CCW_VERTEX_LEAVES_CH for KT" + event.getTriangle().id + " infinite_vtx_idx " + event.getRelevantEdge());
		KineticTriangle tri = event.getTriangle();
		int infiniteVertexIndex = event.getRelevantEdge(); // Index of the infinite vertex
		int edgeToFlip = (infiniteVertexIndex + 2) % 3; // Edge CW of infinite vertex
		double time = event.getTime();

		// Basic checks
		WavefrontVertex infVtx = tri.getVertex(infiniteVertexIndex);
		if (infVtx == null || !infVtx.isInfinite()) {
			System.err
					.println("Error: CCW_VERTEX_LEAVES_CH event triangle KT" + tri.id + " vertex " + infiniteVertexIndex + " is not infinite. Event: " + event);
			updateEventsFor(Collections.singletonList(tri));
			return;
		}
		if (tri.isConstrained(edgeToFlip)) {
			System.err.println("Error: Attempting to flip constrained edge " + edgeToFlip + " in CCW_VERTEX_LEAVES_CH for KT" + tri.id + ". Event: " + event);
			updateEventsFor(Collections.singletonList(tri));
			return;
		}
		KineticTriangle neighbor = tri.getNeighbor(edgeToFlip);
		if (neighbor == null || neighbor.isDead()) {
			System.err.println("Warning: Flip neighbor null or dead for CCW_VERTEX_LEAVES_CH event " + event + ". Skipping.");
			return;
		}
		int neighborEdgeIndex = neighbor.indexOf(tri);
		WavefrontVertex neighborOppositeVtx = neighbor.getVertex(neighborEdgeIndex);
		if (neighborOppositeVtx == null || !neighborOppositeVtx.isInfinite()) {
			System.err.println(
					"Error: Neighbor KT" + neighbor.id + " in CCW_VERTEX_LEAVES_CH event vertex opposite KT" + tri.id + " is not infinite. Event: " + event);
			updateEventsFor(Arrays.asList(tri, neighbor)); // Update both just in case
			return;
		}

		// Geometric check: Are the 3 finite vertices collinear at this time?
		WavefrontVertex v_cw = tri.getVertex((infiniteVertexIndex + 2) % 3); // Vertex being flipped over
		WavefrontVertex v_ccw = tri.getVertex((infiniteVertexIndex + 1) % 3); // Vertex staying on CH edge
		WavefrontVertex v_neighbor_opp = neighbor.getVertex(neighborEdgeIndex); // == infVtx
		WavefrontVertex v_neighbor_other = neighbor.getVertex((neighborEdgeIndex + 1) % 3); // Other finite vertex on CH

		Coordinate p_cw = v_cw.getPositionAt(time);
		Coordinate p_ccw = v_ccw.getPositionAt(time);
		Coordinate p_neighbor_other = v_neighbor_other.getPositionAt(time);

		if (Math.abs(GeoUtils.orientation(p_ccw, p_cw, p_neighbor_other)) > SurfConstants.ZERO_NT * 10) { // Allow some tolerance
			System.err.printf(
					"Warning: CCW_VERTEX_LEAVES_CH event %s scheduled, but vertices V%d, V%d, V%d are not collinear (orientation = %.3e) at t=%.9f. Skipping flip.%n",
					event, v_ccw.id, v_cw.id, v_neighbor_other.id, GeoUtils.orientation(p_ccw, p_cw, p_neighbor_other), time);
			updateEventsFor(Arrays.asList(tri, neighbor));
			return;
		}

		// 1. Perform the flip
		System.out.println("    Performing flip on KT" + tri.id + " edge " + edgeToFlip);
		kineticTriangulation.flipEdge(tri, edgeToFlip);

		// 2. Update events
		List<KineticTriangle> trianglesToUpdate = Arrays.asList(tri, neighbor);
		updateEventsFor(trianglesToUpdate);
		System.out.println("    Completed CCW_VERTEX_LEAVES_CH (Flip).");
	}

	/**
	 * Mirrors C++
	 * KineticTriangulation::handle_face_with_infintely_fast_opposing_vertex
	 */
	private void handleFaceWithInfinitelyFastOpposingVertex(CollapseSpec event) {
		System.out.println("  Handling FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING for KT" + event.getTriangle().id);
		// TODO: Port complex logic from C++
		// - Find the infinitely fast vertex 'v' (type OPPOSING)
		// - Find the most CW triangle 't' around 'v' using incident_faces_iterator
		// logic
		// - Flip away all spokes incident to 'v' in the fan, updating 't' if necessary
		// - Identify vertices v_cw, v_ccw adjacent to 'v' in the final 't' (which now
		// must have constraints)
		// - Determine which edge collapses (v-vcw or v-vccw) based on distance or if
		// one adjacent vertex is also fast. Handle tie (spoke_collapse true).
		// - If spoke_collapse: Stop all 3 vertices (v, v_cw, v_ccw), mark edges dead,
		// DCEL link v to both, handle neighbor across the spoke, mark t dead.
		// - If not spoke_collapse: Identify the vertex 'o' to keep. Stop 'v' and the
		// other vertex ('o' stops normally). DCEL link v to the other. Call
		// handleConstraintCollapsePart2 for 't' and the collapsed edge index.
		// - Update events for triangles modified by flips.
		fail("handleFaceWithInfinitelyFastOpposingVertex logic not implemented.");
	}

	/**
	 * Mirrors C++
	 * KineticTriangulation::handle_face_with_infintely_fast_weighted_vertex
	 */
	private void handleFaceWithInfinitelyFastWeightedVertex(CollapseSpec event) {
		System.out.println("  Handling FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED for KT" + event.getTriangle().id + " edge " + event.getRelevantEdge());
		// TODO: Port complex logic from C++
		// - Identify the fast vertex 'v_fast' (type WEIGHTED) based on
		// event.relevantEdge.
		// - Identify the winning_edge (the one specified in event.relevantEdge) and
		// losing_edge incident to v_fast.
		// - Find the most CW triangle 't' around 'v_fast' (similar to OPPOSING case,
		// may involve walking across edges).
		// - Flip away spokes incident to 'v_fast' if the triangle is bounded. Special
		// handling if unbounded (check C++).
		// - Stop 'v_fast' and the vertex 'o' at the other end of the losing_edge.
		// - DCEL link 'v_fast' to 'o' across the losing edge.
		// - Call handleConstraintCollapsePart2 for 't' and the index corresponding to
		// the losing_edge.
		// - Update events for triangles modified by flips.
		fail("handleFaceWithInfinitelyFastWeightedVertex logic not implemented.");
	}

	// --- Helper Methods ---

	/**
	 * Helper to handle the second part of a spoke collapse for a given triangle.
	 * Determines if the triangle collapses completely or just merges vertices.
	 * Mirrors C++ KineticTriangulation::do_spoke_collapse_part2.
	 *
	 * @param t         The triangle experiencing the collapse on edgeIndex.
	 * @param edgeIndex The index of the collapsed spoke (where neighbor was
	 *                  removed).
	 * @param time      The current simulation time.
	 */
	private void handleSpokeCollapsePart2(KineticTriangle t, int edgeIndex, double time) {
		System.out.println("    Helper: Spoke Collapse Part 2 for KT" + t.id + " edge " + edgeIndex);
		if (t.isDead()) {
			System.out.println("      KT" + t.id + " already dead. Skipping.");
			return;
		}

		WavefrontVertex va = t.getVertex((edgeIndex + 1) % 3); // CCW vertex of spoke
		WavefrontVertex vb = t.getVertex((edgeIndex + 2) % 3); // CW vertex of spoke
		WavefrontVertex v = t.getVertex(edgeIndex); // Vertex opposite spoke

		// Vertices va and vb should have been stopped by the caller
		// (handleSpokeCollapse)
		if (!va.hasStopped() || !vb.hasStopped()) {
			System.err.printf("Error in handleSpokeCollapsePart2: Vertices V%d or V%d not stopped for KT%d edge %d.%n", va.id, vb.id, t.id, edgeIndex);
			// Attempt to stop them now? Risky.
			Coordinate stopPos = va.getPositionAt(time);
			va.stop(time, stopPos);
			vb.stop(time, stopPos);
			// return; // Or try to continue?
		}
		// Check if vertex v is also coincident (triangle collapses fully)
		Coordinate posV = v.getPositionAt(time);
		Coordinate posA = va.getPosStop(); // Use the definitive stop position

		// Need special handling for infinite/fast vertices? C++ does this.
		boolean v_collapses_too = !v.isInfinite() && (posV.distanceSq(posA) < SurfConstants.ZERO_NT_SQ);
		boolean v_infinitely_fast = v.getInfiniteSpeed() != WavefrontVertex.InfiniteSpeedType.NONE;

		if (v_collapses_too || (v_infinitely_fast && !t.hasNeighbor(edgeIndex + 1) && !t.hasNeighbor(edgeIndex + 2))) { // Case from C++ for fully constrained
																														// fast vertex
			System.out.println("      Triangle KT" + t.id + " collapses completely.");
			t.markDead(); // Mark dying
			if (!v.hasStopped()) {
				System.out.printf("      Stopping vertex V%d at (%.6f, %.6f)%n", v.id, posA.x, posA.y);
				v.stop(time, posA); // Stop v at the common point
			}

			// Process remaining neighbors/constraints
			for (int i = 1; i <= 2; i++) {
				int currentEdgeIndex = (edgeIndex + i) % 3;
				KineticTriangle neighbor = t.getNeighbor(currentEdgeIndex);
				WavefrontEdge wfe = t.getWavefront(currentEdgeIndex);

				if (neighbor != null) {
					System.out.println("      Edge " + currentEdgeIndex + ": Spoke to KT" + neighbor.id);
					if (!neighbor.isDead()) {
						int idxInN = neighbor.indexOf(t);
						System.out.println("        Removing neighbor link from KT" + neighbor.id + " edge " + idxInN);
						neighbor.setNeighbor(idxInN, null);
						// Trigger collapse for neighbor
						System.out.println("        Triggering spoke collapse part 2 for neighbor KT" + neighbor.id + " edge " + idxInN);
						handleSpokeCollapsePart2(neighbor, idxInN, time);
					} else {
						System.out.println("        Neighbor KT" + neighbor.id + " already dead.");
					}
				} else if (wfe != null) {
					System.out.println("      Edge " + currentEdgeIndex + ": Constraint WFE" + wfe.id);
					if (!wfe.isDead()) {
						wfe.markDead();
						// DCEL Link head-to-head
						WavefrontVertex edgeV_ccw = t.getVertex((currentEdgeIndex + 1) % 3);
						WavefrontVertex edgeV_cw = t.getVertex((currentEdgeIndex + 2) % 3);
						System.out.printf("        DCEL Link H2H: V%d (side 0) <-> V%d (side 1)%n", edgeV_cw.id, edgeV_ccw.id);
						edgeV_cw.linkHeadToHead(0, edgeV_ccw, 1);
					} else {
						System.out.println("        Constraint WFE" + wfe.id + " already dead.");
					}
				}
			}
			// No new vertex created, only handle neighbors and mark dead.
			// Event queue updates handled recursively by neighbor calls.
		} else {
			// Triangle does not collapse completely, treat as constraint collapse on the
			// spoke
			System.out.println("      Triangle KT" + t.id + " does not collapse completely. Treating as constraint collapse on spoke edge " + edgeIndex);
			handleConstraintCollapsePart2(t, edgeIndex, time);
		}
	}

	/**
	 * Helper to handle the second part of a constraint collapse (or spoke collapse
	 * treated as such). Creates the new merged vertex, updates
	 * neighbors/constraints, and updates the event queue. Mirrors C++
	 * KineticTriangulation::do_constraint_collapse_part2.
	 *
	 * @param t         The triangle being removed.
	 * @param edgeIndex The index of the collapsing constraint or spoke.
	 * @param time      The current simulation time.
	 */
	private void handleConstraintCollapsePart2(KineticTriangle t, int edgeIndex, double time) {
		System.out.println("    Helper: Constraint Collapse Part 2 for KT" + t.id + " edge " + edgeIndex);
		if (t.isDead()) {
			System.out.println("      KT" + t.id + " already dead. Skipping.");
			return;
		}

		WavefrontVertex va = t.getVertex((edgeIndex + 1) % 3); // CCW vertex
		WavefrontVertex vb = t.getVertex((edgeIndex + 2) % 3); // CW vertex

		if (!va.hasStopped() || !vb.hasStopped()) {
			System.err.println("Error in handleConstraintCollapsePart2: Vertices not stopped for KT" + t.id);
			Coordinate stopPos = va.getPositionAt(time); // Attempt recovery
			va.stop(time, stopPos);
			vb.stop(time, stopPos);
		}
		Coordinate stopPos = va.getPosStop(); // Use the definitive stop position

		System.out.println("      Marking KT" + t.id + " dead.");
		t.markDead(); // Mark dying

		// 1. Transfer constraints away from the dying triangle t
		System.out.println("      Transferring constraints away from KT" + t.id);
		List<KineticTriangle> neighborsModifiedByConstraintMove = kineticTriangulation.moveConstraintsToNeighbors(t, edgeIndex);

		// 2. Identify incident edges for the new vertex
		WavefrontEdge edgeIncidentA = va.getIncidentEdge(0); // Edge CCW of va (if exists)
		WavefrontEdge edgeIncidentB = vb.getIncidentEdge(1); // Edge CW of vb (if exists)
		System.out.println("      Incident edges for new vertex: WFE" + (edgeIncidentA != null ? edgeIncidentA.id : "null") + " (A side), WFE"
				+ (edgeIncidentB != null ? edgeIncidentB.id : "null") + " (B side)");
		if (edgeIncidentA == null || edgeIncidentB == null) {
			System.err.println("Error: Cannot find incident edges for merging V" + va.id + " and V" + vb.id);
			// This indicates a problem at the boundary or with DCEL links.
			// Cannot create new vertex. Maybe just stop here?
			return;
		}

		// 3. Create the new merged vertex
		System.out.println("      Creating merged vertex.");
		WavefrontVertex newVertex = kineticTriangulation.createAndAddVertex(stopPos, // pos_zero (needs back-calculation if velocity non-zero)
				stopPos, time, edgeIncidentA, edgeIncidentB);
		System.out.println("      Created merged vertex: V" + newVertex.id);

		// 4. DCEL Link: Link old vertices' outer paths to the new vertex
		// Link tail of A's previous vertex to newVertex (side 0)
		if (va.getPrevVertex(0) != null) {
			System.out.printf("      DCEL Link H2T: V%d (side 0) -> V%d%n", va.getPrevVertex(0).id, newVertex.id);
			va.getPrevVertex(0).linkHeadToTail(0, newVertex);
		} else {
			System.out.println("      V" + va.id + " has no prev vertex on side 0 (boundary?).");
		}
		// Link newVertex (side 1) to head of B's next vertex
		if (vb.getNextVertex(1) != null) {
			System.out.printf("      DCEL Link H2T: V%d (side 1) -> V%d%n", newVertex.id, vb.getNextVertex(1).id);
			newVertex.linkHeadToTail(1, vb.getNextVertex(1));
		} else {
			System.out.println("      V" + vb.id + " has no next vertex on side 1 (boundary?).");
		}

		// 5. Update incident triangles around the OLD vertices (va, vb) to point to the
		// NEW vertex
		System.out.println("      Updating triangles incident to old V" + va.id + " and V" + vb.id);
		Set<KineticTriangle> trianglesToUpdate = new HashSet<>();
		// Add neighbors potentially modified by constraint transfer earlier
		trianglesToUpdate.addAll(neighborsModifiedByConstraintMove);

		Collection<KineticTriangle> incidentToVa = kineticTriangulation.findIncidentTriangles(va);
		for (KineticTriangle incidentT : incidentToVa) {
			if (incidentT == t || incidentT.isDead())
				continue;
			int idx = incidentT.indexOf(va);
			if (idx != -1) {
				System.out.println("        Updating KT" + incidentT.id + " vertex " + idx + " from V" + va.id + " to V" + newVertex.id);
				incidentT.setVertex(idx, newVertex);
				trianglesToUpdate.add(incidentT);
			}
		}
		Collection<KineticTriangle> incidentToVb = kineticTriangulation.findIncidentTriangles(vb);
		for (KineticTriangle incidentT : incidentToVb) {
			if (incidentT == t || incidentT.isDead())
				continue;
			int idx = incidentT.indexOf(vb);
			if (idx != -1) {
				System.out.println("        Updating KT" + incidentT.id + " vertex " + idx + " from V" + vb.id + " to V" + newVertex.id);
				incidentT.setVertex(idx, newVertex);
				trianglesToUpdate.add(incidentT);
			}
		}
		System.out.println("      Total triangles needing event update: " + trianglesToUpdate.size());

		// 6. Update neighbor links for triangles that bordered t along spokes
		KineticTriangle neighborCW = t.getNeighbor((edgeIndex + 1) % 3); // Neighbor opposite vb
		KineticTriangle neighborCCW = t.getNeighbor((edgeIndex + 2) % 3); // Neighbor opposite va

		if (neighborCW != null && neighborCCW != null) {
			System.out.println("      Linking neighbors KT" + neighborCW.id + " and KT" + neighborCCW.id);
			int idxInNCW = neighborCW.indexOf(t);
			int idxInNCCW = neighborCCW.indexOf(t);
			if (idxInNCW != -1)
				neighborCW.setNeighbor(idxInNCW, neighborCCW);
			if (idxInNCCW != -1)
				neighborCCW.setNeighbor(idxInNCCW, neighborCW);
			// Vertex updates already handled in step 5
		} else {
			System.out.println("      One or both neighbors across spokes were null (likely constraints handled earlier).");
		}

		// 7. Mark the original collapsing wavefront edge dead (if it was a constraint)
		WavefrontEdge collapsingEdge = t.getWavefront(edgeIndex); // Get original WFE before marking t dead
		if (collapsingEdge != null && !collapsingEdge.isDead()) {
			System.out.println("      Marking collapsing constraint WFE" + collapsingEdge.id + " dead.");
			collapsingEdge.markDead();
		} else {
			System.out.println("      Collapsing edge was a spoke or already dead.");
		}

		// 8. Update Event Queue for all affected triangles
		updateEventsFor(trianglesToUpdate);

		// Increment C++ stats counters
		avgTrianglesPerEdgeEventSum += trianglesToUpdate.size(); // Approx count
		avgTrianglesPerEdgeEventCtr++;
		maxTrianglesPerEdgeEvent = Math.max(maxTrianglesPerEdgeEvent, trianglesToUpdate.size());
	}

	/** Helper to invalidate and recompute/add events for a list of triangles. */
	private void updateEventsFor(Collection<KineticTriangle> triangles) {
		// Use a Set to avoid updating the same triangle multiple times if it was
		// affected in different ways
		Set<KineticTriangle> uniqueTriangles = new HashSet<>();
		for (KineticTriangle t : triangles) {
			if (t != null && !t.isDead()) {
				uniqueTriangles.add(t);
			}
		}

		if (uniqueTriangles.isEmpty()) {
			return;
		}

		System.out.print("    Updating events for KTs: ");
		for (KineticTriangle t : uniqueTriangles) {
			System.out.print(t.id + " ");
			// getCollapseSpec recalculates if needed
			CollapseSpec newSpec = t.getCollapseSpec(currentTime);
			eventQueue.addOrUpdate(newSpec);
		}
		System.out.println();
	}

	/**
	 * Helper to update events for newly created elements during a step. Mirrors C++
	 * process_check_refinement_queue somewhat.
	 */
	private void updateEventsForNewElements() {
		// Gather all triangles adjacent to newly created vertices and edges.
		Set<KineticTriangle> trianglesToUpdate = new HashSet<>();

		for (WavefrontVertex newV : kineticTriangulation.getNewVerticesThisStep()) {
			if (newV.isDead())
				continue; // Should not happen ideally
			trianglesToUpdate.addAll(kineticTriangulation.findIncidentTriangles(newV));
		}

		for (WavefrontEdge newE : kineticTriangulation.getNewEdgesThisStep()) {
			if (newE.isDead())
				continue;
			KineticTriangle incidentT = newE.getIncidentTriangle();
			if (incidentT != null && !incidentT.isDead()) {
				trianglesToUpdate.add(incidentT);
				// Also need the triangle on the *other* side if this edge becomes a boundary
				// This requires more sophisticated neighbor finding in KT or assuming the edge
				// creation handles it.
			}
		}

		if (!trianglesToUpdate.isEmpty()) {
			System.out.println("  Updating events for " + trianglesToUpdate.size() + " triangles due to new elements created this step.");
			updateEventsFor(trianglesToUpdate);
		}
	}

	/**
	 * Finalization: Builds the final DCEL structure. Mirrors C++
	 * WavefrontPropagator::finalize
	 */
	private void finalizePropagation() {
		if (finalized) {
			return;
		}
		System.out.println("Finalizing propagation steps...");
		// Call the method on KineticTriangulation that builds the DCEL from the final
		// state.
		// This corresponds to C++
		// KineticTriangulation::create_remaining_skeleton_dcel()
		try {
			kineticTriangulation.createRemainingSkeletonDCEL();
			this.finalized = true;
			this.currentComponent = -1; // Reset component tracking
			// Final update for stats
			updateEventTimingStats(-1.0, -1); // Use dummy values to trigger final calculation
			System.out.println("  DCEL creation complete. Propagation finalized.");
		} catch (Exception e) {
			System.err.println("FATAL ERROR during DCEL finalization: " + e.getMessage());
			e.printStackTrace();
			// Mark as finalized anyway to prevent further processing attempts
			this.finalized = true;
		}
	}

	/** Mirrors C++ KineticTriangulation::update_event_timing_stats */
	private void updateEventTimingStats(double now, int component) {
		// Use a small tolerance for comparing double times
		if (Math.abs(now - lastEventTime) > TIME_COMPARISON_TOLERANCE || component != lastEventComponent) {
			if (eventCounter > 1) { // Don't count the very first "event" start
				maxEventsPerTime = Math.max(maxEventsPerTime, eventsPerCurrentEventTime);
				avgEventsPerTimeSum += eventsPerCurrentEventTime;
				avgEventsPerTimeCtr++;
			}
			lastEventTime = now;
			lastEventComponent = component;
			eventsPerCurrentEventTime = (now >= 0) ? 1 : 0; // Reset counter (handle dummy call from finalize)
		} else {
			eventsPerCurrentEventTime++;
		}
	}

	private void printStatistics() {
		System.out.println("\n--- Propagation Statistics ---");
		System.out.println("Total Events Processed: " + eventCounter);
		System.out.println("Final Simulation Time: " + currentTime);

		System.out.println("\nEvent Type Distribution:");
		long totalEventsLogged = eventTypeCounter[CollapseType.UNDEFINED.ordinal()];
		for (CollapseType type : CollapseType.values()) {
			if (eventTypeCounter[type.ordinal()] > 0 && type != CollapseType.UNDEFINED) {
				System.out.printf("  %-40s: %8d (%.2f%%)%n", type, eventTypeCounter[type.ordinal()],
						(double) 100 * eventTypeCounter[type.ordinal()] / totalEventsLogged);
			}
		}

		System.out.println("\nEvent Concurrency:");
		System.out.println("  Max Events at Same Time Step: " + maxEventsPerTime);
		if (avgEventsPerTimeCtr > 0) {
			System.out.printf("  Avg Events at Same Time Step: %.2f%n", (double) avgEventsPerTimeSum / avgEventsPerTimeCtr);
		}

		System.out.println("\nTriangle Updates per Event Type:");
		if (avgTrianglesPerEdgeEventCtr > 0) {
			System.out.println("  Edge Collapse Events: ");
			System.out.println("    Max Affected Tris: " + maxTrianglesPerEdgeEvent);
			System.out.printf("    Avg Affected Tris: %.2f%n", (double) avgTrianglesPerEdgeEventSum / avgTrianglesPerEdgeEventCtr);
		}
		if (avgTrianglesPerSplitEventCtr > 0) {
			System.out.println("  Split Events: ");
			System.out.println("    Max Affected Tris: " + maxTrianglesPerSplitEvent);
			System.out.printf("    Avg Affected Tris: %.2f%n", (double) avgTrianglesPerSplitEventSum / avgTrianglesPerSplitEventCtr);
		}
		System.out.println("----------------------------\n");
	}

	// --- Getters ---
	public double getCurrentTime() {
		return currentTime;
	}

	public long getEventCounter() {
		return eventCounter;
	}

	public boolean isFinalized() {
		return finalized;
	}

	public SkeletonDCEL getSkeleton() {
		if (!finalized) {
			System.err.println("Warning: Requesting skeleton before finalization.");
			// Optionally finalize automatically:
			// finalizePropagation();
		}
		return kineticTriangulation.getSkeleton();
	}

	// --- Placeholder for unimplemented logic ---
	private void fail(String message) {
		System.err.println("TODO/FAIL: " + message);
		// In a real implementation, you might throw an UnsupportedOperationException
		// throw new UnsupportedOperationException("Not implemented yet: " + message);
	}
}