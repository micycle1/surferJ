package com.github.micycle1.surferj.kinetics;

import static com.github.micycle1.surferj.TriangulationUtils.ccw;
import static com.github.micycle1.surferj.TriangulationUtils.cw;
import static com.github.micycle1.surferj.TriangulationUtils.mod3;

import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.lang3.tuple.Pair;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.math.Vector2D;

import com.github.micycle1.surferj.SurfConstants;
import com.github.micycle1.surferj.collapse.CollapseSpec;
import com.github.micycle1.surferj.collapse.CollapseType;
import com.github.micycle1.surferj.wavefront.Event;

/**
 * Handles the processing of kinetic events for a KineticTriangulation. This
 * class encapsulates the logic for modifying the triangulation state based on
 * different event types (collapses, splits, flips).
 */
public class KineticEventHandler {
	private static final Logger LOGGER = Logger.getLogger(KineticEventHandler.class.getName());

	private final KineticTriangulation kt;

	public KineticEventHandler(KineticTriangulation kineticTriangulation) {
		this.kt = Objects.requireNonNull(kineticTriangulation, "KineticTriangulation cannot be null");
	}

	// --- Main Event DisgetPositionAtcher ---

	public void handleEvent(Event event) {
		final double time = event.getTime();
		if (Double.isNaN(time)) {
			LOGGER.severe("Attempting to handle event with NaN time: " + event);
			// Potentially throw an exception or return early
			return;
		}

		// Basic check against processing events out of order (optional)
		// if (time < kt.lastEventTime - SurfConstants.TIME_TOL) { // Use tolerance
		// LOGGER.warning("Handling event potentially out of order: current_time=" +
		// kt.lastEventTime + ", event_time=" + time);
		// }

		// Increment counters
		kt.eventTypeCounter[CollapseType.UNDEFINED.ordinal()]++; // Assuming UNDEFINED is a valid index
		if (event.getType().ordinal() < kt.eventTypeCounter.length) {
			kt.eventTypeCounter[event.getType().ordinal()]++;
		} else {
			LOGGER.warning("Event type ordinal out of bounds for counter array: " + event.getType());
		}
		updateEventTimingStats(time);

		// Get the triangle *reference* from the triangulation list using the ID from
		// the event
		// This ensures we operate on the potentially updated triangle object in the
		// list,
		// not just the potentially stale one inside the event object itself.
		KineticTriangle triangle = kt.getTriangles().get((int) event.getTriangle().getId()); // Assuming ID fits int and is used as index/lookup key
		if (triangle == null || triangle.isDead() || triangle.isDying()) {
			LOGGER.log(Level.FINE, "Skipping event for dead/dying triangle " + event.getTriangle().getId() + ": " + event);
			return; // Skip event if triangle is already dead/dying
		}
		// Create a new Event object referencing the *live* triangle from the list
		// This ensures handlers work with the correct triangle instance.
		Event liveEvent = new Event(triangle, time); // Recalculate spec? No, use original event data but with live triangle ref
		// It's better to pass the original event's data but ensure triangle reference
		// is live
		// Let's pass the original event, but handlers should use
		// event.getTriangle().getId() to look up live triangle if needed.

		LOGGER.log(Level.FINER, "Handling event: {0}", event);

		switch (event.getType()) {
			case TRIANGLE_COLLAPSE :
				handleTriangleCollapseEvent(event);
				break;
			case CONSTRAINT_COLLAPSE :
				handleConstraintEvent(event);
				break;
			case FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING :
				handleFaceWithInfinitelyFastOpposingVertex(event);
				break;
			case FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED :
				handleFaceWithInfinitelyFastWeightedVertex(event);
				break;
			case SPOKE_COLLAPSE :
				handleSpokeCollapseEvent(event);
				break;
			case SPLIT_OR_FLIP_REFINE :
				handleSplitOrFlipRefineEvent(event);
				break;
			case VERTEX_MOVES_OVER_SPOKE :
				handleVertexMovesOverSpokeEvent(event);
				break;
			case CCW_VERTEX_LEAVES_CH :
				handleCcwVertexLeavesChEvent(event);
				break;
			case INVALID_EVENT :
			case UNDEFINED :
			case NEVER :
				LOGGER.severe("Should not receive event to handle: " + event);
				// Consider throwing an exception
				break;
			default :
				LOGGER.severe("Unexpected event type: " + event);
				// Consider throwing an exception
				break;
		}

		// Process refinement queue after handling the event
		processCheckRefinementQueue(time);

		// Post-event validation (optional, can be expensive)
		// kt.assertValid(event.getComponent(), time);

	}

	// --- Statistics ---
	private void updateEventTimingStats(double now) {
		// Using package-private access or getters/setters for kt fields
		if (Math.abs(now - kt.lastEventTime) > SurfConstants.TIME_TOL) { // Use tolerance
			kt.maxEventsPerTime = Math.max(kt.maxEventsPerTime, kt.eventsPerCurrentEventTime);
			kt.avgEventsPerTimeSum += kt.eventsPerCurrentEventTime;
			kt.avgEventsPerTimeCtr++;
			kt.lastEventTime = now;
			kt.eventsPerCurrentEventTime = 1;
		} else {
			kt.eventsPerCurrentEventTime++;
		}
	}

	// --- Specific Event Handlers (Port C++ logic into these methods) ---

	private void handleTriangleCollapseEvent(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		final double time = event.getTime();
		LOGGER.log(Level.FINEST, "Handling TRIANGLE_COLLAPSE for T{0} @ {1}", new Object[] { t.getId(), time });

		t.markDying();

		for (int i = 0; i < 3; ++i) {
			WavefrontVertex v = t.getVertex(i);
			if (v != null) {
				v.stop(time);
			}
		}
		for (int i = 0; i < 3; ++i) {
			if (t.isConstrained(i)) {
				WavefrontEdge w = t.getWavefront(i);
				if (w != null) {
					w.markDead();
				}
				// update prev/next links for DCEL representation
				WavefrontVertex cwV = t.getVertex(cw(i));
				WavefrontVertex ccwV = t.getVertex(ccw(i));
				if (cwV != null && ccwV != null) {
					cwV.setNextVertex(0, ccwV, false); // Assuming side 0 is CW direction link
				}
			} else {
				KineticTriangle n = t.getNeighbor(i);
				if (n != null) {
					int idxInN = n.indexOfNeighbor(t);
					if (idxInN != -1) {
						t.setNeighborRaw(i, null); // Break link from t
						n.setNeighborRaw(idxInN, null); // Break link from n
						if (!n.isDying()) { // Avoid recursion on already dying triangles
							doSpokeCollapsePart2(n, idxInN, time);
						}
					} else {
						LOGGER.warning("Neighbor inconsistency during TRIANGLE_COLLAPSE for T" + t.getId() + " and neighbor T" + n.getId());
					}
				}
			}
		}
		kt.getQueue().needsDropping(t);
	}

	private void handleConstraintEvent(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		final double time = event.getTime();
		final int edgeIdx = event.getRelevantEdge();
		LOGGER.log(Level.FINEST, "Handling CONSTRAINT_COLLAPSE for T{0} edge {1} @ {2}", new Object[] { t.getId(), edgeIdx, time });

		WavefrontVertex va = t.getVertex(ccw(edgeIdx));
		WavefrontVertex vb = t.getVertex(cw(edgeIdx));
		if (va != null) {
			va.stop(time);
		}
		if (vb != null) {
			vb.stop(time);
		}

		// Update prev/next links
		if (va != null && vb != null) {
			va.setNextVertex(1, vb, false); // Assuming side 1 is CCW direction link
		}

		doConstraintCollapsePart2(t, edgeIdx, time);
	}

	private void handleSpokeCollapseEvent(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		final double time = event.getTime();
		final int edgeIdx = event.getRelevantEdge();
		LOGGER.log(Level.FINEST, "Handling SPOKE_COLLAPSE for T{0} edge {1} @ {2}", new Object[] { t.getId(), edgeIdx, time });

		WavefrontVertex va = t.getVertex(ccw(edgeIdx));
		WavefrontVertex vb = t.getVertex(cw(edgeIdx));
		if (va != null) {
			va.stop(time);
		}
		if (vb != null) {
			vb.stop(time);
		}

		KineticTriangle n = t.getNeighbor(edgeIdx);
		if (n != null) {
			int idxInN = n.indexOfNeighbor(t);
			if (idxInN != -1) {
				t.setNeighborRaw(edgeIdx, null);
				n.setNeighborRaw(idxInN, null);
				if (!n.isDying()) { // Avoid recursion on dying triangles
					doSpokeCollapsePart2(n, idxInN, time);
				}
			} else {
				LOGGER.warning("Neighbor inconsistency during SPOKE_COLLAPSE for T" + t.getId() + " and neighbor T" + n.getId());
			}
		} else {
			LOGGER.warning("Missing neighbor during SPOKE_COLLAPSE for T" + t.getId() + " edge " + edgeIdx);
		}

		// Note: C++ links t and n before calling part2. We do it above.
		// The primary triangle 't' also needs part2 processing.
		doSpokeCollapsePart2(t, edgeIdx, time);

		// C++ comment mentions DCEL linking handled by part2 calls.
	}

	private void handleSplitEvent(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		final double time = event.getTime();
		final int edgeIdx = event.getRelevantEdge();
		LOGGER.log(Level.FINEST, "Handling SPLIT for T{0} edge {1} @ {2}", new Object[] { t.getId(), edgeIdx, time });

		WavefrontVertex v = t.getVertex(edgeIdx);
		WavefrontEdge edge = t.getWavefront(edgeIdx);
		WavefrontVertex va = t.getVertex(ccw(edgeIdx)); // Endpoint of edge
		WavefrontVertex vb = t.getVertex(cw(edgeIdx)); // Endpoint of edge

		if (v == null || edge == null || va == null || vb == null) {
			LOGGER.severe("Null component found during SPLIT event for T" + t.getId());
			return; // Cannot proceed
		}
		final WavefrontEdge edgeA = v.getIncidentEdge(1); // Edge towards va? Check convention
		final WavefrontEdge edgeB = v.getIncidentEdge(0); // Edge towards vb? Check convention
		if (edgeA == null || edgeB == null) {
			LOGGER.severe("Null incident edge found on vertex V" + v.id + " during SPLIT event for T" + t.getId());
			return;
		}

		// TODO: Port C++ assertion checks for orientation

		Coordinate pos = v.getPositionAt(time);
		v.stop(time, pos); // Stop vertex v at calculated position

		// Split the wavefront edge
		Pair<WavefrontEdge, WavefrontEdge> newEdges = edge.split(kt.getWavefrontEdges(), pos); // Assumes split returns Pair/Tuple
		WavefrontEdge nea = newEdges.getLeft(); // New edge part towards original va
		WavefrontEdge neb = newEdges.getRight(); // New edge part towards original vb

		// Create new wavefront vertices
		WavefrontVertex nva = WavefrontVertex.makeVertex(pos, time, nea, edgeA, true, kt.getVertices()); // New vertex using edgeA
		WavefrontVertex nvb = WavefrontVertex.makeVertex(pos, time, edgeB, neb, true, kt.getVertices()); // New vertex using edgeB

		// Update edge endpoint references
		if (nea != null && nva != null) {
			WavefrontVertex neaV0 = nea.getVertex(0); // Should be original va
			if (neaV0 != null) {
				neaV0.setIncidentEdge(1, nea); // Update original va
			}
			nea.setVertexRaw(1, nva); // Set vertex 1 of new edge nea
		}
		if (neb != null && nvb != null) {
			WavefrontVertex nebV1 = neb.getVertex(1); // Should be original vb
			if (nebV1 != null) {
				nebV1.setIncidentEdge(0, neb); // Update original vb
			}
			neb.setVertexRaw(0, nvb); // Set vertex 0 of new edge neb
		}

		// Update vertex links (DCEL structure)
		v.setNextVertex(0, nvb);
		v.setNextVertex(1, nva);
		if (nva != null && nvb != null) {
			nva.linkTailToTail(nvb);
		}

		t.markDying();
		long affectedTriangles = 0;

		// Iterate around vertex v and update vertex pointers to nva/nvb
		// This requires the AroundVertexIterator or equivalent logic
		KineticTriangulation.AroundVertexIterator iter = kt.incidentFacesIterator(t, edgeIdx);
		KineticTriangulation.AroundVertexIterator end = kt.incidentFacesEnd();

		// Iterate CCW (towards edgeA side -> nva)
		KineticTriangulation.AroundVertexIterator iterCCW = new KineticTriangulation.AroundVertexIterator(iter.t(), iter.vInTIdx()); // Copy start
		iterCCW.walkCcw(); // Move one step CCW
		while (!iterCCW.isEnd()) {
			KineticTriangle currentTri = iterCCW.t();
			if (currentTri != null && !currentTri.isDying()) {
				currentTri.setVertex(iterCCW.vInTIdx(), nva);
				modified(currentTri, false);
				affectedTriangles++;
			} else {
				break; // Stop if null or dying
			}
			iterCCW.walkCcw();
		}

		// Iterate CW (towards edgeB side -> nvb)
		KineticTriangulation.AroundVertexIterator iterCW = new KineticTriangulation.AroundVertexIterator(iter.t(), iter.vInTIdx()); // Copy start
		iterCW.walkCw(); // Move one step CW
		while (!iterCW.isEnd()) {
			KineticTriangle currentTri = iterCW.t();
			if (currentTri != null && !currentTri.isDying()) {
				currentTri.setVertex(iterCW.vInTIdx(), nvb);
				modified(currentTri, false);
				affectedTriangles++;
			} else {
				break; // Stop if null or dying
			}
			iterCW.walkCw();
		}

		// Update wavefronts in neighbors across original edge 'edge'
		KineticTriangle na = t.getNeighbor(cw(edgeIdx)); // Neighbor towards vb
		KineticTriangle nb = t.getNeighbor(ccw(edgeIdx)); // Neighbor towards va
		if (na != null) {
			int idxInNa = na.indexOfNeighbor(t);
			if (idxInNa != -1) {
				na.setWavefront(idxInNa, neb); // Should use new edge neb
			}
		}
		if (nb != null) {
			int idxInNb = nb.indexOfNeighbor(t);
			if (idxInNb != -1) {
				nb.setWavefront(idxInNb, nea); // Should use new edge nea
			}
		}

		// Update stats
		kt.maxTrianglesPerSplitEvent = Math.max(kt.maxTrianglesPerSplitEvent, affectedTriangles);
		kt.avgTrianglesPerSplitEventSum += affectedTriangles;
		kt.avgTrianglesPerSplitEventCtr++;

		kt.getQueue().needsDropping(t);
	}

	private void handleSplitOrFlipRefineEvent(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		final double time = event.getTime();
		final int edgeIdx = event.getRelevantEdge();
		LOGGER.log(Level.FINEST, "Handling SPLIT_OR_FLIP_REFINE for T{0} edge {1} @ {2}", new Object[] { t.getId(), edgeIdx, time });

		WavefrontVertex v = t.getVertex(edgeIdx);
		WavefrontVertex va = t.getVertex(ccw(edgeIdx));
		WavefrontVertex vb = t.getVertex(cw(edgeIdx));
		if (v == null || va == null || vb == null) {
			return;
		}

		Coordinate pos = v.getPositionAt(time);

		Coordinate posa = va.getPositionAt(time);
		Coordinate posb = vb.getPositionAt(time);
		// Segment s = new Segment(posa, posb); // Need Segment class/logic

		// Use tolerance for comparisons
		double sqDistA = pos.distanceSq(posa);
		double sqDistB = pos.distanceSq(posb);
		double sqDistAB = posa.distanceSq(posb);

		// Check collinearity and position relative to segment AB
		// Simplified check: If distance(A,V) + distance(V,B) approx equals
		// distance(A,B) -> V is on segment AB
		boolean onSegment;

		if (Orientation.index(pos, posa, posb) == Orientation.COLLINEAR) {
			// If collinear, check bounds
			double dot = Vector2D.create(posb).subtract(Vector2D.create(posa)).dot(Vector2D.create(pos).subtract(Vector2D.create(posa)));
			onSegment = dot >= -SurfConstants.ZERO_TOL_SQ && dot <= sqDistAB + SurfConstants.ZERO_TOL_SQ;
		} else {
			onSegment = false;
		}

		if (!onSegment) {
			// Case 1: Outside segment -> Refine to VERTEX_MOVES_OVER_SPOKE
			LOGGER.log(Level.FINEST, "Refining SPLIT_OR_FLIP to VERTEX_MOVES_OVER_SPOKE for T{0}", t.getId());
			double longestSpokeSq;
			int flipEdge;
			if (sqDistA > sqDistB) { // Closer to B
				longestSpokeSq = sqDistA;
				flipEdge = cw(edgeIdx); // Edge opposite A
			} else { // Closer to A or equidistant
				longestSpokeSq = sqDistB;
				flipEdge = ccw(edgeIdx); // Edge opposite B
			}
			// Refine spec: Create new spec and call refineCollapseSpec
			CollapseSpec refinedSpec = new CollapseSpec(CollapseType.VERTEX_MOVES_OVER_SPOKE, time, t, flipEdge, Math.sqrt(longestSpokeSq), t.getComponent()); // Need
																																								// sqrt
																																								// for
																																								// C++
																																								// secondary
																																								// key?
																																								// Check
																																								// usage.
			t.refineCollapseSpec(refinedSpec); // Update internal spec
			kt.getQueue().needsUpdate(t, true); // Mark for queue update

		} else if (sqDistA < SurfConstants.ZERO_TOL_SQ || pos.distanceSq(posa) < SurfConstants.ZERO_TOL_SQ) {
			// Case 2: At endpoint A -> Refine to SPOKE_COLLAPSE
			LOGGER.log(Level.FINEST, "Refining SPLIT_OR_FLIP to SPOKE_COLLAPSE (at A) for T{0}", t.getId());
			CollapseSpec refinedSpec = new CollapseSpec(CollapseType.SPOKE_COLLAPSE, time, t, cw(edgeIdx), t.getComponent()); // Edge
																																// opposite
																																// A
																																// collapses
			t.refineCollapseSpec(refinedSpec);
			kt.getQueue().needsUpdate(t, true);

		} else if (sqDistB < SurfConstants.ZERO_TOL_SQ || pos.distanceSq(posb) < SurfConstants.ZERO_TOL_SQ) {
			// Case 2: At endpoint B -> Refine to SPOKE_COLLAPSE
			LOGGER.log(Level.FINEST, "Refining SPLIT_OR_FLIP to SPOKE_COLLAPSE (at B) for T{0}", t.getId());
			CollapseSpec refinedSpec = new CollapseSpec(CollapseType.SPOKE_COLLAPSE, time, t, ccw(edgeIdx), t.getComponent()); // Edge
																																// opposite
																																// B
																																// collapses
			t.refineCollapseSpec(refinedSpec);
			kt.getQueue().needsUpdate(t, true);

		} else {
			// Case 3: On interior -> Handle as true SPLIT event
			LOGGER.log(Level.FINEST, "Executing true SPLIT for T{0}", t.getId());
			handleSplitEvent(event);
		}
	}

	private void handleVertexMovesOverSpokeEvent(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		LOGGER.log(Level.FINEST, "Handling VERTEX_MOVES_OVER_SPOKE for T{0}", t.getId());
		doFlipEvent(event.getTime(), t, event.getRelevantEdge());
	}

	private void handleCcwVertexLeavesChEvent(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		LOGGER.log(Level.FINEST, "Handling CCW_VERTEX_LEAVES_CH for T{0}", t.getId());
		int idx = event.getRelevantEdge(); // Finite edge index
		// Flip the *other* finite edge (the one CW from the infinite vertex)
		doFlip(t, cw(idx), event.getTime(), true); // Allow collinear flip
	}

	// --- Complex Handlers for Infinite Speed Cases (Require careful porting) ---

	private void handleFaceWithInfinitelyFastOpposingVertex(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		final double time = event.getTime();
		LOGGER.log(Level.FINEST, "Handling FACE_HAS_INFINITELY_FAST_OPPOSING for T{0} @ {1}", new Object[] { t.getId(), time });

		// ... Port the complex logic from C++ ...
		// Includes: checking constraints, stopping vertices, finding fast vertex index,
		// potentially flipping away spokes using AroundVertexIterator,
		// determining collapse edge, calling doConstraintCollapsePart2 or handling
		// spoke collapse.
		LOGGER.warning("Handler for FACE_HAS_INFINITELY_FAST_OPPOSING not fully ported.");

	}

	private void handleFaceWithInfinitelyFastWeightedVertex(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		final double time = event.getTime();
		LOGGER.log(Level.FINEST, "Handling FACE_HAS_INFINITELY_FAST_WEIGHTED for T{0} @ {1}", new Object[] { t.getId(), time });

		// ... Port the complex logic from C++ ...
		// Includes: finding fastest vertex, finding most_cw_triangle, flipping spokes,
		// finding losing edge, stopping vertices, updating links, calling
		// doConstraintCollapsePart2.
		LOGGER.warning("Handler for FACE_HAS_INFINITELY_FAST_WEIGHTED not fully ported.");
	}

	// --- Helper Methods for Event Handling ---

	/**
	 * Handles the second part of a spoke collapse. The two vertices forming the
	 * spoke (va, vb) have already been stopped. This determines if the third vertex
	 * (v) also collapses or if it forms a new constraint. Called recursively for
	 * neighbors.
	 */
	private void doSpokeCollapsePart2(KineticTriangle t, int edgeIdx, double time) {
		if (t == null || t.isDying()) {
			return; // Skip if already processed/dying
		}

		WavefrontVertex v = t.getVertex(edgeIdx);
		WavefrontVertex va = t.getVertex(ccw(edgeIdx)); // Already stopped
		WavefrontVertex vb = t.getVertex(cw(edgeIdx)); // Already stopped

		if (v == null || va == null || vb == null) {
			LOGGER.severe("Null vertex during doSpokeCollapsePart2 for T" + t.getId());
			return;
		}

		// Check if triangle collapses completely: Does v end up at the same point as
		// va/vb?
		// Need a reliable way to get stop position, using va's stop position assuming
		// they match.
		Coordinate stopPos = va.getPosStop(); // Assuming getPosStop() exists
		if (stopPos == null) {
			LOGGER.severe("Cannot get stop position during doSpokeCollapsePart2 for T" + t.getId());
			// Fallback or error - maybe use getPositionAt(time)? Needs careful check.
			stopPos = va.getPositionAt(time); // Use position at time as fallback
		}

		boolean fullCollapse = false;
		if (!v.isInfinite()) {
			Coordinate vPosAtTime = v.getPositionAt(time);
			// Use tolerance for position check
			fullCollapse = vPosAtTime != null && vPosAtTime.distanceSq(stopPos) < SurfConstants.ZERO_TOL_SQ;
			if (fullCollapse) {
				v.stop(time, stopPos); // Ensure v is stopped at the collapse point
			}
		}

		if (fullCollapse) {
			LOGGER.log(Level.FINEST, "--> Full collapse detected in doSpokeCollapsePart2 for T{0}", t.getId());
			t.markDying();
			// Handle neighbors recursively or update DCEL links
			for (int i = 1; i <= 2; ++i) {
				int currentEdge = mod3(edgeIdx + i);
				KineticTriangle n = t.getNeighbor(currentEdge);
				if (n != null) {
					int idxInN = n.indexOfNeighbor(t);
					if (idxInN != -1) {
						t.setNeighborRaw(currentEdge, null);
						n.setNeighborRaw(idxInN, null);
						if (!n.isDying()) { // Avoid infinite recursion
							doSpokeCollapsePart2(n, idxInN, time); // Recurse
						}
					} else {
						LOGGER.warning("Neighbor inconsistency during full spoke collapse part2 for T" + t.getId() + " and neighbor T" + n.getId());
					}
				} else { // Must be a constraint
					// Update DCEL vertex links for the boundary edge
					WavefrontVertex cwV = t.getVertex(cw(currentEdge));
					WavefrontVertex ccwV = t.getVertex(ccw(currentEdge));
					if (cwV != null && ccwV != null) {
						// Assuming side 1 links CCW neighbour, side 0 links CW
						WavefrontVertex outerV = (cwV == v) ? ccwV : cwV; // The vertex that wasn't v
						WavefrontVertex collapsingV = (outerV == cwV) ? ccwV : cwV; // The vertex that was v
						// Need to link the two outer vertices directly? C++ logic unclear here.
						// C++: t.vertices[ccw(edge)]->set_next_vertex(1, t.vertices[cw(edge)], false);
						// Let's assume we link ccwV -> cwV using side 1 (CCW link)
						ccwV.setNextVertex(1, cwV, false);
					}

					WavefrontEdge w = t.getWavefront(currentEdge);
					if (w != null) {
						w.markDead();
					}
				}
			}
			kt.getQueue().needsDropping(t);

		} else {
			LOGGER.log(Level.FINEST, "--> Spoke turning into constraint in doSpokeCollapsePart2 for T{0}", t.getId());
			// Handle infinite speed case (from C++) - check if v is infinite and both
			// neighbors are constraints
			boolean dropTriangle = false;
			if (v.getInfiniteSpeed() != WavefrontVertex.InfiniteSpeedType.NONE) {
				int i1 = mod3(edgeIdx + 1);
				int i2 = mod3(edgeIdx + 2);
				if (t.getNeighbor(i1) == null && t.getNeighbor(i2) == null) { // Both sides constrained?
					LOGGER.log(Level.FINEST, "--> Dropping triangle T{0} due to infinite speed vertex and constraints", t.getId());
					t.markDying();
					v.stop(time, stopPos); // Stop v at the collapse point
					// Update DCEL links across the constraints
					for (int i = 1; i <= 2; ++i) {
						int currentEdge = mod3(edgeIdx + i);
						WavefrontVertex cwV = t.getVertex(cw(currentEdge));
						WavefrontVertex ccwV = t.getVertex(ccw(currentEdge));
						if (cwV != null && ccwV != null) {
							ccwV.setNextVertex(1, cwV, false); // Link across constraint
						}
						WavefrontEdge w = t.getWavefront(currentEdge);
						if (w != null) {
							w.markDead();
						}
					}
					dropTriangle = true;
					kt.getQueue().needsDropping(t);
				}
			}

			if (!dropTriangle) {
				// Triangle does not collapse completely, handle as constraint formation
				doConstraintCollapsePart2(t, edgeIdx, time);
			}
		}
	}

	/**
	 * Handles the second part of a constraint collapse OR a spoke collapse that
	 * turns into a constraint. Vertex va and vb are stopped at the same point.
	 * Creates a new vertex v, updates neighbors and triangle links.
	 */
	private void doConstraintCollapsePart2(KineticTriangle t, int edgeIdx, double time) {
		if (t == null || t.isDying()) {
			return; // Skip if already dying
		}

		WavefrontVertex va = t.getVertex(ccw(edgeIdx));
		WavefrontVertex vb = t.getVertex(cw(edgeIdx));
		if (va == null || vb == null) {
			LOGGER.severe("Null vertex during doConstraintCollapsePart2 for T" + t.getId());
			return;
		}

		// Stop position should be consistent
		Coordinate pos = va.getPosStop();
		if (pos == null) {
			pos = va.getPositionAt(time); // Fallback if stop pos isn't set yet
		}
		if (pos == null) {
			LOGGER.severe("Cannot determine collapse position for T" + t.getId());
			return; // Cannot proceed
		}

		t.markDying();
		moveConstraintsToNeighbor(t, edgeIdx); // Transfer constraints away from t

		WavefrontEdge edgeA = va.getIncidentEdge(0); // Edge leaving va CW? Check convention
		WavefrontEdge edgeB = vb.getIncidentEdge(1); // Edge leaving vb CCW? Check convention
		if (edgeA == null || edgeB == null) {
			LOGGER.severe("Cannot find incident edges for new vertex in doConstraintCollapsePart2 for T" + t.getId());
			// This might happen if va/vb were part of a triangle collapse handled earlier
			// Mark t for dropping and return
			WavefrontEdge w = t.getWavefront(edgeIdx);
			if (w != null) {
				w.markDead(); // Mark original edge dead if it existed
			}
			kt.getQueue().needsDropping(t);
			return;
		}

		// Create the new vertex
		WavefrontVertex v = WavefrontVertex.makeVertex(pos, time, edgeA, edgeB, false, kt.getVertices());
		if (v == null) {
			LOGGER.severe("Failed to create new vertex in doConstraintCollapsePart2 for T" + t.getId());
			return;
		}

		// Link old vertices to the new one
		va.setNextVertex(0, v); // Link side 0 of va
		vb.setNextVertex(1, v); // Link side 1 of vb

		// Update vertex references in the fan of triangles around the collapse point
		long affectedTriangles = 0;
		// Iterate CCW from va
		KineticTriangulation.AroundVertexIterator iterCCW = kt.incidentFacesIterator(t, ccw(edgeIdx));
		KineticTriangulation.AroundVertexIterator end = kt.incidentFacesEnd();
		iterCCW.walkCcw(); // Start one step away
		while (!iterCCW.isEnd()) {
			KineticTriangle currentTri = iterCCW.t();
			if (currentTri == null || currentTri.isDying()) {
				break; // Stop if we hit null or dying
			}
			currentTri.setVertex(iterCCW.vInTIdx(), v);
			modified(currentTri, false);
			affectedTriangles++;
			iterCCW.walkCcw();
		}
		// Iterate CW from vb
		KineticTriangulation.AroundVertexIterator iterCW = kt.incidentFacesIterator(t, cw(edgeIdx));
		iterCW.walkCw(); // Start one step away
		while (!iterCW.isEnd()) {
			KineticTriangle currentTri = iterCW.t();
			if (currentTri == null || currentTri.isDying()) {
				break; // Stop if we hit null or dying
			}
			currentTri.setVertex(iterCW.vInTIdx(), v);
			modified(currentTri, false);
			affectedTriangles++;
			iterCW.walkCw();
		}

		// Update neighbor links (bypass t)
		KineticTriangle na = t.getNeighbor(cw(edgeIdx)); // Neighbor opposite va
		KineticTriangle nb = t.getNeighbor(ccw(edgeIdx)); // Neighbor opposite vb
		if (na != null) {
			int idxInNa = na.indexOfNeighbor(t);
			if (idxInNa != -1) {
				na.setNeighborRaw(idxInNa, nb); // Link na -> nb
			}
		}
		if (nb != null) {
			int idxInNb = nb.indexOfNeighbor(t);
			if (idxInNb != -1) {
				nb.setNeighborRaw(idxInNb, na); // Link nb -> na
			}
		}

		// Mark original wavefront edge (if it was a constraint) dead
		WavefrontEdge w = t.getWavefront(edgeIdx);
		if (w != null) {
			w.markDead();
		}

		// Update stats
		kt.maxTrianglesPerEdgeEvent = Math.max(kt.maxTrianglesPerEdgeEvent, affectedTriangles);
		kt.avgTrianglesPerEdgeEventSum += affectedTriangles;
		kt.avgTrianglesPerEdgeEventCtr++;

		kt.getQueue().needsDropping(t);
	}

	/**
	 * Transfers constraints from triangle t adjacent to edgeIdx to neighbors.
	 * Assumes t is collapsing/dying.
	 */
	private static void moveConstraintsToNeighbor(KineticTriangle t, int edgeIdx) {
		// C++ Precondition: Not both other edges are constrained
		assert !(t.isConstrained(cw(edgeIdx)) && t.isConstrained(ccw(edgeIdx)));

		for (int i = 1; i <= 2; ++i) { // Check edges CW and CCW from edgeIdx
			int constrainedEdge = mod3(edgeIdx + i);
			if (t.isConstrained(constrainedEdge)) {
				// Find the neighbor *opposite* the constrained edge
				int neighborEdge = mod3(edgeIdx + (3 - i)); // The other non-collapsing edge index
				KineticTriangle n = t.getNeighbor(neighborEdge);
				if (n != null) {
					int nidx = n.indexOfNeighbor(t); // Index in neighbor n corresponding to edge towards t
					if (nidx != -1) {
						// Assert that neighbor n does not already have a constraint opposite t
						// C++: CGAL_assertion(! n->is_constrained(nidx) );
						assert !n.isConstrained(nidx) : "Neighbor already constrained where constraint should be moved";
						n.moveConstraintFrom(nidx, t, constrainedEdge); // Move constraint data
					} else {
						LOGGER.warning("Cannot find index of dying triangle T" + t.getId() + " in neighbor T" + n.getId() + " during moveConstraint");
					}
				} else {
					LOGGER.warning("Cannot find neighbor opposite constrained edge " + constrainedEdge + " of dying triangle T" + t.getId());
				}
			}
		}
	}

	/**
	 * Marks triangle t as modified, schedules for refinement check and queue
	 * update.
	 */
	private void modified(KineticTriangle t, boolean front) {
		if (t == null || t.isDying() || t.isDead()) { // Don't modify dying/dead triangles
			return;
		}

		putOnCheckRefinement(t, front);

		if (kt.getQueue() != null) {
			kt.getQueue().needsUpdate(t, false); // Mark for priority queue update
		}

		// Handle unbounded case: notify neighbor across the finite edge
		if (t.isUnbounded()) {
			int infIdx = t.getInfiniteVertexIndex();
			if (infIdx != -1) {
				// The neighbor sharing the finite edge opposite the infinite vertex also needs
				// update
				KineticTriangle n = t.getNeighbor(infIdx); // Neighbor across finite edge
				if (n != null && !n.isDying() && !n.isDead()) { // Check neighbor state
					n.invalidateCollapseSpec(); // Invalidate its current spec
					if (kt.getQueue() != null) {
						kt.getQueue().needsUpdate(n, false); // Mark neighbor for queue update
					}
				}
			}
		}
	}

	// --- Refinement Methods ---

	/** Puts triangle on refinement queue if not already present. */
	private void putOnCheckRefinement(KineticTriangle t, boolean front) {
		if (t == null) {
			return;
		}
		int tid = (int) t.getId(); // Assuming ID fits int
		// Check bounds for BitSet access
		if (tid < 0) {
			LOGGER.warning("Invalid triangle ID for refinement check: " + tid);
			return;
		}
		// Ensure BitSet is large enough (should happen in initialize or when adding
		// triangles)
		// if (tid >= kt.getTidxInCheckRefinement().size()) { /* resize logic or error
		// */ }

		if (kt.getTidxInCheckRefinement().get(tid)) {
			return; // Already in queue
		}

		kt.getTidxInCheckRefinement().set(tid);
		if (front) {
			kt.getCheckRefinement().addFirst(t);
		} else {
			kt.getCheckRefinement().addLast(t);
		}
	}

	/** Removes and returns the next triangle from the refinement queue. */
	private KineticTriangle checkRefinementPop() {
		if (kt.getCheckRefinement().isEmpty()) {
			return null;
		}
		KineticTriangle t = kt.getCheckRefinement().removeFirst();
		if (t != null) {
			int tid = (int) t.getId();
			if (tid >= 0) { // Check bounds
				kt.getTidxInCheckRefinement().clear(tid);
			}
		}
		return t;
	}

	/** Processes all triangles currently in the refinement queue. */
	private void processCheckRefinementQueue(double time) {
		while (!kt.getCheckRefinement().isEmpty()) {
			KineticTriangle t = checkRefinementPop();
			if (t != null && !t.isDying() && !t.isDead()) { // Check state before refining
				refineTriangulation(t, time);
			}
		}
	}

	/** Performs initial refinement pass after triangulation setup. */
	void refineTriangulationInitial() {
		LOGGER.fine("Performing initial triangulation refinement...");
		for (KineticTriangle t : kt.getTriangles()) {
			if (t != null && !t.isDying() && !t.isDead()) {
				// Put all initial triangles on the queue for checking
				putOnCheckRefinement(t, false);
			}
		}
		processCheckRefinementQueue(SurfConstants.CORE_ZERO); // Refine at time 0
		LOGGER.fine("Initial refinement complete.");
	}

	/**
	 * Checks local Delaunay condition for edge opposite vertex `reflexIdx` in
	 * triangle `t`. If the edge is flippable (based on reflex vertices and
	 * geometry), performs the flip. C++ version has complex logic based on number
	 * of reflex vertices in the quad.
	 */
	private void refineTriangulation(KineticTriangle t, double time) {
		if (t.isUnbounded()) {
			return; // Don't refine unbounded triangles this way
		}

		// Simplified logic based on C++ case 1 (one reflex vertex in triangle t)
		// Find reflex vertex in t (if exactly one)
		int reflexIdx = -1;
		int reflexCount = 0;
		for (int i = 0; i < 3; i++) {
			WavefrontVertex v = t.getVertex(i);
			// Check if v exists and is reflex (not convex or straight)
			if (v != null && !v.isConvexOrStraight()) { // Need isConvexOrStraight(time)
				reflexIdx = i;
				reflexCount++;
			}
		}

		if (reflexCount == 1) {
			// Found exactly one reflex vertex at index reflexIdx
			int edgeToFlip = reflexIdx; // Edge opposite the reflex vertex
			if (t.isConstrained(edgeToFlip)) {
				return; // Cannot flip constraint
			}

			KineticTriangle n = t.getNeighbor(edgeToFlip);
			// Conditions from C++ to *not* flip:
			if (n == null || n.isDying() || n.isDead() || n.isUnbounded()) {
				return;
			}

			int idxInN = n.indexOfNeighbor(t);
			if (idxInN == -1) {
				return; // Should not happen
			}

			WavefrontVertex v = t.getVertex(reflexIdx);
			WavefrontVertex va = t.getVertex(ccw(reflexIdx));
			WavefrontVertex vb = t.getVertex(cw(reflexIdx));
			WavefrontVertex o = n.getVertex(idxInN); // Vertex in n opposite the shared edge
			if (v == null || va == null || vb == null || o == null) {
				return; // Need all vertices
			}

			// C++ has further checks based on 'straight' vertices at corners of quad (va,
			// vb, o)
			// Porting these requires robust is_reflex_or_straight checks and vertex linking
			// Example check (simplified): if va is straight and is endpoint of
			// constraints... skip flip
			// boolean skipFlip = false;
			// if (va.isReflexOrStraight(time) && (t.isConstrained(ccw(reflexIdx)) &&
			// n.isConstrained(cw(idxInN))) ) skipFlip = true;
			// if (vb.isReflexOrStraight(time) && (t.isConstrained(cw(reflexIdx)) &&
			// n.isConstrained(ccw(idxInN))) ) skipFlip = true;
			// etc...

			// Geometric check (inCircle or equivalent for Delaunay)
			// For non-Delaunay refinement, this might be different. C++ doesn't show
			// explicit inCircle.
			// It seems to flip *only* if the quad has exactly one reflex vertex (at t's
			// reflexIdx).
			// Let's assume the C++ logic implies a flip is needed/beneficial here unless
			// skipped by straight checks.
			// We need the vertex type checks ported first.

			// Placeholder: Assume flip should happen if not skipped by
			// constraints/unbounded/straight checks
			boolean shouldFlip = true; // Assume true unless specific checks implemented above say no

			if (shouldFlip) {
				LOGGER.log(Level.FINEST, "Refinement: Flipping edge {0} of T{1} (neighbor T{2})", new Object[] { edgeToFlip, t.getId(), n.getId() });
				doFlip(t, edgeToFlip, time, false); // Perform the flip
			}
		}
		// C++ cases for 2 or 3 reflex vertices are omitted or NOPs in the #ifdef block
	}

	// --- Flip Operations ---

	/** Performs the flip, marks triangles modified. */
	private void doFlip(KineticTriangle t, int edgeIdx, double time, boolean allowCollinear) {
		KineticTriangle n = t.getNeighbor(edgeIdx);
		if (n == null || n.isDying() || n.isDead()) {
			LOGGER.warning("Attempted to flip edge " + edgeIdx + " of T" + t.getId() + " but neighbor is null or dead.");
			return;
		}
		doRawFlip(t, edgeIdx, time, allowCollinear);
		modified(t, true); // Put t on front of refinement queue
		modified(n, false); // Put n on back
	}

	/** Handles the geometry check for a flip event and performs the flip. */
	private void doFlipEvent(double time, KineticTriangle t, int edgeIdx) {
		// C++ performs some squared length calculations - unclear purpose here, maybe
		// assertions?
		// Simply perform the flip based on the event data.
		doFlip(t, edgeIdx, time, false); // Assume non-collinear flip for events? Check C++ usage.
	}

	/** Performs the actual pointer manipulation for a flip. */
	private void doRawFlip(KineticTriangle t, int edgeIdx, double time, boolean allowCollinear) {
		KineticTriangle n = t.getNeighbor(edgeIdx);
		if (n == null) {
			return; // Should have been checked by caller
		}
		int nidx = n.indexOfNeighbor(t);
		if (nidx == -1) {
			return; // Should have been checked
		}

		WavefrontVertex v = t.getVertex(edgeIdx);
		WavefrontVertex v1 = t.getVertex(ccw(edgeIdx));
		WavefrontVertex v2 = t.getVertex(cw(edgeIdx));
		WavefrontVertex o = n.getVertex(nidx);

		// TODO: Port C++ assertion logic using geometry checks at 'time' if needed
		// boolean isUnboundedFlip = (v1 != null && v1.isInfinite()) || (o != null &&
		// o.isInfinite());
		// if (isUnboundedFlip) { ... } else { ... }

		// Perform the combinatorial flip
		t.doRawFlip(edgeIdx); // Assumes KineticTriangle.doRawFlip handles pointers
		// TODO call modified()?
	}

}