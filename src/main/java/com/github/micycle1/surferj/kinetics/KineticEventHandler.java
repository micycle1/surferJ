package com.github.micycle1.surferj.kinetics;

import static com.github.micycle1.surferj.TriangulationUtils.ccw;
import static com.github.micycle1.surferj.TriangulationUtils.cw;
import static com.github.micycle1.surferj.TriangulationUtils.mod3;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import org.apache.commons.lang3.tuple.Pair;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.github.micycle1.surferj.SurfConstants;
import com.github.micycle1.surferj.TriangulationUtils;
import com.github.micycle1.surferj.collapse.CollapseSpec;
import com.github.micycle1.surferj.collapse.CollapseType;
import com.github.micycle1.surferj.kinetics.WavefrontVertex.InfiniteSpeedType;
import com.github.micycle1.surferj.wavefront.Event;

/**
 * Handles the processing of kinetic events for a KineticTriangulation. This
 * class encapsulates the logic for modifying the triangulation state based on
 * different event types (collapses, splits, flips).
 */
public class KineticEventHandler {

	// NOTE this replaces C++ handle_event() method on KineticTriangulation
	// a standalone class for code cleanness

	private static final Logger LOGGER = LoggerFactory.getLogger(KineticEventHandler.class);

	private final KineticTriangulation kt;

	public KineticEventHandler(KineticTriangulation kineticTriangulation) {
		this.kt = Objects.requireNonNull(kineticTriangulation, "KineticTriangulation cannot be null");
	}

	public void handleEvent(Event event) {
		final double time = event.getTime();
//		LOGGER.
		if (Double.isNaN(time)) {
			LOGGER.error("Attempting to handle event with NaN time: " + event);
			return;
		}

		// Re-grab the live triangle up front and replace the event’s pointer
		long tid = event.getTriangle().getId();
		KineticTriangle liveT = kt.getTriangles().get((int) tid);
		if (liveT == null || liveT.isDead() || liveT.isDying()) {
			LOGGER.debug("Skipping stale/dead triangle " + tid + ": " + event);
			return;
		}
		// Rebuild event so all handlers see the live triangle instance
		event = new Event(liveT, time); // NOTE

		// Increment counters exactly as in CGAL
		kt.eventTypeCounter[CollapseType.UNDEFINED.ordinal()]++;
		kt.eventTypeCounter[event.getType().ordinal()]++;

		updateEventTimingStats(time);

		LOGGER.debug("Handling event: " + event);
		switch (event.getType()) {
			case TRIANGLE_COLLAPSE :
				handleTriangleCollapseEvent(event);
				break;
			case CONSTRAINT_COLLAPSE :
				handleConstraintEvent(event);
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
			case FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING :
				handleFaceWithInfinitelyFastOpposingVertex(event);
				break;
			case FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED :
				handleFaceWithInfinitelyFastWeightedVertex(event);
				break;
			default :
				LOGGER.error("Unexpected event type: " + event);
				return;
		}

		// Process any triangles queued for local Delaunay‐refinement
		processCheckRefinementQueue(time);

		// NOTE (Optional) expensive debug‐only validity check
		kt.assertValid(liveT.getComponent(), time);
	}

	// --- Statistics ---
	private void updateEventTimingStats(double now) {
		if (now != kt.lastEventTime) {
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
		LOGGER.debug("Handling TRIANGLE_COLLAPSE for T{} @ {}", new Object[] { t.getId(), time });

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
						LOGGER.warn("Neighbor inconsistency during TRIANGLE_COLLAPSE for T" + t.getId() + " and neighbor T" + n.getId());
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
		LOGGER.debug("Handling CONSTRAINT_COLLAPSE for T{} edge {} @ {}", new Object[] { t.getId(), edgeIdx, time });

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

		LOGGER.trace("Incident edges: {}. {}", va.getIncidentEdge(0), va.getIncidentEdge(1));
		LOGGER.trace("Incident edges: {}. {}", vb.getIncidentEdge(0), vb.getIncidentEdge(1));

		doConstraintCollapsePart2(t, edgeIdx, time);
	}

	private void handleSpokeCollapseEvent(Event event) {
		KineticTriangle t = event.getTriangle();
		int edgeIdx = event.getRelevantEdge();
		double time = event.getTime();

		LOGGER.debug("SPOKE_COLLAPSE T" + t.getId() + " edge " + edgeIdx + " @ " + time);
		// stop the two spoke endpoints
		WavefrontVertex va = t.getVertex(ccw(edgeIdx));
		WavefrontVertex vb = t.getVertex(cw(edgeIdx));
		va.stop(time);
		vb.stop(time);

		// break link to neighbor
		KineticTriangle n = t.getNeighbor(edgeIdx);
		if (n != null) {
			int idxInN = n.indexOfNeighbor(t);
			t.setNeighborRaw(edgeIdx, null);
			n.setNeighborRaw(idxInN, null);
		}

		// 2a) first handle collapse on the primary triangle
		doSpokeCollapsePart2(t, edgeIdx, time);

		// 2b) now handle collapse on the neighbor
		if (n != null && !n.isDying()) {
			int idxInN = n.indexOfNeighbor(t); // should still be correct
			doSpokeCollapsePart2(n, idxInN, time);
		}
	}

	private void handleSplitEvent(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		final double time = event.getTime();
		final int edgeIdx = event.getRelevantEdge();
		LOGGER.debug("Handling SPLIT for T{} edge {} @ {}", new Object[] { t.getId(), edgeIdx, time });

		WavefrontVertex v = t.getVertex(edgeIdx);
		WavefrontEdge edge = t.getWavefront(edgeIdx);
		WavefrontVertex va = t.getVertex(ccw(edgeIdx)); // Endpoint of edge
		WavefrontVertex vb = t.getVertex(cw(edgeIdx)); // Endpoint of edge

		if (v == null || edge == null || va == null || vb == null) {
			LOGGER.error("Null component found during SPLIT event for T" + t.getId());
			return; // Cannot proceed
		}
		final WavefrontEdge edgeA = v.getIncidentEdge(1); // Edge towards va? Check convention
		final WavefrontEdge edgeB = v.getIncidentEdge(0); // Edge towards vb? Check convention
		if (edgeA == null || edgeB == null) {
			LOGGER.error("Null incident edge found on vertex V" + v.id + " during SPLIT event for T" + t.getId());
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
		AroundVertexIterator iter = kt.incidentFacesIterator(t, edgeIdx);
		AroundVertexIterator end = kt.incidentFacesEnd();

		// Iterate CCW (towards edgeA side -> nva)
		AroundVertexIterator iterCCW = new AroundVertexIterator(iter.t(), iter.vInTIdx()); // Copy start
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
		AroundVertexIterator iterCW = new AroundVertexIterator(iter.t(), iter.vInTIdx()); // Copy start
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
		int tid = (int) event.getTriangle().getId();
		KineticTriangle t = kt.getTriangles().get(tid);
//		t = event.getTriangle();

		if (t == null || t.isDying()) {
			return;
		}
		double time = event.getTime();
		int ei = event.getRelevantEdge();

		LOGGER.debug("Handling SPLIT_OR_FLIP_REFINE for T{} edge {} @ {}", new Object[] { t.getId(), ei, time });

		WavefrontVertex v = t.getVertex(ei);
		WavefrontVertex va = t.getVertex(ccw(ei));
		WavefrontVertex vb = t.getVertex(cw(ei));
		if (v == null || va == null || vb == null) {
			return;
		}

		Coordinate pos = v.getPositionAt(time);
		Coordinate posa = va.getPositionAt(time);
		Coordinate posb = vb.getPositionAt(time);

		// Use supporting line to test collinearity+between:

		boolean onSegment = TriangulationUtils.pointLiesOnSegment(posa, posb, pos);

		if (!onSegment) {
			// refine into VERTEX_MOVES_OVER_SPOKE
			LOGGER.debug("Refining SPLIT_OR_FLIP to VERTEX_MOVES_OVER_SPOKE for T{}", t.getId());
			double d2a = pos.distanceSq(posa), d2b = pos.distanceSq(posb);
			double longestSpoke = Math.sqrt(Math.max(d2a, d2b));
			int flipEdge = (d2a > d2b) ? cw(ei) : ccw(ei);
			CollapseSpec cs = new CollapseSpec(CollapseType.VERTEX_MOVES_OVER_SPOKE, time, t, flipEdge, longestSpoke, t.getComponent());
			LOGGER.debug("tID: {}, spec: {}", t.getId(), t.getCollapseSpec(time).toString());
			t.refineCollapseSpec(cs);
			kt.getQueue().needsUpdate(t, true);
		} else if (pos.distanceSq(posa) < SurfConstants.ZERO_TOL_SQ) {
			// refine into SPOKE_COLLAPSE at A
			LOGGER.debug("Refining SPLIT_OR_FLIP to SPOKE_COLLAPSE (at A) for T{}", t.getId());
			CollapseSpec cs = new CollapseSpec(CollapseType.SPOKE_COLLAPSE, time, t, cw(ei), t.getComponent());
			t.refineCollapseSpec(cs);
			kt.getQueue().needsUpdate(t, true);
		} else if (pos.distanceSq(posb) < SurfConstants.ZERO_TOL_SQ) {
			// refine into SPOKE_COLLAPSE at B
			LOGGER.debug("Refining SPLIT_OR_FLIP to SPOKE_COLLAPSE (at B) for T{}", t.getId());
			CollapseSpec cs = new CollapseSpec(CollapseType.SPOKE_COLLAPSE, time, t, ccw(ei), t.getComponent());
			t.refineCollapseSpec(cs);
			kt.getQueue().needsUpdate(t, true);
		} else {
			// true split
			LOGGER.debug("Executing true SPLIT for T{}", t.getId());
			handleSplitEvent(event);
		}
	}

	private void handleVertexMovesOverSpokeEvent(Event event) {
		KineticTriangle t = kt.getTriangles().get((int) event.getTriangle().getId()); // Get live triangle
		if (t == null || t.isDying()) {
			return; // Check again
		}
		LOGGER.debug("Handling VERTEX_MOVES_OVER_SPOKE for T{}", t.getId());
		doFlipEvent(event.getTime(), t, event.getRelevantEdge());
	}

	private void handleCcwVertexLeavesChEvent(Event event) {
		KineticTriangle t = event.getTriangle();
		int idx = event.getRelevantEdge();
		double time = event.getTime();

		LOGGER.debug("CCW_VERTEX_LEAVES_CH T" + t.getId() + " edge " + idx + " @ " + time);
		// CGAL calls do_flip(…, allow_collinear=false)
		doFlip(t, cw(idx), time, false);
	}

	private void handleFaceWithInfinitelyFastOpposingVertex(Event event) {
		// NOTE complex original impl -- maybe check this port?
		KineticTriangle t = event.getTriangle();
		double time = event.getTime();

		// 1) Count constraints & infinite‐speed vertices
		int numConstrained = 0, numFast = 0;
		for (int i = 0; i < 3; ++i) {
			if (t.isConstrained(i)) {
				++numConstrained;
			}
			if (t.getVertex(i).getInfiniteSpeed() != InfiniteSpeedType.NONE) {
				++numFast;
			}
		}

		// 2) CASE A: all three edges are constrained → full triangle collapse
		if (numConstrained == 3) {
			// mark dying
			t.markDying();

			// stop the unique *finite*‐speed vertex, record its stop‐position p
			Coordinate p = null;
			for (int i = 0; i < 3; ++i) {
				WavefrontVertex vi = t.getVertex(i);
				if (vi.getInfiniteSpeed() == InfiniteSpeedType.NONE) {
					vi.stop(time);
					p = vi.getPosStop();
					break;
				}
			}
			// stop the two infinite‐speed vertices at the same point p
			for (int i = 0; i < 3; ++i) {
				WavefrontVertex vi = t.getVertex(i);
				if (vi.getInfiniteSpeed() != InfiniteSpeedType.NONE) {
					vi.stop(time, p);
				}
			}
			// kill all three wavefront edges and stitch DCEL‐links around the boundary
			for (int i = 0; i < 3; ++i) {
				WavefrontEdge w = t.getWavefront(i);
				if (w != null) {
					w.markDead();
				}
				// link cw(i) → ccw(i) on side 0 (CW link)
				WavefrontVertex cw = t.getVertex(cw(i));
				WavefrontVertex ccw = t.getVertex(ccw(i));
				if (cw != null && ccw != null) {
					cw.setNextVertex(0, ccw, false);
				}
			}
			// finally remove triangle from the queue
			kt.getQueue().needsDropping(t);

			kt.assertValid(t.getComponent(), time);
			return;
		}

		// 3) CASE B: exactly one vertex v has infinite speed opposing a constrained
		// face
		int fastIdx = t.getInfiniteSpeedOpposingVertexIndex();
		WavefrontVertex v = t.getVertex(fastIdx);

		// 3a) walk around v to the most‐CW triangle, flipping spokes until you hit a
		// constraint
		List<KineticTriangle> flippedNeighbors = new ArrayList<>();
		AroundVertexIterator avi = kt.incidentFacesIterator(t, fastIdx);
		AroundVertexIterator mostCw = avi.mostCw();
		while (true) {
			int e = ccw(mostCw.vInTIdx());
			if (mostCw.t().isConstrained(e)) {
				break;
			}
			// do a raw flip (allow collinear) on that spoke
			KineticTriangle flipT = mostCw.t();
			doRawFlip(flipT, e, time, true);
			// remember its neighbor for re‐validation
			KineticTriangle neighbor = flipT.getNeighbor(e);
			if (neighbor != null) {
				flippedNeighbors.add(neighbor);
			}
			// advance to the new most‐CW after flip
			mostCw = new AroundVertexIterator(flipT, mostCw.vInTIdx()).mostCw();
		}

		// 3b) Now at the collapsing triangle t2
		KineticTriangle t2 = mostCw.t();
		int vidx2 = mostCw.vInTIdx();
		WavefrontVertex vCw = t2.getVertex(cw(vidx2));
		WavefrontVertex vCcw = t2.getVertex(ccw(vidx2));

		// 3c) Decide: spoke‐collapse vs. constraint‐collapse
		int collapseEdge = -1;
		boolean doSpokeCollapse = false;
		// if both cw and ccw are infinite → spoke collapse
		if (vCw.getInfiniteSpeed() != InfiniteSpeedType.NONE && vCcw.getInfiniteSpeed() != InfiniteSpeedType.NONE) {
			doSpokeCollapse = true;
		}
		// if exactly one neighbor infinite → collapse the opposite constraint
		else if (vCw.getInfiniteSpeed() != InfiniteSpeedType.NONE) {
			collapseEdge = cw(vidx2);
		} else if (vCcw.getInfiniteSpeed() != InfiniteSpeedType.NONE) {
			collapseEdge = ccw(vidx2);
		}
		// else both finite → pick by shortest distance, tie→spoke
		else {
			Coordinate posV = v.getPositionAt(time);
			double dCW = posV.distanceSq(vCw.getPositionAt(time));
			double dCCW = posV.distanceSq(vCcw.getPositionAt(time));
			if (dCW < dCCW) {
				collapseEdge = ccw(vidx2);
			} else if (dCCW < dCW) {
				collapseEdge = cw(vidx2);
			} else {
				doSpokeCollapse = true;
			}
		}

		// 4) Perform the chosen collapse
		if (doSpokeCollapse) {
			// exactly as in handleSpokeCollapseEvent but on t2:
			t2.markDying();
			vCw.stop(time);
			vCcw.stop(time);
			Coordinate stopPoint = vCw.getPosStop();
			v.stop(time, stopPoint);

			// kill the two wavefronts
			WavefrontEdge w1 = t2.getWavefront(cw(vidx2));
			if (w1 != null) {
				w1.markDead();
			}
			WavefrontEdge w2 = t2.getWavefront(ccw(vidx2));
			if (w2 != null) {
				w2.markDead();
			}

			// unlink t2 ↔ neighbor
			KineticTriangle n2 = t2.getNeighbor(vidx2);
			int n2idx = (n2 == null ? -1 : n2.indexOfNeighbor(t2));
			t2.setNeighborRaw(vidx2, null);
			if (n2 != null && n2idx >= 0) {
				n2.setNeighborRaw(n2idx, null);
			}

			// stitch v → its two ends
			v.setNextVertex(0, vCw, false);
			v.setNextVertex(1, vCcw, false);

			// recurse onto neighbor
			if (n2 != null && !n2.isDying()) {
				doSpokeCollapsePart2(n2, n2idx, time);
			}
			kt.getQueue().needsDropping(t2);
		} else {
			// constraint‐collapse on edge=collapseEdge
			WavefrontVertex o = (vCw.getInfiniteSpeed() != InfiniteSpeedType.NONE ? vCcw : vCw);
			o.stop(time);
			Coordinate oPos = o.getPosStop();
			v.stop(time, oPos);

			// link v → o
			int which = (o == vCcw ? 1 : 0);
			v.setNextVertex(which, o, false);

			doConstraintCollapsePart2(t2, collapseEdge, time);
		}

		// 5) re‐validate any triangles we flipped above
		for (KineticTriangle nt : flippedNeighbors) {
			if (!nt.isDead()) {
				modified(nt, false);
			}
		}

		kt.assertValid(t.getComponent(), time);
	}

	private void handleFaceWithInfinitelyFastWeightedVertex(Event event) {
		KineticTriangle t = event.getTriangle();
		double time = event.getTime();

		// 1) Identify which endpoint has WEIGHTED infinite speed
		// and which is the “opposite” wavefront.
		WavefrontVertex vw, vo; // vw = weighted‐infinite, vo = other
		boolean weightedIsCcw;
		int edgeIdx = event.getRelevantEdge(); // edge opposite the infinite‐weighted

		WavefrontVertex vCcw = t.getVertex(ccw(edgeIdx));
		WavefrontVertex vCw = t.getVertex(cw(edgeIdx));

		if (vCcw.getInfiniteSpeed() == InfiniteSpeedType.WEIGHTED) {
			vw = vCcw;
			vo = vCw;
			weightedIsCcw = true;
		} else {
			vw = vCw;
			vo = vCcw;
			weightedIsCcw = false;
		}

		int vIdxInT = weightedIsCcw ? ccw(edgeIdx) : cw(edgeIdx);
		AroundVertexIterator mw = kt.incidentFacesIterator(t, vIdxInT).mostCw();

		// if apex (edgeIdx‐vertex) is infinite speed, do a single flip …
		WavefrontVertex apex = t.getVertex(edgeIdx);
		if (apex.getInfiniteSpeed() != InfiniteSpeedType.NONE) {
			int nidx = ccw(mw.vInTIdx());
			KineticTriangle nbr = mw.t().getNeighbor(nidx);
			doRawFlip(mw.t(), nidx, time, true);
			if (nbr != null) {
				modified(nbr, false);
			}
		} else {
			// general “walk‐and‐flip” loop
			while (true) {
				int nidx = ccw(mw.vInTIdx());
				if (mw.t().isConstrained(nidx)) {
					break;
				}
				KineticTriangle tri = mw.t();
				Coordinate pa = tri.getVertex(ccw(mw.vInTIdx())).getPositionAt(time);
				Coordinate pb = tri.getVertex(cw(mw.vInTIdx())).getPositionAt(time);

				if (tri != t) {
					// decide whether to advance or flip
					KineticTriangle n2 = tri.getNeighbor(cw(mw.vInTIdx()));
					Coordinate pc = n2.getVertex(n2.indexOfNeighbor(tri)).getPositionAt(time);
					if (Orientation.index(pc, pa, pb) != Orientation.RIGHT) {
						// advance one step CW around v
						mw.walkCw(); // <<<<< was mw++ in C++
						continue;
					}
				}

				// perform the flip
				doRawFlip(tri, nidx, time, true);
				KineticTriangle post = tri.getNeighbor(nidx);
				if (post != null) {
					modified(post, false);
					// stay on this triangle and re‐check
				}
			}
		}

		// 2) Identify the “losing” edge and other vertex
		// losingEdge is the non‐weighted side
		WavefrontEdge losingEdge = vw.getIncidentEdge(weightedIsCcw ? 0 : 1);
		WavefrontVertex o = losingEdge.getVertex(weightedIsCcw ? 1 : 0);

		// 3) stop them in correct order
		if (o.getInfiniteSpeed() != InfiniteSpeedType.NONE) {
			o.stop(time, o.getPosStop()); // snap back to start‐pos
		} else {
			o.stop(time);
		}
		// weighted vertex stops at o’s position
		vw.stop(time, o.getPosStop());

		// link weighted→other in DCEL
		int side = weightedIsCcw ? 1 : 0;
		vw.setNextVertex(side, o, false);

		// 4) now do a constraint collapse across the losing edge
		KineticTriangle collapseTri = mw.t();
		int collapseIdx = collapseTri.indexOfNeighbor(t);
		// in C++ they call: do_constraint_collapse_part2(*most_cw_triangle,
		// most_cw_triangle->index(losing_edge), time);
		doConstraintCollapsePart2(collapseTri, collapseIdx, time);

		kt.assertValid(t.getComponent(), time);
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
			LOGGER.error("Null vertex during doSpokeCollapsePart2 for T" + t.getId());
			return;
		}

		// Check if triangle collapses completely: Does v end up at the same point as
		// va/vb?
		// Need a reliable way to get stop position, using va's stop position assuming
		// they match.
		Coordinate stopPos = va.getPosStop(); // Assuming getPosStop() exists
		if (stopPos == null) {
			LOGGER.error("Cannot get stop position during doSpokeCollapsePart2 for T" + t.getId());
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
			LOGGER.debug("--> Full collapse detected in doSpokeCollapsePart2 for T{}", t.getId());
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
						LOGGER.warn("Neighbor inconsistency during full spoke collapse part2 for T" + t.getId() + " and neighbor T" + n.getId());
					}
				} else { // Must be a constraint
					// Update DCEL vertex links for the boundary edge
					WavefrontVertex cwV = t.getVertex(cw(currentEdge));
					WavefrontVertex ccwV = t.getVertex(ccw(currentEdge));
					if (cwV != null && ccwV != null) {
						// Assuming side 1 links CCW neighbour, side 0 links CW
						WavefrontVertex outerV = (cwV == v) ? ccwV : cwV; // The vertex that wasn't v
						WavefrontVertex collapsingV = (outerV == cwV) ? ccwV : cwV; // The vertex that was v
						// TODO Need to link the two outer vertices directly? C++ logic unclear here.
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
			LOGGER.debug("--> Spoke turning into constraint in doSpokeCollapsePart2 for T{}", t.getId());
			// Handle infinite speed case (from C++) - check if v is infinite and both
			// neighbors are constraints
			boolean dropTriangle = false;
			if (v.getInfiniteSpeed() != WavefrontVertex.InfiniteSpeedType.NONE) {
				int i1 = mod3(edgeIdx + 1);
				int i2 = mod3(edgeIdx + 2);
				if (t.getNeighbor(i1) == null && t.getNeighbor(i2) == null) { // Both sides constrained?
					LOGGER.debug("--> Dropping triangle T{} due to infinite speed vertex and constraints", t.getId());
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
			LOGGER.error("Null vertex during doConstraintCollapsePart2 for T" + t.getId());
			return;
		}

		// Stop position should be consistent
		Coordinate pos = va.getPosStop();
		if (pos == null) {
			pos = va.getPositionAt(time); // Fallback if stop pos isn't set yet
		}
		if (pos == null) {
			LOGGER.error("Cannot determine collapse position for T" + t.getId());
			return; // Cannot proceed
		}

		t.markDying();
		moveConstraintsToNeighbor(t, edgeIdx); // Transfer constraints away from t

		WavefrontEdge edgeA = va.getIncidentEdge(0); // Edge leaving va CW? Check convention
		WavefrontEdge edgeB = vb.getIncidentEdge(1); // Edge leaving vb CCW? Check convention
		if (edgeA == null || edgeB == null) {
			LOGGER.error("Cannot find incident edges for new vertex in doConstraintCollapsePart2 for T" + t.getId());
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
			LOGGER.error("Failed to create new vertex in doConstraintCollapsePart2 for T" + t.getId());
			return;
		}

		if (edgeA != null) {
			// edgeA was incoming edge for va (V36).
			// The new vertex v replaces va as the target endpoint (vertex 1) for edgeA.
			edgeA.setVertexRaw(1, v); // Update WE48's vertex 1 to be V50
		}
		if (edgeB != null) {
			// edgeB was outgoing edge for vb (V37).
			// The new vertex v replaces vb as the source endpoint (vertex 0) for edgeB.
			edgeB.setVertexRaw(0, v); // Update WE11's vertex 0 to be V50
		}

		// Link old vertices to the new one
		va.setNextVertex(0, v); // Link side 0 of va
		vb.setNextVertex(1, v); // Link side 1 of vb

		// Update vertex references in the fan of triangles around the collapse point
		long affectedTriangles = 0;
		// Iterate CCW from va
		  // --- Iterate CCW from va ---
	    // Start iterator at the dying triangle t, positioned at va's index
	    AroundVertexIterator iterCCW = kt.incidentFacesIterator(t, ccw(edgeIdx));
	    iterCCW.walkCcw(); // *** STEP AWAY FIRST *** to the first living neighbor CCW

	    while (!iterCCW.isEnd()) {
	        KineticTriangle currentTri = iterCCW.t();
	        if (currentTri == null) { // Should not happen if topology is good
	            LOGGER.error("Iterator CCW found null triangle unexpectedly.");
	            break;
	        }
	         if (currentTri.isDying()) { // Should not encounter dying triangle after first step
	             LOGGER.warn("Iterator CCW encountered dying triangle {} after initial step.", currentTri.getId());
	             // Decide whether to break or try to continue cautiously
	             break; // Breaking is safer
	         }

	        // We are now on a living triangle adjacent to the new vertex v
	        currentTri.setVertex(iterCCW.vInTIdx(), v); // Set vertex to v
	        modified(currentTri, false); // Invalidate/notify
	        affectedTriangles++;

	        // Check before walking: can we walk further?
	        KineticTriangle next = iterCCW.nextTriangleCcw(); // Peek ahead
	        if (next == null) { // Stop if we hit a constraint boundary
	            break;
	        }

	        iterCCW.walkCcw(); // Walk to the next one
	    }

	    // --- Iterate CW from vb ---
	    // Start iterator at the dying triangle t, positioned at vb's index
	    AroundVertexIterator iterCW = kt.incidentFacesIterator(t, cw(edgeIdx));
	    iterCW.walkCw(); // *** STEP AWAY FIRST *** to the first living neighbor CW

	    while (!iterCW.isEnd()) {
	       KineticTriangle currentTri = iterCW.t();
	        if (currentTri == null) {
	            LOGGER.error("Iterator CW found null triangle unexpectedly.");
	            break;
	        }
	         if (currentTri.isDying()) {
	             LOGGER.warn("Iterator CW encountered dying triangle {} after initial step.", currentTri.getId());
	             break;
	         }

	        // We are now on a living triangle adjacent to the new vertex v
	        currentTri.setVertex(iterCW.vInTIdx(), v); // Set vertex to v
	        modified(currentTri, false); // Invalidate/notify
	        affectedTriangles++;

	        // Check before walking
	        KineticTriangle next = iterCW.nextTriangleCw(); // Peek ahead
	        if (next == null) { // Stop if we hit a constraint boundary
	            break;
	        }

	        iterCW.walkCw(); // Walk to the next one
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
						LOGGER.warn("Cannot find index of dying triangle T" + t.getId() + " in neighbor T" + n.getId() + " during moveConstraint");
					}
				} else {
					LOGGER.warn("Cannot find neighbor opposite constrained edge " + constrainedEdge + " of dying triangle T" + t.getId());
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
			LOGGER.warn("Invalid triangle ID for refinement check: " + tid);
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
		LOGGER.debug("Performing initial triangulation refinement...");
		for (KineticTriangle t : kt.getTriangles()) {
			if (t != null && !t.isDying() && !t.isDead()) {
				// Put all initial triangles on the queue for checking
				putOnCheckRefinement(t, false);
			}
		}
		processCheckRefinementQueue(SurfConstants.CORE_ZERO); // Refine at time 0
		LOGGER.debug("Initial refinement complete.");
	}

	static boolean refine = false;

	/**
	 * Checks local Delaunay condition for edge opposite vertex `reflexIdx` in
	 * triangle `t`. If the edge is flippable (based on reflex vertices and
	 * geometry), performs the flip. C++ version has complex logic based on number
	 * of reflex vertices in the quad.
	 */
	private void refineTriangulation(KineticTriangle t, double time) {
		if (refine) {
			// NOTE under #ifdef REFINE_TRIANGULATION
			if (t.isUnbounded()) {
				return;
			}

			// count reflex‐or‐straight
			int reflexCount = 0, reflexIdx = -1;
			for (int i = 0; i < 3; i++) {
				WavefrontVertex v = t.getVertex(i);
				if (v != null && !v.isConvexOrStraight()) {
					reflexIdx = i;
					reflexCount++;
				}
			}
			if (reflexCount != 1) {
				return;
			}

			int edgeToFlip = reflexIdx;
			if (t.isConstrained(edgeToFlip)) {
				return;
			}

			KineticTriangle n = t.getNeighbor(edgeToFlip);
			if (n == null || n.isDying() || n.isDead() || n.isUnbounded()) {
				return;
			}
			int idxInN = n.indexOfNeighbor(t);
			if (idxInN < 0) {
				return;
			}

			// now port CGAL’s “straight‐corner” skips:
			WavefrontVertex va = t.getVertex(ccw(reflexIdx));
			WavefrontVertex vb = t.getVertex(cw(reflexIdx));
			WavefrontVertex o = n.getVertex(idxInN);

			if ((va.isConvexOrStraight() && t.isConstrained(ccw(reflexIdx)) && n.isConstrained(cw(idxInN)))
					|| (vb.isConvexOrStraight() && t.isConstrained(cw(reflexIdx)) && n.isConstrained(ccw(idxInN)))
					|| (o.isConvexOrStraight() && n.isConstrained(idxInN))) {
				return;
			}

			// (Optional) you could add an inCircle test here if you want exact-Delaunay
			// flips.

			doFlip(t, edgeToFlip, time, /* allowCollinear= */false);
		}
	}

	// --- Flip Operations ---

	/** Performs the flip, marks triangles modified. */
	private void doFlip(KineticTriangle t, int edgeIdx, double time, boolean allowCollinear) {
		KineticTriangle n = t.getNeighbor(edgeIdx);
		if (n == null || n.isDying() || n.isDead()) {
			LOGGER.warn("Attempted to flip edge " + edgeIdx + " of T" + t.getId() + " but neighbor is null or dead.");
			return;
		}
		doRawFlip(t, edgeIdx, time, allowCollinear);
		modified(t, true); // Put t on front of refinement queue
		modified(n, false); // Put n on back
	}

	private void doFlipEvent(double time, KineticTriangle t, int edgeIdx) {
		// replicate CGAL's squared‐distance calculations (for debug/assert only)
		WavefrontVertex v0 = t.getVertex(edgeIdx);
		WavefrontVertex v1 = t.getVertex(ccw(edgeIdx));
		WavefrontVertex v2 = t.getVertex(cw(edgeIdx));
		Coordinate p0 = v0.getPositionAt(time), p1 = v1.getPositionAt(time), p2 = v2.getPositionAt(time);

		double l01 = p0.distanceSq(p1), l12 = p1.distanceSq(p2), l20 = p2.distanceSq(p0);
		// (You can log or assert if you like)

		// now perform the real flip
		doFlip(t, edgeIdx, time, /* allowCollinear= */false);
	}

	private void doRawFlip(KineticTriangle t, int edgeIdx, double time, boolean allowCollinear) {
		KineticTriangle n = t.getNeighbor(edgeIdx);
		if (n == null) {
			return;
		}
		int nidx = n.indexOfNeighbor(t);

		// fetch the four involved vertices
		WavefrontVertex v = t.getVertex(edgeIdx);
		WavefrontVertex v1 = t.getVertex(ccw(edgeIdx));
		WavefrontVertex v2 = t.getVertex(cw(edgeIdx));
		WavefrontVertex o = n.getVertex(nidx);

		// CGAL‐style orientation checks:
		if (v1.isInfinite() || o.isInfinite()) {
			// unbounded flip: require collinearity
			assert Orientation.index(v.getPositionAt(time), o.getPositionAt(time), v2.getPositionAt(time)) == Orientation.COLLINEAR;
		} else {
			// bounded: require non‐right‐turn
			assert Orientation.index(v1.getPositionAt(time), v2.getPositionAt(time), v.getPositionAt(time)) != Orientation.RIGHT;
			if (!allowCollinear) {
				assert Orientation.index(v1.getPositionAt(time), v2.getPositionAt(time), o.getPositionAt(time)) != Orientation.RIGHT;
			}
		}

		// finally do the pointer swap
		t.doRawFlip(edgeIdx);
	}

}