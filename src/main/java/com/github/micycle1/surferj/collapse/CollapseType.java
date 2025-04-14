package com.github.micycle1.surferj.collapse;

public enum CollapseType {
	UNDEFINED, FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING, // Degenerate: Parallel opposing edges same speed
	TRIANGLE_COLLAPSE, // All 3 vertices/edges meet simultaneously
	CONSTRAINT_COLLAPSE, // A constraint edge shrinks to zero length
	SPOKE_COLLAPSE, // An internal edge shrinks to zero length
	SPLIT_OR_FLIP_REFINE, // Vertex hits constraint's supporting line (needs refinement)
	FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED, // Degenerate: Parallel adjacent edges different speeds
	VERTEX_MOVES_OVER_SPOKE, // Flip event for internal edges
	CCW_VERTEX_LEAVES_CH, // Unbounded triangle event
	INVALID_EVENT, // Calculation error or unexpected state
	NEVER // No collapse predicted in the future
}