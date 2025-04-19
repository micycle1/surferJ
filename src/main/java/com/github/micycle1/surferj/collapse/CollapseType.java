package com.github.micycle1.surferj.collapse;

public enum CollapseType {
	/**
	 * Undefined collapse type.
	 */
	UNDEFINED,

	/**
	 * This triangle has a vertex which is between parallel, opposing wavefront
	 * elements that have crashed into each other and their intersection is now a
	 * line segment.
	 */
	FACE_HAS_INFINITELY_FAST_VERTEX_OPPOSING,

	/**
	 * All 3 vertices/edges meet simultaneously.
	 */
	TRIANGLE_COLLAPSE,

	/**
	 * A constraint edge shrinks to zero length.
	 */
	CONSTRAINT_COLLAPSE,

	/**
	 * Two non-incident vertices become incident, splitting the wavefront here.
	 */
	SPOKE_COLLAPSE,

	/**
	 * Vertex moves onto supporting line of constraint, can refine event type when
	 * it comes to it.
	 */
	SPLIT_OR_FLIP_REFINE,

	/**
	 * This triangle has a vertex which is between parallel adjacent wavefront
	 * elements that have different weights but move in the same direction.
	 */
	FACE_HAS_INFINITELY_FAST_VERTEX_WEIGHTED,

	/**
	 * Vertex moves into spoke (triangulation edge interior), flip event.
	 */
	VERTEX_MOVES_OVER_SPOKE,

	/**
	 * The ccw vertex of the infinite vertex in an unbounded triangle leaves the
	 * convex hull of the wavefront polygon.
	 */
	CCW_VERTEX_LEAVES_CH,

	/**
	 * The triangle will collapse at this time, but we should never see this as
	 * prior events should have rebuilt the triangulation in some way. If this is
	 * the next event, something went wrong.
	 */
	INVALID_EVENT,

	/**
	 * No collapse predicted in the future.
	 */
	NEVER
}