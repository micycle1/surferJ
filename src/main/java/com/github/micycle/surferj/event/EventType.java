package com.github.micycle.surferj.event;

public enum EventType {
    NONE,               // No event / Placeholder
    EDGE_COLLAPSE,      // Two vertices of a wavefront edge meet
    TRIANGLE_COLLAPSE,  // Three vertices of a triangle meet simultaneously (or spoke collapse simplified)
    FLIP_EVENT,         // Topological change due to vertex crossing an edge
    // --- More complex events (Deferred) ---
    // SPOKE_COLLAPSE,    // Specific case of triangle collapse where one edge shrinks to zero first
    // SPLIT_EVENT,       // Vertex hits a non-incident wavefront edge
    // SPLIT_OR_FLIP,     // Intermediate event type needing refinement
    // INFINITE_SPEED_OPPOSING, // Vertex between parallel opposing edges
    // INFINITE_SPEED_WEIGHTED // Vertex between parallel same-direction edges with different weights
}
