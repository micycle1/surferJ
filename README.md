# surferj
surferJ

### Core Concepts & Mapping:
1. Input: C++ uses BGL/GraphML -> BasicInput. Java will use JTS Polygon.
2. Initial Triangulation: C++ uses CGAL Constrained_Delaunay_triangulation_2. Java will use TinFour IncrementalTin.
3. Kinetic Data Structures: The core concepts (WavefrontVertex, WavefrontEdge, KineticTriangle, KineticTriangulation) need direct Java counterparts.
4. Geometry Primitives: C++ uses CGAL Point_2, Vector_2, Segment_2, Line_2, NT (Number Type). Java will use JTS Coordinate, Vector2D, LineSegment, potentially custom Line class, and double for coordinates/time initially. TinFour primitives (Vertex, Edge, Triangle) will be used for the initial triangulation structure.
5. Events & Queue: C++ uses CollapseSpec, Event, EventQueue (custom heap). Java will use CollapseEvent (enum + class), PriorityQueue.
6. Output: C++ uses SkeletonDCEL. Java will initially use a simpler structure, like List<LineSegment> for skeleton arcs.
7. Numerics: C++ uses CGAL's exact predicates/constructions. Java will start with double and JTS/TinFour robustness features. Precision issues might need later attention (e.g., using BigDecimal or specialized libraries if double proves insufficient).


### Simplified Scope:
- Single, simple input polygon (no holes, non-self-intersecting).
- Uniform edge weights (weight = 1.0).
- No complex beveling (degree-1 vertices, reflex vertices needing >1 split). We'll handle basic convex/reflex vertices.
- Focus on Edge Collapse and Triangle Collapse events first, then basic Flip events. Split events are complex and deferred.
- Output is a list of skeleton line segments.

### How to Use & Next Steps:

1. Dependencies: Add JTS (e.g., org.locationtech.jts:jts-core) and TinFour (org.tinfour:tinfour) to your project (e.g., Maven/Gradle).
2. Main Class: Create a main class that reads a Polygon (e.g., from WKT using WKTReader), creates StraightSkeletonGenerator, calls generate(), and processes the SkeletonOutput.
3. Testing: Create test cases with simple polygons (square, rectangle, L-shape, convex polygon) and verify the output visually or programmatically if possible.
4. Refinement & Debugging: This is a complex algorithm. Expect significant debugging.
    - Visualization: Add debug code to output the state of the kinetic triangulation at different times (e.g., write triangle/vertex positions to WKT or a simple graphics format).
    - Logging: Use a proper logging framework (Log4j, SLF4j) to trace event handling and calculations.
    - Velocity Calculation: The GeometryUtil.calculateVelocity is crucial and needs careful testing and potential refinement, especially regarding edge orientations and handling collinear cases correctly.
    - Connectivity: The establishConnectivity and event handling (especially handleEdgeCollapse topology updates) are complex points of failure. Ensure neighbor and vertex pointers are updated correctly.
    - Flip Logic: The handleFlipEvent needs a full implementation of the topological updates.
    - Precision: Monitor for issues arising from double precision. Compare results with known outputs or the C++ version if possible.
  
5. Extensibility:
    - Weights: Modify KineticEdge and velocity calculations to handle non-uniform weights.
    - Beveling: Implement the logic from the C++ code for splitting vertices (degree-1, reflex) during initialization (KineticTriangulation.create_bevels). This involves adding new triangles and edges dynamically.
    - Split Events: Implement the calculation and handling for split events (vertex hitting non-incident edge). This requires geometric tests (point-line distance over time).
    - Output: Implement a more robust DCEL-like structure for the output if needed.
    - Holes/Multi-Polygons: Adapt the initial triangulation and event handling for more complex inputs.