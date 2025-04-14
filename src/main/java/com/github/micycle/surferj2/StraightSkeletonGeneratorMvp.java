package com.github.micycle.surferj2;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.triangulate.quadedge.QuadEdge;
import org.locationtech.jts.triangulate.quadedge.QuadEdgeSubdivision;
import java.util.*;

/**
 * Generates a Straight Skeleton for a JTS Polygon (MVP Implementation).
 *
 * <p>This is a simplified port of concepts found in triangulation-based
 * straight skeleton algorithms like 'surfer2'.</p>
 *
 * <p><b>MVP Limitations:</b></p>
 * <ul>
 *   <li>Uses standard double-precision floating-point arithmetic (potential robustness issues).</li>
 *   <li>Primarily handles edge collapse events (two vertices meeting).</li>
 *   <li>May not correctly handle triangle collapses, split events (reflex vertices), or complex topologies.</li>
 *   <li>Likely works best for convex or simple concave polygons.</li>
 *   <li>Output skeleton might be incomplete or incorrect for complex inputs.</li>
 *   <li>No beveling or advanced offset features.</li>
 * </ul>
 *
 * <p>All classes are defined within this single file for simplicity.</p>
 */
public class StraightSkeletonGeneratorMvp {

    private static final GeometryFactory gf = new GeometryFactory();
    private static final double EPSILON = 1e-8; // Tolerance for floating-point comparisons

    // -------------------------------------------------------------------------
    // Data Structures
    // -------------------------------------------------------------------------

    /** Represents a moving vertex in the wavefront. (Corresponds to C++ WavefrontVertex) */
    private static class LAVertex { // LAV = Live Arc Vertex
        private static int nextId = 0;
        final int id = nextId++;

        final Coordinate posZero; // Position at time t=0 (initial intersection)
        final Vector2D velocity; // Velocity vector (dx/dt, dy/dt)
        final double timeStart;   // Time this vertex was created

        Coordinate posStart; // Actual position when created (might be != posZero if timeStart > 0)

        // Incident wavefront edges (conceptually left/ccw and right/cw)
        // Edge 0 is the "left" edge when looking from the vertex outwards
        // Edge 1 is the "right" edge when looking from the vertex outwards
        LAEdge edge0; // Left / CCW edge
        LAEdge edge1; // Right / CW edge

        boolean stopped = false;
        double timeStop = Double.POSITIVE_INFINITY;
        Coordinate posStop = null;

        // --- Skeleton Structure Tracking ---
        // Pointers to vertices earlier/later in time along the skeleton arcs
        // We track two sides corresponding to the two faces the skeleton arc separates
        LAVertex prevVertex0 = null; // Vertex this one grew *from* on side 0
        LAVertex nextVertex0 = null; // Vertex that grew *from* this one on side 0
        LAVertex prevVertex1 = null; // Vertex this one grew *from* on side 1
        LAVertex nextVertex1 = null; // Vertex that grew *from* this one on side 1

        LAVertex(Coordinate posZero, Vector2D velocity, double timeStart, Coordinate posStart) {
            this.posZero = posZero != null ? posZero.copy() : null; // Can be null for initial infinite vertex concept
            this.velocity = velocity;
            this.timeStart = timeStart;
            this.posStart = posStart != null ? posStart.copy() : null;
        }

        Coordinate p_at(double time) {
            if (stopped && time >= timeStop) {
                return posStop;
            }
            if (posZero == null || velocity == null) { // Should only be for the conceptual infinite vertex placeholder
                return null;
            }
            // p(t) = p(0) + v*t
            return new Coordinate(
                    posZero.x + velocity.getX() * time,
                    posZero.y + velocity.getY() * time);
        }

        void stop(double time, Coordinate stopPos) {
            if (!stopped) {
                this.stopped = true;
                this.timeStop = time;
                this.posStop = stopPos.copy();

                // Check if degenerate (stopped immediately at creation)
                // if (Math.abs(timeStop - timeStart) < EPSILON && posStart.distance(posStop) < EPSILON) {
                // Degeneracy check might need refinement
                // }
            }
        }

        // Link this vertex (which just stopped) to the newly created one
        void setNext(LAVertex next, int side) {
             if (side == 0) {
                 if (this.nextVertex0 != null && this.nextVertex0 != next)
                    System.err.printf("Warning: Overwriting nextVertex0 for LAV %d%n", id);
                 this.nextVertex0 = next;
                 if (next != null) {
                    if (next.prevVertex0 != null && next.prevVertex0 != this)
                        System.err.printf("Warning: Overwriting prevVertex0 for LAV %d%n", next.id);
                    next.prevVertex0 = this;
                 }
             } else {
                  if (this.nextVertex1 != null && this.nextVertex1 != next)
                    System.err.printf("Warning: Overwriting nextVertex1 for LAV %d%n", id);
                 this.nextVertex1 = next;
                  if (next != null) {
                     if (next.prevVertex1 != null && next.prevVertex1 != this)
                        System.err.printf("Warning: Overwriting prevVertex1 for LAV %d%n", next.id);
                    next.prevVertex1 = this;
                  }
             }
        }


        @Override
        public String toString() {
            return "LAV" + id + "[" + (posStart != null ? formatCoord(posStart) : "null") +
                   (stopped ? " -> " + formatCoord(posStop) + "@" + String.format("%.2f", timeStop) : "") + "]";
        }

        @Override public int hashCode() { return id; }
        @Override public boolean equals(Object obj) {
            if (this == obj) return true;
            if (obj == null || getClass() != obj.getClass()) return false;
            LAVertex other = (LAVertex) obj;
            return id == other.id;
        }
    }

    /** Represents a moving edge in the wavefront. (Corresponds to C++ WavefrontEdge) */
    private static class LAEdge { // LA = Live Arc Edge
        private static int nextId = 0;
        final int id = nextId++;

        final LineSegment initialSegment; // The original polygon edge
        final SupportingLine supportingLine;
        LAVertex vStart; // Start vertex (moving)
        LAVertex vEnd; // End vertex (moving)
        KineticTriangle triangle; // The triangle "behind" this edge
        boolean isProcessed = false; // If edge event has been handled

        LAEdge(LineSegment initialSegment, SupportingLine supportingLine, KineticTriangle triangle) {
            this.initialSegment = initialSegment; // Can be null for internal edges
            this.supportingLine = supportingLine;
            this.triangle = triangle;
        }

        // Calculate the time when this edge collapses (its vertices meet)
        // Returns Double.POSITIVE_INFINITY if they don't meet or move apart/parallel
        EdgeCollapseResult calculateCollapseTime(double currentTime) {
            if (vStart == null || vEnd == null || vStart.stopped || vEnd.stopped) {
                return new EdgeCollapseResult(Double.POSITIVE_INFINITY, null);
            }

            Vector2D vA = vStart.velocity;
            Vector2D vB = vEnd.velocity;
            Coordinate pA0 = vStart.posZero;
            Coordinate pB0 = vEnd.posZero;

            // Relative velocity and position
            Vector2D vRel = vB.subtract(vA); // Velocity of B relative to A
            Coordinate pRel0 = new Coordinate(pB0.x - pA0.x, pB0.y - pA0.y); // Position of B relative to A at t=0

            // We are looking for time t > currentTime such that || pA(t) - pB(t) ||^2 = 0
            // || (pA0 + vA*t) - (pB0 + vB*t) ||^2 = 0
            // || (pA0 - pB0) + (vA - vB)*t ||^2 = 0
            // || -pRel0 - vRel*t ||^2 = 0
            // || pRel0 + vRel*t ||^2 = 0
            // Let p = pRel0, v = vRel. We want ||p + v*t||^2 = 0
            // (px + vx*t)^2 + (py + vy*t)^2 = 0
            // (px^2 + 2*px*vx*t + vx^2*t^2) + (py^2 + 2*py*vy*t + vy^2*t^2) = 0
            // (vx^2 + vy^2)*t^2 + 2*(px*vx + py*vy)*t + (px^2 + py^2) = 0
            // This is a quadratic equation At^2 + Bt + C = 0
            // A = ||vRel||^2 = vRel.dot(vRel)
            // B = 2 * pRel0.dot(vRel)
            // C = ||pRel0||^2 = pRel0.dot(pRel0)

            double A = vRel.dot(vRel);
            double B = 2 * (pRel0.x * vRel.getX() + pRel0.y * vRel.getY());
            double C = pRel0.x * pRel0.x + pRel0.y * pRel0.y;

            if (Math.abs(A) < EPSILON) {
                // If A is near zero, relative velocity is near zero (parallel movement).
                // If B is also near zero, they are essentially stationary relative to each other.
                // If C is near zero, they started coincident -> collapse time = start time
                // If C is non-zero, they started apart and remain so -> never collapse.
                 if (Math.abs(B) < EPSILON && Math.abs(C) < EPSILON) {
                     // Started coincident and moving parallel/stationary - effectively collapse at start
                     // This case is tricky, depends on how start times are handled. For simplicity, treat as never.
                     return new EdgeCollapseResult(Double.POSITIVE_INFINITY, null);
                 } else if (Math.abs(B) < EPSILON) { // A=0, B=0, C!=0 -> Never
                     return new EdgeCollapseResult(Double.POSITIVE_INFINITY, null);
                 } else {
                     // Linear equation Bt + C = 0 -> t = -C / B
                     double t = -C / B;
                      if (t >= currentTime - EPSILON) {
                          Coordinate collapsePoint = vStart.p_at(t); // Or vEnd.p_at(t)
                          return new EdgeCollapseResult(t, collapsePoint);
                      } else {
                           return new EdgeCollapseResult(Double.POSITIVE_INFINITY, null); // Collapse was in the past
                      }
                 }

            }

            // Quadratic formula: t = [-B +/- sqrt(B^2 - 4AC)] / 2A
            double discriminant = B * B - 4 * A * C;

            if (discriminant < -EPSILON) {
                // No real roots -> vertices never meet exactly (they pass by)
                return new EdgeCollapseResult(Double.POSITIVE_INFINITY, null);
            }

            // Clamp small negative discriminant to 0 due to potential float errors
             if (discriminant < 0) discriminant = 0;

            double sqrtDiscriminant = Math.sqrt(discriminant);
            double t1 = (-B - sqrtDiscriminant) / (2 * A);
            double t2 = (-B + sqrtDiscriminant) / (2 * A);

            // We want the smallest time t >= currentTime
            double collapseTime = Double.POSITIVE_INFINITY;
            if (t1 >= currentTime - EPSILON && t1 < collapseTime) {
                collapseTime = t1;
            }
            if (t2 >= currentTime - EPSILON && t2 < collapseTime) {
                collapseTime = t2;
            }


            if (collapseTime == Double.POSITIVE_INFINITY) {
                 return new EdgeCollapseResult(Double.POSITIVE_INFINITY, null); // Collapsed in past or never
            } else {
                 // Calculate collapse point
                 Coordinate collapsePoint = vStart.p_at(collapseTime); // Or vEnd.p_at(collapseTime)
                 return new EdgeCollapseResult(collapseTime, collapsePoint);
            }
        }

        @Override
        public String toString() {
            return "LAE" + id + "(" + vStart + " -> " + vEnd + ")";
        }
        @Override public int hashCode() { return id; }
         @Override public boolean equals(Object obj) {
            if (this == obj) return true;
            if (obj == null || getClass() != obj.getClass()) return false;
            LAEdge other = (LAEdge) obj;
            return id == other.id;
        }
    }

     /** Stores result of edge collapse calculation */
    private static class EdgeCollapseResult {
        final double time;
        final Coordinate point;

        EdgeCollapseResult(double time, Coordinate point) {
            this.time = time;
            this.point = point;
        }
    }

    /** Represents the supporting line of a wavefront edge. (Corresponds to C++ WavefrontSupportingLine) */
    private static class SupportingLine {
        final LineSegment segment; // Original defining segment (p0->p1)
        Vector2D normal;    // Unit normal vector pointing "outward" (offset direction)
        final double weight;      // Speed of offset (usually 1.0)

        SupportingLine(LineSegment segment, double weight) {
            this.segment = segment;
            this.weight = weight;

            // Calculate outward normal (assuming CCW polygon vertices)
            // Vector tangent = p1 - p0
            double dx = segment.p1.x - segment.p0.x;
            double dy = segment.p1.y - segment.p0.y;
            // Normal is (-dy, dx) for CCW outward normal
            Vector2D n = new Vector2D(-dy, dx);
            this.normal = n.normalize(); // Ensure unit vector
            if (this.normal == null) {
                // Handle degenerate segment if necessary
                System.err.println("Warning: Degenerate segment used for SupportingLine");
                 this.normal = new Vector2D(0, 0); // Or throw error
            }
        }
    }

    /** Represents a triangle in the kinetic triangulation. (Corresponds to C++ KineticTriangle) */
    private static class KineticTriangle {
         private static int nextId = 0;
         final int id = nextId++;
        // Vertices (indices into a main vertex list, or direct references)
        // These are *not* the LAVertex moving points, but the fixed triangle corners
        // For MVP, we might not need explicit vertices if we rely on QuadEdge info
        // Coordinate v0, v1, v2; // Or references to QuadEdge vertices

        // Edges (references to LAEdge for constrained edges, or null/internal marker)
        LAEdge edge0, edge1, edge2; // Edges opposite v0, v1, v2
        KineticTriangle neighbor0, neighbor1, neighbor2; // Neighbors opposite v0, v1, v2

        boolean isProcessed = false; // If triangle event has been handled

        // Store reference to the original QuadEdge face/triangle if needed
        Object originalTriangleRef; // Could be QuadEdge triangle or JTS Geometry

        @Override public String toString() { return "KT" + id; }
        @Override public int hashCode() { return id; }
         @Override public boolean equals(Object obj) {
             if (this == obj) return true;
             if (obj == null || getClass() != obj.getClass()) return false;
             KineticTriangle other = (KineticTriangle) obj;
             return id == other.id;
         }

         LAEdge getConstrainedEdge() {
             if (edge0 != null && edge0.initialSegment != null) return edge0;
             if (edge1 != null && edge1.initialSegment != null) return edge1;
             if (edge2 != null && edge2.initialSegment != null) return edge2;
             return null; // Should not happen for initial triangles attached to polygon boundary
         }
    }

    /** Event occurring at a specific time involving a triangle or edge. (Corresponds to C++ Event) */
    private static class Event implements Comparable<Event> {
        final double time;
        final EventType type;
        final KineticTriangle triangle; // Triangle primarily involved (can be null for edge events in some designs)
        final LAEdge edge; // Edge involved (relevant for edge events)
        final Coordinate collapsePoint; // Location where event occurs

        // Constructor for Edge Collapse Event
        Event(double time, KineticTriangle triangle, LAEdge edge, Coordinate collapsePoint) {
            this.time = time;
            this.type = EventType.EDGE_COLLAPSE;
            this.triangle = triangle;
            this.edge = edge;
            this.collapsePoint = collapsePoint;
        }

        // Add constructors for other event types if/when implemented

        @Override
        public int compareTo(Event other) {
            int timeCompare = Double.compare(this.time, other.time);
            if (timeCompare != 0) {
                return timeCompare;
            }
            // Add tie-breaking rules if necessary (e.g., based on event type prio)
            return Integer.compare(this.hashCode(), other.hashCode()); // Basic tie-breaker
        }

        @Override
        public String toString() {
            return String.format("Event[%.2f, %s, %s]", time, type, (edge != null ? edge : triangle));
        }
    }

    private enum EventType {
        EDGE_COLLAPSE, // Two LAVertices meet, collapsing an LAEdge
        // TRIANGLE_COLLAPSE, // Three LAVertices meet simultaneously (often follows edge collapses)
        // SPLIT_EVENT // A LAVertex hits an opposite LAEdge (requires reflex vertex handling)
        // Add other types as needed
    }

    /** Simple structure to hold skeleton edges before final geometry creation */
    private static class SkeletonArc {
        final LAVertex startVertex;
        final LAVertex endVertex; // Can be null if ray goes to infinity

        SkeletonArc(LAVertex start, LAVertex end) {
            this.startVertex = start;
            this.endVertex = end;
        }
    }


    // -------------------------------------------------------------------------
    // Helper Classes and Methods
    // -------------------------------------------------------------------------

    private static String formatCoord(Coordinate c) {
        if (c == null) return "null";
        return String.format("(%.2f, %.2f)", c.x, c.y);
    }

    // Basic 2D Vector utility
    private static class Vector2D {
        final double x, y;

        Vector2D(double x, double y) { this.x = x; this.y = y; }
        Vector2D(Coordinate start, Coordinate end) { this.x = end.x - start.x; this.y = end.y - start.y; }

        double getX() { return x; }
        double getY() { return y; }
        double length() { return Math.sqrt(x * x + y * y); }

        Vector2D normalize() {
            double len = length();
            if (len < EPSILON) return null; // Avoid division by zero
            return new Vector2D(x / len, y / len);
        }

        Vector2D subtract(Vector2D other) { return new Vector2D(x - other.x, y - other.y); }
        Vector2D add(Vector2D other) { return new Vector2D(x + other.x, y + other.y); }
        Vector2D multiply(double scalar) { return new Vector2D(x * scalar, y * scalar); }
        double dot(Vector2D other) { return x * other.x + y * other.y; }
        double dot(Coordinate other) { return x * other.x + y * other.y; }

        // Perpendicular vector (-y, x)
        Vector2D perpendicular() { return new Vector2D(-y, x); }

         @Override public String toString() { return String.format("<%.2f, %.2f>", x, y); }
    }

    /** Calculates the intersection point of two lines defined by segments. */
    private static Coordinate lineIntersection(LineSegment s1, LineSegment s2) {
        return s1.lineIntersection(s2);
         // Robust intersection might be needed for production code
    }

    /** Calculate velocity of a vertex formed by intersection of offsets of lines a and b. */
    private static Vector2D calculateVelocity(SupportingLine lineA, SupportingLine lineB, Coordinate intersectionP0) {
        LineSegment segA = lineA.segment;
        LineSegment segB = lineB.segment;
        Vector2D normalA = lineA.normal.multiply(lineA.weight);
        Vector2D normalB = lineB.normal.multiply(lineB.weight);

        Coordinate pA0_offset = new Coordinate(segA.p0.x + normalA.x, segA.p0.y + normalA.y);
        Coordinate pA1_offset = new Coordinate(segA.p1.x + normalA.x, segA.p1.y + normalA.y);
        Coordinate pB0_offset = new Coordinate(segB.p0.x + normalB.x, segB.p0.y + normalB.y);
        Coordinate pB1_offset = new Coordinate(segB.p1.x + normalB.x, segB.p1.y + normalB.y);

        LineSegment offsetSegA = new LineSegment(pA0_offset, pA1_offset);
        LineSegment offsetSegB = new LineSegment(pB0_offset, pB1_offset);

        Coordinate intersectionP1 = lineIntersection(offsetSegA, offsetSegB);

        if (intersectionP0 == null || intersectionP1 == null) {
            Vector2D normA = lineA.normal;
            Vector2D normB = lineB.normal;
            if (normA == null || normB == null) {
                System.err.println("Warning: Degenerate supporting line in velocity calculation.");
                return new Vector2D(0, 0);
            }
            double cross = normA.x * normB.y - normA.y * normB.x;
            if (Math.abs(cross) < EPSILON) {
                double dot = normA.dot(normB);
                if (dot > 0 && Math.abs(lineA.weight - lineB.weight) < EPSILON) {
                    return normA.multiply(lineA.weight);
                }
                System.out.println("Parallel lines detected. Using normal velocity.");
                return normA.multiply(lineA.weight);
            }
            System.err.println("Warning: Invalid intersection in velocity calculation.");
            return new Vector2D(0, 0);
        }

        return new Vector2D(intersectionP0, intersectionP1);
    }
    
    private Vector2D calculateBisectorVelocity(LAEdge edgeA, LAEdge edgeB, Coordinate collapsePoint, double currentTime) {
        if (edgeA == null || edgeB == null || edgeA.supportingLine == null || edgeB.supportingLine == null) {
            System.err.println("Warning: Null edges or supporting lines in calculateBisectorVelocity.");
            return new Vector2D(0, 0);
        }

        // Get the vertices’ positions at currentTime
        LAVertex vA = edgeA.vStart != null && !edgeA.vStart.stopped ? edgeA.vStart : edgeA.vEnd;
        LAVertex vB = edgeB.vEnd != null && !edgeB.vEnd.stopped ? edgeB.vEnd : edgeB.vStart;
        if (vA == null || vB == null) {
            System.err.println("Warning: Could not determine active vertices for edges LAE" + edgeA.id + ", LAE" + edgeB.id);
            return new Vector2D(0, 0);
        }

        // Compute edge directions using vertex positions at collapse time
        Coordinate pA = vA.p_at(currentTime);
        Coordinate pB = vB.p_at(currentTime);
        if (pA == null || pB == null) {
            System.err.println("Warning: Null vertex positions at T=" + currentTime);
            return new Vector2D(0, 0);
        }

        // Approximate edge directions (since vertices are collapsing, use adjacent vertices or velocities)
        Vector2D dirA, dirB;
        if (edgeA.vStart != vA && edgeA.vStart != null) {
            Coordinate pOtherA = edgeA.vStart.p_at(currentTime);
            dirA = pOtherA != null ? new Vector2D(pA, pOtherA).normalize() : vA.velocity.normalize();
        } else {
            Coordinate pOtherA = edgeA.vEnd != null ? edgeA.vEnd.p_at(currentTime) : null;
            dirA = pOtherA != null ? new Vector2D(pOtherA, pA).normalize() : vA.velocity.normalize();
        }
        if (edgeB.vEnd != vB && edgeB.vEnd != null) {
            Coordinate pOtherB = edgeB.vEnd.p_at(currentTime);
            dirB = pOtherB != null ? new Vector2D(pOtherB, pB).normalize() : vB.velocity.normalize();
        } else {
            Coordinate pOtherB = edgeB.vStart != null ? edgeB.vStart.p_at(currentTime) : null;
            dirB = pOtherB != null ? new Vector2D(pB, pOtherB).normalize() : vB.velocity.normalize();
        }

        if (dirA == null || dirB == null) {
            System.err.println("Warning: Could not compute edge directions for LAE" + edgeA.id + ", LAE" + edgeB.id);
            return new Vector2D(0, 0);
        }

        // Compute the angle bisector (average of unit vectors for interior angle)
        Vector2D bisector = dirA.add(dirB).normalize();
        if (bisector == null) {
            System.err.println("Warning: Zero bisector vector. Falling back to normal average.");
            bisector = new Vector2D(0, 0);
        }

        // Ensure bisector points inward (toward shrinking polygon interior)
        // Compute the normal at collapse point (average of supporting line normals)
        Vector2D avgNormal = edgeA.supportingLine.normal.add(edgeB.supportingLine.normal).normalize();
        if (avgNormal != null && bisector.dot(avgNormal) < 0) {
            bisector = bisector.multiply(-1); // Flip to point inward
        }

        // Scale by offset speed (assume weight = 1.0 for simplicity)
        double speed = edgeA.supportingLine.weight; // Should equal edgeB’s weight in MVP
        return bisector.multiply(speed);
    }

    // -------------------------------------------------------------------------
    // Main Algorithm Logic
    // -------------------------------------------------------------------------

    private final Polygon inputPolygon;
    private final List<LAVertex> liveVertices = new ArrayList<>();
    private final List<LAEdge> liveEdges = new ArrayList<>();
    private final List<KineticTriangle> kineticTriangles = new ArrayList<>();
    private final PriorityQueue<Event> eventQueue = new PriorityQueue<>();
    private final List<SkeletonArc> skeletonArcs = new ArrayList<>();
    private double currentTime = 0.0;

    // Map JTS vertices/edges/triangles to internal structures
    private final Map<Coordinate, LAVertex> vertexMap = new HashMap<>();
    private final Map<LineSegment, LAEdge> edgeMap = new HashMap<>(); // Might need better key (edge ID?)
    private final Map<Object, KineticTriangle> triangleMap = new HashMap<>(); // Key is QuadEdge triangle/face object


    public StraightSkeletonGeneratorMvp(Polygon polygon) {
        if (polygon == null || polygon.isEmpty() || !polygon.isValid() || !(polygon.getExteriorRing() instanceof LinearRing)) {
            throw new IllegalArgumentException("Input must be a valid, non-empty Polygon with a LinearRing shell.");
        }
         // Ensure CCW orientation for exterior, CW for interior (JTS standard is typically CCW for exterior)
         // For simplicity, we assume valid input orientation for now.
         // polygon.normalize(); // May change coordinate order
        this.inputPolygon = polygon;
    }

    /**
     * Generates the straight skeleton.
     * @return A Geometry (likely MultiLineString) representing the skeleton, or null on failure.
     */
    public Geometry generate() {
        try {
            initialize();
            propagateWavefront();
            return buildSkeletonGeometry();
        } catch (Exception e) {
            System.err.println("Error during skeleton generation: " + e.getMessage());
            e.printStackTrace();
            return null;
        }
    }

    private void initialize() {
        System.out.println("Initializing Straight Skeleton...");
        vertexMap.clear();
        edgeMap.clear();
        triangleMap.clear();
        liveVertices.clear();
        liveEdges.clear();
        kineticTriangles.clear();
        eventQueue.clear();
        LAVertex.nextId = 0;
        LAEdge.nextId = 0;
        KineticTriangle.nextId = 0;

        // Get polygon exterior ring coordinates
        LineString shell = inputPolygon.getExteriorRing();
        Coordinate[] coords = shell.getCoordinates(); // Includes closing point

        // Ensure CCW orientation
//        if (!CGAlgorithms.isCCW(coords)) {
//            System.err.println("Warning: Input polygon is not CCW. Reversing coordinates.");
//            coords = Arrays.copyOf(coords, coords.length);
//            Collections.reverse(Arrays.asList(coords).subList(0, coords.length - 1));
//        }

        // Create LAVertices and LAEdges for boundary
        LAVertex[] vertices = new LAVertex[coords.length - 1];
        LAEdge[] edges = new LAEdge[coords.length - 1];

        for (int i = 0; i < coords.length - 1; i++) {
            Coordinate p0 = coords[i];
            Coordinate p1 = coords[(i + 1) % (coords.length - 1)];
            LineSegment segment = new LineSegment(p0, p1);

            // Previous segment for velocity calculation
            Coordinate p_prev = coords[i == 0 ? coords.length - 2 : i - 1];
            LineSegment prevSegment = new LineSegment(p_prev, p0);

            // Create SupportingLines
            SupportingLine line = new SupportingLine(segment, 1.0);
            SupportingLine prevLine = new SupportingLine(prevSegment, 1.0);

            // Create or get LAVertex at p0
            LAVertex currentVertex = vertexMap.get(p0);
            if (currentVertex == null) {
                Vector2D velocity = calculateVelocity(prevLine, line, p0);
                currentVertex = new LAVertex(p0, velocity, 0.0, p0);
                vertexMap.put(p0, currentVertex);
                liveVertices.add(currentVertex);
                System.out.println("Created LAV" + currentVertex.id + " at " + formatCoord(p0) + " vel=" + velocity);
            }
            vertices[i] = currentVertex;

            // Create KineticTriangle (placeholder for boundary edge)
            KineticTriangle kt = new KineticTriangle();
            kineticTriangles.add(kt);
            // triangleMap.put(segment, kt); // Use segment as key for now

            // Create LAEdge
            LAEdge edge = new LAEdge(segment, line, kt);
            edge.vStart = currentVertex;
            edges[i] = edge;
            liveEdges.add(edge);
            edgeMap.put(segment, edge);
            kt.edge0 = edge; // Assign to triangle
            System.out.println("Created LAE" + edge.id + " for segment " + formatCoord(p0) + "->" + formatCoord(p1));
        }

        // Assign vEnd and link edges to vertices
        for (int i = 0; i < edges.length; i++) {
            edges[i].vEnd = vertices[(i + 1) % vertices.length];
            vertices[i].edge1 = edges[i]; // Right edge
            vertices[i].edge0 = edges[(i - 1 + edges.length) % edges.length]; // Left edge
        }

        // Verify connectivity
        for (LAEdge edge : liveEdges) {
            if (edge.vStart == null || edge.vEnd == null) {
                System.err.println("Warning: Boundary edge LAE" + edge.id + " missing vertex assignments.");
            }
        }

        // Create initial edge collapse events
        for (LAEdge edge : liveEdges) {
            if (edge.isProcessed || edge.vStart == null || edge.vEnd == null) continue;
            EdgeCollapseResult result = edge.calculateCollapseTime(currentTime);
            if (result.time < Double.POSITIVE_INFINITY && result.point != null) {
                System.out.println("Initial Event for " + edge + ": Collapse at T=" + String.format("%.2f", result.time) + " P=" + formatCoord(result.point));
                eventQueue.add(new Event(result.time, edge.triangle, edge, result.point));
            }
        }
        System.out.println("Initialization complete. " + eventQueue.size() + " initial events.");
    }

    // Placeholder: Find the first KT associated with a JTS edge segment
    private KineticTriangle findTriangleContainingEdge(LineSegment segment, QuadEdgeSubdivision subdivision) {
         // This requires iterating subdivision edges, finding one matching segment coords,
         // then finding its associated KT via the triangleMap. Robust matching needed.
         // Simple version: Iterate KTs and check if segment matches one of their boundary edges.
         for (Map.Entry<Object, KineticTriangle> entry : triangleMap.entrySet()) {
             QuadEdge qe = (QuadEdge) entry.getKey(); // Assuming key is the representative edge
             Coordinate p0 = qe.orig().getCoordinate();
             Coordinate p1 = qe.dest().getCoordinate();
             LineSegment triSeg = new LineSegment(p0, p1);
             if (areSegmentsEqual(segment, triSeg)) return entry.getValue();

             qe = qe.lNext();
             p0 = qe.orig().getCoordinate();
             p1 = qe.dest().getCoordinate();
             triSeg = new LineSegment(p0, p1);
              if (areSegmentsEqual(segment, triSeg)) return entry.getValue();

             qe = qe.lNext();
              p0 = qe.orig().getCoordinate();
             p1 = qe.dest().getCoordinate();
             triSeg = new LineSegment(p0, p1);
              if (areSegmentsEqual(segment, triSeg)) return entry.getValue();
         }
        return null;
    }

     // Placeholder: Assign edge to correct slot in triangle (edge0, edge1, or edge2)
     private void assignEdgeToTriangle(KineticTriangle kt, LAEdge laEdge, LineSegment segment) {
         // Requires knowing the triangle's vertices and checking which edge is opposite which vertex.
         // This info isn't fully populated in the MVP init yet.
         // Simple assignment for now, assuming only one constrained edge per initial triangle.
         if (kt.edge0 == null) kt.edge0 = laEdge;
         else if (kt.edge1 == null) kt.edge1 = laEdge;
         else if (kt.edge2 == null) kt.edge2 = laEdge;
         else System.err.println("Warning: Triangle KT"+kt.id+" already has 3 edges assigned.");
         laEdge.triangle = kt;
     }

      // Placeholder: Find LAEdge matching a segment
     private LAEdge findLAEdgeForSegment(LineSegment segment) {
          // Requires robust key or iteration
          for(LAEdge edge : liveEdges) {
              if (edge.initialSegment != null && areSegmentsEqual(segment, edge.initialSegment)) {
                  return edge;
              }
          }
          return null;
     }

     // Helper for segment comparison (ignoring direction)
     private boolean areSegmentsEqual(LineSegment s1, LineSegment s2) {
         return (s1.p0.equals2D(s2.p0, EPSILON) && s1.p1.equals2D(s2.p1, EPSILON)) ||
                (s1.p0.equals2D(s2.p1, EPSILON) && s1.p1.equals2D(s2.p0, EPSILON));
     }


    private void propagateWavefront() {
        System.out.println("\nPropagating wavefront...");
        int eventCount = 0;
        while (!eventQueue.isEmpty()) {
            Event event = eventQueue.poll();

            if (event.time < currentTime - EPSILON) {
                 System.err.println("Warning: Event time "+event.time+" is before current time "+currentTime+". Skipping.");
                continue; // Time paradox! Should not happen with correct logic.
            }

            // Check if event is obsolete (triangle/edge already processed)
            if ((event.edge != null && event.edge.isProcessed) ||
                (event.triangle != null && event.triangle.isProcessed)) {
                System.out.println("Skipping obsolete event: " + event);
                continue;
            }

            currentTime = event.time;
            eventCount++;
            System.out.printf("--- Event %d @ T=%.2f: %s ---\n", eventCount, currentTime, event);

            switch (event.type) {
                case EDGE_COLLAPSE:
                    handleEdgeCollapse(event);
                    break;
                // Add cases for other event types here
                default:
                    System.err.println("Warning: Unhandled event type: " + event.type);
                    break;
            }
             // Basic sanity check for infinite loops
             if (eventCount > kineticTriangles.size() * 5 + liveEdges.size() * 5) { // Heuristic limit
                 System.err.println("Error: Exceeded maximum event count, potential infinite loop.");
                 return;
             }
        }
         System.out.println("Wavefront propagation finished after " + eventCount + " events. CurrentTime=" + String.format("%.2f", currentTime));
    }

    private void handleEdgeCollapse(Event event) {
        LAEdge collapsingEdge = event.edge;
        if (collapsingEdge == null || collapsingEdge.isProcessed) {
            System.out.println("Skipping processed or null edge collapse event at T=" + event.time);
            return;
        }

        LAVertex vStart = collapsingEdge.vStart;
        LAVertex vEnd = collapsingEdge.vEnd;
        Coordinate collapsePoint = event.collapsePoint;

        // Validate state
        if (vStart == null || vEnd == null || collapsePoint == null) {
            System.err.println("Error: Invalid edge collapse for LAE" + collapsingEdge.id + ": null vertex or collapse point.");
            collapsingEdge.isProcessed = true;
            return;
        }
        if (vStart.stopped || vEnd.stopped) {
            System.out.println("Skipping edge collapse LAE" + collapsingEdge.id + ": vertices LAV" + vStart.id + ", LAV" + vEnd.id + " already stopped.");
            vStart.stop(currentTime, collapsePoint);
            vEnd.stop(currentTime, collapsePoint);
            collapsingEdge.isProcessed = true;
            return;
        }

        System.out.println("Handling EDGE_COLLAPSE for " + collapsingEdge + " at " + formatCoord(collapsePoint));

        // Mark as processed
        collapsingEdge.isProcessed = true;
        if (collapsingEdge.triangle != null) {
            collapsingEdge.triangle.isProcessed = true;
        }

        // Stop vertices
        vStart.stop(currentTime, collapsePoint);
        vEnd.stop(currentTime, collapsePoint);
        System.out.println("  Stopped " + vStart + " at " + formatCoord(vStart.posStop));
        System.out.println("  Stopped " + vEnd + " at " + formatCoord(vEnd.posStop));

        // Find continuing edges
        LAEdge edgeA = (vStart.edge0 == collapsingEdge) ? vStart.edge1 : vStart.edge0;
        LAEdge edgeB = (vEnd.edge1 == collapsingEdge) ? vEnd.edge0 : vEnd.edge1;
        int sideA = (edgeA == vStart.edge1) ? 1 : 0;
        int sideB = (edgeB == vEnd.edge0) ? 0 : 1;

        System.out.println("  Continuing edge A: " + (edgeA != null ? "LAE" + edgeA.id + " (side " + sideA + ")" : "null"));
        System.out.println("  Continuing edge B: " + (edgeB != null ? "LAE" + edgeB.id + " (side " + sideB + ")" : "null"));

        // Handle termination cases
        if (edgeA == null || edgeB == null || edgeA == edgeB || edgeA.isProcessed || edgeB.isProcessed) {
            System.out.println("  Terminating at collapse point.");
            skeletonArcs.add(new SkeletonArc(vStart, null));
            skeletonArcs.add(new SkeletonArc(vEnd, null));
            if (edgeA != null) edgeA.isProcessed = true;
            if (edgeB != null) edgeB.isProcessed = true;
            return;
        }

        // Calculate new vertex velocity
     // Calculate new vertex velocity
        Vector2D newVelocity = calculateBisectorVelocity(edgeA, edgeB, collapsePoint, currentTime);
        if (newVelocity == null || newVelocity.length() < EPSILON) {
            System.out.println("  Zero or null velocity at collapse point. Terminating.");
            skeletonArcs.add(new SkeletonArc(vStart, null));
            skeletonArcs.add(new SkeletonArc(vEnd, null));
            edgeA.isProcessed = true;
            edgeB.isProcessed = true;
            return;
        }

        // Create new vertex
        LAVertex newVertex = new LAVertex(collapsePoint, newVelocity, currentTime, collapsePoint);
        liveVertices.add(newVertex);
        newVertex.edge0 = edgeA;
        newVertex.edge1 = edgeB;
        System.out.println("  Created " + newVertex + " vel=" + newVelocity);

        // Update continuing edges
        updateEdgeEndpoint(edgeA, vStart, newVertex);
        updateEdgeEndpoint(edgeB, vEnd, newVertex);

        // Create skeleton arcs
        vStart.setNext(newVertex, sideA);
        vEnd.setNext(newVertex, sideB);
        skeletonArcs.add(new SkeletonArc(vStart, newVertex));
        skeletonArcs.add(new SkeletonArc(vEnd, newVertex));
        System.out.println("  Added arc: LAV" + vStart.id + " -> LAV" + newVertex.id);
        System.out.println("  Added arc: LAV" + vEnd.id + " -> LAV" + newVertex.id);

        // Queue new events
        calculateAndAddEvent(edgeA, currentTime);
        calculateAndAddEvent(edgeB, currentTime);
    }
    
 // Helper to get the *other* vertex of an edge
    private LAVertex getOtherVertex(LAEdge edge, LAVertex vertex) {
        if (edge == null || vertex == null) return null;
        if (edge.vStart == vertex) return edge.vEnd;
        if (edge.vEnd == vertex) return edge.vStart;
        // Should not happen if vertex is indeed part of the edge
        System.err.println("Error: Vertex LAV"+vertex.id+" not found on edge LAE"+edge.id+" in getOtherVertex.");
        return null;
    }
    
    // Helper to update edge endpoint robustly
    private void updateEdgeEndpoint(LAEdge edge, LAVertex oldVertex, LAVertex newVertex) {
        if (edge == null || oldVertex == null || newVertex == null) {
            System.err.printf("Error: Null parameter in updateEdgeEndpoint (edge=%s, old=%s, new=%s)%n",
                edge == null ? "null" : "LAE" + edge.id,
                oldVertex == null ? "null" : "LAV" + oldVertex.id,
                newVertex == null ? "null" : "LAV" + newVertex.id);
            return;
        }
        if (edge.isProcessed) {
            System.out.println("Skipping update for processed edge LAE" + edge.id);
            return;
        }
        if (edge.vStart == oldVertex) {
            edge.vStart = newVertex;
            System.out.println("    Updated LAE" + edge.id + ".vStart to LAV" + newVertex.id);
        } else if (edge.vEnd == oldVertex) {
            edge.vEnd = newVertex;
            System.out.println("    Updated LAE" + edge.id + ".vEnd to LAV" + newVertex.id);
        } else {
            System.err.printf("Error: Edge LAE%d not connected to LAV%d%n", edge.id, oldVertex.id);
        }
    }
    
 // Helper to find the third vertex of a triangle given the edge's vertices
    private LAVertex findThirdVertex(KineticTriangle kt, LAVertex v1, LAVertex v2) {
        if (kt == null) return null;
        // This requires the KineticTriangle to actually store references to its LAVertices,
        // which the current MVP initialization doesn't fully establish.
        // Placeholder logic:
        if (kt.edge0 != null && kt.edge0.vStart != v1 && kt.edge0.vStart != v2) return kt.edge0.vStart;
        if (kt.edge0 != null && kt.edge0.vEnd != v1 && kt.edge0.vEnd != v2) return kt.edge0.vEnd;
        if (kt.edge1 != null && kt.edge1.vStart != v1 && kt.edge1.vStart != v2) return kt.edge1.vStart;
        if (kt.edge1 != null && kt.edge1.vEnd != v1 && kt.edge1.vEnd != v2) return kt.edge1.vEnd;
        if (kt.edge2 != null && kt.edge2.vStart != v1 && kt.edge2.vStart != v2) return kt.edge2.vStart;
        if (kt.edge2 != null && kt.edge2.vEnd != v1 && kt.edge2.vEnd != v2) return kt.edge2.vEnd;
        // This needs to be implemented based on how KTs store vertex info.
        System.err.println("Warning: findThirdVertex placeholder could not find third vertex for KT"+kt.id);
        return null;
    }

     // Placeholder: Find neighbor across a specific edge
     private KineticTriangle findNeighborTriangle(KineticTriangle kt, LAEdge edge) {
         // Requires robust neighbor info in KineticTriangle, which is missing in MVP init.
         // Return null for now.
         return null;
     }

     // Helper to calculate and add/update event for an edge
     private void calculateAndAddEvent(LAEdge edge, double currentTime) {
         if (edge == null || edge.isProcessed) return;

         // In a full implementation, we would remove any existing event for this edge first.
         // Java's PriorityQueue doesn't support easy removal/update.
         // Workaround: Add new event, handle obsolete events when polled.

         EdgeCollapseResult result = edge.calculateCollapseTime(currentTime);
         if (result.time < Double.POSITIVE_INFINITY) {
              System.out.println("  Queueing new event for "+edge+": Collapse at T="+String.format("%.2f", result.time)+" P="+formatCoord(result.point));
             eventQueue.add(new Event(result.time, edge.triangle, edge, result.point));
         }
     }


    private Geometry buildSkeletonGeometry() {
        System.out.println("\nBuilding final skeleton geometry...");
        List<LineString> lines = new ArrayList<>();
        Set<LAVertex> processedStarts = new HashSet<>(); // Avoid duplicate lines if structure is complex

        for (SkeletonArc arc : skeletonArcs) {
             LAVertex start = arc.startVertex;
             LAVertex end = arc.endVertex;

             if (start == null || start.posStart == null) {
                 System.err.println("Warning: Skipping arc with null start vertex/position.");
                 continue;
             }
             if (processedStarts.contains(start)) {
                 // This indicates a potential issue in skeleton arc tracking or definition,
                 // but for MVP we might see cases where a vertex is a start point twice
                 // if linking isn't perfect. Allow for now.
                 // System.err.println("Warning: Start vertex LAV"+start.id+" already processed for an arc.");
                 // continue;
             }
             processedStarts.add(start);


             Coordinate startCoord = start.posStart;
             Coordinate endCoord;

             if (end != null) {
                 // Finite arc: goes from start creation point to end creation point (collapse point)
                 endCoord = end.posStart; // End vertex starts where this arc ends
                 // Alternative: use stop positions?
                 // startCoord = start.posStart;
                 // endCoord = start.posStop; // Arc goes from start to where start stopped
                 // if (!end.posStart.equals2D(start.posStop, EPSILON)) {
                 //    System.err.printf("Warning: Arc end mismatch: end.posStart %s != start.posStop %s%n", formatCoord(end.posStart), formatCoord(start.posStop));
                 // }
             } else {
                 // Arc goes to infinity (or boundary) - represented by where the start vertex stopped
                  if (!start.stopped || start.posStop == null) {
                       System.err.println("Warning: Cannot draw arc for non-stopped vertex LAV"+start.id+" with null end.");
                       continue; // Skip arc if start didn't stop properly
                  }
                 endCoord = start.posStop;
             }


             // Avoid creating zero-length lines
             if (startCoord.distance(endCoord) > EPSILON) {
                 LineString line = gf.createLineString(new Coordinate[]{startCoord.copy(), endCoord.copy()});
                 lines.add(line);
                  System.out.println("  Created LineString: "+formatCoord(startCoord)+" -> "+formatCoord(endCoord));
             } else {
                   System.out.println("  Skipping zero-length arc for LAV"+start.id);
             }
        }

         if (lines.isEmpty()) {
             System.out.println("No skeleton lines generated.");
             return gf.createMultiLineString(null); // Return empty geometry
         }

        System.out.println("Generated " + lines.size() + " skeleton segments.");
        return gf.createMultiLineString(lines.toArray(new LineString[0]));
    }


    // -------------------------------------------------------------------------
    // Example Usage
    // -------------------------------------------------------------------------
    public static void main(String[] args) {
        System.out.println("Straight Skeleton Generator MVP - Test");

        // Test 1: Simple Triangle
        Coordinate[] triangleCoords = new Coordinate[]{
            new Coordinate(0, 0),
            new Coordinate(10, 0),
            new Coordinate(5, 8.66), // Approx equilateral
            new Coordinate(0, 0)
        };
        Polygon triangle = gf.createPolygon(triangleCoords);
        System.out.println("\nTest 1: Triangle Polygon: " + triangle);
        StraightSkeletonGeneratorMvp generator1 = new StraightSkeletonGeneratorMvp(triangle);
        Geometry skeleton1 = generator1.generate();
        System.out.println("Skeleton: " + (skeleton1 != null ? skeleton1.toText() : "Failed"));

        // Test 2: Square
        Coordinate[] squareCoords = new Coordinate[]{
            new Coordinate(0, 0),
            new Coordinate(10, 0),
            new Coordinate(10, 10),
            new Coordinate(0, 10),
            new Coordinate(0, 0)
        };
        Polygon square = gf.createPolygon(squareCoords);
        System.out.println("\nTest 2: Square Polygon: " + square);
        StraightSkeletonGeneratorMvp generator2 = new StraightSkeletonGeneratorMvp(square);
        Geometry skeleton2 = generator2.generate();
        System.out.println("Skeleton: " + (skeleton2 != null ? skeleton2.toText() : "Failed"));

        // Test 3: L-Shape
        Coordinate[] lShapeCoords = new Coordinate[]{
            new Coordinate(0, 0),
            new Coordinate(5, 0),
            new Coordinate(5, 5),
            new Coordinate(10, 5),
            new Coordinate(10, 10),
            new Coordinate(0, 10),
            new Coordinate(0, 0)
        };
        Polygon lShape = gf.createPolygon(lShapeCoords);
        System.out.println("\nTest 3: L-Shape Polygon: " + lShape);
        StraightSkeletonGeneratorMvp generator3 = new StraightSkeletonGeneratorMvp(lShape);
        Geometry skeleton3 = generator3.generate();
        System.out.println("Skeleton: " + (skeleton3 != null ? skeleton3.toText() : "Failed"));
        System.out.println("(Note: L-shape may produce incomplete skeleton due to missing split events)");
    }
}