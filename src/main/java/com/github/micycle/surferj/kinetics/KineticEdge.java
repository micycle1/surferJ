package com.github.micycle.surferj.kinetics;

import org.locationtech.jts.geom.LineSegment;
import org.tinfour.common.IQuadEdge; // Link back to original constraint

import com.github.micycle.surferj.event.CollapseEvent;
import com.github.micycle.surferj.event.EventType;
import com.github.micycle.surferj.util.Constants;
import com.github.micycle.surferj.util.GeometryUtil;

import java.util.Objects;

public class KineticEdge {
    private static int counter = 0;
    public final int id;

    private final LineSegment segmentZero; // Segment at time t=0
    private final double weight;
    private KineticVertex[] vertices = new KineticVertex[2]; // Endpoints: vertices[0] = start, vertices[1] = end
    private KineticTriangle incidentTriangle; // The triangle "behind" this edge
    private boolean stopped = false;
     public final IQuadEdge originalTinEdge; // Link back

    public KineticEdge(KineticVertex v0, KineticVertex v1, double weight, IQuadEdge originalTinEdge) {
        this.id = counter++;
        this.vertices[0] = Objects.requireNonNull(v0);
        this.vertices[1] = Objects.requireNonNull(v1);
        this.segmentZero = new LineSegment(v0.getPosZero(), v1.getPosZero());
        this.weight = weight;
         this.originalTinEdge = originalTinEdge; // Can be null if synthetic

        // Initial edge setup (might be called later)
        v0.setIncidentEdges(null, this); // Assume this edge is CCW for v0 initially
        v1.setIncidentEdges(this, null); // Assume this edge is CW for v1 initially
    }

    public LineSegment getSegmentAt(double time) {
        if (stopped) {
            // Maybe return the final segment? Or null?
             return new LineSegment(vertices[0].getPositionAt(vertices[0].getTimeStop()),
                                    vertices[1].getPositionAt(vertices[1].getTimeStop()));
        }
        return new LineSegment(vertices[0].getPositionAt(time), vertices[1].getPositionAt(time));
    }

    public KineticVertex getVertex(int index) {
        return vertices[index];
    }

     public double getWeight() { return weight; }
    public KineticTriangle getIncidentTriangle() { return incidentTriangle; }
    public void setIncidentTriangle(KineticTriangle t) { this.incidentTriangle = t; }
     public IQuadEdge getOriginalTinEdge() { return originalTinEdge; }


    // Called by KineticVertex to ensure references are correct
    public void updateVertexRef(KineticVertex vertex) {
        if (vertices[0] == vertex) {
            // Vertex is start, ensure its CCW edge is this
            if (vertex.getIncidentEdge(1) != this) vertex.setIncidentEdges(vertex.getIncidentEdge(0), this);
        } else if (vertices[1] == vertex) {
             // Vertex is end, ensure its CW edge is this
             if (vertex.getIncidentEdge(0) != this) vertex.setIncidentEdges(this, vertex.getIncidentEdge(1));
        } else {
             System.err.println("Error: Vertex " + vertex.id + " claims incidence to edge " + id + " but is not an endpoint.");
        }
    }


    public CollapseEvent computeCollapseEvent(double timeNow) {
        if (stopped || vertices[0].hasStopped() || vertices[1].hasStopped()) {
            return new CollapseEvent(EventType.NONE, Double.POSITIVE_INFINITY, null); // No future event
        }

        KineticVertex v0 = vertices[0];
        KineticVertex v1 = vertices[1];

        // Relative position and velocity
        double dx = v1.posZero.x - v0.posZero.x;
        double dy = v1.posZero.y - v0.posZero.y;
        double dvx = v1.velocity.getX() - v0.velocity.getX();
        double dvy = v1.velocity.getY() - v0.velocity.getY();

        // Edge collapses when squared distance is 0: || P1(t) - P0(t) ||^2 = 0
        // (x1(0)+vx1*t - x0(0)-vx0*t)^2 + (y1(0)+vy1*t - y0(0)-vy0*t)^2 = 0
        // (dx + dvx*t)^2 + (dy + dvy*t)^2 = 0
        // (dvx^2 + dvy^2) * t^2 + 2*(dx*dvx + dy*dvy) * t + (dx^2 + dy^2) = 0
        // This is a quadratic equation At^2 + Bt + C = 0

        double a = dvx * dvx + dvy * dvy; // Squared relative speed
        double b = 2 * (dx * dvx + dy * dvy);
        double c = dx * dx + dy * dy;     // Initial squared distance

        if (Math.abs(a) < Constants.EPSILON) {
            // Relative velocity is zero or near zero
            if (Math.abs(c) < Constants.EPSILON) {
                // Already collapsed (or very close) and not moving relatively
                 // Need to check if moving apart *exactly* collinearly - complex
                 // For now, treat as potentially collapsing now if close
                 if (getSegmentAt(timeNow).getLength() < Constants.EPSILON) {
                      System.out.println("Edge " + id + " appears parallel and coincident at time " + timeNow);
                      // This might be an "ALWAYS" case in C++, needs careful handling.
                      // Return an event slightly in the future to avoid infinite loops? Or handle specifically.
                      // Returning NONE for now, as this state shouldn't persist if velocities calculated right.
                       return new CollapseEvent(EventType.NONE, Double.POSITIVE_INFINITY, incidentTriangle);
                 } else {
                     return new CollapseEvent(EventType.NONE, Double.POSITIVE_INFINITY, incidentTriangle); // Parallel, distinct
                 }
            } else {
                 // Parallel, moving at same speed, never collapse
                 return new CollapseEvent(EventType.NONE, Double.POSITIVE_INFINITY, incidentTriangle);
            }
        }

        // Solve At^2 + Bt + C = 0 for t >= timeNow
        double collapseTime = GeometryUtil.solveQuadratic(a, b, c);

        // Ensure the collapse is not in the past relative to timeNow
        if (collapseTime < timeNow - Constants.EPSILON) {
             // Check if B < 0, indicating they were moving towards each other in the past
             if (b < -Constants.EPSILON) {
                 // They collided in the past and are now moving apart
                 return new CollapseEvent(EventType.NONE, Double.POSITIVE_INFINITY, incidentTriangle); // PAST event type
             } else {
                  // Roots are in the past, but they are currently moving apart or parallel
                  return new CollapseEvent(EventType.NONE, Double.POSITIVE_INFINITY, incidentTriangle);
             }
        }

         // If collapse time is very close to now, clamp it to now?
         if (Math.abs(collapseTime - timeNow) < Constants.EPSILON) {
             collapseTime = timeNow;
         }


        if (collapseTime != Double.POSITIVE_INFINITY) {
            return new CollapseEvent(EventType.EDGE_COLLAPSE, collapseTime, incidentTriangle);
        } else {
            return new CollapseEvent(EventType.NONE, Double.POSITIVE_INFINITY, incidentTriangle);
        }
    }

     public void stop() {
         this.stopped = true;
         // Optionally stop vertices if they aren't already? Depends on event logic.
     }
     public boolean hasStopped() { return stopped; }


    @Override
    public String toString() {
        return "KE" + id + "(" + vertices[0].id + "," + vertices[1].id + ")";
    }
     @Override
     public int hashCode() { return id; }

     @Override
     public boolean equals(Object obj) {
         if (this == obj) return true;
         if (obj == null || getClass() != obj.getClass()) return false;
         KineticEdge other = (KineticEdge) obj;
         return id == other.id;
     }
}