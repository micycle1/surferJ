package com.github.micycle.surferj.kinetics;

import org.locationtech.jts.geom.Coordinate;
import org.tinfour.common.SimpleTriangle; // Link back

import com.github.micycle.surferj.event.CollapseEvent;
import com.github.micycle.surferj.event.EventType;
import com.github.micycle.surferj.util.Constants;
import com.github.micycle.surferj.util.GeometryUtil;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class KineticTriangle {
    private static int counter = 0;
    public final int id;

    KineticVertex[] vertices = new KineticVertex[3]; // v0, v1, v2 in CCW order
    private KineticEdge[] edges = new KineticEdge[3]; // edge opposite vertex i (e.g., edges[0] is edge v1-v2)
    private KineticTriangle[] neighbors = new KineticTriangle[3]; // neighbor opposite vertex i
    private boolean stopped = false;
    public final SimpleTriangle originalTinTriangle; // Link back

    public KineticTriangle(KineticVertex v0, KineticVertex v1, KineticVertex v2, SimpleTriangle originalTinTriangle) {
        this.id = counter++;
        // Ensure CCW order
        if (GeometryUtil.orientation(v0.posZero, v1.posZero, v2.posZero) < 0) {
            this.vertices[0] = Objects.requireNonNull(v0);
            this.vertices[1] = Objects.requireNonNull(v2); // Swap v1 and v2
            this.vertices[2] = Objects.requireNonNull(v1);
        } else {
            this.vertices[0] = Objects.requireNonNull(v0);
            this.vertices[1] = Objects.requireNonNull(v1);
            this.vertices[2] = Objects.requireNonNull(v2);
        }
         this.originalTinTriangle = originalTinTriangle;
    }

    public KineticVertex getVertex(int index) {
        return vertices[index];
    }

    public KineticEdge getEdge(int index) {
        return edges[index];
    }

    public KineticTriangle getNeighbor(int index) {
        return neighbors[index];
    }

     // Set edge opposite vertex i
     public void setEdge(int index, KineticEdge edge) {
         edges[index] = edge;
         if (edge != null) {
             edge.setIncidentTriangle(this);
         }
     }

     // Set neighbor opposite vertex i
     public void setNeighbor(int index, KineticTriangle neighbor) {
         neighbors[index] = neighbor;
     }

     public int getVertexIndex(KineticVertex v) {
         for (int i = 0; i < 3; i++) {
             if (vertices[i].equals(v)) return i;
         }
         return -1;
     }
      public int getNeighborIndex(KineticTriangle n) {
          for (int i = 0; i < 3; i++) {
              if (neighbors[i] != null && neighbors[i].equals(n)) return i;
          }
          return -1;
      }

    public CollapseEvent computeCollapseEvent(double timeNow) {
        if (stopped) return new CollapseEvent(EventType.NONE, Double.POSITIVE_INFINITY, this);

        List<CollapseEvent> potentialEvents = new ArrayList<>();

        // --- 1. Edge Collapse Events ---
        for (int i = 0; i < 3; i++) {
            if (edges[i] != null) { // Is it a wavefront edge?
                CollapseEvent edgeEvent = edges[i].computeCollapseEvent(timeNow);
                 if (edgeEvent.getType() == EventType.EDGE_COLLAPSE) {
                     // Add edge index info? C++ version stores this.
                     // We might need it for handling. Let's store it in the event.
                     edgeEvent.setEdgeIndex(i); // Custom field in CollapseEvent needed
                     potentialEvents.add(edgeEvent);
                 }
            }
        }

        // --- 2. Triangle Collapse Event ---
        // Check when the triangle area becomes zero (vertices become collinear)
        // Area(t) = 0.5 * | (x1(t)-x0(t))*(y2(t)-y0(t)) - (x2(t)-x0(t))*(y1(t)-y0(t)) |
        // Area = 0 when (x1-x0)(y2-y0) - (x2-x0)(y1-y0) = 0
        // Substitute x_i(t) = x_i(0) + vx_i*t, y_i(t) = y_i(0) + vy_i*t
        // This results in a quadratic equation in t: At^2 + Bt + C = 0
        KineticVertex v0 = vertices[0];
        KineticVertex v1 = vertices[1];
        KineticVertex v2 = vertices[2];

        double x0=v0.posZero.x, y0=v0.posZero.y, vx0=v0.velocity.getX(), vy0=v0.velocity.getY();
        double x1=v1.posZero.x, y1=v1.posZero.y, vx1=v1.velocity.getX(), vy1=v1.velocity.getY();
        double x2=v2.posZero.x, y2=v2.posZero.y, vx2=v2.velocity.getX(), vy2=v2.velocity.getY();

        // Coefficients for (x1-x0)(y2-y0) - (x2-x0)(y1-y0) = 0
        // A = (vx1-vx0)(vy2-vy0) - (vx2-vx0)(vy1-vy0)
        // B = (x1-x0)(vy2-vy0) + (vx1-vx0)(y2-y0) - (x2-x0)(vy1-vy0) - (vx2-vx0)(y1-y0)
        // C = (x1-x0)(y2-y0) - (x2-x0)(y1-y0)  ( == 2 * signed area at t=0 )

        double A = (vx1 - vx0) * (vy2 - vy0) - (vx2 - vx0) * (vy1 - vy0);
        double B = (x1 - x0) * (vy2 - vy0) + (vx1 - vx0) * (y2 - y0)
                 - (x2 - x0) * (vy1 - vy0) - (vx2 - vx0) * (y1 - y0);
        double C = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

        double collapseTime = GeometryUtil.solveQuadratic(A, B, C);

         // Check if valid collapse (not in past, and shrinking)
         if (collapseTime >= timeNow - Constants.EPSILON && collapseTime != Double.POSITIVE_INFINITY) {
             // Check derivative at collapseTime: 2At + B
             double derivative = 2 * A * collapseTime + B;
             // We need the signed area to be decreasing (or zero derivative if A=0)
             // The sign of the area at timeNow is sign(At^2+Bt+C) evaluated at timeNow
             double areaNow = A*timeNow*timeNow + B*timeNow + C;

             // If area is positive now, derivative at collapse should be negative (or zero)
             // If area is negative now (shouldn't happen if always CCW), derivative should be positive (or zero)
             // Simplified: If A != 0, the sign of A determines parabola opening.
             // If A > 0 (opens up), we want the *first* root (smaller t), derivative should be negative.
             // If A < 0 (opens down), we want the *second* root (larger t), derivative should be positive.
             // If A = 0 (linear), derivative is B. Collapse happens if B has opposite sign to C (area at t=0).

             boolean acceptCollapse = false;
             if (Math.abs(A) > Constants.EPSILON) {
                 // Quadratic case
                  if (A > 0 && derivative <= Constants.EPSILON) acceptCollapse = true; // Opens up, want first root (decreasing area)
                  if (A < 0 && derivative >= -Constants.EPSILON) acceptCollapse = true; // Opens down, want second root (decreasing area)
             } else {
                 // Linear case (A approx 0)
                 if (Math.abs(B) > Constants.EPSILON) {
                    // If B > 0, area increases. Collapse only if C (area at t=0) was negative.
                    // If B < 0, area decreases. Collapse only if C was positive.
                    // Since we assume C >= 0 initially, only accept if B < 0.
                    if (B < -Constants.EPSILON) acceptCollapse = true;
                 } else {
                     // A=0, B=0. Area is constant. Collapse only if C=0 (always collinear).
                      if (Math.abs(C) < Constants.EPSILON) acceptCollapse = true; // Already collapsed
                 }
             }


             if (acceptCollapse) {
                  // Distinguish between triangle collapse and flip/spoke collapse based on edge lengths at collapse time
                  Coordinate p0 = v0.getPositionAt(collapseTime);
                  Coordinate p1 = v1.getPositionAt(collapseTime);
                  Coordinate p2 = v2.getPositionAt(collapseTime);
                  double d01 = p0.distance(p1);
                  double d12 = p1.distance(p2);
                  double d20 = p2.distance(p0);

                  // Check if one edge length is zero (spoke collapse) or all are zero (triangle collapse)
                  // Need tolerance
                   int zeroEdges = 0;
                   int zeroEdgeIndex = -1;
                   if (d12 < Constants.EPSILON) { zeroEdges++; zeroEdgeIndex = 0; } // Edge opposite v0
                   if (d20 < Constants.EPSILON) { zeroEdges++; zeroEdgeIndex = 1; } // Edge opposite v1
                   if (d01 < Constants.EPSILON) { zeroEdges++; zeroEdgeIndex = 2; } // Edge opposite v2


                   if (zeroEdges >= 2) { // All points coincide
                       potentialEvents.add(new CollapseEvent(EventType.TRIANGLE_COLLAPSE, collapseTime, this));
                   } else if (zeroEdges == 1) {
                        // Spoke collapse - the edge opposite zeroEdgeIndex collapsed
                        // In C++, this is SPOKE_COLLAPSE with the *non-collapsing* edge index
                        // Let's use TRIANGLE_COLLAPSE for now, refinement needed.
                         potentialEvents.add(new CollapseEvent(EventType.TRIANGLE_COLLAPSE, collapseTime, this));
                        // potentialEvents.add(new CollapseEvent(EventType.SPOKE_COLLAPSE, collapseTime, this, (zeroEdgeIndex + 1) % 3)); // Example: Pass index of a non-zero edge
                   } else {
                       // Flip event - Collinear, but no edges collapsed
                       // In C++, this is VERTEX_MOVES_OVER_SPOKE with the index of the edge being crossed.
                       // The edge being crossed is the one *opposite* the vertex that moved.
                       // Needs determining which vertex moved "through". For now, use generic FLIP.
                       potentialEvents.add(new CollapseEvent(EventType.FLIP_EVENT, collapseTime, this));
                   }
             }
         }


        // --- 3. Split/Refine Events (Deferred) ---
        // Would involve checking vertex against non-incident wavefront edges.

        // --- Find minimum time event ---
        CollapseEvent minEvent = new CollapseEvent(EventType.NONE, Double.POSITIVE_INFINITY, this);
        for (CollapseEvent event : potentialEvents) {
            if (event.compareTo(minEvent) < 0) {
                minEvent = event;
            }
        }

        // If multiple events happen at the *exact* same time, C++ code prioritizes.
        // EDGE_COLLAPSE < TRIANGLE_COLLAPSE < FLIP_EVENT etc.
        // Our compareTo in CollapseEvent handles this.

        return minEvent;
    }

     public void stop() {
         this.stopped = true;
         // Stop associated edges/vertices? Depends on event logic.
     }
     public boolean hasStopped() { return stopped; }

    @Override
    public String toString() {
        return "KT" + id + "(" + vertices[0].id + "," + vertices[1].id + "," + vertices[2].id + ")";
    }
     @Override
     public int hashCode() { return id; }

     @Override
     public boolean equals(Object obj) {
         if (this == obj) return true;
         if (obj == null || getClass() != obj.getClass()) return false;
         KineticTriangle other = (KineticTriangle) obj;
         return id == other.id;
     }
}