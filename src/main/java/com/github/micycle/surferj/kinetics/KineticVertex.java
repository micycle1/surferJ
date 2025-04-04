package com.github.micycle.surferj.kinetics;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.math.Vector2D;
import org.tinfour.common.Vertex; // Keep reference to original

import com.github.micycle.surferj.util.Constants;
import com.github.micycle.surferj.util.GeometryUtil;

import java.util.Objects;

public class KineticVertex {
    private static int counter = 0;
    public final int id;

    public final Coordinate posZero; // Position at time t=0
    public final Vector2D velocity;
    public final double timeStart;
    public final Vertex originalTinVertex; // Link back to TinFour vertex

    private double timeStop = Double.POSITIVE_INFINITY;
    private Coordinate posStop = null;
    private boolean stopped = false;

    // Incident wavefront edges (order matters: edge[0] is CW, edge[1] is CCW)
    private KineticEdge[] incidentEdges = new KineticEdge[2];

    // Skeleton structure links (simplified)
    // In a full DCEL, these would be half-edges
    private KineticVertex nextVertex = null;
    private KineticVertex prevVertex = null;


    public KineticVertex(Coordinate posZero, Vector2D velocity, double timeStart, Vertex originalTinVertex) {
        this.id = counter++;
        this.posZero = Objects.requireNonNull(posZero);
        this.velocity = Objects.requireNonNull(velocity);
        this.timeStart = timeStart;
        this.originalTinVertex = originalTinVertex; // Can be null if synthetic
    }

    // Static factory for initial vertices
    public static KineticVertex createInitial(Coordinate pos, KineticEdge edgeA, KineticEdge edgeB, Vertex tinVertex) {
         // Calculate velocity based on incident edges
         // Note: Edges might be null during initial setup before full connectivity
         LineSegment segA = (edgeA != null) ? edgeA.getSegmentAt(0) : null;
         LineSegment segB = (edgeB != null) ? edgeB.getSegmentAt(0) : null;
         double weightA = (edgeA != null) ? edgeA.getWeight() : Constants.DEFAULT_WEIGHT;
         double weightB = (edgeB != null) ? edgeB.getWeight() : Constants.DEFAULT_WEIGHT;

         Vector2D vel = Vector2D.create(0,0); // Default
         if (segA != null && segB != null) {
             vel = GeometryUtil.calculateVelocity(pos, segA, segB, weightA, weightB);
         } else {
             System.err.println("Warning: Creating initial vertex " + counter + " with incomplete edge info.");
             // This might happen at polygon corners initially. Velocity might need recalculation later.
         }

         KineticVertex kv = new KineticVertex(pos, vel, 0.0, tinVertex);
         kv.setIncidentEdges(edgeA, edgeB); // Set edges after creation
         return kv;
    }
     // Static factory for event vertices
     public static KineticVertex createEventVertex(Coordinate pos, double time, KineticEdge edgeA, KineticEdge edgeB) {
         // Velocity needs careful calculation based on the *new* configuration
         // For simplicity, we might copy velocity from predecessors or recalculate
         LineSegment segA = edgeA.getSegmentAt(0); // Use t=0 segment for consistent velocity calc
         LineSegment segB = edgeB.getSegmentAt(0);
         Vector2D vel = GeometryUtil.calculateVelocity(pos, segA, segB, edgeA.getWeight(), edgeB.getWeight());

         // Calculate posZero based on pos, time, and velocity
         Coordinate posZero = new Coordinate(pos.x - vel.getX() * time, pos.y - vel.getY() * time);

         KineticVertex kv = new KineticVertex(posZero, vel, time, null); // No original Tin vertex
         kv.setIncidentEdges(edgeA, edgeB);
         return kv;
     }


    public Coordinate getPositionAt(double time) {
        if (stopped && time >= timeStop) {
            return posStop;
        }
        if (time < timeStart) {
             System.err.println("Warning: Querying position before start time for vertex " + id);
             // Extrapolate backwards - might be needed for some calculations
        }
        return new Coordinate(
            posZero.x + velocity.getX() * time,
            posZero.y + velocity.getY() * time
        );
    }

    public void stop(double time) {
        if (!stopped) {
            this.timeStop = time;
            this.posStop = getPositionAt(time);
            this.stopped = true;
            // TODO: Update skeleton structure (e.g., mark associated arc as segment)
        }
    }

    public boolean hasStopped() {
        return stopped;
    }
     public double getTimeStop() { return timeStop; }
     public Coordinate getPosStop() { return posStop; }
     public Coordinate getPosZero() { return posZero; }
     public Vector2D getVelocity() { return velocity; }
      public Vertex getOriginalTinVertex() { return originalTinVertex; }

     public void setIncidentEdges(KineticEdge edgeA, KineticEdge edgeB) {
         // Need consistent orientation, e.g., edgeA is CW, edgeB is CCW
         // This requires careful handling during triangulation setup
         this.incidentEdges[0] = edgeA; // Assume edgeA is CW for now
         this.incidentEdges[1] = edgeB; // Assume edgeB is CCW for now
         if (edgeA != null) edgeA.updateVertexRef(this);
         if (edgeB != null) edgeB.updateVertexRef(this);
     }

     public KineticEdge getIncidentEdge(int index) { return incidentEdges[index]; }

      // --- Skeleton Linking ---
      public void setNext(KineticVertex next) { this.nextVertex = next;}
      public void setPrev(KineticVertex prev) { this.prevVertex = prev;}
      public KineticVertex getNext() { return nextVertex; }
      public KineticVertex getPrev() { return prevVertex; }

    @Override
    public String toString() {
        return "KV" + id + "(" + String.format("%.2f", posZero.x) + "," + String.format("%.2f", posZero.y) + ")";
    }

    @Override
    public int hashCode() { return id; }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        KineticVertex other = (KineticVertex) obj;
        return id == other.id;
    }
}