package com.github.micycle.surferj.kinetics;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.math.Vector2D;
import org.tinfour.common.IQuadEdge;
import org.tinfour.common.SimpleTriangle;
import org.tinfour.common.Vertex;

import com.github.micycle.surferj.SkeletonOutput;
import com.github.micycle.surferj.event.CollapseEvent;
import com.github.micycle.surferj.event.EventQueue;
import com.github.micycle.surferj.event.EventType;
import com.github.micycle.surferj.triangulation.InitialTriangulator;
import com.github.micycle.surferj.util.Constants;
import com.github.micycle.surferj.util.GeometryUtil;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

public class KineticTriangulation {

    private final InitialTriangulator initialTriangulator;
    private List<KineticVertex> vertices;
    private List<KineticEdge> edges;
    private List<KineticTriangle> triangles;
    private Map<SimpleTriangle, KineticTriangle> tinToKineticMap;
    private EventQueue eventQueue; // Needs to be set externally

    public KineticTriangulation(InitialTriangulator initialTriangulator) {
        this.initialTriangulator = Objects.requireNonNull(initialTriangulator);
        this.vertices = new ArrayList<>();
        this.edges = new ArrayList<>();
        this.triangles = new ArrayList<>();
        this.tinToKineticMap = new HashMap<>();
    }

    public void initialize(EventQueue queue) {
        this.eventQueue = Objects.requireNonNull(queue);
        createKineticVertices();
        createKineticEdges();
        createKineticTriangles();
        establishConnectivity();
        calculateInitialEvents();
    }

    private void createKineticVertices() {
        Map<Vertex, KineticVertex> vertexMap = new HashMap<>();
        for (Vertex v : initialTriangulator.getTin().getVertices()) {
            // Initial velocity is tricky - needs neighboring edges. Calculate later?
            // For now, create with zero velocity, update in establishConnectivity
            KineticVertex kv = new KineticVertex(v.getCoordinate(), Vector2D.create(0,0), 0.0, v);
            vertices.add(kv);
            vertexMap.put(v, kv);
        }
         // Store map for lookup during edge/triangle creation
         // (This assumes vertices list maps 1:1 with getVertices result, which is typical)
         System.out.println("Created " + vertices.size() + " initial kinetic vertices.");
    }

    private Map<Vertex, KineticVertex> getVertexMap() {
         Map<Vertex, KineticVertex> vertexMap = new HashMap<>();
         for (KineticVertex kv : vertices) {
             if (kv.getOriginalTinVertex() != null) {
                 vertexMap.put(kv.getOriginalTinVertex(), kv);
             }
         }
         return vertexMap;
    }


    private void createKineticEdges() {
         Map<Vertex, KineticVertex> vertexMap = getVertexMap();
         Map<IQuadEdge, KineticEdge> edgeMap = new HashMap<>(); // Avoid duplicates

        for (IQuadEdge edge : initialTriangulator.getTin().getEdges()) {
             if (edge.isConstraint()) {
                  // Only create one KineticEdge per constraint pair (edge/twin)
                  if (edgeMap.containsKey(edge) || edgeMap.containsKey(edge.getTwin())) {
                      continue;
                  }

                 Vertex vA = edge.getA();
                 Vertex vB = edge.getB();
                 KineticVertex kvA = vertexMap.get(vA);
                 KineticVertex kvB = vertexMap.get(vB);

                 if (kvA != null && kvB != null) {
                     double weight = initialTriangulator.getConstraintWeight(edge);
                     KineticEdge ke = new KineticEdge(kvA, kvB, weight, edge);
                     edges.add(ke);
                     edgeMap.put(edge, ke);
                     edgeMap.put(edge.getTwin(), ke); // Map twin as well
                 } else {
                      System.err.println("Error: Could not find KineticVertices for constraint edge " + edge);
                 }
             }
         }
         System.out.println("Created " + edges.size() + " initial kinetic edges.");
    }

     private Map<IQuadEdge, KineticEdge> getEdgeMap() {
         Map<IQuadEdge, KineticEdge> edgeMap = new HashMap<>();
         for(KineticEdge ke : edges) {
             if (ke.getOriginalTinEdge() != null) {
                 edgeMap.put(ke.getOriginalTinEdge(), ke);
                 edgeMap.put(ke.getOriginalTinEdge().getTwin(), ke); // Map twin too
             }
         }
         return edgeMap;
     }

    private void createKineticTriangles() {
        Map<Vertex, KineticVertex> vertexMap = getVertexMap();
        for (SimpleTriangle tinTri : initialTriangulator.getTriangles()) {
            Vertex v0 = tinTri.getVertexA();
            Vertex v1 = tinTri.getVertexB();
            Vertex v2 = tinTri.getVertexC();

            KineticVertex kv0 = vertexMap.get(v0);
            KineticVertex kv1 = vertexMap.get(v1);
            KineticVertex kv2 = vertexMap.get(v2);

            if (kv0 != null && kv1 != null && kv2 != null) {
                KineticTriangle kt = new KineticTriangle(kv0, kv1, kv2, tinTri);
                triangles.add(kt);
                tinToKineticMap.put(tinTri, kt);
            } else {
                 System.err.println("Error: Could not find KineticVertices for triangle " + tinTri);
            }
        }
        System.out.println("Created " + triangles.size() + " initial kinetic triangles.");
    }

    private void establishConnectivity() {
         Map<Vertex, KineticVertex> vertexMap = getVertexMap();
         Map<IQuadEdge, KineticEdge> edgeMap = getEdgeMap(); // Map includes twins

        // 1. Link Triangles to Edges and Neighbors
        for (KineticTriangle kt : triangles) {
            SimpleTriangle tinTri = kt.originalTinTriangle;
            for (int i = 0; i < 3; i++) {
                IQuadEdge tinEdge = tinTri.getEdge(i); // Edge opposite vertex i
                KineticEdge ke = edgeMap.get(tinEdge); // Find corresponding KineticEdge
                kt.setEdge(i, ke); // Can be null if not a constraint

                // Find neighbor triangle
                IQuadEdge twin = tinEdge.getTwin();
                SimpleTriangle neighborTinTri = initialTriangulator.getTin().getTriangle(twin); // Get triangle attached to twin
                if (neighborTinTri != null) {
                    KineticTriangle neighborKt = tinToKineticMap.get(neighborTinTri);
                    kt.setNeighbor(i, neighborKt);
                }
            }
        }

         // 2. Link Vertices to Edges and Calculate Initial Velocities
         for (KineticVertex kv : vertices) {
             Vertex tinV = kv.originalTinVertex;
             if (tinV == null) continue;

             List<IQuadEdge> edgesAround = tinV.getEdges();
             KineticEdge edgeA = null; // CW edge
             KineticEdge edgeB = null; // CCW edge

             // Find the two *constraint* edges incident to this vertex
             // This logic is simplified and might need refinement for complex cases
             List<KineticEdge> incidentConstraints = new ArrayList<>();
             for (IQuadEdge te : edgesAround) {
                 if (te.isConstraint()) {
                     KineticEdge ke = edgeMap.get(te);
                     if (ke != null && !incidentConstraints.contains(ke)) {
                         incidentConstraints.add(ke);
                     }
                 }
             }

             if (incidentConstraints.size() == 2) {
                  // Determine CW/CCW order (relative to vertex)
                  // This requires careful geometric check or knowledge from TIN structure
                  edgeA = incidentConstraints.get(0); // Placeholder
                  edgeB = incidentConstraints.get(1); // Placeholder
                  // TODO: Implement robust CW/CCW ordering check
                  kv.setIncidentEdges(edgeA, edgeB);

                  // Now calculate velocity
                  Vector2D vel = GeometryUtil.calculateVelocity(kv.posZero,
                                                               edgeA.getSegmentAt(0),
                                                               edgeB.getSegmentAt(0),
                                                               edgeA.getWeight(), edgeB.getWeight());
                   // Recreate vertex with correct velocity? Or update? Update is easier.
                   // Need to make velocity mutable in KineticVertex or use a builder pattern.
                   // For now, let's assume KineticVertex constructor was called *after* edges were known
                   // (requires restructuring initialization) OR we update it here (requires mutable field).
                   // HACK: Create a new vertex and replace - less efficient but avoids mutable fields now.
                   // vertices.set(vertices.indexOf(kv), new KineticVertex(kv.posZero, vel, kv.timeStart, kv.originalTinVertex));
                   // System.out.println("WARN: Velocity update mechanism is basic.");

             } else if (incidentConstraints.size() > 2) {
                  System.err.println("Warning: Vertex " + kv.id + " has > 2 incident constraints. Beveling needed (not implemented).");
             } else {
                 System.err.println("Warning: Vertex " + kv.id + " has < 2 incident constraints. Boundary vertex?");
                  // Handle boundary cases - velocity might be based on one edge + angle? Simplified for now.
             }
         }
         System.out.println("Established connectivity.");
    }

    private void calculateInitialEvents() {
        for (KineticTriangle kt : triangles) {
            CollapseEvent event = kt.computeCollapseEvent(0.0);
            if (event.getType() != EventType.NONE) {
                eventQueue.add(event);
            }
        }
         System.out.println("Calculated initial events. Queue size: " + eventQueue.size());
    }

    // --- Event Handling ---

    public void handleEvent(CollapseEvent event) {
        System.out.println("Handling event: " + event);
        KineticTriangle triangle = event.getTriangle();
        if (triangle.hasStopped()) {
            System.out.println("  Triangle already stopped. Skipping.");
            return; // Already processed
        }

        double time = event.getTime();

        switch (event.getType()) {
            case EDGE_COLLAPSE:
                handleEdgeCollapse(event, time);
                break;
            case TRIANGLE_COLLAPSE:
                 handleTriangleCollapse(event, time);
                break;
            case FLIP_EVENT:
                 handleFlipEvent(event, time);
                break;
            // Add cases for SPOKE_COLLAPSE, SPLIT_OR_FLIP etc. when implemented
            default:
                System.err.println("Unhandled event type: " + event.getType());
        }

        // Invalidate events of affected neighbors (and the triangle itself if it wasn't stopped)
        updateAffectedEvents(triangle, time);
    }

    private void handleEdgeCollapse(CollapseEvent event, double time) {
        KineticTriangle t = event.getTriangle();
        int edgeIndex = event.getEdgeIndex(); // Assumes CollapseEvent stores this
        KineticEdge edge = t.getEdge(edgeIndex);

        if (edge == null || edge.hasStopped()) return;

        System.out.println("  Edge Collapse: Triangle " + t.id + ", Edge opposite vertex " + edgeIndex);

        KineticVertex v0 = edge.getVertex(0);
        KineticVertex v1 = edge.getVertex(1);

        // Stop the collapsing edge and its vertices
        edge.stop();
        v0.stop(time);
        v1.stop(time);

        // Create the new vertex at the collision point
        Coordinate collisionPoint = v0.getPosStop(); // Should be same as v1.getPosStop()
        KineticEdge edgeA = v0.getIncidentEdge(0); // Edge CW from v0 (not the collapsing one)
        KineticEdge edgeB = v1.getIncidentEdge(1); // Edge CCW from v1 (not the collapsing one)

        // Check if edgeA and edgeB are valid (might be null at boundary)
        if (edgeA == null || edgeB == null) {
            System.err.println("Error: Null edge encountered during edge collapse vertex creation at vertex " + v0.id + " or " + v1.id);
             // This indicates a boundary issue or incomplete connectivity handling
             t.stop(); // Stop the triangle as a fallback
             return;
        }


        KineticVertex newVertex = KineticVertex.createEventVertex(collisionPoint, time, edgeA, edgeB);
        vertices.add(newVertex);

        // Update topology:
        // The triangle t is removed.
        t.stop();

        // Update neighbor pointers and vertex references in adjacent triangles
        // This involves iterating around v0 and v1 before the collapse, finding the
        // triangles that shared v0/v1, and updating their vertex pointers to newVertex.
        // Also need to update neighbor links across the collapsed edge's "gap".
        // This is complex and requires careful iteration/pointer updates.

        // Simplified update (may leave dangling refs):
        KineticTriangle neighbor = t.getNeighbor(edgeIndex); // Neighbor across the non-wavefront edge
        if (neighbor != null && !neighbor.hasStopped()) {
             int neighborEdgeIndex = neighbor.getNeighborIndex(t);
             if (neighborEdgeIndex != -1) {
                 neighbor.setNeighbor(neighborEdgeIndex, null); // Remove link to collapsed triangle
                 // TODO: Update neighbor's vertex refs if needed
                 eventQueue.update(neighbor, time); // Mark neighbor for event recalculation
             }
        }
         // Update triangles adjacent to v0 and v1
         updateVertexReferences(v0, newVertex, time);
         updateVertexReferences(v1, newVertex, time);

         // Link new vertex in skeleton structure (simplified)
         v0.setNext(newVertex);
         newVertex.setPrev(v0);
         v1.setPrev(newVertex);
         newVertex.setNext(v1);

        System.out.println("  Created new vertex: " + newVertex.id);
    }

     private void handleTriangleCollapse(CollapseEvent event, double time) {
         KineticTriangle t = event.getTriangle();
         System.out.println("  Triangle Collapse: Triangle " + t.id);

         if (t.hasStopped()) return;
         t.stop();

         // Stop all vertices and edges
         KineticVertex v0=null, v1=null, v2=null; // For linking
         for (int i = 0; i < 3; i++) {
             KineticVertex v = t.getVertex(i);
             if (!v.hasStopped()) v.stop(time);
             if (i==0) v0=v; else if (i==1) v1=v; else v2=v;

             KineticEdge edge = t.getEdge(i);
             if (edge != null && !edge.hasStopped()) edge.stop();
         }

         // Update neighbors
         for (int i = 0; i < 3; i++) {
              KineticTriangle neighbor = t.getNeighbor(i);
              if (neighbor != null && !neighbor.hasStopped()) {
                  int neighborEdgeIndex = neighbor.getNeighborIndex(t);
                  if (neighborEdgeIndex != -1) {
                      neighbor.setNeighbor(neighborEdgeIndex, null);
                      eventQueue.update(neighbor, time);
                  }
              }
         }
          // Link vertices in skeleton structure (assuming they meet at one point)
          if(v0 != null && v1 != null && v2 != null) {
              v0.setNext(v1); v1.setPrev(v0);
              v1.setNext(v2); v2.setPrev(v1);
              v2.setNext(v0); v0.setPrev(v2); // Form a small loop at the collapse point
          }
     }

      private void handleFlipEvent(CollapseEvent event, double time) {
          KineticTriangle t1 = event.getTriangle();
          System.out.println("  Flip Event: Triangle " + t1.id);

          // A flip involves two triangles, t1 and t2, sharing an edge e (which must NOT be a constraint)
          // Find the edge e to flip. In C++, the event often stores the index of the vertex *opposite* the edge to flip.
          // Let's assume the event tells us which edge index *within t1* is the flip edge.
          int edgeIndexInT1 = event.getEdgeIndex(); // Need to add this to CollapseEvent based on calculation logic
          if (t1.getEdge(edgeIndexInT1) != null) {
              System.err.println("Error: Trying to flip a constraint edge in T" + t1.id);
              return;
          }
          KineticTriangle t2 = t1.getNeighbor(edgeIndexInT1);
          if (t2 == null || t2.hasStopped() || t1.hasStopped()) {
              System.out.println("  Cannot flip: neighbor missing or stopped.");
              return; // Cannot flip with boundary or stopped triangle
          }
          int edgeIndexInT2 = t2.getNeighborIndex(t1);
          if (edgeIndexInT2 == -1) {
               System.err.println("Error: Neighbor T" + t2.id + " does not point back to T" + t1.id);
               return;
          }

          // Vertices involved:
          // t1: v0, v1, v2 (CCW) -> edgeIndexInT1 is opposite v[edgeIndexInT1]
          // t2: v3, v2, v1 (CCW) -> where edge is v1-v2. edgeIndexInT2 is opposite v3.
          KineticVertex v_opp_t1 = t1.getVertex(edgeIndexInT1); // Vertex in t1 opposite the shared edge
          KineticVertex v_shared1 = t1.getVertex((edgeIndexInT1 + 1) % 3); // Shared vertex 1
          KineticVertex v_shared2 = t1.getVertex((edgeIndexInT1 + 2) % 3); // Shared vertex 2
          KineticVertex v_opp_t2 = t2.getVertex(edgeIndexInT2); // Vertex in t2 opposite the shared edge

           // --- Perform the flip ---
           // 1. Update vertices of the triangles
           t1.vertices[(edgeIndexInT1 + 1) % 3] = v_opp_t2; // Replace shared1 with v_opp_t2 in t1
           t2.vertices[(edgeIndexInT2 + 1) % 3] = v_opp_t1; // Replace shared1 with v_opp_t1 in t2
           // Ensure CCW order might need checks/swaps after vertex change

           // 2. Update edges and neighbors - this is complex
           // The old shared edge is gone. A new shared edge (v_opp_t1 - v_opp_t2) exists.
           // Edges/neighbors opposite the vertices need reassignment between t1 and t2.

           // Simplified: Just mark triangles for event recalculation
           System.out.println("  Performed flip between T" + t1.id + " and T" + t2.id + " (Connectivity update simplified)");
           eventQueue.update(t1, time);
           eventQueue.update(t2, time);
           // A full implementation requires careful reassignment of all 4 edge/neighbor pointers
           // for t1 and t2 based on the new vertex configuration.
      }


     // Helper to update vertex references in triangles around a collapsed vertex
     private void updateVertexReferences(KineticVertex oldVertex, KineticVertex newVertex, double time) {
         // This requires iterating through triangles incident to oldVertex.
         // The KineticVertex itself doesn't store incident triangles.
         // We need to iterate through the *global* list of active triangles.
         for (KineticTriangle tri : triangles) {
             if (!tri.hasStopped()) {
                 for (int i = 0; i < 3; i++) {
                     if (tri.getVertex(i) == oldVertex) {
                         tri.vertices[i] = newVertex; // Direct update (needs vertices to be non-final)
                         eventQueue.update(tri, time); // Mark for update
                         break; // Assume vertex appears only once per triangle
                     }
                 }
             }
         }
     }

      // Helper to mark affected triangles for event updates
      private void updateAffectedEvents(KineticTriangle sourceTriangle, double time) {
          if (!sourceTriangle.hasStopped()) {
              eventQueue.update(sourceTriangle, time);
          }
          for (int i = 0; i < 3; i++) {
              KineticTriangle neighbor = sourceTriangle.getNeighbor(i);
              if (neighbor != null && !neighbor.hasStopped()) {
                  eventQueue.update(neighbor, time);
              }
          }
      }

    // --- Getters ---
    public List<KineticTriangle> getActiveTriangles() {
        List<KineticTriangle> active = new ArrayList<>();
        for(KineticTriangle t : triangles) {
            if (!t.hasStopped()) active.add(t);
        }
        return active;
    }
     public List<KineticVertex> getVertices() { return vertices; }
     public List<KineticEdge> getEdges() { return edges; }


     // --- Output Generation ---
     public SkeletonOutput generateSkeletonOutput() {
         List<LineSegment> arcs = new ArrayList<>();
         // Iterate through vertices and construct segments/rays from prev/next links
         List<KineticVertex> processedStarts = new ArrayList<>();

         for (KineticVertex v : vertices) {
             if (v.hasStopped() && v.getNext() != null && !processedStarts.contains(v)) {
                  // Follow the chain from v
                  KineticVertex current = v;
                  while(current != null && current.getNext() != null && current.getNext() != v) { // Follow chain until stop or loop
                      KineticVertex next = current.getNext();
                      if (current.hasStopped() && next.hasStopped()) {
                          // Check if points are distinct enough to draw
                          if (current.getPosStop().distance(next.getPosStop()) > Constants.EPSILON) {
                             arcs.add(new LineSegment(current.getPosStop(), next.getPosStop()));
                             processedStarts.add(current); // Mark segment start as processed
                          }
                      } else if (current.hasStopped() && !next.hasStopped()) {
                          // Ray from current's stop point along next's velocity
                          // Need to decide on ray length or representation
                          System.out.println("Ray generation needed (not implemented)");
                      }
                      current = next;
                      if (processedStarts.contains(current)) break; // Avoid reprocessing part of a chain
                  }
             }
         }
          System.out.println("Generated " + arcs.size() + " skeleton arcs (simplified).");
         return new SkeletonOutput(arcs);
     }

}