package com.github.micycle.surferj.event;

import java.util.Objects;

import com.github.micycle.surferj.kinetics.KineticTriangle;

public class CollapseEvent implements Comparable<CollapseEvent> {

    private final EventType type;
    private final double time;
    private final KineticTriangle triangle; // The triangle detecting the event
    private int edgeIndex = -1; // Relevant edge index (for EDGE_COLLAPSE, FLIP_EVENT etc.)

    public CollapseEvent(EventType type, double time, KineticTriangle triangle) {
        this.type = Objects.requireNonNull(type);
        this.time = time;
        this.triangle = Objects.requireNonNull(triangle); // Can be null for NONE event? No, needs source.
    }

     // Constructor for events needing an edge index
     public CollapseEvent(EventType type, double time, KineticTriangle triangle, int edgeIndex) {
         this(type, time, triangle);
         this.edgeIndex = edgeIndex;
     }


    public EventType getType() {
        return type;
    }

    public double getTime() {
        return time;
    }

    public KineticTriangle getTriangle() {
        return triangle;
    }

     public int getEdgeIndex() { return edgeIndex; }
     public void setEdgeIndex(int index) { this.edgeIndex = index; } // If needed post-creation


     @Override
     public int compareTo(CollapseEvent other) {
         // 0. Component ID (if implemented)
         // int componentCompare = Integer.compare(this.componentId, other.componentId);
         // if (componentCompare != 0) return componentCompare;

         // 1. Time
         int timeCompare = Double.compare(this.time, other.time);
         if (timeCompare != 0) {
             return timeCompare;
         }

         // 2. Event Type Priority (Lower ordinal = HIGHER priority)
         int priorityThis = getEventTypePriority(this.type);
         int priorityOther = getEventTypePriority(other.type);
         int typePriorityCompare = Integer.compare(priorityThis, priorityOther);
         if (typePriorityCompare != 0) {
             return typePriorityCompare;
         }

         // 3. Flip Event Tie-breaking (if both are flips)
         if (this.type == EventType.FLIP_EVENT && other.type == EventType.FLIP_EVENT) {
              // Requires storing the length^2 of the edge being flipped (the spoke)
              // Higher length = HIGHER priority (so compare other.flipEdgeLengthSq to this.flipEdgeLengthSq)
        	 
        	 // NOTE UNCOMMENT BELOW!
//              int flipCompare = Double.compare(other.getFlipEdgeLengthSquared(), this.getFlipEdgeLengthSquared()); // Needs modification to CollapseEvent
//              if (flipCompare != 0) {
//                  return flipCompare;
//              }
         }

         // 4. Triangle ID (Final deterministic tie-breaking)
         return Integer.compare(this.triangle.id, other.triangle.id);
     }

     // Helper method for type priority (Lower value = higher priority)
     private int getEventTypePriority(EventType type) {
         switch (type) {
             // Highest priority - Infinite speed cases (add these enums if needed)
             // case INFINITE_SPEED_OPPOSING: return 0;
             // case INFINITE_SPEED_WEIGHTED: return 1;
             case EDGE_COLLAPSE:       return 10; // Real event
             case TRIANGLE_COLLAPSE:   return 11; // Real event (includes spoke collapse for now)
             // case SPOKE_COLLAPSE:   return 12; // Real event
             // case SPLIT_EVENT:      return 13; // Real event
             case FLIP_EVENT:          return 20; // Internal event
             // case CCW_VERTEX_LEAVES_CH: return 21; // Internal/Boundary event
             case NONE:                return 99; // Lowest priority
             default:                  return 50; // Other/Unknown
         }
     }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        CollapseEvent that = (CollapseEvent) o;
        // Use compareTo for equality check considering potential floating point issues
        return compareTo(that) == 0;
         // Direct comparison might fail due to floating point:
         // return Double.compare(that.time, time) == 0 && type == that.type && Objects.equals(triangle, that.triangle);
    }

    @Override
    public int hashCode() {
        // Hash code should be consistent with equals
        // Hashing doubles is tricky, use Long.doubleToLongBits
        return Objects.hash(type, Double.doubleToLongBits(time), triangle.id); // Use ID for triangle hash
    }

    @Override
    public String toString() {
        return "Event{" +
               "type=" + type +
               ", time=" + String.format("%.5f", time) +
               ", triangle=T" + triangle.id +
                (edgeIndex != -1 ? ", edgeIdx=" + edgeIndex : "") +
               '}';
    }
}