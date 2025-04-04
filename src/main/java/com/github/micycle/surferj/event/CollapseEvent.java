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
        // Primary sort: time
        int timeCompare = Double.compare(this.time, other.time);
        if (timeCompare != 0) {
            return timeCompare;
        }

        // Secondary sort: event type priority (lower enum ordinal = higher priority)
        // Example: Edge collapse happens "before" triangle collapse if times are equal
        int typeCompare = this.type.ordinal() - other.type.ordinal();
        if (typeCompare != 0) {
            return typeCompare;
        }

        // Tertiary sort: triangle ID (for deterministic tie-breaking)
        return Integer.compare(this.triangle.id, other.triangle.id);
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