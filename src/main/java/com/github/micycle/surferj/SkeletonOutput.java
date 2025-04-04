package com.github.micycle.surferj;

import org.locationtech.jts.geom.LineSegment;
import java.util.List;

public class SkeletonOutput {
    private final List<LineSegment> skeletonArcs;

    public SkeletonOutput(List<LineSegment> skeletonArcs) {
        this.skeletonArcs = List.copyOf(skeletonArcs); // Make immutable
    }

    public List<LineSegment> getSkeletonArcs() {
        return skeletonArcs;
    }

    // TODO: Add methods for outputting in different formats (WKT, GeoJSON, etc.)
}
