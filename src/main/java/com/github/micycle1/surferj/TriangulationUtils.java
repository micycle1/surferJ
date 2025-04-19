package com.github.micycle1.surferj;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;

public class TriangulationUtils {

	// moving two steps forward (equivalent to one step backward in a 3-element
	// cycle)

	// faster to use tables?
//	private static final int[] CW = { 2, 0, 1 };
//	private static final int[] CCW = { 1, 2, 0 };

	public static int cw(int i) {
		return (i + 2) % 3;
	}

	// move one step forward (CCW)
	public static int ccw(int i) {
		return (i + 1) % 3;
	}

	public static int mod3(int i) {
		return ((i % 3) + 3) % 3;
	}

	public static boolean pointLiesOnSegment(Coordinate p1, Coordinate p2, Coordinate q) {
		if (Orientation.index(p1, p2, q) == Orientation.COLLINEAR) {
			// If collinear, check bounds
			if (p1.x != p2.x) { // Non-vertical segment
				// Using min/max to handle segments in either direction
				return (Math.min(p1.x, p2.x) <= q.x) && (q.x <= Math.max(p1.x, p2.x));
			} else { // Vertical segment
				return (Math.min(p1.y, p2.y) <= q.y) && (q.y <= Math.max(p1.y, p2.y));
			}
		} else {
			return false;
		}
	}
}