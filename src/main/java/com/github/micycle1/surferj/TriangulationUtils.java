package com.github.micycle1.surferj;

public class TriangulationUtils {

	// moving two steps forward (equivalent to one step backward in a 3-element
	// cycle)
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
}