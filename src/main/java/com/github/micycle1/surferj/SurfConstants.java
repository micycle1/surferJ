package com.github.micycle1.surferj;

public class SurfConstants {

	public static final double ZERO_NT = 1e-12;
	public static final double ZERO_NT_SQ = 1e-12;
	/**
	 * Set to -1 for slower edge wins. Changing this during the propagation will
	 * yield undefined behaviour as the event queue is not correct then.
	 */
	public static final double FASTER_EDGE_WINS_IN_COLLINEAR_CASES = 1.0;
	public static final double ZERO_AREA_SQ = 1e-10;
	// threshold by which polynomial coeffcieitns is considered exactly 0.
	public static final double ZERO_POLYNOMIAL_COEFF = 1e-10;
	public static final double ZERO_POLYNOMIAL_VALUE = 1e-10;
	public static final double ZERO_NORM_SQ = 1e-12;
	public static final double ZERO_DIST = 1e-10;
	public static final double ZERO_DIST_SQ = 1e-12;
	public static final double ZERO_SPEED = 1e-10;
	public static final double ZERO_SPEED_SQ = 1e-10;
	public static final double ZERO_AREA = 1e-12; // an appropriate tolerance for cross product (area)
	public static final double VERTICAL_SLOPE_TOL = 1e-12;
	public static final double TIME_TOL = 1e-12;
	public static final double ZERO_DET = 1e-12;
	public static final double ZERO_WEIGHT_DIFF = 1e-12;
}
