package com.github.micycle1.surferj.collapse;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.github.micycle1.surferj.wavefront.WavefrontPropagator;

public class QuadraticSolver {

	private static final Logger LOGGER = LoggerFactory.getLogger(QuadraticSolver.class);

	// Tolerance for discriminant check and coefficient checks
	private static final double SOLVER_TOLERANCE = 1e-12; // Use a small tolerance

	/**
	 * Solves the quadratic equation p(x) = ax^2 + bx + c = 0. Returns ALL real
	 * roots, sorted in ascending order. Does NOT filter for non-negativity or time
	 * constraints.
	 *
	 * @param p The Polynomial (expected degree 2, but handles lower degrees)
	 * @return A sorted array of real roots (0, 1, or 2 elements).
	 */
	public static double[] solve(Polynomial p) {
		// Handle lower degrees explicitly if needed, though the quadratic formula
		// might degenerate correctly if 'a' is zero. However, the Polynomial
		// constructor ensures degree is set, so we rely on that.

		if (p.getDegree() < 2) {
			if (p.getDegree() == 1) { // bx + c = 0
				if (Math.abs(p.b) > SOLVER_TOLERANCE) {
					return new double[] { -p.c / p.b };
				} else {
					// 0x + c = 0. No single root unless c=0 (infinite roots).
					// Returning empty is safest for a solver context.
					// The calling function handles the c=0 case.
					return new double[0];
				}
			} else { // Degree 0: c = 0
				// Equation is just c = 0. No variable 'x'.
				// No specific root value. Caller handles this.
				return new double[0];
			}
		}

		// --- Degree 2 logic ---
		double a = p.a;
		double b = p.b;
		double c = p.c;

		// Should not happen if degree is 2, but check anyway
		if (Math.abs(a) <= SOLVER_TOLERANCE) {
			LOGGER.warn("Solver called with degree 2 polynomial, but 'a' coeff is near zero: {}", p);
			// Fallback to linear solver logic
			if (Math.abs(b) > SOLVER_TOLERANCE) {
				return new double[] { -c / b };
			} else {
				return new double[0];
			}
		}

		double discriminant = b * b - 4 * a * c;

		if (discriminant < -SOLVER_TOLERANCE) {
			// No real roots
			return new double[0];
		} else if (Math.abs(discriminant) <= SOLVER_TOLERANCE) {
			// One real root (discriminant is effectively zero)
			// Use stable formula, but -b/(2a) is standard
			double root = -b / (2 * a);
			return new double[] { root };
		} else {
			// Two distinct real roots
			double sqrtDiscriminant = Math.sqrt(discriminant);
			// Use the more stable quadratic formula variant to avoid cancellation errors
			double root1, root2;
			if (b >= 0) {
				root1 = (-b - sqrtDiscriminant) / (2 * a);
				root2 = (2 * c) / (-b - sqrtDiscriminant);
			} else {
				root1 = (2 * c) / (-b + sqrtDiscriminant);
				root2 = (-b + sqrtDiscriminant) / (2 * a);
			}

			// Sort roots
			if (root1 <= root2) {
				return new double[] { root1, root2 };
			} else {
				return new double[] { root2, root1 };
			}
		}
	}
}