package com.github.micycle1.surferj.collapse;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

// --- Quadratic Solver ---
public class QuadraticSolver {
	/**
	 * Solves the quadratic equation ax^2 + bx + c = 0.
	 *
	 * @param p The quadratic polynomial.
	 * @return A double array: - empty if no real roots. - {root} if one real root
	 *         (multiplicity 2). - {root1, root2} if two distinct real roots, sorted
	 *         root1 <= root2. Returns roots >= -tolerance.
	 */
	public static double[] solve(Polynomial p) {
		double tolerance = 1e-9; // Tolerance for floating point comparisons
		if (p.getDegree() < 2) {
			if (p.getDegree() == 1) { // bx + c = 0
				if (Math.abs(p.b) > tolerance) {
					double root = -p.c / p.b;
					if (root >= -tolerance) {
						return new double[] { Math.max(0.0, root) }; // Clamp small negatives
					} else {
						return new double[0]; // Root is negative
					}
				} else { // 0x + c = 0
					return (Math.abs(p.c) <= tolerance) ? new double[] { 0.0 } : new double[0]; // Or infinite roots? For events, 0 is fine.
				}
			} else { // c = 0
				return (Math.abs(p.c) <= tolerance) ? new double[] { 0.0 } : new double[0];
			}
		}

		double a = p.a;
		double b = p.b;
		double c = p.c;
		double discriminant = b * b - 4 * a * c;

		if (discriminant < -tolerance) { // Allow small negatives due to precision
			return new double[0]; // No real roots
		} else if (Math.abs(discriminant) <= tolerance) { // One real root (multiplicity 2)
			double root = -b / (2 * a);
			if (root >= -tolerance) {
				return new double[] { Math.max(0.0, root) };
			} else {
				return new double[0]; // Root is negative
			}
		} else { // Two distinct real roots
			double sqrtDiscriminant = Math.sqrt(discriminant);
			double root1 = (-b - sqrtDiscriminant) / (2 * a);
			double root2 = (-b + sqrtDiscriminant) / (2 * a);

			// Filter roots less than -tolerance and clamp small negatives
			List<Double> validRoots = new ArrayList<>();
			if (root1 >= -tolerance) {
				validRoots.add(Math.max(0.0, root1));
			}
			if (root2 >= -tolerance) {
				validRoots.add(Math.max(0.0, root2));
			}

			// Sort and return
			Collections.sort(validRoots);
			return validRoots.stream().mapToDouble(d -> d).toArray();
		}
	}
}