package com.github.micycle1.surferj.collapse;

// --- Polynomial Class (Simplified for Quadratic) ---
public class Polynomial {
	// Coefficients for ax^2 + bx + c
	public final double a, b, c;
	final int degree;

	// Represent ax^2 + bx + c
	public Polynomial(double a, double b, double c) {
		// Determine actual degree based on non-zero coefficients
		if (Math.abs(a) > 1e-12) {
			this.a = a;
			this.b = b;
			this.c = c;
			this.degree = 2;
		} else if (Math.abs(b) > 1e-12) {
			this.a = 0;
			this.b = b;
			this.c = c;
			this.degree = 1;
		} else {
			this.a = 0;
			this.b = 0;
			this.c = c;
			this.degree = 0;
		}
	}

	// Represent bx + c
	public Polynomial(double b, double c) {
		this(0, b, c);
	}

	// Represent c
	public Polynomial(double c) {
		this(0, 0, c);
	}

	public double evaluate(double t) {
		return a * t * t + b * t + c;
	}

	public Polynomial differentiate() {
		// Derivative of ax^2 + bx + c is 2ax + b
		return new Polynomial(2 * a, b); // Creates a linear polynomial
	}

	public int getDegree() {
		return degree;
	}

	// Sign of the leading non-zero coefficient
	public int sign() {
		if (degree == 2) {
			return Double.compare(a, 0.0);
		}
		if (degree == 1) {
			return Double.compare(b, 0.0);
		}
		return Double.compare(c, 0.0);
	}

	@Override
	public String toString() {
		return String.format("%.4fx² + %.4fx + %.4f (deg %d)", a, b, c, degree);
	}
}