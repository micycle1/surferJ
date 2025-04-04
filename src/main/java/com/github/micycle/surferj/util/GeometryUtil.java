package com.github.micycle.surferj.util;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.math.Vector2D;

public class GeometryUtil {

	/**
	 * Calculates the velocity vector for a vertex based on its two incident
	 * wavefront edges, assuming the vertex is the common endpoint. The velocity
	 * vector V satisfies V . nA = weightA and V . nB = weightB, where nA and nB are
	 * the unit normals of the edges pointing outwards from the polygon boundary.
	 *
	 * @param vertexPos The position of the vertex (at t=0 for initial calculation).
	 * @param edgeA     The first incident edge segment. Assumed orientation: p0 ->
	 *                  vertexPos.
	 * @param edgeB     The second incident edge segment. Assumed orientation:
	 *                  vertexPos -> p1.
	 * @param weightA   The weight (speed) associated with edgeA.
	 * @param weightB   The weight (speed) associated with edgeB.
	 * @return The calculated velocity vector. Handles collinear cases. Returns zero
	 *         vector for opposing edges.
	 */
	public static Vector2D calculateVelocity(Coordinate vertexPos, LineSegment edgeA, LineSegment edgeB, double weightA, double weightB) {
		// Ensure edges are valid
		if (edgeA.getLength() < Constants.EPSILON || edgeB.getLength() < Constants.EPSILON) {
			System.err.println("Warning: Zero-length edge provided to calculateVelocity.");
			return Vector2D.create(0, 0);
		}

		// 1. Calculate outward unit normals (assuming CCW polygon boundary)
		// Normal = CCW 90-degree rotation of the edge vector along the boundary.

		// Vector along edge A (p0 -> vertexPos)
		Vector2D vecA_edge = Vector2D.create(vertexPos.x - edgeA.p0.x, vertexPos.y - edgeA.p0.y);
		Vector2D normA = vecA_edge.rotate(Math.PI / 2.0).normalize(); // Outward normal for A

		// Vector along edge B (vertexPos -> p1)
		Vector2D vecB_edge = Vector2D.create(edgeB.p1.x - vertexPos.x, edgeB.p1.y - vertexPos.y);
		Vector2D normB = vecB_edge.rotate(Math.PI / 2.0).normalize(); // Outward normal for B

		// Check for NaN result from normalization (if edge vector was zero length,
		// though checked above)
		if (Double.isNaN(normA.getX()) || Double.isNaN(normB.getX())) {
			System.err.println("Error: Normalization failed in calculateVelocity (unexpected zero vector).");
			return Vector2D.create(0, 0);
		}

		// 2. Solve the system V.nA = wA, V.nB = wB for V=(Vx, Vy)
		// Vx*normA.x + Vy*normA.y = wA
		// Vx*normB.x + Vy*normB.y = wB

		double ax = normA.getX();
		double ay = normA.getY();
		double bx = normB.getX();
		double by = normB.getY();

		// Determinant D = ax*by - ay*bx ( = normA x normB, related to sin(angle between
		// normals))
		double D = ax * by - ay * bx;

		// 3. Handle Collinear Cases (D == 0)
		if (Math.abs(D) < Constants.EPSILON) {
			// Normals are parallel. Edges are collinear.
			// Check if normals are same or opposite direction.
			if (normA.distance(normB) < Constants.EPSILON) {
				// Same direction (e.g., straight vertex)
				if (Math.abs(weightA - weightB) > Constants.EPSILON) {
					// WEIGHTED infinite speed case
					System.err.println("Warning: Collinear edges with different weights detected (Weighted Infinite Speed).");
					// Return velocity based on faster edge
					return (weightA > weightB) ? normA.multiply(weightA) : normB.multiply(weightB);
				} else {
					// Same direction, same weight. Velocity is simply normal * weight.
					return normA.multiply(weightA);
				}
			} else if (normA.distance(normB.multiply(-1.0)) < Constants.EPSILON) {
				// Opposite direction (e.g., reflex 360 or cusp)
				System.err.println("Warning: Opposing collinear edges detected (Opposing Infinite Speed).");
				// This case requires specific event handling. Return zero velocity for now.
				return Vector2D.create(0, 0);
			} else {
				// Should not happen if D is near zero unless vectors are zero
				System.err.println("Error: Near-zero determinant in calculateVelocity but normals neither same nor opposite.");
				return Vector2D.create(0, 0);
			}
		}

		// 4. Solve using Cramer's Rule for non-collinear case
		double vx = (weightA * by - weightB * ay) / D;
		double vy = (weightB * ax - weightA * bx) / D; // Note: Original formula had swapped terms for Vy

		return Vector2D.create(vx, vy);
	}

	// ... (rest of GeometryUtil methods: lineIntersection, solveQuadratic,
	// orientation, angleBisector)
	// Make sure solveQuadratic is the corrected version from the previous step.
	public static Coordinate lineIntersection(LineSegment l1, LineSegment l2) {
		double x1 = l1.p0.x, y1 = l1.p0.y;
		double x2 = l1.p1.x, y2 = l1.p1.y;
		double x3 = l2.p0.x, y3 = l2.p0.y;
		double x4 = l2.p1.x, y4 = l2.p1.y;

		double denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

		if (Math.abs(denominator) < Constants.EPSILON) {
			return null; // Parallel or coincident
		}

		double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denominator;

		double intersectX = x1 + t * (x2 - x1);
		double intersectY = y1 + t * (y2 - y1);

		return new Coordinate(intersectX, intersectY);
	}

	public static double solveQuadratic(double a, double b, double c) {
		if (Math.abs(a) < Constants.EPSILON) { // Linear equation bx + c = 0
			if (Math.abs(b) < Constants.EPSILON) {
				return (Math.abs(c) < Constants.EPSILON) ? 0.0 : Double.POSITIVE_INFINITY;
			}
			double root = -c / b;
			return (root >= -Constants.EPSILON) ? Math.max(0.0, root) : Double.POSITIVE_INFINITY;
		}

		double discriminant = b * b - 4 * a * c;
		if (discriminant < -Constants.EPSILON) {
			return Double.POSITIVE_INFINITY; // No real roots
		}
		if (discriminant < 0)
			discriminant = 0; // Clamp near-zero negative

		double sqrtDiscriminant = Math.sqrt(discriminant);
		double root1 = (-b + sqrtDiscriminant) / (2 * a);
		double root2 = (-b - sqrtDiscriminant) / (2 * a);

		double minPositiveRoot = Double.POSITIVE_INFINITY;
		if (root1 >= -Constants.EPSILON) {
			minPositiveRoot = Math.max(0.0, root1);
		}
		if (root2 >= -Constants.EPSILON) {
			double clampedRoot2 = Math.max(0.0, root2);
			minPositiveRoot = Math.min(minPositiveRoot, clampedRoot2);
		}

		return minPositiveRoot;
	}

	// Calculate the bisector of the angle between two vectors (normalized)
	// Note: This simple addition works for the bisector direction but might need
	// adjustments depending on the exact definition needed (e.g., handling 180
	// degrees)
	public static Vector2D angleBisector(Vector2D v1, Vector2D v2) {
		Vector2D n1 = v1.normalize();
		Vector2D n2 = v2.normalize();
		// Ensure vectors are valid before adding
		if (Double.isNaN(n1.getX()) || Double.isNaN(n2.getX())) {
			System.err.println("Warning: NaN detected in angleBisector input normalization.");
			// Handle degenerate case, e.g., return one of the vectors or zero
			return Double.isNaN(n1.getX()) ? (Double.isNaN(n2.getX()) ? Vector2D.create(0, 0) : n2) : n1;
		}
		// Check if vectors are opposite
		if (n1.distance(n2.multiply(-1.0)) < Constants.EPSILON) {
			System.err.println("Warning: Trying to bisect opposing vectors.");
			// Bisector is ambiguous. Return a perpendicular vector or handle as error.
			// Returning a perpendicular for now.
			return n1.rotate(Math.PI / 2.0);
		}

		return n1.add(n2).normalize(); // Simple average for bisector direction
	}

	// Determine orientation: > 0 CCW, < 0 CW, = 0 Collinear
	public static double orientation(Coordinate p, Coordinate q, Coordinate r) {
		// Using JTS robust orientation is recommended if precision becomes an issue
		// return Orientation.index(p, q, r); // Returns -1, 0, 1
		// Manual calculation for double return:
		return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
	}
}