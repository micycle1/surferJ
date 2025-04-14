package com.github.micycle.surferj2.kinetics;

import java.util.concurrent.atomic.AtomicLong;

import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.math.Vector2D;

import com.github.micycle.surferj2.SurfConstants;

public class WavefrontVertex {
	private static final AtomicLong idCounter = new AtomicLong(0);
	public final long id;
	public final Coordinate initialPosition; // Corresponds to pos_zero/pos_start
	public boolean isInfinite; // NOTE messy? with also having InfiniteSpeedType?

	// Kinetic properties
	private VertexAngle angle;
	private InfiniteSpeedType infiniteSpeed;
	private Coordinate velocity;

	// References to the two incident wavefront edges forming this vertex
	// Order matters: edge0 is CCW, edge1 is CW when looking from inside the
	// skeleton towards the vertex
	private WavefrontEdge edge0; // Left/CCW in C++ incident_wavefront_edges[0] (a)
	private WavefrontEdge edge1; // Right/CW in C++ incident_wavefront_edges[1] (b)

	private boolean hasStopped = false;
	private double timeStop = Double.NaN;
	private Coordinate posStop = null;

	// Links for DCEL structure built during propagation
	private WavefrontVertex nextVertex0 = null; // Next vertex along edge0 (CCW) boundary path
	private WavefrontVertex prevVertex0 = null; // Previous vertex along edge0 boundary path
	private WavefrontVertex nextVertex1 = null; // Next vertex along edge1 (CW) boundary path
	private WavefrontVertex prevVertex1 = null; // Previous vertex along edge1 boundary path

	// Placeholder for the infinite vertex
	public static final WavefrontVertex INFINITE_VERTEX = new WavefrontVertex();

	public enum VertexAngle {
		LEFT_TURN, // CONVEX
		RIGHT_TURN, // REFLEX
		COLLINEAR // STRAIGHT
	}

	public WavefrontEdge getWavefront(int id) {
		if (id == 0) {
			return edge0;
		} else if (id == 1) {
			return edge1;
		}
		System.err.println("Tried to get wavefront not in [0,1]!");
		return null;
	}

	/**
	 * This wavefront vertex is either between parallel, opposing wavefront elements
	 * that have crashed into each other and become collinear, or it is between
	 * neighboring wavefront edges that have become collinear yet have different
	 * weights.
	 */
	public enum InfiniteSpeedType {
		NONE, OPPOSING, WEIGHTED
	}

	private WavefrontVertex() {
		this.id = -1; // Special ID for infinite
		this.initialPosition = null;
		this.isInfinite = true;
		this.velocity = new Coordinate(0, 0);
		this.angle = VertexAngle.COLLINEAR;
		this.infiniteSpeed = InfiniteSpeedType.NONE;
	}

	// for tests
	public static WavefrontVertex makeInfinite(InfiniteSpeedType t) {
		WavefrontVertex v = new WavefrontVertex(new Coordinate(0, 0));
		v.infiniteSpeed = t;
		v.isInfinite = true;
		return v;
	}

	public WavefrontVertex(double x, double y, double vx, double vy) {
		// NOTE constructor for tests
		this.id = idCounter.incrementAndGet();
		this.initialPosition = new Coordinate(x, y);
		this.isInfinite = false;
		this.velocity = new Coordinate(vx, vy);
		this.angle = VertexAngle.COLLINEAR;
		this.infiniteSpeed = InfiniteSpeedType.NONE;
	}

	/**
	 * NOTE Unlike the main C++ constructor, these constructors do not perform any
	 * geometric calculations (orientation, get_infinite_speed_type,
	 * compute_velocity) based on incident edges because they don't receive them.
	 */
	public WavefrontVertex(Coordinate initialPosition) {
		this.id = idCounter.incrementAndGet();
		this.initialPosition = initialPosition;
		this.isInfinite = false;
		this.velocity = new Coordinate(0, 0);
		this.angle = VertexAngle.COLLINEAR;
		this.infiniteSpeed = InfiniteSpeedType.NONE;
	}

	public void setIncidentEdge(int index, WavefrontEdge edge) {
		if (isInfinite) {
			// Store reference even for infinite? Check C++ intent. Maybe not needed?
			if (index == 0) {
				this.edge0 = edge;
			} else if (index == 1) {
				this.edge1 = edge;
			} else {
				throw new IndexOutOfBoundsException("");
			}
			return;
		}

		// Only store the reference for finite vertices too
		if (index == 0) {
			this.edge0 = edge;
		} else if (index == 1) {
			this.edge1 = edge;
		} else {
			throw new IndexOutOfBoundsException("Index must be 0 or 1");
		}
	}

	/**
	 * Sets the incident wavefront edges for this vertex and recalculates its
	 * geometric properties (angle, infinite speed type, velocity). This method
	 * assumes it's called when both edges are known or intended to be set
	 * definitively (potentially null if on boundary without an edge).
	 *
	 * @param edge0 The edge arriving at this vertex (following CCW boundary).
	 * @param edge1 The edge leaving this vertex (following CCW boundary).
	 */
	public void setIncidentEdges(WavefrontEdge edge0, WavefrontEdge edge1) {
		if (this.isInfinite) {
			// Store references for topology if needed, but geometry is fixed.
			this.edge0 = edge0;
			this.edge1 = edge1;
			// Potentially throw exception if caller tries to set edges on infinite?
			// Or just return silently. For now, store and return.
			System.err.println("Warning: setIncidentEdges called on infinite vertex V" + id + ". Geometry not recalculated.");
			return;
		}

		// Store the references first, so helper methods can access them via `this`
		this.edge0 = edge0;
		this.edge1 = edge1;

		// Now recalculate geometry based on the stored edges
		// Use a try-catch block for robustness during calculations
		try {
			if (this.edge0 != null && this.edge1 != null) {
				// --- Optional: Perform Sanity Checks ---
				// Check if edges actually contain this vertex at the correct ends
				// (This helps catch logical errors in how edges are passed)
				if (!checkEdgeConsistency()) {
					// Log error and set defaults if consistency check fails
					System.err.println("Error: Edge consistency check failed for V" + id + ". Setting default geometry.");
					setDefaultGeometry();
					return; // Exit after setting defaults
				}
				// --- End Sanity Checks ---

				// Calculate geometric properties using helper methods
				calculateAngle(); // Calculates and sets this.angle
				calculateInfiniteSpeedType(); // Calculates and sets this.infiniteSpeed (needs this.angle)
				calculateVelocity(); // Calculates and sets this.velocity (needs this.angle, this.infiniteSpeed)

				// Optional: Log calculated values
				// System.out.println("V" + id + " geometry calculated: Angle=" + this.angle
				// + ", Speed=" + this.infiniteSpeed + ", Vel=" + this.velocity);

			} else {
				// One or both edges are null, set default geometry
				setDefaultGeometry();
			}
		} catch (Exception e) {
			// Catch potential errors during calculation (e.g., null pointers if checks
			// fail, math errors)
			System.err.println("Error during setIncidentEdges geometry calculation for V" + id + ": " + e.getMessage());
			e.printStackTrace(); // Or log properly
			setDefaultGeometry(); // Fallback to safe defaults on error
		}
	}

	public Coordinate getVelocity() {
		return velocity;
	}

	public double getVx() {
		return velocity.x;
	}

	public double getVy() {
		return velocity.y;
	}

	// Getter for infinite speed type
	public InfiniteSpeedType getInfiniteSpeed() {
		return infiniteSpeed;
	}

	public Coordinate getInitialPosition() {
		return initialPosition;
	}

	public boolean isInfinite() {
		return isInfinite;
	}

	/**
	 * Checks if the angle formed by the incident wavefront edges is convex (left
	 * turn) or straight (collinear). Conventionally, the infinite vertex is
	 * considered straight/convex.
	 *
	 * @return true if the vertex angle is CONVEX or STRAIGHT, false if REFLEX.
	 */
	public boolean isConvexOrStraight() {
		// The infinite vertex doesn't have a geometric angle in the same sense,
		// but by convention in geometric algorithms, it often doesn't participate
		// in reflex checks, so consider it "not reflex".
		if (this.isInfinite) {
			return true;
		}
		if (this.angle == null) {
			// This should not happen after constructor initialization
			System.err.println("Error: Vertex angle not initialized for V" + id + ". Returning default true.");
			return true;
		}
		return this.angle != VertexAngle.RIGHT_TURN;
	}

	/**
	 * Checks if the angle formed by the incident wavefront edges is reflex (right
	 * turn) or straight (collinear). Conventionally, the infinite vertex is
	 * considered straight/convex.
	 *
	 * @return true if the vertex angle is REFLEX or STRAIGHT, false if CONVEX.
	 */
	public boolean isReflexOrStraight() {
		// Similar to isConvexOrStraight, the infinite vertex is typically not
		// considered reflex.
		if (this.isInfinite) {
			return true; // Consistent with C++ (angle != CONVEX) default for null edges
		}
		if (this.angle == null) {
			// This should not happen after constructor initialization
			System.err.println("Error: Vertex angle not initialized for V" + id + ". Returning default true.");
			return true;
		}
		return this.angle != VertexAngle.LEFT_TURN;
	}

	public VertexAngle getAngle() {
		return angle;
	}

	// Compute position at time t (assuming motion starts at t=0)
	public Coordinate getPositionAt(double time) {
		if (isInfinite || infiniteSpeed != InfiniteSpeedType.NONE) {
			return initialPosition; // Does not move
		}
		return new Coordinate(initialPosition.getX() + velocity.getX() * time, initialPosition.getY() + velocity.getY() * time);
	}

	public boolean hasStopped() {
		return hasStopped;
	}

	public double getTimeStop() {
		return timeStop;
	}

	public Coordinate getPosStop() {
		return posStop;
	}

	public WavefrontEdge getIncidentEdge(int index) {
		if (index == 0)
			return edge0;
		if (index == 1)
			return edge1;
		throw new IndexOutOfBoundsException("Index must be 0 or 1");
	}

	/** Stops the vertex's motion at a specific time and position. */
	public void stop(double time, Coordinate position) {
		if (hasStopped) {
			// Allow stopping again if time and position are consistent (within tolerance)?
			if (Math.abs(time - this.timeStop) > 1e-9 || !position.equals2D(this.posStop, 1e-9)) {
				System.err.println("Warning: Vertex " + this.id + " already stopped at " + this.timeStop + "/" + this.posStop + ", attempting to stop again at "
						+ time + "/" + position);
				// Optionally throw exception or just ignore? Ignoring for now.
				return;
			}
		}
		if (isInfinite) {
			throw new IllegalStateException("Cannot stop the infinite vertex.");
		}
		this.timeStop = time;
		this.posStop = position; // Should calculate based on velocity if not provided? Assume provided for now.
		this.hasStopped = true;
	}

	// Overload for calculating position (Requires velocity implementation later)
	public void stop(double time) {
		if (isInfinite)
			throw new IllegalStateException("Cannot stop the infinite vertex.");
		// TODO: Calculate stop position based on initialPosition, velocity, and time
		// Coordinate calculatedPos = getPositionAt(time); // Need getPositionAt with
		// velocity
		// stop(time, calculatedPos);
		stop(time, this.initialPosition); // Placeholder: stops at initial position
		System.err.println("Warning: Vertex.stop(time) called without position - stopping at initialPosition. Velocity needed.");
	}

	// --- DCEL Linking Methods ---
	public void setNextVertex(int side, WavefrontVertex next) {
		if (side == 0)
			this.nextVertex0 = next;
		else if (side == 1)
			this.nextVertex1 = next;
		else
			throw new IndexOutOfBoundsException("Side must be 0 or 1");
	}

	public void setPrevVertex(int side, WavefrontVertex prev) {
		if (side == 0)
			this.prevVertex0 = prev;
		else if (side == 1)
			this.prevVertex1 = prev;
		else
			throw new IndexOutOfBoundsException("Side must be 0 or 1");
	}

	public WavefrontVertex getNextVertex(int side) {
		return side == 0 ? nextVertex0 : nextVertex1;
	}

	public WavefrontVertex getPrevVertex(int side) {
		return side == 0 ? prevVertex0 : prevVertex1;
	}

	/**
	 * Links the tail of this vertex (prev pointers) to the tail of another vertex.
	 * Used after splits.
	 */
	public void linkTailToTail(WavefrontVertex other) {
		if (this.prevVertex0 != null || other.prevVertex1 != null) {
			throw new IllegalStateException("Cannot link tails: previous pointers already set.");
		}
		this.prevVertex0 = other;
		other.prevVertex1 = this;
		System.out.println("DCEL Link T2T: " + this + "[0] <- " + other + "[1]");
	}

	/**
	 * Links the head of this vertex (next[side]) to the tail of another vertex
	 * (prev[side]).
	 */
	public void linkHeadToTail(int side, WavefrontVertex next) {
		if (side != 0 && side != 1)
			throw new IndexOutOfBoundsException("Side must be 0 or 1");
		if (getNextVertex(side) != null || next.getPrevVertex(side) != null) {
			System.err.println("Warning: Overwriting DCEL H2T links between " + this + " and " + next + " on side " + side);
			// throw new IllegalStateException("Cannot link head-to-tail: pointers already
			// set.");
		}
		setNextVertex(side, next);
		next.setPrevVertex(side, this);
		System.out.println("DCEL Link H2T: " + this + "[" + side + "] -> " + next + "[" + side + "]");
	}

	/**
	 * Links the head of this vertex (next[thisSide]) to the head of another vertex
	 * (next[otherSide]).
	 */
	public void linkHeadToHead(int thisSide, WavefrontVertex other, int otherSide) {
		if ((thisSide != 0 && thisSide != 1) || (otherSide != 0 && otherSide != 1))
			throw new IndexOutOfBoundsException("Side must be 0 or 1");
		if (getNextVertex(thisSide) != null || other.getNextVertex(otherSide) != null) {
			System.err.println("Warning: Overwriting DCEL H2H links between " + this + " and " + other);
			// throw new IllegalStateException("Cannot link head-to-head: pointers already
			// set.");
		}
		setNextVertex(thisSide, other);
		other.setNextVertex(otherSide, this);
		System.out.println("DCEL Link H2H: " + this + "[" + thisSide + "] <-> " + other + "[" + otherSide + "]");
	}

	/**
	 * (re)calculate vertex geometry once edges are added. This is not needed in C++
	 * version since edges are created first. In Java version, edges can be created
	 * afterwards.
	 */
	 public void recalculateGeometry() {
		// NOTE does this duplicate code?
		// Check if both edges are now known
		if (this.edge0 != null && this.edge1 != null) {
			// --- Logic copied/moved from the original setIncidentEdges(e0, e1) ---
			try { // Add try-catch for robustness during calculation
					// Step 1: Compute angle between edge directions
					// **** Use NORMAL DIRECTION consistently ****
				Vector2D dir0 = this.edge0.getSupportingLine().getNormalDirection();
				Vector2D dir1 = this.edge1.getSupportingLine().getNormalDirection();

				// Check for degenerate edges before calculation
				if (dir0.lengthSquared() < SurfConstants.ZERO_NORM_SQ || dir1.lengthSquared() < SurfConstants.ZERO_NORM_SQ) {
					System.err.println("Warning: Degenerate edge geometry detected for V" + id + ". Setting default angle/velocity.");
					setDefaultGeometry(); // Use the same defaults as the else block
					return;
				}

				final double det = dir0.getX() * dir1.getY() - dir0.getY() * dir1.getX();
				if (det > SurfConstants.ZERO_DET) { // Use tolerance
					this.angle = VertexAngle.LEFT_TURN;
				} else if (det < -SurfConstants.ZERO_DET) { // Use tolerance
					this.angle = VertexAngle.RIGHT_TURN;
				} else {
					this.angle = VertexAngle.COLLINEAR;
				}

				// Step 2: Determine infinite speed type (if applicable)
				// ... (use this.edge0, this.edge1) ...
				calculateInfiniteSpeedType(); // Extracted to helper?

				// Step 3: Compute velocity (if applicable)
				// ... (use this.edge0, this.edge1) ...
				calculateVelocity(); // Extracted to helper?

			} catch (Exception e) {
				System.err.println("Error during geometry recalculation for V" + id + ": " + e.getMessage());
				e.printStackTrace(); // Or log properly
				setDefaultGeometry(); // Fallback to safe defaults
			}
			// --- End of copied logic ---

		} else {
			// Default state if one or both edges are null
			setDefaultGeometry();
		}
	}

	private void setDefaultGeometry() {
		this.angle = VertexAngle.COLLINEAR; // Defaulting to COLLINEAR might cause issues like this test failure. Consider
											// VertexAngle.UNDEFINED?
		this.infiniteSpeed = InfiniteSpeedType.NONE;
		this.velocity = new Coordinate(0, 0); // Safe default velocity
	}

	/**
	 * Helper to check if the stored edge0 and edge1 correctly reference this
	 * vertex. Assumes edge0 is incoming CCW (ends at this vertex) and edge1 is
	 * outgoing CCW (starts at this vertex).
	 * 
	 * @return true if consistent, false otherwise.
	 */
	private boolean checkEdgeConsistency() {
		// Assumes edge0 and edge1 are non-null when called
		boolean edge0_ok = (this.edge0.getVertex(1) == this); // edge0 should end here
		boolean edge1_ok = (this.edge1.getVertex(0) == this); // edge1 should start here

		if (!edge0_ok) {
			System.err.println("Consistency Error: V" + id + " is not vertex 1 of its edge0 (" + this.edge0 + ")");
		}
		if (!edge1_ok) {
			System.err.println("Consistency Error: V" + id + " is not vertex 0 of its edge1 (" + this.edge1 + ")");
		}
		return edge0_ok && edge1_ok;
	}

	private void calculateAngle() {
		// Ensure both edges and their vertices are available
		if (edge0 == null || edge1 == null || edge0.getVertex(0) == null || edge0.getVertex(1) == null || // Check vertices of edge0
				edge1.getVertex(0) == null || edge1.getVertex(1) == null) { // Check vertices of edge1

			// Default to COLLINEAR if edge info is incomplete (matches C++ default)
			// Could also consider an UNDEFINED state if preferred
			this.angle = VertexAngle.COLLINEAR;
			return;
		}

		try {
			// Identify the three relevant points:
			// p0: The vertex preceding this one along edge0 (CCW incoming)
			// p1: This vertex (the common endpoint)
			// p2: The vertex following this one along edge1 (CCW outgoing)

			Coordinate p1 = this.initialPosition; // The vertex itself

			// Find p0: It's the "other" vertex of edge0
			WavefrontVertex v_p0 = (edge0.getVertex(0) == this) ? edge0.getVertex(1) : edge0.getVertex(0);
			if (v_p0 == null) {
				System.err.println("Warning: Could not find preceding vertex p0 for angle calculation at V" + id);
				this.angle = VertexAngle.COLLINEAR;
				return;
			}
			Coordinate p0 = v_p0.initialPosition;

			// Find p2: It's the "other" vertex of edge1
			WavefrontVertex v_p2 = (edge1.getVertex(0) == this) ? edge1.getVertex(1) : edge1.getVertex(0);
			if (v_p2 == null) {
				System.err.println("Warning: Could not find succeeding vertex p2 for angle calculation at V" + id);
				this.angle = VertexAngle.COLLINEAR;
				return;
			}
			Coordinate p2 = v_p2.initialPosition;

			// Calculate orientation using JTS: Orientation.index(p0, p1, p2)
			// p0 = incoming point, p1 = vertex, p2 = outgoing point
			int orientationIndex = Orientation.index(p0, p1, p2);

			// Map JTS orientation index to VertexAngle enum
			switch (orientationIndex) {
				case Orientation.COUNTERCLOCKWISE : // CCW turn from vector p0->p1 to p1->p2
					this.angle = VertexAngle.LEFT_TURN; // Corresponds to CGAL::LEFT_TURN
					break;
				case Orientation.CLOCKWISE : // CW turn from vector p0->p1 to p1->p2
					this.angle = VertexAngle.RIGHT_TURN; // Corresponds to CGAL::RIGHT_TURN
					break;
				case Orientation.COLLINEAR :
					this.angle = VertexAngle.COLLINEAR; // Corresponds to CGAL::COLLINEAR (STRAIGHT)
					break;
				default :
					// Should not happen with Orientation.index
					System.err.println("Warning: Unexpected JTS Orientation index: " + orientationIndex + " at V" + id);
					this.angle = VertexAngle.COLLINEAR; // Default on unexpected result
					break;
			}
		} catch (Exception e) {
			System.err.println("Error during JTS Orientation calculation for V" + id + ": " + e.getMessage());
			e.printStackTrace();
			this.angle = VertexAngle.COLLINEAR; // Default on error
		}
	}

	// Optional helper methods extracted from above
	private void calculateInfiniteSpeedType() {
		if (this.angle == VertexAngle.COLLINEAR) {
			// Use unit normals and weights as before
			Vector2D normal1 = this.edge1.getSupportingLine().getUnitNormal(); // Requires getUnitNormal()
			Vector2D dir0 = this.edge0.getSupportingLine().getNormalDirection();
			// Ensure unitNormal calculation is safe (handles zero length dir)
			if (normal1 == null || dir0 == null) { // Check if normals are valid
				this.infiniteSpeed = InfiniteSpeedType.NONE;
				return;
			}
			double orient = dir0.getX() * normal1.getY() - dir0.getY() * normal1.getX();
			if (orient > SurfConstants.ZERO_DET) { // Tolerance
				this.infiniteSpeed = InfiniteSpeedType.OPPOSING;
			} else if (Math.abs(this.edge0.getWeight() - this.edge1.getWeight()) > SurfConstants.ZERO_WEIGHT_DIFF) { // Tolerance
				this.infiniteSpeed = InfiniteSpeedType.WEIGHTED;
			} else {
				this.infiniteSpeed = InfiniteSpeedType.NONE;
			}
		} else {
			this.infiniteSpeed = InfiniteSpeedType.NONE;
		}
	}

	private void calculateVelocity() {
		if (this.infiniteSpeed == InfiniteSpeedType.NONE) {
			if (this.angle != VertexAngle.COLLINEAR) {
				WavefrontSupportingLine line0 = this.edge0.getSupportingLine();
				WavefrontSupportingLine line1 = this.edge1.getSupportingLine();
				Coordinate intersect = line0.lineAtOne().lineIntersection(line1.lineAtOne()); // Intersection at t=1
				if (intersect != null) {
					// Velocity is vector from initial pos to intersection at t=1
					this.velocity = new Coordinate(intersect.getX() - initialPosition.getX(), intersect.getY() - initialPosition.getY());
				} else {
					// Parallel lines? Should have been caught by COLLINEAR? Error.
					System.err.println("Error: No intersection for non-collinear vertex V" + id + " at time 1. Setting zero velocity.");
					this.velocity = new Coordinate(0, 0);
				}
			} else { // COLLINEAR and NONE infinite speed (implies same direction, same weight)
				// Velocity is along the normal (weighted normal?) - check C++ intent
				// Using unit normal * weight seems consistent with wavefront propagation
				Vector2D weightedNormal = this.edge0.getSupportingLine().getNormal(); // Get weighted normal
				if (weightedNormal == null) { // Handle degenerate case
					System.err.println("Error: Cannot get weighted normal for collinear vertex V" + id + ". Setting zero velocity.");
					this.velocity = new Coordinate(0, 0);
				} else {
					// Velocity is w * unit_normal. Since getNormal() = w * unit_normal, this is
					// just getNormal()
					this.velocity = new Coordinate(weightedNormal.getX(), weightedNormal.getY());
				}
			}
		} else { // Infinite speed
			this.velocity = new Coordinate(0, 0); // Finite velocity is zero
		}
	}

	@Override
	public String toString() {
		if (isInfinite)
			return "WV(Inf)";
		return "WV" + id + "(" + initialPosition.getX() + "," + initialPosition.getY() + ")";
	}
	// equals/hashCode based on ID for identity semantics if needed,
	// but using object identity is often sufficient.
	// Coordinate equality is based on value.
}