package com.github.micycle1.surferj.kinetics;

import java.util.List;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicLong;

import org.locationtech.jts.algorithm.LineIntersector;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.math.Vector2D;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.github.micycle1.surferj.SurfConstants;
import com.github.micycle1.surferj.collapse.Polynomial;
import com.github.micycle1.surferj.collapse.QuadraticSolver;
import com.github.micycle1.surferj.wavefront.WavefrontPropagator;

public class WavefrontVertex {

	private static final Logger LOGGER = LoggerFactory.getLogger(WavefrontVertex.class);

	private static final AtomicLong idCounter = new AtomicLong(0);
	private static final Coordinate ORIGIN = new Coordinate(0, 0);
	public final long id;
	public final Coordinate posStart; // Corresponds to pos_zero/pos_start
	public Coordinate posZero; // Corresponds to pos_zero/pos_start  NOTE required?
	public boolean isInfinite; // NOTE messy? with also having InfiniteSpeedType?
	public boolean isDegenerate = false;

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
	private double timeStart = 0;
	private Coordinate posStop = null;
	
	private boolean isInitial = false;

	// Links for DCEL structure built during propagation
	WavefrontVertex[] nextVertex = new WavefrontVertex[2];
	WavefrontVertex[] prevVertex = new WavefrontVertex[2];

	private Polynomial px;
	private Polynomial py;

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
		LOGGER.error("Tried to get wavefront not in [0,1]!");
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
		this.posStart = null;
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
		this.posStart = new Coordinate(x, y);
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
		this.posStart = initialPosition;
		this.isInfinite = false;
		this.velocity = new Coordinate(0, 0);
		this.angle = VertexAngle.COLLINEAR;
		this.infiniteSpeed = InfiniteSpeedType.NONE;
	}

	// --- Constructor used by makeVertex ---
	/**
	 * Private constructor to create a vertex with known initial position and
	 * potentially stopping information if created during an event.
	 *
	 * @param initialPos Initial position (at time=0).
	 * @param stopPos    Position where the vertex stops (often the creation
	 *                   position).
	 * @param stopTime   Time when the vertex stops (often the creation time).
	 * @param edge0      Incident edge 0.
	 * @param edge1      Incident edge 1.
	 */
	private WavefrontVertex(Coordinate posZero, Coordinate posStart, double timeStart, WavefrontEdge edge0, WavefrontEdge edge1) {
	    this.id = idCounter.incrementAndGet();
	    this.posZero = Objects.requireNonNull(posZero, "posZero cannot be null"); // Store pos_zero equivalent
	    this.posStart = Objects.requireNonNull(posStart, "posStart cannot be null"); // Store pos_start equivalent
	    this.timeStart = timeStart; // Store time_start equivalent

	    // is_initial should be false for vertices created via makeVertex
	    this.isInitial = false;
//	    this.isBeveling = false; // makeVertex doesn't handle beveling directly
	    this.isInfinite = false; // makeVertex is for finite vertices

	    // Set incident edges AND CALCULATE angle, infiniteSpeed, velocity
	    // Make sure setIncidentEdges uses posZero for geometric velocity calculation
        this.edge0 = edge0;
        this.edge1 = edge1;
        // actually Initialize geometry to defaults, will be calculated later
        setDefaultGeometry(); // Sets angle=COLLINEAR, infiniteSpeed=NONE, velocity=(0,0)

	    // --- Initialize state correctly ---
	    this.hasStopped = false; // Vertex just started, not stopped
	    this.timeStop = Double.NaN;
	    this.posStop = null;
	    this.isDegenerate = false; // Cannot be degenerate at creation time unless timeStart=0 and posStart=posZero

	    // --- C++ degeneracy check at creation (only relevant if timeStart == 0) ---
	    if (Math.abs(this.timeStart) < SurfConstants.TIME_TOL && this.posStart.equals2D(this.posZero, SurfConstants.ZERO_DIST_SQ)) {
	        LOGGER.warn("Vertex {} created degenerate (timeStart=0, posStart=posZero)", this.id);
	        // C++ doesn't explicitly set is_degenerate_ here, it happens in stop()
	        // but good to be aware.
	    }

	    // Initialize links
	    this.nextVertex[0] = this.nextVertex[1] = null;
	    this.prevVertex[0] = this.prevVertex[1] = null;

	    // Calculate polynomial representations AFTER velocity is known
	    if (this.infiniteSpeed == InfiniteSpeedType.NONE && this.velocity != null) {
	        this.px = new Polynomial(this.velocity.getX(), this.posZero.x); // Assuming Polynomial(slope, intercept)
	        this.py = new Polynomial(this.velocity.getY(), this.posZero.y);
	    } else {
	        // For infinite speed, polynomials are constant at the creation point?
	        // Or maybe just represent initial position? C++ needs review here.
	        // Let's make them constant at the *initial* position (posZero).
	         this.px = new Polynomial(0, this.posZero.x);
	         this.py = new Polynomial(0, this.posZero.y);
	         LOGGER.warn("Vertex {} created with infinite speed {}. Polynomial representation might be simplified.", this.id, this.infiniteSpeed);
	    }
	}

	/**
	 * Creates a new WavefrontVertex instance resulting from an event (like split or
	 * constraint collapse) occurring at a specific time and position. Calculates
	 * the initial position (time=0) by back-propagating.
	 *
	 * @param pos       The known position of the vertex at 'time'.
	 * @param time      The time at which the vertex is created/known to be at
	 *                  'pos'.
	 * @param edgeA     The first incident wavefront edge.
	 * @param edgeB     The second incident wavefront edge.
	 * @param fromSplit Flag indicating if this vertex results from a split event
	 *                  (C++ uses this for assertions, maybe logging in Java).
	 * @return The newly created WavefrontVertex.
	 */
	public static WavefrontVertex makeVertex(Coordinate pos, double time, WavefrontEdge edgeA, WavefrontEdge edgeB, boolean fromSplit,
			// TODO check this method
			List<WavefrontVertex> vertices) {
		if (Double.isNaN(time) || time < -SurfConstants.TIME_TOL) { // Use tolerance
			throw new IllegalArgumentException("Invalid creation time: " + time);
		}
		// Clamp time >= 0
//		double creationTime = Math.max(0.0, time);

		// C++ Assertion check
		if (!fromSplit) {
			// LOG(INFO) or DEBUG check if needed. C++ has empty block here.
		}

		Coordinate posZero; // The initial position (t=0) we need to calculate
		WavefrontSupportingLine lineA = edgeA.getSupportingLine(); // Assuming getter exists
		WavefrontSupportingLine lineB = edgeB.getSupportingLine(); // Assuming getter exists

		// Compute intersection of the supporting lines
		var intersectionPoint = lineA.lineIntersection(lineB); // Assuming method exists
		int intersection_type = LineIntersector.NO_INTERSECTION; // NONE
		if (intersectionPoint != null) {
			intersection_type = LineIntersector.POINT_INTERSECTION;
		} else {
			boolean collinear = lineA.collinear(lineB);
			if (collinear) {
				intersection_type = LineIntersector.COLLINEAR_INTERSECTION;
			}
		}

		switch (intersection_type) {
			case LineIntersector.COLLINEAR_INTERSECTION : // Lines are identical
				// Back-propagate using velocity. Need to compute velocity first.
				// Velocity depends on the angle type at the *intersection point*.
				// We need to determine angle type (CONVEX, REFLEX, STRAIGHT) based on edge
				// directions relative to each other.
//				VertexAngle angleAtIntersection = calculateAngle(edgeA, edgeB);
				Coordinate v = calculateVelocity(edgeA, edgeB, VertexAngle.COLLINEAR, InfiniteSpeedType.NONE, ORIGIN);
				var velocity = Vector2D.create(v);
				// pos_zero = pos - time * velocity;
				posZero = Vector2D.create(pos).subtract(velocity.multiply(time)).toCoordinate();
				break;

			case LineIntersector.POINT_INTERSECTION : // Lines intersect at a single point (at t=0)
				// The intersection point *is* the initial position (t=0).
				posZero = intersectionPoint;
				// Optional: Assert that pos ≈ intersectionPoint + time * velocity
				// Vec2d checkVel = computeVelocity(intersectionPoint, lineA, lineB,
				// determineAngleType(lineA, lineB));
				// Vec2d expectedPos = intersectionPoint.add(checkVel.scale(creationTime));
				// assert pos.distanceSq(expectedPos) < SurfConstants.ZERO_TOL_SQ : "Position
				// mismatch after back-propagation";
				break;

			case LineIntersector.NO_INTERSECTION : // Lines are parallel (never intersect)
				// This case might indicate an issue or a special scenario (e.g., straight
				// vertex)
				// C++ sets pos_zero = pos. This implies velocity is zero or lines shouldn't be
				// parallel?
				// If lines are parallel, the vertex shouldn't exist unless it's infinitely far
				// away or angle is straight.
				// Let's follow C++ for now, but this might need review.
//                 LOGGER.warning("Parallel lines encountered in makeVertex. Setting posZero = current pos. Verify logic.");
				LOGGER.warn("Parallel lines encountered in makeVertex. Setting posZero = current pos. Verify logic.");
				posZero = pos;
				break;

			default :
				throw new IllegalStateException("Unhandled LineIntersectionType: " + intersection_type);
		}

		// Create the vertex using the calculated posZero and the provided stopping info
		// (pos, time)
		WavefrontVertex v = new WavefrontVertex(posZero, pos, time, edgeA, edgeB);

		// Optional: Immediately link the new vertex to the edges (if constructor
		// doesn't)
		// edgeA.setVertex(?, v); // Need to know which end corresponds to edgeA/edgeB
		// edgeB.setVertex(?, v);

		vertices.add(v);
		return v;
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
			LOGGER.warn("setIncidentEdges called on infinite vertex V" + id + ". Geometry not recalculated.");
			return;
		}

		// Store the references first, so helper methods can access them via `this`
		this.edge0 = edge0;
		this.edge1 = edge1;

		// Now recalculate geometry based on the stored edges
		recalculateGeometry();
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
		return posStart;
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
			LOGGER.error("Vertex angle not initialized for V" + id + ". Returning default true.");
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
			LOGGER.error("Vertex angle not initialized for V" + id + ". Returning default true.");
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
			return posStart; // Does not move
		}
		// posStart.add(velocity.multiply(t))
		return new Coordinate(posStart.getX() + velocity.getX() * time, posStart.getY() + velocity.getY() * time);
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

	// incident_wavefront_edge(int) in C++
	public WavefrontEdge getIncidentEdge(int index) {
		if (index == 0) {
			return edge0;
		}
		if (index == 1) {
			return edge1;
		}
		throw new IndexOutOfBoundsException("Index must be 0 or 1");
	}

    // Existing stop method that takes position (likely okay)
    public void stop(double time, Coordinate position) {
        if (hasStopped) {
            return;
        }
        if (isInfinite) {
            throw new IllegalStateException("Cannot stop the infinite vertex.");
        }
         // Check for infinite speed - C++ uses this version ONLY for infinite speed
         if (getInfiniteSpeed() == InfiniteSpeedType.NONE) {
             LOGGER.warn("Vertex " + this.id + ": stop(time, position) called but vertex does not have infinite speed.");
             // Decide how to handle: maybe call the other stop method? Or log and proceed?
             // Let's proceed for now, assuming the caller knows best, but log it.
         }

        this.timeStop = time;
        this.posStop = position;
        this.hasStopped = true;

        // Check for degeneracy (C++ logic)
        if (this.posStart != null && this.posStop != null && this.posStop.equals2D(this.posStart, SurfConstants.ZERO_DIST_SQ)) {
             this.isDegenerate = true;
             // Optional: Add C++ assertion check: assert timeStop == timeStart within tolerance?
             if (Math.abs(this.timeStop - this.timeStart) > SurfConstants.TIME_TOL) {
                  LOGGER.warn("Vertex {} degenerate stop: timeStop {} != timeStart {}", this.id, this.timeStop, this.timeStart);
             }
        } else {
             this.isDegenerate = false;
        }
    }

    public void stop(double time) {
        if (hasStopped) {
             // ... (existing warning/check logic - check against calculated position?) ...
             // Recalculate intended stop position to compare with existing posStop
             Coordinate calculatedPos = getPositionAt(time);
             if (Math.abs(time - this.timeStop) > SurfConstants.TIME_TOL ||
                 (this.posStop != null && calculatedPos != null && !this.posStop.equals2D(calculatedPos, SurfConstants.ZERO_DIST_SQ))) {
                 LOGGER.warn("Vertex " + this.id + " already stopped at " + this.timeStop + "/" + this.posStop + ", attempting to stop again at " + time + " (calculated pos: " + calculatedPos + ")");
             }
            return;
        }
        if (isInfinite) {
            throw new IllegalStateException("Cannot stop the infinite vertex.");
        }
        // Check for infinite speed - C++ uses this version ONLY for NON-infinite speed
        if (getInfiniteSpeed() != InfiniteSpeedType.NONE) {
             LOGGER.error("Vertex " + this.id + ": stop(time) called but vertex has infinite speed type " + getInfiniteSpeed() + ". Use stop(time, position).");
             // Decide how to handle: throw exception? Or try to proceed?
             // Throwing is safer as the state is unexpected.
             throw new IllegalStateException("stop(time) called on vertex with infinite speed.");
        }

        this.timeStop = time;
        // *** Calculate the stopping position ***
        this.posStop = getPositionAt(time); // Assumes getPositionAt(time) is correct
        this.hasStopped = true;

         // Check for calculation failure
         if (this.posStop == null) {
             LOGGER.error("Vertex " + this.id + ": Failed to calculate stop position at time " + time);
             this.hasStopped = false; // Revert state
             this.timeStop = Double.NaN;
             // Optionally throw an exception
             throw new IllegalStateException("Failed to calculate stop position.");
         }


        // Check for degeneracy (C++ logic)
        if (this.posStart != null && this.posStop.equals2D(this.posStart, SurfConstants.ZERO_DIST_SQ)) {
             this.isDegenerate = true;
             // Optional: Add C++ assertion check: assert timeStop == timeStart within tolerance?
             if (Math.abs(this.timeStop - this.timeStart) > SurfConstants.TIME_TOL) {
                  LOGGER.warn("Vertex {} degenerate stop: timeStop {} != timeStart {}", this.id, this.timeStop, this.timeStart);
             }
        } else {
            this.isDegenerate = false;
        }
         LOGGER.debug("Vertex {} stopped via stop(time) at {} / {}", this.id, this.timeStop, this.posStop);
    }

	public void setNextVertex(int side, WavefrontVertex next) {
		setNextVertex(side, next, true);
	}

	public void setNextVertex(int side, WavefrontVertex next, boolean headToTail) {
		nextVertex[side] = next;
		if (headToTail) {
			next.prevVertex[side] = this;
		} else {
			next.nextVertex[1 - side] = this;
		}
	}

	/**
	 * Links the tail of this vertex (prev pointers) to the tail of another vertex.
	 * Used after splits.
	 */
	public void linkTailToTail(WavefrontVertex other) {
		prevVertex[0] = other;
		other.prevVertex[1] = this;
		LOGGER.debug("DCEL Link T2T: " + this + "[0] <- " + other + "[1]");
	}

	/**
	 * (re)calculate vertex geometry once edges are added. This is not needed in C++
	 * version since edges are created first. In Java version, edges can be created
	 * afterwards.
	 */
	public void recalculateGeometry() {
		// Check if both edges are now known
		if (this.edge0 != null && this.edge1 != null) {
			// --- Logic copied/moved from the original setIncidentEdges(e0, e1) ---
			try { // Add try-catch for robustness during calculation

				// --- Optional: Perform Sanity Checks ---
				// Check if edges actually contain this vertex at the correct ends
				// (This helps catch logical errors in how edges are passed)
				if (!checkEdgeConsistency()) {
					// Log error and set defaults if consistency check fails
					LOGGER.error("Edge consistency check failed for V" + id + ". Setting default geometry.");
					setDefaultGeometry();
					return; // Exit after setting defaults
				}
				// --- End Sanity Checks ---

				// Step 1: Compute angle between edge directions
				// **** Use NORMAL DIRECTION consistently ****
				Vector2D dir0 = this.edge0.getSupportingLine().getNormalDirection();
				Vector2D dir1 = this.edge1.getSupportingLine().getNormalDirection();

				// Check for degenerate edges before calculation
				if (dir0.lengthSquared() < SurfConstants.ZERO_NORM_SQ || dir1.lengthSquared() < SurfConstants.ZERO_NORM_SQ) {
					LOGGER.warn("Degenerate edge geometry detected for V" + id + ". Setting default angle/velocity.");
					setDefaultGeometry(); // Use the same defaults as the else block
					return;
				}

				// Calculate geometric properties using helper methods
				angle = calculateAngle(edge0, edge1); // Calculates and sets this.angle
				infiniteSpeed = calculateInfiniteSpeedType(edge0, edge1, angle); // Calculates and sets this.infiniteSpeed (needs this.angle)
				velocity = calculateVelocity(edge0, edge1, angle, infiniteSpeed, posStart);

			} catch (Exception e) {
				LOGGER.error("Error during geometry recalculation for V" + id + ": " + e.getMessage());
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
			LOGGER.error("Consistency V" + id + " is not vertex 1 of its edge0 (" + this.edge0 + ")");
		}
		if (!edge1_ok) {
			LOGGER.error("Consistency V" + id + " is not vertex 0 of its edge1 (" + this.edge1 + ")");
		}
		return edge0_ok && edge1_ok;
	}

	static VertexAngle calculateAngle(WavefrontEdge edge0, WavefrontEdge edge1) {
		// Ensure both edges and their vertices are available
		if (edge0 == null || edge1 == null || edge0.getVertex(0) == null || edge0.getVertex(1) == null || // Check vertices of edge0
				edge1.getVertex(0) == null || edge1.getVertex(1) == null) { // Check vertices of edge1
			return VertexAngle.COLLINEAR;
		}

		try {
			// Find the common vertex between edge0 and edge1
			WavefrontVertex commonVertex = null;
			if (edge0.getVertex(0) == edge1.getVertex(0)) {
				commonVertex = edge0.getVertex(0);
			} else if (edge0.getVertex(0) == edge1.getVertex(1)) {
				commonVertex = edge0.getVertex(0);
			} else if (edge0.getVertex(1) == edge1.getVertex(0)) {
				commonVertex = edge0.getVertex(1);
			} else if (edge0.getVertex(1) == edge1.getVertex(1)) {
				commonVertex = edge0.getVertex(1);
			}

			if (commonVertex == null) {
				LOGGER.warn("No common vertex found between edges");
				return VertexAngle.COLLINEAR;
			}

			Coordinate p1 = commonVertex.posStart; // The common vertex

			// Find p0: The "other" vertex of edge0
			WavefrontVertex v_p0 = (edge0.getVertex(0) == commonVertex) ? edge0.getVertex(1) : edge0.getVertex(0);
			if (v_p0 == null) {
				LOGGER.warn("Could not find preceding vertex p0 for angle calculation");
				return VertexAngle.COLLINEAR;
			}
			Coordinate p0 = v_p0.posStart;

			// Find p2: The "other" vertex of edge1
			WavefrontVertex v_p2 = (edge1.getVertex(0) == commonVertex) ? edge1.getVertex(1) : edge1.getVertex(0);
			if (v_p2 == null) {
				LOGGER.warn("Could not find succeeding vertex p2 for angle calculation");
				return VertexAngle.COLLINEAR;
			}
			Coordinate p2 = v_p2.posStart;

			// Calculate orientation using JTS: Orientation.index(p0, p1, p2)
			int orientationIndex = Orientation.index(p0, p1, p2);

			// Map JTS orientation index to VertexAngle enum
			switch (orientationIndex) {
				case Orientation.COUNTERCLOCKWISE :
					return VertexAngle.LEFT_TURN;
				case Orientation.CLOCKWISE :
					return VertexAngle.RIGHT_TURN;
				case Orientation.COLLINEAR :
					return VertexAngle.COLLINEAR;
				default :
					LOGGER.warn("Unexpected JTS Orientation index: " + orientationIndex);
					return VertexAngle.COLLINEAR;
			}
		} catch (Exception e) {
			LOGGER.error("Error during JTS Orientation calculation: " + e.getMessage());
			e.printStackTrace();
			return VertexAngle.COLLINEAR;
		}
	}

	static InfiniteSpeedType calculateInfiniteSpeedType(WavefrontEdge edge0, WavefrontEdge edge1, VertexAngle angle) {
		if (angle == VertexAngle.COLLINEAR) {
			// Use unit normals and weights as before
			Vector2D normal1 = edge1.getSupportingLine().getUnitNormal();
			Vector2D dir0 = edge0.getSupportingLine().getNormalDirection();
			// Ensure unitNormal calculation is safe (handles zero length dir)
			if (normal1 == null || dir0 == null) { // Check if normals are valid
				return InfiniteSpeedType.NONE;
			}
			double orient = dir0.getX() * normal1.getY() - dir0.getY() * normal1.getX();
			if (orient > SurfConstants.ZERO_DET) { // Tolerance
				return InfiniteSpeedType.OPPOSING;
			} else if (Math.abs(edge0.getWeight() - edge1.getWeight()) > SurfConstants.ZERO_WEIGHT_DIFF) { // Tolerance
				return InfiniteSpeedType.WEIGHTED;
			} else {
				return InfiniteSpeedType.NONE;
			}
		} else {
			return InfiniteSpeedType.NONE;
		}
	}

	static Coordinate calculateVelocity(WavefrontEdge edge0, WavefrontEdge edge1, VertexAngle angle, InfiniteSpeedType infiniteSpeed,
			Coordinate initialPosition) {
		if (infiniteSpeed == InfiniteSpeedType.NONE) {
			if (angle != VertexAngle.COLLINEAR) {
				WavefrontSupportingLine line0 = edge0.getSupportingLine();
				WavefrontSupportingLine line1 = edge1.getSupportingLine();
				Coordinate intersect = line0.lineAtOne().lineIntersection(line1.lineAtOne()); // Intersection at t=1
				if (intersect != null) {
					// Velocity is vector from initial pos to intersection at t=1
					return new Coordinate(intersect.getX() - initialPosition.getX(), intersect.getY() - initialPosition.getY());
				} else {
					// Parallel lines? Should have been caught by COLLINEAR? Error.
					LOGGER.error("No intersection for non-collinear vertex at time 1. Setting zero velocity.");
					return new Coordinate(0, 0);
				}
			} else { // COLLINEAR and NONE infinite speed (implies same direction, same weight)
				// Velocity is along the normal (weighted normal?) - check C++ intent
				Vector2D weightedNormal = edge0.getSupportingLine().getNormal(); // Get weighted normal
				if (weightedNormal == null) { // Handle degenerate case
					LOGGER.error("Cannot get weighted normal for collinear vertex. Setting zero velocity.");
					return new Coordinate(0, 0);
				} else {
					// Velocity is w * unit_normal. Since getNormal() = w * unit_normal, this is
					// just getNormal()
					return new Coordinate(weightedNormal.getX(), weightedNormal.getY());
				}
			}
		} else { // Infinite speed
			return new Coordinate(0, 0); // Finite velocity is zero
		}
	}

	@Override
	public String toString() {
		if (isInfinite) {
			return "WV(Inf)";
		}
		return "WV" + id + "(" + posStart.getX() + "," + posStart.getY() + ")";
	}
	// equals/hashCode based on ID for identity semantics if needed,
	// but using object identity is often sufficient.
	// Coordinate equality is based on value.
}