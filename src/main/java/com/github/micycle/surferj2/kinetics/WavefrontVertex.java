package com.github.micycle.surferj2.kinetics;

import java.util.concurrent.atomic.AtomicLong;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.math.Vector2D;

// --- Wavefront Vertex ---
public class WavefrontVertex {
	private static final AtomicLong idCounter = new AtomicLong(0);
	public final long id;
	public final Coordinate initialPosition; // Corresponds to pos_zero/pos_start
	public final boolean isInfinite; // NOTE messy? with also having InfiniteSpeedType?

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
		if (id == 0 ) {
			return edge0;
		}
		else if (id==1) {
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

	public WavefrontVertex(double x, double y, double vx, double vy) {
		// NOTE constructor for tests
		this.id = idCounter.incrementAndGet();
		this.initialPosition = new Coordinate(x, y);
		this.isInfinite = false;
		this.velocity = new Coordinate(vx, vy);
		this.angle = VertexAngle.COLLINEAR;
		this.infiniteSpeed = InfiniteSpeedType.NONE;
	}

	public WavefrontVertex(Coordinate initialPosition) {
		this.id = idCounter.incrementAndGet();
		this.initialPosition = initialPosition;
		this.isInfinite = false;
		this.velocity = new Coordinate(0, 0);
		this.angle = VertexAngle.COLLINEAR;
		this.infiniteSpeed = InfiniteSpeedType.NONE;
	}

	public void setIncidentEdges(WavefrontEdge edge0, WavefrontEdge edge1) {
		if (this.isInfinite)
			throw new IllegalStateException("Cannot set edges for infinite vertex");

		// Basic check: the edges should share this vertex's coordinate
		Coordinate sharedCoord = null;
		if (edge0 != null && edge1 != null) {
			if (edge0.getVertex(0) == this)
				sharedCoord = edge0.getVertex(0).initialPosition;
			else if (edge0.getVertex(1) == this)
				sharedCoord = edge0.getVertex(1).initialPosition;

			if (sharedCoord == null || !sharedCoord.equals(this.initialPosition)) {
				throw new IllegalArgumentException(
						"Edges do not consistently share this vertex coordinate " + this.initialPosition + " edge0: " + edge0 + " edge1: " + edge1);
			}
			if ((edge1.getVertex(0) != this && edge1.getVertex(1) != this)) {
				throw new IllegalArgumentException("Edge1 does not contain this vertex coordinate " + this.initialPosition + " edge1: " + edge1);
			}

			// Compute kinetic properties
			// Step 1: Compute angle between edge directions
			Vector2D dir0 = edge0.getSupportingLine().getNormalDirection();
			Vector2D dir1 = edge1.getSupportingLine().getNormal();
			double det = dir0.getX() * dir1.getY() - dir0.getY() * dir1.getX();
			if (det > 1e-12) {
				angle = VertexAngle.LEFT_TURN;
			} else if (det < -1e-12) {
				angle = VertexAngle.RIGHT_TURN;
			} else {
				angle = VertexAngle.COLLINEAR;
			}

			// Step 2: Determine infinite speed type
			if (angle == VertexAngle.COLLINEAR) {
				Vector2D normal1 = edge1.getSupportingLine().getUnitNormal(); // NOTE unit normal
				double orient = dir0.getX() * normal1.getY() - dir0.getY() * normal1.getX();
				if (orient > 1e-12) {
					infiniteSpeed = InfiniteSpeedType.OPPOSING;
				} else if (Math.abs(edge0.getWeight() - edge1.getWeight()) > 1e-12) {
					infiniteSpeed = InfiniteSpeedType.WEIGHTED;
				} else {
					infiniteSpeed = InfiniteSpeedType.NONE;
				}
			} else {
				infiniteSpeed = InfiniteSpeedType.NONE;
			}

			// Step 3: Compute velocity
			if (infiniteSpeed == InfiniteSpeedType.NONE) {
				if (angle != VertexAngle.COLLINEAR) {
					// Intersect supporting lines at time 1
					WavefrontSupportingLine line0 = edge0.getSupportingLine();
					WavefrontSupportingLine line1 = edge1.getSupportingLine();
					Coordinate intersect = line0.lineIntersection(line1);
					if (intersect != null) {
						velocity = new Coordinate(intersect.getX() - initialPosition.getX(), intersect.getY() - initialPosition.getY());
					} else {
						throw new IllegalStateException("No intersection for non-collinear vertex at time 1");
					}
				} else {
					// Collinear with same direction and weight
					Vector2D normal = edge0.getSupportingLine().getUnitNormal();
					velocity = new Coordinate(normal.getX(), normal.getY());
				}
			} else {
				velocity = new Coordinate(0, 0);
			}
		} else {
			angle = VertexAngle.COLLINEAR;
			infiniteSpeed = InfiniteSpeedType.NONE;
			velocity = new Coordinate(0, 0);
		}
		this.edge0 = edge0;
		this.edge1 = edge1;
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
	 * Checks if the angle formed by the incident wavefront edges is convex (left turn)
	 * or straight (collinear). Conventionally, the infinite vertex is considered straight/convex.
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
	 * Checks if the angle formed by the incident wavefront edges is reflex (right turn)
	 * or straight (collinear). Conventionally, the infinite vertex is considered straight/convex.
	 *
	 * @return true if the vertex angle is REFLEX or STRAIGHT, false if CONVEX.
	 */
	public boolean isReflexOrStraight() {
	    // Similar to isConvexOrStraight, the infinite vertex is typically not considered reflex.
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

	public void setIncidentEdge(int index, WavefrontEdge edge) {
		if (isInfinite)
			throw new IllegalStateException("Cannot set edges for infinite vertex");
		if (index == 0)
			this.edge0 = edge;
		else if (index == 1)
			this.edge1 = edge;
		else
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