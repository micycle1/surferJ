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

	// Placeholder for the infinite vertex
	public static final WavefrontVertex INFINITE_VERTEX = new WavefrontVertex();
	
	public enum VertexAngle {
        LEFT_TURN, RIGHT_TURN, COLLINEAR
    }

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
		// NOTE for tests
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
            Vector2D dir0 = edge0.getSupportingLine().getDirection();
            Vector2D dir1 = edge1.getSupportingLine().getDirection();
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
                        velocity = new Coordinate(
                            intersect.getX() - initialPosition.getX(),
                            intersect.getY() - initialPosition.getY()
                        );
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
    
    public VertexAngle getAngle() {
    	return angle;
    }

    // Compute position at time t (assuming motion starts at t=0)
    public Coordinate getPositionAt(double time) {
        if (isInfinite || infiniteSpeed != InfiniteSpeedType.NONE) {
            return initialPosition; // Does not move
        }
        return new Coordinate(
            initialPosition.getX() + velocity.getX() * time,
            initialPosition.getY() + velocity.getY() * time
        );
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