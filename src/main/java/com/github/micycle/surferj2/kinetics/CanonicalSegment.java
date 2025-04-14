package com.github.micycle.surferj2.kinetics;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineSegment;
import java.util.Objects;

/**
 * Represents a LineSegment with a canonical representation (sorted endpoints)
 * suitable for use as a Map key or in a Set.
 * 
 * USED FOR HASHING. Useless, can't i use unique edges from the triangulation?
 */
public class CanonicalSegment {
	private final Coordinate p0;
	private final Coordinate p1;

	// original coords (so that LineSegment can be created from this having original
	// orientation).
	private final Coordinate p0o;
	private final Coordinate p1o;

	public CanonicalSegment(Coordinate c1, Coordinate c2) {
		this.p0o = c1;
		this.p1o = c2;

		if (c1.compareTo(c2) <= 0) {
			this.p0 = c1;
			this.p1 = c2;
		} else {
			this.p0 = c2;
			this.p1 = c1;
		}
	}

	public CanonicalSegment(LineSegment segment) {
		this(segment.p0, segment.p1);
	}

	public Coordinate getP0() {
		return p0;
	}

	public Coordinate getP1() {
		return p1;
	}

	public LineSegment getSegment() {
		return new LineSegment(p0o, p1o);
	}

	@Override
	public boolean equals(Object o) {
		if (this == o)
			return true;
		if (o == null || getClass() != o.getClass())
			return false;
		CanonicalSegment that = (CanonicalSegment) o;
		// JTS Coordinate equals handles coordinate equality
		return Objects.equals(p0, that.p0) && Objects.equals(p1, that.p1);
	}

	@Override
	public int hashCode() {
		// JTS Coordinate hashCode is value-based
		return Objects.hash(p0, p1);
	}

	@Override
	public String toString() {
		return "CanonicalSegment{" + p0 + " -> " + p1 + '}';
	}
}
