package com.github.micycle1.surferj.wavefront;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import org.locationtech.jts.geom.Polygon;

import com.github.micycle1.surferj.kinetics.KineticTriangulation;
import com.github.micycle1.surferj.kinetics.WavefrontEdge;

/**
 * Main class orchestrating the straight skeleton computation. It holds the
 * input data, the kinetic triangulation, the resulting skeleton DCEL, and the
 * wavefront propagator that drives the simulation.
 */
public class SkeletonStructure {

	// Make fields accessible within the package or via getters as needed
	final Polygon input;
	final List<WavefrontEdge> wavefrontEdges;
	private final KineticTriangulation kt;
	final SkeletonDCEL skeleton; // The (initially empty) skeleton structure
	final WavefrontPropagator wp;
	
	public SkeletonStructure() {
		// NOTE for tests
		this.input = null;
		this.wavefrontEdges = new ArrayList<>(); // Initialize empty list
		this.skeleton = new SkeletonDCEL(); // Initialize placeholder DCEL
		this.kt = new KineticTriangulation(this.skeleton); // Pass DCEL reference
		this.wp = new WavefrontPropagator(this); // Propagator needs reference back
	}

	/**
	 * Constructs a SkeletonStructure with the given input geometry.
	 *
	 * @param input The input polygon data.
	 */
	public SkeletonStructure(Polygon input) {
		// NOTE: C++ uses std::forward/move. Java copy/assignment is different.
		// Assume Polygon is either immutable or defensively copied if needed.
		this.input = Objects.requireNonNull(input, "Input cannot be null");
		this.wavefrontEdges = new ArrayList<>(); // Initialize empty list
		this.skeleton = new SkeletonDCEL(); // Initialize placeholder DCEL
		this.kt = new KineticTriangulation(this.skeleton); // Pass DCEL reference
		this.wp = new WavefrontPropagator(this); // Propagator needs reference back
	}

	/**
	 * Initializes the kinetic triangulation based on the input data and sets up the
	 * event queue in the wavefront propagator.
	 *
	 * @param restrictComponent If >= 0, restrict processing to only this input
	 *                          component. -1 processes all components.
	 */
	public void initialise(int restrictComponent) {
		// 1. Initialize the Kinetic Triangulation
		kt.initialise(input); // , wavefrontEdges, restrictComponent

		// 2. Setup the Event Queue using the initialized triangulation
		wp.setupQueue(kt); // This also calls kt.setQueue(eq)
	}

	/**
	 * Gets the initial input data.
	 *
	 * @return The Polygon object.
	 */
	public Polygon getInput() {
		return input;
	}

	/**
	 * Gets the kinetic triangulation managing the wavefront state.
	 *
	 * @return The KineticTriangulation object.
	 */
	public KineticTriangulation getKt() {
		return kt;
	}

	/**
	 * Gets the skeleton DCEL structure where the final result is stored.
	 *
	 * @return The SkeletonDCEL object (initially a placeholder).
	 */
	public SkeletonDCEL getSkeleton() {
		return skeleton; // Directly return the member
		// Or return kt.getSkeleton() if KT holds the primary reference
	}

	/**
	 * Gets the wavefront propagator instance that drives the simulation.
	 *
	 * @return The WavefrontPropagator object.
	 */
	public WavefrontPropagator getWp() {
		return wp;
	}
}