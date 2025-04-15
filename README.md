# surferJ

This is an ongoing Java port of [surfer2](https://github.com/cgalab/surfer2), a highly robust, triangulation-based straight skeleton algorithm – the most powerful implementation for straight skeletons out there.

However this is not an easy undertaking... it's conceptually difficult, the algorithm details are complex, the data structures are nasty and the original code is academic C++. On top of that, I haven't actually been able to get the original code built and running (the dependencies are nuisance). Hopefully I'll get there...

## Status

#### Done:

The **static** structure of the kinetic triangulation at t=0, along with when kinetic triangles collapse (solving quadratic equations), the movement of wavefront vertices and edges, and the various "events" between vertices and edges (such as when vertex crashes into an opposing wavefront edge). These parts are also well-tested, deriving assumptions about intended behavior from the reference implementation (since no tests are given there).

#### TODO:

The building blocks are in place so the big next step is the wavefront **propagation** itself. This entails stepping through the simulation, deciding which event to handle next and updating topology as triangles collapse. Events must be processed in their correct order.

I also might need to revist the triangulation initalisation (whcih currently comprises constrained triangles only): 
> First, create a constrained triangulation of the convex hull of the input graph. Then, for every edge e on the convex hull add one unbounded triangle, that is, one triangle with e as one edge and two infinite edges going outwards. All these edges are thought to meet at infinity
