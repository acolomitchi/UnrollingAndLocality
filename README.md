# UnrollingAndLocality
Experiments in C++ on template loop unrolling and data locality effects on performance 

==In brief==

Some experiments on the effects on performance of:
1. data representation influence on CPU cache(s) data locality - 
   I felt the need of seeing some numbers after watching 
   [this part](https://youtu.be/fHNmRkzxHWs?t=2081) of 
   Chandler Carruth's 
   "Efficiency with Algorithms, Performance with Data Structures" presentation
   at CPPcon 2014; and
2. explicit loop unrolling by using of C++ templates.

**The playground**: have different representations of a `N-dimensional point` and
a collection thereof and compute all the (distinct) distances between these
points. Study (roughly) the execution time for the above in the context of a
largish set of points - so far I used `10000` points, which amounts to about
5 millions distinct distances to be calculated.

The 'loop unrolling' technique is applied to the individual distance point
calculations, while the 'CPU caches data locality' is explored by storing
the N-dim point collection either as:

1. a vector of structs - all coordinates of a point are locally packed); vs 
2. a set of array of coordinates, each array containing the coordinate
   pertaining to one axis/dimension for all points.

(similar to [row-major vs column major](https://en.wikipedia.org/wiki/Row-major_order))

Don't expect deep-digging into the guts of the CPU's caches, not even detailed
dissection (by instrumentation) of how much to CPU time different parts of the
application contributes - the analysis only looks at the overall execution
time.<br>In fact, as it stands right now (late Oct 2016), not much 
analysis is available, only some preliminary observations and 
*a stable piece of code*. In regards with the analysis, I do intend to 
look a bit more into it if the conditions allow in the future.

[to be continued]
