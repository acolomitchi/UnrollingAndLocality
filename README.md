# UnrollingAndLocality
Experiments in C++ on template loop unwinding and data locality effects on performance 

## In brief

Some experiments on the effects on performance of:
1. data representation influence on CPU cache(s) data locality - 
   I felt the need of seeing some numbers after watching 
   [this part](https://youtu.be/fHNmRkzxHWs?t=2081) of 
   Chandler Carruth's 
   "Efficiency with Algorithms, Performance with Data Structures" presentation
   at CPPcon 2014; and
2. explicit loop unwinding by using of C++ templates.

**The playground**: have different representations of a `N-dimensional point` and
a collection thereof and compute all the (distinct) distances between these
points. Study (roughly) the execution time for the above in the context of a
largish set of points - so far I used `10000` points, which amounts to about
5 millions distinct distances to be calculated.

The 'loop unwinding' technique is applied to the individual distance point
calculations, while the 'CPU caches data locality' is explored by storing
the N-dim point collection either as:

1. a vector of structs - all coordinates of a point are locally packed; vs 
2. a set of arrays of coordinates, each array containing the coordinates
   pertaining to one axis/dimension for all points.

(see [row-major vs column major](https://en.wikipedia.org/wiki/Row-major_order))

Don't expect deep-digging into the guts of the CPU's caches, not even detailed
dissection (by instrumentation) of how much to CPU time different parts of the
application contributes - the analysis only looks at the overall execution
time.<br>In fact, as it stands right now (late Oct 2016), not much 
analysis is available, only some preliminary observations and 
*a stable piece of code*. In regards with the analysis, I do intend to 
look a bit more into it if the conditions allow in the future.

<h2 id='look-unwinding'>The compile-time loop unwinding - the technique</h2>

Suppose you have a method that computes the squared difference between the 
coordinates of two N-dim points along one of the axis/dimension (pseudo-code):

```
sq_diff(index, point a, point b) {
   diff=a[index]-b[index];
   return diff*diff;
}
```
The iterative way of computing the euclidian distance is the trivial:
```
const ndim=...; // number of dimensions

euclid_dist(point a, point b) {
    sq_dist=0.0;
    for(i in [0, ndim) ) {
      sq_dist+=sq_diff(i, a, b);
    }
    return sqrt(sq_dist);
}
```

A (runtime) way of recursively computing the same, would be on the line of:

```
const ndim=...; // number of dimensions

function accumulator(index, stepsToGo, point a, point b) {
  if(stepsToGo == 0) { // safety recursion stop
    return 0.0;
  }
  else if(stepsToGo == 1) { 
    // formally, this isn't necessary, but why make an empty step?
    return sq_diff(index, a, b);
  }
  return sq_diff(index, a, b)+accumulator(index+1, stepsToGo-1, a, b);
}

unwound_euclid_dist(point a, point b) {
  return sqrt(accumulator(0, ndim, a, b);
}
```

While the recursive solution will be clearly sub-optimal if applied runtime
(there will be an overhead of calling the `accumulator` function which adds to 
the call of `sq_diff`), if one can determine the compiler to inline both the
`sq_diff` and the `accumulator` functions then the net result is the 
'compile time loop unwinding':

```
unwound_euclid_dist(point a, point b) {
  return sqrt((a[0]-b[0])^2 + (a[1]-b[1])^2 + ... + (a[ndim-1]-b[ndim-1])^2);
}
```

Using templates to implement the recursion and hoping that the C++ compiler 
will be coerced to inline by optimization flags is exactly the
`compile time loop unwinding` is about. 

Before jumping into the code

A note: examining the code above, one can really see that the only
requirement against a `point` type is to provide the 
'access by index operator', so any structure which implements
`operator[]` will make a sufficiently-good point representation
(including an array of an arithmetic type to hold the coordinates)

Second, the `accumulator` function can be made generic enough to work
not only the `sq_diff` but other functions as well - one will need just to pass
the `sq_diff` (or other function) as a parameter.

```
 accumulator(index, stepsToGo, operation, point a, point b) {
   // ... recursion stop guard
   return operation(index, a, b)+accumulator(index+1, stepsToGo-1, a, b);
 }
 unwound_euclid_dist(point a, point b) {
   return sqrt(accumulator(0, ndim, sq_diff, a, b));
 } 
 // computes step of Manhattan distance along the `index` dimension
 abs_diff(index, point a, point b) {
   return abs(a[i]-b[i]);
 }
 unwound_manhattan_dist(point a, point b) {
   return accumulator(0, ndim, abs_diff, a, b);
 }
 ```
Going one step further (and that's the last one), the `operation` parameter
can be thought as a function operating on  coordinate of `point a` 
and having `point b` as sort of a 'context', so that the most general concept
of the `accumulator` could be expressed by a signature like

```
accumulator(begin_index, len, input_array, operation, operation_context)
```
with

* `begin_index - offset into the `array_argument` to start applying the
   `operation` and accumulate the result
* `len` - the number of steps to go, starting from the `begin_index`
* `input_array` - any structure supporting the `[index]` operator and
   acting as input
* `operation` - the transformation to apply to each element of the
  `input_array` before accumulating the results of this transformation;
* `operation_context` - any additional information required for the `operation`
  to transform the input. In the most general case, the type of this
  parameter can be anything, it just so happens that in the case of `sq_diff`, 
  this context is  made from the second `point` (the one to calculate the
  distance to).

So, now, the code... (btw, it's located in the `utils\utils.hpp` header).
There's the need to express the `accumulator` and `sq_diff` functions in
a way that will.. well... compel the compiler to:

1. inline those functions when any optimisations are considered
2. make the runtime dependency of those functions limited to _only_ what will
   be runtime specific: namely the two points.

Therefore:
* express them in the form of templates, taking the dimensions and indexed/offset
  as template parameters
* rely on the 'template definition recursion' rather than runtime recursion
  (that's _compile time_ loop unwinding, right?)
* as C++(11) does not allow partial specialization of function templates,
  we'll need to deal with class templates and implement the runtime functionality
  as methods of these templates 
 
```C++

template <
  std::size_t N, // number of steps to go in the recursive definotion
  typename Result, // the result type of the accumulation
  typename container_type, // the input type. Expected to define the [] operator
  std::size_t Offset, // offset from the start 
  template <
     // result type of the transformation
     typename R_type,
     // which element of the input to use in transformation 
     std::size_t OffsetIx, 
     // the type of additional info
     typename context_type, 
     // other types on which the transformation depends 
     typename...Other
  > class op, // the `operation` type. Expecte to implement an `apply` static method
  typename ctx_type, // the type of additional info that to be passed when 
                     // specializing the `op` 
  typename... opArgs // other types to pass at `op` specialization
> struct rec_accumulator {
  
  inline static Result accumulate(
      const container_type& src, const ctx_type& context
  ) {
    return
        op<Result, Offset, ctx_type, opArgs...>::apply(src, context) // one step...
      + rec_accumulator< // ... added to the rest
          N-1, Result, 
          container_type, Offset+1, 
          op, ctx_type, opArgs...
        >::accumulate(
          src, context
        )
    ;
  }
};

// template specialization for the case the number of steps is 1
// (help the compiler recurse one step less)
template <
  typename Result,
  typename container_type,
  std::size_t Offset,
  template <
     typename R_type, std::size_t OffsetIx, 
     typename context_type, typename...Other
  > class op,
  typename ctx_type,
  typename... opArgs
> 
struct rec_accumulator<1, Result, container_type, Offset, op, ctx_type, opArgs...> {
  
  inline static Result accumulate(
      const container_type& c, 
      const ctx_type& context
  ) {
    // one step without adding "the rest", as "the rest" will be zero
    return op<Result, Offset, ctx_type, opArgs...>::apply(c, context);
  }
};

// just for safety, in case one calls into rec_accumulator with a N of zero 
template <
  typename Result,
  typename container_type,
  std::size_t Offset,
  template <
     typename R_type, std::size_t OffsetIx, 
     typename context_type, typename...Other
  > class op,
  typename ctx_type,
  typename... opArgs
> 
struct rec_accumulator<0, Result, container_type, Offset, op, ctx_type, opArgs...> {
  inline static Result accumulate(
      const container_type& /* c */, 
      const ctx_type& /* context */
  ) {
    return static_cast<Result>(0.0);
  }
};
```

<h2 id='data-locality'>CPU caches data locality</h2>
Well... not exactly. Yes, choosing different storage layouts will influence
how the data is brought into CPU caches, but... no, this is not actually an
proper experiment to actually _determine_ the CPU cache hits/misses rates and/or
out-of-order execution.

[To be continued]

