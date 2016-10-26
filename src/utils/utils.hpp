/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   utils.hpp
 * Author: acolomitchi
 *
 * Created on 24 October 2016, 10:00 PM
 */

#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstddef>
#include <cmath>

namespace distspctrm {

namespace detail {

template <
  std::size_t N,
  typename Result,
  typename container_type,
  std::size_t Offset,
  template <
     typename R_type, std::size_t OffsetIx, 
     typename context_type, typename...Other
  > class op,
  typename ctx_type,
  typename... opArgs
> struct rec_accumulator {
  
  inline static Result accumulate(
      const container_type& src, const ctx_type& context
  ) {
    return
        op<Result, Offset, ctx_type, opArgs...>::apply(src, context)
      + rec_accumulator<N-1, Result, container_type, Offset+1, op, ctx_type, opArgs...>::accumulate(
          src, context
        )
    ;
  }
};

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
  
template <typename Result, std::size_t CoordIx, class RefPointT>
struct sq_coord_diff_to {
  
  template <typename OtherPointT> inline static Result apply(
    const OtherPointT& src, const RefPointT& reference
  ) {
    Result diff=static_cast<Result>(src[CoordIx]-reference[CoordIx]);
    return diff*diff;
  } 
};


} // namespace detail

template <typename Result, typename Point0, typename Point1, size_t N>
Result unrolled_euclid_dist(const Point0& p0, const Point1& p1) 
{
  Result val=detail::rec_accumulator<
    N, Result, Point1, 0,
    detail::sq_coord_diff_to, Point0
  >::accumulate(p1, p0);
  
  return static_cast<Result>(std::sqrt(val)); 
}

template <typename R, typename Point0, typename Point1, std::size_t N>
R euclid_dist(const Point0 & p0, const Point1 & p1) {
  R val=static_cast<R>(0.0);
  for(size_t i=0; i<N; i++) {
    R diff=p0[i]-p1[i];
    val+=diff*diff;
  }
  return static_cast<R>(std::sqrt(val)); 
}


}
#endif /* UTILS_HPP */

