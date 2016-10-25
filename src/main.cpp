/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: acolomitchi
 *
 * Created on 25 October 2016, 5:46 PM
 */

#include <cassert>

#include <chrono>
#include <iostream>

#include <random>

#include <vector>
#include <signal.h>

#include "utils/utils.hpp"
#include "model/model.hpp"

namespace distspctrm {


template <typename R, std::size_t nDims>
void initUniformRndPoints(std::vector<ndpoint<R,nDims>>& dest, std::size_t numPoints=1000) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<R> dis(0, 1);
    
  dest.clear();
  ndpoint<R,nDims> val;
  for(size_t i=0; i<numPoints; i++) {
    for(size_t j=0; j<nDims; j++) {
      val[j]=dis(gen);
    }
    dest.push_back(val);
  }
}

template <typename R, std::size_t nDims, std::size_t numSlots>
void consume_iterated_inc(
  const std::vector<ndpoint<R,nDims>>& src, 
  dist_histogram<R, numSlots>& dest
) {
  using point_type=ndpoint<R,nDims>;
  std::size_t len=src.size();
  if(len > 0) {
    for(std::size_t i=0; i<len-1; i++) {
      // const ndpoint<R,nDims>& p0=src[i];
      ndpoint<R,nDims> p0(src[i]);
      for(std::size_t j=i+1; j<len; j++) {
        const ndpoint<R,nDims>& p1=src[j];
        dest.consume( 
          euclid_dist<R, point_type, point_type,nDims>(p0, p1) 
        );
      }
    }
  }
}

template <typename R, std::size_t nDims, std::size_t numSlots>
void consume_unrolled_inc(
  const std::vector<ndpoint<R,nDims>>& src, 
  dist_histogram<R, numSlots>& dest
) {
  using point_type=ndpoint<R,nDims>;
  std::size_t len=src.size();
  if( len>0 ) {
    for(std::size_t i=0; i<len-1; i++) {
      // const ndpoint<R,nDims>& p0=src[i];
      ndpoint<R,nDims> p0(src[i]);
      for(std::size_t j=i+1; j<len; j++) {
        const ndpoint<R,nDims>& p1=src[j];
        dest.consume( 
          unrolled_euclid_dist<R,point_type, point_type, nDims>(p0, p1) 
        );
      }
    }
  }
}

template <typename R, std::size_t nDims, std::size_t numSlots>
void consume_iterated_inc(
  const ParallelPointStorage<R,nDims>& src, 
  dist_histogram<R, numSlots>& dest
) {
  using point_type0=ndpoint<R,nDims>; 
  using point_type1=typename ParallelPointStorage<R,nDims>::PointIt;
  std::size_t len=src.size();
  if( len>0 ) {
    for(point_type1 p0=src.begin(); p0<len-1; ++p0) {
      point_type0 p0Copy(p0); 
      point_type1 p1=p0;
      for(++p1; p1<len; ++p1) {
        dest.consume( 
          euclid_dist<R, point_type0, point_type1, nDims>(p0Copy, p1) 
        );
      }
    }
  }
}

template <typename R, std::size_t nDims, std::size_t numSlots>
void consume_unrolled_inc(
  const ParallelPointStorage<R,nDims>& src, 
  dist_histogram<R, numSlots>& dest
) {
  using point_type0=ndpoint<R,nDims>; 
  using point_type1=typename ParallelPointStorage<R,nDims>::PointIt;
  std::size_t len=src.size();
  if( len>0 ) {
    for(point_type1 p0=src.begin(); p0<len-1; ++p0) {
      point_type0 p0Copy(p0);
      point_type1 p1=p0;
      for(++p1; p1<len; ++p1) {
        dest.consume( 
          unrolled_euclid_dist<R, point_type0, point_type1, nDims>(p0Copy, p1) 
        );
      }
    }
  }
}

} // namespace distspctrm

int main(int argc, char *argv[]) {
  
  using namespace distspctrm;
  
  namespace ch=std::chrono;
  
  using ptype=ndpoint<double,16>;
  std::vector<ptype> points;
  initUniformRndPoints(points, 10000);
  
  dist_histogram<double,4> consumer_iterated;
  dist_histogram<double,4> consumer_unrolled;
  
  ParallelPointStorage<ptype::elem_t,ptype::DIM> parallelPoints;
  for(auto p : points) {
    parallelPoints.add(p);
  }
  std::cout << "Init done" << std::endl;
  for(size_t i=0; i<32; i++) {
    consumer_iterated.clear();
    consumer_unrolled.clear();
 
    ch::high_resolution_clock::time_point t1, t2;
    t1 = ch::high_resolution_clock::now();
    consume_iterated_inc(points, consumer_iterated);
    t2 = ch::high_resolution_clock::now();

    auto duration = ch::duration_cast<ch::microseconds>( t2 - t1 ).count();

    std::cout << "+Iterated vec (us): " << duration;

    t1 = ch::high_resolution_clock::now();
    consume_unrolled_inc(points, consumer_unrolled);
    t2 = ch::high_resolution_clock::now();
    duration = ch::duration_cast<ch::microseconds>( t2 - t1 ).count();

    std::cout << " Unrolled vec (us): " << duration << std::endl;

    t1 = ch::high_resolution_clock::now();
    consume_iterated_inc(parallelPoints, consumer_iterated);
    t2 = ch::high_resolution_clock::now();

    duration = ch::duration_cast<ch::microseconds>( t2 - t1 ).count();

    std::cout << "+Iterated par (us): " << duration;

    t1 = ch::high_resolution_clock::now();
    consume_unrolled_inc(parallelPoints, consumer_unrolled);
    t2 = ch::high_resolution_clock::now();
    duration = ch::duration_cast<ch::microseconds>( t2 - t1 ).count();

    std::cout << " Unrolled par (us): " << duration << std::endl;

  }
}