/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   model.hpp
 * Author: acolomitchi
 *
 * Created on 24 October 2016, 8:02 PM
 */

#ifndef MODEL_HPP
#define MODEL_HPP

#include <cstdlib>
#include <exception>

#include <vector>

#include <type_traits>

#include "../utils/utils.hpp"

namespace distspctrm {

template <typename R, std::size_t N> class ndpoint 
{
public:
  using elem_t=
    typename std::enable_if<std::is_floating_point<R>::value, R>::type;
  
  static constexpr std::size_t DIM=N;
  
  ndpoint() = default;
  
  template <typename... coordt> ndpoint(coordt... x) : elems_ {static_cast<R>(x)...} {
  }
  ndpoint(const ndpoint& other) : elems_() {
    *this=other;
  }
  template <typename PointType> ndpoint(const PointType& other) : elems_() {
    *this = other;
  }
  
  ndpoint& operator=(const ndpoint& other) {
    for(size_t i=0; i<N; i++) {
      this->elems_[i]=other.elems_[i];
    }
    return *this;
  }
  
  template <typename PointT> ndpoint& operator=(const PointT& other) {
    for(size_t i=0; i<N; i++) {
      this->elems_[i]=other[i];
    }
  }
  
  const R& operator[](std::size_t i) const { return this->elems_[i]; }

  R& operator[](std::size_t i) { return this->elems_[i]; }
  
protected:
  R elems_[N];
};


template <typename R, std::size_t ndim>
class ParallelPointStorage {
public:

  class PointIt {
  private:
    const ParallelPointStorage& owner_;
    std::size_t           index_;
    bool                  validatedOk_;

  protected:
    PointIt(const ParallelPointStorage& owner, std::size_t index) 
      : owner_(owner), index_(index), validatedOk_(false)
    { }
  public:
    
    PointIt(const PointIt& other) 
      : owner_(other.owner_), index_(other.index_), validatedOk_(false)
    { }
    
    PointIt& operator=(const PointIt& other) {
      this->owner_=other.owner_;
      this->index_=other.index_;
      this->validatedOk_=other.validatedOk_;
    }

    const R& operator[](std::size_t coordIx) const {
//      if( !this->validatedOk_ && ! this->valid() ) {
//        throw std::out_of_range("");
//      }
      return this->owner_.storage_[coordIx][this->index_];
    }
    
    operator std::size_t() const {
      return this->index_;
    }
    
    operator ndpoint<R, ndim>() const {
      return ndpoint<R,ndim>(*this);
    }
    
    
    bool valid() const {
      bool ret=this->index_<this->owner_.size();
      return ret;
    }
    
    PointIt& operator+=(std::size_t delta) {
      this->index_+=delta;
      // this->validatedOk_=this->valid();
      return *this;
    }
    PointIt& operator-=(std::size_t delta) {
      this->index_=(delta>this->index_ ? 0 : this->index_ - delta);
      // this->validatedOk_=this->valid();
      return *this;
    }
    PointIt& operator++() {
      this->index_++;
      // this->validatedOk_=this->valid();
      return *this;
    }
    PointIt& operator--() {
      if(this->valid() && this->index_>0 ) {
        this->index_--;
      }
      else if( ! this->valid() && this->owner_.size() ){
        this->index_=this->owner_.size()-1;
      }
      // this->validatedOk_=this->valid();
      return *this;
    }
    
    PointIt& slide(std::size_t to) {
      this->index_ = to;
      // this->validatedOk_=this->valid();
      return *this;
    }
    
    bool operator==(const PointIt& other) const {
      bool ret=(&this->owner_ == &other.owner_);
      if(ret) {
        ret = 
            this->valid() 
          ? (other.valid() && this->index_==other.index_)
          : ! other.valid()
        ;
      }
      return ret;
    }
    
    bool operator != (const PointIt& other) const {
      return ! (*this==other);
    }
    
    friend PointIt operator+(const PointIt& left, std::size_t delta) {
      PointIt ret(left);
      ret+=(delta);
      return ret;
    }
    friend PointIt operator-(const PointIt& left, std::size_t delta) {
      PointIt ret(left);
      ret-=(delta);
      return ret;
    }
    
    friend class ParallelPointStorage;
  };
private:
  void ensureCap(std::size_t newCap) {
    if(newCap<8) {
      newCap=8;
    }
    std::size_t cap=this->cap_;
    if(0==this->cap_) {
      for(std::size_t i=0; i<ndim; i++) {
        this->storage_[i]=static_cast<arr_type>(
          std::calloc(newCap, sizeof(elem_type))
        );
        if( ! this->storage_[i]) {
          throw std::bad_alloc();
        }
      }
      this->cap_=newCap;
    }
    else if(this->cap_ < newCap) {
      if(newCap>2048) {
        newCap=((newCap+1023)/1024)*1024;
      }
      else {
        std::size_t cap=this->cap_;
        while(cap<newCap) {
          cap+=(cap>>1); // 3/2 of oldCap
        }
        newCap=cap;
      }
      for(std::size_t i=0; i<ndim; i++) {
        arr_type newArr=static_cast<arr_type>(
          std::realloc(this->storage_[i], newCap*sizeof(elem_type))
        );
        if( ! newArr ) {
          throw std::bad_alloc();
        }
        this->storage_[i]=newArr;
      }
      this->cap_=newCap;
    }
  }
public:
  ParallelPointStorage(std::size_t initialCap=8) : storage_(), cap_(0), count_(0) {
    this->ensureCap(initialCap);
  }
  
  ~ParallelPointStorage() {
    for(std::size_t i=0; i<ndim; i++) {
      if(this->storage_[i]) {
        std::free(this->storage_[i]);
        this->storage_[i]=NULL;
      }
    }
  }
  
  template<class Point>
  void add(const Point& p) {
    if(this->count_>=this->cap_) {
      this->ensureCap(this->count_+1);
    }
    for(std::size_t i=0; i<ndim; i++) {
      this->storage_[i][this->count_]=p[i];
    }
    this->count_++;
  }
  
  PointIt operator[](std::size_t i) {
    return PointIt(*this, i);
  }
  
  const PointIt operator[](std::size_t i) const {
    return PointIt(*this, i);
  }
  
  std::size_t size() const {
    return this->count_;
  }
  
  PointIt begin() const {
    return PointIt(*this, 0);
  }
  
  PointIt end() const {
    return PointIt(*this, static_cast<std::size_t>(-1));
  }
protected:
  using elem_type=typename std::remove_cv<R>::type;
  using arr_type=
    typename std::enable_if<
      std::is_floating_point<R>::value,
      elem_type
    >::type *
  ;
  arr_type storage_[ndim];
  std::size_t cap_, count_;
};

template <typename R, std::size_t numDivs=4>
class dist_histogram {
public:
  dist_histogram() : counts_()
  {
    this->clear();
//    for(std::size_t i=0; i<numDivs; i++) {
//      this->counts_[i]=0;
//    }
  }

  void clear() {
    for(std::size_t i=0; i<numDivs; i++) {
      this->counts_[i]=0;
    }
  }
  
  void consume(R distance) {
    const R (&limits)[numDivs]=thresholds();
    for(size_t i=numDivs-1; i<numDivs; i--) {
      if(distance>limits[i]) {
        this->counts_[i]++;
        break;
      }
    }
  }
  
  
protected:
  std::size_t counts_[numDivs];
private:
  static R limits_[numDivs];
  static const decltype(limits_)& thresholds() {
    static bool inited=false;
    if( ! inited) {
      for(size_t i=0; i<numDivs; i++) {
        limits_[i]=static_cast<R>(i)/numDivs;
      }
      inited = true;
    }
    return limits_;
  }
};

template <typename R, std::size_t numDivs> 
R dist_histogram<R,numDivs>::limits_[numDivs];

} // namespace distspctrm


#endif /* MODEL_HPP */

