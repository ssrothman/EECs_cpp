#ifndef VECND_H
#define VECND_H

#include <vector>
#include "combinatorics.h"
#include "iterating.h"
#include <algorithm>

template <class T>
class vecND{
  public:

    vecND(unsigned N, unsigned dim) : vec_(intPow(N, dim)),
                                      dim_(dim), N_(N) {}
    vecND(unsigned N, unsigned dim, T fillval) : vec_(intPow(N, dim), fillval),
                                                 dim_(dim), N_(N) {}

    template <typename I>
    T& at(const std::vector<I>& ord){
      if(ord.size() != dim_){
        throw std::invalid_argument("Trying to index vecND with the wrong dimensionality");
      }
      return vec_.at(idx_(ord));
    }

  private:
    std::vector<T> vec_;
    const unsigned dim_, N_;

    //NB assumes ord is the right size
    //sorts a copy of ord
    template <typename I>
    size_t idx_(const std::vector<I>& ord){
      std::vector<I> sortOrd(ord);
      std::sort(sortOrd.begin(), sortOrd.end());

      size_t result=0;
      for(size_t i=0; i<dim_; ++i){
        result += sortOrd[i] * intPow<size_t>(N_, dim_-i-1);
      }
      return result;
    }
};


#endif

