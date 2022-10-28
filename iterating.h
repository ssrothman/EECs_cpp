#ifndef ITERATING_H
#define ITERATING_H

#include <vector>
#include <iostream>

template <typename T>
void printOrd(const std::vector<T> ord) {
  std::cout << "(";
  unsigned int i = 0;
  for (i = 0; i < ord.size() - 1; ++i) {
    std::cout << ord[i] << ", ";
  }
  std::cout << ord[ord.size() - 1] << ")";
}

template <typename T, typename R>
void iterate(const T dims, std::vector<R>& ordinates, const T nParts) {
  // iterate over dimensions in reverse...
  for (int dim = dims - 1; dim >= 0; --dim) {
    if (ordinates[dim] < nParts - dims + dim) {
      ++ordinates[dim];
      for (int d2 = dim + 1; d2 < dims; ++d2) {
        ordinates[d2] = ordinates[d2 - 1] + 1;
      }
      return;
    }
    if (dim > 0) {
      ordinates[dim] = ordinates[dim - 1] + 1;
    } else {
      ordinates[dim] = 0;
    }
  }
}

template <typename T, typename R>
void iterate_wdiag(const T dims, std::vector<R>& ordinates, const T nParts) {
  // iterate over dimensions in reverse...
  for (int dim = dims - 1; dim >= 0; --dim) {
    if (ordinates[dim] < nParts-1) {
      ++ordinates[dim];
      for (int d2 = dim + 1; d2 < dims; ++d2) {
        ordinates[d2] = ordinates[d2 - 1];
      }
      return;
    }
    if (dim > 0) {
      ordinates[dim] = ordinates[dim - 1];
    } else {
      ordinates[dim] = 0;
    }
  }
}



template <typename T>
void iterate_all(const T dims, std::vector<T>& ordinates, const T nParts);

size_t getIndex(const std::vector<int>& ordinates, const int nPart);

template <typename T>
void iterate_all(const T dims, std::vector<T>& ordinates, const T nParts) {
  // iterate over dimensions in reverse...
  for (int dim = dims - 1; dim >= 0; --dim) {
    if (ordinates[dim] < nParts-1) {
      ++ordinates[dim];
      break;
    } else {
      ordinates[dim] = 0;
    }
  }
}

#endif
