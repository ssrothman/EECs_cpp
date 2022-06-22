#include "iterating.h"

#include <iostream>

void printOrd(const std::vector<int> ord) {
  std::cout << "(";
  unsigned int i = 0;
  for (i = 0; i < ord.size() - 1; ++i) {
    std::cout << ord[i] << ", ";
  }
  std::cout << ord[ord.size() - 1] << ")";
}

void iterate(const int dims, std::vector<int>& ordinates, const int nParts) {
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

size_t getIndex(const std::vector<int>& ordinates, const int nPart) {
  size_t result = 0;
  size_t partPow = 1;
  for (size_t dim = 0; dim < ordinates.size(); ++dim) {
    result += ordinates[dim] * partPow;
    partPow *= nPart;
  }
  return result;
}