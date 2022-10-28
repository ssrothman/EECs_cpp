#include "iterating.h"

size_t getIndex(const std::vector<int>& ordinates, const int nPart) {
  size_t result = 0;
  size_t partPow = 1;
  for (size_t dim = 0; dim < ordinates.size(); ++dim) {
    result += ordinates[dim] * partPow;
    partPow *= nPart;
  }
  return result;
}
