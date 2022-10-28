#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include <vector>

typedef std::vector<std::vector<std::vector<int>>> comp_t;
typedef std::vector<std::vector<int>> factor_t;

int fact(int n);

void fillCompositions(const int n, comp_t& out);

void fillSymFactors(const int N, const comp_t& compositions, factor_t& out);

template <typename T, typename R>
T intPow(const T a, const R b) {
  //mine, very stupid
  //should be upgraded to to square multiply
  //also, assumes b > 0
  T result = 1;
  for (unsigned i = 0; i < (unsigned)b; ++i) {
    result *= a;
  }
  return result;
}

size_t choose(int n, int k);

size_t simp(const size_t nPart, const size_t order);

#endif
