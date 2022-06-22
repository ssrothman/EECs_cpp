#include "combinatorics.h"

#include <stdio.h>
#include <math.h>
#include <algorithm>

int fact(int n) {
  int result = 1;
  for (int i = 2; i < n + 1; ++i) {
    result *= i;
  }
  return result;
}

void fillCompositions(const int n, comp_t& out) {
  out.clear();
  for (int i = 0; i < n; ++i) {
    std::vector<std::vector<int>> next;
    out.push_back(next);
  }

  int a[n + 1];
  int k, y, x, l;
  for (int q = 0; q < n + 1; ++q) {
    a[q] = 0;
  }
  k = 1;
  y = n - 1;
  while (k != 0) {
    x = a[k - 1] + 1;
    k -= 1;
    while (2 * x <= y) {
      a[k] = x;
      y -= x;
      k += 1;
    }
    l = k + 1;
    while (x <= y) {
      a[k] = x;
      a[l] = y;
      //yield a[:k + 2];
      do {
        std::vector<int> next;
        next.reserve(k + 1);
        for (int q = 0; q < k + 2; ++q) {
          next.push_back(a[q]);
        }
        out[k + 1].push_back(next);
      } while (std::next_permutation(a, a + k + 2));
      x += 1;
      y -= 1;
    }
    a[k] = x + y;
    y = x + y - 1;
    //yield a[:k + 1];
    do {
      std::vector<int> next;
      next.reserve(k);
      for (int q = 0; q < k + 1; ++q) {
        next.push_back(a[q]);
      }
      out[k].push_back(next);
    } while (std::next_permutation(a, a + k + 1));
  }
}

void fillSymFactors(const int N, const comp_t& compositions, factor_t& out) {
  out.clear();
  out.reserve(N);
  int factN = fact(N);
  int nextFactor;
  for (int i = 0; i < N; ++i) {  //for each possible length of composition
    std::vector<int> next;
    next.reserve(compositions[i].size());
    for (size_t j = 0; j < compositions[i].size(); ++j) {  //for each composition of that length
      nextFactor = factN;
      for (size_t k = 0; k < compositions[i][j].size(); ++k) {  //for each int in that composition
        nextFactor /= fact(compositions[i][j][k]);
      }  //end loop over individual composition
      next.push_back(nextFactor);
    }  //end loop over compositions of a given length
    out.push_back(next);
  }  //end loop over composition length
}

size_t intPow(int a, int b) {
  //mine, very stupid
  //should be upgraded to to square multiply
  size_t result = 1;
  for (int i = 0; i < b; ++i) {
    result *= a;
  }
  return result;
}

float intPow(float a, int b) {
  //mine, very stupid
  //should be upgraded to to square multiply
  float result = 1;
  for (int i = 0; i < b; ++i) {
    result *= a;
  }
  return result;
}

size_t choose(int n, int k) {
  if (k > n)
    return 0;
  if (k * 2 > n)
    k = n - k;
  if (k == 0)
    return 1;

  size_t result = n;
  for (int i = 2; i <= k; ++i) {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}