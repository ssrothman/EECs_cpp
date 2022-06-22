#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include <vector>

typedef std::vector<std::vector<std::vector<int>>> comp_t;
typedef std::vector<std::vector<int>> factor_t;

int fact(int n);

void fillCompositions(const int n, comp_t& out);

void fillSymFactors(const int N, const comp_t& compositions, factor_t& out);

size_t intPow(int a, int b);

float intPow(float a, int b);

size_t choose(int n, int k);

#endif