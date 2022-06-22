#ifndef ITERATING_H
#define ITERATING_H

#include <vector>

void printOrd(const std::vector<int> ord);

void iterate(const int dims, std::vector<int>& ordinates, const int nParts);

size_t getIndex(const std::vector<int>& ordinates, const int nPart);

#endif