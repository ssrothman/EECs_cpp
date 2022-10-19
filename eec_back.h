#ifndef EEC_H
#define EEC_H

//#define VERBOSE

#include <vector>
#include <tuple>
#include <iostream>
#include <array>
#include <math.h>
#include <stdint.h>
#include <algorithm>

typedef uint16_t idx_t;
typedef double coord_t;
typedef std::tuple<idx_t, idx_t> pair;
typedef std::vector<std::vector<std::vector<idx_t>>> comp_t;
typedef std::vector<std::vector<idx_t>> factor_t;

void eec_onejet(double* jet, int nPart, int nFeat, int N, double* dRs, int nDR, double* wts, int nWT, int maxL);

void eec(double* jets,
         int nPartTot,
         int nFeat,
         int* jetIdxs,
         int nJets,
         int N,
         double* dRs,
         int nDRTot,
         double* wts,
         int nWTTot,
         int* dRIdxs,
         int nDRIdxs,
         int maxL);

size_t choose(idx_t n, idx_t k);

#endif
