#include "iterating.h"
#include "combinatorics.h"
#include "deltaR.h"
#include "eec.h"

#include <iostream>
#include <algorithm>

#ifdef TESTPROJ 
#define TEST
#endif

#ifdef TESTFULL
#define TEST
#endif

float getWt(const float* const pt,
            const int nPart,
            const int N,
            const std::vector<int>& ord,
            const int M,
            const comp_t& compositions,
            const factor_t& symFactors) {
  float result = 0;
  for (size_t i = 0; i < compositions[M - 1].size(); ++i) {  //for each composition
    float nextWt = symFactors[M - 1][i];
    for (size_t j = 0; j < compositions[M - 1][i].size(); ++j) {  //for each element
      nextWt *= intPow(pt[ord[j]], compositions[M - 1][i][j]);
    }  //end for each element
    result += nextWt;
  }  //end for each composition

  return result;
}

float getWt_nonIRC(const float* const pt,
                   const int nPart,
                   const int i, const int j,
                   const int p1, const int p2){
  return intPow(pt[i], p1) * intPow(pt[j], p2);
}

void projectedMway(const float* const pt,
                   const float* const eta,
                   const float* const phi,
                   const int nPart,

                   const int N,
                   const int M,

                   const std::vector<float>& dRs,
                   std::vector<float>& wts,

                   const comp_t compositions,
                   const factor_t symFactors,

                   const std::vector<int>* const cache,
                   const int L,
                   std::vector<int>* newCache) {
  /*
     * compute the weight contribution from 
     * M-way combinations of particles
     * 
     * pt: particle transverse momenta, normalized to the jet pt
     * eta: particle eta
     * phi: particle eta
     * nPart: number of particles
     * 
     * N: correlator order
     * M: number of particles in combination
     * 
     * dRs: precomputed vector of delta R values
     * wts: vector of weights. 1-1 correspondence with dRs
     *     will be modified 
     * 
     * compositions: vector of compositions
     * syFactors: vector of symmetry factors 
     *      associated with each composistion 
     * 
     * cache: cached L-way max delta Rs
     * L: number of particles in cache combinations
     * newCache: to be filled with new cache values if not nullptr
     */

  if (newCache)
    newCache->resize(intPow(nPart, M));

  size_t maxIter = choose(nPart, M);
  std::vector<int> ord(M);  //which M-fold combination?
  for (int i = 0; i < M; ++i) {
    ord[i] = i;
  }
  //loop over all M-way combinations of nPart elements
  for (size_t iter = 0; iter < maxIter; ++iter) {
    int bestIdx = 0;
    float bestDR = -1;

    if (M == 2) {
      bestIdx = iter;
    } else {
      std::vector<int> ordL(L);  //subcombinations from cache
      for (int i = 0; i < L; ++i) {
        ordL[i] = i;
      }
      size_t maxIterL = choose(M, L);

      //iterate over L-way combinations of M elements
      //to find the subcombination with the largest cached dR
      for (size_t iterL = 0; iterL < maxIterL; ++iterL) {
        std::vector<int> sub(L);  //actual subcombination
        for (int i = 0; i < L; ++i) {
          sub[i] = ord[ordL[i]];
        }
        int idx = (*cache)[getIndex(sub, nPart)];
        if (dRs[idx] > bestDR) {
          bestIdx = idx;
          bestDR = dRs[idx];
        }
        iterate(L, ordL, M);
      }  //end iterate over L-way combinations of M elements
    }
    float newWt = getWt(pt, nPart, N, ord, M, compositions, symFactors);
    wts[bestIdx] += newWt;
    if (newCache)
      (*newCache)[getIndex(ord, nPart)] = bestIdx;
    iterate(M, ord, nPart);
  }  //end iterate over M-way combinations of nPart elements
}

void resolved3way(const float* const pt,
                  const float* const eta,
                  const float* const phi,
                  const int nPart,

                  const int N,

                  std::vector<float>& dR1,
                  std::vector<float>& dR2,
                  std::vector<float>& dR3,
                  std::vector<float>& wts,

                  const comp_t compositions,
                  const factor_t symFactors,

                  const std::vector<int> cache) {
  /*
     * compute the fully-resolved weights and dRs from 
     * 3-way combinations of particles
     * 
     * pt: particle transverse momenta, normalized to the jet pt
     * eta: particle eta
     * phi: particle eta
     * nPart: number of particles
     * 
     * N: correlator order. Obviously must be >= 3
     *      * 
     * dR1: to be appended with longest side dRs
     * dR2: to be appended with second longest side dRs
     * dR3: to be appended with third longest side dRs
     * wts: to be appended with corresponding weights
     * 
     * compositions: vector of compositions
     * syFactors: vector of symmetry factors 
     *      associated with each composistion 
     * 
     * cache: cached 2-way dR locations
     */

  float R1, R2, R3;

  if(N<3){
    std::cerr << "Can't call resolved3way() with N<3" << std::endl;
    return;
  }

  size_t maxIter = choose(nPart, 3);
  std::vector<int> ord = {0, 1, 2};
  for (size_t iter = 0; iter < maxIter; ++iter) {
    std::vector<int> sub1 = {ord[0], ord[1]};
    std::vector<int> sub2 = {ord[0], ord[2]};
    std::vector<int> sub3 = {ord[1], ord[2]};

    R1 = dR1[cache[getIndex(sub1, nPart)]];
    R2 = dR1[cache[getIndex(sub2, nPart)]];
    R3 = dR1[cache[getIndex(sub3, nPart)]];

    std::vector<float> Rs = {R1, R2, R3};
    std::sort(Rs.begin(), Rs.end(), std::greater<float>());

    dR1.push_back(Rs[0]);
    dR2.push_back(Rs[1]);
    dR3.push_back(Rs[2]);
    wts.push_back(getWt(pt, nPart, N, ord, 3, compositions, symFactors));
    iterate(3, ord, nPart);
  }
}

void resolved4way(const float* const pt,
                  const float* const eta,
                  const float* const phi,
                  const int nPart,

                  const int N,

                  std::vector<float>& dR1,
                  std::vector<float>& dR2,
                  std::vector<float>& dR3,
                  std::vector<float>& dR4,
                  std::vector<float>& dR5,
                  std::vector<float>& dR6,
                  std::vector<float>& wts,

                  const comp_t compositions,
                  const factor_t symFactors,

                  const std::vector<int> cache) {
  /*
     * compute the fully-resolved weights and dRs from 
     * 4-way combinations of particles
     * 
     * pt: particle transverse momenta, normalized to the jet pt
     * eta: particle eta
     * phi: particle eta
     * nPart: number of particles
     * 
     * N: correlator order. Obviously must be >= 4
     *      * 
     * dR1: to be appended with longest side dRs
     * dR2: to be appended with second longest side dRs
     * dR3: to be appended with third longest side dRs
     * dR4: to be appended with fourth longest side dRs
     * dR5: to be appended with fifth longest side dRs
     * dR6: to be appended with sixth longest side dRs
     * wts: to be appended with corresponding weights
     * 
     * compositions: vector of compositions
     * syFactors: vector of symmetry factors 
     *      associated with each composistion 
     * 
     * cache: cached 2-way dR locations
     */

  float R1, R2, R3, R4, R5, R6;

  if(N<4){
    std::cerr << "Can't call resolved4way() with N<4" << std::endl;
    return;
  }

  size_t maxIter = choose(nPart, 4);
  std::vector<int> ord = {0, 1, 2, 3};
  for (size_t iter = 0; iter < maxIter; ++iter) {
    std::vector<int> sub1 = {ord[0], ord[1]};
    std::vector<int> sub2 = {ord[0], ord[2]};
    std::vector<int> sub3 = {ord[0], ord[3]};
    std::vector<int> sub4 = {ord[1], ord[2]};
    std::vector<int> sub5 = {ord[1], ord[3]};
    std::vector<int> sub6 = {ord[2], ord[3]};

    R1 = dR1[cache[getIndex(sub1, nPart)]];
    R2 = dR1[cache[getIndex(sub2, nPart)]];
    R3 = dR1[cache[getIndex(sub3, nPart)]];
    R4 = dR1[cache[getIndex(sub4, nPart)]];
    R5 = dR1[cache[getIndex(sub5, nPart)]];
    R6 = dR1[cache[getIndex(sub6, nPart)]];

    std::vector<float> Rs = {R1, R2, R3, R4, R5, R6};
    std::sort(Rs.begin(), Rs.end(), std::greater<float>());

    dR1.push_back(Rs[0]);
    dR2.push_back(Rs[1]);
    dR3.push_back(Rs[2]);
    dR4.push_back(Rs[3]);
    dR5.push_back(Rs[4]);
    dR6.push_back(Rs[5]);

    wts.push_back(getWt(pt, nPart, N, ord, 4, compositions, symFactors));
    iterate(4, ord, nPart);
  }
}

void projectedEEC(const float* const pt,
                  const float* const eta,
                  const float* const phi,
                  const int nPart,
                  const int N,
                  const int maxL,
                  std::vector<float>& dRs, 
                  std::vector <float>& wts) { 
  /*
   * EEC, projected onto shortest side
   * 
   * pt: particle pt, normalized to jet pt
   * eta: particle eta
   * phi: particle  phi,
   * 
   * N: correaltor order
   * 
   * maxL: max cache order to allow
   * 
   * dRs: to be filled with dR values
   * wts: to be filled with weight values
   */

  //fill dR array
  dRs.clear();
  fillDR(eta, phi, nPart, dRs);
 
  //initialize zero weights
  wts.clear();
  wts.resize(dRs.size());

  //precompute compositions and symmetry factors
  comp_t compositions;
  fillCompositions(N, compositions);

  factor_t symFactors;
  fillSymFactors(N, compositions, symFactors);

  std::vector<int> cache(0, 0);
  std::vector<int> newCache(0, 0);

  for (int M = 2; M <= N; ++M) {
    if (M == 2) {
      projectedMway(pt, eta, phi, 
                    nPart, 
                    N, M, 
                    dRs, wts, 
                    compositions, symFactors, 
                    nullptr, 0, &newCache);
    } else if (M <= maxL) {
      std::swap(cache, newCache);
      projectedMway(pt, eta, phi, 
                    nPart, 
                    N, M, 
                    dRs, wts, 
                    compositions, symFactors, 
                    &cache, M - 1, &newCache);
    } else {
      projectedMway(pt, eta, phi, 
                    nPart, 
                    N, M, 
                    dRs, wts, 
                    compositions, symFactors, 
                    &newCache, maxL, nullptr);
    }
  }
}

void EECnonIRC(const float* const pt,
               const float* const eta,
               const float* const phi,
               const int nPart,
               const int p1, const int p2,
               std::vector<float>& dRs,
               std::vector<float>& wts){
  /*
   * EEC, with non IRC-safe energy weighting
   * Currently only support 2-way EEC
   *
   * pt: particle pt, normalized to jet pt,
   * eta: particle eta
   * phi: particle phi
   * nPart: number of particles
   *
   * p1, p2: powers in correlator weighting (p1 = p2 = 1 is IRC-safe case)
   * 
   * dRs: to be filled with dR values
   * wts: to be filled with weight values
   */
  //fill dR array
  dRs.clear();
  fillDR(eta, phi, nPart, dRs);
 
  //initialize zero weights
  wts.clear();
  wts.resize(dRs.size());

  //which pair?
  std::vector<int> ord(2);
  size_t maxIter = choose(nPart, 2);

  //initialize 0th pair
  ord[0] = 0;
  ord[1] = 1;

  for(size_t iter=0; iter<maxIter; ++iter){//for each pair
    wts[iter] += getWt_nonIRC(pt, nPart, ord[0], ord[1], p1, p2);
    if(p1 != p2){ //if powers not symmetric, add symmetric term
      wts[iter] += getWt_nonIRC(pt, nPart, ord[0], ord[1], p2, p1);
    }
    iterate(2, ord, nPart);
  } //end for each pair
}

void full3ptEEC(const float* const pt,
                const float* const eta,
                const float* const phi,
                const int nPart,
                std::vector<float>& dR1,
                std::vector<float>& dR2,
                std::vector<float>& dR3,
                std::vector<float>& wts) {
  /*
   * fully resolved 3-point correaltor
   * 
   * pt: particle pt, normalized to jet pt
   * eta: particle eta
   * phi: particle phi,
   * 
   * dR1: to be filled with longest side dR
   * dR2: to be filled with second longest side dR
   * dR3: to be filled with third longest side dR
   * wts: to be filled with weight values
   */

#ifdef TEST
  std::cout << "top of full 3-point EEC" << std::endl << std::endl;
#endif

  dR1.clear();
  dR2.clear();
  dR3.clear();

  //prime dR array with 2-point values
  fillDR(eta, phi, nPart, dR1);

  //initialize zero weights
  wts.clear();
  wts.resize(dR1.size());

  //precompute compositions and symmetry factors
  comp_t compositions;
  fillCompositions(3, compositions);

  factor_t symFactors;
  fillSymFactors(3, compositions, symFactors);

#ifdef TEST
  std::cout << "composition\tsymmetry factor" << std::endl;
  for (unsigned i = 0; i < symFactors.size(); ++i) {
    for (unsigned j = 0; j < symFactors[i].size(); ++j) {
      printOrd(compositions[i][j]);
      std::cout << " " << symFactors[i][j] << std::endl;
    }
  }
  std::cout << std::endl;
#endif

  std::vector<int> cache(0, 0);

  //start by computing 2-point component
  //already fully-resolved, so use projected function
  projectedMway(pt, eta, phi, nPart, 3, 2, dR1, wts, compositions, symFactors, nullptr, 0, &cache);
  int n2 = choose(nPart, 2);
  dR2.insert(dR2.end(), n2, -1.);
  dR3.insert(dR3.end(), n2, -1.);

  //append 3-point component
  resolved3way(pt, eta, phi, nPart, 3, dR1, dR2, dR3, wts, compositions, symFactors, cache);
}

void full4ptEEC(const float* const pt,
                const float* const eta,
                const float* const phi,
                const int nPart,
                std::vector<float>& dR1,
                std::vector<float>& dR2,
                std::vector<float>& dR3,
                std::vector<float>& dR4,
                std::vector<float>& dR5,
                std::vector<float>& dR6,
                std::vector<float>& wts) {
  /*
   * fully resolved 4-point correaltor
   * 
   * pt: particle pt, normalized to jet pt
   * eta: particle eta
   * phi: particle phi,
   * 
   * dR1: to be filled with longest side dR
   * dR2: to be filled with second longest side dR
   * dR3: to be filled with third longest side dR
   * dR4: to be filled with fourth longest side dR
   * dR5: to be filled with fifth longest side dR
   * dR6: to be filled with sixth longest side dR
   * wts: to be filled with weight values
   */


  dR1.clear();
  dR2.clear();
  dR3.clear();
  dR4.clear();
  dR5.clear();
  dR6.clear();

  //print dR array with 2-point values
  fillDR(eta, phi, nPart, dR1);

  //initialize zero weights
  wts.clear();
  wts.resize(dR1.size());

  //precompute compositions and symmetry factors
  comp_t compositions;
  fillCompositions(4, compositions);

  factor_t symFactors;
  fillSymFactors(4, compositions, symFactors);

  std::vector<int> cache(0, 0);

  //start by computing 2-point component
  //already fully-resolved, so use projected function
  projectedMway(pt, eta, phi, nPart, 4, 2, dR1, wts, compositions, symFactors, nullptr, 0, &cache);
  int n2 = choose(nPart, 2);
  dR2.insert(dR2.end(), n2, -1.);
  dR3.insert(dR3.end(), n2, -1.);
  dR4.insert(dR4.end(), n2, -1.);
  dR5.insert(dR5.end(), n2, -1.);
  dR6.insert(dR6.end(), n2, -1.);

  resolved3way(pt, eta, phi, nPart, 4, dR1, dR2, dR3, wts, compositions, symFactors, cache);
  int n3 = choose(nPart, 3);
  dR4.insert(dR4.end(), n3, -1.);
  dR5.insert(dR5.end(), n3, -1.);
  dR6.insert(dR6.end(), n3, -1.);

  resolved4way(pt, eta, phi, nPart, 4, dR1, dR2, dR3, dR4, dR5, dR6, wts, compositions, symFactors, cache);
}

#ifdef TESTPROJ
int main() {
  float pT[] = {1., 2., 0.5, 2., 3.};
  float eta[] = {0., 0.1, 0.4, 1.0, 0.4};
  float phi[] = {0., 0.2, 0.4, 0.0, -0.5};

  float totalWT =0 ;

  std::vector<float> dRs, wts;
  projectedEEC(pT, eta, phi, 5, 4, 10, dRs, wts);

  std::cout << "dR\twt" << std::endl;
  size_t maxIter = dRs.size(), iter = 0;
  std::vector<int> ord = {0, 1};
  for (iter = 0; iter < maxIter; ++iter) {
    printOrd(ord);
    std::cout << ": " << dRs[iter] << "\t" << wts[iter] << std::endl;
    iterate(2, ord, 5);
    totalWT += wts[iter];
  }
  std::cout << "TOTAL WEIGHT: " << totalWT << std::endl;
}
#endif

#ifdef TESTFULL
int main() {
  float pT[] = {1., 2., 0.5, 2., 3.};
  float eta[] = {0., 0.1, 0.4, 1.0, 0.4};
  float phi[] = {0., 0.2, 0.4, 0.0, -0.5};

  std::vector<float> dR1, dR2, dR3, dR4, dR5, dR6, wts;

  float totalWT=0;

  full4ptEEC(pT, eta, phi, 5, dR1, dR2, dR3, dR4, dR5, dR6, wts);
  
  printf("Combo\t\t\tdR1\tdR2\tdR3\tdR4\tdR5\tdR6\twt\n");
  size_t maxIter = choose(5,2), iter = 0;
  std::vector<int> ord = {0, 1};
  for (iter = 0; iter < maxIter; ++iter) {
    printOrd(ord);
    printf("\t\t\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n",dR1[iter], dR2[iter], dR3[iter], dR4[iter], dR5[iter], dR6[iter], wts[iter]);
    iterate(2, ord, 5);
    totalWT += wts[iter];
  }

  maxIter = choose(5,3) + choose(5,2);
  std::vector<int> ord2 = {0, 1, 2};
  for (iter = choose(5,2); iter < maxIter; ++iter) {
    printOrd(ord2);
    printf("\t\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n",dR1[iter], dR2[iter], dR3[iter], dR4[iter], dR5[iter], dR6[iter], wts[iter]);
    iterate(3, ord2, 5);
    totalWT += wts[iter];
  }

  maxIter = choose(5,4) + choose(5,3) + choose(5,2);
  std::vector<int> ord3 = {0, 1, 2, 3};
  for (iter = choose(5, 3) + choose(5,2); iter < maxIter; ++iter) {
    printOrd(ord3);
    printf("\t\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n",dR1[iter], dR2[iter], dR3[iter], dR4[iter], dR5[iter], dR6[iter], wts[iter]);
    iterate(4, ord3, 5);
    totalWT += wts[iter];
  }


  std::cout << "TOTAL WEIGHT: " << totalWT << std::endl;
}
#endif
