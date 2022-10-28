#include "iterating.h"
#include "combinatorics.h"
#include "deltaR.h"
#include "eec.h"
#include "vecND.h"

#include <iostream>
#include <algorithm>
#include <array>

#define TESTPROJ

#ifdef TESTPROJ 
#define TEST
#endif

#ifdef TESTFULL
#define TEST
#endif


double getWt(const double* const pt,
            const int nPart,
            const std::vector<int>& ord,
            const int M,
            const comp_t& compositions,
            const factor_t& symFactors,
            const unsigned int N,
            std::vector<std::vector<std::vector<double>>>* coefs=nullptr,
            vecND<double>* tuplewts=nullptr,
            vecND<int>* tupleiDR=nullptr,
            const int iDR=0) {
  double result = 0;
  for (size_t i = 0; i < compositions[M - 1].size(); ++i) {  //for each composition
    double nextWt = symFactors[M - 1][i];
    std::vector<int> fullTuple;
    for (size_t j = 0; j < compositions[M - 1][i].size(); ++j) {  //for each element
      nextWt *= intPow<double>(pt[ord[j]], compositions[M - 1][i][j]);
      for(int copy=0; copy<compositions[M-1][i][j]; ++copy){
        fullTuple.push_back(ord[j]);
      }
    }  //end for each element
    //
    //for each element, but only if we need to fill the coefs array
    for(unsigned int iord=0; coefs && iord<ord.size(); ++iord){ 
      //we overcount M times
      //coefs->at(compositions[M-1][i][iord]-1).at(ord[iord]).at(iDR) 
      //  += nextWt/(M*intPow(pt[ord[iord]], compositions[M-1][i][iord]));

      coefs->at(0).at(ord[iord]).at(iDR) 
        += nextWt/(M);
    } //end for each element, coefs edition

    if(tuplewts){
      tupleiDR->at(fullTuple) = iDR;
      tuplewts->at(fullTuple) = nextWt;
    }

    result += nextWt;
  }  //end for each composition


  return result;
}

double getWt_nonIRC(const double* const pt,
                   const int nPart,
                   const int i, const int j,
                   const int p1, const int p2){
  return intPow<double>(pt[i], p1) * intPow<double>(pt[j], p2);
}

void projectedMway(const double* const pt,
                   const double* const eta,
                   const double* const phi,
                   const int nPart,

                   const int M,

                   const std::vector<double>& dRs,
                   std::vector<double>& wts,

                   const comp_t compositions,
                   const factor_t symFactors,

                   const std::vector<int>* const cache,
                   const int L,
                   std::vector<int>* newCache,
                   const unsigned int N,
                   
                   std::vector<std::vector<std::vector<double>>>* coefs=nullptr,
                   vecND<double>* tuplewts=nullptr,
                   vecND<int>* tupleiDR=nullptr
                   ) {
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
    double bestDR = -1;

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
        iterate<int, int>(L, ordL, M);
      }  //end iterate over L-way combinations of M elements
    }
    double newWt = getWt(pt, nPart, ord, M, compositions, symFactors, N, coefs, tuplewts, tupleiDR, bestIdx);
    wts[bestIdx] += newWt;
    if (newCache)
      (*newCache)[getIndex(ord, nPart)] = bestIdx;
    iterate<int, int>(M, ord, nPart);
  }  //end iterate over M-way combinations of nPart elements
}

void resolved3way(const double* const pt,
                  const double* const eta,
                  const double* const phi,
                  const int nPart,

                  std::vector<double>& dR1,
                  std::vector<double>& dR2,
                  std::vector<double>& dR3,
                  std::vector<double>& wts,

                  const comp_t compositions,
                  const factor_t symFactors,

                  const std::vector<int> cache,

                  const unsigned int N) {
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

  double R1, R2, R3;

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

    std::vector<double> Rs = {R1, R2, R3};
    std::sort(Rs.begin(), Rs.end(), std::greater<double>());

    dR1.push_back(Rs[0]);
    dR2.push_back(Rs[1]);
    dR3.push_back(Rs[2]);
    wts.push_back(getWt(pt, nPart, ord, 3, compositions, symFactors, N, nullptr, 0));
    iterate<int,int>(3, ord, nPart);
  }
}

void resolved4way(const double* const pt,
                  const double* const eta,
                  const double* const phi,
                  const int nPart,

                  std::vector<double>& dR1,
                  std::vector<double>& dR2,
                  std::vector<double>& dR3,
                  std::vector<double>& dR4,
                  std::vector<double>& dR5,
                  std::vector<double>& dR6,
                  std::vector<double>& wts,

                  const comp_t compositions,
                  const factor_t symFactors,

                  const std::vector<int> cache,
                  
                  const unsigned int N) {
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

  double R1, R2, R3, R4, R5, R6;

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

    std::vector<double> Rs = {R1, R2, R3, R4, R5, R6};
    std::sort(Rs.begin(), Rs.end(), std::greater<double>());

    dR1.push_back(Rs[0]);
    dR2.push_back(Rs[1]);
    dR3.push_back(Rs[2]);
    dR4.push_back(Rs[3]);
    dR5.push_back(Rs[4]);
    dR6.push_back(Rs[5]);

    wts.push_back(getWt(pt, nPart, ord, 4, compositions, symFactors, N, nullptr, 0));
    iterate<int,int>(4, ord, nPart);
  }
}

void projectedEEC(const double* const pt,
                  const double* const eta,
                  const double* const phi,
                  const int nPart,
                  const unsigned int maxL,
                  std::vector<double>& dRs, 
                  std::vector <double>& wts,
                  const unsigned int N,
                  std::vector<std::vector<std::vector<double>>>* coefs,
                  vecND<double>* tuplewts=nullptr,
                  vecND<int>* tupleiDR=nullptr
                  ) { 
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
   *
   * coefs: if not nullptr, to be filled with coefficients with 
   *  which each power of each particle contributes. 
   *  ie $EEC(\Delta R) = \sum_i^{N_{part}} \sum_\alpha \text{coefs}\left[\alpha, i, \Delta R\right] p_T(i)^\alpha$
   */

  //fill dR array
  dRs.clear();
  fillDR(eta, phi, nPart, dRs);
 
  //initialize zero weights
  wts.clear();
  wts.resize(dRs.size());

  //initialize empty coefs array
  if(coefs){
    coefs->clear();
    coefs->resize(N-1);
    for(unsigned i=0; i<coefs->size(); ++i){
      coefs->at(i).resize(nPart);
      for(signed j=0; j<nPart; ++j){
        coefs->at(i).at(j).resize(dRs.size());
        for(unsigned k=0; k<dRs.size(); ++k){
          coefs->at(i).at(j).at(k) = 0;
        }
      }
    }
  }

  //precompute compositions and symmetry factors
  comp_t compositions;
  fillCompositions(N, compositions);

  factor_t symFactors;
  fillSymFactors(N, compositions, symFactors);

  std::vector<int> cache(0, 0);
  std::vector<int> newCache(0, 0);

  for (unsigned int M = 2; M <= N; ++M) {
    if (M == 2) {
      projectedMway(pt, eta, phi, 
                    nPart, 
                    M, 
                    dRs, wts, 
                    compositions, symFactors, 
                    nullptr, 0, &newCache,
                    N,
                    coefs, tuplewts, tupleiDR);
    } else if (M <= maxL) {
      std::swap(cache, newCache);
      projectedMway(pt, eta, phi, 
                    nPart, 
                    M, 
                    dRs, wts, 
                    compositions, symFactors, 
                    &cache, M - 1, &newCache,
                    N,
                    coefs, tuplewts, tupleiDR);
    } else {
      projectedMway(pt, eta, phi, 
                    nPart, 
                    M, 
                    dRs, wts, 
                    compositions, symFactors, 
                    &newCache, maxL, nullptr,
                    N,
                    coefs, tuplewts, tupleiDR);
    }
  }
  
  //do the zero-dR terms
  //ie the configuration where you pick the same particle every time
  //what's the symmetry factor? 1
  dRs.push_back(0);
  std::vector<int> zeroOrd(N, 0);
  double zerowt = 0;
  for(int iPart=0; iPart<nPart; ++iPart){
    for(unsigned i=0; i<N; ++i){
      zeroOrd[i] = iPart;
    }
    double nextwt = intPow(pt[iPart], N);
    zerowt += nextwt;
    if(tuplewts){
      tuplewts->at(zeroOrd) = nextwt;
      tupleiDR->at(zeroOrd) = dRs.size()-1;
    }
  }
  wts.push_back(zerowt);
}



void EECnonIRC(const double* const pt,
               const double* const eta,
               const double* const phi,
               const int nPart,
               const int p1, const int p2,
               std::vector<double>& dRs,
               std::vector<double>& wts){
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
    iterate<int, int>(2, ord, nPart);
  } //end for each pair
}

void full3ptEEC(const double* const pt,
                const double* const eta,
                const double* const phi,
                const int nPart,
                std::vector<double>& dR1,
                std::vector<double>& dR2,
                std::vector<double>& dR3,
                std::vector<double>& wts) {
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
  projectedMway(pt, eta, phi, nPart, 2, dR1, wts, compositions, symFactors, nullptr, 0, &cache, 3);
  int n2 = choose(nPart, 2);
  dR2.insert(dR2.end(), n2, -1.);
  dR3.insert(dR3.end(), n2, -1.);

  //append 3-point component
  resolved3way(pt, eta, phi, nPart, dR1, dR2, dR3, wts, compositions, symFactors, cache, 3);
}

void full4ptEEC(const double* const pt,
                const double* const eta,
                const double* const phi,
                const int nPart,
                std::vector<double>& dR1,
                std::vector<double>& dR2,
                std::vector<double>& dR3,
                std::vector<double>& dR4,
                std::vector<double>& dR5,
                std::vector<double>& dR6,
                std::vector<double>& wts) {
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
  projectedMway(pt, eta, phi, nPart, 2, dR1, wts, compositions, symFactors, nullptr, 0, &cache, 4);
  int n2 = choose(nPart, 2);
  dR2.insert(dR2.end(), n2, -1.);
  dR3.insert(dR3.end(), n2, -1.);
  dR4.insert(dR4.end(), n2, -1.);
  dR5.insert(dR5.end(), n2, -1.);
  dR6.insert(dR6.end(), n2, -1.);

  resolved3way(pt, eta, phi, nPart, dR1, dR2, dR3, wts, compositions, symFactors, cache, 4);
  int n3 = choose(nPart, 3);
  dR4.insert(dR4.end(), n3, -1.);
  dR5.insert(dR5.end(), n3, -1.);
  dR6.insert(dR6.end(), n3, -1.);

  resolved4way(pt, eta, phi, nPart, dR1, dR2, dR3, dR4, dR5, dR6, wts, compositions, symFactors, cache, 4);
}

#ifdef TESTPROJ

#define ORDER 4

int main() {
  double pT[] = {1., 2., 0.5, 2., 3.};
  double eta[] = {0., 0.1, 0.4, 1.0, 0.4};
  double phi[] = {0., 0.2, 0.4, 0.0, -0.5};

  double totalWT =0 ;

  std::vector<double> dRs, wts;
  std::vector<std::vector<std::vector<double>>> coefs;
  std::vector<std::vector<tuple_t>> tuples;
  vecND<double> tuplewts(5, ORDER, 0);
  vecND<int> tupleiDR(5, ORDER, -1);

  projectedEEC(pT, eta, phi, 5, 10, dRs, wts, ORDER, &coefs, &tuplewts, &tupleiDR);

  std::cout << "ord\tdR\twt" << std::endl;
  size_t maxIter = dRs.size(), iter = 0;
  std::vector<int> ord = {0, 1};
  for (iter = 0; iter < maxIter; ++iter) {
    printOrd(ord);
    std::cout << ": " << dRs[iter] << "\t " << wts[iter] << std::endl;
    std::cout << "\t= ";
    for(unsigned i=0; i<5; ++i){
      if(coefs[0][i][iter]>0){
        printf("%0.2f ", coefs[0][i][iter]);
        printf("+ ");
      } 
    }
    std::cout << std::endl;
    iterate<int, int>(2, ord, 5);
    totalWT += wts[iter];
  }

  std::cout << "Tuple contributions are:" << std::endl;

  maxIter = intPow(5, ORDER);
  std::vector<int> ord2(ORDER, 0);
  for(unsigned iter=0; iter<maxIter; ++iter){
    if(tuplewts.at(ord2) > 0){
      printOrd(ord2);
      printf(" -> dR %d, wt %0.3f\n", tupleiDR.at(ord2), tuplewts.at(ord2));
    }
    iterate_all(ORDER, ord2, 5);
  }
  std::cout << "TOTAL WEIGHT: " << totalWT << std::endl;

  std::vector<int> ord3(ORDER, 0);
  size_t maxIter3 = choose(ORDER + 5 - 1, ORDER);
  for(size_t iter=0; iter<maxIter3; ++iter){
    printOrd(ord3);
    printf("\n");
    iterate_wdiag(ORDER, ord3, 5);
  }
}
#endif

#ifdef TESTFULL
int main() {
  double pT[] = {1., 2., 0.5, 2., 3.};
  double eta[] = {0., 0.1, 0.4, 1.0, 0.4};
  double phi[] = {0., 0.2, 0.4, 0.0, -0.5};

  std::vector<double> dR1, dR2, dR3, dR4, dR5, dR6, wts;

  double totalWT=0;

  full4ptEEC(pT, eta, phi, 5, dR1, dR2, dR3, dR4, dR5, dR6, wts);
  
  printf("Combo\t\t\tdR1\tdR2\tdR3\tdR4\tdR5\tdR6\twt\n");
  size_t maxIter = choose(5,2), iter = 0;
  std::vector<int> ord = {0, 1};
  for (iter = 0; iter < maxIter; ++iter) {
    printOrd(ord);
    printf("\t\t\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n",dR1[iter], dR2[iter], dR3[iter], dR4[iter], dR5[iter], dR6[iter], wts[iter]);
    iterate<int, int>(2, ord, 5);
    totalWT += wts[iter];
  }

  maxIter = choose(5,3) + choose(5,2);
  std::vector<int> ord2 = {0, 1, 2};
  for (iter = choose(5,2); iter < maxIter; ++iter) {
    printOrd(ord2);
    printf("\t\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n",dR1[iter], dR2[iter], dR3[iter], dR4[iter], dR5[iter], dR6[iter], wts[iter]);
    iterate<int ,int>(3, ord2, 5);
    totalWT += wts[iter];
  }

  maxIter = choose(5,4) + choose(5,3) + choose(5,2);
  std::vector<int> ord3 = {0, 1, 2, 3};
  for (iter = choose(5, 3) + choose(5,2); iter < maxIter; ++iter) {
    printOrd(ord3);
    printf("\t\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n",dR1[iter], dR2[iter], dR3[iter], dR4[iter], dR5[iter], dR6[iter], wts[iter]);
    iterate<int,int>(4, ord3, 5);
    totalWT += wts[iter];
  }



  std::cout << "TOTAL WEIGHT: " << totalWT << std::endl;
}
#endif
