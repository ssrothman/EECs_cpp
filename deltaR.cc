#include <math.h>

#include "deltaR.h"

float dR2(float eta1, float phi1, float eta2, float phi2) {
  /*
    * Compute delta R^2 between (eta1, phi1) and (eta2, phi2)
    */
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  if (dphi > M_PI) {
    dphi = 2 * M_PI - dphi;
  } else if (dphi < -M_PI) {
    dphi = 2 * M_PI + dphi;  //correct up to sign
  }
  return deta * deta + dphi * dphi;
}

void fillDR(const float* const eta, const float* const phi, const int nPart, std::vector<float>& dRs) {
  /*
   *  Compute the pairwise delta r^2 between all the particles in the jet
   *  
   *  nPart: number of jet constituents
   *  dRs: array to fill with dR^2 values
   */
  dRs.clear();
  int i, j;
  for (i = 0; i < nPart - 1; ++i) {
    for (j = i + 1; j < nPart; ++j) {
      dRs.push_back(sqrtf(dR2(eta[i], phi[i], eta[j], phi[j])));
    }
  }
}