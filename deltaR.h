#ifndef DELTA_R_H
#define DELTA_R_H

#include <vector>

float dR2(float eta1, float phi1, float eta2, float phi2);

void fillDR(const float* const eta, const float* const phi, const int nPart, std::vector<float>& dRs);

#endif