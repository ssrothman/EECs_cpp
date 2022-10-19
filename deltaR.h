#ifndef DELTA_R_H
#define DELTA_R_H

#include <vector>

double dR2(double eta1, double phi1, double eta2, double phi2);

void fillDR(const double* const eta, const double* const phi, const int nPart, std::vector<double>& dRs);

#endif
