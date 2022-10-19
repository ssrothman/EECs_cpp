#ifndef EEC_H
#define EEC_H

#include <vector>

void projectedEEC(const double* const pt,
                  const double* const eta,
                  const double* const phi,
                  const int nPart,
                  const unsigned int maxL,
                  std::vector<double>& dRs, 
                  std::vector <double>& wts,
                  const unsigned int N,
                  std::vector<std::vector<std::vector<double>>>* coefs); 

void full3ptEEC(const double* const pt,
                const double* const eta,
                const double* const phi,
                const int nPart,
                std::vector<double>& dR1,
                std::vector<double>& dR2,
                std::vector<double>& dR3,
                std::vector<double>& wts);

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
                std::vector<double>& wts);

void EECnonIRC(const double* const pt,
               const double* const eta,
               const double* const phi,
               const int nPart,
               const int p1, const int p2,
               std::vector<double>& dRs,
               std::vector<double>& wts);

#endif
