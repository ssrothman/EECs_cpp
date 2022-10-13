#ifndef EEC_H
#define EEC_H

#include<array>

template <unsigned int N>
void projectedEEC(const float* const pt,
                  const float* const eta,
                  const float* const phi,
                  const int nPart,
                  const int maxL,
                  std::vector<float>& dRs,
                  std::vector<float>& wts,
                  std::vector<std::vector<std::array<int, N>>>* tuples=nullptr);

void full3ptEEC(const float* const pt,
                const float* const eta,
                const float* const phi,
                const int nPart,
                std::vector<float>& dR1,
                std::vector<float>& dR2,
                std::vector<float>& dR3,
                std::vector<float>& wts);

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
                std::vector<float>& wts);

void EECnonIRC(const float* const pt,
               const float* const eta,
               const float* const phi,
               const int nPart,
               const int p1, const int p2,
               std::vector<float>& dRs,
               std::vector<float>& wts);

#endif
