#pragma once
#include <vector>

namespace MAT {

// XXX Pass as reference, const?? + overload w. template ?? XXX
int invertMatrix(const std::vector<std::vector<double>> &inmat,
                 std::vector<std::vector<double>> &outmat);

int invertMatrix(std::vector<std::vector<double>> &inmat);

int invertMatrix(const std::vector<std::vector<float>> &inmat,
                 std::vector<std::vector<float>> &outmat);

int invertMatrix(std::vector<std::vector<float>> &inmat);

double calcDeterminant(const std::vector<std::vector<double>> &inmat);

double calcDeterminant(const std::vector<std::vector<float>> &inmat);

int linsolve(const std::vector<std::vector<double>> &inmat,
             const std::vector<double> &invec, std::vector<double> &outvec);

} // namespace MAT
