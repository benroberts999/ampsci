#pragma once
#include "Physics/AtomData.hpp"
#include "Physics/NuclearData.hpp"
#include <string>

namespace AtomData {

void printData(const Nuclear::Isotope &nuc);

int parse_A(const std::string &A_str, int z = 0);

//==============================================================================
void printConstants();

//==============================================================================
void periodicTable(std::string z_str, std::string a_str);

} // namespace AtomData
