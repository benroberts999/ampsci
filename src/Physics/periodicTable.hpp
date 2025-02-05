#pragma once
#include "Physics/AtomData.hpp"
#include "Physics/NuclearData.hpp"
#include <string>

namespace AtomData {

void printData(const Nuclear::Isotope &nuc);

int parse_A(const std::string &A_str, int z = 0);

//==============================================================================
//! Writes some physical constants to screen
void printConstants();

//! Performs basic unit conversions (writes to screen)
void conversions(double number, const std::string &unit);

//==============================================================================
//! Prints a basic periodic table, and write isotope info ro screen
void periodicTable(std::string z_str, std::string a_str);

} // namespace AtomData
