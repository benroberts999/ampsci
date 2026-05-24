#pragma once
#include "IO/InputBlockLegacy.hpp"
#include "Wavefunction/Wavefunction.hpp"

//! Calculates wavefunction and runs optional modules
Wavefunction ampsci(const IO::InputBlockLegacy &input);