// #define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include "IO/ChronoTimer.hpp"
#include "catch2/catch.hpp"
#include "version/version.hpp"
#include <iostream>

int main(int argc, char *argv[]) {

  IO::ChronoTimer timer{"Catch2 tests"};

  std::cout << "ampsci tests (Catch2):\n";
  std::cout << "AMPSCI v: " << version::version() << '\n';
  std::cout << "Libraries:\n" << version::libraries() << '\n';
  std::cout << "Compiled: " << version::compiled() << '\n';

  int result = Catch::Session().run(argc, argv);

  return result;
}