#include <iostream>

namespace qip {

void progbar(int i, int max, int length = 50) {
  const int len = (length - 2) / 2; // total bar length is 2*len+2
  const int current = int(len * double(i) / double(max));
  std::cout << "[";
  for (auto j = 0; j < current; ++j) {
    std::cout << "**";
  }
  for (auto j = current; j < len; ++j) {
    std::cout << "  ";
  }
  std::cout << "]   \r" << std::flush;
}

void progbar50(int i, int max) {
  const int len = 49;
  const int prev = int(len * double(i) / double(max) + 1.0e-6);
  const int current = int(len * double(i + 1) / double(max) + 1.0e-6);
  // this ensures progbar line always exactly 50 chars long
  for (int j = prev; j < current; ++j)
    std::cout << "=" << std::flush;
  if (current == len)
    std::cout << "|\n" << std::flush;
}

} // namespace qip
