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

} // namespace qip
