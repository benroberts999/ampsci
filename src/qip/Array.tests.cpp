#include "qip/Array.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <cassert>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>

TEST_CASE("qip::ArrayView", "[qip][Array][ArrayView][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "qip::ArrayView\n";

  // Tests arrayview, which tests stride iterator

  const std::vector<int> v1{1, 2,  3,  4,  5,  6,  7,  8,
                            9, 10, 11, 12, 13, 14, 15, 16};

  for (std::size_t stride = 1; stride <= v1.size() + 1; ++stride) {

    // nb: when stride = v1.size() + 1 => size=0
    // Should still work, though cannot access elements (undefined behaviour)

    const std::size_t size = v1.size() / stride;
    auto view_1 = qip::ArrayView(v1.data(), v1.size() / stride, stride);
    const auto view_2 = qip::ArrayView(v1.data(), v1.size() / stride, stride);

    // We modify this vector!
    std::vector<int> vb = v1;
    auto view_b = qip::ArrayView(vb.data(), vb.size() / stride, stride);

    REQUIRE(view_1.size() == size);
    REQUIRE(view_2.size() == size);
    REQUIRE(view_1.data() == v1.data());

    // Test element access

    for (std::size_t i = 0; i < view_1.size(); ++i) {
      REQUIRE(view_1[i] == v1[i * stride]);
      REQUIRE(view_1.at(i) == v1[i * stride]);
      REQUIRE(view_1(i) == v1[i * stride]);
      REQUIRE(view_2[i] == v1[i * stride]);
      REQUIRE(view_2.at(i) == v1[i * stride]);
      REQUIRE(view_2(i) == v1[i * stride]);
    }
    if (size != 0) {
      REQUIRE(view_1.front() == v1[0]);
      REQUIRE(view_2.front() == v1[0]);
      REQUIRE(view_b.front() == v1[0]);
      REQUIRE(view_1.back() == v1[(size - 1) * stride]);
      REQUIRE(view_2.back() == v1[(size - 1) * stride]);
      REQUIRE(view_b.back() == v1[(size - 1) * stride]);
    }

    // Test element modification
    for (std::size_t i = 0; i < view_1.size(); ++i) {
      REQUIRE(view_b[i] == vb[i * stride]);
      view_b[i] *= 2;
      view_b.at(i) *= 2;
      view_b(i) *= 2;
      REQUIRE(view_b[i] == 2 * 2 * 2 * view_1[i]);
      // Directly check that underlying vector was actually modified
      REQUIRE(view_b[i] == vb[i * stride]);
      REQUIRE(vb[i * stride] == 2 * 2 * 2 * v1[i * stride]);
    }
    if (size != 0) {
      REQUIRE(view_b.front() == 2 * 2 * 2 * view_1.front());
      REQUIRE(view_b.back() == 2 * 2 * 2 * view_1.back());
      view_b.front() += 5;
      REQUIRE(view_b.front() == 2 * 2 * 2 * view_1.front() + 5);
      // if size=1, front _is_ back, don't increment again
      if (size != 1)
        view_b.back() += 5;
      REQUIRE(view_b.back() == 2 * 2 * 2 * view_1.back() + 5);
    }

    // Check iterator access
    {
      std::size_t count = 0;
      for (auto it = view_1.begin(); it != view_1.end(); ++it) {
        REQUIRE(*it == v1[count]);
        count += stride;
      }
    }

    // ... for const functions
    {
      std::size_t count = 0;
      for (auto it = view_1.cbegin(); it != view_1.cend(); ++it) {
        REQUIRE(*it == v1[count]);
        count += stride;
      }
    }
    // ... for const view
    {
      std::size_t count = 0;
      for (auto it = view_2.cbegin(); it != view_2.cend(); ++it) {
        REQUIRE(*it == v1[count]);
        count += stride;
      }
    }

    // Same, use post-increment
    {
      std::size_t count = 0;
      for (auto it = view_1.begin(); it != view_1.end(); it++) {
        REQUIRE(*it == v1[count]);
        count += stride;
      }
    }

    // Same, +=
    {
      std::size_t count = 0;
      for (auto it = view_1.begin(); it != view_1.end(); it += 1) {
        REQUIRE(*it == v1[count]);
        count += stride;
      }
    }

    // Check iterator access - backwards manually
    {
      std::size_t count = (size - 1) * stride;
      for (auto it = view_1.end() - 1; it != view_1.begin() - 1; --it) {
        REQUIRE(*it == v1[count]);
        count -= stride;
      }
    }

    // Check iterator access - backwards reverse iterator
    {
      std::size_t count = (size - 1) * stride;
      for (auto it = view_1.rbegin(); it != view_1.rend(); ++it) {
        REQUIRE(*it == v1[count]);
        count -= stride;
      }
    }

    std::size_t count = 0;
    std::for_each(view_1.begin(), view_1.end(), [&count, v1, stride](auto a) {
      REQUIRE(a == v1[count]);
      count += stride;
    });
  }

  // test a few algorithms, using default stride of 1
  auto view_1 = qip::ArrayView(v1.data(), v1.size());
  const auto view_2 = qip::ArrayView(v1.data(), v1.size());

  const auto expected_sum = std::accumulate(v1.begin(), v1.end(), 0);

  REQUIRE(expected_sum == std::accumulate(view_1.begin(), view_1.end(), 0));
  REQUIRE(expected_sum == std::accumulate(view_1.cbegin(), view_1.cend(), 0));
  REQUIRE(expected_sum == std::accumulate(view_1.rbegin(), view_1.rend(), 0));
  REQUIRE(expected_sum == std::accumulate(view_1.crbegin(), view_1.crend(), 0));
  REQUIRE(expected_sum == std::accumulate(view_2.cbegin(), view_2.cend(), 0));
  REQUIRE(expected_sum == std::accumulate(view_2.crbegin(), view_2.crend(), 0));

  auto view_3 = view_1;
  const auto view_4 = view_1;
  auto view_5 = view_2;
  const auto view_6 = view_2;
  REQUIRE(expected_sum == std::accumulate(view_3.begin(), view_3.end(), 0));
  REQUIRE(expected_sum == std::accumulate(view_4.cbegin(), view_4.cend(), 0));
  REQUIRE(expected_sum == std::accumulate(view_5.cbegin(), view_5.cend(), 0));
  // Can call non-const member function, but can't actually modify underlying
  REQUIRE(expected_sum == std::accumulate(view_5.begin(), view_5.end(), 0));
  // This should not compile: (and it doesn't)
  // view_5(0) *= 2;
  REQUIRE(expected_sum == std::accumulate(view_6.cbegin(), view_6.cend(), 0));
}

//==============================================================================

TEST_CASE("qip::Array", "[qip][Array][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "qip::Array\n";

  qip::Array<int> IntArray(4, 2, 3);

  for (const auto &el : IntArray) {
    REQUIRE(el == 0);
  }

  int x = 1;
  for (std::size_t i = 0; i < IntArray.size(0); ++i) {
    for (std::size_t j = 0; j < IntArray.size(1); ++j) {
      for (std::size_t k = 0; k < IntArray.size(2); ++k) {
        IntArray.at(i, j, k) = x;
        REQUIRE(IntArray.at(i, j, k) == x);
        REQUIRE(IntArray.at(i, j, k) == IntArray(i, j, k));
        ++x;
      }
    }
  }

  x = 1;
  for (std::size_t i = 0; i < IntArray.vector().size(); ++i) {
    REQUIRE(IntArray.vector().at(i) == x);
    ++x;
  }

  x = 1;
  for (std::size_t i = 0; i < IntArray.size(0); ++i) {
    for (std::size_t j = 0; j < IntArray.size(1); ++j) {
      for (std::size_t k = 0; k < IntArray.size(2); ++k) {

        IntArray.at(i, j, k) += int(i * j * k);
        REQUIRE(IntArray.at(i, j, k) == x + int(i * j * k));

        IntArray.at(i, j, k) -= 2 * int(i * j * k);
        REQUIRE(IntArray.at(i, j, k) == x - int(i * j * k));

        ++x;
      }
    }
  }

  const auto IntArray2 = IntArray;
  for (std::size_t i = 0; i < IntArray.vector().size(); ++i) {
    REQUIRE(IntArray2.vector().at(i) == IntArray.vector().at(i));
    REQUIRE(IntArray2.vector().at(i) == IntArray.data()[i]);
    REQUIRE(IntArray.vector().at(i) == IntArray2.data()[i]);
  }

  auto a0 = IntArray;
  const auto a1 = IntArray2;
  auto a3 = (1 + 2 * IntArray2 + IntArray2 - 2) * 3 / 2;
  auto a4 = IntArray2 * IntArray2;
  IntArray *= 2;
  auto a5 = IntArray2 - IntArray2;
  auto a6 = IntArray2 / IntArray2;

  for (std::size_t i = 0; i < IntArray.size(0); ++i) {
    for (std::size_t j = 0; j < IntArray.size(1); ++j) {
      for (std::size_t k = 0; k < IntArray.size(2); ++k) {
        REQUIRE(a3(i, j, k) == (1 + 2 * a1(i, j, k) + a1(i, j, k) - 2) * 3 / 2);
        REQUIRE(a4(i, j, k) == a1(i, j, k) * a1(i, j, k));
        REQUIRE(IntArray(i, j, k) == IntArray2(i, j, k) * 2);
        REQUIRE(a5(i, j, k) == 0);
        REQUIRE(a6(i, j, k) == 1);
      }
    }
  }

  {
    std::size_t count = 0;
    for (auto it = a0.begin(); it != a0.end(); ++it) {
      REQUIRE(*it == a0.vector()[count]);
      ++count;
    }
  }
  {
    std::size_t count = 0;
    for (auto it = a0.cbegin(); it != a0.cend(); ++it) {
      REQUIRE(*it == a0.vector()[count]);
      ++count;
    }
  }
  {
    std::size_t count = 0;
    for (auto it = a1.cbegin(); it != a1.cend(); ++it) {
      REQUIRE(*it == a1.vector()[count]);
      ++count;
    }
  }

  {
    std::size_t count = a0.size() - 1;
    for (auto it = a0.rbegin(); it != a0.rend(); ++it) {
      REQUIRE(*it == a0.vector()[count]);
      --count;
    }
  }
  {
    std::size_t count = a0.size() - 1;
    for (auto it = a0.crbegin(); it != a0.crend(); ++it) {
      REQUIRE(*it == a0.vector()[count]);
      --count;
    }
  }
  {
    std::size_t count = a1.size() - 1;
    for (auto it = a1.crbegin(); it != a1.crend(); ++it) {
      REQUIRE(*it == a1.vector()[count]);
      --count;
    }
  }

  a0.resize(0);
  REQUIRE(a0.size() == 0);
  REQUIRE(a0.dimensions() == 1);
  REQUIRE(a0.vector().size() == 0);

  a0.resize(3, 6, 2, 2);
  REQUIRE(a0.size() == 3 * 6 * 2 * 2);
  REQUIRE(a0.dimensions() == 4);
  REQUIRE(a0.vector().size() == 3 * 6 * 2 * 2);

  REQUIRE(a0.shape() == std::vector<std::size_t>{3, 6, 2, 2});

  for (const auto &el : a0) {
    REQUIRE(el == 0);
  }

  {
    qip::Array<double> dArray(2, 2);
    dArray += 2.0;
    dArray = dArray * 3.0;
    for (const auto &el : dArray) {
      REQUIRE(el == Approx(2.0 * 3.0));
    }
  }

  {
    qip::Array<std::complex<double>> cArray(2, 2);
    cArray += std::complex{2.0, 3.0};
    cArray = cArray * std::complex{3.0, -1.0};
    const auto tmp = std::complex{2.0, 3.0} * std::complex{3.0, -1.0};
    for (const auto &el : cArray) {
      REQUIRE(el.real() == Approx(tmp.real()));
      REQUIRE(el.imag() == Approx(tmp.imag()));
    }
  }

  // Fine, as long as you don't call * or /
  {
    qip::Array<std::string> sArray(2, 2);
    sArray += "a b c";
    for (const auto &el : sArray) {
      REQUIRE(el == "a b c");
    }
  }

  //----------------------------------------------------
  qip::Array<int> iArray(5, 7);

  int count = 10;
  for (std::size_t i = 0; i < iArray.size(0); ++i) {
    for (std::size_t j = 0; j < iArray.size(1); ++j) {
      iArray(i, j) = int(i + j) + count;
      ++count;
    }
  }
  const auto iArray2 = iArray;

  for (std::size_t row = 0; row < iArray.size(0); row++) {
    auto r = iArray.row(row);
    auto r2 = iArray2.row(row);
    for (std::size_t i = 0; i < r.size(); ++i) {
      REQUIRE(r(i) == iArray(row, i));
      REQUIRE(r2(i) == iArray2(row, i));
      r(i) *= 2;
      REQUIRE(r(i) == iArray(row, i));
      REQUIRE(r(i) == 2 * r2(i));
      REQUIRE(iArray(row, i) == 2 * iArray2(row, i));
    }
  }

  for (std::size_t col = 0; col < iArray.size(1); col++) {
    auto r = iArray.col(col);
    auto r2 = iArray2.col(col);
    for (std::size_t i = 0; i < r.size(); ++i) {
      REQUIRE(r(i) == iArray(i, col));
      REQUIRE(r2(i) == iArray2(i, col));
      r(i) *= 2;
      REQUIRE(r(i) == iArray(i, col));
      // 4x, because was 2x in the row loops
      REQUIRE(r(i) == 4 * r2(i));
      REQUIRE(iArray(i, col) == 4 * iArray2(i, col));
    }
  }

  //----------------------------------------------------
  // Test NDrange
  {
    const auto indexes = qip::NDrange(2, 1, 3);

    std::vector<std::array<std::size_t, 3>> expected{
        {0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {1, 0, 0}, {1, 0, 1}, {1, 0, 2}};

    REQUIRE(indexes == expected);
  }
  {
    const auto indexes = qip::NDrange(2, 2, 2, 2);

    std::vector<std::array<std::size_t, 4>> expected{
        {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 0, 1, 1},
        {0, 1, 0, 0}, {0, 1, 0, 1}, {0, 1, 1, 0}, {0, 1, 1, 1},
        {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 0, 1, 0}, {1, 0, 1, 1},
        {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 0}, {1, 1, 1, 1},
    };

    REQUIRE(indexes == expected);
  }
}
