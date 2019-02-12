#include "ExponentialGrid.h"
#include <cmath> //log, pow
#include <vector>

ExpGrid::ExpGrid(int in_N, double in_min, double in_max)
// constructor, includes some safety features/checks
{
  m_N = in_N;
  m_min = in_min;
  m_max = in_max;
  if (m_min > m_max) {
    // Do I want this? Or allow 'backwards' grid?
    m_min = in_max;
    m_max = in_min;
  }
  if (in_N <= 1) {
    // do I want to do this? Or throw an error?
    m_max = in_min;
    m_min = in_min;
    m_N = 1;
  }
  m_dxonx = log(m_max / m_min) / (m_N - 1);
  for (int i = 0; i < m_N; i++) {
    x_array.push_back(f_x(i));
    dxdi_array.push_back(f_x(i));
  }
}

int ExpGrid::N() const { return m_N; }
double ExpGrid::min() const { return m_min; }
double ExpGrid::max() const { return m_max; }
double ExpGrid::dxonx() const { return m_dxonx; }

const std::vector<double> &ExpGrid::x_grid() const { return x_array; }
const std::vector<double> &ExpGrid::dxdt_grid() const { return dxdi_array; }
double ExpGrid::dt() const { return m_dxonx; } // ONLY for exp....

double ExpGrid::x(int i) const {
  return x_array[i]; // no bounds-checking, don't waste time
}

double ExpGrid::dxdi(int i) const { return dxdi_array[i]; }

int ExpGrid::findNextIndex(double x) const
// Returns index correspoding to given value
// Note: finds NEXT LARGEST grid point (greater then or equal to.)
// Note: this is slow - don't rely on it inside big loops!
{
  double tmp = (m_N - 1) * log(x / m_min) / log(m_max / m_min);
  int i = (int)ceil(tmp);
  return i;
}

int ExpGrid::findNearestIndex(double x) const
// Note: this is slow - don't rely on it inside big loops!
// Create array if needed!
// Returns index correspoding to given value
// Uses 'round' - finds closest index!
{
  double tmp = (m_N - 1) * log(x / m_min) / log(m_max / m_min);
  int i = (int)round(tmp);
  return i;
}

double ExpGrid::f_x(int i) const {
  // Calculates value at given grid point
  if (m_N == 1)
    return m_min;
  double y = double(i) / (m_N - 1);
  return m_min * pow(m_max / m_min, y);
}

double ExpGrid::f_dxdi(int i) const { return m_dxonx * x(i); }
