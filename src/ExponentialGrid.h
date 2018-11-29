#pragma once
#include <vector>
/*
Simple struct to hold/define logarithmic/exponential grids.
In future, should be able to create options for other grid types!
*/

//******************************************************************************
class ExpGrid
/*
Simple struct to hold/define logarithmic/exponential grids.
Including Jacobian [dxdi := dx/di]
Note: dxonx := (dx/di)/x -- this is a constant! Used often
*/
{
public:


public:
  ExpGrid(int in_N, double in_min, double in_max);

  int N() const; //number of points in grid
  double min() const;
  double max() const;
  double dxonx() const; //:= (dx/di)/x -- this is a constant! Used often

  double x(int i) const;
  double dxdi(int i) const;
  int findNextIndex(double x) const;
  int findNearestIndex(double x) const;

private:
  int m_N; //number of points in grid
  double m_min;
  double m_max;
  double m_dxonx; //:= (dx/di)/x -- this is a constant! Used often
  std::vector<double> x_array;
  std::vector<double> dxdi_array;

private:
  double f_x(int i) const;
  double f_dxdi(int i) const;
};
