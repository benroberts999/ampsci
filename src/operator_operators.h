#pragma once
#include <vector>

namespace operator{

  enum Operator {unity, r, gamma0, gamma5, dr, dr2};

  void operate(std::vector<double> &f, std::vector<double> &g, Operator v);
  void operate(std::vector<double> &f, std::vector<double> &g,
    const std::vector<double> &v);

  // void operate_r(std::vector<double> &f, std::vector<double> &g);
  // //Needs acess to grid! As does d/dr etc.!! Best way?
  void operate_gamma0(std::vector<double> &f, std::vector<double> &g);
  void operate_gamma5(std::vector<double> &f, std::vector<double> &g);

}//namespace


void operate(std::vector<double> &f, std::vector<double> &g, Operator v)
{

  switch(v){
    case unity  : break;
    case gamma0 : operate_gamma0(f,g); break;
    case gamma5 : operate_gamma5(f,g); break;
    default : std::cout<<"\nError18 in operator.cpp: unknown operator\n";
  }

}

void operate(std::vector<double> &f, std::vector<double> &g,
  const std::vector<double> &v)
{
  for(auto i=0u; i<f.size(); i++){
    f[i] *= v[i];
    g[i] *= v[i];
  }
}

void operate_gamma0(std::vector<double> &f, std::vector<double> &g){
  for(auto i=0u; i<f.size(); i++)
    g[i] = -g[i];
}

void operate_gamma5(std::vector<double> &f, std::vector<double> &g){
  for(auto i=0u; i<f.size(); i++){
    double fi = f[i];
    f[i] = g[i];
    g[i] = fi;
  }
}

// void operate_gamma(std::vector<double> &f, std::vector<double> &g,
//   const GammaMatrix &gm){
//
//   double a=gm.el[0][0];
//   double b=gm.el[0][1];
//   double c=gm.el[1][0];
//   double d=gm.el[1][1];
//
//   for(auto i=0u; i<f.size(); i++){
//     double fi = f[i];
//     f[i] = a*fi + b*g[i];
//     g[i] = c*fi + d*g[i];
//   }
// }
