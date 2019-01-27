

enum class GridType {loglinear, logarithmic, linear};

//******************************************************************************
struct Grid{

  Grid();

  const double r0;
  const double rmax;
  const int ngp;
  const double du;

  std::vector<double> r;
  std::vector<double> drdu;
  std::vector<double> drduor;


private:
  void exponentialRadialGrid(int ngp_in, double r0, double rmax);
  void logLinearRadialGrid(int ngp_in, double r0, double rmax, double b=4.);
  void logLinearRadialGrid(double in_h, double r0, double rmax, double b=4.);

  double calcdu(GridType gridtype, double in_r0, double in_rmax, int in_ngp,
    double b);

};


//******************************************************************************
double calcdu(GridType gridtype, double in_r0, double in_rmax, int in_ngp,
  double b)
{
  switch(gridtype){
    case GridType::loglinear:
      return (in_rmax-in_r0+b*log(in_rmax/in_r0))/(in_ngp-1);
  }
}

//******************************************************************************
Grid::Grid(GridType gridtype, double in_r0, double in_rmax, int in_ngp,
  double b)
  : r0(in_r0), rmax(in_rmax), ngp(in_ngp),
  du(calcdu(gridtype, in_r0, in_rmax, in_ngp, b))
{

  r.clear();
  drdu.clear();
  drduor.clear();

  r.reserve(ngp);
  drdu.reserve(ngp);
  drduor.reserve(ngp);

}



void Grid::form_loglinear_grid()
{

  //initial points:
  r.push_back(r0);
  drduor.push_back(1./(b+r0));
  drdu.push_back(drduor[0]*r0);

  /*
  r(u[i]) = r0 ...... put equations!
  */

  // Use method from Dzuba code to calculate r grid
  double u = r0 + b*log(r0); //t is linear/uniform grid
  for (int i=1; i<ngp; i++) {
    u += du;
    double r_tmp = r[i-1];
    //Integrate dr/dt to find r:
    double delta_r = 1.;
    int ii = 0; //to count number of iterations
    while (delta_r > (r0*1.e-15)){
      double delta_u = u - (r_tmp + b*log(r_tmp));
      double drdu_tmp = r_tmp/(r_tmp + b);
      delta_r = delta_u * drdu_tmp;
      r_tmp += delta_r;
      ii++;
      if(ii>30) break; //usually converges in ~ 2 or 3 steps!
    }
    r.push_back(r_tmp);
    drduor.push_back(1./(b+r_tmp));
    drdu.push_back(drduor[i]*r_tmp);
  }
}

/*
// Lattice coordinate x = rmin + ih,
//      y = x + beta * log(rmin)
double y = rmin + beta*log(rmin) + h*double(i);
double r, rold;

if(beta <= y)
    r = y;
else
    r = exp(y/beta);

// Solve using Newton's method f(r) == 0, where
//    f(r) = r + beta*log(r/rmin) - x
//         = r + beta*log(r) - y
do
{   rold = r;
    r = r - (r + beta*log(r) - y)/(1. + beta/r);
    if(r <= 0.0)
        r = rold/2.0;
} while(fabs(r - rold)/r > 1.e-13);
*/
