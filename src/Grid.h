

enum class GridType {loglinear, logarithmic, linear};

/*

Use the keyword 'static' to declare the method:

static int MyMethod( int * a, int * b );
Then you can call the method without an instance like so:

int one = 1;
int two = 2;

MyClass::MyMethod( &two, &one );
'static' methods are functions which only use the class as a namespace, and do not require an instance.

The, grid MUST be called with ngp only (no h)!
But I can use h to generate ngp outside class if I want!


*/



//******************************************************************************
struct Grid{

  Grid(double in_r0, double in_rmax, int in_ngp,
    GridType gridtype = GridType::loglinear, double in_b=4.);

  // const GridType gridtype;
  const double r0;
  const double rmax;
  const int ngp;
  const double du;

  std::vector<double> r;
  std::vector<double> drdu;
  std::vector<double> drduor;

  void printGrid();

private:
  void formLogLinearGrid();
  void formLogGrid();
  void formLinearGrid();

  double calcdu(GridType gridtype, double in_r0, double in_rmax, int in_ngp,
    double b);
};

//******************************************************************************
Grid::Grid(double in_r0, double in_rmax, int in_ngp, GridType gridtype,
  double b)
  : r0(in_r0), rmax(in_rmax), ngp(in_ngp),
  du(calcdu(in_ngp, gridtype, b))
{

  r.clear();
  drdu.clear();
  drduor.clear();

  r.reserve(ngp);
  drdu.reserve(ngp);
  drduor.reserve(ngp);

  switch(gridtype){
    case GridType::loglinear : formLogLinearGrid(); break;
    case GridType::logarithmic : formLogGrid(); break;
    case GridType::loglinear : logLinearRadialGrid(); break;
    default : std::cerr<<"\n FAIL 49 in Grid: no grid type?\n";
  }

}




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
void Grid::form_loglinear_grid()
/*
Roughly, grid is logarithmically spaced below r=b, and linear above.
Definition:
  dr/du = r/(b+r)
  => u_0 = r0 + b*ln(r0)
  du = (in_rmax-in_r0+b*log(in_rmax/in_r0))/(in_ngp-1)
  du is constant (step-size for uniformly spaced grid)
Typically (and by default), b = 4 (unit units/bohr radius)
*/
{

  //initial points:
  r.push_back(r0);
  drduor.push_back(1./(b+r0));
  drdu.push_back(drduor[0]*r0);

  // Use iterative method from Dzuba code to calculate r grid
  double u = r0 + b*log(r0);
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
