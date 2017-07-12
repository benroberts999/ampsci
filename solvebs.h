#ifndef _SOLVEBS_H
#define _SOLVEBS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include "funs.h"
#include "params.h"

#endif

int adamsmoulton(double *p, double *q, double *v, int ka, double &en, int ni, int nf);

int outint(double *p, double *q, double *v, int Z, int ka, double &en, int ctp);

int inint(double *p, double *q, double *v, int Z, int ka, double &en, int ctp, int pinf);

int solveDBS(double *p, double *q, double *v, int Z, int n, int ka, double &en, int &pinf, int &its, double &eps);

int AMcoefs(double *mia, double &mid, double &miaa);

int OIcoefs(double (*oie)[amo2], double *oia, double &oid);
