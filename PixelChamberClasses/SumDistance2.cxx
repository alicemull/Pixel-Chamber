#include <math.h>
#include <stdio.h>
#include <math.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::max
#include <stdlib.h>     /* abs */
#include <Math/Vector3D.h>
#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include <TVirtualFitter.h>
#include "TMatrixTUtils.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TRandom2.h"
#include "TLegend.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <TGraph2DErrors.h>
#include <TClonesArray.h>
#include <TPolyLine3D.h>
#include <TPolyLine.h>
#include <vector>
#include "TH3.h"
#include "TH2.h"
#include <TCanvas.h>
#include "TMatrixD.h"

using namespace ROOT::Math;

#include "SumDistance2.h"

//*********CONSTRUCTOR*********//

SumDistance2::SumDistance2(TGraph2DErrors *graph)
{
  fGraph = graph;
}

//*********METHOD & OPERATOR*********//

double SumDistance2::distance2(double x,double y,double z, double ex, double ey, double ez, const double *p)
{

  double x0 = p[0];
  double y0 = p[2];
  double z0 = p[4];
  double vy = p[1];
  double vz = p[3];
  double ux = p[5];

  double d2;
  double ym=y0+vy*(x-x0);
  double by = vy;
  //by=0.;
  double zm=z0+vz*(x-x0);
  double bz = vz;
  //bz=0.;

  double zm2=z0+vz/vy*(y-y0);
  double bz2 = vz/vy;

  // // d2 = (y-ym)*(y-ym)/(ey*ey+by*by*ex*ex)+(z-zm)*(z-zm)/(ez*ez+bz*bz*ex*ex);
  // if((x-x0) == 0.)
  // d2 = (z-zm2)*(z-zm2)/(ez*ez+bz2*bz2*ey*ey);
  // else
  //   d2 = (y-ym)*(y-ym)/(ey*ey+by*by*ex*ex) + (z-zm)*(z-zm)/(ez*ez+bz*bz*ex*ex);

  // if(ux != 0.) d2 = (y-ym)*(y-ym)/(ey*ey+by*by*ex*ex) + (z-zm)*(z-zm)/(ez*ez+bz*bz*ex*ex);
  // else
  //   d2 = (z-zm2)*(z-zm2)/(ez*ez+bz2*bz2*ey*ey);
  if(ux != 0.) d2 = (y-ym)*(y-ym)/(ey*ey) + (z-zm)*(z-zm)/(ez*ez);
  else
    d2 = (z-zm2)*(z-zm2)/(ez*ez);

  return d2;


  //*******************************************//
}

double SumDistance2::operator() (const double *par)
{
  //Operator that uses the function distance2 to set lines and errors parameters
  //and makes the sum of chi2 for all tracks considered
  bool first = true;
  assert(fGraph != 0);
  double * x = fGraph->GetX();
  double * y = fGraph->GetY();
  double * z = fGraph->GetZ();
  double * ex = fGraph->GetEX();
  double * ey = fGraph->GetEY();
  double * ez = fGraph->GetEZ();
  int npoints = fGraph->GetN();
  double sum = 0;
  for (int i  = 0; i < npoints; ++i)
  {
    double d = distance2(x[i],y[i],z[i],ex[i],ey[i],ez[i],par);
    sum += d;
  }
  first = false;
  return sum;
}
