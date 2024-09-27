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

#include "Datapoint.h"

//*********CONSTRUCTOR & DESTRUCTOR*********//

Datapoint::Datapoint() :
TObject()
{
  fx = 0;
  fy = 0;
  fz = 0;
  fi = -1;
  fj = -1;
  fk = -1;
  fused = false;
  fneighbours = 0;
}

Datapoint::Datapoint(double x, double y, double z, int i, int j, int k) :
TObject()
{
  fx = x;
  fy = y;
  fz = z;
  fi = i;
  fj = j;
  fk = k;
  fused = false;
  fneighbours = 0;
}

Datapoint::Datapoint(Datapoint &p) :
TObject(p),
fx(p.fx),
fy(p.fy),
fz(p.fz),
fi(p.fi),
fj(p.fj),
fk(p.fk),
fused(p.fused),
fneighbours(p.fneighbours){}

//*********METHODS*********//

bool Datapoint::isNoise(int NB_MIN, int NB_MAX)
{
  //Condition about neighbours to consider points as noise
  return fneighbours < NB_MIN || fneighbours > NB_MAX;
}

bool Datapoint::isValid(int NB_MIN, int NB_MAX)
{
  return !fused && !isNoise(NB_MIN,NB_MAX);
}

int Datapoint::distance(Datapoint *p2)
{
  //Distance to consider a point as neighbour (considering idices)
  int i = abs(fi - p2->GetI());
  int j = abs(fj - p2->GetJ());
  int k = abs(fk - p2->GetK());

  return std::max(std::max(i, j), k);
}

double Datapoint::distance3D(Datapoint *p2)
{
  //Distance to consider a point as neighbour (considering idices)
  double x = (fx - p2->GetX());
  double y = (fy - p2->GetY());
  double z = (fz - p2->GetZ());
  double distance = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
  return distance;
}
