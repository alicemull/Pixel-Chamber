#ifndef SUMDISTANCE2_H
#define SUMDISTANCE2_H

class SumDistance2
{
public:

  //*********CONSTRUCTOR AND DESTRUCTOR*********//

  SumDistance2(TGraph2DErrors *graph);
  ~SumDistance2(){};

  //*********METHOD & OPERATOR*********//

  double distance2(double x,double y,double z, double ex, double ey, double ez, const double *p);
  double operator() (const double *par);
private:
  TGraph2DErrors *fGraph;

};

#endif
