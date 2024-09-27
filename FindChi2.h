#ifndef FINDCHI2_H
#define FINDCHI2_H

class FindChi2
{
public:

  //*********CONSTRUCTOR AND DESTRUCTOR*********//

  FindChi2(TClonesArray *tracks, double xm);
  FindChi2(TClonesArray *tracks, double xm, double *WT);

  ~FindChi2(){};
  void Clean();
  double * GetWT(){return fWTupdated;}
  //*********METHOD & OPERATOR*********//

  double Chi2(int itr,double *q, TMatrixDSym *Cov, const double *p);
  double operator() (const double *parameters);
private:
  TClonesArray *fTracklets;
  double fxm;
  double fWT[200];
  double fWTupdated[200];
};

#endif
