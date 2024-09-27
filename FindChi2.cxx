#include "math.h"
#include "stdio.h"
#include "math.h"
#include "iostream"     // std::cout
#include "algorithm"    // std::max
#include "stdlib.h"     /* abs */
#include "Math/Vector3D.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "TVirtualFitter.h"
#include "TMatrixTUtils.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TRandom2.h"
#include "TLegend.h"
#include "iostream"
#include "string"
#include "fstream"
#include "sstream"
#include "TGraph2DErrors.h"
#include "TClonesArray.h"
#include "TPolyLine3D.h"
#include "TPolyLine.h"
#include "TH3.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TDecompChol.h"

using namespace ROOT::Math;

#include "Datapoint.h"
#include "Cluster.h"
#include "FindChi2.h"

//*********CONSTRUCTOR*********//

FindChi2::FindChi2(TClonesArray *tracks, double xm)
{
  fTracklets = tracks;
  fxm = xm;
  for(int i=0; i<fTracklets->GetEntries(); i++){
    fWT[i] = -9999.;
    fWTupdated[i] = -9999.;
  }
  //fWT = new double(fTracklets->GetEntries());
}

FindChi2::FindChi2(TClonesArray *tracks, double xm, double *WT)
{
  fTracklets = tracks;
  fxm = xm;
  for(int i=0; i<fTracklets->GetEntries(); i++){
    fWT[i] = WT[i];
    fWTupdated[i] = -9999.;
  }
  //fWT = new double(fTracklets->GetEntries());
}

//*********METHOD & OPERATOR*********//

double FindChi2::Chi2(int itr, double *q, TMatrixDSym *Cov, const double *p)
{

  //  printf("first fWT[%d]=%f\n",itr,fWT[itr]);
  double qh[2] = {q[0]-p[1],q[1]-p[2]};
  // printf("q[0]=%f q[1]=%f\n",q[0],q[1]);
  // printf("p[1]=%f p[2]=%f\n",p[1],p[2]);
  // printf("qh[0]=%f qh[1]=%f\n",qh[0],qh[1]);
  //Cov->Print();
  TMatrixD *QH = (TMatrixD*) new TMatrixD(2,1);
  for(int i=0; i<2; i++)(*QH)[i][0]=qh[i];
  TMatrixD *tQH = (TMatrixD*) new TMatrixD(1,2);
  tQH->Transpose(*QH);
  TMatrixD *mchi = (TMatrixD*) new TMatrixD(1,1);
  *mchi = (*tQH) * (*Cov) * (*QH);
  double chi2 = (*mchi)[0][0];

  // recalculate weights
//  printf("chi2=%f fWT[itr]=%f CT=%f\n",chi2,fWT[itr],fxm);
  if(chi2 < fxm)
    fWTupdated[itr] = pow(1-chi2/pow(fxm,2),2);
  else
    fWTupdated[itr] = 0.;

    chi2 *= fWT[itr];
  //  printf("fWTupdated[%d]=%f\n",itr,fWTupdated[itr]);

  //printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fxm= %f chi = %f\n",fxm,chi);
  // double chi = pow((y0-y0_model),2)/pow(sy0,2)+pow((z0-z0_model),2)/pow(sz0,2)+
  //   pow((tx-p[3]),2)/pow(stx,2)+pow((ty-p[4]),2)/pow(sty,2)+pow((tz-p[5]),2)/pow(stz,2);
  return chi2;
}

double FindChi2::operator()(const double *parameters)
{
  //Operator that uses the function Chi2 to set lines and errors parameters
  //and makes the sum of chi2 for all tracks considered
  int n_tracks = fTracklets->GetEntries();
  const double *parLine;
  const double *errorsLine;

  double CHISquare=0;
  int incr=0;
  double par[3];
  double q[2];
  double D;
  for(int i = 0; i<n_tracks; i++)
  {
    Cluster *track = (Cluster*)(fTracklets->At(i));
  //  printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! itr = %d ClusterID=%d ntr=%d\n",i,track->GetClusterID(),n_tracks);

    //printf("parameters[0]=%f\n",parameters[0]);
    parLine = track->GetParFit();
    errorsLine = track->GetParFitErrors();
    double alpha = parLine[3]/parLine[1];
    double y0 = parLine[2];
    double beta = parLine[5]/parLine[1];
    double z0 = parLine[4];
    double DeltaX = parameters[0] - parLine[0];
    //    if(i == 0) D=parameters[0] - parLine[0];
    D=parameters[0] - parLine[0];
    // printf("parameters[0]=%f parLine[0]=%f D=%f\n",parameters[0],parLine[0],D);
    // printf("parLine[0]=%f parLine[2]=%f parLine[4]=%f\n",parLine[0],parLine[2],parLine[4]);
    // printf("parLine[1]=%f parLine[3]=%f parLine[5]=%f\n",parLine[1],parLine[3],parLine[5]);

    //printf("parLine[0]=%f parLine[1]=%f parLine[2]=%f parLine[3]=%f parLine[4]=%f\n",parLine[0],parLine[1],parLine[2],parLine[3],parLine[4]);

    q[0] = y0 + DeltaX * alpha;
    q[1] = z0 + DeltaX * beta;

    TMatrixDSym *C = (TMatrixDSym*) track->GetCovarianceMatrix();
    // printf("original covariance matrix\n");
    // C->Print();

    TMatrixDSym *CovQ = (TMatrixDSym*) new TMatrixDSym(2);

    (*CovQ)[0][0] = (*C)[0][0]*pow(D,2) + 2.*(*C)[0][1]*D + (*C)[1][1];
    (*CovQ)[0][1] = (*C)[0][2]*pow(D,2) + (*C)[0][3]*D + (*C)[1][2]*D + (*C)[1][3];
    (*CovQ)[1][0] = (*CovQ)[0][1];
    (*CovQ)[1][1] = (*C)[2][2]*pow(D,2) + 2.*(*C)[2][3]*D + (*C)[3][3];

    TDecompChol *Decomp = new TDecompChol(*CovQ);
    Decomp->Invert(*CovQ);
    
    par[0]=parameters[0];
    par[1]=parameters[1];
    par[2]=parameters[2];
    // par[3]=parameters[3+2*incr];
    // par[4]=parameters[3+2*incr+1];

    double C2 = Chi2(i,q, CovQ, par);
    //double C2 = Chi2(q, C, par);

    CHISquare += C2;
    incr++;

  }
  return CHISquare;
}

void FindChi2::Clean()
{
  if(fTracklets) fTracklets->Clear("C");
}
