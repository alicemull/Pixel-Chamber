#include "math.h"
#include "stdio.h"
#include "math.h"
#include "iostream"
#include "algorithm"
#include "stdlib.h"
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
#include "TDecompChol.h"
#include "TDecompLU.h"

using namespace ROOT::Math;

#include "Datapoint.h"
#include "FindChi2.h"
#include "SumDistance2.h"
#include "Cluster.h"

//*********CONSTRUCTOR & DESTRUCTOR*********//

Cluster::Cluster() :
TObject(),
fClusterID(0),
fndf(1),
fChi2(-1),
fReducedChi2(-1.),
fVx(-16.),
fVy(-16.),
fVz(-16.),
fEVx(0.),
fEVy(0.),
fEVz(0.),
fExPt1(0),
fExPt2(0),
fPDG_Code(0),
fPx(0),
fPy(0),
fPz(0),
fETot(0),
fPoints(0x0),
fGraph(0x0),
// fCovarianceMatrix(0x0),
// fInverseCovarianceMatrix(0x0),
fCovMatrixStatus(-1.),
fDecompositionFlag(kFALSE),
ft0(0),
ftf(0)
{
  fPoints = new TClonesArray("Datapoint",300);
  fGraph = new TGraph2DErrors();
  for(int i=0; i<6; i++)
  {
    fParFit[i] = 0.;
    fParFitErrors[i] = 0.;
  }
  for(int i=0; i<500; i++)
  {
    for(int j=0; j<11; j++)
    {
      fSecondaryVertices[i][j]=-16.00;
    }
  }
  fCovarianceMatrix = new TMatrixDSym(4);
  fInverseCovarianceMatrix = new TMatrixDSym(4);

}

Cluster::Cluster(int ClusterID) :
TObject()
{
  fClusterID = ClusterID;
  fPoints = new TClonesArray("Datapoint",300);
  fGraph = new TGraph2DErrors();
  fCovarianceMatrix = new TMatrixDSym(4);
  fInverseCovarianceMatrix = new TMatrixDSym(4);
  fCovMatrixStatus=-1.;
  fDecompositionFlag=kFALSE;
  fndf = 1.;
  fChi2 = -1.;
  fReducedChi2 = -1.;
  fVx=-16.;
  fVy=-16.;
  fVz=-16.;
  fEVx = 0.;
  fEVy = 0.;
  fEVz = 0.;
  fPx = 0;
  fPy = 0;
  fPz = 0;
  fETot = 0;
  fExPt1 = 0;
  fExPt2 = 0;
  fPDG_Code=0;
  for(int i=0; i<6; i++)
  {
    fParFit[i] = 0.;
    fParFitErrors[i] = 0.;
  }
  for(int i=0; i<500; i++)
  {
    for(int j=0; j<11; j++)
    {
      fSecondaryVertices[i][j]=-16.00;
    }
  }
  ft0=0.;
  ftf=0.;
}

Cluster::Cluster(Cluster &c) :
TObject(c),
fClusterID(c.fClusterID),
fndf(c.fndf),
fChi2(c.fChi2),
fReducedChi2(c.fReducedChi2),
fVx(c.fVx),
fVy(c.fVy),
fVz(c.fVz),
fEVx(c.fEVx),
fEVy(c.fEVy),
fEVz(c.fEVz),
fExPt1(c.fExPt1),
fExPt2(c.fExPt2),
fPx(c.fPx),
fPy(c.fPy),
fPz(c.fPz),
fETot(c.fETot),
fPDG_Code(c.fPDG_Code),
fCovMatrixStatus(c.fCovMatrixStatus),
fDecompositionFlag(c.fDecompositionFlag),
ft0(c.ft0),
ftf(c.ftf){

  fGraph = new TGraph2DErrors();
  fCovarianceMatrix = (TMatrixDSym*) c.fCovarianceMatrix->Clone();
  fInverseCovarianceMatrix = (TMatrixDSym*) c.fInverseCovarianceMatrix->Clone();

  fPoints = new TClonesArray("Datapoint",300);
  TClonesArray *fnp =  c.GetPointsObj();
  for(int i=0; i<fnp->GetEntries();i++) {
    Datapoint *p = (Datapoint*) fnp->At(i);
    this->AddDataPoint(p);
  }

  for (Int_t i=0; i<6; i++) {
    fParFit[i] = c.fParFit[i];
    fParFitErrors[i] = c.fParFitErrors[i];
  }
  for(int i=0; i<500; i++)
  {
    for(int j=0; j<11; j++)
    {
      fSecondaryVertices[i][j]=c.fSecondaryVertices[i][j];
    }
  }
}

Cluster::~Cluster()
{
  fPoints->Clear("C");
  delete fPoints;
  delete fGraph;
  delete fCovarianceMatrix;
  delete fInverseCovarianceMatrix;
  // delete[] fSecondaryVertices;
  
}

void Cluster::PrintCluster()
{
  printf("\n");
  printf("********* Printing cluster info *********\n");
  printf("\n");
  printf("fClusterID=%d\n",fClusterID);
  Datapoint *p;
  printf("fPoints->GetEntries()=%d\n",fPoints->GetEntries());

  int npoints = fGraph->GetN();
  printf("graph npoints = %d\n",npoints);
  printf("fndf=%d\n",fndf);
  printf("fChi2=%f\n",fChi2);
  printf("fReducedChi2=%f\n",fReducedChi2);
  printf("fParFit[0]=%f\n",fParFit[0]);
  printf("fParFit[1]=%f\n",fParFit[1]);
  printf("fParFit[2]=%f\n",fParFit[2]);
  printf("fParFit[3]=%f\n",fParFit[3]);
  printf("fParFit[4]=%f\n",fParFit[4]);
  printf("fParFit[5]=%f\n",fParFit[5]);
  printf("fParFitErrors[0]=%f\n",fParFitErrors[0]);
  printf("fParFitErrors[1]=%f\n",fParFitErrors[1]);
  printf("fParFitErrors[2]=%f\n",fParFitErrors[2]);
  printf("fParFitErrors[3]=%f\n",fParFitErrors[3]);
  printf("fParFitErrors[4]=%f\n",fParFitErrors[4]);
  printf("fParFitErrors[5]=%f\n",fParFitErrors[5]);
  printf("fExPt1=%d\n",fExPt1);
  printf("fExPt2=%d\n",fExPt2);
  printf("ft0=%f\n",ft0);
  printf("ftf=%f\n",ftf);
  printf("fVx=%f\n",fVx);
  printf("fVy=%f\n",fVy);
  printf("fVz=%f\n",fVz);
  printf("fEVx=%f\n",fEVx);
  printf("fEVy=%f\n",fEVy);
  printf("fEVz=%f\n",fEVz);
  printf("fPDG_Code=%d\n",fPDG_Code);
  printf("fPx=%f\n",fPx);
  printf("fPy=%f\n",fPy);
  printf("fPz=%f\n",fPz);
  printf("fETot=%f\n",fETot);

  printf("Covariance matrix\n");
  fCovarianceMatrix->Print();
  printf("Inverse covariance matrix\n");
  fInverseCovarianceMatrix->Print();

}
//*********SOME USEFULL METHODS*********//

void Cluster::AddDataPoint(Datapoint *p)
{
  TClonesArray &arrTemp  = *fPoints;
  new (arrTemp[arrTemp.GetEntriesFast()]) Datapoint(*p);
}

void Cluster::RemoveDataPointAndCompress(int ip)
{
  fPoints->RemoveAt(ip);
  fPoints->Compress();
}

void Cluster::SortCluster()
{
  TClonesArray &arrTemp  = *fPoints;
  //pts->Sort();

  int index=fPoints->GetEntries()-1;
  printf("pts->GetEntries()=%d\n", fPoints->GetEntries());
  for(int l=0; l<index; l++)
  {
    for(int m=0; m<index-l; m++)
    {
      // {
      Datapoint *p1 = (Datapoint*)(fPoints->At(m));
      Datapoint *p2 = (Datapoint*)(fPoints->At(m+1));

      if(p1->GetI()>p2->GetI())
      {
        Datapoint*p0=(Datapoint*)p1->Clone();
        fPoints->RemoveAt(m);
        new (arrTemp[m]) Datapoint(*p2);
        fPoints->RemoveAt(m+1);
        new (arrTemp[m+1]) Datapoint(*p0);
      }
      if(p1->GetI()==p2->GetI()&&p1->GetJ()>p2->GetJ())
      {
        Datapoint*p0=(Datapoint*)p1->Clone();
        fPoints->RemoveAt(m);
        new (arrTemp[m]) Datapoint(*p2);
        fPoints->RemoveAt(m+1);
        new (arrTemp[m+1]) Datapoint(*p0);
      }
      if(p1->GetI()==p2->GetI()&&p1->GetJ()==p2->GetJ()&&p1->GetK()>p2->GetK())
      {
        Datapoint*p0=(Datapoint*)p1->Clone();
        fPoints->RemoveAt(m);
        new (arrTemp[m]) Datapoint(*p2);
        fPoints->RemoveAt(m+1);
        new (arrTemp[m+1]) Datapoint(*p0);
      }
      // }

    }
  }
  //printf("sort->GetEntries()=%d\n", sort->GetEntries());
}

int Cluster::findFarestPoint(int j, double *xg, double *yg, double *zg, int np)
{
  //Method to find the farest point from one point. This is used to find the extreme points of a cluster
  int ifar = np-1;
  double maxdst = 0;
  for(int i=0; i<np;i++)
  {
    double dst = sqrt((xg[i]-xg[j])*(xg[i]-xg[j])+(yg[i]-yg[j])*(yg[i]-yg[j])+(zg[i]-zg[j])*(zg[i]-zg[j]));
    if (dst > maxdst)
    {
      ifar = i;
      maxdst = dst;
    }
  }
  return ifar;
}

bool Cluster::isNeighbour(Datapoint *p, int NB_MIN, int NB_MAX)
{
  //Method to find neighbours using methods of Datapoint
  if (!p->isValid(NB_MIN,NB_MAX))
  return false;
  //Loop on all points
  for (int n = 0; n < fPoints->GetEntries(); n++)
  {
    Datapoint *p0 = (Datapoint*)(fPoints->At(n));
    //To consider a point as neighbour the distance (in terms of indices) has to be smaller than EPS
    if (p0->distance(p) <= EPS)
    return true;
  }
  return false;
}

//*********FIT METHODS*********//

void Cluster::makeGraph(bool kError)
{
  //SumDistance2 operator needs a graph as input so every cluster is in a TGraph
  int ip = 0;
  double ex, ey, ez;
  //first errors are the semi pixel pitch for the first fit
  //second errors are the pixel pitch over sqrt(12) for the final fit
  if(!kError)
  {
    ex=0.01462;
    ey=0.025;
    ez=0.01344;
  }
  else
  {
    ex=0.02924/sqrt(12.);
    ey=0.050/sqrt(12.);
    ez=0.02688/sqrt(12.);
  }
  if (fGraph) delete fGraph;
  fGraph = new TGraph2DErrors();
  for(int m = 0; m < fPoints->GetEntries(); m++)
  {
    Datapoint *pl = (Datapoint*)(fPoints->At(m));
    // printf("pl->GetX()=%f, pl->GetY()=%f\n",pl->GetX(), pl->GetY());
    fGraph->SetPoint(ip, pl->GetX(), pl->GetY(), pl->GetZ());
    fGraph->SetPointError(ip, ex, ey, ez);
    ip++;
  }
}

void Cluster::findExtremes()
{
  // finds extreme points of a cluster
  double *xg = fGraph->GetX();
  double *yg = fGraph->GetY();
  double *zg = fGraph->GetZ();
  int np = fGraph->GetN();
  //
  fExPt1 = findFarestPoint( 0, xg, yg, zg, np);
  fExPt2 = findFarestPoint( fExPt1, xg, yg, zg, np);
}

void Cluster::fitCluster(double xstart)
{
  // Root fitter: by default minuit algorithm is used
  ROOT::Fit::Fitter  fitter;
  // make the functor object
  SumDistance2 sdist(fGraph);
  // set the function
  ROOT::Math::Functor fcn(sdist,6);
  double *xg = fGraph->GetX();
  double *yg = fGraph->GetY();
  double *zg = fGraph->GetZ();
  int np = fGraph->GetN();
  double maxX=-999., maxY=-999., maxZ=-999.;
  // finds extreme points of a cluster
  findExtremes();
  // set initial parameter values
  double ux = xg[fExPt2] - xg[fExPt1];
  // double uy = (yg[fExPt2] - yg[fExPt1]) / (ux + 0.000001);
  // double uz = (zg[fExPt2] - zg[fExPt1]) / (ux + 0.000001);
  double uy = (yg[fExPt2] - yg[fExPt1]);
  double uz = (zg[fExPt2] - zg[fExPt1]);
  double pStart[6];
  // printf("fExPt1=%d, fExPt2=%d\n", fExPt1, fExPt2);
  // printf("xg[fExPt1]=%f, xg[fExPt2]=%f\n",xg[fExPt1], xg[fExPt2] );
  // printf("yg[fExPt1]=%f, yg[fExPt2]=%f\n",yg[fExPt1], yg[fExPt2] );
  // printf("zg[fExPt1]=%f, zg[fExPt2]=%f\n",zg[fExPt1], zg[fExPt2] );
  // printf("xstart=%f ystart=%f zstart=%f\n",xg[fExPt1],yg[fExPt1],zg[fExPt1]);
  // printf("ux=%f uy=%f uz=%f\n",ux,uy,uz);
  if(fParFit[1] == 0. & fParFit[3] == 0.)
  {
    if(xstart == -999.) pStart[0] = xg[fExPt1];
    else
    pStart[0] = xstart;
    //    pStart[0] = 5;
    pStart[1] = uy;
    pStart[2] = yg[fExPt1];
    pStart[3] = uz;
    pStart[4] = zg[fExPt1];
    pStart[5] = ux;
  }
  else
  {

    if(xstart == -999.) pStart[0] = xg[fExPt1];
    else
    pStart[0] = xstart;
    //pStart[0] = 3.164934;
    pStart[1] = uy;
    pStart[2] = yg[fExPt1];
    pStart[3] = uz;
    pStart[4] = zg[fExPt1];
    pStart[5] = ux;
  }
  // for(int i=0; i<6; i++)printf("pStart=%f\n",pStart[i]);

  fitter.Config().SetMinimizer("Minuit2","Migrad");
  if(ux != 0.) fitter.SetFCN(fcn,pStart,np*2.,1);
  else
  fitter.SetFCN(fcn,pStart,np,1);

  //  fitter.SetFCN(fcn,pStart,np*3.,1);
  for (int i = 0; i < 5; ++i)
  {
    fitter.Config().ParSettings(i).SetStepSize(0.001);
  }
  fitter.Config().ParSettings(0).Fix();
  fitter.Config().ParSettings(5).Fix();

  bool ok = fitter.FitFCN();
  // Taking fit results and saving as private members of the class
  const ROOT::Fit::FitResult & result = fitter.Result();
  double chi2 = result.Chi2();
  double ndf = result.Ndf();
  double reducedChi2 = chi2/ndf;
  double norm = 1.;
  for(int cov_i=1; cov_i<5; cov_i++)
  {
    for(int cov_j=1; cov_j<5; cov_j++)
    {
      (*fCovarianceMatrix)[cov_i-1][cov_j-1]=result.CovMatrix(cov_i,cov_j);
      (*fInverseCovarianceMatrix)[cov_i-1][cov_j-1]=result.CovMatrix(cov_i,cov_j);
    }
  }
  fCovMatrixStatus = result.CovMatrixStatus();
  //  fCovarianceMatrix->Print();
  //fInverseCovarianceMatrix=(TMatrixDSym*)fCovarianceMatrix->Clone();

  TDecompLU lu(*fInverseCovarianceMatrix);
  fDecompositionFlag=lu.Decompose();
  if(!lu.Decompose())
  {
    for(int cov_i=0; cov_i<4; cov_i++)
    {
      for(int cov_j=0; cov_j<4; cov_j++)
      {
        (*fInverseCovarianceMatrix)[cov_i][cov_j]=0.;
      }
    }
  }else{
    TDecompChol *Decomp = new TDecompChol(*fInverseCovarianceMatrix);
    Decomp->Invert(*fInverseCovarianceMatrix);
  }

  fChi2 = chi2*norm;
  fndf = ndf*norm;
  fReducedChi2 = reducedChi2*norm;

  const double *ParFit;
  ParFit = result.GetParams();

  const double *ParFitErrors;
  ParFitErrors = result.GetErrors();

  double vx = 1. / sqrt(1. + ParFit[1]*ParFit[1]+ParFit[3]*ParFit[3]);
  double vy = ParFit[1] * vx;
  double vz = ParFit[3] * vx;

  double evx = sqrt(pow(ParFit[1]*pow(vx,3),2)*pow(ParFitErrors[1],2) + pow(ParFit[3]*pow(vx,3),2)*pow(ParFitErrors[3],2));
  double evy = sqrt(pow(vx,2)*pow(ParFitErrors[1],2) + pow(ParFit[1],2)*pow(evx,2));
  double evz = sqrt(pow(vx,2)*pow(ParFitErrors[3],2) + pow(ParFit[3],2)*pow(evx,2));

  fParFit[0] = norm*ParFit[0];
  fParFit[1] = vx;
  fParFit[2] = norm*ParFit[2];
  fParFit[3] = vy;
  fParFit[4] = norm*ParFit[4];
  fParFit[5] = vz;

  fParFitErrors[0] = norm*ParFitErrors[0];
  fParFitErrors[1] = evx;
  fParFitErrors[2] = norm*ParFitErrors[2];
  fParFitErrors[3] = evy;
  fParFitErrors[4] = norm*ParFitErrors[4];
  fParFitErrors[5] = evz;

  // printf("chi2=%f ndf=%f\n",chi2,ndf);
  // printf("CovMatrixStatus=%d\n",fCovMatrixStatus);
  // printf(" npoints=%d\n",np);
  // printf("x0=%f+-%f alfa=%f+-%f y0=%f+-%f beta=%f+-%f z0=%f+-%f\n",ParFit[0],ParFitErrors[0],ParFit[1],ParFitErrors[1],ParFit[2],ParFitErrors[2],ParFit[3],ParFitErrors[3],ParFit[4],ParFitErrors[4]);
  //   printf("CovMatrixStatus = %d\n",fCovMatrixStatus);
  if(fCovMatrixStatus < 2) {
    printf("!!!!!!!!!!!!!!!!!!!!! Covariance matrix not positive definite !!!!!!!\n");
    // printf("chi2=%f ndf=%f\n",chi2,ndf);
    // printf("x0=%f+-%f alfa=%f+-%f y0=%f+-%f beta=%f+-%f z0=%f+-%f\n",ParFit[0],ParFitErrors[0],ParFit[1],ParFitErrors[1],ParFit[2],ParFitErrors[2],ParFit[3],ParFitErrors[3],ParFit[4],ParFitErrors[4]);
    // printf("CovMatrixStatus = %d\n",fCovMatrixStatus);
    // fCovarianceMatrix->Print();
    // printf("parameters\n");
    // printf("%15.12f\n",fParFit[0]);
    // printf("%15.12f\n",fParFit[1]);
    // printf("%15.12f\n",fParFit[2]);
    // printf("%15.12f\n",fParFit[3]);
    // printf("%15.12f\n",fParFit[4]);
    // printf("%15.12f\n",fParFit[5]);
    //
    // printf("errors\n");
    // printf("%15.12f\n",fParFitErrors[0]);
    // printf("%15.12f\n",fParFitErrors[1]);
    // printf("%15.12f\n",fParFitErrors[2]);
    // printf("%15.12f\n",fParFitErrors[3]);
    // printf("%15.12f\n",fParFitErrors[4]);
    // printf("%15.12f\n",fParFitErrors[5]);
    //
    // printf("Extremes\n");
    // printf("%d\n",fExPt1);
    // printf("%d\n",fExPt2);
    //
    // printf("xg\n");
    // for(int ii=0;ii<np;ii++){
    //   printf("%f,\n",xg[ii]);
    // }
    // printf("yg\n");
    // for(int ii=0;ii<np;ii++){
    //   printf("%f,\n",yg[ii]);
    // }
    // printf("zg\n");
    // for(int ii=0;ii<np;ii++){
    //   printf("%f,\n",zg[ii]);
    // }
  }

  int n = 10000;

  double ParLine[6]={xg[fExPt1],fParFit[1],yg[fExPt1],fParFit[3],zg[fExPt1],fParFit[5]};
  for (int i = 0; i <np;++i)
  {
    if(fParFit[1]!=0.){
      ft0=(xg[fExPt1]-fParFit[0])/fParFit[1];
      ftf=(xg[fExPt2]-fParFit[0])/fParFit[1];
    } else if(fParFit[3]!=0.){
      ft0=(yg[fExPt1]-fParFit[2])/fParFit[3];
      ftf=(yg[fExPt2]-fParFit[2])/fParFit[3];
    } else if(fParFit[5]!=0.){
      ft0=(zg[fExPt1]-fParFit[4])/fParFit[5];
      ftf=(zg[fExPt2]-fParFit[4])/fParFit[5];
    }
    double t = ft0 + (ftf-ft0)*i/(np-1.);
    double x,y,z;
    line(t,fParFit,x,y,z);
  }

}


// double Cluster::KalmanFilter()
void Cluster::KalmanFilter()
{
  // printf("@@@@@Starting with Cluster %d with %d points@@@@@\n", this->GetClusterID(), this->GetNumPoints());
  this->makeGraph(1);
  this->fitCluster();
  TClonesArray *pts = (TClonesArray*)(this->GetPointsObj());
  int NbPts=pts->GetEntries()-1;
  const double *parFit;
  parFit = this->GetParFit();
  TMatrixD *CovMatrix = new TMatrixD(4,4);
  TMatrixD *CovMatrix_k1 = new TMatrixD(4,4);
  TMatrixD *NewCovMatrix_k1 = new TMatrixD(4,4);
  TMatrixD *CovMatrixI = new TMatrixD(4,4);
  TMatrixD *CovMatrixI_k1 = new TMatrixD(4,4);
  TMatrixD *NewCovMatrixI_k1 = new TMatrixD(4,4);

  Datapoint *Hit_k = new Datapoint();
  Datapoint *Hit_k1 = new Datapoint();
  TMatrixD *SV = new TMatrixD(4,1);//State Vector
  (*SV)[0][0]=parFit[3]/parFit[1];
  (*SV)[1][0]=parFit[2];
  (*SV)[2][0]=parFit[5]/parFit[1];
  (*SV)[3][0]=parFit[4];
  TMatrixD *V=new TMatrixD(2,2);
  (*V)[0][0]=0.050/sqrt(12.);
  (*V)[1][1]=0.02688/sqrt(12.);
  TMatrixD *VI=new TMatrixD(2,2);
  VI=(TMatrixD*)V->Clone();
  double DetV = VI->Determinant();
  VI->Invert(&DetV);
  TMatrixD *MS = new TMatrixD(4,4);//Multiple Scattering matrix

  TMatrixD *p = new TMatrixD(4,1);//Parameters Vector
  TMatrixD *p_start = new TMatrixD(4,1);//Parameters Vector
  TMatrixD *p_k1 = new TMatrixD(4,1);//Updated Parameters Vector
  TMatrixD *f_k1 = new TMatrixD(4,4);//Extrapolation Matrix
  TMatrixD *pTilde_k1 = new TMatrixD(4,1);//Extrapolated Parameters Vector
  TMatrixD *D = new TMatrixD(4,4);//Trasformation matrix for covariance matrix
  TMatrixD *Dt = new TMatrixD(4,4);
  TMatrixD *H = new TMatrixD(2,4);
  TMatrixD *dChi2 = new TMatrixD(1,1);
  double L=0.;
  double L_r=93.7; //mm
  int z=1;
  double c=2.998e8;
  double P = TMath::Sqrt(pow(this->GetPx(),2)+ pow(this->GetPy(),2)+ pow(this->GetPz(),2));
  // printf("cluster %d has %f GeV of momentum\n", this->GetClusterID(), P);
  double theta_proj=0.;
  double chi2=0.;
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<4; j++)
    {
      (*H)[i][j]=0;
      if(i==0 && j==1)(*H)[i][j]=1;
      if(i==1 && j==3)(*H)[i][j]=1;
      // else(*H)[i][j]=0;
    }
  }
  // printf("Matrix HT\n");
  TMatrixD *Ht = new TMatrixD(4,2);
  Ht->Transpose(*H);
  TMatrixD *m=new TMatrixD(2,1);
  TMatrixD *r=new TMatrixD(2,1);
  TMatrixD *rt=new TMatrixD(1,2);
  for(int k=NbPts; k>0; k--)
  {
    // printf("#####k=%d\n",k);
    // printf("Sono qui 1\n");
    Hit_k=(Datapoint*)(pts->At(k));
    Hit_k1=(Datapoint*)(pts->At(k-1));
    L=Hit_k->distance3D(Hit_k1);
    theta_proj=pow((13.6/(P))*TMath::Sqrt(L/L_r)*(1+0.0038*TMath::Log(L/L_r)),2);
    (*m)[0][0]=Hit_k1->GetY();
    (*m)[1][0]=Hit_k1->GetZ();
    //printf("imposto *p\n");
    if(k==NbPts)
    {
      //STARTING VALUES FOR COV MATRIX (C_K) AND PTILDE_K
      //printf("primo if\n");
      (*p)[0][0]=parFit[3]/parFit[1];
      (*p)[1][0]=parFit[2];
      (*p)[2][0]=parFit[5]/parFit[1];
      (*p)[3][0]=parFit[4];
      (*p_start)[0][0]=parFit[3]/parFit[1];
      (*p_start)[1][0]=parFit[2];
      (*p_start)[2][0]=parFit[5]/parFit[1];
      (*p_start)[3][0]=parFit[4];

      // printf("Sono qui 2\n");
      // printf("Matrix Mat\n");
      TMatrixDSym* Mat=(TMatrixDSym*)this->GetCovarianceMatrix();
      // printf("Mat->Determinant()=%f\n", Mat->Determinant());
      if(Mat->Determinant()==0)break;
      // printf("Matrix MatI\n");
      TMatrixDSym* MatI=(TMatrixDSym*)this->GetInverseCovarianceMatrix();
      // printf("trackl->GetCovMatrixStatus()=%f\n", trackl->GetCovMatrixStatus());
      // printf("trackm->GetCovMatrixStatus()=%f\n", trackm->GetCovMatrixStatus());
      // if(trackm->GetDecompositionFlag()==false)printf("Ciaone\n" );

      for(int i=0; i<4; i++)
      {
        for(int j=0; j<4; j++)
        {
          (*CovMatrix)[i][j]=(*Mat)[i][j];
          (*CovMatrixI)[i][j]=(*MatI)[i][j];
        }
      }
      // printf("Starting Covariance matrix\n");
      // Mat->Print();
    }
    else
    {
      //printf("secondo if\n");
      //UPDATES VALUES FOR COV MATRIX (C_K) AND PTILDE_K
      // printf("Matrix CovMatrix\n");
      CovMatrix=(TMatrixD*)NewCovMatrix_k1->Clone();
      // printf("Matrix CovMatrixI\n");
      CovMatrixI=(TMatrixD*)CovMatrix->Clone();
      double DetCovMatrixI = CovMatrixI->Determinant();
      // printf("Sono qui 3\n");
      CovMatrixI->Invert(&DetCovMatrixI);
      // printf("Matrix p_k1\n");
      p = (TMatrixD*)p_k1->Clone();
    }
    // printf("*p\n");
    // p->Print();
    // printf("*CovMatrix\n");
    // CovMatrix->Print();
    for(int i=0; i<4; i++)
    {
      for(int j=0; j<4; j++)
      {
        //F_K+1
        (*f_k1)[i][j]=0;
        if(i==j)(*f_k1)[i][j]=1;
        if(i==1 && j==0)(*f_k1)[i][j]=Hit_k1->GetX()-Hit_k->GetX();
        if(i==3 && j==2)(*f_k1)[i][j]=Hit_k1->GetX()-Hit_k->GetX();

      }
    }
    //printf("Hit_k1->GetX()-Hit_k->GetX()=%f\n",Hit_k1->GetX()-Hit_k->GetX());
    //printf("imposto *pTilde_k1\n");
    (*pTilde_k1)=(*f_k1)*(*p);
    //printf("imposto *D\n");
    D->Transpose(*f_k1);
    Dt->Transpose(*D);
    //printf("imposto *MS\n");
    for(int i = 0; i<4; i++)
    {
      for(int j = 0; j<4; j++)
      {
        (*MS)[i][j]=0;
        //MULTIPLE SCATTERING MATRIX
        if(i==0 && j==0)(*MS)[i][j]=theta_proj*(1+pow((*pTilde_k1)[0][0],2))*(1+pow((*pTilde_k1)[0][0],2)+pow((*pTilde_k1)[2][0],2));
        if(i==2 && j==2)(*MS)[i][j]=theta_proj*(1+pow((*pTilde_k1)[2][0],2))*(1+pow((*pTilde_k1)[0][0],2)+pow((*pTilde_k1)[2][0],2));
        if(i==0 && j==2)(*MS)[i][j]=theta_proj*(*pTilde_k1)[0][0]*(*pTilde_k1)[2][0]*(1+pow((*pTilde_k1)[0][0],2)+pow((*pTilde_k1)[2][0],2));
        if(i==2 && j==0)(*MS)[i][j]=theta_proj*(*pTilde_k1)[0][0]*(*pTilde_k1)[2][0]*(1+pow((*pTilde_k1)[0][0],2)+pow((*pTilde_k1)[2][0],2));
      }
    }
    //MS[0][2]=MS[2][0];
    //COVARIANCE MATRIX FOR K-1 THAT TAKES IN ACCOUNT MS_K-1
    // printf("Matrix CovMatrix_k1\n");
    (*CovMatrix_k1) = (*D)*(*CovMatrix)*(*Dt)+(*MS);
    // printf("CovMatrix\n");
    // CovMatrix->Print();
    // printf("D\n");
    // D->Print();
    // printf("MS\n");
    // MS->Print();
    // printf("CovMatrix_k1\n");
    // CovMatrix_k1->Print();
    CovMatrixI_k1=(TMatrixD*)CovMatrix_k1->Clone();
    double detCovMatrixI_k1 = CovMatrixI_k1->Determinant();
    // printf("Sono qui 4\n");
    // printf("Matrix CovMatrixI_k1\n");
    CovMatrixI_k1->Invert(&detCovMatrixI_k1);
    //printf("imposto *prodMatrix\n");
    // printf("Matrix prodMatrix\n");
    TMatrixD *prodMatrix = new TMatrixD(4,4);
    (*prodMatrix)=(*CovMatrixI_k1)+(*Ht)*(*VI)*(*H);
    double det=prodMatrix->Determinant();
    // printf("Sono qui 5\n");
    // printf("Matrix prodMatrixI\n");
    prodMatrix->Invert(&det);

    //VALUE P_K-1 THAT WILL UPDATE THE *P_K
    (*p_k1)=(*prodMatrix)*((*Ht)*(*VI)*(*m)+(*CovMatrixI_k1)*(*pTilde_k1));
    // printf("p_k1\n");
    // p_k1->Print();
    // printf("Matrix NewCovMatrixI_k1\n");
    (*NewCovMatrixI_k1)=(*CovMatrixI_k1)+(*Ht)*(*VI)*(*H);
    // printf("NewCovMatrix_k1\n");
    // NewCovMatrix_k1->Print();
    //VALUE OF COVARIANCE MATRIX_K-1 THAT WILL UPDATE THE COVMATRIX_K
    NewCovMatrix_k1=(TMatrixD*)NewCovMatrixI_k1->Clone();
    double detNewCovMatrix_k1 = NewCovMatrix_k1->Determinant();
    // printf("Sono qui 6\n");
    // printf("Matrix NewCovMatrixI_k1\n");
    NewCovMatrix_k1->Invert(&detNewCovMatrix_k1);
    //printf("imposto *r\n");
    //RESIDUAL VECTOR
    (*r)=(*m)-(*H)*(*p_k1);
    rt->Transpose(*r);
    //printf("imposto *Diff\n");
    TMatrixD *Diff = new TMatrixD(4,1);
    TMatrixD *Difft = new TMatrixD(1,4);
    (*Diff)=(*p_k1)-(*pTilde_k1);
    Difft->Transpose(*Diff);
    //DELTACHI2
    (*dChi2)=(*rt)*(*VI)*(*r)+(*Difft)*(*NewCovMatrixI_k1)*(*Diff);

    chi2 += (*dChi2)[0][0];
  }
  //SAVING FIT PARAMETERS
  for(int i=0; i<4; i++)
  {
    for(int j=0; j<4; j++)
    {
      (*fCovarianceMatrix)[i][j]=(*NewCovMatrix_k1)[i][j];
      (*fInverseCovarianceMatrix)[i][j]=(*NewCovMatrixI_k1)[i][j];
    }
  }

  // Datapoint* p;
  // for(int i=0; i<this->GetNumPoints(); i++)
  // {
  //   p=(Datapoint*)fPoints->At(i);
  //   if()
  // }
  fParFit[1]=1./(1+(*p_k1)[0][0]+(*p_k1)[2][0]);
  fParFit[2]=(*p_k1)[1][0];
  fParFit[3]=(*p_k1)[0][0]*fParFit[1];
  fParFit[4]=(*p_k1)[3][0];
  fParFit[5]=(*p_k1)[2][0]*fParFit[1];
  fParFitErrors[1]=sqrt(pow((*p_k1)[0][0]*pow(fParFit[1],3),2)*pow((*NewCovMatrix_k1)[0][0],2) + pow((*p_k1)[2][0]*pow(fParFit[1],3),2)*pow((*NewCovMatrix_k1)[2][2],2));
  fParFitErrors[2]=(*NewCovMatrix_k1)[1][1];
  fParFitErrors[3]=sqrt(pow(fParFit[1],2)*pow((*NewCovMatrix_k1)[0][0],2) + pow((*p_k1)[0][0],2)*pow(fParFitErrors[1],2));
  fParFitErrors[4]=(*NewCovMatrix_k1)[3][3];
  fParFitErrors[5]=sqrt(pow(fParFit[1],2)*pow((*NewCovMatrix_k1)[2][2],2) + pow((*p_k1)[2][0],2)*pow(fParFitErrors[1],2));
  // p_start->Print();
  // printf("last *p\n");
  // p_k1->Print();
  // printf("@@@@@Cluster %d ended filtering@@@@@\n", this->GetClusterID());
  delete CovMatrix;
  delete CovMatrix_k1;
  delete NewCovMatrix_k1;
  delete CovMatrixI;
  delete CovMatrixI_k1;
  delete NewCovMatrixI_k1;
  delete SV;//State Vector
  delete V;
  delete VI;
  delete MS;//Multiple Scattering matrix
  delete p;//Parameters Vector
  delete p_start;//Parameters Vector
  delete p_k1;//Updated Parameters Vector
  delete f_k1;//Extrapolation Matrix
  delete pTilde_k1;//Extrapolated Parameters Vector
  delete D;//Trasformation matrix for covariance matrix
  delete Dt;
  delete H;
  delete dChi2;
  delete Ht;
  delete m;
  delete r;
  delete rt;
  if(isnan(chi2))chi2=999999.;
  fChi2=chi2;
  fndf=(2*this->GetNumPoints()-4);
  fReducedChi2=chi2/(2*this->GetNumPoints()-4);
  // return chi2;

}

//*********COMPATIBILITY AND MERGE METHODS*********//
double Cluster::KalmanToMerge(Cluster* trackm)
{
  // this->makeGraph(1);
  // trackl->fitCluster();
  // trackm->makeGraph(1);
  // trackm->fitCluster();
  double Chi2_kalman=9999.;
  TClonesArray* pointsl = this->GetPointsObj();
  // printf("this momentum = %f\n", CalculatePTot(this->GetPx(), this->GetPy(), this->GetPz()));
  // printf("this has %d points\n", this->GetNumPoints());
  TClonesArray* pointsm = trackm->GetPointsObj();
  // printf("trackm momentum = %f\n", CalculatePTot(trackm->GetPx(), trackm->GetPy(), trackm->GetPz()));
  // printf("trackm has %d points\n", trackm->GetNumPoints());
  Cluster *c0=new Cluster();
  for(int i=0; i<pointsl->GetEntries(); i++)
  {
    Datapoint* p= (Datapoint*)pointsl->At(i);
    c0->AddDataPoint(p);
  }
  for(int i=0; i<pointsm->GetEntries(); i++)
  {
    Datapoint* p= (Datapoint*)pointsm->At(i);
    c0->AddDataPoint(p);
  }
  // printf("c0 has %d points\n", c0->GetNumPoints());
  // c0->makeGraph(1);
  // c0->fitCluster();
  c0->SetPx(this->GetPx());
  c0->SetPy(this->GetPy());
  c0->SetPz(this->GetPz());

  c0->KalmanFilter();
  Chi2_kalman=c0->GetReducedChi2();
  return Chi2_kalman;

}

double Cluster::KalmanToMerge2(Cluster* trackm)
{
  this->KalmanFilter();
}

double Cluster::DataPointDistanceFromLine(Datapoint *p)
{
  // Distance between a point and a fit line
  // Used in isCompatible and isCompatible2 and in CheckNoisePoints
  double x = p->GetX();
  double y = p->GetY();
  double z = p->GetZ();

  XYZVector xp(x,y,z);
  XYZVector x0(fParFit[0], fParFit[2], fParFit[4]);
  //another point of line (for t = 1)
  XYZVector x1(fParFit[0] + fParFit[1], fParFit[2] + fParFit[3], fParFit[4] + fParFit[5] );

  XYZVector u = (x1-x0).Unit();
  double d = sqrt( ((xp-x0).Cross(u)).X()*((xp-x0).Cross(u)).X() + ((xp-x0).Cross(u)).Y()*((xp-x0).Cross(u)).Y() +
  ((xp-x0).Cross(u)).Z()*((xp-x0).Cross(u)).Z() );
  return d;

}

bool Cluster::CheckIfInsideCluster(TClonesArray *fnp, int i_start, int i_end)
{
  // Distance between extreme points of clusters
  // Used in isCompatible and isCompatible2

  bool IsInside = kTRUE;

  double NbPoints1 = fPoints->GetEntries();
  Datapoint *p1_start = (Datapoint*) fPoints->At(fExPt1);
  Datapoint *p1_end = (Datapoint*) fPoints->At(fExPt2);

  double ClusterLength = std::sqrt((p1_start->GetX()-p1_end->GetX()) * (p1_start->GetX()-p1_end->GetX())+
  (p1_start->GetY()-p1_end->GetY()) * (p1_start->GetY()-p1_end->GetY())+
  (p1_start->GetZ()-p1_end->GetZ()) * (p1_start->GetZ()-p1_end->GetZ()));

  double NbPoints2 = fnp->GetEntries();
  Datapoint *p2_start = (Datapoint*) fnp->At(i_start);

  double D1 = 0.;
  double D2 = 0.;

  D1 =  std::sqrt((p1_start->GetX()-p2_start->GetX()) * (p1_start->GetX()-p2_start->GetX())+
  (p1_start->GetY()-p2_start->GetY()) * (p1_start->GetY()-p2_start->GetY())+
  (p1_start->GetZ()-p2_start->GetZ()) * (p1_start->GetZ()-p2_start->GetZ()));

  D2 =  std::sqrt((p1_end->GetX()-p2_start->GetX()) * (p1_end->GetX()-p2_start->GetX())+
  (p1_end->GetY()-p2_start->GetY()) * (p1_end->GetY()-p2_start->GetY())+
  (p1_end->GetZ()-p2_start->GetZ()) * (p1_end->GetZ()-p2_start->GetZ()));

  if(D1 > ClusterLength | D2 > ClusterLength) IsInside = kFALSE;

  return IsInside;

}

bool Cluster::CheckIfInsideCluster(Datapoint *p)
{
  // Distance between extreme points of clusters
  // Used in isCompatible and isCompatible2

  bool IsInside = kTRUE;

  Datapoint *p1_start = (Datapoint*) fPoints->At(fExPt1);
  Datapoint *p1_end = (Datapoint*) fPoints->At(fExPt2);

  double ClusterLength = std::sqrt((p1_start->GetX()-p1_end->GetX()) * (p1_start->GetX()-p1_end->GetX())+
  (p1_start->GetY()-p1_end->GetY()) * (p1_start->GetY()-p1_end->GetY())+
  (p1_start->GetZ()-p1_end->GetZ()) * (p1_start->GetZ()-p1_end->GetZ()));

  double D1 = 0.;
  double D2 = 0.;

  D1 =  std::sqrt((p1_start->GetX()-p->GetX()) * (p1_start->GetX()-p->GetX())+
  (p1_start->GetY()-p->GetY()) * (p1_start->GetY()-p->GetY())+
  (p1_start->GetZ()-p->GetZ()) * (p1_start->GetZ()-p->GetZ()));

  D2 =  std::sqrt((p1_end->GetX()-p->GetX()) * (p1_end->GetX()-p->GetX())+
  (p1_end->GetY()-p->GetY()) * (p1_end->GetY()-p->GetY())+
  (p1_end->GetZ()-p->GetZ()) * (p1_end->GetZ()-p->GetZ()));

  if(D1 > ClusterLength | D2 > ClusterLength) IsInside = kFALSE;

  return IsInside;

}

double Cluster::ClusterPointsDistance(TClonesArray *fnp, int i_start, int i_end)
{
  // Distance between extreme points of clusters
  // Used in isCompatible and isCompatible2
  double NbPoints1 = fPoints->GetEntries();
  Datapoint *p1_start = (Datapoint*) fPoints->At(fExPt1);
  Datapoint *p1_end = (Datapoint*) fPoints->At(fExPt2);
  double NbPoints2 = fnp->GetEntries();
  Datapoint *p2_start = (Datapoint*) fnp->At(i_start);
  Datapoint *p2_end = (Datapoint*) fnp->At(i_end);
  double D_ss = 0.;
  double D_se = 0.;
  double D_es = 0.;
  double D_ee = 0.;

  D_ss = std::sqrt((p1_start->GetX()-p2_start->GetX()) * (p1_start->GetX()-p2_start->GetX())+
  (p1_start->GetY()-p2_start->GetY()) * (p1_start->GetY()-p2_start->GetY())+
  (p1_start->GetZ()-p2_start->GetZ()) * (p1_start->GetZ()-p2_start->GetZ()));

  D_se = std::sqrt((p1_start->GetX()-p2_end->GetX()) * (p1_start->GetX()-p2_end->GetX())+
  (p1_start->GetY()-p2_end->GetY()) * (p1_start->GetY()-p2_end->GetY())+
  (p1_start->GetZ()-p2_end->GetZ()) * (p1_start->GetZ()-p2_end->GetZ()));

  D_es =  std::sqrt((p1_end->GetX()-p2_start->GetX()) * (p1_end->GetX()-p2_start->GetX())+
  (p1_end->GetY()-p2_start->GetY()) * (p1_end->GetY()-p2_start->GetY())+
  (p1_end->GetZ()-p2_start->GetZ()) * (p1_end->GetZ()-p2_start->GetZ()));

  D_ee = std::sqrt((p1_end->GetX()-p2_end->GetX()) * (p1_end->GetX()-p2_end->GetX())+
  (p1_end->GetY()-p2_end->GetY()) * (p1_end->GetY()-p2_end->GetY())+
  (p1_end->GetZ()-p2_end->GetZ()) * (p1_end->GetZ()-p2_end->GetZ()));

  return std::min(std::min(D_ss, D_ee), std::min(D_se, D_es));
}

double Cluster::ClusterDistance(Cluster *cluster)
{
  TClonesArray *pts1 = (TClonesArray *)this->GetPointsObj();
  TClonesArray *pts2 = (TClonesArray *)cluster->GetPointsObj();

  double MinDistance = 1000.;
  double distance = 0.;
  for(int i1 = 0; i1<pts1->GetEntries(); i1++)
  {
    Datapoint *p1 = (Datapoint*)pts1->At(i1);
    for(int i2 = 0; i2<pts2->GetEntries(); i2++)
    {
      Datapoint *p2 = (Datapoint*)pts2->At(i2);
      distance=p1->distance3D(p2);
      if(distance < MinDistance) MinDistance = distance;
    }
  }
  return MinDistance;

}

double Cluster::ClusterDistance(Datapoint *p)
{
  TClonesArray *pts1 = (TClonesArray *)this->GetPointsObj();

  double MinDistance = 1000.;
  double distance = 0.;
  for(int i1 = 0; i1<pts1->GetEntries(); i1++)
  {
    Datapoint *p1 = (Datapoint*)pts1->At(i1);
    distance=p1->distance3D(p);
    if(distance < MinDistance) MinDistance = distance;
  }
  return MinDistance;

}

bool Cluster::isCompatible(Cluster *cluster/*, Double_t min_cos, Double_t min_dist_cl*/)
{
  // Used in PixelChamberEvent's method mergeCompatibleClusters
  // Check on the compatibility of two clusters using the following conditions:
  // scalar product between vectors made by director cosines of 2 fit lines has to be bigger than MIN_COS:
  // this means that 2 lines are parallels
  // the result of ClusterPointsDistance has to be smaller than 1:
  // this means that the 2 clusters have to be contiguos
  // this->makeGraph(1);
  // this->fitCluster();
  // cluster->makeGraph(1);
  // cluster->fitCluster();
  const double *ParFit = cluster->GetParFit();
  const double *ParFitErrors = cluster->GetParFitErrors();
  double cosX = fParFit[1];
  double cosY = fParFit[3];
  double cosZ = fParFit[5];
  double cosX1 = ParFit[1];
  double cosY1 = ParFit[3];
  double cosZ1 = ParFit[5];
  double EcosX=fParFitErrors[1];
  double EcosY=fParFitErrors[3];
  double EcosZ=fParFitErrors[5];
  double EcosX1=ParFitErrors[1];
  double EcosY1=ParFitErrors[3];
  double EcosZ1=ParFitErrors[5];

  double e1 = sqrt(EcosX*EcosX+EcosY*EcosY+EcosZ*EcosZ);
  double e2 = sqrt(EcosX1*EcosX1+EcosY1*EcosY1+EcosZ1*EcosZ1);

  double w;
  if(e1<e2) w = e1/e2;
  else
  w = e2/e1;

  XYZVector u(cosX, cosY, cosZ);
  XYZVector u1(cosX1, cosY1, cosZ1);
  double scalarproduct = u.Dot(u1);

  double ErrorScalarProduct = sqrt(pow(cosX1*EcosX,2)+pow(cosX*EcosX1,2)+pow(cosY1*EcosY,2)+pow(cosY*EcosY1,2)+pow(cosZ1*EcosZ,2)+pow(cosZ*EcosZ1,2));
  TClonesArray *fnp = (TClonesArray *) cluster->GetPointsObj();

  Datapoint *p1 = (Datapoint*) fnp->At(cluster->GetfExPt1());
  double d1 = this->DataPointDistanceFromLine(p1);
  Datapoint *p2 = (Datapoint*) fnp->At(cluster->GetfExPt2());
  double d2 = this->DataPointDistanceFromLine(p2);
  double d  = min(d1,d2);
  w=2.;
  //
  if(this->GetReducedChi2()<2.5
  && cluster->GetReducedChi2()<2.5
  && fabs(scalarproduct) >= min(MIN_COS,1-ErrorScalarProduct*w)
  && d <= MIN_DIST
  && this->ClusterPointsDistance(fnp, cluster->GetfExPt1(), cluster->GetfExPt2())<=MIN_DIST_CL
  && this->CheckIfInsideCluster(fnp, cluster->GetfExPt1(), cluster->GetfExPt2())==kFALSE){
    Cluster *c0=new Cluster();
    for(int i=0; i<fPoints->GetEntries();i++){
      Datapoint *p=(Datapoint *) fPoints->At(i);
      c0->AddDataPoint(p);
    }
    Cluster *c1=new Cluster();
    for(int i=0; i<fnp->GetEntries();i++){
      Datapoint *p=(Datapoint *) fnp->At(i);
      c1->AddDataPoint(p);
    }

    c0->merge(c1);
    c0->makeGraph(1);
    c0->fitCluster();
    double chi2max = (fPoints->GetEntries() * this->GetReducedChi2() + fnp->GetEntries() * cluster->GetReducedChi2())/
    (fPoints->GetEntries()+fnp->GetEntries());
    if(c0->GetReducedChi2()<1.5){
      return true;
    }
  } else if ( this->GetReducedChi2()<2.5
  && cluster->GetReducedChi2()<2.5
  && fabs(scalarproduct) >= min(MIN_COS,1-ErrorScalarProduct*w)
  && d <= MIN_DIST
  && this->CheckIfInsideCluster(fnp, cluster->GetfExPt1(), cluster->GetfExPt2())==kTRUE){
    //
    Cluster *c0=new Cluster();
    for(int i=0; i<fPoints->GetEntries();i++){
      Datapoint *p=(Datapoint *) fPoints->At(i);
      c0->AddDataPoint(p);
    }
    Cluster *c1=new Cluster();
    for(int i=0; i<fnp->GetEntries();i++){
      Datapoint *p=(Datapoint *) fnp->At(i);
      c1->AddDataPoint(p);
    }

    c0->merge(c1);
    c0->makeGraph(1);
    c0->fitCluster();
    double chi2max = (fPoints->GetEntries() * this->GetReducedChi2() + fnp->GetEntries() * cluster->GetReducedChi2())/
    (fPoints->GetEntries()+fnp->GetEntries());
    if(c0->GetReducedChi2()<1.5){
      return true;
    }
  }

  return false;
}

bool Cluster::isCompatibleKalman(Cluster *cluster/*, Double_t min_cos, Double_t min_dist_cl*/)
{
  // Used in PixelChamberEvent's method mergeCompatibleClusters
  // Check on the compatibility of two clusters using the following conditions:
  // scalar product between vectors made by director cosines of 2 fit lines has to be bigger than MIN_COS:
  // this means that 2 lines are parallels
  // the result of ClusterPointsDistance has to be smaller than 1:
  // this means that the 2 clusters have to be contiguos
  this->KalmanFilter();
  cluster->KalmanFilter();
  const double *ParFit = cluster->GetParFit();
  const double *ParFitErrors = cluster->GetParFitErrors();
  double cosX = fParFit[1];
  double cosY = fParFit[3];
  double cosZ = fParFit[5];
  double cosX1 = ParFit[1];
  double cosY1 = ParFit[3];
  double cosZ1 = ParFit[5];
  double EcosX=fParFitErrors[1];
  double EcosY=fParFitErrors[3];
  double EcosZ=fParFitErrors[5];
  double EcosX1=ParFitErrors[1];
  double EcosY1=ParFitErrors[3];
  double EcosZ1=ParFitErrors[5];

  double e1 = sqrt(EcosX*EcosX+EcosY*EcosY+EcosZ*EcosZ);
  double e2 = sqrt(EcosX1*EcosX1+EcosY1*EcosY1+EcosZ1*EcosZ1);

  double w;
  if(e1<e2) w = e1/e2;
  else
  w = e2/e1;

  XYZVector u(cosX, cosY, cosZ);
  XYZVector u1(cosX1, cosY1, cosZ1);
  double scalarproduct = u.Dot(u1);

  double ErrorScalarProduct = sqrt(pow(cosX1*EcosX,2)+pow(cosX*EcosX1,2)+pow(cosY1*EcosY,2)+pow(cosY*EcosY1,2)+pow(cosZ1*EcosZ,2)+pow(cosZ*EcosZ1,2));
  TClonesArray *fnp = (TClonesArray *) cluster->GetPointsObj();

  Datapoint *p1 = (Datapoint*) fnp->At(cluster->GetfExPt1());
  double d1 = this->DataPointDistanceFromLine(p1);
  Datapoint *p2 = (Datapoint*) fnp->At(cluster->GetfExPt2());
  double d2 = this->DataPointDistanceFromLine(p2);
  double d  = min(d1,d2);
  w=2.;
  //
  if(this->GetReducedChi2()<2.5
  && cluster->GetReducedChi2()<2.5
  && fabs(scalarproduct) >= min(MIN_COS,1-ErrorScalarProduct*w)
  && d <= MIN_DIST
  && this->ClusterPointsDistance(fnp, cluster->GetfExPt1(), cluster->GetfExPt2())<=MIN_DIST_CL
  && this->CheckIfInsideCluster(fnp, cluster->GetfExPt1(), cluster->GetfExPt2())==kFALSE)
  {

    double chiKalman=this->KalmanToMerge(cluster);
    if(chiKalman<2.5){
      return true;
    }
  } else if ( this->GetReducedChi2()<2.5
  && cluster->GetReducedChi2()<2.5
  && fabs(scalarproduct) >= min(MIN_COS,1-ErrorScalarProduct*w)
  && d <= MIN_DIST
  && this->CheckIfInsideCluster(fnp, cluster->GetfExPt1(), cluster->GetfExPt2())==kTRUE)
  {

    double chiKalman=this->KalmanToMerge(cluster);
    if(chiKalman<2.5){
      return true;
    }
  }

  return false;
}


bool Cluster::isCompatible2(Cluster *cluster)
{
  // Used in PixelChamberEvent's method mergeCompatibleClusters2
  // Check on the compatibility of points and fit lines:
  // the distance between datapoints and fit lines has to be smaller than 0.07
  // this method is used to check the fraction of compatible clusters
  bool Compatibility = kFALSE;
  if(this->GetReducedChi2()<2.5 && cluster->GetReducedChi2()<2.5){
    double FracCompatible = 0.;
    TClonesArray *fnp1 = (TClonesArray *)this->GetPointsObj();
    TClonesArray *fnp2 = (TClonesArray *)cluster->GetPointsObj();
    if(fnp1->GetEntries() >fnp2->GetEntries()){
      Datapoint *p;
      double d = 0.;
      int nPCompatible = 0;
      for(int m = 0; m < fnp2->GetEntries(); m++) {
        p = (Datapoint*) fnp2->At(m);
        d = this->DataPointDistanceFromLine(p);
        if(d < MIN_DIST) nPCompatible++;
      }
      if(FracCompatible > FRAC_COMPATIBLE) Compatibility = kTRUE;
      return Compatibility;
    } else if(fnp2->GetEntries() >fnp1->GetEntries()){
      Datapoint *p;
      double d = 0.;
      int nPCompatible = 0;
      for(int m = 0; m < fnp1->GetEntries(); m++) {
        p = (Datapoint*) fnp1->At(m);
        d = cluster->DataPointDistanceFromLine(p);
        if(d < MIN_DIST) nPCompatible++;
      }
      if(FracCompatible > FRAC_COMPATIBLE) Compatibility = kTRUE;
      return Compatibility;
    }
  } else
  return Compatibility;
}

void Cluster::merge(Cluster *cluster)
{
  //
  TClonesArray *pts = (TClonesArray *)cluster->GetPointsObj();
  Datapoint *p;
  for( int l = 0; l<pts->GetEntries(); l++)
  {
    p = (Datapoint*) pts->At(l);
    this->AddDataPoint(p);
  }
}

//*********PRIVATE METHOD*********//

void Cluster::line(double t, const double *p, double &x, double &y, double &z)
{
  x = p[0] + p[1]*t;
  y = p[2] + p[3]*t;
  z = p[4] + p[5]*t;
}

double Cluster::Chi2TrackToVertex(double Vx, double Vy, double Vz, TMatrixD *CovVert){

  // printf("vertex covariance matrix in cluster method\n");
  // CovVert->Print();

  double alpha = fParFit[3]/fParFit[1];
  double y0    = fParFit[2];
  double beta  = fParFit[5]/fParFit[1];
  double z0    = fParFit[4];

  double DeltaX = Vx - fParFit[0];

  TMatrixDSym *VComb = new TMatrixDSym(2);

  //   with dr = { y0 + alpha * dX - Vy,  z0 + betha * dX - Vz}, dX = Vx-x0
  // and VComb ={
  // sig^2_Vy + V11+ dX^2*V00 + alpha^2*sig^2_Vx + 2*dX*V10,
  // sig_Vy_Vz + dX^2*V13 + 2*alpha*betha*sig^2_vx
  //    sym.matrix                                                                            ,
  // sig^2_Vz + V33 + dX^2*V22 + betha^2*sig^2_Vx + 2*dX*V23
  // }


  (*VComb)[0][0] = (*CovVert)[1][1] +
  (*fCovarianceMatrix)[1][1] +
  pow(DeltaX,2)* (*fCovarianceMatrix)[0][0] +
  pow(alpha,2)*(*CovVert)[0][0] +
  2*DeltaX*(*fCovarianceMatrix)[1][0] -
  2.*alpha*(*CovVert)[0][1];

  (*VComb)[0][1] = (*CovVert)[1][2] +
  pow(DeltaX,2)*(*fCovarianceMatrix)[0][2] +
  alpha*beta*(*CovVert)[0][0] +
  DeltaX * (*fCovarianceMatrix)[0][3] +
  DeltaX * (*fCovarianceMatrix)[1][2] +
  (*fCovarianceMatrix)[1][3] -
  alpha * (*CovVert)[0][2] -
  beta * (*CovVert)[0][1];

  (*VComb)[1][0] = (*VComb)[0][1];

  (*VComb)[1][1] = (*CovVert)[2][2] +
  (*fCovarianceMatrix)[3][3] +
  pow(DeltaX,2)*(*fCovarianceMatrix)[2][2]+
  pow(beta,2)*(*CovVert)[0][0]+
  2.*DeltaX*(*fCovarianceMatrix)[2][3]-
  2.*beta*(*CovVert)[0][2];

  double qh[2] = {y0 + alpha * DeltaX - Vy,z0 + beta * DeltaX - Vz};

  TDecompChol *Decomp = new TDecompChol(*VComb);
  Decomp->Invert(*VComb);

  TMatrixD *QH = (TMatrixD*) new TMatrixD(2,1);
  for(int i=0; i<2; i++)(*QH)[i][0]=qh[i];
  TMatrixD *tQH = (TMatrixD*) new TMatrixD(1,2);
  tQH->Transpose(*QH);
  TMatrixD *mchi = (TMatrixD*) new TMatrixD(1,1);
  *mchi = (*tQH) * (*VComb) * (*QH);
  double chi = (*mchi)[0][0];
  // printf("chi=%f\n",chi);

  return chi;

}
