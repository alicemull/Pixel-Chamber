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
#include "cluster.h"

using namespace ROOT::Math;

// #define NB_MIN				2
// #define NB_MAX				4
#define EPS				1
#define MIN_TRACK_POINTS	        10
#define MAX_CHI2			2.5
#define MIN_COS				0.9998
//#define MIN_COS				0.9993
#define MIN_DIST			0.05

Datapoint::Datapoint(double x, double y, double z, int i, int j, int k)
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

Datapoint::Datapoint()
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

Datapoint::Datapoint(Datapoint &p) : fx(p.fx),
				     fy(p.fy),
				     fz(p.fz),
				     fi(p.fi),
				     fj(p.fj),
				     fk(p.fk),
				     fused(p.fused),
				     fneighbours(p.fneighbours){}
 

int Datapoint::distance(Datapoint *p2)
{
	//
	int i = abs(fi - p2->getI());
	int j = abs(fj - p2->getJ());
	int k = abs(fk - p2->getK());
	//
	return std::max(std::max(i, j), k);
}

bool Datapoint::isNoise(int NB_MIN, int NB_MAX){
return fneighbours < NB_MIN || fneighbours > NB_MAX;
}

bool Datapoint::isValid(int NB_MIN, int NB_MAX){
  //  printf("inside isvalid: fused = %d fneighbours = %d isNoise(NB_MIN,NB_MAX) = %d\n",fused,fneighbours,isNoise(NB_MIN,NB_MAX));
  return !fused && !isNoise(NB_MIN,NB_MAX);
}


//************************************************************************//

Cluster::Cluster()
{
	fClusterID = 0;
	// fPoints = new TObjArray(100);	
	// fPoints->SetOwner(kTRUE);
	fPoints = new TClonesArray("Datapoint",2500);	
	fGraph = new TGraph2DErrors();
	fndf = 1.;
	fChi2 = -1.;
	fReducedChi2 = -1.;
	fParFit = new double(6);
	fParFitErrors = new double(6);

	for(int i=0; i<6; i++){
	  fParFit[i] = 0.;
	  fParFitErrors[i] = 0.;
	}

	int n = 10000;
	
	fPline3D = new TPolyLine3D(n);
	flXZ = new TPolyLine(n);
	flXY = new TPolyLine(n);
	flYZ = new TPolyLine(n);
	
	ft0=0.;
	ftf=0.;
}

Cluster::Cluster(int ClusterID)
{
	fClusterID = ClusterID;
	// fPoints = new TObjArray(100);
	// fPoints->SetOwner(kTRUE);
	fPoints = new TClonesArray("Datapoint",2500);	
	fGraph = new TGraph2DErrors();
	fndf = 1.;
	fChi2 = -1.;
	fReducedChi2 = -1.;
	fParFit = new double[6];
	fParFitErrors = new double[6];

	for(int i=0; i<6; i++){
	  fParFit[i] = 0.;
	  fParFitErrors[i] = 0.;
	}
	int n = 10000;
	
	fPline3D = new TPolyLine3D(n);
	flXZ = new TPolyLine(n);
	flXY = new TPolyLine(n);
	flYZ = new TPolyLine(n);

	ft0=0.;
	ftf=0.;
	//	printf("allocating cluster n. %d\n",fClusterID);
}
	      

Cluster::Cluster(Cluster &c) : fClusterID(c.fClusterID),
			       fPoints(c.fPoints),
			       fPline3D(c.fPline3D),
			       fGraph(c.fGraph),
			       fndf(c.fndf),
			       fChi2(c.fChi2),
			       fReducedChi2(c.fReducedChi2),
			       fParFit(c.fParFit),
			       fParFitErrors(c.fParFitErrors),
			       flXZ(c.flXZ),
			       flXY(c.flXY),
			       flYZ(c.flYZ),
			       ft0(c.ft0),
			       ftf(c.ftf){}

Cluster::~Cluster()
{
  //  printf("deallocating cluster n. %d\n",fClusterID);

	//
  if(fPoints) fPoints->Clear("C");
  if(fPline3D) delete fPline3D;
  if(fGraph) delete fGraph;
  if(fParFit) delete fParFit;
  if(fParFitErrors) delete fParFitErrors;
  if(flXY) delete flXY;
  if(flXZ) delete flXZ;
  if(flYZ) delete flYZ;
}

void Cluster::AddDataPoint(Datapoint *p){ 

  TClonesArray &arrTemp  = *fPoints; 
  new (arrTemp[arrTemp.GetEntriesFast()]) Datapoint(*p);

}

bool Cluster::isNeighbour(Datapoint *p, int NB_MIN, int NB_MAX)
{
	//
  if (!p->isValid(NB_MIN,NB_MAX))
		return false;
	//
	for (int n = 0; n < fPoints->GetEntries(); n++)
	{
		//
		Datapoint *p0 = (Datapoint*)(fPoints->At(n));
		//
		if (p0->distance(p) <= EPS)
			return true;
	}
	//
	return false;
}

void Cluster::makeGraph()
{
	int ip = 0;
	double ex=0.015, ey=0.025, ez=0.015;
	//
	if (fGraph) delete fGraph;
	// //
	fGraph = new TGraph2DErrors();
	// printf("fPoints->GetEntries()=%d\n",fPoints->GetEntries());
    for(int m = 0; m < fPoints->GetEntries(); m++)
	{
		Datapoint *pl = (Datapoint*)(fPoints->At(m));
		fGraph->SetPoint(ip, pl->getX(), pl->getY(), pl->getZ());
		fGraph->SetPointError(ip, ex, ey, ez);
		ip++;
	}
}

void Cluster::line(double t, const double *p, double &x, double &y, double &z)
{
	x = p[0] + p[1]*t;
    y = p[2] + p[3]*t;
    z = p[4] + p[5]*t;
}

void Cluster::fitCluster()
{
  double p0[6] = {-1,-4,-0.01, -0.05,-0.1, -0.5}; //line parameters
   	ROOT::Fit::Fitter  fitter;
 
 
    // make the functor object
   	SumDistance2 sdist(fGraph);
   	ROOT::Math::Functor fcn(sdist,6); // set the function and the initial parameter values
   	double * xg = fGraph->GetX();
   	double * yg = fGraph->GetY();
   	double * zg = fGraph->GetZ();
   	int np = fGraph->GetN();
	double minX=999., minY=999., minZ=999.;
	double maxX=-999., maxY=-999., maxZ=-999.;

	for(int i=0; i<np;i++) {
	  if(xg[i] < minX) minX = xg[i];
	  if(yg[i] < minY) minY = yg[i];
	  if(zg[i] < minZ) minZ = zg[i];
	  if(xg[i] > maxX) maxX = xg[i];
	  if(yg[i] > maxY) maxY = yg[i];
	  if(zg[i] > maxZ) maxZ = zg[i];
	  //	  printf("graph x=%f y=%f z=%f\n",xg[i],yg[i],zg[i]);
	}

	double ux = xg[np-1] - xg[0];
	double uy = yg[np-1] - yg[0];
	double uz = zg[np-1] - zg[0];
	// double ux = maxX - minX;
	// double uy = maxY - minY;
	// double uz = maxZ - minZ;
	double norm = sqrt(ux*ux+uy*uy+uz*uz);
	ux /= norm;
	uy /= norm;
	uz /= norm;

	// double pStart[6] = {xg[0],-6.3,yg[0], -0.08,zg[0], -0.4};
	printf("fit starting point x=%f, y=%f z=%f\n",minX,minY,minZ);

	double pStart[6]; 
	if(fParFit[1] == 0. & fParFit[3] == 0. & fParFit[5] == 0.) {
	  pStart[0] = xg[np-1];
	  pStart[1] = ux;
	  pStart[2] = yg[np-1];
	  pStart[3] = uy;
	  pStart[4] = zg[np-1];
	  pStart[5] = uz;

	} else {
	  printf("second round\n");
	  pStart[0] = xg[np-1];
	  pStart[1] = fParFit[1];
	  pStart[2] = yg[np-1];
	  pStart[3] = fParFit[3];
	  pStart[4] = zg[np-1];
	  pStart[5] = fParFit[5];
	}

	printf("np = %d\n",np);
   	fitter.SetFCN(fcn,pStart,np,1);
   	for (int i = 0; i < 6; ++i) {
	  fitter.Config().ParSettings(i).SetStepSize(0.01);
	} 
	//	fitter.Config().ParSettings(0).SetLimits(-16.,16.);

	bool ok = fitter.FitFCN();

	const ROOT::Fit::FitResult & result = fitter.Result();
 
    //std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
    //result.Print(std::cout);

    double chi2 = result.Chi2();
    double ndf = result.Ndf();
    double reducedChi2 = chi2/ndf;
    
    //std::cout << "reduced chi2 = " << reducedChi2 << std::endl;
    norm = 1.;

    fChi2 = chi2*norm;
    fndf = ndf*norm;
    fReducedChi2 = reducedChi2*norm;
   
    const double *ParFit;
    ParFit = result.GetParams();
 
    const double *ParFitErrors;
    ParFitErrors = result.GetErrors();

    for(int i=0;i<6; i++) {
      fParFit[i] = norm*ParFit[i];
      fParFitErrors[i] = norm*ParFitErrors[i];
    } 

    printf("fit x0=%f y0=%f z0=%f\n",fParFit[0],fParFit[2],fParFit[4]);
    printf("fit cosX=%f cosY=%f cosZ=%f\n",fParFit[1],fParFit[3],fParFit[5]);
    printf("fit chi2=%f ndf=%d\n",fChi2,fndf);

    // draw the fitted line
    //    double t0 = 0;
    //    double dt = 10;
    int n = 10000;

    if (fPline3D) delete fPline3D;
    if (flXY) delete flXY;
    if (flXZ) delete flXZ;
    if (flYZ) delete flYZ;
    fPline3D = new TPolyLine3D(n);
    flXZ = new TPolyLine(n);
    flXY = new TPolyLine(n);
    flYZ = new TPolyLine(n);

    printf("Poly line initial coordinates x=%f y=%f z=%f\n",xg[0],yg[0],zg[0]);
    printf("Poly line final coordinates x=%f y=%f z=%f\n",xg[np-1],yg[np-1],zg[np-1]);

    double ParLine[6]={minX,fParFit[1],minY,fParFit[3],minZ,fParFit[5]};
    for (int i = 0; i <np;++i) 
    {
      	// double t0=(xg[0]-fParFit[0])/fParFit[1];
      	// double tf=(xg[np-1]-fParFit[0])/fParFit[1];
     	ft0=(minX-fParFit[0])/fParFit[1];
      	ftf=(maxX-fParFit[0])/fParFit[1];
        //       double t = t0+ dt*i/np;
        double t = ft0 + (ftf-ft0)*i/(np-1.);
        double x,y,z;
	// ParLine[0] = minX;
	// ParLine[2] = minY;
	// ParLine[4] = minZ;
	// ParLine[1] = fParFit[1];
	// ParLine[3] = fParFit[3];
	// ParLine[5] = fParFit[5];
        line(t,fParFit,x,y,z);
        //       printf("x0=%f y0=%f z0=%f\n",parFit[0],parFit[2],parFit[4]);
        //       printf("x=%f y=%f z=%f\n",x,y,z);
        fPline3D->SetPoint(i,x,y,z);
        flXZ->SetPoint(i,x,z);
        flXY->SetPoint(i,x,y);
        flYZ->SetPoint(i,y,z);
    }
}

double Cluster::pLineDistance(const double *ParFit)
{
	//
	XYZVector u(fParFit[1], fParFit[3], fParFit[5]);
	XYZVector u1(ParFit[1], ParFit[3], ParFit[5]);
	
	XYZVector P(fParFit[0], fParFit[2], fParFit[4]);
	XYZVector P1(ParFit[0], ParFit[2], ParFit[4]);
	
	XYZVector P_diff = P-P1;

	XYZVector e = u.Cross(u1);
	XYZVector e_norm = e.Unit();

	return P_diff.Dot(e_norm);

}

bool Cluster::isCompatible(Cluster *cluster)
{
	//
	const double *ParFit = cluster->getParFit();
	const double *ParFitErrors = cluster->getParFitErrors();
	double norm = sqrt(fParFit[1]*fParFit[1]+fParFit[3]*fParFit[3]+fParFit[5]*fParFit[5]);
	double norm1 = sqrt(ParFit[1]*ParFit[1]+ParFit[3]*ParFit[3]+ParFit[5]*ParFit[5]);
	double cosX = fParFit[1]/norm;
	double cosY = fParFit[3]/norm;
	double cosZ = fParFit[5]/norm;
	double cosX1 = ParFit[1]/norm1;
	double cosY1 = ParFit[3]/norm1;
	double cosZ1 = ParFit[5]/norm1;
	// double cosX = fParFit[1];
	// double cosY = fParFit[3];
	// double cosZ = fParFit[5];
	// double cosX1 = ParFit[1];
	// double cosY1 = ParFit[3];
	// double cosZ1 = ParFit[5];
	XYZVector u(cosX, cosY, cosZ);
	XYZVector u1(cosX1, cosY1, cosZ1);
	double scalarproduct = u.Dot(u1);
	printf("scalar product = %f\n",scalarproduct);

	double diffCosX = fabs(cosX - cosX1);
	double diffCosY = fabs(cosY - cosY1);
	double diffCosZ = fabs(cosZ - cosZ1);
	printf("First Cluster cosX=%f cosY=%f cosZ=%f\n",cosX,cosY,cosZ);
	printf("Second Cluster cosX=%f cosY=%f cosZ=%f\n",cosX1,cosY1,cosZ1);

	// if( diffCosX < MIN_COS && diffCosY < MIN_COS && 
	// 	diffCosZ < MIN_COS && pLineDistance(ParFit) <= MIN_DIST)
		if( fabs(scalarproduct) > MIN_COS && pLineDistance(ParFit) <= MIN_DIST)
			return true;
	//
	return false;
}

/*if( fabs(fParFit[1]-ParFit[1]) < MIN_COS && fabs(fParFit[3]-ParFit[3]) < MIN_COS && 
		fabs(fParFit[5]-ParFit[5]) < MIN_COS && pLineDistance(ParFit) < 0.05)
	*/

void Cluster::merge(Cluster *cluster)
{
	//
	TClonesArray *pts = cluster->getPointsObj();
	//	pts->SetOwner();
	//
	Datapoint *p;
	for( int l = 0; l<pts->GetEntries(); l++){
	  p = (Datapoint*) pts->At(l);
	  this->AddDataPoint(p);
	}
}

void Cluster::ViewClusterFit()
{
	TH2F *h_XZ=new TH2F("h_coord_XZ","Coordinates of Pixel Centre", 512*2 , -7.68*2, 7.68*2, 512 , -7.68, 7.68);
	TH2F *h_XY=new TH2F("h_coord_XY","Coordinates of Pixel Centre", 512*2 , -7.68*2, 7.68*2, 10 , -0.25, 0.25);
	TH2F *h_YZ=new TH2F("h_coord_YZ","Coordinates of Pixel Centre", 10 , -0.25, 0.25, 512 , -7.68, 7.68);
	for (int i= 0; i<fPoints->GetEntries(); i++)
	{
		Datapoint *p = (Datapoint*)(fPoints->At(i));
		double x = p->getX();
		double y = p->getY();
		double z = p->getZ();
		h_XZ->Fill(x,z);
		h_XY->Fill(x,y);
		h_YZ->Fill(y,z);
	}
	TCanvas *c1 =new TCanvas ("c1","ZX projection",1000,1200);
      c1->cd();
      h_XZ->Draw("colz");
      flXZ->Draw("same");
      
    TCanvas *c2 =new TCanvas ("c2","XY projection",1000,1200);
      c2->cd();
      h_XY->Draw("colz");
      flXY->Draw("same");

    TCanvas *c3 =new TCanvas ("c3","YZ projection",1000,1200);
      c3->cd();
      h_YZ->Draw("colz");
      flYZ->Draw("same");

}

void Cluster::RemoveDataPointAndCompress(int ip){
  fPoints->RemoveAt(ip);
  fPoints->Compress();
}

double Cluster::DataPointDistanceFromLine(Datapoint *p){

  double x = p->getX();
  double y = p->getY();
  double z = p->getZ();

  XYZVector xp(x,y,z);
  XYZVector x0(fParFit[0], fParFit[2], fParFit[4]);
  XYZVector x1(fParFit[0] + fParFit[1], fParFit[2] + fParFit[3], fParFit[4] + fParFit[5] ); //another point of line (for t = 1)

  XYZVector u = (x1-x0).Unit();
  double d = sqrt( ((xp-x0).Cross(u)).X()*((xp-x0).Cross(u)).X() + ((xp-x0).Cross(u)).Y()*((xp-x0).Cross(u)).Y() +
    ((xp-x0).Cross(u)).Z()*((xp-x0).Cross(u)).Z() );
   return d;

}

double * Cluster::getStartLine(){
 
  double *xyz = new double[3]; 
  line(ft0,fParFit,xyz[0],xyz[1],xyz[2]);
  return xyz;
}

double * Cluster::getEndLine(){
 
  double *xyz = new double[3]; 
  line(ftf,fParFit,xyz[0],xyz[1],xyz[2]);
  return xyz;
}


//************************************************************************//
SumDistance2::SumDistance2(TGraph2DErrors *graph)
{
	fGraph = graph;
}

//SumDistance2::SumDistance2(SumDistance2 &s) : fGraph(s.fGraph){}

double SumDistance2::distance2(double x,double y,double z, double ex, double ey, double ez, const double *p)
{
	XYZVector xp(x,y,z);
    XYZVector x0(p[0], p[2], p[4]);
    XYZVector x1(p[0] + p[1], p[2] + p[3], p[4] + p[5] ); //another point of line (for t = 1)
    XYZVector u = (x1-x0).Unit();
    double d2 = ((xp-x0).Cross(u)).X()*((xp-x0).Cross(u)).X() / ex/ex + ((xp-x0).Cross(u)).Y()*((xp-x0).Cross(u)).Y() / ey/ey +
 ((xp-x0).Cross(u)).Z()*((xp-x0).Cross(u)).Z() / ez/ez;
   return d2;
}

double SumDistance2::operator() (const double *par)
{
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
   /* if (first) 
    {
        std::cout << "Total Initial distance square = " << sum << std::endl;
    }*/
    first = false;
       //       std::cout << " distance square = " << sum << std::endl;
    return sum;
}

//************************************************************************//

PixelChamberEvent::PixelChamberEvent()
{
	// fPoints = new TObjArray(2000);
	// fPoints->SetOwner(kTRUE);
	fPoints = new TClonesArray("Datapoint",2500);	
	fEventID = 0;
	// fTracklets = new TObjArray(100);
	// fTracklets->SetOwner(kTRUE);
	fTracklets = new TClonesArray("Cluster",2500);	
	fNoise = new Cluster(0);
	fClustersNumber = 0;
}

PixelChamberEvent::PixelChamberEvent(int EventID)
{
	// fPoints = new TObjArray(100);
	// fPoints->SetOwner(kTRUE);
	fPoints = new TClonesArray("Datapoint",2500);	
	fEventID = EventID;
	// fTracklets = new TObjArray(100);
	// fTracklets->SetOwner(kTRUE);
	fTracklets = new TClonesArray("Cluster",2500);	
	fNoise = new Cluster(0);
	fClustersNumber = 0;
	//	printf("allocating pixel chamber event n. %d\n",fEventID);

}

PixelChamberEvent::PixelChamberEvent(PixelChamberEvent &e) : fPoints(e.fPoints),
							     fEventID(e.fEventID),
							     fTracklets(e.fTracklets),
							     fNoise(e.fNoise),
							     fClustersNumber(e.fClustersNumber){}

PixelChamberEvent::~PixelChamberEvent()
{
  //  printf("deallocating pixel chamber event n. %d\n",fEventID);

  // if(fPoints) delete fPoints;
  // if(fTracklets) delete fTracklets;
  // if(fNoise) delete fNoise;
  fPoints->Clear("C");
  fTracklets->Clear("C");
  delete fNoise;

}

void PixelChamberEvent::AddDataPoint(Datapoint *p){ 
  //  printf("adding data point\n");
  TClonesArray &arrTemp  = *fPoints; 
  new (arrTemp[arrTemp.GetEntriesFast()]) Datapoint(*p);

}

void PixelChamberEvent::AddTracklet(Cluster *c){ 
  //  printf("adding tracklet\n");
  TClonesArray &arrTemp  = *fTracklets; 
  new (arrTemp[arrTemp.GetEntriesFast()]) Cluster(*c);

}

int PixelChamberEvent::findClusters(int NB_MIN, int NB_MAX)
{
	// count neighbours
  printf("fPoints->GetEntries()=%d\n",fPoints->GetEntries());
  Datapoint *pl, *pm, *pn;
  for (int l = 0; l < fPoints->GetEntries(); l++)
    {
      for (int m = 0; m < fPoints->GetEntries(); m++)
  	{
  	  pl = (Datapoint*)(fPoints->At(l));
  	  pm = (Datapoint*)(fPoints->At(m));				
  	  if (l != m && (pl->distance(pm) <= EPS)){
  	    pl->incNeighbours();
	    //printf("increase neighbours\n");
	  }
  	}
    }
  
  // find tracks
  printf("now find tracks\n");
  Cluster *track;

  TClonesArray *pts;
  int tid = 0;
  // printf("fPoints address = %p\n",fPoints);
  // printf("fTracklets address = %p\n",fTracklets);
  // printf("event address = %p\n",this);

  for (int l = 0; l < fPoints->GetEntries(); l++)
    {
      // skip used point or noise
      pl = (Datapoint*)(fPoints->At(l));
      if (pl->isValid(NB_MIN,NB_MAX))
  	{
  	  // new track
	  //	  printf("new track\n");
  	  tid++;
  	  track = new Cluster( tid );
	  //	  printf("track address = %p\n",track);
  	  track->AddDataPoint(pl);
  	  // scans all points from l+1
  	  int m = l + 1;
  	  //
  	  while(m < fPoints->GetEntries())
  	    {
  	      // check if point m is neighbour of current track
  	      pm = (Datapoint*)(fPoints->At(m));				
	      //printf("check if point is valid = %d\n",pm->isValid(NB_MIN,NB_MAX));
	      //printf("check if point is neighbour = %d\n",track->isNeighbour(pm,NB_MIN,NB_MAX));
  	      if (track->isNeighbour(pm,NB_MIN,NB_MAX))
  		{
  		  // add point
		  //		  printf("add point to track cluster\n");
  		  track->AddDataPoint(pm);
  		  // mark point as used
  		  pm->setUsed(true);
  		  // reset loop
  		  m = l + 1;
  		}
  	      else
  		// next point
  		m++;
  	    }
  	  //
  	  pts = track->getPointsObj();
	  //	  printf("pts->GetEntries()=%d\n",pts->GetEntries());
  	  if (pts->GetEntries() >= MIN_TRACK_POINTS){
	    //  	    printf("adding tracklet\n");
  	    this->AddTracklet(track);
  	  } else
  	    {
	      //  	      printf("adding to noise cluster\n");
  	      for (int n = 0; n < pts->GetEntries(); n++){
  		pn = (Datapoint*)(pts->At(n));
  		fNoise->AddDataPoint(pn);
  	      }	
  	    }
	  
  	}
    }
  
  // noise
  for (int l = 0; l < fPoints->GetEntries(); l++)
    {
      pn = (Datapoint*)(fPoints->At(l));
      if (pn->isNoise(NB_MIN,NB_MAX))
  	fNoise->AddDataPoint(pn);
    }
  
  //
  printf("fTracklets->GetEntries()=%d\n",fTracklets->GetEntries());
  return fTracklets->GetEntries();
}

void PixelChamberEvent::fitClusters(bool kReject)
{
	//
	int l = 0;
	const double *parfit;
	Cluster *track[fTracklets->GetEntries()];
	printf("sono qui\n");
	for(int i=0;i<fTracklets->GetEntries();i++) track[i] = new Cluster(i);

	printf("number of tracklets = %d\n", fTracklets->GetEntries());
	int nentries = fTracklets->GetEntries();
	while (l < fTracklets->GetEntries())
	{
	  //	  track[l] = (Cluster*)(fTracklets->At(l));
	  printf("l=%d\n",l);
	  track[l] = (Cluster*)(fTracklets->At(l));
		TClonesArray *pts = track[l]->getPointsObj();
		//		pts->SetOwner();
		printf("tracklet n. %d has %d points\n", l,pts->GetEntries());
		// make graph
		track[l]->makeGraph();
		printf("created graph of tracklet = %d\n",l);
		// fit
		track[l]->fitCluster();	
		// fTracklets->RemoveAt(l);
		// fTracklets->AddAt(track[l],l);

		parfit = track[l]->getParFit();
		printf("@@@@@@@@@@@@ cluster %d parfit[1]=%f parfit[3]=%f partfit[5]=%f\n",l,parfit[1],parfit[3],parfit[5]);
		// check chi2
		if (track[l]->getReducedChi2()>MAX_CHI2 && kReject)
		{
			// moves back bad track to noise
			fNoise->merge(track[l]);
			//
			fTracklets->Remove(track[l]);
			fTracklets->Compress();
		}
		else
			l++;
	}
	
	//  parfit = track[0]->getParFit();
	//  printf("###### parfit[1]=%f parfit[3]=%f partfit[5]=%f\n",parfit[1],parfit[3],parfit[5]);
	// fTracklets->Clear();
	// for(l=0;l<nentries;l++) fTracklets->Add(track[l]);

	// Cluster *t = (Cluster*)(fTracklets->At(1));
	// parfit = t->getParFit();
	// printf("parfit[1]=%f parfit[3]=%f partfit[5]=%f\n",parfit[1],parfit[3],parfit[5]);
}


void PixelChamberEvent::mergeCompatibleClusters()
{
	//
  printf("find compatibility of %d tracklets\n",fTracklets->GetEntries());
  for (int l = 0; l < fTracklets->GetEntries(); l++) {
    Cluster *trackl = (Cluster*)(fTracklets->At(l));
    printf("cluster n. %d has trackid = %d\n",l,trackl->getClusterID());
  }
	  
	for (int l = 0; l < fTracklets->GetEntries(); l++)
	{
		int m = l + 1;
		while (m < fTracklets->GetEntries())
		{
			//
			Cluster *trackl = (Cluster*)(fTracklets->At(l));
			Cluster *trackm = (Cluster*)(fTracklets->At(m));
			
			//
			printf("now comparing cluster n. %d with cluster n. %d\n",l,m);
			printf("now comparing cluster n. %d with cluster n. %d\n",trackl->getClusterID(),trackm->getClusterID());
			if (trackl->isCompatible(trackm))
			{
				//
				trackl->merge(trackm);
				// make graph
				// trackl->makeGraph();
				// // fit
				// trackl->fitCluster();
				fTracklets->RemoveAt(m);
				fTracklets->Compress();

			}
			else
				m++;
		}
	}
	fClustersNumber = fTracklets->GetEntries();
	std::cout<< "Number of clusters found = " << fClustersNumber << std::endl;
}

void PixelChamberEvent::CheckNoisePoints(){

  printf("find compatibility of noise points with tracklets\n");

  TClonesArray *fnp = fNoise->getPointsObj();
  Datapoint *p;
  for (int l = 0; l < fTracklets->GetEntries(); l++) {
   Cluster *trackl = (Cluster*)(fTracklets->At(l));
   int m = 0;
   double d, d1, d2;
   double *xyz;
   while(m < fnp->GetEntries()) {
      p = (Datapoint*) fnp->At(m);
      d = trackl->DataPointDistanceFromLine(p);
      xyz = trackl->getStartLine();
      d1 = sqrt( (p->getX()-xyz[0])*(p->getX()-xyz[0]) + (p->getY()-xyz[1])*(p->getY()-xyz[1]) + (p->getZ()-xyz[2])*(p->getZ()-xyz[2]) );
      xyz = trackl->getEndLine();
      d2 = sqrt( (p->getX()-xyz[0])*(p->getX()-xyz[0]) + (p->getY()-xyz[1])*(p->getY()-xyz[1]) + (p->getZ()-xyz[2])*(p->getZ()-xyz[2]) );

      if(d < 0.07){
	//      if(d < 0.07 && (d1 < 0.07 || d2 < 0.07)){
	trackl->AddDataPoint(p);
	fNoise->RemoveDataPointAndCompress(m);
      } else
	m++;
    }
    
  }

  this->fitClusters(1);

  int m = 0;
  if(fPoints) delete fPoints;
  fPoints = new TClonesArray("Datapoint",2500);
  printf("fnp->GetEntries()=%d\n",fnp->GetEntries());
  while(m < fnp->GetEntries()) {
    p = (Datapoint*) fnp->At(m);
    printf("p->getUsed()=%d\n",p->getUsed());
    p->setNeighbours(0);
    this->AddDataPoint(p);
    m++;
    printf("placing back noise cluster point %d \n",m);
  }
  
  for(int i=0; i<fnp->GetEntries(); i++) fNoise->RemoveDataPointAndCompress(i);

  this->findClusters(2,15);


}


void PixelChamberEvent::ViewClusters()
{
  int nb = fTracklets->GetEntries();
  printf("nb = %d\n",nb);
  TH3F *h_array[nb];
  TH2F *h_array_XZ[nb];
  TH2F *h_array_XY[nb];
  TH2F *h_array_YZ[nb];
  
  TH3F *h_coord_n = new TH3F("h_coord_n","Clusters_noise", 512*2 , -7.68*2, 7.68*2, 10 , -0.25, 0.25, 512 , -7.68, 7.68);
  TH2F *h_coord_nXZ = new TH2F("h_coord_nXZ","XZ projection", 512*2 , -7.68*2, 7.68*2, 512 , -7.68, 7.68);
  TH2F *h_coord_nXY = new TH2F("h_coord_nXY","XY projection", 512*2 , -7.68*2, 7.68*2, 10 , -0.25, 0.25);
  TH2F *h_coord_nYZ = new TH2F("h_coord_nYZ","YZ projection", 10 , -0.25, 0.25, 512 , -7.68, 7.68);
  TPolyLine *FLXZ[nb], *FLXY[nb], *FLYZ[nb];
  double ReducedChi2[nb];

  for(int i=0; i< nb; i++)
    {
      h_array[i]=new TH3F(Form("h_coord%d",i),"Clusters", 512*2 , -7.68*2, 7.68*2, 10 , -0.25, 0.25, 512 , -7.68, 7.68);
      h_array_XZ[i]=new TH2F(Form("h_coord_XZ%d",i),"XZ projection", 512*2 , -7.68*2, 7.68*2, 512 , -7.68, 7.68);
      h_array_XY[i]=new TH2F(Form("h_coord_XY%d",i),"XY projection", 512*2 , -7.68*2, 7.68*2, 10 , -0.25, 0.25);
      h_array_YZ[i]=new TH2F(Form("h_coord_YZ%d",i),"YZ projection", 10 , -0.25, 0.25, 512 , -7.68, 7.68);
      
      Cluster *track = (Cluster*)(fTracklets->At(i));
      ReducedChi2[i] = track->getReducedChi2();
      FLXZ[i] = (TPolyLine*) track->getPlineXZObj();
      FLXY[i] = (TPolyLine*) track->getPlineXYObj();
      FLYZ[i] = (TPolyLine*) track->getPlineYZObj();
      TClonesArray *pts = (TClonesArray*)(track->getPointsObj());
      //		pts->SetOwner();
      for(int j=0; j<pts->GetEntries(); j++)
	{
	  Datapoint *p = (Datapoint*)(pts->At(j));
	  double x = p->getX();
	  double y = p->getY();
	  double z = p->getZ();
	  h_array[i]->Fill(x,y,z);
	  h_array_XZ[i]->Fill(x,z);
	  h_array_XY[i]->Fill(x,y);
	  h_array_YZ[i]->Fill(y,z);
	}
    }
  
  TClonesArray *n_pts = (TClonesArray*)(fNoise->getPointsObj());
  //	n_pts->SetOwner();
  for(int i=0; i<n_pts->GetEntries(); i++)
    {
      Datapoint *n_p = (Datapoint*)(n_pts->At(i));
      double x_n = n_p->getX();
      double y_n = n_p->getY();
      double z_n = n_p->getZ();
      h_coord_n->Fill(x_n,y_n,z_n);
      h_coord_nXZ->Fill(x_n,z_n);
      h_coord_nXY->Fill(x_n,y_n);
      h_coord_nYZ->Fill(y_n,z_n);
    }
  
  //printf("plotting canvas 1\n");
  TCanvas *c1 =new TCanvas("c1","",1000,1000);
  c1->cd();
  h_coord_n->GetXaxis()->SetTitle("x[mm]");
  h_coord_n->GetXaxis()->SetTitle("y[mm]");
  h_coord_n->GetZaxis()->SetTitle("z[mm]");
  h_coord_n->SetMarkerStyle(2);
  h_coord_n->SetMarkerColor(16);
  h_coord_n->Draw();
  printf("sono qui\n");
  for(int i=0; i<nb; i++){
    h_array[i]->SetMarkerStyle(20+(i%4));
    h_array[i]->SetMarkerSize(1);
    h_array[i]->SetMarkerColor(1+(i%9));
    h_array[i]->Draw("same");    
  }

  TCanvas *c2 =new TCanvas ("c2","",1000,1200);
  c2->cd();
    h_coord_nXZ->GetXaxis()->SetTitle("x[mm]");
    h_coord_nXZ->GetYaxis()->SetTitle("z[mm]");
    h_coord_nXZ->SetMarkerStyle(2);
    h_coord_nXZ->SetMarkerColor(16);
    h_coord_nXZ->Draw("");
    //printf("print lines only if chi2/ndf < %f\n",MAX_CHI2);
    for(int i=0; i<nb; i++){
    	h_array_XZ[i]->SetMarkerStyle(20+(i%4));
    	h_array_XZ[i]->SetMarkerSize(1);
    	h_array_XZ[i]->SetMarkerColor(1+(i%9));
    	h_array_XZ[i]->Draw("same");
	printf("ReducedChi2[%d] = %f\n",i,ReducedChi2[i]);
	//  	if(ReducedChi2[i] < MAX_CHI2) FLXZ[i]->Draw("same");
  	FLXZ[i]->Draw("same");
   }

    TCanvas *c3 =new TCanvas ("c3","",1000,1200);
    c3->cd();
    h_coord_nXY->GetXaxis()->SetTitle("x[mm]");
    h_coord_nXY->GetYaxis()->SetTitle("y[mm]");
    h_coord_nXY->SetMarkerStyle(2);
  	h_coord_nXY->SetMarkerColor(16);
  	h_coord_nXY->Draw();
    for(int i=0; i<nb; i++){
    	h_array_XY[i]->SetMarkerStyle(20+(i%4));
    	h_array_XY[i]->SetMarkerSize(1);
    	h_array_XY[i]->SetMarkerColor(1+(i%9));
    	h_array_XY[i]->Draw("same");
  	FLXY[i]->Draw("same");
    }

    TCanvas *c4 =new TCanvas ("c4","",1000,1200);
    c4->cd();
    h_coord_nYZ->GetXaxis()->SetTitle("y[mm]");
    h_coord_nYZ->GetYaxis()->SetTitle("z[mm]");
    h_coord_nYZ->SetMarkerStyle(2);
  	h_coord_nYZ->SetMarkerColor(16);
  	h_coord_nYZ->Draw();
    for(int i=0; i<nb; i++){
    	h_array_YZ[i]->SetMarkerStyle(20+(i%4));
    	h_array_YZ[i]->SetMarkerSize(1);
    	h_array_YZ[i]->SetMarkerColor(1+(i%9));
    	h_array_YZ[i]->Draw("same");
	FLYZ[i]->Draw("same");
    }


}
