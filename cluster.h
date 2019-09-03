#ifndef DATAPOINT_H
#define DATAPOINT_H
#include <TGraph2DErrors.h>
#include <TClonesArray.h>
#include <TPolyLine3D.h>		
#include <TPolyLine.h>		
#include <vector>
#include "TH3.h"
#include "TH2.h"
#include <TCanvas.h>

#define NBMIN				2
#define NBMAX				5

class Datapoint : public TObject
{
	public:
		
		//costruttori e distruttore 
		Datapoint(double x, double y, double z, int i, int j, int k);
		Datapoint();
		Datapoint(Datapoint &p); 
		~Datapoint(){}
		//metodi
   
		void setUsed(bool used){fused = used;}
		void setNeighbours(int neighbours){fneighbours = neighbours;}
		double getX(){return fx;}
		double getY(){return fy;}
		double getZ(){return fz;}
		int getI(){return fi;}
		int getJ(){return fj;}
		int getK(){return fk;}
		bool getUsed(){return fused;}
		int getNeighbours(){return fneighbours;}
		void incNeighbours(){fneighbours++;}
		bool isNoise(int NB_MIN=NBMIN, int NB_MAX=NBMAX);
		int distance(Datapoint *p2);
		bool isValid(int NB_MIN=NBMIN, int NB_MAX=NBMAX);

		
	private:
		//dati della classe
		double fx;
		double fy;
		double fz;
		int fi;
		int fj;
		int fk;
		bool fused;
		int fneighbours;

		ClassDef(Datapoint,1); 

};

#endif

//************************************************************************//
#ifndef CLUSTER_H
#define CLUSTER_H
class Cluster : public TObject
{
	public:
		
		//costruttori e distruttore 
		Cluster();
		Cluster(int ClusterID);
		Cluster(Cluster &c); 
		~Cluster();
		//metodi
		void AddDataPoint(Datapoint *p);
		int getClusterID(){return fClusterID;}
		TClonesArray * getPointsObj(){return fPoints;}
		TPolyLine3D * getPline3DObj(){return fPline3D;}
		TPolyLine * getPlineXZObj(){return flXZ;}
		TPolyLine * getPlineXYObj(){return flXY;}
		TPolyLine * getPlineYZObj(){return flYZ;}
		TGraph2DErrors * getGraphObj(){return fGraph;}
		int getNdf(){return fndf;}
		double getChi2(){return fChi2;}
		double getReducedChi2(){return fReducedChi2;}
		double *getParFit(){return fParFit;}
		double *getParFitErrors(){return fParFitErrors;}
		bool isNeighbour(Datapoint *p, int NB_MIN=NBMIN, int NB_MAX=NBMAX);
		void makeGraph();
		void fitCluster();
		bool isCompatible(Cluster *cluster);
		double isCompatible2(Cluster *cluster);
		void merge(Cluster *cluster);
		void ViewClusterFit();
		void RemoveDataPointAndCompress(int ip);
		double DataPointDistanceFromLine(Datapoint *p);
		double *getStartLine();
		double *getEndLine();


	private:
		//
		double pLineDistance(const double *ParFit);
		void line(double t, const double *p, double &x, double &y, double &z);
		//dati della classe
		int fClusterID;
		//TObjArray *fPoints;
		TClonesArray *fPoints;
		TPolyLine3D *fPline3D;
		TGraph2DErrors *fGraph;
		int fndf;
		double fChi2;
		double fReducedChi2;
		double *fParFit;
		double *fParFitErrors;
		TPolyLine *flXZ;
		TPolyLine *flXY;
		TPolyLine *flYZ;
		double ft0;
		double ftf;
		ClassDef(Cluster,1); 

};

#endif

//************************************************************************//

#ifndef SUMDISTANCE2_H
#define SUMDISTANCE2_H

class SumDistance2
{
	public:
		SumDistance2(TGraph2DErrors *graph);
		//		SumDistance2(SumDistance2 &s); 
		~SumDistance2(){}
		double distance2(double x,double y,double z, double ex, double ey, double ez, const double *p);
		double operator() (const double *par);
	private:
		TGraph2DErrors *fGraph;

		//		ClassDef(SumDistance2,1); 

};

#endif

//************************************************************************//

#ifndef PIXELCHAMBEREVENT_H
#define PIXELCHAMBEREVENT_H

class PixelChamberEvent : public TObject
{
	public:
		
		//costruttori e distruttore 
		PixelChamberEvent();
		PixelChamberEvent(int EventID);
		PixelChamberEvent(PixelChamberEvent &e); 
		~PixelChamberEvent();
		//metodi
		void AddDataPoint(Datapoint *p);
		void AddTracklet(Cluster *c);
		void setEventID(int eventID){fEventID = eventID;}
		int getEventID(){return fEventID;}
		TClonesArray *getPointsObj(){return fPoints;}
		TClonesArray *getTrackletsObj(){return fTracklets;}
		Cluster *getNoiseObj(){return fNoise;}
		int getClustersNumber(){return fClustersNumber;}
		int findClusters(int NB_MIN=NBMIN, int NB_MAX=NBMAX);
		void fitClusters(bool kReject=kTRUE);
		void mergeCompatibleClusters();
		void mergeCompatibleClusters2();
		void ViewClusters();
		void CheckNoisePoints();

	private:
		//dati della classe
		TClonesArray *fPoints;
		int fEventID;
		TClonesArray *fTracklets;
		Cluster *fNoise;
		int fClustersNumber;

		ClassDef(PixelChamberEvent,1); 

};

#endif
