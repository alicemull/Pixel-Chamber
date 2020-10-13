#ifndef PIXELCHAMBEREVENT_H
#define PIXELCHAMBEREVENT_H

class PixelChamberEvent : public TObject
{
public:

  //*********CONSTRUCTOR & DESTRUCTOR*********//

  PixelChamberEvent();
  PixelChamberEvent(int EventID);
  PixelChamberEvent(PixelChamberEvent &e);
  ~PixelChamberEvent();

  //*********GETTERS*********//

  int GetEventID(){return fEventID;}
  int GetClustersNumber(){return fClustersNumber;}
  TClonesArray *GetTrackletsObj(){return fTracklets;}
  // TClonesArray *GetTrackletsObj();
  TClonesArray *GetPointsObj(){return fPoints;}
  Cluster *GetNoiseObj(){return fNoise;}
  int GetProtonTrackID(){return fProtonID;}
  Cluster *GetProtonTrack(){return fProtonTrack;}
  double GetChiVert(){return fchiVert;}
  int GetNDFVert(){return fNDFVert;}
  double GetXVert(){return fXVert;}
  double GetYVert(){return fYVert;}
  double GetZVert(){return fZVert;}
  void SetXVert(double vx){fXVert = vx;}
  void SetYVert(double vy){fYVert = vy;}
  void SetZVert(double vz){fZVert = vz;}
  double GetErrorXVert(){return fEXVert;}
  double GetErrorYVert(){return fEYVert;}
  double GetErrorZVert(){return fEZVert;}
  void setEventID(int eventID){fEventID = eventID;}
  int GetNParVert(){return fnpar;}
  TMatrixD *GetCovVert(){return fCovVert;}

  //********************************//
  void ClearTracklets(){fTracklets->Clear("C");}
  //*********SOME USEFULL METHODS*********//
  //void SetTracklets(Cluster *cluster);
  void SetTracklets(TClonesArray *Tracks);//{fTracklets=Tracks;}
  void AddDataPoint(Datapoint *p);
  void AddTracklet(Cluster *c);
  void quickSortPointsByI(int l, int r);

  //*********METHODS FOR CLUSTERS*********//

  int findClusters(int NB_MIN=NBMIN, int NB_MAX=NBMAX);
  void ViewClusters();
  void write_cluster_coordinates(int n);

  //*********PROTON METHODS*********//

  Cluster *FindProton();
  void viewProton();

  //*********FIT AND MERGE METHODS*********//

  void fitClusters(bool kReject=kTRUE, bool kError = kFALSE, double xstart=-999.);
  void mergeCompatibleClusters();
  void mergeCompatibleClusters2();
  void attachNoisePoints();
  void moveNoiseToPoints();

  //*************WRITE TTREE**********************//
  void Write_ttree(int TreeID, int Event_ID);
  //*********FINDING AND FITTING VERTICES*********//

  double FitVertex(TClonesArray *Tracklets, int i, double xm);
  void FindPrimaryVertex();
  void FindSecondaryVertex();
  void write_vertex_info();

private:
  int fEventID;
  int fClustersNumber;
  TClonesArray *fTracklets;
  TClonesArray *fPoints;
  Cluster *fNoise;
  int fProtonID;
  Cluster *fProtonTrack;
  double fchiVert;
  int fNDFVert;
  double fXVert;
  double fYVert;
  double fZVert;
  double fEXVert;
  double fEYVert;
  double fEZVert;
  TMatrixD *fCovVert;

  int *qsort;
  int fnpar;
  ClassDef(PixelChamberEvent,1);

};

#endif
