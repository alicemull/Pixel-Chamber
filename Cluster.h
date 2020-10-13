#ifndef CLUSTER_H
#define CLUSTER_H
class Cluster : public TObject
{
public:

  //*********CONSTRUCTOR & DESTRUCTOR*********//

  Cluster();
  Cluster(int ClusterID);
  Cluster(Cluster &c);
  ~Cluster();

  //*********GETTERS*********//

  int GetClusterID(){return fClusterID;}
  TClonesArray * GetPointsObj(){return fPoints;}
  /*TPolyLine3D * GetPline3DObj(){return fPline3D;}
  TPolyLine * GetPlineXZObj(){return flXZ;}
  TPolyLine * GetPlineXYObj(){return flXY;}
  TPolyLine * GetPlineYZObj(){return flYZ;}*/
  TGraph2DErrors * GetGraphObj(){return fGraph;}
  int GetNdf(){return fndf;}
  double GetChi2(){return fChi2;}
  double GetReducedChi2(){return fReducedChi2;}
  double *GetParFit(){return fParFit;}
  double *GetParFitErrors(){return fParFitErrors;}
  // double **GetSecondaryVerticesArray(){return fSecondaryVertices;}
  double GetSecondaryVerticesArray(int row, int column){return fSecondaryVertices[row][column];}
  int GetfExPt1() {return fExPt1;}
  int GetfExPt2() {return fExPt2;}
  int GetNumPoints(){return fPoints->GetEntries();}
  double GetVertex_x(){return fVx;}
  double GetVertex_y(){return fVy;}
  double GetVertex_z(){return fVz;}
  double GetErrorVertex_x(){return fEVx;}
  double GetErrorVertex_y(){return fEVy;}
  double GetErrorVertex_z(){return fEVz;}
  int GetPDG_Code(){return fPDG_Code;}
  double GetPx(){return fPx;}
  double GetPy(){return fPy;}
  double GetPz(){return fPz;}
  double GetETot(){return fETot;}
  TMatrixDSym *GetCovarianceMatrix(){return fCovarianceMatrix;}
  TMatrixDSym *GetInverseCovarianceMatrix(){return fInverseCovarianceMatrix;}
  double GetCovMatrixStatus(){return fCovMatrixStatus;}
  bool GetDecompositionFlag(){return fDecompositionFlag;}
  //*************SETTERS******************//
  void SetClusterID(int ClusterID){fClusterID=ClusterID;}
  void SetVertex_x(double vx){fVx = vx;}
  void SetVertex_y(double vy){fVy = vy;}
  void SetVertex_z(double vz){fVz = vz;}
  void SetErrorVertex_x(double evx){fEVx = evx;}
  void SetErrorVertex_y(double evy){fEVy = evy;}
  void SetErrorVertex_z(double evz){fEVz = evz;}
  void SetPDG_Code(int PDG){fPDG_Code = PDG;}
  void SetPx(double Px){fPx=Px;}
  void SetPy(double Py){fPy=Py;}
  void SetPz(double Pz){fPz=Pz;}
  void SetETot(double Etot){fETot=Etot;}

  void SetSecondaryVerticesArray(int row, int column, double value)
  {
        fSecondaryVertices[row][column]=value;
  }
  // void SetSecondaryVerticesArray(double clusterID2, double Vx, double Vy, double Vz, double RedChi2, double covmatrix00, double covmatrix01, double covmatrix02, double covmatrix11, double covmatrix12, double covmatrix22);
  //*********SOME USEFULL METHODS*********//

  void AddDataPoint(Datapoint *p);
  void RemoveDataPointAndCompress(int ip);
  int findFarestPoint(int j, double *xg, double *yg, double *zg, int np);
  bool isNeighbour(Datapoint *p, int NB_MIN=NBMIN, int NB_MAX=NBMAX);
  void PrintCluster();

  //*********FIT METHODS*********//

  void makeGraph(bool KError);
  void findExtremes();
  void fitCluster(double xstart=-999.);
  double Chi2TrackToVertex(double Vx, double Vy, double Vz, TMatrixD *CovVert);
  void SortCluster();
  // double KalmanFilter();
  void KalmanFilter();

  //*********COMPATIBILITY AND MERGE METHODS*********//

  double KalmanToMerge(Cluster* trackm);
  double KalmanToMerge2(Cluster* trackm);
  double DataPointDistanceFromLine(Datapoint *p);
  double ClusterPointsDistance(TClonesArray *fnp, int i_start, int i_end );
  bool isCompatible(Cluster *cluster/*, Double_t min_cos, Double_t min_dist_cl*/);
  bool isCompatibleKalman(Cluster *cluster);
  bool isCompatible2(Cluster *cluster);
  bool CheckIfInsideCluster(TClonesArray *fnp, int i_start, int i_end);
  bool CheckIfInsideCluster(Datapoint *p);
  void merge(Cluster *cluster);
  double ClusterDistance(Cluster *cluster);
  double ClusterDistance(Datapoint *p);

private:
  int fClusterID;
  TClonesArray *fPoints;
  /*TPolyLine3D *fPline3D;
  TPolyLine *flXZ;
  TPolyLine *flXY;
  TPolyLine *flYZ;*/
  TGraph2DErrors *fGraph;
  int fndf;
  double fChi2;
  double fReducedChi2;
  double fParFit[6];
  double fParFitErrors[6];
  double fSecondaryVertices[500][11];
  // double
  int fExPt1;
  int fExPt2;
  double ft0;
  double ftf;
  double fVx;
  double fVy;
  double fVz;
  double fEVx;
  double fEVy;
  double fEVz;
  int fPDG_Code;
  double fPx;
  double fPy;
  double fPz;
  double fETot;
  TMatrixDSym *fCovarianceMatrix;
  TMatrixDSym *fInverseCovarianceMatrix;
  bool fDecompositionFlag;
  int fCovMatrixStatus;

  //*********PRIVATE METHOD*********//

  void line(double t, const double *p, double &x, double &y, double &z);

  ClassDef(Cluster,1);

};

#endif
