#define NBMIN				2
#define NBMAX				3

#define N     3

#define EPS				1
#define MIN_TRACK_POINTS	        10
#define MAX_CHI2			2
#define MIN_COS				0.9993
//#define MIN_DIST			0.07
#define MIN_DIST			0.07
//#define MIN_DIST_CL			0.15
#define MIN_DIST_CL			0.3
#define STACK                           216
#define FRAC_COMPATIBLE 0.95
#define MAX_PROB 0.0000001

#ifndef DATAPOINT_H
#define DATAPOINT_H

class Datapoint : public TObject
{
public:

  //*********CONSTRUCTOR & DESTRUCTOR*********//

  Datapoint(double x, double y, double z, int i, int j, int k);
  Datapoint();
  Datapoint(Datapoint &p);
  ~Datapoint(){}

  //*********GetTERS*********//

  double GetX(){return fx;}
  double GetY(){return fy;}
  double GetZ(){return fz;}
  int GetI(){return fi;}
  int GetJ(){return fj;}
  int GetK(){return fk;}
  bool GetUsed(){return fused;}
  int GetNeighbours(){return fneighbours;}

  //*********METHODS*********//

  void incNeighbours(){fneighbours++;}
  void setUsed(bool used){fused = used;}
  void setNeighbours(int neighbours){fneighbours = neighbours;}
  bool isNoise(int NB_MIN=NBMIN, int NB_MAX=NBMAX);
  bool isValid(int NB_MIN=NBMIN, int NB_MAX=NBMAX);
  int distance(Datapoint *p2);
  double distance3D(Datapoint *p2);

private:
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
