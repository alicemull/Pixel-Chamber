#include "math.h"
#include "stdio.h"
#include "math.h"
#include "iostream"     // std::cout
#include "algorithm"    // std::max
#include "stdlib.h"     /* abs */
#include "Math/Vector3D.h"
#include "Math/Vector2D.h"
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
#include "TFile.h"
#include "TTree.h"
#include "TDecompLU.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"
#include "TSystem.h"

using namespace ROOT::Math;

#include "Datapoint.h"
#include "Cluster.h"
#include "SumDistance2.h"
#include "FindChi2.h"
#include "PixelChamberEvent.h"

//*********CONSTRUCTOR & DESTRUCTOR*********//

PixelChamberEvent::PixelChamberEvent() :
TObject(),
fPoints(0x0),
fEventID(0),
fTracklets(0x0),
fNoise(0x0),
fProtonTrack(0x0),
fClustersNumber(0),
fchiVert(0),
fXVert(0),
fYVert(0),
fZVert(0),
fEXVert(0),
fEYVert(0),
fEZVert(0),
fNDFVert(0),

fnpar(0)
{
  fPoints = new TClonesArray("Datapoint",10000);
  fTracklets = new TClonesArray("Cluster",50);
  fNoise = new Cluster(0);
  fProtonTrack = new Cluster(0);
  fCovVert = new TMatrixD(3,3);
}

PixelChamberEvent::PixelChamberEvent(int EventID) :
TObject()
{
  fPoints = new TClonesArray("Datapoint",10000);
  fEventID = EventID;
  fTracklets = new TClonesArray("Cluster",50);
  fNoise = new Cluster(0);
  fProtonTrack = new Cluster(0);
  fClustersNumber = 0;
  fchiVert = -1.;
  fXVert=-16.;
  fYVert=-16.;
  fZVert=-16.;
  fEXVert=999.;
  fEYVert=-999.;
  fEZVert=-999.;
  fNDFVert = 1;
  fnpar = -1;
  fCovVert = new TMatrixD(3,3);
}

PixelChamberEvent::PixelChamberEvent(PixelChamberEvent &e) :
TObject(e),
fPoints(e.fPoints),
fEventID(e.fEventID),
fTracklets(e.fTracklets),
fNoise(e.fNoise),
fProtonTrack(e.fProtonTrack),
fClustersNumber(e.fClustersNumber),
fXVert(e.fXVert),
fYVert(e.fYVert),
fZVert(e.fZVert),
fEXVert(e.fEXVert),
fEYVert(e.fEYVert),
fEZVert(e.fEZVert),
fchiVert(e.fchiVert),
fNDFVert(e.fNDFVert),
fnpar(e.fnpar){
  fCovVert = (TMatrixD*) e.fCovVert->Clone();
}


PixelChamberEvent::~PixelChamberEvent()
{
  fPoints->Clear("C");
  delete fPoints;
  fTracklets->Clear("C");
  delete fTracklets;
  delete fNoise;
  delete fCovVert;
  //  delete fProtonTrack;
}

//*********SOME USEFULL METHODS*********//
//void PixelChamberEvent::SetTracklets(Cluster *cluster)

// TClonesArray *PixelChamberEvent::GetTrackletsObj()
// {
//   printf("Getter fTracklets->GetEntries()=%d\n", fTracklets->GetEntries());
//   return fTracklets;
// }
void PixelChamberEvent::SetTracklets(TClonesArray *Tracklets)
{
  // TClonesArray &arrTemp  = *fTracklets;
  // new (arrTemp[arrTemp.GetEntriesFast()]) Cluster(*cluster);
  //fClustersNumber=arrTemp.GetEntriesFast();
  fTracklets= (TClonesArray *) Tracklets->Clone();
  printf("fTracklets=%d\n",fTracklets->GetEntries());
  fClustersNumber=fTracklets->GetEntries();
  printf("fClustersNumber=%d\n",fClustersNumber);
}
void PixelChamberEvent::AddDataPoint(Datapoint *p)
{
  TClonesArray &arrTemp  = *fPoints;
  new (arrTemp[arrTemp.GetEntriesFast()]) Datapoint(*p);
}

void PixelChamberEvent::AddTracklet(Cluster *c)
{
  TClonesArray &arrTemp  = *fTracklets;
  new (arrTemp[arrTemp.GetEntriesFast()]) Cluster(*c);
}

void PixelChamberEvent::quickSortPointsByI(int l, int r)
{
  // Quick sorting algorithm: it orderes point for consecutive I
  Datapoint *qsort_pi, *qsort_pj, *qsort_pp;
  int i = l,j,p;
  while (i<r)
  {
    i = l;
    j = r;
    p = qsort[(l+r)/2];
    qsort_pp = (Datapoint*)(fPoints->At(p));
    while (i<=j)
    {
      qsort_pi = (Datapoint*)(fPoints->At(qsort[i]));
      while (qsort_pi->GetI() < qsort_pp->GetI())
      {
        i++;
        qsort_pi = (Datapoint*)(fPoints->At(qsort[i]));
      }
      qsort_pj = (Datapoint*)(fPoints->At(qsort[j]));
      while (qsort_pj->GetI() > qsort_pp->GetI())
      {
        j--;
        qsort_pj = (Datapoint*)(fPoints->At(qsort[j]));
      }
      if (i<=j)
      {
        // swap
        int k = qsort[i];
        qsort[i] = qsort[j];
        qsort[j] = k;
        i++;
        j--;
      }
    }
    if(l<j)
    quickSortPointsByI(l, j);
    l = i;
  }
}

//*********METHODS FOR CLUSTERS*********//

int PixelChamberEvent::findClusters(int NB_MIN, int NB_MAX)
{
  // Count neighbours
  //printf("fPoints->GetEntries()=%d\n",fPoints->GetEntries());
  Datapoint *pl, *pm, *pn;
  // Array for sorting in I
  qsort = new int[fPoints->GetEntries()];
  for (int l = 0; l < fPoints->GetEntries(); l++)
  qsort[l] = l;
  // Sorting
  quickSortPointsByI(0, fPoints->GetEntries()-1);
  // for(int m=0; m<fPoints->GetEntries(); m++)
  // {
  //   Datapoint* p= (Datapoint*)fPoints->At(qsort[m]);
  //   printf("i=%d\n",p->GetI() );
  // }
  // Search neighbour, limiting search where I<=EPS
  for (int l = 0; l < fPoints->GetEntries(); l++)
  {
    pl = (Datapoint*)(fPoints->At(qsort[l]));
    // Search backwards
    for (int m = l-1; m >= 0; m--)
    {
      pm = (Datapoint*)(fPoints->At(qsort[m]));
      if (pl->distance(pm) <= EPS)
      {
        pl->incNeighbours();
      }
      // Stop if I > EPS
      if (fabs(pl->GetI() - pm->GetI()) > EPS)
      break;
    }
    // Search forwards
    for (int m = l+1; m < fPoints->GetEntries(); m++)
    {
      pm = (Datapoint*)(fPoints->At(qsort[m]));
      if (pl->distance(pm) <= EPS)
      {
        pl->incNeighbours();
      }
      // Stop if I > EPS
      if (fabs(pl->GetI() - pm->GetI()) > EPS)
      break;
    }
  }
  // Find tracks
  //printf("now find tracks\n");
  Cluster *track;

  TClonesArray *pts = new TClonesArray("Datapoint");
  int tid = 0;
  for (int l = 0; l < fPoints->GetEntries(); l++)
  {
    pl = (Datapoint*)(fPoints->At(l));
    if (pl->isValid(NB_MIN,NB_MAX))
    {
      tid++;
      track = new Cluster( tid );
      track->AddDataPoint(pl);
      // Mark point as used
      pl->setUsed(true);
      // Scans all points from l+1
      int m = l + 1;
      while(m < fPoints->GetEntries())
      {
        // Check if point m is neighbour of current track
        pm = (Datapoint*)(fPoints->At(m));
        if (track->isNeighbour(pm,NB_MIN,NB_MAX))
        {
          // Add point
          track->AddDataPoint(pm);
          // Mark point as used
          pm->setUsed(true);
          // Reset loop
          m = l + 1;
        }
        else
        // Next point
        m++;
      }
      pts = track->GetPointsObj();
      // Check on the number of points in a cluster
      if (pts->GetEntries() >= MIN_TRACK_POINTS)
      {
        this->AddTracklet(track);
      }
      else
      {
        for (int n = 0; n < pts->GetEntries(); n++)
        {
          pn = (Datapoint*)(pts->At(n));
          fNoise->AddDataPoint(pn);
        }
      }
    }
  }
  // Noise
  for (int l = 0; l < fPoints->GetEntries(); l++)
  {
    pn = (Datapoint*)(fPoints->At(l));
    if (pn->isNoise(NB_MIN,NB_MAX))
    fNoise->AddDataPoint(pn);
  }
  // free memory
  // this->ViewClusters();
  delete qsort;
  pts->Clear("C");
  // printf("fNoise->GetNumPoints()=%d\n",fNoise->GetNumPoints());
  // printf("fTracklets->GetEntries()=%d\n",fTracklets->GetEntries());
  return fTracklets->GetEntries();

}

void PixelChamberEvent::ViewClusters()
{
  // Visualization method
  int nb = fTracklets->GetEntries();
  //printf("nb = %d\n",nb);
  //TH3F *h_array[nb];
  TH2F *h_array_XZ[nb];
  TH2F *h_array_XY[nb];
  TH2F *h_array_YZ[nb];

  // TH3F *h_coord_n = new TH3F("h_coord_n","Clusters_noise", 512*2 , -7.68*2, 7.68*2, 512 , -7.68, 7.68 , STACK , -STACK*0.05/2., STACK*0.05/2.);
  TH2F *h_coord_nXZ = new TH2F("h_coord_nXZ","XZ projection", 512*2 , -14.95626, 14.95626, 512 , -6.86784, 6.86784);
  TH2F *h_coord_nXY = new TH2F("h_coord_nXY","XY projection", 512*2 , -14.95626, 14.95626, STACK , -5.4, 5.35);
  TH2F *h_coord_nYZ = new TH2F("h_coord_nYZ","YZ projection", STACK ,-5.4, 5.35, 512 , -6.86784, 6.86784);
  //  TPolyLine *FLXZ[nb], *FLXY[nb], *FLYZ[nb];
  double ReducedChi2[nb];
  double Chi2[nb];
  double Ndf[nb];


  // nb=10;
  for(int i=0; i< nb; i++)
  {
    // h_array[i]=new TH3F(Form("h_coord%d",i),"Clusters", 512*2 , -7.68*2, 7.68*2, 512 , -7.68, 7.68, STACK , -STACK*0.05/2., STACK*0.05/2.);
    h_array_XZ[i]=new TH2F(Form("h_coord_XZ%d",i),"XZ projection", 512*2 , -14.95626, 14.95626, 512 , -6.86784, 6.86784);
    h_array_XY[i]=new TH2F(Form("h_coord_XY%d",i),"XY projection", 512*2 , -14.95626, 14.95626, STACK , -5.4, 5.35);
    h_array_YZ[i]=new TH2F(Form("h_coord_YZ%d",i),"YZ projection", STACK , -5.4, 5.35, 512 , -6.86784, 6.86784);

    Cluster *track = (Cluster*)(fTracklets->At(i));
    ReducedChi2[i] = track->GetReducedChi2();
    Chi2[i] = track->GetChi2();
    Ndf[i] = track->GetNdf();
    // FLXZ[i] = (TPolyLine*) track->GetPlineXZObj();
    // FLXY[i] = (TPolyLine*) track->GetPlineXYObj();
    // FLYZ[i] = (TPolyLine*) track->GetPlineYZObj();
    TClonesArray *pts = (TClonesArray*)(track->GetPointsObj());
    for(int j=0; j<pts->GetEntries(); j++)
    {
      //if(Chi2[i]/Ndf[i] >2)continue;
      Datapoint *p = (Datapoint*)(pts->At(j));
      double x = p->GetX();
      double y = p->GetY();
      double z = p->GetZ();
      // printf("x=%f\n",x);
      // h_array[i]->Fill(x,y,z);
      h_array_XZ[i]->Fill(x,z);
      h_array_XY[i]->Fill(x,y);
      h_array_YZ[i]->Fill(y,z);
    }
  }

  TClonesArray *n_pts = (TClonesArray*)(fNoise->GetPointsObj());
  for(int i=0; i<n_pts->GetEntries(); i++)
  {
    Datapoint *n_p = (Datapoint*)(n_pts->At(i));
    double x_n = n_p->GetX();
    double y_n = n_p->GetY();
    double z_n = n_p->GetZ();
    //printf("x=%f\n",x_n);
    // h_coord_n->Fill(x_n,y_n,z_n);
    h_coord_nXZ->Fill(x_n,z_n);
    h_coord_nXY->Fill(x_n,y_n);
    h_coord_nYZ->Fill(y_n,z_n);
  }

  auto legend = new TLegend(0.1,0.1,0.25,0.9);
  TCanvas *c2 =new TCanvas ("c2","xz",1000,1200);
  h_coord_nXZ->GetXaxis()->SetTitle("x[mm]");
  h_coord_nXZ->GetYaxis()->SetTitle("z[mm]");
  h_coord_nXZ->SetMarkerStyle(2);
  h_coord_nXZ->SetMarkerColor(16);
  h_coord_nXZ->Draw("");
  for(int i=0; i<nb; i++)
  {
    h_array_XZ[i]->SetMarkerStyle(20+(i%4));
    h_array_XZ[i]->SetMarkerSize(1);
    h_array_XZ[i]->SetMarkerColor(1+(i%9));
    h_array_XZ[i]->Draw("same");
    //printf("h_array_XZ[i]->GetEntries()=%f\n",h_array_XZ[i]->GetEntries());
    //printf("ReducedChi2[%d] = %f \n Chi2 / Ndf = %f / %f \n",i,ReducedChi2[i], Chi2[i], Ndf[i]);

    //  FLXZ[i]->Draw("same");
    // legend->AddEntry(h_array_XZ[i],Form(" %d",i),"lep");
    // legend->Draw("same");
  }

  // // auto legend1 = new TLegend(0.1,0.1,0.25,0.9);
  // TCanvas *c3 =new TCanvas ("c3","",1000,1200);
  // c3->cd();
  // h_coord_nXY->GetXaxis()->SetTitle("x[mm]");
  // h_coord_nXY->GetYaxis()->SetTitle("y[mm]");
  // h_coord_nXY->SetMarkerStyle(2);
  // h_coord_nXY->SetMarkerColor(16);
  // h_coord_nXY->Draw();
  // for(int i=0; i<nb; i++){
  //   h_array_XY[i]->SetMarkerStyle(20+(i%4));
  //   h_array_XY[i]->SetMarkerSize(1);
  //   h_array_XY[i]->SetMarkerColor(1+(i%9));
  //   h_array_XY[i]->Draw("same");
  //   //FLXY[i]->Draw("same");
  //   // legend1->AddEntry(h_array_XY[i],Form(" %d",i),"lep");
  //   // legend1->Draw("same");
  // }
  //
  // TCanvas *c4 =new TCanvas ("c4","",1000,1200);
  // c4->cd();
  // h_coord_nYZ->GetXaxis()->SetTitle("y[mm]");
  // h_coord_nYZ->GetYaxis()->SetTitle("z[mm]");
  // h_coord_nYZ->SetMarkerStyle(2);
  // h_coord_nYZ->SetMarkerColor(16);
  // h_coord_nYZ->Draw();
  // for(int i=0; i<nb; i++){
  //   h_array_YZ[i]->SetMarkerStyle(20+(i%4));
  //   h_array_YZ[i]->SetMarkerSize(1);
  //   h_array_YZ[i]->SetMarkerColor(1+(i%9));
  //   h_array_YZ[i]->Draw("same");
  //   //FLYZ[i]->Draw("same");
  // }

}

void PixelChamberEvent::write_cluster_coordinates(int n)
{
  // Method to save coordinates of a particular cluster in a .csv file
  ofstream file_coord;
  file_coord.open(Form("Coordinates_cluster_qsort_%d.csv",n));
  Cluster *cluster = (Cluster*)fTracklets->At(n);
  TClonesArray *pts = cluster->GetPointsObj();
  int n_pts = pts->GetEntries();
  Datapoint *p;
  //printf("\n/*************/\nCluster %d has %d points\n/*************/\n",n, n_pts);
  for(int i=0; i<n_pts; i++)
  {
    p = (Datapoint*) pts->At(i);
    file_coord << p->GetX() << ", " << p->GetY() << ", " << p->GetZ() << std::endl;
  }
  file_coord.close();
}
//*********PROTON METHODS*********//

Cluster *PixelChamberEvent::FindProton()
{
  // Finds proton track looking at tracks at the very beginning of the detector
  int n_tracks = fTracklets->GetEntries();
  int i_p = 0;
  int iFirst = 9999;
  fProtonID = -1;
  //  while(fProtonID == -1.) {
  for(int i = 0; i<n_tracks; i++)
  {
    Cluster *track = (Cluster*)(fTracklets->At(i));
    if(track->GetReducedChi2()>1.5) continue;
    TClonesArray *pts = (TClonesArray*)track->GetPointsObj();
    Datapoint *p = (Datapoint*) pts->At(0);
    i_p = p->GetI();
    if(i_p < iFirst)
    {
      fProtonTrack = track;
      //printf("******proton i=%d**********\n", i);
      fProtonID= i;
      iFirst = i_p;
    }
    //	if(fProtonID != -1) break;
  }
  //    iFirst++;
  //  }
  // printf("iFirst = %d\n",iFirst);
  // printf("proton track has %d points\n",fProtonTrack->GetNumPoints());

  const double *parProton;
  parProton = fProtonTrack -> GetParFit();
  return fProtonTrack;
}

void PixelChamberEvent::viewProton()
{
  // Visualization method to draw proton track
  TH2F *h_coord_XZ = new TH2F("h_coord_XZ","XZ projection", 512*2 , -7.68*2, 7.68*2, 512 , -7.68, 7.68);
  TH2F *h_coord_XY = new TH2F("h_coord_XY","XY projection", 512*2 , -7.68*2, 7.68*2, STACK , -STACK*0.05/2., STACK*0.05/2.);
  TH2F *h_coord_YZ = new TH2F("h_coord_YZ","YZ projection", STACK , -STACK*0.05/2., STACK*0.05/2., 512 , -7.68, 7.68);

  TClonesArray *pts = (TClonesArray*)fProtonTrack->GetPointsObj();
  const double *parLine;
  parLine = fProtonTrack->GetParFit();
  for(int j=0; j<pts->GetEntries(); j++)
  {
    Datapoint *p = (Datapoint*)(pts->At(j));
    double x = p->GetX();
    double y = p->GetY();
    double z = p->GetZ();
    h_coord_XZ->Fill(x,z);
    h_coord_XY->Fill(x,y);
    h_coord_YZ->Fill(y,z);
  }
  TCanvas *c1 = new TCanvas ("c1","",1000,1200);
  c1->Divide(3);
  c1->cd(1);
  h_coord_XZ->Draw();
  c1->cd(2);
  h_coord_XY->Draw();
  c1->cd(3);
  h_coord_YZ->Draw();

}

//*********FIT AND MERGE METHODS*********//

void PixelChamberEvent::fitClusters(bool kReject, bool kError, double xstart)
{
  // Fitt clusters using Claster class' method fitCluster
  int l = 0;
  const double *parfit;
  const double *parfitError;
  Cluster *track;

  //printf("number of tracklets = %d\n", fTracklets->GetEntries());
  while (l < fTracklets->GetEntries())
  {
    // printf("!!!!!!!!!!!!!! \n");
    // printf("fitting track %d\n",l);
    track = (Cluster*)(fTracklets->At(l));
    TClonesArray *pts = track->GetPointsObj();
    // Make graph
    track->makeGraph(kError);
    // Fit
    track->fitCluster(xstart);

    parfit = track->GetParFit();
    parfitError = track->GetParFitErrors();

    if(l==0){
      XYZVector x0(parfit[0], parfit[2], parfit[4]);
      XYZVector x1(parfit[0] + parfit[1], parfit[2] + parfit[3], parfit[4] + parfit[5] ); //another point of line (for t = 1)
      XYZVector u = (x1-x0).Unit();
      double dx, dy, dz, d;
      double dmax = sqrt(0.02924/2.*0.02924/2.+0.025*0.025+0.02688/2.*0.02688/2.);
      Datapoint *p;
      for(int i=0; i<pts->GetEntries(); i++){
        p = (Datapoint *) pts->At(i);
        XYZVector xp(p->GetX(),p->GetY(),p->GetZ());

        dx=((xp-x0).Cross(u)).X();
        dy=((xp-x0).Cross(u)).Y();
        dz=((xp-x0).Cross(u)).Z();
        d=sqrt(dx*dx+dy*dy+dz*dz);
      }
    }
    // Check on fit Reduced chi2; if it is bigger than MAX_CHI2 it is put inside the noise cluster
    if (track->GetReducedChi2()>MAX_CHI2 && kReject)
    {
      fNoise->merge(track);
      fTracklets->Remove(track);
      fTracklets->Compress();
    }
    else
    l++;
  }
}

void PixelChamberEvent::mergeCompatibleClusters()
{
  //Merge tracks that satisfy the "isCompatible"'s conditions
  //printf("find compatibility of %d tracklets (first merge algorithm)\n",fTracklets->GetEntries());
  for (int l = 0; l < fTracklets->GetEntries(); l++)
  {
    int m = l + 1;
    while (m < fTracklets->GetEntries())
    {
      Cluster *trackl = (Cluster*)(fTracklets->At(l));
      Cluster *trackm = (Cluster*)(fTracklets->At(m));
      if (trackl->isCompatible(trackm/*, MIN_COS, MIN_DIST_CL*/))
      {
        trackl->merge(trackm);
        // New fit of merged tracks
        trackl->makeGraph(1);
        trackl->fitCluster();
        fTracklets->RemoveAt(m);
        fTracklets->Compress();
        m = l + 1;
      }
      else
      m++;
    }

  }

  fClustersNumber = fTracklets->GetEntries();
  //std::cout<< "Number of clusters found = " << fTracklets->GetEntries() << std::endl;
}


void PixelChamberEvent::mergeCompatibleClusters2()
{
  //Merge tracks if reduced Chi2 of tracks is smaller than MAX_CHI2 and "isCompatible2"'s conditions are satisfied
  //printf("find compatibility of %d tracklets (second algorithm)\n",fTracklets->GetEntries());
  for (int l = 0; l < fTracklets->GetEntries(); l++)
  {
    // Loop to compare all clusters
    int m=l+1;
    while (m < fTracklets->GetEntries())
    {
      Cluster *trackl = (Cluster*)(fTracklets->At(l));
      Cluster *trackm = (Cluster*)(fTracklets->At(m));
      TClonesArray *pts = (TClonesArray*)trackm->GetPointsObj();
      if (trackl->isCompatible2(trackm))      {
        trackl->merge(trackm);
        // Fit
        trackl->makeGraph(0);
        trackl->fitCluster();
        fTracklets->RemoveAt(m);
        fTracklets->Compress();
        m = l+1;
      }
      else
      m++;
    }
  }
  fClustersNumber = fTracklets->GetEntries();
  std::cout<< "Number of clusters found = " << fClustersNumber << std::endl;
}

void PixelChamberEvent::attachNoisePoints()
{
  // Checks if there are noise points compatibles with clusters
  //printf("#####1 fNoise->GetNumPoints() = %d 1#####\n", fNoise->GetNumPoints());
  TClonesArray *fnp = fNoise->GetPointsObj(); //Non sto creando un oggetto ma prendendo il puntatore ai punti di fNoise
  Datapoint *p;

  for (int l = 0; l < fTracklets->GetEntries(); l++)
  {
    // Loop to all tracks
    Cluster *trackl = (Cluster*)(fTracklets->At(l));
    int m = 0;
    double d, d1, d2;
    while(m < fnp->GetEntries())
    {
      // Calculate distances between noise points and clusters
      p = (Datapoint*) fnp->At(m);
      // Check distance between points and cluster points
      d = trackl->DataPointDistanceFromLine(p);
      TGraph2DErrors *Graph = trackl->GetGraphObj();
      double *xg = Graph->GetX();
      double *yg = Graph->GetY();
      double *zg = Graph->GetZ();
      int np = Graph->GetN();

      // Check the distance between the extreme points of clusters and noise points
      int ex1 = trackl->GetfExPt1();
      int ex2 = trackl->GetfExPt2();
      double x_start = xg[ex1];
      double y_start = yg[ex1];
      double z_start = zg[ex1];
      double x_end = xg[ex2];
      double y_end = yg[ex2];
      double z_end = zg[ex2];

      d1 = sqrt( (p->GetX()-x_start)*(p->GetX()-x_start) + (p->GetY()-y_start)*(p->GetY()-y_start) + (p->GetZ()-z_start)*(p->GetZ()-z_start) );
      d2 = sqrt( (p->GetX()-x_end)*(p->GetX()-x_end) + (p->GetY()-y_end)*(p->GetY()-y_end) + (p->GetZ()-z_end)*(p->GetZ()-z_end) );

      if(d<0.04 && trackl->CheckIfInsideCluster(p) && trackl->GetReducedChi2()<2.5)
      {
        p->setUsed(true);
        trackl->AddDataPoint(p);
        fNoise->RemoveDataPointAndCompress(m);
        // ricalcola estremi
        trackl->makeGraph(0);
        trackl->findExtremes();
        //
        m = 0;
      }
      else
      m++;
    }
  }
}

void PixelChamberEvent::moveNoiseToPoints()
{
  // Checks if there are noise points compatibles with clusters
  TClonesArray *fnp = fNoise->GetPointsObj(); //Non sto creando un oggetto ma prendendo il puntatore ai punti di fNoise
  Datapoint *p;

  int m = 0;
  // Deleting and updating points
  if(fPoints) delete fPoints;
  fPoints = new TClonesArray("Datapoint",10000);
  while(m < fnp->GetEntries())
  {
    p = (Datapoint*) fnp->At(m);
    p->setNeighbours(0);
    this->AddDataPoint(p);
    m++;
  }
  fnp->Clear("C");//l'ho aggiunto io ora. Ed è giusto cosi!!!! deve svuotare fNoise. DM
}

//*********WRITE INFORMATIONS FOR VERTEX FINDER**********//
void PixelChamberEvent::Write_ttree(int TTreeID, int Event_ID)
{
  TFile *fOut = new TFile(Form("ClustersInfo_Tree%d.root", TTreeID), "RECREATE");
  TClonesArray *Clusters = new TClonesArray("Cluster",200);
  TClonesArray &ToWriteTemp = *Clusters;

  int n_tracks = fTracklets->GetEntries();
  const double *parFit;
  const double *parFitErrors;
  int TTree_ID, Ev_ID;
  TTree *ClustersTree = new TTree("ClustersInfo", "Info");
  ClustersTree->Branch("TreeID", &TTree_ID, "TTree_ID/I");
  ClustersTree->Branch("EvID", &Ev_ID, "Ev_ID/I");
  ClustersTree->Branch("clusters", &Clusters);

  // printf("writing tree\n");
  // printf("n_tracks=%d\n",n_tracks);
  int iv=0;
  Cluster *c0;
  Cluster *cluster;
  Cluster &c = *cluster;
  for(int i = 0; i<n_tracks; i++)
  {
    cluster = (Cluster*) fTracklets->At(i);
    if(cluster->GetVertex_x()>-20.)new(ToWriteTemp[ToWriteTemp.GetEntries()])(Cluster)(*cluster);

  }
  //printf ("Clusters array to be written: %p\n",Clusters);
  TTree_ID = TTreeID;
  Ev_ID = Event_ID;
  ClustersTree->Fill();
  fOut->Write();
  fOut->Close();
}
//*********FINDING AND FITTING VERTICES*********//

double PixelChamberEvent::FitVertex(TClonesArray *Tracklets, int n_tracks, double CT)
{
  // printf("FitVertex invoked with n_tracks=%d and CT=%f\n",n_tracks,CT);
  //printf("xm=%f\n",xm);
  // This method uses the class FindChi2 and the root fitter
  // It's a new fit for tracks asking that new fit lines pass in the same point
  Cluster *VertexClusters=new Cluster();
  // Number of parmeters
  //  int n = 3+3*n_tracks;
  //int n = 3+2*n_tracks;
  int n = 3;
  double *p = new double[n];
  // Parameters initialization

  const double *parLine;

  //Questo viene settato come seed dall'operazione in find vertex con il setter.
  //Nel caso del vertice primario è l'ultimo punto del protone
  p[0] = fXVert;
  p[1] = fYVert;
  p[2] = fZVert;
  //  printf("p[0]=%f p[1]=%f p[2]=%f\n",p[0],p[1],p[2]);
  ROOT::Fit::Fitter  fitter1;

  double *WT;
  FindChi2 chiPrevIteration(Tracklets, CT);
  chiPrevIteration.operator()(p);
  WT =  chiPrevIteration.GetWT();
  // for(int i = 0; i<n_tracks; i++)
  //      printf("*WT[%d]=%f\n",i,*(WT+i));

  // force weight of proton track to 1
  //  int ProtonID = this->GetProtonTrackID();
  // // printf("proton id inside fit vertex=%d\n",ProtonID);
  //  *(WT+ProtonID) = 1.;

  FindChi2 chi(Tracklets,CT,WT);

  ROOT::Math::Functor fcn(chi, n);
  fitter1.Config().SetMinimizer("Minuit2","Migrad");
  fitter1.SetFCN(fcn,p,n_tracks*2,1);
  //  fitter1.SetFCN(fcn,p,n_tracks*5,1);
  for (int j = 0; j < n; ++j)
  {
    fitter1.Config().ParSettings(j).SetStepSize(0.001);
    fitter1.Config().ParSettings(j).SetLimits(-15.,15.);

  }
  // for(int i = 0; i<n_tracks; i++)
  // {
  //   fitter1.Config().ParSettings(3+2*i).Fix();
  //   fitter1.Config().ParSettings(3+2*i+1).Fix();
  //   //    p[3+2*i+2]=parLine[1];
  // }

  // Saving fit results
  bool ok = fitter1.FitFCN();

  const ROOT::Fit::FitResult & result1 = fitter1.Result();

  double chi2 = result1.Chi2();
  // fNDFVert = result1.Ndf();
  fNDFVert = 2. * n_tracks - 3;
  //  printf("result1.Chi2()=%f result1.Ndf() = %d\n",result1.Chi2(),result1.Ndf());
  fnpar = result1.NPar();
  int CovMatrixStatus = result1.CovMatrixStatus();

  const double *ParFit1 ;
  ParFit1 = result1.GetParams();
  const double *ParFitErrors1;
  ParFitErrors1 = result1.GetErrors();
  fchiVert=chi2;

  //  for(int i=0;i<n; i++) //printf("i= %d; par = %f; err = %f \n", i, ParFit1[i], ParFitErrors1[i]);
  fXVert=ParFit1[0];
  fYVert=ParFit1[1];
  fZVert=ParFit1[2];
  fEXVert=ParFitErrors1[0];
  fEYVert=ParFitErrors1[1];
  fEZVert=ParFitErrors1[2];

  int nVT=0;
  chi.operator()(ParFit1);
  WT =  chi.GetWT();
  // *(WT+ProtonID) = 1.;
  for(int i = 0; i<n_tracks; i++){
    //printf("*WT[%d]=%f\n",i,*(WT+i));
    VertexClusters = (Cluster*)Tracklets->At(i);
    VertexClusters->SetVertex_x(-16.);
    VertexClusters->SetVertex_y(-16.);
    VertexClusters->SetVertex_z(-16.);
    if(*(WT+i) > 0.) {
      nVT++;
      VertexClusters->SetVertex_x(fXVert);
      VertexClusters->SetVertex_y(fYVert);
      VertexClusters->SetVertex_z(fZVert);
    }
  }

  //  printf("%d associated track to vertex CT=%f\n",nVT,CT);
  TMatrixD *CovVert = new TMatrixD(3,3);
  for(int cov_i=0; cov_i<3; cov_i++)
  {
    for(int cov_j=0; cov_j<3; cov_j++)
    {
      (*CovVert)[cov_i][cov_j]=result1.CovMatrix(cov_i,cov_j);
    }
  }
  fCovVert=(TMatrixD*)CovVert->Clone();
  // printf("vertex covariance matrix\n");
  // CovVert->Print();
  // Calculate track to vertex chi2
  double Chi2Track=0;
  for(int i = 0; i<n_tracks; i++)
  {
    VertexClusters = (Cluster*)Tracklets->At(i);
    if(*(WT+i) > 0.) Chi2Track += VertexClusters->Chi2TrackToVertex(fXVert, fYVert, fZVert,CovVert);
  }
  fNDFVert = 2. * n_tracks;
  fchiVert = Chi2Track;
  // VertexClusters = (Cluster*)Tracklets->At(n_tracks-1);
  // Chi2Track = VertexClusters->Chi2TrackToVertex(fXVert, fYVert, fZVert,CovVert);

  //  printf("fchiVert=%f fNDFVert=%d\n",fchiVert,fNDFVert);
  //printf("xvert=%f yvert=%f zvert=%f\n",fXVert,fYVert,fZVert);
  delete[] p;
  //  chi.Clean();
  return Chi2Track;
}

void PixelChamberEvent::FindPrimaryVertex()
{
  // Finds Vertex starting making the fit (VertexFit) of tracks starting from proton track
  TClonesArray *PrimaryVertClusters = new TClonesArray("Cluster");
  TClonesArray &ClustVertTemp = *PrimaryVertClusters;

  // Gets proton track and adds it to the TClones array containing tracks to be fitted (VertexFit)

  int n_tracks = fTracklets->GetEntries();
  for(int i = 0; i<n_tracks; i++){
    Cluster *tr = (Cluster*) fTracklets->At(i);
    tr->SetClusterID(i);
  }

  printf("****** Vertexing first stage ******\n");

  Cluster *ProtonTrack = (Cluster*)this->FindProton();
  int ProtonID = this->GetProtonTrackID();
  const double *ParFitProton = ProtonTrack->GetParFit();
  double x0 = ParFitProton[0];
  double y0 = ParFitProton[2];
  double z0 = ParFitProton[4];

  // seed
  this->SetXVert(x0);
  this->SetYVert(y0);
  this->SetZVert(z0);
  new(ClustVertTemp[ClustVertTemp.GetEntries()])(Cluster)(*ProtonTrack);

  const double *parLine;
  const double *errorsLine;

  //************//create pool of tracks for vertexing
  //************//
  Cluster *ClusterToAdd;

  for(int i = 0; i<n_tracks; i++)
  {
    // Loop on all tracks
    ClusterToAdd = (Cluster*) fTracklets->At(i);
    Cluster *c0 = new Cluster(*ClusterToAdd);

    int n_pts = ClusterToAdd->GetNumPoints();
    //TMatrixD *Cov0 = (TMatrixD*) c0->GetCovarianceMatrix();
    TMatrixDSym *SubMatrix = (TMatrixDSym*) c0->GetCovarianceMatrix();

    bool Flag = c0->GetDecompositionFlag();
    if (!Flag) {
      cout << "Decomposition failed, matrix singular ?" << endl;
      //cout << "condition number = " << = Cov0->GetCondition() << endl;
    }
    double DetCTrack = SubMatrix->Determinant();
    if(c0->GetClusterID()!=ProtonID & DetCTrack>0.& Flag==true &  c0->GetReducedChi2()<1.5 && n_pts>50)
    new(ClustVertTemp[ClustVertTemp.GetEntries()])(Cluster)(*c0);
  }

  double CT[15] = {1000000.,100000.,10000.,3000.,200.,100.,20.,10.,9.,8.,7.,6.,5.,4.,3.};
  int nIterations = 30;
  for(int i = 0; i<7; i++)
  {
    this->FitVertex(PrimaryVertClusters, PrimaryVertClusters->GetEntries(),CT[i]);
  }
  //
  printf("PrimaryVertClusters->GetEntries()=%d\n",PrimaryVertClusters->GetEntries());
  for(int i = 0; i<PrimaryVertClusters->GetEntries(); i++)
  {
    Cluster *tr1 = (Cluster*) PrimaryVertClusters->At(i);
    //printf("vertex track %d\n",tr1->GetClusterID());
    for(int j = 0; j<fTracklets->GetEntries(); j++)
    {
      Cluster *tr = (Cluster*) fTracklets->At(j);
      if(tr1->GetClusterID()==tr->GetClusterID())
      {
        tr->SetVertex_x(tr1->GetVertex_x());
        tr->SetVertex_y(tr1->GetVertex_y());
        tr->SetVertex_z(tr1->GetVertex_z());
        tr->SetErrorVertex_x(this->GetErrorXVert());
        tr->SetErrorVertex_y(this->GetErrorYVert());
        tr->SetErrorVertex_z(this->GetErrorZVert());
      }
    }
  }
  printf("Primary Vertex fit: final ch2/ndf=%f\n",this->GetChiVert()/this->GetNDFVert());
  printf("Vx=%f+-%f Vy=%f+-%f Vz=%f+-%f\n",this->GetXVert(),this->GetErrorXVert(),this->GetYVert(),this->GetErrorYVert(),this->GetZVert(),this->GetErrorZVert());
}

void PixelChamberEvent::FindSecondaryVertex()
{
  printf("Starting secondary vertexing\n");
  // Finds Vertex starting making the fit (VertexFit) of tracks starting from proton track
  TClonesArray *SecondaryVertClusters = new TClonesArray("Cluster");
  TClonesArray &ClustVertTemp = *SecondaryVertClusters;

  // Gets proton track and adds it to the TClones array containing tracks to be fitted (VertexFit)

  int n_tracks = fTracklets->GetEntries();
  // for(int i = 0; i<n_tracks; i++){
  //   Cluster *tr = (Cluster*) fTracklets->At(i);
  //   tr->SetClusterID(i);
  // }

  printf("****** Vertexing first stage ******\n");

  double primaryVertX=this->GetXVert();
  double primaryVertY=this->GetYVert();
  double primaryVertZ=this->GetZVert();
  // double vertices[n_tracks][n_tracks][12];
  if(this->GetChiVert()/this->GetNDFVert()<2.5)
  {
    Datapoint* Primary_vertex=new Datapoint(primaryVertX, primaryVertY, primaryVertZ, 0, 0, 0);

    const double *parLine;
    const double *errorsLine;

    Cluster *ProtonTrack = (Cluster*)this->FindProton();
    int ProtonID = this->GetProtonTrackID();
    //************//create pool of tracks for vertexing
    //************//
    Cluster *ClusterToAdd;
    int n=0;
    int leftTracks=0;
    // for(int i=0; i<n_tracks; i++)
    // {
    //   Cluster *c2 = (Cluster*) fTracklets->At(i);
    //   if(c2->GetReducedChi2()<2.5|c2->GetNumPoints()>50|c2->GetVertex_x()==-16.)
    //   leftTracks++;
    // }
    // printf("leftTracks=%d\n",leftTracks);
    while(n<n_tracks)
    {

      double vertices[n_tracks][11];
      Cluster *c1 = (Cluster*) fTracklets->At(n);
      TMatrixDSym *SubMatrix1 = (TMatrixDSym*) c1->GetCovarianceMatrix();
      double DetCTrack1 = SubMatrix1->Determinant();
      int n_pts1 = c1->GetNumPoints();
      // printf("c1->GetClusterID()=%d, DetCTrack1=%f, c1->GetReducedChi2()=%f, n_pts1=%d, c1->GetVertex_x()=%f\n",c1->GetClusterID(), DetCTrack1, c1->GetReducedChi2(), n_pts1, c1->GetVertex_x());
      int ID=c1->GetClusterID();
      // if(ID==ProtonID|c1->GetVertex_x()!=-16.|c1->GetReducedChi2()>=2.5|n_pts1<50)
      if(ID==ProtonID|c1->GetVertex_x()!=-16.|c1->GetReducedChi2()>=2.5|n_pts1<50)
      {
        n++;
        continue;
      }
      // printf("leftTracks=%d\n",leftTracks);

      // ProcInfo_t *info = new ProcInfo_t;
      // gSystem->GetProcInfo(info);
      // Long_t memo =info->fMemResident;
      // printf ("Iteration %d -- Resident memory used by root: %d\n", n, memo);
      // delete info;
      int ext1_c1 = c1->GetfExPt1();
      int ext2_c1 = c1->GetfExPt2();
      TClonesArray *pts_c1=(TClonesArray*)c1->GetPointsObj();
      Datapoint* pt1_c1 = (Datapoint*)pts_c1->At(ext1_c1);
      Datapoint* pt2_c1 = (Datapoint*)pts_c1->At(ext2_c1);

      double distance1_c1 = pt1_c1->distance3D(Primary_vertex);
      double distance2_c1 = pt2_c1->distance3D(Primary_vertex);

      if(c1->GetClusterID()!=ProtonID && c1->GetReducedChi2()<2.5 && n_pts1>50 && c1->GetVertex_x()==-16.)
      {
        new(ClustVertTemp[ClustVertTemp.GetEntries()])(Cluster)(*c1);
        ID=c1->GetClusterID();
      }

      for(int i = 0; i<n_tracks; i++)
      {
        // ProcInfo_t *info = new ProcInfo_t;
        // gSystem->GetProcInfo(info);
        // Long_t memo =info->fMemResident;
        // printf ("Iteration %d -- Resident memory used by root: %d\n", i, memo);
        // delete info;
        if(abs(distance1_c1)<=abs(distance2_c1))
        {
          this->SetXVert(pt1_c1->GetX());
          this->SetYVert(pt1_c1->GetY());
          this->SetZVert(pt1_c1->GetZ());
          // if(c1->GetClusterID()==24)printf("this->GetXVert()=%f\n", this->GetXVert());
        }

        if(abs(distance1_c1)>abs(distance2_c1))
        {
          this->SetXVert(pt2_c1->GetX());
          this->SetYVert(pt2_c1->GetY());
          this->SetZVert(pt2_c1->GetZ());
          // if(c1->GetClusterID()==24)printf("this->GetXVert()=%f\n", this->GetXVert());
        }
        // printf("this->GetXVert()=%f\n", this->GetXVert());
        ClusterToAdd = (Cluster*) fTracklets->At(i);
        int n_pts = ClusterToAdd->GetNumPoints();
        if(ClusterToAdd->GetReducedChi2()>=2.5|n_pts<50|ClusterToAdd->GetVertex_x()!=-16.)continue;
        //TMatrixD *Cov0 = (TMatrixD*) c0->GetCovarianceMatrix();
        Cluster *c0 = new Cluster(*ClusterToAdd);

        if(c0->GetClusterID()!=ID && c0->GetClusterID()!=ProtonID &&  c0->GetReducedChi2()<2.5 && n_pts>50 && c0->GetVertex_x()==-16.)
        {
            new(ClustVertTemp[ClustVertTemp.GetEntries()])(Cluster)(*c0);
            // printf("OK!\n");
        }
        if(SecondaryVertClusters->GetEntries()==2)
        {
          double CT[15] = {1000000.,100000.,10000.,3000.,200.,100.,20.,10.,9.,8.,7.,6.,5.,4.,3.};
          int nIterations = 30;
          // printf("OK!\n");
          for(int i = 0; i<7; i++)
          {
            this->FitVertex(SecondaryVertClusters, SecondaryVertClusters->GetEntries(),CT[i]);
          }

          vertices[i][0]=ClusterToAdd->GetClusterID();
          vertices[i][1]=this->GetXVert();
          vertices[i][2]=this->GetYVert();
          vertices[i][3]=this->GetZVert();
          vertices[i][4]=this->GetChiVert()/this->GetNDFVert();
          vertices[i][5]=(*fCovVert)[0][0];
          vertices[i][6]=(*fCovVert)[0][1];
          vertices[i][7]=(*fCovVert)[0][2];
          vertices[i][8]=(*fCovVert)[1][1];
          vertices[i][9]=(*fCovVert)[1][2];
          vertices[i][10]=(*fCovVert)[2][2];

          // for(int k=0; k<11; k++)printf("TrackID=%d, vertices[%d][%d]=%f\n",c1->GetClusterID(),i,k,vertices[i][k]);

          for(int h=0; h<11; h++)
          {
            c1->SetSecondaryVerticesArray(i,h,vertices[i][h]);
          }


          if(this->GetChiVert()/this->GetNDFVert()>=2.5&&SecondaryVertClusters->GetEntries()!=2)
          {
            // for(int i = 0; i<SecondaryVertClusters->GetEntries(); i++)
            // {
            //   Cluster *tr1 = (Cluster*) SecondaryVertClusters->At(i);
            //   for(int j = 0; j<fTracklets->GetEntries(); j++)
            //   {
            //     Cluster *tr = (Cluster*) fTracklets->At(j);
            //   }
            // }
            SecondaryVertClusters->RemoveAt(1);
            SecondaryVertClusters->Compress();
            continue;
          }
        }

        else{
          for(int h=0; h<11; h++)
          {
            vertices[i][h]=-16.00;
          }
        }

        for(int i=0; i<SecondaryVertClusters->GetEntries(); i++)
        {
          Cluster*track=(Cluster*)SecondaryVertClusters->At(i);

          // printf("Tracks associated to vertex: %d\n", track->GetClusterID());
          if(i>0)
          {
            SecondaryVertClusters->RemoveAt(1);
            SecondaryVertClusters->Compress();
          }
        }
        if(SecondaryVertClusters->GetEntries()<2)continue;

      }
      SecondaryVertClusters->Clear("C");
      n++;
      }
      // delete ClusterToAdd;
    }
    printf("Ciaone\n");
}


void PixelChamberEvent::write_vertex_info()
{
  // Method to save vertex info in a .csv file
  double XV = this->GetXVert();
  double eXV = this->GetErrorXVert();
  double YV = this->GetYVert();
  double eYV = this->GetErrorYVert();
  double ZV = this->GetZVert();
  double eZV = this->GetErrorZVert();
  double Chi2V = this->GetChiVert();
  int NDFV = this->GetNDFVert();
  int NPAR = this->GetNParVert();
  int eventID = this->GetEventID();
  //printf("\n/*************/\n/*************/ \n Event %d Vertex: \n x = %f +- %f \n y = %f +- %f \n z = %f +- %f \n ReducedChi2 = %f \n/*************/\n/*************/\n", eventID, XV, eXV, YV, eYV, ZV, eZV, Chi2V/NDFV);
}
