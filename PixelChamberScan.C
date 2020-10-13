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
#include "TNtuple.h"
#include "TFile.h"
#include "TSystem.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"

using namespace ROOT::Math;
using namespace std;

#include "Datapoint.h"
#include "SumDistance2.h"
#include "Cluster.h"
#include "FindChi2.h"
#include "PixelChamberEvent.h"

void PixelChamberScan(Int_t myevent=0, Int_t myevent2 = 0, Int_t TreeNumber = 0)
{

  PixelChamberEvent *event;
  //TFile *file1 = new TFile(Form("/gr3-str/PixelChamber/Geant/TTree/NoCharm/new_gauss5000/Proton%d/inelasticTree_216Planes1tree.root",TreeNumber));
  //TFile *file1 = new TFile("/Users/alicemulliri/Desktop/PixelChamber/TTree/DPlus/newMergeDPlus.root");
  //  TFile *file1 = new TFile(Form("/gr3-str/PixelChamber/Geant/TTree/NoCharm/new_gauss5000/Proton%d/newMergeD0250Ev.root",TreeNumber));
  // TFile *file1 = new TFile(Form("/gr3-str/PixelChamber/Geant/TTree/NoCharm/new_gauss5000/Proton%d/inelasticTree_216Planes_hit1tree.root",TreeNumber));
  TFile *file1 = new TFile(Form("/Volumes/AliceM_SSD/Varie/PixelChamberAll/backup_pixelchamber/DPlus/newMergeDPlus.root",TreeNumber));
  TNtuple *ntuple = (TNtuple*) file1->Get("Interest_D");
  //Takes Ntuple
  Int_t  EventID, ic, jc, kc, ParticlePDGCode;
  Double_t  x, y, z, TrackID_1, TrackID_2, TrackID_3, Mom_x, Mom_y, Mom_z;
  //  ntuple->SetBranchAddress("EventID_D",&EventID);
  ntuple->SetBranchAddress("EventID_D",&EventID);
  ntuple->SetBranchAddress("i",&ic);
  ntuple->SetBranchAddress("j",&jc);
  ntuple->SetBranchAddress("k",&kc);
  //  ntuple->SetBranchAddress("ParticlePDGCode_D",&ParticlePDGCode);
  ntuple->SetBranchAddress("ParticlePDGCode",&ParticlePDGCode);
  ntuple->SetBranchAddress("x",&x);
  ntuple->SetBranchAddress("y",&y);
  ntuple->SetBranchAddress("z",&z);
  ntuple->SetBranchAddress("TrackID",&TrackID_1);
  /*ntuple->SetBranchAddress("TrackID_1",&TrackID_1);
  ntuple->SetBranchAddress("TrackID_2",&TrackID_2);
  ntuple->SetBranchAddress("TrackID_3",&TrackID_3);*/
  ntuple->SetBranchAddress("px", &Mom_x);
  ntuple->SetBranchAddress("py", &Mom_y);
  ntuple->SetBranchAddress("pz", &Mom_z);

  Int_t nentries = ntuple->GetEntries();

  //**********************************************************************************//
  TFile *fOut = new TFile(Form("ClustersDPlus%dNoMerge.root",TreeNumber), "RECREATE");
  TClonesArray *Clusters = new TClonesArray("Cluster");
  //TClonesArray &ToWriteTemp = *Clusters;
  int TTree_ID, Ev_ID;
  double Chi2Vert;
  TTree *ClustersTree = new TTree("ClustersInfo", "Info");
  ClustersTree->Branch("TreeID", &TTree_ID, "TTree_ID/I");
  ClustersTree->Branch("EvID", &Ev_ID, "Ev_ID/I");
  ClustersTree->Branch("Chi2Vert", &Chi2Vert, "Chi2Vert/D");
  ClustersTree->Branch("clusters", &Clusters);
  Cluster *c0 = new Cluster();
  Cluster *cluster;
  //**********************************************************************************//
  int l_event = 0;
  int l_event1=0;
  Datapoint *pt;
  Clusters->BypassStreamer();
  printf("start reconstruction\n");

  for (int ievent = myevent; ievent<myevent2; ievent++)
  {
    // New event
    l_event1=l_event;
    TClonesArray &ToWriteTemp = *Clusters;
    event = new PixelChamberEvent(ievent);
    printf("reading event %d\n",ievent);
    printf("nentries = %d\n",nentries);
    //
    for (int l=l_event1; l<nentries; l++)
    {
      ntuple->GetEvent(l);
      if (EventID == ievent){
        pt = new Datapoint(x, y, z, ic, jc, kc);
        event->AddDataPoint(pt);
        l_event=l;
      } else if(EventID==ievent+1) break;
    }

    printf("Cluster finding of event %d\n",ievent);
    event->findClusters();
    event->fitClusters(0,1);	//***************
    event->mergeCompatibleClusters();

    int id = event->GetEventID();
    std::cout << id << std::endl;

    //
    //
    event->fitClusters(0,1);//***************
    //      event->attachNoisePoints();//*************
    event->moveNoiseToPoints();//*************
    event->findClusters(2,4);
    event->moveNoiseToPoints();//*************
    event->findClusters(2,6);
    event->moveNoiseToPoints();//*************
    //event->ViewClusters();

    //      event->findClusters(2,8);

    //  event->fitClusters(0,1);
    TClonesArray *EvTracks = event->GetTrackletsObj();

    printf("fClusters entries = %d\n", EvTracks->GetEntries());
    for(int ntracks=0; ntracks<EvTracks->GetEntries(); ntracks++)
    {
      cluster = (Cluster*) EvTracks->At(ntracks);
      c0 = new Cluster(*cluster);
      new(ToWriteTemp[ToWriteTemp.GetEntries()])(Cluster)(*c0);
    }

    TTree_ID=TreeNumber;
    Ev_ID = event->GetEventID();
    Chi2Vert = event->GetChiVert()/event->GetNDFVert();
    ClustersTree->Fill();
    printf("EvID = %d, Ev_ID = %d\n",event->GetEventID(), Ev_ID);
    printf("TreeNumber = %d,   TTree_ID = %d\n",TreeNumber,   TTree_ID);
    Clusters->Clear("C");

    delete event;

    ProcInfo_t *info = new ProcInfo_t;
    gSystem->GetProcInfo(info);
    Long_t memo =info->fMemResident;
    printf ("Iteration %d -- Resident memory used by root: %ld\n", event->GetEventID(), memo);
    delete info;
  }
  fOut->Write();
  fOut->Close();
}
