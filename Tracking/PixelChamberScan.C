// ConsoleApplication2.cpp : Defines the entry point for the console application.
//
#include <math.h>
#include <stdio.h>  
#include <math.h>
#include <TEfficiency.h>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPad.h"
#include "cluster.h"
#include "TClonesArray.h"

void PixelChamberScan(int myevent=0)
{
  //  gObjectTable->Print();
	//
  PixelChamberEvent *event[1000];
  // TObjArray *events = new TObjArray(5000);
  // events->SetOwner();

  
  //
  TFile *file1 = new TFile("B2_Digit_5000ev.root");
  TNtuple *ntupleDigit = (TNtuple*) file1->Get("Interest_D"); 
  //
  
  //Takes Ntuple
  int EventID_D, ic, jc, kc, ParticlePDGCode;
  double x, y, z, TrackID_1, TrackID_2, TrackID_3;
 
  //
  ntupleDigit->SetBranchAddress("EventID_D",&EventID_D); 	
  ntupleDigit->SetBranchAddress("i",&ic); 
  ntupleDigit->SetBranchAddress("j",&jc); 
  ntupleDigit->SetBranchAddress("k",&kc);
  ntupleDigit->SetBranchAddress("ParticlePDGCode_D",&ParticlePDGCode); 
  ntupleDigit->SetBranchAddress("x",&x); 
  ntupleDigit->SetBranchAddress("y",&y); 
  ntupleDigit->SetBranchAddress("z",&z); 
  ntupleDigit->SetBranchAddress("TrackID_1",&TrackID_1); 
  ntupleDigit->SetBranchAddress("TrackID_2",&TrackID_2); 
  ntupleDigit->SetBranchAddress("TrackID_3",&TrackID_3); 
  int nentries = ntupleDigit->GetEntries();
  TH1F *h_cluster_n[4];
  int nMC[4]={0};
  
  for (int i=myevent; i<myevent+1; i++) 
    {
      h_cluster_n[i] = new TH1F(Form("h_cluster_n%d",i),"cluster_n",50,0.,49.);
      
      for (int k=0; k<nentries; k++) 
  	{ 
  	  ntupleDigit->GetEvent(k); 	 
	  
  	  if(EventID_D==i)
  	    {
  	      h_cluster_n[i]->Fill(TrackID_1);
  	    }
  	}
    }
  for (int i=myevent; i<myevent+1; i++) 
    {
      for (int j = 0; j <= h_cluster_n[i]->GetNbinsX(); j++)
  	{
  	  if (h_cluster_n[i]->GetBinContent(j))
  	    nMC[i]++;
  	}
      std::cout << "Event: " << i << " Number of MC clusters = "<< nMC[i] << std::endl;
    }
  
  
	
  // reads firsts 4 events
  //PixelChamberEvent *event;
  Datapoint *pt;
    for (int m = myevent; m<myevent+1; m++)
    {
		// New event
      printf("reading event %d\n",m);
      event[m] = new PixelChamberEvent(m);
      //event = new PixelChamberEvent(m);
      printf("nentries = %d\n",nentries);
      //
      
      for (int l=0; l<nentries; l++)				
	{
	  //
	  ntupleDigit->GetEntry(l); 
	  //printf("EventID_D = %d\n",EventID_D);
	  //
	  if (EventID_D == m){	
	    printf("adding point n. %d x=%f y=%f z=%f i=%d j=%d k=%d\n",l,x,y,z,ic,jc,kc);
	    pt = new Datapoint(x, y, z, ic, jc, kc);
	    event[m]->AddDataPoint(pt);
	    //event->AddDataPoint(pt);
	    //events->Add(event);
	  } else if(EventID_D==m+1) break;
	}
      
    }
    printf("sono qui\n");
    TClonesArray *p = event[myevent]->getPointsObj();
    //    p->SetOwner();
    printf("There are %d  data points in event %d\n", p->GetEntries(), event[myevent]->getEventID());
    //   	find clusters, fit & merge
    printf("start reconstruction\n");
	for (int ievent=myevent; ievent<myevent+1; ievent++)
	{
		//
		printf("Cluster finding of event %d\n",ievent);
		event[ievent]->findClusters();
		//event = (PixelChamberEvent*)(events->At(ievent));
		//event->findClusters();
	  
		int id = event[ievent]->getEventID();
		//int id = event->getEventID();
		std::cout << id << std::endl;
		//		
		event[ievent]->fitClusters();	
		//event->fitClusters();	
		// //
		event[ievent]->mergeCompatibleClusters();
		//event->mergeCompatibleClusters();

		event[ievent]->fitClusters(1);	
		//event->fitClusters();	
		
		event[ievent]->CheckNoisePoints();
		event[ievent]->mergeCompatibleClusters();
		event[ievent]->fitClusters(0);	

		// event[ievent]->CheckNoisePoints();
		//event->CheckNoisePoints();

		//		event[ievent]->fitClusters();	
		//		event->fitClusters();	

		event[ievent]->ViewClusters();
		//event->ViewClusters();

	}

	//
	//	If vis is needed
	  //	   event = (PixelChamberEvent*)(events->At(3));
	//	   event[myevent]->ViewClusters();

    // TObjArray *track = (TObjArray*)(event->getTrackletsObj());
    // Cluster *cluster = (Cluster*)(track->At(2));
    //cluster->ViewClusterFit();
	
	//
	//delete events;
}

