#include <TString.h>
#include <TROOT.h>

void loadLibs(){

  gROOT->LoadMacro("/Volumes/AliceM_SSD/AllEnv/EnvNewCode/PixelChamber/Datapoint.cxx+");
  gROOT->LoadMacro("/Volumes/AliceM_SSD/AllEnv/EnvNewCode/PixelChamber/SumDistance2.cxx+");
  gROOT->LoadMacro("/Volumes/AliceM_SSD/AllEnv/EnvNewCode/PixelChamber/FindChi2.cxx+");
  gROOT->LoadMacro("/Volumes/AliceM_SSD/AllEnv/EnvNewCode/PixelChamber/Cluster.cxx+");
  gROOT->LoadMacro("/Volumes/AliceM_SSD/AllEnv/EnvNewCode/PixelChamber/PixelChamberEvent.cxx+");

}
