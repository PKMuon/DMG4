#include "globals.hh"

#include <iostream>

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "ALP.hh"
#include "DarkZ.hh"

#include "Randomize.hh"


int main() {

  double MA = 0.017;
  double SigmaNorm = 1.;

  //double EThresh = 2.*MA; // for full shape
  //double EThresh = 35.; // for sensitivity calculations
  double EThresh = 1.; // for shape studies
  //double EThresh = 2000.; // to turn off A emissions

  DarkMatter* myDarkMatter = new DarkPhotons(MA, EThresh, SigmaNorm); // Initialize by default for Pb with eps=0.0001
  //DarkMatter* myDarkMatter = new DarkZ(MA, EThresh, SigmaNorm); // Initialize by default for Pb with eps=0.0001
  //DarkMatter* myDarkMatter = new ALP(MA, EThresh, SigmaNorm); // Initialize by default for Pb with eps=0.0001
  myDarkMatter->PrepareTable();

  double ekin = 100.;
  double angles[2];

  std::cout << std::endl;
  std::cout << "Test of the DarkMatter package: DM emission simulation, energy = " << ekin << " GeV, mass = " << MA << " GeV" << std::endl;
  std::cout << std::endl;

  std::ofstream fout("result.d");

  int NTry=100;
  int ITry;
  for(int i=0; i<NTry; i++) {

    ITry = myDarkMatter->Emission(ekin, 11.35, 1.);
    //double XAcc = myDarkMatter->SimulateEmission(ekin, angles);            // used in invisible mode for electrons
    //double XAcc = myDarkMatter->SimulateEmissionWithAngle(ekin, angles);   // one-step sampling, for electrons, not used by default
    double XAcc = myDarkMatter->SimulateEmissionWithAngle2(ekin, angles);    // two-step sampling, used by default for electrons if decays are enabled
    //double XAcc = myDarkMatter->SimulateEmissionWithAngle3(ekin, angles);  // one-step sampling of Theta, for ALP
    //double XAcc = myDarkMatter->SimulateEmissionByMuon(ekin, angles);      // one-step sampling (2-dim), for muons
    //double XAcc = myDarkMatter->SimulateEmissionByMuon2(ekin, angles);     // two-step sampling, used by default for muons

    if(XAcc > 0.0000001) {
      std::cout << "Emission simulated, X = " << XAcc << " Theta = " << angles[0] << std::endl;
      fout << XAcc << " " << angles[0] << std::endl;
    }

  }
  (void)ITry; // to avoid warning

  std::cout << std::endl;
  std::cout << "Cross section in pb for eps=0.0001 cs = " << myDarkMatter->GetAccumulatedProbability()/((double)NTry) << std::endl;

  return 0;
}
