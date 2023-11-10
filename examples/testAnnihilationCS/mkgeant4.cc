#include "globals.hh"

#include "G4ios.hh"

#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "DarkZAnnihilation.hh"

#include "Randomize.hh"

// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"

int main() {

  //G4double EThresh = 2.*MA; // for full shape
  //G4double EThresh = 35.; // for sensitivity calculations
  G4double EThresh = 1.; // for shape studies
  //G4double EThresh = 2000.; // to turn off A emissions

  // Set parameters to be passed to DarkMatterAnnihilation class
  G4double MA = 0.220;
  G4double coupling = 1e-3;
  G4double SigmaNorm = 1.;
  G4double ANuclPb = 207.;
  G4double ZNuclPb = 82.;
  G4double DensityPb = 11.35;
  G4int IDecayIn = 1;
  G4int IBranchingIn = 0;
  G4double rIn = 1./3.;
  G4double alphaD = 1./3.;

  // Initialize DarkMatterAnnihilation instance
  DarkMatterAnnihilation* myDarkMatter = new DarkZAnnihilation(MA, EThresh, SigmaNorm, ANuclPb, ZNuclPb, DensityPb, coupling, IDecayIn, rIn, alphaD, IBranchingIn);

  myDarkMatter->PrepareTable();

  // Set energy range and steps
  G4double Emax = 100.;
  G4double Emin = EThresh;
  unsigned int nSteps = 1000;
  G4double Ediff = (Emax - Emin)/nSteps;

  // Define output ROOT files with plots
  TFile* hOutputFile = new TFile("result.root", "RECREATE");
  TGraph *hSigma = new TGraph(nSteps);
  //TH1D *hSigmaAE = new TH1D("hSigmaAE","hSigma Atomic Effects", nSteps+1,Emin,Emax);

  G4double width = myDarkMatter->Width();

  G4cout << G4endl;
  G4cout << "Test of the DarkMatter package: Resonant Annihilation production simulation, coupling = " << coupling << ", mass = " << MA << " GeV" << G4endl;
  G4cout << G4endl;

  for(unsigned int i=0; i<nSteps; i++) {
    double ekin = Emin + Ediff * i;
    if (myDarkMatter->EmissionAllowed(ekin, DensityPb)) {
      double totalCS = myDarkMatter->GetSigmaTot(ekin);
      double preF = myDarkMatter->PreFactor(ekin);
      double BW = myDarkMatter->BreitWignerDenominator(ekin);
      //std::cout << ekin << ", prefactor: " << preF << ", 1/BW: " << 1./BW << ", total CS: " << totalCS << std::endl;
      //double totalCS_AE = myDarkMatter->GetSigmaTotAtomicEffectsOneShell(ekin);
      hSigma->SetPoint(i, ekin, totalCS);
      //hSigmaAE->Fill(totalCS_AE);
    }
  }

  G4cout << G4endl;
  G4cout << "Cross section in pb for eps=0.0001 cs = " << myDarkMatter->GetAccumulatedProbability() << G4endl;

  TCanvas *c=new TCanvas("c","c");
    
  hSigma->SetLineColor(1);
  hSigma->SetTitle(Form("Total CS for resonant annihilation production with #Gamma = %3.2e MeV", width*1.e3));
  hSigma->Draw("AC*");
  c->SetLogy();
  c->Write();
  //hSigmaAE->SetLineColor(2);
  //hSigmaAE->Draw("SAMES");

  hOutputFile->Write();

  return 0;
}
