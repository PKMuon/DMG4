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
  G4double MA = 0.25;
  G4double coupling = 1e-5;
  G4double SigmaNorm = 1.;
  G4double ANuclPb = 207.;
  G4double ZNuclPb = 82.;
  G4double DensityPb = 11.35;
  G4int IDecayIn = 1;
  G4int IBranchingIn = 0; // 0-9: Lmu-Ltau, 10-19: B-L
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
  double ekin;

  // Define output ROOT files with plots
  TFile* hOutputFile = new TFile("result.root", "RECREATE");
  TGraph *gSigma = new TGraph(nSteps);
  TGraph *gPreFactor = new TGraph(nSteps);
  TGraph *gBW = new TGraph(nSteps);
  TH1D *hAngle = new TH1D("hAngle","Angular distribution; #eta; nevts [-]", 100,-1,1);
  //TH1D *gSigmaAE = new TH1D("gSigmaAE","gSigma Atomic Effects", nSteps+1,Emin,Emax);

  G4double width = myDarkMatter->Width();

  // ---------------------------------------------------------------
  // Test of total cross-section as a function of the primary energy
  // ---------------------------------------------------------------

  G4cout << G4endl;
  G4cout << "Test of the DarkMatter package: Resonant Annihilation production simulation, coupling = " << coupling << ", mass = " << MA << " GeV" << G4endl;
  G4cout << G4endl;

  for(unsigned int i=0; i<nSteps; i++) {
    ekin = Emin + Ediff * i;
    if (myDarkMatter->EmissionAllowed(ekin, DensityPb)) {
      double totalCS = myDarkMatter->GetSigmaTot(ekin);
      double preF = myDarkMatter->PreFactor(ekin);
      double BW = myDarkMatter->BreitWignerDenominator(ekin);
      //std::cout << Form("%3.2f GeV:   Prefactor = %5.2e [pb]  --  1/BW = %3.2e [GeV^4]  --  total CS = %5.2e [pb]\n", ekin, preF, 1./BW, totalCS);
      //double totalCS_AE = myDarkMatter->GetSigmaTotAtomicEffectsOneShell(ekin);
      gSigma->SetPoint(i, ekin, totalCS);
      gPreFactor->SetPoint(i, ekin, preF);
      gBW->SetPoint(i, ekin, BW);
      //gSigmaAE->Fill(totalCS_AE);
    }
  }

  // ---------------------------------------------------------
  // Test sampling of cross-section at peak E0 = MA*MA / 2*m_e
  // ---------------------------------------------------------

  ekin = MA*MA/(2*5.11e-4);

  G4cout << "Testing sampling for resonant annihilation at E = " << ekin << " GeV, for coupling = " << coupling << ", mass = " << MA << " GeV" << G4endl;

  int NTry=1000;
  int ITry;
  double angle;
  for(int i=0; i<NTry; i++) {

    ITry = myDarkMatter->Emission(ekin, DensityPb, 1.);
    angle = myDarkMatter->SimulateEmissionResonant(ekin);
    hAngle->Fill(angle);

    G4cout << "Emission simulated, Theta = " << angle << G4endl;
    G4cout << angle << G4endl;

  }
  (void)ITry; // to avoid warning


  G4cout << G4endl;
  G4cout << "Cross section in pb for eps=0.0001 cs = " << myDarkMatter->GetAccumulatedProbability() << G4endl;

  // ------------------------------
  // Save plots to ROOT output file
  // ------------------------------

  TCanvas *c=new TCanvas("c","c");

  // Save CS graph
  c->SetName("totalCS");
  c->SetGrid();
  gSigma->SetTitle(Form("Total CS for MA = %3.2e MeV and coupling = %3.2e with #Gamma = %3.2e MeV", MA*1.e3, coupling, width*1.e3));
  gSigma->SetMarkerStyle(21);
  gSigma->SetMarkerColor(kBlue);
  gSigma->Draw("ACP");
  gSigma->GetXaxis()->SetTitle("E_{primary} [GeV]");
  gSigma->GetYaxis()->SetTitle("Cross-section [pb]");
  c->SetLogy();
  c->Write();

  // Save prefactor graph
  c->Clear();
  c->SetName("CSprefactor");
  c->SetGrid();
  gPreFactor->SetTitle(Form("Prefactor of resonant CS for MA = %3.2e MeV and coupling = %3.2e with #Gamma = %3.2e MeV", MA*1.e3, coupling, width*1.e3));
  gPreFactor->SetMarkerStyle(21);
  gPreFactor->SetMarkerColor(kBlue);
  gPreFactor->Draw("ACP");
  gPreFactor->GetXaxis()->SetTitle("E_{primary} [GeV]");
  gPreFactor->GetYaxis()->SetTitle("CS prefactor [pb*GeV^{4}]");
  c->SetLogy();
  c->Write();

  // Save Breit-Wigner denominator graph
  c->Clear();
  c->SetName("BWdenominator");
  c->SetGrid();
  gBW->SetTitle(Form("Breit-Wigner denominator for resonant annihilation production with #Gamma = %3.2e MeV", width*1.e3));
  gBW->SetMarkerStyle(21);
  gBW->SetMarkerColor(kBlue);
  gBW->Draw("ACP");
  gBW->GetXaxis()->SetTitle("E_{primary} [GeV]");
  gBW->GetYaxis()->SetTitle("Denominator [GeV^{-4}]");
  c->SetLogy();
  c->Write();

  // Save angular distribution histogram
  c->Clear();
  c->SetName("angle");
  c->SetGrid();
  hAngle->Draw("hist");
  c->SetLogy();
  c->Write();


  hOutputFile->Write();

  return 0;
}
