#include "globals.hh"

#include "G4ios.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkScalars.hh"
#include "ALP.hh"
#include "DarkLFCScalars.hh"

#include "Randomize.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLorentzVector.h"


int main() {

  G4double MA = Mtau;// - 0.7*Mmu; // GeV, for mu N -> tau N + phi
  std::cout << "Mphi: " << MA*1000 << " MeV" << std::endl;
  //  G4double MA = Mmu - Mel; // GeV, for mu N -> e N + phi
  G4double SigmaNorm = 1.;

  //G4double EThresh = 2.*MA; // for full shape
  //G4double EThresh = 35.; // for sensitivity calculations
  G4double EThresh = 1.; // for shape studies
  //G4double EThresh = 2000.; // to turn off A emissions

  DarkLFCScalars* myDarkMatter = new DarkLFCScalars(MA, EThresh, SigmaNorm); // Initialize by default for Pb with eps=0.0001
  myDarkMatter->PrepareTable();

  double ekin = 160;// GeV
  double angles[2]; // [0]: theta, [1]: phi

  TLorentzVector pMu(ekin, 0, 0, ekin); // incoming muon 4-mom
  
  G4cout << G4endl;
  G4cout << "Test of the DarkMatter package: DM emission simulation, energy = " << ekin << " GeV, mass = " << MA << " GeV" << G4endl;
  G4cout << G4endl;


  TFile * output = new TFile("output.root", "RECREATE");
  TH1F * hEphi = new TH1F("hEphi", "E_{#phi}", 100, 0, 200);
  TH1F * hEtau = new TH1F("hEtau", "E_{#tau}", 100, 0, 200);
  TH1F * hEtauL = new TH1F("hEtauL", "E_{#tau L}", 100, 0, 200);
  TH1F * hThetaphi = new TH1F("hThetaphi", "#Theta_{#phi} Lab", 100, 0, 0.15); 
  
  int NTry=10000;
  int ITry;
  for(int i=0; i<NTry; i++) {

    ITry = myDarkMatter->Emission(ekin, 11.35, 1.);
    //double XAcc = myDarkMatter->SimulateEmission(ekin, angles);            // used in invisible mode for electrons
    //double XAcc = myDarkMatter->SimulateEmissionWithAngle(ekin, angles);   // one-step sampling, for electrons, not used by default
    //    double XAcc = myDarkMatter->SimulateEmissionWithAngle2(ekin, angles);    // two-step sampling, used by default for electrons if decays are enabled
    //double XAcc = myDarkMatter->SimulateEmissionByMuon(ekin, angles);      // one-step sampling (2-dim), for muons
    //double XAcc = myDarkMatter->SimulateEmissionByMuon2(ekin, angles);     // two-step sampling, used by default for muons
    double XAcc = myDarkMatter->SimulateEmissionVector(ekin, angles);    // Vector sampling using DSDXDTheta
    if(XAcc > 0.0000001) {
      G4cout << "Emission simulated, X = " << XAcc << " Theta = " << angles[0] << G4endl;
    }
    double theta = angles[0];
    double phi = angles[1];
    // Calculate Tau 4-mom from p_tau = p_mu - p_phi (approximation)
    double pPhiMom = TMath::Sqrt(XAcc*ekin*XAcc*ekin - MA*MA);
    TLorentzVector pPhi(XAcc*ekin, pPhiMom*TMath::Sin(theta)*TMath::Cos(phi), pPhiMom*TMath::Sin(theta)*TMath::Sin(phi), pPhiMom*TMath::Cos(theta) );
    TLorentzVector pTau = pMu - pPhi;
    hEtau->Fill((1-XAcc)*ekin);
    hEtauL->Fill(pTau.E());
    hEphi->Fill(XAcc*ekin);
    hThetaphi->Fill(theta);

  }
  (void)ITry; // to avoid warning

  hEphi->Write();
  hThetaphi->Write();
  hEtau->Write();
  hEtauL->Write();
  output->Close();
  
  G4cout << G4endl;
  G4cout << "Cross section in pb for eps=0.0001 cs = " << myDarkMatter->GetAccumulatedProbability()/((double)NTry) << G4endl;

  return 0;
}
