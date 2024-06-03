#include "globals.hh"

#include "G4ios.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "DarkZAnnihilation.hh"
#include "AnnihilationStepLimiter.hh"

#include "Randomize.hh"

// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"


//An exponential model
G4double GetOneRandomEleEnergy(G4double B){
  G4double u=G4UniformRand();
  G4double x=-B*log(u);
  return x;
}

std::vector<G4double> SimulateElectronEnergies(const G4Element *elm,G4int is){

  std::vector<G4double> v;
  G4double B=elm->GetAtomicShell(is); //binding energy, positive, in g4 units

  const G4int nEne=100;
  for (G4int ii=0;ii<nEne;ii++){
    G4double x=GetOneRandomEleEnergy(B); //in G4 units
    v.push_back(x/GeV); //convert it in GeV
  }
  return v;
}


int main() {

  //G4double EThresh = 2.*MA; // for full shape
  //G4double EThresh = 35.; // for sensitivity calculations
  G4double EThresh = 1.; // for shape studies
  //G4double EThresh = 2000.; // to turn off A emissions

  // Set parameters to be passed to DarkMatterAnnihilation class
  G4double MA = 0.25; // GeV
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
  auto *dmLimiterProc=new AnnihilationStepLimiter(myDarkMatter,"StepLimiterAnnihilation");
  dmLimiterProc->SetFactor(2);

  myDarkMatter->PrepareTable();

  // Set energy range and steps
  G4double Emax = 100.;
  G4double Emin = EThresh;
  G4double me = CLHEP::electron_mass_c2/GeV;
  G4double Eres=(MA*MA-2*me*me)/(2*me);
  unsigned int nSteps = 1000;
  G4double Ediff = (Emax - Emin)/nSteps;
  double etot;

  // Define output ROOT files with plots
  TFile* hOutputFile = new TFile("result.root", "RECREATE");
  TGraph *gSigma = new TGraph(nSteps);
  TGraph *gSigmaAE = new TGraph(nSteps);
  TGraph *gSigmaAEFull = new TGraph(nSteps);
  TGraph *gPreFactor = new TGraph(nSteps);
  TGraph *gBW = new TGraph(nSteps);
  TGraph *gEnergyLoss = new TGraph(nSteps);
  TH1D *hAngle = new TH1D("hAngle","Angular distribution; #eta; nevts [-]", 100,-1,1);

  G4double width = myDarkMatter->Width()*GeV;

  // ---------------------------------------------------------------
  // Test of total cross-section as a function of the primary energy
  // ---------------------------------------------------------------

  G4cout << G4endl;
  G4cout << "Test of the DarkMatter package: Resonant Annihilation production simulation, coupling = " << coupling << ", mass = " << MA << " GeV" << G4endl;
  G4cout << G4endl;

  G4NistManager* nistManager = G4NistManager::Instance();

  // Retrieve the lead element from the NIST database
  const G4Element* elm = nistManager->FindOrBuildElement("Pb");
  G4int Z=(G4int)(elm->GetZ());

  // Set up maps for atomic effects calculations
  std::map<G4int, std::map<G4int,std::vector<G4double> > > shellElectronEnergies;
  std::map<G4int, std::map<G4int,G4int> > shellElectronZ;
  /* Compute the electron energies for this material*/
  if (shellElectronEnergies.find(Z)==shellElectronEnergies.end()){
    for (G4int is=0;is<elm->GetNbOfAtomicShells();is++){
      auto v=SimulateElectronEnergies(elm,is);
      shellElectronZ[Z][is]=elm->GetNbOfShellElectrons(is);
      shellElectronEnergies[Z][is]=v;
    }
  }

  for(unsigned int i=0; i<nSteps; i++) {
    etot = Emin + Ediff * i;
    if (myDarkMatter->EmissionAllowed(etot, DensityPb)) {
      double preF = myDarkMatter->PreFactor(etot);
      double BW = myDarkMatter->BreitWignerDenominator(etot);
      double totalCS = myDarkMatter->GetSigmaTot(etot);

      // Calculate CS with atomic effects
      double totalCSAtomicEffects = 0.;
      double totalCSAtomicEffectsFull = 0.;
      for (int is = 0; is < shellElectronZ[Z].size(); is++) {
        int ZeleShell = shellElectronZ[Z].at(is);
        const std::vector<double>& eneShell = shellElectronEnergies[Z].at(is);

        double sigmaShell = myDarkMatter->GetSigmaTotAtomicEffectsOneShell(etot, ZeleShell, eneShell);
        double sigmaShellFull = myDarkMatter->GetSigmaTotAtomicEffectsOneShellFull(etot, ZeleShell, eneShell);

        totalCSAtomicEffects += sigmaShell;
        totalCSAtomicEffectsFull += sigmaShellFull;
      }
      double energyLoss = dmLimiterProc->GetMaxEloss(etot*GeV)/GeV;

      //G4cout << Form("E = %3.2f GeV:   Prefactor = %5.2e [pb]  --  1/BW = %3.2e [GeV^4]  --  total CS = %5.2e [pb]", etot, preF, 1./BW, totalCS) << G4endl;
      //G4cout << Form("E = %3.2f [GeV]  --  Energy Loss = %3.2e [GeV]  --  Eres = %3.2e [GeV]", etot, energyLoss, Eres) << G4endl;

      // Save info in TGraph
      gSigma->SetPoint(i, etot, totalCS);
      gPreFactor->SetPoint(i, etot, preF);
      gBW->SetPoint(i, etot, BW);
      double ediff = (etot==Eres)? 1 : energyLoss/fabs(Eres-etot);
      gEnergyLoss->SetPoint(i, etot, ediff);
      gSigmaAE->SetPoint(i, etot, totalCSAtomicEffects);
      gSigmaAEFull->SetPoint(i, etot, totalCSAtomicEffectsFull);
    }
  }

  // ---------------------------------------------------------
  // Test sampling of cross-section at peak E0 = MA*MA / 2*m_e
  // ---------------------------------------------------------

  G4cout << "Testing sampling for resonant annihilation at E = " << Eres << " GeV, for coupling = " << coupling << ", mass = " << MA << " GeV" << G4endl;

  int NTry=1000;
  int ITry;
  double angle;
  for(int i=0; i<NTry; i++) {

    ITry = myDarkMatter->Emission(Eres, DensityPb, 1.);
    angle = myDarkMatter->SimulateEmissionResonant(Eres);
    hAngle->Fill(angle);

    G4cout << "Emission simulated, Theta = " << angle << G4endl;

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
  TLegend* legend = new TLegend();
  gSigma->SetTitle(Form("Total CS for MA = %3.2e MeV and coupling = %3.2e with #Gamma = %3.2e MeV", MA*1.e3, coupling, width));
  gSigma->SetMarkerStyle(21);
  gSigma->SetMarkerColor(kBlack);
  legend->AddEntry(gSigma, "CS without atomic effects", "lp");
  gSigma->Draw("ACP");
  gSigma->GetXaxis()->SetTitle("E_{primary} [GeV]");
  gSigma->GetYaxis()->SetTitle("Cross-section [pb]");
  gSigmaAE->SetTitle(Form("Total CS (with atomic effects for MA = %3.2e MeV and coupling = %3.2e with #Gamma = %3.2e MeV", MA*1.e3, coupling, width));
  legend->AddEntry(gSigmaAE, "CS with atomic effects in BW denominator", "lp");
  gSigmaAE->SetMarkerStyle(21);
  gSigmaAE->SetMarkerColor(kBlue);
  gSigmaAE->Draw("SAME CP");
  legend->AddEntry(gSigmaAEFull, "CS with atomic effects in full expression", "lp");
  gSigmaAEFull->SetMarkerStyle(21);
  gSigmaAEFull->SetMarkerColor(kGreen);
  gSigmaAEFull->Draw("SAME CP");
  legend->Draw();
  c->Update();
  c->SetLogy();
  c->Write();

  // Save prefactor graph
  c->Clear();
  c->SetName("CSprefactor");
  c->SetGrid();
  gPreFactor->SetTitle(Form("Prefactor of resonant CS for MA = %3.2e MeV and coupling = %3.2e with #Gamma = %3.2e MeV", MA*1.e3, coupling, width));
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
  gBW->SetTitle(Form("Breit-Wigner denominator for resonant annihilation production with #Gamma = %3.2e MeV", width));
  gBW->SetMarkerStyle(21);
  gBW->SetMarkerColor(kBlue);
  gBW->Draw("ACP");
  gBW->GetXaxis()->SetTitle("E_{primary} [GeV]");
  gBW->GetYaxis()->SetTitle("Denominator [GeV^{-4}]");
  c->SetLogy();
  c->Write();

  // Save Energy Loss graph
  c->Clear();
  c->SetName("EnergyLoss");
  c->SetGrid();
  gEnergyLoss->SetTitle(Form("Max energy loss from step limiter #frac{#deltaE_{max}}{E-E_{res}} for E_{res} = %3.2e GeV and #Gamma = %3.2e MeV", Eres*1.e3, width));
  gEnergyLoss->SetMarkerStyle(21);
  gEnergyLoss->SetMarkerColor(kBlue);
  gEnergyLoss->Draw("ACP");
  gEnergyLoss->GetXaxis()->SetTitle("E_{primary} [GeV]");
  gEnergyLoss->GetYaxis()->SetTitle("Diff");
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
