#define Mel 5.109989461E-04 // electron mass in GeV
#define Mmu 0.1056583745 // muon mass in GeV
#define Mtau 1.77686     // tau mass in GeV
#define alphaEW 1./137.
#define APb 207
#define ZPb 82
#define rhoPb 11.35
#define GeV2pb 2.56819E-09


std::map<TString, TH1F*> getHistograms(std::string filename="output.txt") {

  ifstream inFile(filename.c_str());
  if (!inFile.is_open()) {
    throw std::runtime_error("Error: Unable to open the file.");
  }

  std::map<TString, TH1F*> histograms;
  histograms["Fractional Energy"] = new TH1F("X", "Histogram of fractional energy; X [-]; nevts [-]", 100, 0.901, 1.001);
  histograms["#theta Angle"]      = new TH1F("theta", "Histogram of #theta angle; Angle [rad]; nevts [-]", 200, 0., 0.01);
  histograms["#phi Angle"]        = new TH1F("phi", "Histogram of #phi angle; Angle [rad]; nevts [-]", 100, 0., 2*TMath::Pi());

  float X, theta, phi;
  char comma;
  while (inFile >> X >> comma >> theta >> comma >> phi) {
    histograms["Fractional Energy"]->Fill(X);
    histograms["#theta Angle"]->Fill(theta);
    histograms["#phi Angle"]->Fill(phi);
  }
  inFile.close();

  for (const auto& [key, histo] : histograms){
    histo->SetLineWidth(2);
    histo->SetMarkerStyle(20);
  }

  return histograms;
}


double formFactor(double t) {
  double a = 111*TMath::Power(ZPb, -1./3.)/Mel;
  double a2 = a*a;
  double d = 0.164*TMath::Power(APb, -2./3.);
  return ZPb * ((a2*t) / (1.+a2*t)) * (1./(1.+t/d));
}

double diffCS(double theta, double g=1., double DMMass=0.1, double Egamma=50.) {
  double E2 = Egamma*Egamma;
  double t = E2*theta*theta + 0.25*DMMass*DMMass*DMMass*DMMass/E2;
  double ft = formFactor(t);
  double delta = DMMass*DMMass/(2*E2);
  return alphaEW*g*g/4.*ft*ft*theta*theta*theta/(theta*theta+delta*delta)/(theta*theta+delta*delta)*GeV2pb;
}

double totCS(double g=1., double DMMass=0.1, double Egamma=50.) {
  double E2 = Egamma*Egamma;
  double a = 111*TMath::Power(ZPb, -1./3.)/Mel;
  double a2 = a*a;
  double d = 0.164*TMath::Power(APb, -2./3.);
  double tmin = DMMass*DMMass*DMMass*DMMass/(4*E2);
  return alphaEW*g*g/8.*ZPb*ZPb*(TMath::Log(d/(1/a2+tmin))-2)*GeV2pb;
}

double thetaMAX(double DMMass=0.1, double Egamma=50.) {
  double E2 = Egamma*Egamma;
  double a = 111*TMath::Power(ZPb, -1./3.)/Mel;
  double a2 = a*a;
  double tmin = DMMass*DMMass*DMMass*DMMass/(4*E2);
  return 1/(a*Egamma)*TMath::Sqrt(3*(1+a2*tmin));
}



void plotDiffCS() {

  // Define histograms using data from output.txt
  std::map<TString, TH1F*> histograms = getHistograms("output.txt");

  // Take values used in mkgeant4!
  double g = 1.e-4;
  double DMMass = 0.1; // GeV
  double Egamma = 50.; // GeV

  // Define function that describes diff CS
  auto fun = new TF1("diffcs", Form("diffCS(x, %f, %f, %f)", g, DMMass, Egamma), 0., 1.e-2);
  auto maximumCS = diffCS(thetaMAX());

  // Print some information on the CS
  std::cout << "Value of theta_max = " << thetaMAX() << std::endl;
  std::cout << "Cross-section value at theta_max = " << maximumCS << std::endl;
  std::cout << "Total Cross-section = " << totCS(g, DMMass, Egamma) << std::endl;

  // Rescale histogram to represent pb units
  // Scale with totCS using the parameters from mkgeant.cc and the bin width
  auto factor = totCS(g, DMMass, Egamma)/histograms["#theta Angle"]->Integral()/histograms["#theta Angle"]->GetBinWidth(1);
  std::cout << factor << std::endl;
  histograms["#theta Angle"]->Scale(factor);

  // Create Canvas to plot all related histograms
  TCanvas* canvas = new TCanvas("canvas", "ALP distributions", 1600,900);
  canvas->Divide(3,1);

  // Draw Fractional Energy plot
  canvas->cd(1);
  histograms["Fractional Energy"]->Draw("E");
  gPad->SetLogy();
  gPad->SetLogx();

  // Draw Theta angle distribution
  canvas->cd(2);
  histograms["#theta Angle"]->Draw("E");
  histograms["#theta Angle"]->GetYaxis()->SetTitleOffset(1.5);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  fun->Draw("same");
  gPad->SetLogy();
  //gPad->SetLogx();

  auto legend = new TLegend();
  legend->AddEntry(histograms["#theta Angle"], "Sample from DMG4");
  legend->AddEntry(fun, "Reference PDF");
  legend->Draw();

  // Draw Phi angle distribution
  canvas->cd(3);
  histograms["#phi Angle"]->Draw();
  gPad->SetLogy();

  // Save .pdf
  canvas->SaveAs("ALP.pdf");

}
