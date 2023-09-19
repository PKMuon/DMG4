void histresults()
{
  double xval, yval;
  ifstream infile;

  TFile* hOutputFile = new TFile("result.root", "RECREATE");

  TH1D X("X", "X of process", 100, 1., 0.);               // lower limit > upper limit : choose limits automatically 
  TH1D Theta("Theta", "Theta of process", 100, 1., 0.);   // lower limit > upper limit : choose limits automatically

  // read file
  infile.open("result.d");
  int nread=0;
  if (infile.is_open()) {
    while (!infile.eof()) {
      infile >> xval >> yval;
      if(!infile.eof()) {
        nread++;
        X.Fill(xval);
        Theta.Fill(yval);
      }
    }
  }
  infile.close();
  cout << "read from file: " << nread << endl;
  hOutputFile->Write();
}
