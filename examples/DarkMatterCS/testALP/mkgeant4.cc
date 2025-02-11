#include <math.h>
#include "globals.hh"
#include "DarkMatter.hh"
#include "ALP.hh"
#include <iostream>
#include <fstream>
#include <string>


int main() {

  double E0, massTested, EThresh, sigmaNorm, ANucl, ZNucl, Dens, Epsil;
  int DecayType = 1;
  unsigned int nPoints = 100;

  EThresh = 0.12; // for total CS testing applicable to mass range mA<1MeV
  massTested = 0.1; 
  sigmaNorm = 1.e9;
  ANucl = 207;
  ZNucl = 82;
  Dens = 11.35;
  Epsil = 1.e-4;
  E0 = 50.;

  // Specify the file name
  std::string filename = "output.txt";

  // Open the file for writing (creates the file if it doesn't exist)
  std::ofstream outputFile(filename);

  // Check if the file is opened successfully
  if (!outputFile.is_open()) {
    std::cerr << "Error: Unable to open the file." << std::endl;
    return 1; // Return an error code
  }

  // Loop over number of points needed
  double XAcc, angles[2];
  for(int jj=0; jj<nPoints; jj++) {
    DarkMatter* myDarkMatter = new ALP(massTested, EThresh, sigmaNorm, ANucl, ZNucl, Dens, Epsil, DecayType);  // Initialize DM by default for Pb with eps=0.0001
    XAcc = myDarkMatter->SimulateEmissionWithAngle3(E0, angles);                                               // Get calculated ETL cs from DMG4       
    outputFile << XAcc << "," << angles[0] << "," << angles[1] << std::endl;
  }

  // Close the file
  outputFile.close();

  std::cout << "File '" << filename << "' has been written successfully." << std::endl;


  return 0;
}
