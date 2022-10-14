#include "globals.hh"

#include "G4ios.hh"

#include "DarkMatter.hh"
#include "DarkMassSpin2.hh"


#include "Randomize.hh"


int main() {

  G4double MA = 0.017;
  G4double SigmaNorm = 1.;

  //G4double EThresh = 2.*MA; // for full shape
  //G4double EThresh = 35.; // for sensitivity calculations
  G4double EThresh = 1.; // for shape studies
  //G4double EThresh = 2000.; // to turn off A emissions

  double E0Tab[] = {1,2,5,10,15,25,50,150,200}; 

  double TCSWM[] = { 10.43451525586200, 15.69935095150541, 24.75787553520542
                   , 33.26810768234646, 38.90819192662305, 46.6922349817397
                   , 58.42391896338646, 79.65344183434661, 85.72892252870289 };

  DarkMatter* myDarkMatter = new DarkMassSpin2( MA, EThresh, SigmaNorm ); 
  std::cout << "E0" << "\t" << "TCSDMG4" << "\t" << "TCSWM" <<  std::endl;
  for(int i = 0; i < 9; i++){
    std::cout << E0Tab[i] << "\t" 
              << myDarkMatter->TotalCrossSectionCalc(E0Tab[i]) <<  "\t"
              << TCSWM[i]
              <<  std::endl;
  }
  return 0;
}
