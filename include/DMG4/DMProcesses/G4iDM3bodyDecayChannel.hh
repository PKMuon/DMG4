#ifndef G4iDM3bodyDecayChannel_h
#define G4iDM3bodyDecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4iDM3bodyDecayChannel :public G4VDecayChannel
{
 // Class Description
 //   This class describes kinematics of the three body decay X2->X1e+e-


 public: // With Description
  //Constructors 
  G4iDM3bodyDecayChannel(const G4String& theParentName,
                           G4double theBR,
                           const G4String& theDaughterName,
                           const G4String& theLeptonName,
                           const G4String& theAntiLeptonName);
  //  Destructor
  virtual ~G4iDM3bodyDecayChannel();

  protected:
    // Copy constructor and assignment operator
      G4iDM3bodyDecayChannel(const G4iDM3bodyDecayChannel &);
      G4iDM3bodyDecayChannel & operator=(const G4iDM3bodyDecayChannel &);

  private:
      double DiffWidth(double, double, double, double, double, double);
      G4iDM3bodyDecayChannel();

  public: // With Description
     virtual G4DecayProducts *DecayIt(G4double);     

  private:
     enum{idDaughter=0, idLepton=1, idAntiLepton=2}; 

};  


#endif
