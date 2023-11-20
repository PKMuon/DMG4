#include "DarkMatterPhysics.hh"
#include "DarkMatterParametersFactory.hh"

#include "G4SystemOfUnits.hh"


// BiasSigmaFactor Invisible mode Vector EThresh=35
// 900.  9.e12
// 16.7  1.e9
//  5.   1.75e8
//  2.   3.5e7
// 0.5   1.2e7
// 0.1   3.6e6
// 0.01  1.55e6
// 0.002 1.38e6

// BiasSigmaFactor Invisible mode Scalar EThresh=35
// 16.7  2.3e9

// BiasSigmaFactor Visible mode Vector EThresh=18
// 16.7  3.4e8



//A.C. the default version does nothing with ptr
/*
 *
 * SYSTEM OF UNITS
 *
 * All entered quantities must be accompained by their unit of measurement.
 *
 * (ANucl is the atomic number)
 *
 */
bool DarkMatterPhysics::DarkMatterPhysicsConfigure(void *ptr)
{
    //call an instance of the class
     DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();

     DMpar->RegisterNewParam("BiasSigmaFactor0", 1.e9);
     DMpar->RegisterNewParam("EThresh", 35.*GeV); // for sensitivity calculations invisible mode
     //G4double EThresh = 18.; // for sensitivity calculations visible mode
     //G4double EThresh = 1.; // for shape studies
     //G4double EThresh = 2000.; // to turn off A emissions

     //select particle type and details
     DMpar->RegisterNewParam("DMProcessType", 1.); // 1 - 4: Brem. process for Vector, Scalar, Axial, Pseudoscalar, 5 - spin 2, 21 - ALP
                                                   // 31 - ZPrime (muon beams), 11 - 14: Annihilation
     DMpar->RegisterNewParam("DMMass", 0.0167*GeV);
     DMpar->RegisterNewParam("Epsilon", 0.0001);

     // Initialize for Pb
     DMpar->RegisterNewParam("ANucl"      ,207.   );
     DMpar->RegisterNewParam("ZNucl"      ,82.    );
     DMpar->RegisterNewParam("Density"    ,11.35 *(g/cm3) );

     // Initialize for W
   //  DMpar->RegisterNewParam("ANucl"   ,184.   );
   //  DMpar->RegisterNewParam("ZNucl"   ,74.    );
   //  DMpar->RegisterNewParam("Density" ,19.25 *(g/cm3)s );

     /* Comments to "DecayType"
      0: invisible (no decays simulated), 1: visible, 2: visible with constraints.
      *
      For annihilation, the process reads e+ e- -> R -> ff, where "R" is the intermediate resonance, and ff the final state particles.
      Setting this to "0" means that the code accounts for the total cross-section for e+e- --> R --> ff, but only produces the R in the final state,
      with a fully invisible energy (the missing energy is properly accounted for, since the e+ is killed).
      Setting this to "1" means that the code will produce the "ff" final state, according to Branching type.
     */

     DMpar->RegisterNewParam("DecayType", 0.);

     /* Comments to "BranchingType"
      * Optional to make different decay table; 0 : invisible decays or SM lepton decays, depending on DecayType; 1 : B - L model; 2 : Semivisible: Inelastic DM; 3 : Semivisible: Dirac Inelastic DM    
      * For annihilation e+e- --> R --> ff the meaning is DIFFERENT!

       CASE Dark Photon, Dark Scalar, Dark Axial Vector, Dark Pseudo Scalar

       0: default: fermionic ff final state
       1: scalar ff final state
       2: asymmetric fermionic DM final state produced (iDM)
       3: Dirac DM final state produced (i2DM)

       CASE Dark Z'

       0: Lmu-Ltau vanilla model,  nu-nu final state - sum nuMu + nuTau (default)
       1: Lmu-Ltau DM model,       scalar DM final state
       10: B-L vanilla model,      nu-nu final state
       11: B-L DM model,           scalar DM final state
      */
     //DMpar->RegisterNewParam("BranchingType", 0.);

     /* additional parameters for annihilation (if absent the default ones will be used)
      * RDM -> mChi/mA ratio
      * AlphaD -> value of alphaDark
      * */
     //DMpar->RegisterNewParam("RDM", 1./3.);
     //DMpar->RegisterNewParam("AlphaD", 0.5);

     /* additional parameters for semivisible DM, in addition to above parameters for annihilation (if absent the default ones will be used)*/
     //DMpar->RegisterNewParam("Ffactor", 0.4);

     /* additional parameter for Dirac DM, describing the dark fermion mixing (if absent the default ones will be used)*/
     //DMpar->RegisterNewParam("IDMTheta", 1.e-3);

     /* Additional parameters to handle narrow width resonances in annihilation
      * dEmaxPerStep -> the maximum energy loss per step in the material(s) where the annihilation is allowed
      */
     //DMpar->RegisterNewParam("dEmaxPerStep",5*MeV)

  return true;
}
