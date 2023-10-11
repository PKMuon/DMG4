/*
 * DarkMatterAnnihilation.hh
 *
 *  Created on: Sep, 2023
 *      Author: bbantoob
 */

#ifndef INCLUDE_DARKMATTERANNIHILATION_HH_
#define INCLUDE_DARKMATTERANNIHILATION_HH_

class DarkMatter;

class DarkMatterAnnihilation : public DarkMatter {
  friend class DarkMassSpin2Annihilation;
  friend class DarkPhotonsAnnihilation;
  friend class DarkScalarsAnnihilation;
  friend class DarkPseudoScalarsAnnihilation;
  friend class DarkAxialsAnnihilation;

  public:


    DarkMatterAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn = 1., double ANuclIn = 207., double ZNuclIn = 82., double DensityIn = 11.35, double epsilIn = 0.0001, int IDecayIn = 0,double rIn=1./3,double alphaD=0.5,int IBranchingIn=0, double fIn=0.1);

    virtual ~DarkMatterAnnihilation()=0;

    // methods to be implemented in each derived class
    virtual double PreFactor(double E0) = 0;
    virtual double Width() = 0;
    virtual double AngularDistributionResonant(double eta,double E0) = 0;

    //for resonant production e+ e- --> R --> f f, this function returns the cosine of the angle of the f in the CM frame.
    double SimulateEmissionResonant(double E0); //E0 in GeV
    double TotalCrossSectionCalc(double E0);
    double GetTotalCrossSectionMax();

  private:

    int iBranchingType;
    double r,f;
    double alphaD;
    double mChi;
    double mChi1,mChi2,deltaMchi;
};

#endif /* INCLUDE_DARKMATTERANNIHILATION_HH_ */
