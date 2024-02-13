#pragma once

#include "DarkMatter.hh"

#include <map>
#include <vector>

/*
 * DarkMatterAnnihilation.hh
 *
 *  Created on: Sep, 2023
 *      Author: bbantoob
 */

class DarkMatterAnnihilation : public DarkMatter
{
  public:


    DarkMatterAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn = 1., double ANuclIn = 207., double ZNuclIn = 82., double DensityIn = 11.35, double epsilIn = 0.0001, int IDecayIn = 0,double rIn=1./3,double alphaD=0.5,int IBranchingIn=0, double fIn=0.1);

    virtual ~DarkMatterAnnihilation()=0;

    // methods to be implemented in each derived class
    virtual double PreFactor(double s) = 0; //s is the e+e- invariant mass
    virtual double Width() = 0;
    virtual double AngularDistributionResonant(double eta,double E0) = 0;

    virtual bool EmissionAllowed(double E0, double DensityMat);// E0 in GeV, density in g/cm3
    virtual double sMin(); //minimum value for s allowing the reaction
    virtual double BreitWignerDenominator(double E0);

    virtual double q(double s); //return the magnitude of the 3-momentum of the daughter particles in the Cm frame


    //for resonant production e+ e- --> R --> f f, this function returns the cosine of the angle of the f in the CM frame.
    double SimulateEmissionResonant(double E0); //E0 in GeV
    double TotalCrossSectionCalc(double E0);
    double GetTotalCrossSectionMax();

    double GetSigmaTotAtomicEffects(double E0,const std::map<int, int>& Z,const std::map<int, std::vector<double>>& ene);
    double GetTotalCrossSectionMaxAtomicEffects(const std::map<int,int>& Z,const std::map<int,std::vector<double>>& ene);

    double GetSigmaTotAtomicEffectsOneShell(double E0,int Zshell,const std::vector<double>& eneShell);
    double GetSigmaTotAtomicEffectsOneShellFull(double E0,int Zshell,const std::vector<double>& eneShell);

  protected:

    int iBranchingType;
    double r,f;
    double alphaD;
    double mChi;
    double mChi1,mChi2,deltaMchi;

    double funIntegral(long double x,long double eta2,long double delta2);
};

