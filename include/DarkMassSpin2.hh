#pragma once

#include "DarkMatter.hh"

class DarkMassSpin2 : public DarkMatter {
  public:
    
    DarkMassSpin2( double MAIn, double EThreshIn, double SigmaNormIn=1.
                 , double ANuclIn=207., double ZNuclIn=82.
                 , double DensityIn=11.35, double epsilIn=0.0001
                 , int IDecayIn=0 );

    virtual ~DarkMassSpin2();

    virtual double TotalCrossSectionCalc( double E0 ); 
    virtual double GetSigmaTot( double E0 );
    virtual double CrossSectionDSDX( double XEv, double E0 );
    virtual double CrossSectionDSDXDTheta( double XEv, double Theta, double E0 );
    virtual double CrossSectionDSDXDU(double XEv, double UThetaEv, double E0); 
    virtual double Width();

    double TotalCrossSectionCalc_WW( double E0 );
    double CrossSectionDSDX_WW( double XEv, double E0 );
    double CrossSectionDSDXDTheta_WW( double XEv, double Theta, double E0 );

  private:

    //int IApprox;
    //int IMethodTotalCS;
    double tMax;
    double ThetaMax;
    double mInit;
    //double PsiMax;
};
