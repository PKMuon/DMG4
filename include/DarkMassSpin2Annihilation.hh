class DarkMatter;

class DarkMassSpin2Annihilation : public DarkMatter
{

  public:

    DarkMassSpin2Annihilation( double MAIn, double EThreshIn
                             , double SigmaNormIn = 1.0, double ANuclIn = 207.0
                             , double ZNuclIn = 82.0, double DensityIn = 11.35
                             , double epsilIn = 0.0001, int IDecayIn = 0
                             , double rIn = 1.0/3.0, double alphaD = 0.5
                             , int IBranchingIn = 0, double fIn = 0.1 );
    virtual ~DarkMassSpin2Annihilation();

    virtual double TotalCrossSectionCalc(double E0);
    virtual double GetSigmaTot(double E0);
    virtual bool EmissionAllowed(double E0, double DensityMat); 
    virtual double CrossSectionDSDX(double Xev, double E0);
    virtual double CrossSectionDSDXDU(double Xev, double UThetaEv, double E0);
    virtual double Width();
    //virtual void SetMA(double MAIn);

  private:

         int iBranchingType;

         double r,f;
         double alphaD;
         double mChi;
         double mChi1, mChi2, deltaMchi;
};
