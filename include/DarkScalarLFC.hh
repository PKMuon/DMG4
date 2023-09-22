class DarkScalarLFC : public DarkMatter
{

  public:

    DarkScalarLFC(double MAIn, double EThreshIn, double SigmaNormIn=1., double ANuclIn=207., double ZNuclIn=82., double DensityIn=11.35,
                double epsilIn=0.0001, int IDecayIn=0, int ParentPDGIDIn=-9999);

    virtual ~DarkScalarLFC();

    virtual double TotalCrossSectionCalc(double E0); // wrapper for computing total cs using functions below
    double TotalCrossSectionCalc_ETL(double E0);     // from dsdx IWW with KFactors from ETL
    double TotalCrossSectionCalc_IWW(double E0);     // from dsdx IWW
    double TotalCrossSectionCalc_WW(double E0);      // from dsdx WW
    double TotalCrossSectionCalc_WW2(double E0);     // from dsdxdtheta WW
    double TotalCrossSectionCalc_WW3(double E0);     // from dsdxdpsi WW
    virtual double GetSigmaTot(double E0);                                          // calls function GetSigmaTot0 in parent class
    virtual double CrossSectionDSDX(double XEv, double E0);                         // wrapper for CrossSectionDSDX_WW and CrossSectionDSDX_IWW
    virtual double CrossSectionDSDX_WW(double XEv, double E0);                      // 
    virtual double CrossSectionDSDX_IWW(double XEv, double E0);                     // not including photon flux

    virtual double CrossSectionDSDXDU(double XEv, double UThetaEv, double E0);      // UTheta is 0.5theta^2.
    double CrossSectionDSDXDU_WW(double XEv, double UThetaEv, double E0);
    double CrossSectionDSDXDU_IWW(double XEv, double UThetaEv, double E0);      // UTheta is 0.5theta^2. not including photon flux

    virtual double CrossSectionDSDXDPSI(double XEv, double auxpsi, double E0);      // wrapper for _WW and _IWW. 
    virtual double CrossSectionDSDXDPSI_IWW(double XEv, double auxpsi, double E0);  // not including photon flux and h^2 alpha^2 / 2pi
    virtual double CrossSectionDSDXDPSI_WW(double XEv, double auxpsi, double E0);

    virtual double CrossSectionDSDXDTheta(double XEv, double Theta, double E0);     // WW. this is dxdtheta, not dxdcostheta
    virtual double Width();
    double chiIWW(double E0);

    // for plotting purposes, not sampling
    double CrossSectionDSDX_WW_withprefactor(double XEv, double E0);                      // 
    double CrossSectionDSDX_IWW_withprefactor(double XEv, double E0);                     // not including flux and h^2 alpha^2 / 4pi
    double CrossSectionDSDPSI_IWW_noprefactor(double Psi, double E0);

    // setters and getters (CAREFUL, only use BEFORE calling PrepareTable)
    void SetIApprox(int IApproxIn) {IApprox = IApproxIn;}
    void SetIMethodTotalCS(int IMethodTotalCSIn) {IMethodTotalCS = IMethodTotalCSIn;}
    void SetThetaMax(double ThetaMaxIn) {ThetaMax = ThetaMaxIn;}
    int GetIApprox() {return IApprox;}
    int GetIMethodTotalCS() {return IMethodTotalCS;}
    double GetThetaMax() {return ThetaMax;}
    int GetParentPDGID() {return ParentPDGID;}
  private:

    int IApprox;
    int IMethodTotalCS;
    double tMax;
    double ThetaMax;
    double PsiMax;
    double MChild; // mass of scattered final state lepton
};
