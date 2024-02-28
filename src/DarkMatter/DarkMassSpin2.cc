#include <iostream>

#include "DarkMassSpin2.hh"
#include "Utils.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_integration.h>

#include <cmath>

#define LIM_ITER_QAGS 25000
#define REL_ERR_QAGS 1.0e-8
#define REL_ERR_QAGS_MIN 1.0e-2
#define ABS_ERR_QAGS 0.0
#define REL_ERR_QAGS_STEP 1.2


// Auxiliary structures and functions:


// Additional structure holding E0 value to be passed into GSL callbacks.
struct BoundParmsDs {
  DarkMassSpin2 * this_;
  double E0;
};


// Additional structure holding E0, xEv values to be passed into GSL callbacks.
struct BoundParmsDsDx {
  DarkMassSpin2 * this_;
  double E0;
  double xEv;
};


// A callback wrapping function for single differential cross-section as 
// function fraction of energy
static double _DarkMassSpin2DsDx(double x, void* parms_) {
  // Geting parameters
  BoundParmsDsDx* parms = reinterpret_cast< BoundParmsDsDx* >( parms_ );
  // Forward invocation to target method
  return parms->this_->CrossSectionDSDXDTheta( parms->xEv, x, parms->E0 );
}

// A callback wrapping function for total differential cross-section as 
// function energy
static double _DarkMassSpin2Ds(double x, void* parms_) {
  // Geting parameters
  BoundParmsDs* parms = reinterpret_cast< BoundParmsDs* >( parms_ );
  // Forward invocation to target method
  return parms->this_->CrossSectionDSDX( x, parms->E0 );
}


// Class methods:  -------------------------------------------------------------

DarkMassSpin2::DarkMassSpin2( double MAIn, double EThreshIn
                            , double SigmaNormIn, double ANuclIn
                            , double ZNuclIn, double DensityIn
                            , double epsilIn, int IDecayIn 
                            ) : DarkMatter( MAIn, EThreshIn, SigmaNormIn
                                          , ANuclIn, ZNuclIn, DensityIn
                                          , epsilIn, IDecayIn){
  // Setting id for type particle 
  DMType = 5;
  ParentPDGID = 11;
  DaughterPDGID = 11;

  PrepareVariables();

  mInit = Mel;          // mass of initial particle
  tMax =         100.;  // tmax initial; tmax = E0*E0 will be taken
  ThetaMax =      0.1;  // Max. angle of radiation
}


DarkMassSpin2::~DarkMassSpin2(  ){ ; }


// Double differential cross-section as a function of angle and energy fraction of radiated particle
double DarkMassSpin2::CrossSectionDSDXDTheta( double XEv, double Theta, double E0 )
{
  return CrossSectionDSDXDTheta_WW( XEv, Theta, E0 );
}


// Single differential cross-section as a function of energy fraction
double DarkMassSpin2::CrossSectionDSDX( double XEv, double E0 )
{
  return CrossSectionDSDX_WW( XEv, E0 );
}


// Total cross-section
double DarkMassSpin2::TotalCrossSectionCalc( double E0 ){
  return TotalCrossSectionCalc_WW( E0 );
}


// Double differential cross-section as a function of angle and fraction 
// of energy of emitted particle in WW approximation. Formula for the double 
// differential cross-section from 2210.00751 ( mass of lepton is zero )
double DarkMassSpin2::CrossSectionDSDXDTheta_WW( double XEv, double Theta, double E0 )
{
  // Checking correct the low limit of fraction of energy
  if( XEv*E0 <= MA ){ return 0.0; }
  // Setting momentums of photon flux.
  // Typical square of inverse of screening and nuclear radius
  double aa = 111.0*pow(ZNucl,-1.0/3.0)/Mel, d = 0.164*pow(ANucl,-2.0/3.0);
  // Transmitted momentum, GeV^2
  double ta = pow(1.0/aa, 2.0), td = d;                              
  // Expression for lepton Mandelstam variable u in WW approx. shifting on mass
  double Uxth = E0*E0 * Theta*Theta * XEv 
              + MA*MA * (1.0 - XEv)/XEv + mInit*mInit * XEv;
  // Limits for photon flux
  double tmin = Uxth*Uxth / ( 4.0 * E0*E0 * (1.0 - XEv)*(1.0 - XEv) )
       , tmax = tMax;
  // Checking correct limits
  if ( tmax < tmin ){ return 0.0; }
  // Mandelstam variable in WW approx.
  double tMnd = - XEv * Uxth / (1.0 - XEv) + MA*MA
       , sMnd =   Uxth / (1.0 - XEv) + mInit*mInit
       , uMnd = - Uxth + mInit*mInit;
  // Matrix element 
  double
  MatElwhoutMass = ((uMnd + tMnd)*(uMnd + tMnd) + (uMnd - MA*MA)*(uMnd - MA*MA)) 
                 * (4.0 * uMnd * (uMnd + tMnd - MA*MA) - MA*MA * tMnd)
                 / (4.0 * tMnd * (sMnd -  mInit*mInit) * (uMnd -  mInit*mInit));
  // Photon flux with using Tsai's form-factor
  double ChiWWAnalytical = 2.0*( td - ta );
         ChiWWAnalytical -= ( ta + td + 2.0*tmin ) 
                            * std::log( (td + tmin ) / ( ta + tmin ) );
         ChiWWAnalytical *= ZNucl*ZNucl * td*td / pow( (ta - td), 3 );
  // Calc result 
  double
  Prefactor = 4.0 * E0*E0 *alphaEW*alphaEW  * epsilBench*epsilBench 
            * sin(Theta) / (1.0 - XEv) * sqrt( XEv*XEv - ( MA*MA ) / ( E0*E0 ) )
            / ( 8.0 * M_PI * (sMnd -  mInit*mInit)*(sMnd -  mInit*mInit) );
  return Prefactor * ChiWWAnalytical * MatElwhoutMass;
}


// Single differential cross-section as a function of fraction of energy 
// of emitted particle in WW approximation. Uses GSL integration method QAGS
// and double differential cross-section
double DarkMassSpin2::CrossSectionDSDX_WW( double XEv, double E0 )
{
  // Checking correct the low limit of fraction of energy
  if( E0 < 2.0 * MA ) { return 0.0; }
  // Setting upper limit of number of subintervals and relative error
  double limInterv = LIM_ITER_QAGS, relerr = REL_ERR_QAGS, abserr = ABS_ERR_QAGS;
  // Allocation of space for the work of the integrator
  gsl_integration_workspace* w = gsl_integration_workspace_alloc( limInterv );
  double res, err, thMin = 0.0, thMax = ThetaMax;
  // Setting function of integrand with parameters
  gsl_function F;
  BoundParmsDsDx parms = { this, E0, XEv };
  F.function = &_DarkMassSpin2DsDx;
  F.params = &parms;
  // Starting integration with increase relative error in case fall
  int status = 1;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  while(status) {
    status = gsl_integration_qags( &F, thMin, thMax
                                 , abserr, relerr, limInterv, w, &res, &err );
    relerr *= REL_ERR_QAGS_STEP;
  }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free( w ); 
  return res;
}

// Total cross-section in WW approximation. Uses GSL integration method QAGS and
// single differential cross-section in pBarn
double DarkMassSpin2::TotalCrossSectionCalc_WW( double E0 )
{
  // Checking correct the low limit of fraction of energy
  if( E0 < 2.0 * MA ) { return 0.0; }
  // Setting upper limit of number of subintervals and relative error
  double limInterv = LIM_ITER_QAGS
       , relerr = REL_ERR_QAGS, abserr = ABS_ERR_QAGS;
  // Allocation of space for the work of the integrator
  gsl_integration_workspace* w = gsl_integration_workspace_alloc( limInterv );
  double res, err, xMin = MA/E0, xMax = 1.0 - mInit/E0;
  // Setting function of integrand with parameters
  gsl_function F;
  BoundParmsDs parms = { this, E0 };
  F.function = &_DarkMassSpin2Ds;
  F.params = &parms;
  // Starting integration with increase relative error in case fall
  int status = 1;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  while(status) {
    status = gsl_integration_qags( &F, xMin, xMax
                                 , abserr, relerr, limInterv, w, &res, &err );
    relerr *= REL_ERR_QAGS_STEP;
  }
  gsl_set_error_handler(old_handler);
  gsl_integration_workspace_free( w ); 
  return GeVtoPb * res;
}


double DarkMassSpin2::GetSigmaTot( double E0 )
{
  return GetSigmaTot0( E0 );
}


double DarkMassSpin2::CrossSectionDSDXDU( double XEv, double UThetaEv, double E0 )
{
  (void) XEv; (void) UThetaEv; (void) E0;
  return 0;
}


// Decay width into a pair of fermions
double DarkMassSpin2::Width()
{
  double rIn =  mInit / MA;
  return  ( 1.0/(160.0*M_PI) ) * MA*MA*MA * epsil*epsil
        * pow( 1.0 - 4.0 * rIn*rIn, 3.0/2.0 ) * ( 1.0 + (8.0/3.0) * rIn*rIn );
}
