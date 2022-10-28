// -----------------------------------------------------------------------
//      Implementation of the inelastic DM 3-body decay channel X2->X1l+l-
//      based on the G4DalitzDecayChannel.cc class
//      
//      Author: Martina Mongillo 
//      Date: 9 August 2022
// -----------------------------------------------------------------------

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4iDM3bodyDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4RandomDirection.hh"

G4iDM3bodyDecayChannel::G4iDM3bodyDecayChannel()
  :G4VDecayChannel()
{
}

G4iDM3bodyDecayChannel::G4iDM3bodyDecayChannel(
   const G4String& theParentName,
   G4double theBR,
   const G4String& theDaughterName,
   const G4String& theLeptonName,
   const G4String& theAntiLeptonName)
   :G4VDecayChannel("iDM Decay",1)
{
  // set names for daughter particles
  SetParent(theParentName);
  SetBR(theBR);
  SetNumberOfDaughters(3);
  SetDaughter(idDaughter, theDaughterName);
  SetDaughter(idLepton, theLeptonName);
  SetDaughter(idAntiLepton, theAntiLeptonName);
}

G4iDM3bodyDecayChannel::~G4iDM3bodyDecayChannel()
{
}

G4iDM3bodyDecayChannel::G4iDM3bodyDecayChannel(const G4iDM3bodyDecayChannel &right):
  G4VDecayChannel(right)
{
}

G4iDM3bodyDecayChannel & G4iDM3bodyDecayChannel::operator=(const G4iDM3bodyDecayChannel & right)
{
  if (this != &right) { 
    kinematics_name = right.kinematics_name;
    verboseLevel = right.verboseLevel;
    rbranch = right.rbranch;

    // copy parent name
    parent_name = new G4String(*right.parent_name);

    // clear daughters_name array
    ClearDaughtersName();

    // recreate array
    numberOfDaughters = right.numberOfDaughters;
    if ( numberOfDaughters >0 ) {
      if (daughters_name !=0) ClearDaughtersName();
      daughters_name = new G4String*[numberOfDaughters];
      //copy daughters name
      for (G4int index=0; index < numberOfDaughters; index++) {
          daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
  }
  return *this;
}


//Define spin-averaged matrix element of the X2->X1l+l- 3-body decay
//==================================================================
//s1=(p_lep1+p_lep2)^2 -> invariant mass squared l+l- system
//s2=(p_lep2+p_chi1)^2 -> invariant mass squared X1l+ system

double G4iDM3bodyDecayChannel::DiffWidth(double s1, double s2, double mchi1, double mchi2, double ma, double mlep)
{
  //Define Dalitz plot boundaries (kinematically allowed region)
  double bound_mid=pow(mchi1,2.)+pow(mlep,2.)+(pow(mchi2,2.)-s1-pow(mchi1,2.))/2.;
  double bound_var=(std::sqrt(pow(s1,2.)+pow(mchi2,4.)+pow(mchi1,4.)-2.*s1*pow(mchi2,2.)-2.*s1*pow(mchi1,2.)-2.*pow(mchi2,2.)*pow(mchi1,2.))*std::sqrt(pow(s1,2.)-4.*s1*pow(mlep,2.)))/(2.*s1);
  double upper_boundary=bound_mid+bound_var;
  double lower_boundary=bound_mid-bound_var;
  //Define matrix element
  double denominator=pow((s1-ma*ma),2.);
  double numerator=-(s1*s1-pow(mchi2,2.)*(s1+2.*s2-2.*pow(mchi1,2.))+2.*mchi1*mchi2*(s1+2.*pow(mlep,2.))+s1*(2.*s2-pow(mchi1,2.))+2.*(pow(s2,2.)-s2*(pow(mchi1,2.)+2.*pow(mlep,2.))+pow(mlep,4.)));
  //Set negative matrix element outside boundary 
  if (denominator == 0 || s2 > upper_boundary || s2 < lower_boundary) return -1.;
  return numerator/denominator;
}


G4DecayProducts *G4iDM3bodyDecayChannel::DecayIt(G4double) 
{
#ifdef G4VERBOSE 
  if (GetVerboseLevel()>1) G4cout << "G4iDM3bodyDecayChannel::DecayIt ";
#endif 
  CheckAndFillParent();
  CheckAndFillDaughters();

  //Get masses of involved particles
  G4double parentmass = G4MT_parent->GetPDGMass();
  G4double leptonmass = G4MT_daughters[idLepton]->GetPDGMass();
  G4double daughtermass= G4MT_daughters[idDaughter]->GetPDGMass();
  G4ParticleDefinition * aprime = particletable->FindParticle(5500322);
  G4double aprimemass = aprime->GetPDGMass();
  
  
  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle(G4MT_parent, dummy, 0.0);
 

  //Matrix element sampling
  //=======================

  //Bounds of the sampled variables                                                                                                                      
  G4double s1min  = 4*pow(leptonmass,2.);
  G4double s1max  = pow(parentmass,2.)+pow(daughtermass,2.)-2*parentmass*daughtermass;
  G4double s2min  = pow(daughtermass,2.)+pow(leptonmass,2.)+2*daughtermass*leptonmass;
  G4double s2max  = pow(parentmass,2.)+pow(leptonmass,2.)-2*parentmass*leptonmass;

  //Maximum of the distribution is at (4*mlep^2,1/2(m_chi1^2+m_chi2^2-2mlep^2)) for s1min<s1<s1max, m_chi2>0, m_chi2/2<m_chi1<m_chi2, 0<mlep<1/2(m_chi2-m_chi1), ma>m_chi2+m_chi1
  G4double peak_ycoord = (pow(daughtermass,2.)+pow(parentmass,2.)-2*pow(leptonmass,2.))/2.;

  //Sampling of the distribution to get s1,s2
  G4double M_val, M_max, s1, s2, f;
  do {
    s1 = G4UniformRand()*(s1max-s1min) + s1min;
    s2 = G4UniformRand()*(s2max-s2min) + s2min;
    M_max=DiffWidth(s1min, peak_ycoord, daughtermass, parentmass, aprimemass, leptonmass);
    M_val=G4UniformRand()*M_max;
    f=DiffWidth(s1,s2, daughtermass, parentmass, aprimemass, leptonmass);
  } while (M_val > f); 
  
  //Kinematics calculation
  //======================

  //Calculate X1 momentum in CMS frame (restframe of X2)
  G4double Pchi1 = G4PhaseSpaceDecayChannel::Pmx(parentmass, daughtermass, std::sqrt(s1));
  G4double Echi1 = std::sqrt(Pchi1*Pchi1 + daughtermass*daughtermass);
  G4double Ekinchi1=Echi1-daughtermass;
  //Sample direction (isotropic)
  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta)); 
  G4double phi  = twopi*G4UniformRand()*rad;
  G4ThreeVector X1direction(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);
  
  //Create G4DynamicParticle for X1                                                                       
  G4DynamicParticle * chi1particle = new G4DynamicParticle(G4MT_daughters[idDaughter] , X1direction, Ekinchi1);
  G4LorentzVector pchi1_CMS = chi1particle->Get4Momentum();
 
  //Calcurate beta of (l+ l-) system - boost of Jackson frame in which l+ and l- are back to back (rest frame of A'*)
  G4double beta = -1.0*Pchi1/(parentmass-Echi1);
 
  //Boost X1 to this frame
  G4LorentzVector pchi1_JF = pchi1_CMS;
  pchi1_JF.boost(-1.0*X1direction.x()*beta, -1.0*X1direction.y()*beta, -1.0*X1direction.z()*beta);
  G4double Echi1_JF = pchi1_JF.e();
  G4ThreeVector pchi1_JF_3v=pchi1_JF.vect();
  G4ThreeVector X1direction_JF=pchi1_JF_3v.unit();

  //Calculate momenta of leptons in the rest frame of (l+ l-) system
  G4double Plepton=G4PhaseSpaceDecayChannel::Pmx(std::sqrt(s1),leptonmass,leptonmass);
  G4double Elepton = std::sqrt(Plepton*Plepton + leptonmass*leptonmass);
  
  //Calculate angle between X1 and lep in the JF frame, get random orthogonal axis to the plane, rotate to get lep direction
  G4double alpha=std::acos((-s2+pow(daughtermass,2.)+pow(leptonmass,2.)+2*Elepton*Echi1_JF)/(2.0*Plepton*pchi1_JF_3v.mag()));
  G4ThreeVector rand=G4RandomDirection();
  G4ThreeVector rot_ax=rand-(rand.dot(X1direction_JF))*X1direction_JF; 
  rot_ax=rot_ax.unit();
  G4ThreeVector ldirection_JF=X1direction_JF;
  ldirection_JF.rotate(alpha,rot_ax);
  
  //Create G4DynamicParticle for leptons  in the rest frame of (l+ l-)system                                                 
  G4DynamicParticle * leptonparticle = new G4DynamicParticle(G4MT_daughters[idLepton], -1.0*ldirection_JF, Elepton-leptonmass);
  G4DynamicParticle * antileptonparticle= new G4DynamicParticle(G4MT_daughters[idAntiLepton], ldirection_JF, Elepton-leptonmass);
  
  //Boost leptons to CMS frame (rest frame of the parent X2) 
  G4LorentzVector pl1 = leptonparticle->Get4Momentum();
  pl1.boost(X1direction.x()*beta, X1direction.y()*beta, X1direction.z()*beta);
  leptonparticle->Set4Momentum(pl1);
  G4LorentzVector pl2 = antileptonparticle->Get4Momentum();
  pl2.boost(X1direction.x()*beta, X1direction.y()*beta, X1direction.z()*beta);
  antileptonparticle->Set4Momentum(pl2);

  //Create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;
  products->PushProducts(chi1particle);
  products->PushProducts(leptonparticle);
  products->PushProducts(antileptonparticle);
  
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "G4iDM3bodyDecayChannel::DecayIt ";
     G4cout << "  create decay products in rest frame " <<G4endl;
     products->DumpInfo();
  }
#endif
  return products;
}





