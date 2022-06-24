#include "DMProcessAnnihilation.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"

#include "DMParticleAPrime.hh"

#include "DMParticleChi.hh"
#include "DMParticleChi1.hh"
#include "DMParticleChi2.hh"


#include "G4ProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"


#include "DarkMatterParametersFactory.hh"



DMProcessAnnihilation::DMProcessAnnihilation(DarkMatter* DarkMatterPointerIn, G4ParticleDefinition* theDMParticlePtrIn,
                                             G4double BiasSigmaFactorIn)
: G4VDiscreteProcess( "DMProcessAnnihilation", fUserDefined ),  // fElectromagnetic
  myDarkMatter(DarkMatterPointerIn),
  theDMParticlePtr(theDMParticlePtrIn),
  BiasSigmaFactor(BiasSigmaFactorIn),
  DMpar(0),
  iBranchingType(0),
  mChi(0),
  mChi1(0),
  mChi2(0)
{
  SetProcessSubType( 1 ); //fBremsstrahlung? // TODO: verify this

  DMpar = DarkMatterParametersFactory::GetInstance();
  if (DMpar){
      iBranchingType=(int)DMpar->GetRegisteredParam("BranchingType");

      if (iBranchingType==0){
          mChi=DMpar->GetRegisteredParam("RDM")*myDarkMatter->GetMA(); //A.C. by default, DarkMatter units are GeV
      }
      else{
          double r = DMpar->GetRegisteredParam("RDM", 1. / 3);
          mChi1 = myDarkMatter->GetMA() * r;
          mChi2 = (1. + DMpar->GetRegisteredParam("Ffactor")) * mChi1;
      }
  }
}

G4bool DMProcessAnnihilation::IsApplicable(const G4ParticleDefinition & pDef)
{
  return ("e+" == pDef.GetParticleName());
}

G4double DMProcessAnnihilation::GetMeanFreePath( const G4Track& aTrack,
                                                 G4double, /*previousStepSize*/
                                                 G4ForceCondition* /*condition*/ )
{
  G4double DensityMat = aTrack.GetMaterial()->GetDensity()/(g/cm3);
  G4double ekin = aTrack.GetKineticEnergy()/GeV;


  if( myDarkMatter->EmissionAllowed(ekin, DensityMat) ) {

    G4double CrossSection = myDarkMatter->GetSigmaTot(ekin); //A.C. by DarkMatter definition, this is in picobarn
    CrossSection *= picobarn;

    //The DarkMatter classes compute the cross section for eps = epsilBench. Here, we revert back to epsilon
    CrossSection *= (myDarkMatter->Getepsil()* myDarkMatter->Getepsil())/(myDarkMatter->GetepsilBench()* myDarkMatter->GetepsilBench());
    CrossSection /= myDarkMatter->GetSigmaNorm();


    G4double n = aTrack.GetMaterial()->GetTotNbOfAtomsPerVolume(); //The annihilation cross section already contains a multiplicative factor "Z".
    G4double XMeanFreePath = 1./(n*CrossSection);

    XMeanFreePath /= BiasSigmaFactor;


    return XMeanFreePath;

  }
  return DBL_MAX;
}

G4VParticleChange* DMProcessAnnihilation::PostStepDoIt( const G4Track& aTrack,
                                                        const G4Step & aStep )
{
  const G4double incidentE = aTrack.GetKineticEnergy();
  //const G4double DMMass = theDMParticleAPrimePtr->GetPDGMass();
  G4ThreeVector incidentDir = aTrack.GetMomentumDirection();

  //For the e+e- --> R --> ff process, the only relevant variable is the cosine of the angle of the f in the CM frame.

  //Do not simulate the decay.
  if(myDarkMatter->Decay()==0) {
      G4ThreeVector DMDirection=incidentDir;
      G4double DME = incidentE;
      G4double DMM = myDarkMatter->GetMA()*GeV; //A.C. Dark Matter units are, by default, GeV
      G4double DMKinE = incidentE - DMM;



      G4DynamicParticle* movingDM = new G4DynamicParticle( theDMParticlePtr,
                                                            DMDirection,
                                                            DMKinE );
       aParticleChange.Initialize( aTrack );

       // Set DM:
       aParticleChange.SetNumberOfSecondaries( 1 );
       aParticleChange.AddSecondary( movingDM );
       // Kill projectile:
       aParticleChange.ProposeEnergy( 0. );
       aParticleChange.ProposeTrackStatus( fStopAndKill ) ;

       std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding()
                 << " emitted by " << aTrack.GetDefinition()->GetParticleName()
                 << " with energy = " << incidentE/GeV << " GeV, DM energy = " << DME/GeV <<" GeV "<< std::endl;

       return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  else{ //simulate the decay e+e- -->A' -->ff

      if( myDarkMatter->GetDMType() == 1){ //dark photon

          //1: Get the cosine of the final state f in the CM frame.
          G4double DMeta_CM=myDarkMatter->SimulateEmissionResonant(incidentE);
          //2: Get the phi angle of the final state f in the CM frame
          G4double DMphi_CM=G4UniformRand()*2*CLHEP::pi;
          //3: Define the Lorentz Vector of the CM (e+ + e-)
          G4LorentzVector vCM(aTrack.GetMomentum(),aTrack.GetTotalEnergy()+CLHEP::electron_mass_c2);

          G4double ss=2*(incidentE)*CLHEP::electron_mass_c2+2*CLHEP::electron_mass_c2*CLHEP::electron_mass_c2;
          G4double Pcm,Px_cm,Py_cm,Pz_cm;
          switch(iBranchingType){
              case 0: //fermionic and scalar elastic
              case 1:
              {
                  G4double Ecm=sqrt(ss)/2;
                  Pcm=sqrt(Ecm*Ecm-mChi*mChi);

                  Px_cm=Pcm*sqrt(1-DMeta_CM*DMeta_CM)*sin(DMphi_CM);
                  Py_cm=Pcm*sqrt(1-DMeta_CM*DMeta_CM)*cos(DMphi_CM);
                  Pz_cm=Pcm*DMeta_CM;

                  //The two lorentz vectors in the CM frame
                  G4LorentzVector v1(Px_cm,Py_cm,Pz_cm,Ecm);
                  G4LorentzVector v2(-Px_cm,-Py_cm,-Pz_cm,Ecm);

                  v1.boost(vCM.boostVector());
                  v2.boost(vCM.boostVector());

                  G4DynamicParticle* movingDM1 = new G4DynamicParticle(DMParticleChi::Definition(),v1.vect());
                  G4DynamicParticle* movingDM2 = new G4DynamicParticle(DMParticleChi::Definition(),v2.vect());

                  aParticleChange.Initialize( aTrack );

                  // Set DM:
                  aParticleChange.SetNumberOfSecondaries( 2 );
                  aParticleChange.AddSecondary( movingDM1);
                  aParticleChange.AddSecondary( movingDM2);

                  // Kill projectile:
                  aParticleChange.ProposeEnergy( 0. );
                  aParticleChange.ProposeTrackStatus( fStopAndKill ) ;

                  std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding()
                                  << " emitted by " << aTrack.GetDefinition()->GetParticleName()
                                  << " with energy = " << incidentE/GeV << "GeV, DM energy = " << incidentE/GeV <<" GeV "<< std::endl;

                  std::cout <<"Daugther1: "<<DMParticleChi::Definition()->GetPDGEncoding()<<" energy= "<<v1.e()/GeV <<"GeV "<<std::endl;
                  std::cout <<"Daugther2: "<<DMParticleChi::Definition()->GetPDGEncoding()<<" energy= "<<v2.e()/GeV <<"GeV "<<std::endl;



                  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

                  break;
              }
              case 2: //asymmetric DM
              {

                  G4double Ecm1=(ss+mChi1*mChi1-mChi2*mChi2)/(2*sqrt(ss));
                  G4double Ecm2=(ss-mChi1*mChi1+mChi2*mChi2)/(2*sqrt(ss));

                  Pcm=sqrt(Ecm1*Ecm1-mChi1*mChi1);
                  Px_cm=Pcm*sqrt(1-DMeta_CM*DMeta_CM)*sin(DMphi_CM);
                  Py_cm=Pcm*sqrt(1-DMeta_CM*DMeta_CM)*cos(DMphi_CM);
                  Pz_cm=Pcm*DMeta_CM;

                  //The two lorentz vectors in the CM frame
                  G4LorentzVector v1(Px_cm,Py_cm,Pz_cm,Ecm1);
                  G4LorentzVector v2(-Px_cm,-Py_cm,-Pz_cm,Ecm2);

                  v1.boost(vCM.boostVector());
                  v2.boost(vCM.boostVector());

                  G4DynamicParticle* movingDM1 = new G4DynamicParticle(DMParticleChi1::Definition(),v1.vect());
                  G4DynamicParticle* movingDM2 = new G4DynamicParticle(DMParticleChi2::Definition(),v2.vect());

                  aParticleChange.Initialize( aTrack );

                  // Set DM:
                  aParticleChange.SetNumberOfSecondaries( 2 );
                  aParticleChange.AddSecondary( movingDM1);
                  aParticleChange.AddSecondary( movingDM2);

                  // Kill projectile:
                  aParticleChange.ProposeEnergy( 0. );
                  aParticleChange.ProposeTrackStatus( fStopAndKill ) ;

                  std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding()
                            << " emitted by " << aTrack.GetDefinition()->GetParticleName()
                            << " with energy = " << incidentE/GeV << "GeV, DM total energy = " <<incidentE/GeV << std::endl;

                  std::cout <<"Daugther1: "<<DMParticleChi1::Definition()->GetPDGEncoding()<<" energy= "<<v1.e()/GeV <<"GeV "<<std::endl;
                  std::cout <<"Daugther2: "<<DMParticleChi2::Definition()->GetPDGEncoding()<<" energy= "<<v2.e()/GeV <<"GeV "<<std::endl;


                  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

                  break;
              }
          }
      }//end DMtype==1
      else{
          G4cout<<"DMProcessAnnihilation::PostStepDoIt DMType== "<< myDarkMatter->GetDMType()<<G4endl;
          G4cout<<"The decay to final state LDM particles is not yet implemented "<<G4endl;
          //revert to Decay()==0a
          G4ThreeVector DMDirection=incidentDir;
          G4double DME = incidentE;
          G4double DMM = myDarkMatter->GetMA()*GeV;
          G4double DMKinE = incidentE - DMM;
          G4DynamicParticle* movingDM = new G4DynamicParticle( theDMParticlePtr,
                                                                      DMDirection,
                                                                      DMKinE );
          aParticleChange.Initialize( aTrack );

          // Set DM:
          aParticleChange.SetNumberOfSecondaries( 1 );
          aParticleChange.AddSecondary( movingDM );
          // Kill projectile:
          aParticleChange.ProposeEnergy( 0. );
          aParticleChange.ProposeTrackStatus( fStopAndKill ) ;

          std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding()
                    << " emitted by " << aTrack.GetDefinition()->GetParticleName()
                    << " with energy = " << incidentE/GeV << "GeV, DM total energy = " <<incidentE/GeV << std::endl;

          return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
  }
}
