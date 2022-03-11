#include "DarkMatterPhysicsMessenger.hh"
#include "CommandsFactory.hh"

#include "DarkMatterPhysics.hh"

// STD Library
#include <iostream>

DarkMatterPhysicsMessenger::DarkMatterPhysicsMessenger(std::string Name)
  : fName(Name)
{
  DarkMatterCmd = new CommandsFactory("DarkMatterCmd");
  DarkMatterCmd->SetInputFile(fName);

  DarkMatterCmd->SetCmdDirectory("/DarkMatter/");
  DarkMatterCmd->AddNewCmdWithADouble("/DarkMatter/BiasSigmaFactor0", 1.e8);
  DarkMatterCmd->AddNewCmdWithADouble("/DarkMatter/DMMass", 0.0167);
  DarkMatterCmd->AddNewCmdWithADouble("/DarkMatter/EThresh", 35.0);
  DarkMatterCmd->AddNewCmdWithAString("/DarkMatter/DMParticle", "DarkPhotons");
  DarkMatterCmd->AddNewCmdWithADouble("/DarkMatter/SigmaNorm", 1.);
  DarkMatterCmd->AddNewCmdWithADouble("/DarkMatter/ANucl", 207.0);
  DarkMatterCmd->AddNewCmdWithADouble("/DarkMatter/ZNucl", 82.0);
  DarkMatterCmd->AddNewCmdWithADouble("/DarkMatter/Density", 11.35);
  DarkMatterCmd->AddNewCmdWithADouble("/DarkMatter/Epsilon", 0.0001);
  DarkMatterCmd->AddNewCmdWithADouble("/DarkMatter/Decay", 0);

  // gather the commands values
  DarkMatterCmd->ParseInputFile();
}

DarkMatterPhysicsMessenger::~DarkMatterPhysicsMessenger()
{
  delete DarkMatterCmd;
}

void DarkMatterPhysicsMessenger::SetPhysicsList(DarkMatterPhysics* PhysList)
{
  fPhysList = PhysList;
}

void DarkMatterPhysicsMessenger::ApplyCmd()
{
  fPhysList->SetBiasSigmaFactor0( 
      DarkMatterCmd->GetCmdValueWithADouble("/DarkMatter/BiasSigmaFactor0") 
      );
  fPhysList->SetEThresh( 
      DarkMatterCmd->GetCmdValueWithADouble("/DarkMatter/EThresh") 
      );
  fPhysList->SetDMMass( 
      DarkMatterCmd->GetCmdValueWithADouble("/DarkMatter/DMMass") 
      );
  fPhysList->SetDMParticle( 
      DarkMatterCmd->GetCmdValueWithAString("/DarkMatter/DMParticle") 
      );
  fPhysList->SetSigmaNorm( 
      DarkMatterCmd->GetCmdValueWithADouble("/DarkMatter/SigmaNorm") 
      );
  fPhysList->SetANucl( 
      DarkMatterCmd->GetCmdValueWithADouble("/DarkMatter/ANucl") 
      );
  fPhysList->SetZNucl( 
      DarkMatterCmd->GetCmdValueWithADouble("/DarkMatter/ZNucl") 
      );
  fPhysList->SetDensity( 
      DarkMatterCmd->GetCmdValueWithADouble("/DarkMatter/Density") 
      );
  fPhysList->SetEpsilon( 
      DarkMatterCmd->GetCmdValueWithADouble("/DarkMatter/Epsilon") 
      );
  fPhysList->SetDecay( 
      DarkMatterCmd->GetCmdValueWithADouble("/DarkMatter/Decay") 
      );
}

