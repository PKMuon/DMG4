#include "DarkMatterParametersFactory.hh"
#include "DarkMatterParametersRegistry.hh"

// STD Library
#include <iostream>
#include <utility>
#include <stdlib.h>

DarkMatterParametersFactory* DarkMatterParametersFactory::fFactoryInstance = nullptr;

DarkMatterParametersFactory::DarkMatterParametersFactory(std::string factoryName)
{
  fFactoryName = factoryName;
  std::cout << std::endl;
  std::cout << "WARNING: the name DarkMatterParametersFactory is obsolete, will be removed soon. Please use DarkMatterParametersRegistry instead"
            << std::endl;
  std::cout << std::endl;
}

DarkMatterParametersFactory::~DarkMatterParametersFactory()
{
  // ...
}

void DarkMatterParametersFactory::RegisterNewParam(std::string name, double value)
{
  std::cout << std::endl;
  std::cout << "WARNING: the name DarkMatterParametersFactory is obsolete, will be removed soon. Please use DarkMatterParametersRegistry instead"
            << std::endl;
  std::cout << std::endl;
  DarkMatterParametersRegistry* DMreg = DarkMatterParametersRegistry::GetInstance();
  DMreg->RegisterNewParam(name, value);

#if 0
  std::pair<std::map<std::string, double>::iterator, bool> insert 
    = fFactoryParams.insert(std::make_pair(name, value));
  if (!insert.second) {
    std::cerr << "Name  " << name << "  already booked, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
#endif
}

void DarkMatterParametersFactory::UpdateRegisteredParam(std::string name, double& value)
{
  std::cout << std::endl;
  std::cout << "WARNING: the name DarkMatterParametersFactory is obsolete, will be removed soon. Please use DarkMatterParametersRegistry instead"
            << std::endl;
  std::cout << std::endl;
  DarkMatterParametersRegistry* DMreg = DarkMatterParametersRegistry::GetInstance();
  DMreg->UpdateRegisteredParam(name, value);

#if 0
  std::map<std::string, double>::iterator it = fFactoryParams.find(name);
  if (it != fFactoryParams.end()) {
    it->second = value;
  }
  else {
    std::cerr << "Name " << name << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
#endif
}

double DarkMatterParametersFactory::GetRegisteredParam(std::string name)
{
  std::cout << std::endl;
  std::cout << "WARNING: the name DarkMatterParametersFactory is obsolete, will be removed soon. Please use DarkMatterParametersRegistry instead"
            << std::endl;
  std::cout << std::endl;
  DarkMatterParametersRegistry* DMreg = DarkMatterParametersRegistry::GetInstance();
  return DMreg->GetRegisteredParam(name);

#if 0
  std::map<std::string, double>::iterator it = fFactoryParams.find(name);
  if (it != fFactoryParams.end()) {
    return it->second;
  }
  else {
    std::cerr << "Name " << name << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
#endif
}

double DarkMatterParametersFactory::GetRegisteredParam(std::string name, double def)
{
  std::cout << std::endl;
  std::cout << "WARNING: the name DarkMatterParametersFactory is obsolete, will be removed soon. Please use DarkMatterParametersRegistry instead"
            << std::endl;
  std::cout << std::endl;
  DarkMatterParametersRegistry* DMreg = DarkMatterParametersRegistry::GetInstance();
  return DMreg->GetRegisteredParam(name, def);

#if 0
  std::map<std::string, double>::iterator it = fFactoryParams.find(name);
  if (it != fFactoryParams.end()) {
    return it->second;
  }
  else {
    std::cerr << "Name " << name << " not found, using default parameter:" << def << "\n";
    RegisterNewParam(name, def);
    return def;
  }
#endif
}

void DarkMatterParametersFactory::PrintListOfRegisteredParams()
{
  DarkMatterParametersRegistry* DMreg = DarkMatterParametersRegistry::GetInstance();
  DMreg->PrintListOfRegisteredParams();

#if 0
  std::map<std::string, double>::iterator it;
  for (it = fFactoryParams.begin(); it != fFactoryParams.end(); it++) {
    std::cout << it->first << "=" << it->second << std::endl;
  }
#endif
}

void DarkMatterParametersFactory::DumpToFile(std::ofstream* outFile)
{
  DarkMatterParametersRegistry* DMreg = DarkMatterParametersRegistry::GetInstance();
  DMreg->DumpToFile(outFile);

#if 0
  std::map<std::string, double>::iterator it;
  for (it = fFactoryParams.begin(); it != fFactoryParams.end(); it++) {
    *outFile << it->first << "=" << it->second << std::endl;
  }
#endif
}


bool DarkMatterParametersFactory::ExistsRegisteredParam(std::string name)
{
  DarkMatterParametersRegistry* DMreg = DarkMatterParametersRegistry::GetInstance();
  return DMreg->ExistsRegisteredParam(name);

#if 0
  std::map<std::string, double>::iterator it = fFactoryParams.find(name);
  return (it != fFactoryParams.end());
#endif
}


DarkMatterParametersFactory* DarkMatterParametersFactory::GetInstance()
{
  if (!fFactoryInstance) fFactoryInstance = new DarkMatterParametersFactory("DMParameters");
  return fFactoryInstance;
}
