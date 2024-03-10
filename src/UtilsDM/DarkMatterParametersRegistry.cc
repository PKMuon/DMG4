#include "DarkMatterParametersRegistry.hh"

// STD Library
#include <iostream>
#include <utility>
#include <stdlib.h>

DarkMatterParametersRegistry* DarkMatterParametersRegistry::fRegistryInstance = nullptr;

DarkMatterParametersRegistry::DarkMatterParametersRegistry(std::string registryName)
{
  fRegistryName = registryName;
}

DarkMatterParametersRegistry::~DarkMatterParametersRegistry()
{
  // ...
}

void DarkMatterParametersRegistry::RegisterNewParam(std::string name, double value)
{
  std::pair<std::map<std::string, double>::iterator, bool> insert 
    = fRegistryParams.insert(std::make_pair(name, value));
  if (!insert.second) {
    std::cerr << "Name  " << name << "  already booked, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void DarkMatterParametersRegistry::UpdateRegisteredParam(std::string name, double& value)
{
  std::map<std::string, double>::iterator it = fRegistryParams.find(name);
  if (it != fRegistryParams.end()) {
    it->second = value;
  }
  else {
    std::cerr << "Name " << name << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

double DarkMatterParametersRegistry::GetRegisteredParam(std::string name)
{
  std::map<std::string, double>::iterator it = fRegistryParams.find(name);
  if (it != fRegistryParams.end()) {
    return it->second;
  }
  else {
    std::cerr << "Name " << name << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

double DarkMatterParametersRegistry::GetRegisteredParam(std::string name, double def)
{
  std::map<std::string, double>::iterator it = fRegistryParams.find(name);
  if (it != fRegistryParams.end()) {
    return it->second;
  }
  else {
    std::cerr << "Name " << name << " not found, using default parameter:" << def << "\n";
    RegisterNewParam(name, def);
    return def;
  }
}

void DarkMatterParametersRegistry::PrintListOfRegisteredParams()
{
  std::map<std::string, double>::iterator it;
  for (it = fRegistryParams.begin(); it != fRegistryParams.end(); it++) {
    std::cout << it->first << "=" << it->second << std::endl;
  }
}

void DarkMatterParametersRegistry::DumpToFile(std::ofstream* outFile)
{
  std::map<std::string, double>::iterator it;
  for (it = fRegistryParams.begin(); it != fRegistryParams.end(); it++) {
    *outFile << it->first << "=" << it->second << std::endl;
  }
}


bool DarkMatterParametersRegistry::ExistsRegisteredParam(std::string name){

    std::map<std::string, double>::iterator it = fRegistryParams.find(name);
    return (it != fRegistryParams.end());
}


DarkMatterParametersRegistry* DarkMatterParametersRegistry::GetInstance()
{
  if (!fRegistryInstance) fRegistryInstance = new DarkMatterParametersRegistry("DMParameters");
  return fRegistryInstance;
}
