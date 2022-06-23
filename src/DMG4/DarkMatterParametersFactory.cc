#include "DarkMatterParametersFactory.hh"

// STD Library
#include <iostream>
#include <utility>
#include <stdlib.h>

DarkMatterParametersFactory* DarkMatterParametersFactory::fFactoryInstance = nullptr;

DarkMatterParametersFactory::DarkMatterParametersFactory(std::string factoryName)
{
  fFactoryName = factoryName;
}

DarkMatterParametersFactory::~DarkMatterParametersFactory()
{
  // ...
}

void DarkMatterParametersFactory::RegisterNewParam(std::string name, double value)
{
  std::pair<std::map<std::string, double>::iterator, bool> insert 
    = fFactoryParams.insert(std::make_pair(name, value));
  if (!insert.second) {
    std::cerr << "Name  " << name << "  already booked, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void DarkMatterParametersFactory::UpdateRegisteredParam(std::string name, double& value)
{
  std::map<std::string, double>::iterator it = fFactoryParams.find(name);
  if (it != fFactoryParams.end()) {
    it->second = value;
  }
  else {
    std::cerr << "Name " << name << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

double DarkMatterParametersFactory::GetRegisteredParam(std::string name)
{
  std::map<std::string, double>::iterator it = fFactoryParams.find(name);
  if (it != fFactoryParams.end()) {
    return it->second;
  }
  else {
    std::cerr << "Name " << name << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

double DarkMatterParametersFactory::GetRegisteredParam(std::string name, double def)
{
  std::map<std::string, double>::iterator it = fFactoryParams.find(name);
  if (it != fFactoryParams.end()) {
    return it->second;
  }
  else {
    std::cerr << "Name " << name << " not found, using default parameter:" << def << "\n";
    RegisterNewParam(name, def);
    return def;
  }
}

void DarkMatterParametersFactory::PrintListOfRegisteredParams()
{
  std::map<std::string, double>::iterator it;
  for (it = fFactoryParams.begin(); it != fFactoryParams.end(); it++) {
    std::cout << it->first << "=" << it->second << std::endl;
  }
}

void DarkMatterParametersFactory::DumpToFile(std::ofstream* outFile)
{
  std::map<std::string, double>::iterator it;
  for (it = fFactoryParams.begin(); it != fFactoryParams.end(); it++) {
    *outFile << it->first << "=" << it->second << std::endl;
  }
}


bool DarkMatterParametersFactory::ExistsRegisteredParam(std::string name){

    std::map<std::string, double>::iterator it = fFactoryParams.find(name);
    return (it != fFactoryParams.end());
}


DarkMatterParametersFactory* DarkMatterParametersFactory::GetInstance()
{
  if (!fFactoryInstance) fFactoryInstance = new DarkMatterParametersFactory("DMParameters");
  return fFactoryInstance;
}
