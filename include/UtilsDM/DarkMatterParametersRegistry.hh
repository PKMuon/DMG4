#pragma once

// STD Library
#include <string>
#include <map>
#include <fstream>

class DarkMatterParametersRegistry
{
  public:
    DarkMatterParametersRegistry(std::string);
    ~DarkMatterParametersRegistry();

    static DarkMatterParametersRegistry* GetInstance();

    void RegisterNewParam(std::string, double);
    void UpdateRegisteredParam(std::string, double&);
    double GetRegisteredParam(std::string);
    double GetRegisteredParam(std::string, double);
    void PrintListOfRegisteredParams();
    void DumpToFile(std::ofstream*);

    bool ExistsRegisteredParam(std::string);

    void SetRegistryName(std::string registryName) {fRegistryName = registryName;}

  private:
    static DarkMatterParametersRegistry* fRegistryInstance;

    std::string fRegistryName;
    std::map<std::string, double> fRegistryParams;
};
