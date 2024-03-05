#pragma once

// STD Library
#include <string>
#include <map>
#include <fstream>

class DarkMatterParametersFactory
{
  public:
    DarkMatterParametersFactory(std::string);
    ~DarkMatterParametersFactory();

    static DarkMatterParametersFactory* GetInstance();

    void RegisterNewParam(std::string, double);
    void UpdateRegisteredParam(std::string, double&);
    double GetRegisteredParam(std::string);
    double GetRegisteredParam(std::string, double);
    void PrintListOfRegisteredParams();
    void DumpToFile(std::ofstream*);

    bool ExistsRegisteredParam(std::string);

    void SetFactoryName(std::string factoryName) {fFactoryName = factoryName;}

  private:
    static DarkMatterParametersFactory* fFactoryInstance;

    std::string fFactoryName;
    std::map<std::string, double> fFactoryParams;
};

