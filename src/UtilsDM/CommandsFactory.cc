#include "CommandsFactory.hh"

// STD Library
#include <iostream>
#include <utility>
#include <stdlib.h>
#include <exception>


CommandsFactory* CommandsFactory::fCmdFactory = nullptr;

CommandsFactory::CommandsFactory(std::string Name)
  : fName(Name)
{
  // ...
}

CommandsFactory::~CommandsFactory()
{
  // ...
}

void CommandsFactory::AddNewCmdWithAnInt(std::string CmdName, int value)
{
  std::pair<kType::CmdInt::iterator, bool> insert 
    = fCmdsWithAnInt.insert(std::make_pair(CmdName, value));
  if (!insert.second) {
    std::cerr << "Command  " << CmdName << "  already booked, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void CommandsFactory::AddNewCmdWithADouble(std::string CmdName, double value)
{
  std::pair<kType::CmdDouble::iterator, bool> insert 
    = fCmdsWithADouble.insert(std::make_pair(CmdName, value));
  if (!insert.second) {
    std::cerr << "Command  " << CmdName << "  already booked, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void CommandsFactory::AddNewCmdWithAString(std::string CmdName, std::string value)
{
  std::pair<kType::CmdString::iterator, bool> insert 
    = fCmdsWithAString.insert(std::make_pair(CmdName, value));
  if (!insert.second) {
    std::cerr << "Command  " << CmdName << "  already booked, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void CommandsFactory::SetCmdValueWithAnInt(std::string CmdName, int& value)
{
  kType::CmdInt::iterator it = fCmdsWithAnInt.find(CmdName);
  if (it != fCmdsWithAnInt.end()) {
    it->second = value;
  }
  else {
    std::cerr << "Command " << CmdName << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void CommandsFactory::SetCmdValueWithADouble(std::string CmdName, double& value)
{
  kType::CmdDouble::iterator it = fCmdsWithADouble.find(CmdName);
  if (it != fCmdsWithADouble.end()) {
    it->second = value;
  }
  else {
    std::cerr << "Command " << CmdName << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void CommandsFactory::SetCmdValueWithAString(std::string CmdName, std::string& value)
{
  kType::CmdString::iterator it = fCmdsWithAString.find(CmdName);
  if (it != fCmdsWithAString.end()) {
    it->second = value;
  }
  else {
    std::cerr << "Command " << CmdName << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

int CommandsFactory::GetCmdValueWithAnInt(std::string CmdName)
{
  kType::CmdInt::iterator it = fCmdsWithAnInt.find(CmdName);
  if (it != fCmdsWithAnInt.end()) {
    return it->second;
  }
  else {
    std::cerr << "Command " << CmdName << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

double CommandsFactory::GetCmdValueWithADouble(std::string CmdName)
{
  kType::CmdDouble::iterator it = fCmdsWithADouble.find(CmdName);
  if (it != fCmdsWithADouble.end()) {
    return it->second;
  }
  else {
    std::cerr << "Command " << CmdName << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

std::string CommandsFactory::GetCmdValueWithAString(std::string CmdName)
{
  kType::CmdString::iterator it = fCmdsWithAString.find(CmdName);
  if (it != fCmdsWithAString.end()) {
    return it->second;
  }
  else {
    std::cerr << "Command " << CmdName << " not found, exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
}

CommandsFactory* CommandsFactory::GetInstance()
{
  if (!fCmdFactory) fCmdFactory = new CommandsFactory("Factory");
  return fCmdFactory;
}

void CommandsFactory::SetInputFile(std::string Name)
{
  fInputFile.open(Name, std::ios_base::in);
}

void CommandsFactory::SetCmdDirectory(std::string CmdDirectory)
{
  fCmdDirectory = CmdDirectory;
}

void CommandsFactory::ParseInputFile()
{
  double CmdDoubleValue;
  std::string CmdStringValue;

  // loop over the file
  std::string str;
  while (getline(fInputFile, str)) {
    if (str.find(fCmdDirectory) != std::string::npos) { // contains directory keyword
      CmdDoubleValue = 0.0;
      std::string CmdName = str.substr(0, str.find_first_of(' '));
      CmdStringValue = str.substr(str.find_last_of(' ')+1, str.length());

      // check whether the string is a double or not
      std::string::size_type st;
      try {
        CmdDoubleValue = std::stod(CmdStringValue, &st);
        this->SetCmdValueWithADouble(CmdName, CmdDoubleValue);
      }
      catch(std::exception& err) {
        this->SetCmdValueWithAString(CmdName, CmdStringValue);
      }

    }
  }
}

