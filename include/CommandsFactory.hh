#ifndef CommandsFactory_hh
#define CommandsFactory_hh 1

// STD Library
#include <string>
#include <fstream>
#include <map>

namespace kType{
  typedef std::map<std::string, int> CmdInt;
  typedef std::map<std::string, double> CmdDouble;
  typedef std::map<std::string, std::string> CmdString;
}

class CommandsFactory{

  public:
    CommandsFactory(std::string);
    ~CommandsFactory();

    // add methods
    void AddNewCmdWithAnInt(std::string, int);
    void AddNewCmdWithADouble(std::string, double);
    void AddNewCmdWithAString(std::string, std::string);

    // set methods
    void SetCmdValueWithAnInt(std::string, int&);
    void SetCmdValueWithADouble(std::string, double&);
    void SetCmdValueWithAString(std::string, std::string&);

    // parsing
    void SetInputFile(std::string);
    void SetCmdDirectory(std::string);
    void ParseInputFile();
    std::string GetInputFileName() {return fName;}
    std::ifstream& GetInputFile() {return fInputFile;}

    // get methods
    int         GetCmdValueWithAnInt(std::string);
    double      GetCmdValueWithADouble(std::string);
    std::string GetCmdValueWithAString(std::string);

    static CommandsFactory* GetInstance();

  private:

    std::string fName;

    std::ifstream fInputFile;

    std::string fCmdDirectory;
    kType::CmdInt    fCmdsWithAnInt;
    kType::CmdDouble fCmdsWithADouble;
    kType::CmdString fCmdsWithAString;

    static CommandsFactory* fCmdFactory;

};

#endif // CommandsFactory_hh
