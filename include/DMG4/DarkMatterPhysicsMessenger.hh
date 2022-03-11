#ifndef DarkMatterPhysicsMessenger_hh
#define DarkMatterPhysicsMessenger_hh 1

#include <string>

class CommandsFactory;
class DarkMatterPhysics;

class DarkMatterPhysicsMessenger{
  
  public:
    DarkMatterPhysicsMessenger(std::string);
    ~DarkMatterPhysicsMessenger();

    void SetPhysicsList(DarkMatterPhysics*);
    void ApplyCmd();

    CommandsFactory* GetCmdFactory() {return DarkMatterCmd;}

  private:

    CommandsFactory* DarkMatterCmd;
    DarkMatterPhysics* fPhysList;

    std::string fName;

};

#endif // DarkMatterPhysicsMessenger_hh
