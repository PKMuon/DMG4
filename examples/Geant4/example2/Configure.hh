#pragma once

#include <map>
#include <string>
#include <cmath>

#include <G4Types.hh>

void print_usage( std::ostream & os, const char * appName );
int parseArguments( int argc , char ** argv, char *& wCardFile );

// This class is used to propagate input parameters such as epsilon 
// and mass provided in batch mode 
// in the subsequent Dark Matter and DMG4 classes

class PhysicsParametersSet : public std::map<std::string, double> {
  private:
    static PhysicsParametersSet * _self;
  public:
    static PhysicsParametersSet & self();
    G4double get(const std::string & name, double default_=NAN) const;
};  // class PhysicsParametersSet
