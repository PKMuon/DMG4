#include "Configure.hh"
#include <iostream>
#include <map>
#include <string>
#include <unistd.h>
#include <cmath>

#include <G4ios.hh>


PhysicsParametersSet * PhysicsParametersSet::_self = nullptr;

PhysicsParametersSet & PhysicsParametersSet::self() {
  if( ! _self ) {
    _self = new PhysicsParametersSet();
  }
  return *_self;
}


G4double PhysicsParametersSet::get(const std::string & name, double default_) const {
  std::map<std::string, double>::const_iterator it = find(name);
  if( end() != it ) {
    return it->second;
  }
  if( std::isnan(default_) ) {
    char ebf[128];
    snprintf( ebf, sizeof(ebf), "Value \"%s\" is not set."
        , name.c_str() );
    throw std::runtime_error(ebf);
  }
  return default_;
}


void print_usage( std::ostream & os, const char * appName ) {
  os << "Usage:" << std::endl
    << "  $ " << appName
    << "-p <param name>=<param value>" << std::endl
    ;
}


int parseArguments( int argc , char ** argv, char *& wCardFile ) {
  int opt;
  while(-1 != (opt = getopt(argc, argv, "hp:w:"))) {
    switch(opt) {
      case 'p' : {
                   const std::string expr( optarg );
                   size_t pDelim = expr.find( '=' );
                   if( std::string::npos == pDelim || 3 >= expr.size() ) {
                     G4cerr << "Unable to interpret option \"" << optarg
                       << "\"." << std::endl;
                     return -1;
                   }
                   std::string key = expr.substr( 0, pDelim )
                     , value = expr.substr( pDelim + 1, std::string::npos );
                   PhysicsParametersSet::self().emplace( key, atof(value.c_str()) );
                 } break;
      case 'w' : {
                   wCardFile = optarg;
                 } break;
      case 'h' : {
                   print_usage( std::cout, argv[0] );
                   return 1;
                 } break;
      default:
                 // Do nothing and rely on default values
                 std::cout << "INFO: No input parameters detected!" << std::endl;
                 return 1;
    }
  }
  if( !wCardFile ) {
    std::cerr << ".mac file is not set. Use -w optarg." << std::endl;
    return -1;
  }
  return 1;
}
