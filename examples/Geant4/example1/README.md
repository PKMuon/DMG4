# EXAMPLE 1 (Basic Geant4 example)

This G4 application uses default configuration in DMG4/DarkMatterPhysicsConfigure.cc. If you want to change it, copy the file here or try examples/Geant4/example2

Usage: type the following commands:
```bash
  runconfigure
  make
  source config.csh     #  (or "source config.sh" for sh, bash shells, this is to be done only once)
  ./mkgeant4 cardg.g4mac
```

## OUTPUT

This application outputs a file `Report.txt` that saves only the number of produced DM particles. 
