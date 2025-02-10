# Test ALP CS

This is example allows to display the differential cross-section as a function of the angle in the case of ALP production. 

Usage:
```bash
  ../clean        # Remove older CMake files
  runconfigure
  make
  source config.sh
  ./mkgeant4
```

## OUTPUT

The application prints out the results from the sampling and writes to `output.txt`. The entries of this file correspond to the obtained fractional energy, the polar angle theta and the azimuthal angle phi in the center-of-mass frame:

```
<Frac. Energy> <Theta> <Phi>
```

