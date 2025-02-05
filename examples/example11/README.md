# EXAMPLE 11 (Geant4)

This example features the production of DM in a simple MC simulation using Geant4. The detector
consists of a Pb block, towards which 100 GeV electrons are shot. The beam definition is given
through the cardg.g4mac parameters. 

Usage: type the following commands:
```cpp
  runconfigure
  make
  source config.csh  #  (or "source config.sh" for sh, bash shells, this is to be done only once)
  mkgeant4 -w cardg.g4mac -p<Param>=<Value>
```

See the definition of all possible parameters in DarkPhysicsConfigure.cc.
Adjust mass, coupling and bias factors as needed.

## OUTPUT

The output consists of a text file named `Report.txt`. It contains all the information saved in case
of DM production. The structure of the output is defined in EventInfo and is as follows:

```
...

EVENT
<Event ID> <DM Emission ID>
DMINFO <Nr. of DM Processes>
<Prod. Vertex X> <Prod. Vertex Y> <Prod. Vertex Z>
<Decay Vertex X> <Decay Vertex Y> <Decay Vertex Z>
<Parent    PID> <Parent    E> <Parent    Mom. X> <Parent    Mom. Y> <Parent    Mom. Z>
<Daugther1 PID> <Daugther1 E> <Daugther1 Mom. X> <Daugther1 Mom. Y> <Daugther1 Mom. Z>
<Daugther2 PID> <Daugther2 E> <Daugther2 Mom. X> <Daugther2 Mom. Y> <Daugther2 Mom. Z>
...
ENDEVENT

...
```

## USAGE WITH DarkMatterPhysicsConfigure.cc

Here are a few examples for DM definitions:

- $`A^{\prime} \rightarrow`$ invis
  `./mkgeant4 -w cardg_test.g4mac -pDMMass=0.0167 -pEpsilon=0.001 -pBiasSigmaFactor0=1.e9 -pDMProcess=1 -pDecayType=0 -pBranchingType=0`

- $`B-L`$ $`Z^{\prime}`$ (production via Dark Bremsstrahlung):
  `./mkgeant4 -w cardg_test.g4mac -pDMMass=0.0167 -pEpsilon=0.001 -pBiasSigmaFactor0=1.e9 -pDMProcess=1 -pDecayType=1 -pBranchingType=1`

- $`B-L`$ $`Z^{\prime}`$ (production via resonant $`e^{-}e^{+}`$ annihilation):
  `./mkgeant4 -w cardg_test.g4mac -pDMMass=0.0167 -pEpsilon=0.001 -pBiasSigmaFactor0=1.e9 -pDMProcess=16 -pDecayType=1 -pBranchingType=10`

