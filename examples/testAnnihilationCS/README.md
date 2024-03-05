# Resonant Annihilation Cross-Section Test

## Overview

This repository contains an example script to test the implementation of the resonant annihilation
production cross-section in the DMG4 package.

## Usage

Compile the example:

```bash
../clean
./runconfigure_roothisto
make
```

Run the example script:

```bash
./mkgeant
```

The output is a ROOT file containing plots for the *total resonant annihilation* cross-section as 
a function of the energy as well as for the *angular sampling* for the maximum cross-section at
$`E_{e^{+}} = m_{Z'}^2/(2\times m_{e^{-}})`$.
