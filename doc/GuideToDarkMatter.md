# Simulating Dark Matter with DMG4

## DarkMatterPhysicsConfigure

The following tables provide an overview of the different Dark Matter (DM) particles and 
processes implemented in DMG4, according to the production channel.

For given values of `DMProcessType`, `DecayType` and `BranchingType` in the 
`DarkMatterPhysicsConfigure()` the following particles and processes are
initialized:

### Bremstrahlung production

| Description |`DMParticle` | `DMProcess` | `DMProcessType` | `DecayType` | `BranchingType`  | 
| ------ | ------ | ------ | ------ | ------ | ------ |
| $`A^{\prime}\rightarrow \mathrm{invis}`$ (no decay simulated) | `DMParticleAPrime` | `DMProcessDMBrem` | 1 | 0 | any |
| $`A^{\prime}\rightarrow e^{-}e^{+}`$ | `DMParticleAPrime` | `DMProcessDMBrem` | 1 | 1 | 0 |
| Inelastic DM: <br>$`A^{\prime}\rightarrow \chi_{1}\chi_{2}`$ | `DMParticleAPrime` | `DMProcessDMBrem` | 1 | 2 | 2 |
| Dirac inelastic DM: <br>$`A^{\prime}\rightarrow \chi_{1}\chi_{1}`$<br>$`A^{\prime}\rightarrow \chi_{2}\chi_{2}`$<br>$`A^{\prime}\rightarrow \chi_{1}\chi_{2}`$ | `DMParticleAPrime` | `DMProcessDMBrem` | 1 | 2 | 3 |
| $`S\rightarrow \mathrm{invis}`$ | `DMParticleScalar` | `DMProcessDMBrem` | 2 | 0 | any |
| $`S\rightarrow e^{-}e^{+}`$ | `DMParticleScalar` | `DMProcessDMBrem` | 2 | 1 | any |
| $`A\rightarrow \mathrm{invis}`$ | `DMParticleAxial` | `DMProcessDMBrem` | 3 | 0 | any |
| $`A\rightarrow e^{-}e^{+}`$ | `DMParticleAxial` | `DMProcessDMBrem` | 3 | 1 | any |
| $`P\rightarrow \mathrm{invis}`$ | `DMParticlePseudoScalar` | `DMProcessDMBrem` | 4 | 0 | any |
| $`P\rightarrow e^{-}e^{+}`$ | `DMParticlePseudoScalar` | `DMProcessDMBrem` | 4 | 1 | any |
| Spin 2 particle:<br>$`X\rightarrow \mathrm{invis}`$ | `DMParticleAPrime` | `DMProcessDMBrem` | 5 | 0 | any |

### Annihilation production

| Description |`DMParticle` | `DMProcess` | `DMProcessType` | `DecayType` | `BranchingType`  | 
| ------ | ------ | ------ | ------ | ------ | ------ |
| $`S\rightarrow \mathrm{invis}`$ | `DMParticleScalar` | `DMProcessAnnihilation` | 12 | 0 | any |
| $`S\rightarrow e^{-}e^{+}`$ | `DMParticleScalar` | `DMProcessAnnihilation` | 12 | 1 | any |
| $`A\rightarrow \mathrm{invis}`$ | `DMParticleAxial` | `DMProcessAnnihilation` | 13 | 0 | any |
| $`A\rightarrow e^{-}e^{+}`$ | `DMParticleAxial` | `DMProcessAnnihilation` | 13 | 1 | any |
| $`P\rightarrow \mathrm{invis}`$ | `DMParticlePseudoScalar` | `DMProcessAnnihilation` | 14 | 0 | any |
| $`P\rightarrow e^{-}e^{+}`$ | `DMParticlePseudoScalar` | `DMProcessAnnihilation` | 14 | 1 | any |
| Spin 2 particle:<br>$`X\rightarrow \mathrm{invis}`$ | `DMParticleAPrime` | `DMProcessAnnihilation` | 15 | 0 | any |

### Primakoff production

| Description |`DMParticle` | `DMProcess` | `DMProcessType` | `DecayType` | `BranchingType`  | 
| ------ | ------ | ------ | ------ | ------ | ------ |
| $`a\rightarrow \mathrm{invis}`$ | `DMParticleALP` | `DMProcessPrimakoffALP` | 21 | 0 | any |
| $`a\rightarrow \gamma \gamma`$ | `DMParticleALP` | `DMProcessPrimakoffALP` | 21 | 1 | any |
