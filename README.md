# persisters

Supplementary code for the preprint "Persistence as an optimal hedging strategy" [available on bioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.19.883645v3)

## Code

### Requirements
* MATLAB

### Files

* The file `Persisters.m` is a script run to produce results in the main document. To reproduce results in the main document, only this file needs to be modified.
  * To switch environment, modify line 22. For example, `Env = Envs{3}` corresponds to environment 3, the Poisson environment.
* The file `Parameters.m` stores non-environment specific parameters.
* The file `StateEquations.m` stores non-environment specific state equations (namely, the state equations for n and theta).
* The `Environments` directory contains environment specific code for each of the following environments:
  1. Constant (`Env_1_Constant.m`)
  2. Monod (`Env_2_Monod.m`)
  3. Poisson (`Env_3_Poisson.m`)
  4. Ornstein-Uhlenbeck (`Env_4_OrnsteinUhlenbeck.m`) (results in supporting material only)
  5. Duffing (`Env_5_Duffing.m`) (results in supporting material only)
* The `Environment` directory also contains the file `EnvironmentSeed.m` which is called to return the random number generator (RNG) seed used to reproduce results in the paper.
* The `HJB` directory contains code used to solve the HJB PDE for the persister problem. 
  * The file `HJB_Persisters.m` is called to solve the PDE
  * The file `HJB_Forward_Persisters.m` is called to solve the SDE forward, coupling the control to the solution of the PDE
  * The file `HJB_CreateGrid.m` is called to create the spatial mesh
  * The file `BoundaryExtension.m` is called to approximate derivatives on the boundary through linear interpolation.
  * Details on the solution technique, including the boundary approximation method, are provided in the supporting material document.
