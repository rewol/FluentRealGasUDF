# FluentRealGasUDF
A User Defined Function for Real-Gas Thermodynamics for FLUENT CFD simulations for 1-step CH4-O2 global mechanism.

AUTHOR          : REFİK ALPER TUNCER
DATE / REVISION : FEBRUARY 2018 / R0


DESCRIPTION     : This UDF file includes PENG-ROBINSON real-gas thermodynamic models and transport models (viscosity, thermal                               conductivity and mass diffusivity) with pressure-correction. It is purpose is to be utilized in high-pressure                             combustion applications in FLUENT simulations.

CONTENTS        : Peng-Robinson CEOS is utilized to compute mixture density. Thermodynamic properties are calculated with departure                         functions and ideal properties for the thermodynamic properties are calculated from NASA coefficients. Also, for                           CRYOGENIC TEMPERATURES, a VOLUME-CORRECTION method, proposed  by ABUDOUR (Fluid Phase Equilib. 349 37–55) is                               employed. As for viscosity and thermal conductivity, CHUNG' s method is employed, since it is very common and                             accurate in high-pressure conditions. For mass diffusivity, TAKAHASHI' s method is utilized.

REVISION NOTES  :


