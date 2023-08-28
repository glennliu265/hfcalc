hfcalc

Scripts to calculate heat flux feedback (HFF).

# Introduction and Notes



# Organization
See legend below

## Main
### Preprocessing
    - [preproc_ncep] 				: (G) Prep NCEP Reanalysis for HFF calculations
    - [preproc_ERA20C]                          : (G) Prep ERA20C Reanalysis for HFF calculations
    - [preproc_CESMPIC]                         : (C) Re-process CESM1 output to match output of [preproc_ncep]
    - [preproc_CESM1_LENS]   			: (C) Prep CESM1 TS + Heat Flux data for calculations (landice mask + QNET + Ensemble Average + regridding)
    - [preproc_damping_LENS] 			: (G) Prep CMIP5 LENS for calculations (regrid + landice mask + ensemble average + QNET)
### ENSO and HFF Calculation
    - [calc_enso_general]    			: (G) For all LENS datasets: (1) preprocess data (2) Compute ENSO Index (3) Regress and Remove ENSO from variables (4) Compute Heat Flux Feedback
    - [calc_ENSO_PIC]                   	: (C) Calculate ENSO index in CESM1 PiControl
    - [calc_HF_func] 				: (C) Old script that first turned HFF calculations into functions. Was used in stochastic model project.
### Postprocessing and HFF Combination
    - [prep_HF]      				: (C) Old script that preps output from [calc_HF_func] for stochastic model input.
    - [combine_damping_lens] 			: (G) Combine separate HFF components (general CMIP5 LENS version of [combine_damping_CESM1LE]

## Analysis
    - [viz_hfdamping_CESM1PIC]     		: (C) Visualize HFF for CESM1 Preindustrial Control
    - [viz_hf_comparison]          		: (G) Compare HFF in different datasets computed with [calc_enso_general]
    - [compare_damping_components] 		: (C) Compare components of HFF across CESM1 scenarios (historical, rcp85)
    - [compare_damping_CESM1LE]         	: (C) Compare HFF between Historical and RCP85 of CESM1
    - [compare_damping_lens]            	: (G) Compare HFF across CMIP5 LENs, adapted from [compare_damping_CESM1LE]
    - [calc_autocorrelation_SLAB_PIC]    	: (C) Compute autocorrelation of regionally-averaged SSTs with ENSO removed. Stochastic model project?
    - [viz_hfdamping_CESM1LE] 			: (C) Visualize HFF in CESM1, script version of [viz_hfdamping_CESM1LE.ipynb] 
    - [crop_lens_atm]              		: (V) Crop atmospheric variables for vertical gradient analysis
    - [viz_atmvar_verticaldiff]         	: (V) Visualize output of [crop_lens_atm]
    - [check_UBOT]        			: (V) Compare lower level of U with other variables

## MATLAB Scripts
Old versions of the sript written for MATLAB

## Scrap
Includes older versions of scripts, etc
    - [remove_ENSO]  			: (C) Remove ENSO component of a variable using regression. Works with CESM1 PiControl (SLAB)
    - [compare_hff_calc_method] 	: (D) Compare output of [calc_enso_general] with original MATLAB scripts
    - [hfdamping_viz]                   : (C) Old script for visualizing HFF output with significance
    - [test_lag_monwin]                 : (D) Debugging script for month windows in lag calculations
    - [viz_ENSO_PIC]                    : (C) Visualize ENSO calculated with [calc_ENSO_PIC]
    - [viz_hfdamping_CESM1LE.ipynb]     : (C) Visualize HFF in CESM1 (see non python notebook version in Analysis)

## Legend
(C) --> Scripts originally written for analysis in CESM1 Large Ensemble, and for comparing historical + RCP85 periods
(G) --> Scripts adapted for HFF calculations over general CMIP5 large ensembles
(V) --> Additional analysis on vertical gradients to understand heat flux estimates
(D) --> Debugging

