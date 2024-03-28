# stress-strain-analysis
Code for data processing and stress-strain analysis for measured strain signals of two microbubble populations. 

**SS_Processing_cropping_strain_CN**

Determines the cross-correlation of all selected strain signals with a simulated strain. The mean delay determines the shift of all signals for further processing.


**SS_Processing_errorminimization_lsqnonlin_CN**

Data after bandpass filtering (not included in this repository) are loaded and processed using script. Strain signals recorded of single microbubbles are selected based on their signal-to-noise ratio and average strain value. Each included strain signal is filtered and its time derivatives are determined in the frequency domain. This script determines the nondimensional viscoelastic pressure contribution as a function of strain and strain rate. The elastic contribution is extracted and the surface tension as a function of strain is determined.

A fitting procedure is used to determine the initial state: initial radius, initial surface tension and relative delay of the microbubble strain w.r.t. the driving ultrasound pulse. To this end, the native function 'lsqnonlin' is used, which uses 100 random starting points within the relevant parameter ranges. The median values corresponding to the lowest 10% residual norms are used as initial state values.

Settings for the analysis are stored in structure "S".


**Functions**

Functions used in the processing are named "SS_func_*.m'
  SS_func_strain_filtering.m          applies a double-band-pass filter around the fundamental frequency and the second harmonic
  SS_func_strain_derivatives.m        determines the time derivatives of the strain; first and second order, in the frequency domain
  SS_func_surface_tension.m           determines the surface tension from the strain with a given initial state
  SS_func_surface_tension_Segers.m    determines the surface tension according to a parametric curve of DPPC:DPPE-PEG5k, 9:1 mol% bubbles, from Segers et al Soft Matter 2018 "High-precision                                       acoustic measurements of the nonlinear dilatational elasticity of phospholipid coated monodisperse microbubbles"
  SS_func_error_lsqnonlin.m           determines the surface tension by first determining the viscoelastic contribution (elastic+viscous contribution of the bubble shell) following the                                            stress-strain analysis, then uses an initial radius, initial surface tension, and delay to determine the surface tension, and the error compared to the                                       parametric curve is calculated.

  lininterp1.m                        faster linear interpolation, see reference below.
  Jeffrey Wu (2024). Faster linear Interpolation (https://www.mathworks.com/matlabcentral/fileexchange/28376-faster-linear-interpolation), MATLAB Central File Exchange.


**SS_analysis_error_minimization_result_combined_CN**

Combines the analyzed data per microbubble to determine the median surface tension and elasticity over time, and determine the distributions of the initial surface tension and initial elasticity.


**SS_Processing_viscous_PCA_CN**

Combines the analyzed data per microbubble to determine the viscous contribution per microbubble as a function of strain and strain rate. Also performs a principal component analysis (PCA) to determine the most dominant feature per point in the strain and strain rate space.
