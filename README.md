# ADRESSTools: A suite of data analysis tools for photoelectron spectroscopy at ADRESS

A compilation of MATLAB scripts used for the analysis of soft X-ray angle-resolved photoemission spectroscopy (SX-ARPES) experiments that give direct access to the electronic band-structure of a material. Designed to be directly compatible with the data format of SX-ARPES experiments at the ADRESS beamline, at the Swiss Light Source (SLS) in the Paul Scherrer Institute (PSI), but can be generalised to other data formats if required. The scripts can also be used for perform x-ray photoelectron spectroscopy (PES) curve fitting, using a variety of background shapes.

## Installation  
1. Download the *ADRESSTools* repository.
2. Open MATLAB and use *Set Path* in the *Home* tab to add the *ADRESSTools* repository and all its sub-folders into its saved search paths.
3. Make sure you also use *Set Path* to add the repository / folder that contains all of your data to be loaded in.
4. In the *PESTools* folder, there is an *Examples* folder, which contains many ARPES / XPS data processing and curve fitting templates / examples that can be used. You should use this as a starting point.

## Materials Database Tools
**(1) Material Properties Database (MPD)**: This is a local MATLAB database that compiles the most useful physical, electronic, optical and magnetic material properties of elements / compounds. The data is taken from a range of sources, where the 'average' values are used for parameters that had more than 1 unique value. Accessible in MATLAB via 'get_mpd_props()'.
![000_MPD](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/000_MPD.png)

**(2) Electron Inelastic Mean Free Path Database (eIMFPD)**: This is a local MATLAB database that compiles the optical data from the NIST Electron Inelastic-Mean-Free-Path Database (http://dx.doi.org/10.18434/T48C78) so that the results can be easily called and accessed within MATLAB. This is accessible in MATLAB via 'get_eimfpd_props()'. Furthermore, predictive eIMFP formulas's are also available using the (1) Universal, (2) TPP-2M and (3) S1 & S2 formalisms.
![000_eIMFPD](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/000_eIMFPD.png)

**(3) Photoionisation Cross-Section and Asymmetry Database (PIXSAD)**:
This is a local MATLAB database that compiles the photoionisation cross-section and asymmetry parameter data of the elements. This is useful when modelling the total photoelectron intensity that originates from a given layer of a sample when performing ARPES / XPS experiments. Accessible in MATLAB via 'get_pixsad_props()'.
![000_PIXSAD](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/000_PIXSAD.png)

**(4) Photoionisation Energy and Fluorescence Database (PIEFD)**:
This is a local MATLAB database that compiles the photoionisation energy and fluorescence yield data of the elements. Accessible in MATLAB via 'get_piefd_props()'.
![000_PIEFD](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/000_PIEFD.png)

**(5) Physics Constants Database**:
This is a MATLAB function that loads in many physics constants that can be used for data processing, or modelling of physical systems. Accessible in MATLAB via 'physics_constants()'.


## Photoelectron Spectroscopy (PES) Analysis Tools
**PES Multilayer Modelling**:
We define function that allows the user simulate the photoelectron intensity contributions of an N-layered sample as a function of photon energy. This can be used to model PES data acquired versus photon energy. The function allows the user to define the material type of each layer, the thickness, the core-level being probed and then determines the emitted photoelectron intensity of each layer using the Beer-Lambert law. All the material, photoionisation and core-level parameters are extracted from the Materials Database Tools.
![020_PES_Intensity_Modelling](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/020_PES_Intensity_Modelling.png)

**PES Model Curves**:
A whole suite of PES model curves are defined and available for the user to use within MATLAB. From Gaussian, Lorentzian, Voigt and Doniach-Sunjic curve shapes. Furthermore, the Fermi-Dirac Distribution (FDD) is also defined, as well as its Gaussian broadened and linear background form. Combining all of these together, the fundamental curve-shapes for PES curve fitting to XPS / ARPES data is possible. Functions for background subtracting PES data according to a polynomial or Shirley background is also available.
![003_Curve_Shapes](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/003_Curve_Shapes.png)
![003_PES_Model_Curves](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/003_PES_Model_Curves.png)
![003_PES_Backgrounds](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/003_PES_Backgrounds.png)

**PES Curve Fitting**:
The PESTools functions offer a host of different line-shapes and background that you can use for fitting PES spectra. The available functions are; (i) For PES curves: Gaussian, Lorentzian, Voigt, asymmetrical Voigt and Doniach-Sunjic functions are available; (ii) For Backgrounds: Linear, Shirley, Offset Shirley and linear step functions are available. A host of miscellaneous tools are also available; (i) peak location finder, (ii) full-width at half-maximum (FWHM) finder, (iii) Fermi-Dirac Distribution (FDD) fitter for binding energy referencing XPS data. The general objective of fitting PES data is to define an initial photoelectron spectrum that provides a reasonable first guess and then run an optimisation procedure to hone in to the best fit solution. Here, the general procedures are listed below;  
(1) Define the shape and properties of the N primary (+ spin-orbit if required) components to be fitted  
(2) Define the type of background to be used  
(3) Preview the initial conditions of the XPS curve fit model vs data (with 'pes2ncurve_view_init()')  
(4) If you are happy with the initial guess, run the optimisation algorithm (with 'pes2ncurve_solver()')  
(5) Extract the best fit variables of the converged fit and plot the final fit with residuals (with 'pes2ncurve_view_fit()')
![022_PES_Curve_Fitting](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/022_PES_Curve_Fitting.png)


## Angle-Resolved Photoelectron Spectroscopy (ARPES) Analysis Tools
**ARPES Data Processing**:
A suite of tools are available to process and analyse ARPES data, along with several miscellaneous tools that can crop, filter or correct the variables and data associated with the ARPES spectra. The 1st stage is to align the binding energies to the valence band maximum (or Fermi-edge), the 2nd is to normalise the intensity of the ARPES spectra over the scan parameter and the 3rd and final step is to convert the angles into wave-vectors.  
![013_ARPES_Data_Processing](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/013_ARPES_Data_Processing.png)

**ARPES Model Curves**:
A whole suite of ARPES model curves are defined and available for the user to use within MATLAB. These are 2D parabolic dispersions that have Gaussian spot profiles, whose FWHM along the x- and y-axes can be defined. Furthermore, the Fermi-Dirac Distribution (FDD) is also defined, as well as its Gaussian broadened and linear background form. Combining these, parabolic ARPES curves can be simulated and used for: (i) accurately fitting N quantum well states near the Fermi-edge or (ii) estimating the VBM position by fitting etc.
![003_ARPES_Model_Curves](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/003_ARPES_Model_Curves.png)

**ARPES State Fitting**:
Allows the user to walk through all of the fitting steps to parabolic features within ARPES spectra. The 1st stage is to load in the ARPES data, and perform any pre-fitting corrections (energy re-alignment, resampling, filtering and cropping). The 2nd stage is to crop around the region of interest and also some reference background with the same density of points. The option for polynomial background subtraction is also available using the reference background. The 3rd stage is to the perform the fitting to the parabolic feature using several different approaches and routines. The Fermi-Dirac distribution can be defined and used within the fitting constraints. A parabolic approximation can be set, whose parabolic trail can be used as a constraint to the fitting parameters to allow the solution to converge easier. Various constraints can also be placed on EDC and MDC cuts through the region of interest, but through trial-and-error, the best solution can be found.  
![013_ARPES_Curve_Fitting](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/013_ARPES_Curve_Fitting.png)

**kF Analysis**:
Allows the user to determine how Kf varies as a function of either the photon energy or tilt angle. This is important when it comes to finding out whether a state is behaving with 2D or 3D characteristics. Can also be used to track how the width of a feature changes as a function of a scan parameter.  
![005_ARPESkf](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/005_ARPESkf.png)

**IsoE Analysis**:
Allows the user to load in single or multiple processed ARPES data-files so that a consistent iso-energetic map can be created even across multiple data-files (if scans were segmented, rather than a single whole scan). The 2D Brillouin-Zone slice can be defined and overlaid on top of the iso-energetic slice. Multiple corrections can also be performed; filtering, background subtraction, origin, shear and scale corrections. There is also the option to determine the area over some threshold, so that the Luttinger Area can be determined to estimate the number density of electrons that are occupying the states.  
![005_ARPESisoe](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/005_ARPESisoe.png)

**BZ Navigation**:
Allows the user to define any crystal geometry and to determine the Wigner-Seitz cell in real-space and the First Brilluoin-Zone in reciprocal space. Additionally, by defining the crystal planes (100), (010) or (001), planar cuts through the stacked Brillouin-Zone can be determined, which is important in the interpretation of ARPES spectra. By imposing the ARPES geometry and experimental variables, the user can define the photon energy being used and extract the ARPES cut/slice being probed through the 2D Brillouin-Zone. This is very useful when it comes to finding what photon energies are needed to probe certain high-symmetry lines or valleys.
![005_BZNavi](ADRESSTools_Dionysus/PESTools_PCC/0_ReadMeImages/005_BZNavi.png)


## Miscellaneous Tools 
**ADRESS Job File Generators**:
Allows the user to generate custom job-files in the form of a text-file that is compatible to use in the 'Restore Tab' button on 'SmartGUI' at the ADRESS beamline. This allows the user to very quickly create a list of ARPES / PES scans to be taken as a function of photon energy, time or position, where each scan can be referenced to a known Fermi-edge, if desired. List of available job file generators:  
- 'ADRESS_JobFileGenerator_CLvsTIME()'  
- 'ADRESS_JobFileGenerator_CLvsHV()'  
- 'ADRESS_JobFileGenerator_CLvsPOS()'  


**Data Loaders**:
A MATLAB data loader for the ARPES / PES data acquired at the ADRESS, SIS and PEARL beamlines in the Swiss Light Source (SLS) are all available. 


## MATLAB Version control  
MATLAB version:   2022a  
MATLAB add-ons (recommended): Database Toolbox, Image Processing Toolbox, Global Optimization Toolbox, Optimization Toolbox, Curve Fitting Toolbox, Parallel Processing Toolbox, Statistics and Machine Learning Toolbox, Signal Processing Toolbox.


## Authors
MaterialsDatabase & PESTools Scripts:  
**Dr. Procopios Constantinou**,  
Swiss Light Source (SLS),  
Paul Scherrer Institute (PSI),  
email: procopios.constantinou@psi.ch

MATools Scripts:  
**Dr. Vladimir Strocov**,  
Swiss Light Source (SLS),  
Paul Scherrer Institute (PSI),  
email: vladimir.strocov@psi.ch

## License  
This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.

--PCC, October 2022
