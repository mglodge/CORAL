# CORAL
CORAL (Comparison Of Radiative AnaLyses) calculates the optical properties of particulate aerosols for use in atmospheric models. It is designed to compare three different optical scattering theories:

1) Mie scattering - assumes that particles are spherical
2) MMF (Modified Mean Field theory) - calculates the fractal properties of an input file, and uses these to calculate scattering properties
3) DDA (Discrete Dipole Approximation) - uses the exact shape file to determine the scattering properties.

The main outputs are the optical efficiencies $Q_{ext}$, $Q_{sca}$, $Q_{abs}$, and $g=<\cos(\theta)>$, as well as the associated cross-sections in $\mu m^2$, for a given wavelength range, aggregate radius range, shape-type, and refractive index profile.

All three methods (including details of equations) used in the code are explained in detail in the two papers available within the root folder:

- <b>PAPER 1:</b> "Aerosols are not Spherical Cows: Using Discrete Dipole Approximation to Model the Properties of Fractal Particles" (2023), Lodge et al. This paper includes a thorough explanation of CORAL and the methods used to compare Mie, MMF and DDA (original LDR method)
      
- <b>PAPER 2:</b> "MANTA-Ray: Supercharging Speeds for Calculating the Optical Properties of Fractal Aggregates in the Long-Wavelength Limit" (2024), Lodge et al. In section 2.2.2 of this paper, we describe an update to the DDA calculation in CORAL, which allows the user choice of LDR (the original DDA method) or the Filtered Coupled Dipole method of Yurkin et al. (2010).

- <b>NOTE</b>: Several small manuscript typos in the FCD theory were identified in Paper 2 by Maxim Yurkin and Clément Argentin, and a corrected version of the paper has now been uploaded here. In addition, there is a sign error in the FCD polarisability term (Eq 7+8), following from an error in the FCD implementation within DDSCAT. Fortunately the difference is negligible for the small-particle regime explored in this paper, and so the results presented in the manuscript are unchanged, but the code has now been corrected. Many thanks to Maxim Yurkin and Clément Argentin for finding these errors.

We are incredibly grateful to Kouji Adachi, Peter Buseck and Serena Chung for allowing the use of the shape.dat files included in the templates folder, which were obtained using the methods outlined in their paper ["Shapes of soot aerosol particles and implications for their effects on climate"](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009JD012868) (2010), JGR.

# Running the Code

1) Download all of the files into a directory of your choice -- these are already set up with default values to run as demonstration. There are 6 key files in this folder:

    - CORAL.c (the main code which we will compile)
    - shape.dat (the shape file that describes the aerosol)
    - parameters.txt (sets the parameters for the scattering models)
    - refractive_index_data.txt (describes the wavelength dependence of refractive index)
    - STAG.py (allows users to view the shapes interactively in 3D)
    - Shape Templates (this includes the three shapes used in Lodge et al. (2023))

2) To compile, simply navigate to the directory from within the terminal window, and type: 

         gcc -o CORAL CORAL.c -fopenmp -lm -O3

The flag -fopenmp enables OpenMP and parallel processing, which rapidly speed up the code, and -O3 gives a further increased speed boost. OpenMP is usually necessary 
to enable calculations at a reasonable speed, but users can remove the '#include omp.h' header from within CORAL.c, and remove the '-fopenmp' flag from above if they 
do not want to use it. Further details about OpenMP are [here.](https://www.openmp.org/)

3) To run the code, simply type:

       ./CORAL

# Parameters file

The parameters file details the range of wavelengths to study and the size of the particles etc. There are several options for the user to choose, and these are detailed below.

<b>IMPORTANT NOTE 1: the formatting has to be kept <i>exactly</i> the same as the original file i.e. if you change "shapeswitch= 1" to "shapeswitch=1" (if you delete the space after the equals sign), the code will not work. </b>

<b>IMPORTANT NOTE 2: CORAL v2.0 includes an option to use the filtered coupled dipole method, and so the old parameter files from v1.0 can no longer be used. Be sure that parameter files use the updated structure below.</b>

    shapeswitch: set this to 1 to use shape data file as an input, or set to 0 to make our own spheres with requested N (see below)
    marbleswitch: set to 1 to use Biconjugate Method with Stabilisation method of solving linear equations, or 0 to use standard conjugate gradient method. In our experience, the latter is faster.
    filtered_coupled_dipole_switch: set to 1 to use filtered coupled dipole method from Yurkin et al. (2010), or 0 to use LDR method from Draine (1994)
    accuracy: this sets the accuracy required for three consecutive iterations of the equation solver to agree and determine that they have converged. We recommend a value of 0.00001, representing (0.001% agreement) 
    aerosol_radius_min: set the minimum radius of aerosols to be examined, in um
    aerosol_radius_max= set the maximum radius of aerosols to be examined, in um
    n_radii: sets the number of radii to be analysed within the range above (max-min)
    radius_increment_type: set equal to 0 for linear-spaced incrememnts, or 1 for log-spaced increments
    requested_N: chooses a rough number of dipoles to create pseudospheres out of, if making our own spheres using shapeswitch=0 (rather than an input file)
    speed_test: Set to 0 to ignore, or set to 1 to record the times that various parts in the code reach. Used for debugging/speed benchmrking.
    
    /* wavelength range */
    wavelength_min: the minimum wavelength of radiation to be examined, in um
    wavelength_max: the maximum wavelength of radiation to be examined, in um
    num_wavelengths: sets the number of wavelengths to be analysed within the range above (max-min)
    wavelength_increment_type= 1  set equal to 0 for linear-spaced incrememnts, or 1 for log-spaced increments
    
    /* resolution that determines the number of angles that will be used */
    orientational_average_resolution: determines how many angles to fire incoming radiation from (see copy of Fig. C1 in Lodge_et_al_2023.pdf). For example, a value of 0 corresponds to 12 equally-spaced angles, a value of 1 corresponds to 42 angles, a value of 2 corresponds to 92 angles and so on... 
    
    threadnumber: choose the number of threads to distribute the parallelised calculations across
    
    N_monomers: set the number of monomers for the shape that you are analysing (often you can judge this by eye using STAG)
    MMF_method: select the cut-off function for the structure factor calculation a value of 0 gives a gaussian cutoff and a value of 1 gives a fractal cutoff. We use the fractal cutoff as default to match OPTOOL. For more details, see Eq. (19) of Tazaki, R. and Tanaka, H., 2018. Light scattering by fractal dust aggregates. II. opacity and asymmetry parameter. The Astrophysical Journal, 860(1), p.79.

Fig C1 from Lodge et al. (2023), with example angles for incoming radiation marked:

  - "orientational_average_resolution= 0" would assess 12 angles
  - "orientational_average_resolution= 1" would assess 42 angles
  - "orientational_average_resolution= 2" would assess 92 angles
  - "orientational_average_resolution= 5" would assess 362 angles, etc.

![angle_positions](https://github.com/mglodge/CORAL/assets/126600438/d866dd7f-af26-4723-9f41-f09e3f71db68)

# Shape templates

<b>IMPORTANT NOTE: When changing shape files, remember to state how many monomers the aggregate is composed of as the 'N_monomers' variable in "parameters.txt".</b>

Three demonstration shape files used in Lodge et al. (2023) are included, but any DDSCAT-style shape file will work.

   - compact cluster (N_monomers = 21)
   - linear branched fractal (N_monomers = 76)
   - elongated cluster (N_monomers = 15)

These three shapes are shown (in order, left to right) below. For more detail about how N_monomers is determined, please read section 3.3 of Lodge et al. (2023).

![realistic_fractals](https://github.com/mglodge/CORAL/assets/126600438/9398e6be-4253-4d91-9fed-8b1accdbd2b7)

# STAG

To view your 3D model in STAG (assuming that python is intalled on your device), simply type the following into the terminal:

    python STAG.py

This should produce a rotatable interactive figure, similar to the shapes shown above.

# Output files

If left at the default settings, CORAL will now calculate scattering properties for IR radiation of wavelength 1.5um incident on a linear branched fractal of radius 0.5 um, 
calculated over an average of 12 different directions and two perpendicular polarisation states. Several output files will appear:

MAIN OUTPUTS: 

Each of the following files finds $Q_{ext}$, $Q_{sca}$, $Q_{abs}$, and $g=<\cos(\theta)>$, as well as the associated cross-sections in $\mu m^2$, for each of the Mie, MMF and DDA models, at each particle radius and wavelength.

   - mie_results.txt
   - mmf_results.txt
   - DDA_results.txt

OTHER FILES:

The other files may be of use to those that wish to do further detailed analysis.

   - angles.txt (the list of incident angles that the EM radiation came from, in spherical polar coordinates defined by ISO 80000-2:2019)
   - ddscat_shape.dat (creates a DDSCAT-type shape.dat file in the case that the user chooses to create pseudosphere (see 'shapeswitch' in parameter file)
   - dipole_coords.txt (creates a list of coordinates for STAG.py to read and display the shape in 3D)
   - LDR_criteria.txt (ignore this)
   - MMF_key_parameters.txt (outputs fractal parameters $d_f$, $k_0$, $N$ and $R_g$ for a given shape type for single-wavelength and particle size calculations).
   - qtable.txt (lists results in exactly the same format as DDSCAT for easier comparison)
   - refractive_index_data.refrind (outputs the refractive index in the format required by VIRGA for further calculations - see folder of the same name below)
   - results_angle_averages_y.txt (angle averages of the optical properties for each wavelength assessed, y-polarisation (0,1,0) only)
   - results_angle_averages_z.txt (angle averages of the optical properties for each wavelength assessed, z-polarisation (0,0,1) only)
   - results_angles_y.txt (optical properties for every angle for each wavelength assessed, y-polarisation (0,1,0) only)
   - results_angles_z.txt (optical properties for every angle for each wavelength assessed, z-polarisation (0,0,1) only)
   - VIRGA (outputs optical properties for all three models in the format required by VIRGA for further calculations/use in cloud models; see [VIRGA](https://github.com/natashabatalha/virga) for more details.

# Citations/Disclaimer

If CORAL is useful to your research, please cite Lodge et al. (2023). You are free to use and modify this code, but please attribute credit to the author (Matt Lodge) and also the authors below.

I would like to thank the authors of the following two codes, which were really well-written benchmark models. The detailed comments 
throughout were also indredibly helpful for debugging our own models, and we are grateful for their attention to detail:

    OPTOOL: Dominik, C., Min, M. & Tazaki, R. 2021, Optool, 1.9, Astrophysics Source Code Library, ascl:2104.010

    DDSCAT: Draine, B.T., & Flatau, P.J., "Discrete dipole approximation for scattering calculations", J. Opt. Soc. Am. A, 11, 1491-1499 (1994)

For the FCD implementation of DDA, please also cite M. Yurkin:

      Yurkin, M.A., Min, M. and Hoekstra, A.G., 2010. Application of the discrete dipole approximation to very large refractive indices: Filtered coupled dipoles revived. Physical Review E—Statistical, Nonlinear, and Soft Matter Physics, 82(3), p.036703.

Although every care has been taken to benchmark the code and test it under a wide range of conditions (see attached paper for details), 
this code is provided without a warantee of any kind. The author cannot be held liable or accept any resposibility for any claim, loss 
or damages that arise in connection with the use of this code.
