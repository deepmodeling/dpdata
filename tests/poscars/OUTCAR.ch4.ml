 vasp.6.3.0 20Jan22 (build May 08 2022 00:08:36) gamma-only                     
  
 executed on             LinuxIFC date 2022.05.12  14:11:53
 running   16 mpi-ranks, with    1 threads/rank
 distrk:  each k-point on   16 cores,    1 groups
 distr:  one band on NCORE=   4 cores,    4 groups


--------------------------------------------------------------------------------------------------------


 INCAR:
   SYSTEM = CH4
   NCORE = 4
   KGAMMA = .TRUE.
   KSPACING = 1.0
   ENCUT = 520.0
   PREC = Normal
   ISTART = 0
   LWAVE = .FALSE.
   LCHARG = .FALSE.
   ISMEAR = 0
   SIGMA = 0.05
   EDIFF = 1e-5
   LREAL = .FALSE.
   NELM = 200
   NELMIN = 4
   IBRION = 0
   MDALGO = 3
   ISIF = 3
   ALGO = Fast
   ISYM = 0
   TEBEG = 10
   TEEND = 10
   NSW = 10
   POTIM = 1
   LANGEVIN_GAMMA = 1 1
   LANGEVIN_GAMMA_L = 1
   PMASS = 10
   PSTRESS = 1.0
   NWRITE = 2
   NBLOCK = 1
   ML_LMLFF = .TRUE.
   ML_ISTART = 0
   ML_MB = 2000
   ML_WTOTEN = 10
   ML_WTIFOR = 10
   RANDOM_SEED = 283862281                0                0

 POTCAR:   PAW_GGA H 07Jul1998                    
 POTCAR:   PAW_GGA C 05Jan2001                    
 POTCAR:   PAW_GGA H 07Jul1998                    
   VRHFIN =H: ultrasoft test                                                                        
   LEXCH  = 91                                                                                      
   EATOM  =    12.5313 eV,     .9210 Ry                                                             
                                                                                                    
   TITEL  = PAW_GGA H 07Jul1998                                                                     
   LULTRA =        F    use ultrasoft PP ?                                                          
   IUNSCR =        0    unscreen: 0-lin 1-nonlin 2-no                                               
   RPACOR =     .000    partial core radius                                                         
   POMASS =    1.000; ZVAL   =    1.000    mass and valenz                                          
   RCORE  =    1.100    outmost cutoff radius                                                       
   RWIGS  =     .700; RWIGS  =     .370    wigner-seitz radius (au A)                               
   ENMAX  =  250.000; ENMIN  =  200.000 eV                                                          
   RCLOC  =     .701    cutoff for local pot                                                        
   LCOR   =        T    correct aug charges                                                         
   LPAW   =        T    paw PP                                                                      
   EAUG   =  400.000                                                                                
   RMAX   =    2.174    core radius for proj-oper                                                   
   RAUG   =    1.200    factor for augmentation sphere                                              
   RDEP   =    1.112    core radius for depl-charge                                                 
   QCUT   =   -5.749; QGAM   =   11.498    optimization parameters                                  
                                                                                                    
   Description                                                                                      
     l     E      TYP  RCUT    TYP  RCUT                                                            
     0   .000     23  1.100                                                                         
     0   .500     23  1.100                                                                         
     1  -.300     23  1.100                                                                         
  local pseudopotential read in
  atomic valenz-charges read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
    PAW grid and wavefunctions read in
 
   number of l-projection  operators is LMAX  =           3
   number of lm-projection operators is LMMAX =           5
 
 POTCAR:   PAW_GGA C 05Jan2001                    
   VRHFIN =C: s2p2                                                                                  
   LEXCH  = 91                                                                                      
   EATOM  =   147.4688 eV,   10.8386 Ry                                                             
                                                                                                    
   TITEL  = PAW_GGA C 05Jan2001                                                                     
   LULTRA =        F    use ultrasoft PP ?                                                          
   IUNSCR =        0    unscreen: 0-lin 1-nonlin 2-no                                               
   RPACOR =     .000    partial core radius                                                         
   POMASS =   12.011; ZVAL   =    4.000    mass and valenz                                          
   RCORE  =    1.500    outmost cutoff radius                                                       
   RWIGS  =    1.630; RWIGS  =     .863    wigner-seitz radius (au A)                               
   ENMAX  =  400.000; ENMIN  =  300.000 eV                                                          
   ICORE  =        2    local potential                                                             
   LCOR   =        T    correct aug charges                                                         
   LPAW   =        T    paw PP                                                                      
   EAUG   =  644.873                                                                                
   DEXC   =     .000                                                                                
   RMAX   =    2.266    core radius for proj-oper                                                   
   RAUG   =    1.300    factor for augmentation sphere                                              
   RDEP   =    1.501    radius for radial grids                                                     
   RDEPT  =    1.300    core radius for aug-charge                                                  
   QCUT   =   -5.516; QGAM   =   11.033    optimization parameters                                  
                                                                                                    
   Description                                                                                      
     l     E      TYP  RCUT    TYP  RCUT                                                            
     0   .000     23  1.200                                                                         
     0   .000     23  1.200                                                                         
     1   .000     23  1.500                                                                         
     1  2.500     23  1.500                                                                         
     2   .000      7  1.500                                                                         
  local pseudopotential read in
  atomic valenz-charges read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
    PAW grid and wavefunctions read in
 
   number of l-projection  operators is LMAX  =           4
   number of lm-projection operators is LMMAX =           8
 
 PAW_GGA H 07Jul1998                    :
 energy of atom  1       EATOM=  -12.5313
 kinetic energy error for atom=    0.0014 (will be added to EATOM!!)
 PAW_GGA C 05Jan2001                    :
 energy of atom  2       EATOM= -147.4688
 kinetic energy error for atom=    0.0071 (will be added to EATOM!!)
 
 
 POSCAR: POSCAR file written by OVITO
  positions in direct lattice
  velocities in cartesian coordinates

  MD-specific parameters
                MDALGO =   3
        LANGEVIN_GAMMA =     1.000    1.000
      LANGEVIN_GAMMA_L =     1.000
                 CNEXP =     9.000   14.000
 exchange correlation table for  LEXCH =        7
   RHO(1)=    0.500       N(1)  =     2000
   RHO(2)=  100.500       N(2)  =     4000
 


--------------------------------------------------------------------------------------------------------


 ion  position               nearest neighbor table
   1  0.538  0.407  0.361-   5 1.10
   2  0.395  0.480  0.438-   5 1.10
   3  0.552  0.565  0.443-   5 1.10
   4  0.528  0.416  0.539-   5 1.10
   5  0.503  0.467  0.445-   1 1.10   4 1.10   3 1.10   2 1.10
 

IMPORTANT INFORMATION: All symmetrisations will be switched off!
NOSYMM: (Re-)initialisation of all symmetry stuff for point group C_1.


----------------------------------------------------------------------------------------

                                     Primitive cell                                     

  volume of cell :    1000.0000

  direct lattice vectors                    reciprocal lattice vectors
    10.000000000  0.000000000  0.000000000     0.100000000  0.000000000  0.000000000
     0.000000000 10.000000000  0.000000000     0.000000000  0.100000000  0.000000000
     0.000000000  0.000000000 10.000000000     0.000000000  0.000000000  0.100000000

  length of vectors
    10.000000000 10.000000000 10.000000000     0.100000000  0.100000000  0.100000000

  position of ions in fractional coordinates (direct lattice)
     0.538154339  0.406860801  0.360573014
     0.394539656  0.480320569  0.438468840
     0.552092428  0.565450285  0.442708738
     0.528185305  0.416414755  0.539182657
     0.503250593  0.467255160  0.445232341

  ion indices of the primitive-cell ions
   primitive index   ion index
                 1           1
                 2           2
                 3           3
                 4           4
                 5           5

----------------------------------------------------------------------------------------

 
 -----------------------------------------------------------------------------
|                                                                             |
|           W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!           |
|           W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!           |
|           W    W  A    A  R    R  N N  N  II  N N  N  G       !!!           |
|           W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !            |
|           WW  WW  A    A  R   R   N   NN  II  N   NN  G    G                |
|           W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!           |
|                                                                             |
|     The requested file  could not be found or opened for reading            |
|     k-point information. Automatic k-point generation is used as a          |
|     fallback, which may lead to unwanted results.                           |
|                                                                             |
 -----------------------------------------------------------------------------

 

Automatic generation of k-mesh.
 Grid dimensions derived from KSPACING:
 generate k-points for:    1    1    1

 Generating k-lattice:

  Cartesian coordinates                     Fractional coordinates (reciprocal lattice)
     0.100000000  0.000000000  0.000000000     1.000000000  0.000000000  0.000000000
     0.000000000  0.100000000  0.000000000     0.000000000  1.000000000  0.000000000
     0.000000000  0.000000000  0.100000000     0.000000000  0.000000000  1.000000000

  Length of vectors
     0.100000000  0.100000000  0.100000000

  Shift w.r.t. Gamma in fractional coordinates (k-lattice)
     0.000000000  0.000000000  0.000000000

 
 Subroutine IBZKPT returns following result:
 ===========================================
 
 Found      1 irreducible k-points:
 
 Following reciprocal coordinates:
            Coordinates               Weight
  0.000000  0.000000  0.000000      1.000000
 
 Following cartesian coordinates:
            Coordinates               Weight
  0.000000  0.000000  0.000000      1.000000
 


--------------------------------------------------------------------------------------------------------




 Dimension of arrays:
   k-points           NKPTS =      1   k-points in BZ     NKDIM =      1   number of bands    NBANDS=      8
   number of dos      NEDOS =    301   number of ions     NIONS =      5
   non local maximal  LDIM  =      4   non local SUM 2l+1 LMDIM =      8
   total plane-waves  NPLWV = 175616
   max r-space proj   IRMAX =      1   max aug-charges    IRDMAX=   2272
   dimension x,y,z NGX =    56 NGY =   56 NGZ =   56
   dimension x,y,z NGXF=   112 NGYF=  112 NGZF=  112
   support grid    NGXF=   112 NGYF=  112 NGZF=  112
   ions per type =               4   1
   NGX,Y,Z   is equivalent  to a cutoff of   9.31,  9.31,  9.31 a.u.
   NGXF,Y,Z  is equivalent  to a cutoff of  18.62, 18.62, 18.62 a.u.

 SYSTEM =  CH4                                     
 POSCAR =  POSCAR file written by OVITO            

 Startparameter for this run:
   NWRITE =      2    write-flag & timer
   PREC   = normal    normal or accurate (medium, high low for compatibility)
   ISTART =      0    job   : 0-new  1-cont  2-samecut
   ICHARG =      2    charge: 1-file 2-atom 10-const
   ISPIN  =      1    spin polarized calculation?
   LNONCOLLINEAR =      F non collinear calculations
   LSORBIT =      F    spin-orbit coupling
   INIWAV =      1    electr: 0-lowe 1-rand  2-diag
   LASPH  =      F    aspherical Exc in radial PAW
 Electronic Relaxation 1
   ENCUT  =  520.0 eV  38.22 Ry    6.18 a.u.  18.59 18.59 18.59*2*pi/ulx,y,z
   ENINI  =  520.0     initial cutoff
   ENAUG  =  644.9 eV  augmentation charge cutoff
   NELM   =    200;   NELMIN=  4; NELMDL= -5     # of ELM steps 
   EDIFF  = 0.1E-04   stopping-criterion for ELM
   LREAL  =      F    real-space projection
   NLSPLINE    = F    spline interpolate recip. space projectors
   LCOMPAT=      F    compatible to vasp.4.4
   GGA_COMPAT  = T    GGA compatible to vasp.4.4-vasp.4.6
   LMAXPAW     = -100 max onsite density
   LMAXMIX     =    2 max onsite mixed and CHGCAR
   VOSKOWN=      0    Vosko Wilk Nusair interpolation
   ROPT   =    0.00000   0.00000
 Ionic relaxation
   EDIFFG = 0.1E-03   stopping-criterion for IOM
   NSW    =     10    number of steps for IOM
   NBLOCK =      1;   KBLOCK =     10    inner block; outer block 
   IBRION =      0    ionic relax: 0-MD 1-quasi-New 2-CG
   NFREE  =      0    steps in history (QN), initial steepest desc. (CG)
   ISIF   =      3    stress and relaxation
   IWAVPR =     11    prediction:  0-non 1-charg 2-wave 3-comb
   ISYM   =      0    0-nonsym 1-usesym 2-fastsym
   LCORR  =      T    Harris-Foulkes like correction to forces

   POTIM  = 1.0000    time-step for ionic-motion
   TEIN   =    0.0    initial temperature
   TEBEG  =   10.0;   TEEND  =  10.0 temperature during run
   SMASS  =  -3.00    Nose mass-parameter (am)
   estimated Nose-frequenzy (Omega)   =  0.10E-29 period in steps = 0.63E+46 mass=  -0.229E-26a.u.
   SCALEE = 1.0000    scale energy and forces
   NPACO  =    256;   APACO  = 10.0  distance and # of slots for P.C.
   PSTRESS=    1.0 pullay stress

  Mass of Ions in am
   POMASS =   1.00 12.01
  Ionic Valenz
   ZVAL   =   1.00  4.00
  Atomic Wigner-Seitz radii
   RWIGS  =  -1.00 -1.00
  virtual crystal weights 
   VCA    =   1.00  1.00
   NELECT =       8.0000    total number of electrons
   NUPDOWN=      -1.0000    fix difference up-down

 DOS related values:
   EMIN   =  10.00;   EMAX   =-10.00  energy-range for DOS
   EFERMI =   0.00
   ISMEAR =     0;   SIGMA  =   0.05  broadening in eV -4-tet -1-fermi 0-gaus

 Electronic relaxation 2 (details)
   IALGO  =     68    algorithm
   LDIAG  =      T    sub-space diagonalisation (order eigenvalues)
   LSUBROT=      F    optimize rotation matrix (better conditioning)
   TURBO    =      0    0=normal 1=particle mesh
   IRESTART =      0    0=no restart 2=restart with 2 vectors
   NREBOOT  =      0    no. of reboots
   NMIN     =      0    reboot dimension
   EREF     =   0.00    reference energy to select bands
   IMIX   =      4    mixing-type and parameters
     AMIX     =   0.40;   BMIX     =  1.00
     AMIX_MAG =   1.60;   BMIX_MAG =  1.00
     AMIN     =   0.10
     WC   =   100.;   INIMIX=   1;  MIXPRE=   1;  MAXMIX= -45

 Intra band minimization:
   WEIMIN = 0.0010     energy-eigenvalue tresh-hold
   EBREAK =  0.31E-06  absolut break condition
   DEPER  =   0.30     relativ break condition  

   TIME   =   0.40     timestep for ELM

  volume/ion in A,a.u.               =     200.00      1349.67
  Fermi-wavevector in a.u.,A,eV,Ry     =   0.327420  0.618734  1.458594  0.107204
  Thomas-Fermi vector in A             =   1.220131
 
 Write flags
   LWAVE        =      F    write WAVECAR
   LDOWNSAMPLE  =      F    k-point downsampling of WAVECAR
   LCHARG       =      F    write CHGCAR
   LVTOT        =      F    write LOCPOT, total local potential
   LVHAR        =      F    write LOCPOT, Hartree potential only
   LELF         =      F    write electronic localiz. function (ELF)
   LORBIT       =      0    0 simple, 1 ext, 2 COOP (PROOUT), +10 PAW based schemes


 Dipole corrections
   LMONO  =      F    monopole corrections only (constant potential shift)
   LDIPOL =      F    correct potential (dipole corrections)
   IDIPOL =      0    1-x, 2-y, 3-z, 4-all directions 
   EPSILON=  1.0000000 bulk dielectric constant

 Exchange correlation treatment:
   GGA     =    --    GGA type
   LEXCH   =     7    internal setting for exchange type
   LIBXC   =     F    Libxc                    
   VOSKOWN =     0    Vosko Wilk Nusair interpolation
   LHFCALC =     F    Hartree Fock is set to
   LHFONE  =     F    Hartree Fock one center treatment
   AEXX    =    0.0000 exact exchange contribution

 Linear response parameters
   LEPSILON=     F    determine dielectric tensor
   LRPA    =     F    only Hartree local field effects (RPA)
   LNABLA  =     F    use nabla operator in PAW spheres
   LVEL    =     F    velocity operator in full k-point grid
   CSHIFT  =0.1000    complex shift for real part using Kramers Kronig
   OMEGAMAX=  -1.0    maximum frequency
   DEG_THRESHOLD= 0.2000000E-02 threshold for treating states as degnerate
   RTIME   =   -0.100 relaxation time in fs
  (WPLASMAI=    0.000 imaginary part of plasma frequency in eV, 0.658/RTIME)
   DFIELD  = 0.0000000 0.0000000 0.0000000 field for delta impulse in time
 
  Optional k-point grid parameters
   LKPOINTS_OPT  =     F    use optional k-point grid
   KPOINTS_OPT_MODE=     1    mode for optional k-point grid
 
 Orbital magnetization related:
   ORBITALMAG=     F  switch on orbital magnetization
   LCHIMAG   =     F  perturbation theory with respect to B field
   DQ        =  0.001000  dq finite difference perturbation B field
   LLRAUG    =     F  two centre corrections for induced B field



--------------------------------------------------------------------------------------------------------


 molecular dynamics for ions
   using a microcanonical ensemble
 charge density and potential will be updated during run
 non-spin polarized calculation
 RMM-DIIS sequential band-by-band and
  variant of blocked Davidson during initial phase
 perform sub-space diagonalisation
    before iterative eigenvector-optimisation
 modified Broyden-mixing scheme, WC =      100.0
 initial mixing is a Kerker type mixing with AMIX =  0.4000 and BMIX =      1.0000
 Hartree-type preconditioning will be used
 using additional bands            4
 reciprocal scheme for non local part
 calculate Harris-corrections to forces 
   (improved forces if not selfconsistent)
 use gradient corrections 
 use of overlap-Matrix (Vanderbilt PP)
 Gauss-broadening in eV      SIGMA  =   0.05


--------------------------------------------------------------------------------------------------------


  energy-cutoff  :      520.00
  volume of cell :     1000.00
      direct lattice vectors                 reciprocal lattice vectors
    10.000000000  0.000000000  0.000000000     0.100000000  0.000000000  0.000000000
     0.000000000 10.000000000  0.000000000     0.000000000  0.100000000  0.000000000
     0.000000000  0.000000000 10.000000000     0.000000000  0.000000000  0.100000000

  length of vectors
    10.000000000 10.000000000 10.000000000     0.100000000  0.100000000  0.100000000


 
 k-points in units of 2pi/SCALE and weight: read from INCAR                         
   0.00000000  0.00000000  0.00000000       1.000
 
 k-points in reciprocal lattice and weights: read from INCAR                         
   0.00000000  0.00000000  0.00000000       1.000
 
 position of ions in fractional coordinates (direct lattice) 
   0.53815434  0.40686080  0.36057301
   0.39453966  0.48032057  0.43846884
   0.55209243  0.56545029  0.44270874
   0.52818530  0.41641476  0.53918266
   0.50325059  0.46725516  0.44523234
 
 position of ions in cartesian coordinates  (Angst):
   5.38154339  4.06860801  3.60573014
   3.94539656  4.80320569  4.38468840
   5.52092428  5.65450285  4.42708738
   5.28185305  4.16414755  5.39182657
   5.03250593  4.67255160  4.45232341
 


--------------------------------------------------------------------------------------------------------


 use parallel FFT for orbitals z direction half grid
 k-point   1 :   0.0000 0.0000 0.0000  plane waves:   13481

 maximum and minimum number of plane-waves per node :      3374     3363

 maximum number of plane-waves:     13481
 maximum index in each direction: 
   IXMAX=   18   IYMAX=   18   IZMAX=   18
   IXMIN=  -18   IYMIN=  -18   IZMIN=    0


 parallel 3D FFT for wavefunctions:
    minimum data exchange during FFTs selected (reduces bandwidth)
 parallel 3D FFT for charge:
    minimum data exchange during FFTs selected (reduces bandwidth)


 total amount of memory used by VASP MPI-rank0    41513. kBytes
=======================================================================

   base      :      30000. kBytes
   nonl-proj :        699. kBytes
   fftplans  :       3956. kBytes
   grid      :       6748. kBytes
   one-center:          3. kBytes
   wavefun   :        107. kBytes
 
     INWAV:  cpu time      0.0000: real time      0.0001
 Broyden mixing: mesh for mixing (old mesh)
   NGX = 37   NGY = 37   NGZ = 37
  (NGX  =112   NGY  =112   NGZ  =112)
  gives a total of  50653 points

 initial charge density was supplied:
 charge density of overlapping atoms calculated
 number of electron       8.0000000 magnetization 
 keeping initial charge density in first step


--------------------------------------------------------------------------------------------------------


 Maximum index for augmentation-charges           90 (set IRDMAX)


--------------------------------------------------------------------------------------------------------


 First call to EWALD:  gamma=   0.177
 Maximum number of real-space cells 3x 3x 3
 Maximum number of reciprocal cells 3x 3x 3

    FEWALD:  cpu time      0.0013: real time      0.0015
    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:      0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  in kB       0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  external pressure =       -1.00 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =      0.00 kB
  Total+kin.     0.000       0.000       0.000       0.000       0.000       0.000
  volume of cell :     1000.00
      direct lattice vectors                 reciprocal lattice vectors
    10.000000000  0.000000000  0.000000000     0.100000000  0.000000000  0.000000000
     0.000000000 10.000000000  0.000000000     0.000000000  0.100000000  0.000000000
     0.000000000  0.000000000 10.000000000     0.000000000  0.000000000  0.100000000

  length of vectors
    10.000000000 10.000000000 10.000000000     0.100000000  0.100000000  0.100000000


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38154      4.06861      3.60573         0.000000      0.000000      0.000000
      3.94540      4.80321      4.38469         0.000000      0.000000      0.000000
      5.52092      5.65450      4.42709         0.000000      0.000000      0.000000
      5.28185      4.16415      5.39183         0.000000      0.000000      0.000000
      5.03251      4.67255      4.45232         0.000000      0.000000      0.000000
 -----------------------------------------------------------------------------------
    total drift:                                0.000000      0.000000      0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =         0.00000000 eV

  ML energy  without entropy=        0.00000000  ML energy(sigma->0) =        0.00000000

  enthalpy is ML TOTEN    =         0.62415064 eV   P V=        0.62415064



--------------------------------------- Iteration      1(   1)  ---------------------------------------


    POTLOK:  cpu time      0.0594: real time      0.0680
    SETDIJ:  cpu time      0.0020: real time      0.0044
     EDDAV:  cpu time      0.0231: real time      0.0316
       DOS:  cpu time      0.0009: real time      0.0027
    --------------------------------------------
      LOOP:  cpu time      0.0854: real time      0.1066

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.6999765E+02  (-0.1857416E+03)
 number of electron       8.0000000 magnetization 
 augmentation part        8.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -244.12240337
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        22.64847752
  PAW double counting   =        41.11719426      -42.35612335
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -33.79259455
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =        69.99764688 eV

  energy without entropy =       69.99764688  energy(sigma->0) =       69.99764688


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   2)  ---------------------------------------


     EDDAV:  cpu time      0.0261: real time      0.0261
       DOS:  cpu time      0.0001: real time      0.0001
    --------------------------------------------
      LOOP:  cpu time      0.0262: real time      0.0262

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.6232328E+02  (-0.6232328E+02)
 number of electron       8.0000000 magnetization 
 augmentation part        8.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -244.12240337
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        22.64847752
  PAW double counting   =        41.11719426      -42.35612335
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -96.11587390
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =         7.67436753 eV

  energy without entropy =        7.67436753  energy(sigma->0) =        7.67436753


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   3)  ---------------------------------------


     EDDAV:  cpu time      0.0146: real time      0.0148
       DOS:  cpu time      0.0001: real time      0.0001
    --------------------------------------------
      LOOP:  cpu time      0.0147: real time      0.0149

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.3181075E+02  (-0.3181075E+02)
 number of electron       8.0000000 magnetization 
 augmentation part        8.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -244.12240337
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        22.64847752
  PAW double counting   =        41.11719426      -42.35612335
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -127.92662209
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.13638066 eV

  energy without entropy =      -24.13638066  energy(sigma->0) =      -24.13638066


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   4)  ---------------------------------------


     EDDAV:  cpu time      0.0145: real time      0.0147
       DOS:  cpu time      0.0001: real time      0.0001
    --------------------------------------------
      LOOP:  cpu time      0.0146: real time      0.0148

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.3413468E+01  (-0.3413468E+01)
 number of electron       8.0000000 magnetization 
 augmentation part        8.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -244.12240337
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        22.64847752
  PAW double counting   =        41.11719426      -42.35612335
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -131.34009031
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -27.54984888 eV

  energy without entropy =      -27.54984888  energy(sigma->0) =      -27.54984888


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   5)  ---------------------------------------


     EDDAV:  cpu time      0.0247: real time      0.0248
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0064: real time      0.0081
    MIXING:  cpu time      0.0041: real time      0.0059
    --------------------------------------------
      LOOP:  cpu time      0.0354: real time      0.0389

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.1233522E+00  (-0.1233522E+00)
 number of electron       8.0000007 magnetization 
 augmentation part        0.2323676 magnetization 

 Broyden mixing:
  rms(total) = 0.84942E+00    rms(broyden)= 0.84931E+00
  rms(prec ) = 0.12300E+01
  weight for this iteration     100.00

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -244.12240337
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        22.64847752
  PAW double counting   =        41.11719426      -42.35612335
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -131.46344255
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -27.67320112 eV

  energy without entropy =      -27.67320112  energy(sigma->0) =      -27.67320112


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   6)  ---------------------------------------


    POTLOK:  cpu time      0.0581: real time      0.0582
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0124: real time      0.0219
  RMM-DIIS:  cpu time      0.0139: real time      0.0147
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0053: real time      0.0053
    MIXING:  cpu time      0.0027: real time      0.0035
    --------------------------------------------
      LOOP:  cpu time      0.0936: real time      0.1048

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.2944242E+01  (-0.4748563E+00)
 number of electron       8.0000007 magnetization 
 augmentation part        0.2052833 magnetization 

 Broyden mixing:
  rms(total) = 0.37023E+00    rms(broyden)= 0.37021E+00
  rms(prec ) = 0.51924E+00
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.3828
  1.3828

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -263.92806877
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        23.86234804
  PAW double counting   =        69.00369650      -70.45380803
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -109.71622293
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.72895883 eV

  energy without entropy =      -24.72895883  energy(sigma->0) =      -24.72895883


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   7)  ---------------------------------------


    POTLOK:  cpu time      0.0593: real time      0.0595
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0034: real time      0.0034
  RMM-DIIS:  cpu time      0.0132: real time      0.0132
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0018: real time      0.0018
    --------------------------------------------
      LOOP:  cpu time      0.0839: real time      0.0841

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.4217318E+00  (-0.1875670E+00)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1908460 magnetization 

 Broyden mixing:
  rms(total) = 0.20477E+00    rms(broyden)= 0.20472E+00
  rms(prec ) = 0.25962E+00
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.0981
  1.7106  2.4856

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -276.64447254
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        24.76761791
  PAW double counting   =        86.49593432      -88.01099629
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -97.41840673
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.30722698 eV

  energy without entropy =      -24.30722698  energy(sigma->0) =      -24.30722698


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   8)  ---------------------------------------


    POTLOK:  cpu time      0.0556: real time      0.0557
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0131: real time      0.0131
    ORTHCH:  cpu time      0.0003: real time      0.0003
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0049: real time      0.0049
    MIXING:  cpu time      0.0021: real time      0.0021
    --------------------------------------------
      LOOP:  cpu time      0.0801: real time      0.0803

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.1493866E+00  (-0.4073743E-01)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1942004 magnetization 

 Broyden mixing:
  rms(total) = 0.80159E-01    rms(broyden)= 0.80158E-01
  rms(prec ) = 0.11958E+00
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.5868
  2.4490  0.8688  1.4426

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -279.85430287
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.22242630
  PAW double counting   =        86.90836895      -88.31856446
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -94.61886465
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.15784036 eV

  energy without entropy =      -24.15784036  energy(sigma->0) =      -24.15784036


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   9)  ---------------------------------------


    POTLOK:  cpu time      0.0559: real time      0.0560
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0129: real time      0.0130
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0049: real time      0.0049
    MIXING:  cpu time      0.0020: real time      0.0020
    --------------------------------------------
      LOOP:  cpu time      0.0803: real time      0.0804

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.2753904E-01  (-0.3356277E-02)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1933423 magnetization 

 Broyden mixing:
  rms(total) = 0.49990E-01    rms(broyden)= 0.49990E-01
  rms(prec ) = 0.76800E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.5906
  1.8279  1.8279  1.3533  1.3533

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -281.34977031
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.29664039
  PAW double counting   =        88.88104200      -90.31680163
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -93.14450813
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.13030132 eV

  energy without entropy =      -24.13030132  energy(sigma->0) =      -24.13030132


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  10)  ---------------------------------------


    POTLOK:  cpu time      0.0556: real time      0.0560
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0130: real time      0.0130
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0049: real time      0.0049
    MIXING:  cpu time      0.0023: real time      0.0023
    --------------------------------------------
      LOOP:  cpu time      0.0802: real time      0.0806

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.5752872E-02  (-0.2159009E-01)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1896610 magnetization 

 Broyden mixing:
  rms(total) = 0.64119E-01    rms(broyden)= 0.64115E-01
  rms(prec ) = 0.85062E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.7914
  2.9369  2.4806  0.9636  1.2878  1.2878

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -283.96361568
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.35709251
  PAW double counting   =        90.61238336      -92.09095073
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.55406001
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.13605419 eV

  energy without entropy =      -24.13605419  energy(sigma->0) =      -24.13605419


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  11)  ---------------------------------------


    POTLOK:  cpu time      0.0529: real time      0.0531
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0132: real time      0.0132
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0049: real time      0.0049
    MIXING:  cpu time      0.0023: real time      0.0023
    --------------------------------------------
      LOOP:  cpu time      0.0778: real time      0.0779

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.2374416E-01  (-0.2168259E-02)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1904599 magnetization 

 Broyden mixing:
  rms(total) = 0.13585E-01    rms(broyden)= 0.13585E-01
  rms(prec ) = 0.19949E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.8785
  3.9503  2.1794  1.6771  0.9410  1.2617  1.2617

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -284.36569701
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.40602284
  PAW double counting   =        89.08931980      -90.54081604
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.20423597
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.11231002 eV

  energy without entropy =      -24.11231002  energy(sigma->0) =      -24.11231002


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  12)  ---------------------------------------


    POTLOK:  cpu time      0.0505: real time      0.0506
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0129: real time      0.0129
    ORTHCH:  cpu time      0.0003: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0049: real time      0.0049
    MIXING:  cpu time      0.0023: real time      0.0023
    --------------------------------------------
      LOOP:  cpu time      0.0752: real time      0.0753

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.5616810E-02  (-0.3038805E-03)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1903692 magnetization 

 Broyden mixing:
  rms(total) = 0.91223E-02    rms(broyden)= 0.91222E-02
  rms(prec ) = 0.15608E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.1046
  5.0953  2.8604  2.1934  1.3197  1.3197  0.9717  0.9717

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -285.08904445
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.44199507
  PAW double counting   =        89.62931534      -91.08054492
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -89.52274424
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.11792683 eV

  energy without entropy =      -24.11792683  energy(sigma->0) =      -24.11792683


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  13)  ---------------------------------------


    POTLOK:  cpu time      0.0492: real time      0.0494
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0128: real time      0.0128
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0024: real time      0.0024
    --------------------------------------------
      LOOP:  cpu time      0.0740: real time      0.0742

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.6372986E-02  (-0.9410624E-03)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1907689 magnetization 

 Broyden mixing:
  rms(total) = 0.43509E-02    rms(broyden)= 0.43500E-02
  rms(prec ) = 0.60487E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.1515
  4.9398  3.4568  2.1253  2.1253  1.3091  1.3091  0.9935  0.9531

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -284.52848187
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.39819219
  PAW double counting   =        89.09346772      -90.53922434
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.05134988
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12429982 eV

  energy without entropy =      -24.12429982  energy(sigma->0) =      -24.12429982


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  14)  ---------------------------------------


    POTLOK:  cpu time      0.0493: real time      0.0494
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0034: real time      0.0034
  RMM-DIIS:  cpu time      0.0129: real time      0.0129
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0026: real time      0.0026
    --------------------------------------------
      LOOP:  cpu time      0.0743: real time      0.0744

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.3192683E-02  (-0.7624847E-03)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1912804 magnetization 

 Broyden mixing:
  rms(total) = 0.13606E-01    rms(broyden)= 0.13605E-01
  rms(prec ) = 0.19348E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.3626
  6.5374  4.3140  2.5610  2.1865  1.3075  1.3075  1.1135  1.0034  0.9324

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -284.15389486
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.37442709
  PAW double counting   =        88.56511515      -90.00909297
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.40714328
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12749250 eV

  energy without entropy =      -24.12749250  energy(sigma->0) =      -24.12749250


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  15)  ---------------------------------------


    POTLOK:  cpu time      0.0483: real time      0.0484
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0128: real time      0.0128
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0049: real time      0.0049
    MIXING:  cpu time      0.0028: real time      0.0028
    --------------------------------------------
      LOOP:  cpu time      0.0734: real time      0.0734

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.7275167E-03  (-0.3948155E-03)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1908414 magnetization 

 Broyden mixing:
  rms(total) = 0.17697E-02    rms(broyden)= 0.17695E-02
  rms(prec ) = 0.24022E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.3085
  6.9182  4.1461  2.5704  1.9486  1.9486  1.3118  1.3118  0.9419  0.9419  1.0452

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -284.68181353
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.39948386
  PAW double counting   =        89.28729910      -90.73777251
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -89.89705826
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12676499 eV

  energy without entropy =      -24.12676499  energy(sigma->0) =      -24.12676499


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  16)  ---------------------------------------


    POTLOK:  cpu time      0.0480: real time      0.0492
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0130: real time      0.0130
    ORTHCH:  cpu time      0.0003: real time      0.0003
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0048: real time      0.0048
    MIXING:  cpu time      0.0031: real time      0.0031
    --------------------------------------------
      LOOP:  cpu time      0.0735: real time      0.0747

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.1791412E-03  (-0.2065426E-04)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1908472 magnetization 

 Broyden mixing:
  rms(total) = 0.23969E-03    rms(broyden)= 0.23940E-03
  rms(prec ) = 0.53105E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.4806
  7.7572  5.1870  2.8325  2.8325  2.0944  1.3182  1.3182  1.0757  0.9696  0.9507
  0.9507

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -284.64325100
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.39810953
  PAW double counting   =        89.20067227      -90.65026209
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -89.93530920
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12694413 eV

  energy without entropy =      -24.12694413  energy(sigma->0) =      -24.12694413


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  17)  ---------------------------------------


    POTLOK:  cpu time      0.0496: real time      0.0498
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0129: real time      0.0129
    ORTHCH:  cpu time      0.0003: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0030: real time      0.0030
    --------------------------------------------
      LOOP:  cpu time      0.0750: real time      0.0752

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.1247330E-03  (-0.5936357E-05)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1908195 magnetization 

 Broyden mixing:
  rms(total) = 0.81358E-03    rms(broyden)= 0.81346E-03
  rms(prec ) = 0.10332E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.4194
  8.0154  5.1811  3.4343  2.5600  2.0070  1.3182  1.3182  1.3149  1.0597  0.9315
  0.9461  0.9461

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -284.63349018
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.39651980
  PAW double counting   =        89.19533399      -90.64480882
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -89.94372002
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12706886 eV

  energy without entropy =      -24.12706886  energy(sigma->0) =      -24.12706886


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  18)  ---------------------------------------


    POTLOK:  cpu time      0.0508: real time      0.0510
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0104: real time      0.0104
    ORTHCH:  cpu time      0.0003: real time      0.0003
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0032: real time      0.0032
    --------------------------------------------
      LOOP:  cpu time      0.0739: real time      0.0740

 eigenvalue-minimisations  :    12
 total energy-change (2. order) :-0.2014560E-04  (-0.7436772E-06)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1908158 magnetization 

 Broyden mixing:
  rms(total) = 0.51177E-03    rms(broyden)= 0.51174E-03
  rms(prec ) = 0.66915E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.5131
  8.3362  5.3697  3.7299  2.8780  2.5903  2.0916  1.3143  1.3143  1.1182  1.0752
  0.9535  0.9494  0.9494

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -284.63831349
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.39701494
  PAW double counting   =        89.19004721      -90.63941212
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -89.93952191
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12708901 eV

  energy without entropy =      -24.12708901  energy(sigma->0) =      -24.12708901


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  19)  ---------------------------------------


    POTLOK:  cpu time      0.0532: real time      0.0533
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0104: real time      0.0104
    ORTHCH:  cpu time      0.0003: real time      0.0003
       DOS:  cpu time      0.0001: real time      0.0001
    --------------------------------------------
      LOOP:  cpu time      0.0681: real time      0.0682

 eigenvalue-minimisations  :    12
 total energy-change (2. order) :-0.9011489E-05  (-0.5561172E-06)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1908158 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23919328
  Ewald energy   TEWEN  =       128.68270354
  -Hartree energ DENC   =      -284.62523194
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.39657989
  PAW double counting   =        89.17924264      -90.62885322
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -89.95193175
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12709802 eV

  energy without entropy =      -24.12709802  energy(sigma->0) =      -24.12709802


--------------------------------------------------------------------------------------------------------




 average (electrostatic) potential at core
  the test charge radii are     0.5201  0.6991
  (the norm of the test charge is              1.0000)
       1 -41.5589       2 -41.5565       3 -41.5570       4 -41.5581       5 -59.0029
 
 
 
 E-fermi :  -9.1045     XC(G=0):  -0.6163     alpha+bet : -0.2111

 Fermi energy:        -9.1044786619

 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -16.9228      2.00000
      2      -9.3517      2.00000
      3      -9.3512      2.00000
      4      -9.3503      2.00000
      5      -0.5421      0.00000
      6       1.0562      0.00000
      7       1.0670      0.00000
      8       1.0822      0.00000


--------------------------------------------------------------------------------------------------------


 soft charge-density along one line, spin component           1
         0         1         2         3         4         5         6         7         8         9
 total charge-density along one line
 
 pseudopotential strength for first ion, spin component:           1
 -2.357  -0.011  -0.001  -0.002   0.001
 -0.011   0.050  -0.001  -0.002   0.001
 -0.001  -0.001  -0.339   0.002  -0.001
 -0.002  -0.002   0.002  -0.338  -0.001
  0.001   0.001  -0.001  -0.001  -0.340
 total augmentation occupancy for first ion, spin component:           1
  1.738  -0.485   0.166   0.233  -0.096
 -0.485   0.185  -0.053  -0.074   0.030
  0.166  -0.053   0.023   0.020  -0.008
  0.233  -0.074   0.020   0.036  -0.011
 -0.096   0.030  -0.008  -0.011   0.014


------------------------ aborting loop because EDIFF is reached ----------------------------------------


    CHARGE:  cpu time      0.0049: real time      0.0049
    FORLOC:  cpu time      0.0037: real time      0.0037
    FORNL :  cpu time      0.0019: real time      0.0033
    STRESS:  cpu time      0.0221: real time      0.0231
    FORHAR:  cpu time      0.0134: real time      0.0141
    MIXING:  cpu time      0.0037: real time      0.0037
    OFIELD:  cpu time      0.0001: real time      0.0001

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Alpha Z     0.23919     0.23919     0.23919
  Ewald      42.88200    42.89495    42.90570    -0.00850     0.00232    -0.00338
  Hartree    94.87523    94.87841    94.88063    -0.00495     0.00061    -0.00128
  E(xc)     -27.05613   -27.05607   -27.05607    -0.00000     0.00001    -0.00004
  Local    -196.91996  -196.93284  -196.94306     0.01270    -0.00232     0.00407
  n-local   -15.55573   -15.55569   -15.55571     0.00001    -0.00001    -0.00006
  augment     0.19051     0.19083     0.19060     0.00003    -0.00006    -0.00005
  Kinetic   101.17913   101.17641   101.17459     0.00035    -0.00054     0.00055
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  -------------------------------------------------------------------------------------
  Total      -0.16576    -0.16480    -0.16412    -0.00035     0.00000    -0.00019
  in kB      -0.26557    -0.26404    -0.26295    -0.00057     0.00000    -0.00031
  external pressure =       -1.26 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =     -0.26 kB
  Total+kin.    -0.266      -0.264      -0.263      -0.001       0.000      -0.000

 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
  energy-cutoff  :      520.00
  volume of cell :     1000.00
      direct lattice vectors                 reciprocal lattice vectors
    10.000000000  0.000000000  0.000000000     0.100000000  0.000000000  0.000000000
     0.000000000 10.000000000  0.000000000     0.000000000  0.100000000  0.000000000
     0.000000000  0.000000000 10.000000000     0.000000000  0.000000000  0.100000000

  length of vectors
    10.000000000 10.000000000 10.000000000     0.100000000  0.100000000  0.100000000


 FORCES acting on ions
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction
 -----------------------------------------------------------------------------------------------
   -.168E+02 0.291E+02 0.408E+02   0.186E+02 -.321E+02 -.450E+02   -.177E+01 0.306E+01 0.429E+01   0.352E-03 -.646E-03 -.928E-03
   0.524E+02 -.630E+01 0.326E+01   -.578E+02 0.695E+01 -.360E+01   0.551E+01 -.662E+00 0.342E+00   -.114E-02 0.154E-03 -.739E-04
   -.235E+02 -.473E+02 0.122E+01   0.260E+02 0.522E+02 -.134E+01   -.247E+01 -.497E+01 0.128E+00   0.486E-03 0.106E-02 -.290E-04
   -.120E+02 0.245E+02 -.453E+02   0.133E+02 -.270E+02 0.500E+02   -.126E+01 0.257E+01 -.476E+01   0.245E-03 -.539E-03 0.102E-02
   0.102E-01 -.566E-02 -.177E-02   -.123E-01 0.661E-02 0.255E-02   0.201E-03 0.179E-02 -.108E-02   -.218E-03 0.131E-03 -.388E-04
 -----------------------------------------------------------------------------------------------
   -.115E-02 -.108E-02 0.806E-03   -.356E-14 0.712E-14 -.711E-14   0.279E-02 0.586E-03 -.711E-03   -.271E-03 0.156E-03 -.481E-04
 
 
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      5.38154      4.06861      3.60573        -0.020999      0.037724      0.054091
      3.94540      4.80321      4.38469         0.071295     -0.008772      0.003936
      5.52092      5.65450      4.42709        -0.032029     -0.063643      0.001782
      5.28185      4.16415      5.39183        -0.014802      0.031482     -0.059424
      5.03251      4.67255      4.45232        -0.002097      0.002870     -0.000339
 -----------------------------------------------------------------------------------
    total drift:                                0.001368     -0.000338      0.000047


--------------------------------------------------------------------------------------------------------



  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -24.12709802 eV

  energy  without entropy=      -24.12709802  energy(sigma->0) =      -24.12709802
  enthalpy is  TOTEN    =       -23.50294738 eV   P V=        0.62415064



--------------------------------------------------------------------------------------------------------


    POTLOK:  cpu time      0.0551: real time      0.0561


--------------------------------------------------------------------------------------------------------


    OFIELD:  cpu time      0.0000: real time      0.0000
           RANDOM_SEED =         283862281                0                0
           RANDOM_SEED =         283862281               66                0
   IONSTEP:  cpu time      0.0451: real time      0.0604

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -24.127098  see above
  kinetic energy EKIN   =         0.000027
  kin. lattice  EKIN_LAT=         0.003516  (temperature    3.43 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -24.123556 eV

  maximum distance moved by ions :      0.78E-04

    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:      0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  in kB       0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  external pressure =       -1.00 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =      0.00 kB
  Total+kin.     0.000       0.000       0.000       0.000       0.000      -0.000
  volume of cell :     1000.12
      direct lattice vectors                 reciprocal lattice vectors
    10.001524767 -0.002342671  0.000853163     0.099984755  0.000000000  0.000000000
     0.000000000  9.999300091  0.000429911     0.000023425  0.100007000  0.000000000
     0.000000000  0.000000000 10.000423975    -0.000008531 -0.000004299  0.099995760

  length of vectors
    10.001525078  9.999300101 10.000423975     0.099984755  0.100007002  0.099995761


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38218      4.06754      3.60708         0.000000      0.000000      0.000000
      3.94666      4.80192      4.38551         0.000000      0.000000      0.000000
      5.52131      5.65218      4.42807         0.000000      0.000000      0.000000
      5.28249      4.16291      5.39221         0.000000      0.000000      0.000000
      5.03327      4.67105      4.45322         0.000000      0.000000      0.000000
 -----------------------------------------------------------------------------------
    total drift:                                0.000000      0.000000      0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =         0.00000000 eV

  ML energy  without entropy=        0.00000000  ML energy(sigma->0) =        0.00000000

  enthalpy is ML TOTEN    =         0.62422858 eV   P V=        0.62422858

    WAVPRE:  cpu time      0.0086: real time      0.0088
    FEWALD:  cpu time      0.0007: real time      0.0007
    GENKIN:  cpu time      0.0005: real time      0.0006
    ORTHCH:  cpu time      0.0010: real time      0.0010
     LOOP+:  cpu time      1.4332: real time      1.5063


--------------------------------------- Iteration      2(   1)  ---------------------------------------


    POTLOK:  cpu time      0.0526: real time      0.0531
    SETDIJ:  cpu time      0.0009: real time      0.0009
     EDDAV:  cpu time      0.0266: real time      0.0266
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0049: real time      0.0049
    MIXING:  cpu time      0.0023: real time      0.0022
    --------------------------------------------
      LOOP:  cpu time      0.0875: real time      0.0879

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.1706435E-03  (-0.1165652E-03)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1911029 magnetization 

 Broyden mixing:
  rms(total) = 0.23045E-02    rms(broyden)= 0.23039E-02
  rms(prec ) = 0.31603E-02
  weight for this iteration     100.00

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23916342
  Ewald energy   TEWEN  =       128.83911758
  -Hartree energ DENC   =      -284.70882986
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.40030945
  PAW double counting   =        89.17992019      -90.62940751
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.02873246
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12725965 eV

  energy without entropy =      -24.12725965  energy(sigma->0) =      -24.12725965


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      2(   2)  ---------------------------------------


    POTLOK:  cpu time      0.0545: real time      0.0546
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0123: real time      0.0123
    ORTHCH:  cpu time      0.0003: real time      0.0003
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0051: real time      0.0051
    MIXING:  cpu time      0.0017: real time      0.0017
    --------------------------------------------
      LOOP:  cpu time      0.0781: real time      0.0782

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.1740768E-04  (-0.4814525E-05)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1910459 magnetization 

 Broyden mixing:
  rms(total) = 0.83309E-03    rms(broyden)= 0.83304E-03
  rms(prec ) = 0.10968E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.5109
  1.5109

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23916342
  Ewald energy   TEWEN  =       128.83911758
  -Hartree energ DENC   =      -284.75195558
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.40263145
  PAW double counting   =        89.26612019      -90.71674417
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -89.98677466
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12724224 eV

  energy without entropy =      -24.12724224  energy(sigma->0) =      -24.12724224


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      2(   3)  ---------------------------------------


    POTLOK:  cpu time      0.0549: real time      0.0551
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0034: real time      0.0034
  RMM-DIIS:  cpu time      0.0121: real time      0.0121
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0049: real time      0.0049
    MIXING:  cpu time      0.0019: real time      0.0019
    --------------------------------------------
      LOOP:  cpu time      0.0784: real time      0.0785

 eigenvalue-minimisations  :    15
 total energy-change (2. order) :-0.6968962E-06  (-0.3198050E-05)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1910017 magnetization 

 Broyden mixing:
  rms(total) = 0.36355E-03    rms(broyden)= 0.36333E-03
  rms(prec ) = 0.47551E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.2250
  1.2250  1.2250

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23916342
  Ewald energy   TEWEN  =       128.83911758
  -Hartree energ DENC   =      -284.77226636
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.40378940
  PAW double counting   =        89.31370581      -90.76512092
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -89.96683140
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12724294 eV

  energy without entropy =      -24.12724294  energy(sigma->0) =      -24.12724294


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      2(   4)  ---------------------------------------


    POTLOK:  cpu time      0.0553: real time      0.0554
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0034
  RMM-DIIS:  cpu time      0.0104: real time      0.0104
    ORTHCH:  cpu time      0.0003: real time      0.0003
       DOS:  cpu time      0.0001: real time      0.0001
    --------------------------------------------
      LOOP:  cpu time      0.0702: real time      0.0704

 eigenvalue-minimisations  :    12
 total energy-change (2. order) : 0.1863032E-06  (-0.1750405E-06)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1910017 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23916342
  Ewald energy   TEWEN  =       128.83911758
  -Hartree energ DENC   =      -284.76495037
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.40346872
  PAW double counting   =        89.30544624      -90.75683607
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -89.97385181
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12724275 eV

  energy without entropy =      -24.12724275  energy(sigma->0) =      -24.12724275


--------------------------------------------------------------------------------------------------------




 average (electrostatic) potential at core
  the test charge radii are     0.5201  0.6991
  (the norm of the test charge is              1.0000)
       1 -41.5622       2 -41.5589       3 -41.5625       4 -41.5612       5 -58.9962
 
 
 
 E-fermi :  -9.1113     XC(G=0):  -0.6146     alpha+bet : -0.2111

 Fermi energy:        -9.1113244552

 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -16.9271      2.00000
      2      -9.3541      2.00000
      3      -9.3537      2.00000
      4      -9.3520      2.00000
      5      -0.5426      0.00000
      6       1.0435      0.00000
      7       1.0580      0.00000
      8       1.0626      0.00000


--------------------------------------------------------------------------------------------------------


 soft charge-density along one line, spin component           1
         0         1         2         3         4         5         6         7         8         9
 total charge-density along one line
 
 pseudopotential strength for first ion, spin component:           1
 -2.357  -0.011  -0.001  -0.002   0.001
 -0.011   0.050  -0.001  -0.002   0.001
 -0.001  -0.001  -0.339   0.002  -0.001
 -0.002  -0.002   0.002  -0.338  -0.001
  0.001   0.001  -0.001  -0.001  -0.340
 total augmentation occupancy for first ion, spin component:           1
  1.739  -0.486   0.166   0.233  -0.096
 -0.486   0.185  -0.053  -0.074   0.030
  0.166  -0.053   0.023   0.020  -0.008
  0.233  -0.074   0.020   0.037  -0.011
 -0.096   0.030  -0.008  -0.011   0.014


------------------------ aborting loop because EDIFF is reached ----------------------------------------


    CHARGE:  cpu time      0.0051: real time      0.0051
    FORLOC:  cpu time      0.0037: real time      0.0037
    FORNL :  cpu time      0.0013: real time      0.0013
    STRESS:  cpu time      0.0214: real time      0.0214
    FORHAR:  cpu time      0.0130: real time      0.0130
    MIXING:  cpu time      0.0021: real time      0.0021
    OFIELD:  cpu time      0.0000: real time      0.0000

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Alpha Z     0.23916     0.23916     0.23916
  Ewald      42.91451    42.96761    42.95695    -0.01990     0.00257     0.00040
  Hartree    94.90593    94.93242    94.92459    -0.01349     0.00076     0.00110
  E(xc)     -27.06317   -27.06326   -27.06315    -0.00016     0.00001     0.00000
  Local    -196.97779  -197.04696  -197.03093     0.03332    -0.00258    -0.00190
  n-local   -15.56929   -15.57035   -15.56939    -0.00050    -0.00010     0.00002
  augment     0.19029     0.19050     0.19035     0.00001    -0.00006    -0.00004
  Kinetic   101.21889   101.22059   101.21464     0.00583     0.00000    -0.00069
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  -------------------------------------------------------------------------------------
  Total      -0.14146    -0.13029    -0.13777     0.00511     0.00061    -0.00111
  in kB      -0.22662    -0.20872    -0.22071     0.00818     0.00097    -0.00178
  external pressure =       -1.22 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =     -0.22 kB
  Total+kin.    -0.227      -0.209      -0.221       0.008       0.001      -0.002

 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
  energy-cutoff  :      520.00
  volume of cell :     1000.12
      direct lattice vectors                 reciprocal lattice vectors
    10.001524767 -0.002342671  0.000853163     0.099984755  0.000000000  0.000000000
     0.000000000  9.999300091  0.000429911     0.000023425  0.100007000  0.000000000
     0.000000000  0.000000000 10.000423975    -0.000008531 -0.000004299  0.099995760

  length of vectors
    10.001525078  9.999300101 10.000423975     0.099984755  0.100007002  0.099995761


 FORCES acting on ions
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction
 -----------------------------------------------------------------------------------------------
   -.168E+02 0.291E+02 0.408E+02   0.186E+02 -.322E+02 -.451E+02   -.177E+01 0.306E+01 0.429E+01   0.319E-05 0.284E-03 0.313E-03
   0.524E+02 -.631E+01 0.327E+01   -.578E+02 0.697E+01 -.361E+01   0.551E+01 -.664E+00 0.343E+00   0.458E-03 0.131E-03 -.780E-04
   -.236E+02 -.473E+02 0.121E+01   0.260E+02 0.523E+02 -.134E+01   -.248E+01 -.498E+01 0.128E+00   -.894E-04 -.361E-03 0.260E-04
   -.120E+02 0.245E+02 -.453E+02   0.133E+02 -.271E+02 0.500E+02   -.126E+01 0.258E+01 -.476E+01   0.162E-03 0.316E-03 -.402E-03
   0.352E-01 0.260E-01 -.662E-02   -.424E-01 -.295E-01 0.829E-02   -.203E-02 0.794E-03 -.103E-02   -.175E-03 -.233E-03 -.163E-03
 -----------------------------------------------------------------------------------------------
   0.530E-02 0.576E-02 0.233E-03   -.178E-14 0.354E-14 -.709E-14   -.456E-02 -.637E-02 0.235E-03   0.358E-03 0.137E-03 -.304E-03
 
 
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      5.38218      4.06754      3.60708        -0.014768      0.025269      0.037224
      3.94666      4.80192      4.38551         0.054453     -0.008200      0.003241
      5.52131      5.65218      4.42807        -0.018044     -0.036240      0.001120
      5.28249      4.16291      5.39221        -0.011155      0.021607     -0.041895
      5.03327      4.67105      4.45322        -0.009395     -0.002909      0.000474
 -----------------------------------------------------------------------------------
    total drift:                                0.001091     -0.000473      0.000164


--------------------------------------------------------------------------------------------------------



  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -24.12724275 eV

  energy  without entropy=      -24.12724275  energy(sigma->0) =      -24.12724275
  enthalpy is  TOTEN    =       -23.50301418 eV   P V=        0.62422858

 d Force = 0.1633910E-03[ 0.132E-03, 0.195E-03]  d Energy = 0.6679337E-04 0.966E-04
 d Force =-0.1617756E+00[-0.162E+00,-0.162E+00]  d Ewald  =-0.1564140E+00-0.536E-02


--------------------------------------------------------------------------------------------------------


    POTLOK:  cpu time      0.0588: real time      0.0588


--------------------------------------------------------------------------------------------------------


    OFIELD:  cpu time      0.0000: real time      0.0000
           RANDOM_SEED =         283862281              114                0
   IONSTEP:  cpu time      0.1698: real time      0.1843

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -24.127243  see above
  kinetic energy EKIN   =         0.000200
  kin. lattice  EKIN_LAT=         0.003501  (temperature    3.58 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -24.123542 eV

  maximum distance moved by ions :      0.14E-03

    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:     -0.02269    -0.00091    -0.02248     0.01050     0.00429    -0.00320
  in kB      -0.03635    -0.00145    -0.03601     0.01682     0.00687    -0.00513
  external pressure =       -1.02 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =     -0.02 kB
  Total+kin.    -0.036      -0.001      -0.036       0.017       0.007      -0.005
  volume of cell :     1000.23
      direct lattice vectors                 reciprocal lattice vectors
    10.003006767 -0.004715163  0.001714120     0.099969941 -0.000000000 -0.000000000
     0.000000000  9.998558676  0.000846627     0.000047144  0.100014415  0.000000000
     0.000000000 -0.000000000 10.000757138    -0.000017139 -0.000008467  0.099992429

  length of vectors
    10.003008025  9.998558712 10.000757138     0.099969941  0.100014426  0.099992431


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38273      4.06687      3.60888        -0.005703      0.001764      0.008718
      3.94848      4.80068      4.38628         0.017890     -0.005978      0.002622
      5.52145      5.64931      4.42908         0.001066      0.003984      0.001151
      5.28312      4.16171      5.39214        -0.007233      0.005835     -0.015892
      5.03397      4.66952      4.45406        -0.006021     -0.005604      0.003401
 -----------------------------------------------------------------------------------
    total drift:                                0.000000      0.000000     -0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =       -24.12738815 eV

  ML energy  without entropy=      -24.12738815  ML energy(sigma->0) =      -24.12738815

  enthalpy is ML TOTEN    =       -23.50309257 eV   P V=        0.62429558

    WAVPRE:  cpu time      0.0096: real time      0.0098
    FEWALD:  cpu time      0.0002: real time      0.0002
    GENKIN:  cpu time      0.0005: real time      0.0005
    ORTHCH:  cpu time      0.0010: real time      0.0010
     LOOP+:  cpu time      0.6044: real time      0.6211


--------------------------------------- Iteration      3(   1)  ---------------------------------------


    POTLOK:  cpu time      0.0564: real time      0.0565
    SETDIJ:  cpu time      0.0009: real time      0.0009
     EDDAV:  cpu time      0.0267: real time      0.0267
       DOS:  cpu time      0.0005: real time      0.0005
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0019: real time      0.0019
    --------------------------------------------
      LOOP:  cpu time      0.0915: real time      0.0916

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.1746511E-03  (-0.3450564E-03)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1915077 magnetization 

 Broyden mixing:
  rms(total) = 0.38422E-02    rms(broyden)= 0.38412E-02
  rms(prec ) = 0.52817E-02
  weight for this iteration     100.00

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23913775
  Ewald energy   TEWEN  =       129.12002510
  -Hartree energ DENC   =      -284.91442329
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.41084863
  PAW double counting   =        89.30894760      -90.76019412
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.11295880
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12741759 eV

  energy without entropy =      -24.12741759  energy(sigma->0) =      -24.12741759


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      3(   2)  ---------------------------------------


    POTLOK:  cpu time      0.0573: real time      0.0575
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0039: real time      0.0039
  RMM-DIIS:  cpu time      0.0136: real time      0.0136
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0017: real time      0.0017
    --------------------------------------------
      LOOP:  cpu time      0.0828: real time      0.0829

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.4828479E-04  (-0.1444238E-04)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1914134 magnetization 

 Broyden mixing:
  rms(total) = 0.13377E-02    rms(broyden)= 0.13376E-02
  rms(prec ) = 0.17650E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.4868
  1.4868

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23913775
  Ewald energy   TEWEN  =       129.12002510
  -Hartree energ DENC   =      -284.98784245
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.41478583
  PAW double counting   =        89.46501121      -90.91652449
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.04316178
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12736931 eV

  energy without entropy =      -24.12736931  energy(sigma->0) =      -24.12736931


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      3(   3)  ---------------------------------------


    POTLOK:  cpu time      0.0571: real time      0.0572
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0035: real time      0.0035
  RMM-DIIS:  cpu time      0.0126: real time      0.0126
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0019: real time      0.0019
    --------------------------------------------
      LOOP:  cpu time      0.0812: real time      0.0814

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.2422001E-05  (-0.9368841E-05)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1913363 magnetization 

 Broyden mixing:
  rms(total) = 0.60998E-03    rms(broyden)= 0.60958E-03
  rms(prec ) = 0.77351E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.0041
  1.3123  0.6960

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23913775
  Ewald energy   TEWEN  =       129.12002510
  -Hartree energ DENC   =      -285.02230627
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.41677469
  PAW double counting   =        89.55329552      -91.00469810
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.01079995
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12737173 eV

  energy without entropy =      -24.12737173  energy(sigma->0) =      -24.12737173


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      3(   4)  ---------------------------------------


    POTLOK:  cpu time      0.0589: real time      0.0594
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0098: real time      0.0100
    ORTHCH:  cpu time      0.0003: real time      0.0761
       DOS:  cpu time      0.0001: real time      0.0001
    --------------------------------------------
      LOOP:  cpu time      0.0733: real time      0.1497

 eigenvalue-minimisations  :    12
 total energy-change (2. order) : 0.2510364E-06  (-0.4260472E-06)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1913363 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23913775
  Ewald energy   TEWEN  =       129.12002510
  -Hartree energ DENC   =      -285.01556747
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.41650856
  PAW double counting   =        89.54721031      -90.99833546
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.01754980
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12737148 eV

  energy without entropy =      -24.12737148  energy(sigma->0) =      -24.12737148


--------------------------------------------------------------------------------------------------------




 average (electrostatic) potential at core
  the test charge radii are     0.5201  0.6991
  (the norm of the test charge is              1.0000)
       1 -41.5730       2 -41.5669       3 -41.5745       4 -41.5692       5 -58.9867
 
 
 
 E-fermi :  -9.1185     XC(G=0):  -0.6228     alpha+bet : -0.2111

 Fermi energy:        -9.1185412810

 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -16.9357      2.00000
      2      -9.3595      2.00000
      3      -9.3573      2.00000
      4      -9.3563      2.00000
      5      -0.5440      0.00000
      6       0.9905      0.00000
      7       1.0451      0.00000
      8       1.0571      0.00000


--------------------------------------------------------------------------------------------------------


 soft charge-density along one line, spin component           1
         0         1         2         3         4         5         6         7         8         9
 total charge-density along one line
 
 pseudopotential strength for first ion, spin component:           1
 -2.358  -0.011  -0.001  -0.002   0.001
 -0.011   0.050  -0.001  -0.002   0.001
 -0.001  -0.001  -0.339   0.002  -0.001
 -0.002  -0.002   0.002  -0.338  -0.001
  0.001   0.001  -0.001  -0.001  -0.340
 total augmentation occupancy for first ion, spin component:           1
  1.742  -0.486   0.167   0.234  -0.097
 -0.486   0.185  -0.053  -0.074   0.031
  0.167  -0.053   0.023   0.020  -0.008
  0.234  -0.074   0.020   0.037  -0.011
 -0.097   0.031  -0.008  -0.011   0.014


------------------------ aborting loop because EDIFF is reached ----------------------------------------


    CHARGE:  cpu time      0.0058: real time      0.0058
    FORLOC:  cpu time      0.0043: real time      0.0043
    FORNL :  cpu time      0.0014: real time      0.0014
    STRESS:  cpu time      0.0214: real time      0.0214
    FORHAR:  cpu time      0.0129: real time      0.0129
    MIXING:  cpu time      0.0021: real time      0.0021
    OFIELD:  cpu time      0.0000: real time      0.0000

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Alpha Z     0.23914     0.23914     0.23914
  Ewald      42.99968    43.08225    43.03805    -0.04447     0.00385     0.00308
  Hartree    94.97529    95.02367    95.00474    -0.02586    -0.00149     0.00457
  E(xc)     -27.07690   -27.07707   -27.07697    -0.00034    -0.00004     0.00004
  Local    -197.11680  -197.23261  -197.18136     0.06895     0.00068    -0.00913
  n-local   -15.59867   -15.60095   -15.59872    -0.00114    -0.00063     0.00034
  augment     0.18991     0.19003     0.19001     0.00001    -0.00009    -0.00002
  Kinetic   101.29657   101.30477   101.29624     0.01522     0.00306    -0.00305
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  -------------------------------------------------------------------------------------
  Total      -0.09178    -0.07077    -0.08889     0.01238     0.00533    -0.00418
  in kB      -0.14701    -0.11335    -0.14239     0.01982     0.00854    -0.00669
  external pressure =       -1.13 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =     -0.13 kB
  Total+kin.    -0.147      -0.113      -0.142       0.020       0.009      -0.007

 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
  energy-cutoff  :      520.00
  volume of cell :     1000.23
      direct lattice vectors                 reciprocal lattice vectors
    10.003006767 -0.004715163  0.001714120     0.099969941 -0.000000000 -0.000000000
     0.000000000  9.998558676  0.000846627     0.000047144  0.100014415  0.000000000
     0.000000000 -0.000000000 10.000757138    -0.000017139 -0.000008467  0.099992429

  length of vectors
    10.003008025  9.998558712 10.000757138     0.099969941  0.100014426  0.099992431


 FORCES acting on ions
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction
 -----------------------------------------------------------------------------------------------
   -.169E+02 0.292E+02 0.409E+02   0.186E+02 -.322E+02 -.452E+02   -.178E+01 0.307E+01 0.431E+01   -.429E-03 0.101E-02 0.133E-02
   0.524E+02 -.634E+01 0.327E+01   -.579E+02 0.700E+01 -.362E+01   0.553E+01 -.669E+00 0.345E+00   0.178E-02 -.434E-04 -.231E-04
   -.236E+02 -.474E+02 0.121E+01   0.261E+02 0.524E+02 -.134E+01   -.249E+01 -.500E+01 0.128E+00   -.677E-03 -.156E-02 0.845E-06
   -.120E+02 0.245E+02 -.453E+02   0.133E+02 -.271E+02 0.501E+02   -.127E+01 0.258E+01 -.478E+01   -.154E-03 0.101E-02 -.164E-02
   0.727E-01 0.652E-01 -.306E-01   -.845E-01 -.745E-01 0.366E-01   -.130E-02 -.450E-03 0.168E-02   -.143E-03 -.293E-03 -.287E-04
 -----------------------------------------------------------------------------------------------
   0.115E-01 0.143E-01 -.731E-02   0.534E-14 -.357E-14 0.713E-14   -.108E-01 -.148E-01 0.783E-02   0.378E-03 0.121E-03 -.361E-03
 
 
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      5.38273      4.06687      3.60888        -0.003244      0.002252      0.005300
      3.94848      4.80068      4.38628         0.019887     -0.006329      0.001733
      5.52145      5.64931      4.42908         0.003680      0.006508      0.000486
      5.28312      4.16171      5.39214        -0.005964      0.007259     -0.014990
      5.03397      4.66952      4.45406        -0.013299     -0.010050      0.007630
 -----------------------------------------------------------------------------------
    total drift:                                0.001061     -0.000360      0.000160


--------------------------------------------------------------------------------------------------------



  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -24.12737148 eV

  energy  without entropy=      -24.12737148  energy(sigma->0) =      -24.12737148
  enthalpy is  TOTEN    =       -23.50307590 eV   P V=        0.62429558

 d Force = 0.1365306E-03[ 0.403E-04, 0.233E-03]  d Energy = 0.6172069E-04 0.748E-04
 d Force =-0.2855262E+00[-0.286E+00,-0.285E+00]  d Ewald  =-0.2809075E+00-0.462E-02


--------------------------------------------------------------------------------------------------------


    POTLOK:  cpu time      0.0562: real time      0.0572


--------------------------------------------------------------------------------------------------------


    OFIELD:  cpu time      0.0000: real time      0.0000
           RANDOM_SEED =         283862281              162                0
   IONSTEP:  cpu time      0.1630: real time      0.1634

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -24.127371  see above
  kinetic energy EKIN   =         0.000354
  kin. lattice  EKIN_LAT=         0.003480  (temperature    3.71 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -24.123537 eV

  maximum distance moved by ions :      0.15E-03

    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:      0.04463     0.07078     0.03938     0.01440     0.00696    -0.00472
  in kB       0.07148     0.11337     0.06307     0.02307     0.01115    -0.00756
  external pressure =       -0.92 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =      0.08 kB
  Total+kin.     0.072       0.114       0.063       0.023       0.011      -0.008
  volume of cell :     1000.33
      direct lattice vectors                 reciprocal lattice vectors
    10.004437802 -0.007098076  0.002622609     0.099955642 -0.000000000  0.000000000
     0.000000000  9.997809556  0.001204539     0.000070965  0.100021909  0.000000000
    -0.000000000 -0.000000000 10.001063584    -0.000026220 -0.000012047  0.099989365

  length of vectors
    10.004440664  9.997809629 10.001063584     0.099955642  0.100021934  0.099989369


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38317      4.06633      3.61060         0.008312     -0.025622     -0.025738
      3.95055      4.79954      4.38702        -0.027423     -0.002673      0.000867
      5.52148      5.64655      4.42994         0.024088      0.049406      0.000984
      5.28370      4.16073      5.39195         0.001108     -0.014532      0.018921
      5.03467      4.66795      4.45488        -0.006085     -0.006580      0.004965
 -----------------------------------------------------------------------------------
    total drift:                                0.000000     -0.000000      0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =       -24.12729692 eV

  ML energy  without entropy=      -24.12729692  ML energy(sigma->0) =      -24.12729692

  enthalpy is ML TOTEN    =       -23.50293968 eV   P V=        0.62435724

     LOOP+:  cpu time      0.5988: real time      0.7244
           RANDOM_SEED =         283862281              210                0
   IONSTEP:  cpu time      0.0006: real time      0.0006

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -23.502940  see above
  kinetic energy EKIN   =         0.000324
  kin. lattice  EKIN_LAT=         0.003419  (temperature    3.62 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -23.499197 eV

    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:      0.09499     0.12555     0.09187     0.01588     0.00662    -0.00475
  in kB       0.15213     0.20107     0.14713     0.02544     0.01061    -0.00760
  external pressure =       -0.83 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =      0.17 kB
  Total+kin.     0.152       0.201       0.147       0.025       0.011      -0.008
  volume of cell :     1000.41
      direct lattice vectors                 reciprocal lattice vectors
    10.005811992 -0.009476043  0.003546260     0.099941914 -0.000000000 -0.000000000
     0.000000000  9.996997221  0.001586820     0.000094734  0.100030037  0.000000000
     0.000000000 -0.000000000 10.001339909    -0.000035452 -0.000015871  0.099986603

  length of vectors
    10.005817108  9.996997347 10.001339909     0.099941914  0.100030082  0.099986610


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38366      4.06531      3.61238         0.018774     -0.046009     -0.053197
      3.95233      4.79835      4.38784        -0.061552      0.000611     -0.001105
      5.52169      5.64406      4.43077         0.040478      0.083410      0.000202
      5.28436      4.15971      5.39182         0.008457     -0.032109      0.050141
      5.03534      4.66638      4.45578        -0.006157     -0.005903      0.003959
 -----------------------------------------------------------------------------------
    total drift:                               -0.000000      0.000000     -0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =       -24.12706616 eV

  ML energy  without entropy=      -24.12706616  ML energy(sigma->0) =      -24.12706616

  enthalpy is ML TOTEN    =       -23.50265664 eV   P V=        0.62440952

     LOOP+:  cpu time      0.0015: real time      0.0017
           RANDOM_SEED =         283862281              258                0
   IONSTEP:  cpu time      0.0006: real time      0.0006

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -23.502657  see above
  kinetic energy EKIN   =         0.000141
  kin. lattice  EKIN_LAT=         0.003394  (temperature    3.42 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -23.499121 eV

    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:      0.11219     0.14285     0.11266     0.01421     0.00782    -0.00544
  in kB       0.17966     0.22876     0.18042     0.02276     0.01253    -0.00872
  external pressure =       -0.80 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =      0.20 kB
  Total+kin.     0.180       0.229       0.181       0.023       0.013      -0.009
  volume of cell :     1000.48
      direct lattice vectors                 reciprocal lattice vectors
    10.007124992 -0.011886170  0.004436715     0.099928801 -0.000000000  0.000000000
     0.000000000  9.996125727  0.001863572     0.000118823  0.100038758  0.000000000
    -0.000000000 -0.000000000 10.001547736    -0.000044351 -0.000018640  0.099984525

  length of vectors
    10.007133034  9.996125901 10.001547736     0.099928801  0.100038828  0.099984537


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38426      4.06385      3.61363         0.023426     -0.054326     -0.064816
      3.95338      4.79697      4.38836        -0.073490      0.002314     -0.001376
      5.52206      5.64234      4.43155         0.044640      0.092764      0.000430
      5.28515      4.15831      5.39216         0.011124     -0.037935      0.061930
      5.03597      4.66473      4.45658        -0.005701     -0.002818      0.003831
 -----------------------------------------------------------------------------------
    total drift:                               -0.000000     -0.000000     -0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =       -24.12695084 eV

  ML energy  without entropy=      -24.12695084  ML energy(sigma->0) =      -24.12695084

  enthalpy is ML TOTEN    =       -23.50250084 eV   P V=        0.62444999

    WAVPRE:  cpu time      0.0096: real time      0.0098
    FEWALD:  cpu time      0.0002: real time      0.0002
    GENKIN:  cpu time      0.0005: real time      0.0005
    ORTHCH:  cpu time      0.0007: real time      0.0007
     LOOP+:  cpu time      0.0133: real time      0.0137


--------------------------------------- Iteration      6(   1)  ---------------------------------------


    POTLOK:  cpu time      0.0562: real time      0.0564
    SETDIJ:  cpu time      0.0009: real time      0.0009
     EDDAV:  cpu time      0.0203: real time      0.0203
       DOS:  cpu time      0.0007: real time      0.0007
    --------------------------------------------
      LOOP:  cpu time      0.0782: real time      0.0783

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.1199621E-03  (-0.1932763E-02)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1913366 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23907862
  Ewald energy   TEWEN  =       129.77161038
  -Hartree energ DENC   =      -285.35065823
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.43314619
  PAW double counting   =        89.54024632      -90.99063737
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.35123721
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12725177 eV

  energy without entropy =      -24.12725177  energy(sigma->0) =      -24.12725177


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      6(   2)  ---------------------------------------


     EDDAV:  cpu time      0.0251: real time      0.0253
       DOS:  cpu time      0.0001: real time      0.0001
    --------------------------------------------
      LOOP:  cpu time      0.0252: real time      0.0254

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.6033886E-05  (-0.6033886E-05)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1913366 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23907862
  Ewald energy   TEWEN  =       129.77161038
  -Hartree energ DENC   =      -285.35065823
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.43314619
  PAW double counting   =        89.54024632      -90.99063737
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.35124325
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12725780 eV

  energy without entropy =      -24.12725780  energy(sigma->0) =      -24.12725780


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      6(   3)  ---------------------------------------


     EDDAV:  cpu time      0.0246: real time      0.0247
       DOS:  cpu time      0.0001: real time      0.0001
    --------------------------------------------
      LOOP:  cpu time      0.0247: real time      0.0248

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.1278295E-08  (-0.1278273E-08)
 number of electron       8.0000007 magnetization 
 augmentation part        0.1913366 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23907862
  Ewald energy   TEWEN  =       129.77161038
  -Hartree energ DENC   =      -285.35065823
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.43314619
  PAW double counting   =        89.54024632      -90.99063737
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.35124325
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12725780 eV

  energy without entropy =      -24.12725780  energy(sigma->0) =      -24.12725780


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      6(   4)  ---------------------------------------


     EDDAV:  cpu time      0.0146: real time      0.0148
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0023: real time      0.0023
    --------------------------------------------
      LOOP:  cpu time      0.0220: real time      0.0222

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.2842171E-12  (-0.3268497E-12)
 number of electron       8.0000006 magnetization 
 augmentation part        0.1925236 magnetization 

 Broyden mixing:
  rms(total) = 0.93376E-02    rms(broyden)= 0.93351E-02
  rms(prec ) = 0.12840E-01
  weight for this iteration     100.00

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23907862
  Ewald energy   TEWEN  =       129.77161038
  -Hartree energ DENC   =      -285.35065823
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.43314619
  PAW double counting   =        89.54024632      -90.99063737
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.35124325
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12725780 eV

  energy without entropy =      -24.12725780  energy(sigma->0) =      -24.12725780


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      6(   5)  ---------------------------------------


    POTLOK:  cpu time      0.0575: real time      0.0577
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0053: real time      0.0053
  RMM-DIIS:  cpu time      0.0140: real time      0.0140
    ORTHCH:  cpu time      0.0008: real time      0.0008
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0051: real time      0.0051
    MIXING:  cpu time      0.0017: real time      0.0017
    --------------------------------------------
      LOOP:  cpu time      0.0854: real time      0.0855

 eigenvalue-minimisations  :    16
 total energy-change (2. order) : 0.2949021E-03  (-0.8794435E-04)
 number of electron       8.0000006 magnetization 
 augmentation part        0.1922784 magnetization 

 Broyden mixing:
  rms(total) = 0.30884E-02    rms(broyden)= 0.30882E-02
  rms(prec ) = 0.41050E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.4563
  1.4563

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23907862
  Ewald energy   TEWEN  =       129.77161038
  -Hartree energ DENC   =      -285.52790995
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.44255722
  PAW double counting   =        89.91225353      -91.36517099
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.18058123
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12696290 eV

  energy without entropy =      -24.12696290  energy(sigma->0) =      -24.12696290


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      6(   6)  ---------------------------------------


    POTLOK:  cpu time      0.0577: real time      0.0578
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0133: real time      0.0133
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    CHARGE:  cpu time      0.0050: real time      0.0050
    MIXING:  cpu time      0.0017: real time      0.0017
    --------------------------------------------
      LOOP:  cpu time      0.0822: real time      0.0823

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.1251783E-04  (-0.5293058E-04)
 number of electron       8.0000006 magnetization 
 augmentation part        0.1921152 magnetization 

 Broyden mixing:
  rms(total) = 0.14286E-02    rms(broyden)= 0.14275E-02
  rms(prec ) = 0.18066E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   0.8500
  1.4071  0.2929

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23907862
  Ewald energy   TEWEN  =       129.77161038
  -Hartree energ DENC   =      -285.60703382
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.44702809
  PAW double counting   =        90.11216482      -91.56507542
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.10594763
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12697542 eV

  energy without entropy =      -24.12697542  energy(sigma->0) =      -24.12697542


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      6(   7)  ---------------------------------------


    POTLOK:  cpu time      0.0577: real time      0.0578
    SETDIJ:  cpu time      0.0009: real time      0.0009
    EDDIAG:  cpu time      0.0033: real time      0.0033
  RMM-DIIS:  cpu time      0.0129: real time      0.0129
    ORTHCH:  cpu time      0.0002: real time      0.0002
       DOS:  cpu time      0.0001: real time      0.0001
    --------------------------------------------
      LOOP:  cpu time      0.0751: real time      0.0752

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.1325224E-05  (-0.2970894E-05)
 number of electron       8.0000006 magnetization 
 augmentation part        0.1921152 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.23907862
  Ewald energy   TEWEN  =       129.77161038
  -Hartree energ DENC   =      -285.59869592
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        25.44671036
  PAW double counting   =        90.10531907      -91.55788882
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -90.11430998
  atomic energy  EATOM  =       197.58119954
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -24.12697674 eV

  energy without entropy =      -24.12697674  energy(sigma->0) =      -24.12697674


--------------------------------------------------------------------------------------------------------




 average (electrostatic) potential at core
  the test charge radii are     0.5201  0.6991
  (the norm of the test charge is              1.0000)
       1 -41.5951       2 -41.5855       3 -41.5946       4 -41.5893       5 -58.9628
 
 
 
 E-fermi :  -9.0713     XC(G=0):  -0.6467     alpha+bet : -0.2110

 Fermi energy:        -9.0713031286

 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -16.9548      2.00000
      2      -9.3737      2.00000
      3      -9.3672      2.00000
      4      -9.3614      2.00000
      5      -0.5503      0.00000
      6       0.9390      0.00000
      7       1.0355      0.00000
      8       1.0484      0.00000


--------------------------------------------------------------------------------------------------------


 soft charge-density along one line, spin component           1
         0         1         2         3         4         5         6         7         8         9
 total charge-density along one line
 
 pseudopotential strength for first ion, spin component:           1
 -2.358  -0.011  -0.001  -0.002   0.001
 -0.011   0.050  -0.001  -0.002   0.001
 -0.001  -0.001  -0.339   0.002  -0.001
 -0.002  -0.002   0.002  -0.338  -0.001
  0.001   0.001  -0.001  -0.001  -0.340
 total augmentation occupancy for first ion, spin component:           1
  1.748  -0.487   0.168   0.236  -0.098
 -0.487   0.186  -0.053  -0.075   0.031
  0.168  -0.053   0.023   0.020  -0.008
  0.236  -0.075   0.020   0.037  -0.012
 -0.098   0.031  -0.008  -0.012   0.014


------------------------ aborting loop because EDIFF is reached ----------------------------------------


    CHARGE:  cpu time      0.0051: real time      0.0051
    FORLOC:  cpu time      0.0037: real time      0.0037
    FORNL :  cpu time      0.0014: real time      0.0014
    STRESS:  cpu time      0.0214: real time      0.0214
    FORHAR:  cpu time      0.0127: real time      0.0129
    MIXING:  cpu time      0.0021: real time      0.0020
    OFIELD:  cpu time      0.0000: real time      0.0000

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Alpha Z     0.23908     0.23908     0.23908
  Ewald      43.16907    43.35547    43.24702    -0.16439    -0.00729     0.04687
  Hartree    95.11615    95.25745    95.18234    -0.06239    -0.00387     0.01945
  E(xc)     -27.10865   -27.10891   -27.10894    -0.00111    -0.00010     0.00034
  Local    -197.39092  -197.70358  -197.53102     0.20353     0.01233    -0.06063
  n-local   -15.66511   -15.66738   -15.66556    -0.00114    -0.00096     0.00058
  augment     0.18897     0.18898     0.18899     0.00020    -0.00007    -0.00007
  Kinetic   101.48298   101.48419   101.48693     0.04225     0.00750    -0.01441
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  -------------------------------------------------------------------------------------
  Total       0.03157     0.04530     0.03883     0.01695     0.00754    -0.00786
  in kB       0.05056     0.07254     0.06219     0.02715     0.01208    -0.01259
  external pressure =       -0.94 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =      0.06 kB
  Total+kin.     0.051       0.073       0.062       0.027       0.012      -0.013

 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
  energy-cutoff  :      520.00
  volume of cell :     1000.48
      direct lattice vectors                 reciprocal lattice vectors
    10.007124992 -0.011886170  0.004436715     0.099928801 -0.000000000  0.000000000
     0.000000000  9.996125727  0.001863572     0.000118823  0.100038758  0.000000000
    -0.000000000 -0.000000000 10.001547736    -0.000044351 -0.000018640  0.099984525

  length of vectors
    10.007133034  9.996125901 10.001547736     0.099928801  0.100038828  0.099984537


 FORCES acting on ions
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction
 -----------------------------------------------------------------------------------------------
   -.169E+02 0.292E+02 0.410E+02   0.188E+02 -.324E+02 -.454E+02   -.179E+01 0.309E+01 0.434E+01   -.155E-02 0.323E-02 0.447E-02
   0.526E+02 -.642E+01 0.331E+01   -.582E+02 0.710E+01 -.366E+01   0.556E+01 -.682E+00 0.351E+00   0.561E-02 -.354E-03 -.281E-04
   -.237E+02 -.476E+02 0.121E+01   0.262E+02 0.527E+02 -.135E+01   -.250E+01 -.504E+01 0.129E+00   -.246E-02 -.578E-02 0.201E-03
   -.121E+02 0.246E+02 -.455E+02   0.134E+02 -.273E+02 0.503E+02   -.128E+01 0.260E+01 -.481E+01   -.662E-03 0.281E-02 -.518E-02
   0.116E+00 0.125E+00 -.666E-01   -.130E+00 -.135E+00 0.757E-01   0.715E-02 0.156E-01 -.305E-02   -.438E-03 -.106E-02 -.408E-03
 -----------------------------------------------------------------------------------------------
   0.100E-01 0.800E-02 -.809E-02   0.888E-14 0.000E+00 0.000E+00   -.980E-02 -.679E-02 0.913E-02   0.509E-03 -.115E-02 -.951E-03
 
 
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      5.38426      4.06385      3.61363         0.024419     -0.049633     -0.069040
      3.95338      4.79697      4.38836        -0.071515     -0.003306     -0.000696
      5.52206      5.64234      4.43155         0.046633      0.083639     -0.001086
      5.28515      4.15831      5.39216         0.009183     -0.034893      0.065236
      5.03597      4.66473      4.45658        -0.007990      0.004255      0.005674
 -----------------------------------------------------------------------------------
    total drift:                                0.000730      0.000064      0.000088


--------------------------------------------------------------------------------------------------------



  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -24.12697674 eV

  energy  without entropy=      -24.12697674  energy(sigma->0) =      -24.12697674
  enthalpy is  TOTEN    =       -23.50252675 eV   P V=        0.62444999

 d Force =-0.4086392E-03[-0.933E-03, 0.116E-03]  d Energy =-0.5491456E-03 0.141E-03
 d Force =-0.6623187E+00[-0.664E+00,-0.661E+00]  d Ewald  =-0.6515853E+00-0.107E-01


--------------------------------------------------------------------------------------------------------


    POTLOK:  cpu time      0.0563: real time      0.0563


--------------------------------------------------------------------------------------------------------


    OFIELD:  cpu time      0.0000: real time      0.0000
           RANDOM_SEED =         283862281              306                0
   IONSTEP:  cpu time      0.1635: real time      0.1641

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -24.126977  see above
  kinetic energy EKIN   =         0.000043
  kin. lattice  EKIN_LAT=         0.003370  (temperature    3.30 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -24.123564 eV

  maximum distance moved by ions :      0.96E-04

    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:      0.08931     0.11524     0.09970     0.00852     0.00418    -0.00346
  in kB       0.14301     0.18453     0.15966     0.01364     0.00669    -0.00554
  external pressure =       -0.84 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =      0.16 kB
  Total+kin.     0.143       0.185       0.160       0.014       0.007      -0.006
  volume of cell :     1000.53
      direct lattice vectors                 reciprocal lattice vectors
    10.008368905 -0.014284955  0.005356384     0.099916381 -0.000000000  0.000000000
     0.000000000  9.995201345  0.002089127     0.000142799  0.100048010  0.000000000
    -0.000000000 -0.000000000 10.001694428    -0.000053540 -0.000020898  0.099983059

  length of vectors
    10.008380533  9.995201563 10.001694428     0.099916381  0.100048111  0.099983075


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38502      4.06171      3.61434         0.019134     -0.044730     -0.055248
      3.95360      4.79533      4.38886        -0.058857      0.002524     -0.001598
      5.52270      5.64145      4.43227         0.033670      0.072506     -0.000146
      5.28601      4.15677      5.39290         0.010140     -0.033710      0.056667
      5.03649      4.66302      4.45732        -0.004087      0.003410      0.000325
 -----------------------------------------------------------------------------------
    total drift:                                0.000000     -0.000000     -0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =       -24.12708815 eV

  ML energy  without entropy=      -24.12708815  ML energy(sigma->0) =      -24.12708815

  enthalpy is ML TOTEN    =       -23.50260913 eV   P V=        0.62447902

     LOOP+:  cpu time      0.6619: real time      0.6662
           RANDOM_SEED =         283862281              354                0
   IONSTEP:  cpu time      0.0007: real time      0.0013

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -23.502609  see above
  kinetic energy EKIN   =         0.000179
  kin. lattice  EKIN_LAT=         0.003334  (temperature    3.40 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -23.499096 eV

    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:      0.03455     0.05217     0.05630     0.00044     0.00065    -0.00144
  in kB       0.05532     0.08353     0.09016     0.00071     0.00104    -0.00230
  external pressure =       -0.92 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =      0.08 kB
  Total+kin.     0.056       0.084       0.090       0.001       0.001      -0.002
  volume of cell :     1000.55
      direct lattice vectors                 reciprocal lattice vectors
    10.009538389 -0.016644227  0.006272819     0.099904707  0.000000000 -0.000000000
    -0.000000000  9.994202986  0.002391195     0.000166380  0.100058004  0.000000000
     0.000000000 -0.000000000 10.001788318    -0.000062697 -0.000023922  0.099982120

  length of vectors
    10.009554193  9.994203272 10.001788318     0.099904707  0.100058142  0.099982143


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38575      4.05929      3.61446         0.008579     -0.022966     -0.029710
      3.95304      4.79376      4.38957        -0.022923      0.000856     -0.000680
      5.52384      5.64127      4.43321         0.011309      0.029103     -0.000282
      5.28688      4.15481      5.39431         0.004818     -0.019167      0.033309
      5.03695      4.66126      4.45814        -0.001784      0.012174     -0.002637
 -----------------------------------------------------------------------------------
    total drift:                                0.000000      0.000000      0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =       -24.12731626 eV

  ML energy  without entropy=      -24.12731626  ML energy(sigma->0) =      -24.12731626

  enthalpy is ML TOTEN    =       -23.50282079 eV   P V=        0.62449547

     LOOP+:  cpu time      0.0016: real time      0.0024
           RANDOM_SEED =         283862281              402                0
   IONSTEP:  cpu time      0.0006: real time      0.0030

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -23.502821  see above
  kinetic energy EKIN   =         0.000433
  kin. lattice  EKIN_LAT=         0.003358  (temperature    3.67 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -23.499029 eV

    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:     -0.04264    -0.02738    -0.01300    -0.00314    -0.00257     0.00019
  in kB      -0.06828    -0.04385    -0.02081    -0.00504    -0.00412     0.00030
  external pressure =       -1.04 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =     -0.04 kB
  Total+kin.    -0.068      -0.043      -0.020      -0.005      -0.004       0.000
  volume of cell :     1000.57
      direct lattice vectors                 reciprocal lattice vectors
    10.010652335 -0.019025302  0.007189362     0.099893590  0.000000000 -0.000000000
    -0.000000000  9.993136707  0.002664900     0.000190181  0.100068680  0.000000000
     0.000000000 -0.000000000 10.001883736    -0.000071854 -0.000026662  0.099981166

  length of vectors
    10.010672996  9.993137062 10.001883736     0.099893590  0.100068861  0.099981196


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38675      4.05671      3.61409        -0.008107      0.006771      0.009126
      3.95196      4.79233      4.39031         0.028755     -0.004185      0.001455
      5.52498      5.64124      4.43418        -0.014858     -0.021888     -0.000071
      5.28791      4.15241      5.39608        -0.005395      0.002477     -0.005537
      5.03734      4.65950      4.45892        -0.000396      0.016825     -0.004974
 -----------------------------------------------------------------------------------
    total drift:                                0.000000     -0.000000     -0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =       -24.12735379 eV

  ML energy  without entropy=      -24.12735379  ML energy(sigma->0) =      -24.12735379

  enthalpy is ML TOTEN    =       -23.50284950 eV   P V=        0.62450429

     LOOP+:  cpu time      0.0017: real time      0.0042
           RANDOM_SEED =         283862281              450                0
   IONSTEP:  cpu time      0.0006: real time      0.0006

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -23.502849  see above
  kinetic energy EKIN   =         0.000509
  kin. lattice  EKIN_LAT=         0.003431  (temperature    3.81 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -23.498909 eV

    OFIELD:  cpu time      0.0000: real time      0.0000

  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Total:     -0.10733    -0.09866    -0.07502    -0.00639    -0.00814     0.00334
  in kB      -0.17186    -0.15799    -0.12013    -0.01023    -0.01303     0.00534
  external pressure =       -1.15 kB  Pullay stress =        1.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =     -0.15 kB
  Total+kin.    -0.171      -0.157      -0.120      -0.010      -0.013       0.005
  volume of cell :     1000.58
      direct lattice vectors                 reciprocal lattice vectors
    10.011776419 -0.021496755  0.008079854     0.099882374  0.000000000 -0.000000000
    -0.000000000  9.992079789  0.002921446     0.000214885  0.100079265  0.000000000
     0.000000000 -0.000000000 10.001923990    -0.000080751 -0.000029232  0.099980764

  length of vectors
    10.011802757  9.992080216 10.001923990     0.099882374  0.100079496  0.099980801


  POSITION                                       TOTAL-FORCE (eV/Angst) (ML)
 -----------------------------------------------------------------------------------
      5.38787      4.05405      3.61382        -0.023390      0.034980      0.045718
      3.95119      4.79079      4.39094         0.071382     -0.008284      0.002535
      5.52595      5.64098      4.43523        -0.037376     -0.067152     -0.000477
      5.28871      4.14990      5.39763        -0.013280      0.020513     -0.038585
      5.03772      4.65769      4.45972         0.002665      0.019943     -0.009192
 -----------------------------------------------------------------------------------
    total drift:                               -0.000000     -0.000000      0.000000


--------------------------------------------------------------------------------------------------------



  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy ML TOTEN  =       -24.12712228 eV

  ML energy  without entropy=      -24.12712228  ML energy(sigma->0) =      -24.12712228

  enthalpy is ML TOTEN    =       -23.50261141 eV   P V=        0.62451087

     LOOP+:  cpu time      0.0015: real time      0.0016
           RANDOM_SEED =         283862281              498                0
   IONSTEP:  cpu time      0.0006: real time      0.0006

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -23.502611  see above
  kinetic energy EKIN   =         0.000323
  kin. lattice  EKIN_LAT=         0.003518  (temperature    3.71 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -23.498771 eV


 mean value of Nose-termostat <S>:     1.000 mean value of <T> :     3.564
 mean temperature <T/S>/<1/S>  :     3.564

     LOOP+:  cpu time      0.0021: real time      0.0075

 total amount of memory used by VASP MPI-rank0    41513. kBytes
=======================================================================

   base      :      30000. kBytes
   nonl-proj :        699. kBytes
   fftplans  :       3956. kBytes
   grid      :       6748. kBytes
   one-center:          3. kBytes
   wavefun   :        107. kBytes
 
  
  
 General timing and accounting informations for this job:
 ========================================================
  
                  Total CPU time used (sec):        3.965
                            User time (sec):        3.483
                          System time (sec):        0.482
                         Elapsed time (sec):        4.234
  
                   Maximum memory used (kb):      470384.
                   Average memory used (kb):          N/A
  
                          Minor page faults:        43208
                          Major page faults:            0
                 Voluntary context switches:         4416
