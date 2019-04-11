 vasp.5.4.4.18Apr17-6-g9f103f2a35 (build Sep 18 2018 16:57:57) complex          
  
 executed on             LinuxIFC date 2019.04.10  04:05:36
 running on    1 total cores
 distrk:  each k-point on    1 cores,    1 groups
 distr:  one band on NCORES_PER_BAND=   1 cores,    1 groups


--------------------------------------------------------------------------------------------------------


 INCAR:
 POTCAR:    PAW_PBE O 08Apr2002                   
 POTCAR:    PAW_PBE H 15Jun2001                   
 POTCAR:    PAW_PBE O 08Apr2002                   
   VRHFIN =O: s2p4                                                              
   LEXCH  = PE                                                                  
   EATOM  =   432.3788 eV,   31.7789 Ry                                         
                                                                                
   TITEL  = PAW_PBE O 08Apr2002                                                 
   LULTRA =        F    use ultrasoft PP ?                                      
   IUNSCR =        1    unscreen: 0-lin 1-nonlin 2-no                           
   RPACOR =    1.200    partial core radius                                     
   POMASS =   16.000; ZVAL   =    6.000    mass and valenz                      
   RCORE  =    1.520    outmost cutoff radius                                   
   RWIGS  =    1.550; RWIGS  =    0.820    wigner-seitz radius (au A)           
   ENMAX  =  400.000; ENMIN  =  300.000 eV                                      
   ICORE  =        2    local potential                                         
   LCOR   =        T    correct aug charges                                     
   LPAW   =        T    paw PP                                                  
   EAUG   =  605.392                                                            
   DEXC   =    0.000                                                            
   RMAX   =    1.553    core radius for proj-oper                               
   RAUG   =    1.300    factor for augmentation sphere                          
   RDEP   =    1.550    radius for radial grids                                 
   RDEPT  =    1.329    core radius for aug-charge                              
                                                                                
   Atomic configuration                                                         
    4 entries                                                                   
     n  l   j            E        occ.                                          
     1  0  0.50      -514.6923   2.0000                                         
     2  0  0.50       -23.9615   2.0000                                         
     2  1  0.50        -9.0305   4.0000                                         
     3  2  1.50        -9.5241   0.0000                                         
   Description                                                                  
     l       E           TYP  RCUT    TYP  RCUT                                 
     0    -23.9615318     23  1.200                                             
     0     -9.5240782     23  1.200                                             
     1     -9.0304911     23  1.520                                             
     1      8.1634956     23  1.520                                             
     2     -9.5240782      7  1.500                                             
  local pseudopotential read in
  partial core-charges read in
  partial kinetic energy density read in
  kinetic energy density of atom read in
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
 
 POTCAR:    PAW_PBE H 15Jun2001                   
   VRHFIN =H: ultrasoft test                                                    
   LEXCH  = PE                                                                  
   EATOM  =    12.4884 eV,    0.9179 Ry                                         
                                                                                
   TITEL  = PAW_PBE H 15Jun2001                                                 
   LULTRA =        F    use ultrasoft PP ?                                      
   IUNSCR =        0    unscreen: 0-lin 1-nonlin 2-no                           
   RPACOR =    0.000    partial core radius                                     
   POMASS =    1.000; ZVAL   =    1.000    mass and valenz                      
   RCORE  =    1.100    outmost cutoff radius                                   
   RWIGS  =    0.700; RWIGS  =    0.370    wigner-seitz radius (au A)           
   ENMAX  =  250.000; ENMIN  =  200.000 eV                                      
   RCLOC  =    0.701    cutoff for local pot                                    
   LCOR   =        T    correct aug charges                                     
   LPAW   =        T    paw PP                                                  
   EAUG   =  400.000                                                            
   RMAX   =    1.123    core radius for proj-oper                               
   RAUG   =    1.200    factor for augmentation sphere                          
   RDEP   =    1.112    radius for radial grids                                 
   RDEPT  =    0.926    core radius for aug-charge                              
                                                                                
   Atomic configuration                                                         
    2 entries                                                                   
     n  l   j            E        occ.                                          
     1  0  0.50        -6.4927   1.0000                                         
     2  1  0.50        -3.4015   0.0000                                         
   Description                                                                  
     l       E           TYP  RCUT    TYP  RCUT                                 
     0     -6.4927494     23  1.100                                             
     0      6.8029130     23  1.100                                             
     1     -4.0817478     23  1.100                                             
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
 

 ----------------------------------------------------------------------------- 
|                                                                             |
|  ADVICE TO THIS USER RUNNING 'VASP/VAMP'   (HEAR YOUR MASTER'S VOICE ...):  |
|                                                                             |
|      You have a (more or less) 'small supercell' and for smaller cells      |
|      it is recommended  to use the reciprocal-space projection scheme!      |
|      The real space optimization is not  efficient for small cells and it   |
|      is also less accurate ...                                              |
|      Therefore set LREAL=.FALSE. in the  INCAR file                         |
|                                                                             |
 ----------------------------------------------------------------------------- 

 Optimization of the real space projectors (new method)

 maximal supplied QI-value         = 24.76
 optimisation between [QCUT,QGAM] = [ 14.36, 28.96] = [ 57.73,234.92] Ry 
 Optimized for a Real-space Cutoff    0.98 Angstroem

   l    n(q)    QCUT    max X(q) W(low)/X(q) W(high)/X(q)  e(spline) 
   0      8    14.358    20.381    0.72E-04    0.23E-03    0.89E-07
   0      8    14.358    15.268    0.76E-04    0.24E-03    0.10E-06
   1      8    14.358     5.964    0.17E-03    0.15E-03    0.14E-06
   1      8    14.358     5.382    0.15E-03    0.12E-03    0.12E-06
 Optimization of the real space projectors (new method)

 maximal supplied QI-value         = 34.20
 optimisation between [QCUT,QGAM] = [ 14.37, 28.73] = [ 57.79,231.16] Ry 
 Optimized for a Real-space Cutoff    0.95 Angstroem

   l    n(q)    QCUT    max X(q) W(low)/X(q) W(high)/X(q)  e(spline) 
   0      8    14.366    19.460    0.22E-03    0.25E-03    0.13E-06
   0      8    14.366    12.209    0.21E-03    0.23E-03    0.12E-06
   1      8    14.366     4.655    0.29E-04    0.14E-03    0.19E-06
  PAW_PBE O 08Apr2002                   :
 energy of atom  1       EATOM= -432.3788
 kinetic energy error for atom=    0.0041 (will be added to EATOM!!)
  PAW_PBE H 15Jun2001                   :
 energy of atom  2       EATOM=  -12.4884
 kinetic energy error for atom=    0.0002 (will be added to EATOM!!)
 
 
 POSCAR: POSCAR file written by OVITO            
  positions in direct lattice
  No initial velocities read in
 exchange correlation table for  LEXCH =        8
   RHO(1)=    0.500       N(1)  =     2000
   RHO(2)=  100.500       N(2)  =     4000
 


--------------------------------------------------------------------------------------------------------


 ion  position               nearest neighbor table
   1  0.428  0.424  0.520-   3 0.99   4 1.00
   2  0.230  0.628  0.113-   5 1.00   6 1.00
   3  0.458  0.352  0.458-   1 0.99
   4  0.389  0.384  0.603-   1 1.00
   5  0.137  0.626  0.150-   2 1.00
   6  0.231  0.589  0.021-   2 1.00
 

IMPORTANT INFORMATION: All symmetrisations will be switched off!
NOSYMM: (Re-)initialisation of all symmetry stuff for point group C_1.

 
 

Automatic generation of k-mesh.
 generate k-points for:    1    1    1
Space group operators:
 irot       det(A)        alpha          n_x          n_y          n_z        tau_x        tau_y        tau_z
    1     1.000000     0.000000     1.000000     0.000000     0.000000     0.000000     0.000000     0.000000
 
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
   k-points           NKPTS =      1   k-points in BZ     NKDIM =      1   number of bands    NBANDS=     12
   number of dos      NEDOS =    301   number of ions     NIONS =      6
   non local maximal  LDIM  =      4   non local SUM 2l+1 LMDIM =      8
   total plane-waves  NPLWV = 884736
   max r-space proj   IRMAX =   3497   max aug-charges    IRDMAX=  11677
   dimension x,y,z NGX =    96 NGY =   96 NGZ =   96
   dimension x,y,z NGXF=   192 NGYF=  192 NGZF=  192
   support grid    NGXF=   192 NGYF=  192 NGZF=  192
   ions per type =               2   4
   NGX,Y,Z   is equivalent  to a cutoff of  15.96, 15.96, 15.96 a.u.
   NGXF,Y,Z  is equivalent  to a cutoff of  31.92, 31.92, 31.92 a.u.

 SYSTEM =  unknown system                          
 POSCAR =  POSCAR file written by OVITO            

 Startparameter for this run:
   NWRITE =      2    write-flag & timer
   PREC   = a         normal or accurate (medium, high low for compatibility)
   ISTART =      0    job   : 0-new  1-cont  2-samecut
   ICHARG =      2    charge: 1-file 2-atom 10-const
   ISPIN  =      1    spin polarized calculation?
   LNONCOLLINEAR =      F non collinear calculations
   LSORBIT =      F    spin-orbit coupling
   INIWAV =      1    electr: 0-lowe 1-rand  2-diag
   LASPH  =      T    aspherical Exc in radial PAW
   METAGGA=      F    non-selfconsistent MetaGGA calc.

 Electronic Relaxation 1
   ENCUT  =  800.0 eV  58.80 Ry    7.67 a.u.  23.06 23.06 23.06*2*pi/ulx,y,z
   ENINI  =  800.0     initial cutoff
   ENAUG  =  605.4 eV  augmentation charge cutoff
   NELM   =     60;   NELMIN=  4; NELMDL=  0     # of ELM steps 
   EDIFF  = 0.1E-05   stopping-criterion for ELM
   LREAL  =      T    real-space projection
   NLSPLINE    = F    spline interpolate recip. space projectors
   LCOMPAT=      F    compatible to vasp.4.4
   GGA_COMPAT  = T    GGA compatible to vasp.4.4-vasp.4.6
   LMAXPAW     = -100 max onsite density
   LMAXMIX     =    2 max onsite mixed and CHGCAR
   VOSKOWN=      0    Vosko Wilk Nusair interpolation
   ROPT   =   -0.00025  -0.00025
 Ionic relaxation
   EDIFFG = 0.1E-04   stopping-criterion for IOM
   NSW    =      3    number of steps for IOM
   NBLOCK =      1;   KBLOCK =      3    inner block; outer block 
   IBRION =      0    ionic relax: 0-MD 1-quasi-New 2-CG
   NFREE  =      0    steps in history (QN), initial steepest desc. (CG)
   ISIF   =      2    stress and relaxation
   IWAVPR =     12    prediction:  0-non 1-charg 2-wave 3-comb
   ISYM   =      0    0-nonsym 1-usesym 2-fastsym
   LCORR  =      T    Harris-Foulkes like correction to forces

   POTIM  = 1.0000    time-step for ionic-motion
   TEIN   =    0.0    initial temperature
   TEBEG  =    0.0;   TEEND  =   0.0 temperature during run
   SMASS  =  -3.00    Nose mass-parameter (am)
   estimated Nose-frequenzy (Omega)   =  0.10E-29 period in steps =****** mass=  -0.229E-26a.u.
   SCALEE = 1.0000    scale energy and forces
   NPACO  =    256;   APACO  = 16.0  distance and # of slots for P.C.
   PSTRESS=    0.0 pullay stress

  Mass of Ions in am
   POMASS =  16.00  1.00
  Ionic Valenz
   ZVAL   =   6.00  1.00
  Atomic Wigner-Seitz radii
   RWIGS  =  -1.00 -1.00
  virtual crystal weights 
   VCA    =   1.00  1.00
   NELECT =      16.0000    total number of electrons
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
   EBREAK =  0.21E-07  absolut break condition
   DEPER  =   0.30     relativ break condition  

   TIME   =   0.40     timestep for ELM

  volume/ion in A,a.u.               =     166.67      1124.72
  Fermi-wavevector in a.u.,A,eV,Ry     =   0.412523  0.779555  2.315374  0.170175
  Thomas-Fermi vector in A             =   1.369550
 
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
   LEXCH   =     8    internal setting for exchange type
   VOSKOWN=      0    Vosko Wilk Nusair interpolation
   LHFCALC =     F    Hartree Fock is set to
   LHFONE  =     F    Hartree Fock one center treatment
   AEXX    =    0.0000 exact exchange contribution

 Linear response parameters
   LEPSILON=     F    determine dielectric tensor
   LRPA    =     F    only Hartree local field effects (RPA)
   LNABLA  =     F    use nabla operator in PAW spheres
   LVEL    =     F    velocity operator in full k-point grid
   LINTERFAST=   F  fast interpolation
   KINTER  =     0    interpolate to denser k-point grid
   CSHIFT  =0.1000    complex shift for real part using Kramers Kronig
   OMEGAMAX=  -1.0    maximum frequency
   DEG_THRESHOLD= 0.2000000E-02 threshold for treating states as degnerate
   RTIME   =   -0.100 relaxation time in fs
  (WPLASMAI=    0.000 imaginary part of plasma frequency in eV, 0.658/RTIME)
   DFIELD  = 0.0000000 0.0000000 0.0000000 field for delta impulse in time
 
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
 real space projection scheme for non local part
 use partial core corrections
 calculate Harris-corrections to forces 
   (improved forces if not selfconsistent)
 use gradient corrections 
 use of overlap-Matrix (Vanderbilt PP)
 Gauss-broadening in eV      SIGMA  =   0.05


--------------------------------------------------------------------------------------------------------


  energy-cutoff  :      800.00
  volume of cell :     1000.00
      direct lattice vectors                 reciprocal lattice vectors
    10.000000000  0.000000000  0.000000000     0.100000000  0.000114090 -0.001410404
    -0.011409000 10.000000000  0.000000000     0.000000000  0.100000000  0.000595569
     0.141108300 -0.059556900 10.000000000     0.000000000  0.000000000  0.100000000

  length of vectors
    10.000000000 10.000006508 10.001172860     0.100010011  0.100001773  0.100000000


 
 k-points in units of 2pi/SCALE and weight: read from INCAR                         
   0.00000000  0.00000000  0.00000000       1.000
 
 k-points in reciprocal lattice and weights: read from INCAR                         
   0.00000000  0.00000000  0.00000000       1.000
 
 position of ions in fractional coordinates (direct lattice) 
   0.42800000  0.42400000  0.52000000
   0.23000000  0.62800000  0.11300000
   0.45800000  0.35200000  0.45800000
   0.38900000  0.38400000  0.60300000
   0.13700000  0.62600000  0.15000000
   0.23100000  0.58900000  0.02100000
 
 position of ions in cartesian coordinates  (Angst):
   4.34853890  4.20903041  5.20000000
   2.30878039  6.27327007  1.13000000
   4.64061163  3.49272294  4.58000000
   3.97070725  3.80408719  6.03000000
   1.38402421  6.25106646  1.50000000
   2.30624337  5.88874931  0.21000000
 


--------------------------------------------------------------------------------------------------------


 k-point  1 :   0.0000 0.0000 0.0000  plane waves:   51381

 maximum and minimum number of plane-waves per node :     51381    51381

 maximum number of plane-waves:     51381
 maximum index in each direction: 
   IXMAX=   23   IYMAX=   23   IZMAX=   23
   IXMIN=  -23   IYMIN=  -23   IZMIN=  -23


 serial   3D FFT for wavefunctions
 parallel 3D FFT for charge:
    minimum data exchange during FFTs selected (reduces bandwidth)


 total amount of memory used by VASP MPI-rank0   719393. kBytes
=======================================================================

   base      :      30000. kBytes
   nonlr-proj:       1360. kBytes
   fftplans  :     243296. kBytes
   grid      :     434847. kBytes
   one-center:         18. kBytes
   wavefun   :       9872. kBytes
 
     INWAV:  cpu time    0.0000: real time    0.0000
 Broyden mixing: mesh for mixing (old mesh)
   NGX = 47   NGY = 47   NGZ = 47
  (NGX  =192   NGY  =192   NGZ  =192)
  gives a total of 103823 points

 initial charge density was supplied:
 charge density of overlapping atoms calculated
 number of electron      16.0000000 magnetization 
 keeping initial charge density in first step


--------------------------------------------------------------------------------------------------------


 Maximum index for non-local projection operator         3346
 Maximum index for augmentation-charges        11137 (set IRDMAX)


--------------------------------------------------------------------------------------------------------


 First call to EWALD:  gamma=   0.177
 Maximum number of real-space cells 3x 3x 3
 Maximum number of reciprocal cells 3x 3x 3

    FEWALD:  cpu time    0.0012: real time    0.0012


----------------------------------------- Iteration    1(   1)  ---------------------------------------


    POTLOK:  cpu time    3.0759: real time    3.0767
    SETDIJ:  cpu time    0.2569: real time    0.2570
     EDDAV:  cpu time    0.0922: real time    0.0922
       DOS:  cpu time    0.0001: real time    0.0001
    --------------------------------------------
      LOOP:  cpu time    3.4265: real time    3.4299

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.1910465E+03  (-0.4893166E+03)
 number of electron      16.0000000 magnetization 
 augmentation part       16.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -796.12263672
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        64.16083839
  PAW double counting   =       701.63506352     -705.26791599
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =       -57.83011466
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       191.04653371 eV

  energy without entropy =      191.04653371  energy(sigma->0) =      191.04653371


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(   2)  ---------------------------------------


     EDDAV:  cpu time    0.0815: real time    0.0815
       DOS:  cpu time    0.0000: real time    0.0000
    --------------------------------------------
      LOOP:  cpu time    0.0821: real time    0.0857

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.1270414E+03  (-0.1090457E+03)
 number of electron      16.0000000 magnetization 
 augmentation part       16.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -796.12263672
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        64.16083839
  PAW double counting   =       701.63506352     -705.26791599
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -184.87156302
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =        64.00508535 eV

  energy without entropy =       64.00508535  energy(sigma->0) =       64.00508535


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(   3)  ---------------------------------------


     EDDAV:  cpu time    0.0822: real time    0.0874
       DOS:  cpu time    0.0000: real time    0.0000
    --------------------------------------------
      LOOP:  cpu time    0.0826: real time    0.0878

 eigenvalue-minimisations  :    36
 total energy-change (2. order) :-0.9054589E+02  (-0.9028999E+02)
 number of electron      16.0000000 magnetization 
 augmentation part       16.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -796.12263672
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        64.16083839
  PAW double counting   =       701.63506352     -705.26791599
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -275.41745229
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -26.54080392 eV

  energy without entropy =      -26.54080392  energy(sigma->0) =      -26.54080392


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(   4)  ---------------------------------------


     EDDAV:  cpu time    0.0757: real time    0.0757
       DOS:  cpu time    0.0000: real time    0.0000
    --------------------------------------------
      LOOP:  cpu time    0.0762: real time    0.0762

 eigenvalue-minimisations  :    28
 total energy-change (2. order) :-0.5761346E+01  (-0.5758539E+01)
 number of electron      16.0000000 magnetization 
 augmentation part       16.0000000 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -796.12263672
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        64.16083839
  PAW double counting   =       701.63506352     -705.26791599
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -281.17879852
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -32.30215016 eV

  energy without entropy =      -32.30215016  energy(sigma->0) =      -32.30215016


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(   5)  ---------------------------------------


     EDDAV:  cpu time    0.0801: real time    0.0824
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2212: real time    0.2213
    MIXING:  cpu time    0.0968: real time    0.0968
    --------------------------------------------
      LOOP:  cpu time    0.3987: real time    0.4011

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.4119188E-01  (-0.4118307E-01)
 number of electron      16.0000002 magnetization 
 augmentation part        1.6558953 magnetization 

 Broyden mixing:
  rms(total) = 0.10981E+01    rms(broyden)= 0.10980E+01
  rms(prec ) = 0.14872E+01
  weight for this iteration     100.00

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -796.12263672
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        64.16083839
  PAW double counting   =       701.63506352     -705.26791599
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -281.21999040
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -32.34334203 eV

  energy without entropy =      -32.34334203  energy(sigma->0) =      -32.34334203


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(   6)  ---------------------------------------


    POTLOK:  cpu time    2.9395: real time    2.9430
    SETDIJ:  cpu time    0.2563: real time    0.2565
    EDDIAG:  cpu time    0.0303: real time    0.0303
  RMM-DIIS:  cpu time    0.0546: real time    0.0546
    ORTHCH:  cpu time    0.0102: real time    0.0102
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2002: real time    0.2002
    MIXING:  cpu time    0.0969: real time    0.0969
    --------------------------------------------
      LOOP:  cpu time    3.5888: real time    3.5954

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.3637436E+01  (-0.1514950E+01)
 number of electron      16.0000001 magnetization 
 augmentation part        1.4051504 magnetization 

 Broyden mixing:
  rms(total) = 0.44609E+00    rms(broyden)= 0.44607E+00
  rms(prec ) = 0.54163E+00
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   0.9523
  0.9523

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -843.30387941
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        66.76391762
  PAW double counting   =       880.91709234     -885.25738070
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -232.29695492
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.70590591 eV

  energy without entropy =      -28.70590591  energy(sigma->0) =      -28.70590591


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(   7)  ---------------------------------------


    POTLOK:  cpu time    2.9233: real time    2.9240
    SETDIJ:  cpu time    0.2519: real time    0.2536
    EDDIAG:  cpu time    0.0283: real time    0.0283
  RMM-DIIS:  cpu time    0.0552: real time    0.0552
    ORTHCH:  cpu time    0.0102: real time    0.0102
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1986: real time    0.1988
    MIXING:  cpu time    0.0994: real time    0.0994
    --------------------------------------------
      LOOP:  cpu time    3.5679: real time    3.5759

 eigenvalue-minimisations  :    26
 total energy-change (2. order) : 0.1932297E+00  (-0.7927167E-01)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3651844 magnetization 

 Broyden mixing:
  rms(total) = 0.29730E+00    rms(broyden)= 0.29730E+00
  rms(prec ) = 0.35429E+00
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.6726
  1.2046  2.1406

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -852.30152685
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.24163390
  PAW double counting   =       968.99589907     -973.29683957
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -223.62314191
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.51267621 eV

  energy without entropy =      -28.51267621  energy(sigma->0) =      -28.51267621


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(   8)  ---------------------------------------


    POTLOK:  cpu time    2.8995: real time    2.9002
    SETDIJ:  cpu time    0.2523: real time    0.2523
    EDDIAG:  cpu time    0.0281: real time    0.0281
  RMM-DIIS:  cpu time    0.0533: real time    0.0533
    ORTHCH:  cpu time    0.0104: real time    0.0104
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1989: real time    0.2016
    MIXING:  cpu time    0.1019: real time    0.1020
    --------------------------------------------
      LOOP:  cpu time    3.5455: real time    3.5582

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.1371898E+00  (-0.3377750E-01)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3922697 magnetization 

 Broyden mixing:
  rms(total) = 0.83917E-01    rms(broyden)= 0.83915E-01
  rms(prec ) = 0.12897E+00
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.5251
  2.4707  1.0523  1.0523

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -855.92978989
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.42385433
  PAW double counting   =      1065.53034689    -1069.63559436
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -220.23560252
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.37548639 eV

  energy without entropy =      -28.37548639  energy(sigma->0) =      -28.37548639


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(   9)  ---------------------------------------


    POTLOK:  cpu time    2.8932: real time    2.8939
    SETDIJ:  cpu time    0.2515: real time    0.2516
    EDDIAG:  cpu time    0.0265: real time    0.0265
  RMM-DIIS:  cpu time    0.0509: real time    0.0509
    ORTHCH:  cpu time    0.0098: real time    0.0098
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1970: real time    0.1972
    MIXING:  cpu time    0.1041: real time    0.1041
    --------------------------------------------
      LOOP:  cpu time    3.5341: real time    3.5394

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.2936671E-01  (-0.1555725E-01)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3687098 magnetization 

 Broyden mixing:
  rms(total) = 0.26359E-01    rms(broyden)= 0.26358E-01
  rms(prec ) = 0.52331E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.8276
  2.6430  2.6430  1.0122  1.0122

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -862.10052074
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.74930083
  PAW double counting   =      1084.01965559    -1088.18860120
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.29725332
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.34611968 eV

  energy without entropy =      -28.34611968  energy(sigma->0) =      -28.34611968


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  10)  ---------------------------------------


    POTLOK:  cpu time    2.9093: real time    2.9102
    SETDIJ:  cpu time    0.2527: real time    0.2527
    EDDIAG:  cpu time    0.0266: real time    0.0266
  RMM-DIIS:  cpu time    0.0549: real time    0.0549
    ORTHCH:  cpu time    0.0098: real time    0.0098
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2034: real time    0.2034
    MIXING:  cpu time    0.1079: real time    0.1079
    --------------------------------------------
      LOOP:  cpu time    3.5658: real time    3.5781

 eigenvalue-minimisations  :    26
 total energy-change (2. order) :-0.2226041E-01  (-0.4359631E-02)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3588466 magnetization 

 Broyden mixing:
  rms(total) = 0.44676E-01    rms(broyden)= 0.44676E-01
  rms(prec ) = 0.58176E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.7242
  3.1260  2.4256  1.0139  1.0139  1.0416

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -864.34269587
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.79301251
  PAW double counting   =      1060.14710643    -1064.33153975
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.10556258
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.36838009 eV

  energy without entropy =      -28.36838009  energy(sigma->0) =      -28.36838009


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  11)  ---------------------------------------


    POTLOK:  cpu time    2.8423: real time    2.8437
    SETDIJ:  cpu time    0.2525: real time    0.2544
    EDDIAG:  cpu time    0.0272: real time    0.0272
  RMM-DIIS:  cpu time    0.0512: real time    0.0513
    ORTHCH:  cpu time    0.0096: real time    0.0096
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1950: real time    0.1951
    MIXING:  cpu time    0.1097: real time    0.1097
    --------------------------------------------
      LOOP:  cpu time    3.4887: real time    3.5000

 eigenvalue-minimisations  :    24
 total energy-change (2. order) :-0.1667795E-02  (-0.2842950E-02)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3687021 magnetization 

 Broyden mixing:
  rms(total) = 0.99433E-02    rms(broyden)= 0.99424E-02
  rms(prec ) = 0.18060E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.0304
  4.2185  2.6667  2.1826  1.0496  1.0496  1.0153

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -862.99284124
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.70332864
  PAW double counting   =      1061.49507687    -1065.64156863
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -213.40534269
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.37004789 eV

  energy without entropy =      -28.37004789  energy(sigma->0) =      -28.37004789


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  12)  ---------------------------------------


    POTLOK:  cpu time    2.7550: real time    2.7557
    SETDIJ:  cpu time    0.2527: real time    0.2550
    EDDIAG:  cpu time    0.0268: real time    0.0268
  RMM-DIIS:  cpu time    0.0513: real time    0.0513
    ORTHCH:  cpu time    0.0097: real time    0.0099
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1969: real time    0.1970
    MIXING:  cpu time    0.1121: real time    0.1121
    --------------------------------------------
      LOOP:  cpu time    3.4055: real time    3.4142

 eigenvalue-minimisations  :    24
 total energy-change (2. order) :-0.1043649E-01  (-0.3541362E-03)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3699881 magnetization 

 Broyden mixing:
  rms(total) = 0.99041E-02    rms(broyden)= 0.99041E-02
  rms(prec ) = 0.14274E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.0156
  4.5190  2.8008  2.4900  1.1388  1.0282  1.0663  1.0663

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.02068799
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.67683584
  PAW double counting   =      1062.01362580    -1066.15927834
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -213.36227886
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38048438 eV

  energy without entropy =      -28.38048438  energy(sigma->0) =      -28.38048438


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  13)  ---------------------------------------


    POTLOK:  cpu time    2.8211: real time    2.8219
    SETDIJ:  cpu time    0.2525: real time    0.2525
    EDDIAG:  cpu time    0.0268: real time    0.0268
  RMM-DIIS:  cpu time    0.0526: real time    0.0526
    ORTHCH:  cpu time    0.0098: real time    0.0099
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2007: real time    0.2008
    MIXING:  cpu time    0.1153: real time    0.1153
    --------------------------------------------
      LOOP:  cpu time    3.4799: real time    3.4876

 eigenvalue-minimisations  :    24
 total energy-change (2. order) :-0.2221504E-02  (-0.1412282E-03)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3681773 magnetization 

 Broyden mixing:
  rms(total) = 0.28751E-02    rms(broyden)= 0.28751E-02
  rms(prec ) = 0.53233E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.2084
  5.9118  3.0799  2.5556  2.0408  1.0680  1.0680  1.0282  0.9151

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.46092575
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.69609028
  PAW double counting   =      1064.06582480    -1068.21582563
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.93916874
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38270588 eV

  energy without entropy =      -28.38270588  energy(sigma->0) =      -28.38270588


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  14)  ---------------------------------------


    POTLOK:  cpu time    2.7183: real time    2.7190
    SETDIJ:  cpu time    0.2522: real time    0.2569
    EDDIAG:  cpu time    0.0265: real time    0.0265
  RMM-DIIS:  cpu time    0.0543: real time    0.0543
    ORTHCH:  cpu time    0.0096: real time    0.0096
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1958: real time    0.1958
    MIXING:  cpu time    0.1182: real time    0.1182
    --------------------------------------------
      LOOP:  cpu time    3.3759: real time    3.3844

 eigenvalue-minimisations  :    26
 total energy-change (2. order) :-0.2422282E-02  (-0.7992793E-04)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3671981 magnetization 

 Broyden mixing:
  rms(total) = 0.13901E-02    rms(broyden)= 0.13900E-02
  rms(prec ) = 0.24688E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.3391
  6.5906  3.7212  2.5234  2.5234  1.5796  1.0727  1.0727  1.0191  0.9491

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.59232358
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.69928625
  PAW double counting   =      1063.96123117    -1068.11176626
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.81285491
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38512817 eV

  energy without entropy =      -28.38512817  energy(sigma->0) =      -28.38512817


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  15)  ---------------------------------------


    POTLOK:  cpu time    2.7420: real time    2.7427
    SETDIJ:  cpu time    0.2518: real time    0.2518
    EDDIAG:  cpu time    0.0264: real time    0.0265
  RMM-DIIS:  cpu time    0.0510: real time    0.0510
    ORTHCH:  cpu time    0.0097: real time    0.0097
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1967: real time    0.1967
    MIXING:  cpu time    0.1253: real time    0.1254
    --------------------------------------------
      LOOP:  cpu time    3.4038: real time    3.4062

 eigenvalue-minimisations  :    24
 total energy-change (2. order) :-0.8270306E-03  (-0.1739563E-04)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3668864 magnetization 

 Broyden mixing:
  rms(total) = 0.21404E-02    rms(broyden)= 0.21403E-02
  rms(prec ) = 0.28094E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.3756
  7.3413  4.1647  2.7379  2.4587  2.0164  1.0672  1.0672  1.0158  1.0158  0.8714

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.58211573
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.69804759
  PAW double counting   =      1063.37165485    -1067.52296615
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.82187492
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38595520 eV

  energy without entropy =      -28.38595520  energy(sigma->0) =      -28.38595520


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  16)  ---------------------------------------


    POTLOK:  cpu time    2.8218: real time    2.8225
    SETDIJ:  cpu time    0.2541: real time    0.2550
    EDDIAG:  cpu time    0.0267: real time    0.0267
  RMM-DIIS:  cpu time    0.0509: real time    0.0509
    ORTHCH:  cpu time    0.0097: real time    0.0097
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1974: real time    0.1974
    MIXING:  cpu time    0.1246: real time    0.1246
    --------------------------------------------
      LOOP:  cpu time    3.4861: real time    3.4912

 eigenvalue-minimisations  :    24
 total energy-change (2. order) :-0.1635372E-03  (-0.8609277E-05)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3673480 magnetization 

 Broyden mixing:
  rms(total) = 0.27160E-03    rms(broyden)= 0.27154E-03
  rms(prec ) = 0.46293E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.4604
  7.8696  4.7872  2.8985  2.3702  2.3702  1.7675  1.0665  1.0665  1.0186  0.9247
  0.9247

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.49640326
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.69289760
  PAW double counting   =      1063.24654029    -1067.39635964
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.90409289
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38611873 eV

  energy without entropy =      -28.38611873  energy(sigma->0) =      -28.38611873


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  17)  ---------------------------------------


    POTLOK:  cpu time    2.8014: real time    2.8023
    SETDIJ:  cpu time    0.2503: real time    0.2503
    EDDIAG:  cpu time    0.0264: real time    0.0264
  RMM-DIIS:  cpu time    0.0487: real time    0.0487
    ORTHCH:  cpu time    0.0097: real time    0.0097
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1957: real time    0.1958
    MIXING:  cpu time    0.1295: real time    0.1295
    --------------------------------------------
      LOOP:  cpu time    3.4627: real time    3.4667

 eigenvalue-minimisations  :    22
 total energy-change (2. order) :-0.7783725E-04  (-0.1038845E-05)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3674069 magnetization 

 Broyden mixing:
  rms(total) = 0.24223E-03    rms(broyden)= 0.24222E-03
  rms(prec ) = 0.34042E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.5409
  8.4098  5.2725  3.4087  2.7079  2.4275  1.9629  1.0649  1.0649  1.2535  1.0296
  0.9445  0.9445

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.49545982
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.69253900
  PAW double counting   =      1063.27757726    -1067.42740388
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.90474829
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38619657 eV

  energy without entropy =      -28.38619657  energy(sigma->0) =      -28.38619657


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  18)  ---------------------------------------


    POTLOK:  cpu time    2.8524: real time    2.8532
    SETDIJ:  cpu time    0.2543: real time    0.2545
    EDDIAG:  cpu time    0.0264: real time    0.0264
  RMM-DIIS:  cpu time    0.0457: real time    0.0457
    ORTHCH:  cpu time    0.0097: real time    0.0097
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2013: real time    0.2013
    MIXING:  cpu time    0.1325: real time    0.1325
    --------------------------------------------
      LOOP:  cpu time    3.5233: real time    3.5292

 eigenvalue-minimisations  :    19
 total energy-change (2. order) :-0.2251838E-04  (-0.3818942E-06)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3674063 magnetization 

 Broyden mixing:
  rms(total) = 0.17070E-03    rms(broyden)= 0.17070E-03
  rms(prec ) = 0.23427E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.5265
  8.5312  5.3854  3.6257  2.5898  2.5898  2.2236  1.8139  1.0660  1.0660  0.9248
  0.9657  1.0309  1.0309

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.49991904
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.69274234
  PAW double counting   =      1063.33745782    -1067.48734229
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.90045708
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38621909 eV

  energy without entropy =      -28.38621909  energy(sigma->0) =      -28.38621909


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  19)  ---------------------------------------


    POTLOK:  cpu time    2.9464: real time    2.9497
    SETDIJ:  cpu time    0.2529: real time    0.2530
    EDDIAG:  cpu time    0.0279: real time    0.0279
  RMM-DIIS:  cpu time    0.0440: real time    0.0440
    ORTHCH:  cpu time    0.0102: real time    0.0102
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2040: real time    0.2040
    MIXING:  cpu time    0.1382: real time    0.1383
    --------------------------------------------
      LOOP:  cpu time    3.6246: real time    3.6317

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.5438037E-05  (-0.7792679E-07)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3673849 magnetization 

 Broyden mixing:
  rms(total) = 0.70572E-04    rms(broyden)= 0.70572E-04
  rms(prec ) = 0.89229E-04
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.5120
  8.5975  5.4679  3.9917  2.7571  2.7571  2.4641  1.9388  1.0656  1.0656  1.1860
  1.0509  0.9864  0.9300  0.9094

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.50821217
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.69318789
  PAW double counting   =      1063.35438045    -1067.50433405
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.89254581
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38622453 eV

  energy without entropy =      -28.38622453  energy(sigma->0) =      -28.38622453


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  20)  ---------------------------------------


    POTLOK:  cpu time    2.8871: real time    2.8883
    SETDIJ:  cpu time    0.2516: real time    0.2516
    EDDIAG:  cpu time    0.0277: real time    0.0277
  RMM-DIIS:  cpu time    0.0427: real time    0.0428
    ORTHCH:  cpu time    0.0102: real time    0.0102
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1983: real time    0.2009
    MIXING:  cpu time    0.1398: real time    0.1399
    --------------------------------------------
      LOOP:  cpu time    3.5583: real time    3.5666

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.1323549E-05  (-0.6673245E-07)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3673553 magnetization 

 Broyden mixing:
  rms(total) = 0.73789E-04    rms(broyden)= 0.73787E-04
  rms(prec ) = 0.96084E-04
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   2.4769
  8.6072  5.5378  4.0716  3.0134  2.5385  2.5385  2.0145  1.6907  1.0664  1.0664
  1.0790  1.0790  0.9371  0.9371  0.9766

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.51302974
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.69346062
  PAW double counting   =      1063.35079228    -1067.50085892
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.88788925
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38622585 eV

  energy without entropy =      -28.38622585  energy(sigma->0) =      -28.38622585


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    1(  21)  ---------------------------------------


    POTLOK:  cpu time    2.9326: real time    2.9333
    SETDIJ:  cpu time    0.2521: real time    0.2521
    EDDIAG:  cpu time    0.0274: real time    0.0274
  RMM-DIIS:  cpu time    0.0434: real time    0.0434
    ORTHCH:  cpu time    0.0100: real time    0.0100
       DOS:  cpu time    0.0000: real time    0.0000
    --------------------------------------------
      LOOP:  cpu time    3.2664: real time    3.2707

 eigenvalue-minimisations  :    12
 total energy-change (2. order) :-0.3901877E-06  (-0.1066007E-07)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3673553 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        68.87223147
  -Hartree energ DENC   =      -863.50988411
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.69329751
  PAW double counting   =      1063.34603123    -1067.49607372
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -212.89089632
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.38622624 eV

  energy without entropy =      -28.38622624  energy(sigma->0) =      -28.38622624


--------------------------------------------------------------------------------------------------------




 average (electrostatic) potential at core
  the test charge radii are     0.7215  0.5201
  (the norm of the test charge is              1.0000)
       1 -80.2716       2 -80.2675       3 -44.1621       4 -44.1278       5 -44.1080
       6 -44.1197
 
 
 
 E-fermi :  -6.5950     XC(G=0):  -1.0103     alpha+bet : -0.2602


 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -24.7178      2.00000
      2     -24.6963      2.00000
      3     -12.8876      2.00000
      4     -12.8419      2.00000
      5      -8.7555      2.00000
      6      -8.7412      2.00000
      7      -6.9128      2.00000
      8      -6.8922      2.00000
      9      -1.2066      0.00000
     10      -0.4645      0.00000
     11       0.4428      0.00000
     12       0.8350      0.00000


--------------------------------------------------------------------------------------------------------


 soft charge-density along one line, spin component           1
         0         1         2         3         4         5         6         7         8         9
 total charge-density along one line
 
 pseudopotential strength for first ion, spin component:           1
 13.767 -16.875   0.094  -0.017   0.006  -0.115   0.020  -0.008
-16.875  20.716  -0.120   0.021  -0.008   0.147  -0.026   0.010
  0.094  -0.120 -10.421  -0.018   0.009  12.902   0.023  -0.012
 -0.017   0.021  -0.018 -10.478   0.077   0.023  12.977  -0.102
  0.006  -0.008   0.009   0.077 -10.347  -0.012  -0.102  12.803
 -0.115   0.147  12.902   0.023  -0.012 -15.898  -0.031   0.016
  0.020  -0.026   0.023  12.977  -0.102  -0.031 -15.998   0.135
 -0.008   0.010  -0.012  -0.102  12.803   0.016   0.135 -15.768
 total augmentation occupancy for first ion, spin component:           1
  2.665   0.410  -0.377   0.065  -0.025  -0.149   0.026  -0.010
  0.410   0.113  -0.383   0.069  -0.027  -0.065   0.012  -0.005
 -0.377  -0.383   2.253   0.031  -0.015   0.313   0.029  -0.015
  0.065   0.069   0.031   2.339  -0.098   0.029   0.408  -0.126
 -0.025  -0.027  -0.015  -0.098   2.171  -0.015  -0.126   0.193
 -0.149  -0.065   0.313   0.029  -0.015   0.047   0.008  -0.004
  0.026   0.012   0.029   0.408  -0.126   0.008   0.077  -0.032
 -0.010  -0.005  -0.015  -0.126   0.193  -0.004  -0.032   0.023


------------------------ aborting loop because EDIFF is reached ----------------------------------------


    CHARGE:  cpu time    0.2170: real time    0.2172
    FORLOC:  cpu time    0.3375: real time    0.3375
    FORNL :  cpu time    0.0408: real time    0.0408
    STRESS:  cpu time    0.3686: real time    0.3687
    FORCOR:  cpu time    3.0609: real time    3.0617
    FORHAR:  cpu time    0.9477: real time    0.9481
    MIXING:  cpu time    0.1451: real time    0.1451
    OFIELD:  cpu time    0.0001: real time    0.0001

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Alpha Z     0.89701     0.89701     0.89701
  Ewald     -37.78103   -62.17243   168.82554   -25.04698    24.58938   -60.21377
  Hartree   246.06036   242.91145   374.53304   -19.50123   -10.77994     4.77240
  E(xc)     -72.04402   -72.15038   -71.64995    -0.01804     0.18563    -0.34030
  Local    -397.16945  -376.25165  -706.20338    43.57899    -1.35440    33.03151
  n-local   -60.60872   -60.88565   -59.86597    -0.08060     0.28357    -0.64219
  augment    14.19230    14.61316    12.44467     0.03181    -0.84087     1.49832
  Kinetic   305.10655   312.44204   277.99182     1.25728   -12.69392    23.39258
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  -------------------------------------------------------------------------------------
  Total      -1.34701    -0.59646    -3.02722     0.22123    -0.61055     1.49855
  in kB      -2.15814    -0.95564    -4.85014     0.35445    -0.97821     2.40095
  external pressure =       -2.65 kB  Pullay stress =        0.00 kB

  kinetic pressure (ideal gas correction) =      0.00 kB
  total pressure  =     -2.65 kB
  Total+kin.    -2.158      -0.956      -4.850       0.354      -0.978       2.401

 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
  energy-cutoff  :      800.00
  volume of cell :     1000.00
      direct lattice vectors                 reciprocal lattice vectors
    10.000000000  0.000000000  0.000000000     0.100000000  0.000114090 -0.001410404
    -0.011409000 10.000000000  0.000000000     0.000000000  0.100000000  0.000595569
     0.141108300 -0.059556900 10.000000000     0.000000000  0.000000000  0.100000000

  length of vectors
    10.000000000 10.000006508 10.001172860     0.100010011  0.100001773  0.100000000


 FORCES acting on ions
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction
 -----------------------------------------------------------------------------------------------
   -.133E+02 -.440E+02 0.707E+00   0.164E+02 0.891E+02 -.876E+01   -.326E+01 -.457E+02 0.841E+01   0.835E-05 0.421E-03 -.103E-03
   -.346E+02 -.289E+02 -.170E+02   0.715E+02 0.451E+02 0.390E+02   -.375E+02 -.164E+02 -.223E+02   0.318E-03 0.620E-04 0.215E-03
   -.277E+02 0.611E+02 0.537E+02   0.296E+02 -.660E+02 -.579E+02   -.239E+01 0.535E+01 0.509E+01   -.461E-04 0.140E-03 0.801E-04
   0.314E+02 0.329E+02 -.714E+02   -.338E+02 -.356E+02 0.766E+02   0.293E+01 0.275E+01 -.642E+01   0.533E-04 0.899E-04 -.131E-03
   0.792E+02 -.820E+00 -.327E+02   -.851E+02 0.646E+00 0.350E+02   0.695E+01 0.251E-01 -.313E+01   0.163E-03 0.528E-06 -.400E-04
   -.148E+01 0.308E+02 0.781E+02   0.139E+01 -.332E+02 -.839E+02   -.323E+00 0.287E+01 0.702E+01   0.215E-04 0.554E-04 0.160E-03
 -----------------------------------------------------------------------------------------------
   0.336E+02 0.511E+02 0.114E+02   -.195E-13 0.000E+00 0.142E-13   -.336E+02 -.511E+02 -.114E+02   0.519E-03 0.768E-03 0.180E-03
 
 
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      4.34854      4.20903      5.20000        -0.141986     -0.573729      0.356515
      2.30878      6.27327      1.13000        -0.603569     -0.289019     -0.411629
      4.64061      3.49272      4.58000        -0.446503      0.453279      0.921972
      3.97071      3.80409      6.03000         0.584218      0.121041     -1.267066
      1.38402      6.25107      1.50000         1.014160     -0.149409     -0.794823
      2.30624      5.88875      0.21000        -0.412170      0.435811      1.203659
 -----------------------------------------------------------------------------------
    total drift:                               -0.005850     -0.002026      0.008629


--------------------------------------------------------------------------------------------------------



  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -28.38622624 eV

  energy  without entropy=      -28.38622624  energy(sigma->0) =      -28.38622624
 


--------------------------------------------------------------------------------------------------------


    POTLOK:  cpu time    3.1955: real time    3.1968


--------------------------------------------------------------------------------------------------------


           RANDOM_SEED =         273936164                0                0
           RANDOM_SEED =         273936164               36                0
   IONSTEP:  cpu time    0.0020: real time    0.0045

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -28.386226  see above
  kinetic energy EKIN   =         0.008168
  kin. lattice  EKIN_LAT=         0.000000  (temperature   12.64 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -28.378059 eV

  maximum distance moved by ions :      0.14E-02

    WAVPRE:  cpu time    1.0518: real time    1.0649
    FEWALD:  cpu time    0.0004: real time    0.0004
    ORTHCH:  cpu time    0.0180: real time    0.0180
     LOOP+:  cpu time   69.5910: real time   69.7456


----------------------------------------- Iteration    2(   1)  ---------------------------------------


    POTLOK:  cpu time    3.0137: real time    3.0147
    SETDIJ:  cpu time    0.2524: real time    0.2524
     EDDAV:  cpu time    0.0835: real time    0.0835
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2006: real time    0.2008
    MIXING:  cpu time    0.0943: real time    0.0943
    --------------------------------------------
      LOOP:  cpu time    3.6456: real time    3.6475

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.5387530E-01  (-0.2415222E-01)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3806988 magnetization 

 Broyden mixing:
  rms(total) = 0.29143E-01    rms(broyden)= 0.29142E-01
  rms(prec ) = 0.35627E-01
  weight for this iteration     100.00

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -865.87808935
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.80887378
  PAW double counting   =      1063.34382299    -1067.49383233
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -215.14192616
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.44010115 eV

  energy without entropy =      -28.44010115  energy(sigma->0) =      -28.44010115


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(   2)  ---------------------------------------


    POTLOK:  cpu time    2.9865: real time    2.9873
    SETDIJ:  cpu time    0.2514: real time    0.2514
    EDDIAG:  cpu time    0.0271: real time    0.0271
  RMM-DIIS:  cpu time    0.0517: real time    0.0517
    ORTHCH:  cpu time    0.0098: real time    0.0099
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1935: real time    0.1935
    MIXING:  cpu time    0.0964: real time    0.0965
    --------------------------------------------
      LOOP:  cpu time    3.6173: real time    3.6215

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.1255625E-02  (-0.7445992E-03)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3786525 magnetization 

 Broyden mixing:
  rms(total) = 0.13118E-01    rms(broyden)= 0.13117E-01
  rms(prec ) = 0.15202E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.4738
  1.4738

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -866.77746243
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.85949281
  PAW double counting   =      1069.65002259    -1073.81483896
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.27710947
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43884553 eV

  energy without entropy =      -28.43884553  energy(sigma->0) =      -28.43884553


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(   3)  ---------------------------------------


    POTLOK:  cpu time    3.0259: real time    3.0267
    SETDIJ:  cpu time    0.2544: real time    0.2545
    EDDIAG:  cpu time    0.0271: real time    0.0271
  RMM-DIIS:  cpu time    0.0530: real time    0.0531
    ORTHCH:  cpu time    0.0100: real time    0.0100
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2091: real time    0.2091
    MIXING:  cpu time    0.1012: real time    0.1012
    --------------------------------------------
      LOOP:  cpu time    3.6816: real time    3.6861

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.1203537E-03  (-0.1714170E-03)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3765009 magnetization 

 Broyden mixing:
  rms(total) = 0.60859E-02    rms(broyden)= 0.60858E-02
  rms(prec ) = 0.70364E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.5869
  1.5869  1.5869

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -867.15347103
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.87997734
  PAW double counting   =      1074.53843662    -1078.71295176
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -213.91176627
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43872518 eV

  energy without entropy =      -28.43872518  energy(sigma->0) =      -28.43872518


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(   4)  ---------------------------------------


    POTLOK:  cpu time    2.9837: real time    2.9846
    SETDIJ:  cpu time    0.2511: real time    0.2511
    EDDIAG:  cpu time    0.0269: real time    0.0269
  RMM-DIIS:  cpu time    0.0515: real time    0.0515
    ORTHCH:  cpu time    0.0100: real time    0.0100
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1957: real time    0.1991
    MIXING:  cpu time    0.1006: real time    0.1006
    --------------------------------------------
      LOOP:  cpu time    3.6205: real time    3.6306

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.2060560E-04  (-0.4697320E-04)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3777278 magnetization 

 Broyden mixing:
  rms(total) = 0.34673E-02    rms(broyden)= 0.34673E-02
  rms(prec ) = 0.42101E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.5033
  2.4224  0.9617  1.1258

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -866.90641366
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.86587347
  PAW double counting   =      1075.71946701    -1079.89116515
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.14751617
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43870457 eV

  energy without entropy =      -28.43870457  energy(sigma->0) =      -28.43870457


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(   5)  ---------------------------------------


    POTLOK:  cpu time    2.9176: real time    2.9185
    SETDIJ:  cpu time    0.2533: real time    0.2533
    EDDIAG:  cpu time    0.0266: real time    0.0266
  RMM-DIIS:  cpu time    0.0509: real time    0.0509
    ORTHCH:  cpu time    0.0097: real time    0.0097
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1968: real time    0.1968
    MIXING:  cpu time    0.1039: real time    0.1040
    --------------------------------------------
      LOOP:  cpu time    3.5599: real time    3.5642

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.1085067E-04  (-0.4423497E-05)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3774431 magnetization 

 Broyden mixing:
  rms(total) = 0.10704E-02    rms(broyden)= 0.10704E-02
  rms(prec ) = 0.13948E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.5596
  2.3616  1.4605  1.4605  0.9558

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -866.99229218
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.87039168
  PAW double counting   =      1077.24164905    -1081.41645049
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.06304170
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43869372 eV

  energy without entropy =      -28.43869372  energy(sigma->0) =      -28.43869372


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(   6)  ---------------------------------------


    POTLOK:  cpu time    3.0205: real time    3.0215
    SETDIJ:  cpu time    0.2549: real time    0.2550
    EDDIAG:  cpu time    0.0266: real time    0.0266
  RMM-DIIS:  cpu time    0.0522: real time    0.0522
    ORTHCH:  cpu time    0.0097: real time    0.0097
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1992: real time    0.1993
    MIXING:  cpu time    0.1073: real time    0.1073
    --------------------------------------------
      LOOP:  cpu time    3.6714: real time    3.6754

 eigenvalue-minimisations  :    24
 total energy-change (2. order) :-0.2373588E-04  (-0.9466176E-05)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3768386 magnetization 

 Broyden mixing:
  rms(total) = 0.15018E-02    rms(broyden)= 0.15018E-02
  rms(prec ) = 0.19199E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.5250
  2.5268  2.1219  1.0609  0.9578  0.9578

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -867.06246967
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.87421797
  PAW double counting   =      1077.39489230    -1081.57131017
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -213.99509780
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43871746 eV

  energy without entropy =      -28.43871746  energy(sigma->0) =      -28.43871746


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(   7)  ---------------------------------------


    POTLOK:  cpu time    2.9219: real time    2.9226
    SETDIJ:  cpu time    0.2516: real time    0.2517
    EDDIAG:  cpu time    0.0263: real time    0.0264
  RMM-DIIS:  cpu time    0.0509: real time    0.0509
    ORTHCH:  cpu time    0.0097: real time    0.0097
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1936: real time    0.1936
    MIXING:  cpu time    0.1098: real time    0.1099
    --------------------------------------------
      LOOP:  cpu time    3.5648: real time    3.5689

 eigenvalue-minimisations  :    24
 total energy-change (2. order) :-0.3371101E-05  (-0.1906500E-05)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3770898 magnetization 

 Broyden mixing:
  rms(total) = 0.27232E-03    rms(broyden)= 0.27231E-03
  rms(prec ) = 0.39188E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.6406
  2.7862  2.1680  1.8776  0.9941  1.0088  1.0088

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -867.00687766
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.87106464
  PAW double counting   =      1077.38245099    -1081.55806234
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.04834638
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43872083 eV

  energy without entropy =      -28.43872083  energy(sigma->0) =      -28.43872083


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(   8)  ---------------------------------------


    POTLOK:  cpu time    2.9373: real time    2.9405
    SETDIJ:  cpu time    0.2533: real time    0.2535
    EDDIAG:  cpu time    0.0277: real time    0.0277
  RMM-DIIS:  cpu time    0.0476: real time    0.0476
    ORTHCH:  cpu time    0.0103: real time    0.0103
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1976: real time    0.1976
    MIXING:  cpu time    0.1125: real time    0.1126
    --------------------------------------------
      LOOP:  cpu time    3.5872: real time    3.5939

 eigenvalue-minimisations  :    19
 total energy-change (2. order) :-0.1264246E-04  (-0.2907987E-06)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3771518 magnetization 

 Broyden mixing:
  rms(total) = 0.58566E-04    rms(broyden)= 0.58553E-04
  rms(prec ) = 0.10893E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.6468
  3.1306  2.4200  2.0371  1.0043  1.0043  1.0047  0.9263

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -866.99626030
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.87045672
  PAW double counting   =      1077.36592863    -1081.54137048
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.05853797
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43873347 eV

  energy without entropy =      -28.43873347  energy(sigma->0) =      -28.43873347


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(   9)  ---------------------------------------


    POTLOK:  cpu time    3.0067: real time    3.0073
    SETDIJ:  cpu time    0.2518: real time    0.2518
    EDDIAG:  cpu time    0.0279: real time    0.0279
  RMM-DIIS:  cpu time    0.0447: real time    0.0447
    ORTHCH:  cpu time    0.0104: real time    0.0104
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2027: real time    0.2055
    MIXING:  cpu time    0.1160: real time    0.1160
    --------------------------------------------
      LOOP:  cpu time    3.6611: real time    3.6672

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.5852328E-05  (-0.1773891E-06)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3772001 magnetization 

 Broyden mixing:
  rms(total) = 0.15794E-03    rms(broyden)= 0.15794E-03
  rms(prec ) = 0.21165E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.6256
  3.4966  2.4762  2.0532  0.9985  0.9985  1.1653  0.9080  0.9080

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -866.98913970
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.87005298
  PAW double counting   =      1077.34664919    -1081.52193624
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.06541548
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43873932 eV

  energy without entropy =      -28.43873932  energy(sigma->0) =      -28.43873932


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(  10)  ---------------------------------------


    POTLOK:  cpu time    2.9512: real time    2.9521
    SETDIJ:  cpu time    0.2525: real time    0.2526
    EDDIAG:  cpu time    0.0273: real time    0.0273
  RMM-DIIS:  cpu time    0.0427: real time    0.0427
    ORTHCH:  cpu time    0.0099: real time    0.0099
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1939: real time    0.1939
    MIXING:  cpu time    0.1172: real time    0.1172
    --------------------------------------------
      LOOP:  cpu time    3.5959: real time    3.6013

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.1933179E-05  (-0.4088900E-07)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3772096 magnetization 

 Broyden mixing:
  rms(total) = 0.16818E-03    rms(broyden)= 0.16818E-03
  rms(prec ) = 0.22168E-03
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.7877
  4.3163  2.6998  2.0625  2.0625  1.0446  1.0446  0.9797  0.9396  0.9396

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -866.99074010
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.87013419
  PAW double counting   =      1077.34974543    -1081.52504100
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.06388970
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43874125 eV

  energy without entropy =      -28.43874125  energy(sigma->0) =      -28.43874125


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(  11)  ---------------------------------------


    POTLOK:  cpu time    3.0018: real time    3.0026
    SETDIJ:  cpu time    0.2550: real time    0.2551
    EDDIAG:  cpu time    0.0272: real time    0.0272
  RMM-DIIS:  cpu time    0.0439: real time    0.0440
    ORTHCH:  cpu time    0.0101: real time    0.0101
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2108: real time    0.2108
    MIXING:  cpu time    0.1249: real time    0.1249
    --------------------------------------------
      LOOP:  cpu time    3.6747: real time    3.6798

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.2331584E-05  (-0.5528636E-07)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3771833 magnetization 

 Broyden mixing:
  rms(total) = 0.31122E-04    rms(broyden)= 0.31120E-04
  rms(prec ) = 0.45444E-04
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.7791
  4.6497  2.7450  2.1550  2.1550  1.1645  1.1645  0.9643  0.9643  0.9907  0.8377

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -866.99908767
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.87057744
  PAW double counting   =      1077.36146830    -1081.53687868
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.05587290
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43874359 eV

  energy without entropy =      -28.43874359  energy(sigma->0) =      -28.43874359


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    2(  12)  ---------------------------------------


    POTLOK:  cpu time    3.0825: real time    3.0834
    SETDIJ:  cpu time    0.2512: real time    0.2512
    EDDIAG:  cpu time    0.0274: real time    0.0274
  RMM-DIIS:  cpu time    0.0426: real time    0.0426
    ORTHCH:  cpu time    0.0100: real time    0.0100
       DOS:  cpu time    0.0000: real time    0.0000
    --------------------------------------------
      LOOP:  cpu time    3.4146: real time    3.4197

 eigenvalue-minimisations  :    14
 total energy-change (2. order) :-0.4715274E-06  (-0.1483489E-07)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3771833 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        73.32198224
  -Hartree energ DENC   =      -867.00003177
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        67.87061708
  PAW double counting   =      1077.36177636    -1081.53720451
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -214.05495115
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.43874406 eV

  energy without entropy =      -28.43874406  energy(sigma->0) =      -28.43874406


--------------------------------------------------------------------------------------------------------




 average (electrostatic) potential at core
  the test charge radii are     0.7215  0.5201
  (the norm of the test charge is              1.0000)
       1 -80.2844       2 -80.2802       3 -44.2788       4 -44.2719       5 -44.2440
       6 -44.2600
 
 
 
 E-fermi :  -6.6366     XC(G=0):  -1.0177     alpha+bet : -0.2602


 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -24.8475      2.00000
      2     -24.8332      2.00000
      3     -12.9761      2.00000
      4     -12.9389      2.00000
      5      -8.7929      2.00000
      6      -8.7795      2.00000
      7      -6.9397      2.00000
      8      -6.9197      2.00000
      9      -1.1725      0.00000
     10      -0.4176      0.00000
     11       0.4670      0.00000
     12       0.8212      0.00000


--------------------------------------------------------------------------------------------------------


 soft charge-density along one line, spin component           1
         0         1         2         3         4         5         6         7         8         9
 total charge-density along one line
 
 pseudopotential strength for first ion, spin component:           1
 13.768 -16.877   0.097  -0.018   0.007  -0.119   0.021  -0.008
-16.877  20.718  -0.123   0.022  -0.009   0.151  -0.027   0.011
  0.097  -0.123 -10.422  -0.018   0.009  12.902   0.023  -0.012
 -0.018   0.022  -0.018 -10.480   0.079   0.023  12.979  -0.104
  0.007  -0.009   0.009   0.079 -10.346  -0.012  -0.104  12.802
 -0.119   0.151  12.902   0.023  -0.012 -15.898  -0.031   0.016
  0.021  -0.027   0.023  12.979  -0.104  -0.031 -16.000   0.138
 -0.008   0.011  -0.012  -0.104  12.802   0.016   0.138 -15.765
 total augmentation occupancy for first ion, spin component:           1
  2.691   0.424  -0.393   0.070  -0.028  -0.155   0.028  -0.011
  0.424   0.118  -0.390   0.072  -0.028  -0.067   0.012  -0.005
 -0.393  -0.390   2.275   0.030  -0.015   0.320   0.029  -0.015
  0.070   0.072   0.030   2.365  -0.108   0.029   0.415  -0.128
 -0.028  -0.028  -0.015  -0.108   2.181  -0.015  -0.128   0.197
 -0.155  -0.067   0.320   0.029  -0.015   0.049   0.008  -0.004
  0.028   0.012   0.029   0.415  -0.128   0.008   0.079  -0.032
 -0.011  -0.005  -0.015  -0.128   0.197  -0.004  -0.032   0.024


------------------------ aborting loop because EDIFF is reached ----------------------------------------


    CHARGE:  cpu time    0.1993: real time    0.1993
    FORLOC:  cpu time    0.3315: real time    0.3315
    FORNL :  cpu time    0.0400: real time    0.0400
    STRESS:  cpu time    0.3673: real time    0.3680
    FORCOR:  cpu time    3.0354: real time    3.0363
    FORHAR:  cpu time    0.9673: real time    0.9675
    MIXING:  cpu time    0.1257: real time    0.1257
    OFIELD:  cpu time    0.0001: real time    0.0001

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Alpha Z     0.89701     0.89701     0.89701
  Ewald     -36.52973   -60.71387   170.56543   -24.69297    24.74526   -60.44698
  Hartree   247.22496   244.24254   375.53437   -19.42250   -10.89331     5.02128
  E(xc)     -72.22903   -72.33214   -71.83960    -0.01593     0.18414    -0.33597
  Local    -399.34434  -378.90236  -708.21069    43.19539    -1.21178    32.61947
  n-local   -61.04957   -61.31459   -60.39544    -0.10453     0.24499    -0.57267
  augment    14.24222    14.66455    12.48551     0.03241    -0.84428     1.50591
  Kinetic   306.01383   313.21218   279.12978     1.19721   -12.61173    23.18496
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  -------------------------------------------------------------------------------------
  Total      -0.77466    -0.24667    -1.83364     0.18908    -0.38672     0.97600
  in kB      -1.24113    -0.39521    -2.93782     0.30294    -0.61959     1.56373
  external pressure =       -1.52 kB  Pullay stress =        0.00 kB

  kinetic pressure (ideal gas correction) =      0.03 kB
  total pressure  =     -1.49 kB
  Total+kin.    -1.214      -0.388      -2.867       0.296      -0.606       1.526

 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
  energy-cutoff  :      800.00
  volume of cell :     1000.00
      direct lattice vectors                 reciprocal lattice vectors
    10.000000000  0.000000000  0.000000000     0.100000000  0.000114090 -0.001410404
    -0.011409000 10.000000000  0.000000000     0.000000000  0.100000000  0.000595569
     0.141108300 -0.059556900 10.000000000     0.000000000  0.000000000  0.100000000

  length of vectors
    10.000000000 10.000006508 10.001172860     0.100010011  0.100001773  0.100000000


 FORCES acting on ions
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction
 -----------------------------------------------------------------------------------------------
   -.136E+02 -.458E+02 0.124E+01   0.168E+02 0.919E+02 -.963E+01   -.334E+01 -.462E+02 0.858E+01   -.286E-04 -.880E-04 -.633E-05
   -.361E+02 -.296E+02 -.180E+02   0.739E+02 0.462E+02 0.405E+02   -.380E+02 -.166E+02 -.226E+02   -.362E-04 0.822E-04 -.129E-05
   -.279E+02 0.622E+02 0.543E+02   0.301E+02 -.676E+02 -.589E+02   -.246E+01 0.558E+01 0.525E+01   0.999E-05 -.474E-04 -.187E-04
   0.319E+02 0.338E+02 -.726E+02   -.346E+02 -.368E+02 0.785E+02   0.306E+01 0.293E+01 -.671E+01   -.243E-04 -.344E-04 0.344E-04
   0.806E+02 -.676E+00 -.330E+02   -.873E+02 0.495E+00 0.357E+02   0.728E+01 0.456E-01 -.323E+01   -.475E-04 0.955E-05 0.184E-04
   -.121E+01 0.314E+02 0.795E+02   0.115E+01 -.342E+02 -.861E+02   -.291E+00 0.301E+01 0.735E+01   -.893E-05 -.143E-04 -.438E-04
 -----------------------------------------------------------------------------------------------
   0.337E+02 0.513E+02 0.114E+02   0.266E-13 0.711E-14 0.000E+00   -.337E+02 -.513E+02 -.114E+02   -.136E-03 -.925E-04 -.173E-04
 
 
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      4.34845      4.20869      5.20021        -0.067101     -0.121912      0.184825
      2.30842      6.27310      1.12975        -0.197981     -0.103058     -0.153232
      4.63631      3.49709      4.58889        -0.296092      0.174308      0.601445
      3.97636      3.80524      6.01777         0.358273     -0.051564     -0.772878
      1.39380      6.24962      1.49232         0.548030     -0.135375     -0.548307
      2.30227      5.89294      0.22162        -0.352073      0.235595      0.698292
 -----------------------------------------------------------------------------------
    total drift:                               -0.006943     -0.002007      0.010144


--------------------------------------------------------------------------------------------------------



  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -28.43874406 eV

  energy  without entropy=      -28.43874406  energy(sigma->0) =      -28.43874406
 
 d Force = 0.5232045E-01[ 0.393E-01, 0.654E-01]  d Energy = 0.5251782E-01-0.197E-03
 d Force =-0.4450050E+01[-0.450E+01,-0.440E+01]  d Ewald  =-0.4449751E+01-0.300E-03


--------------------------------------------------------------------------------------------------------


    POTLOK:  cpu time    3.2829: real time    3.2845


--------------------------------------------------------------------------------------------------------


           RANDOM_SEED =         273936164               36                0
   IONSTEP:  cpu time    0.0002: real time    0.0002

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -28.438744  see above
  kinetic energy EKIN   =         0.055334
  kin. lattice  EKIN_LAT=         0.000000  (temperature   85.62 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -28.383410 eV

  maximum distance moved by ions :      0.22E-02

    WAVPRE:  cpu time    1.0658: real time    1.0721
    FEWALD:  cpu time    0.0003: real time    0.0003
    ORTHCH:  cpu time    0.0179: real time    0.0179
 Prediction of Wavefunctions ALPHA= 1.604 BETA= 0.000
     LOOP+:  cpu time   52.9897: real time   53.0655


----------------------------------------- Iteration    3(   1)  ---------------------------------------


    POTLOK:  cpu time    3.0409: real time    3.0417
    SETDIJ:  cpu time    0.2545: real time    0.2545
     EDDAV:  cpu time    0.0821: real time    0.0821
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1974: real time    0.1974
    MIXING:  cpu time    0.0946: real time    0.0946
    --------------------------------------------
      LOOP:  cpu time    3.6704: real time    3.6713

 eigenvalue-minimisations  :    32
 total energy-change (2. order) :-0.2717900E-01  (-0.7480506E-02)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3904361 magnetization 

 Broyden mixing:
  rms(total) = 0.70174E-02    rms(broyden)= 0.70167E-02
  rms(prec ) = 0.82993E-02
  weight for this iteration     100.00

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        80.52437211
  -Hartree energ DENC   =      -872.56999640
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        68.15559307
  PAW double counting   =      1100.24489929    -1104.46162064
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -215.95823771
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.46592258 eV

  energy without entropy =      -28.46592258  energy(sigma->0) =      -28.46592258


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    3(   2)  ---------------------------------------


    POTLOK:  cpu time    3.0020: real time    3.0029
    SETDIJ:  cpu time    0.2512: real time    0.2512
    EDDIAG:  cpu time    0.0270: real time    0.0270
  RMM-DIIS:  cpu time    0.0517: real time    0.0517
    ORTHCH:  cpu time    0.0099: real time    0.0099
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1965: real time    0.1994
    MIXING:  cpu time    0.0979: real time    0.0980
    --------------------------------------------
      LOOP:  cpu time    3.6373: real time    3.6448

 eigenvalue-minimisations  :    24
 total energy-change (2. order) :-0.7227433E-03  (-0.8177834E-03)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3935597 magnetization 

 Broyden mixing:
  rms(total) = 0.51890E-02    rms(broyden)= 0.51888E-02
  rms(prec ) = 0.62744E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   0.5837
  0.5837

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        80.52437211
  -Hartree energ DENC   =      -872.39621604
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        68.14655627
  PAW double counting   =      1100.06153301    -1104.27306012
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -216.12889824
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.46664533 eV

  energy without entropy =      -28.46664533  energy(sigma->0) =      -28.46664533


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    3(   3)  ---------------------------------------


    POTLOK:  cpu time    3.0309: real time    3.0318
    SETDIJ:  cpu time    0.2533: real time    0.2533
    EDDIAG:  cpu time    0.0264: real time    0.0264
  RMM-DIIS:  cpu time    0.0498: real time    0.0507
    ORTHCH:  cpu time    0.0097: real time    0.0097
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.2004: real time    0.2005
    MIXING:  cpu time    0.0994: real time    0.0995
    --------------------------------------------
      LOOP:  cpu time    3.6709: real time    3.6757

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.3483602E-04  (-0.1796221E-04)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3936067 magnetization 

 Broyden mixing:
  rms(total) = 0.41864E-02    rms(broyden)= 0.41864E-02
  rms(prec ) = 0.50459E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.0125
  1.0125  1.0125

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        80.52437211
  -Hartree energ DENC   =      -872.48234254
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        68.15119877
  PAW double counting   =      1100.21940867    -1104.43269672
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -216.04561846
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.46661049 eV

  energy without entropy =      -28.46661049  energy(sigma->0) =      -28.46661049


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    3(   4)  ---------------------------------------


    POTLOK:  cpu time    3.0148: real time    3.0157
    SETDIJ:  cpu time    0.2518: real time    0.2519
    EDDIAG:  cpu time    0.0264: real time    0.0264
  RMM-DIIS:  cpu time    0.0507: real time    0.0508
    ORTHCH:  cpu time    0.0098: real time    0.0098
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1971: real time    0.1971
    MIXING:  cpu time    0.1019: real time    0.1019
    --------------------------------------------
      LOOP:  cpu time    3.6536: real time    3.6579

 eigenvalue-minimisations  :    24
 total energy-change (2. order) : 0.2768832E-04  (-0.1309207E-04)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3929724 magnetization 

 Broyden mixing:
  rms(total) = 0.13407E-02    rms(broyden)= 0.13406E-02
  rms(prec ) = 0.14893E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   0.9522
  0.7647  1.0459  1.0459

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        80.52437211
  -Hartree energ DENC   =      -872.59244148
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        68.15685414
  PAW double counting   =      1100.43097020    -1104.64662883
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -215.93877663
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.46658280 eV

  energy without entropy =      -28.46658280  energy(sigma->0) =      -28.46658280


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    3(   5)  ---------------------------------------


    POTLOK:  cpu time    3.0138: real time    3.0146
    SETDIJ:  cpu time    0.2514: real time    0.2515
    EDDIAG:  cpu time    0.0265: real time    0.0265
  RMM-DIIS:  cpu time    0.0418: real time    0.0418
    ORTHCH:  cpu time    0.0097: real time    0.0097
       DOS:  cpu time    0.0000: real time    0.0000
    CHARGE:  cpu time    0.1984: real time    0.1984
    MIXING:  cpu time    0.1060: real time    0.1061
    --------------------------------------------
      LOOP:  cpu time    3.6485: real time    3.6528

 eigenvalue-minimisations  :    16
 total energy-change (2. order) :-0.1189557E-05  (-0.1985530E-06)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3930050 magnetization 

 Broyden mixing:
  rms(total) = 0.13023E-02    rms(broyden)= 0.13023E-02
  rms(prec ) = 0.14781E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.4416
  2.4460  1.1621  1.1621  0.9962

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        80.52437211
  -Hartree energ DENC   =      -872.57377807
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        68.15575096
  PAW double counting   =      1100.40986142    -1104.62524039
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -215.95661771
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.46658399 eV

  energy without entropy =      -28.46658399  energy(sigma->0) =      -28.46658399


--------------------------------------------------------------------------------------------------------




----------------------------------------- Iteration    3(   6)  ---------------------------------------


    POTLOK:  cpu time    3.0712: real time    3.0753
    SETDIJ:  cpu time    0.2531: real time    0.2531
    EDDIAG:  cpu time    0.0277: real time    0.0277
  RMM-DIIS:  cpu time    0.0513: real time    0.0513
    ORTHCH:  cpu time    0.0104: real time    0.0104
       DOS:  cpu time    0.0000: real time    0.0000
    --------------------------------------------
      LOOP:  cpu time    3.4146: real time    3.4231

 eigenvalue-minimisations  :    22
 total energy-change (2. order) :-0.6965331E-06  (-0.5215791E-06)
 number of electron      16.0000001 magnetization 
 augmentation part        1.3930050 magnetization 

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         0.89700515
  Ewald energy   TEWEN  =        80.52437211
  -Hartree energ DENC   =      -872.55780043
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        68.15437980
  PAW double counting   =      1100.43491193    -1104.65062337
  entropy T*S    EENTRO =        -0.00000000
  eigenvalues    EBANDS =      -215.97089241
  atomic energy  EATOM  =       914.70206254
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       -28.46658469 eV

  energy without entropy =      -28.46658469  energy(sigma->0) =      -28.46658469


--------------------------------------------------------------------------------------------------------




 average (electrostatic) potential at core
  the test charge radii are     0.7215  0.5201
  (the norm of the test charge is              1.0000)
       1 -80.3126       2 -80.3085       3 -44.4789       4 -44.5270       5 -44.4828
       6 -44.5072
 
 
 
 E-fermi :  -6.6864     XC(G=0):  -1.0129     alpha+bet : -0.2602


 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -25.0631      2.00000
      2     -25.0613      2.00000
      3     -13.1093      2.00000
      4     -13.0873      2.00000
      5      -8.8621      2.00000
      6      -8.8504      2.00000
      7      -6.9854      2.00000
      8      -6.9664      2.00000
      9      -1.1216      0.00000
     10      -0.3509      0.00000
     11       0.5109      0.00000
     12       0.8393      0.00000


--------------------------------------------------------------------------------------------------------


 soft charge-density along one line, spin component           1
         0         1         2         3         4         5         6         7         8         9
 total charge-density along one line
 
 pseudopotential strength for first ion, spin component:           1
 13.771 -16.880   0.102  -0.019   0.008  -0.124   0.023  -0.009
-16.880  20.722  -0.129   0.024  -0.010   0.159  -0.030   0.012
  0.102  -0.129 -10.426  -0.017   0.009  12.907   0.023  -0.012
 -0.019   0.024  -0.017 -10.485   0.082   0.023  12.985  -0.108
  0.008  -0.010   0.009   0.082 -10.346  -0.012  -0.108  12.801
 -0.124   0.159  12.907   0.023  -0.012 -15.903  -0.030   0.016
  0.023  -0.030   0.023  12.985  -0.108  -0.030 -16.006   0.143
 -0.009   0.012  -0.012  -0.108  12.801   0.016   0.143 -15.762
 total augmentation occupancy for first ion, spin component:           1
  2.734   0.446  -0.418   0.079  -0.032  -0.165   0.031  -0.013
  0.446   0.125  -0.402   0.075  -0.030  -0.071   0.013  -0.005
 -0.418  -0.402   2.310   0.028  -0.014   0.331   0.028  -0.015
  0.079   0.075   0.028   2.405  -0.122   0.028   0.427  -0.132
 -0.032  -0.030  -0.014  -0.122   2.196  -0.015  -0.132   0.201
 -0.165  -0.071   0.331   0.028  -0.015   0.052   0.008  -0.004
  0.031   0.013   0.028   0.427  -0.132   0.008   0.082  -0.033
 -0.013  -0.005  -0.015  -0.132   0.201  -0.004  -0.033   0.025


------------------------ aborting loop because EDIFF is reached ----------------------------------------


    CHARGE:  cpu time    0.2018: real time    0.2020
    FORLOC:  cpu time    0.3341: real time    0.3341
    FORNL :  cpu time    0.0408: real time    0.0416
    STRESS:  cpu time    0.3644: real time    0.3644
    FORCOR:  cpu time    3.0297: real time    3.0304
    FORHAR:  cpu time    0.9527: real time    0.9529
    MIXING:  cpu time    0.1063: real time    0.1064
    OFIELD:  cpu time    0.0001: real time    0.0001

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Alpha Z     0.89701     0.89701     0.89701
  Ewald     -34.47019   -58.09440   173.08881   -23.96382    24.89538   -60.52901
  Hartree   249.11427   246.43744   377.10392   -19.27383   -11.08852     5.46685
  E(xc)     -72.52558   -72.62230   -72.14548    -0.01196     0.18118    -0.32734
  Local    -402.84565  -383.37806  -711.11954    42.45177    -0.89071    31.67941
  n-local   -61.76435   -62.01361   -61.24250    -0.14616     0.18748    -0.46849
  augment    14.32066    14.74312    12.54922     0.03001    -0.84978     1.51725
  Kinetic   307.50914   314.41944   281.09234     1.05867   -12.45048    22.76107
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  -------------------------------------------------------------------------------------
  Total       0.23532     0.38865     0.22378     0.14467    -0.01545     0.09974
  in kB       0.37702     0.62269     0.35853     0.23179    -0.02476     0.15980
  external pressure =        0.45 kB  Pullay stress =        0.00 kB

  kinetic pressure (ideal gas correction) =      0.09 kB
  total pressure  =      0.54 kB
  Total+kin.     0.447       0.638       0.543       0.211       0.013       0.058

 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
  energy-cutoff  :      800.00
  volume of cell :     1000.00
      direct lattice vectors                 reciprocal lattice vectors
    10.000000000  0.000000000  0.000000000     0.100000000  0.000114090 -0.001410404
    -0.011409000 10.000000000  0.000000000     0.000000000  0.100000000  0.000595569
     0.141108300 -0.059556900 10.000000000     0.000000000  0.000000000  0.100000000

  length of vectors
    10.000000000 10.000006508 10.001172860     0.100010011  0.100001773  0.100000000


 FORCES acting on ions
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction
 -----------------------------------------------------------------------------------------------
   -.140E+02 -.490E+02 0.214E+01   0.175E+02 0.968E+02 -.111E+02   -.348E+01 -.471E+02 0.885E+01   -.659E-03 -.106E-01 0.178E-02
   -.388E+02 -.308E+02 -.196E+02   0.781E+02 0.480E+02 0.430E+02   -.388E+02 -.170E+02 -.231E+02   -.906E-02 -.382E-02 -.531E-02
   -.283E+02 0.640E+02 0.551E+02   0.308E+02 -.702E+02 -.606E+02   -.258E+01 0.598E+01 0.550E+01   0.819E-03 -.291E-02 -.167E-02
   0.327E+02 0.354E+02 -.744E+02   -.361E+02 -.390E+02 0.817E+02   0.328E+01 0.325E+01 -.719E+01   -.116E-02 -.199E-02 0.254E-02
   0.829E+02 -.388E+00 -.334E+02   -.910E+02 0.189E+00 0.367E+02   0.783E+01 0.876E-01 -.339E+01   -.350E-02 -.334E-03 0.735E-03
   -.647E+00 0.324E+02 0.817E+02   0.619E+00 -.358E+02 -.897E+02   -.221E+00 0.325E+01 0.789E+01   -.692E-03 -.145E-02 -.318E-02
 -----------------------------------------------------------------------------------------------
   0.339E+02 0.515E+02 0.115E+02   -.178E-14 -.711E-14 0.000E+00   -.339E+02 -.515E+02 -.114E+02   -.143E-01 -.211E-01 -.511E-02
 
 
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      4.34833      4.20827      5.20054         0.072413      0.667359     -0.134111
      2.30794      6.27287      1.12940         0.521522      0.226802      0.305916
      4.62915      3.50315      4.60357        -0.036438     -0.306247      0.046743
      3.98547      3.80589      5.99808        -0.041300     -0.361604      0.101786
      1.40887      6.24687      1.47935        -0.275308     -0.111298     -0.113423
      2.29491      5.89940      0.23997        -0.249497     -0.118610     -0.197687
 -----------------------------------------------------------------------------------
    total drift:                               -0.008608     -0.003597      0.009222


--------------------------------------------------------------------------------------------------------



  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =       -28.46658469 eV

  energy  without entropy=      -28.46658469  energy(sigma->0) =      -28.46658469
 
 d Force = 0.2711380E-01[-0.919E-02, 0.634E-01]  d Energy = 0.2784063E-01-0.727E-03
 d Force =-0.7203530E+01[-0.734E+01,-0.707E+01]  d Ewald  =-0.7202390E+01-0.114E-02


--------------------------------------------------------------------------------------------------------


    POTLOK:  cpu time    3.3558: real time    3.3575


--------------------------------------------------------------------------------------------------------


           RANDOM_SEED =         273936164               36                0
   IONSTEP:  cpu time    0.0002: real time    0.0002

  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)
  ---------------------------------------------------
% ion-electron   TOTEN  =       -28.466585  see above
  kinetic energy EKIN   =         0.080047
  kin. lattice  EKIN_LAT=         0.000000  (temperature  123.85 K)
  nose potential ES     =         0.000000
  nose kinetic   EPS    =         0.000000
  ---------------------------------------------------
  total energy   ETOTAL =       -28.386537 eV

  maximum distance moved by ions :      0.21E-02


 mean value of Nose-termostat <S>:     1.000 mean value of <T> :    74.036
 mean temperature <T/S>/<1/S>  :    74.036

    WAVPRE:  cpu time    1.0674: real time    1.0751
    FEWALD:  cpu time    0.0003: real time    0.0003
    ORTHCH:  cpu time    0.0183: real time    0.0183
 Prediction of Wavefunctions ALPHA= 3.790 BETA=-4.569
    POTLOK:  cpu time    3.1859: real time    3.1867
    EDDIAG:  cpu time    0.0277: real time    0.0277
     LOOP+:  cpu time   34.6434: real time   34.6905
    4ORBIT:  cpu time    0.0000: real time    0.0000

 total amount of memory used by VASP MPI-rank0   739137. kBytes
=======================================================================

   base      :      30000. kBytes
   nonlr-proj:       1360. kBytes
   fftplans  :     243296. kBytes
   grid      :     434847. kBytes
   one-center:         18. kBytes
   wavefun   :      29616. kBytes
 
  
  
 General timing and accounting informations for this job:
 ========================================================
  
                  Total CPU time used (sec):      162.932
                            User time (sec):      155.084
                          System time (sec):        7.848
                         Elapsed time (sec):      163.263
  
                   Maximum memory used (kb):     1660936.
                   Average memory used (kb):           0.
  
                          Minor page faults:       223955
                          Major page faults:            0
                 Voluntary context switches:         1984
