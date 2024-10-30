INPUT_PARAMETERS
suffix                          spin
calculation                     relax
ecutwfc                         100
scf_thr                         1.0e-3
scf_nmax                        200
smearing_method                 gauss
smearing_sigma                  0.01
#ks_solver                       genelpa
basis_type                      pw
symmetry                        0
noncolin			1
nspin                           4

onsite_radius 3

relax_nmax 3

force_thr 0.00001

sc_mag_switch                   1
decay_grad_switch               1
sc_thr                          1e-4
nsc                             100
nsc_min                         2
alpha_trial                     0.01
sccut                           3



#md_nstep 3
#md_type                nve
#md_dt                  1
#md_tfirst              300

