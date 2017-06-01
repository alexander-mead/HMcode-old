# HMcode

This code is produces the non-linear matter power spectrum using the halomodel approach described in Mead et al. (2015; http://arxiv.org/abs/1505.07833). Appendix B of that paper details the methods for doing the calculation. It also now includes some small updates from Mead et al. (2016; http://arxiv.org/abs/1602.02154).

If you use this work, or this code, I would be very grateful if you were to cite the original paper, and the updated paper if you use results from that. The code itself can also be cited: http://ascl.net/1508.001

The code should compile with any fortran compiler, and does not need to be pointed to any libraries. I use 'ifort' and compile with '>ifort HMcode.f90' but it also works with gfortran.

When it starts, the code fills up some arrays for the wavenumber 'k' and redshift 'z' values for which the user wants the power spectrum. The code then calls 'assign_cosmology', which sets the cosmological parameters - if these need to be changed then it should be done in this subroutine. The code uses the Eistenstein + Hu (1998) approximation for the linear power spectrum, but anything (e.g., a CAMB linear spectrum) could be wired in quite easily. The code then calls 'initialise_cosmology', which normalises the power spectrum to the correct sigma_8 and fills up look-up tables that contain the linear growth function and sigma(R) for later use.

The option 'imead' can be changed in the code. If imead=0 the code performs the standard halo model calculation (Dv=200, dc=1.686, Sheth & Tormen (1999) mass function, Bullock (2001) concentration-mass relation; although it neglects the standard bias term in the two-halo term, because this is not important for the matter spectrum), whereas imead=1 performs the accurate calculation detailed in Appendix B of http://arxiv.org/abs/1505.07833 

The code then loops through 'z' (outer loop) and 'k' (inner loop) producing power at each redshift and wave number. The ordering of loops (z then k) is important because for each new redshift the code needs to call 'halomod_init' to fill up redshift-dependent look-up tables for the halomodel calculation. These look-up tables contain various redshift-dependent halo properties, such as mass, radius, nu, concentration etc. which are used in the one-halo integral.

Once these look-up tables have been filled the halomodel integral can then be carried out. This calculation is done by the routine 'halomod', which calls 'p_1h' and 'p_2h' to evaluate 1- and 2-halo terms and then uses these to compute the full power spectrum. The power spectrum at each k and z is then added to an array which is printed out to power.dat (k, pow(z1), pow(z2), ...) when the code finishes.

In testing I was able to get 16 redshifts, with 200 k-points, in 0.4 seconds (using ifort with -O3 optimisation). 

The file 'plot.p' is a gnuplot script to plot the output. It can be run using "gnuplot> load 'plot.p'".

Please let me know if you need any help running the code. Or if you have any questions whatsoever.

HMcode is now also integrated within CAMB (https://github.com/cmbant/CAMB).

Alexander Mead
(mead@phas.ubc.ca)

UPDATE July 7th, 2015
One user reported crashes that occured for quite extreme cosmological models (n_s < 0.5, sig8 < 0.3 z>5). I have fixed this rather crudely by adding IF statements that catch problems (which manifest themsevles as extreme parameter values). The physical reason for these problems is that models with these odd cosmological parameters have R_nl<<1 Mpc and so have very few haloes. Some of the routines I was using previously had assumed that R_nl would not be so tiny.

UPDATE January 23rd, 2016
Updated the code a little so that it no longer calculates a range in nu and simply converts a mass range into a nu range to do the integral. The mass range is fixed from haloes of 1 Solar mass to 10^16 Solar masses, it is difficult to imagine an application of this code where this halo mass range would not be sufficient. This further helps when considering strange cosmological models at high redshift that suffer from a lack of haloes, for these models doing a nu to M inversion occasionally reached incredibly tiny halo masses that contribute negligbly to the power spectrum on cosmological scales due to their tiny sizes.

UPDATE February 4th, 2016
Included updates from Mead et al. (2016) including parameterising the two-halo damping term in terms of f(sigma_v) and slightly recalibrated values of alpha and f_damp. Now works for w(a)CDM models, where w(a)=w_0+(1.-a)*w_a.

UPDATE August 2nd, 2016
Small update to the README and some very minor fixes in the code to facilitate direct comparisons with other halomodel power codes.

######

Adding in a CAMB linear P(k)

Given the differences between CAMB and Eisenstein + Hu (1998) one might wish to make HMcode read in a linear CAMB power spectrum and work with that instead (or any other tabulated power spectrum). This is fine, and is what I did when comparing to Cosmic Emu in the paper (where the difference between CAMB and Eisenstein + Hu *was* important) but there is a subtle thing to bear in mind:

The halo-model calculation does integrals that are technically defined over the entire k range from 0 to infinity. In practice contributions from both very small and very large scales are suppressed, but the integration routines still need to know about the power at these scales sometimes, otherwise they may go bonkers. Obviously this is a problem given that one usually has a tabulated linear power spectrum defined on some finite range in k. The way I dealt with this was to read in the linear spectrum but then extrapolate if the code called 'p_lin' for k values below the minimum, or above the maximum, using physically motivated extrapolations. For example you know that the linear power spectrum is a power-law down to k->0 (\Delta^2 \propto k^{3+n}) and the high-k part can be approximated as \Delta^2 \propto k^{n-1}log(k)^2. 


