MODULE cosdef

  TYPE cosmology
     !Contains only things that do not need to be recalculated with each new z
     REAL :: om_m, om_b, om_v, om_c, h, n, sig8, w, wa
     REAL :: A
     REAL, ALLOCATABLE :: r_sigma(:), sigma1d(:)
     REAL, ALLOCATABLE :: growth(:), a_growth(:)
     REAL, ALLOCATABLE :: ktab(:), tktab(:), pktab(:)
  END TYPE cosmology

  TYPE tables
     !Stuff that needs to be recalculated for each new z
     REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
     REAL :: sigv, sigv100, c3, knl, rnl, neff, sig8z
     INTEGER :: n
  END TYPE tables

END MODULE cosdef

PROGRAM HMcode
  USE cosdef
  IMPLICIT NONE
  REAL :: z
  REAL :: p1h, p2h, pfull, plin
  REAL, ALLOCATABLE :: k(:), ztab(:), ptab(:,:)
  INTEGER :: i, j, nk, nz
  INTEGER :: ihm, imead
  REAL :: kmin, kmax, zmin, zmax
  REAL, PARAMETER :: pi=3.141592654
  TYPE(cosmology) :: cosi
  TYPE(tables) :: lut
  LOGICAL :: lexist
  CHARACTER(len=64) :: input, output

  !HMcode developed by Alexander Mead
  !If you use this in your work please cite the original paper: http://arxiv.org/abs/1505.07833
  !and consider citing the source code at ASCL: http://ascl.net/1508.001

  !ihm
  !0 - Non-verbose
  !1 - Verbose
  ihm=1

  !imead
  !0 - Do the standard halo model calculation (Dv=200, dc=1.686, Sheth & Tormen (199) mass function, Bullock (2001) c(M)'
  !1 - Do the accurate calculation detailed in 1505.07833 with updates in Mead et al. (2016)
  imead=1

  WRITE(*,*)
  WRITE(*,*) 'Welcome to HMcode'
  WRITE(*,*) '================='
  WRITE(*,*)
  IF(imead==0) WRITE(*,*) 'Doing standard calculation'
  IF(imead==1) WRITE(*,*) 'Doing accurate calculation'
  WRITE(*,*)

  !Set number of k points and k range (log spaced)
  nk=200
  kmin=0.001
  kmax=1.e4
  CALL fill_table(kmin,kmax,k,nk,1)

  WRITE(*,*) 'k min:', kmin
  WRITE(*,*) 'k max:', kmax
  WRITE(*,*) 'number of k:', nk
  WRITE(*,*)

  !Set the number of redshifts and range (linearly spaced)
  nz=16
  zmin=0.
  zmax=4.
  CALL fill_table(zmin,zmax,ztab,nz,0)

  WRITE(*,*) 'z min:', zmin
  WRITE(*,*) 'z max:', zmax
  WRITE(*,*) 'number of z:', nz
  WRITE(*,*)
  
  !Fill table for output power
  ALLOCATE(ptab(nz,nk))

  !Assigns the cosmological model
  CALL assign_cosmology(cosi)

  !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
  CALL initialise_cosmology(cosi)

  !Ignore this, only useful for bug tests
  !DO
  !CALL random_cosmology(cosi)
  
  CALL write_cosmology(cosi)

  !Loop over redshifts
  DO j=1,nz

     !Sets the current redshift from the table
     z=ztab(j)

     !Initiliasation for the halomodel calcualtion
     CALL halomod_init(z,lut,cosi)

     !Loop over k values
     DO i=1,SIZE(k)

        plin=p_lin(k(i),z,cosi)

        CALL halomod(k(i),z,p1h,p2h,pfull,plin,lut,cosi)

        ptab(j,i)=pfull

     END DO

     IF(j==1) THEN
        WRITE(*,fmt='(A5,A7)') 'i', 'z'
        WRITE(*,fmt='(A13)') '   ============'
     END IF
     WRITE(*,fmt='(I5,F8.3)') j, ztab(j)

  END DO
  WRITE(*,*)

  output='power.dat'
  WRITE(*,fmt='(A19,A10)') 'Writing output to:', TRIM(output)
  WRITE(*,*)
  WRITE(*,*) 'The top row of the file contains the redshifts (the first entry is hashes - #####)'
  WRITE(*,*) 'Subsequent rows contain ''k'' and then the halo-model power for each redshift'
  OPEN(7,file=output)
  DO i=0,nk
     IF(i==0) THEN
        WRITE(7,fmt='(A20,40F20.10)') '#####', (ztab(j), j=1,nz)
     ELSE
        WRITE(7,fmt='(F20.10,40F20.10)') k(i), (ptab(j,i), j=1,nz)
     END IF
  END DO
  CLOSE(7)
  WRITE(*,*) 'Done'
  WRITE(*,*)

  !Ignore this, only useful for bug tests
  !END DO

CONTAINS

  FUNCTION Delta_v(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: Delta_v
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Virialised overdensity
    IF(imead==0) THEN
       Delta_v=200.
    ELSE IF(imead==1) THEN
       Delta_v=418.*(omega_m(z,cosm)**(-0.352))       
    END IF
    
  END FUNCTION Delta_v

  FUNCTION delta_c(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: delta_c
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Linear collapse density
    IF(imead==0) THEN
       delta_c=1.686
    ELSE IF(imead==1) THEN
       delta_c=1.59+0.0314*log(sigma_cb(8.,z,cosm))
    END IF

    !Nakamura & Suto (1997) fitting formula for LCDM
    delta_c=delta_c*(1.+0.0123*log10(omega_m(z,cosm)))

  END FUNCTION delta_c

  FUNCTION eta(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: eta
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(imead==0) THEN
       eta=0.
    ELSE IF(imead==1) THEN
       !The first parameter here is 'eta_0' in Mead et al. (2015; arXiv 1505.07833)
       eta=0.603-0.3*(sigma_cb(8.,z,cosm))
    END IF

  END FUNCTION eta
 
  FUNCTION kstar(lut,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: kstar
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut

    IF(imead==0) THEN
       !Set to zero for the standard Poisson one-halo term
       kstar=0.
    ELSE IF(imead==1) THEN
       !One-halo cut-off wavenumber
       kstar=0.584*(lut%sigv)**(-1.)
    END IF

  END FUNCTION kstar

  FUNCTION As(cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: As
    TYPE(cosmology), INTENT(IN) :: cosm

    !Halo concentration pre-factor
    IF(imead==0) THEN
       !Set to 4 for the standard Bullock value
       As=4.
    ELSE IF(imead==1) THEN
    !This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833)
       As=3.13
    END IF

  END FUNCTION As

  FUNCTION fdamp(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL ::fdamp
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Linear theory damping factor
    IF(imead==0) THEN
       !Set to 0 for the standard linear theory two halo term
       fdamp=0.
    ELSE IF(imead==1) THEN
       !fdamp=0.188*sigma_cb(8.,z,cosm)**4.29
       fdamp=0.0095*lut%sigv100**1.37
    END IF

    !Catches extreme values of fdamp
    IF(fdamp<1.e-3) fdamp=1.e-3
    IF(fdamp>0.99)  fdamp=0.99

  END FUNCTION fdamp

  FUNCTION alpha(lut,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: alpha
    TYPE(tables), INTENT(IN) :: lut
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(imead==0) THEN
       !Set to 1 for the standard halo model addition of one- and two-halo terms
       alpha=1.
    ELSE IF(imead==1) THEN
       !This uses the top-hat defined neff
       !alpha=2.93*(1.77**lut%neff)
       alpha=3.24*1.85**lut%neff
    END IF

    !Catches values of alpha that are crazy
    IF(alpha>2.)  alpha=2.
    IF(alpha<0.5) alpha=0.5

  END FUNCTION alpha

  SUBROUTINE write_parameters(z,lut,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    !TYPE(parameters), INTENT(IN) :: p
    TYPE(tables), INTENT(IN) :: lut

    !This subroutine writes out the physical parameters at some redshift 
    !(e.g. Delta_v) rather than the model parameters

    WRITE(*,*) 'Parameters at your redshift'
    WRITE(*,*) '==========================='
    WRITE(*,fmt='(A10,F10.5)') 'z:', z
    WRITE(*,fmt='(A10,F10.5)') 'Dv:', Delta_v(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'dc:', delta_c(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'eta:', eta(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'k*:', kstar(lut,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'A:', As(cosm)
    WRITE(*,fmt='(A10,F10.5)') 'fdamp:', fdamp(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'alpha:', alpha(lut,cosm)
    WRITE(*,*)

  END SUBROUTINE write_parameters

  FUNCTION r_nl(lut)

    USE cosdef
    TYPE(tables), INTENT(IN) :: lut
    REAL :: r_nl

    !Calculates k_nl as 1/R where nu(R)=1.

    IF(lut%nu(1)>1.) THEN
       !This catches some very strange values
       r_nl=lut%rr(1)
    ELSE
       r_nl=exp(find(log(1.),log(lut%nu),log(lut%rr),3,3))
    END IF

  END FUNCTION r_nl

  SUBROUTINE halomod(k,z,p1h,p2h,pfull,plin,lut,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL, INTENT(OUT) :: p1h, p2h, pfull
    REAL, INTENT(IN) :: plin, k, z
    REAL :: alp
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut

    !Calls expressions for one- and two-halo terms and then combines
    !to form the full power spectrum
    IF(k==0.) THEN
       p1h=0.
       p2h=0.
    ELSE
       p1h=p_1h(k,z,lut,cosm)
       p2h=p_2h(k,z,plin,lut,cosm)
    END IF

    alp=alpha(lut,cosm)
    pfull=(p2h**alp+p1h**alp)**(1./alp)

  END SUBROUTINE halomod

  SUBROUTINE fill_table(min,max,arr,n,ilog)

    IMPLICIT NONE
    INTEGER :: i
    REAL, INTENT(IN) :: min, max
    REAL :: a, b
    REAL, ALLOCATABLE :: arr(:)
    INTEGER, INTENT(IN) :: ilog, n

    !Fills array 'arr' in equally spaced intervals
    !ilog=0 does linear spacing
    !ilog=1 does log spacing

    IF(ALLOCATED(arr)) DEALLOCATE(arr)

    ALLOCATE(arr(n))

    arr=0.

    IF(ilog==0) THEN
       a=min
       b=max
    ELSE IF(ilog==1) THEN
       a=log(min)
       b=log(max)
    END IF

    IF(n==1) THEN
       arr(1)=a
    ELSE IF(n>1) THEN
       DO i=1,n
          arr(i)=a+(b-a)*float(i-1)/float(n-1)
       END DO
    END IF

    IF(ilog==1) arr=exp(arr)

  END SUBROUTINE fill_table

  SUBROUTINE write_cosmology(cosm)

    USE cosdef
    IMPLICIT NONE
    TYPE(cosmology) :: cosm

    IF(ihm==1) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_m:', cosm%om_m
    IF(ihm==1) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_b:', cosm%om_b
    IF(ihm==1) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_c:', cosm%om_c
    IF(ihm==1) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_v:', cosm%om_v
    IF(ihm==1) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'h:', cosm%h
    IF(ihm==1) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'w_0:', cosm%w
    IF(ihm==1) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'w_a:', cosm%wa
    IF(ihm==1) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'sig8:', cosm%sig8
    IF(ihm==1) WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'n:', cosm%n   
    IF(ihm==1) WRITE(*,*)

  END SUBROUTINE write_cosmology

  SUBROUTINE assign_cosmology(cosm)

    USE cosdef
    IMPLICIT NONE
    CHARACTER(len=64) :: input
    TYPE(cosmology) :: cosm
    LOGICAL :: lexist

    cosm%om_m=0.3
    cosm%om_v=1.-cosm%om_m
    cosm%om_b=0.05
    cosm%om_c=cosm%om_m-cosm%om_b
    cosm%h=0.7
    cosm%w=-1.
    cosm%sig8=0.8
    cosm%n=0.96
    cosm%wa=0.

  END SUBROUTINE assign_cosmology

  SUBROUTINE initialise_cosmology(cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: sigi
    TYPE(cosmology) :: cosm

    CALL fill_growtab(cosm)

    cosm%A=1.

    sigi=sigma(8.,0.,cosm)

    IF(ihm==1) WRITE(*,*) 'INITIALISE: Initial sigma_8:', sigi

    cosm%A=cosm%sig8/sigi    

    sigi=sigma(8.,0.,cosm)

    IF(ihm==1) THEN
       WRITE(*,*) 'INITIALISE: Normalisation factor:', cosm%A
       WRITE(*,*) 'INITIALISE: Target sigma_8:', cosm%sig8
       WRITE(*,*) 'INITIALISE: Final sigma_8 (calculated):', sigi
       WRITE(*,*) 'INITIALISE: Complete'
       WRITE(*,*)
    END IF

    !Fill tables of r vs. sigma(r)
    CALL fill_sigtab(cosm)

  END SUBROUTINE initialise_cosmology

  SUBROUTINE random_cosmology(cosm)

    USE cosdef
    TYPE(cosmology) :: cosm
    REAL :: om_m_min, om_m_max, om_b_min, om_b_max, n_min, n_max
    REAL :: w_min, w_max, h_min, h_max, sig8_min, sig8_max, wa_min, wa_max

    !Needs to be set to normalise P_lin
    cosm%A=1.

    om_m_min=0.1
    om_m_max=1.
    cosm%om_m=ran(om_m_min,om_m_max)

    cosm%om_v=1.-cosm%om_m

    om_b_min=0.005
    om_b_max=MIN(0.095,cosm%om_m)
    cosm%om_b=ran(om_b_min,om_b_max)

    cosm%om_c=cosm%om_m-cosm%om_b

    n_min=0.5
    n_max=1.5
    cosm%n=ran(n_min,n_max)

    h_min=0.4
    h_max=1.2
    cosm%h=ran(h_min,h_max)

    w_min=-1.5
    w_max=-0.5
    cosm%w=ran(w_min,w_max)

    wa_min=-1.
    wa_max=-cosm%w*0.8
    cosm%wa=ran(wa_min,wa_max)

    sig8_min=0.2
    sig8_max=1.5
    cosm%sig8=ran(sig8_min,sig8_max)

  END SUBROUTINE random_cosmology

  SUBROUTINE RNG_set

    IMPLICIT NONE
    INTEGER :: int, timearray(3)
    REAL :: rand

    WRITE(*,*) 'Initialising RNG'
    !This fills the time array using the system clock!
    !If called within the same second the numbers will be identical!
    CALL itime(timeArray)
    !This then initialises the generator!
    int=rand(timeArray(1)+timeArray(2)+timeArray(3))
    WRITE(*,*) 'RNG set'
    WRITE(*,*)

  END SUBROUTINE RNG_set

  FUNCTION ran(x1,x2)

    IMPLICIT NONE
    REAL :: rand, ran
    REAL :: x1,x2

    !Generates a random number in the interval x1->x2 with uniform probability
    !rand is some inbuilt function!
    ran=x1+(x2-x1)*(rand(0))

  END FUNCTION ran

  SUBROUTINE allocate_LUT(lut)

    USE cosdef
    TYPE(tables) :: lut
    INTEGER :: n

    !Allocates memory for the look-up tables
    n=lut%n

    ALLOCATE(lut%zc(n),lut%m(n),lut%c(n),lut%rv(n))
    ALLOCATE(lut%nu(n),lut%rr(n),lut%sigf(n),lut%sig(n))

    lut%zc=0.
    lut%m=0.
    lut%c=0.
    lut%rv=0.
    lut%nu=0.
    lut%rr=0.
    lut%sigf=0.
    lut%sig=0.

  END SUBROUTINE allocate_LUT

  SUBROUTINE deallocate_LUT(lut)

    USE cosdef
    TYPE(tables) :: lut

    !Deallocates look-up tables
    DEALLOCATE(lut%zc,lut%m,lut%c,lut%rv,lut%nu,lut%rr,lut%sigf,lut%sig)

  END SUBROUTINE deallocate_LUT

  SUBROUTINE halomod_init(z,lut,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    INTEGER :: i, imin, imax, n
    REAL :: rin, dr, Dv, dc, f, m, mmin, mmax, nu, r, sig
    !REAL, ALLOCATABLE :: rg_up(:), rg_dn(:), nu_up(:), nu_dn(:), mg_up(:), mg_dn(:)
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut

    !Halo-model initialisation routine
    !The computes other tables necessary for the one-halo integral

    !Find value of sigma_v
    lut%sigv=sqrt(dispint(z,cosm))
    lut%sigv100=sigma_v(100.,z,cosm)
    lut%sig8z=sigma(8.,z,cosm)

    !n=500
    !ALLOCATE(rg_up(n),rg_dn(n),nu_up(n),nu_dn(n),mg_up(n),mg_dn(n))

    IF(ihm==1) WRITE(*,*) 'HALOMOD: Filling look-up tables'
    IF(ihm==1) WRITE(*,*) 'HALOMOD: Tables being filled at redshift:', z
      
    IF(ihm==1) WRITE(*,*) 'HALOMOD: sigv [Mpc/h]:', lut%sigv
    IF(ihm==1) WRITE(*,*) 'HALOMOD: sigv100 [Mpc/h]:', lut%sigv100
    IF(ihm==1) WRITE(*,*) 'HALOMOD: sig8(z):', lut%sig8z

    IF(ALLOCATED(lut%rr)) CALL deallocate_LUT(lut)

    n=256

    lut%n=n
    CALL allocate_lut(lut)

    !Mass range for halo model calculation
    mmin=1.e0
    mmax=1.e16

    dc=delta_c(z,cosm)

    DO i=1,n 

       m=exp(log(mmin)+log(mmax/mmin)*float(i-1)/float(n-1))
       r=radius_m(m,cosm)
       sig=sigma_cb(r,z,cosm)
       nu=dc/sig

       lut%m(i)=m
       lut%rr(i)=r
       lut%sig(i)=sig
       lut%nu(i)=nu

    END DO

    IF(ihm==1) WRITE(*,*) 'HALOMOD: m, r, nu, sig tables filled'

    !Fills up a table for sigma(fM) for Bullock c(m) relation
    !This is the f=0.01 parameter in the Bullock realtion sigma(fM,z)
    f=0.01**(1./3.)
    DO i=1,lut%n
       lut%sigf(i)=sigma_cb(lut%rr(i)*f,z,cosm)
    END DO
    IF(ihm==1) WRITE(*,*) 'HALOMOD: sigf tables filled'  

    !Fill virial radius table using real radius table
    Dv=Delta_v(z,cosm)
    lut%rv=lut%rr/(Dv**(1./3.))

    IF(ihm==1) WRITE(*,*) 'HALOMOD: rv tables filled'  
    IF(ihm==1) WRITE(*,*) 'HALOMOD: nu min:', lut%nu(1)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: nu max:', lut%nu(lut%n)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: R_v min [Mpc/h]:', lut%rv(1)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: R_v max [Mpc/h]:', lut%rv(lut%n)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: M min [Msun/h]:', lut%m(1)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: M max [Msun/h]:', lut%m(lut%n)

    !Find non-linear radius and scale
    lut%rnl=r_nl(lut)
    lut%knl=1./lut%rnl

    IF(ihm==1) WRITE(*,*) 'HALOMOD: r_nl [Mpc/h]:', lut%rnl
    IF(ihm==1) WRITE(*,*) 'HALOMOD: k_nl [h/Mpc]:', lut%knl

    lut%neff=neff(lut,cosm)

    IF(ihm==1) WRITE(*,*) 'HALOMOD: n_eff:', lut%neff

    CALL conc_bull(z,cosm,lut)

    IF(ihm==1) WRITE(*,*) 'HALOMOD: c tables filled'
    IF(ihm==1) WRITE(*,*) 'HALOMOD: c min [Msun/h]:', lut%c(lut%n)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: c max [Msun/h]:', lut%c(1)
    IF(ihm==1) WRITE(*,*) 'HALOMOD: Done'

    IF(ihm==1) WRITE(*,*)

    if(ihm==1) CALL write_parameters(z,lut,cosi)

    ihm=0

  END SUBROUTINE halomod_init

  PURE FUNCTION radius_m(m,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: radius_m
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, PARAMETER :: pi=3.141592654

    radius_m=(3.*m/(4.*pi*cosmic_density(cosm)))**(1./3.)

  END FUNCTION radius_m

  FUNCTION neff(lut,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: neff
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut

    !Numerical differentiation to find effective index at collapse
    neff=-3.-derivative_table(log(lut%rnl),log(lut%rr),log(lut%sig**2.),3,3)

    !For some bizarre cosmologies r_nl is very small, so almost no collapse has occured
    !In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
    IF(neff<cosm%n-4.) neff=cosm%n-4.
    IF(neff>cosm%n)    neff=cosm%n

  END FUNCTION neff

  SUBROUTINE conc_bull(z,cosm,lut)

    USE cosdef
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology) :: cosm, cos_lcdm
    TYPE(tables) :: lut
    REAL :: A, zinf, ainf, zf, g_lcdm, g_wcdm, w
    INTEGER :: i

    !Calculates the Bullock et al. (2001) mass-concentration relation

    A=As(cosm)

    !Fill the collapse z look-up table
    CALL zcoll_bull(z,cosm,lut)

    !Fill the concentration look-up table
    DO i=1,lut%n

       zf=lut%zc(i)
       lut%c(i)=A*(1.+zf)/(1.+z)

       !Dolag2004 prescription for adding DE dependence
       !IF((cosm%w .NE. -1.) .OR. (cosm%wa .NE. 0)) THEN

       zinf=10.

       g_wcdm=grow(zinf,cosm)

       !Make a LCDM cosmology
       cos_lcdm=cosm
       DEALLOCATE(cos_lcdm%growth)
       DEALLOCATE(cos_lcdm%a_growth)
       cos_lcdm%w=-1.
       cos_lcdm%wa=0.

       ainf=1./(1.+zinf)

       !Needs to use grow_int explicitly in case tabulated values are stored
       g_lcdm=grow_int(ainf,0.001,cos_lcdm)

       !Changed this to a power of 1.5, which produces more accurate results for extreme DE
       lut%c(i)=lut%c(i)*((g_wcdm/g_lcdm)**1.5)

       !END IF

    END DO

  END SUBROUTINE conc_bull

  SUBROUTINE zcoll_bull(z,cosm,lut)

    USE cosdef
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut
    REAL :: dc
    REAL :: amin, amax, af, zf, RHS, a, growz
    REAL, ALLOCATABLE :: af_tab(:), grow_tab(:)
    INTEGER :: i, j, ntab

    !This fills up the halo collapse redshift table as per Bullock relations       

    ntab=SIZE(cosm%growth)
    ALLOCATE(af_tab(ntab),grow_tab(ntab))

    af_tab=cosm%a_growth
    grow_tab=cosm%growth

    !Do numerical inversion
    DO i=1,lut%n

       !I don't think this is really consistent with dc varying as a function of z
       !But the change will be very small
       dc=delta_c(z,cosm)

       RHS=dc*grow(z,cosm)/lut%sigf(i)

       a=1./(1.+z)
       growz=find(a,af_tab,grow_tab,3,3)

       IF(RHS>growz) THEN
          zf=z
       ELSE
          af=find(RHS,grow_tab,af_tab,3,3)
          zf=-1.+1./af
       END IF

       lut%zc(i)=zf

    END DO

    DEALLOCATE(af_tab,grow_tab)

  END SUBROUTINE zcoll_bull

  FUNCTION mass_r(r,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: mass_r, r
    TYPE(cosmology) :: cosm
    REAL, PARAMETER :: pi=3.141592654

    !Relation between mean cosmological mass and radius

    mass_r=(4.*pi/3.)*cosmic_density(cosm)*(r**3.)

  END FUNCTION mass_r

  PURE FUNCTION cosmic_density(cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: cosmic_density
    TYPE(cosmology), INTENT(IN) :: cosm

    !In Msun per Mpc^3 with h factors included. The constant does this.
    cosmic_density=(2.775e11)*cosm%om_m

  END FUNCTION cosmic_density

  FUNCTION Tk(k,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: Tk, k
    TYPE(cosmology) :: cosm

    !Transfer function
    Tk=Tk_eh(k,cosm)

  END FUNCTION Tk

!!$  FUNCTION find_tk(k,cosm)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: find_tk
!!$    REAL :: kmin, kmax
!!$    REAL, INTENT(IN) :: k
!!$    INTEGER :: n
!!$    TYPE(cosmology) :: cosm
!!$
!!$    !Look-up and interpolation for T(k)
!!$
!!$    n=SIZE(cosm%ktab)
!!$    kmin=cosm%ktab(1)
!!$    kmax=cosm%ktab(n)
!!$
!!$    IF(k<kmin) THEN
!!$       !For k<<keq Tk=1.
!!$       find_tk=1.
!!$    ELSE IF(k>kmax) THEN
!!$       !Do some interpolation here based on knowledge of things at high k
!!$       find_tk=cosm%tktab(n)*(log(k)/log(kmax))*((k/kmax)**(-2.))
!!$    ELSE
!!$       !Otherwise use the standard find algorithm
!!$       find_tk=exp(find(log(k),log(cosm%ktab),log(cosm%tktab),3,3))
!!$    END IF
!!$
!!$  END FUNCTION find_tk

!!$  FUNCTION find_pk(k,cosm)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: find_pk
!!$    REAL :: kmax
!!$    REAL, INTENT(IN) :: k
!!$    INTEGER :: n
!!$    TYPE(cosmology) :: cosm
!!$
!!$    !Look-up and interpolation for P(k)
!!$
!!$    n=SIZE(cosm%ktab)
!!$    kmax=cosm%ktab(n)
!!$
!!$    IF(k>kmax) THEN
!!$       !Do some interpolation here based on knowledge of things at high k
!!$       find_pk=cosm%pktab(n)*((log(k)/log(kmax))**2.)*((k/kmax)**(cosm%n-1.))
!!$    ELSE
!!$       !Otherwise use the standard find algorithm
!!$       find_pk=exp(find(log(k),log(cosm%ktab),log(cosm%pktab),3,3))
!!$    END IF
!!$
!!$  END FUNCTION find_pk

  FUNCTION Tk_eh(yy,cosm)

    ! the astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])
    ! remember I use k/h, whereas they use pure k, om_m is cdm + baryons

    USE cosdef
    IMPLICIT NONE

    REAL :: Tk_eh
    REAL, INTENT(IN) :: yy
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL*8 :: rk, e, thet, b1, b2, zd, ze, rd, re, rke, s, rks
    REAL*8 :: q
    REAL*8 :: y, g, ab
    REAL*8 :: a1, a2, ac
    REAL*8 :: bc
    REAL*8 :: f, fac
    REAL*8 :: c1, c2, tc
    REAL*8 :: bb, bn, ss, tb
    REAL*8 :: om_m, om_b, h

    om_m=cosm%om_m
    om_b=cosm%om_b
    h=cosm%h

    rk=yy*h

    e=exp(1.)

    thet=2.728/2.7
    b1=0.313*(om_m*h*h)**(-0.419)*(1+0.607*(om_m*h*h)**0.674)
    b2=0.238*(om_m*h*h)**0.223
    zd=1291.*(1+b1*(om_b*h*h)**b2)*(om_m*h*h)**0.251/(1.+0.659*(om_m*h*h)**0.828)
    ze=2.50e4*om_m*h*h/thet**4.
    rd=31500.*om_b*h*h/thet**4./zd
    re=31500.*om_b*h*h/thet**4./ze
    rke=7.46e-2*om_m*h*h/thet**2.
    s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1+sqrt(re)))
    rks=1.6*( (om_b*h*h)**0.52 ) * ( (om_m*h*h)**0.73 ) * (1.+(10.4*om_m*h*h)**(-0.95))

    q=rk/13.41/rke

    y=(1.+ze)/(1.+zd)
    g=y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)))
    ab=g*2.07*rke*s/(1.+rd)**(0.75)

    a1=(46.9*om_m*h*h)**0.670*(1+(32.1*om_m*h*h)**(-0.532))
    a2=(12.0*om_m*h*h)**0.424*(1+(45.0*om_m*h*h)**(-0.582))
    ac=(a1**(-om_b/om_m)) * (a2**(-(om_b/om_m)**3.))

    b1=0.944/(1+(458.*om_m*h*h)**(-0.708))
    b2=(0.395*om_m*h*h)**(-0.0266)
    bc=1./(1.+b1*((1.-om_b/om_m)**b2-1.))

    f=1./(1.+(rk*s/5.4)**4.)

    c1=14.2 + 386./(1.+69.9*q**1.08)
    c2=14.2/ac + 386./(1.+69.9*q**1.08)
    tc=f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q*q) +(1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q*q)

    bb=0.5+(om_b/om_m) + (3.-2.*om_b/om_m)*sqrt((17.2*om_m*h*h)**2.+1.)
    bn=8.41*(om_m*h*h)**0.435
    ss=s/(1.+(bn/rk/s)**3.)**(1./3.)
    tb=log(e+1.8*q)/(log(e+1.8*q)+c1*q*q)/(1+(rk*s/5.2)**2.)
    IF((rk/rks**1.4)>7.) THEN
       fac=0.
    ELSE
       fac=exp(-(rk/rks)**1.4)
    END IF
    tb=(tb+ab*fac/(1.+(bb/rk/s)**3.))*sin(rk*ss)/rk/ss

    tk_eh=(om_b/om_m)*tb+(1-om_b/om_m)*tc

  END FUNCTION TK_EH

  FUNCTION p_lin(k,z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: p_lin
    REAL, INTENT (IN) :: k, z
    TYPE(cosmology), INTENT(IN) :: cosm

    !This gives the linear power spectrum for the model in question
    !P(k) should have been previously normalised so as to get the amplitude 'A' correct

    IF(k==0.) THEN
       !If p_lin happens to be foolishly called for 0 mode (which should never happen, but might in integrals)
       p_lin=0.
    ELSE IF(k>1.e8) THEN
       !Avoids some issues if p_lin is called for very (absurdly) high k values
       !For some reason crashes can occur if this is the case
       p_lin=0.
!    ELSE IF(cosm%itk==5) THEN
!       !itk==5 means P(k) has been taken as an input file
!       !In this case use the input P(k) file
!       p_lin=(cosm%A**2.)*(grow(z,cosm)**2.)*find_pk(k,cosm)
    ELSE
       !In this case look for the transfer function
       p_lin=(cosm%A**2.)*(grow(z,cosm)**2.)*(Tk(k,cosm)**2.)*(k**(cosm%n+3.))
    END IF

  END FUNCTION p_lin

  FUNCTION p_2h(k,z,plin,lut,cosm)

    USE cosdef
    REAL :: p_2h
    REAL, INTENT(IN) :: k, plin
    REAL :: sigv, frac
    REAL, INTENT(IN) :: z
    TYPE(tables), INTENT(IN) :: lut
    TYPE(cosmology), INTENT(IN) :: cosm

    !Produces the 'two-halo' power

    sigv=lut%sigv
    frac=fdamp(z,cosm)

    IF(frac==0.) THEN
       p_2h=plin
    ELSE
       p_2h=plin*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2.)
    END IF

    !For some extreme cosmologies frac>1. so this must be added to prevent p_2h<0.
    IF(p_2h<0.) p_2h=0.

  END FUNCTION p_2h

  FUNCTION p_1h(k,z,lut,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: p_1h
    REAL, INTENT(IN) :: k, z
    TYPE(tables), INTENT(IN) :: lut
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: Dv, g, fac, et, ks, wk
    REAL, ALLOCATABLE :: integrand(:)
    REAL :: sum
    INTEGER :: i
    REAL, PARAMETER :: pi=3.141592654

    !Does the one-halo power integral

    ALLOCATE(integrand(lut%n))
    integrand=0.
    
    !Only call eta once
    et=eta(z,cosm)

    !Calculates the value of the integrand at all nu values!
    DO i=1,lut%n
       g=gnu(lut%nu(i))
       wk=win(k*(lut%nu(i)**et),lut%rv(i),lut%c(i))
       integrand(i)=(lut%rv(i)**3.)*g*(wk**2.)
    END DO
    
    !Carries out the integration
    sum=inttab(lut%nu,REAL(integrand),1)
    
    DEALLOCATE(integrand)

    Dv=Delta_v(z,cosm)

    !These are just numerical factors from the 1-halo integral in terms of nu!
    p_1h=sum*2.*Dv*(k**3.)/(3.*pi)

    !Damping of the 1-halo term at very large scales
    ks=kstar(lut,cosm)

    !Prevents problems if k/ks is very large

    IF(ks>0.) THEN

       IF((k/ks)**2.>7.) THEN
          fac=0.
       ELSE
          fac=exp(-((k/ks)**2.))
       END IF

       p_1h=p_1h*(1.-fac)

    END IF

  END FUNCTION p_1h

  SUBROUTINE fill_sigtab(cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: rmin, rmax
    REAL, ALLOCATABLE :: rtab(:), sigtab(:)
    REAL :: r, sig
    INTEGER :: i, nsig
    TYPE(cosmology) :: cosm

    !This fills up tables of r vs. sigma(r) across a range in r!
    !It is used only in look-up for further calculations of sigma(r) and not otherwise!
    !and prevents a large number of calls to the sigint functions
    !rmin and rmax need to be decided in advance and are chosen such that
    !R vs. sigma(R) is approximately power-law below and above these values of R   
    !This wouldn't be appropriate for models with a linear spectrum cut-off (e.g. WDM)

    !These must be not allocated before sigma calculations otherwise when sigma(r) is called
    !otherwise sigma(R) looks for the result in the tables
    IF(ALLOCATED(cosm%r_sigma)) DEALLOCATE(cosm%r_sigma)
    IF(ALLOCATED(cosm%sigma1d)) DEALLOCATE(cosm%sigma1d)   

    !These values of 'r' work fine for any power spectrum of cosmological importance
    !Having nsig as a 2** number is most efficient for the look-up routines
    rmin=1e-4
    rmax=1e3
    nsig=64

    IF(ihm==1) WRITE(*,*) 'SIGTAB: Filling sigma interpolation table'
    IF(ihm==1) WRITE(*,*) 'SIGTAB: Rmin:', rmin
    IF(ihm==1) WRITE(*,*) 'SIGTAB: Rmax:', rmax
    IF(ihm==1) WRITE(*,*) 'SIGTAB: Values:', nsig

    ALLOCATE(rtab(nsig),sigtab(nsig))

    DO i=1,nsig

       !Equally spaced r in log
       r=exp(log(rmin)+log(rmax/rmin)*float(i-1)/float(nsig-1))

       sig=sigma(r,0.,cosm)

       rtab(i)=r
       sigtab(i)=sig

    END DO

    !Must be allocated after the sigtab calulation above
    ALLOCATE(cosm%r_sigma(nsig),cosm%sigma1d(nsig))

    cosm%r_sigma=rtab
    cosm%sigma1d=sigtab

    DEALLOCATE(rtab,sigtab)

    IF(ihm==1) WRITE(*,*) 'SIGTAB: Done'
    IF(ihm==1) WRITE(*,*)

  END SUBROUTINE fill_sigtab

  FUNCTION sigma(r,z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: sigma
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(r>=1.e-2) THEN

       sigma=sigint0(r,z,cosm)

    ELSE IF(r<1.e-2) THEN

       sigma=sqrt(sigint1(r,z,cosm)+sigint2(r,z,cosm))

    END IF

  END FUNCTION sigma

  FUNCTION sigma_v(R,z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: sigma_v
    REAL, INTENT(IN) :: z, R
    REAL*8 :: sum
    REAL :: alpha_min, alpha_max, alpha
    INTEGER :: l, nmax
    REAL :: dtheta, k, theta, oldsum, acc
    REAL, PARAMETER :: pi=3.141592654
    INTEGER :: i, j, n, ninit, jmax
    TYPE(cosmology), INTENT(IN) :: cosm

    acc=0.001
    ninit=100
    jmax=30

    alpha=1.65

    DO j=1,jmax

       n=ninit*(2**(j-1))

       sum=0.d0
       dtheta=1./float(n)

       DO i=2,n-1

          !theta converts integral to 0->1 range
          !Values at the end points are 0 so removed for convenience
          theta=float(i-1)/float(n-1)
          k=(-1.+1./theta)/r**alpha
          sum=sum+p_lin(k,z,cosm)*(wk_tophat(k*r)**2.)/((k**2.)*theta*(1.-theta))

       END DO

       sum=sum*dtheta

       IF(j>1 .AND. ABS(-1.+sum/oldsum)<acc) THEN
          !Convert from sigma_v^2 and to 1D dispersion
          sigma_v=sqrt(sum/3.)
          EXIT
       ELSE
          oldsum=sum
       END IF

    END DO

  END FUNCTION sigma_v

  FUNCTION sigma_cb(r,z,cosm)

    USE cosdef
    REAL :: sigma_cb
    REAL, INTENT(IN) :: r, z
    REAL :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    !Finds sigma_cold from look-up tables
    
    sigma_cb=grow(z,cosm)*exp(find(log(r),log(cosm%r_sigma),log(cosm%sigma1d),3,3))

  END FUNCTION sigma_cb

  FUNCTION wk_tophat(x)

    IMPLICIT NONE
    REAL :: wk_tophat, x

    !The normlaised Fourier Transform of a top-hat
    !Taylor expansion used for low |x| to avoid cancellation problems

    IF(x<0.01) THEN
       wk_tophat=1.-(x**2.)/10.
    ELSE
       wk_tophat=3.*(sin(x)-x*cos(x))/(x**3.)
    END IF

  END FUNCTION wk_tophat

  FUNCTION inttab(x,y,iorder)

    IMPLICIT NONE
    REAL :: inttab
    REAL, INTENT(IN) :: x(:), y(:)
    REAL :: a, b, c, d, h
    REAL :: q1, q2, q3, qi, qf
    REAL :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
    REAL*8 :: sum
    INTEGER :: i, n, i1, i2, i3, i4
    INTEGER, INTENT(IN) :: iorder

    !Routine to integrate tables of data using the trapezium rule
    !Can either use linear, quadratic or cubic methods

    n=SIZE(x)

    IF(n .NE. SIZE(y)) STOP 'Tables must be of the same length'

    sum=0.d0

    IF(iorder==1) THEN

       !Sums over all Trapezia (a+b)*h/2
       DO i=1,n-1
          a=y(i+1)
          b=y(i)
          h=x(i+1)-x(i)
          sum=sum+(a+b)*h/2.d0
       END DO

    ELSE IF(iorder==2) THEN

       DO i=1,n-2

          x1=x(i)
          x2=x(i+1)
          x3=x(i+2)

          y1=y(i)
          y2=y(i+1)
          y3=y(i+2)

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          q1=a*(x1**3.)/3.+b*(x1**2.)/2.+c*x1
          q2=a*(x2**3.)/3.+b*(x2**2.)/2.+c*x2
          q3=a*(x3**3.)/3.+b*(x3**2.)/2.+c*x3

          !Takes value for first and last sections but averages over sections where you
          !have two independent estimates of the area
          IF(n==3) THEN
             sum=sum+q3-q1
          ELSE IF(i==1) THEN
             sum=sum+(q2-q1)+(q3-q2)/2.d0
          ELSE IF(i==n-2) THEN
             sum=sum+(q2-q1)/2.d0+(q3-q2)
          ELSE
             sum=sum+(q3-q1)/2.
          END IF

       END DO

    ELSE IF(iorder==3) THEN

       DO i=1,n-1

          !First choose the integers used for defining cubics for each section
          !First and last are different because the section does not lie in the *middle* of a cubic

          IF(i==1) THEN

             i1=1
             i2=2
             i3=3
             i4=4

          ELSE IF(i==n-1) THEN

             i1=n-3
             i2=n-2
             i3=n-1
             i4=n

          ELSE

             i1=i-1
             i2=i
             i3=i+1
             i4=i+2

          END IF

          x1=x(i1)
          x2=x(i2)
          x3=x(i3)
          x4=x(i4)

          y1=y(i1)
          y2=y(i2)
          y3=y(i3)
          y4=y(i4)

          CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

          !These are the limits of the particular section of integral
          xi=x(i)
          xf=x(i+1)

          qi=a*(xi**4.)/4.+b*(xi**3.)/3.+c*(xi**2.)/2.+d*xi
          qf=a*(xf**4.)/4.+b*(xf**3.)/3.+c*(xf**2.)/2.+d*xf

          sum=sum+qf-qi

       END DO

    END IF

    inttab=sum

  END FUNCTION inttab

  FUNCTION sigma_integrand(t,R,f,z,cosm)

    USE cosdef
    REAL :: sigma_integrand
    REAL, INTENT(IN) :: t, R, z
    REAL :: k, y, w_hat
    TYPE(cosmology), INTENT(IN) :: cosm

    INTERFACE
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    !Integrand to the sigma integral in terms of t. Defined by k=(1/t-1)/f(R) where f(R) is *any* function

    IF(t==0.) THEN
       !t=0 corresponds to k=infintiy when W(kR)=0.
       sigma_integrand=0.
    ELSE IF(t==1.) THEN
       !t=1 corresponds to k=0. when P(k)=0.
       sigma_integrand=0.
    ELSE
       !f(R) can be *any* function of R here to improve integration speed
       k=(-1.+1./t)/f(R)
       y=k*R
       w_hat=wk_tophat(y)
       sigma_integrand=p_lin(k,z,cosm)*(w_hat**2.)/(t*(1.-t))
    END IF

  END FUNCTION sigma_integrand

  FUNCTION f_rapid(r)

    IMPLICIT NONE
    REAL :: f_rapid
    REAL, INTENT(IN) :: r
    REAL :: alpha

    !This is the 'rapidising' function to increase integration speed
    !for sigma(R). Found by trial-and-error

    IF(r>1.e-2) THEN
       !alpha 0.3-0.5 works well
       alpha=0.5
    ELSE
       !If alpha=1 this goes tits up
       !alpha 0.7-0.9 works well
       alpha=0.8
    END IF

    f_rapid=r**alpha

  END FUNCTION f_rapid

  FUNCTION sigint0(r,z,cosm)

    !Integrates between a and b until desired accuracy is reached!

    USE cosdef
    IMPLICIT NONE
    REAL, INTENT(IN) :: r, z
    INTEGER :: i, j, jmax
    REAL :: sigint0, acc, dx
    INTEGER :: ninit, n
    REAL :: x
    REAL*8 :: sum1, sum2
    TYPE(cosmology), INTENT(IN) :: cosm

    acc=0.001

    sum1=0.d0
    sum2=0.d0

    ninit=50
    jmax=20
    
    DO j=1,jmax

       n=ninit*2**(j-1)

       !Avoids the end-points where the integrand is 0 anyway
       DO i=2,n-1

          !x is defined on the interval 0 -> 1
          x=float(i-1)/float(n-1)

          sum2=sum2+sigma_integrand(x,r,f_rapid,z,cosm)

       END DO

       dx=1./float(n-1)
       sum2=sum2*dx
       sum2=sqrt(sum2)

       IF(j==1) THEN
          sum1=sum2
       ELSE IF(ABS(-1.+sum2/sum1)<acc) THEN
          sigint0=sum2
          EXIT
       ELSE IF(j==jmax) THEN
          WRITE(*,*)
          WRITE(*,*) 'SIGINT: r:', r
          WRITE(*,*) 'SIGINT: Integration timed out'
          WRITE(*,*)
          STOP
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION sigint0

  FUNCTION sigint1(r,z,cosm)

    !Integrates between a and b until desired accuracy is reached!

    USE cosdef
    IMPLICIT NONE
    REAL, INTENT(IN) :: r, z
    INTEGER :: i, j, jmax
    REAL :: sigint1, acc, dx
    INTEGER :: ninit, n
    REAL :: x, fac, xmin, xmax, k
    REAL*8 :: sum1, sum2
    TYPE(cosmology), INTENT(IN) :: cosm

    acc=0.001

    sum1=0.d0
    sum2=0.d0

    ninit=50
    jmax=20

    xmin=r/(r+r**.5)
    xmax=1.
    
    DO j=1,jmax

       n=ninit*2**(j-1)

       !Avoids the end-point where the integrand is 0 anyway
       DO i=1,n-1

          x=xmin+(xmax-xmin)*float(i-1)/float(n-1)

          IF(i==1 .OR. i==n) THEN
             fac=0.5
          ELSE
             fac=1.
          END IF

          k=(-1.+1./x)/r**.5
          sum2=sum2+fac*p_lin(k,z,cosm)*(wk_tophat(k*r)**2.)/(x*(1.-x))

       END DO

       dx=(xmax-xmin)/float(n-1)
       sum2=sum2*dx

       IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
          sigint1=sum2
          EXIT
       ELSE IF(j==jmax) THEN
          WRITE(*,*)
          WRITE(*,*) 'SIGINT1: r:', r
          WRITE(*,*) 'SIGINT1: Integration timed out'
          WRITE(*,*)
          STOP
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION sigint1

  FUNCTION sigint2(r,z,cosm)

    !Integrates between a and b until desired accuracy is reached!

    USE cosdef
    IMPLICIT NONE
    REAL, INTENT(IN) :: r, z
    INTEGER :: i, j, jmax
    REAL :: sigint2, acc, dx
    INTEGER :: ninit, n
    REAL :: x, fac, xmin, xmax, A
    REAL*8 :: sum1, sum2
    TYPE(cosmology), INTENT(IN) :: cosm

    acc=0.001

    sum1=0.d0
    sum2=0.d0

    ninit=50
    jmax=20

    !How far to go out in 1/r units for integral
    A=10.

    xmin=1./r
    xmax=A/r
    
    DO j=1,jmax

       n=ninit*2**(j-1)

       DO i=1,n

          x=xmin+(xmax-xmin)*float(i-1)/float(n-1)

          IF(i==1 .OR. i==n) THEN
             fac=0.5
          ELSE
             fac=1.
          END IF          

          !Integrate linearly in k for the rapidly oscillating part
          sum2=sum2+fac*p_lin(x,z,cosm)*(wk_tophat(x*r)**2.)/x

       END DO

       dx=(xmax-xmin)/float(n-1)
       sum2=sum2*dx

       IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
          sigint2=sum2
          EXIT
       ELSE IF(j==jmax) THEN
          WRITE(*,*)
          WRITE(*,*) 'SIGINT2: r:', r
          WRITE(*,*) 'SIGINT2: Integration timed out'
          WRITE(*,*)
          STOP
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION sigint2

  FUNCTION win(k,rv,c)

    IMPLICIT NONE
    REAL :: win, k, rv, c

    !Calls the analytic Fourier Transform of the NFW profile
    win=winnfw(k,rv,c)

    !Correct for the case of disasters (a bit sloppy, not sure if this is ever used)
    IF(win>1.) win=1.
    IF(win<0.) win=0.

  END FUNCTION win

  FUNCTION winnfw(k,rv,c)

    IMPLICIT NONE
    REAL :: winnfw
    REAL, INTENT(IN) :: k, rv, c
    REAL :: si1, si2, ci1, ci2, ks
    REAL :: p1, p2, p3

    !The analytic Fourier Transform of the NFW

    ks=k*rv/c

    si1=si(ks)
    si2=si((1.+c)*ks)
    ci1=ci(ks)
    ci2=ci((1.+c)*ks)

    p1=cos(ks)*(ci2-ci1)
    p2=sin(ks)*(si2-si1)
    p3=sin(ks*c)/(ks*(1.+c))

    winnfw=p1+p2-p3
    winnfw=winnfw/mass(c)

  END FUNCTION winnfw

  FUNCTION mass(c)

    !This calculates the (normalised) mass of a halo of concentration c
    !This is mass with some factors missing

    IMPLICIT NONE
    REAL :: mass, c

    mass=log(1.+c)-c/(1.+c)

  END FUNCTION mass

  FUNCTION gnu(nu)

    IMPLICIT NONE
    REAL :: gnu, nu

    !Mass function

    gnu=gst(nu)

  END FUNCTION gnu

  FUNCTION gst(nu)

    !Sheth Tormen mass function!
    !Note I use nu=dc/sigma(M) and this Sheth & Tormen (1999) use nu=(dc/sigma)^2
    !This accounts for some small differences

    IMPLICIT NONE
    REAL :: nu, gst
    REAL :: p, q

    p=0.3
    q=0.707

    gst=0.21616*(1.+((q*nu*nu)**(-p)))*exp(-q*nu*nu/2.)

  END FUNCTION gst

!!$  FUNCTION hubble2(z,cosm)
!!$
!!$    !This calculates the dimensionless squared hubble parameter at redshift z!
!!$    !and it ignores contributions from radiation (not accurate at high z, but consistent with simulations)!
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: hubble2
!!$    REAL, INTENT(IN) :: z
!!$    REAL :: om_m, om_v, w
!!$    TYPE(cosmology), INTENT(IN) :: cosm
!!$
!!$    om_m=cosm%om_m
!!$    om_v=cosm%om_v
!!$    w=cosm%w
!!$
!!$    hubble2=(om_m*(1.+z)**3.)+om_v*((1.+z)**(3.*(1.+w)))+((1.-om_m-om_v)*(1.+z)**2.)
!!$
!!$  END FUNCTION hubble2

  FUNCTION hubble2(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: hubble2
    REAL, INTENT(IN) :: z
    REAL :: om_m, om_v, a
    TYPE(cosmology), INTENT(IN) :: cosm

    !Calculates Hubble^2 in units such that H^2(z=0)=1.
    om_m=cosm%om_m
    om_v=cosm%om_v
    a=1./(1.+z)
    hubble2=(om_m*(1.+z)**3.)+om_v*x_de(a,cosm)+((1.-om_m-om_v)*(1.+z)**2.)

  END FUNCTION hubble2

  FUNCTION omega_m(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: omega_m
    REAL, INTENT(IN) :: z
    REAL :: om_m
    TYPE(cosmology), INTENT(IN) :: cosm

    !This calculates Omega_m variations with z!
    om_m=cosm%om_m
    omega_m=(om_m*(1.+z)**3.)/hubble2(z,cosm)

  END FUNCTION omega_m

!!$  FUNCTION omega_v(z,cosm)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: omega_v
!!$    REAL, INTENT(IN) :: z
!!$    REAL :: om_v, w
!!$    TYPE(cosmology), INTENT(IN) :: cosm
!!$
!!$    !This calculates Omega_v variations with z for any w
!!$    om_v=cosm%om_v
!!$    w=cosm%w
!!$
!!$    omega_v=om_v*((1.+z)**(3.*(1.+w)))/hubble2(z,cosm)
!!$
!!$  END FUNCTION omega_v

  FUNCTION omega_v(z,cosm)
    
    USE cosdef
    IMPLICIT NONE
    REAL :: omega_v
    REAL, INTENT(IN) :: z
    REAL :: om_v, a
    TYPE(cosmology), INTENT(IN) :: cosm

    !This calculates omega_v (really omega_de) variations with z!
    om_v=cosm%om_v
    a=1./(1.+z)
    omega_v=om_v*x_de(a,cosm)/hubble2(z,cosm)

  END FUNCTION omega_v

  FUNCTION grow(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: grow
    REAL, INTENT(IN) :: z
    REAL :: a, acc
    TYPE(cosmology), INTENT(IN) :: cosm

    !Scale-independent growth function at z=0
    IF(z==0.) THEN
       grow=1.
    ELSE
       a=1./(1.+z)
       grow=find(a,cosm%a_growth,cosm%growth,3,3)
    END IF

  END FUNCTION grow

  FUNCTION grow_int(b,acc,cosm)

    !Integrates between a and b with nint points until desired accuracy is reached!
    USE cosdef
    IMPLICIT NONE
    INTEGER :: i, j, jmax
    REAL :: grow_int, a, b, acc, dx
    INTEGER :: nint
    REAL :: x, fac, func, gam
    REAL*8 :: sum1, sum2
    TYPE(cosmology) :: cosm

    sum1=0.d0
    sum2=0.d0

    jmax=20

    a=1.
    
    IF(cosm%wa .NE. 0.) STOP 'This does not work for w(a)'

    IF(a==b) THEN

       grow_int=1.

    ELSE

       DO j=1,jmax

          nint=10.*(2.**j)

          DO i=1,nint

             x=a+(b-a)*((float(i)-1)/(float(nint)-1))

             IF(i==1 .OR. i==nint) THEN
                !multiple of 1 for beginning and end and multiple of 2 for middle points!
                fac=1.
             ELSE
                fac=2.
             END IF

             !Insert function here!
             IF(cosm%w<-1.) THEN
                gam=0.55+0.02*(1.+cosm%w)
             ELSE IF(cosm%w>-1) THEN
                gam=0.55+0.05*(1.+cosm%w)
             ELSE
                gam=0.55
             END IF

             func=(omega_m(-1.+1./x,cosm)**gam)/x

             sum2=sum2+fac*func

          END DO

          dx=((b-a)/(float(nint)-1.))
          sum2=sum2*dx/2.

          IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
             grow_int=exp(sum2)
             EXIT
          ELSE
             sum1=sum2
             sum2=0.d0
          END IF

       END DO

    END IF

  END FUNCTION grow_int

  FUNCTION dispint(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: dispint
    REAL, INTENT(IN) :: z
    REAL*8 :: sum
    REAL :: dtheta, k, theta, oldsum, acc
    REAL, PARAMETER :: pi=3.141592654
    INTEGER :: i, j, n, ninit, jmax
    TYPE(cosmology), INTENT(IN) :: cosm

    acc=0.001
    ninit=100
    jmax=30

    DO j=1,jmax

       n=ninit*(2**(j-1))

       sum=0.d0
       dtheta=1./float(n)

       DO i=2,n-1

          !theta converts integral to 0->1 range
          !Values at the end points are 0 so removed for convenience
          theta=float(i-1)/float(n-1)
          k=(-1.+1./theta)          
          sum=sum+((1.+k)**2.)*p_lin(k,z,cosm)/(k**3.)!((k**3.)*theta**2.)

       END DO

       sum=sum*dtheta

       IF(j>1 .AND. ABS(-1.+sum/oldsum)<acc) THEN  
          dispint=sum/3.
          EXIT
       ELSE
          oldsum=sum
       END IF

    END DO

  END FUNCTION dispint

  FUNCTION si(x)

    IMPLICIT NONE
    REAL :: si, x
    REAL*8 :: x2, y, f, g, si8
    REAL*8, PARAMETER :: pi=3.1415926535897932384626433d0

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.

    IF(ABS(x)<=4.) THEN

       x2=x*x

       si8 = x*(1.d0+x2*(-4.54393409816329991d-2+x2*(1.15457225751016682d-3&
            +x2*(-1.41018536821330254d-5+x2*(9.43280809438713025d-8+x2*(-3.53201978997168357d-10&
            +x2*(7.08240282274875911d-13+x2*(-6.05338212010422477d-16))))))))/ &
            (1.+x2*(1.01162145739225565d-2 +x2*(4.99175116169755106d-5+&
            x2*(1.55654986308745614d-7+x2*(3.28067571055789734d-10+x2*(4.5049097575386581d-13&
            +x2*(3.21107051193712168d-16)))))))

       si=si8

    ELSE IF(ABS(x)>4.) THEN

       y=1.d0/(x*x)

       f = (1.d0 + y*(7.44437068161936700618d2 + y*(1.96396372895146869801d5 +&
            y*(2.37750310125431834034d7 +y*(1.43073403821274636888d9 + y*(4.33736238870432522765d10 &
            + y*(6.40533830574022022911d11 + y*(4.20968180571076940208d12 + &
            y*(1.00795182980368574617d13 + y*(4.94816688199951963482d12 +&
            y*(-4.94701168645415959931d11)))))))))))/ (x*(1. +y*(7.46437068161927678031d2 +&
            y*(1.97865247031583951450d5 +y*(2.41535670165126845144d7 + &
            y*(1.47478952192985464958d9 + y*(4.58595115847765779830d10 +&
            y*(7.08501308149515401563d11 + y*(5.06084464593475076774d12 + &
            y*(1.43468549171581016479d13 + y*(1.11535493509914254097d13)))))))))))


       g = y*(1.d0 + y*(8.1359520115168615d2 + y*(2.35239181626478200d5 + &
            y*(3.12557570795778731d7 + y*(2.06297595146763354d9 + y*(6.83052205423625007d10 +&
            y*(1.09049528450362786d12 + y*(7.57664583257834349d12 +y*(1.81004487464664575d13 +&
            y*(6.43291613143049485d12 +y*(-1.36517137670871689d12)))))))))))/&
            (1. + y*(8.19595201151451564d2 +y*(2.40036752835578777d5 + y*(3.26026661647090822d7 &
            + y*(2.23355543278099360d9 + y*(7.87465017341829930d10 + y*(1.39866710696414565d12 &
            + y*(1.17164723371736605d13 + y*(4.01839087307656620d13 +y*(3.99653257887490811d13))))))))))

       si=pi/2.d0-f*cos(x)-g*sin(x)

    END IF

  END FUNCTION si

  FUNCTION ci(x)

    IMPLICIT NONE
    REAL :: ci, x
    REAL*8 :: x2, y, f, g, ci8
    REAL*8, PARAMETER :: em_const=0.577215664901532861d0

    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.

    IF(ABS(x)<=4.) THEN

       x2=x*x

       ci8=em_const+log(x)+x2*(-0.25d0+x2*(7.51851524438898291d-3+x2*(-1.27528342240267686d-4&
            +x2*(1.05297363846239184d-6+x2*(-4.68889508144848019d-9+x2*(1.06480802891189243d-11&
            +x2*(-9.93728488857585407d-15)))))))/ (1.+x2*(1.1592605689110735d-2+&
            x2*(6.72126800814254432d-5+x2*(2.55533277086129636d-7+x2*(6.97071295760958946d-10+&
            x2*(1.38536352772778619d-12+x2*(1.89106054713059759d-15+x2*(1.39759616731376855d-18))))))))

       ci=ci8

    ELSE IF(ABS(x)>4.) THEN

       y=1./(x*x) 

       f = (1.d0 + y*(7.44437068161936700618d2 + y*(1.96396372895146869801d5 + &
            y*(2.37750310125431834034d7 +y*(1.43073403821274636888d9 + y*(4.33736238870432522765d10&
            + y*(6.40533830574022022911d11 + y*(4.20968180571076940208d12 + y*(1.00795182980368574617d13&
            + y*(4.94816688199951963482d12 +y*(-4.94701168645415959931d11)))))))))))/&
            (x*(1. +y*(7.46437068161927678031d2 +y*(1.97865247031583951450d5 +&
            y*(2.41535670165126845144d7 + y*(1.47478952192985464958d9 + &
            y*(4.58595115847765779830d10 +y*(7.08501308149515401563d11 + y*(5.06084464593475076774d12 &
            + y*(1.43468549171581016479d13 + y*(1.11535493509914254097d13)))))))))))   

       g = y*(1.d0 + y*(8.1359520115168615d2 + y*(2.35239181626478200d5 + y*(3.12557570795778731d7&
            + y*(2.06297595146763354d9 + y*(6.83052205423625007d10 +&
            y*(1.09049528450362786d12 + y*(7.57664583257834349d12 +&
            y*(1.81004487464664575d13 + y*(6.43291613143049485d12 +y*(-1.36517137670871689d12)))))))))))&
            / (1. + y*(8.19595201151451564d2 +y*(2.40036752835578777d5 +&
            y*(3.26026661647090822d7 + y*(2.23355543278099360d9 + y*(7.87465017341829930d10 &
            + y*(1.39866710696414565d12 + y*(1.17164723371736605d13 + y*(4.01839087307656620d13 +y*(3.99653257887490811d13))))))))))

       ci=f*sin(x)-g*cos(x)

    END IF

  END FUNCTION ci

!!$  FUNCTION derivative_table(x,xin,yin,iorder,imeth)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: derivative_table
!!$    REAL, INTENT(IN) :: x, xin(:), yin(:)
!!$    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
!!$    REAL :: a, b, c, d
!!$    REAL :: x1, x2, x3, x4
!!$    REAL :: y1, y2, y3, y4
!!$    INTEGER :: i, n
!!$    INTEGER, INTENT(IN) :: imeth, iorder
!!$    INTEGER :: maxorder, maxmethod
!!$
!!$    !Finds the derivative f'(x) given tables x, f(x)
!!$
!!$    !This version interpolates if the value is off either end of the array!
!!$    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
!!$    !Results from the interpolation!
!!$
!!$    !imeth = 1 => find x in xtab by crudely searching
!!$    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
!!$    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))
!!$
!!$    !iorder = 1 => linear interpolation
!!$    !iorder = 2 => quadratic interpolation
!!$    !iorder = 3 => cubic interpolation
!!$
!!$    n=SIZE(xtab)
!!$
!!$    maxorder=3
!!$    maxmethod=3
!!$
!!$    n=SIZE(xin)
!!$    IF(n .NE. SIZE(yin)) STOP 'FIND: Tables not of the same size'
!!$    ALLOCATE(xtab(n),ytab(n))
!!$
!!$    xtab=xin
!!$    ytab=yin
!!$
!!$    IF(xtab(1)>xtab(n)) THEN
!!$       !Reverse the arrays in this case
!!$       CALL reverse(xtab)
!!$       CALL reverse(ytab)
!!$    END IF
!!$
!!$    IF(iorder<1) STOP 'FIND: find order not specified correctly'
!!$    IF(iorder>maxorder) STOP 'FIND: find order not specified correctly'
!!$    IF(imeth<1) STOP 'FIND: Method of finding within a table not specified correctly'
!!$    IF(imeth>maxmethod) STOP 'FIND: Method of finding within a table not specified correctly'
!!$
!!$    IF(iorder==1) THEN
!!$
!!$       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'
!!$
!!$       IF(x<=xtab(2)) THEN
!!$
!!$          x2=xtab(2)
!!$          x1=xtab(1)
!!$
!!$          y2=ytab(2)
!!$          y1=ytab(1)
!!$
!!$       ELSE IF (x>=xtab(n-1)) THEN
!!$
!!$          x2=xtab(n)
!!$          x1=xtab(n-1)
!!$
!!$          y2=ytab(n)
!!$          y1=ytab(n-1)
!!$
!!$       ELSE
!!$
!!$          IF(imeth==1) i=search_int(x,xtab)
!!$          IF(imeth==2) i=linear_table_integer(x,xtab)
!!$          IF(imeth==3) i=int_split(x,xtab)
!!$
!!$          x2=xtab(i+1)
!!$          x1=xtab(i)
!!$
!!$          y2=ytab(i+1)
!!$          y1=ytab(i)
!!$
!!$       END IF
!!$
!!$       CALL fit_line(a,b,x1,y1,x2,y2)
!!$       derivative_table=a
!!$
!!$    ELSE IF(iorder==2) THEN
!!$
!!$       IF(n<3) STOP 'FIND_QUADRATIC: Not enough points in your table'
!!$
!!$       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN
!!$
!!$          IF(x<=xtab(2)) THEN
!!$
!!$             x3=xtab(3)
!!$             x2=xtab(2)
!!$             x1=xtab(1)
!!$
!!$             y3=ytab(3)
!!$             y2=ytab(2)
!!$             y1=ytab(1)
!!$
!!$          ELSE IF (x>=xtab(n-1)) THEN
!!$
!!$             x3=xtab(n)
!!$             x2=xtab(n-1)
!!$             x1=xtab(n-2)
!!$
!!$             y3=ytab(n)
!!$             y2=ytab(n-1)
!!$             y1=ytab(n-2)
!!$
!!$          END IF
!!$
!!$          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
!!$
!!$          derivative_table=2.*a*x+b
!!$
!!$       ELSE
!!$
!!$          IF(imeth==1) i=search_int(x,xtab)
!!$          IF(imeth==2) i=linear_table_integer(x,xtab)
!!$          IF(imeth==3) i=int_split(x,xtab)
!!$
!!$          x1=xtab(i-1)
!!$          x2=xtab(i)
!!$          x3=xtab(i+1)
!!$          x4=xtab(i+2)
!!$
!!$          y1=ytab(i-1)
!!$          y2=ytab(i)
!!$          y3=ytab(i+1)
!!$          y4=ytab(i+2)
!!$
!!$          !In this case take the average of two separate quadratic spline values
!!$
!!$          derivative_table=0.
!!$
!!$          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
!!$          derivative_table=derivative_table+(2.*a*x+b)/2.
!!$
!!$          CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
!!$          derivative_table=derivative_table+(2.*a*x+b)/2.
!!$
!!$       END IF
!!$
!!$    ELSE IF(iorder==3) THEN
!!$
!!$       IF(n<4) STOP 'FIND_CUBIC: Not enough points in your table'
!!$
!!$       IF(x<=xtab(3)) THEN
!!$
!!$          x4=xtab(4)
!!$          x3=xtab(3)
!!$          x2=xtab(2)
!!$          x1=xtab(1)
!!$
!!$          y4=ytab(4)
!!$          y3=ytab(3)
!!$          y2=ytab(2)
!!$          y1=ytab(1)
!!$
!!$       ELSE IF (x>=xtab(n-2)) THEN
!!$
!!$          x4=xtab(n)
!!$          x3=xtab(n-1)
!!$          x2=xtab(n-2)
!!$          x1=xtab(n-3)
!!$
!!$          y4=ytab(n)
!!$          y3=ytab(n-1)
!!$          y2=ytab(n-2)
!!$          y1=ytab(n-3)
!!$
!!$       ELSE
!!$
!!$          IF(imeth==1) i=search_int(x,xtab)
!!$          IF(imeth==2) i=linear_table_integer(x,xtab)
!!$          IF(imeth==3) i=int_split(x,xtab)
!!$
!!$          x1=xtab(i-1)
!!$          x2=xtab(i)
!!$          x3=xtab(i+1)
!!$          x4=xtab(i+2)
!!$
!!$          y1=ytab(i-1)
!!$          y2=ytab(i)
!!$          y3=ytab(i+1)
!!$          y4=ytab(i+2)
!!$
!!$       END IF
!!$
!!$       CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
!!$       derivative_table=3.*a*(x**2.)+2.*b*x+c
!!$
!!$    END IF
!!$
!!$  END FUNCTION derivative_table

   FUNCTION derivative_table(x,xin,yin,iorder,imeth)

    IMPLICIT NONE
    REAL :: derivative_table
    REAL, INTENT(IN) :: x, xin(:), yin(:)
    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i, n
    INTEGER, INTENT(IN) :: imeth, iorder
    INTEGER :: maxorder, maxmethod

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !imeth = 1 => find x in xtab by crudely searching
    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    n=SIZE(xtab)

    maxorder=3
    maxmethod=3

    n=SIZE(xin)
    IF(n .NE. SIZE(yin)) STOP 'FIND: Tables not of the same size'
    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       !Reverse the arrays in this case
       CALL reverse(xtab)
       CALL reverse(ytab)
    END IF

    IF(iorder<1) STOP 'FIND: find order not specified correctly'
    IF(iorder>maxorder) STOP 'FIND: find order not specified correctly'
    IF(imeth<1) STOP 'FIND: Method of finding within a table not specified correctly'
    IF(imeth>maxmethod) STOP 'FIND: Method of finding within a table not specified correctly'

!    IF(xtab(1)>xtab(n)) STOP 'DERIVATIVE_TABLE: x table in wrong order'
!    IF(n .NE. SIZE(ytab)) STOP 'DERIVATIVE_TABLE: Tables not of the same size'
!    IF(iorder<1) STOP 'DERIVATIVE_TABLE: find order not specified correctly'
!    IF(iorder>maxorder) STOP 'DERIVATIVE_TABLE: find order not specified correctly'
!    IF(imeth<1) STOP 'DERIVATIVE_TABLE: Method of finding within a table not specified correctly'
!    IF(imeth>maxmethod) STOP 'DERIVATIVE_TABLE: Method of finding within a table not specified correctly'

    IF(iorder==1) THEN

       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          x2=xtab(2)
          x1=xtab(1)

          y2=ytab(2)
          y1=ytab(1)

       ELSE IF (x>=xtab(n-1)) THEN

          x2=xtab(n)
          x1=xtab(n-1)

          y2=ytab(n)
          y1=ytab(n-1)

       ELSE

          i=table_integer(x,xtab,imeth)
          !IF(imeth==1) i=search_int(x,xtab)
          !IF(imeth==2) i=linear_table_integer(x,xtab)
          !IF(imeth==3) i=int_split(x,xtab)

          x2=xtab(i+1)
          x1=xtab(i)

          y2=ytab(i+1)
          y1=ytab(i)

       END IF

       CALL fit_line(a,b,x1,y1,x2,y2)
       derivative_table=a

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'FIND_QUADRATIC: Not enough points in your table'

       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

          IF(x<=xtab(2)) THEN

             x3=xtab(3)
             x2=xtab(2)
             x1=xtab(1)

             y3=ytab(3)
             y2=ytab(2)
             y1=ytab(1)

          ELSE IF (x>=xtab(n-1)) THEN

             x3=xtab(n)
             x2=xtab(n-1)
             x1=xtab(n-2)

             y3=ytab(n)
             y2=ytab(n-1)
             y1=ytab(n-2)

          END IF

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          derivative_table=2.*a*x+b

       ELSE

          i=table_integer(x,xtab,imeth)
          !IF(imeth==1) i=search_int(x,xtab)
          !IF(imeth==2) i=linear_table_integer(x,xtab)
          !IF(imeth==3) i=int_split(x,xtab)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          !In this case take the average of two separate quadratic spline values

          derivative_table=0.

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
          derivative_table=derivative_table+(2.*a*x+b)/2.

          CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
          derivative_table=derivative_table+(2.*a*x+b)/2.

       END IF

    ELSE IF(iorder==3) THEN

       IF(n<4) STOP 'FIND_CUBIC: Not enough points in your table'

       IF(x<=xtab(3)) THEN

          x4=xtab(4)
          x3=xtab(3)
          x2=xtab(2)
          x1=xtab(1)

          y4=ytab(4)
          y3=ytab(3)
          y2=ytab(2)
          y1=ytab(1)

       ELSE IF (x>=xtab(n-2)) THEN

          x4=xtab(n)
          x3=xtab(n-1)
          x2=xtab(n-2)
          x1=xtab(n-3)

          y4=ytab(n)
          y3=ytab(n-1)
          y2=ytab(n-2)
          y1=ytab(n-3)

       ELSE

          i=table_integer(x,xtab,imeth)
          !IF(imeth==1) i=search_int(x,xtab)
          !IF(imeth==2) i=linear_table_integer(x,xtab)
          !IF(imeth==3) i=int_split(x,xtab)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
       derivative_table=3.*a*(x**2.)+2.*b*x+c

    END IF

  END FUNCTION derivative_table
  
!!$  FUNCTION find(x,xin,yin,iorder,imeth)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: find
!!$    REAL, INTENT(IN) :: x, xin(:), yin(:)
!!$    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
!!$    REAL :: a, b, c, d
!!$    REAL :: x1, x2, x3, x4
!!$    REAL :: y1, y2, y3, y4
!!$    INTEGER :: i, n
!!$    INTEGER, INTENT(IN) :: imeth, iorder
!!$    INTEGER :: maxorder, maxmethod
!!$
!!$    !Interpolation routine.
!!$
!!$    !This version interpolates if the value is off either end of the array!
!!$    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
!!$    !Results from the interpolation!
!!$
!!$    !If the value required is off the table edge the interpolation is always linear
!!$
!!$    !imeth = 1 => find x in xtab by crudely searching from x(1) to x(n)
!!$    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
!!$    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))
!!$
!!$    !iorder = 1 => linear interpolation
!!$    !iorder = 2 => quadratic interpolation
!!$    !iorder = 3 => cubic interpolation
!!$
!!$    maxorder=3
!!$    maxmethod=3
!!$
!!$    n=SIZE(xin)
!!$    IF(n .NE. SIZE(yin)) STOP 'FIND: Tables not of the same size'
!!$    ALLOCATE(xtab(n),ytab(n))
!!$
!!$    xtab=xin
!!$    ytab=yin
!!$
!!$    IF(xtab(1)>xtab(n)) THEN
!!$       !Reverse the arrays in this case
!!$       CALL reverse(xtab)
!!$       CALL reverse(ytab)
!!$    END IF
!!$
!!$    IF(iorder<1) STOP 'FIND: find order not specified correctly'
!!$    IF(iorder>maxorder) STOP 'FIND: find order not specified correctly'
!!$    IF(imeth<1) STOP 'FIND: Method of finding within a table not specified correctly'
!!$    IF(imeth>maxmethod) STOP 'FIND: Method of finding within a table not specified correctly'
!!$
!!$    IF(x<xtab(1)) THEN
!!$
!!$       x1=xtab(1)
!!$       x2=xtab(2)
!!$
!!$       y1=ytab(1)
!!$       y2=ytab(2)
!!$
!!$       CALL fit_line(a,b,x1,y1,x2,y2)
!!$       find=a*x+b
!!$
!!$    ELSE IF(x>xtab(n)) THEN
!!$
!!$       x1=xtab(n-1)
!!$       x2=xtab(n)
!!$
!!$       y1=ytab(n-1)
!!$       y2=ytab(n)
!!$
!!$       CALL fit_line(a,b,x1,y1,x2,y2)
!!$       find=a*x+b
!!$
!!$    ELSE IF(iorder==1) THEN
!!$
!!$       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'
!!$
!!$       IF(x<=xtab(2)) THEN
!!$
!!$          x1=xtab(1)
!!$          x2=xtab(2)
!!$
!!$          y1=ytab(1)
!!$          y2=ytab(2)
!!$
!!$       ELSE IF (x>=xtab(n-1)) THEN
!!$
!!$          x1=xtab(n-1)
!!$          x2=xtab(n)
!!$
!!$          y1=ytab(n-1)
!!$          y2=ytab(n)
!!$
!!$       ELSE
!!$
!!$          IF(imeth==1) i=search_int(x,xtab)
!!$          IF(imeth==2) i=linear_table_integer(x,xtab)
!!$          IF(imeth==3) i=int_split(x,xtab)
!!$
!!$          x1=xtab(i)
!!$          x2=xtab(i+1)
!!$
!!$          y1=ytab(i)
!!$          y2=ytab(i+1)
!!$
!!$       END IF
!!$
!!$       CALL fit_line(a,b,x1,y1,x2,y2)
!!$       find=a*x+b
!!$
!!$    ELSE IF(iorder==2) THEN
!!$
!!$       IF(n<3) STOP 'FIND: Not enough points in your table'
!!$
!!$       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN
!!$
!!$          IF(x<=xtab(2)) THEN
!!$
!!$             x1=xtab(1)
!!$             x2=xtab(2)
!!$             x3=xtab(3)
!!$
!!$             y1=ytab(1)
!!$             y2=ytab(2)
!!$             y3=ytab(3)
!!$
!!$          ELSE IF (x>=xtab(n-1)) THEN
!!$
!!$             x1=xtab(n-2)
!!$             x2=xtab(n-1)
!!$             x3=xtab(n)
!!$
!!$             y1=ytab(n-2)
!!$             y2=ytab(n-1)
!!$             y3=ytab(n)
!!$
!!$          END IF
!!$
!!$          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
!!$
!!$          find=a*(x**2.)+b*x+c
!!$
!!$       ELSE
!!$
!!$          i=table_integer(x,xtab,imeth)
!!$          !IF(imeth==1) i=search_int(x,xtab)
!!$          !IF(imeth==2) i=linear_table_integer(x,xtab)
!!$          !IF(imeth==3) i=int_split(x,xtab)
!!$
!!$          x1=xtab(i-1)
!!$          x2=xtab(i)
!!$          x3=xtab(i+1)
!!$          x4=xtab(i+2)
!!$
!!$          y1=ytab(i-1)
!!$          y2=ytab(i)
!!$          y3=ytab(i+1)
!!$          y4=ytab(i+2)
!!$
!!$          !In this case take the average of two separate quadratic spline values
!!$
!!$          find=0.
!!$
!!$          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
!!$          find=find+(a*(x**2.)+b*x+c)/2.
!!$
!!$          CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
!!$          find=find+(a*(x**2.)+b*x+c)/2.
!!$
!!$       END IF
!!$
!!$    ELSE IF(iorder==3) THEN
!!$
!!$       IF(n<4) STOP 'FIND: Not enough points in your table'
!!$
!!$       IF(x<=xtab(3)) THEN
!!$
!!$          x1=xtab(1)
!!$          x2=xtab(2)
!!$          x3=xtab(3)
!!$          x4=xtab(4)        
!!$
!!$          y1=ytab(1)
!!$          y2=ytab(2)
!!$          y3=ytab(3)
!!$          y4=ytab(4)
!!$
!!$       ELSE IF (x>=xtab(n-2)) THEN
!!$
!!$          x1=xtab(n-3)
!!$          x2=xtab(n-2)
!!$          x3=xtab(n-1)
!!$          x4=xtab(n)
!!$
!!$          y1=ytab(n-3)
!!$          y2=ytab(n-2)
!!$          y3=ytab(n-1)
!!$          y4=ytab(n)
!!$
!!$       ELSE
!!$
!!$          IF(imeth==1) i=search_int(x,xtab)
!!$          IF(imeth==2) i=linear_table_integer(x,xtab)
!!$          IF(imeth==3) i=int_split(x,xtab)
!!$
!!$          x1=xtab(i-1)
!!$          x2=xtab(i)
!!$          x3=xtab(i+1)
!!$          x4=xtab(i+2)
!!$
!!$          y1=ytab(i-1)
!!$          y2=ytab(i)
!!$          y3=ytab(i+1)
!!$          y4=ytab(i+2)
!!$
!!$       END IF
!!$
!!$       CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
!!$       find=a*x**3.+b*x**2.+c*x+d
!!$
!!$    END IF
!!$
!!$  END FUNCTION find

  FUNCTION find(x,xin,yin,iorder,imeth)

    IMPLICIT NONE
    REAL :: find
    REAL, INTENT(IN) :: x, xin(:), yin(:)
    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i, n
    INTEGER, INTENT(IN) :: imeth, iorder
    INTEGER :: maxorder, maxmethod

    !This version interpolates if the value is off either end of the array!
    !Insert 'x, xtab' or ytab as log if this gives better results from the interpolation

    !If the value required is off the table edge the interpolation is always linear

    !imeth = 1 => find x in xtab by crudely searching from x(1) to x(n)
    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    n=SIZE(xin)
    IF(n .NE. SIZE(yin)) STOP 'FIND: Tables not of the same size'
    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       !Reverse the arrays in this case
       CALL reverse(xtab)
       CALL reverse(ytab)
    END IF

    IF(x<xtab(1)) THEN

       !Do a linear interpolation beyond the table boundary

       x1=xtab(1)
       x2=xtab(2)

       y1=ytab(1)
       y2=ytab(2)

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b
       
    ELSE IF(x>xtab(n)) THEN

       !Do a linear interpolation beyond the table boundary
       
       x1=xtab(n-1)
       x2=xtab(n)

       y1=ytab(n-1)
       y2=ytab(n)

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(iorder==1) THEN

       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          x1=xtab(1)
          x2=xtab(2)

          y1=ytab(1)
          y2=ytab(2)

       ELSE IF (x>=xtab(n-1)) THEN

          x1=xtab(n-1)
          x2=xtab(n)

          y1=ytab(n-1)
          y2=ytab(n)

       ELSE

          i=table_integer(x,xtab,imeth)
          
          x1=xtab(i)
          x2=xtab(i+1)

          y1=ytab(i)
          y2=ytab(i+1)

       END IF

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'FIND: Not enough points in your table'

       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

          IF(x<=xtab(2)) THEN

             x1=xtab(1)
             x2=xtab(2)
             x3=xtab(3)

             y1=ytab(1)
             y2=ytab(2)
             y3=ytab(3)

          ELSE IF (x>=xtab(n-1)) THEN

             x1=xtab(n-2)
             x2=xtab(n-1)
             x3=xtab(n)

             y1=ytab(n-2)
             y2=ytab(n-1)
             y3=ytab(n)

          END IF

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          find=a*(x**2.)+b*x+c

       ELSE

          i=table_integer(x,xtab,imeth)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          !In this case take the average of two separate quadratic spline values

          find=0.

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
          find=find+(a*(x**2.)+b*x+c)/2.

          CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
          find=find+(a*(x**2.)+b*x+c)/2.

       END IF

    ELSE IF(iorder==3) THEN

       IF(n<4) STOP 'FIND: Not enough points in your table'

       IF(x<=xtab(3)) THEN

          x1=xtab(1)
          x2=xtab(2)
          x3=xtab(3)
          x4=xtab(4)        

          y1=ytab(1)
          y2=ytab(2)
          y3=ytab(3)
          y4=ytab(4)

       ELSE IF (x>=xtab(n-2)) THEN

          x1=xtab(n-3)
          x2=xtab(n-2)
          x3=xtab(n-1)
          x4=xtab(n)

          y1=ytab(n-3)
          y2=ytab(n-2)
          y3=ytab(n-1)
          y4=ytab(n)

       ELSE

          i=table_integer(x,xtab,imeth)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
       find=a*x**3.+b*x**2.+c*x+d

    ELSE

       STOP 'FIND: Error, interpolation order specified incorrectly'

    END IF

  END FUNCTION find

  FUNCTION table_integer(x,xtab,imeth)

    IMPLICIT NONE
    INTEGER :: table_integer 
    REAL, INTENT(IN) :: x, xtab(:)
    INTEGER, INTENT(IN) :: imeth

    IF(imeth==1) THEN
       table_integer=linear_table_integer(x,xtab)
    ELSE IF(imeth==2) THEN
       table_integer=search_int(x,xtab)
    ELSE IF(imeth==3) THEN
       table_integer=int_split(x,xtab)
    ELSE
       STOP 'TABLE INTEGER: Method specified incorrectly'
    END IF

  END FUNCTION table_integer

  FUNCTION linear_table_integer(x,xtab)

    IMPLICIT NONE
    INTEGER :: linear_table_integer
    REAL, INTENT(IN) :: x, xtab(:)
    INTEGER :: n
    REAL :: x1, x2, xn
    REAL :: acc

    !Returns the integer (table position) below the value of x
    !eg. if x(3)=6. and x(4)=7. and x=6.5 this will return 6
    !Assumes table is organised linearly (care for logs)

    n=SIZE(xtab)
    x1=xtab(1)
    x2=xtab(2)
    xn=xtab(n)

    !Test for linear table
    acc=0.001

    IF(x1>xn) STOP 'LINEAR_TABLE_INTEGER :: table in the wrong order'
    IF(ABS(-1.+float(n-1)*(x2-x1)/(xn-x1))>acc) STOP 'LINEAR_TABLE_INTEGER :: table does not seem to be linear'

    linear_table_integer=1+FLOOR(float(n-1)*(x-x1)/(xn-x1))

  END FUNCTION linear_table_integer

  FUNCTION search_int(x,xtab)

    IMPLICIT NONE
    INTEGER :: search_int
    INTEGER :: i, n
    REAL, INTENT(IN) :: x, xtab(:)

    !Searches for the point in the table brute force.
    !This is usually a stupid thing to do

    n=SIZE(xtab)

    IF(xtab(1)>xtab(n)) STOP 'SEARCH_INT: table in wrong order'

    DO i=1,n
       IF(x>=xtab(i) .AND. x<=xtab(i+1)) EXIT
    END DO

    search_int=i

  END FUNCTION search_int

  FUNCTION int_split(x,xtab)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x, xtab(:)
    INTEGER :: i1, i2, imid, n
    INTEGER :: int_split

    !Finds the position of the value in the table by continually splitting it in half

    n=SIZE(xtab)

    IF(xtab(1)>xtab(n)) STOP 'INT_SPLIT: table in wrong order'

    i1=1
    i2=n

    DO

       imid=NINT((i1+i2)/2.)

       IF(x<xtab(imid)) THEN
          i2=imid
       ELSE
          i1=imid
       END IF

       IF(i2==i1+1) EXIT

    END DO

    int_split=i1

  END FUNCTION int_split

  SUBROUTINE fit_line(a1,a0,x1,y1,x2,y2)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1
    REAL, INTENT(IN) :: x1, y1, x2, y2

    !Given xi, yi i=1,2 fits a line between these points

    a1=(y2-y1)/(x2-x1)
    a0=y1-a1*x1

  END SUBROUTINE fit_line

  SUBROUTINE fit_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1, a2
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3

    !Given xi, yi i=1,2,3 fits a quadratic between these points

    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
    a0=y1-a2*(x1**2.)-a1*x1

  END SUBROUTINE fit_quadratic

  SUBROUTINE fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a, b, c, d
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3, x4, y4
    REAL :: f1, f2, f3

    !Given xi, yi i=1,2,3,4 fits a cubic between these points

    f1=(y4-y1)/((x4-x2)*(x4-x1)*(x4-x3))
    f2=(y3-y1)/((x3-x2)*(x3-x1)*(x4-x3))
    f3=(y2-y1)/((x2-x1)*(x4-x3))*(1./(x4-x2)-1./(x3-x2))

    a=f1-f2-f3

    f1=(y3-y1)/((x3-x2)*(x3-x1))
    f2=(y2-y1)/((x2-x1)*(x3-x2))
    f3=a*(x3+x2+x1)

    b=f1-f2-f3

    f1=(y4-y1)/(x4-x1)
    f2=a*(x4**2.+x4*x1+x1**2.)
    f3=b*(x4+x1)

    c=f1-f2-f3

    d=y1-a*x1**3.-b*x1**2.-c*x1

  END SUBROUTINE fit_cubic

  SUBROUTINE reverse(arry)

    IMPLICIT NONE
    INTEGER :: n, i
    REAL, ALLOCATABLE :: hold(:)
    REAL :: arry(:)

    !This reverses the contents of 'arry'

    n=SIZE(arry)

    ALLOCATE(hold(n))

    hold=arry

    DO i=1,n
       arry(i)=hold(n-i+1)
    END DO

    DEALLOCATE(hold)

  END SUBROUTINE reverse

  FUNCTION file_length(file_name)

    IMPLICIT NONE
    CHARACTER(len=64) :: file_name
    INTEGER ::n, file_length
    REAL :: data

    !Finds the length of a file

    OPEN(7,file=file_name)
    n=0
    DO
       n=n+1
       READ(7,*, end=301) data
    END DO

    !301 is just the label to jump to when the end of the file is reached

301 CLOSE(7)  

    file_length=n-1

  END FUNCTION file_length

  SUBROUTINE fill_growtab(cosm)

    USE cosdef
    IMPLICIT NONE
    TYPE(cosmology) :: cosm
    INTEGER :: i, n
    REAL :: a, norm
    REAL, ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
    REAL :: ainit, amax, dinit, vinit
    REAL :: acc

    !The calculation should start at a z when Om_m(z)=1., so that the assumption
    !of starting in the g\propto a growing mode is valid (this will not work for early DE)
    ainit=0.001
    !Final should be a=1. unless considering models in the future
    amax=1.

    !These set the initial conditions to be the Om_m=1. growing mode
    dinit=ainit
    vinit=1.

    !Overall accuracy for the ODE solver
    acc=0.001

    IF(ihm==1) WRITE(*,*) 'GROWTH: Solving growth equation'
    CALL ode_growth(d_tab,v_tab,a_tab,0.,ainit,amax,dinit,vinit,acc,3,cosm)
    IF(ihm==1) WRITE(*,*) 'GROWTH: ODE done'

    !Normalise so that g(z=0)=1
    norm=find(1.,a_tab,d_tab,3,3)
    IF(ihm==1) WRITE(*,*) 'GROWTH: unnormalised g(a=1):', norm
    d_tab=d_tab/norm
    IF(ihm==1) WRITE(*,*)

    !This downsamples the tables that come out of the ODE solver (which can be a bit long)
    !Could use some table-interpolation routine here to save time
    IF(ALLOCATED(cosm%a_growth)) DEALLOCATE(cosm%a_growth,cosm%growth)
    n=64
    ALLOCATE(cosm%a_growth(n),cosm%growth(n))
    DO i=1,n
       a=ainit+(amax-ainit)*float(i-1)/float(n-1)
       cosm%a_growth(i)=a
       cosm%growth(i)=find(a,a_tab,d_tab,3,3)
    END DO  

  END SUBROUTINE fill_growtab

  SUBROUTINE ode_growth(x,v,t,kk,ti,tf,xi,vi,acc,imeth,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: xi, ti, tf, dt, acc, vi, x4, v4, t4, kk
    REAL :: kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4
    REAL*8, ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)
    REAL, ALLOCATABLE :: x(:), v(:), t(:)
    INTEGER :: i, j, k, n, np, ifail, kn, imeth
!    REAL, EXTERNAL :: fx, fv
    TYPE(cosmology) :: cosm

    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values
    !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    !fx is what x' is equal to
    !fv is what v' is equal to
    !acc is the desired accuracy across the entire solution
    !imeth selects method

    IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(v)) DEALLOCATE(v)
    IF(ALLOCATED(t)) DEALLOCATE(t)

    DO j=1,30

       n=100*(2**(j-1))
       n=n+1  

       ALLOCATE(x8(n),t8(n),v8(n))

       x8=0.d0
       t8=0.d0
       v8=0.d0

       dt=(tf-ti)/float(n-1)

       x8(1)=xi
       v8(1)=vi
       t8(1)=ti

       ifail=0

       DO i=1,n-1

          x4=real(x8(i))
          v4=real(v8(i))
          t4=real(t8(i))

          IF(imeth==1) THEN

             !Crude method!
             v8(i+1)=v8(i)+fv(x4,v4,kk,t4,cosm)*dt
             x8(i+1)=x8(i)+fd(x4,v4,kk,t4,cosm)*dt
             t8(i+1)=t8(i)+dt

          ELSE IF(imeth==2) THEN

             !Mid-point method!
             kx1=dt*fd(x4,v4,kk,t4,cosm)
             kv1=dt*fv(x4,v4,kk,t4,cosm)
             kx2=dt*fd(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
             kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)

             v8(i+1)=v8(i)+kv2*dt
             x8(i+1)=x8(i)+kx2*dt
             t8(i+1)=t8(i)+dt

          ELSE IF(imeth==3) THEN

             !4th order Runge-Kutta method (fucking fast)!
             kx1=dt*fd(x4,v4,kk,t4,cosm)
             kv1=dt*fv(x4,v4,kk,t4,cosm)
             kx2=dt*fd(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
             kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
             kx3=dt*fd(x4+kx2/2.,v4+kv2/2.,kk,t4+dt/2.,cosm)
             kv3=dt*fv(x4+kx2/2.,v4+kv2/2.,kk,t4+dt/2.,cosm)
             kx4=dt*fd(x4+kx3,v4+kv3,kk,t4+dt,cosm)
             kv4=dt*fv(x4+kx3,v4+kv3,kk,t4+dt,cosm)

             x8(i+1)=x8(i)+(kx1+(2.*kx2)+(2.*kx3)+kx4)/6.
             v8(i+1)=v8(i)+(kv1+(2.*kv2)+(2.*kv3)+kv4)/6.
             t8(i+1)=t8(i)+dt

          END IF

       END DO

       IF(j==1) ifail=1

       IF(j .NE. 1) THEN

          np=1+(n-1)/2

          DO k=1,1+(n-1)/2

             kn=2*k-1

             IF(ifail==0) THEN

                IF(xh(k)>acc .AND. x8(kn)>acc .AND. (ABS(xh(k)/x8(kn))-1.)>acc) ifail=1
                IF(vh(k)>acc .AND. v8(kn)>acc .AND. (ABS(vh(k)/v8(kn))-1.)>acc) ifail=1

                IF(ifail==1) THEN
                   DEALLOCATE(xh,th,vh)
                   EXIT
                END IF

             END IF
          END DO

       END IF

       IF(ifail==0) THEN
          ALLOCATE(x(n),t(n),v(n))
          x=x8
          v=v8
          t=t8
          EXIT
       END IF

       ALLOCATE(xh(n),th(n),vh(n))
       xh=x8
       vh=v8
       th=t8
       DEALLOCATE(x8,t8,v8)

    END DO

  END SUBROUTINE ode_growth

  FUNCTION fv(d,v,k,a,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: fv
    REAL, INTENT(IN) :: d, v, k, a
    REAL :: f1, f2, z
    TYPE(cosmology), INTENT(IN) :: cosm

    !Needed for growth function solution
    !This is the fv in \ddot{\delta}=fv

    z=-1.+(1./a)

    f1=3.*omega_m(z,cosm)*d/(2.*(a**2.))
    f2=(2.+AH(z,cosm)/hubble2(z,cosm))*(v/a)

    fv=f1-f2

  END FUNCTION fv

  FUNCTION fd(d,v,k,a,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: fd
    REAL, INTENT(IN) :: d, v, k, a
    TYPE(cosmology), INTENT(IN) :: cosm
    
    !Needed for growth function solution
    !This is the fd in \dot{\delta}=fd
    
    fd=v

  END FUNCTION fd

  FUNCTION AH(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: AH
    REAL, INTENT(IN) :: z
    REAL :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    !\ddot{a}/a

    a=1./(1.+z)

    AH=cosm%om_m*(a**(-3.))+cosm%om_v*(1.+3.*w_de(a,cosm))*x_de(a,cosm)

    AH=-AH/2.

  END FUNCTION AH

  FUNCTION x_de(a,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: x_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    !The time evolution for Om_w for w(a) DE models
    x_de=(a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))

  END FUNCTION x_de

  FUNCTION w_de(a,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: w_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    !w(a) for DE models
    w_de=cosm%w+(1.-a)*cosm%wa
    
  END FUNCTION w_de

END PROGRAM HMcode


