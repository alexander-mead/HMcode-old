MODULE cosdef

  TYPE cosmology
     !Contains only things that do not need to be recalculated with each new z
     REAL :: om_m, om_b, om_v, om_c, h, n, sig8, w, wa
     REAL :: A, pkz
     REAL :: eta0, Abary
     REAL, ALLOCATABLE :: k_plin(:), plin(:)
     REAL, ALLOCATABLE :: r_sigma(:), sigma(:)
     REAL, ALLOCATABLE :: a_grow(:), grow(:)
     INTEGER :: nsig, ng, nk
     LOGICAL :: itab
  END TYPE cosmology

  TYPE tables
     !Stuff that needs to be recalculated for each new z
     REAL, ALLOCATABLE :: c(:), rv(:), nu(:), sig(:), zc(:), m(:), rr(:), sigf(:)
     REAL :: sigv, sigv100, c3, knl, rnl, neff, sig8z
     REAL :: gmin, gmax, gbmin, gbmax
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
  LOGICAL :: verbose
  REAL :: kmin, kmax, zmin, zmax, zin
  TYPE(cosmology) :: cosm
  TYPE(tables) :: lut
  CHARACTER(len=256) :: infile, outfile, redshift
  LOGICAL :: lexist, itab

  !Constants
  REAL, PARAMETER :: pi=3.141592654

  !Accuracy parameter
  REAL, PARAMETER :: acc=1e-4

  !ihm parameter
  !1 - Do the accurate calculation detailed in Mead et al. (2015; 1505.07833) with updates from Mead et al. (2016; 1602.02154)
  !2 - Standard halo-model calculation (Dv=200, dc=1.686) with linear two-halo term'
  !3 - Standard halo-model calculation (Dv=200, dc=1.686) with full two-halo term'
  INTEGER, PARAMETER :: ihm=1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !2016/09/19 - Changed subroutines so that none assume an array size
  !2017/06/01 - Added comments, made recompatible with ifort (thanks Dipak Munshi)
  !2017/06/15 - Removed bug in E&H (1997) fitting function that created small spikes in linear P(k) (thanks David Copeland)
  !2017/08/?? - Increased integration routine speed
  !2017/09/27 - Added baryon models explicitly
  !2018/01/18 - Added capacity for input CAMB linear P(k)
  !2018/02/04 - Added two-halo bias integral

  !HMcode developed by Alexander Mead
  !If you use this in your work please cite the original paper: http://arxiv.org/abs/1505.07833
  !Also maybe cite the update: http://arxiv.org/abs/1602.02154
  !Also consider citing the source code at ASCL: http://ascl.net/1508.001

  !verbose variable (cannot be parameter because value is changed in code)
  verbose=.TRUE.

  CALL get_command_argument(1,infile)
  IF(infile=='') THEN
     itab=.FALSE.
  ELSE
     INQUIRE(file=infile,exist=lexist)
     IF(lexist .EQV. .FALSE.) STOP 'This input file does not exist'
     itab=.TRUE.
     CALL get_command_argument(2,redshift)
     IF(redshift=='') STOP 'You need to supply a redshift for the input P(k)'
     READ(redshift,*) zin
  END IF
 
  WRITE(*,*)
  WRITE(*,*) 'Welcome to HMcode'
  WRITE(*,*) '================='
  WRITE(*,*)
  IF(ihm==1) THEN
     WRITE(*,*) 'HMcode: Doing accurate calculation'
  ELSE IF(ihm==2) THEN
     WRITE(*,*) 'HMcode: Doing standard calculation with linear two-halo term'
  ELSE IF(ihm==3) THEN
     WRITE(*,*) 'HMcode: Doing standard calculation with full two-halo term'
  ELSE
     STOP 'HMcode: Error, ihm specified incorrectly'
  END IF
  WRITE(*,*)

  !Set number of k points and k range (log spaced)
  nk=200
  kmin=1e-3
  kmax=1e4
  CALL fill_table(log(kmin),log(kmax),k,nk)
  k=exp(k)

  WRITE(*,*) 'HMcode: k min:', kmin
  WRITE(*,*) 'HMcode: k max:', kmax
  WRITE(*,*) 'HMcode: number of k:', nk
  WRITE(*,*)

  !Set the number of redshifts and range (linearly spaced)
  nz=16
  zmin=0.
  zmax=4.
  CALL fill_table(zmin,zmax,ztab,nz)

  WRITE(*,*) 'HMcode: z min:', zmin
  WRITE(*,*) 'HMcode: z max:', zmax
  WRITE(*,*) 'HMcode: number of z:', nz
  WRITE(*,*)

  !Fill table for output power
  ALLOCATE(ptab(nz,nk))

  !Assigns the cosmological model
  CALL assign_cosmology(cosm)

  !Read in linear P(k) if provided
  IF(itab) THEN
     CALL read_CAMB_Pk(cosm%k_plin,cosm%plin,cosm%nk,infile)
     cosm%itab=.TRUE. !Change this variable to true
     cosm%pkz=zin
  END IF

  !Normalises power spectrum (via sigma_8) and fills sigma(R) look-up tables
  CALL initialise_cosmology(cosm)

!!$  !Ignore this, useful only for bug tests
!!$  CALL RNG_set(0)
!!$  DO
!!$  CALL random_cosmology(cosm)

  CALL write_cosmology(cosm)

  !Loop over redshifts
  DO j=1,nz

     !Sets the current redshift from the table
     z=ztab(j)

     !Initiliasation for the halomodel calcualtion
     CALL halomod_init(z,lut,cosm)

     !Loop over k values
     DO i=1,nk

        plin=p_lin(k(i),z,cosm)

        CALL halomod(k(i),z,p1h,p2h,pfull,plin,lut,cosm)

        ptab(j,i)=pfull

     END DO

     IF(j==1) THEN
        WRITE(*,fmt='(A5,A7)') 'i', 'z'
        WRITE(*,fmt='(A13)') '   ============'
     END IF
     WRITE(*,fmt='(I5,F8.3)') j, ztab(j)

  END DO
  WRITE(*,fmt='(A13)') '   ============'
  WRITE(*,*)

  outfile='power.dat'
  WRITE(*,*) 'HMcode: Writing output to: ', TRIM(outfile)
  WRITE(*,*) 'HMcode: The first entry in the file is hashes - #####'
  WRITE(*,*) 'HMcode: The remainder of the top row contains the redshifts'
  WRITE(*,*) 'HMcode: The first column containts ''k'' after the hashes'
  WRITE(*,*) 'HMcode: The rows then contain the power at that ''k'' for each redshift'
  OPEN(7,file=outfile)
  DO i=0,nk
     IF(i==0) THEN
        WRITE(7,fmt='(A20,40F20.10)') '#####', (ztab(j), j=1,nz)
     ELSE
        WRITE(7,fmt='(F20.10,40F20.10)') k(i), (ptab(j,i), j=1,nz)
     END IF
  END DO
  CLOSE(7)
  WRITE(*,*) 'HMcode: Done'
  WRITE(*,*)

!!$  !Ignore this, only useful for bug tests
!!$  END DO

CONTAINS

  SUBROUTINE read_CAMB_Pk(k,p,n,infile)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: infile
    REAL, ALLOCATABLE, INTENT(OUT) :: k(:), p(:)
    INTEGER, INTENT(OUT) :: n
    INTEGER :: i
    
    n=file_length(infile)
    n=n-1
    WRITE(*,*) 'READ_CAMB_PK: CAMB file: ', TRIM(infile)
    WRITE(*,*) 'READ_CAMB_PK: Number of points:', n
    WRITE(*,*)

    ALLOCATE(k(n),p(n))

    OPEN(7,file=infile)
    DO i=0,n
       IF(i==0) THEN
          READ(7,*)
       ELSE
          READ(7,*) k(i), p(i)
       END IF
    END DO
    CLOSE(7)

    !Convert from P(k) -> Delta^2(k)
    p=4.*pi*p*(k**3)/(2.*pi)**3

  END SUBROUTINE read_CAMB_Pk

  FUNCTION Delta_v(z,cosm)

    !Virialised overdensity
    IMPLICIT NONE
    REAL :: Delta_v
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(ihm==2 .OR. ihm==3) THEN
       Delta_v=200.
    ELSE IF(ihm==1) THEN
       Delta_v=418.*(omega_m(z,cosm)**(-0.352))
    ELSE
       STOP 'Error, ihm defined incorrectly'
    END IF

  END FUNCTION Delta_v

  FUNCTION delta_c(z,cosm)

    !Linear collapse density
    IMPLICIT NONE
    REAL :: delta_c
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(ihm==2 .OR. ihm==3) THEN
       delta_c=1.686
    ELSE IF(ihm==1) THEN
       delta_c=1.59+0.0314*log(sigma_cb(8.,z,cosm))
    ELSE
       STOP 'Error, ihm defined incorrectly'
    END IF

    !Nakamura & Suto (1997) fitting formula for LCDM
    delta_c=delta_c*(1.+0.0123*log10(omega_m(z,cosm)))

  END FUNCTION delta_c

  FUNCTION eta(z,cosm)

    IMPLICIT NONE
    REAL :: eta
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(ihm==2 .OR. ihm==3) THEN
       eta=0.
    ELSE IF(ihm==1) THEN
       !The first parameter here is 'eta_0' in Mead et al. (2015; arXiv 1505.07833)
       !eta=0.603-0.3*(sigma_cb(8.,z,cosm))
       eta=cosm%eta0-0.3*(sigma_cb(8.,z,cosm))
    ELSE
       STOP 'Error, ihm defined incorrectly'
    END IF

  END FUNCTION eta

  FUNCTION kstar(lut)

    IMPLICIT NONE
    REAL :: kstar
    !TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut

    IF(ihm==2 .OR. ihm==3) THEN
       !Set to zero for the standard Poisson one-halo term
       kstar=0.
    ELSE IF(ihm==1) THEN
       !One-halo cut-off wavenumber
       kstar=0.584*(lut%sigv)**(-1.)
    ELSE
       STOP 'Error, ihm defined incorrectly'
    END IF

  END FUNCTION kstar

  FUNCTION As(cosm)

    IMPLICIT NONE
    REAL :: As
    TYPE(cosmology), INTENT(IN) :: cosm

    !Halo concentration pre-factor
    IF(ihm==2 .OR. ihm==3) THEN
       !Set to 4 for the standard Bullock value
       As=4.
    ELSE IF(ihm==1) THEN
       !This is the 'A' halo-concentration parameter in Mead et al. (2015; arXiv 1505.07833)
       As=cosm%Abary
    ELSE
       STOP 'Error, ihm defined incorrectly'
    END IF

  END FUNCTION As

  FUNCTION fdamp(lut)

    IMPLICIT NONE
    REAL ::fdamp
    !REAL, INTENT(IN) :: z
    !TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut

    !Linear theory damping factor
    IF(ihm==2 .OR. ihm==3) THEN
       !Set to 0 for the standard linear theory two halo term
       fdamp=0.
    ELSE IF(ihm==1) THEN
       !fdamp=0.188*sigma_cb(8.,z,cosm)**4.29
       fdamp=0.0095*lut%sigv100**1.37
       !Catches extreme values of fdamp that occur for ridiculous cosmologies
       IF(fdamp<1.e-3) fdamp=0.
       IF(fdamp>0.99)  fdamp=0.99
    ELSE
       STOP 'Error, ihm defined incorrectly'
    END IF

  END FUNCTION fdamp

  FUNCTION alpha(lut)

    IMPLICIT NONE
    REAL :: alpha
    TYPE(tables), INTENT(IN) :: lut

    IF(ihm==2 .OR. ihm==3) THEN
       !Set to 1 for the standard halo model addition of one- and two-halo terms
       alpha=1.
    ELSE IF(ihm==1) THEN
       !This uses the top-hat defined neff
       alpha=3.24*1.85**lut%neff
    ELSE
       STOP 'Error, ihm defined incorrectly'
    END IF

    !Catches values of alpha that are crazy
    IF(alpha>2.)  alpha=2.
    IF(alpha<0.5) alpha=0.5

  END FUNCTION alpha

  SUBROUTINE write_parameters(z,lut,cosm)

    !This subroutine writes out the physical parameters at some redshift 
    !(e.g. Delta_v) rather than the model parameters
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut    

    WRITE(*,*) 'HMcode: Parameters'
    WRITE(*,*) '==================='
    WRITE(*,fmt='(A10,F10.5)') 'z:', z
    WRITE(*,fmt='(A10,F10.5)') 'Dv:', Delta_v(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'dc:', delta_c(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'eta:', eta(z,cosm)
    WRITE(*,fmt='(A10,F10.5)') 'k*:', kstar(lut)
    WRITE(*,fmt='(A10,F10.5)') 'A:', As(cosm)
    WRITE(*,fmt='(A10,F10.5)') 'fdamp:', fdamp(lut)
    WRITE(*,fmt='(A10,F10.5)') 'alpha:', alpha(lut)
    WRITE(*,*) '==================='
    WRITE(*,*)

  END SUBROUTINE write_parameters

  FUNCTION r_nl(lut)

    !Calculates k_nl as 1/R where nu(R)=1.
    IMPLICIT NONE
    TYPE(tables), INTENT(IN) :: lut
    REAL :: r_nl    

    IF(lut%nu(1)>1.) THEN
       !This catches some very strange values
       r_nl=lut%rr(1)
    ELSE
       r_nl=exp(find(log(1.),log(lut%nu),log(lut%rr),lut%n,3,3,2))
    END IF

  END FUNCTION r_nl

  SUBROUTINE halomod(k,z,p1h,p2h,pfull,plin,lut,cosm)

    !Calls expressions for one- and two-halo terms and then combines
    !to form the full power spectrum
    IMPLICIT NONE
    REAL, INTENT(OUT) :: p1h, p2h, pfull
    REAL, INTENT(IN) :: plin, k, z
    TYPE(cosmology), INTENT(IN) :: cosm
    TYPE(tables), INTENT(IN) :: lut
    REAL :: alp, rv, c, nu, et
    REAL :: wk(lut%n)
    INTEGER :: i

    IF(k==0.) THEN
       p1h=0.
       p2h=0.
    ELSE

       !Only call eta once
       et=eta(z,cosm)

       DO i=1,lut%n
          nu=lut%nu(i)
          rv=lut%rv(i)
          c=lut%c(i)
          wk(i)=win(k*nu**et,rv,c)
       END DO
       
       p1h=p_1h(wk,k,lut,cosm)
       p2h=p_2h(wk,k,plin,lut,cosm)
       
    END IF

    alp=alpha(lut)
    pfull=(p2h**alp+p1h**alp)**(1./alp)

  END SUBROUTINE halomod

  SUBROUTINE fill_table(min,max,arr,n)

    !Fills array 'arr' in equally spaced intervals
    !I'm not sure if inputting an array like this is okay
    IMPLICIT NONE
    INTEGER :: i
    REAL, INTENT(IN) :: min, max
    REAL, ALLOCATABLE :: arr(:)
    INTEGER, INTENT(IN) :: n

    !Allocate the array, and deallocate it if it is full
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))
    arr=0.

    IF(n==1) THEN
       arr(1)=min
    ELSE IF(n>1) THEN
       DO i=1,n
          arr(i)=min+(max-min)*float(i-1)/float(n-1)
       END DO
    END IF

  END SUBROUTINE fill_table

  SUBROUTINE fill_table8(min,max,arr,n)

    !Fills array 'arr' in equally spaced intervals
    !I'm not sure if inputting an array like this is okay
    IMPLICIT NONE
    INTEGER :: i
    REAL*8, INTENT(IN) :: min, max
    REAL*8, ALLOCATABLE :: arr(:)
    INTEGER, INTENT(IN) :: n

    !Allocate the array, and deallocate it if it is full
    IF(ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(n))
    arr=0.

    IF(n==1) THEN
       arr(1)=min
    ELSE IF(n>1) THEN
       DO i=1,n
          arr(i)=min+(max-min)*float(i-1)/float(n-1)
       END DO
    END IF

  END SUBROUTINE fill_table8

  SUBROUTINE write_cosmology(cosm)

    !Writes the cosmological parameters to the screen
    IMPLICIT NONE
    TYPE(cosmology) :: cosm

    IF(verbose) THEN
       WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_m:', cosm%om_m
       WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_b:', cosm%om_b
       WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_c:', cosm%om_c
       WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'Omega_v:', cosm%om_v
       WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'h:', cosm%h
       WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'w_0:', cosm%w
       WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'w_a:', cosm%wa
       WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'sig8:', cosm%sig8
       WRITE(*,fmt='(A11,A10,F10.5)') 'COSMOLOGY:', 'n:', cosm%n   
       WRITE(*,*)
    END IF

  END SUBROUTINE write_cosmology

  SUBROUTINE assign_cosmology(cosm)

    !Assigns the cosmological parameters - edit this to change cosmology
    IMPLICIT NONE
    !LOGICAL, INTENT(IN) :: itab
    TYPE(cosmology), INTENT(OUT) :: cosm

    !Read in linear P(k) or not; set to false by default
    cosm%itab=.FALSE.

    !Standard cosmological parameters
    cosm%om_m=0.3
    cosm%om_v=1.-cosm%om_m
    cosm%om_b=0.05
    cosm%om_c=cosm%om_m-cosm%om_b
    cosm%h=0.7
    cosm%w=-1.
    cosm%sig8=0.8
    cosm%n=0.96
    cosm%wa=0.

    !Baryon feedback parameters
    cosm%Abary=3.13 !DMONLY
    !cosm%Abary=3.12 !REF
    !cosm%Abary=2.41 !DBLIM
    !cosm%Abary=2.11 !AGN
    !cosm%Abary=2.75

    !Enfore one-parameter baryon model (see Mead et al. 2015)
    cosm%eta0=0.98-0.12*cosm%Abary

  END SUBROUTINE assign_cosmology

  SUBROUTINE initialise_cosmology(cosm)

    !Does some initialisation steps that are necessary for each cosmological model
    !1 - Normalises the power spectrum to have correct sigma8
    !2 - Fills look-up tables of R vs. Sigma(R)
    IMPLICIT NONE
    REAL :: sigi
    TYPE(cosmology) :: cosm

    CALL fill_growtab(cosm)

    IF(cosm%itab) THEN

       !'Grow' the input power to z=0
       IF(verbose) WRITE(*,*) 'INITIALISE: Growing input P(k) to z=0'
       cosm%plin=cosm%plin*(grow(0.,cosm)/grow(cosm%pkz,cosm))**2
       cosm%plin=log(cosm%plin)
       cosm%k_plin=log(cosm%k_plin)
       sigi=sigma(8.,0.,cosm)
       cosm%sig8=sigi
       IF(verbose) WRITE(*,*) 'INITIALISE: sigma_8:', sigi
       IF(verbose) WRITE(*,*) 'INITIALISE: Done'
       IF(verbose) WRITE(*,*)

    ELSE

       !Initially set amplitude: A=1
       cosm%A=1.

       !Calculate sigma_8 under the assumption of A=1
       sigi=sigma(8.,0.,cosm)

       IF(verbose) WRITE(*,*) 'INITIALISE: Initial sigma_8:', sigi

       !Change the amplitude to ensure that the power will be normalised correctly at z=0
       cosm%A=cosm%sig8/sigi    

       !Check that the normalisation is correct
       sigi=sigma(8.,0.,cosm)

       IF(verbose) THEN
          WRITE(*,*) 'INITIALISE: Normalisation factor:', cosm%A
          WRITE(*,*) 'INITIALISE: Target sigma_8:', cosm%sig8
          WRITE(*,*) 'INITIALISE: Final sigma_8 (calculated):', sigi
          WRITE(*,*) 'INITIALISE: Complete'
          WRITE(*,*)
       END IF 

    END IF

    !Fill tables of r vs. sigma(r)
    CALL fill_sigtab(cosm)

  END SUBROUTINE initialise_cosmology

!!$  SUBROUTINE random_cosmology(cosm)
!!$
!!$    !Makes a 'random' cosmological model - good for testing
!!$    IMPLICIT NONE
!!$    TYPE(cosmology) :: cosm
!!$    REAL :: om_m_min, om_m_max, om_b_min, om_b_max, n_min, n_max
!!$    REAL :: w_min, w_max, h_min, h_max, sig8_min, sig8_max, wa_min, wa_max
!!$
!!$    !Needs to be set to normalise P_lin
!!$    cosm%A=1.
!!$
!!$    om_m_min=0.1
!!$    om_m_max=1.
!!$    cosm%om_m=uniform(om_m_min,om_m_max)
!!$
!!$    cosm%om_v=1.-cosm%om_m
!!$
!!$    om_b_min=0.005
!!$    om_b_max=MIN(0.095,cosm%om_m)
!!$    cosm%om_b=uniform(om_b_min,om_b_max)
!!$
!!$    cosm%om_c=cosm%om_m-cosm%om_b
!!$
!!$    n_min=0.5
!!$    n_max=1.5
!!$    cosm%n=uniform(n_min,n_max)
!!$
!!$    h_min=0.4
!!$    h_max=1.2
!!$    cosm%h=uniform(h_min,h_max)
!!$
!!$    w_min=-1.5
!!$    w_max=-0.5
!!$    cosm%w=uniform(w_min,w_max)
!!$
!!$    wa_min=-1.
!!$    wa_max=-cosm%w*0.8
!!$    cosm%wa=uniform(wa_min,wa_max)
!!$
!!$    sig8_min=0.2
!!$    sig8_max=1.5
!!$    cosm%sig8=uniform(sig8_min,sig8_max)
!!$
!!$  END SUBROUTINE random_cosmology
!!$
!!$  SUBROUTINE RNG_set(seed)
!!$
!!$    !Seeds the RNG using the system clock so that it is different each time
!!$    IMPLICIT NONE
!!$    INTEGER :: int, timearray(3)
!!$    REAL :: rand !Needs to be defined for ifort (thanks Dipak Munshi)
!!$    INTEGER, INTENT(IN) :: seed
!!$
!!$    WRITE(*,*) 'Initialising RNG'
!!$
!!$    IF(seed==0) THEN
!!$       !This fills the time array using the system clock!
!!$       !If called within the same second the numbers will be identical!
!!$       CALL itime(timeArray)
!!$       !This then initialises the generator!
!!$       int=FLOOR(rand(timeArray(1)+timeArray(2)+timeArray(3)))
!!$    ELSE
!!$       int=FLOOR(rand(seed))
!!$    END IF
!!$    WRITE(*,*) 'RNG set'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE RNG_set

  FUNCTION uniform(x1,x2)

    !Produces a uniform random number between x1 and x2
    IMPLICIT NONE
    REAL :: uniform
    REAL :: rand !Needs to be defined for ifort (thanks Dipak Munshi)
    REAL, INTENT(IN) :: x1,x2

    !Rand is some inbuilt function
    uniform=x1+(x2-x1)*(rand(0))

  END FUNCTION uniform

  SUBROUTINE allocate_LUT(lut)

    !Allocates memory for the look-up tables
    IMPLICIT NONE
    TYPE(tables) :: lut
    INTEGER :: n

    !The number of points in the look-up tables
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

    !Deallocates memory for the look-up tables
    IMPLICIT NONE
    TYPE(tables) :: lut
    
    DEALLOCATE(lut%zc,lut%m,lut%c,lut%rv,lut%nu,lut%rr,lut%sigf,lut%sig)

  END SUBROUTINE deallocate_LUT

  SUBROUTINE halomod_init(z,lut,cosm)

    !Halo-model initialisation routine
    !The computes other tables necessary for the one-halo integral
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    INTEGER :: i
    REAL :: Dv, dc, f, m, nu, r, sig
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut
    REAL, PARAMETER :: mmin=1e2 !Minimum mass for integration
    REAL, PARAMETER :: mmax=1e18 !Maximum mass for integration
    INTEGER, PARAMETER :: n=256 !Number of points for integration
    REAL, PARAMETER :: large_nu=10. !A large value of nu(M)
    
    !Find value of sigma_v
    lut%sigv=sqrt(dispint(0.,z,cosm)/3.)
    lut%sigv100=sqrt(dispint(100.,z,cosm)/3.)
    lut%sig8z=sigma(8.,z,cosm)

    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD: Filling look-up tables'
       WRITE(*,*) 'HALOMOD: Tables being filled at redshift:', z
       WRITE(*,*) 'HALOMOD: sigv [Mpc/h]:', lut%sigv
       WRITE(*,*) 'HALOMOD: sigv100 [Mpc/h]:', lut%sigv100
       WRITE(*,*) 'HALOMOD: sig8(z):', lut%sig8z
    END IF

    IF(ALLOCATED(lut%rr)) CALL deallocate_LUT(lut)

    lut%n=n
    CALL allocate_lut(lut)

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

    IF(verbose) WRITE(*,*) 'HALOMOD: m, r, nu, sig tables filled'

    !Fills up a table for sigma(fM) for Bullock c(m) relation
    !This is the f=0.01 parameter in the Bullock realtion sigma(fM,z)
    f=0.01**(1./3.)
    DO i=1,lut%n
       lut%sigf(i)=sigma_cb(lut%rr(i)*f,z,cosm)
    END DO
    IF(verbose) WRITE(*,*) 'HALOMOD: sigf tables filled'  

    !Fill virial radius table using real radius table
    Dv=Delta_v(z,cosm)
    lut%rv=lut%rr/(Dv**(1./3.))

    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD: rv tables filled'  
       WRITE(*,*) 'HALOMOD: nu min:', lut%nu(1)
       WRITE(*,*) 'HALOMOD: nu max:', lut%nu(lut%n)
       WRITE(*,*) 'HALOMOD: R_v min [Mpc/h]:', lut%rv(1)
       WRITE(*,*) 'HALOMOD: R_v max [Mpc/h]:', lut%rv(lut%n)
       WRITE(*,*) 'HALOMOD: M min [Msun/h]:', lut%m(1)
       WRITE(*,*) 'HALOMOD: M max [Msun/h]:', lut%m(lut%n)
    END IF

    lut%gmin=1.-integrate(lut%nu(1),large_nu,gnu,acc,3)
    lut%gmax=integrate(lut%nu(lut%n),large_nu,gnu,acc,3)
    lut%gbmin=1.-integrate(lut%nu(1),large_nu,gnubnu,acc,3)
    lut%gbmax=integrate(lut%nu(lut%n),large_nu,gnubnu,acc,3)
    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD: missing g(nu) at low end:', REAL(lut%gmin)
       WRITE(*,*) 'HALOMOD: missing g(nu) at high end:', REAL(lut%gmax)
       WRITE(*,*) 'HALOMOD: missing g(nu)b(nu) at low end:', REAL(lut%gbmin)
       WRITE(*,*) 'HALOMOD: missing g(nu)b(nu) at high end:', REAL(lut%gbmax)
    END IF

    !Find non-linear radius and scale
    lut%rnl=r_nl(lut)
    lut%knl=1./lut%rnl

    IF(verbose) WRITE(*,*) 'HALOMOD: r_nl [Mpc/h]:', lut%rnl
    IF(verbose) WRITE(*,*) 'HALOMOD: k_nl [h/Mpc]:', lut%knl

    lut%neff=neff(lut,cosm)

    IF(verbose) WRITE(*,*) 'HALOMOD: n_eff:', lut%neff

    CALL conc_bull(z,cosm,lut)

    IF(verbose) THEN
       WRITE(*,*) 'HALOMOD: c tables filled'
       WRITE(*,*) 'HALOMOD: c min:', lut%c(lut%n)
       WRITE(*,*) 'HALOMOD: c max:', lut%c(1)
       WRITE(*,*) 'HALOMOD: Done'
       WRITE(*,*)
       CALL write_parameters(z,lut,cosm)
    END IF

    verbose=.FALSE.

  END SUBROUTINE halomod_init

  PURE FUNCTION radius_m(m,cosm)

    !Finds the comoving radius that encloses mass 'M' in the homogeneous universe
    IMPLICIT NONE
    REAL :: radius_m
    REAL, INTENT(IN) :: m
    TYPE(cosmology), INTENT(IN) :: cosm
    !REAL, PARAMETER :: pi=3.141592654

    radius_m=(3.*m/(4.*pi*cosmic_density(cosm)))**(1./3.)

  END FUNCTION radius_m

  FUNCTION neff(lut,cosm)

    !Calculates the effective spectral index at the non-linear scale
    IMPLICIT NONE
    REAL :: neff
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut

    !Numerical differentiation to find effective index at collapse
    neff=-3.-derivative_table(log(lut%rnl),log(lut%rr),log(lut%sig**2.),lut%n,3,3)

    !For some bizarre cosmologies r_nl is very small, so almost no collapse has occured
    !In this case the n_eff calculation goes mad and needs to be fixed using this fudge.
    IF(neff<cosm%n-4.) neff=cosm%n-4.
    IF(neff>cosm%n)    neff=cosm%n

  END FUNCTION neff

  SUBROUTINE conc_bull(z,cosm,lut)

    !Calculates the Bullock et al. (2001) mass-concentration relation
    !Note that there are two c(M) relations in that paper and this routine computes the more complex one
    !i.e., not the one that is a simple power law in mass and z
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology) :: cosm, cos_lcdm
    TYPE(tables) :: lut
    REAL :: A, zinf, ainf, zf, g_lcdm, g_wcdm
    INTEGER :: i

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
       DEALLOCATE(cos_lcdm%grow)
       DEALLOCATE(cos_lcdm%a_grow)
       cos_lcdm%w=-1.
       cos_lcdm%wa=0.
       cos_lcdm%om_v=1.-cosm%om_m !Added this so that 'making a LCDM cosmology' works for curved models.

       ainf=1./(1.+zinf)

       !Needs to use grow_int explicitly in case tabulated values are stored
       !g_lcdm=grow_int(ainf,cos_lcdm)
       g_lcdm=growint(ainf,cos_lcdm)

       !Changed this to a power of 1.5, which produces more accurate results for extreme DE
       lut%c(i)=lut%c(i)*((g_wcdm/g_lcdm)**1.5)

       !END IF

    END DO

  END SUBROUTINE conc_bull

  SUBROUTINE zcoll_bull(z,cosm,lut)

    !This fills up the halo collapse redshift table as per Bullock relations     
    IMPLICIT NONE
    REAL, INTENT(IN) :: z
    TYPE(cosmology) :: cosm
    TYPE(tables) :: lut
    REAL :: dc
    REAL :: af, zf, RHS, a, growz
    REAL, ALLOCATABLE :: af_tab(:), grow_tab(:)
    INTEGER :: i, ntab      

    ntab=SIZE(cosm%grow)
    ALLOCATE(af_tab(ntab),grow_tab(ntab))

    af_tab=cosm%a_grow
    grow_tab=cosm%grow

    !Do numerical inversion
    DO i=1,lut%n

       !I don't think this is really consistent with dc varying as a function of z
       !But the change will be very small
       dc=delta_c(z,cosm)

       RHS=dc*grow(z,cosm)/lut%sigf(i)

       a=1./(1.+z)
       growz=find(a,af_tab,grow_tab,cosm%ng,3,3,2)

       IF(RHS>growz) THEN
          zf=z
       ELSE
          af=find(RHS,grow_tab,af_tab,cosm%ng,3,3,2)
          zf=-1.+1./af
       END IF

       lut%zc(i)=zf

    END DO

    DEALLOCATE(af_tab,grow_tab)

  END SUBROUTINE zcoll_bull

  FUNCTION mass_r(r,cosm)

    !Calculates the mass enclosed by comoving radius 'r' in the homogeneous universe
    IMPLICIT NONE
    REAL :: mass_r, r
    TYPE(cosmology) :: cosm

    !Relation between mean cosmological mass and radius
    mass_r=(4.*pi/3.)*cosmic_density(cosm)*(r**3.)

  END FUNCTION mass_r

  PURE FUNCTION cosmic_density(cosm)

    !Calculates the comoving cosmological mass density (independent of time obviously)
    IMPLICIT NONE
    REAL :: cosmic_density
    TYPE(cosmology), INTENT(IN) :: cosm

    !In Msun per Mpc^3 with h factors included. The constant does this.
    cosmic_density=(2.775e11)*cosm%om_m

  END FUNCTION cosmic_density

  FUNCTION Tk(k,cosm)

    !Get the transfer function
    IMPLICIT NONE
    REAL :: Tk, k
    TYPE(cosmology) :: cosm

    !Transfer function
    Tk=Tk_eh(k,cosm)

  END FUNCTION Tk

  FUNCTION Tk_eh(yy,cosm)

    !This routine was originally written by John Peacock
    !The astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])
    !remember I use k/h, whereas they use pure k, om_m is cdm + baryons
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
    rd=31500.*om_b*h*h/thet**4./zd !Should this be 1+zd (Steven Murray enquirey)?
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
    !IF((rk/rks**1.4)>10.) THEN
    !   fac=0.
    !ELSE
    !Removed this IF statement as it produced a discontinuity in P_lin(k) as cosmology
    !was varied - thanks David Copeland for pointing this out
    fac=exp(-(rk/rks)**1.4)
    !END IF
    tb=(tb+ab*fac/(1.+(bb/rk/s)**3.))*sin(rk*ss)/rk/ss

    tk_eh=real((om_b/om_m)*tb+(1-om_b/om_m)*tc)

  END FUNCTION TK_EH

  FUNCTION p_lin(k,z,cosm)

    !This gives the linear power spectrum for the model in question
    !P(k) should have been previously normalised so as to get the amplitude 'A' correct
    IMPLICIT NONE
    REAL :: p_lin
    REAL, INTENT (IN) :: k, z
    TYPE(cosmology), INTENT(IN) :: cosm

    IF(k==0.) THEN
       !If p_lin happens to be foolishly called for 0 mode (which should never happen, but might in integrals)
       p_lin=0.
    ELSE IF(k>1.e8) THEN
       !Avoids some issues if p_lin is called for very (absurdly) high k values
       !For some reason crashes can occur if this is the case
       p_lin=0.
    ELSE IF(cosm%itab) THEN
       !Get the linear power from the table, using interpolation
       p_lin=(grow(z,cosm)**2)*exp(find(log(k),cosm%k_plin,cosm%plin,cosm%nk,3,3,2))
    ELSE
       !In this case look for the transfer function
       p_lin=(cosm%A**2)*(grow(z,cosm)**2)*(Tk(k,cosm)**2)*(k**(cosm%n+3.))
    END IF

  END FUNCTION p_lin

  FUNCTION p_2h(wk,k,plin,lut,cosm)

    !Produces the 'two-halo' power
    IMPLICIT NONE
    REAL :: p_2h
    REAL, INTENT(IN) :: k, plin!, wk(lut%n)
    TYPE(tables), INTENT(IN) :: lut
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, INTENT(IN) :: wk(lut%n)
    REAL :: sigv, frac!, rhom
    REAL, ALLOCATABLE :: integrand(:)
    REAL :: nu, m, m0, w
    REAL :: sum, crap
    INTEGER :: i

    INTEGER, PARAMETER :: ip2h=2

    !rhom=comoving_matter_density(cosm)

    !Stop compile-time warnings
    crap=cosm%A

    IF(ihm==3) THEN

       ALLOCATE(integrand(lut%n))

       DO i=1,lut%n

          m=lut%m(i)
          nu=lut%nu(i)
          w=wk(i)

          !Linear bias term
          integrand(i)=gnu(nu)*bnu(nu)*w

       END DO

       !Evaluate these integrals from the tabled values
       sum=integrate_table(lut%nu,integrand,lut%n,1,lut%n,3)

       IF(ip2h==0) THEN
          !Do nothing in this case
       ELSE IF(ip2h==1) THEN
          !Add on the value of integral b(nu)*g(nu) assuming w=1
          sum=sum+lut%gbmin
       ELSE IF(ip2h==2) THEN
          !Put the missing part of the integrand as a delta function at nu1
          m0=lut%m(1)
          w=wk(1)
          sum=sum+lut%gbmin*w
       ELSE
          STOP 'P_2h: Error, ip2h not specified correctly'
       END IF

       p_2h=plin*sum**2

    ELSE IF(ihm==1 .OR. ihm==2) THEN

       sigv=lut%sigv
       frac=fdamp(lut)

       IF(frac==0.) THEN
          p_2h=plin
       ELSE
          p_2h=plin*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2)
       END IF

       !For some extreme cosmologies frac>1. so this must be added to prevent p_2h<0.
       IF(p_2h<0.) p_2h=0.

    ELSE

       STOP 'P_2h: Error, ihm specified incorrectly'

    END IF

  END FUNCTION p_2h

!!$  FUNCTION p_2h(k,z,plin,lut,cosm)
!!$
!!$    !Produces the 'two-halo' power
!!$    IMPLICIT NONE
!!$    REAL :: p_2h
!!$    REAL, INTENT(IN) :: k, plin
!!$    REAL, INTENT(IN) :: z
!!$    TYPE(tables), INTENT(IN) :: lut
!!$    TYPE(cosmology), INTENT(IN) :: cosm
!!$    REAL :: sigv, frac
!!$    REAL :: crap
!!$
!!$    !Stop compile-time warnings
!!$    crap=z
!!$    crap=cosm%A
!!$
!!$    sigv=lut%sigv
!!$    frac=fdamp(lut)
!!$
!!$    IF(frac==0.) THEN
!!$       p_2h=plin
!!$    ELSE
!!$       p_2h=plin*(1.-frac*(tanh(k*sigv/sqrt(ABS(frac))))**2)
!!$    END IF
!!$
!!$    !For some extreme cosmologies frac>1. so this must be added to prevent p_2h<0.
!!$    IF(p_2h<0.) p_2h=0.
!!$
!!$  END FUNCTION p_2h

  FUNCTION p_1h(wk,k,lut,cosm)

    !Does the one-halo power integral
    IMPLICIT NONE
    REAL :: p_1h
    REAL, INTENT(IN) :: k
    TYPE(tables), INTENT(IN) :: lut
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, INTENT(IN) :: wk(lut%n)
    REAL :: g, fac, ks, w, m
    REAL, ALLOCATABLE :: integrand(:)
    INTEGER :: i

    ALLOCATE(integrand(lut%n))

    !Only call eta once
    !et=eta(z,cosm)

    !Calculates the value of the integrand at all nu values!
    DO i=1,lut%n
       g=gnu(lut%nu(i))
       m=lut%m(i)
       !wk=win(k*(lut%nu(i)**et),lut%rv(i),lut%c(i))
       w=wk(i)
       integrand(i)=m*g*w**2
    END DO

    !Carries out the integration
    !Important to use basic trapezium rule because the integrand is messy
    p_1h=integrate_table(lut%nu,REAL(integrand),lut%n,1,lut%n,1)*(4.*pi)*(k**3)/(cosmic_density(cosm)*(2.*pi)**3)    

    DEALLOCATE(integrand)

    !Damping of the 1-halo term at very large scales
    ks=kstar(lut)

    !Prevents problems if k/ks is very large
    IF(ks>0.) THEN
       IF((k/ks)**2>7.) THEN
          fac=0.
       ELSE
          fac=exp(-((k/ks)**2))
       END IF
       p_1h=p_1h*(1.-fac)
    END IF

  END FUNCTION p_1h

  SUBROUTINE fill_sigtab(cosm)

    !This fills up tables of r vs. sigma(r) across a range in r!
    !It is used only in look-up for further calculations of sigma(r) and not otherwise!
    !and prevents a large number of calls to the sigint functions    
    IMPLICIT NONE
    REAL :: rmin, rmax
    REAL, ALLOCATABLE :: rtab(:), sigtab(:)
    REAL :: r, sig
    INTEGER :: i
    INTEGER, PARAMETER :: nsig=64 !Number of entries for sigma(R) tables
    TYPE(cosmology) :: cosm

    !These must be not allocated before sigma calculations otherwise when sigma(r) is called
    !otherwise sigma(R) looks for the result in the tables
    IF(ALLOCATED(cosm%r_sigma)) DEALLOCATE(cosm%r_sigma)
    IF(ALLOCATED(cosm%sigma)) DEALLOCATE(cosm%sigma)   

    !These values of 'r' work fine for any power spectrum of cosmological importance
    !Having nsig as a 2** number is most efficient for the look-up routines
    !rmin and rmax need to be decided in advance and are chosen such that
    !R vs. sigma(R) is a power-law below and above these values of R   
    rmin=1e-4
    rmax=1e3
    cosm%nsig=nsig

    IF(verbose) WRITE(*,*) 'SIGTAB: Filling sigma interpolation table'
    IF(verbose) WRITE(*,*) 'SIGTAB: Rmin:', rmin
    IF(verbose) WRITE(*,*) 'SIGTAB: Rmax:', rmax
    IF(verbose) WRITE(*,*) 'SIGTAB: Values:', nsig

    ALLOCATE(rtab(nsig),sigtab(nsig))

    DO i=1,nsig

       !Equally spaced r in log
       r=exp(log(rmin)+log(rmax/rmin)*float(i-1)/float(nsig-1))

       sig=sigma(r,0.,cosm)

       rtab(i)=r
       sigtab(i)=sig

    END DO

    !Must be allocated after the sigtab calulation above
    ALLOCATE(cosm%r_sigma(nsig),cosm%sigma(nsig))

    cosm%r_sigma=log(rtab)
    cosm%sigma=log(sigtab)

    DEALLOCATE(rtab,sigtab)

    IF(verbose) WRITE(*,*) 'SIGTAB: Done'
    IF(verbose) WRITE(*,*)

  END SUBROUTINE fill_sigtab

  FUNCTION sigma_cb(r,z,cosm)

    !Finds sigma_cold from look-up tables
    !In this version of HMcode sigma_cold=sigma
    !This is no longer true if massive neutrinos are considered
    IMPLICIT NONE
    REAL :: sigma_cb
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm

    sigma_cb=grow(z,cosm)*exp(find(log(r),cosm%r_sigma,cosm%sigma,cosm%nsig,3,3,2))

  END FUNCTION sigma_cb

  FUNCTION wk_tophat(x)

    !The normlaised Fourier Transform of a top-hat
    !Taylor expansion used for low |x| to avoid cancellation problems
    IMPLICIT NONE
    REAL :: wk_tophat
    REAL, INTENT(IN) :: x

    IF(x<0.01) THEN
       wk_tophat=1.-(x**2.)/10.
    ELSE
       wk_tophat=3.*(sin(x)-x*cos(x))/(x**3.)
    END IF

  END FUNCTION wk_tophat

  FUNCTION integrate(a,b,f,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: integrate
    REAL, INTENT(IN) :: a, b, acc
    INTEGER, INTENT(IN) :: iorder
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    DOUBLE PRECISION :: sum_n, sum_2n, sum_new, sum_old
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    INTERFACE
       FUNCTION f(x)
         REAL :: f
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       integrate=0.

    ELSE

       !Set the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_old=0.d0
       sum_new=0.d0

       DO j=1,jmax
          
          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             !The first go is just the trapezium of the end points
             f1=f(a)
             f2=f(b)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=f(x)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'INTEGRATE: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             !integrate=REAL(sum_new)
             !WRITE(*,*) 'INTEGRATE: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'INTEGRATE: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       integrate=REAL(sum_new)

    END IF

  END FUNCTION integrate

  FUNCTION integrate_table(x,y,n,n1,n2,iorder)

    !Integrates tables y(x)dx
    IMPLICIT NONE
    REAL :: integrate_table
    INTEGER, INTENT(IN) :: n, n1, n2
    REAL, INTENT(IN) :: x(n), y(n)
    REAL :: a, b, c, d, h
    REAL :: q1, q2, q3, qi, qf
    REAL :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
    DOUBLE PRECISION :: sum
    INTEGER :: i, i1, i2, i3, i4
    INTEGER, INTENT(IN) :: iorder

    sum=0.d0

    !I think if n1=n2 then the result will just be zero anyway
    !IF(n2<=n1) STOP 'INTEGRATE_TABLE: Error n2 must be greater than n1'

    IF(iorder==1) THEN

       !Sums over all Trapezia (a+b)*h/2
       DO i=n1,n2-1
          a=y(i+1)
          b=y(i)
          h=x(i+1)-x(i)
          sum=sum+(a+b)*h/2.d0
       END DO

    ELSE IF(iorder==2) THEN

       DO i=n1,n2-2

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

       DO i=n1,n2-1

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

    ELSE

       STOP 'INTEGRATE_TABLE: Error, order not specified correctly'

    END IF

    integrate_table=REAL(sum)

  END FUNCTION integrate_table

  FUNCTION sigma(r,z,cosm)

    !USE cosdef
    IMPLICIT NONE
    REAL :: sigma
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm
    INTEGER, PARAMETER :: iorder=3
    REAL, PARAMETER :: rsplit=1e-2

    IF(r>=rsplit) THEN
       sigma=sqrt(sigint0(r,z,cosm,acc,iorder))
    ELSE IF(r<rsplit) THEN
       sigma=sqrt(sigint1(r,z,cosm,acc,iorder)+sigint2(r,z,cosm,acc,iorder))
    ELSE
       STOP 'SIGMA: Error, something went wrong'
    END IF

  END FUNCTION sigma

  FUNCTION sigma_integrand(k,R,z,cosm)

    !The integrand for the sigma(R) integrals
    IMPLICIT NONE
    REAL :: sigma_integrand
    REAL, INTENT(IN) :: k, R, z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: y, w_hat

    IF(k==0.) THEN
       sigma_integrand=0.
    ELSE
       y=k*R
       w_hat=wk_tophat(y)
       sigma_integrand=p_lin(k,z,cosm)*(w_hat**2)/k
    END IF

  END FUNCTION sigma_integrand

  FUNCTION sigma_integrand_transformed(t,R,f,z,cosm)

    !The integrand for the sigma(R) integrals
    IMPLICIT NONE
    REAL :: sigma_integrand_transformed
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
       sigma_integrand_transformed=0.
    ELSE IF(t==1.) THEN
       !t=1 corresponds to k=0. when P(k)=0.
       sigma_integrand_transformed=0.
    ELSE
       !f(R) can be *any* function of R here to improve integration speed
       k=(-1.+1./t)/f(R)
       y=k*R
       w_hat=wk_tophat(y)
       sigma_integrand_transformed=p_lin(k,z,cosm)*(w_hat**2)/(t*(1.-t))
    END IF

  END FUNCTION sigma_integrand_transformed

  FUNCTION sigint0(r,z,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigint0
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: iorder
    REAL :: a, b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    REAL*8 :: sum_n, sum_2n, sum_new, sum_old
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    a=0.d0
    b=1.d0

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       sigint0=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_new=0.d0
       sum_old=0.d0

       DO j=1,jmax
          
          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             !The first go is just the trapezium of the end points
             f1=sigma_integrand_transformed(a,r,f0_rapid,z,cosm)
             f2=sigma_integrand_transformed(b,r,f0_rapid,z,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=sigma_integrand_transformed(x,r,f0_rapid,z,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'SIGINT0: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'SIGINT0: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       sigint0=REAL(sum_new)

    END IF

  END FUNCTION sigint0

  FUNCTION f0_rapid(r)

    !This is the 'rapidising' function to increase integration speed
    !for sigma(R). Found by trial-and-error
    IMPLICIT NONE
    REAL :: f0_rapid
    REAL, INTENT(IN) :: r
    REAL :: alpha
    REAL, PARAMETER :: rsplit=1e-2

    IF(r>rsplit) THEN
       !alpha 0.3-0.5 works well
       alpha=0.5
    ELSE
       !If alpha=1 this goes tits up
       !alpha 0.7-0.9 works well
       alpha=0.8
    END IF

    f0_rapid=r**alpha

  END FUNCTION f0_rapid

  FUNCTION sigint1(r,z,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigint1
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: iorder
    REAL :: a, b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    REAL*8 :: sum_n, sum_2n, sum_new, sum_old
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30

    a=r/(r+r**.5)
    b=1.d0

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       sigint1=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_new=0.d0
       sum_old=0.d0

       DO j=1,jmax
          
          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             !The first go is just the trapezium of the end points
             f1=sigma_integrand_transformed(a,r,f1_rapid,z,cosm)
             f2=sigma_integrand_transformed(b,r,f1_rapid,z,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=sigma_integrand_transformed(x,r,f1_rapid,z,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'SIGINT1: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence             
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'SIGINT1: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       sigint1=REAL(sum_new)

    END IF

  END FUNCTION sigint1

  FUNCTION f1_rapid(r)

    !This is the 'rapidising' function to increase integration speed
    !for sigma(R). Found by trial-and-error
    IMPLICIT NONE
    REAL :: f1_rapid
    REAL, INTENT(IN) :: r
    REAL, PARAMETER :: alpha=0.5  

    f1_rapid=r**alpha

  END FUNCTION f1_rapid

  FUNCTION sigint2(r,z,cosm,acc,iorder)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: sigint2
    REAL, INTENT(IN) :: r, z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL, INTENT(IN) :: acc
    INTEGER, INTENT(IN) :: iorder
    REAL :: a, b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    REAL*8 :: sum_n, sum_2n, sum_new, sum_old
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    REAL, PARAMETER :: C=10. !How far to go out in 1/r units for integral
    
    a=1./r
    b=C/r

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       sigint2=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_new=0.d0
       sum_old=0.d0

       DO j=1,jmax
          
          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             !The first go is just the trapezium of the end points
             f1=sigma_integrand(a,r,z,cosm)
             f2=sigma_integrand(b,r,z,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=sigma_integrand(x,r,z,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'SIGINT2: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             !WRITE(*,*) 'INTEGRATE_STORE: Nint:', n
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'SIGINT2: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       sigint2=REAL(sum_new)

    END IF

  END FUNCTION sigint2

  FUNCTION win(k,rv,c)

    !Calculates the Fourier Transform of the halo density profile
    !Normalised such that W(k=0)=1
    IMPLICIT NONE
    REAL :: win
    REAL, INTENT(IN) :: k, rv, c

    !Calls the analytic Fourier Transform of the NFW profile
    win=winnfw(k,rv,c)

    !Correct for the case of disasters (a bit sloppy, not sure if this is ever used)
    IF(win>1.) win=1.
    IF(win<0.) win=0.

  END FUNCTION win

  FUNCTION winnfw(k,rv,c)

    !The magical analytic Fourier Transform of the NFW profile
    IMPLICIT NONE
    REAL :: winnfw
    REAL, INTENT(IN) :: k, rv, c
    REAL :: si1, si2, ci1, ci2, ks
    REAL :: p1, p2, p3

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
    !This is mass with some factors missing (4*pi, rs, ...)
    IMPLICIT NONE
    REAL :: mass
    REAL, INTENT(IN) :: c

    mass=log(1.+c)-c/(1.+c)

  END FUNCTION mass

  FUNCTION bnu(nu)

    !Halo bias as a function of nu
    IMPLICIT NONE
    REAL :: bnu
    REAL, INTENT(IN) :: nu

    bnu=bst(nu)

  END FUNCTION bnu

  FUNCTION bst(nu)

    !Sheth Tormen bias
    IMPLICIT NONE
    REAL :: bst
    REAL, INTENT(IN) :: nu

    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707
    REAL, PARAMETER :: dc=1.686

    bst=1.+(q*(nu**2)-1.+2.*p/(1.+(q*nu**2)**p))/dc

  END FUNCTION bst

  FUNCTION gnu(nu)

    !Mass function
    IMPLICIT NONE
    REAL :: gnu
    REAL, INTENT(IN) :: nu

    gnu=gst(nu)

  END FUNCTION gnu

  FUNCTION gst(nu)

    !Sheth Tormen mass function!
    !Note I use nu=dc/sigma(M) and this Sheth & Tormen (1999) use nu=(dc/sigma)^2
    !This accounts for small differences in the equations
    IMPLICIT NONE
    REAL :: gst
    REAL, INTENT(IN) :: nu

    REAL, PARAMETER :: p=0.3
    REAL, PARAMETER :: q=0.707
    REAL, PARAMETER :: A=0.21616

    gst=A*(1.+((q*nu*nu)**(-p)))*exp(-q*nu*nu/2.)

  END FUNCTION gst

  FUNCTION gnubnu(nu)

    !g(nu) times b(nu)
    IMPLICIT NONE
    REAL :: gnubnu
    REAL, INTENT(IN) :: nu

    gnubnu=gnu(nu)*bnu(nu)

  END FUNCTION gnubnu
  
  FUNCTION omega_m(z,cosm)

    !This calculates Omega_m variations with z!
    IMPLICIT NONE
    REAL :: omega_m
    REAL, INTENT(IN) :: z
    REAL :: om_m
    TYPE(cosmology), INTENT(IN) :: cosm
    
    om_m=cosm%om_m
    omega_m=(om_m*(1.+z)**3.)/H2(z,cosm)

  END FUNCTION omega_m

  FUNCTION grow(z,cosm)

    !Scale-independent growth function at z=0
    IMPLICIT NONE
    REAL :: grow
    REAL, INTENT(IN) :: z
    REAL :: a
    TYPE(cosmology), INTENT(IN) :: cosm
    
    IF(z==0.) THEN
       grow=1.
    ELSE
       a=1./(1.+z)
       grow=find(a,cosm%a_grow,cosm%grow,cosm%ng,3,3,2)
    END IF

  END FUNCTION grow

  FUNCTION growint(a,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: growint
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    REAL*8 :: sum_n, sum_2n, sum_new, sum_old
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: iorder=3   

    !Integration range for integration parameter
    !Note a -> 1
    b=1.d0

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       growint=exp(0.)

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_new=0.d0
       sum_old=0.d0

       DO j=1,jmax
          
          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             !The first go is just the trapezium of the end points
             f1=growint_integrand(a,cosm)
             f2=growint_integrand(b,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=growint_integrand(x,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'GROWINT: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'GROWINT: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       growint=REAL(exp(sum_new))

    END IF

  END FUNCTION growint

  FUNCTION growint_integrand(a,cosm)

    !Integrand for the approximate growth integral
    IMPLICIT NONE
    REAL :: growint_integrand
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: gam
    
    IF(cosm%w<-1.) THEN
       gam=0.55+0.02*(1.+cosm%w)
    ELSE IF(cosm%w>-1) THEN
       gam=0.55+0.05*(1.+cosm%w)
    ELSE
       gam=0.55
    END IF

    !Note the minus sign here
    growint_integrand=-(Omega_m(-1.+1./a,cosm)**gam)/a
    
  END FUNCTION growint_integrand

  FUNCTION dispint(R,z,cosm)

    !Integrates between a and b until desired accuracy is reached
    !Stores information to reduce function calls
    IMPLICIT NONE
    REAL :: dispint
    REAL, INTENT(IN) :: z, R
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: a, b
    INTEGER :: i, j
    INTEGER :: n
    REAL :: x, dx
    REAL :: f1, f2, fx
    REAL*8 :: sum_n, sum_2n, sum_new, sum_old
    INTEGER, PARAMETER :: jmin=5
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: iorder=3   

    !Integration range for integration parameter
    !Note 0 -> infinity in k has changed to 0 -> 1 in x
    a=0.d0
    b=1.d0

    IF(a==b) THEN

       !Fix the answer to zero if the integration limits are identical
       dispint=0.

    ELSE

       !Reset the sum variable for the integration
       sum_2n=0.d0
       sum_n=0.d0
       sum_new=0.d0
       sum_old=0.d0

       DO j=1,jmax
          
          !Note, you need this to be 1+2**n for some integer n
          !j=1 n=2; j=2 n=3; j=3 n=5; j=4 n=9; ...'
          n=1+2**(j-1)

          !Calculate the dx interval for this value of 'n'
          dx=(b-a)/REAL(n-1)

          IF(j==1) THEN
             
             !The first go is just the trapezium of the end points
             f1=dispint_integrand(a,R,z,cosm)
             f2=dispint_integrand(b,R,z,cosm)
             sum_2n=0.5d0*(f1+f2)*dx
             sum_new=sum_2n
             
          ELSE

             !Loop over only new even points to add these to the integral
             DO i=2,n,2
                x=a+(b-a)*REAL(i-1)/REAL(n-1)
                fx=dispint_integrand(x,R,z,cosm)
                sum_2n=sum_2n+fx
             END DO

             !Now create the total using the old and new parts
             sum_2n=sum_n/2.d0+sum_2n*dx

             !Now calculate the new sum depending on the integration order
             IF(iorder==1) THEN  
                sum_new=sum_2n
             ELSE IF(iorder==3) THEN         
                sum_new=(4.d0*sum_2n-sum_n)/3.d0 !This is Simpson's rule and cancels error
             ELSE
                STOP 'DISPINT: Error, iorder specified incorrectly'
             END IF

          END IF

          IF((j>=jmin) .AND. (ABS(-1.d0+sum_new/sum_old)<acc)) THEN
             !jmin avoids spurious early convergence
             EXIT
          ELSE IF(j==jmax) THEN
             STOP 'DISPINT: Integration timed out'
          ELSE
             !Integral has not converged so store old sums and reset sum variables
             sum_old=sum_new
             sum_n=sum_2n
             sum_2n=0.d0
          END IF

       END DO

       dispint=REAL(sum_new)

    END IF

  END FUNCTION dispint

  FUNCTION dispint_integrand(theta,R,z,cosm)

    !This is the integrand for the velocity dispersion integral
    IMPLICIT NONE
    REAL :: dispint_integrand
    REAL, INTENT(IN) :: theta, z, R
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: k
    REAL, PARAMETER :: alpha=1.65 !Speeds up integral for large 'R'
    REAL, PARAMETER :: Rsplit=10. !Value to impliment speed up

    !Note that I have not included the speed up alpha and Rsplit
    !The choice of alpha=1.65 seemed to work well for R=100.
    !Rsplit=10 is thoughlessly chosen (only because 100.>10.)
    !Including this seems to make things slower (faster integration but slower IF statements?)

    IF(theta==0.d0 .OR. theta==1.d0) THEN
       dispint_integrand=0.d0
    ELSE
       !IF(r>Rsplit) THEN
       !   k=(-1.d0+1.d0/theta)/r**alpha
       !ELSE
       k=(-1.+1./theta)
       !END IF
       dispint_integrand=(p_lin(k,z,cosm)/k**2)*(wk_tophat(k*r)**2)/(theta*(1.-theta))
    END IF
    
  END FUNCTION dispint_integrand

  FUNCTION Si(x)

    !Sine integral function
    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
    IMPLICIT NONE
    REAL :: Si
    REAL, INTENT(IN) :: x
    REAL*8 :: x2, y, f, g, si8
    REAL*8, PARAMETER :: pi8=3.1415926535897932384626433d0

    IF(ABS(x)<=4.) THEN

       x2=x*x

       si8 = x*(1.d0+x2*(-4.54393409816329991d-2+x2*(1.15457225751016682d-3&
            +x2*(-1.41018536821330254d-5+x2*(9.43280809438713025d-8+x2*(-3.53201978997168357d-10&
            +x2*(7.08240282274875911d-13+x2*(-6.05338212010422477d-16))))))))/ &
            (1.+x2*(1.01162145739225565d-2 +x2*(4.99175116169755106d-5+&
            x2*(1.55654986308745614d-7+x2*(3.28067571055789734d-10+x2*(4.5049097575386581d-13&
            +x2*(3.21107051193712168d-16)))))))

       si=real(si8)

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

       si=real(pi8/2.d0-f*cos(x)-g*sin(x))

    ELSE

       STOP 'ERROR: Si, something went wrong'

    END IF

  END FUNCTION Si

  FUNCTION Ci(x)

    !Cosine integral function
    !Expansions for high and low x thieved from Wikipedia, two different expansions for above and below 4.
    IMPLICIT NONE
    REAL :: Ci
    REAL, INTENT(IN) :: x
    REAL*8 :: x2, y, f, g, ci8
    REAL*8, PARAMETER :: em_const=0.577215664901532861d0    

    IF(ABS(x)<=4.) THEN

       x2=x*x

       ci8=em_const+log(x)+x2*(-0.25d0+x2*(7.51851524438898291d-3+x2*(-1.27528342240267686d-4&
            +x2*(1.05297363846239184d-6+x2*(-4.68889508144848019d-9+x2*(1.06480802891189243d-11&
            +x2*(-9.93728488857585407d-15)))))))/ (1.+x2*(1.1592605689110735d-2+&
            x2*(6.72126800814254432d-5+x2*(2.55533277086129636d-7+x2*(6.97071295760958946d-10+&
            x2*(1.38536352772778619d-12+x2*(1.89106054713059759d-15+x2*(1.39759616731376855d-18))))))))

       ci=real(ci8)

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

       ci=real(f*sin(x)-g*cos(x))

    ELSE

       STOP 'ERROR: Ci, something went wrong'

    END IF

  END FUNCTION Ci

  FUNCTION derivative_table(x,xin,yin,n,iorder,imeth)

    !Given two arrays x and y such that y=y(x) this uses interpolation to calculate the derivative y'(x_i) at position x_i
    IMPLICIT NONE
    REAL :: derivative_table
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xin(n), yin(n)
    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i
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

    maxorder=3
    maxmethod=3

    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       !Reverse the arrays in this case
       CALL reverse(xtab,n)
       CALL reverse(ytab,n)
    END IF

    IF(iorder<1) STOP 'DERIVATIVE_TABLE: find order not specified correctly'
    IF(iorder>maxorder) STOP 'DERIVATIVE_TABLE: find order not specified correctly'
    IF(imeth<1) STOP 'DERIVATIVE_TABLE: Method of finding within a table not specified correctly'
    IF(imeth>maxmethod) STOP 'DERIVATIVE_TABLE: Method of finding within a table not specified correctly'

    IF(iorder==1) THEN

       IF(n<2) STOP 'DERIVATIVE_TABLE: Not enough points in your table for linear interpolation'

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

          i=table_integer(x,xtab,n,imeth)

          x2=xtab(i+1)
          x1=xtab(i)

          y2=ytab(i+1)
          y1=ytab(i)

       END IF

       CALL fit_line(a,b,x1,y1,x2,y2)
       derivative_table=a

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'DERIVATIVE_TABLE: Not enough points in your table'

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

          i=table_integer(x,xtab,n,imeth)
          
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

       IF(n<4) STOP 'DERIVATIVE_TABLE: Not enough points in your table'

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

          i=table_integer(x,xtab,n,imeth)

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

    ELSE

       STOP 'DERIVATIVE_TABLE: Error, iorder not specified correctly'

    END IF

  END FUNCTION derivative_table

  FUNCTION find(x,xin,yin,n,iorder,ifind,imeth)

    !Given two arrays x and y this routine interpolates to find the y_i value at position x_i
    IMPLICIT NONE
    REAL :: find
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xin(n), yin(n)
    REAL, ALLOCATABLE ::  xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i
    INTEGER, INTENT(IN) :: imeth, iorder, ifind

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !If the value required is off the table edge the interpolation is always linear

    !imeth = 1 => find x in xtab by crudely searching from x(1) to x(n)
    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    ALLOCATE(xtab(n),ytab(n))

    xtab=xin
    ytab=yin

    IF(xtab(1)>xtab(n)) THEN
       !Reverse the arrays in this case
       CALL reverse(xtab,n)
       CALL reverse(ytab,n)
    END IF

    IF(x<xtab(1)) THEN

       !Do a linear interpolation beyond the table boundary

       x1=xtab(1)
       x2=xtab(2)

       y1=ytab(1)
       y2=ytab(2)

       IF(imeth==1) THEN
          CALL fit_line(a,b,x1,y1,x2,y2)
          find=a*x+b
       ELSE IF(imeth==2) THEN
          find=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
       ELSE
          STOP 'FIND: Error, method not specified correctly'
       END IF

    ELSE IF(x>xtab(n)) THEN

       !Do a linear interpolation beyond the table boundary

       x1=xtab(n-1)
       x2=xtab(n)

       y1=ytab(n-1)
       y2=ytab(n)

       IF(imeth==1) THEN
          CALL fit_line(a,b,x1,y1,x2,y2)
          find=a*x+b
       ELSE IF(imeth==2) THEN
          find=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
       ELSE
          STOP 'FIND: Error, method not specified correctly'
       END IF

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

          i=table_integer(x,xtab,n,ifind)

          x1=xtab(i)
          x2=xtab(i+1)

          y1=ytab(i)
          y2=ytab(i+1)

       END IF

       IF(imeth==1) THEN
          CALL fit_line(a,b,x1,y1,x2,y2)
          find=a*x+b
       ELSE IF(imeth==2) THEN
          find=Lagrange_polynomial(x,1,(/x1,x2/),(/y1,y2/))
       ELSE
          STOP 'FIND: Error, method not specified correctly'
       END IF

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

          IF(imeth==1) THEN
             CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
             find=a*(x**2.)+b*x+c
          ELSE IF(imeth==2) THEN
             find=Lagrange_polynomial(x,2,(/x1,x2,x3/),(/y1,y2,y3/))
          ELSE
             STOP 'FIND: Error, method not specified correctly'
          END IF

       ELSE

          i=table_integer(x,xtab,n,ifind)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          IF(imeth==1) THEN
             !In this case take the average of two separate quadratic spline values
             CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
             find=(a*(x**2.)+b*x+c)/2.
             CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
             find=find+(a*(x**2.)+b*x+c)/2.
          ELSE IF(imeth==2) THEN
             !In this case take the average of two quadratic Lagrange polynomials
             find=(Lagrange_polynomial(x,2,(/x1,x2,x3/),(/y1,y2,y3/))+Lagrange_polynomial(x,2,(/x2,x3,x4/),(/y2,y3,y4/)))/2.
          ELSE
             STOP 'FIND: Error, method not specified correctly'
          END IF

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

          i=table_integer(x,xtab,n,ifind)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       IF(imeth==1) THEN
          CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
          find=a*x**3.+b*x**2.+c*x+d
       ELSE IF(imeth==2) THEN
          find=Lagrange_polynomial(x,3,(/x1,x2,x3,x4/),(/y1,y2,y3,y4/))
       ELSE
          STOP 'FIND: Error, method not specified correctly'
       END IF

    ELSE

       STOP 'FIND: Error, interpolation order specified incorrectly'

    END IF

  END FUNCTION find

  FUNCTION Lagrange_polynomial(x,n,xv,yv)

    !Computes the result of the nth order Lagrange polynomial at point x, L(x)
    IMPLICIT NONE
    REAL :: Lagrange_polynomial
    REAL, INTENT(IN) :: x, xv(n+1), yv(n+1)
    REAL :: l(n+1)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j

    !Initialise variables, one for sum and one for multiplication
    Lagrange_polynomial=0.
    l=1.

    !Loops to find the polynomials, one is a sum and one is a multiple
    DO i=0,n
       DO j=0,n
          IF(i .NE. j) l(i+1)=l(i+1)*(x-xv(j+1))/(xv(i+1)-xv(j+1))
       END DO
       Lagrange_polynomial=Lagrange_polynomial+l(i+1)*yv(i+1)
    END DO

  END FUNCTION Lagrange_polynomial

  FUNCTION table_integer(x,xtab,n,imeth)

    !Chooses between ways to find the integer location below some value in an array
    IMPLICIT NONE
    INTEGER :: table_integer
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xtab(n)
    INTEGER, INTENT(IN) :: imeth

    IF(imeth==1) THEN
       table_integer=linear_table_integer(x,xtab,n)
    ELSE IF(imeth==2) THEN
       table_integer=search_int(x,xtab,n)
    ELSE IF(imeth==3) THEN
       table_integer=int_split(x,xtab,n)
    ELSE
       STOP 'TABLE INTEGER: Method specified incorrectly'
    END IF

  END FUNCTION table_integer

  FUNCTION linear_table_integer(x,xtab,n)

    !Assuming the table is exactly linear this gives you the integer position
    IMPLICIT NONE
    INTEGER :: linear_table_integer
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xtab(n)
    REAL :: x1, x2, xn

    !Returns the integer (table position) below the value of x
    !eg. if x(3)=6. and x(4)=7. and x=6.5 this will return 6
    !Assumes table is organised linearly (care for logs)

    !n=SIZE(xtab)
    x1=xtab(1)
    x2=xtab(2)
    xn=xtab(n)   

    IF(x1>xn) STOP 'LINEAR_TABLE_INTEGER :: table in the wrong order'
    IF(ABS(-1.+float(n-1)*(x2-x1)/(xn-x1))>acc) STOP 'LINEAR_TABLE_INTEGER :: table does not seem to be linear'

    linear_table_integer=1+FLOOR(float(n-1)*(x-x1)/(xn-x1))

  END FUNCTION linear_table_integer

  FUNCTION search_int(x,xtab,n)

    !Does a stupid search through the table from beginning to end to find integer
    IMPLICIT NONE
    INTEGER :: search_int
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xtab(n)
    INTEGER :: i

    IF(xtab(1)>xtab(n)) STOP 'SEARCH_INT: table in wrong order'

    DO i=1,n
       IF(x>=xtab(i) .AND. x<=xtab(i+1)) EXIT
    END DO

    search_int=i

  END FUNCTION search_int

  FUNCTION int_split(x,xtab,n)

    !Finds the position of the value in the table by continually splitting it in half
    IMPLICIT NONE
    INTEGER :: int_split
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: x, xtab(n)
    INTEGER :: i1, i2, imid

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

    !Given xi, yi i=1,2 fits a line between these points
    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1
    REAL, INTENT(IN) :: x1, y1, x2, y2   

    a1=(y2-y1)/(x2-x1)
    a0=y1-a1*x1

  END SUBROUTINE fit_line

  SUBROUTINE fit_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)

    !Given xi, yi i=1,2,3 fits a quadratic between these points
    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1, a2
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3   

    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
    a0=y1-a2*(x1**2.)-a1*x1

  END SUBROUTINE fit_quadratic

  SUBROUTINE fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

    !Given xi, yi i=1,2,3,4 fits a cubic between these points
    IMPLICIT NONE
    REAL, INTENT(OUT) :: a, b, c, d
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3, x4, y4
    REAL :: f1, f2, f3    

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

  SUBROUTINE reverse(arry,n)

    !This reverses the contents of arry!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(INOUT) :: arry(n)
    INTEGER :: i
    REAL, ALLOCATABLE :: hold(:) 

    ALLOCATE(hold(n))

    hold=arry

    DO i=1,n
       arry(i)=hold(n-i+1)
    END DO

    DEALLOCATE(hold)

  END SUBROUTINE reverse

  FUNCTION file_length(file_name)

    !Finds the length of a file
    IMPLICIT NONE
    CHARACTER(len=64) :: file_name
    INTEGER ::n, file_length  

    OPEN(7,file=file_name)
    n=0
    DO
       n=n+1
       READ(7,*, end=301)
    END DO

    !301 is just the label to jump to when the end of the file is reached

301 CLOSE(7)  

    file_length=n-1

  END FUNCTION file_length

  SUBROUTINE fill_growtab(cosm)

    !Fills a table of the growth function vs. a   
    IMPLICIT NONE
    TYPE(cosmology) :: cosm
    INTEGER :: i
    REAL :: a, norm
    REAL, ALLOCATABLE :: d_tab(:), v_tab(:), a_tab(:)
    REAL :: ainit, amax, dinit, vinit
    INTEGER, PARAMETER :: n=64 !Number of entries for growth tables

    !The calculation should start at a z when Om_m(z)=1., so that the assumption
    !of starting in the g\propto a growing mode is valid (this will not work for early DE)
    ainit=0.001
    !Final should be a=1. unless considering models in the future
    amax=1.

    !These set the initial conditions to be the Om_m=1. growing mode
    dinit=ainit
    vinit=1.    

    IF(verbose) WRITE(*,*) 'GROWTH: Solving growth equation'
    CALL ode_growth(d_tab,v_tab,a_tab,0.,ainit,amax,dinit,vinit,acc,3,cosm)
    IF(verbose) WRITE(*,*) 'GROWTH: ODE done'

    !Normalise so that g(z=0)=1
    norm=find(1.,a_tab,d_tab,SIZE(a_tab),3,3,2)
    IF(verbose) WRITE(*,*) 'GROWTH: unnormalised g(a=1):', norm
    d_tab=d_tab/norm

    !This downsamples the tables that come out of the ODE solver (which can be a bit long)
    !Could use some table-interpolation routine here to save time
    IF(ALLOCATED(cosm%a_grow)) DEALLOCATE(cosm%a_grow,cosm%grow)
    cosm%ng=n

    IF(verbose) WRITE(*,*) 'GROWTH: Entries in table:', n

    ALLOCATE(cosm%a_grow(n),cosm%grow(n))
    DO i=1,n
       a=ainit+(amax-ainit)*float(i-1)/float(n-1)
       cosm%a_grow(i)=a
       cosm%grow(i)=find(a,a_tab,d_tab,SIZE(a_tab),3,3,2)
    END DO

    IF(verbose) WRITE(*,*) 'GROWTH: Done'
    IF(verbose) WRITE(*,*)

  END SUBROUTINE fill_growtab

  SUBROUTINE ode_growth(x,v,t,kk,ti,tf,xi,vi,acc,imeth,cosm)

    !Solves 2nd order ODE x''(t) from ti to tf and writes out array of x, v, t values 
    IMPLICIT NONE
    REAL :: xi, ti, tf, dt, acc, vi, x4, v4, t4, kk
    REAL :: kx1, kx2, kx3, kx4, kv1, kv2, kv3, kv4
    REAL*8, ALLOCATABLE :: x8(:), t8(:), v8(:), xh(:), th(:), vh(:)
    REAL, ALLOCATABLE :: x(:), v(:), t(:)
    INTEGER :: i, j, k, n, np, ifail, kn, imeth
    TYPE(cosmology) :: cosm
    INTEGER, PARAMETER :: jmax=30
    INTEGER, PARAMETER :: ninit=100

    !xi and vi are the initial values of x and v (i.e. x(ti), v(ti))
    !fx is what x' is equal to
    !fv is what v' is equal to
    !acc is the desired accuracy across the entire solution
    !imeth selects method

    IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(v)) DEALLOCATE(v)
    IF(ALLOCATED(t)) DEALLOCATE(t)

    DO j=1,jmax

       !Set the number of points for the forward integration
       n=ninit*(2**(j-1))
       n=n+1  

       !Allocate arrays
       ALLOCATE(x8(n),t8(n),v8(n))

       !Set the arrays to initialy be zeroes (is this neceseary?)
       x8=0.d0
       t8=0.d0
       v8=0.d0

       !Set the intial conditions at the intial time
       x8(1)=xi
       v8(1)=vi

       !Fill up a table for the time values
       CALL fill_table8(DBLE(ti),DBLE(tf),t8,n)

       !Set the time interval
       dt=(tf-ti)/float(n-1)

       !Intially fix this to zero. It will change to 1 if method is a 'failure'
       ifail=0

       DO i=1,n-1

          x4=real(x8(i))
          v4=real(v8(i))
          t4=real(t8(i))

          IF(imeth==1) THEN

             !Crude method
             kx1=dt*fd(x4,v4,kk,t4,cosm)
             kv1=dt*fv(x4,v4,kk,t4,cosm)

             x8(i+1)=x8(i)+kx1
             v8(i+1)=v8(i)+kv1
                  
          ELSE IF(imeth==2) THEN

             !Mid-point method
             !2017/06/18 - There was a bug in this part before. Luckily it was not used. Thanks Dipak Munshi.
             kx1=dt*fd(x4,v4,kk,t4,cosm)
             kv1=dt*fv(x4,v4,kk,t4,cosm)
             kx2=dt*fd(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)
             kv2=dt*fv(x4+kx1/2.,v4+kv1/2.,kk,t4+dt/2.,cosm)

             x8(i+1)=x8(i)+kx2
             v8(i+1)=v8(i)+kv2
             
          ELSE IF(imeth==3) THEN

             !4th order Runge-Kutta method (fast!)
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

          END IF

          !t8(i+1)=t8(i)+dt

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
          x=real(x8)
          v=real(v8)
          t=real(t8)
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

    !This is the fv in delta''=fv
    !Needed for growth function solution
    IMPLICIT NONE
    REAL :: fv
    REAL, INTENT(IN) :: d, v, k, a
    REAL :: f1, f2, z
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: crap

    crap=k !Stop compile-time warnings

    z=-1.+(1./a)

    f1=3.*omega_m(z,cosm)*d/(2.*(a**2.))
    f2=(2.+AH(z,cosm)/H2(z,cosm))*(v/a)

    fv=f1-f2

  END FUNCTION fv

  FUNCTION fd(d,v,k,a,cosm)

    !This is the fd in delta'=fd
    !Needed for growth function solution
    IMPLICIT NONE
    REAL :: fd
    REAL, INTENT(IN) :: d, v, k, a
    TYPE(cosmology), INTENT(IN) :: cosm
    REAL :: crap

    !Stop compile-time warnings
    crap=a
    crap=cosm%A
    crap=d
    crap=k

    fd=v

  END FUNCTION fd

  FUNCTION H2(z,cosm)

    !Calculates Hubble^2 in units such that H^2(z=0)=1.
    IMPLICIT NONE
    REAL :: H2
    REAL, INTENT(IN) :: z
    REAL :: om_m, om_v, a
    TYPE(cosmology), INTENT(IN) :: cosm
    
    om_m=cosm%om_m
    om_v=cosm%om_v
    a=1./(1.+z)
    H2=(om_m*(1.+z)**3.)+om_v*X_de(a,cosm)+((1.-om_m-om_v)*(1.+z)**2.)

  END FUNCTION H2

  FUNCTION AH(z,cosm)

    !Acceleration function - \ddot{a}/a
    IMPLICIT NONE
    REAL :: AH
    REAL, INTENT(IN) :: z
    REAL :: a
    TYPE(cosmology), INTENT(IN) :: cosm

    a=1./(1.+z)

    AH=cosm%om_m*(a**(-3.))+cosm%om_v*(1.+3.*w_de(a,cosm))*X_de(a,cosm)

    AH=-AH/2.

  END FUNCTION AH

  FUNCTION X_de(a,cosm)

    !The time evolution for Om_w for w(a) DE models
    IMPLICIT NONE
    REAL :: X_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm
    
    X_de=(a**(-3.*(1.+cosm%w+cosm%wa)))*exp(-3.*cosm%wa*(1.-a))

  END FUNCTION X_de

  FUNCTION w_de(a,cosm)

    !w(a) for DE models
    IMPLICIT NONE
    REAL :: w_de
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(IN) :: cosm
    
    w_de=cosm%w+(1.-a)*cosm%wa

  END FUNCTION w_de

END PROGRAM HMcode
