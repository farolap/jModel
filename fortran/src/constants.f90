module constants
  !======================================================
  !-
  !	Module containing parameters needed for several 
  !       other routines
  !-
  !-
  !   Pedro Peixoto - 03-03-2011
  !======================================================
  implicit none
  save 
  public

  !---------------------------------------------------
  !Kind attributions
  !---------------------------------------------------

  integer, parameter :: i2  = selected_real_kind(2,20)
  integer, parameter :: i4  = selected_real_kind(6,20)
  integer, parameter :: i8  = selected_real_kind(14,40)

  integer, parameter :: r4  = selected_real_kind(6,37)  
  integer, parameter :: r8  = selected_real_kind(12,100)
  integer, parameter :: r16 = max(r8,selected_real_kind(27,2400))

  !---------------------------------------------------
  ! General Parameters
  !---------------------------------------------------
  logical :: nlin=.true., &
      saveenergy=.false., &
      vectinv=.true.

  !Pi 
!  real(r8), parameter :: pi   = 4._r8* 0.78539816339744830961566 !datan (1._r8)
  real(r8), parameter :: pi   = 3.141592653589793d0 !datan (1._r8)
  real(r8), parameter :: pi2  = 2._r8*pi
  real(r8), parameter :: pio2 = pi/2._r8
  real(r8), parameter :: piby2 = pi*0.5_r8
  real(r8), parameter :: pio4 = pi/4._r8

  !Degrees to radians coversion (multiply to obtain conversion)
  real(r8), parameter :: deg2rad = pi / 180._r8

  !Radians to Degrees coversion (multiply to obtain conversion)
  real(r8), parameter :: rad2deg = 1._r8/deg2rad

  !Very small real number (aprox 1.10e-7)
  !real(r8), parameter :: eps = epsilon(1.)
  real(r8), parameter :: eps = epsilon(1.)

  !Very very small real number (aprox 1.10e-16)
  real(r8), parameter :: eps2 = epsilon(pi)
  
  !---------------------------------------------------
  ! Physical Parameters
  !---------------------------------------------------

  real*8 :: Hm     = 1d4;
!  real*8, parameter :: Hm     = 1d1;

  ! Earth mean radius (meters)
  real(r8), parameter :: rad     = 6.37122d6
!  real*8, parameter :: rad     = 1d0

  real(r8), parameter :: unitspharea    = 4._r8*pi

  !Gravitational accerlaration of the Earth (m/s^2)
!  real*8, parameter :: g        = 1d1
  real*8, parameter :: g        = 9.81d0

  ! Angular velocity of the Earth (rot/s)
  real (r8), parameter :: omega   = 2*pi/86400
  real (r8), parameter :: rotatn   = 7.292e-5_r8

  !Days to seconds
  real (r8), parameter :: day2sec = 86400_r8
  real (r8), parameter :: sec2day = 1._r8/86400_r8

  ! Dry air gas constant [J/(kg K)]
  real(r8), parameter :: rdry     = 287.   

  ! Dry air spec heat at const P [J/(kg K)]
  real(r8), parameter :: cp       = 1004.  

  ! Dry air spec heat at const vol [J/(kg K)]
  real(r8), parameter :: cv       = 717.              

  ! Water vapor gas constant [J/(kg K)]
  real(r8), parameter :: rvap     = 461.               

  ! Reference pressure [Pa]
  real(r8), parameter :: p00      = 1.e5            

  ! 0 Celsius temperature [K]
  real(r8), parameter :: t00      = 273.15         

  !---------------------------------------------------
  ! Global constant variables
  !---------------------------------------------------
  !Flag for verbose output
  logical :: showonscreen=.false., &
             tiledarea=.false.

  !Simulation to be done
  integer (i4) :: simulcase

  !---------------------------------------------------
  ! Density Interpolation File Parameters
  !---------------------------------------------------

  integer, parameter :: n_lat = 4*180+1
  integer, parameter :: n_lon = 4*360+1
  real (kind=8), parameter :: latmin = -pio2
  real (kind=8), parameter :: latmax = pio2
  real (kind=8), parameter :: lonmin = -pi
  real (kind=8), parameter :: lonmax = pi

  !---------------------------------------------------
  ! Altitude Parameters
  !---------------------------------------------------
  integer, parameter :: nlat_alt = 4*180+1 
  integer, parameter :: nlon_alt = 4*360+1
  character (len=60), parameter::  altdir = "altitude/"

end module constants
