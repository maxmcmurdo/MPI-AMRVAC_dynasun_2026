module mod_usr

  use mod_uawsom
  use mod_global_parameters
  use mod_random

  implicit none

  double precision :: usr_grav, Hr
  double precision :: stefboltz,heatunit,B0,theta,SRadius,kx,ly,bQ0, Ttop, rhohc, Heatscale, hc, dya
  double precision, save :: avgrho0 = 1.9d5, avgrho = 1.9d5, q
  double precision, allocatable :: pa(:),ra(:), ya(:)
  double precision, allocatable, dimension(:), save :: rdist, rho, p, vr
  integer, save :: fsize
  integer, parameter :: jmax = 10000000

  ! The pseudo-random number generator
  integer :: n = 12
  type(rng_t) :: rng
  integer, parameter :: i8 = selected_int_kind(18)
  integer(i8)        :: seed(16334132)

  integer, parameter :: st_maxmodes=1000
  !OU variance corresponding to decay time and energy input rate
  double precision,save :: st_OUvar
 
  !Number of modes
  integer, save :: st_nmodes

  double precision ,save, dimension(3,st_maxmodes) :: st_mode, st_aka, st_akb
  double precision,save, dimension(6*st_maxmodes) ::st_OUphases
  double precision,save, dimension(st_maxmodes) ::st_ampl


  logical,save :: st_useStir, st_computeDt

  logical,save :: reproducible=.false.,saveReproducible=.true.
  integer,save :: randomSaveUnit

  double precision,save :: st_decay = 0.05
  double precision,save :: st_energy = 2.d2
  double precision,save :: st_stirmin = 1.00
  double precision,save :: st_stirmax = 8.00
  double precision,save :: st_spectform = 0.0
  double precision,save :: st_solweight = 0.3
  double precision,save :: st_solweightnorm = 1.00
  integer,save :: st_freq = 10

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_1.75D")

    unit_length        = 6.961d10                                     ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_init_vector_potential=>initvecpot_usr
    usr_set_B0          => specialset_B0
    usr_set_J0          => specialset_J0
    usr_special_bc      => specialbound_usr
    usr_aux_output      => specialvar_output
    usr_process_global  => usrprocess_global
    usr_add_aux_names   => specialvarnames_output 

    call uawsom_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    integer :: j, stat, ixO^L
    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=1.d-4/heatunit ! background heating power density
    stefboltz = const_sigma/heatunit/unit_length*unit_temperature**4
    Hr = 1.55d7/unit_length ! 155 km 
    Heatscale = 5.d9
    ! a coronal height at 10 Mm
    hc=1.d9/unit_length
    ! the number density at the coronal height
    rhohc=6.2d8/unit_numberdensity

    dya = 6*hc/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    SRadius=6.961d10/unit_length ! Solar radius
    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic

    fsize=0
    open(unit=100, file='rho.xyz', status = 'old', action = 'read')
    do
        read(100,*,iostat=stat)
        if (stat /= 0) exit
        fsize=fsize+1
    end do
    close(unit=100)
    allocate(rdist(fsize))
    allocate(rho(fsize))
    allocate(p(fsize))
    allocate(vr(fsize))

    open(unit=100, file='rho.xyz', status = 'old', action = 'read')
    open(unit=200, file='p.xyz', status = 'old', action = 'read')
    open(unit=300, file='vr.xyz', status = 'old', action = 'read')
    
    do j=1,fsize
       read(100,*) rdist(j), rho(j)
       read(200,*) rdist(j), p(j)
       read(300,*) rdist(j), vr(j)
    end do
    close(unit=100)
    close(unit=200)
    close(unit=300)
 
    !call Stir_init(restart=.false.)

  end subroutine initglobaldata_usr

  subroutine inithdstatic
    use mod_solar_atmosphere
    use mod_global_parameters
    ! initialize the table in a vertical line through the global domain
    integer :: j,na,ibc
    double precision, allocatable :: Ta(:),gg(:)
    double precision :: rpho,Tpho,wtra,res,rhob,pb,htra,htrc,Ttr,Ttrc,Fc,invT,kappa
    integer :: simple_temperature_curve

    simple_temperature_curve=2

    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))

    select case(simple_temperature_curve)
    case(0)
      rpho=2.d14/unit_numberdensity ! number density at the bottom of height table
      Tpho=6.161d3/unit_temperature ! temperature of chromosphere
      Ttop=1.5d6/unit_temperature ! estimated temperature in the top
      htra=2.d8/unit_length ! height of initial transition region
      wtra=2.7d7/unit_length ! width of initial transition region 
      Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
      Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
      kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      do j=1,jmax
         ya(j)=(dble(j)-0.5d0)*dya
         if(ya(j)>htra) then
           Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
         else
           Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-2.7d7/unit_length)/wtra)+1.d0)
         endif
         gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      enddo
      !! solution of hydrostatic equation 
      ra(1)=rpho
      pa(1)=rpho*Tpho
      invT=gg(1)/Ta(1)
      invT=0.d0
      do j=2,jmax
         invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
         pa(j)=pa(1)*dexp(invT*dya)
         ra(j)=pa(j)/Ta(j)
      end do
    case(1)
      do j=1,jmax
         ! get height table
         ya(j)=(dble(j)-0.5d0)*dya
         ! get gravity table
         gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      enddo
      ! get density and pressure table in hydrostatic state with a preset temperature table
      call get_atm_para(ya,ra,pa,gg,jmax,'AL-C7',hc,rhohc)
    case(2)
     rpho=2.d14/unit_numberdensity ! number density at the bottom of height table
     Tpho=6.161d3/unit_temperature ! temperature of chromosphere
     htra=2.d8/unit_length ! height of initial transition region
     htrc=0.002! height of initial transition region
     Ttrc=(0.75d0*0.0058d0**4*(exp(-htrc/Hr)+&
                0.71d0-0.133d0/(1.d0 + 0.15d0*exp(-htrc/Hr)**0.73d0)**17.4d0))**0.25d0
     Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
     kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
     Ttop=1.5d6/unit_temperature ! estimated temperature in the top
     Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
     wtra=2.7d7/unit_length ! width of initial transition region 
     do j=1,jmax
      ya(j)=(dble(j)-0.5d0)*dya
      if (ya(j) < htrc) then
        Ta(j) = (0.75d0*0.0058d0**4*(exp(-ya(j)/Hr)+&
                0.71d0-0.133d0/(1.d0 + 0.15d0*exp(-ya(j)/Hr)**0.73d0)**17.4d0))**0.25d0
      else if (ya(j)>htra) then
               Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
           else
               Ta(j)=Ttrc+0.5d0*(Ttop-Ttrc)*(tanh((ya(j)-htra-2.7d7/unit_length)/wtra)+1.d0)
      endif
      gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
     enddo
     ra(1)=rpho
     pa(1)=rpho*Ta(1)
     invT=gg(1)/Ta(1)
     invT=0.d0
     do j=2,jmax
      invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
      pa(j)=pa(1)*dexp(invT*dya)
      ra(j)=pa(j)/Ta(j)
     end do

    end select

    deallocate(gg,Ta)

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    ! initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, allocatable :: Ta(:),gg(:)
    double precision :: rpho,Tpho,wtra,res,rhob,pb,htra,htrc,Ttr,Ttrc,Fc,invT,kappa,delta
    integer :: j,ibc
    double precision, parameter :: pi = 3.141592653589793d0

    double precision :: trans, Hsc, Tmp, rbg, xdist
    integer :: ix^D, na, ind

    logical, save :: first = .true.

    if (mype==0 .and. first) then
      write(*,*)'Unit time (s): ',unit_time
      write(*,*)'Unit velocity (km/s):',unit_velocity*1.d-5
      write(*,*)'Unit heat: ',heatunit
      write(*,*)'Unit pressure: ',unit_pressure
      write(*,*)'Unit density: ',unit_density
      write(*,*)'Unit B field (G):',unit_magneticfield
      first = .false.
    endif

    !trans = 1.d0 + 6*hc    
    !Tmp = 1.5d0
    !rbg = 1.d2
    !Hsc = Tmp/(2.74d4*unit_length/unit_velocity**2)
     
    !{do ix^DB=ixOmin^DB,ixOmax^DB\}
    !  if (x(ix^D,1) .le. trans) then
    !    na = floor((x(ix^D,1) - xprobmin1)/dya + 0.5d0) 
    !    res = x(ix^D,1) - xprobmin1 - (dble(na) - 0.5d0)*dya
    !    w(ix^D,rho_) = ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
    !    w(ix^D,p_) = pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na)) 
    !  endif
    !{end do\}
   
    !where(x(ixO^S,1) >= trans)
    !  w(ixO^S,rho_) = minval(ra)*dlog(x(ixO^S,1) - trans + 1.d0 + 0.1)**(-2.5)/(dlog(1.1d0)**(-2.5))
    !  w(ixO^S,p_) = minval(pa)*dlog(x(ixO^S,1) - trans + 1.d0 + 0.1)**(-2.5)/(dlog(1.1d0)**(-2.5))
    !  w(ixO^S,mom(1)) = 0.1 * (x(ixO^S,1) - trans)
    !endwhere

    w(ixO^S,mom(:)) = zero
     
    !pre-tabulated data
    do ix1 = ixOmin1,ixOmax1
       xdist = x(ix1,1)
       call findindex(xdist, fsize, rdist, ind)
       w(ix1^%1ixO^S,rho_) = rho(ind)
       w(ix1^%1ixO^S,p_) = p(ind)*(2.d0/3.d0)
       w(ix1^%1ixO^S,mom(1)) = vr(ind)
    end do
 
    if(B0field) then
      w(ixO^S,mag(:))=zero
    else if(stagger_grid) then
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)
      call uawsom_face_to_center(ixO^L,block)
      !w(ixO^S,mag(3))= 0.d0
    else
      w(ixO^S,mag(1))= B0*x(ixO^S,1)**(-2.0)
      w(ixO^S,mag(2))= 0.d0
      w(ixO^S,mag(3))= 0.d0
    endif
     
    ! wkplus/minus are kink and wAplus/minus are Alfven waves 
 
    delta = 3e9/unit_length

    !do ix1 = ixOmin1,ixOmax1
    !  if (x(ix1,1) < 1 + delta) then 
    !    w(ix1,wkminus_) = 1.67d-2*(0.5d0*cos(pi*(x(ix1,1)-1.0d0)/delta)+0.5d0)
    !  else
    !    w(ix1,wkminus_) = 0.0d0
    !  end if
    !end do

    !do ix1 = ixOmin1,ixOmax1
    !  if (x(ix1,1) > xprobmax1 - delta) then
    !    w(ix1,wkplus_) = 0.005d0*1.67d-2*(0.5d0*cos(pi*(x(ix1,1)-xprobmax1)/delta)+0.5d0) 
    !  else 
    !    w(ix1,wkplus_) = 0.d0
    !  end if
    !end do

    !> Luka found this value stabilises the atmosphere with radiative losses on
    !do ix1 = ixOmin1,ixOmax1
    !  if (x(ix1,1) < 1 + delta) then 
    !    w(ix1,wkminus_) = 1.23d-2*(0.5d0*cos(pi*(x(ix1,1)-1.0d0)/delta)+0.5d0)
    !    w(ix1,wAminus_) = 1.23d-2*(0.5d0*cos(pi*(x(ix1,1)-1.0d0)/delta)+0.5d0)
    !  else
    !    w(ix1,wkminus_) = 0.0d0
    !    w(ix1,wAminus_) = 0.0d0
    !  end if
    !end do

    do ix1 = ixOmin1,ixOmax1
      if (x(ix1,1) < 1 + delta) then 
        !w(ix1,wkplus_) = 1.67d-2*0.4d0*(cos(pi*(x(ix1,1)-1.0d0)/delta)+1.5d0)  
        w(ix1,wAminus_) = 1.67d-2*0.4d0*(cos(pi*(x(ix1,1)-1.0d0)/delta)+1.5d0)  
      else
        !w(ix1,wkplus_) = 1.67d-2*0.2d0 !> usually zero but now has energy everywhere initially
        w(ix1,wAminus_) = 1.67d-2*0.2d0 !> usually zero but now has energy everywhere initially
      end if
    end do

    !do ix1 = ixOmin1,ixOmax1
    !  if (x(ix1,1) > xprobmax1 - delta) then
    !    w(ix1,wAplus_) = 1.d-7*1.67d-2*0.4d0*(cos(pi*(x(ix1,1)-xprobmax1)/delta)+1.5d0) 
    !  else 
    !    w(ix1,wAplus_) = 1.d-7*1.67d-2*0.2d0
    !  end if
    !end do

    !do ix1 = ixOmin1,ixOmax1
    !  if (x(ix1,1) < 1 + delta) then 
    !    w(ix1,wkminus_) = 5.d0*1.23d-2*0.4d0*(cos(pi*(x(ix1,1)-1.0d0)/delta)+1.5d0) 
    !    w(ix1,wAminus_) = 5.d0*1.23d-2*0.4d0*(cos(pi*(x(ix1,1)-1.0d0)/delta)+1.5d0) 
    !  else
    !    w(ix1,wkminus_) = 5.d0*1.23d-2*0.2d0 !> usually zero but now has energy everywhere initially
    !    w(ix1,wAminus_) = 5.d0*1.23d-2*0.2d0 !> usually zero but now has energy everywhere initially
    !  end if
    !end do

    !w(ixO^S,wkminus_) = 1.d-6
    !w(ixO^S,wkplus_) = 0.d0!1.d-8 + 1.d-6*exp(-(x(ixO^S,1)-1.1435d0)**2.d0/0.05d0**2.d0)
    !w(ixO^S,wAplus_) = 1.d-8 + 1.d-6*exp(-(x(ixO^S,1)-1.1435d0)**2.d0/0.05d0**2.d0)
    
    call uawsom_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine findindex(xdist, fsize, rdist, ind)
    integer, intent(in) :: fsize
    double precision, intent(in) :: xdist
    double precision, intent(in) :: rdist(fsize)
    integer, intent(out) :: ind
    
    integer :: i
    double precision :: least

    least = xprobmax1
    ind = 1
    
    do i = 1, fsize
       if (least .gt. abs(xdist - rdist(i))) then
         least = abs(xdist - rdist(i))
         ind = i
       endif
    end do

  end subroutine findindex 

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    !if (idir==3) then
    !  A(ixC^S) = -B0/dtan(xC(ixC^S,2))/xC(ixC^S,1)
    !else
    !  A(ixC^S) = 0.d0
    !end if

  end subroutine initvecpot_usr


  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! In 3D BC govern: rho, m1, m2, m3, e, b1, b2, b3, wA+/-, wk+/-
    ! In 1D BC govern: rho, m1, e, b1, wA+/-, wk+/- (8 items)
    ! Currently, there is no wA+/- or wk+, so need below to cover: rho, m1, e, b1, wk-

    ! TVD has conditions for rho based off his dT/dz = 0 conversion to d(Eth/dens*(gamma-1))/dz = 0. We can do this by setting P by extrapolating
    ! the HSE using FD, then we can assume isothermality by saying T = P/rho, by setting rho as if the ghost cells maintain isothermality.

    ! Momentum must be set to 0 at the bottom boundary 
    ! Momentum must be set with a zero gradient at top boundary allows mass to flow out of the box

    ! Energy is set to not include the wave energy at the bottom boundary 
    ! Energy is set to include the wave energy at the top boundary 

    ! Magnetic field b1 is kept constant in all cases and does not change throughout the simulation so needs no BC (I think).

    ! TODO: figure out why Norbert had:       w(ixOmin1,mom(:)) = w(ixOmin1-1,mom(:))
    !                                         w(ixOmax1,mom(:)) = w(ixOmin1-2,mom(:))
    ! for the top boundary.

    use mod_global_parameters
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    double precision :: bQgrid(ixI^S), wave_ref_fraction !(ixI^S)
    
    double precision :: pth(ixI^S),Pin(ixI^S),Tmp(ixI^S),h1,h2,h3,h4,Hh1,Hh2
    double precision :: Ti,Temp, xdist, ggrid(ixI^S), gradvA_int(ixI^S), vA_int(ixI^S)
    
    integer:: ixIM^L, ix^D, idir
   
    integer, save :: ind1, ind2, iter = 1
 
    logical, save :: first = .true.

    select case(iB)
    
    case(1)

      ixIMmin1=ixOmax1+1;ixIMmax1=ixOmax1+4;
      
      !do ix1=ixOmin1,ixOmax1
      ! do idir=1,ndir
      !  w(ix1^%1ixO^S,mom(idir)) = w(ixOmax1+1^%1ixO^S,mom(idir)) ! zero gradient boundary condition at LH boundary
      ! end do
      !end do

      w(ixOmax1^%1ixO^S,mag(1)) = w(ixOmax1+1^%1ixO^S,mag(1)) 
      w(ixOmin1^%1ixO^S,mag(1)) = w(ixOmax1+1^%1ixO^S,mag(1))

      w(ixOmax1^%1ixO^S,mom(1)) = w(ixOmax1+1^%1ixO^S,mom(1)) 
      w(ixOmin1^%1ixO^S,mom(1)) = w(ixOmax1+1^%1ixO^S,mom(1))

      w(ixOmin1^%1ixO^S,wkminus_) = w(ixOmax1+1^%1ixO^S,wkminus_) ! Need wkminus in the ghost cells such that energy is calculated properly
      w(ixOmax1^%1ixO^S,wkminus_) = w(ixOmax1+1^%1ixO^S,wkminus_)

      w(ixOmin1^%1ixO^S,wkplus_) = w(ixOmax1+1^%1ixO^S,wkplus_)  ! wk+ open BBC to allow outflow (equally could set to zero)
      w(ixOmax1^%1ixO^S,wkplus_) = w(ixOmax1+1^%1ixO^S,wkplus_)

      w(ixOmin1^%1ixO^S,wAminus_) = 1.67d-2 !w(ixOmax1+1^%1ixO^S,wAminus_) ! Need wAminus in the ghost cells such that energy is calculated properly
      w(ixOmax1^%1ixO^S,wAminus_) = 1.67d-2 !w(ixOmax1+1^%1ixO^S,wAminus_)

      w(ixOmin1^%1ixO^S,wAplus_) = w(ixOmax1+1^%1ixO^S,wAplus_)
      w(ixOmax1^%1ixO^S,wAplus_) = w(ixOmax1+1^%1ixO^S,wAplus_)
       
      call uawsom_get_pthermal(w,x,ixI^L,ixIM^L,pth)

      !Tmp(ixOmax1^%1ixO^S) = pth(ixOmax1+1^%1ixO^S)/(w(ixOmax1+1^%1ixO^S,rho_))   ! isothermal Temp into ghost cells
      h1 = x(ixOmax1+1,1) - x(ixOmax1,1) ! grid spacing
      
      ! Pressure update with constant gravity. Uses HSE to extrapolate back into the ghost cells the correct pressure
      ! - sign might appear wrong way round but gravity is positive so swaps round +-
      Pin(ixOmax1^%1ixO^S) = pth(2+ixOmax1^%1ixO^S) - 2*h1*usr_grav*w(ixOmax1+1^%1ixO^S,rho_)

      ! Old density update assumes isothermality in the boundary
      !w(ixOmax1^%1ixO^S,rho_) = Pin(ixOmax1^%1ixO^S)/(Tmp(ixOmax1^%1ixO^S))
      !w(ixOmax1^%1ixO^S,rho_) = w(ixOmax1+1^%1ixO^S,rho_)
      w(ixOmax1^%1ixO^S,rho_) = w(ixOmax1+2^%1ixO^S,rho_) - 2*h1*usr_grav*w(ixOmax1+1^%1ixO^S,rho_)**2.d0/pth(ixOmax1+1^%1ixO^S)

      Tmp(ixOmax1^%1ixO^S) = Pin(ixOmax1^%1ixO^S)/(w(ixOmax1^%1ixO^S,rho_))

      ! Energy in the ghost cell ignores contribution from wave energy TODO: Figure out why this is commom and what it means
      w(ixOmax1^%1ixO^S,e_) = Pin(ixOmax1^%1ixO^S)/(uawsom_gamma-1) + &
                   0.5d0*(w(ixOmax1^%1ixO^S,mag(1))**2 + w(ixOmax1^%1ixO^S,mag(2))**2 + w(ixOmax1^%1ixO^S,mag(3))**2) + &
                   0.5d0*(w(ixOmax1^%1ixO^S,mom(1))**2 + w(ixOmax1^%1ixO^S,mom(2))**2 + w(ixOmax1^%1ixO^S,mom(3))**2)/w(ixOmax1^%1ixO^S,rho_) + &
                   w(ixOmax1^%1ixO^S,wkplus_) + w(ixOmax1^%1ixO^S,wkminus_) + &
                   w(ixOmax1^%1ixO^S,wAplus_) + w(ixOmax1^%1ixO^S,wAminus_) 

      ! Arbitrary again when the grid spacing is equal
      h1 = x(ixOmax1,1) - x(ixOmin1,1)
      
      ! Pressure update with constant gravity extends the HSE into the ghost cells again
      Pin(ixOmin1^%1ixO^S) = pth(ixOmax1+1^%1ixO^S) - 2*h1*usr_grav*w(ixOmax1^%1ixO^S,rho_)

      ! Temp calculated to use for density later
      !Tmp(ixOmin1^%1ixO^S) = Pin(ixOmax1^%1ixO^S)/(w(ixOmax1^%1ixO^S,rho_))
      
      ! Old density update based on Norberts commented out work 
      !w(ixOmin1^%1ixO^S,rho_) = Pin(ixOmin1^%1ixO^S)/Tmp(ixOmin1^%1ixO^S)
      !w(ixOmin1^%1ixO^S,rho_) = w(ixOmax1+1^%1ixO^S,rho_)
      w(ixOmin1^%1ixO^S,rho_) = w(ixOmax1+1^%1ixO^S,rho_) - 2*h1*usr_grav*w(ixOmax1^%1ixO^S,rho_)**2.d0/Pin(ixOmax1^%1ixO^S)

      Tmp(ixOmin1^%1ixO^S) = Pin(ixOmin1^%1ixO^S)/(w(ixOmin1^%1ixO^S,rho_))

      ! Energy in ghost cell, same todo as before
      w(ixOmin1^%1ixO^S,e_) = Pin(ixOmin1^%1ixO^S)/(uawsom_gamma-1) + &
                   0.5d0*(w(ixOmin1^%1ixO^S,mag(1))**2 + w(ixOmin1^%1ixO^S,mag(2))**2 + w(ixOmin1^%1ixO^S,mag(3))**2) + &
                   0.5d0*(w(ixOmin1^%1ixO^S,mom(1))**2 + w(ixOmin1^%1ixO^S,mom(2))**2 + w(ixOmin1^%1ixO^S,mom(3))**2)/w(ixOmin1^%1ixO^S,rho_) + &
                   w(ixOmin1^%1ixO^S,wkplus_) + w(ixOmin1^%1ixO^S,wkminus_) + &
                   w(ixOmin1^%1ixO^S,wAplus_) + w(ixOmin1^%1ixO^S,wAminus_)

      !if (mod(it,1000)==0 .and. mype==0) then 
      !  write(*,*)"B bot min ghost = ", w(ixOmin1,mag(1))
      !  write(*,*)"B bot min ghost = ", w(ixOmin1,mag(1))
      !  write(*,*)"B bot max cell  = ", w(ixOmax1+1,mag(1))
      !end if
     
    case (2) ! Upper radial boundary

     if (iprob==0) then

      !pre-tabulated data
       if (first) then
         xdist = x(ixOmin1,1)
         call findindex(xdist, fsize, rdist, ind1)
         xdist = x(ixOmax1,1)
         call findindex(xdist, fsize, rdist, ind2)
         first = .false.
       endif
        w(ixOmin1^%1ixO^S,rho_) = rho(ind1)
        w(ixOmin1^%1ixO^S,p_) = p(ind1)
        w(ixOmin1^%1ixO^S,mom(1)) = vr(ind1)
        w(ixOmax1^%1ixO^S,rho_) = rho(ind2)
        w(ixOmax1^%1ixO^S,p_) = p(ind2)
        w(ixOmax1^%1ixO^S,mom(1)) = vr(ind2)
      
       w(ixOmin1,mom(2:3)) = w(ixOmin1-1,mom(2:3))
       w(ixOmax1,mom(2:3)) = w(ixOmin1-1,mom(2:3))
      
      call uawsom_to_conserved(ixI^L,ixO^L,w,x)
     else if (iprob==1) then
      ixIMmin1=ixOmin1-2;ixIMmax1=ixOmin1-1;

      !w(ixOmax1^%1ixO^S,mag(1)) = w(ixOmax1+1^%1ixO^S,mag(1)) 
      !w(ixOmin1^%1ixO^S,mag(1)) = w(ixOmax1+1^%1ixO^S,mag(1))

      w(ixOmin1^%1ixO^S,mom(1)) = w(ixOmin1-1^%1ixO^S,mom(1))
      w(ixOmax1^%1ixO^S,mom(1)) = w(ixOmin1-1^%1ixO^S,mom(1))

      w(ixOmin1^%1ixO^S,wkminus_) = w(ixOmin1-1^%1ixO^S,wkminus_)
      w(ixOmax1^%1ixO^S,wkminus_) = w(ixOmin1-1^%1ixO^S,wkminus_)

      w(ixOmin1^%1ixO^S,wkplus_) = w(ixOmin1-1^%1ixO^S,wkplus_)
      w(ixOmax1^%1ixO^S,wkplus_) = w(ixOmin1-1^%1ixO^S,wkplus_)

      !if (mod(it,10000) == 0 .and. mype == 7) then
      !  write(*,*) 'w(ixOmin1^%1ixO^S,wkplus_) = ', w(ixOmin1^%1ixO^S,wkplus_)
      !  write(*,*)'w(ixOmax1^%1ixO^S,wkplus_) = ', w(ixOmax1^%1ixO^S,wkplus_)
      !  write(*,*)'w(ixOmax1+1^%1ixO^S,wkplus_) = ', w(ixOmin1-1^%1ixO^S,wkplus_)
      !end if

      w(ixOmin1^%1ixO^S,wAminus_) = w(ixOmin1-1^%1ixO^S,wAminus_)
      w(ixOmax1^%1ixO^S,wAminus_) = w(ixOmin1-1^%1ixO^S,wAminus_)

      w(ixOmin1^%1ixO^S,wAplus_) = 0.d0! 5.d0*4.d-4
      w(ixOmax1^%1ixO^S,wAplus_) = 0.d0! 5.d0*4.d-4 

      !vA_int(ixOmin1-10:ixOmax1) = dsqrt((20.d0/x(ixOmin1-10:ixOmax1,1)**2.d0)/w(ixOmin1-10:ixOmax1,rho_))
      !call gradient(vA_int,ixI^L,ixO^L,ndim,gradvA_int)
      !do ix1 = ixOmin1,ixOmax1
      !  gradvA_int(ix1) = gradvA_int(ixOmin1-1)
      !end do

      !if (mod(it,5000)==0) then
      !  write(*,*) "vA_int     = ", vA_int(ixOmin1-10:ixOmax1)
      !  write(*,*) "gradvA_int = ", gradvA_int(ixOmin1-10:ixOmax1)
      !  write(*,*) "ref term: ", ((w(ixOmin1-1,mom(1)) + vA_int(ixOmin1-1))/vA_int(ixOmin1-1))*gradvA_int(ixOmin1-1)*0.1d0*w(ixOmin1-1^%1ixO^S,wAminus_)
      !end if

      !w(ixOmin1^%1ixO^S,wAplus_) = ((w(ixOmin1-1,mom(1)) + vA_int(ixOmin1-1))/vA_int(ixOmin1-1))*gradvA_int(ixOmin1-1)*0.1d0*w(ixOmin1-1^%1ixO^S,wAminus_)   
      !w(ixOmax1^%1ixO^S,wAplus_) = ((w(ixOmin1-1,mom(1)) + vA_int(ixOmin1-1))/vA_int(ixOmin1-1))*gradvA_int(ixOmin1-1)*0.1d0*w(ixOmin1-1^%1ixO^S,wAminus_)    

      call uawsom_get_pthermal(w,x,ixI^L,ixIM^L,pth)

      ! Calculate Temp in the domain for use in updating density 
      !Tmp(ixOmin1^%1ixO^S) = pth(ixOmin1-1^%1ixO^S)/w(ixOmin1-1^%1ixO^S,rho_)

      h1 = x(ixOmax1,1) - x(ixOmin1,1)
       
      ! Pressure calculated by extending HSE into ghost cells 
      !Pin(ixOmin1^%1ixO^S) = pth(ixOmin1-2^%1ixO^S) + usr_grav*(SRadius/x(ixOmin1-1^%1ixO^S,1))**2*(x(ixOmin1^%1ixO^S,1)-x(ixOmin1-2^%1ixO^S,1))*w(ixOmin1-1^%1ixO^S,rho_)
      Pin(ixOmin1^%1ixO^S) = pth(ixOmin1-2^%1ixO^S) + usr_grav*(SRadius/x(ixOmin1-1^%1ixO^S,1))**2*2.d0*h1*w(ixOmin1-1^%1ixO^S,rho_)
  
      ! Old density update- Assumes isothermality at the boundary
      !w(ixOmin1^%1ixO^S,rho_) = Pin(ixOmin1^%1ixO^S)/Tmp(ixOmin1^%1ixO^S)
      !w(ixOmin1^%1ixO^S,rho_) = w(ixOmin1-1^%1ixO^S,rho_)
      w(ixOmin1^%1ixO^S,rho_) = w(ixOmin1-2^%1ixO^S,rho_) + 2*h1*usr_grav*(SRadius/x(ixOmin1-1^%1ixO^S,1))**2*w(ixOmin1-1^%1ixO^S,rho_)**2.d0/pth(ixOmin1-1^%1ixO^S)

      Tmp(ixOmin1^%1ixO^S) = Pin(ixOmin1^%1ixO^S)/w(ixOmin1^%1ixO^S,rho_)
      !Tmp(ixOmin1^%1ixO^S) = pth(ixOmin1-1^%1ixO^S)/w(ixOmin1-1^%1ixO^S,rho_)

      ! Energy in the ghost cells does include wave energy TODO: Why now ok here and not earlier?    
      w(ixOmin1^%1ixO^S,e_) = Pin(ixOmin1^%1ixO^S)/(uawsom_gamma-1) + &
                   0.5d0*(w(ixOmin1^%1ixO^S,mag(1))**2 + w(ixOmin1^%1ixO^S,mag(2))**2 + w(ixOmin1^%1ixO^S,mag(3))**2) + &
                   0.5d0*(w(ixOmin1^%1ixO^S,mom(1))**2 + w(ixOmin1^%1ixO^S,mom(2))**2 + w(ixOmin1^%1ixO^S,mom(3))**2)/w(ixOmin1^%1ixO^S,rho_) + &
                   w(ixOmin1^%1ixO^S,wkplus_) + w(ixOmin1^%1ixO^S,wkminus_) + &
                   w(ixOmin1^%1ixO^S,wAplus_) + w(ixOmin1^%1ixO^S,wAminus_)

      ! TODO: Check this, I think in Norberts original this was only done once and not twice like you have now.
      ! 10/12/24 - Checked NM original sim, this was only done once, since temp assumed isothermal I guess it does not 'need' to be done twice
      
      !Pin(ixOmax1^%1ixO^S) = pth(ixOmin1-1^%1ixO^S) + usr_grav*(SRadius/x(ixOmin1^%1ixO^S,1))**2*(x(ixOmax1^%1ixO^S,1)-x(ixOmin1-1^%1ixO^S,1))*w(ixOmin1^%1ixO^S,rho_)
      Pin(ixOmax1^%1ixO^S) = pth(ixOmin1-1^%1ixO^S) + usr_grav*(SRadius/x(ixOmin1-1^%1ixO^S,1))**2*(x(ixOmax1^%1ixO^S,1)-x(ixOmin1-1^%1ixO^S,1))*w(ixOmin1^%1ixO^S,rho_)
      
      ! Now w(ixOmin1^%1ixO^S,rho_) = w(ixImax1-1^%1ixO^S,rho_) hence the ixO notation above in Pin
      
      !Tmp(ixOmax1^%1ixO^S) = Pin(ixOmin1^%1ixO^S)/w(ixOmin1^%1ixO^S,rho_)
      
      ! Old density update
      !w(ixOmax1^%1ixO^S,rho_) = Pin(ixOmax1^%1ixO^S)/Tmp(ixOmax1^%1ixO^S)
      !w(ixOmax1^%1ixO^S,rho_) = w(ixOmin1-1^%1ixO^S,rho_)
      w(ixOmax1^%1ixO^S,rho_) = w(ixOmin1-1^%1ixO^S,rho_) + 2.d0*h1*usr_grav*(SRadius/x(ixOmin1-1^%1ixO^S,1))**2*w(ixOmin1^%1ixO^S,rho_)**2.d0/Pin(ixOmin1^%1ixO^S)

      Tmp(ixOmax1^%1ixO^S) = Pin(ixOmax1^%1ixO^S)/w(ixOmax1^%1ixO^S,rho_)
      !Tmp(ixOmax1^%1ixO^S) = Tmp(ixOmin1^%1ixO^S)
      !Tmp(ixOmax1^%1ixO^S) = Pin(ixOmin1^%1ixO^S)/w(ixOmin1^%1ixO^S,rho_)
      !Tmp(ixOmax1^%1ixO^S) = Tmp(ixOmin1^%1ixO^S)

      !if (mod(it,10000) == 0 .and. mype==7) then
      !  write(*,*) 'T vals 1st cell, 1st ghost, top ghost'
      !  write(*,*) '= ', pth(ixOmin1-1^%1ixO^S)/w(ixOmin1-1^%1ixO^S,rho_), Tmp(ixOmin1^%1ixO^S), Tmp(ixOmax1^%1ixO^S)
      !end if
      
      w(ixOmax1^%1ixO^S,e_) = Pin(ixOmax1^%1ixO^S)/(uawsom_gamma-1) + &
                   0.5d0*(w(ixOmax1^%1ixO^S,mag(1))**2 + w(ixOmax1^%1ixO^S,mag(2))**2 + w(ixOmax1^%1ixO^S,mag(3))**2) + &
                   0.5d0*(w(ixOmax1^%1ixO^S,mom(1))**2 + w(ixOmax1^%1ixO^S,mom(2))**2 + w(ixOmax1^%1ixO^S,mom(3))**2)/w(ixOmax1^%1ixO^S,rho_) + &
                   w(ixOmax1^%1ixO^S,wkplus_) + w(ixOmax1^%1ixO^S,wkminus_) + &
                   w(ixOmax1^%1ixO^S,wAplus_) + w(ixOmax1^%1ixO^S,wAminus_) 

      
      !if (mod(it,1000)==0 .and. mype==7) then 
      !  write(*,*)"e top min cell  = ", w(ixOmin1-1,e_)
      !  write(*,*)"e top min ghost = ", w(ixOmin1,e_)
      !  write(*,*)"e top max ghost = ", w(ixOmax1,e_)
      !end if

      endif
    case default
       call mpistop("Special boundary is not defined for this region")
    end select
    
  end subroutine specialbound_usr



  subroutine usrprocess_global(iit,qt)
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt

    double precision :: davgrho, avgrho_pe
    integer:: iigrid, igrid, count
   
    avgrho_pe = 0.d0
    davgrho = 0.d0
    count = 0
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       block=>ps(igrid)
       if (ps(igrid)%is_physical_boundary(1)) then
         count = count + 1
         call usrprocess_grid(ixG^LL,ixM^LL,qt,ps(igrid)%w,ps(igrid)%x,avgrho_pe)
         davgrho = davgrho + avgrho_pe
       endif
    end do
    if (count/=0) davgrho = davgrho/count
    call MPI_ALLREDUCE(davgrho,avgrho,1,MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)

  end subroutine usrprocess_global
 
  subroutine usrprocess_grid(ixI^L,ixO^L,qt,w,x,avgrho_pe)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
   
    double precision :: avgrho_pe

    avgrho_pe = w(3^%1ixO^S,rho_)
    
  end subroutine usrprocess_grid
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,1)=ggrid(ixO^S)

  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/x(ixO^S,1))**2

  end subroutine

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ix^D
    double precision :: pth(ixO^S)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)

    !if (mod(it,10000)==0 .and. mype==0) then
    !  write(*,*) "ixOmin1 = ", ixOmin1, ", ixImin1 = ", ixImin1
    !  write(*,*) "ixOmax1 = ", ixOmax1, ", ixImax1 = ", ixImax1  
    !end if
    
    !> Time dependent injection at the base 
    !do ix1 = ixOmin1,ixOmin1+1
    !  if (x(ix1,1) < 1+3.0d-5) then
    !    w(ix1,wkminus_) = 4.d-4*(0.5d0*(1-tanh((qt-0.015d0)/0.001d0))) 
    !    !w(ix1,wAminus_) = 2.d-4
    !  end if 
    !end do

    !> Constant injection at the base 
    !if (mod(it,1000)==0.and.mype==0)then
    !  do ix1 = ixOmin1,ixOmin1+1
    !    if (x(ix1,1) < 1+2.0d-4) then
    !      write(*,*) 'Wkminus = ', w(ix1,wkminus_)
    !    end if
    !  end do
    !end if

    !do ix1 = ixOmin1,ixOmin1+1
    !  if (x(ix1,1) < 1+2.0d-4) then
    !    w(ix1,wkminus_) = 1.67d-2 !+ w(ix1,wkplus_)
    !    w(ix1,wAminus_) = 1.67d-2 !+ w(ix1,wAplus_) !> Turn on or off Re-Reflected waves 
    !    w(ix1,wAminus_) = 5.d0*1.23d-2 + w(ix1,wAplus_)
    !    w(ix1,wkminus_) = 4.d-4 !+ w(ix1,wkplus_)
    !    w(ix1,wAminus_) = 2.d-4 + w(ix1,wAplus_)
    !    w(ix1,wkplus_) = 1.67d-2
    !  end if 
    !end do

    !do ix1 = ixOmax1, ixOmax1+1
    !  if (x(ix1,1) > xprobmax1-3.0d-4) then
    !    w(ix1,wkplus_) = 0.005d0*1.67d-2
    !  end if
    !end do

    !do ix1 = ixOmax1, ixOmax1+1
    !  if (x(ix1,1) > xprobmax1-3.0d-4) then
    !    w(ix1,wAplus_) = 1.d-7*1.67d-2 + 0.01d0*w(ix1,wAminus_)
    !    w(ix1,wAplus_) = 0.d0
    !  end if
    !end do

    !> data_velocity_reduction_layer_0_100_domain_SuperGaussian simulation
    !do ix1 = ixOmin1, ixOmax1
    !  if (x(ix1,1) > 1.1d0) then
    !    w(ix1,mom(1)) = w(ix1,mom(1))*(1 - exp(-((x(ix1,1) - 1.12d0)/0.01d0)**8.d0))
    !  end if
    !end do

    !w(ixO^S,wkminus_) = w(ixO^S,wkminus_) + 0.01d0*qdt*exp(-(x(ixO^S,1)-1.0d0)**2/(0.002d0)**2)

    !w(ixO^S,wkminus_) = w(ixO^S,wkminus_) + 10.0d0*qdt*exp(-(x(ixO^S,1)-1.05d0)**2/(0.02d0)**2)*(0.5d0*(1-tanh((qt-0.002d0)/(0.01d0)**2)))
    !w(ixO^S,wkminus_) = w(ixO^S,wkminus_) + qdt*(0.5d0*(1-tanh(qt-3.d0)))*0.0125d0*exp(-(x(ixO^S,1)-1.d0)**2.d0/0.05d0**2.d0)
    !w(ixO^S,wkminus_) = w(ixO^S,wkminus_) + qdt*1.d0*exp(-(x(ixO^S,1)-1.0d0)**2/(0.002d0)**2)
    !w(ixO^S,wkminus_) = w(ixO^S,wkminus_) + 0.5d0*(1+tanh((qt-2.d0)/0.1d0))*qdt*exp(-(x(ixO^S,1)-1.02d0)**2/(0.02d0)**2)
    !w(ixO^S,wkminus_) = w(ixO^S,wkminus_) + qdt*1.d0*exp(-(x(ixO^S,1)-xprobmin1)/(Heatscale/unit_length))

    !w(ixO^S,wAminus_) = w(ixO^S,wAminus_) + 0.0d0*qdt*1.d0*exp(-(x(ixO^S,1)-1.03d0)**2/(0.02d0)**2)
    !w(ixO^S,wAminus_) = w(ixO^S,wAminus_) + qdt*2.d0*exp(-(x(ixO^S,1)-1.005d0)**2/(0.04d0)**2)
    !w(ixO^S,wAminus_) = w(ixO^S,wAminus_) + qdt*1.d0*exp(-(x(ixO^S,1)-xprobmin1)/(Heatscale/unit_length))

    !> Velocity rewrite layer for inflows from top boundary 
    !do ix1 = ixOmin1, ixOmax1
    !  if (x(ix1,1)>1.14d0) then
    !    w(ix1^%1ixO^S,mom(1)) = 0.d0
    !  end if
    !end do

    !do ix1 = ixOmin1, ixOmax1
    !  if (x(ix1,1)>1.14d0) then
    !    w(ix1^%1ixO^S,mom(1)) = 0.5d0*(1+tanh((qt-2.d0)/0.1))*w(ix1^%1ixO^S,mom(1))
    !  end if
    !end do

  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    use mod_random
  ! calculate background heating bQ
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S), pth(ixI^S), Tmp(ixO^S), tau(ixO^S)
    double precision, save :: qt0 = 0.0, fqt = 0.0, rr(50)

    !integer :: i
    !logical, save :: firstqt=.true.
 
    !if (floor(500*fqt+1) .ne. floor(500*qt+1)) then
    !   firstqt = .true.
    !endif
    
    !fqt = qt
    
    !if (firstqt) then
    ! call rng%unif_01_vec(rr)
    ! rr = rr*2.d0+1.d0
    ! firstqt = .false.
    ! qt0 = qt
    !endif

    bQgrid(ixO^S) = 0.5d0*bQ0*dexp(-(x(ixO^S,1)-xprobmin1)/(Heatscale/unit_length))*0.5d0*(1-tanh((qt-2.5d0)/0.25d0))
    !bQgrid(ixO^S) = 0.5d0*bQ0*dexp(-(x(ixO^S,1)-xprobmin1)/(0.01d0))
    
    !call uawsom_get_pthermal(w,x,ixI^L,ixI^L,pth)
    !Tmp(ixO^S) = pth(ixO^S)/w(ixO^S,rho_)
    
    !where(Tmp(ixO^S) < 1.d4/unit_temperature) 
    !  tau(ixO^S) = exp(-(x(ixO^S,1)-1.d0)/Hr)
    !  bQgrid(ixO^S) = 5.d-2*(1.d0 - exp(-tau(ixO^S)/0.1d0))*4.d0*tau(ixO^S)/Hr*stefboltz*(0.75d0*0.0058d0**4*(tau(ixO^S) +&
    !                  0.71d0-0.133d0/(1.d0 + 0.15d0*tau(ixO^S)**0.73d0)**17.4d0) - Tmp(ixO^S)**4.d0)
    !endwhere
 
    !bQgrid(ixO^S)=bQ0*dexp(-(x(ixO^S,1)-xprobmin1)/(Heatscale/unit_length))*&
    !                  0.5d0*(1.d0-tanh(qt-3.d0))
    
    !do i=1,50
    !  bQgrid(ixO^S) = bQgrid(ixO^S) + 1.d-2*dexp(-(x(ixO^S,1)-xprobmin1))*dexp(-((x(ixO^S,1)-rr(i))**2)/(0.01d0)**2)*&
    !                  exp(-(qt-qt0-(2.d-3*(rr(i)-1.25d0)))**2/4.d-6) 
    !end do

  end subroutine getbQ

  subroutine get_zeta(w,x,ixI^L,ixO^L,zeta)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: zeta(ixI^S)
    double precision :: zeta0 = 5.0d0
    
    !zeta(ixI^S) = zeta0*exp(-(x(ixI^S,1)-xprobmin1)/5.d0)
    zeta(ixI^S) = (zeta0-1.d0)*exp(-(x(ixI^S,1)-xprobmin1)/5.d0)+1.d0
    !where(zeta(ixI^S) < 1.d0)
    !  zeta(ixI^S) = 1.d0
    !end where

  end subroutine get_zeta

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_radiative_cooling
    use mod_uawsom_phys
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                :: ggrid(ixI^S)

    double precision :: pth(ixI^S),B2(ixI^S),tmp2(ixI^S),dRdT(ixI^S)
    double precision :: ens(ixI^S),divb(ixI^S),wlocal(ixI^S,1:nw)
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir), zeta(ixI^S), radius(ixI^S), Lperp_AW(ixI^S), Lperp(ixI^S), Gamma_plus(ixI^S), Gamma_minus(ixI^S)
    double precision :: Te(ixI^S),tco_local
    double precision, dimension(ixI^S,1:ndim) :: gradT, bunitvec
    integer :: idirmin,idir,ix^D
    logical :: lrlt(ixI^S)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call uawsom_get_pthermal(wlocal,x,ixI^L,ixI^L,pth)
    Te(ixI^S)=pth(ixI^S)/w(ixI^S,rho_)
    w(ixO^S,nw+1)=Te(ixO^S)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0)
      else
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))
      endif
    end do
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
      
    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+2)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))

    ! output divB1
    call get_divb(wlocal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)
    ! output heating rate
    call getbQ(ens,ixI^L,ixO^L,global_time,wlocal,x)
    w(ixO^S,nw+5)=ens(ixO^S)
    ! store the cooling rate 
    if(uawsom_radiative_cooling)call getvar_cooling(ixI^L,ixO^L,wlocal,x,ens,rc_fl)
    w(ixO^S,nw+6)=ens(ixO^S)

    ! store current
    call get_current(wlocal,ixI^L,ixO^L,idirmin,curlvec)
    do idir=1,ndir
      w(ixO^S,nw+6+idir)=curlvec(ixO^S,idir)
    end do

    w(ixO^S,nw+10) = w(ixO^S,mag(1)) !> perturbation: always zero in 1D as are the others
    w(ixO^S,nw+11) = w(ixO^S,mag(2))
    w(ixO^S,nw+12) = w(ixO^S,mag(3))
    
    w(ixO^S,nw+13) = w(ixO^S,mom(1))/w(ixO^S,rho_)
    w(ixO^S,nw+14) = w(ixO^S,mom(2))/w(ixO^S,rho_)
    w(ixO^S,nw+15) = w(ixO^S,mom(3))/w(ixO^S,rho_)

    !call getggrav(ggrid,ixI^L,ixO^L,x)

    Lperp_AW(ixO^S) = (1/6.961d10) * 1.5d5 * 1.0d2 * 1.0d2 * (Te(ixO^S) / Btotal(ixO^S,1))**0.5d0 
    Gamma_plus(ixO^S) = (2.0d0 / Lperp_AW(ixO^S)) * (w(ixO^S, wAminus_)/w(ixO^S,rho_))**0.5d0
    Gamma_minus(ixO^S) = (2.0d0 / Lperp_AW(ixO^S)) * (w(ixO^S, wAplus_)/w(ixO^S,rho_))**0.5d0

    w(ixO^S,nw+16) = Gamma_plus(ixO^S)*w(ixO^S,wAplus_) + Gamma_minus(ixO^S)*w(ixO^S,wAminus_)

    !ff = 0.1d0
    
    call get_zeta(w,x,ixI^L,ixO^L,zeta)

    radius(ixO^S) = 1.d8/unit_length * ((Busr/unit_magneticfield)/Btotal(ixO^S,1))**0.5d0 !Radius = R_0*(B0/B)^0.5 TVD 2025 paper uses R_0 = 1Mm (1e8 cm)

    Lperp(ixO^S) = (zeta(ixO^S) + 1.d0 - ff)**(3.d0/2.d0)/(1.d0 - ff**(5.d0/2.d0))/&
                   (zeta(ixO^S) - 1.d0)*3.1622776*(ff*dpi)**0.5d0*radius(ixO^S) 
                
    !if (mype==0 .and. mod(it,10000)==0) then  
    !  write(*,*) "Lperp = ", Lperp(ixOmin1:ixOmin1+10)
    !end if

    w(ixO^S,nw+17) = w(ixO^S,wkplus_)**(3.d0/2.d0)/(w(ixO^S,rho_)*(1+ff*zeta(ixO^S)-ff)**(-1.d0))**0.5d0/Lperp(ixO^S) + w(ixO^S,wkminus_)**(3.d0/2.d0)/(w(ixO^S,rho_)*(1+ff*zeta(ixO^S)-ff)**(-1.d0))**0.5d0/Lperp(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
      
    varnames='Te Alfv divB beta bQ rad j1 j2 j3 br bphi bth vr vphi vth Qaw Qk'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a time-independent background magnetic field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,1)=B0*x(ixO^S,1)**(-2.d0)
    wB0(ixO^S,2)=0.d0
    wB0(ixO^S,3)=0.d0

  end subroutine specialset_B0
  
  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)

    wJ0(ixO^S,1)= 0.d0
    wJ0(ixO^S,2)= 0.d0
    wJ0(ixO^S,3)= 0.d0

  end subroutine specialset_J0
end module mod_usr
