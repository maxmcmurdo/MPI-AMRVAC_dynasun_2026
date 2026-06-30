! 03/08/2022 - Copying from mhd/solar_atmosphere_2.5D.t
! 05/10/2022 - Setting the B field splitting with 'linedpowel'

module mod_usr
  use mod_uawsom
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav, trelax
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,dya  !bQ0,
  double precision :: y0, Bpo, Bzo
  double precision, allocatable :: pa(:),ra(:)
  integer, parameter :: jmax=8000

  ! Storing additional variables in the dat file
  integer :: Temp_, b1t_, b2t_, b3t_, j3_
  integer :: Alfv_, vsound_, divB_, beta_
  integer :: bQ_, lQ_, rad_, thcond_
  integer :: Qk_, QAW_ 

contains

  subroutine usr_init()

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_set_B0          => specialset_B0
    usr_modify_output   => set_output_vars

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9 !cm-3,cm-3

    call set_coordinate_system("Cartesian_2.5D")
    call uawsom_activate()

    Temp_ = var_set_extravar("Te", "Te")
    b1t_ = var_set_extravar("b1t", "b1t")
    b2t_ = var_set_extravar("b2t", "b2t")
    b3t_ = var_set_extravar("b3t", "b3t")
    j3_ = var_set_extravar("j3", "j3")

    Alfv_ = var_set_extravar("Alfv", "Alfv")
    vsound_ = var_set_extravar("vsound", "vsound")
    divB_ = var_set_extravar("divB", "divB")
    beta_ = var_set_extravar("beta", "beta")

    bQ_ = var_set_extravar("bQ", "bQ")
    lQ_ = var_set_extravar("lQ", "lQ")
    rad_ = var_set_extravar("rad", "rad")
    thcond_ = var_set_extravar("thcond", "thcond")
    Qk_ = var_set_extravar("Qk", "Qk")
    QAW_ = var_set_extravar("QAW", "QAW")

  end subroutine usr_init

  subroutine initglobaldata_usr()
    heatunit=unit_pressure/unit_time !3.697693390805347E-003 erg*cm-3/s,erg*cm-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    ! bQ0=1.d-4/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) !cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    Bpo = B0     !/(dsqrt(2.d0)*(dexp(-1.d0)-dexp(-3.d0)))
    Bzo = B0     !/dsqrt(2.d0)
    y0 = -0.4d0 !=4Mm as in Zhang+2019
    kx=dpi/(xprobmax1-xprobmin1)
    SRadius=69.61d0 ! Solar radius
    trelax=90.d0

    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic

    if(mype .eq. 0) then
      print*, 'unit_density = ', unit_density
      print*, 'unit_pressure = ', unit_pressure
      print*, 'unit_velocity = ', unit_velocity
      print*, 'unit_magneticfield = ', unit_magneticfield
      print*, 'unit_time = ', unit_time
      print*, 'usr_grav = ', usr_grav
    end if

  end subroutine initglobaldata_usr

  subroutine inithdstatic
    use mod_solar_atmosphere
    ! initialize the table in a vertical line through the global domain
    integer :: j,na,ibc
    double precision, allocatable :: Ta(:),gg(:),ya(:)
    double precision :: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa
    double precision :: rhohc,hc

    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))

    rpho=1.151d15/unit_numberdensity !number density at the bottom of height table
    Tpho=8.d3/unit_temperature ! temperature of chromosphere
    Ttop=1.5d6/unit_temperature ! estimated temperature in the top
    htra=0.2d0 ! height of initial transition region
    wtra=0.02d0 ! width of initial transition region
    Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
    Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**&
       3
    do j=1,jmax
       ya(j)=(dble(j)-0.5d0)*dya-gzone
       if(ya(j)>htra) then
         Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
       else
         Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
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
    deallocate(ya,gg,Ta)

    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    ! initialize one grid
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: MAX0wB02(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        MIN0wB02(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision :: res
    integer :: ix1,ix2,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D solar atmosphere'
      endif
      first=.false.
    endif
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
        na=floor((x(ix1,ix2,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix1,ix2,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix1,ix2,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix1,ix2,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    end do
    end do
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(:))=zero
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))=zero

    MAX0wB02(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = max(0.0d0,&
       -Bpo*dsin(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1))*dexp(-kx*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2)-y0)) + Bpo*dsin(3*kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1))*dexp(-3*kx*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)-y0)))
    MIN0wB02(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = min(0.0d0,&
       -Bpo*dsin(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1))*dexp(-kx*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2)-y0)) + Bpo*dsin(3*kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1))*dexp(-3*kx*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)-y0)))

    !w(ixO^S,wkminus_) = 1.d-6 + 5.0d-1*(0.0d0+1.0d0*MAX0wB02(ixO^S))*dexp(-(x(ixO^S,2)-(xprobmin2))) !1+0.1Max seems fine !so is 0.02+1Max
    !w(ixO^S,wkplus_) =  1.d-6 - 1.0d-1*(-0.0d0+1.0d0*MIN0wB02(ixO^S))*dexp(-(x(ixO^S,2)-(xprobmin2)))

    !w(ixO^S,wAminus_) = 1.d-6 + 5.0d-1*(0.0d0+1.0d0*MAX0wB02(ixO^S))*dexp(-(x(ixO^S,2)-(xprobmin2))) 
    !w(ixO^S,wAplus_) =  1.d-6 - 1.0d-1*(-0.0d0+1.0d0*MIN0wB02(ixO^S))*dexp(-(x(ixO^S,2)-(xprobmin2)))

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,wkminus_) = 1.d-6 + &
       0.7d0*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2)-(xprobmin2)))*merge(1.d0, 0.d0, MAX0wB02(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) > 0.d0)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,wkplus_) = 1.d-6 + &
       0.7d0*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2)-(xprobmin2)))*merge(1.d0, 0.d0, MIN0wB02(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) < 0.d0)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,wAminus_) = 1.d-6 + &
       0.7d0*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2)-(xprobmin2)))*merge(1.d0, 0.d0, MAX0wB02(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) > 0.d0)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,wAplus_) = 1.d-6 + &
       0.7d0*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2)-(xprobmin2)))*merge(1.d0, 0.d0, MIN0wB02(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) < 0.d0)

    call uawsom_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x)

  end subroutine initonegrid_usr

  subroutine specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x,wB0)
  ! Here add a steady (time-independent) potential or
  ! linear force-free background field
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)

    wB0(ixImin1:ixImax1,ixImin2:ixImax2,1)=Bpo*dcos(kx*x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1))*dexp(-kx*(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       2)-y0)) - Bpo*dcos(3*kx*x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1))*dexp(-3*kx*(x(ixImin1:ixImax1,ixImin2:ixImax2,2)-y0))
    wB0(ixImin1:ixImax1,ixImin2:ixImax2,2)=-Bpo*dsin(kx*x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1))*dexp(-kx*(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       2)-y0)) + Bpo*dsin(3*kx*x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1))*dexp(-3*kx*(x(ixImin1:ixImax1,ixImin2:ixImax2,2)-y0))
    wB0(ixImin1:ixImax1,ixImin2:ixImax2,3)=Bzo

  end subroutine specialset_B0

  subroutine specialbound_usr(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, iB, ixImin1,&
       ixImin2,ixImax1,ixImax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2),ggrid(ixImin1:ixImax1,&
       ixImin2:ixImax2),invT(ixImin1:ixImax1,ixImin2:ixImax2),&
        MAX0wB02G(ixImin1:ixImax1,ixImin2:ixImax2),MIN0wB02G(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: Q(ixImin1:ixImax1,ixImin2:ixImax2),Qp(ixImin1:ixImax1,&
       ixImin2:ixImax2),zeta(ixImin1:ixImax1,ixImin2:ixImax2)
    integer          :: ix1,ix2,ixOsmin1,ixOsmin2,ixOsmax1,ixOsmax2,ixCmin1,&
       ixCmin2,ixCmax1,ixCmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2,jxOmin1,jxOmin2,&
       jxOmax1,jxOmax2,idir
    double precision :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    MAX0wB02G(ixImin1:ixImax1,ixImin2:ixImax2) = max(0.0d0,&
       -Bpo*dsin(kx*x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1))*dexp(-kx*(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       2)-y0)) + Bpo*dsin(3*kx*x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1))*dexp(-3*kx*(x(ixImin1:ixImax1,ixImin2:ixImax2,2)-y0)))
    MIN0wB02G(ixImin1:ixImax1,ixImin2:ixImax2) = min(0.0d0,&
       -Bpo*dsin(kx*x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1))*dexp(-kx*(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       2)-y0)) + Bpo*dsin(3*kx*x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1))*dexp(-3*kx*(x(ixImin1:ixImax1,ixImin2:ixImax2,2)-y0)))

    ! iB	integer indicating direction of boundary
    select case(iB)
    case(3)
      !! Fixed zero velocity:
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=-w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))/w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do

      do ix1=ixOmin1,ixOmax1
        w(ix1,ixOmin2,wkminus_) = 0.7d0*merge(1.d0, 0.d0, MAX0wB02G(ix1,&
           ixOmin2) > 0.d0) !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2)) 
        w(ix1,ixOmax2,wkminus_) = 0.7d0*merge(1.d0, 0.d0, MAX0wB02G(ix1,&
           ixOmin2) > 0.d0) !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2)) 

        w(ix1,ixOmin2,wkplus_) =  0.7d0*merge(1.d0, 0.d0, MIN0wB02G(ix1,&
           ixOmin2) < 0.d0) !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))
        w(ix1,ixOmax2,wkplus_) =  0.7d0*merge(1.d0, 0.d0, MIN0wB02G(ix1,&
           ixOmin2) < 0.d0) !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))

        w(ix1,ixOmin2,wAminus_) = 0.7d0*merge(1.d0, 0.d0, MAX0wB02G(ix1,&
           ixOmin2) > 0.d0) !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2)) 
        w(ix1,ixOmax2,wAminus_) = 0.7d0*merge(1.d0, 0.d0, MAX0wB02G(ix1,&
           ixOmin2) > 0.d0) !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2)) 
        
        w(ix1,ixOmin2,wAplus_) =  0.7d0*merge(1.d0, 0.d0, MIN0wB02G(ix1,&
           ixOmin2) < 0.d0) !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))
        w(ix1,ixOmax2,wAplus_) =  0.7d0*merge(1.d0, 0.d0, MIN0wB02G(ix1,&
           ixOmin2) < 0.d0) !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))

      end do

      !! fixed b1 b2 b3
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))=0.d0

      !! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
      enddo
      call uawsom_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)
    case(4)
      ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOsmin1,&
         ixOsmin2,ixOsmax1,ixOsmax2,pth)
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmax2;
      call getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOsmin1,ixOsmin2,&
         ixOsmax1,ixOsmax2,x)

      !> Fill pth, rho ghost layers according to gravity stratification:
      invT(ixOmin1:ixOmax1,ixOmin2-1)=w(ixOmin1:ixOmax1,ixOmin2-1,&
         rho_)/pth(ixOmin1:ixOmax1,ixOmin2-1)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin1:ixOmax1,ixOmin2-1)=tmp(ixOmin1:ixOmax1,&
           ixOmin2-1)+0.5d0*(ggrid(ixOmin1:ixOmax1,ix2)+ggrid(ixOmin1:ixOmax1,&
           ix2-1))*invT(ixOmin1:ixOmax1,ixOmin2-1)
        w(ixOmin1:ixOmax1,ix2,p_)=pth(ixOmin1:ixOmax1,&
           ixOmin2-1)*dexp(tmp(ixOmin1:ixOmax1,ixOmin2-1)*dxlevel(2))
        w(ixOmin1:ixOmax1,ix2,rho_)=w(ixOmin1:ixOmax1,ix2,&
           p_)*invT(ixOmin1:ixOmax1,ixOmin2-1)
      enddo

      !> Fixed zero velocity:
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir)) =-w(ixOmin1:ixOmax1,&
           ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))/w(ixOmin1:ixOmax1,&
           ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do

      !> Magnetic field:
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* (-w(ixOmin1:ixOmax1,ix2-2,&
           mag(:))+4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
      enddo
      call uawsom_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)

    case default
       call mpistop("Special boundary is not defined for this region")

    end select   

  end subroutine specialbound_usr

  subroutine gravity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,wCT,x,gravity_field)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
       ixImin2:ixImax2,ndim)

    double precision                :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2)

    gravity_field=0.d0
    call getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,x)
    gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=ggrid(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine gravity

  subroutine getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2)

    ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=usr_grav*(SRadius/(SRadius+&
       x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))**2
  end subroutine

  subroutine special_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    use mod_global_parameters
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: ix1,ix2,var
    
    double precision :: lQgrid(ixImin1:ixImax1,ixImin2:ixImax2),&
       bQgrid(ixImin1:ixImax1,ixImin2:ixImax2), zeta(ixImin1:ixImax1,&
       ixImin2:ixImax2), dx_inj, MAX0wB02(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        MIN0wB02(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    ! add global background heating bQ
    !call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    !w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    ! add steady localized heating
    !call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
    !w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)

  end subroutine special_source

  subroutine get_zeta(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,zeta)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
       x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) :: zeta(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: zeta0 = 5.d0 !> Changed the density contrast denominator to be scaled with R_suns not unit_lengths, here: 10Mm 
    
    zeta(ixImin1:ixImax1,ixImin2:ixImax2) = &
       (zeta0-xprobmin2)*exp(-(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       2)-xprobmin2)/(5.d0*69.61d0))+1.d0 !> MAX: Open field MAX TODO: Should we have a better representation for this expansion with distance along a field line rather than just height? (for 2D+)
    !zeta(ixI^S) = (zeta0-1.d0)*exp(-((xprobmax1-xprobmin1)/dpi * sin(dpi * (x(ixI^S,1)-xprobmin1)/(xprobmax1-xprobmin1)))/5.d0)+1.d0 !> MAX: 1D Loops
  
  end subroutine get_zeta

  subroutine getbQ(bQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,w,x)
    use mod_global_parameters
  ! calculate background heating bQ
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: bQgrid(ixImin1:ixImax1,ixImin2:ixImax2),bQ0

    !> temporal distribution
    !> initially a little bit higher to prevent condensation
    if(qt .lt. 10) then
      bQ0=4.d-4/heatunit
    else if(qt .lt. 20) then
      bQ0=3.d-4/heatunit
    else if(qt .lt. 30) then
      bQ0=2.d-4/heatunit
    else
      bQ0=1.d-4/heatunit
    end if

    ! spatial distribution
    bQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=bQ0*dexp(-x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)/5.d0)

  end subroutine getbQ

  subroutine getlQ(lQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,w,x)
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: ix1,ix2,i
    double precision :: lQgrid(ixImin1:ixImax1,ixImin2:ixImax2),lQ0,tramp
    double precision :: yh,lambdah,sigma2,xr,xl

    lQ0 = 2.d-2/heatunit
    tramp = 10.d0

    yh = 0.4d0
    lambdah = 0.25d0
    sigma2 = 0.2d0

    xr = 4.2d0
    xl = -xr

    if(qt .le. trelax) then
      lQgrid=zero
    else if(qt .le. trelax+tramp) then
      lQgrid=dsin((qt-trelax)/tramp*dpi/2.d0)
    else
      lQgrid=one
    end if

    where (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) .lt. yh)
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = lQ0*lQgrid(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*(dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1)-xr)**2/sigma2) + dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1)-xl)**2/sigma2))
    elsewhere
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = lQ0*lQgrid(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2)-yh)**2/lambdah)* (dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1)-xr)**2/sigma2) + dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1)-xl)**2/sigma2))
    end where

  end subroutine getlQ

  subroutine special_refine_grid(igrid,level,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,w,x,refine,coarsen)
    integer, intent(in) :: igrid, level, ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: qt, w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if(any(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) .ge. 5.d0)) then
      refine=1
      coarsen=-1
    end if

  end subroutine special_refine_grid

  subroutine set_output_vars(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,w,x)
    use mod_radiative_cooling
    use mod_thermal_conduction
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: wlocal(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
       Te(ixImin1:ixImax1,ixImin2:ixImax2),pth(ixImin1:ixImax1,&
       ixImin2:ixImax2),B2(ixImin1:ixImax1,ixImin2:ixImax2),&
       Btotal(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
    double precision :: Bmag(ixImin1:ixImax1,ixImin2:ixImax2),&
       pmag(ixImin1:ixImax1,ixImin2:ixImax2),divb(ixImin1:ixImax1,&
       ixImin2:ixImax2),curlvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
    double precision :: rc(ixImin1:ixImax1,ixImin2:ixImax2),tc(ixImin1:ixImax1,&
       ixImin2:ixImax2),ens(ixImin1:ixImax1,ixImin2:ixImax2),&
       loc_heat(ixImin1:ixImax1,ixImin2:ixImax2)
    integer          :: idir,idirmin
    double precision :: zeta(ixImin1:ixImax1,ixImin2:ixImax2),&
        radius(ixImin1:ixImax1,ixImin2:ixImax2), Lperp_AW(ixImin1:ixImax1,&
       ixImin2:ixImax2), Lperp(ixImin1:ixImax1,ixImin2:ixImax2),&
        Gamma_plus(ixImin1:ixImax1,ixImin2:ixImax2),&
        Gamma_minus(ixImin1:ixImax1,ixImin2:ixImax2)

    wlocal(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    !   output temperature
    call uawsom_get_pthermal(wlocal,x,ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,&
       ixImin2,ixImax1,ixImax2,pth)
    Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Temp_)=Te(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixImin1:ixImax1,ixImin2:ixImax2,idir)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,mag(idir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
           idir,0)
      else
        Btotal(ixImin1:ixImax1,ixImin2:ixImax2,idir)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,mag(idir))
      endif
    end do
    !   store total magnetic field
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1t_)=Btotal(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2t_)=Btotal(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3t_)=Btotal(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,3)

    !   B^2
    B2(ixImin1:ixImax1,ixImin2:ixImax2)=sum((Btotal(ixImin1:ixImax1,&
       ixImin2:ixImax2,:))**2,dim=ndim+1)
    Bmag(ixImin1:ixImax1,ixImin2:ixImax2)=sqrt(B2(ixImin1:ixImax1,&
       ixImin2:ixImax2))
    pmag(ixImin1:ixImax1,ixImin2:ixImax2)=B2(ixImin1:ixImax1,&
       ixImin2:ixImax2)/2.0d0

    !   store current
    call get_current(wlocal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,idirmin,curlvec)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,j3_)=curlvec(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,3)

    !    output Alfven wave speed B/sqrt(rho)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Alfv_)=dsqrt(B2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

    !    output the sound speed sqrt(gamm p/rho)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,vsound_)= sqrt(uawsom_gamma * &
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_))

    !   output divB1
    call get_divb(wlocal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,divb)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,divB_)=divb(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    ! output the plasma beta p*2/B**2
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,beta_)=pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*two/B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !  store the cooling rate
    if(uawsom_radiative_cooling) call getvar_cooling(ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wlocal,x,rc,rc_fl)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rad_)=rc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    ! output thermal conduction TC
    !tc_fl%e_=1
    !call uawsom_sts_set_source_tc_mhd(ixI^L,ixO^L,wlocal,x,tc,.false.,1.d0,1,1,tc_fl)
    !w(ixO^S,thcond_)=tc(ixO^S)
    !tc_fl%e_=e_

    !    output heating rate
    call getbQ(ens,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,global_time,wlocal,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bQ_)=ens(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    !   output local heating
    call getlQ(loc_heat,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,global_time,wlocal,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,lQ_)=loc_heat(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    !   output kink wave heating rate
    call get_zeta(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,zeta)
    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1.d8/unit_length * &
       ((Busr/unit_magneticfield)/dsqrt((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       b1t_)**2.d0 + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       b2t_)**2.d0 + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3t_)**2.d0)))**0.5d0 !Radius = R_0*(B0/B)0.5,R_0*(B0/B)0.5 TVD 2025 paper uses R_0 = 1Mm (1e8 cm)
    Lperp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (zeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + 1.d0 - ff)**(3.d0/2.d0)/(1.d0 - &
       ff**(5.d0/2.d0))/(zeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - 1.d0)*3.1622776*(ff*dpi)**&
       0.5d0*radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Qk_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       wkplus_)**(3.d0/2.d0)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)*(1+ff*zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-ff)**(-&
       1.d0))**0.5d0/Lperp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) +w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,wkminus_)**(3.d0/2.d0)/(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)*(1+ff*zeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-ff)**(-1.d0))**0.5d0/Lperp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)  
                  

    !   output Alfven wave heating rate
    Lperp_AW(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (1/unit_length*10**2) * 1.5d5 &
       * 1.0d2 * (Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2) / Btotal(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1))**0.5d0 
    Gamma_plus(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (2.0d0 / &
       Lperp_AW(ixOmin1:ixOmax1,ixOmin2:ixOmax2)) * (w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, wAminus_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_))**0.5d0
    Gamma_minus(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (2.0d0 / &
       Lperp_AW(ixOmin1:ixOmax1,ixOmin2:ixOmax2)) * (w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, wAplus_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_))**0.5d0
    
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,QAW_) = Gamma_plus(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       wAplus_) + Gamma_minus(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,wAminus_)

  end subroutine set_output_vars
end module mod_usr
