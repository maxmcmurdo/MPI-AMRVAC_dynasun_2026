!> Magnetic bipolar field
module mod_usr
  use mod_uawsom
  implicit none
  integer, save :: nx1,nx2,nxbc1,nxbc2,nxbc3
  integer, parameter :: jmax=8000
  double precision, allocatable :: pbc(:),rbc(:)
  ! 1D solar atmosphere table for pressure, density, and height
  double precision :: pa(jmax),ra(jmax),ya(jmax)
  double precision :: usr_grav,SRadius,rhob,Tiso,dr,gzone,bQ0
  double precision :: q_para,d_para,L_para, charge1_x(3), charge2_x(3),&
      charge1, charge2

contains

  !==============================================================================
  ! Purpose: to include global parameters, set user methods, set coordinate 
  !          system and activate physics module.
  !==============================================================================
  subroutine usr_init()

    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm-3,cm-3,cm-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => specialsource
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_init_vector_potential=>initvecpot_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_set_electric_field => driven_electric_field
    usr_set_B0          => specialset_B0
    usr_set_J0          => specialset_J0

    call set_coordinate_system("Cartesian_3D")
    call uawsom_activate()

  end subroutine usr_init

  !==============================================================================
  ! Purpose: to initialize user public parameters and reset global parameters.
  !          Input data are also read here.
  !==============================================================================
  subroutine initglobaldata_usr()
    double precision :: x0, y0, z0, r0, h, theta, L, Bperp
    integer :: ixp, nphalf

    !unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    !unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    !unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    !unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    !unit_time          = unit_length/unit_velocity                     ! 85.8746159942810 s 
    if(.not.uawsom_energy) then
      ! bottom density
      rhob=2.d0
      ! isothermal uniform temperature
      Tiso= uawsom_adiab
    end if

    bQ0=1.d-4/unit_pressure*unit_time !3.697693390805347E-003 erg*cm-3/s,erg*cm-3/s,erg*cm-3/s

    gzone=2.d8/unit_length
    ! cell size in 1D solar atmosphere table
    dr=(2.d0*gzone+xprobmax3-xprobmin3)/dble(jmax)
    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius

    !q_para=7.d19/(unit_magneticfield*unit_length**2) ! strength and sign of magnetic charges
    q_para=Busr*2.10025d18/(unit_magneticfield*unit_length**2) !strength and sign of magnetic charges
    d_para=1.d9/unit_length ! depth of magnetic charges
    L_para=1.5d9/unit_length ! half distance between magnetic charges

    charge1=-q_para
    charge1_x(1)=-L_para
    charge1_x(2)=0.d0
    charge1_x(3)=-d_para

    charge2=q_para
    charge2_x(1)=L_para
    charge2_x(2)=0.d0
    charge2_x(3)=-d_para

    if(uawsom_energy) call inithdstatic

  end subroutine initglobaldata_usr

  !> initialize solar atmosphere table in a vertical line through the global domain
  subroutine inithdstatic
    use mod_global_parameters

    double precision :: Ta(jmax),gg(jmax)
    double precision :: rpho,Tpho,Ttop,htra,wtra,ftra,Ttr,Fc,k_para
    double precision :: res,pb,rhob,invT
    integer :: j,na,ibc,btlevel

    rpho=1.151d15/unit_numberdensity ! number density at the bottom relaxla
    Tpho=8.d3/unit_temperature ! temperature of chromosphere
    Ttop=1.5d6/unit_temperature ! estimated temperature in the top
    htra=2.d8/unit_length ! height of initial transition region
    wtra=2.d7/unit_length ! width of initial transition region 
    htra=2.d8/unit_length ! height of initial transition region
    wtra=0.2d8/unit_length! width of initial transition region 
    Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
    Fc=2.d5/unit_pressure/unit_velocity  ! constant thermal conduction flux
    ftra=wtra*atanh(2.d0*(Ttr-Tpho)/(Ttop-Tpho)-1.d0)
    ! Spitzer thermal conductivity with cgs units
    k_para=8.d-7*unit_temperature&
       **3.5d0/unit_length/unit_density/unit_velocity**3 
    !! set T distribution with height
    do j=1,jmax
       ya(j)=(dble(j)-0.5d0)*dr-gzone
       if(ya(j)>htra) then
         Ta(j)=(3.5d0*Fc/k_para*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
       else
         Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra+ftra)/wtra)+1.d0)
       endif
       gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    enddo
    !! solution of hydrostatic equation 
    ra(1)=rpho
    pa(1)=rpho*Tpho
    !invT=0.d0
    !do j=2,jmax
    !   invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
    !   pa(j)=pa(1)*dexp(invT*dr)
    !   ra(j)=pa(j)/Ta(j)
    !end do
    do j=2,jmax
       pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dr*(gg(j)+&
          gg(j-1))/Ta(j)/4.d0)
       ra(j)=pa(j)/Ta(j)
    end do
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dr+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dr
    rhob=ra(na)+res/dr*(ra(na+1)-ra(na))
    pb=pa(na)+res/dr*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    btlevel=refine_max_level
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(3,btlevel)*(dble(nghostcells-ibc+1)-0.5d0))/dr+0.5d0)
      res=gzone-dx(3,btlevel)*(dble(nghostcells-ibc+&
         1)-0.5d0)-(dble(na)-0.5d0)*dr
      rbc(ibc)=ra(na)+res/dr*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dr*(pa(na+1)-pa(na))
    end do
    
    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'maxTa',Ta(jmax)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: Bfr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    double precision :: res
    integer :: ix1,ix2,ix3,na
    logical, save :: first=.true.
    double precision :: MAX0wB03(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), MIN0wB03(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if(first)then
      if(mype==0) then
      write(*,*)'Stable solar atmosphere with a bipolar magnetic field'
    endif
      first=.false.
    endif

    if(B0field) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))=zero
    else 
      if(stagger_grid) then
        call b_from_vector_potential(block%ixGsmin1,block%ixGsmin2,&
           block%ixGsmin3,block%ixGsmax1,block%ixGsmax2,block%ixGsmax3,ixImin1,&
           ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
           ixOmax1,ixOmax2,ixOmax3,block%ws,x)
        call uawsom_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3,block)
      else
        call bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,A,Bfr)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(:))=Bfr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)
      end if
    end if
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:))=0.d0

    if(uawsom_energy) then
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
         na=floor((x(ix1,ix2,ix3,3)-xprobmin3+gzone)/dr+0.5d0)
         res=x(ix1,ix2,ix3,3)-xprobmin3+gzone-(dble(na)-0.5d0)*dr
         w(ix1,ix2,ix3,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+&
            1)-ra(na))
         w(ix1,ix2,ix3,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+&
            1)-pa(na))
      end do
      end do
      end do
    else if(uawsom_adiab/=0) then
      ! isothermal
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*(1.d0/SRadius-&
         1.d0/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)+SRadius)))
    else
      ! zero beta
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)=sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)
    end if

    MAX0wB03(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = max(0.0d0,&
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3)))
    MIN0wB03(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = min(0.0d0,&
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(3)))
    
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_) = 1.d-7 + 1.d-1*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)-(xprobmin3))/1.0d0)*merge(1.d0, 0.d0,&
        MAX0wB03(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) > 0.d0)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_) = 1.d-7 + 1.d-1*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)-(xprobmin3))/1.0d0)*merge(1.d0, 0.d0,&
        MIN0wB03(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) < 0.d0)

    !w(ixO^S,wAminus_) = 1.d-6 + 1.4d0*dexp(-(x(ixO^S,3)-(xprobmin3))/0.1d0)*merge(1.d0, 0.d0, MAX0wB03(ixO^S) > 0.d0)
    !w(ixO^S,wAplus_) = 1.d-6 + 1.4d0*dexp(-(x(ixO^S,3)-(xprobmin3))/0.1d0)*merge(1.d0, 0.d0, MIN0wB03(ixO^S) < 0.d0)

    call uawsom_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir
    double precision, intent(in)       :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)      :: A(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    ! vector potential
    double precision :: Avec1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)

    call bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,xC,Avec1)

    if (idir==3) then
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=Avec1(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,3)
    else if(idir==2) then 
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=Avec1(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,2)
    else
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=Avec1(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)
    end if

  end subroutine initvecpot_usr

  subroutine bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,A,Bbp)

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! vector potential
    double precision, intent(out)   :: A(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! magnetic field
    double precision, optional, intent(out)   :: Bbp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)

    double precision :: tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       f1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),f2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)=0.d0
    f1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)-&
       charge1_x(1))/(sqrt((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)-charge1_x(1))**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)-charge1_x(2))**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)-charge1_x(3))**2)*((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)-charge1_x(2))**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)-charge1_x(3))**2))
    f2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)-&
       charge2_x(1))/(sqrt((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)-charge2_x(1))**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)-charge2_x(2))**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)-charge2_x(3))**2)*((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)-charge2_x(2))**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)-charge2_x(3))**2))
    A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)=charge1*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)-charge1_x(3))*f1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+charge2*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)-charge2_x(3))*f2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)=-charge1*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)-charge1_x(2))*f1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-charge2*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)-charge2_x(2))*f2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if(present(Bbp)) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=-sqrt((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)-charge1_x(1))**2+(x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)-charge1_x(2))**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)-charge1_x(3))**2)**3
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)-charge1_x(1))/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)-charge1_x(2))/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)-charge1_x(3))/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)-charge2_x(1))**2+(x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)-charge2_x(2))**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)-charge2_x(3))**2)**3
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)=Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)-charge2_x(1))/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)=Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)-charge2_x(2))/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)=Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)-charge2_x(3))/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)=q_para*Bbp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)
    end if

  end subroutine bipolar_field

  ! allow user to change inductive electric field, especially for boundary driven applications
  subroutine driven_electric_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fE,s)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt,qdt
    type(state)                        :: s
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3)

    integer :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3

    ! fix Bz at bottom boundary
    if(s%is_physical_boundary(5)) then
      ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,nghostcells,1:2)=0.d0
    end if

  end subroutine driven_electric_field

  subroutine specialbound_usr(qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters
    integer, intent(in) :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       pth(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),Qp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: xlen1,xlen2,xlen3,dxb1,dxb2,dxb3,startpos1,startpos2,&
       startpos3,coeffrho
    integer :: idir,ix1,ix2,ix3,ixIMmin1,ixIMmin2,ixIMmin3,ixIMmax1,ixIMmax2,&
       ixIMmax3,ixOsmin1,ixOsmin2,ixOsmin3,ixOsmax1,ixOsmax2,ixOsmax3,jxOmin1,&
       jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3
    double precision :: A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: Bfr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    double precision :: MIN0wB03G(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), MAX0wB03G(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(uawsom_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       psi_)=0.d0
    select case(iB)
     case(1)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(1))/w(ixOmax1+nghostcells:ixOmax1+1:-1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(2))/w(ixOmax1+nghostcells:ixOmax1+1:-1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(3))=-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(3))/w(ixOmax1+nghostcells:ixOmax1+1:-1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
           ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
           ixOsmin3=ixOmin3-kr(3,idir);
           do ix1=ixOsmax1,ixOsmin1,-1
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
             !       (-block%ws(ix1+2^%1ixOs^S,idir)&
             !   +4.d0*block%ws(ix1+1^%1ixOs^S,idir))
             ! 2nd order one-sided equal gradient extrapolation
             block%ws(ix1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
                idir) = third*( block%ws(ix1+3,ixOsmin2:ixOsmax2,&
                ixOsmin3:ixOsmax3,idir)-5.d0*block%ws(ix1+2,ixOsmin2:ixOsmax2,&
                ixOsmin3:ixOsmax3,idir)+7.d0*block%ws(ix1+1,ixOsmin2:ixOsmax2,&
                ixOsmin3:ixOsmax3,idir))
           end do
         end do
         ixOsmin1=ixOmin1-kr(1,1);ixOsmin2=ixOmin2-kr(1,2)
         ixOsmin3=ixOmin3-kr(1,3);ixOsmax1=ixOmax1-kr(1,1)
         ixOsmax2=ixOmax2-kr(1,2);ixOsmax3=ixOmax3-kr(1,3);
         jxOmin1=ixOmin1+nghostcells*kr(1,1)
         jxOmin2=ixOmin2+nghostcells*kr(1,2)
         jxOmin3=ixOmin3+nghostcells*kr(1,3)
         jxOmax1=ixOmax1+nghostcells*kr(1,1)
         jxOmax2=ixOmax2+nghostcells*kr(1,2)
         jxOmax3=ixOmax3+nghostcells*kr(1,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            1)=zero
         do ix1=ixOsmax1,ixOsmin1,-1
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ix1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,1)=Qp(ix1+1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)*block%dvolume(ix1+1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%surfaceC(ix1,&
              ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,1)
         end do
         call uawsom_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)
       else
         do ix1=ixOmax1,ixOmin1,-1
           !! 2nd order accuacy zero gradient extrapolation
           !w(ix1^%1ixO^S,mag(:))=third* &
           !           (-w(ix1+2^%1ixO^S,mag(:)) &
           !      +4.0d0*w(ix1+1^%1ixO^S,mag(:)))
           ! 3rd order one-sided equal gradient extrapolation
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)) = third*(w(ix1+3,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))-5.d0*w(ix1+2,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))+7.d0*w(ix1+1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)))
         end do
       end if
       if(uawsom_energy) then
         ixIMmin1=ixOmin1;ixIMmin2=ixOmin2;ixIMmin3=ixOmin3;ixIMmax1=ixOmax1
         ixIMmax2=ixOmax2;ixIMmax3=ixOmax3;
         ixIMmin1=ixOmax1+1;ixIMmax1=ixOmax1+nghostcells;
         call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
            ixImax3,ixIMmin1,ixIMmin2,ixIMmin3,ixIMmax1,ixIMmax2,ixIMmax3,pth)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,rho_)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            p_)=pth(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       else if(uawsom_adiab==0) then
         ! zero beta
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(:))**2,dim=ndim+1)
       else
         ! isothermal
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*(1.d0/SRadius-&
            1.d0/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            3)+SRadius)))
       end if
       call uawsom_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case(2)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(1))/w(ixOmin1-1:ixOmin1-nghostcells:-1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(2))/w(ixOmin1-1:ixOmin1-nghostcells:-1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(3))=-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(3))/w(ixOmin1-1:ixOmin1-nghostcells:-1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==1) cycle
           ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
           ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
           ixOsmin3=ixOmin3-kr(3,idir);
           do ix1=ixOsmin1,ixOsmax1
             !block%ws(ix1^%1ixOs^S,idir) = 1.d0/3.d0*&
             !       (-block%ws(ix1-2^%1ixOs^S,idir)&
             !   +4.d0*block%ws(ix1-1^%1ixOs^S,idir))
             ! 2nd order one-sided equal gradient extrapolation
             block%ws(ix1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
                idir) = third*( block%ws(ix1-3,ixOsmin2:ixOsmax2,&
                ixOsmin3:ixOsmax3,idir)-5.d0*block%ws(ix1-2,ixOsmin2:ixOsmax2,&
                ixOsmin3:ixOsmax3,idir)+7.d0*block%ws(ix1-1,ixOsmin2:ixOsmax2,&
                ixOsmin3:ixOsmax3,idir))
           end do
         end do
         ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmin3=ixOmin3;ixOsmax1=ixOmax1
         ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
         jxOmin1=ixOmin1-nghostcells*kr(1,1)
         jxOmin2=ixOmin2-nghostcells*kr(1,2)
         jxOmin3=ixOmin3-nghostcells*kr(1,3)
         jxOmax1=ixOmax1-nghostcells*kr(1,1)
         jxOmax2=ixOmax2-nghostcells*kr(1,2)
         jxOmax3=ixOmax3-nghostcells*kr(1,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            1)=zero
         do ix1=ixOsmin1,ixOsmax1
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ix1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,1)=-Qp(ix1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)*block%dvolume(ix1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%surfaceC(ix1,&
              ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,1)
         end do
         call uawsom_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)
       else
         do ix1=ixOmin1,ixOmax1
           !! 2nd order accuacy zero gradient extrapolation
           !w(ix1^%1ixO^S,mag(:))=third* &
           !           (-w(ix1-2^%1ixO^S,mag(:)) &
           !      +4.0d0*w(ix1-1^%1ixO^S,mag(:)))
           ! 3rd order one-sided equal gradient extrapolation
           w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)) = third*(w(ix1-3,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))-5.d0*w(ix1-2,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))+7.d0*w(ix1-1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)))
         end do
       end if
       if(uawsom_energy) then
         ixIMmin1=ixOmin1;ixIMmin2=ixOmin2;ixIMmin3=ixOmin3;ixIMmax1=ixOmax1
         ixIMmax2=ixOmax2;ixIMmax3=ixOmax3;
         ixIMmin1=ixOmin1-nghostcells;ixIMmax1=ixOmin1-1;
         call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
            ixImax3,ixIMmin1,ixIMmin2,ixIMmin3,ixIMmax1,ixIMmax2,ixIMmax3,pth)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,rho_)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            p_)=pth(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       else if(uawsom_adiab==0) then
         ! zero beta
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(:))**2,dim=ndim+1)
       else
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*(1.d0/SRadius-&
            1.d0/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            3)+SRadius)))
       end if
       call uawsom_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case(3)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,&
          ixOmin3:ixOmax3,mom(1))/w(ixOmin1:ixOmax1,&
          ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,&
          ixOmin3:ixOmax3,mom(2))/w(ixOmin1:ixOmax1,&
          ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,&
          ixOmin3:ixOmax3,mom(3))/w(ixOmin1:ixOmax1,&
          ixOmax2+nghostcells:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
           ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
           ixOsmin3=ixOmin3-kr(3,idir);
           do ix2=ixOsmax2,ixOsmin2,-1
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
             !       (-block%ws(ix2+2^%2ixOs^S,idir)&
             !   +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
             ! 2nd order one-sided equal gradient extrapolation
             block%ws(ixOsmin1:ixOsmax1,ix2,ixOsmin3:ixOsmax3,&
                idir) = third*( block%ws(ixOsmin1:ixOsmax1,ix2+3,&
                ixOsmin3:ixOsmax3,idir)-5.d0*block%ws(ixOsmin1:ixOsmax1,ix2+2,&
                ixOsmin3:ixOsmax3,idir)+7.d0*block%ws(ixOsmin1:ixOsmax1,ix2+1,&
                ixOsmin3:ixOsmax3,idir))
           end do
         end do
         ixOsmin1=ixOmin1-kr(2,1);ixOsmin2=ixOmin2-kr(2,2)
         ixOsmin3=ixOmin3-kr(2,3);ixOsmax1=ixOmax1-kr(2,1)
         ixOsmax2=ixOmax2-kr(2,2);ixOsmax3=ixOmax3-kr(2,3);
         jxOmin1=ixOmin1+nghostcells*kr(2,1)
         jxOmin2=ixOmin2+nghostcells*kr(2,2)
         jxOmin3=ixOmin3+nghostcells*kr(2,3)
         jxOmax1=ixOmax1+nghostcells*kr(2,1)
         jxOmax2=ixOmax2+nghostcells*kr(2,2)
         jxOmax3=ixOmax3+nghostcells*kr(2,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            2)=zero
         do ix2=ixOsmax2,ixOsmin2,-1
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ixOsmin1:ixOsmax1,ix2,ixOsmin3:ixOsmax3,&
              2)=Qp(ixOmin1:ixOmax1,ix2+1,&
              ixOmin3:ixOmax3)*block%dvolume(ixOmin1:ixOmax1,ix2+1,&
              ixOmin3:ixOmax3)/block%surfaceC(ixOsmin1:ixOsmax1,ix2,&
              ixOsmin3:ixOsmax3,2)
         end do
         call uawsom_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)
       else
         do ix2=ixOmax2,ixOmin2,-1
           !! 2nd order accuacy zero gradient extrapolation
           !w(ix2^%2ixO^S,mag(:))=third* &
           !           (-w(ix2+2^%2ixO^S,mag(:)) &
           !      +4.0d0*w(ix2+1^%2ixO^S,mag(:)))
           ! 3rd order one-sided equal gradient extrapolation
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
              mag(:)) = third*(w(ixOmin1:ixOmax1,ix2+3,ixOmin3:ixOmax3,&
              mag(:))-5.d0*w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,&
              mag(:))+7.d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,mag(:)))
         end do
       end if
       if(uawsom_energy) then
         ixIMmin1=ixOmin1;ixIMmin2=ixOmin2;ixIMmin3=ixOmin3;ixIMmax1=ixOmax1
         ixIMmax2=ixOmax2;ixIMmax3=ixOmax3;
         ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+nghostcells;
         call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
            ixImax3,ixIMmin1,ixIMmin2,ixIMmin3,ixIMmax1,ixIMmax2,ixIMmax3,pth)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            p_)=pth(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,&
            ixOmin3:ixOmax3)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,&
            ixOmin3:ixOmax3,rho_)
       else if(uawsom_adiab==0) then
         ! zero beta
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(:))**2,dim=ndim+1)
       else
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*(1.d0/SRadius-&
            1.d0/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            3)+SRadius)))
       end if
       call uawsom_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case(4)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,&
          ixOmin3:ixOmax3,mom(1))/w(ixOmin1:ixOmax1,&
          ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,&
          ixOmin3:ixOmax3,mom(2))/w(ixOmin1:ixOmax1,&
          ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(3))=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,&
          ixOmin3:ixOmax3,mom(3))/w(ixOmin1:ixOmax1,&
          ixOmin2-1:ixOmin2-nghostcells:-1,ixOmin3:ixOmax3,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==2) cycle
           ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
           ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
           ixOsmin3=ixOmin3-kr(3,idir);
           do ix2=ixOsmin2,ixOsmax2
             !block%ws(ix2^%2ixOs^S,idir) = 1.d0/3.d0*&
             !       (-block%ws(ix2-2^%2ixOs^S,idir)&
             !   +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
             ! 2nd order one-sided equal gradient extrapolation
             block%ws(ixOsmin1:ixOsmax1,ix2,ixOsmin3:ixOsmax3,&
                idir) = third*( block%ws(ixOsmin1:ixOsmax1,ix2-3,&
                ixOsmin3:ixOsmax3,idir)-5.d0*block%ws(ixOsmin1:ixOsmax1,ix2-2,&
                ixOsmin3:ixOsmax3,idir)+7.d0*block%ws(ixOsmin1:ixOsmax1,ix2-1,&
                ixOsmin3:ixOsmax3,idir))
           end do
         end do
         ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmin3=ixOmin3;ixOsmax1=ixOmax1
         ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
         jxOmin1=ixOmin1-nghostcells*kr(2,1)
         jxOmin2=ixOmin2-nghostcells*kr(2,2)
         jxOmin3=ixOmin3-nghostcells*kr(2,3)
         jxOmax1=ixOmax1-nghostcells*kr(2,1)
         jxOmax2=ixOmax2-nghostcells*kr(2,2)
         jxOmax3=ixOmax3-nghostcells*kr(2,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            2)=zero
         do ix2=ixOsmin2,ixOsmax2
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ixOsmin1:ixOsmax1,ix2,ixOsmin3:ixOsmax3,&
              2)=-Qp(ixOmin1:ixOmax1,ix2,&
              ixOmin3:ixOmax3)*block%dvolume(ixOmin1:ixOmax1,ix2,&
              ixOmin3:ixOmax3)/block%surfaceC(ixOsmin1:ixOsmax1,ix2,&
              ixOsmin3:ixOsmax3,2)
         end do
         call uawsom_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)
       else
         do ix2=ixOmax2,ixOmin2,-1
           !! 2nd order accuacy zero gradient extrapolation
           !w(ix2^%2ixO^S,mag(:))=third* &
           !           (-w(ix2-2^%2ixO^S,mag(:)) &
           !      +4.0d0*w(ix2-1^%2ixO^S,mag(:)))
           ! 3rd order one-sided equal gradient extrapolation
           w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
              mag(:)) = third*(w(ixOmin1:ixOmax1,ix2-3,ixOmin3:ixOmax3,&
              mag(:))-5.d0*w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,&
              mag(:))+7.d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,mag(:)))
         end do
       end if
       if(uawsom_energy) then
         ixIMmin1=ixOmin1;ixIMmin2=ixOmin2;ixIMmin3=ixOmin3;ixIMmax1=ixOmax1
         ixIMmax2=ixOmax2;ixIMmax3=ixOmax3;
         ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmin2-1;
         call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
            ixImax3,ixIMmin1,ixIMmin2,ixIMmin3,ixIMmax1,ixIMmax2,ixIMmax3,pth)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            p_)=pth(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,&
            ixOmin3:ixOmax3)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,&
            ixOmin3:ixOmax3,rho_)
       else if(uawsom_adiab==0) then
         ! zero beta
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(:))**2,dim=ndim+1)
       else
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*(1.d0/SRadius-&
            1.d0/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            3)+SRadius)))
       end if
       call uawsom_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case(5) 
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))/w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmax3+nghostcells:ixOmax3+1:-1,mom(2))/w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(3))=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmax3+nghostcells:ixOmax3+1:-1,mom(3))/w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmax3+nghostcells:ixOmax3+1:-1,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==3) cycle
             ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
             ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
             ixOsmin3=ixOmin3-kr(3,idir);
             do ix3=ixOsmax3,ixOsmin3,-1
               ! 3rd order one-sided equal gradient extrapolation
               block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ix3,&
                  idir) = 1.d0/11.d0*( -2.d0*block%ws(ixOsmin1:ixOsmax1,&
                  ixOsmin2:ixOsmax2,ix3+4,&
                  idir)+11.d0*block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                  ix3+3,idir)-27.d0*block%ws(ixOsmin1:ixOsmax1,&
                  ixOsmin2:ixOsmax2,ix3+2,&
                  idir)+29.d0*block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                  ix3+1,idir))
             end do
         end do
         ixOsmin1=ixOmin1-kr(3,1);ixOsmin2=ixOmin2-kr(3,2)
         ixOsmin3=ixOmin3-kr(3,3);ixOsmax1=ixOmax1-kr(3,1)
         ixOsmax2=ixOmax2-kr(3,2);ixOsmax3=ixOmax3-kr(3,3);
         jxOmin1=ixOmin1+nghostcells*kr(3,1)
         jxOmin2=ixOmin2+nghostcells*kr(3,2)
         jxOmin3=ixOmin3+nghostcells*kr(3,3)
         jxOmax1=ixOmax1+nghostcells*kr(3,1)
         jxOmax2=ixOmax2+nghostcells*kr(3,2)
         jxOmax3=ixOmax3+nghostcells*kr(3,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            3)=zero
         do ix3=ixOsmax3,ixOsmin3,-1
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ix3,&
              3)=Qp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ix3+1)*block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ix3+1)/block%surfaceC(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ix3,3)
         end do
         call uawsom_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)
       else
         !! 4th order accuacy zero gradient extrapolation
         !do ix3=ixOmax3,ixOmin3,-1
         !  w(ix3^%3ixO^S,mag(:)) = 0.04d0*&
         !     ( -3.d0*w(ix3+4^%3ixO^S,mag(:))&
         !      +16.d0*w(ix3+3^%3ixO^S,mag(:))&
         !      -36.d0*w(ix3+2^%3ixO^S,mag(:))&
         !      +48.d0*w(ix3+1^%3ixO^S,mag(:)))
         !end do
         do ix3=ixOmax3,ixOmin3,-1
           ! 3rd order one-sided equal gradient extrapolation
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
              mag(:)) = 1.d0/11.d0*( -2.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ix3+4,mag(:))+11.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,&
              mag(:))-27.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,&
              mag(:))+29.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(:)))
         end do
       end if
       if(uawsom_energy) then
         do ix3=ixOmin3,ixOmax3
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,rho_)=rbc(ix3)
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,p_)=pbc(ix3)
         enddo
       else if(uawsom_adiab==0) then
         ! zero beta
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(:))**2,dim=ndim+1)
       else
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=rhob*dexp(usr_grav*SRadius**2/Tiso*(1.d0/SRadius-&
            1.d0/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            3)+SRadius)))
       end if

       MAX0wB03G(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = max(0.0d0,&
          w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(3)))
       MIN0wB03G(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = min(0.0d0,&
          w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(3)))

       do ix2=ixOmin2,ixOmax2
        do ix1=ixOmin1,ixOmax1
         w(ix1,ix2,ixOmin3,wkminus_) = 1.0d-1*merge(1.d0, 0.d0, MAX0wB03G(ix1,&
            ix2,ixOmin3) > 0.d0) !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2)) 
         w(ix1,ix2,ixOmax3,wkminus_) = 1.0d-1*merge(1.d0, 0.d0, MAX0wB03G(ix1,&
            ix2,ixOmin3) > 0.d0) !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2)) 

         w(ix1,ix2,ixOmin3,wkplus_) =  1.0d-1*merge(1.d0, 0.d0, MIN0wB03G(ix1,&
            ix2,ixOmin3) < 0.d0) !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))
         w(ix1,ix2,ixOmax3,wkplus_) =  1.0d-1*merge(1.d0, 0.d0, MIN0wB03G(ix1,&
            ix2,ixOmin3) < 0.d0) !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))

         w(ix1,ix2,ixOmin3,wAminus_) = 0.0d-1*merge(1.d0, 0.d0, MAX0wB03G(ix1,&
            ix2,ixOmin3) > 0.d0) !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2)) 
         w(ix1,ix2,ixOmax3,wAminus_) = 0.0d-1*merge(1.d0, 0.d0, MAX0wB03G(ix1,&
            ix2,ixOmin3) > 0.d0) !+ 5.0d-1*(0.0d0+1.0d0*MAX0wB02G(ix1,ixOmax2)) 
        
         w(ix1,ix2,ixOmin3,wAplus_) =  0.0d-1*merge(1.d0, 0.d0, MIN0wB03G(ix1,&
            ix2,ixOmin3) < 0.d0) !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))
         w(ix1,ix2,ixOmax3,wAplus_) =  0.0d-1*merge(1.d0, 0.d0, MIN0wB03G(ix1,&
            ix2,ixOmin3) < 0.d0) !- 1.0d-1*(-0.0d0+1.0d0*MIN0wB02G(ix1,ixOmax2))
        end do
       end do
       
       call bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,A,Bfr)
       !do ix1 = ixOmin1,ixOmax1
       ! do ix2 = ixOmin2,ixOmax2
       !  do ix3 = ixOmin3,ixOmax3
       !    if (Bfr(ix1,ix2,ix3,3) < 0.0d0) then
       !      w(ix1,ix2,ix3,wkplus_) = 1.0d-1
       !      w(ix1,ix2,ix3,wkminus_) = w(ix1,ix2,ixOmax3+1,wkminus_)
       !    elseif (Bfr(ix1,ix2,ix3,3) > 0.0d0) then
       !      w(ix1,ix2,ix3,wkminus_) = 1.0d-1
       !      w(ix1,ix2,ix3,wkplus_) = w(ix1,ix2,ixOmax3+1,wkplus_)
       !    end if
       !  end do
       ! end do
       !end do

       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAminus_) = 0.0d0
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAplus_) = 0.0d0

       !w(ixO^S,wkminus_) = 0.0d0
       !w(ixO^S,wkplus_) = 0.0d0

       call uawsom_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case(6)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3-1:ixOmin3-nghostcells:-1,mom(1))/w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3-1:ixOmin3-nghostcells:-1,mom(2))/w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3-1:ixOmin3-nghostcells:-1,mom(3))/w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-nghostcells:-1,rho_)
       if(stagger_grid) then
         do idir=1,nws
           if(idir==3) cycle
           ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
           ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir)
           ixOsmin3=ixOmin3-kr(3,idir);
           do ix3=ixOsmin3,ixOsmax3
             ! 3rd order one-sided equal gradient extrapolation
             block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ix3,&
                idir) = 1.d0/3.d0*( block%ws(ixOsmin1:ixOsmax1,&
                ixOsmin2:ixOsmax2,ix3-3,idir)-5.d0*block%ws(ixOsmin1:ixOsmax1,&
                ixOsmin2:ixOsmax2,ix3-2,idir)+7.d0*block%ws(ixOsmin1:ixOsmax1,&
                ixOsmin2:ixOsmax2,ix3-1,idir))
             ! 3rd order one-sided zero gradient extrapolation
             !block%ws(ix3^%3ixOs^S,idir) = 1.d0/11.d0*&
             !   (+2.d0*block%ws(ix3-3^%3ixOs^S,idir)&
             !    -9.d0*block%ws(ix3-2^%3ixOs^S,idir)&
             !   +18.d0*block%ws(ix3-1^%3ixOs^S,idir))
           end do
         end do
         ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmin3=ixOmin3;ixOsmax1=ixOmax1
         ixOsmax2=ixOmax2;ixOsmax3=ixOmax3;
         jxOmin1=ixOmin1-nghostcells*kr(3,1)
         jxOmin2=ixOmin2-nghostcells*kr(3,2)
         jxOmin3=ixOmin3-nghostcells*kr(3,3)
         jxOmax1=ixOmax1-nghostcells*kr(3,1)
         jxOmax2=ixOmax2-nghostcells*kr(3,2)
         jxOmax3=ixOmax3-nghostcells*kr(3,3);
         block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ixOsmin3:ixOsmax3,&
            3)=zero
         do ix3=ixOsmin3,ixOsmax3
           call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Qp)
           block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ix3,&
              3)=-Qp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ix3)*block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ix3)/block%surfaceC(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,ix3,3)
         end do
         call uawsom_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
            ixOmax3,block)
       else
         do ix3=ixOmin3,ixOmax3
           ! 3rd order one-sided equal gradient extrapolation
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
              mag(:)) = 1.d0/3.d0*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-3,&
              mag(:))-5.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,&
              mag(:))+7.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,mag(:)))
           !! 2nd order one-sided zero gradient extrapolation
           !w(ix3^%3ixO^S,mag(:)) = 1.d0/3.d0*&
           !       (-w(ix3-2^%3ixO^S,mag(:))&
           !   +4.d0*w(ix3-1^%3ixO^S,mag(:)))
         end do
       end if
       if(uawsom_energy) then
         ixIMmin1=ixOmin1;ixIMmin2=ixOmin2;ixIMmin3=ixOmin3;ixIMmax1=ixOmax1
         ixIMmax2=ixOmax2;ixIMmax3=ixOmax3;
         ixIMmin3=ixOmin3-1;ixIMmax3=ixOmax3;
         call getggrav(tmp1,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixIMmin1,ixIMmin2,ixIMmin3,ixIMmax1,ixIMmax2,ixIMmax3,x)
         ixIMmin3=ixOmin3-1;ixIMmax3=ixOmin3-1;
         call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
            ixImax3,ixIMmin1,ixIMmin2,ixIMmin3,ixIMmax1,ixIMmax2,ixIMmax3,pth)
         pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1)=pth(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3-1)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3-1,rho_)
         where(pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1)<1.618d0)
           pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1)=1.618d0
         end where
         tmp2=0.d0
         do ix3=ixOmin3,ixOmax3
           tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3)=tmp2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3)+0.5d0*(tmp1(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ix3-1)+tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ix3))/pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1)
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,rho_)=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3-1,rho_)*dexp(tmp2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3)*dxlevel(3))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,p_)=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ix3,rho_)*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3-1)
         enddo
       else if(uawsom_adiab==0) then
         ! zero beta
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)=sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(:))**2,dim=ndim+1)
       else
         coeffrho=usr_grav*SRadius**2/Tiso
         do ix3=ixOmin3,ixOmax3
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,rho_)=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3-1,rho_)*dexp(coeffrho*(1.d0/(SRadius+&
              x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1,&
              3))-1.d0/(SRadius+x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,3))))
         enddo
       end if
       call uawsom_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine getggrav(ggrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=usr_grav*(SRadius/(SRadius+x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)))**2
  end subroutine

  !==============================================================================
  ! Purpose: get gravity field
  !==============================================================================
  subroutine gravity(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    double precision                :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    gravity_field=0.d0
    call getggrav(ggrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x)
    gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)=ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
  end subroutine gravity

  subroutine specialsource(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,wCT,qt,w,&
     x)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim), wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: lQgrid(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),bQgrid(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    !! add global background heating bQ
    call getbQ(bQgrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qtC,wCT,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)+qdt*bQgrid(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    !! add localized heating lQ
    !  call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
    !  w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
    !endif

  end subroutine specialsource

  subroutine getbQ(bQgrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    !!calculate background heating bQ, Mok 2016 ApJ
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(out) :: bQgrid(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: Bmag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir),&
       cvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: curlo,curhi
    integer :: idims,idir

    bQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=bQ0*dexp(-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)/6.d0)
    !if(B0field) then
    !  Bmag(ixI^S)=dsqrt(sum((w(ixI^S,mag(:))+block%B0(ixI^S,:,0))**2,dim=ndim+1))
    !  do idir=1,ndir
    !    bvec(ixI^S,idir)=(w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0))/Bmag(ixI^S)
    !  end do
    !else
    !  Bmag(ixI^S)=dsqrt(sum(w(ixI^S,mag(:))**2,dim=ndim+1))
    !  do idir=1,ndir
    !    bvec(ixI^S,idir)=w(ixI^S,mag(idir))/Bmag(ixI^S)
    !  end do
    !endif
    !cvec=0.d0
    !! calculate local curvature of magnetic field
    !do idims=1,ndim
    !  call gradient(bvec(ixI^S,1),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,1)=cvec(ixO^S,1)+bvec(ixO^S,idims)*tmp(ixO^S)
    !  call gradient(bvec(ixI^S,2),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,2)=cvec(ixO^S,2)+bvec(ixO^S,idims)*tmp(ixO^S)
    !  call gradient(bvec(ixI^S,3),ixI^L,ixO^L,idims,tmp) 
    !  cvec(ixO^S,3)=cvec(ixO^S,3)+bvec(ixO^S,idims)*tmp(ixO^S)
    !end do 
    !tmp(ixO^S)=dsqrt(sum(cvec(ixO^S,:)**2,dim=ndim+1))
    !! set lower and upper limit for curvature
    !curlo=0.1d0 ! 1/10 
    !curhi=2.d0 ! 1/0.5
    !where(tmp(ixO^S)<curlo)
    !  tmp(ixO^S)=curlo
    !elsewhere(tmp(ixO^S)>curhi)
    !  tmp(ixO^S)=curhi
    !end where
    !
    !bQgrid(ixO^S)=bQ0*Bmag(ixO^S)**1.75d0*w(ixO^S,rho_)**0.125d0*tmp(ixO^S)**0.75d0
  
  end subroutine getbQ

  !==============================================================================
  ! Purpose: Enforce additional refinement or coarsening. One can use the
  !          coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !==============================================================================
  subroutine special_refine_grid(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x,&
     refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qt, w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (block%is_physical_boundary(5)) then
      refine=1
      coarsen=-1
    endif
  end subroutine special_refine_grid

  !==============================================================================
  ! Purpose: 
  !   this subroutine can be used in convert, to add auxiliary variables to the
  !   converted output file, for further analysis using tecplot, paraview, ....
  !   these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  !   the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  !   corresponding normalization values (default value 1)
  !==============================================================================
  subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),dip(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       B2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir) :: Btotal,qvec,curlvec
    integer                            :: ix1,ix2,ix3,idirmin,idims,idir,jdir,&
       kdir

    ! Btotal & B^2
    if(B0field) then
      Btotal(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir,b0i)
    else
      Btotal(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndir))
    end if
    B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sum((Btotal(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+1)=dsqrt(B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_))
    ! output divB1
    call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+2)=divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    ! output the plasma beta p*2/B**2
    call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
    where(B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/=0.d0)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         nw+3)=2.d0*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    else where
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+3)=0.d0
    end where
    ! store current
    call curlvector(Btotal,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         nw+3+idir)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)
    end do
    ! calculate Lorentz force
    qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)=zero
    do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
      if(lvc(idir,jdir,kdir)/=0)then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mag(kdir))
        if(lvc(idir,jdir,kdir)==1)then
          qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        else
          qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        endif
      endif
    enddo; enddo; enddo
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         nw+3+ndir+idir)=qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)
    end do
    ! find magnetic dips
    !dip=0.d0
    !do idir=1,ndir
    !  call gradient(w(ixI^S,mag(3)),ixI^L,ixO^L,idir,tmp)
    !  dip(ixO^S)=dip(ixO^S)+w(ixO^S,b0_+idir)*tmp(ixO^S)
    !end do
    !where(dabs(w(ixO^S,mag(3)))<0.08d0 .and. dip(ixO^S)>=0.d0)
    !  w(ixO^S,nw+8)=1.d0
    !elsewhere
    !  w(ixO^S,nw+8)=0.d0
    !end where
  end subroutine specialvar_output

  subroutine specialset_B0(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)

    double precision :: A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)

    call bipolar_field(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,A,wB0)

  end subroutine specialset_B0

  subroutine specialset_J0(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,wJ0)
  ! Here add a time-independent background current density 
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: wJ0(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:ndir)

    wJ0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)=0.d0

  end subroutine specialset_J0

  !==============================================================================
  ! Purpose: names for special variable output
  !==============================================================================
  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Alfv divB beta j1 j2 j3 L1 L2 L3'
  end subroutine specialvarnames_output

end module mod_usr
