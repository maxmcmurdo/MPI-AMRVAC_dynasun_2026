module mod_usr
  use mod_uawsom
  implicit none
  double precision :: Bbnd,rhobnd,Tbnd,gamma0,gamma2,gamma4,usr_grav

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    uawsom_gamma=5.0d0/3.0d0
    uawsom_eta=zero

    call set_coordinate_system("spherical_2.5D")
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source => specialsource
    usr_refine_grid     => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_set_B0          => specialset_B0
    usr_add_aux_names   => specialvarnames_output

    call uawsom_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

    ! we want to have 1.1 Gauss at the Solar surface at the equator
    Bbnd = 2.40858d0 !2.40858d0
    ! we want a number density of n=1.0E8 cm-3 at the inner boundary
    rhobnd = 2.0d0    !/14.0d0*3.0d0
    ! we want a temperature of 1.5E6 K at the inner boundary
    Tbnd = 2.493976d0      ! before it was 2.493976d0
    ! solar differential rotation:
    gamma0 = 0.0194796d0/2.0d0/dpi
    gamma2 = -0.0024349d0/2.0d0/dpi
    gamma4 = -0.0031306d0/2.0d0/dpi

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: lat(ixI^S), Br_arc, Bth_arc, deltath, alpha, B0, arcade, shift_arc, shift_dip

    double precision :: MAX0wB02(ixI^S), MIN0wB02(ixI^S) 
    logical, save :: first=.true.
    integer :: ix1,ix2
  
    if (first .and. mype==0) then
       write(*,*)'Running solar wind'
       first=.false.
    endif

    arcade=0

    ! n = 1e8 cm-3 at the inner boundary and decreases like r^-2
    w(ixO^S,rho_)=rhobnd/x(ixO^S,1)**2
    ! vr is 10 km/s everywhere
    w(ixO^S,mom(1))  = 2.0d0*w(ixO^S,rho_)
    w(ixO^S,mom(2))  = 0.0d0
    ! at the inner boundary: differential rotation of the Sun
    ! vphi decreases like r^-2
    w(ixO^S,mom(3))  = 0.0d0
 !   w(ixO^S,rho_)*x(ixO^S,1)*dsin(x(ixO^S,2))* &
 !       (gamma0+gamma2*dcos(x(ixO^S,2))**2+gamma4*dcos(x(ixO^S,2))**4)/x(ixO^S,1)**2


    shift_dip = 0.0d0
    
    do ix1=ixOmin1,ixOmax1
       do ix2=ixOmin2,ixOmax2  
             w(ix1,ix2,mag(1))  = Bbnd*2.0d0*dcos(x(ix1,ix2,2)-shift_dip)/x(ix1,ix2,1)**3    
             w(ix1,ix2,mag(2))  = Bbnd*dsin(x(ix1,ix2,2)-shift_dip)/x(ix1,ix2,1)**3
        enddo
    enddo

    w(ixO^S,mag(3))  = 0.0d0

    do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        MAX0wB02(ix1,ix2) = max(0.0d0,w(ix1,ix2,mag(1)))
        MIN0wB02(ix1,ix2) = min(0.0d0,w(ix1,ix2,mag(1)))
      end do
    end do

    w(ixO^S,wkminus_) = 1.d-8 + 1.67d-5*dexp(-(x(ixO^S,1)-(xprobmin1))/0.1d0)*merge(1.d0, 0.d0, MAX0wB02(ixO^S) > 0.d0) !> small background everywhere to avoid steep gradients
    w(ixO^S,wkplus_)  = 1.d-8 + 1.67d-5*dexp(-(x(ixO^S,1)-(xprobmin1))/0.1d0)*merge(1.d0, 0.d0, MIN0wB02(ixO^S) < 0.d0)

    w(ixO^S,wAminus_) = 0.0d0 !1.d-6 + 1.67d-2*dexp(-(x(ixO^S,1)-(xprobmin1)))  !*merge(1.d0, 0.d0, MAX0wB02(ixO^S) > 0.d0)
    w(ixO^S,wAplus_)  = 0.0d0 !1.d-6  + 1.67d-2*dexp(-(x(ixO^S,1)-(xprobmin1))) !*merge(1.d0, 0.d0, MIN0wB02(ixO^S) < 0.d0)
    
    !pressure decreases like 1/r^2 (T is thus fixed)
    w(ixO^S,e_)   = Tbnd*rhobnd/x(ixO^S,1)**2/(uawsom_gamma-one) &
                  +0.5d0*(w(ixO^S,mom(1))**2+w(ixO^S,mom(2))**2+w(ixO^S,mom(3))**2)/w(ixO^S,rho_) &
                  +0.5d0*(w(ixO^S,mag(1))**2+w(ixO^S,mag(2))**2+w(ixO^S,mag(3))**2) &
                  +w(ixO^S,wkplus_) + w(ixO^S,wkminus_) &
                  +w(ixO^S,wAplus_) + w(ixO^S,wAminus_) !> Added wave energy to total energy calculation

    if(uawsom_glm) w(ixO^S,psi_)=0.0d0

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S)
    double precision :: omega, LM, Pg1, Pg2, Ti, Pin, latitude

    double precision :: th, lat, Br_arc, deltath, alpha, B0, arcade, shift_arc, Bth_arc, shift_dip
    double precision :: shear1, v0_max1, deltat1, colat_shear1, shift_shear1, v_shear1, delth_shear1, time_stop1  
    double precision :: shear2, v0_max2, deltat2, colat_shear2, shift_shear2, v_shear2, delth_shear2, time_stop2    

    integer :: ix1,ix2,ixA^L

    select case (iB)
    case (1)

       !!!!!!!!!!!!!!!!!!
       ! inner boundary !
       !!!!!!!!!!!!!!!!!!

       ! Loop through ghost cells along boundary (in the theta direction)
       do ix2=ixOmin2,ixOmax2

          arcade = 0
          shear1  = 0
          shear2  = 0

          ! The density is fixed at the inner boundary
          w(ixOmax1,ix2,rho_) = 2.0d0*rhobnd-w(ixOmax1+1,ix2,rho_)
          w(ixOmin1,ix2,rho_) = 4.0d0*rhobnd-3.0d0*w(ixOmax1+1,ix2,rho_)

   !       omega=gamma0+gamma2*(dcos(x(ixOmax1,ix2,2)))**2+gamma4*(dcos(x(ixOmax1,ix2,2)))**4
          
   !       w(ixOmax1,ix2,mom(3)) = x(ixOmax1,ix2,1)*dsin(x(ixOmax1,ix2,2))*w(ixOmax1,ix2,rho_)&
   !                *(2.0d0*(omega+(v_shear1+v_shear2)/dsin(x(ixOmax1,ix2,2)))-w(ixOmax1+1,ix2,mom(3))/(x(ixOmax1+1,ix2,1)*dsin(x(ixOmax1+1,ix2,2))*w(ixOmax1+1,ix2,rho_)))
   !       w(ixOmin1,ix2,mom(3)) = x(ixOmin1,ix2,1)*dsin(x(ixOmin1,ix2,2))*w(ixOmin1,ix2,rho_)&
   !               *(4.0d0*(omega+(v_shear1+v_shear2)/dsin(x(ixOmax1,ix2,2)))-3.0d0*w(ixOmax1+1,ix2,mom(3))/(x(ixOmax1+1,ix2,1)*dsin(x(ixOmax1+1,ix2,2))*w(ixOmax1+1,ix2,rho_)))

          shift_dip = 0.0d0

          !!! this is where I defined the extra magnetic field for arcades

          !!Br_arc  = 0

          w(ixOmax1,ix2,mag(2)) = w(ixOmax1+1,ix2,mag(2))*(x(ixOmax1+1,ix2,1)**3)/(x(ixOmax1,ix2,1)**3)
          w(ixOmin1,ix2,mag(2)) = w(ixOmax1+1,ix2,mag(2))*(x(ixOmax1+1,ix2,1)**3)/(x(ixOmin1,ix2,1)**3)

          ! r^2*Br is fixed at the inner boundary (this is how we force a dipole field for the Sun)
          w(ixOmax1,ix2,mag(1)) = (2.0d0*(2.0d0*Bbnd*dcos(x(ixOmax1,ix2,2)-shift_dip))-w(ixOmax1+1,ix2,mag(1))*x(ixOmax1+1,ix2,1)**2)/x(ixOmax1,ix2,1)**2
          w(ixOmin1,ix2,mag(1)) = (4.0d0*(2.0d0*Bbnd*dcos(x(ixOmin1,ix2,2)-shift_dip))-3.0d0*w(ixOmax1+1,ix2,mag(1))*x(ixOmax1+1,ix2,1)**2)/x(ixOmin1,ix2,1)**2

          !> MAX: Set injection of wave energy that propagates into the domain: merge ensures B0 dependence respected          
          w(ixOmin1,ix2,wkminus_) = 1.67d-5*merge(1.d0, 0.d0, w(ixOmax1+1,ix2,mag(1)) > 0.d0) 
          w(ixOmax1,ix2,wkminus_) = 1.67d-5*merge(1.d0, 0.d0, w(ixOmax1+1,ix2,mag(1)) > 0.d0)
          w(ixOmin1,ix2,wkplus_)  = 1.67d-5*merge(1.d0, 0.d0, w(ixOmax1+1,ix2,mag(1)) < 0.d0)
          w(ixOmax1,ix2,wkplus_)  = 1.67d-5*merge(1.d0, 0.d0, w(ixOmax1+1,ix2,mag(1)) < 0.d0)

          w(ixOmin1,ix2,wAminus_) = 0.0d0 !1.67d-2*merge(1.d0, 0.d0, w(ixOmax1+1,ix2,mag(1)) > 0.d0)
          w(ixOmax1,ix2,wAminus_) = 0.0d0 !1.67d-2*merge(1.d0, 0.d0, w(ixOmax1+1,ix2,mag(1)) > 0.d0)
          w(ixOmin1,ix2,wAplus_)  = 0.0d0 !1.67d-2*merge(1.d0, 0.d0, w(ixOmax1+1,ix2,mag(1)) < 0.d0)
          w(ixOmax1,ix2,wAplus_)  = 0.0d0 !1.67d-2*merge(1.d0, 0.d0, w(ixOmax1+1,ix2,mag(1)) < 0.d0)

          latitude = 0.5d0*dpi-x(ixOmax1,ix2,2)

          w(ixOmax1,ix2,mom(1)) = w(ixOmax1+1,ix2,mom(1))*(x(ixOmax1+1,ix2,1)/x(ixOmax1,ix2,1))**2
          w(ixOmin1,ix2,mom(1)) = w(ixOmax1,ix2,mom(1))*(x(ixOmax1,ix2,1)/x(ixOmin1,ix2,1))**2
          
          w(ixOmax1,ix2,mom(2)) = -w(ixOmax1+1,ix2,mom(2))
          w(ixOmin1,ix2,mom(2)) = -3*w(ixOmax1+1,ix2,mom(2))
         
          w(ixOmax1,ix2,mag(3)) = w(ixOmax1+1,ix2,mag(3))
          w(ixOmin1,ix2,mag(3)) = w(ixOmax1+1,ix2,mag(3))

          ! T is fixed at the boundary
          Pin = (uawsom_gamma-one)*(w(ixOmax1+1,ix2,e_)&
                   -0.5d0*(w(ixOmax1+1,ix2,mag(1))**2+w(ixOmax1+1,ix2,mag(2))**2+w(ixOmax1+1,ix2,mag(3))**2) &
                   -0.5d0*(w(ixOmax1+1,ix2,mom(1))**2+w(ixOmax1+1,ix2,mom(2))**2+w(ixOmax1+1,ix2,mom(3))**2)/w(ixOmax1+1,ix2,rho_) &
                   -w(ixOmax1+1,ix2,wkplus_) - w(ixOmax1+1,ix2,wkminus_) &
                   -w(ixOmax1+1,ix2,wAplus_) - w(ixOmax1+1,ix2,wAminus_)) 
          Pg1 = w(ixOmax1,ix2,rho_)*(2.0d0*Tbnd-Pin/w(ixOmax1+1,ix2,rho_))
          Pg2 = w(ixOmin1,ix2,rho_)*(4.0d0*Tbnd-3.0d0*Pin/w(ixOmax1+1,ix2,rho_))
          w(ixOmax1,ix2,e_)  = Pg1/(uawsom_gamma-one) &
                             + 0.5d0*(w(ixOmax1,ix2,mag(1))**2+w(ixOmax1,ix2,mag(2))**2+w(ixOmax1,ix2,mag(3))**2) &
                             + 0.5d0*(w(ixOmax1,ix2,mom(1))**2+w(ixOmax1,ix2,mom(2))**2+w(ixOmax1,ix2,mom(3))**2)/w(ixOmax1,ix2,rho_) &
                             + w(ixOmax1,ix2,wkplus_) + w(ixOmax1,ix2,wkminus_) &
                             + w(ixOmax1,ix2,wAplus_) + w(ixOmax1,ix2,wAminus_)
          w(ixOmin1,ix2,e_)  = Pg2/(uawsom_gamma-one) &
                             + 0.5d0*(w(ixOmin1,ix2,mag(1))**2+w(ixOmin1,ix2,mag(2))**2+w(ixOmin1,ix2,mag(3))**2) &
                             + 0.5d0*(w(ixOmin1,ix2,mom(1))**2+w(ixOmin1,ix2,mom(2))**2+w(ixOmin1,ix2,mom(3))**2)/w(ixOmin1,ix2,rho_) &
                             + w(ixOmin1,ix2,wkplus_) + w(ixOmin1,ix2,wkminus_) &
                             + w(ixOmin1,ix2,wAplus_) + w(ixOmin1,ix2,wAminus_)
       enddo


    case (2) ! Upper radial boundary

      !!!!!!!!!!!!!!!!!!!!
      ! Outflow boundary !
      !!!!!!!!!!!!!!!!!!!!

      do ix1=ixOmin1,ixOmax1
        do ix2=ixOmin2,ixOmax2

          ! r^2*rho is continuous
          w(ix1,ix2,rho_)  = w(ixOmin1-1,ix2,rho_)*(x(ixOmin1-1,ix2,1)/x(ix1,ix2,1))**2

          ! r*r*rho*vr is continuous
          w(ix1,ix2,mom(1)) = w(ixOmin1-1,ix2,mom(1))*(x(ixOmin1-1,ix2,1)/x(ix1,ix2,1))**2

          ! rho*vtheta is continuous
          w(ix1,ix2,mom(2)) = w(ixOmin1-1,ix2,mom(2))

          ! r*vphi is continuous
          w(ix1,ix2,mom(3)) = w(ixOmin1-1,ix2,mom(3))/w(ixOmin1-1,ix2,rho_)*x(ixOmin1-1,ix2,1)/x(ix1,ix2,1)*w(ix1,ix2,rho_)

          ! r*r*Br is continuous
          w(ix1,ix2,mag(1)) = w(ixOmin1-1,ix2,mag(1))*(x(ixOmin1-1,ix2,1)/x(ix1,ix2,1))**2

          ! Btheta is continuous
          w(ix1,ix2,mag(2)) = w(ixOmin1-1,ix2,mag(2))

          ! r*Bphi is continuous
          w(ix1,ix2,mag(3)) = w(ixOmin1-1,ix2,mag(3))*x(ixOmin1-1,ix2,1)/x(ix1,ix2,1)

          ! T is continuous
          ! Ti is the temperature in the last inner cell

          !> MAX: Set wave energy to be either zero or outflow depending on B0
          w(ix1,ix2,wkminus_) = merge(w(ixOmin1-1,ix2,wkminus_), 0.d0, w(ixOmin1-1,ix2,mag(1)) > 0.d0)
          w(ix1,ix2,wkplus_)  = merge(w(ixOmin1-1,ix2,wkplus_), 0.d0, w(ixOmin1-1,ix2,mag(1)) < 0.d0)

          w(ix1,ix2,wAminus_) = 0.0d0 !merge(w(ixOmin1-1,ix2,wAplus_), 0.d0, w(ixOmax1+1,ix2,mag(1)) > 0.d0)
          w(ix1,ix2,wAplus_)  = 0.0d0 !merge(w(ixOmin1-1,ix2,wAplus_), 0.d0, w(ixOmax1+1,ix2,mag(1)) < 0.d0)

          Ti=(uawsom_gamma-one)/w(ixOmin1-1,ix2,rho_)*(w(ixOmin1-1,ix2,e_)&
                  - 0.5d0*(w(ixOmin1-1,ix2,mag(1))**2+w(ixOmin1-1,ix2,mag(2))**2+w(ixOmin1-1,ix2,mag(3))**2) &
                  - 0.5d0*(w(ixOmin1-1,ix2,mom(1))**2+w(ixOmin1-1,ix2,mom(2))**2+w(ixOmin1-1,ix2,mom(3))**2)/w(ixOmin1-1,ix2,rho_) &
                  - w(ixOmin1-1,ix2,wkplus_) - w(ixOmin1-1,ix2,wkminus_) &
                  - w(ixOmin1-1,ix2,wAplus_) - w(ixOmin1-1,ix2,wAminus_))
          w(ix1,ix2,e_) = Ti*w(ix1,ix2,rho_)/(uawsom_gamma-one) &
                            + 0.5d0*(w(ix1,ix2,mag(1))**2+w(ix1,ix2,mag(2))**2+w(ix1,ix2,mag(3))**2) &
                            + 0.5d0*(w(ix1,ix2,mom(1))**2+w(ix1,ix2,mom(2))**2+w(ix1,ix2,mom(3))**2)/w(ix1,ix2,rho_) &
                            + w(ix1,ix2,wkplus_) + w(ix1,ix2,wkminus_) &
                            + w(ix1,ix2,wAplus_) + w(ix1,ix2,wAminus_)

        enddo
      enddo

    end select


  end subroutine specialbound_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: current(ixI^S,7-2*ndir:3)
    integer          :: idirmin,idir,jdir,kdir
    double precision :: lr,lthin,lthout,Rin,Rout,c(ixI^S)
    integer :: iw,ix1,ix2
    double precision :: Bt1,Bt2,Bt3,Bt4
    double precision :: S1,S2,S3,S4

    w(ixO^S,nw+1)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    w(ixO^S,nw+2)=w(ixO^S,mom(2))/w(ixO^S,rho_)
    w(ixO^S,nw+3)=w(ixO^S,mom(3))/w(ixO^S,rho_)
    w(ixO^S,nw+4)=w(ixO^S,mag(1))
    w(ixO^S,nw+5)=w(ixO^S,mag(2))
    w(ixO^S,nw+6)=w(ixO^S,mag(3))
    w(ixO^S,nw+7)=x(ixO^S,1)
    w(ixO^S,nw+8)=x(ixO^S,2)

    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+9)= current(ixO^S,1)
    w(ixO^S,nw+10)= current(ixO^S,2)
    w(ixO^S,nw+11)= current(ixO^S,3)
    w(ixO^S,nw+12)=dsqrt(current(ixO^S,1)**2+current(ixO^S,2)**2+current(ixO^S,3)**2)

    lr = x(ixOmax1,ixOmax2,1) - x(ixOmax1-1,ixOmax2,1)

    do ix1=ixOmin1,ixOmax1
      do ix2=ixOmin2,ixOmax2
        Rin = x(ix1,ixOmin2,1) - lr/2.0d0
        Rout = x(ix1,ixOmin2,1) + lr/2.0d0
        lthin = Rin*(x(ix1,ixOmax2,2)-x(ix1,ixOmax2-1,2))
        lthout = Rout*(x(ix1,ixOmax2,2)-x(ix1,ixOmax2-1,2))

        Bt1 = dsqrt((w(ix1,ix2,mag(1))+w(ix1,ix2+1,mag(1)))**2/4.0d0 + (w(ix1,ix2,mag(3))+w(ix1,ix2+1,mag(3)))**2/4.0d0)
        Bt2 = dsqrt((w(ix1,ix2,mag(2))+w(ix1+1,ix2,mag(2)))**2/4.0d0 + (w(ix1,ix2,mag(3))+w(ix1+1,ix2,mag(3)))**2/4.0d0)
        Bt3 = dsqrt((w(ix1,ix2,mag(1))+w(ix1,ix2-1,mag(1)))**2/4.0d0 + (w(ix1,ix2,mag(3))+w(ix1,ix2-1,mag(3)))**2/4.0d0)
        Bt4 = dsqrt((w(ix1,ix2,mag(2))+w(ix1-1,ix2,mag(2)))**2/4.0d0 + (w(ix1,ix2,mag(3))+w(ix1-1,ix2,mag(3)))**2/4.0d0)

        S1 = sign(Bt1*lr,w(ix1,ix2,mag(1)))
        S2 = sign(Bt2*lthout,w(ix1,ix2,mag(2)))
        S3 = sign(Bt3*lr,-w(ix1,ix2,mag(1)))
        S4 = sign(Bt4*lthin,-w(ix1,ix2,mag(2)))

        c(ix1,ix2)=abs((S1+S2+S3+S4))/(Bt1*lr + Bt2*lthout + Bt3*lr + Bt4*lthin)

        w(ixO^S,nw+13) = c(ix1,ix2)

      enddo
    enddo

  end subroutine specialvar_output
  
  
  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
  
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)
    integer :: ix1,ix2
    double precision :: Br_arc, Bth_arc, B0, alpha, deltath, lat(ixI^S), arcade, shift_arc

    wB0(ixO^S,1)= wB0(ixO^S,1)
    wB0(ixO^S,2)= wB0(ixO^S,2)
    wB0(ixO^S,3)= 0.d0

  end subroutine specialset_B0
  

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames
    varnames='vr vth vphi br bth bphi r th j1 j2 j3 jtot c'
  end subroutine specialvarnames_output

  subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: sin2_thc,thc,hscale,TNow,TTarget,q0,TEquatorial,TPolar
    double precision :: Pth,r,th,theta1,theta2,theta3,theta4
    integer :: ix1, ix2

    ! Gravity
    w(ixO^S,e_)  = w(ixO^S,e_)  - qdt*wCT(ixO^S,mom(1))*19.08/((x(ixO^S,1))**2)
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) - qdt*wCT(ixO^S,rho_)*19.08/((x(ixO^S,1))**2)

    ! Heating/cooling

    q0 = 83.69d0
    theta1 = dpi*10.5d0/180.0d0
    theta2 = dpi*75.5d0/180.0d0
    ! 1.5E6 K
    TEquatorial=2.494d0    ! inainte era 2.494d0
    ! 2.625E6
    TPolar= 4.3645d0  ! inainte era 4.3645d0


    {do ix^DB=ixOmin^DB,ixOmax^DB\}
    !if(x(ix^D,1) < 30.0d0)then
       r  = x(ix^D,1)
       th = x(ix^D,2)
       Pth = (uawsom_gamma-one)*(wCT(ix^D,e_)-0.5d0*(sum(wCT(ix^D,mom(:))**2)/wCT(ix^D,rho_)+sum(wCT(ix^D,mag(:))**2))-wCT(ix^D,wkminus_)-wCT(ix^D,wkplus_)-wCT(ix^D,wAminus_)-wCT(ix^D,wAplus_))
     !  sin2_thc = (dsin(theta1)**2)+(dcos(theta1)**2)*(r-1.0d0)/8.0d0
     !  if ( (r >= 7.0d0) .and. (r < 47.0d0) ) then
     !     sin2_thc = (dsin(theta2)**2)+(dcos(theta2)**2)*(r-7.0d0)/40.0d0
     !  else if (r>= 47.0d0) then
     !     sin2_thc = 1.0d0
     !  end if
  !     thc = dasin(dsqrt(sin2_thc))
       hscale = 4.5d0
       TTarget = TEquatorial
  !     if ( (th < thc) .or. (th > (dpi-thc) ) ) then
  !        hscale=4.5d0*(2.0d0-(sin(th)**2)/sin2_thc)
  !        TTarget = TPolar
  !     end if
       TNow = Pth/wCT(ix^D,rho_)
       w(ix^D,e_)  = w(ix^D,e_)  + qdt*q0*wCT(ix^D,rho_)*(TTarget-TNow)*exp(-((r-1.0d0)/hscale)**2)
   ! endif
    {end do\}

  end subroutine specialsource

   subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
  
    double precision :: rmin1,rmin2,rmax,tend
    double precision :: lr,lthin,lthout,Rin,Rout,c(ixG^T)
    integer :: ix1,ix2
    double precision :: Bt1,Bt2,Bt3,Bt4
    double precision :: S1,S2,S3,S4

    lr = x(ixmax1,ixmax2,1) - x(ixmax1-1,ixmax2,1)

    do ix1=ixmin1,ixmax1
      do ix2=ixmin2,ixmax2
        Rin = x(ix1,ixmin2,1) - lr/2.0d0
        Rout = x(ix1,ixmin2,1) + lr/2.0d0
        lthin = Rin*(x(ix1,ixmax2,2)-x(ix1,ixmax2-1,2))
        lthout = Rout*(x(ix1,ixmax2,2)-x(ix1,ixmax2-1,2))

        Bt1 = dsqrt((w(ix1,ix2,mag(1))+w(ix1,ix2+1,mag(1)))**2/4.0d0 + (w(ix1,ix2,mag(3))+w(ix1,ix2+1,mag(3)))**2/4.0d0)
        Bt2 = dsqrt((w(ix1,ix2,mag(2))+w(ix1+1,ix2,mag(2)))**2/4.0d0 + (w(ix1,ix2,mag(3))+w(ix1+1,ix2,mag(3)))**2/4.0d0)
        Bt3 = dsqrt((w(ix1,ix2,mag(1))+w(ix1,ix2-1,mag(1)))**2/4.0d0 + (w(ix1,ix2,mag(3))+w(ix1,ix2-1,mag(3)))**2/4.0d0)
        Bt4 = dsqrt((w(ix1,ix2,mag(2))+w(ix1-1,ix2,mag(2)))**2/4.0d0 + (w(ix1,ix2,mag(3))+w(ix1-1,ix2,mag(3)))**2/4.0d0)

        S1 = sign(Bt1*lr,w(ix1,ix2,mag(1)))
        S2 = sign(Bt2*lthout,w(ix1,ix2,mag(2)))
        S3 = sign(Bt3*lr,-w(ix1,ix2,mag(1)))
        S4 = sign(Bt4*lthin,-w(ix1,ix2,mag(2)))

        c(ix1,ix2)=abs((S1+S2+S3+S4))/(Bt1*lr + Bt2*lthout + Bt3*lr + Bt4*lthin)

      enddo
    enddo

    if ( any(c(ix^S) > 0.02d0) ) then
            refine = 1
            coarsen = -1
    else if ( all(c(ix^S) < 0.01d0) ) then
            coarsen = 1
            refine = 0
    else
            coarsen = 0
            refine = 0
    end if

      !!!  fixed refinement near the sun 
    if ( any (((x(ix^S,1) < 2.5d0)) .and. ((abs(x(ix^S,2)-0.5d0*dpi) < 1.1d0))) ) then        
            refine=1
            coarsen=-1
    endif


    if ( any (((x(ix^S,1) > 2.5d0)) .and. ((x(ix^S,1) < 3.5d0)) .and. ((abs(x(ix^S,2)-0.5d0*dpi) < 0.95d0))) ) then
            refine=1
            coarsen=-1
    endif

    if ( any (((x(ix^S,1) > 3.5d0)) .and. ((x(ix^S,1) < 4.5d0)) .and. ((abs(x(ix^S,2)-0.5d0*dpi) < 0.75d0))) ) then
            refine=1
            coarsen=-1
    endif




    ! COARSEN HIGHER LATITUDES
    if ( (any(abs(x(ix^S,2)-0.5d0*dpi) > 0.5d0)) .AND. (any(x(ix^S,1) > 5.0d0))  ) then
            coarsen=1
            refine=-1
    endif



   end subroutine specialrefine_grid

end module mod_usr
