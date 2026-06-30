!> Module for including rotating frame in hydrodynamics simulations
!> The rotation vector is assumed to be along z direction
!>(both in cylindrical and spherical)

module mod_rotating_frame
  implicit none

  !> Rotation frequency of the frame
  double precision :: omega_frame

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1
  
contains
  !> Read this module's parameters from a file
  subroutine rotating_frame_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    
    namelist /rotating_frame_list/ omega_frame
    
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, rotating_frame_list, end=111)
111    close(unitpar)
    end do
    
  end subroutine rotating_frame_params_read
  
  !> Initialize the module
  subroutine rotating_frame_init()
    use mod_global_parameters
    integer :: nwx,idir
    
    call rotating_frame_params_read(par_files)
    
  end subroutine rotating_frame_init
  
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine rotating_frame_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer                         :: idir

    ! .. local ..

    double precision                :: rotating_terms(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3),frame_omega(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    select case (coordinate)
    case (cylindrical)
       rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) =  omega_frame**2 * x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_) * wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_rho)
       if (phi_ > 0) then
          rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + 2.d0 * omega_frame *wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mom(phi_))
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(r_)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(r_)) + qdt * rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       
       if(phi_>0 .and. .not. angmomfix) then
          rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)   = - two * omega_frame * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mom(r_))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(phi_)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, iw_mom(phi_)) + qdt * &
             rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       end if
    case (spherical)
       !call mpistop("Rotating frame not implemented yet with spherical geometries")
       ! source[mrad] = 2mphi*vframe/r + rho*vframe**2/r

       call rotating_frame_omega(x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,frame_omega)
       rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)**2 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,r_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,iw_rho)
        if (phi_ > 0) then
           rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) = rotating_terms(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3) + 2.d0 * &
              frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,iw_mom(phi_))
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(r_)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(r_)) + qdt * rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)

       
       ! source[mtheta] = cot(theta)*(2mphi*vframe/r + rho*vframe**2/r )
       rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) =  frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)**2 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,r_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,iw_rho)
       if (phi_>0 ) then
          rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + 2.d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, iw_mom(phi_))* frame_omega(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(2)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(2)) + qdt * rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)/ tan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, 2))

       if (phi_>0  .and. .not. angmomfix) then
          ! source[mphi]=-2*mr*vframe/r-2*cot(theta)*mtheta*vframe/r
          rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = -2.d0*frame_omega(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)* wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(r_))- 2.d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, iw_mom(2)) * frame_omega(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)/ tan(x(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, 2))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(3)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(3)) + qdt * rotating_terms(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       end if
      
    case default
       call mpistop("Rotating frame not implemented in this geometries")
    end select
    
  end subroutine rotating_frame_add_source


  subroutine rotating_frame_velocity(x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,frame_vel)
        use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: frame_vel(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    frame_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       r_) * omega_frame
    
    if (coordinate == spherical .and. ndim > 1) then
       frame_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = frame_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,2))
    end if
   

  end subroutine rotating_frame_velocity

    subroutine rotating_frame_omega(x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,frame_omega)
        use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: frame_omega(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) =  omega_frame
    
    if (coordinate == spherical .and. ndim > 1) then
       frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)* dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,2))
    end if
   

  end subroutine rotating_frame_omega
end module mod_rotating_frame
