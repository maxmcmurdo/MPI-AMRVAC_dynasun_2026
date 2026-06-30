!> Module for including gravity in (magneto)hydrodynamics simulations
module mod_gravity
  implicit none

  !> source split or not
  logical :: grav_split= .false.

contains
  !> Read this module's parameters from a file
  subroutine grav_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grav_list/ grav_split

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grav_list, end=111)
111    close(unitpar)
    end do

  end subroutine grav_params_read

  !> Initialize the module
  subroutine gravity_init()
    use mod_global_parameters
    integer :: nwx,idir

    call grav_params_read(par_files)

  end subroutine gravity_init

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine gravity_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
     energy,qsourcesplit,active)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
    integer                         :: idim

    double precision :: gravity_field(ixImin1:ixImax1,ndim)

    if(qsourcesplit .eqv. grav_split) then
      active = .true.

      if (.not. associated(usr_gravity)) then
        write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
        write(*,*) "like the phys_gravity in mod_usr_methods.t"
        call mpistop("gravity_add_source: usr_gravity not defined")
      else
        call usr_gravity(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,gravity_field)
      end if
  
      do idim = 1, ndim
        w(ixOmin1:ixOmax1,iw_mom(idim)) = w(ixOmin1:ixOmax1,&
           iw_mom(idim)) + qdt * gravity_field(ixOmin1:ixOmax1,&
           idim) * wCT(ixOmin1:ixOmax1,iw_rho)
        if(energy) then
          if(iw_equi_rho>0) then
            w(ixOmin1:ixOmax1,iw_e)=w(ixOmin1:ixOmax1,&
               iw_e) + qdt * gravity_field(ixOmin1:ixOmax1,&
               idim) *wCT(ixOmin1:ixOmax1,iw_mom(idim))/(wCT(ixOmin1:ixOmax1,&
               iw_rho)+block%equi_vars(ixOmin1:ixOmax1,iw_equi_rho,&
               0))*wCT(ixOmin1:ixOmax1,iw_rho)

          else
            w(ixOmin1:ixOmax1,iw_e)=w(ixOmin1:ixOmax1,&
               iw_e) + qdt * gravity_field(ixOmin1:ixOmax1,&
               idim) * wCT(ixOmin1:ixOmax1,iw_mom(idim))
          endif
        end if
      end do
    end if

  end subroutine gravity_add_source

  subroutine gravity_get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: dx1, x(ixImin1:ixImax1,1:ndim),&
        w(ixImin1:ixImax1,1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: dxinv(1:ndim), max_grav
    integer                         :: idim

    double precision :: gravity_field(ixImin1:ixImax1,ndim)

    dxinv(1)=one/dx1;

    if(.not. associated(usr_gravity)) then
      write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
      write(*,*) "like the phys_gravity in mod_usr_methods.t"
      call mpistop("gravity_get_dt: usr_gravity not defined")
    else
      call usr_gravity(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,gravity_field)
    end if

    do idim = 1, ndim
      max_grav = maxval(abs(gravity_field(ixOmin1:ixOmax1,idim)))
      max_grav = max(max_grav, epsilon(1.0d0))
      dtnew = min(dtnew, 1.0d0 / sqrt(max_grav * dxinv(idim)))
    end do

  end subroutine gravity_get_dt

end module mod_gravity
