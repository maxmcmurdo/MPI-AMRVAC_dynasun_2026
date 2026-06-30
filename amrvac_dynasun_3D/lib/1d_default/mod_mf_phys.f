!> Magnetofriction module
module mod_mf_phys
  use mod_global_parameters, only: std_len
  implicit none
  private

  !> viscosity coefficient in s cm^-2 for solar corona (Cheung 2012 ApJ)
  double precision, public                :: mf_nu = 1.d-15

  !> maximal limit of magnetofrictional velocity in cm s^-1 (Pomoell 2019 A&A)
  double precision, public                :: mf_vmax = 3.d6

  !> decay scale of frictional velocity 
  double precision, public                :: mf_decay_scale(2*1)=0.d0

  !> Whether particles module is added
  logical, public, protected              :: mf_particles = .false.

  !> Whether GLM-MHD is used
  logical, public, protected              :: mf_glm = .false.

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: mf_glm_alpha = 0.5d0

  !> MHD fourth order
  logical, public, protected              :: mf_4th_order = .false.

  !> set to true if need to record electric field on cell edges
  logical, public, protected              :: mf_record_electric_field = &
     .false.

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the magnetic field
  integer, allocatable, public, protected :: mag(:)

  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> The resistivity
  double precision, public                :: mf_eta = 0.0d0

  !> The hyper-resistivity
  double precision, public                :: mf_eta_hyper = 0.0d0

  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'ct'

  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'average'

  !> Whether divB is computed with a fourth order approximation
  logical, public, protected :: mf_divb_4thorder = .false.

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'

  !> Use a compact way to add resistivity
  logical :: compactres   = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*1)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*1)=0

  !> clean divb in the initial condition
  logical, public, protected :: clean_initial_divb=.false.

  ! DivB cleaning methods
  integer, parameter :: divb_none          = 0
  integer, parameter :: divb_multigrid     = -1
  integer, parameter :: divb_glm           = 1
  integer, parameter :: divb_powel         = 2
  integer, parameter :: divb_janhunen      = 3
  integer, parameter :: divb_linde         = 4
  integer, parameter :: divb_lindejanhunen = 5
  integer, parameter :: divb_lindepowel    = 6
  integer, parameter :: divb_lindeglm      = 7
  integer, parameter :: divb_ct            = 8


  ! Public methods
  public :: mf_phys_init
  public :: mf_get_v
  public :: mf_get_v_idim
  public :: mf_to_conserved
  public :: mf_to_primitive
  public :: mf_face_to_center
  public :: get_divb
  public :: get_current
  public :: get_normalized_divb
  public :: b_from_vector_potential
  public :: mf_mag_en_all
  public :: record_force_free_metrics
  

contains

  !> Read this module's parameters from a file
  subroutine mf_read_params(files)
    use mod_global_parameters
    use mod_particles, only: particles_eta
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mf_list/ mf_nu, mf_vmax, mf_decay_scale, mf_eta, mf_eta_hyper,&
        mf_glm_alpha, mf_particles,particles_eta, mf_record_electric_field,&
       mf_4th_order, typedivbfix, source_split_divb, divbdiff,typedivbdiff,&
        type_ct, compactres, divbwave, He_abundance, SI_unit, Bdip, Bquad,&
        Boct, Busr, clean_initial_divb, boundary_divbfix,&
        boundary_divbfix_skip, mf_divb_4thorder

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mf_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine mf_read_params

  !> Write this module's parameters to a snapsoht
  subroutine mf_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "nu"
    values(1) = mf_nu
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine mf_write_info

  subroutine mf_angmomfix(fC,x,wnew,ixImin1,ixImax1,ixOmin1,ixOmax1,idim)
    use mod_global_parameters
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout)    :: fC(ixImin1:ixImax1,1:nwflux,1:ndim),&
         wnew(ixImin1:ixImax1,1:nw)
    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    integer, intent(in)                :: idim
    integer                            :: hxOmin1,hxOmax1, kxCmin1,kxCmax1, iw
    double precision                   :: inv_volume(ixImin1:ixImax1)

    call mpistop("to do")

  end subroutine mf_angmomfix

  subroutine mf_phys_init()
    use mod_global_parameters
    use mod_physics
    use mod_particles, only: particles_init, particles_eta, particles_etah
    

    integer :: itr, idir

    call mf_read_params(par_files)

    physics_type = "mf"
    phys_energy = .false.
    use_particles=mf_particles

    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    
    case ('glm')
      mf_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_glm
    case ('powel', 'powell')
      type_divb = divb_powel
    case ('janhunen')
      type_divb = divb_janhunen
    case ('linde')
      type_divb = divb_linde
    case ('lindejanhunen')
      type_divb = divb_lindejanhunen
    case ('lindepowel')
      type_divb = divb_lindepowel
    case ('lindeglm')
      mf_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_lindeglm
    case ('ct')
      type_divb = divb_ct
      stagger_grid = .true.
    case default
      call mpistop('Unknown divB fix')
    end select

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! set velocity field as flux variables
    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! set magnetic field as flux variables
    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)

    ! start with magnetic field and skip velocity when update ghostcells
    iwstart=mag(1)
    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=mag(1)

    if (mf_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    ! set number of variables which need update ghostcells
    nwgc=nwflux-ndir

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    ! determine number of stagger variables
    if(stagger_grid) nws=ndim

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    do idir=1,ndir
      if(ndim>1) flux_type(idir,mag(idir))=flux_tvdlf
    end do
    if(mf_glm .and. ndim>1) flux_type(:,psi_)=flux_tvdlf

    phys_get_dt              => mf_get_dt
    phys_get_cmax            => mf_get_cmax
    phys_get_cbounds         => mf_get_cbounds
    phys_get_flux            => mf_get_flux
    phys_get_v               => mf_get_v
    phys_add_source_geom     => mf_add_source_geom
    phys_add_source          => mf_add_source
    phys_to_conserved        => mf_to_conserved
    phys_to_primitive        => mf_to_primitive
    phys_check_params        => mf_check_params
    phys_write_info          => mf_write_info
    phys_angmomfix           => mf_angmomfix
    phys_special_advance     => mf_velocity_update

    if(type_divb==divb_glm) then
      phys_modify_wLR => mf_modify_wLR
    end if

    ! pass to global variable to record electric field
    record_electric_field=mf_record_electric_field

    ! Initialize particles module
    if(mf_particles) then
      call particles_init()
      phys_req_diagonal = .true.
      ! never allow Hall effects in particles when doing magnetofrictional
      particles_etah=0.0d0
      if(mype==0) then
         write(*,*) '*****Using particles:     with mf_eta, mf_eta_hyper :',&
             mf_eta, mf_eta_hyper
         write(*,*) '*****Using particles: particles_eta, particles_etah :',&
             particles_eta, particles_etah
      end if
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_get_ct_velocity => mf_get_ct_velocity
      phys_update_faces => mf_update_faces
      phys_face_to_center => mf_face_to_center
    else if(ndim>1) then
      phys_boundary_adjust => mf_boundary_adjust
    end if

    

    ! derive units from basic units
    call mf_physical_units()

  end subroutine mf_phys_init

  subroutine mf_check_params
    use mod_global_parameters

  end subroutine mf_check_params

  subroutine mf_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0,c_lightspeed
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
      c_lightspeed=c_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi
      c_lightspeed=const_c
    end if
    if(unit_density/=1.d0) then
      unit_numberdensity=unit_density/((1.d0+4.d0*He_abundance)*mp)
    else
      ! unit of numberdensity is independent by default
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
    end if
    if(unit_velocity/=1.d0) then
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+&
         3.d0*He_abundance)*unit_numberdensity*kB)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    else if(unit_pressure/=1.d0) then
      unit_temperature=unit_pressure/((2.d0+&
         3.d0*He_abundance)*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    else if(unit_magneticfield/=1.d0) then
      unit_pressure=unit_magneticfield**2/miu0
      unit_temperature=unit_pressure/((2.d0+&
         3.d0*He_abundance)*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
    else
      ! unit of temperature is independent by default
      unit_pressure=(2.d0+3.d0*He_abundance)&
         *unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    end if
    if(unit_time/=1.d0) then
      unit_length=unit_time*unit_velocity
    else
      ! unit of length is independent by default
      unit_time=unit_length/unit_velocity
    end if
    ! Additional units needed for the particles
    c_norm=c_lightspeed/unit_velocity
    unit_charge=unit_magneticfield*unit_length**2/unit_velocity/miu0
    if (.not. SI_unit) unit_charge = unit_charge*const_c
    unit_mass=unit_density*unit_length**3

    ! get dimensionless mf nu
    mf_nu=mf_nu/unit_time*unit_length**2

    ! get dimensionless maximal mf velocity limit
    mf_vmax=mf_vmax/unit_velocity
 
  end subroutine mf_physical_units

  !> Transform primitive variables into conservative ones
  subroutine mf_to_conserved(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    ! nothing to do for mf
  end subroutine mf_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine mf_to_primitive(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    ! nothing to do for mf
  end subroutine mf_to_primitive

  !> Calculate v vector
  subroutine mf_get_v(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1,nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ndir)

    integer :: idir

    do idir=1,ndir
      v(ixOmin1:ixOmax1,idir) = w(ixOmin1:ixOmax1, mom(idir))
    end do

  end subroutine mf_get_v

  !> Calculate v component
  subroutine mf_get_v_idim(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1)

    v(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1, mom(idim))

  end subroutine mf_get_v_idim

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mf_get_cmax(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in) :: w(ixImin1:ixImax1, nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1)

    cmax(ixOmin1:ixOmax1)=abs(w(ixOmin1:ixOmax1,mom(idim)))+one

  end subroutine mf_get_cmax

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine mf_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImax1,ixOmin1,ixOmax1,&
     idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    use mod_variables

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1, nw),&
        wRC(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1, nw),&
        wRp(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       1:number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,&
       1:number_species)

    double precision, dimension(ixImin1:ixImax1) :: tmp1

    tmp1(ixOmin1:ixOmax1)=0.5d0*(wLC(ixOmin1:ixOmax1,&
       mom(idim))+wRC(ixOmin1:ixOmax1,mom(idim)))
    if(present(cmin)) then
      cmax(ixOmin1:ixOmax1,1)=max(tmp1(ixOmin1:ixOmax1)+one,zero)
      cmin(ixOmin1:ixOmax1,1)=min(tmp1(ixOmin1:ixOmax1)-one,zero)
    else
      cmax(ixOmin1:ixOmax1,1)=abs(tmp1(ixOmin1:ixOmax1))+one
    end if

  end subroutine mf_get_cbounds

  !> prepare velocities for ct methods
  subroutine mf_get_ct_velocity(vcts,wLp,wRp,ixImin1,ixImax1,ixOmin1,ixOmax1,&
     idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)    :: wLp(ixImin1:ixImax1, nw),&
        wRp(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: cmax(ixImin1:ixImax1)
    double precision, intent(in), optional :: cmin(ixImin1:ixImax1)
    type(ct_velocity), intent(inout):: vcts

    integer                         :: idimE,idimN

    ! calculate velocities related to different UCT schemes
    select case(type_ct)
    case('average')
    case('uct_contact')
      if(.not.allocated(vcts%vnorm)) allocate(vcts%vnorm(ixImin1:ixImax1,&
         1:ndim))
      ! get average normal velocity at cell faces
      vcts%vnorm(ixOmin1:ixOmax1,idim)=0.5d0*(wLp(ixOmin1:ixOmax1,&
         mom(idim))+wRp(ixOmin1:ixOmax1,mom(idim)))
    case('uct_hll')
      if(.not.allocated(vcts%vbarC)) then
        allocate(vcts%vbarC(ixImin1:ixImax1,1:ndir,2),&
           vcts%vbarLC(ixImin1:ixImax1,1:ndir,2),vcts%vbarRC(ixImin1:ixImax1,&
           1:ndir,2))
        allocate(vcts%cbarmin(ixImin1:ixImax1,1:ndim),&
           vcts%cbarmax(ixImin1:ixImax1,1:ndim)) 
      end if
      ! Store magnitude of characteristics
      if(present(cmin)) then
        vcts%cbarmin(ixOmin1:ixOmax1,idim)=max(-cmin(ixOmin1:ixOmax1),zero)
        vcts%cbarmax(ixOmin1:ixOmax1,idim)=max( cmax(ixOmin1:ixOmax1),zero)
      else
        vcts%cbarmax(ixOmin1:ixOmax1,idim)=max( cmax(ixOmin1:ixOmax1),zero)
        vcts%cbarmin(ixOmin1:ixOmax1,idim)=vcts%cbarmax(ixOmin1:ixOmax1,idim)
      end if

      idimN=mod(idim,ndir)+1 ! 'Next' direction
      idimE=mod(idim+1,ndir)+1 ! Electric field direction
      ! Store velocities
      vcts%vbarLC(ixOmin1:ixOmax1,idim,1)=wLp(ixOmin1:ixOmax1,mom(idimN))
      vcts%vbarRC(ixOmin1:ixOmax1,idim,1)=wRp(ixOmin1:ixOmax1,mom(idimN))
      vcts%vbarC(ixOmin1:ixOmax1,idim,1)=(vcts%cbarmax(ixOmin1:ixOmax1,&
         idim)*vcts%vbarLC(ixOmin1:ixOmax1,idim,&
         1) +vcts%cbarmin(ixOmin1:ixOmax1,idim)*vcts%vbarRC(ixOmin1:ixOmax1,&
         idim,1))/(vcts%cbarmax(ixOmin1:ixOmax1,&
         idim)+vcts%cbarmin(ixOmin1:ixOmax1,idim))

      vcts%vbarLC(ixOmin1:ixOmax1,idim,2)=wLp(ixOmin1:ixOmax1,mom(idimE))
      vcts%vbarRC(ixOmin1:ixOmax1,idim,2)=wRp(ixOmin1:ixOmax1,mom(idimE))
      vcts%vbarC(ixOmin1:ixOmax1,idim,2)=(vcts%cbarmax(ixOmin1:ixOmax1,&
         idim)*vcts%vbarLC(ixOmin1:ixOmax1,idim,&
         2) +vcts%cbarmin(ixOmin1:ixOmax1,idim)*vcts%vbarRC(ixOmin1:ixOmax1,&
         idim,1))/(vcts%cbarmax(ixOmin1:ixOmax1,&
         idim)+vcts%cbarmin(ixOmin1:ixOmax1,idim))
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine mf_get_ct_velocity

  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine mf_get_p_total(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,p)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out)   :: p(ixImin1:ixImax1)

    p(ixOmin1:ixOmax1) = 0.5d0 * sum(w(ixOmin1:ixOmax1, mag(:))**2,&
        dim=ndim+1)

  end subroutine mf_get_p_total

  !> Calculate fluxes within ixO^L.
  subroutine mf_get_flux(wC,w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,f)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,nwflux)

    integer                      :: idir

    ! flux of velocity field is zero, frictional velocity is given in addsource2
    f(ixOmin1:ixOmax1,mom(:))=0.d0

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mf_glm) then
           f(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,psi_)
        else
           f(ixOmin1:ixOmax1,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,&
           mom(idim))*w(ixOmin1:ixOmax1,mag(idir))-w(ixOmin1:ixOmax1,&
           mag(idim))*w(ixOmin1:ixOmax1,mom(idir))
      end if
    end do

    if (mf_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,mag(idim))
    end if

  end subroutine mf_get_flux

  !> Add global source terms to update frictional velocity on complete domain
  subroutine mf_velocity_update(qt,psa)
    use mod_global_parameters
    use mod_ghostcells_update
    double precision, intent(in) :: qt     !< Current time
    type(state), target :: psa(max_blocks) !< Compute based on this state

    integer :: iigrid,igrid
    logical :: stagger_flag
    logical :: firstmf=.true.

    if(firstmf) then
      ! point bc mpi datatype to partial type for velocity field
      type_send_srl=>type_send_srl_p1
      type_recv_srl=>type_recv_srl_p1
      type_send_r=>type_send_r_p1
      type_recv_r=>type_recv_r_p1
      type_send_p=>type_send_p_p1
      type_recv_p=>type_recv_p_p1
      call create_bc_mpi_datatype(mom(1),ndir)
      firstmf=.false.
    end if

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      block=>psa(igrid)
      call frictional_velocity(psa(igrid)%w,psa(igrid)%x,ixGlo1,ixGhi1,ixMlo1,&
         ixMhi1)
    end do
    !$OMP END PARALLEL DO

    ! only update mf velocity in ghost cells
    stagger_flag=stagger_grid
    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1
    bcphys=.false.
    stagger_grid=.false.
    call getbc(qt,0.d0,psa,mom(1),ndir,.true.)
    bcphys=.true.
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    stagger_grid=stagger_flag

  end subroutine mf_velocity_update

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine mf_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
     qsourcesplit,active,wCTprim)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active
    double precision, intent(in), optional :: wCTprim(ixImin1:ixImax1,1:nw)

    if (.not. qsourcesplit) then

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(mf_eta)>smalldouble)then
        active = .true.
        call add_source_res2(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
      end if

      if (mf_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
      end if

    end if

      

  end subroutine mf_add_source

  subroutine frictional_velocity(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: xmin(ndim),xmax(ndim)
    double precision :: decay(ixOmin1:ixOmax1)
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3),&
       tmp(ixOmin1:ixOmax1)
    integer :: ix1, idirmin,idir,jdir,kdir

    call get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)
    ! extrapolate current for the outmost layer, 
    ! because extrapolation of B at boundaries introduces artificial current

    if(block%is_physical_boundary(2*1-1)) then
      current(ixOmin1,:)=current(ixOmin1+1,:)
    end if
    if(block%is_physical_boundary(2*1)) then
      current(ixOmax1,:)=current(ixOmax1-1,:)
    end if

    w(ixOmin1:ixOmax1,mom(:))=0.d0
    ! calculate Lorentz force
    do idir=1,ndir; do jdir=idirmin,3; do kdir=1,ndir
      if(lvc(idir,jdir,kdir)/=0) then
        if(lvc(idir,jdir,kdir)==1) then
          w(ixOmin1:ixOmax1,mom(idir))=w(ixOmin1:ixOmax1,&
             mom(idir))+current(ixOmin1:ixOmax1,jdir)*w(ixOmin1:ixOmax1,&
             mag(kdir))
        else
          w(ixOmin1:ixOmax1,mom(idir))=w(ixOmin1:ixOmax1,&
             mom(idir))-current(ixOmin1:ixOmax1,jdir)*w(ixOmin1:ixOmax1,&
             mag(kdir))
        end if
      end if
    end do; end do; end do

    tmp(ixOmin1:ixOmax1)=sum(w(ixOmin1:ixOmax1,mag(:))**2,dim=ndim+1)
    ! frictional coefficient
    where(tmp(ixOmin1:ixOmax1)/=0.d0)
      tmp(ixOmin1:ixOmax1)=1.d0/(tmp(ixOmin1:ixOmax1)*mf_nu)
    endwhere

    do idir=1,ndir
      w(ixOmin1:ixOmax1,mom(idir))=w(ixOmin1:ixOmax1,&
         mom(idir))*tmp(ixOmin1:ixOmax1)
    end do

    ! decay frictional velocity near selected boundaries
    xmin(1)=xprobmin1
    xmax(1)=xprobmax1
    decay(ixOmin1:ixOmax1)=1.d0
    do idir=1,ndim
      if(mf_decay_scale(2*idir-1)>0.d0) decay(ixOmin1:ixOmax1)=decay(&
         ixOmin1:ixOmax1)*(1.d0-exp(-(x(ixOmin1:ixOmax1,&
         idir)-xmin(idir))/mf_decay_scale(2*idir-1)))
      if(mf_decay_scale(2*idir)>0.d0) decay(ixOmin1:ixOmax1)=decay(&
         ixOmin1:ixOmax1)*(1.d0-exp((x(ixOmin1:ixOmax1,&
         idir)-xmax(idir))/mf_decay_scale(2*idir)))
    end do

    ! saturate mf velocity at mf_vmax
    tmp(ixOmin1:ixOmax1)=sqrt(sum(w(ixOmin1:ixOmax1,mom(:))**2,&
       dim=ndim+1))/mf_vmax+1.d-12
    tmp(ixOmin1:ixOmax1)=dtanh(tmp(ixOmin1:ixOmax1))/tmp(ixOmin1:ixOmax1)
    do idir=1,ndir
      w(ixOmin1:ixOmax1,mom(idir))=w(ixOmin1:ixOmax1,&
         mom(idir))*tmp(ixOmin1:ixOmax1)*decay(ixOmin1:ixOmax1)
    end do

  end subroutine frictional_velocity

  !> Add resistive source to w within ixO Uses 3 point stencil (1 neighbour) in
  !> each direction, non-conservative. If the fourthorder precompiler flag is
  !> set, uses fourth order central difference for the laplacian. Then the
  !> stencil is 5 (2 neighbours).
  subroutine add_source_res1(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    integer :: ixAmin1,ixAmax1,idir,jdir,kdir,idirmin,idim,jxOmin1,jxOmax1,&
       hxOmin1,hxOmax1,ix
    integer :: lxOmin1,lxOmax1, kxOmin1,kxOmax1

    double precision :: tmp(ixImin1:ixImax1),tmp2(ixImin1:ixImax1)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3),&
       eta(ixImin1:ixImax1)
    double precision :: gradeta(ixImin1:ixImax1,1:ndim), Bf(ixImin1:ixImax1,&
       1:ndir)

    ! Calculating resistive sources involve one extra layer
    if (mf_4th_order) then
      ixAmin1=ixOmin1-2;ixAmax1=ixOmax1+2;
    else
      ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
    end if

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1) call &
       mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)

    if (mf_eta>zero)then
       eta(ixAmin1:ixAmax1)=mf_eta
       gradeta(ixOmin1:ixOmax1,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixImin1,ixImax1,ixAmin1,ixAmax1,&
          idirmin,x,current,eta)
       ! assumes that eta is not function of current?
       do idim=1,ndim
          call gradient(eta,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,tmp)
          gradeta(ixOmin1:ixOmax1,idim)=tmp(ixOmin1:ixOmax1)
       end do
    end if

    Bf(ixImin1:ixImax1,1:ndir)=wCT(ixImin1:ixImax1,mag(1:ndir))

    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (mf_4th_order) then
         tmp(ixOmin1:ixOmax1)=zero
         tmp2(ixImin1:ixImax1)=Bf(ixImin1:ixImax1,idir)
         do idim=1,ndim
            lxOmin1=ixOmin1+2*kr(idim,1);lxOmax1=ixOmax1+2*kr(idim,1);
            jxOmin1=ixOmin1+kr(idim,1);jxOmax1=ixOmax1+kr(idim,1);
            hxOmin1=ixOmin1-kr(idim,1);hxOmax1=ixOmax1-kr(idim,1);
            kxOmin1=ixOmin1-2*kr(idim,1);kxOmax1=ixOmax1-2*kr(idim,1);
            tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)+(-tmp2(lxOmin1:lxOmax1)+&
               16.0d0*tmp2(jxOmin1:jxOmax1)-30.0d0*tmp2(ixOmin1:ixOmax1)+&
               16.0d0*tmp2(hxOmin1:hxOmax1)-tmp2(kxOmin1:kxOmax1)) /(12.0d0 * &
               dxlevel(idim)**2)
         end do
       else
         tmp(ixOmin1:ixOmax1)=zero
         tmp2(ixImin1:ixImax1)=Bf(ixImin1:ixImax1,idir)
         do idim=1,ndim
            jxOmin1=ixOmin1+kr(idim,1);jxOmax1=ixOmax1+kr(idim,1);
            hxOmin1=ixOmin1-kr(idim,1);hxOmax1=ixOmax1-kr(idim,1);
            tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)+&
               (tmp2(jxOmin1:jxOmax1)-2.0d0*tmp2(ixOmin1:ixOmax1)+&
               tmp2(hxOmin1:hxOmax1))/dxlevel(idim)**2
         end do
       end if

       ! Multiply by eta
       tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)*eta(ixOmin1:ixOmax1)

       ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
       if (mf_eta<zero)then
          do jdir=1,ndim; do kdir=idirmin,3
             if (lvc(idir,jdir,kdir)/=0)then
                if (lvc(idir,jdir,kdir)==1)then
                   tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)-&
                      gradeta(ixOmin1:ixOmax1,jdir)*current(ixOmin1:ixOmax1,&
                      kdir)
                else
                   tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)+&
                      gradeta(ixOmin1:ixOmax1,jdir)*current(ixOmin1:ixOmax1,&
                      kdir)
                end if
             end if
          end do; end do
       end if

       ! Add sources related to eta*laplB-grad(eta) x J to B and e
       w(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,&
          mag(idir))+qdt*tmp(ixOmin1:ixOmax1)

    end do ! idir

  end subroutine add_source_res1

  !> Add resistive source to w within ixO
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_res2(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3),&
       eta(ixImin1:ixImax1),curlj(ixImin1:ixImax1,1:3)
    double precision :: tmpvec(ixImin1:ixImax1,1:3)
    integer :: ixAmin1,ixAmax1,idir,idirmin,idirmin1

    ! resistive term already added as an electric field component
    if(stagger_grid .and. ndim==ndir) return

    ixAmin1=ixOmin1-2;ixAmax1=ixOmax1+2;

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1) call &
       mpistop("Error in add_source_res2: Non-conforming input limits")

    ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixImin1,ixImax1,ixAmin1,ixAmax1,idirmin,current)

    if (mf_eta>zero)then
       eta(ixAmin1:ixAmax1)=mf_eta
    else
       call usr_special_resistivity(wCT,ixImin1,ixImax1,ixAmin1,ixAmax1,&
          idirmin,x,current,eta)
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    tmpvec(ixAmin1:ixAmax1,1:ndir)=zero
    do idir=idirmin,3
       tmpvec(ixAmin1:ixAmax1,idir)=current(ixAmin1:ixAmax1,&
          idir)*eta(ixAmin1:ixAmax1)
    end do
    curlj=0.d0
    call curlvector(tmpvec,ixImin1,ixImax1,ixOmin1,ixOmax1,curlj,idirmin1,1,3)
    if(stagger_grid) then
      if(ndim==2.and.ndir==3) then
        ! if 2.5D
        w(ixOmin1:ixOmax1,mag(ndir)) = w(ixOmin1:ixOmax1,&
           mag(ndir))-qdt*curlj(ixOmin1:ixOmax1,ndir)
      end if
    else
      w(ixOmin1:ixOmax1,mag(1:ndir)) = w(ixOmin1:ixOmax1,&
         mag(1:ndir))-qdt*curlj(ixOmin1:ixOmax1,1:ndir)
    end if

  end subroutine add_source_res2

  !> Add Hyper-resistive source to w within ixO
  !> Uses 9 point stencil (4 neighbours) in each direction.
  subroutine add_source_hyperres(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    !.. local ..
    double precision                :: current(ixImin1:ixImax1,7-2*ndir:3)
    double precision                :: tmpvec(ixImin1:ixImax1,1:3),&
       tmpvec2(ixImin1:ixImax1,1:3),tmp(ixImin1:ixImax1),&
       ehyper(ixImin1:ixImax1,1:3)
    integer                         :: ixAmin1,ixAmax1,idir,jdir,kdir,idirmin,&
       idirmin1

    ixAmin1=ixOmin1-3;ixAmax1=ixOmax1+3;
    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1) call &
       mpistop("Error in add_source_hyperres: Non-conforming input limits")

    call get_current(wCT,ixImin1,ixImax1,ixAmin1,ixAmax1,idirmin,current)
    tmpvec(ixAmin1:ixAmax1,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixAmin1:ixAmax1,jdir)=current(ixAmin1:ixAmax1,jdir)
    end do

    ixAmin1=ixOmin1-2;ixAmax1=ixOmax1+2;
    call curlvector(tmpvec,ixImin1,ixImax1,ixAmin1,ixAmax1,tmpvec2,idirmin1,1,&
       3)

    ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
    tmpvec(ixAmin1:ixAmax1,1:ndir)=zero
    call curlvector(tmpvec2,ixImin1,ixImax1,ixAmin1,ixAmax1,tmpvec,idirmin1,1,&
       3)
    ehyper(ixAmin1:ixAmax1,1:ndir) = - tmpvec(ixAmin1:ixAmax1,&
       1:ndir)*mf_eta_hyper

    ixAmin1=ixOmin1;ixAmax1=ixOmax1;
    tmpvec2(ixAmin1:ixAmax1,1:ndir)=zero
    call curlvector(ehyper,ixImin1,ixImax1,ixAmin1,ixAmax1,tmpvec2,idirmin1,1,&
       3)

    do idir=1,ndir
      w(ixOmin1:ixOmax1,mag(idir)) = w(ixOmin1:ixOmax1,&
         mag(idir))-tmpvec2(ixOmin1:ixOmax1,idir)*qdt
    end do

  end subroutine add_source_hyperres

  subroutine add_source_glm(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision:: divb(ixImin1:ixImax1)
    integer          :: idim,idir
    double precision :: gradPsi(ixImin1:ixImax1)

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImax1,ixOmin1,ixOmax1,divb, mf_divb_4thorder)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mf_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,psi_) = abs(mf_glm_alpha)*wCT(ixOmin1:ixOmax1,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixOmin1:ixOmax1,psi_) = dexp(-qdt*cmax_global*mf_glm_alpha/minval(&
           dxlevel(:)))*w(ixOmin1:ixOmax1,psi_)
      else
        w(ixOmin1:ixOmax1,psi_) = dexp(-qdt*cmax_global*mf_glm_alpha/minval(&
           block%ds(ixOmin1:ixOmax1,:),dim=ndim+1))*w(ixOmin1:ixOmax1,psi_)
      end if
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
          call gradient(wCT(ixImin1:ixImax1,psi_),ixImin1,ixImax1,ixOmin1,&
             ixOmax1,idim,gradPsi)
       case("limited")
          call gradientS(wCT(ixImin1:ixImax1,psi_),ixImin1,ixImax1,ixOmin1,&
             ixOmax1,idim,gradPsi)
       end select
    end do

  end subroutine add_source_glm

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision                :: divb(ixImin1:ixImax1),v(ixImin1:ixImax1,&
       1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImax1,ixOmin1,ixOmax1,divb, mf_divb_4thorder)

    ! calculate velocity
    call mf_get_v(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,v)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,&
         mag(idir))-qdt*v(ixOmin1:ixOmax1,idir)*divb(ixOmin1:ixOmax1)
    end do

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision                :: divb(ixImin1:ixImax1),v(ixImin1:ixImax1,&
       1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImax1,ixOmin1,ixOmax1,divb, mf_divb_4thorder)

    ! calculate velocity
    call mf_get_v(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,v)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,&
         mag(idir))-qdt*v(ixOmin1:ixOmax1,idir)*divb(ixOmin1:ixOmax1)
    end do

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    ! Add Linde's divB related sources to wnew within ixO
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt, wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    integer :: idim, idir, ixpmin1,ixpmax1, i1, iside
    double precision :: divb(ixImin1:ixImax1),graddivb(ixImin1:ixImax1)
    logical, dimension(-1:1) :: leveljump

    ! Calculate div B
    ixpmin1=ixOmin1-1;ixpmax1=ixOmax1+1;
    call get_divb(wCT,ixImin1,ixImax1,ixpmin1,ixpmax1,divb, mf_divb_4thorder)

    ! for AMR stability, retreat one cell layer from the boarders of level jump
    do i1=-1,1
      if(i1==0) cycle
      if(neighbor_type(i1,block%igrid)==2 .or. neighbor_type(i1,&
         block%igrid)==4) then
        leveljump(i1)=.true.
      else
        leveljump(i1)=.false.
      end if
    end do

    ixpmin1=ixOmin1;ixpmax1=ixOmax1;
    do idim=1,ndim
      select case(idim)
       case(1)
          do iside=1,2
            i1=kr(1,1)*(2*iside-3);
            if (leveljump(i1)) then
              if (iside==1) then
                ixpmin1=ixOmin1-i1
              else
                ixpmax1=ixOmax1-i1
              end if
            end if
          end do
       
      end select
    end do

    ! Add Linde's diffusive terms
    do idim=1,ndim
       ! Calculate grad_idim(divb)
       select case(typegrad)
       case("central")
         call gradient(divb,ixImin1,ixImax1,ixpmin1,ixpmax1,idim,graddivb)
       case("limited")
         call gradientS(divb,ixImin1,ixImax1,ixpmin1,ixpmax1,idim,graddivb)
       end select

       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab_uniform) then
          graddivb(ixpmin1:ixpmax1)=graddivb(ixpmin1:ixpmax1)*divbdiff/(&
             1.0d0/dxlevel(1)**2)
       else
          graddivb(ixpmin1:ixpmax1)=graddivb(ixpmin1:ixpmax1)*divbdiff &
             /(1.0d0/block%ds(ixpmin1:ixpmax1,1)**2)
       end if

       w(ixpmin1:ixpmax1,mag(idim))=w(ixpmin1:ixpmax1,&
          mag(idim))+graddivb(ixpmin1:ixpmax1)
    end do

  end subroutine add_source_linde

  !> Calculate div B within ixO
  subroutine get_divb(w,ixImin1,ixImax1,ixOmin1,ixOmax1,divb, fourthorder)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(inout) :: divb(ixImin1:ixImax1)
    logical, intent(in), optional   :: fourthorder

    integer                            :: ixCmin1,ixCmax1, idir

    if(stagger_grid) then
      divb(ixOmin1:ixOmax1)=0.d0
      do idir=1,ndim
        ixCmin1=ixOmin1-kr(idir,1);ixCmax1=ixOmax1-kr(idir,1);
        divb(ixOmin1:ixOmax1)=divb(ixOmin1:ixOmax1)+block%ws(ixOmin1:ixOmax1,&
           idir)*block%surfaceC(ixOmin1:ixOmax1,idir)-block%ws(ixCmin1:ixCmax1,&
           idir)*block%surfaceC(ixCmin1:ixCmax1,idir)
      end do
      divb(ixOmin1:ixOmax1)=divb(ixOmin1:ixOmax1)/block%dvolume(&
         ixOmin1:ixOmax1)
    else
      select case(typediv)
      case("central")
        call divvector(w(ixImin1:ixImax1,mag(1:ndir)),ixImin1,ixImax1,ixOmin1,&
           ixOmax1,divb,fourthorder)
      case("limited")
        call divvectorS(w(ixImin1:ixImax1,mag(1:ndir)),ixImin1,ixImax1,ixOmin1,&
           ixOmax1,divb)
      end select
    end if

  end subroutine get_divb

  !> get dimensionless div B = |divB| * volume / area / |B|
  subroutine get_normalized_divb(w,ixImin1,ixImax1,ixOmin1,ixOmax1,divb)

    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)       :: w(ixImin1:ixImax1,1:nw)
    double precision                   :: divb(ixImin1:ixImax1),&
        dsurface(ixImin1:ixImax1)

    double precision :: invB(ixOmin1:ixOmax1)
    integer :: ixAmin1,ixAmax1,idims

    call get_divb(w,ixImin1,ixImax1,ixOmin1,ixOmax1,divb)
    invB(ixOmin1:ixOmax1)=sqrt(mf_mag_en_all(w,ixImin1,ixImax1,ixOmin1,&
       ixOmax1))
    where(invB(ixOmin1:ixOmax1)/=0.d0)
      invB(ixOmin1:ixOmax1)=1.d0/invB(ixOmin1:ixOmax1)
    end where
    if(slab_uniform) then
      divb(ixOmin1:ixOmax1)=0.5d0*abs(divb(ixOmin1:ixOmax1))*invB(&
         ixOmin1:ixOmax1)/sum(1.d0/dxlevel(:))
    else
      ixAmin1=ixOmin1-1;
      ixAmax1=ixOmax1-1;
      dsurface(ixOmin1:ixOmax1)= sum(block%surfaceC(ixOmin1:ixOmax1,:),&
         dim=ndim+1)
      do idims=1,ndim
        ixAmin1=ixOmin1-kr(idims,1);ixAmax1=ixOmax1-kr(idims,1);
        dsurface(ixOmin1:ixOmax1)=dsurface(ixOmin1:ixOmax1)+&
           block%surfaceC(ixAmin1:ixAmax1,idims)
      end do
      divb(ixOmin1:ixOmax1)=abs(divb(ixOmin1:ixOmax1))*invB(&
         ixOmin1:ixOmax1)*block%dvolume(ixOmin1:ixOmax1)/dsurface(&
         ixOmin1:ixOmax1)
    end if

  end subroutine get_normalized_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current,&
     fourthorder)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)  :: ixOmin1,ixOmax1, ixImin1,ixImax1
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw)
    integer, intent(out) :: idirmin
    logical, intent(in), optional :: fourthorder
    integer :: idir, idirmin0

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3)

    idirmin0 = 7-2*ndir

    if(present(fourthorder)) then
      call curlvector(w(ixImin1:ixImax1,mag(1:ndir)),ixImin1,ixImax1,ixOmin1,&
         ixOmax1,current,idirmin,idirmin0,ndir,fourthorder)
    else
      call curlvector(w(ixImin1:ixImax1,mag(1:ndir)),ixImin1,ixImax1,ixOmin1,&
         ixOmax1,current,idirmin,idirmin0,ndir)
    end if

  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine mf_get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx1
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)

    integer                       :: idirmin,idim
    double precision              :: dxarr(ndim)
    double precision              :: current(ixImin1:ixImax1,7-2*ndir:3),&
       eta(ixImin1:ixImax1)

    dtnew = bigdouble

    dxarr(1)=dx1;
    if (mf_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/mf_eta
    else if (mf_eta<zero)then
       call get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)
       call usr_special_resistivity(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,&
          x,current,eta)
       dtnew=bigdouble
       do idim=1,ndim
         if(slab_uniform) then
           dtnew=min(dtnew,dtdiffpar/(smalldouble+&
              maxval(eta(ixOmin1:ixOmax1)/dxarr(idim)**2)))
         else
           dtnew=min(dtnew,dtdiffpar/(smalldouble+&
              maxval(eta(ixOmin1:ixOmax1)/block%ds(ixOmin1:ixOmax1,idim)**2)))
         end if
       end do
    end if

    if(mf_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/mf_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixOmin1:ixOmax1,&
           1:ndim))**4/mf_eta_hyper,dtnew)
      end if
    end if
  end subroutine mf_get_dt

  ! Add geometrical source terms to w
  subroutine mf_add_source_geom(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,1:nw),&
        w(ixImin1:ixImax1,1:nw)

    integer          :: iw,idir
    double precision :: tmp(ixImin1:ixImax1)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    select case (coordinate)
    case (cylindrical)
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
      call mf_get_p_total(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,tmp)
      if(phi_>0) then
        if(.not.stagger_grid) then
          w(ixOmin1:ixOmax1,bphi_)=w(ixOmin1:ixOmax1,&
             bphi_)+qdt/x(ixOmin1:ixOmax1,1)*(wCT(ixOmin1:ixOmax1,&
             bphi_)*wCT(ixOmin1:ixOmax1,mr_) -wCT(ixOmin1:ixOmax1,&
             br_)*wCT(ixOmin1:ixOmax1,mphi_))
        end if
      end if
      if(mf_glm) w(ixOmin1:ixOmax1,br_)=w(ixOmin1:ixOmax1,&
         br_)+qdt*wCT(ixOmin1:ixOmax1,psi_)/x(ixOmin1:ixOmax1,1)
    case (spherical)
       ! b1
       if(mf_glm) then
         w(ixOmin1:ixOmax1,mag(1))=w(ixOmin1:ixOmax1,&
            mag(1))+qdt/x(ixOmin1:ixOmax1,1)*2.0d0*wCT(ixOmin1:ixOmax1,psi_)
       end if

       

       if(ndir==3) then
         ! b3
         if(.not.stagger_grid) then
           tmp(ixOmin1:ixOmax1)=(wCT(ixOmin1:ixOmax1,&
              mom(1))*wCT(ixOmin1:ixOmax1,mag(3)) -wCT(ixOmin1:ixOmax1,&
              mom(3))*wCT(ixOmin1:ixOmax1,mag(1))) 
           w(ixOmin1:ixOmax1,mag(3))=w(ixOmin1:ixOmax1,&
              mag(3))+qdt*tmp(ixOmin1:ixOmax1)/x(ixOmin1:ixOmax1,1)
         end if
       end if
    end select

  end subroutine mf_add_source_geom

  !> Compute 2 times total magnetic energy
  function mf_mag_en_all(w, ixImin1,ixImax1, ixOmin1,ixOmax1) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: mge(ixOmin1:ixOmax1)

    mge = sum(w(ixOmin1:ixOmax1, mag(:))**2, dim=ndim+1)
  end function mf_mag_en_all

  !> Compute full magnetic field by direction
  function mf_mag_i_all(w, ixImin1,ixImax1, ixOmin1,ixOmax1,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1, idir
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: mgf(ixOmin1:ixOmax1)

    mgf = w(ixOmin1:ixOmax1, mag(idir))
  end function mf_mag_i_all

  !> Compute evolving magnetic energy
  function mf_mag_en(w, ixImin1,ixImax1, ixOmin1,ixOmax1) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: mge(ixOmin1:ixOmax1)

    mge = 0.5d0 * sum(w(ixOmin1:ixOmax1, mag(:))**2, dim=ndim+1)
  end function mf_mag_en

  subroutine mf_modify_wLR(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,wLC,wRC,wLp,wRp,&
     s,idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixImin1:ixImax1,1:nw),&
        wRC(ixImin1:ixImax1,1:nw)
    double precision, intent(inout) :: wLp(ixImin1:ixImax1,1:nw),&
        wRp(ixImin1:ixImax1,1:nw)
    type(state)                     :: s
    double precision                :: dB(ixImin1:ixImax1),&
        dPsi(ixImin1:ixImax1)

    ! Solve the Riemann problem for the linear 2x2 system for normal
    ! B-field and GLM_Psi according to Dedner 2002:
    ! This implements eq. (42) in Dedner et al. 2002 JcP 175
    ! Gives the Riemann solution on the interface
    ! for the normal B component and Psi in the GLM-MHD system.
    ! 23/04/2013 Oliver Porth
    dB(ixOmin1:ixOmax1)   = wRp(ixOmin1:ixOmax1,&
       mag(idir)) - wLp(ixOmin1:ixOmax1,mag(idir))
    dPsi(ixOmin1:ixOmax1) = wRp(ixOmin1:ixOmax1,psi_) - wLp(ixOmin1:ixOmax1,&
       psi_)

    wLp(ixOmin1:ixOmax1,mag(idir))   = 0.5d0 * (wRp(ixOmin1:ixOmax1,&
       mag(idir)) + wLp(ixOmin1:ixOmax1,mag(idir))) - 0.5d0/cmax_global * &
       dPsi(ixOmin1:ixOmax1)
    wLp(ixOmin1:ixOmax1,psi_)       = 0.5d0 * (wRp(ixOmin1:ixOmax1,&
       psi_) + wLp(ixOmin1:ixOmax1,psi_)) - 0.5d0*cmax_global * &
       dB(ixOmin1:ixOmax1)

    wRp(ixOmin1:ixOmax1,mag(idir)) = wLp(ixOmin1:ixOmax1,mag(idir))
    wRp(ixOmin1:ixOmax1,psi_) = wLp(ixOmin1:ixOmax1,psi_)

  end subroutine mf_modify_wLR

  subroutine mf_boundary_adjust(igrid,psb)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state), target :: psb(max_blocks)

    integer :: iB, idims, iside, ixOmin1,ixOmax1, i1

    block=>ps(igrid)
    dxlevel(1)=rnode(rpdx1_,igrid);
    do idims=1,ndim
       ! to avoid using as yet unknown corner info in more than 1D, we
       ! fill only interior mesh ranges of the ghost cell ranges at first,
       ! and progressively enlarge the ranges to include corners later
       do iside=1,2
          i1=kr(1,idims)*(2*iside-3);
          if (neighbor_type(i1,igrid)/=1) cycle
          iB=(idims-1)*2+iside
          if(.not.boundary_divbfix(iB)) cycle
          if(any(typeboundary(:,iB)==bc_special)) then
            ! MF nonlinear force-free B field extrapolation and data driven
            ! require normal B of the first ghost cell layer to be untouched by
            ! fixdivB=0 process, set boundary_divbfix_skip(iB)=1 in par file
            select case (idims)
            case (1)
               if (iside==2) then
                  ! maximal boundary
                  ixOmin1=ixGhi1+1-nghostcells+boundary_divbfix_skip(2*1);
                  ixOmax1=ixGhi1;
               else
                  ! minimal boundary
                  ixOmin1=ixGlo1;
                  ixOmax1=ixGlo1-1+nghostcells-boundary_divbfix_skip(2*1-1);
               end if 
            end select
            call fixdivB_boundary(ixGlo1,ixGhi1,ixOmin1,ixOmax1,psb(igrid)%w,&
               psb(igrid)%x,iB)
          end if
       end do
    end do

  end subroutine mf_boundary_adjust

  subroutine fixdivB_boundary(ixGmin1,ixGmax1,ixOmin1,ixOmax1,w,x,iB)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmax1,ixOmin1,ixOmax1,iB
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)

    double precision :: dx1x2,dx1x3,dx2x1,dx2x3,dx3x1,dx3x2
    integer :: ix1,ixFmin1,ixFmax1

    select case(iB)
     case(1)
       ! 2nd order CD for divB=0 to set normal B component better
       
       
     case(2)
       
       
     case(3)
       
       
     case(4)
       
       
     
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  

  subroutine mf_update_faces(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,qdt,wprim,fC,&
     fE,sCT,s,vcts)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixImin1:ixImax1,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixImin1:ixImax1,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,7-2*ndim:3)

    select case(type_ct)
    case('average')
      call update_faces_average(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,qdt,fC,fE,&
         sCT,s)
    case('uct_contact')
      call update_faces_contact(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,qdt,wprim,&
         fC,fE,sCT,s,vcts)
    case('uct_hll')
      call update_faces_hll(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,qdt,fE,sCT,s,&
         vcts)
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine mf_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,qdt,fC,fE,&
     sCT,s)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixImin1:ixImax1,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,7-2*ndim:3)

    integer                            :: hxCmin1,hxCmax1,ixCmin1,ixCmax1,&
       jxCmin1,jxCmax1,ixCmmin1,ixCmmax1
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2
    double precision                   :: circ(ixImin1:ixImax1,1:ndim)
    ! current on cell edges
    double precision :: jce(ixImin1:ixImax1,7-2*ndim:3)

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! if there is resistivity, get eta J
    if(mf_eta/=zero) call get_resistive_electric_field(ixImin1,ixImax1,ixOmin1,&
       ixOmax1,sCT,s,jce)

    do idim1=1,ndim 
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=7-2*ndim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax1=ixOmax1;
            ixCmin1=ixOmin1+kr(idir,1)-1;
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmax1=ixCmax1+kr(idim1,1);
            hxCmin1=ixCmin1+kr(idim2,1);hxCmax1=ixCmax1+kr(idim2,1);
            ! Interpolate to edges
            fE(ixCmin1:ixCmax1,idir)=quarter*(fC(ixCmin1:ixCmax1,iwdim1,&
               idim2)+fC(jxCmin1:jxCmax1,iwdim1,idim2)-fC(ixCmin1:ixCmax1,&
               iwdim2,idim1)-fC(hxCmin1:hxCmax1,iwdim2,idim1))

            ! add current component of electric field at cell edges E=-vxB+eta J
            if(mf_eta/=zero) fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
               idir)+jce(ixCmin1:ixCmax1,idir)

            if(record_electric_field) s%we(ixCmin1:ixCmax1,&
               idir)=fE(ixCmin1:ixCmax1,idir)
            ! times time step and edge length 
            fE(ixCmin1:ixCmax1,idir)=qdt*s%dsC(ixCmin1:ixCmax1,&
               idir)*fE(ixCmin1:ixCmax1,idir)

          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImax1,ixOmin1,ixOmax1,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,1:ndim)=zero

    ! Calculate circulation on each face

    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax1=ixOmax1;
      ixCmin1=ixOmin1-kr(idim1,1);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxCmin1=ixCmin1-kr(idim2,1);hxCmax1=ixCmax1-kr(idim2,1);
            ! Add line integrals in direction idir
            circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,idim1)+lvc(idim1,&
               idim2,idir)*(fE(ixCmin1:ixCmax1,idir)-fE(hxCmin1:hxCmax1,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,&
         idim1)/s%surfaceC(ixCmin1:ixCmax1,idim1)
      ! Time update
      bfaces(ixCmin1:ixCmax1,idim1)=bfaces(ixCmin1:ixCmax1,&
         idim1)-circ(ixCmin1:ixCmax1,idim1)
    end do

    end associate

  end subroutine update_faces_average

  !> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
  subroutine update_faces_contact(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,qdt,wp,fC,&
     fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixImin1:ixImax1,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixImin1:ixImax1,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,7-2*ndim:3)

    double precision                   :: circ(ixImin1:ixImax1,1:ndim)
    ! electric field at cell centers
    double precision                   :: ECC(ixImin1:ixImax1,7-2*ndim:3)
    ! gradient of E at left and right side of a cell face
    double precision                   :: EL(ixImin1:ixImax1),&
       ER(ixImin1:ixImax1)
    ! gradient of E at left and right side of a cell corner
    double precision                   :: ELC(ixImin1:ixImax1),&
       ERC(ixImin1:ixImax1)
    ! current on cell edges
    double precision                   :: jce(ixImin1:ixImax1,7-2*ndim:3)
    integer                            :: hxCmin1,hxCmax1,ixCmin1,ixCmax1,&
       jxCmin1,jxCmax1,ixAmin1,ixAmax1,ixBmin1,ixBmax1
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm)

    ECC=0.d0
    ! Calculate electric field at cell centers
    do idim1=1,ndim; do idim2=1,ndim; do idir=7-2*ndim,3
      if(lvc(idim1,idim2,idir)==1)then
         ECC(ixImin1:ixImax1,idir)=ECC(ixImin1:ixImax1,&
            idir)+wp(ixImin1:ixImax1,mag(idim1))*wp(ixImin1:ixImax1,&
            mom(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixImin1:ixImax1,idir)=ECC(ixImin1:ixImax1,&
            idir)-wp(ixImin1:ixImax1,mag(idim1))*wp(ixImin1:ixImax1,&
            mom(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(mf_eta/=zero) call get_resistive_electric_field(ixImin1,ixImax1,ixOmin1,&
       ixOmax1,sCT,s,jce)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ! evaluate electric field along cell edges according to equation (41)
    do idim1=1,ndim 
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax1=ixOmax1;
            ixCmin1=ixOmin1+kr(idir,1)-1;
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmax1=ixCmax1+kr(idim1,1);
            hxCmin1=ixCmin1+kr(idim2,1);hxCmax1=ixCmax1+kr(idim2,1);
            ! average cell-face electric field to cell edges
            fE(ixCmin1:ixCmax1,idir)=quarter*(fC(ixCmin1:ixCmax1,iwdim1,&
               idim2)+fC(jxCmin1:jxCmax1,iwdim1,idim2)-fC(ixCmin1:ixCmax1,&
               iwdim2,idim1)-fC(hxCmin1:hxCmax1,iwdim2,idim1))

            ! add slope in idim2 direction from equation (50)
            ixAmin1=ixCmin1;
            ixAmax1=ixCmax1+kr(idim1,1);
            EL(ixAmin1:ixAmax1)=fC(ixAmin1:ixAmax1,iwdim1,&
               idim2)-ECC(ixAmin1:ixAmax1,idir)
            hxCmin1=ixAmin1+kr(idim2,1);hxCmax1=ixAmax1+kr(idim2,1);
            ER(ixAmin1:ixAmax1)=fC(ixAmin1:ixAmax1,iwdim1,&
               idim2)-ECC(hxCmin1:hxCmax1,idir)
            where(vnorm(ixCmin1:ixCmax1,idim1)>0.d0)
              ELC(ixCmin1:ixCmax1)=EL(ixCmin1:ixCmax1)
            else where(vnorm(ixCmin1:ixCmax1,idim1)<0.d0)
              ELC(ixCmin1:ixCmax1)=EL(jxCmin1:jxCmax1)
            else where
              ELC(ixCmin1:ixCmax1)=0.5d0*(EL(ixCmin1:ixCmax1)+&
                 EL(jxCmin1:jxCmax1))
            end where
            hxCmin1=ixCmin1+kr(idim2,1);hxCmax1=ixCmax1+kr(idim2,1);
            where(vnorm(hxCmin1:hxCmax1,idim1)>0.d0)
              ERC(ixCmin1:ixCmax1)=ER(ixCmin1:ixCmax1)
            else where(vnorm(hxCmin1:hxCmax1,idim1)<0.d0)
              ERC(ixCmin1:ixCmax1)=ER(jxCmin1:jxCmax1)
            else where
              ERC(ixCmin1:ixCmax1)=0.5d0*(ER(ixCmin1:ixCmax1)+&
                 ER(jxCmin1:jxCmax1))
            end where
            fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
               idir)+0.25d0*(ELC(ixCmin1:ixCmax1)+ERC(ixCmin1:ixCmax1))

            ! add slope in idim1 direction from equation (50)
            jxCmin1=ixCmin1+kr(idim2,1);jxCmax1=ixCmax1+kr(idim2,1);
            ixAmin1=ixCmin1;
            ixAmax1=ixCmax1+kr(idim2,1);
            EL(ixAmin1:ixAmax1)=-fC(ixAmin1:ixAmax1,iwdim2,&
               idim1)-ECC(ixAmin1:ixAmax1,idir)
            hxCmin1=ixAmin1+kr(idim1,1);hxCmax1=ixAmax1+kr(idim1,1);
            ER(ixAmin1:ixAmax1)=-fC(ixAmin1:ixAmax1,iwdim2,&
               idim1)-ECC(hxCmin1:hxCmax1,idir)
            where(vnorm(ixCmin1:ixCmax1,idim2)>0.d0)
              ELC(ixCmin1:ixCmax1)=EL(ixCmin1:ixCmax1)
            else where(vnorm(ixCmin1:ixCmax1,idim2)<0.d0)
              ELC(ixCmin1:ixCmax1)=EL(jxCmin1:jxCmax1)
            else where
              ELC(ixCmin1:ixCmax1)=0.5d0*(EL(ixCmin1:ixCmax1)+&
                 EL(jxCmin1:jxCmax1))
            end where
            hxCmin1=ixCmin1+kr(idim1,1);hxCmax1=ixCmax1+kr(idim1,1);
            where(vnorm(hxCmin1:hxCmax1,idim2)>0.d0)
              ERC(ixCmin1:ixCmax1)=ER(ixCmin1:ixCmax1)
            else where(vnorm(hxCmin1:hxCmax1,idim2)<0.d0)
              ERC(ixCmin1:ixCmax1)=ER(jxCmin1:jxCmax1)
            else where
              ERC(ixCmin1:ixCmax1)=0.5d0*(ER(ixCmin1:ixCmax1)+&
                 ER(jxCmin1:jxCmax1))
            end where
            fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
               idir)+0.25d0*(ELC(ixCmin1:ixCmax1)+ERC(ixCmin1:ixCmax1))

            ! add current component of electric field at cell edges E=-vxB+eta J
            if(mf_eta/=zero) fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
               idir)+jce(ixCmin1:ixCmax1,idir)

            if(record_electric_field) s%we(ixCmin1:ixCmax1,&
               idir)=fE(ixCmin1:ixCmax1,idir)
            ! times time step and edge length 
            fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
               idir)*qdt*s%dsC(ixCmin1:ixCmax1,idir)
            if (.not.slab) then
              where(abs(x(ixCmin1:ixCmax1,r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixCmin1:ixCmax1,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImax1,ixOmin1,ixOmax1,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,1:ndim)=zero

    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax1=ixOmax1;
      ixCmin1=ixOmin1-kr(idim1,1);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxCmin1=ixCmin1-kr(idim2,1);hxCmax1=ixCmax1-kr(idim2,1);
            ! Add line integrals in direction idir
            circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,idim1)+lvc(idim1,&
               idim2,idir)*(fE(ixCmin1:ixCmax1,idir)-fE(hxCmin1:hxCmax1,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixCmin1:ixCmax1,idim1) > &
         1.0d-9*s%dvolume(ixCmin1:ixCmax1))
        circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,&
           idim1)/s%surfaceC(ixCmin1:ixCmax1,idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixCmin1:ixCmax1,idim1)=bfaces(ixCmin1:ixCmax1,&
         idim1)-circ(ixCmin1:ixCmax1,idim1)
    end do

    end associate

  end subroutine update_faces_contact

  !> update faces
  subroutine update_faces_hll(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,qdt,fE,sCT,s,&
     vcts)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)       :: qt, qdt
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,7-2*ndim:3)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts

    double precision                   :: vtilL(ixImin1:ixImax1,2)
    double precision                   :: vtilR(ixImin1:ixImax1,2)
    double precision                   :: btilL(ixImin1:ixImax1,ndim)
    double precision                   :: btilR(ixImin1:ixImax1,ndim)
    double precision                   :: cp(ixImin1:ixImax1,2)
    double precision                   :: cm(ixImin1:ixImax1,2)
    double precision                   :: circ(ixImin1:ixImax1,1:ndim)
    ! current on cell edges
    double precision                   :: jce(ixImin1:ixImax1,7-2*ndim:3)
    integer                            :: hxCmin1,hxCmax1,ixCmin1,ixCmax1,&
       ixCpmin1,ixCpmax1,jxCmin1,jxCmax1,ixCmmin1,ixCmmax1
    integer                            :: idim1,idim2,idir

    associate(bfaces=>s%ws,bfacesCT=>sCT%ws,x=>s%x,vbarC=>vcts%vbarC,&
       cbarmin=>vcts%cbarmin,cbarmax=>vcts%cbarmax)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! Loop over components of electric field

    ! idir: electric field component we need to calculate
    ! idim1: directions in which we already performed the reconstruction
    ! idim2: directions in which we perform the reconstruction

    ! if there is resistivity, get eta J
    if(mf_eta/=zero) call get_resistive_electric_field(ixImin1,ixImax1,ixOmin1,&
       ixOmax1,sCT,s,jce)

    do idir=7-2*ndim,3
      ! Indices
      ! idir: electric field component
      ! idim1: one surface
      ! idim2: the other surface
      ! cyclic permutation: idim1,idim2,idir=1,2,3
      ! Velocity components on the surface
      ! follow cyclic premutations:
      ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

      ixCmax1=ixOmax1;
      ixCmin1=ixOmin1-1+kr(idir,1);

      ! Set indices and directions
      idim1=mod(idir,3)+1
      idim2=mod(idir+1,3)+1

      jxCmin1=ixCmin1+kr(idim1,1);jxCmax1=ixCmax1+kr(idim1,1);
      ixCpmin1=ixCmin1+kr(idim2,1);ixCpmax1=ixCmax1+kr(idim2,1);

      ! Reconstruct transverse transport velocities
      call reconstruct(ixImin1,ixImax1,ixCmin1,ixCmax1,idim2,&
         vbarC(ixImin1:ixImax1,idim1,1),vtilL(ixImin1:ixImax1,2),&
         vtilR(ixImin1:ixImax1,2))

      call reconstruct(ixImin1,ixImax1,ixCmin1,ixCmax1,idim1,&
         vbarC(ixImin1:ixImax1,idim2,2),vtilL(ixImin1:ixImax1,1),&
         vtilR(ixImin1:ixImax1,1))

      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      call reconstruct(ixImin1,ixImax1,ixCmin1,ixCmax1,idim2,&
         bfacesCT(ixImin1:ixImax1,idim1),btilL(ixImin1:ixImax1,idim1),&
         btilR(ixImin1:ixImax1,idim1))

      call reconstruct(ixImin1,ixImax1,ixCmin1,ixCmax1,idim1,&
         bfacesCT(ixImin1:ixImax1,idim2),btilL(ixImin1:ixImax1,idim2),&
         btilR(ixImin1:ixImax1,idim2))

      ! Take the maximum characteristic

      cm(ixCmin1:ixCmax1,1)=max(cbarmin(ixCpmin1:ixCpmax1,idim1),&
         cbarmin(ixCmin1:ixCmax1,idim1))
      cp(ixCmin1:ixCmax1,1)=max(cbarmax(ixCpmin1:ixCpmax1,idim1),&
         cbarmax(ixCmin1:ixCmax1,idim1))

      cm(ixCmin1:ixCmax1,2)=max(cbarmin(jxCmin1:jxCmax1,idim2),&
         cbarmin(ixCmin1:ixCmax1,idim2))
      cp(ixCmin1:ixCmax1,2)=max(cbarmax(jxCmin1:jxCmax1,idim2),&
         cbarmax(ixCmin1:ixCmax1,idim2))
     

      ! Calculate eletric field
      fE(ixCmin1:ixCmax1,idir)=-(cp(ixCmin1:ixCmax1,1)*vtilL(ixCmin1:ixCmax1,&
         1)*btilL(ixCmin1:ixCmax1,idim2) + cm(ixCmin1:ixCmax1,&
         1)*vtilR(ixCmin1:ixCmax1,1)*btilR(ixCmin1:ixCmax1,&
         idim2) - cp(ixCmin1:ixCmax1,1)*cm(ixCmin1:ixCmax1,&
         1)*(btilR(ixCmin1:ixCmax1,idim2)-btilL(ixCmin1:ixCmax1,&
         idim2)))/(cp(ixCmin1:ixCmax1,1)+cm(ixCmin1:ixCmax1,&
         1)) +(cp(ixCmin1:ixCmax1,2)*vtilL(ixCmin1:ixCmax1,&
         2)*btilL(ixCmin1:ixCmax1,idim1) + cm(ixCmin1:ixCmax1,&
         2)*vtilR(ixCmin1:ixCmax1,2)*btilR(ixCmin1:ixCmax1,&
         idim1) - cp(ixCmin1:ixCmax1,2)*cm(ixCmin1:ixCmax1,&
         2)*(btilR(ixCmin1:ixCmax1,idim1)-btilL(ixCmin1:ixCmax1,&
         idim1)))/(cp(ixCmin1:ixCmax1,2)+cm(ixCmin1:ixCmax1,2))

      ! add current component of electric field at cell edges E=-vxB+eta J
      if(mf_eta/=zero) fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
         idir)+jce(ixCmin1:ixCmax1,idir)

      if(record_electric_field) s%we(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
         idir)
      ! times time step and edge length 
      fE(ixCmin1:ixCmax1,idir)=qdt*s%dsC(ixCmin1:ixCmax1,&
         idir)*fE(ixCmin1:ixCmax1,idir)

      if (.not.slab) then
        where(abs(x(ixCmin1:ixCmax1,r_)+half*dxlevel(r_)).lt.1.0d-9)
          fE(ixCmin1:ixCmax1,idir)=zero
        end where
      end if

    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImax1,ixOmin1,ixOmax1,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,1:ndim)=zero

    ! Calculate circulation on each face: interal(fE dot dl)

    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax1=ixOmax1;
      ixCmin1=ixOmin1-kr(idim1,1);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxCmin1=ixCmin1-kr(idim2,1);hxCmax1=ixCmax1-kr(idim2,1);
            ! Add line integrals in direction idir
            circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,idim1)+lvc(idim1,&
               idim2,idir)*(fE(ixCmin1:ixCmax1,idir)-fE(hxCmin1:hxCmax1,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixCmin1:ixCmax1,idim1) > &
         1.0d-9*s%dvolume(ixCmin1:ixCmax1))
        circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,&
           idim1)/s%surfaceC(ixCmin1:ixCmax1,idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,idim1)=zero
      end where
      ! Time update
      bfaces(ixCmin1:ixCmax1,idim1)=bfaces(ixCmin1:ixCmax1,&
         idim1)-circ(ixCmin1:ixCmax1,idim1)
    end do

    end associate
  end subroutine update_faces_hll

  !> calculate eta J at cell edges
  subroutine get_resistive_electric_field(ixImin1,ixImax1,ixOmin1,ixOmax1,sCT,&
     s,jce)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    type(state), intent(in)            :: sCT, s
    ! current on cell edges
    double precision :: jce(ixImin1:ixImax1,7-2*ndim:3)

    ! current on cell centers
    double precision :: jcc(ixImin1:ixImax1,7-2*ndir:3)
    ! location at cell faces
    double precision :: xs(ixGslo1:ixGshi1,1:ndim)
    ! resistivity
    double precision :: eta(ixImin1:ixImax1)
    double precision :: gradi(ixGslo1:ixGshi1)
    integer :: ix1,ixCmin1,ixCmax1,ixAmin1,ixAmax1,ixBmin1,ixBmax1,idir,&
       idirmin,idim1,idim2

    associate(x=>s%x,dx=>s%dx,w=>s%w,wCT=>sCT%w,wCTs=>sCT%ws)
    ! calculate current density at cell edges
    jce=0.d0
    do idim1=1,ndim 
      do idim2=1,ndim
        do idir=7-2*ndim,3
          if (lvc(idim1,idim2,idir)==0) cycle
          ixCmax1=ixOmax1+1;
          ixCmin1=ixOmin1+kr(idir,1)-2;
          ixBmax1=ixCmax1-kr(idir,1)+1;
          ixBmin1=ixCmin1;
          ! current at transverse faces
          xs(ixBmin1:ixBmax1,:)=x(ixBmin1:ixBmax1,:)
          xs(ixBmin1:ixBmax1,idim2)=x(ixBmin1:ixBmax1,&
             idim2)+half*dx(ixBmin1:ixBmax1,idim2)
          call gradientx(wCTs(ixGslo1:ixGshi1,idim2),xs,ixGslo1,ixGshi1,&
             ixCmin1,ixCmax1,idim1,gradi,.false.)
          if (lvc(idim1,idim2,idir)==1) then
            jce(ixCmin1:ixCmax1,idir)=jce(ixCmin1:ixCmax1,&
               idir)+gradi(ixCmin1:ixCmax1)
          else
            jce(ixCmin1:ixCmax1,idir)=jce(ixCmin1:ixCmax1,&
               idir)-gradi(ixCmin1:ixCmax1)
          end if
        end do
      end do
    end do
    ! get resistivity
    if(mf_eta>zero)then
      jce(ixImin1:ixImax1,:)=jce(ixImin1:ixImax1,:)*mf_eta
    else
      ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
      call get_current(wCT,ixImin1,ixImax1,ixAmin1,ixAmax1,idirmin,jcc)
      call usr_special_resistivity(wCT,ixImin1,ixImax1,ixAmin1,ixAmax1,idirmin,&
         x,jcc,eta)
      ! calcuate eta on cell edges
      do idir=7-2*ndim,3
        ixCmax1=ixOmax1;
        ixCmin1=ixOmin1+kr(idir,1)-1;
        jcc(ixCmin1:ixCmax1,idir)=0.d0
       do ix1=0,1
          if( ix1==1 .and. 1==idir ) cycle
          ixAmin1=ixCmin1+ix1;
          ixAmax1=ixCmax1+ix1;
          jcc(ixCmin1:ixCmax1,idir)=jcc(ixCmin1:ixCmax1,&
             idir)+eta(ixAmin1:ixAmax1)
       end do
        jcc(ixCmin1:ixCmax1,idir)=jcc(ixCmin1:ixCmax1,idir)*0.25d0
        jce(ixCmin1:ixCmax1,idir)=jce(ixCmin1:ixCmax1,&
           idir)*jcc(ixCmin1:ixCmax1,idir)
      enddo
    end if

    end associate
  end subroutine get_resistive_electric_field

  !> calculate cell-center values from face-center values
  subroutine mf_face_to_center(ixOmin1,ixOmax1,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixOmin1,ixOmax1
    type(state)                        :: s

    integer                            :: fxOmin1,fxOmax1, gxOmin1,gxOmax1,&
        hxOmin1,hxOmax1, jxOmin1,jxOmax1, kxOmin1,kxOmax1, idim

    associate(w=>s%w, ws=>s%ws)

    ! calculate cell-center values from face-center values in 2nd order
    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the w arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxOmin1=ixOmin1-kr(idim,1);hxOmax1=ixOmax1-kr(idim,1);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      w(ixOmin1:ixOmax1,mag(idim))=half/s%surface(ixOmin1:ixOmax1,&
         idim)*(ws(ixOmin1:ixOmax1,idim)*s%surfaceC(ixOmin1:ixOmax1,&
         idim)+ws(hxOmin1:hxOmax1,idim)*s%surfaceC(hxOmin1:hxOmax1,idim))
    end do

    ! calculate cell-center values from face-center values in 4th order
    !do idim=1,ndim
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);

    !  ! Interpolate to cell barycentre using fourth order central formula
    !  w(ixO^S,mag(idim))=(0.0625d0/s%surface(ixO^S,idim))*&
    !         ( -ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !     +9.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !     +9.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !           -ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) )
    !end do

    ! calculate cell-center values from face-center values in 6th order
    !do idim=1,ndim
    !  fxO^L=ixO^L-3*kr(idim,^D);
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);
    !  kxO^L=ixO^L+2*kr(idim,^D);

    !  ! Interpolate to cell barycentre using sixth order central formula
    !  w(ixO^S,mag(idim))=(0.00390625d0/s%surface(ixO^S,idim))* &
    !     (  +3.0d0*ws(fxO^S,idim)*s%surfaceC(fxO^S,idim) &
    !       -25.0d0*ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !      +150.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !      +150.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !       -25.0d0*ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) &
    !        +3.0d0*ws(kxO^S,idim)*s%surfaceC(kxO^S,idim) )
    !end do

    end associate

  end subroutine mf_face_to_center

  !> calculate magnetic field from vector potential
  subroutine b_from_vector_potential(ixIsmin1,ixIsmax1, ixImin1,ixImax1,&
      ixOmin1,ixOmax1, ws, x)
    use mod_global_parameters
    use mod_constrained_transport

    integer, intent(in)                :: ixIsmin1,ixIsmax1, ixImin1,ixImax1,&
        ixOmin1,ixOmax1
    double precision, intent(inout)    :: ws(ixIsmin1:ixIsmax1,1:nws)
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)

    double precision                   :: Adummy(ixIsmin1:ixIsmax1,1:3)

    call b_from_vector_potentialA(ixIsmin1,ixIsmax1, ixImin1,ixImax1, ixOmin1,&
       ixOmax1, ws, x, Adummy)

  end subroutine b_from_vector_potential

  subroutine record_force_free_metrics()
    use mod_global_parameters

    double precision :: sum_jbb_ipe, sum_j_ipe, sum_l_ipe, f_i_ipe, volume_pe
    double precision :: sum_jbb, sum_j, sum_l, f_i, volume, cw_sin_theta
    integer :: iigrid, igrid, ix1
    integer :: amode, istatus(MPI_STATUS_SIZE)
    integer, save :: fhmf
    character(len=800) :: filename,filehead
    character(len=800) :: line,datastr
    logical :: patchwi(ixGlo1:ixGhi1)
    logical, save :: logmfopened=.false.

    sum_jbb_ipe = 0.d0
    sum_j_ipe = 0.d0
    sum_l_ipe = 0.d0
    f_i_ipe = 0.d0
    volume_pe=0.d0
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      dxlevel(1)=rnode(rpdx1_,igrid);
      call mask_inner(ixGlo1,ixGhi1,ixMlo1,ixMhi1,ps(igrid)%w,ps(igrid)%x,&
         patchwi,volume_pe)
      sum_jbb_ipe = sum_jbb_ipe+integral_grid_mf(ixGlo1,ixGhi1,ixMlo1,ixMhi1,&
         ps(igrid)%w,ps(igrid)%x,1,patchwi)
      sum_j_ipe   = sum_j_ipe+integral_grid_mf(ixGlo1,ixGhi1,ixMlo1,ixMhi1,&
         ps(igrid)%w,ps(igrid)%x,2,patchwi)
      f_i_ipe=f_i_ipe+integral_grid_mf(ixGlo1,ixGhi1,ixMlo1,ixMhi1,ps(igrid)%w,&
         ps(igrid)%x,3,patchwi)
      sum_l_ipe   = sum_l_ipe+integral_grid_mf(ixGlo1,ixGhi1,ixMlo1,ixMhi1,&
         ps(igrid)%w,ps(igrid)%x,4,patchwi)
    end do
    call MPI_REDUCE(sum_jbb_ipe,sum_jbb,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,&
       ierrmpi)
    call MPI_REDUCE(sum_j_ipe,sum_j,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,&
       ierrmpi)
    call MPI_REDUCE(f_i_ipe,f_i,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,&
       ierrmpi)
    call MPI_REDUCE(sum_l_ipe,sum_l,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,&
       ierrmpi)
    call MPI_REDUCE(volume_pe,volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,&
       ierrmpi)

    if(mype==0) then
      ! current- and volume-weighted average of the sine of the angle
      ! between the magnetic field B and the current density J
      cw_sin_theta = sum_jbb/sum_j
      ! volume-weighted average of the absolute value of the fractional
      ! magnetic flux change
      f_i = f_i/volume
      sum_j=sum_j/volume
      sum_l=sum_l/volume
      if(.not.logmfopened) then
        ! generate filename
        write(filename,"(a,a)") TRIM(base_filename), "_mf_metrics.csv"

        ! Delete the log when not doing a restart run
        if (restart_from_file == undefined) then
           open(unit=20,file=trim(filename),status='replace')
           close(20, status='delete')
        end if

        amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        amode=ior(amode,MPI_MODE_APPEND)
        call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,fhmf,&
           ierrmpi)
        logmfopened=.true.
        filehead="  it,time,CW_sin_theta,current,Lorenz_force,f_i"
        if (it == 0) then
          call MPI_FILE_WRITE(fhmf,filehead,len_trim(filehead), MPI_CHARACTER,&
             istatus,ierrmpi)
          call MPI_FILE_WRITE(fhmf,achar(10),1,MPI_CHARACTER,istatus,ierrmpi)
        end if
      end if
      line=''
      write(datastr,'(i6,a)') it,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') global_time,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') cw_sin_theta,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') sum_j,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6,a)') sum_l,','
      line=trim(line)//trim(datastr)
      write(datastr,'(es13.6)') f_i
      line=trim(line)//trim(datastr)//new_line('A')
      call MPI_FILE_WRITE(fhmf,line,len_trim(line),MPI_CHARACTER,istatus,&
         ierrmpi)
      if(.not.time_advance) then
        call MPI_FILE_CLOSE(fhmf,ierrmpi)
        logmfopened=.false.
      end if
    end if

  end subroutine record_force_free_metrics

  subroutine mask_inner(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,patchwi,volume_pe)
    use mod_global_parameters

    integer, intent(in)         :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in):: w(ixImin1:ixImax1,nw),x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(inout) :: volume_pe
    logical, intent(out)        :: patchwi(ixImin1:ixImax1)

    double precision            :: xOmin1,xOmax1
    integer                     :: ix1

    xOmin1 = xprobmin1 + 0.05d0*(xprobmax1-xprobmin1)
    xOmax1 = xprobmax1 - 0.05d0*(xprobmax1-xprobmin1)
    if(slab) then
      xOmin1 = xprobmin1
    else
      xOmin1 = xprobmin1
    end if

    do ix1=ixOmin1,ixOmax1
        if( x(ix1,1) > xOmin1 .and. x(ix1,1) < xOmax1 ) then
          patchwi(ix1)=.true.
          volume_pe=volume_pe+block%dvolume(ix1)
        else
          patchwi(ix1)=.false.
        endif
    end do

  end subroutine mask_inner

  function integral_grid_mf(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,iw,patchwi)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1,iw
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    double precision                   :: w(ixImin1:ixImax1,nw+nwauxio)
    logical, intent(in) :: patchwi(ixImin1:ixImax1)

    double precision, dimension(ixImin1:ixImax1,7-2*ndir:3) :: current
    double precision, dimension(ixImin1:ixImax1,1:ndir) :: bvec,qvec
    double precision :: integral_grid_mf,tmp(ixImin1:ixImax1),bm2
    integer :: ix1,idirmin0,idirmin,idir,jdir,kdir

    integral_grid_mf=0.d0
    select case(iw)
     case(1)
      ! Sum(|JxB|/|B|*dvolume)
      bvec(ixOmin1:ixOmax1,:)=w(ixOmin1:ixOmax1,mag(:))
      call get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)
      ! calculate Lorentz force
      qvec(ixOmin1:ixOmax1,1:ndir)=zero
      do idir=1,ndir; do jdir=idirmin,3; do kdir=1,ndir
         if(lvc(idir,jdir,kdir)/=0)then
           if(lvc(idir,jdir,kdir)==1)then
             qvec(ixOmin1:ixOmax1,idir)=qvec(ixOmin1:ixOmax1,&
                idir)+current(ixOmin1:ixOmax1,jdir)*bvec(ixOmin1:ixOmax1,kdir)
           else
             qvec(ixOmin1:ixOmax1,idir)=qvec(ixOmin1:ixOmax1,&
                idir)-current(ixOmin1:ixOmax1,jdir)*bvec(ixOmin1:ixOmax1,kdir)
           end if
         end if
      end do; end do; end do

      do ix1=ixOmin1,ixOmax1
         if(patchwi(ix1)) then
           bm2=sum(bvec(ix1,:)**2)
           if(bm2/=0.d0) bm2=1.d0/bm2
           integral_grid_mf=integral_grid_mf+sqrt(sum(qvec(ix1,&
              :)**2)*bm2)*block%dvolume(ix1)
         end if
      end do
     case(2)
      ! Sum(|J|*dvolume)
      call get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)
      do ix1=ixOmin1,ixOmax1
         if(patchwi(ix1)) integral_grid_mf=integral_grid_mf+&
            sqrt(sum(current(ix1,:)**2))*block%dvolume(ix1)
      end do
     case(3)
      ! f_i solenoidal property of B: (dvolume |div B|)/(dsurface |B|)
      ! Sum(f_i*dvolume)
      call get_normalized_divb(w,ixImin1,ixImax1,ixOmin1,ixOmax1,tmp)
      do ix1=ixOmin1,ixOmax1
         if(patchwi(ix1)) integral_grid_mf=integral_grid_mf+&
            tmp(ix1)*block%dvolume(ix1)
      end do
     case(4)
      ! Sum(|JxB|*dvolume)
      bvec(ixOmin1:ixOmax1,:)=w(ixOmin1:ixOmax1,mag(:))
      call get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)
      ! calculate Lorentz force
      qvec(ixOmin1:ixOmax1,1:ndir)=zero
      do idir=1,ndir; do jdir=idirmin,3; do kdir=1,ndir
         if(lvc(idir,jdir,kdir)/=0)then
           if(lvc(idir,jdir,kdir)==1)then
             qvec(ixOmin1:ixOmax1,idir)=qvec(ixOmin1:ixOmax1,&
                idir)+current(ixOmin1:ixOmax1,jdir)*bvec(ixOmin1:ixOmax1,kdir)
           else
             qvec(ixOmin1:ixOmax1,idir)=qvec(ixOmin1:ixOmax1,&
                idir)-current(ixOmin1:ixOmax1,jdir)*bvec(ixOmin1:ixOmax1,kdir)
           end if
         end if
      end do; end do; end do

      do ix1=ixOmin1,ixOmax1
         if(patchwi(ix1)) integral_grid_mf=integral_grid_mf+sqrt(sum(qvec(ix1,&
            :)**2))*block%dvolume(ix1)
      end do
    end select
    return
  end function integral_grid_mf

end module mod_mf_phys
