!> Module UAWSOM, based on mhd
module mod_uawsom_phys

#include "amrvac.h"

  use mod_global_parameters, only: std_len, const_c
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  use mod_physics

  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: uawsom_energy = .true.

  !> Whether thermal conduction is used
  logical, public, protected              :: uawsom_thermal_conduction = &
     .false.
  !> type of fluid for thermal conduction
  type(tc_fluid), public, allocatable     :: tc_fl
  type(te_fluid), public, allocatable     :: te_fl_uawsom

  !> Whether radiative cooling is added
  logical, public, protected              :: uawsom_radiative_cooling = &
     .false.
  !> type of fluid for radiative cooling
  type(rc_fluid), public, allocatable     :: rc_fl

  !> Whether viscosity is added
  logical, public, protected              :: uawsom_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: uawsom_gravity = .false.

  !> Whether Hall-uawsom is used
  logical, public, protected              :: uawsom_Hall = .false.

  !> Whether Ambipolar term is used
  logical, public, protected              :: uawsom_ambipolar = .false.

  !> Whether Ambipolar term is implemented using supertimestepping
  logical, public, protected              :: uawsom_ambipolar_sts = .false.

  !> Whether Ambipolar term is implemented explicitly
  logical, public, protected              :: uawsom_ambipolar_exp = .false.

  !> Whether particles module is added
  logical, public, protected              :: uawsom_particles = .false.

  !> Whether magnetofriction is added
  logical, public, protected              :: uawsom_magnetofriction = .false.

  !> Whether GLM-uawsom is used to control div B
  logical, public, protected              :: uawsom_glm = .false.

  !> Whether extended GLM-uawsom is used with additional sources
  logical, public, protected              :: uawsom_glm_extended = .true.

  !> Whether TRAC method is used
  logical, public, protected              :: uawsom_trac = .false.

  !> Which TRAC method is used
  integer, public, protected              :: uawsom_trac_type=1

  !> Height of the mask used in the TRAC method
  double precision, public, protected     :: uawsom_trac_mask = 0.d0

  !> Distance between two adjacent traced magnetic field lines (in finest cell size)
  integer, public, protected              :: uawsom_trac_finegrid=4

  !> Whether auxiliary internal energy is solved
  logical, public, protected              :: uawsom_solve_eaux = .false.

  !> Whether internal energy is solved instead of total energy
  logical, public, protected              :: uawsom_internal_e = .false.

  !TODO this does not work with the splitting: check uawsom_check_w_hde and uawsom_handle_small_values_hde
  !> Whether hydrodynamic energy is solved instead of total energy
  logical, public, protected              :: uawsom_hydrodynamic_e = .false.

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> GLM-uawsom parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: uawsom_glm_alpha = 0.5d0

  !> Whether boris simplified semirelativistic uawsom equations (Gombosi 2002 JCP) are solved
  logical, public, protected              :: uawsom_boris_simplification = &
     .false.

  !> Reduced speed of light for semirelativistic uawsom
  double precision, public, protected     :: uawsom_reduced_c = const_c

  !> Whether CAK radiation line force is activated
  logical, public, protected              :: uawsom_cak_force = .false.

  !> uawsom fourth order
  logical, public, protected              :: uawsom_4th_order = .false.

  !> whether split off equilibrium density
  logical, public :: has_equi_rho0 = .false.
  !> whether split off equilibrium thermal pressure
  logical, public :: has_equi_pe0 = .false.
  logical, public :: uawsom_equi_thermal = .false.

  !> equi vars indices in the state%equi_vars array
  integer, public :: equi_rho0_ = -1
  integer, public :: equi_pe0_ = -1

  !> whether dump full variables (when splitting is used) in a separate dat file
  logical, public, protected              :: uawsom_dump_full_vars = .false.

  !> Number of tracer species
  integer, public, protected              :: uawsom_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> Indices of the magnetic field
  integer, allocatable, public, protected :: mag(:)

  !> Index of wminus
  !integer, public, protected              :: wminus_

  !> Index of wplus
  !integer, public, protected              :: wplus_

  ! -------------------Adding AWs to UAWSOM---------------------------
  ! Max: Requires distinction between kink and Alfven waves

  !> Index of wkminus
  integer, public, protected              :: wkminus_

  !> Index of wkplus
  integer, public, protected              :: wkplus_

  !> Index of wAminus
  integer, public, protected              :: wAminus_

  !> Index of wAplus
  integer, public, protected              :: wAplus_

  ! --------------Back to original mod_uawsom_phys----------------------


  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> Indices of auxiliary internal energy
  integer, public, protected :: eaux_
  integer, public, protected :: paux_

  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_
  integer, public, protected              :: Tweight_

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> The adiabatic index
  double precision, public                :: uawsom_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: uawsom_adiab = 1.0d0

  !> The uawsom resistivity
  double precision, public                :: uawsom_eta = 0.0d0

  !> The uawsom hyper-resistivity
  double precision, public                :: uawsom_eta_hyper = 0.0d0

  !> TODO: what is this?
  double precision, public                :: uawsom_etah = 0.0d0

  !> The uawsom ambipolar coefficient
  double precision, public                :: uawsom_eta_ambi = 0.0d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> The number of waves
  integer :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'linde'

  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'uct_contact'

  !> Whether divB is computed with a fourth order approximation
  logical, public, protected :: uawsom_divb_4thorder = .false.

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

  !> clean initial divB
  logical, public :: clean_initial_divb     = .false.

  !> Base density ratio, zeta0
  double precision, public, protected :: zeta0=10.d0 
  
  !> Filling factor constant
  double precision, public, protected :: ff=0.1d0 

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0
  !> Ionization fraction of H
  !> H_ion_fr = H+/(H+ + H)
  double precision, public, protected  :: H_ion_fr=1d0
  !> Ionization fraction of He
  !> He_ion_fr = (He2+ + He+)/(He2+ + He+ + He)
  double precision, public, protected  :: He_ion_fr=1d0
  !> Ratio of number He2+ / number He+ + He2+
  !> He_ion_fr2 = He2+/(He2+ + He+)
  double precision, public, protected  :: He_ion_fr2=1d0
  ! used for eq of state when it is not defined by units,
  ! the units do not contain terms related to ionization fraction
  ! and it is p = RR * rho * T
  double precision, public, protected  :: RR=1d0
  ! remove the below flag  and assume default value = .false.
  ! when eq state properly implemented everywhere
  ! and not anymore through units
  logical, public, protected :: eq_state_units = .true.

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*1)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*1)=0

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.

  !> Whether an total energy equation is used
  logical :: total_energy = .true.

  !> Whether gravity work is included in energy equation
  logical :: gravity_energy

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

  !> inverse of squared speed of light c0 and reduced speed of light c
  double precision :: inv_squared_c0, inv_squared_c

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

  !define the subroutine interface for the ambipolar mask
  abstract interface

    subroutine mask_subroutine(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,res)
      use mod_global_parameters
      integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in) :: w(ixImin1:ixImax1,1:nw)
      double precision, intent(inout) :: res(ixImin1:ixImax1)
    end subroutine mask_subroutine

    function fun_kin_en(w, ixImin1,ixImax1, ixOmin1,ixOmax1,&
        inv_rho) result(ke)
      use mod_global_parameters, only: nw, ndim,block
      integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
      double precision              :: ke(ixOmin1:ixOmax1)
      double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1)
    end function fun_kin_en

  end interface

  procedure(mask_subroutine), pointer  :: usr_mask_ambipolar => null()
  procedure(sub_get_pthermal), pointer  :: usr_Rfactor => null()
  procedure(sub_convert), pointer      :: uawsom_to_primitive  => null()
  procedure(sub_convert), pointer      :: uawsom_to_conserved  => null()
  procedure(sub_small_values), pointer :: uawsom_handle_small_values => null()
  procedure(sub_get_pthermal), pointer :: uawsom_get_pthermal  => null()
  procedure(sub_get_v), pointer        :: uawsom_get_v         => null()
  procedure(fun_kin_en), pointer       :: uawsom_kin_en        => null()
  ! Public methods
  public :: usr_mask_ambipolar
  public :: usr_Rfactor
  public :: uawsom_phys_init
  public :: uawsom_kin_en
  public :: uawsom_get_pthermal
  public :: uawsom_get_v
  public :: uawsom_get_rho
  public :: uawsom_get_v_idim
  public :: uawsom_to_conserved
  public :: uawsom_to_primitive
  public :: uawsom_get_csound2
  public :: uawsom_e_to_ei
  public :: uawsom_ei_to_e
  public :: uawsom_face_to_center
  public :: get_divb
  public :: get_current
  !> needed  public if we want to use the ambipolar coefficient in the user file
  public :: multiplyAmbiCoef
  public :: get_normalized_divb
  public :: b_from_vector_potential
  public :: uawsom_mag_en_all
  

contains

  !> Read this modules parameters from a file
  subroutine uawsom_read_params(files)
    use mod_global_parameters
    use mod_particles, only: particles_eta, particles_etah
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /uawsom_list/ uawsom_energy, uawsom_n_tracer, uawsom_gamma,&
        uawsom_adiab,uawsom_eta, uawsom_eta_hyper, uawsom_etah,&
        uawsom_eta_ambi, uawsom_glm_alpha, uawsom_glm_extended,&
        uawsom_magnetofriction,uawsom_thermal_conduction,&
        uawsom_radiative_cooling, uawsom_Hall, uawsom_ambipolar,&
        uawsom_ambipolar_sts, uawsom_gravity,uawsom_viscosity,&
        uawsom_4th_order, typedivbfix, source_split_divb, divbdiff,&
       typedivbdiff, type_ct, compactres, divbwave, He_abundance, H_ion_fr,&
        He_ion_fr, He_ion_fr2, eq_state_units, SI_unit, B0field ,&
       uawsom_dump_full_vars,B0field_forcefree, Bdip, Bquad, Boct, Busr,&
        uawsom_particles,particles_eta, particles_etah,has_equi_rho0,&
        has_equi_pe0,uawsom_equi_thermal,boundary_divbfix,&
        boundary_divbfix_skip, uawsom_divb_4thorder,&
       uawsom_boris_simplification, uawsom_reduced_c, clean_initial_divb,&
        uawsom_solve_eaux, uawsom_internal_e, uawsom_hydrodynamic_e,&
        uawsom_trac, uawsom_trac_type, uawsom_trac_mask, uawsom_trac_finegrid,&
        uawsom_cak_force

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, uawsom_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine uawsom_read_params

  !> Write this modules parameters to a snapsoht
  subroutine uawsom_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = uawsom_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine uawsom_write_info

  subroutine uawsom_angmomfix(fC,x,wnew,ixImin1,ixImax1,ixOmin1,ixOmax1,idim)
    use mod_global_parameters
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout)    :: fC(ixImin1:ixImax1,1:nwflux,1:ndim),&
         wnew(ixImin1:ixImax1,1:nw)
    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    integer, intent(in)                :: idim
    integer                            :: hxOmin1,hxOmax1, kxCmin1,kxCmax1, iw
    double precision                   :: inv_volume(ixImin1:ixImax1)

    call mpistop("to do")

  end subroutine uawsom_angmomfix

  subroutine uawsom_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init, particles_eta, particles_etah
    use mod_magnetofriction, only: magnetofriction_init
    use mod_supertimestepping, only: sts_init, add_sts_method,&
       set_conversion_methods_to_head, set_error_handling_to_head
    use mod_cak_force, only: cak_init
    

    integer :: itr, idir

    call uawsom_read_params(par_files)

    if(uawsom_internal_e) then
      if(uawsom_solve_eaux) then
        uawsom_solve_eaux=.false.
        if(mype==0) write(*,*)&
            'WARNING: set uawsom_solve_eaux=F when uawsom_internal_e=T'
      end if
      if(uawsom_hydrodynamic_e) then
        uawsom_hydrodynamic_e=.false.
        if(mype==0) write(*,*)&
            'WARNING: set uawsom_hydrodynamic_e=F when uawsom_internal_e=T'
      end if
    end if


    if(.not. uawsom_energy) then
      if(uawsom_internal_e) then
        uawsom_internal_e=.false.
        if(mype==0) write(*,*)&
            'WARNING: set uawsom_internal_e=F when uawsom_energy=F'
      end if
      if(uawsom_solve_eaux) then
        uawsom_solve_eaux=.false.
        if(mype==0) write(*,*)&
            'WARNING: set uawsom_solve_eaux=F when uawsom_energy=F'
      end if
      if(uawsom_hydrodynamic_e) then
        uawsom_hydrodynamic_e=.false.
        if(mype==0) write(*,*)&
            'WARNING: set uawsom_hydrodynamic_e=F when uawsom_energy=F'
      end if
      if(uawsom_thermal_conduction) then
        uawsom_thermal_conduction=.false.
        if(mype==0) write(*,*)&
            'WARNING: set uawsom_thermal_conduction=F when uawsom_energy=F'
      end if
      if(uawsom_radiative_cooling) then
        uawsom_radiative_cooling=.false.
        if(mype==0) write(*,*)&
            'WARNING: set uawsom_radiative_cooling=F when uawsom_energy=F'
      end if
      if(uawsom_trac) then
        uawsom_trac=.false.
        if(mype==0) write(*,*)&
            'WARNING: set uawsom_trac=F when uawsom_energy=F'
      end if
    end if

    physics_type = "uawsom"
    phys_energy=uawsom_energy
    phys_internal_e=uawsom_internal_e
    phys_solve_eaux=uawsom_solve_eaux
    phys_trac=uawsom_trac
    phys_trac_type=uawsom_trac_type

    phys_gamma = uawsom_gamma

    if(uawsom_energy.and..not.uawsom_internal_e.and.&
       .not.uawsom_hydrodynamic_e) then
      total_energy=.true.
    else
      total_energy=.false.
    end if
    phys_total_energy=total_energy
    if(uawsom_energy) then
      if(uawsom_internal_e) then
        gravity_energy=.false.
      else
        gravity_energy=.true.
      end if
    else
      gravity_energy=.false.
    end if

    
    if(uawsom_trac .and. uawsom_trac_type .gt. 2) then
      uawsom_trac_type=1
      if(mype==0) write(*,*)&
          'WARNING: reset uawsom_trac_type=1 for 1D simulation'
    end if
   
    if(uawsom_trac .and. uawsom_trac_type .le. 4) then
      uawsom_trac_mask=bigdouble
      if(mype==0) write(*,*)&
          'WARNING: set uawsom_trac_mask==bigdouble for global TRAC method'
    end if
    phys_trac_mask=uawsom_trac_mask

    if(uawsom_solve_eaux) prolongprimitive=.true.

    ! set default gamma for polytropic/isothermal process
    use_particles=uawsom_particles
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    
    case ('glm')
      uawsom_glm          = .true.
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
      uawsom_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_lindeglm
    case ('ct')
      type_divb = divb_ct
      stagger_grid = .true.
    case default
      call mpistop('Unknown divB fix')
    end select

    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=1

    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (uawsom_energy) then
      nwwave = 8
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 7
      e_     = -1
      p_     = -1
    end if

    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)

    !add wminus and wplus
    !wplus_ = var_set_wplus()
    !wminus_ = var_set_wminus()

    !add wkminus and wkplus
    wkplus_ = var_set_wkplus()
    wkminus_ = var_set_wkminus()

    !add wAminus and wAplus
    wAplus_ = var_set_wAplus()
    wAminus_ = var_set_wAminus()

    !  set auxiliary internal energy variable
    if(uawsom_energy .and. uawsom_solve_eaux) then
      eaux_ = var_set_internal_energy()
      paux_ = eaux_
    else
      eaux_ = -1
      paux_ = -1
    end if

    if (uawsom_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    allocate(tracer(uawsom_n_tracer))

    ! Set starting index of tracers
    do itr = 1, uawsom_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    ! set cutoff temperature when using the TRAC method, as well as an auxiliary weight
    Tweight_ = -1
    if(uawsom_trac) then
      Tcoff_ = var_set_wextra()
      iw_Tcoff=Tcoff_
      if(uawsom_trac_type .ge. 3) then
        Tweight_ = var_set_wextra()
      endif
    else
      Tcoff_ = -1
    end if

    ! set indices of equi vars and update number_equi_vars
    number_equi_vars = 0
    if(has_equi_rho0) then
      number_equi_vars = number_equi_vars + 1
      equi_rho0_ = number_equi_vars
      iw_equi_rho = equi_rho0_
    endif
    if(has_equi_pe0) then
      number_equi_vars = number_equi_vars + 1
      equi_pe0_ = number_equi_vars
      iw_equi_p = equi_pe0_
    endif
    ! determine number of stagger variables
    if(stagger_grid) nws=ndim

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nwflux))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nwflux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    if(nwflux>mag(ndir)) then
      ! for flux of eaux and tracers, using hll flux
      flux_type(:,mag(ndir)+1:nwflux)=flux_hll
    end if

    if(ndim>1) then
      if(uawsom_glm) then
        flux_type(:,psi_)=flux_special
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_special
        end do
      else
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_tvdlf
        end do
      end if
    end if

    phys_get_dt              => uawsom_get_dt
    phys_get_cmax            => uawsom_get_cmax_origin
    phys_get_a2max           => uawsom_get_a2max
    phys_get_tcutoff         => uawsom_get_tcutoff
    phys_get_H_speed         => uawsom_get_H_speed
    if(has_equi_rho0) then
      phys_get_cbounds         => uawsom_get_cbounds_split_rho
    else
      phys_get_cbounds         => uawsom_get_cbounds
    end if
    if(has_equi_rho0) then
      phys_to_primitive        => uawsom_to_primitive_split_rho
      uawsom_to_primitive         => uawsom_to_primitive_split_rho
      phys_to_conserved        => uawsom_to_conserved_split_rho
      uawsom_to_conserved         => uawsom_to_conserved_split_rho
    else if(uawsom_internal_e) then
      phys_to_primitive        => uawsom_to_primitive_inte
      uawsom_to_primitive         => uawsom_to_primitive_inte
      phys_to_conserved        => uawsom_to_conserved_inte
      uawsom_to_conserved         => uawsom_to_conserved_inte
    else if(uawsom_hydrodynamic_e) then
      phys_to_primitive        => uawsom_to_primitive_hde
      uawsom_to_primitive         => uawsom_to_primitive_hde
      phys_to_conserved        => uawsom_to_conserved_hde
      uawsom_to_conserved         => uawsom_to_conserved_hde
    else
      phys_to_primitive        => uawsom_to_primitive_origin
      uawsom_to_primitive         => uawsom_to_primitive_origin
      phys_to_conserved        => uawsom_to_conserved_origin
      uawsom_to_conserved         => uawsom_to_conserved_origin
    end if
      if(B0field.or.has_equi_rho0.or.has_equi_pe0) then
        phys_get_flux            => uawsom_get_flux_split
      else if(uawsom_hydrodynamic_e) then
        phys_get_flux            => uawsom_get_flux_split
      else
        phys_get_flux            => uawsom_get_flux
      end if
    if(uawsom_boris_simplification) then
      phys_get_v                 => uawsom_get_v_boris
      uawsom_get_v                  => uawsom_get_v_boris
      uawsom_kin_en                 => uawsom_kin_en_boris
    else
      phys_get_v                 => uawsom_get_v_origin
      uawsom_get_v                  => uawsom_get_v_origin
      uawsom_kin_en                 => uawsom_kin_en_origin
    end if
    if(B0field.or.has_equi_rho0) then
      phys_add_source_geom     => uawsom_add_source_geom_split
    else
      phys_add_source_geom     => uawsom_add_source_geom
    end if
    phys_add_source          => uawsom_add_source
    phys_check_params        => uawsom_check_params
    phys_write_info          => uawsom_write_info
    phys_angmomfix           => uawsom_angmomfix
    phys_handle_small_values => uawsom_handle_small_values_origin
    uawsom_handle_small_values  => uawsom_handle_small_values_origin
    phys_check_w             => uawsom_check_w_origin
    phys_get_pthermal        => uawsom_get_pthermal_origin
    uawsom_get_pthermal         => uawsom_get_pthermal_origin
    if(number_equi_vars>0) then
      phys_set_equi_vars => set_equi_vars_grid
    endif

    if(type_divb==divb_glm) then
      phys_modify_wLR => uawsom_modify_wLR
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_get_ct_velocity => uawsom_get_ct_velocity
      phys_update_faces => uawsom_update_faces
      phys_face_to_center => uawsom_face_to_center
      phys_modify_wLR => uawsom_modify_wLR
    else if(ndim>1) then
      phys_boundary_adjust => uawsom_boundary_adjust
    end if

    

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call uawsom_physical_units()

    if(.not. uawsom_energy .and. uawsom_thermal_conduction) then
      call mpistop("thermal conduction needs uawsom_energy=T")
    end if
    if(.not. uawsom_energy .and. uawsom_radiative_cooling) then
      call mpistop("radiative cooling needs uawsom_energy=T")
    end if

    ! resistive uawsom needs diagonal ghost cells
    if(uawsom_eta/=0.d0) phys_req_diagonal = .true.

    ! initialize thermal conduction module
    if (uawsom_thermal_conduction) then
      phys_req_diagonal = .true.

      call sts_init()
      call tc_init_params(uawsom_gamma)

      allocate(tc_fl)
      call tc_get_mhd_params(tc_fl,tc_params_read_uawsom)
      call add_sts_method(uawsom_get_tc_dt_mhd,uawsom_sts_set_source_tc_mhd,e_,&
         1,e_,1,.false.)
      if(phys_internal_e) then
        if(has_equi_pe0 .and. has_equi_rho0) then
          tc_fl%get_temperature_from_conserved => &
             uawsom_get_temperature_from_eint_with_equi
        else
          tc_fl%get_temperature_from_conserved => &
             uawsom_get_temperature_from_eint
        end if
      else if(uawsom_hydrodynamic_e) then
        tc_fl%get_temperature_from_conserved => &
           uawsom_get_temperature_from_hde
      else
        if(has_equi_pe0 .and. has_equi_rho0) then
          tc_fl%get_temperature_from_conserved => &
             uawsom_get_temperature_from_etot_with_equi
        else
          tc_fl%get_temperature_from_conserved => &
             uawsom_get_temperature_from_etot
        end if
      end if
      if(uawsom_solve_eaux) then
        call set_conversion_methods_to_head(uawsom_e_to_ei_aux,&
            uawsom_ei_to_e_aux)
      else if(.not. uawsom_internal_e) then
        call set_conversion_methods_to_head(uawsom_e_to_ei, uawsom_ei_to_e)
      end if
      if(has_equi_pe0 .and. has_equi_rho0) then
        tc_fl%get_temperature_from_eint => &
           uawsom_get_temperature_from_eint_with_equi
        if(uawsom_equi_thermal) then
          tc_fl%has_equi = .true.
          tc_fl%get_temperature_equi => uawsom_get_temperature_equi
          tc_fl%get_rho_equi => uawsom_get_rho_equi
        else
          tc_fl%has_equi = .false.
        endif
      else
        tc_fl%get_temperature_from_eint => uawsom_get_temperature_from_eint
      endif
      call set_error_handling_to_head(uawsom_tc_handle_small_e)
      tc_fl%get_rho => uawsom_get_rho
      tc_fl%e_ = e_
      tc_fl%Tcoff_ = Tcoff_
    end if

    ! Initialize radiative cooling module
    if (uawsom_radiative_cooling) then
      call radiative_cooling_init_params(uawsom_gamma,He_abundance)
      allocate(rc_fl)
      call radiative_cooling_init(rc_fl,rc_params_read)
      rc_fl%get_rho => uawsom_get_rho
      rc_fl%get_pthermal => uawsom_get_pthermal
      rc_fl%e_ = e_
      rc_fl%eaux_ = eaux_
      rc_fl%Tcoff_ = Tcoff_
      if(associated(usr_Rfactor)) then
        rc_fl%get_var_Rfactor => usr_Rfactor
      endif
      if(has_equi_pe0 .and. has_equi_rho0 .and. uawsom_equi_thermal) then
        rc_fl%has_equi = .true.
        rc_fl%get_rho_equi => uawsom_get_rho_equi
        rc_fl%get_pthermal_equi => uawsom_get_pe_equi
      else
        rc_fl%has_equi = .false.
      end if
    end if
    allocate(te_fl_uawsom)
    te_fl_uawsom%get_rho=> uawsom_get_rho
    te_fl_uawsom%get_pthermal=> uawsom_get_pthermal
    te_fl_uawsom%Rfactor = RR
    if(associated(usr_Rfactor)) then
      te_fl_uawsom%get_var_Rfactor => usr_Rfactor
    endif

    ! Initialize viscosity module
    if (uawsom_viscosity) call viscosity_init(phys_wider_stencil,&
       phys_req_diagonal)

    ! Initialize gravity module
    if(uawsom_gravity) then
      call gravity_init()
    end if

    ! Initialize particles module
    if(uawsom_particles) then
      call particles_init()
      if (particles_eta  < zero) particles_eta = uawsom_eta
      if (particles_etah < zero) particles_eta = uawsom_etah
      phys_req_diagonal = .true.
      if(mype==0) then
         write(*,*)&
             '*****Using particles:        with uawsom_eta, uawsom_etah :',&
             uawsom_eta, uawsom_etah
         write(*,*) '*****Using particles: particles_eta, particles_etah :',&
             particles_eta, particles_etah
      end if
    end if

    ! initialize magnetofriction module
    if(uawsom_magnetofriction) then
      phys_req_diagonal = .true.
      call magnetofriction_init()
    end if

    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in uawsom_get_flux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if(uawsom_hall) then
      phys_req_diagonal = .true.
      if(uawsom_4th_order) then
        phys_wider_stencil = 2
      else
        phys_wider_stencil = 1
      end if
    end if

    if(uawsom_ambipolar) then
      phys_req_diagonal = .true.
      if(uawsom_ambipolar_sts) then
        call sts_init()
        if(uawsom_internal_e) then
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,mag(1),&
             ndir,mag(1),ndir,.true.)
        else
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,&
             mom(ndir)+1,mag(ndir)-mom(ndir),mag(1),ndir,.true.)
        end if
      else
        uawsom_ambipolar_exp=.true.
        ! For flux ambipolar term, we need one more reconstructed layer since currents are computed
        ! in uawsom_get_flux: assuming one additional ghost layer (two for FOURTHORDER) was
        ! added in nghostcells.
        if(uawsom_4th_order) then
          phys_wider_stencil = 2
        else
          phys_wider_stencil = 1
        end if
      end if
    end if

    ! Initialize CAK radiation force module
    if (uawsom_cak_force) call cak_init(uawsom_gamma)

  end subroutine uawsom_phys_init



!!start th cond
  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  uawsom_sts_set_source_tc_mhd(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,&
     wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_mhd
    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, igrid, nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,1:nw),&
        w(ixImin1:ixImax1,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_mhd(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,wres,&
       fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine uawsom_sts_set_source_tc_mhd

  function uawsom_get_tc_dt_mhd(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dx1,&
     x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_mhd
 
    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: dx1, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_mhd(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dx1,x,tc_fl)
  end function uawsom_get_tc_dt_mhd

  subroutine uawsom_tc_handle_small_e(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      step)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    integer, intent(in)    :: step
    character(len=140) :: error_msg

    write(error_msg,"(a,i3)") "Thermal conduction step ", step
    call uawsom_handle_small_ei(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,e_,&
       error_msg)
  end subroutine uawsom_tc_handle_small_e

  ! fill in tc_fluid fields from namelist
  subroutine tc_params_read_uawsom(fl)
    use mod_global_parameters, only: unitpar,par_files
    type(tc_fluid), intent(inout) :: fl

    integer                      :: n

    ! list parameters
    logical :: tc_perpendicular=.true.
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0
    double precision :: tc_k_perp=0d0
    character(len=std_len)  :: tc_slope_limiter="MC"

    namelist /tc_list/ tc_perpendicular, tc_saturate, tc_slope_limiter,&
        tc_k_para, tc_k_perp

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, tc_list, end=111)
111     close(unitpar)
    end do

    fl%tc_perpendicular = tc_perpendicular
    fl%tc_saturate = tc_saturate
    fl%tc_k_para = tc_k_para
    fl%tc_k_perp = tc_k_perp
    fl%tc_slope_limiter = tc_slope_limiter
  end subroutine tc_params_read_uawsom
!!end th cond

!!rad cool
  subroutine rc_params_read(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_constants, only: bigdouble
    type(rc_fluid), intent(inout) :: fl
    integer                      :: n
    ! list parameters
    integer :: ncool = 4000
    double precision :: cfrac=0.1d0
  
    !> Name of cooling curve
    character(len=std_len)  :: coolcurve='JCcorona'
  
    !> Name of cooling method
    character(len=std_len)  :: coolmethod='exact'
  
    !> Fixed temperature not lower than tlow
    logical    :: Tfix=.false.
  
    !> Lower limit of temperature
    double precision   :: tlow=bigdouble
  
    !> Add cooling source in a split way (.true.) or un-split way (.false.)
    logical    :: rc_split=.false.


    namelist /rc_list/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix,&
        rc_split

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, rc_list, end=111)
111     close(unitpar)
    end do

    fl%ncool=ncool
    fl%coolcurve=coolcurve
    fl%coolmethod=coolmethod
    fl%tlow=tlow
    fl%Tfix=Tfix
    fl%rc_split=rc_split
    fl%cfrac=cfrac

  end subroutine rc_params_read
!! end rad cool

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid_faces(igrid,x,ixImin1,ixImax1,ixOmin1,ixOmax1)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in) :: igrid, ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)

    double precision :: delx(ixImin1:ixImax1,1:ndim)
    double precision :: xC(ixImin1:ixImax1,1:ndim),xshift1
    integer :: idims, ixCmin1,ixCmax1, hxOmin1,hxOmax1, ix, idims2

    if(slab_uniform)then
      delx(ixImin1:ixImax1,1)=rnode(rpdx1_,igrid)
    else
      ! for all non-cartesian and stretched cartesian coordinates
      delx(ixImin1:ixImax1,1:ndim)=ps(igrid)%dx(ixImin1:ixImax1,1:ndim)
    endif

    do idims=1,ndim
      hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
      if(stagger_grid) then
        ! ct needs all transverse cells
        ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
        ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1);
      else
        ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
        ixCmax1=ixOmax1; ixCmin1=hxOmin1;
      end if
      ! always xshift=0 or 1/2
      xshift1=half*(one-kr(1,idims));
      do idims2=1,ndim
        select case(idims2)
        case(1)
          do ix = ixCmin1,ixCmax1
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ix,1)=x(ix,1)+(half-xshift1)*delx(ix,1)
          end do
        end select
      end do
      call usr_set_equi_vars(ixImin1,ixImax1,ixCmin1,ixCmax1,xC,&
         ps(igrid)%equi_vars(ixImin1:ixImax1,1:number_equi_vars,idims))
    end do

  end subroutine set_equi_vars_grid_faces

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid(igrid)
    use mod_global_parameters
    use mod_usr_methods
  
    integer, intent(in) :: igrid

    !values at the center
    call usr_set_equi_vars(ixGlo1,ixGhi1,ixGlo1,ixGhi1,ps(igrid)%x,&
       ps(igrid)%equi_vars(ixGlo1:ixGhi1,1:number_equi_vars,0))

    !values at the interfaces
    call set_equi_vars_grid_faces(igrid,ps(igrid)%x,ixGlo1,ixGhi1,ixMlo1,&
       ixMhi1)
 
  end subroutine set_equi_vars_grid

  ! w, wnew conserved, add splitted variables back to wnew
  function convert_vars_splitting(ixImin1,ixImax1,ixOmin1,ixOmax1, w, x,&
      nwc) result(wnew)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1, nwc
    double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision   :: wnew(ixOmin1:ixOmax1, 1:nwc)
    double precision   :: rho(ixImin1:ixImax1)

    call  uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,&
       rho(ixImin1:ixImax1))
    wnew(ixOmin1:ixOmax1,rho_) = rho(ixOmin1:ixOmax1)
    wnew(ixOmin1:ixOmax1,mom(:)) =  w(ixOmin1:ixOmax1,mom(:))

    if (B0field) then
      ! add background magnetic field B0 to B
      wnew(ixOmin1:ixOmax1,mag(:))=w(ixOmin1:ixOmax1,&
         mag(:))+block%B0(ixOmin1:ixOmax1,:,0)
    else
      wnew(ixOmin1:ixOmax1,mag(:))=w(ixOmin1:ixOmax1,mag(:))
    end if

    if(uawsom_energy) then
      wnew(ixOmin1:ixOmax1,e_) = w(ixOmin1:ixOmax1,e_)
      if(has_equi_pe0) then
        wnew(ixOmin1:ixOmax1,e_) = wnew(ixOmin1:ixOmax1,&
           e_) + block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,0)* inv_gamma_1
      end if
      if(B0field .and. .not. uawsom_internal_e) then
          wnew(ixOmin1:ixOmax1,e_)=wnew(ixOmin1:ixOmax1,&
             e_)+0.5d0*sum(block%B0(ixOmin1:ixOmax1,:,0)**2,&
             dim=ndim+1) + sum(w(ixOmin1:ixOmax1,&
             mag(:))*block%B0(ixOmin1:ixOmax1,:,0),dim=ndim+1)
      end if
    end if

    !wnew(ixO^S,wminus_) = w(ixO^S,wminus_)
    !wnew(ixO^S,wplus_) = w(ixO^S,wplus_)

    wnew(ixOmin1:ixOmax1,wkminus_) = w(ixOmin1:ixOmax1,wkminus_)
    wnew(ixOmin1:ixOmax1,wkplus_) = w(ixOmin1:ixOmax1,wkplus_)
    wnew(ixOmin1:ixOmax1,wAminus_) = w(ixOmin1:ixOmax1,wAminus_)
    wnew(ixOmin1:ixOmax1,wAplus_) = w(ixOmin1:ixOmax1,wAplus_)


  end function convert_vars_splitting

  subroutine uawsom_check_params
    use mod_global_parameters
    use mod_usr_methods
    use mod_convert, only: add_convert_method

    ! after user parameter setting
    gamma_1=uawsom_gamma-1.d0
    if (.not. uawsom_energy) then
       if (uawsom_gamma <= 0.0d0) call mpistop ("Error: uawsom_gamma <= 0")
       if (uawsom_adiab < 0.0d0) call mpistop ("Error: uawsom_adiab < 0")
       small_pressure = uawsom_adiab*small_density**uawsom_gamma
    else
       if (uawsom_gamma <= 0.0d0 .or. uawsom_gamma == 1.0d0) call mpistop &
          ("Error: uawsom_gamma <= 0 or uawsom_gamma == 1")
       inv_gamma_1=1.d0/gamma_1
       small_e = small_pressure * inv_gamma_1
    end if

    if (number_equi_vars > 0 .and. .not. associated(usr_set_equi_vars)) then
      call mpistop("usr_set_equi_vars has to be implemented in the user file")
    endif
    if(convert .or. autoconvert) then
      if(convert_type .eq. 'dat_generic_mpi') then
        if(uawsom_dump_full_vars) then
          if(mype .eq. 0) print*, " add conversion method: split -> full "
          call add_convert_method(convert_vars_splitting, nw, cons_wnames,&
              "new")
        endif
      endif
    endif
  end subroutine uawsom_check_params

  subroutine uawsom_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0,c_lightspeed
    double precision :: a,b
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
      c_lightspeed=c_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi ! G2 cm2 dyne-1
      c_lightspeed=const_c
    end if
    if(eq_state_units) then
      a = 1d0 + 4d0 * He_abundance
      b = 1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + 1d0)+1d0)
      RR = 1d0
    else
      a = 1d0
      b = 1d0
      RR = (1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + &
         1d0)+1d0))/(1d0 + 4d0 * He_abundance)
    endif
    if(unit_density/=1.d0) then
      unit_numberdensity=unit_density/(a*mp)
    else
      ! unit of numberdensity is independent by default
      unit_density=a*mp*unit_numberdensity
    end if
    if(unit_velocity/=1.d0) then
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    else if(unit_pressure/=1.d0) then
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    else if(unit_magneticfield/=1.d0) then
      unit_pressure=unit_magneticfield**2/miu0
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
    else
      ! unit of temperature is independent by default
      unit_pressure=b*unit_numberdensity*kB*unit_temperature
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

  end subroutine uawsom_physical_units

  subroutine uawsom_check_w_origin(primitive,ixImin1,ixImax1,ixOmin1,ixOmax1,w,&
     flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,nw)
    double precision :: tmp(ixImin1:ixImax1)
    logical, intent(inout) :: flag(ixImin1:ixImax1,1:nw)

    flag=.false.
    if(has_equi_rho0) then
      tmp(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,0)
    else
      tmp(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,rho_)
    endif
    where(tmp(ixOmin1:ixOmax1) < small_density) flag(ixOmin1:ixOmax1,&
       rho_) = .true.

    if(uawsom_energy) then
      if(primitive) then
        tmp(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,e_)
        if(has_equi_pe0) then
          tmp(ixOmin1:ixOmax1) = tmp(ixOmin1:ixOmax1)+&
             block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,0)
        endif
        where(tmp(ixOmin1:ixOmax1) < small_pressure) flag(ixOmin1:ixOmax1,&
           e_) = .true.
      else
        if(uawsom_internal_e) then
          tmp(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,e_)
          if(has_equi_pe0) then
            tmp(ixOmin1:ixOmax1) = tmp(ixOmin1:ixOmax1)+&
               block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,0)*inv_gamma_1
          endif
          where(tmp(ixOmin1:ixOmax1) < small_e) flag(ixOmin1:ixOmax1,&
             e_) = .true.
        else
          tmp(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,e_)-uawsom_kin_en(w,ixImin1,&
             ixImax1,ixOmin1,ixOmax1)-uawsom_mag_en(w,ixImin1,ixImax1,ixOmin1,&
             ixOmax1)-uawsom_wk_en(w,ixImin1,ixImax1,ixOmin1,&
             ixOmax1)-uawsom_wA_en(w,ixImin1,ixImax1,ixOmin1,ixOmax1) !Max: Added AW energy
          if(has_equi_pe0) then
            tmp(ixOmin1:ixOmax1) = tmp(ixOmin1:ixOmax1)+&
               block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,0)*inv_gamma_1
          endif
          where(tmp(ixOmin1:ixOmax1) < small_e) flag(ixOmin1:ixOmax1,&
             e_) = .true.
        end if
      end if
    end if

  end subroutine uawsom_check_w_origin

  subroutine uawsom_check_w_hde(primitive,ixImin1,ixImax1,ixOmin1,ixOmax1,w,&
     flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,nw)
    double precision :: tmp(ixImin1:ixImax1)
    logical, intent(inout) :: flag(ixImin1:ixImax1,1:nw)

    flag=.false.
    where(w(ixOmin1:ixOmax1,rho_) < small_density) flag(ixOmin1:ixOmax1,&
       rho_) = .true.

    if(uawsom_energy) then
      if(primitive) then
        where(w(ixOmin1:ixOmax1,e_) < small_pressure) flag(ixOmin1:ixOmax1,&
           e_) = .true.
      else
        tmp(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,e_)-uawsom_kin_en(w,ixImin1,&
           ixImax1,ixOmin1,ixOmax1)
        where(tmp(ixOmin1:ixOmax1) < small_e) flag(ixOmin1:ixOmax1,&
           e_) = .true.
      end if
    end if

  end subroutine uawsom_check_w_hde

  !> Transform primitive variables into conservative ones
  subroutine uawsom_to_conserved_origin(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    double precision :: inv_gamma2(ixOmin1:ixOmax1)
    integer                         :: idir

    !if (fix_small_values) then
    !  call uawsom_handle_small_values(.true., w, x, ixI^L, ixO^L, 'uawsom_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy and waves
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,&
         p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,mom(:))**2,&
         dim=ndim+1)*w(ixOmin1:ixOmax1,rho_)+uawsom_mag_en(w, ixImin1,ixImax1,&
          ixOmin1,ixOmax1)+uawsom_wk_en(w, ixImin1,ixImax1, ixOmin1,&
         ixOmax1)+uawsom_wA_en(w, ixImin1,ixImax1, ixOmin1,ixOmax1) !Max: Added AW energy 
      if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,eaux_)=w(ixOmin1:ixOmax1,&
         paux_)*inv_gamma_1
    end if
 
    if(uawsom_boris_simplification) then
      inv_gamma2=1.d0+sum(w(ixOmin1:ixOmax1,mag(:))**2,&
         dim=ndim+1)/w(ixOmin1:ixOmax1,rho_)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1, mom(idir)) = inv_gamma2*w(ixOmin1:ixOmax1,&
           rho_)*w(ixOmin1:ixOmax1, mom(idir))
      end do
    else
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1,&
           rho_)*w(ixOmin1:ixOmax1, mom(idir))
      end do
    end if
  end subroutine uawsom_to_conserved_origin

  !> Transform primitive variables into conservative ones
  subroutine uawsom_to_conserved_hde(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    integer                         :: idir

    !if (fix_small_values) then
    !  call uawsom_handle_small_values(.true., w, x, ixI^L, ixO^L, 'uawsom_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy 
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,&
         p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,mom(:))**2,&
         dim=ndim+1)*w(ixOmin1:ixOmax1,rho_) !Max: hde so we do not consider contribution from magnetic waves AW and kink
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
      w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1,&
         rho_)*w(ixOmin1:ixOmax1, mom(idir))
    end do
  end subroutine uawsom_to_conserved_hde

  !> Transform primitive variables into conservative ones
  subroutine uawsom_to_conserved_inte(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    double precision :: inv_gamma2(ixOmin1:ixOmax1)
    integer                         :: idir

    !if (fix_small_values) then
    !  call uawsom_handle_small_values(.true., w, x, ixI^L, ixO^L, 'uawsom_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy ! Max: No kinetic energy or magnetic energy here?
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,p_)*inv_gamma_1
    end if

    if(uawsom_boris_simplification) then
      inv_gamma2=1.d0+sum(w(ixOmin1:ixOmax1,mag(:))**2,&
         dim=ndim+1)/w(ixOmin1:ixOmax1,rho_)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1, mom(idir)) = inv_gamma2*w(ixOmin1:ixOmax1,&
           rho_)*w(ixOmin1:ixOmax1, mom(idir))
      end do
    else
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1,&
           rho_)*w(ixOmin1:ixOmax1, mom(idir))
      end do
    end if
  end subroutine uawsom_to_conserved_inte

  !> Transform primitive variables into conservative ones
  subroutine uawsom_to_conserved_split_rho(ixImin1,ixImax1,ixOmin1,ixOmax1,w,&
     x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
    integer                         :: idir
    double precision                :: rho(ixImin1:ixImax1)

    !if (fix_small_values) then
    !  call uawsom_handle_small_values(.true., w, x, ixI^L, ixO^L, 'uawsom_to_conserved')
    !end if

    rho(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
       rho_) + block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i)
    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(uawsom_energy) then
      if(uawsom_internal_e) then
        w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,p_)*inv_gamma_1
      else
        w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,&
           p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,mom(:))**2,&
           dim=ndim+1)*rho(ixOmin1:ixOmax1)+uawsom_mag_en(w, ixImin1,ixImax1,&
            ixOmin1,ixOmax1)+uawsom_wk_en(w,ixImin1,ixImax1, ixOmin1,&
           ixOmax1)+uawsom_wA_en(w,ixImin1,ixImax1, ixOmin1,ixOmax1) !Max: Added AW energy
        if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,eaux_)=w(ixOmin1:ixOmax1,&
           paux_)*inv_gamma_1
      end if
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixOmin1:ixOmax1, mom(idir)) = rho(ixOmin1:ixOmax1) * &
          w(ixOmin1:ixOmax1, mom(idir))
    end do
  end subroutine uawsom_to_conserved_split_rho

  !> Transform conservative variables into primitive ones
  subroutine uawsom_to_primitive_origin(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1),&
        gamma2(ixOmin1:ixOmax1)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call uawsom_handle_small_values(.false., w, x, ixImin1,ixImax1, ixOmin1,&
         ixOmax1, 'uawsom_to_primitive_origin')
    end if

    inv_rho(ixOmin1:ixOmax1) = 1d0/w(ixOmin1:ixOmax1,rho_)

    ! Calculate pressure = (gamma-1) * (e-ek-eb) ! Max: In Norberts work it was (gamma-1) * (e-ek-eb - E_kink) since he included kink energy, 
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,p_)=gamma_1*(w(ixOmin1:ixOmax1,e_)-uawsom_kin_en(w,&
         ixImin1,ixImax1,ixOmin1,ixOmax1,inv_rho)-uawsom_mag_en(w,ixImin1,&
         ixImax1,ixOmin1,ixOmax1)-uawsom_wk_en(w,ixImin1,ixImax1,ixOmin1,&
         ixOmax1)-uawsom_wA_en(w,ixImin1,ixImax1,ixOmin1,ixOmax1)) !Max: AW energy added
      if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,paux_)=w(ixOmin1:ixOmax1,&
         eaux_)*gamma_1
    end if

    if(uawsom_boris_simplification) then
      gamma2=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,mag(:))**2,&
         dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1,&
             mom(idir))*inv_rho*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1, mom(idir))*inv_rho
      end do
    end if

  end subroutine uawsom_to_primitive_origin

  !> Transform conservative variables into primitive ones
  subroutine uawsom_to_primitive_hde(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call uawsom_handle_small_values(.false., w, x, ixImin1,ixImax1, ixOmin1,&
         ixOmax1, 'uawsom_to_primitive_hde')
    end if

    inv_rho(ixOmin1:ixOmax1) = 1d0/w(ixOmin1:ixOmax1,rho_)

    ! Calculate pressure = (gamma-1) * (e-ek)
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,p_)=gamma_1*(w(ixOmin1:ixOmax1,e_)-uawsom_kin_en(w,&
         ixImin1,ixImax1,ixOmin1,ixOmax1,inv_rho))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1, mom(idir))*inv_rho
    end do

  end subroutine uawsom_to_primitive_hde

  !> Transform conservative variables into primitive ones
  subroutine uawsom_to_primitive_inte(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1),&
        gamma2(ixOmin1:ixOmax1)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call uawsom_handle_small_values(.false., w, x, ixImin1,ixImax1, ixOmin1,&
         ixOmax1, 'uawsom_to_primitive_inte')
    end if

    inv_rho(ixOmin1:ixOmax1) = 1d0/w(ixOmin1:ixOmax1,rho_)

    ! Calculate pressure = (gamma-1) * e_internal
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,p_)=w(ixOmin1:ixOmax1,e_)*gamma_1
    end if

    if(uawsom_boris_simplification) then
      gamma2=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,mag(:))**2,&
         dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1,&
             mom(idir))*inv_rho*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1, mom(idir))*inv_rho
      end do
    end if

  end subroutine uawsom_to_primitive_inte

  !> Transform conservative variables into primitive ones
  subroutine uawsom_to_primitive_split_rho(ixImin1,ixImax1,ixOmin1,ixOmax1,w,&
     x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
    double precision                :: inv_rho(ixOmin1:ixOmax1)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call uawsom_handle_small_values(.false., w, x, ixImin1,ixImax1, ixOmin1,&
         ixOmax1, 'uawsom_to_primitive_split_rho')
    end if

    inv_rho(ixOmin1:ixOmax1) = 1d0/(w(ixOmin1:ixOmax1,&
       rho_) + block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i))

    ! Calculate pressure = (gamma-1) * (e-ek-eb)
    if(uawsom_energy) then
      if(uawsom_internal_e) then
        w(ixOmin1:ixOmax1,p_)=w(ixOmin1:ixOmax1,e_)*gamma_1
      else
        w(ixOmin1:ixOmax1,p_)=gamma_1*(w(ixOmin1:ixOmax1,e_)-uawsom_kin_en(w,&
           ixImin1,ixImax1,ixOmin1,ixOmax1,inv_rho)-uawsom_mag_en(w,ixImin1,&
           ixImax1,ixOmin1,ixOmax1)-uawsom_wk_en(w,ixImin1,ixImax1,ixOmin1,&
           ixOmax1)-uawsom_wA_en(w,ixImin1,ixImax1,ixOmin1,ixOmax1)) !Max: AW energy added
        if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,paux_)=w(ixOmin1:ixOmax1,&
           eaux_)*gamma_1
      end if
    end if
    
    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1, mom(idir))*inv_rho
    end do

  end subroutine uawsom_to_primitive_split_rho

  !> Transform internal energy to total energy
  subroutine uawsom_ei_to_e(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    ! Calculate total energy from internal, kinetic and magnetic energy
    w(ixImin1:ixImax1,e_)=w(ixImin1:ixImax1,e_)+uawsom_kin_en(w,ixImin1,&
       ixImax1,ixImin1,ixImax1)+uawsom_mag_en(w,ixImin1,ixImax1,ixImin1,&
       ixImax1)+uawsom_wk_en(w,ixImin1,ixImax1,ixImin1,ixImax1)+uawsom_wA_en(w,&
       ixImin1,ixImax1,ixImin1,ixImax1) !Max: Added AW energy density

  end subroutine uawsom_ei_to_e

  !> Transform total energy to internal energy
  subroutine uawsom_e_to_ei(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    ! Calculate ei = e - ek - eb
    w(ixImin1:ixImax1,e_)=w(ixImin1:ixImax1,e_)-uawsom_kin_en(w,ixImin1,&
       ixImax1,ixImin1,ixImax1)-uawsom_mag_en(w,ixImin1,ixImax1,ixImin1,&
       ixImax1)-uawsom_wk_en(w,ixImin1,ixImax1,ixImin1,ixImax1)-uawsom_wA_en(w,&
       ixImin1,ixImax1,ixImin1,ixImax1) !Max: Added AW energy

    if(fix_small_values) then
      call uawsom_handle_small_ei(w,x,ixImin1,ixImax1,ixImin1,ixImax1,e_,&
         'uawsom_e_to_ei')
    end if

  end subroutine uawsom_e_to_ei

  !> Transform hydrodynamic energy to internal energy
  subroutine uawsom_e_to_ei_hde(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    ! Calculate ei = e - ek
    w(ixImin1:ixImax1,e_)=w(ixImin1:ixImax1,e_)-uawsom_kin_en(w,ixImin1,&
       ixImax1,ixImin1,ixImax1)

    if(fix_small_values) then
      call uawsom_handle_small_ei(w,x,ixImin1,ixImax1,ixImin1,ixImax1,e_,&
         'uawsom_e_to_ei_hde')
    end if

  end subroutine uawsom_e_to_ei_hde

  !> Update eaux and transform internal energy to total energy
  subroutine uawsom_ei_to_e_aux(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
 
    w(ixImin1:ixImax1,eaux_)=w(ixImin1:ixImax1,e_)
    ! Calculate total energy from internal, kinetic and magnetic energy ! Max: Do we need Kink and AW energy here?
    w(ixImin1:ixImax1,e_)=w(ixImin1:ixImax1,e_)+uawsom_kin_en(w,ixImin1,&
       ixImax1,ixImin1,ixImax1)+uawsom_mag_en(w,ixImin1,ixImax1,ixImin1,&
       ixImax1)

  end subroutine uawsom_ei_to_e_aux

  !> Transform total energy to internal energy via eaux as internal energy
  subroutine uawsom_e_to_ei_aux(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)

    w(ixImin1:ixImax1,e_)=w(ixImin1:ixImax1,eaux_)

  end subroutine uawsom_e_to_ei_aux

  subroutine uawsom_handle_small_values_origin(primitive, w, x, ixImin1,&
     ixImax1, ixOmin1,ixOmax1, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixImin1:ixImax1,1:nw)
    double precision :: tmp2(ixImin1:ixImax1)

    if(small_values_method == "ignore") return

    call phys_check_w(primitive, ixImin1,ixImax1, ixOmin1,ixOmax1, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_rho0) then
          where(flag(ixOmin1:ixOmax1,rho_)) w(ixOmin1:ixOmax1,&
             rho_) = small_density-block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,&
             0)
        else
          where(flag(ixOmin1:ixOmax1,rho_)) w(ixOmin1:ixOmax1,&
             rho_) = small_density
        endif
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixOmin1:ixOmax1,rho_)) w(ixOmin1:ixOmax1,&
                mom(idir)) = 0.0d0
          end if
        end do

        if(uawsom_energy) then
          if(primitive) then
           if(has_equi_pe0) then
            tmp2(ixOmin1:ixOmax1) = small_pressure - &
               block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,0)
           else
            tmp2(ixOmin1:ixOmax1) = small_pressure
           endif
           where(flag(ixOmin1:ixOmax1,e_)) w(ixOmin1:ixOmax1,&
              p_) = tmp2(ixOmin1:ixOmax1)
          else
            ! conserved
            if(has_equi_pe0) then
              tmp2(ixOmin1:ixOmax1) = small_e - &
                 block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,0)*inv_gamma_1
            else
              tmp2(ixOmin1:ixOmax1) = small_e
            endif
            if(uawsom_internal_e) then
              where(flag(ixOmin1:ixOmax1,e_))
                w(ixOmin1:ixOmax1,e_)=tmp2(ixOmin1:ixOmax1)
              end where
            else
              where(flag(ixOmin1:ixOmax1,e_))
                w(ixOmin1:ixOmax1,e_) = tmp2(ixOmin1:ixOmax1)+uawsom_kin_en(w,&
                   ixImin1,ixImax1,ixOmin1,ixOmax1)+uawsom_mag_en(w,ixImin1,&
                   ixImax1,ixOmin1,ixOmax1)+uawsom_wk_en(w,ixImin1,ixImax1,&
                   ixOmin1,ixOmax1)+uawsom_wA_en(w,ixImin1,ixImax1,ixOmin1,&
                   ixOmax1) !Max: Added AW energy
             end where
              if(uawsom_solve_eaux) then
                where(flag(ixOmin1:ixOmax1,e_))
                  w(ixOmin1:ixOmax1,eaux_)=tmp2(ixOmin1:ixOmax1)
                end where
              end if
            end if
          end if
        end if
      case ("average")
        if(primitive)then
          call small_values_average(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x,&
              flag, rho_)
          if(uawsom_energy) then
            call small_values_average(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x,&
                flag, p_)
          end if
        else
          ! do averaging of density
          call small_values_average(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x,&
              flag, rho_)
          if(uawsom_energy) then
             ! do averaging of internal energy
            if(.not.uawsom_internal_e) then
              w(ixImin1:ixImax1,e_)=w(ixImin1:ixImax1,e_)-uawsom_kin_en(w,&
                 ixImin1,ixImax1,ixImin1,ixImax1)-uawsom_mag_en(w,ixImin1,&
                 ixImax1,ixImin1,ixImax1)-uawsom_wk_en(w,ixImin1,ixImax1,&
                 ixImin1,ixImax1)-uawsom_wA_en(w,ixImin1,ixImax1,ixImin1,&
                 ixImax1) !Max: Added AW energy 
            end if
            call small_values_average(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x,&
                flag, e_)
             ! convert back
            if(.not.uawsom_internal_e) then
              w(ixImin1:ixImax1,e_)=w(ixImin1:ixImax1,e_)+uawsom_kin_en(w,&
                 ixImin1,ixImax1,ixImin1,ixImax1)+uawsom_mag_en(w,ixImin1,&
                 ixImax1,ixImin1,ixImax1)+uawsom_wk_en(w,ixImin1,ixImax1,&
                 ixImin1,ixImax1)+uawsom_wA_en(w,ixImin1,ixImax1,ixImin1,&
                 ixImax1) !Max: Added AW energy
            end if
            ! eaux
            if(uawsom_solve_eaux) then
              call small_values_average(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x,&
                  flag, paux_)
            end if
          end if
        endif
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek-eb)
          if(uawsom_energy) then
            if(uawsom_internal_e) then
              w(ixOmin1:ixOmax1,p_)=w(ixOmin1:ixOmax1,e_)*gamma_1
            else
              w(ixOmin1:ixOmax1,p_)=gamma_1*(w(ixOmin1:ixOmax1,&
                 e_)-uawsom_kin_en(w,ixImin1,ixImax1,ixOmin1,&
                 ixOmax1)-uawsom_mag_en(w,ixImin1,ixImax1,ixOmin1,&
                 ixOmax1)-uawsom_wk_en(w,ixImin1,ixImax1,ixOmin1,&
                 ixOmax1)-uawsom_wA_en(w,ixImin1,ixImax1,ixOmin1,ixOmax1)) !Max: AW energy added
              if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,paux_)=w(ixOmin1:ixOmax1,&
                 eaux_)*gamma_1
            end if
          end if
          ! Convert momentum to velocity
          if(has_equi_rho0) then
            tmp2(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
               rho_) + block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,0)
          else
            tmp2(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,rho_)
          endif
          do idir = 1, ndir
             w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1,&
                 mom(idir))/tmp2(ixOmin1:ixOmax1)
          end do
        end if
        call small_values_error(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, flag,&
            subname)
      end select
    end if
  end subroutine uawsom_handle_small_values_origin

  !> Calculate v vector
  subroutine uawsom_get_v_origin(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1,nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ndir)

    double precision :: rho(ixImin1:ixImax1)
    integer :: idir

    call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,rho)

    rho(ixOmin1:ixOmax1)=1.d0/rho(ixOmin1:ixOmax1)
    ! Convert momentum to velocity
    do idir = 1, ndir
       v(ixOmin1:ixOmax1, idir) = w(ixOmin1:ixOmax1,&
           mom(idir))*rho(ixOmin1:ixOmax1)
    end do

  end subroutine uawsom_get_v_origin

  !> Calculate v vector
  subroutine uawsom_get_v_boris(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1,nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ndir)

    double precision              :: rho(ixImin1:ixImax1),&
        gamma2(ixOmin1:ixOmax1)
    integer :: idir

    call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,rho)

    rho(ixOmin1:ixOmax1)=1.d0/rho(ixOmin1:ixOmax1)
    gamma2=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,mag(:))**2,&
       dim=ndim+1)*rho(ixOmin1:ixOmax1)*inv_squared_c)
    ! Convert momentum to velocity
    do idir = 1, ndir
       v(ixOmin1:ixOmax1, idir) = w(ixOmin1:ixOmax1,&
           mom(idir))*rho(ixOmin1:ixOmax1)*gamma2
    end do

  end subroutine uawsom_get_v_boris

  !> Calculate v component
  subroutine uawsom_get_v_idim(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1)

    double precision              :: rho(ixImin1:ixImax1)

    call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,rho)

    if(uawsom_boris_simplification) then
      v(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
          mom(idim)) / rho(ixOmin1:ixOmax1) /(1.d0+sum(w(ixOmin1:ixOmax1,&
         mag(:))**2,dim=ndim+1)/rho(ixOmin1:ixOmax1)*inv_squared_c)
    else
      v(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
          mom(idim)) / rho(ixOmin1:ixOmax1)
    end if

  end subroutine uawsom_get_v_idim

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine uawsom_get_cmax_origin(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
     cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in) :: w(ixImin1:ixImax1, nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1)
    double precision :: vel(ixImin1:ixImax1)

    call uawsom_get_csound(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,cmax)
    call uawsom_get_v_idim(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,vel)

    cmax(ixOmin1:ixOmax1)=abs(vel(ixOmin1:ixOmax1))+cmax(ixOmin1:ixOmax1)

  end subroutine uawsom_get_cmax_origin

  subroutine uawsom_get_a2max(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,a2max)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1, nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(inout) :: a2max(ndim)
    double precision :: a2(ixImin1:ixImax1,ndim,nw)
    integer :: gxOmin1,gxOmax1,hxOmin1,hxOmax1,jxOmin1,jxOmax1,kxOmin1,kxOmax1,&
       i,j

    a2=zero
    do i = 1,ndim
      !> 4th order
      hxOmin1=ixOmin1-kr(i,1);hxOmax1=ixOmax1-kr(i,1);
      gxOmin1=hxOmin1-kr(i,1);gxOmax1=hxOmax1-kr(i,1);
      jxOmin1=ixOmin1+kr(i,1);jxOmax1=ixOmax1+kr(i,1);
      kxOmin1=jxOmin1+kr(i,1);kxOmax1=jxOmax1+kr(i,1);
      a2(ixOmin1:ixOmax1,i,1:nw)=abs(-w(kxOmin1:kxOmax1,&
         1:nw)+16.d0*w(jxOmin1:jxOmax1,1:nw)-30.d0*w(ixOmin1:ixOmax1,&
         1:nw)+16.d0*w(hxOmin1:hxOmax1,1:nw)-w(gxOmin1:gxOmax1,1:nw))
      a2max(i)=maxval(a2(ixOmin1:ixOmax1,i,1:nw))/12.d0/dxlevel(i)**2
    end do
  end subroutine uawsom_get_a2max

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine uawsom_get_tcutoff(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,Tco_local,&
     Tmax_local)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(out) :: Tco_local,Tmax_local

    double precision, parameter :: trac_delta=0.25d0
    double precision :: tmp1(ixImin1:ixImax1),Te(ixImin1:ixImax1),&
       lts(ixImin1:ixImax1)
    double precision, dimension(ixImin1:ixImax1,1:ndir) :: bunitvec
    double precision, dimension(ixImin1:ixImax1,1:ndim) :: gradT
    double precision :: Bdir(ndim)
    double precision :: ltrc,ltrp,altr(ixImin1:ixImax1)
    integer :: idims,jxOmin1,jxOmax1,hxOmin1,hxOmax1,ixA1,ixB1
    integer :: jxPmin1,jxPmax1,hxPmin1,hxPmax1,ixPmin1,ixPmax1
    logical :: lrlt(ixImin1:ixImax1)

    ! reuse lts as rhoc
    call uawsom_get_rho(w,x,ixImin1,ixImax1,ixImin1,ixImax1,lts)
    if(uawsom_internal_e) then
      tmp1(ixImin1:ixImax1)=w(ixImin1:ixImax1,e_)*gamma_1
    else
      call phys_get_pthermal(w,x,ixImin1,ixImax1,ixImin1,ixImax1,tmp1)
    end if
    Te(ixImin1:ixImax1)=tmp1(ixImin1:ixImax1)/lts(ixImin1:ixImax1)
    Tco_local=zero
    Tmax_local=maxval(Te(ixOmin1:ixOmax1))

    
    select case(uawsom_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      block%wextra(ixImin1:ixImax1,Tcoff_)=2.5d5/unit_temperature
    case(1)
      hxOmin1=ixOmin1-1;hxOmax1=ixOmax1-1;
      jxOmin1=ixOmin1+1;jxOmax1=ixOmax1+1;
      lts(ixOmin1:ixOmax1)=0.5d0*abs(Te(jxOmin1:jxOmax1)-&
         Te(hxOmin1:hxOmax1))/Te(ixOmin1:ixOmax1)
      lrlt=.false.
      where(lts(ixOmin1:ixOmax1) > trac_delta)
        lrlt(ixOmin1:ixOmax1)=.true.
      end where
      if(any(lrlt(ixOmin1:ixOmax1))) then
        Tco_local=maxval(Te(ixOmin1:ixOmax1), mask=lrlt(ixOmin1:ixOmax1))
      end if
    case(2)
      !> iijima et al. 2021, LTRAC method
      ltrc=1.5d0
      ltrp=4.d0
      ixPmin1=ixOmin1-1;ixPmax1=ixOmax1+1;
      hxOmin1=ixOmin1-1;hxOmax1=ixOmax1-1;
      jxOmin1=ixOmin1+1;jxOmax1=ixOmax1+1;
      hxPmin1=ixPmin1-1;hxPmax1=ixPmax1-1;
      jxPmin1=ixPmin1+1;jxPmax1=ixPmax1+1;
      lts(ixPmin1:ixPmax1)=0.5d0*abs(Te(jxPmin1:jxPmax1)-&
         Te(hxPmin1:hxPmax1))/Te(ixPmin1:ixPmax1)
      lts(ixPmin1:ixPmax1)=max(one, (exp(lts(ixPmin1:ixPmax1))/ltrc)**ltrp)
      lts(ixOmin1:ixOmax1)=0.25d0*(lts(jxOmin1:jxOmax1)+&
         two*lts(ixOmin1:ixOmax1)+lts(hxOmin1:hxOmax1))
      block%wextra(ixOmin1:ixOmax1,Tcoff_)=Te(ixOmin1:ixOmax1)*lts(&
         ixOmin1:ixOmax1)**0.4d0
    case default
      call mpistop("uawsom_trac_type not allowed for 1D simulation")
    end select
   
    
  end subroutine uawsom_get_tcutoff

  !> get H speed for H-correction to fix the carbuncle problem at grid-aligned shock front
  subroutine uawsom_get_H_speed(wprim,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
     Hspeed)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in)    :: wprim(ixImin1:ixImax1, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out)   :: Hspeed(ixImin1:ixImax1,&
       1:number_species)

    double precision :: csound(ixImin1:ixImax1,ndim),tmp(ixImin1:ixImax1)
    integer :: jxCmin1,jxCmax1, ixCmin1,ixCmax1, ixAmin1,ixAmax1, id, ix1

    Hspeed=0.d0
    ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
    do id=1,ndim
      call uawsom_get_csound_prim(wprim,x,ixImin1,ixImax1,ixAmin1,ixAmax1,id,&
         tmp)
      csound(ixAmin1:ixAmax1,id)=tmp(ixAmin1:ixAmax1)
    end do
    ixCmax1=ixOmax1;
    ixCmin1=ixOmin1+kr(idim,1)-1;
    jxCmax1=ixCmax1+kr(idim,1);
    jxCmin1=ixCmin1+kr(idim,1);
    Hspeed(ixCmin1:ixCmax1,1)=0.5d0*abs(wprim(jxCmin1:jxCmax1,&
       mom(idim))+csound(jxCmin1:jxCmax1,idim)-wprim(ixCmin1:ixCmax1,&
       mom(idim))+csound(ixCmin1:ixCmax1,idim))

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=ixCmax1+kr(id,1);
      ixAmin1=ixCmin1+kr(id,1);
      Hspeed(ixCmin1:ixCmax1,1)=max(Hspeed(ixCmin1:ixCmax1,1),&
         0.5d0*abs(wprim(ixAmin1:ixAmax1,mom(id))+csound(ixAmin1:ixAmax1,&
         id)-wprim(ixCmin1:ixCmax1,mom(id))+csound(ixCmin1:ixCmax1,id)))
      ixAmax1=ixCmax1-kr(id,1);
      ixAmin1=ixCmin1-kr(id,1);
      Hspeed(ixCmin1:ixCmax1,1)=max(Hspeed(ixCmin1:ixCmax1,1),&
         0.5d0*abs(wprim(ixCmin1:ixCmax1,mom(id))+csound(ixCmin1:ixCmax1,&
         id)-wprim(ixAmin1:ixAmax1,mom(id))+csound(ixAmin1:ixAmax1,id)))
    end do

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=jxCmax1+kr(id,1);
      ixAmin1=jxCmin1+kr(id,1);
      Hspeed(ixCmin1:ixCmax1,1)=max(Hspeed(ixCmin1:ixCmax1,1),&
         0.5d0*abs(wprim(ixAmin1:ixAmax1,mom(id))+csound(ixAmin1:ixAmax1,&
         id)-wprim(jxCmin1:jxCmax1,mom(id))+csound(jxCmin1:jxCmax1,id)))
      ixAmax1=jxCmax1-kr(id,1);
      ixAmin1=jxCmin1-kr(id,1);
      Hspeed(ixCmin1:ixCmax1,1)=max(Hspeed(ixCmin1:ixCmax1,1),&
         0.5d0*abs(wprim(jxCmin1:jxCmax1,mom(id))+csound(jxCmin1:jxCmax1,&
         id)-wprim(ixAmin1:ixAmax1,mom(id))+csound(ixAmin1:ixAmax1,id)))
    end do

  end subroutine uawsom_get_H_speed

  !> Estimating bounds for the minimum and maximum signal velocities without split
  subroutine uawsom_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImax1,ixOmin1,&
     ixOmax1,idim,Hspeed,cmax,cmin)
    use mod_global_parameters

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

    double precision :: wmean(ixImin1:ixImax1,nw)
    double precision, dimension(ixImin1:ixImax1) :: umean, dmean, csoundL,&
        csoundR, tmp1,tmp2,tmp3
    integer :: ix1

    select case (boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      tmp1(ixOmin1:ixOmax1)=sqrt(wLp(ixOmin1:ixOmax1,rho_))
      tmp2(ixOmin1:ixOmax1)=sqrt(wRp(ixOmin1:ixOmax1,rho_))
      tmp3(ixOmin1:ixOmax1)=1.d0/(tmp1(ixOmin1:ixOmax1)+tmp2(ixOmin1:ixOmax1))
      umean(ixOmin1:ixOmax1)=(wLp(ixOmin1:ixOmax1,&
         mom(idim))*tmp1(ixOmin1:ixOmax1)+wRp(ixOmin1:ixOmax1,&
         mom(idim))*tmp2(ixOmin1:ixOmax1))*tmp3(ixOmin1:ixOmax1)
      call uawsom_get_csound_prim(wLp,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundL)
      call uawsom_get_csound_prim(wRp,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundR)
      dmean(ixOmin1:ixOmax1)=(tmp1(ixOmin1:ixOmax1)*csoundL(ixOmin1:ixOmax1)**&
         2+tmp2(ixOmin1:ixOmax1)*csoundR(ixOmin1:ixOmax1)**&
         2)*tmp3(ixOmin1:ixOmax1)+0.5d0*tmp1(ixOmin1:ixOmax1)*tmp2(&
         ixOmin1:ixOmax1)*tmp3(ixOmin1:ixOmax1)**2*(wRp(ixOmin1:ixOmax1,&
         mom(idim))-wLp(ixOmin1:ixOmax1,mom(idim)))**2
      dmean(ixOmin1:ixOmax1)=sqrt(dmean(ixOmin1:ixOmax1))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,1)=umean(ixOmin1:ixOmax1)-dmean(ixOmin1:ixOmax1)
        cmax(ixOmin1:ixOmax1,1)=umean(ixOmin1:ixOmax1)+dmean(ixOmin1:ixOmax1)
        if(H_correction) then
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,1)=sign(one,cmin(ix1,1))*max(abs(cmin(ix1,1)),Hspeed(ix1,&
               1))
            cmax(ix1,1)=sign(one,cmax(ix1,1))*max(abs(cmax(ix1,1)),Hspeed(ix1,&
               1))
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,1)=abs(umean(ixOmin1:ixOmax1))+&
           dmean(ixOmin1:ixOmax1)
      end if
    case (2)
      wmean(ixOmin1:ixOmax1,1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,1:nwflux))
      tmp1(ixOmin1:ixOmax1)=wmean(ixOmin1:ixOmax1,&
         mom(idim))/wmean(ixOmin1:ixOmax1,rho_)
      call uawsom_get_csound(wmean,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,1)=max(tmp1(ixOmin1:ixOmax1)+&
           csoundR(ixOmin1:ixOmax1),zero)
        cmin(ixOmin1:ixOmax1,1)=min(tmp1(ixOmin1:ixOmax1)-&
           csoundR(ixOmin1:ixOmax1),zero)
        if(H_correction) then
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,1)=sign(one,cmin(ix1,1))*max(abs(cmin(ix1,1)),Hspeed(ix1,&
               1))
            cmax(ix1,1)=sign(one,cmax(ix1,1))*max(abs(cmax(ix1,1)),Hspeed(ix1,&
               1))
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,1)=abs(tmp1(ixOmin1:ixOmax1))+&
           csoundR(ixOmin1:ixOmax1)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call uawsom_get_csound_prim(wLp,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundL)
      call uawsom_get_csound_prim(wRp,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundR)
      csoundL(ixOmin1:ixOmax1)=max(csoundL(ixOmin1:ixOmax1),&
         csoundR(ixOmin1:ixOmax1))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,1)=min(wLp(ixOmin1:ixOmax1,mom(idim)),&
           wRp(ixOmin1:ixOmax1,mom(idim)))-csoundL(ixOmin1:ixOmax1)
        cmax(ixOmin1:ixOmax1,1)=max(wLp(ixOmin1:ixOmax1,mom(idim)),&
           wRp(ixOmin1:ixOmax1,mom(idim)))+csoundL(ixOmin1:ixOmax1)
        if(H_correction) then
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,1)=sign(one,cmin(ix1,1))*max(abs(cmin(ix1,1)),Hspeed(ix1,&
               1))
            cmax(ix1,1)=sign(one,cmax(ix1,1))*max(abs(cmax(ix1,1)),Hspeed(ix1,&
               1))
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,1)=max(wLp(ixOmin1:ixOmax1,mom(idim)),&
           wRp(ixOmin1:ixOmax1,mom(idim)))+csoundL(ixOmin1:ixOmax1)
      end if
    end select

  end subroutine uawsom_get_cbounds

  !> Estimating bounds for the minimum and maximum signal velocities with rho split
  subroutine uawsom_get_cbounds_split_rho(wLC,wRC,wLp,wRp,x,ixImin1,ixImax1,&
     ixOmin1,ixOmax1,idim,Hspeed,cmax,cmin)
    use mod_global_parameters

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

    double precision :: wmean(ixImin1:ixImax1,nw)
    double precision, dimension(ixImin1:ixImax1) :: umean, dmean, csoundL,&
        csoundR, tmp1,tmp2,tmp3
    integer :: ix1
    double precision :: rho(ixImin1:ixImax1)

    select case (boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      tmp1(ixOmin1:ixOmax1)=sqrt(wLp(ixOmin1:ixOmax1,&
         rho_)+block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i))
      tmp2(ixOmin1:ixOmax1)=sqrt(wRp(ixOmin1:ixOmax1,&
         rho_)+block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i))
      tmp3(ixOmin1:ixOmax1)=1.d0/(tmp1(ixOmin1:ixOmax1)+tmp2(ixOmin1:ixOmax1))
      umean(ixOmin1:ixOmax1)=(wLp(ixOmin1:ixOmax1,&
         mom(idim))*tmp1(ixOmin1:ixOmax1)+wRp(ixOmin1:ixOmax1,&
         mom(idim))*tmp2(ixOmin1:ixOmax1))*tmp3(ixOmin1:ixOmax1)
      call uawsom_get_csound_prim(wLp,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundL)
      call uawsom_get_csound_prim(wRp,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundR)
      dmean(ixOmin1:ixOmax1)=(tmp1(ixOmin1:ixOmax1)*csoundL(ixOmin1:ixOmax1)**&
         2+tmp2(ixOmin1:ixOmax1)*csoundR(ixOmin1:ixOmax1)**&
         2)*tmp3(ixOmin1:ixOmax1)+0.5d0*tmp1(ixOmin1:ixOmax1)*tmp2(&
         ixOmin1:ixOmax1)*tmp3(ixOmin1:ixOmax1)**2*(wRp(ixOmin1:ixOmax1,&
         mom(idim))-wLp(ixOmin1:ixOmax1,mom(idim)))**2
      dmean(ixOmin1:ixOmax1)=sqrt(dmean(ixOmin1:ixOmax1))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,1)=umean(ixOmin1:ixOmax1)-dmean(ixOmin1:ixOmax1)
        cmax(ixOmin1:ixOmax1,1)=umean(ixOmin1:ixOmax1)+dmean(ixOmin1:ixOmax1)
        if(H_correction) then
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,1)=sign(one,cmin(ix1,1))*max(abs(cmin(ix1,1)),Hspeed(ix1,&
               1))
            cmax(ix1,1)=sign(one,cmax(ix1,1))*max(abs(cmax(ix1,1)),Hspeed(ix1,&
               1))
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,1)=abs(umean(ixOmin1:ixOmax1))+&
           dmean(ixOmin1:ixOmax1)
      end if
    case (2)
      wmean(ixOmin1:ixOmax1,1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,1:nwflux))
      tmp1(ixOmin1:ixOmax1)=wmean(ixOmin1:ixOmax1,&
         mom(idim))/(wmean(ixOmin1:ixOmax1,&
         rho_)+block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i))
      call uawsom_get_csound(wmean,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,1)=max(tmp1(ixOmin1:ixOmax1)+&
           csoundR(ixOmin1:ixOmax1),zero)
        cmin(ixOmin1:ixOmax1,1)=min(tmp1(ixOmin1:ixOmax1)-&
           csoundR(ixOmin1:ixOmax1),zero)
        if(H_correction) then
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,1)=sign(one,cmin(ix1,1))*max(abs(cmin(ix1,1)),Hspeed(ix1,&
               1))
            cmax(ix1,1)=sign(one,cmax(ix1,1))*max(abs(cmax(ix1,1)),Hspeed(ix1,&
               1))
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,1)=abs(tmp1(ixOmin1:ixOmax1))+&
           csoundR(ixOmin1:ixOmax1)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call uawsom_get_csound_prim(wLp,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundL)
      call uawsom_get_csound_prim(wRp,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
         csoundR)
      csoundL(ixOmin1:ixOmax1)=max(csoundL(ixOmin1:ixOmax1),&
         csoundR(ixOmin1:ixOmax1))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,1)=min(wLp(ixOmin1:ixOmax1,mom(idim)),&
           wRp(ixOmin1:ixOmax1,mom(idim)))-csoundL(ixOmin1:ixOmax1)
        cmax(ixOmin1:ixOmax1,1)=max(wLp(ixOmin1:ixOmax1,mom(idim)),&
           wRp(ixOmin1:ixOmax1,mom(idim)))+csoundL(ixOmin1:ixOmax1)
        if(H_correction) then
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,1)=sign(one,cmin(ix1,1))*max(abs(cmin(ix1,1)),Hspeed(ix1,&
               1))
            cmax(ix1,1)=sign(one,cmax(ix1,1))*max(abs(cmax(ix1,1)),Hspeed(ix1,&
               1))
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,1)=max(wLp(ixOmin1:ixOmax1,mom(idim)),&
           wRp(ixOmin1:ixOmax1,mom(idim)))+csoundL(ixOmin1:ixOmax1)
      end if
    end select

  end subroutine uawsom_get_cbounds_split_rho

  !> prepare velocities for ct methods
  subroutine uawsom_get_ct_velocity(vcts,wLp,wRp,ixImin1,ixImax1,ixOmin1,&
     ixOmax1,idim,cmax,cmin)
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

  end subroutine uawsom_get_ct_velocity

  !> Calculate fast magnetosonic wave speed
  subroutine uawsom_get_csound(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
     csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in) :: w(ixImin1:ixImax1, nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1)
    double precision :: cfast2(ixImin1:ixImax1), AvMinCs2(ixImin1:ixImax1),&
        b2(ixImin1:ixImax1), kmax
    double precision :: inv_rho(ixOmin1:ixOmax1)

    if(has_equi_rho0) then
      inv_rho(ixOmin1:ixOmax1) = 1d0/(w(ixOmin1:ixOmax1,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i))
    else
      inv_rho(ixOmin1:ixOmax1) = 1d0/w(ixOmin1:ixOmax1,rho_)
    endif

    call uawsom_get_csound2(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,csound)

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1) = uawsom_mag_en_all(w,ixImin1,ixImax1,ixOmin1,ixOmax1)

    cfast2(ixOmin1:ixOmax1)   = b2(ixOmin1:ixOmax1) * &
       inv_rho+csound(ixOmin1:ixOmax1)
    AvMinCs2(ixOmin1:ixOmax1) = cfast2(ixOmin1:ixOmax1)**2-&
       4.0d0*csound(ixOmin1:ixOmax1) * uawsom_mag_i_all(w,ixImin1,ixImax1,&
       ixOmin1,ixOmax1,idim)**2 * inv_rho

    where(AvMinCs2(ixOmin1:ixOmax1)<zero)
       AvMinCs2(ixOmin1:ixOmax1)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1)=sqrt(AvMinCs2(ixOmin1:ixOmax1))

    if (.not. uawsom_Hall) then
       csound(ixOmin1:ixOmax1) = sqrt(half*(cfast2(ixOmin1:ixOmax1)+&
          AvMinCs2(ixOmin1:ixOmax1)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),bigdouble)*half
       csound(ixOmin1:ixOmax1) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1)+&
          AvMinCs2(ixOmin1:ixOmax1))), uawsom_etah * &
          sqrt(b2(ixOmin1:ixOmax1))*inv_rho*kmax)
    end if

  end subroutine uawsom_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine uawsom_get_csound_prim(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
     csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    double precision, intent(in) :: w(ixImin1:ixImax1, nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1)
    double precision :: cfast2(ixImin1:ixImax1), AvMinCs2(ixImin1:ixImax1),&
        b2(ixImin1:ixImax1), kmax
    double precision :: inv_rho(ixOmin1:ixOmax1)
    double precision :: tmp(ixImin1:ixImax1)

    call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,tmp)
    inv_rho(ixOmin1:ixOmax1) = 1d0/tmp(ixOmin1:ixOmax1)


    if(uawsom_energy) then
      if(has_equi_pe0) then
        csound(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
           e_) + block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,b0i)
      else
        csound(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,e_)
      endif
      csound(ixOmin1:ixOmax1)=uawsom_gamma*csound(ixOmin1:ixOmax1)*inv_rho
    else
      csound(ixOmin1:ixOmax1)=uawsom_gamma*uawsom_adiab*tmp(ixOmin1:ixOmax1)**&
         gamma_1
    end if

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1)        = uawsom_mag_en_all(w,ixImin1,ixImax1,ixOmin1,&
       ixOmax1)
    cfast2(ixOmin1:ixOmax1)   = b2(ixOmin1:ixOmax1) * &
       inv_rho+csound(ixOmin1:ixOmax1)
    AvMinCs2(ixOmin1:ixOmax1) = cfast2(ixOmin1:ixOmax1)**2-&
       4.0d0*csound(ixOmin1:ixOmax1) * uawsom_mag_i_all(w,ixImin1,ixImax1,&
       ixOmin1,ixOmax1,idim)**2 * inv_rho

    where(AvMinCs2(ixOmin1:ixOmax1)<zero)
       AvMinCs2(ixOmin1:ixOmax1)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1)=sqrt(AvMinCs2(ixOmin1:ixOmax1))

    if (.not. uawsom_Hall) then
       csound(ixOmin1:ixOmax1) = sqrt(half*(cfast2(ixOmin1:ixOmax1)+&
          AvMinCs2(ixOmin1:ixOmax1)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),bigdouble)*half
       csound(ixOmin1:ixOmax1) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1)+&
          AvMinCs2(ixOmin1:ixOmax1))), uawsom_etah * &
          sqrt(b2(ixOmin1:ixOmax1))*inv_rho*kmax)
    end if

  end subroutine uawsom_get_csound_prim

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine uawsom_get_pthermal_origin(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,&
     pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1)
    integer                      :: iw, ix1

    if(uawsom_energy) then
      if(uawsom_internal_e) then
        pth(ixOmin1:ixOmax1)=gamma_1*w(ixOmin1:ixOmax1,e_)
      else
        pth(ixOmin1:ixOmax1)=gamma_1*(w(ixOmin1:ixOmax1,e_)- uawsom_kin_en(w,&
           ixImin1,ixImax1,ixOmin1,ixOmax1)- uawsom_mag_en(w,ixImin1,ixImax1,&
           ixOmin1,ixOmax1)- uawsom_wk_en(w,ixImin1,ixImax1,ixOmin1,&
           ixOmax1)- uawsom_wA_en(w,ixImin1,ixImax1,ixOmin1,ixOmax1)) !Max: AW energy added
      end if
      if(has_equi_pe0) then
        pth(ixOmin1:ixOmax1) = pth(ixOmin1:ixOmax1) + &
           block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,b0i)
      endif
    else
      call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)
      pth(ixOmin1:ixOmax1)=uawsom_adiab*pth(ixOmin1:ixOmax1)**uawsom_gamma
    end if

    if (fix_small_values) then
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1)<small_pressure) then
            pth(ix1)=small_pressure
         end if
      enddo
    end if

    if (check_small_values) then
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix1),&
              " encountered when call uawsom_get_pthermal"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix1,:)
           write(*,*) "Cell number: ", ix1
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix1,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) sqrt(pth(ix1)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      enddo
    end if

  end subroutine uawsom_get_pthermal_origin

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine uawsom_get_pthermal_hde(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1)
    integer                      :: iw, ix1

    if(uawsom_energy) then
      pth(ixOmin1:ixOmax1)=gamma_1*(w(ixOmin1:ixOmax1,e_)-uawsom_kin_en(w,&
         ixImin1,ixImax1,ixOmin1,ixOmax1)) !Max: No wave energy contribution since we are considering only hydrodynamic 
    else
      pth(ixOmin1:ixOmax1)=uawsom_adiab*w(ixOmin1:ixOmax1,rho_)**uawsom_gamma
    end if

    if (fix_small_values) then
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1)<small_pressure) then
            pth(ix1)=small_pressure
         end if
      enddo
    end if

    if (check_small_values) then
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix1),&
              " encountered when call uawsom_get_pthermal_hde"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix1,:)
           write(*,*) "Cell number: ", ix1
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix1,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) sqrt(pth(ix1)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      enddo
    end if

  end subroutine uawsom_get_pthermal_hde

  !> Calculate temperature=p/rho when in e_ the internal energy is stored
  subroutine uawsom_get_temperature_from_eint(w, x, ixImin1,ixImax1, ixOmin1,&
     ixOmax1, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1)
    res(ixOmin1:ixOmax1) = gamma_1 * w(ixOmin1:ixOmax1, e_) /w(ixOmin1:ixOmax1,&
       rho_)
  end subroutine uawsom_get_temperature_from_eint

  !> Calculate temperature=p/rho when in e_ the total energy is stored
  !> this does not check the values of uawsom_energy and uawsom_internal_e,
  !>  uawsom_energy = .true. and uawsom_internal_e = .false.
  !> also check small_values is avoided
  subroutine uawsom_get_temperature_from_etot(w, x, ixImin1,ixImax1, ixOmin1,&
     ixOmax1, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1)
    res(ixOmin1:ixOmax1)=(gamma_1*(w(ixOmin1:ixOmax1,e_)- uawsom_kin_en(w,&
       ixImin1,ixImax1,ixOmin1,ixOmax1)- uawsom_mag_en(w,ixImin1,ixImax1,&
       ixOmin1,ixOmax1)- uawsom_wk_en(w,ixImin1,ixImax1,ixOmin1,&
       ixOmax1))- uawsom_wA_en(w,ixImin1,ixImax1,ixOmin1,&
       ixOmax1))/w(ixOmin1:ixOmax1,rho_) !Max: added AW energy here. What is res? Res seems to be Temperature? res = P/rho = T
  end subroutine uawsom_get_temperature_from_etot

  !> Calculate temperature from hydrodynamic energy
  subroutine uawsom_get_temperature_from_hde(w, x, ixImin1,ixImax1, ixOmin1,&
     ixOmax1, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1)
    res(ixOmin1:ixOmax1)=gamma_1*(w(ixOmin1:ixOmax1,e_)- uawsom_kin_en(w,&
       ixImin1,ixImax1,ixOmin1,ixOmax1))/w(ixOmin1:ixOmax1,rho_)
  end subroutine uawsom_get_temperature_from_hde

  subroutine uawsom_get_temperature_from_eint_with_equi(w, x, ixImin1,ixImax1,&
      ixOmin1,ixOmax1, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1)
    res(ixOmin1:ixOmax1) = (gamma_1 * w(ixOmin1:ixOmax1,&
        e_) + block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,&
       b0i)) /(w(ixOmin1:ixOmax1,rho_) +block%equi_vars(ixOmin1:ixOmax1,&
       equi_rho0_,b0i))
  end subroutine uawsom_get_temperature_from_eint_with_equi

  subroutine uawsom_get_temperature_equi(w,x, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1)
    res(ixOmin1:ixOmax1)= block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,&
       b0i)/block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i)
  end subroutine uawsom_get_temperature_equi

  subroutine uawsom_get_rho_equi(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1)
    res(ixOmin1:ixOmax1) = block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i)
  end subroutine uawsom_get_rho_equi

  subroutine uawsom_get_pe_equi(w,x, ixImin1,ixImax1, ixOmin1,ixOmax1, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1)
    res(ixOmin1:ixOmax1) = block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,b0i)
  end subroutine uawsom_get_pe_equi

  subroutine uawsom_get_temperature_from_etot_with_equi(w, x, ixImin1,ixImax1,&
      ixOmin1,ixOmax1, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1)
    res(ixOmin1:ixOmax1)=(gamma_1*(w(ixOmin1:ixOmax1,e_)- uawsom_kin_en(w,&
       ixImin1,ixImax1,ixOmin1,ixOmax1)- uawsom_mag_en(w,ixImin1,ixImax1,&
       ixOmin1,ixOmax1)) +  block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,&
       b0i))/(w(ixOmin1:ixOmax1,rho_) +block%equi_vars(ixOmin1:ixOmax1,&
       equi_rho0_,b0i))
            
  end subroutine uawsom_get_temperature_from_etot_with_equi

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine uawsom_get_csound2(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1)
    double precision    :: rho(ixImin1:ixImax1)
    
    call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,rho)
    if(uawsom_energy) then
      call uawsom_get_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,csound2)
      csound2(ixOmin1:ixOmax1)=uawsom_gamma*csound2(ixOmin1:ixOmax1)/rho(&
         ixOmin1:ixOmax1)
    else
      csound2(ixOmin1:ixOmax1)=uawsom_gamma*uawsom_adiab*rho(ixOmin1:ixOmax1)**&
         gamma_1
    end if
  end subroutine uawsom_get_csound2

  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine uawsom_get_p_total(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,p)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out)   :: p(ixImin1:ixImax1)

    call uawsom_get_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,p)

    p(ixOmin1:ixOmax1) = p(ixOmin1:ixOmax1) + 0.5d0 * sum(w(ixOmin1:ixOmax1,&
        mag(:))**2, dim=ndim+1)

  end subroutine uawsom_get_p_total

  !> Calculate fluxes within ixO^L without any splitting
  subroutine uawsom_get_flux(wC,w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,nwflux)

    double precision             :: ptotal(ixOmin1:ixOmax1)
    double precision             :: tmp(ixImin1:ixImax1),&
        zeta(ixImin1:ixImax1)
    double precision             :: vHall(ixImin1:ixImax1,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:,:) :: Jambi, btot
    double precision, allocatable, dimension(:) :: tmp2, tmp3

    ! Get flux of density
    f(ixOmin1:ixOmax1,rho_)=w(ixOmin1:ixOmax1,mom(idim))*w(ixOmin1:ixOmax1,&
       rho_)

    call get_zeta(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,zeta)

    if(uawsom_energy) then 
      ptotal(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,&
         p_)+0.5d0*sum(w(ixOmin1:ixOmax1,mag(:))**2,&
         dim=ndim+1)+(zeta(ixOmin1:ixOmax1)+1.d0)*(w(ixOmin1:ixOmax1,&
         wkplus_) + w(ixOmin1:ixOmax1,wkminus_))/4.d0+(w(ixOmin1:ixOmax1,&
         wAplus_) + w(ixOmin1:ixOmax1,wAminus_))/2.0d0 !Max: Added AW pressure
    else
      ptotal(ixOmin1:ixOmax1)=uawsom_adiab*w(ixOmin1:ixOmax1,&
         rho_)**uawsom_gamma+0.5d0*sum(w(ixOmin1:ixOmax1,mag(:))**2,&
         dim=ndim+1)
    end if

    if (uawsom_Hall) then
      call uawsom_getv_Hall(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,vHall)
    end if

    ! Get flux of tracer
    do iw=1,uawsom_n_tracer
      f(ixOmin1:ixOmax1,tracer(iw))=w(ixOmin1:ixOmax1,&
         mom(idim))*w(ixOmin1:ixOmax1,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixOmin1:ixOmax1,mom(idir))=wC(ixOmin1:ixOmax1,&
           mom(idir))*w(ixOmin1:ixOmax1,mom(idim))+&
           ptotal(ixOmin1:ixOmax1)-w(ixOmin1:ixOmax1,&
           mag(idir))*w(ixOmin1:ixOmax1,mag(idim))
      else
        f(ixOmin1:ixOmax1,mom(idir))=wC(ixOmin1:ixOmax1,&
           mom(idir))*w(ixOmin1:ixOmax1,mom(idim))-w(ixOmin1:ixOmax1,&
           mag(idir))*w(ixOmin1:ixOmax1,mag(idim))
      end if
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k) + Q-wkplus + Q+wkminus ! Max: + Z-wAplus + Z+wAminus for AWs TVD 2024 Eqn. 78
    if(uawsom_energy) then
      if (uawsom_internal_e) then
         f(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,mom(idim))*wC(ixOmin1:ixOmax1,&
            e_)
         ! Max: Adding AW to flux of energy. Note: mu (magnetic permeability) = 1 in cgs so does not appear here           
      else
        f(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,mom(idim))*(wC(ixOmin1:ixOmax1,&
           e_)+ptotal(ixOmin1:ixOmax1))-w(ixOmin1:ixOmax1,&
           mag(idim))*sum(w(ixOmin1:ixOmax1,mag(:))*w(ixOmin1:ixOmax1,mom(:)),&
           dim=ndim+1)+w(ixOmin1:ixOmax1,wkminus_)*(w(ixOmin1:ixOmax1,&
           mom(idim)) + w(ixOmin1:ixOmax1,mag(idim))/(w(ixOmin1:ixOmax1,&
           rho_)*(zeta(ixOmin1:ixOmax1)+1.d0)/2.0d0)**0.5d0)+w(ixOmin1:ixOmax1,&
           wkplus_)*(w(ixOmin1:ixOmax1,mom(idim)) - w(ixOmin1:ixOmax1,&
           mag(idim))/(w(ixOmin1:ixOmax1,&
           rho_)*(zeta(ixOmin1:ixOmax1)+1.d0)/2.0d0)**0.5d0)+w(ixOmin1:ixOmax1,&
           wAminus_)*(w(ixOmin1:ixOmax1,mom(idim)) + w(ixOmin1:ixOmax1,&
           mag(idim))/(w(ixOmin1:ixOmax1,rho_)**0.5d0))+w(ixOmin1:ixOmax1,&
           wAplus_)*(w(ixOmin1:ixOmax1,mom(idim)) - w(ixOmin1:ixOmax1,&
           mag(idim))/(w(ixOmin1:ixOmax1,rho_)**0.5d0))


        if(uawsom_solve_eaux) f(ixOmin1:ixOmax1,eaux_)=w(ixOmin1:ixOmax1,&
           mom(idim))*wC(ixOmin1:ixOmax1,eaux_)
        if(uawsom_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
          f(ixOmin1:ixOmax1,e_) = f(ixOmin1:ixOmax1,&
             e_) + vHall(ixOmin1:ixOmax1,idim) * sum(w(ixOmin1:ixOmax1,&
              mag(:))**2,dim=ndim+1) - w(ixOmin1:ixOmax1,&
             mag(idim)) * sum(vHall(ixOmin1:ixOmax1,:)*w(ixOmin1:ixOmax1,&
             mag(:)),dim=ndim+1)
        end if
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (uawsom_glm) then
           f(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,psi_)
        else
           f(ixOmin1:ixOmax1,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,&
           mom(idim))*w(ixOmin1:ixOmax1,mag(idir))-w(ixOmin1:ixOmax1,&
           mag(idim))*w(ixOmin1:ixOmax1,mom(idir))
        if (uawsom_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          f(ixOmin1:ixOmax1,mag(idir)) = f(ixOmin1:ixOmax1,&
             mag(idir)) - vHall(ixOmin1:ixOmax1,idir)*w(ixOmin1:ixOmax1,&
             mag(idim)) + vHall(ixOmin1:ixOmax1,idim)*w(ixOmin1:ixOmax1,&
             mag(idir))
        end if
      end if
    end do

    ! compute flux of wkminus Q+wkminus
    f(ixOmin1:ixOmax1,wkminus_)= w(ixOmin1:ixOmax1,&
       wkminus_)*(w(ixOmin1:ixOmax1,mom(idim)) + w(ixOmin1:ixOmax1,&
       mag(idim))/(w(ixOmin1:ixOmax1,rho_)*(zeta(ixOmin1:ixOmax1)+&
       1.d0)/2.d0)**0.5d0)

    ! compute flux of wkplus Q-wkplus 
    f(ixOmin1:ixOmax1,wkplus_)= w(ixOmin1:ixOmax1,wkplus_)*(w(ixOmin1:ixOmax1,&
       mom(idim)) - w(ixOmin1:ixOmax1,mag(idim))/(w(ixOmin1:ixOmax1,&
       rho_)*(zeta(ixOmin1:ixOmax1)+1.d0)/2.d0)**0.5d0)
                              
                      
    ! Max: compute flux of wAminus  Z+wAminus
    f(ixOmin1:ixOmax1,wAminus_)= w(ixOmin1:ixOmax1,&
       wAminus_)*(w(ixOmin1:ixOmax1,mom(idim)) + w(ixOmin1:ixOmax1,&
       mag(idim))/(w(ixOmin1:ixOmax1,rho_)**0.5d0))

    ! Max: compute flux of wAplus Z-wAplus
    f(ixOmin1:ixOmax1,wAplus_)= w(ixOmin1:ixOmax1,wAplus_)*(w(ixOmin1:ixOmax1,&
       mom(idim)) - w(ixOmin1:ixOmax1,mag(idim))/(w(ixOmin1:ixOmax1,&
       rho_)**0.5d0))

    if (uawsom_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(uawsom_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * uawsom_eta_ambi/rho**2
      allocate(Jambi(ixImin1:ixImax1,1:3))
      call uawsom_get_Jambi(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Jambi)
      allocate(btot(ixOmin1:ixOmax1,1:3))
      btot(ixOmin1:ixOmax1,1:3) = w(ixOmin1:ixOmax1,mag(1:3))
      allocate(tmp2(ixOmin1:ixOmax1),tmp3(ixOmin1:ixOmax1))
      !tmp2 = Btot^2
      tmp2(ixOmin1:ixOmax1) = sum(btot(ixOmin1:ixOmax1,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixOmin1:ixOmax1) = sum(Jambi(ixOmin1:ixOmax1,&
         :)*btot(ixOmin1:ixOmax1,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,&
             mag(3)) *Jambi(ixOmin1:ixOmax1,2) - w(ixOmin1:ixOmax1,&
             mag(2)) * Jambi(ixOmin1:ixOmax1,3)
          f(ixOmin1:ixOmax1,mag(2))= f(ixOmin1:ixOmax1,&
             mag(2)) - tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             3) + tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,3)
          f(ixOmin1:ixOmax1,mag(3))= f(ixOmin1:ixOmax1,&
             mag(3)) + tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             2) - tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,2)
        case(2)
          tmp(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,&
             mag(1)) *Jambi(ixOmin1:ixOmax1,3) - w(ixOmin1:ixOmax1,&
             mag(3)) * Jambi(ixOmin1:ixOmax1,1)
          f(ixOmin1:ixOmax1,mag(1))= f(ixOmin1:ixOmax1,&
             mag(1)) + tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             3) - tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,3)
          f(ixOmin1:ixOmax1,mag(3))= f(ixOmin1:ixOmax1,&
             mag(3)) - tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             1) + tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,1)
        case(3)
          tmp(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,&
             mag(2)) *Jambi(ixOmin1:ixOmax1,1) - w(ixOmin1:ixOmax1,&
             mag(1)) * Jambi(ixOmin1:ixOmax1,2)
          f(ixOmin1:ixOmax1,mag(1))= f(ixOmin1:ixOmax1,&
             mag(1)) - tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             2) + tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,2)
          f(ixOmin1:ixOmax1,mag(2))= f(ixOmin1:ixOmax1,&
             mag(2)) + tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             1) - tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,1)
      endselect

      if(uawsom_energy .and. .not. uawsom_internal_e) then
        f(ixOmin1:ixOmax1,e_) = f(ixOmin1:ixOmax1,&
           e_) + tmp2(ixOmin1:ixOmax1) *  tmp(ixOmin1:ixOmax1)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

  end subroutine uawsom_get_flux

  !> Calculate fluxes within ixO^L with possible splitting
  subroutine uawsom_get_flux_split(wC,w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idim,&
     f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,nwflux)

    double precision             :: pgas(ixOmin1:ixOmax1),&
        ptotal(ixOmin1:ixOmax1), B(ixOmin1:ixOmax1,1:ndir)
    double precision             :: tmp(ixImin1:ixImax1),&
        zeta(ixImin1:ixImax1)
    double precision             :: vHall(ixImin1:ixImax1,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:,:) :: Jambi, btot
    double precision, allocatable, dimension(:) :: tmp2, tmp3
    double precision :: tmp4(ixOmin1:ixOmax1)


    call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,tmp)
    ! Get flux of density
    f(ixOmin1:ixOmax1,rho_)=w(ixOmin1:ixOmax1,mom(idim))*tmp(ixOmin1:ixOmax1)
    ! pgas is time dependent only
    if(uawsom_energy) then
      pgas(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,p_)
    else
      pgas(ixOmin1:ixOmax1)=uawsom_adiab*tmp(ixOmin1:ixOmax1)**uawsom_gamma
      if(has_equi_pe0) then
        pgas(ixOmin1:ixOmax1)=pgas(ixOmin1:ixOmax1)-&
           block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,b0i)
      endif
    end if

    if (uawsom_Hall) then
      call uawsom_getv_Hall(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,vHall)
    end if

    if(B0field) then
      B(ixOmin1:ixOmax1,1:ndir)=w(ixOmin1:ixOmax1,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,1:ndir,idim)
      pgas(ixOmin1:ixOmax1)=pgas(ixOmin1:ixOmax1)+sum(w(ixOmin1:ixOmax1,&
         mag(:))*block%B0(ixOmin1:ixOmax1,:,idim),dim=ndim+1)
    else
      B(ixOmin1:ixOmax1,1:ndir)=w(ixOmin1:ixOmax1,mag(1:ndir))
    end if

    call get_zeta(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,zeta)

    ptotal(ixOmin1:ixOmax1)=pgas(ixOmin1:ixOmax1)+0.5d0*sum(w(ixOmin1:ixOmax1,&
       mag(:))**2,dim=ndim+1)+(zeta(ixOmin1:ixOmax1)+1.d0)*(w(ixOmin1:ixOmax1,&
       wkplus_) + w(ixOmin1:ixOmax1,wkminus_))/4.d0+(w(ixOmin1:ixOmax1,&
       wAplus_) + w(ixOmin1:ixOmax1,wAminus_))/2.d0 !Max: added AW pressure

    ! Get flux of tracer
    do iw=1,uawsom_n_tracer
      f(ixOmin1:ixOmax1,tracer(iw))=w(ixOmin1:ixOmax1,&
         mom(idim))*w(ixOmin1:ixOmax1,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    if(B0field) then
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,mom(idir))=wC(ixOmin1:ixOmax1,&
             mom(idir))*w(ixOmin1:ixOmax1,&
             mom(idim))+ptotal(ixOmin1:ixOmax1)-w(ixOmin1:ixOmax1,&
             mag(idir))*B(ixOmin1:ixOmax1,idim)-block%B0(ixOmin1:ixOmax1,idir,&
             idim)*w(ixOmin1:ixOmax1,mag(idim))
        else
          f(ixOmin1:ixOmax1,mom(idir))=wC(ixOmin1:ixOmax1,&
             mom(idir))*w(ixOmin1:ixOmax1,mom(idim))-w(ixOmin1:ixOmax1,&
             mag(idir))*B(ixOmin1:ixOmax1,idim)-block%B0(ixOmin1:ixOmax1,idir,&
             idim)*w(ixOmin1:ixOmax1,mag(idim))
        end if
      end do
    else
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,mom(idir))=wC(ixOmin1:ixOmax1,&
             mom(idir))*w(ixOmin1:ixOmax1,&
             mom(idim))+ptotal(ixOmin1:ixOmax1)-w(ixOmin1:ixOmax1,&
             mag(idir))*B(ixOmin1:ixOmax1,idim)
        else
          f(ixOmin1:ixOmax1,mom(idir))=wC(ixOmin1:ixOmax1,&
             mom(idir))*w(ixOmin1:ixOmax1,mom(idim))-w(ixOmin1:ixOmax1,&
             mag(idir))*B(ixOmin1:ixOmax1,idim)
        end if
      end do
    end if

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k) ! Max: Eqn 78 TVD 2024, inside divergence ie div((rho*v^2/2 + p/gamma-1 + B^2/2mu)v - B*v.B/mu)
    !                                       ! Max: Norbert has also added wplus and wminus (Kink waves) as well to this equation-terms in Eqn 78.
    ! Max: This subtraction in the flux calculation a few lines down is because wkplus carries wave energy into the Sun
    if(uawsom_energy) then
      if (uawsom_internal_e) then
         f(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,mom(idim))*wC(ixOmin1:ixOmax1,&
            e_)
      else
         f(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,&
            mom(idim))*(wC(ixOmin1:ixOmax1,&
            e_)+ptotal(ixOmin1:ixOmax1))-B(ixOmin1:ixOmax1,&
            idim)*sum(w(ixOmin1:ixOmax1,mag(:))*w(ixOmin1:ixOmax1,mom(:)),&
            dim=ndim+1)+w(ixOmin1:ixOmax1,wkminus_)*B(ixOmin1:ixOmax1,&
            idim)/((w(ixOmin1:ixOmax1,rho_)*(zeta(ixOmin1:ixOmax1)+&
            1.d0)/2.d0)**0.5d0) -w(ixOmin1:ixOmax1,wkplus_)*B(ixOmin1:ixOmax1,&
            idim)/((w(ixOmin1:ixOmax1,rho_)*(zeta(ixOmin1:ixOmax1)+&
            1.d0)/2.d0)**0.5d0) +w(ixOmin1:ixOmax1,wAminus_)*B(ixOmin1:ixOmax1,&
            idim)/(w(ixOmin1:ixOmax1,rho_)**0.5d0) -w(ixOmin1:ixOmax1,&
            wAplus_)*B(ixOmin1:ixOmax1,idim)/(w(ixOmin1:ixOmax1,rho_)**0.5d0)

        if(uawsom_solve_eaux) f(ixOmin1:ixOmax1,eaux_)=w(ixOmin1:ixOmax1,&
           mom(idim))*wC(ixOmin1:ixOmax1,eaux_)
        if (uawsom_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
           if (uawsom_etah>zero) then
              f(ixOmin1:ixOmax1,e_) = f(ixOmin1:ixOmax1,&
                 e_) + vHall(ixOmin1:ixOmax1,idim) * sum(w(ixOmin1:ixOmax1,&
                  mag(:))*B(ixOmin1:ixOmax1,:),dim=ndim+1) - B(ixOmin1:ixOmax1,&
                 idim) * sum(vHall(ixOmin1:ixOmax1,:)*w(ixOmin1:ixOmax1,&
                 mag(:)),dim=ndim+1)
           end if
        end if
      end if
      if(has_equi_pe0) then
        f(ixOmin1:ixOmax1,e_)=  f(ixOmin1:ixOmax1,e_) + w(ixOmin1:ixOmax1,&
           mom(idim)) * block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,&
           idim) * inv_gamma_1
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (uawsom_glm) then
           f(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,psi_)
        else
           f(ixOmin1:ixOmax1,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,&
           mom(idim))*B(ixOmin1:ixOmax1,idir)-B(ixOmin1:ixOmax1,&
           idim)*w(ixOmin1:ixOmax1,mom(idir))

        if (uawsom_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (uawsom_etah>zero) then
            f(ixOmin1:ixOmax1,mag(idir)) = f(ixOmin1:ixOmax1,&
               mag(idir)) - vHall(ixOmin1:ixOmax1,idir)*B(ixOmin1:ixOmax1,&
               idim) + vHall(ixOmin1:ixOmax1,idim)*B(ixOmin1:ixOmax1,idir)
          end if
        end if

      end if
    end do
     
    ! compute flux of wkminus
    f(ixOmin1:ixOmax1,wkminus_)= w(ixOmin1:ixOmax1,&
       wkminus_)*(w(ixOmin1:ixOmax1,mom(idim)) + B(ixOmin1:ixOmax1,&
       idim)/(w(ixOmin1:ixOmax1,rho_)*(zeta(ixOmin1:ixOmax1)+&
       1.d0)/2.d0)**0.5d0)

    ! compute flux of wkplus
    f(ixOmin1:ixOmax1,wkplus_)= w(ixOmin1:ixOmax1,wkplus_)*(w(ixOmin1:ixOmax1,&
       mom(idim)) - B(ixOmin1:ixOmax1,idim)/(w(ixOmin1:ixOmax1,&
       rho_)*(zeta(ixOmin1:ixOmax1)+1.d0)/2.d0)**0.5d0)

    ! Max: compute flux of wAminus
    f(ixOmin1:ixOmax1,wAminus_)= w(ixOmin1:ixOmax1,&
       wAminus_)*(w(ixOmin1:ixOmax1,mom(idim)) + B(ixOmin1:ixOmax1,&
       idim)/(w(ixOmin1:ixOmax1,rho_)**0.5d0))

    ! Max: compute flux of wAplus
    f(ixOmin1:ixOmax1,wAplus_)= w(ixOmin1:ixOmax1,wAplus_)*(w(ixOmin1:ixOmax1,&
       mom(idim)) - B(ixOmin1:ixOmax1,idim)/(w(ixOmin1:ixOmax1,rho_)**0.5d0))

    if (uawsom_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(uawsom_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * uawsom_eta_ambi/rho**2
      allocate(Jambi(ixImin1:ixImax1,1:3))
      call uawsom_get_Jambi(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Jambi)
      allocate(btot(ixOmin1:ixOmax1,1:3))
      if(B0field) then
        do idir=1,3
          btot(ixOmin1:ixOmax1, idir) = w(ixOmin1:ixOmax1,&
             mag(idir)) + block%B0(ixOmin1:ixOmax1,idir,idim)
        enddo
      else
        btot(ixOmin1:ixOmax1,1:3) = w(ixOmin1:ixOmax1,mag(1:3))
      endif
      allocate(tmp2(ixOmin1:ixOmax1),tmp3(ixOmin1:ixOmax1))
      !tmp2 = Btot^2
      tmp2(ixOmin1:ixOmax1) = sum(btot(ixOmin1:ixOmax1,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixOmin1:ixOmax1) = sum(Jambi(ixOmin1:ixOmax1,&
         :)*btot(ixOmin1:ixOmax1,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,&
             mag(3)) *Jambi(ixOmin1:ixOmax1,2) - w(ixOmin1:ixOmax1,&
             mag(2)) * Jambi(ixOmin1:ixOmax1,3)
          if(B0field) tmp4(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
             mag(2)) * btot(ixOmin1:ixOmax1,3) - w(ixOmin1:ixOmax1,&
             mag(3)) * btot(ixOmin1:ixOmax1,2)
          f(ixOmin1:ixOmax1,mag(2))= f(ixOmin1:ixOmax1,&
             mag(2)) - tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             3) + tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,3)
          f(ixOmin1:ixOmax1,mag(3))= f(ixOmin1:ixOmax1,&
             mag(3)) + tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             2) - tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,2)
        case(2)
          tmp(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,&
             mag(1)) *Jambi(ixOmin1:ixOmax1,3) - w(ixOmin1:ixOmax1,&
             mag(3)) * Jambi(ixOmin1:ixOmax1,1)
          if(B0field) tmp4(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
             mag(3)) * btot(ixOmin1:ixOmax1,1) - w(ixOmin1:ixOmax1,&
             mag(1)) * btot(ixOmin1:ixOmax1,3)
          f(ixOmin1:ixOmax1,mag(1))= f(ixOmin1:ixOmax1,&
             mag(1)) + tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             3) - tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,3)
          f(ixOmin1:ixOmax1,mag(3))= f(ixOmin1:ixOmax1,&
             mag(3)) - tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             1) + tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,1)
        case(3)
          tmp(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,&
             mag(2)) *Jambi(ixOmin1:ixOmax1,1) - w(ixOmin1:ixOmax1,&
             mag(1)) * Jambi(ixOmin1:ixOmax1,2)
          if(B0field) tmp4(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
             mag(1)) * btot(ixOmin1:ixOmax1,2) - w(ixOmin1:ixOmax1,&
             mag(2)) * btot(ixOmin1:ixOmax1,1)
          f(ixOmin1:ixOmax1,mag(1))= f(ixOmin1:ixOmax1,&
             mag(1)) - tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             2) + tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,2)
          f(ixOmin1:ixOmax1,mag(2))= f(ixOmin1:ixOmax1,&
             mag(2)) + tmp2(ixOmin1:ixOmax1) * Jambi(ixOmin1:ixOmax1,&
             1) - tmp3(ixOmin1:ixOmax1) * btot(ixOmin1:ixOmax1,1)
      endselect

      if(uawsom_energy .and. .not. uawsom_internal_e) then
        f(ixOmin1:ixOmax1,e_) = f(ixOmin1:ixOmax1,&
           e_) + tmp2(ixOmin1:ixOmax1) *  tmp(ixOmin1:ixOmax1)
        if(B0field) f(ixOmin1:ixOmax1,e_) = f(ixOmin1:ixOmax1,&
           e_) +  tmp3(ixOmin1:ixOmax1) *  tmp4(ixOmin1:ixOmax1)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

  end subroutine uawsom_get_flux_split

  !> Source terms J.E in internal energy.
  !> For the ambipolar term E = ambiCoef * JxBxB=ambiCoef * B^2(-J_perpB)
  !=> the source term J.E = ambiCoef * B^2 * J_perpB^2 = ambiCoef * JxBxB^2/B^2
  !> ambiCoef is calculated as uawsom_ambi_coef/rho^2,  see also the subroutine uawsom_get_Jambi
  subroutine add_source_ambipolar_internal_energy(qdt,ixImin1,ixImax1,ixOmin1,&
     ixOmax1,wCT,w,x,ie)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision :: tmp(ixImin1:ixImax1)
    double precision :: jxbxb(ixImin1:ixImax1,1:3)

    call uawsom_get_jxbxb(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,jxbxb)
    tmp(ixOmin1:ixOmax1) = sum(jxbxb(ixOmin1:ixOmax1,1:3)**2,&
       dim=ndim+1) / uawsom_mag_en_all(wCT, ixImin1,ixImax1, ixOmin1,ixOmax1)
    call multiplyAmbiCoef(ixImin1,ixImax1,ixOmin1,ixOmax1,tmp,wCT,x)
    w(ixOmin1:ixOmax1,ie)=w(ixOmin1:ixOmax1,ie)+qdt * tmp

  end subroutine add_source_ambipolar_internal_energy

  subroutine uawsom_get_jxbxb(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,res)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out)   :: res(:,:)

    double precision  :: btot(ixImin1:ixImax1,1:3)
    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3)
    double precision :: tmp(ixImin1:ixImax1),b2(ixImin1:ixImax1)

    res=0.d0
    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)
    !!!here we know that current has nonzero values only for components in the range idirmin, 3
 
    if(B0field) then
      do idir=1,3
        btot(ixOmin1:ixOmax1, idir) = w(ixOmin1:ixOmax1,&
           mag(idir)) + block%B0(ixOmin1:ixOmax1,idir,b0i)
      enddo
    else
      btot(ixOmin1:ixOmax1,1:3) = w(ixOmin1:ixOmax1,mag(1:3))
    endif

    tmp(ixOmin1:ixOmax1) = sum(current(ixOmin1:ixOmax1,&
       idirmin:3)*btot(ixOmin1:ixOmax1,idirmin:3),dim=ndim+1) !J.B
    b2(ixOmin1:ixOmax1) = sum(btot(ixOmin1:ixOmax1,1:3)**2,dim=ndim+1) !B2
    do idir=1,idirmin-1
      res(ixOmin1:ixOmax1,idir) = btot(ixOmin1:ixOmax1,&
         idir) * tmp(ixOmin1:ixOmax1)
    enddo
    do idir=idirmin,3
      res(ixOmin1:ixOmax1,idir) = btot(ixOmin1:ixOmax1,&
         idir) * tmp(ixOmin1:ixOmax1) - current(ixOmin1:ixOmax1,&
         idir) * b2(ixOmin1:ixOmax1)
    enddo
  end subroutine uawsom_get_jxbxb

  !> Sets the sources for the ambipolar
  !> this is used for the STS method
  ! The sources are added directly (instead of fluxes as in the explicit)
  !> at the corresponding indices
  !>  store_flux_var is explicitly called for each of the fluxes one by one
  subroutine sts_set_source_ambipolar(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,wres,&
     fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1,igrid,nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,1:nw),&
        w(ixImin1:ixImax1,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step

    double precision, dimension(ixImin1:ixImax1,1:3) :: tmp,ff
    double precision :: fluxall(ixImin1:ixImax1,1:nflux,1:ndim)
    double precision :: fE(ixImin1:ixImax1,7-2*ndim:3)
    double precision  :: btot(ixImin1:ixImax1,1:3),tmp2(ixImin1:ixImax1)
    integer :: i, ixAmin1,ixAmax1, ie_

    ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;

    fluxall=zero

    call uawsom_get_jxbxb(w,x,ixImin1,ixImax1,ixAmin1,ixAmax1,tmp)

    !set electric field in tmp: E=nuA * jxbxb, where nuA=-etaA/rho^2
    do i=1,3
      !tmp(ixA^S,i) = -(uawsom_eta_ambi/w(ixA^S, rho_)**2) * tmp(ixA^S,i)
      call multiplyAmbiCoef(ixImin1,ixImax1,ixAmin1,ixAmax1,&
         tmp(ixImin1:ixImax1,i),w,x)
    enddo

    if(uawsom_energy .and. .not.uawsom_internal_e) then
      !btot should be only mag. pert.
      btot(ixAmin1:ixAmax1,1:3)=0.d0
      !if(B0field) then
      !  do i=1,ndir
      !    btot(ixA^S, i) = w(ixA^S,mag(i)) + block%B0(ixA^S,i,0)
      !  enddo
      !else
        btot(ixAmin1:ixAmax1,1:ndir) = w(ixAmin1:ixAmax1,mag(1:ndir))
      !endif
      call cross_product(ixImin1,ixImax1,ixAmin1,ixAmax1,tmp,btot,ff)
      call get_flux_on_cell_face(ixImin1,ixImax1,ixOmin1,ixOmax1,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,1,&
         1:ndim)=ff(ixImin1:ixImax1,1:ndim)
      !- sign comes from the fact that the flux divergence is a source now
      wres(ixOmin1:ixOmax1,e_)=-tmp2(ixOmin1:ixOmax1)
    endif

    if(stagger_grid) then
      if(ndir>ndim) then
        !!!Bz
        ff(ixAmin1:ixAmax1,1) = tmp(ixAmin1:ixAmax1,2)
        ff(ixAmin1:ixAmax1,2) = -tmp(ixAmin1:ixAmax1,1)
        ff(ixAmin1:ixAmax1,3) = 0.d0
        call get_flux_on_cell_face(ixImin1,ixImax1,ixOmin1,ixOmax1,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,1+ndir,&
           1:ndim)=ff(ixImin1:ixImax1,1:ndim)
        wres(ixOmin1:ixOmax1,mag(ndir))=-tmp2(ixOmin1:ixOmax1)
      end if
      fE=0.d0
      call update_faces_ambipolar(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,tmp,fE,&
         btot)
      ixAmax1=ixOmax1;
      ixAmin1=ixOmin1-1;
      wres(ixAmin1:ixAmax1,mag(1:ndim))=-btot(ixAmin1:ixAmax1,1:ndim)
    else
      !write curl(ele) as the divergence
      !m1={0,ele[[3]],-ele[[2]]}
      !m2={-ele[[3]],0,ele[[1]]}
      !m3={ele[[2]],-ele[[1]],0}

      !!!Bx
      ff(ixAmin1:ixAmax1,1) = 0.d0
      ff(ixAmin1:ixAmax1,2) = tmp(ixAmin1:ixAmax1,3)
      ff(ixAmin1:ixAmax1,3) = -tmp(ixAmin1:ixAmax1,2)
      call get_flux_on_cell_face(ixImin1,ixImax1,ixOmin1,ixOmax1,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,2,&
         1:ndim)=ff(ixImin1:ixImax1,1:ndim)
      !flux divergence is a source now
      wres(ixOmin1:ixOmax1,mag(1))=-tmp2(ixOmin1:ixOmax1)
      !!!By
      ff(ixAmin1:ixAmax1,1) = -tmp(ixAmin1:ixAmax1,3)
      ff(ixAmin1:ixAmax1,2) = 0.d0
      ff(ixAmin1:ixAmax1,3) = tmp(ixAmin1:ixAmax1,1)
      call get_flux_on_cell_face(ixImin1,ixImax1,ixOmin1,ixOmax1,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,3,&
         1:ndim)=ff(ixImin1:ixImax1,1:ndim)
      wres(ixOmin1:ixOmax1,mag(2))=-tmp2(ixOmin1:ixOmax1)

      if(ndir==3) then
        !!!Bz
        ff(ixAmin1:ixAmax1,1) = tmp(ixAmin1:ixAmax1,2)
        ff(ixAmin1:ixAmax1,2) = -tmp(ixAmin1:ixAmax1,1)
        ff(ixAmin1:ixAmax1,3) = 0.d0
        call get_flux_on_cell_face(ixImin1,ixImax1,ixOmin1,ixOmax1,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,1+ndir,&
           1:ndim)=ff(ixImin1:ixImax1,1:ndim)
        wres(ixOmin1:ixOmax1,mag(ndir))=-tmp2(ixOmin1:ixOmax1)
      end if

    end if

    if(fix_conserve_at_step) then
      fluxall=my_dt*fluxall
      call store_flux(igrid,fluxall,1,ndim,nflux)
      if(stagger_grid) then
        call store_edge(igrid,ixImin1,ixImax1,my_dt*fE,1,ndim)
      end if
    end if

  end subroutine sts_set_source_ambipolar

  !> get ambipolar electric field and the integrals around cell faces
  subroutine update_faces_ambipolar(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,ECC,fE,&
     circ)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)       :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    ! amibipolar electric field at cell centers
    double precision, intent(in)       :: ECC(ixImin1:ixImax1,1:3)
    double precision, intent(out)      :: fE(ixImin1:ixImax1,7-2*ndim:3)
    double precision, intent(out)      :: circ(ixImin1:ixImax1,1:ndim)

    integer                            :: hxCmin1,hxCmax1,ixCmin1,ixCmax1,&
       ixAmin1,ixAmax1
    integer                            :: idim1,idim2,idir,ix1

    fE=zero
    ! calcuate ambipolar electric field on cell edges from cell centers
    do idir=7-2*ndim,3
      ixCmax1=ixOmax1;
      ixCmin1=ixOmin1+kr(idir,1)-1;
     do ix1=0,1
        if( ix1==1 .and. 1==idir ) cycle
        ixAmin1=ixCmin1+ix1;
        ixAmax1=ixCmax1+ix1;
        fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,idir)+ECC(ixAmin1:ixAmax1,&
           idir)
     end do
      fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
         idir)*0.25d0*block%dsC(ixCmin1:ixCmax1,idir)
    end do

    ! Calculate circulation on each face to get value of line integral of
    ! electric field in the positive idir direction.
    ixCmax1=ixOmax1;
    ixCmin1=ixOmin1-1;

    circ=zero

    do idim1=1,ndim ! Coordinate perpendicular to face
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxCmin1=ixCmin1-kr(idim2,1);hxCmax1=ixCmax1-kr(idim2,1);
          ! Add line integrals in direction idir
          circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,idim1)+lvc(idim1,&
             idim2,idir)*(fE(ixCmin1:ixCmax1,idir)-fE(hxCmin1:hxCmax1,idir))
        end do
      end do
      circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,&
         idim1)/block%surfaceC(ixCmin1:ixCmax1,idim1)
    end do

  end subroutine update_faces_ambipolar

  !> use cell-center flux to get cell-face flux
  !> and get the source term as the divergence of the flux
  subroutine get_flux_on_cell_face(ixImin1,ixImax1,ixOmin1,ixOmax1,ff,src)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, dimension(:,:), intent(inout) :: ff
    double precision, intent(out) :: src(ixImin1:ixImax1)

    double precision :: ffc(ixImin1:ixImax1,1:ndim)
    double precision :: dxinv(ndim)
    integer :: idims, ix1, ixAmin1,ixAmax1, ixBmin1,ixBmax1, ixCmin1,ixCmax1

    ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
    dxinv=1.d0/dxlevel
    ! cell corner flux in ffc
    ffc=0.d0
    ixCmax1=ixOmax1; ixCmin1=ixOmin1-1;
    do ix1=0,1
      ixBmin1=ixCmin1+ix1;
      ixBmax1=ixCmax1+ix1;
      ffc(ixCmin1:ixCmax1,1:ndim)=ffc(ixCmin1:ixCmax1,&
         1:ndim)+ff(ixBmin1:ixBmax1,1:ndim)
    end do
    ffc(ixCmin1:ixCmax1,1:ndim)=0.5d0**ndim*ffc(ixCmin1:ixCmax1,1:ndim)
    ! flux at cell face
    ff(ixImin1:ixImax1,1:ndim)=0.d0
    do idims=1,ndim
      ixBmin1=ixOmin1-kr(idims,1);ixBmax1=ixOmax1-kr(idims,1);
      ixCmax1=ixOmax1; ixCmin1=ixBmin1;
      do ix1=0,1 
         if( ix1==0 .and. 1==idims ) then
           ixBmin1=ixCmin1-ix1;
           ixBmax1=ixCmax1-ix1;
           ff(ixCmin1:ixCmax1,idims)=ff(ixCmin1:ixCmax1,&
              idims)+ffc(ixBmin1:ixBmax1,idims)
         end if
      end do
      ff(ixCmin1:ixCmax1,idims)=ff(ixCmin1:ixCmax1,idims)*0.5d0**(ndim-1)
    end do
    src=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        ff(ixAmin1:ixAmax1,idims)=dxinv(idims)*ff(ixAmin1:ixAmax1,idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmax1=ixOmax1-kr(idims,1);
        src(ixOmin1:ixOmax1)=src(ixOmin1:ixOmax1)+ff(ixOmin1:ixOmax1,&
           idims)-ff(ixBmin1:ixBmax1,idims)
      end do
    else
      do idims=1,ndim
        ff(ixAmin1:ixAmax1,idims)=ff(ixAmin1:ixAmax1,&
           idims)*block%surfaceC(ixAmin1:ixAmax1,idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmax1=ixOmax1-kr(idims,1);
        src(ixOmin1:ixOmax1)=src(ixOmin1:ixOmax1)+ff(ixOmin1:ixOmax1,&
           idims)-ff(ixBmin1:ixBmax1,idims)
      end do
      src(ixOmin1:ixOmax1)=src(ixOmin1:ixOmax1)/block%dvolume(ixOmin1:ixOmax1)
    end if
  end subroutine get_flux_on_cell_face

  !> Calculates the explicit dt for the ambipokar term
  !> This function is used by both explicit scheme and STS method
  function get_ambipolar_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dx1,&
     x)  result(dtnew)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: dx1, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw)
    double precision :: dtnew

    double precision              :: coef
    double precision              :: dxarr(ndim)
    double precision              :: tmp(ixImin1:ixImax1)

    dxarr(1)=dx1;
    tmp(ixOmin1:ixOmax1) = uawsom_mag_en_all(w, ixImin1,ixImax1, ixOmin1,&
       ixOmax1)
    call multiplyAmbiCoef(ixImin1,ixImax1,ixOmin1,ixOmax1,tmp,w,x)
    coef = maxval(abs(tmp(ixOmin1:ixOmax1)))
    if(coef/=0.d0) then
      coef=1.d0/coef
    else
      coef=bigdouble
    end if
    if(slab_uniform) then
      dtnew=minval(dxarr(1:ndim))**2.0d0*coef
    else
      dtnew=minval(block%ds(ixOmin1:ixOmax1,1:ndim))**2.0d0*coef
    end if

  end function get_ambipolar_dt

  !> multiply res by the ambipolar coefficient
  !> The ambipolar coefficient is calculated as -uawsom_eta_ambi/rho^2
  !> The user may mask its value in the user file
  !> by implemneting usr_mask_ambipolar subroutine
  subroutine multiplyAmbiCoef(ixImin1,ixImax1,ixOmin1,ixOmax1,res,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(inout) :: res(ixImin1:ixImax1)
    double precision :: tmp(ixImin1:ixImax1)
    double precision :: rho(ixImin1:ixImax1)

    call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,rho)
    tmp=0.d0
    tmp(ixOmin1:ixOmax1)=-uawsom_eta_ambi/rho(ixOmin1:ixOmax1)**2
    if (associated(usr_mask_ambipolar)) then
      call usr_mask_ambipolar(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,tmp)
    end if

    res(ixOmin1:ixOmax1) = tmp(ixOmin1:ixOmax1) * res(ixOmin1:ixOmax1)
  end subroutine multiplyAmbiCoef

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine uawsom_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
     qsourcesplit,active,wCTprim)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source
    use mod_cak_force, only: cak_add_source

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active
    double precision, intent(in), optional :: wCTprim(ixImin1:ixImax1,1:nw)

    if (.not. qsourcesplit) then
      if(uawsom_internal_e) then
        ! Source for solving internal energy
        active = .true.
        call internal_energy_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,&
           wCT,w,x,e_)
      else
        if(uawsom_solve_eaux) then
          ! Source for auxiliary internal energy equation
          active = .true.
          call internal_energy_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,&
             wCT,w,x,eaux_)
        endif
        if(has_equi_pe0) then
          active = .true.
          call add_pe0_divv(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
        endif
      endif

      ! Source for B0 splitting
      if (B0field.or.B0field) then
        active = .true.
        call add_source_B0split(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
      end if
    end if

      

    if(uawsom_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,&
         wCT,w,x,qsourcesplit,active, rc_fl)
    end if

    if(uawsom_viscosity) then
      call viscosity_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
         uawsom_energy,qsourcesplit,active)
    end if

    if(uawsom_gravity) then
      call gravity_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
         gravity_energy,qsourcesplit,active)
    end if

    call w_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
       uawsom_energy,qsourcesplit,active)

    if (uawsom_cak_force) then
      call cak_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,&
         uawsom_energy,qsourcesplit,active)
    end if

  end subroutine uawsom_add_source

  subroutine add_pe0_divv(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision                :: v(ixImin1:ixImax1,1:ndir)
    double precision                :: divv(ixImin1:ixImax1)

    call uawsom_get_v(wCT,x,ixImin1,ixImax1,ixImin1,ixImax1,v)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixImin1,ixImax1,ixOmin1,ixOmax1,divv,&
           sixthorder=.true.)
      else
        call divvector(v,ixImin1,ixImax1,ixOmin1,ixOmax1,divv,&
           fourthorder=.true.)
      end if
    else
     call divvector(v,ixImin1,ixImax1,ixOmin1,ixOmax1,divv)
    end if
    w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,&
       e_)-qdt*block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,&
       0)*divv(ixOmin1:ixOmax1)

  end subroutine add_pe0_divv

  subroutine smoothed_data(w,ixImin1,ixImax1,ixOmin1,ixOmax1,smoothed_rho,&
     sum_rho,half_window)
    use mod_global_parameters

    double precision :: w(ixImin1:ixImax1,rho_), smoothed_rho(ixImin1:ixImax1),&
        sum_rho
    integer          :: ix1, ix2, ixImin1,ixImax1, ixOmin1,ixOmax1,&
        half_window

    half_window = 10

    !smoothed_rho(:) = 0.d0
    do ix1 = ixOmin1, ixOmax1
      sum_rho = 0.d0
      do ix2 = max(ixOmin1,ix1 - half_window), min(ixOmax1, ix1 + half_window)
        sum_rho = sum_rho + w(ix2,rho_)
      end do
      smoothed_rho(ix1) = sum_rho/(dble(min(ix1 + half_window,&
          ixOmax1) - max(ixOmin1, ix1 - half_window)))
    end do

  end subroutine smoothed_data

  subroutine w_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x,energy,&
     qsourcesplit,active)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt
    double precision                :: Tmp(ixImin1:ixImax1)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
    double precision                :: v(ixImin1:ixImax1,1:ndir), qt
    double precision                :: divv(ixImin1:ixImax1),&
       Lperp(ixImin1:ixImax1),Lperp_AW(ixImin1:ixImax1),&
       Gamma_minus(ixImin1:ixImax1),Gamma_plus(ixImin1:ixImax1),&
       R(ixImin1:ixImax1),radius(ixImin1:ixImax1),zeta(ixImin1:ixImax1),&
       pth(ixImin1:ixImax1),B(ixImin1:ixImax1,3)
    double precision :: Bsqr(ixImin1:ixImax1),vA(ixImin1:ixImax1),&
       gvA(ixImin1:ixImax1),gradvA(ixImin1:ixImax1),refvA(ixImin1:ixImax1),&
       Bfix(ixImin1:ixImax1)
    integer :: half_window, ix1, ix2
    double precision :: smoothed_rho(ixImin1:ixImax1), sum_rho

    !> may use these later: ,B2(ixI^S),vA(ixI^S),gvA(ixI^S),gradvA(ixI^S),refvA(ixI^S)
  
    if(B0field) then
      B(ixOmin1:ixOmax1,1:ndir)=w(ixOmin1:ixOmax1,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,1:ndir,1)
    else
      B(ixOmin1:ixOmax1,1:ndir)=w(ixOmin1:ixOmax1,mag(1:ndir))
    end if

    !call smoothed_data(w,ixI^L,ixO^L,smoothed_rho,sum_rho,half_window)
    !vA(ixO^S) = B(ixO^S,1)/dsqrt(smoothed_rho(ixO^S))

    !call gradient(vA,ixO^L,ixO^L,ndim,gradvA)

    !if (mod(it,1000)==0) then
    !  write(*,*) 'gradvA = ', gradvA(ixO^S)
    !end if     

    !refvA(ixO^S) = (gradvA(ixO^S)**2)/(vA(ixO^S)**2)

    !if (mype==0 .and. mod(it,1000)==0) then
    !  write(*,*) "refvA with smooth rho = ", refvA(ixO^S)
    !end if
    
    radius(ixOmin1:ixOmax1) = 1.d8/unit_length * (Busr/B(ixOmin1:ixOmax1,&
       1))**0.5d0

    if(qsourcesplit .eqv. .false.) then
      active = .true.
    endif

    call uawsom_get_v(wCT,x,ixImin1,ixImax1,ixImin1,ixImax1,v)
    !call uawsom_get_pthermal_origin(w,x,ixI^L,ixO^L,pth)
    !pth(ixO^S) = pth(ixO^S)/w(ixO^S,rho_)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixImin1,ixImax1,ixOmin1,ixOmax1,divv,&
           sixthorder=.true.)
      else
        call divvector(v,ixImin1,ixImax1,ixOmin1,ixOmax1,divv,&
           fourthorder=.true.)
      end if
    else
     call divvector(v,ixImin1,ixImax1,ixOmin1,ixOmax1,divv)
    end if

    call get_zeta(w,x,ixImin1,ixImax1,ixImin1,ixImax1,zeta)

    !B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    !Bsqr(ixO^S) = B(ixO^S,1)**(2.d0)

    ! output Alfven wave speed B/sqrt(rho)
    !vA(ixI^S) = Bfix(ixI^S)/dsqrt(w(ixI^S,rho_))
    
    !call gradient(vA,ixI^L,ixO^L,ndim,gradvA)

    !refvA(ixO^S) = (gradvA(ixO^S)**2)/(vA(ixO^S)**2)

    !Lperp(ixO^S) = (zeta(ixO^S) + 1.d0 - ff)**(3.d0/2.d0)/(1.d0 - ff**(5.d0/2.d0))/&
    !               (zeta(ixO^S) - 1.d0)*3.1622776*(ff*dpi)**0.5d0*radius(ixO^S) !/min(pth(ixO^S),1.d0) 

    Lperp(ixOmin1:ixOmax1) = ((1/(radius(ixOmin1:ixOmax1)*(ff*dpi)**&
       0.5d0))*((2.d0/5.d0)**0.5d0)*((zeta(ixOmin1:ixOmax1)-1)/2)*((1-&
       ff**2.5d0)/(zeta(ixOmin1:ixOmax1)+1-ff)**1.5d0))**(-1.d0)

    ! Max: This is the last term in Eqn. 78 TVD 2024 contribution due to AWs is zero since mu*rho*alpha**2-1 = 0
    w(ixOmin1:ixOmax1,e_) = w(ixOmin1:ixOmax1,&
       e_)-qdt*(zeta(ixOmin1:ixOmax1)-1.d0)/(zeta(ixOmin1:ixOmax1)+&
       1.d0)*(zeta(ixOmin1:ixOmax1)+1.d0)*(wCT(ixOmin1:ixOmax1,&
       wkplus_) + wCT(ixOmin1:ixOmax1,wkminus_))/4.d0*divv(ixOmin1:ixOmax1) 
 
    ! Max: This is from the RHS of Eqn. 80 with no reflection term as yet
    !w(ixO^S,wkplus_) = w(ixO^S,wkplus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wkplus_)/2.d0 &
    !                  + wCT(ixO^S,wkplus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S)) !+ Reflection_kink_plus

    !w(ixO^S,wkminus_) = w(ixO^S,wkminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wkminus_)/2.d0 &
    !                  + wCT(ixO^S,wkminus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S)) !- Reflection_kink_minus (reverse sign as per Eqn. 80)


    !> Reflection term equal to grad(vA)^2/vA^2 

    !if (mype==0 .and. mod(it,100)==0) then 
    !  write(*,*) "it          = ", it
    !  write(*,*) "min wkminus = ", minval(w(ixO^S,wkminus_))
    !  write(*,*) "min wkplus  = ", minval(w(ixO^S,wkplus_))
    !  write(*,*) "Sum ref     = ", SUM(refvA(ixO^S))
    !  write(*,*) "Ref vA = ", refvA(ixO^S)
    !end if

    !if (it < 100) then
    !  w(ixO^S,wkplus_)  = w(ixO^S,wkplus_)  - qdt*(divv(ixO^S)*wCT(ixO^S,wkplus_)/2.d0  + wCT(ixO^S,wkplus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S))! - refvA(ixO^S)*w(ixO^S,wkminus_) + refvA(ixO^S)*w(ixO^S,wkplus_)) 
    !  w(ixO^S,wkminus_) = w(ixO^S,wkminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wkminus_)/2.d0 + wCT(ixO^S,wkminus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S))!+ refvA(ixO^S)*w(ixO^S,wkminus_) - refvA(ixO^S)*w(ixO^S,wkplus_))
    !else
    !  w(ixO^S,wkplus_)  = w(ixO^S,wkplus_)  - qdt*(divv(ixO^S)*wCT(ixO^S,wkplus_)/2.d0  + wCT(ixO^S,wkplus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S) - refvA(ixO^S)*w(ixO^S,wkminus_) + refvA(ixO^S)*w(ixO^S,wkplus_)) 
    !  w(ixO^S,wkminus_) = w(ixO^S,wkminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wkminus_)/2.d0 + wCT(ixO^S,wkminus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S)+ refvA(ixO^S)*w(ixO^S,wkminus_) - refvA(ixO^S)*w(ixO^S,wkplus_))
    !end if

    !w(ixO^S,wkplus_)  = w(ixO^S,wkplus_)  - qdt*(divv(ixO^S)*wCT(ixO^S,wkplus_)/2.d0  + wCT(ixO^S,wkplus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S))  !- refvA(ixO^S)*w(ixO^S,wkminus_) + refvA(ixO^S)*w(ixO^S,wkplus_)) 
    !w(ixO^S,wkminus_) = w(ixO^S,wkminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wkminus_)/2.d0 + wCT(ixO^S,wkminus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S)) !+ refvA(ixO^S)*w(ixO^S,wkminus_) - refvA(ixO^S)*w(ixO^S,wkplus_))
    
    !> Reflective terms have opposite signs now they are inside the qdt becayse its a -qdt so it swaps them round
    w(ixOmin1:ixOmax1,wkplus_)  = w(ixOmin1:ixOmax1,&
       wkplus_)  - qdt*(divv(ixOmin1:ixOmax1)*wCT(ixOmin1:ixOmax1,&
       wkplus_)/2.d0  + wCT(ixOmin1:ixOmax1,&
       wkplus_)**(3.d0/2.d0)/wCT(ixOmin1:ixOmax1,&
       rho_)**0.5d0/Lperp(ixOmin1:ixOmax1) - 28.14d0*exp(-((x(ixOmin1:ixOmax1,&
       1)-1.04d0)/0.01)**2.d0)*w(ixOmin1:ixOmax1,&
       wkminus_) + 28.14d0*exp(-((x(ixOmin1:ixOmax1,&
       1)-1.04d0)/0.01)**2.d0)*w(ixOmin1:ixOmax1,wkplus_)) 
    w(ixOmin1:ixOmax1,wkminus_) = w(ixOmin1:ixOmax1,&
       wkminus_) - qdt*(divv(ixOmin1:ixOmax1)*wCT(ixOmin1:ixOmax1,&
       wkminus_)/2.d0 + wCT(ixOmin1:ixOmax1,&
       wkminus_)**(3.d0/2.d0)/wCT(ixOmin1:ixOmax1,&
       rho_)**0.5d0/Lperp(ixOmin1:ixOmax1)+ 28.14d0*exp(-((x(ixOmin1:ixOmax1,&
       1)-1.04d0)/0.01)**2.d0)*w(ixOmin1:ixOmax1,&
       wkminus_) - 28.14d0*exp(-((x(ixOmin1:ixOmax1,&
       1)-1.04d0)/0.01)**2.d0)*w(ixOmin1:ixOmax1,wkplus_))
    
    !> Old reflection added on outside the timestep 
    !w(ixO^S,wkplus_)  = w(ixO^S,wkplus_)  - qdt*(divv(ixO^S)*wCT(ixO^S,wkplus_)/2.d0  + wCT(ixO^S,wkplus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S) + 0.2d0*exp(-((x(ixO^S,1)-1.07d0)/0.002)**2.d0)*w(ixO^S,wkminus_) - 0.2d0*exp(-((x(ixO^S,1)-1.07d0)/0.002)**2.d0)*w(ixO^S,wkplus_)) 
    !w(ixO^S,wkminus_) = w(ixO^S,wkminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wkminus_)/2.d0 + wCT(ixO^S,wkminus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S)- 0.2d0*exp(-((x(ixO^S,1)-1.07d0)/0.002)**2.d0)*w(ixO^S,wkminus_) + 0.2d0*exp(-((x(ixO^S,1)-1.07d0)/0.002)**2.d0)*w(ixO^S,wkplus_))

    call uawsom_get_pthermal(w,x,ixImin1,ixImax1,ixImin1,ixImax1,pth)
    Tmp(ixOmin1:ixOmax1) = pth(ixOmin1:ixOmax1)/w(ixOmin1:ixOmax1,rho_)

    !Lperp_AW(ixO^S) = (1/1.d8) * 1.5d5 * 1.0d2 * (Tmp(ixO^S) / B(ixO^S,1))**0.5d0 

    !if (mype==0 .and. mod(it,100000)==0) then
    !  write(*,*) "unit_length = ", unit_length
    !end if

    Lperp_AW(ixOmin1:ixOmax1) = 0.01d0*(20.0d0/B(ixOmin1:ixOmax1,1))**0.5d0 
    Gamma_plus(ixOmin1:ixOmax1) = (2.0d0 / Lperp_AW(ixOmin1:ixOmax1)) * &
       (w(ixOmin1:ixOmax1, wAminus_)/w(ixOmin1:ixOmax1,rho_))**0.5d0
    Gamma_minus(ixOmin1:ixOmax1) = (2.0d0 / Lperp_AW(ixOmin1:ixOmax1)) * &
       (w(ixOmin1:ixOmax1, wAplus_)/w(ixOmin1:ixOmax1,rho_))**0.5d0

    !if (mype==0 .and. mod(it,1000)==0) then
    !  write(*,*) "Max LperpAW = ", maxval(Lperp_AW(ixO^S))
    !  write(*,*) "Min LperpAW = ", minval(Lperp_AW(ixO^S))
    !  write(*,*) "Max Lperp   = ", maxval(Lperp(ixO^S))
    !  write(*,*) "Min Lperp   = ", minval(Lperp(ixO^S))
    !end if

    !if (mype==0 .and. mod(it,1000)==0) then
    !  write(*,*) "refvA(ixO^S)*w(ixO^S,wAminus_) = ", refvA(ixO^S)*w(ixO^S,wAminus_)
    !  write(*,*) "refvA(ixO^S)*w(ixO^S,wAplus_) = ", refvA(ixO^S)*w(ixO^S,wAplus_)
    !  write(*,*) "refvA(ixO^S)*w(ixO^S,wkminus_) = ", refvA(ixO^S)*w(ixO^S,wkminus_)
    !  write(*,*) "refvA(ixO^S)*w(ixO^S,wkplus_) = ", refvA(ixO^S)*w(ixO^S,wkplus_)
    !end if
    
    !> Reflection equal to grad(vA)^2/vA^2

    !if (it < 100) then
    !  w(ixO^S,wAplus_)  = w(ixO^S,wAplus_)  - qdt*(divv(ixO^S)*wCT(ixO^S,wAplus_)/2.0d0 + wCT(ixO^S,wAplus_)*Gamma_plus(ixO^S))!   - refvA(ixO^S)*w(ixO^S,wAminus_) + refvA(ixO^S)*w(ixO^S,wAplus_))
    !  w(ixO^S,wAminus_) = w(ixO^S,wAminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wAminus_)/2.0d0 + wCT(ixO^S,wAminus_)*Gamma_minus(ixO^S))!+ refvA(ixO^S)*w(ixO^S,wAminus_) - refvA(ixO^S)*w(ixO^S,wAplus_))
    !else 
    !  w(ixO^S,wAplus_)  = w(ixO^S,wAplus_)  - qdt*(divv(ixO^S)*wCT(ixO^S,wAplus_)/2.0d0 + wCT(ixO^S,wAplus_)*Gamma_plus(ixO^S)   - refvA(ixO^S)*w(ixO^S,wAminus_) + refvA(ixO^S)*w(ixO^S,wAplus_))
    !  w(ixO^S,wAminus_) = w(ixO^S,wAminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wAminus_)/2.0d0 + wCT(ixO^S,wAminus_)*Gamma_minus(ixO^S)+ refvA(ixO^S)*w(ixO^S,wAminus_) - refvA(ixO^S)*w(ixO^S,wAplus_))
    !end if  

    !w(ixO^S,wAplus_)  = w(ixO^S,wAplus_)  - qdt*(divv(ixO^S)*wCT(ixO^S,wAplus_)/2.0d0 + wCT(ixO^S,wAplus_)*Gamma_plus(ixO^S))    !- refvA(ixO^S)*w(ixO^S,wAminus_) + refvA(ixO^S)*w(ixO^S,wAplus_))
    !w(ixO^S,wAminus_) = w(ixO^S,wAminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wAminus_)/2.0d0 + wCT(ixO^S,wAminus_)*Gamma_minus(ixO^S)) !+ refvA(ixO^S)*w(ixO^S,wAminus_) - refvA(ixO^S)*w(ixO^S,wAplus_))
  
    ! Max: Add AWs for now setting Gammas to zero and no reflection to check the rest of the changes still allow for code to compile 
    !w(ixO^S,wAplus_) = w(ixO^S,wAplus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wAplus_)/2.0d0 + wCT(ixO^S,wAplus_)*Gamma_minus(ixO^S)) ! + Reflection_AW_plus 
    !w(ixO^S,wAminus_) = w(ixO^S,wAminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wAminus_)/2.0d0 + wCT(ixO^S,wAminus_)*Gamma_plus(ixO^S)) ! - Reflection_AW_minus 

    !> Reflection added inside the time step, since it is -qdt we have swapped the signs as to e.g., below where outside the qdt bracket
    w(ixOmin1:ixOmax1,wAplus_)  = w(ixOmin1:ixOmax1,&
       wAplus_)  - qdt*(divv(ixOmin1:ixOmax1)*wCT(ixOmin1:ixOmax1,&
       wAplus_)/2.0d0 + wCT(ixOmin1:ixOmax1,&
       wAplus_)*Gamma_plus(ixOmin1:ixOmax1)   - &
       28.14d0*exp(-((x(ixOmin1:ixOmax1,1)-&
       1.04d0)/0.01)**2.0d0)*w(ixOmin1:ixOmax1,&
       wAminus_) + 28.14d0*exp(-((x(ixOmin1:ixOmax1,&
       1)-1.04d0)/0.01)**2.d0)*w(ixOmin1:ixOmax1,wAplus_))
    w(ixOmin1:ixOmax1,wAminus_) = w(ixOmin1:ixOmax1,&
       wAminus_) - qdt*(divv(ixOmin1:ixOmax1)*wCT(ixOmin1:ixOmax1,&
       wAminus_)/2.0d0 + wCT(ixOmin1:ixOmax1,&
       wAminus_)*Gamma_minus(ixOmin1:ixOmax1)+ &
       28.14d0*exp(-((x(ixOmin1:ixOmax1,1)-&
       1.04d0)/0.01)**2.0d0)*w(ixOmin1:ixOmax1,&
       wAminus_) - 28.14d0*exp(-((x(ixOmin1:ixOmax1,&
       1)-1.04d0)/0.01)**2.d0)*w(ixOmin1:ixOmax1,wAplus_))
    
    !> Old reflection added outside of time step
    !w(ixO^S,wAplus_)  = w(ixO^S,wAplus_)  - qdt*(divv(ixO^S)*wCT(ixO^S,wAplus_)/2.0d0 + wCT(ixO^S,wAplus_)*Gamma_plus(ixO^S)   + 0.2d0*exp(-((x(ixO^S,1)-1.07d0)/0.002)**2.0d0)*w(ixO^S,wAminus_) - 0.2d0*exp(-((x(ixO^S,1)-1.07d0)/0.002)**2.d0)*w(ixO^S,wAplus_))
    !w(ixO^S,wAminus_) = w(ixO^S,wAminus_) - qdt*(divv(ixO^S)*wCT(ixO^S,wAminus_)/2.0d0 + wCT(ixO^S,wAminus_)*Gamma_minus(ixO^S)- 0.2d0*exp(-((x(ixO^S,1)-1.07d0)/0.002)**2.0d0)*w(ixO^S,wAminus_) + 0.2d0*exp(-((x(ixO^S,1)-1.07d0)/0.002)**2.d0)*w(ixO^S,wAplus_))
    
    !if (mype==0 .and. mod(it,4000) == 0) then
    !  write(*,*) "---------------------------------------------------------"
    !  write(*,*) "Max AW heating   = ", maxval(wCT(ixO^S,wAplus_)*Gamma_plus(ixO^S) + wCT(ixO^S,wAminus_)*Gamma_minus(ixO^S))
    !  write(*,*) "Min AW heating   = ", minval(wCT(ixO^S,wAplus_)*Gamma_plus(ixO^S) + wCT(ixO^S,wAminus_)*Gamma_minus(ixO^S))
    !  write(*,*) "Max kink heating = ", maxval(wCT(ixO^S,wkplus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S) + wCT(ixO^S,wkminus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S))
    !  write(*,*) "Min kink heating = ", minval(wCT(ixO^S,wkplus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S) + wCT(ixO^S,wkminus_)**(3.d0/2.d0)/wCT(ixO^S,rho_)**0.5d0/Lperp(ixO^S))
    !  write(*,*) "Kink wave heating MAX components: "
    !  write(*,*) "- wCT(ixO^S,wkplus_)  = ", maxval(wCT(ixO^S,wkplus_))
    !  write(*,*) "- wCT(ixO^S,wkminus_) = ", maxval(wCT(ixO^S,wkminus_))
    !  write(*,*) "- wCT(ixO^S,rho_)     = ", maxval(wCT(ixO^S,rho_))
    !  write(*,*) "- Lperp(ixO^S)        = ", maxval(Lperp(ixO^S))
    !  write(*,*) "Alfven wave heating MAX components: "
    !  write(*,*) "- wCT(ixO^S,wAplus_)  = ", maxval(wCT(ixO^S,wAplus_))
    !  write(*,*) "- wCT(ixO^S,wAminus_) = ", maxval(wCT(ixO^S,wAminus_))
    !  write(*,*) "- Gamma_plus(ixO^S)   = ", maxval(Gamma_plus(ixO^S))
    !  write(*,*) "- Gamma_minus(ixO^S)  = ", maxval(Gamma_minus(ixO^S))
    !end if

  end subroutine w_add_source

  !> Compute the Lorentz force (JxB)
  subroutine get_Lorentz_force(ixImin1,ixImax1,ixOmin1,ixOmax1,w,JxB)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(inout) :: JxB(ixImin1:ixImax1,3)
    double precision                :: a(ixImin1:ixImax1,3), b(ixImin1:ixImax1,&
       3)
    integer                         :: idir, idirmin
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for uawsom
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3)

    b=0.0d0
    do idir = 1, ndir
      b(ixOmin1:ixOmax1, idir) = uawsom_mag_i_all(w, ixImin1,ixImax1, ixOmin1,&
         ixOmax1,idir)
    end do

    ! store J current in a
    call get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)

    a=0.0d0
    do idir=7-2*ndir,3
      a(ixOmin1:ixOmax1,idir)=current(ixOmin1:ixOmax1,idir)
    end do

    call cross_product(ixImin1,ixImax1,ixOmin1,ixOmax1,a,b,JxB)
  end subroutine get_Lorentz_force

  subroutine internal_energy_add_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,&
     wCT,w,x,ie)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision                :: v(ixImin1:ixImax1,1:ndir),&
       divv(ixImin1:ixImax1)

    call uawsom_get_v(wCT,x,ixImin1,ixImax1,ixImin1,ixImax1,v)
    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixImin1,ixImax1,ixOmin1,ixOmax1,divv,&
           sixthorder=.true.)
      else
        call divvector(v,ixImin1,ixImax1,ixOmin1,ixOmax1,divv,&
           fourthorder=.true.)
      end if
    else
     call divvector(v,ixImin1,ixImax1,ixOmin1,ixOmax1,divv)
    end if
    w(ixOmin1:ixOmax1,ie)=w(ixOmin1:ixOmax1,ie)-qdt*wCT(ixOmin1:ixOmax1,&
       ie)*gamma_1*divv(ixOmin1:ixOmax1)
    if(uawsom_ambipolar)then
       call add_source_ambipolar_internal_energy(qdt,ixImin1,ixImax1,ixOmin1,&
          ixOmax1,wCT,w,x,ie)
    end if
    if(fix_small_values) then
      call uawsom_handle_small_ei(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,ie,&
         'internal_energy_add_source')
    end if
  end subroutine internal_energy_add_source

  subroutine uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,rho)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw),x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out) :: rho(ixImin1:ixImax1)

    if(has_equi_rho0) then
      rho(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i)
    else
      rho(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,rho_)
    endif

  end subroutine uawsom_get_rho
  
  subroutine get_zeta(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,zeta)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw),x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out) :: zeta(ixImin1:ixImax1)
    
    zeta(ixImin1:ixImax1) = zeta0*exp(-(x(ixImin1:ixImax1,&
       1)-xprobmin1)/6.957d10*unit_length/5.d0)
    where(zeta(ixImin1:ixImax1) < 1.d0)
      zeta(ixImin1:ixImax1) = 1.d0
    end where

  end subroutine get_zeta

  !> handle small or negative internal energy
  subroutine uawsom_handle_small_ei(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, ie,&
      subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1, ie
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixImin1:ixImax1,1:nw)
    double precision              :: rho(ixImin1:ixImax1)

    flag=.false.
    if(has_equi_pe0) then
      where(w(ixOmin1:ixOmax1,ie)+block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,&
         0)*inv_gamma_1<small_e)flag(ixOmin1:ixOmax1,ie)=.true.
    else
      where(w(ixOmin1:ixOmax1,ie)<small_e) flag(ixOmin1:ixOmax1,ie)=.true.
    endif
    if(any(flag(ixOmin1:ixOmax1,ie))) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_pe0) then
          where(flag(ixOmin1:ixOmax1,ie)) w(ixOmin1:ixOmax1,&
             ie)=small_e - block%equi_vars(ixOmin1:ixOmax1,equi_pe0_,&
             0)*inv_gamma_1
        else
          where(flag(ixOmin1:ixOmax1,ie)) w(ixOmin1:ixOmax1,ie)=small_e
        endif
      case ("average")
        call small_values_average(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x, flag,&
            ie)
      case default
        ! small values error shows primitive variables
        w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,e_)*gamma_1
        call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,rho)
        do idir = 1, ndir
           w(ixOmin1:ixOmax1, mom(idir)) = w(ixOmin1:ixOmax1,&
               mom(idir))/rho(ixOmin1:ixOmax1)
        end do
        call small_values_error(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, flag,&
            subname)
      end select
    end if

  end subroutine uawsom_handle_small_ei

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: a(ixImin1:ixImax1,3), b(ixImin1:ixImax1,3),&
        axb(ixImin1:ixImax1,3)
    integer :: idir

    a=0.d0
    b=0.d0
    ! for force-free field J0xB0 =0
    if(.not.B0field_forcefree) then
      ! store B0 magnetic field in b
      b(ixOmin1:ixOmax1,1:ndir)=block%B0(ixOmin1:ixOmax1,1:ndir,0)

      ! store J0 current in a
      do idir=7-2*ndir,3
        a(ixOmin1:ixOmax1,idir)=block%J0(ixOmin1:ixOmax1,idir)
      end do
      call cross_product(ixImin1,ixImax1,ixOmin1,ixOmax1,a,b,axb)
      axb(ixOmin1:ixOmax1,:)=axb(ixOmin1:ixOmax1,:)*qdt
      ! add J0xB0 source term in momentum equations
      w(ixOmin1:ixOmax1,mom(1:ndir))=w(ixOmin1:ixOmax1,&
         mom(1:ndir))+axb(ixOmin1:ixOmax1,1:ndir)
    end if

    if(total_energy) then
      a=0.d0
      ! for free-free field -(vxB0) dot J0 =0
      b(ixOmin1:ixOmax1,:)=wCT(ixOmin1:ixOmax1,mag(:))
      ! store full magnetic field B0+B1 in b
      if(.not.B0field_forcefree) b(ixOmin1:ixOmax1,:)=b(ixOmin1:ixOmax1,&
         :)+block%B0(ixOmin1:ixOmax1,:,0)
      ! store velocity in a
      call uawsom_get_v(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,&
         a(ixImin1:ixImax1,1:ndir))
      call cross_product(ixImin1,ixImax1,ixOmin1,ixOmax1,a,b,axb)
      axb(ixOmin1:ixOmax1,:)=axb(ixOmin1:ixOmax1,:)*qdt
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,e_)-axb(ixOmin1:ixOmax1,&
           idir)*block%J0(ixOmin1:ixOmax1,idir)
      end do
      if(uawsom_ambipolar) then
        !reuse axb
        call uawsom_get_jxbxb(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,axb)
        ! source J0 * E
        do idir=7-2*ndim,3
          !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
          call multiplyAmbiCoef(ixImin1,ixImax1,ixOmin1,ixOmax1,&
             axb(ixImin1:ixImax1,idir),wCT,x)
          w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,e_)+axb(ixOmin1:ixOmax1,&
             idir)*block%J0(ixOmin1:ixOmax1,idir)
        enddo
      endif
    end if

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,'add_source_B0')

  end subroutine add_source_B0split

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

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for uawsom
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3),&
       eta(ixImin1:ixImax1)
    double precision :: gradeta(ixImin1:ixImax1,1:ndim), Bf(ixImin1:ixImax1,&
       1:ndir)

    ! Calculating resistive sources involve one extra layer
    if (uawsom_4th_order) then
      ixAmin1=ixOmin1-2;ixAmax1=ixOmax1+2;
   else
      ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
    end if

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1) call &
       mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)

    if (uawsom_eta>zero)then
       eta(ixAmin1:ixAmax1)=uawsom_eta
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

    if(B0field) then
      Bf(ixImin1:ixImax1,1:ndir)=wCT(ixImin1:ixImax1,&
         mag(1:ndir))+block%B0(ixImin1:ixImax1,1:ndir,0)
    else
      Bf(ixImin1:ixImax1,1:ndir)=wCT(ixImin1:ixImax1,mag(1:ndir))
    end if

    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (uawsom_4th_order) then
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
       if (uawsom_eta<zero)then
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
       if(total_energy) then
          w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,&
             e_)+qdt*tmp(ixOmin1:ixOmax1)*Bf(ixOmin1:ixOmax1,idir)
       end if
    end do ! idir

    if(uawsom_energy) then
      ! de/dt+=eta*J**2
      tmp(ixOmin1:ixOmax1)=qdt*eta(ixOmin1:ixOmax1)*sum(current(&
         ixOmin1:ixOmax1,:)**2,dim=ndim+1)
      w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,e_)+tmp(ixOmin1:ixOmax1)
      if(uawsom_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixOmin1:ixOmax1,eaux_)=w(ixOmin1:ixOmax1,eaux_)+tmp(ixOmin1:ixOmax1)
      end if
    end if

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,'add_source_res1')

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

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for uawsom
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3),&
       eta(ixImin1:ixImax1),curlj(ixImin1:ixImax1,1:3)
    double precision :: tmpvec(ixImin1:ixImax1,1:3),tmp(ixOmin1:ixOmax1)
    integer :: ixAmin1,ixAmax1,idir,idirmin,idirmin1

    ixAmin1=ixOmin1-2;ixAmax1=ixOmax1+2;

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1) call &
       mpistop("Error in add_source_res2: Non-conforming input limits")

    ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixImin1,ixImax1,ixAmin1,ixAmax1,idirmin,current)

    tmpvec=zero
    if(uawsom_eta>zero)then
      do idir=idirmin,3
        tmpvec(ixAmin1:ixAmax1,idir)=current(ixAmin1:ixAmax1,idir)*uawsom_eta
      end do
    else
      call usr_special_resistivity(wCT,ixImin1,ixImax1,ixAmin1,ixAmax1,idirmin,&
         x,current,eta)
      do idir=idirmin,3
        tmpvec(ixAmin1:ixAmax1,idir)=current(ixAmin1:ixAmax1,&
           idir)*eta(ixAmin1:ixAmax1)
      end do
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
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

    if(uawsom_energy) then
      if(uawsom_eta>zero)then
        tmp(ixOmin1:ixOmax1)=qdt*uawsom_eta*sum(current(ixOmin1:ixOmax1,:)**2,&
           dim=ndim+1)
      else
        tmp(ixOmin1:ixOmax1)=qdt*eta(ixOmin1:ixOmax1)*sum(current(&
           ixOmin1:ixOmax1,:)**2,dim=ndim+1)
      end if
      if(total_energy) then
        ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
        ! de1/dt= eta J^2 - B1 dot curl(eta J)
        w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,&
           e_)+tmp(ixOmin1:ixOmax1)-qdt*sum(wCT(ixOmin1:ixOmax1,&
           mag(1:ndir))*curlj(ixOmin1:ixOmax1,1:ndir),dim=ndim+1)
      else
        ! add eta*J**2 source term in the internal energy equation
        w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,e_)+tmp(ixOmin1:ixOmax1)
      end if
      if(uawsom_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixOmin1:ixOmax1,eaux_)=w(ixOmin1:ixOmax1,eaux_)+tmp(ixOmin1:ixOmax1)
      end if
    end if

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,'add_source_res2')
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
       1:ndir)*uawsom_eta_hyper

    ixAmin1=ixOmin1;ixAmax1=ixOmax1;
    tmpvec2(ixAmin1:ixAmax1,1:ndir)=zero
    call curlvector(ehyper,ixImin1,ixImax1,ixAmin1,ixAmax1,tmpvec2,idirmin1,1,&
       3)

    do idir=1,ndir
      w(ixOmin1:ixOmax1,mag(idir)) = w(ixOmin1:ixOmax1,&
         mag(idir))-tmpvec2(ixOmin1:ixOmax1,idir)*qdt
    end do

    if(total_energy) then
      ! de/dt= +div(B x Ehyper)
      ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
      tmpvec2(ixAmin1:ixAmax1,1:ndir)=zero
      do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
        tmpvec2(ixAmin1:ixAmax1,idir) = tmpvec(ixAmin1:ixAmax1,idir)+ lvc(idir,&
           jdir,kdir)*wCT(ixAmin1:ixAmax1,mag(jdir))*ehyper(ixAmin1:ixAmax1,&
           kdir)
      end do; end do; end do
      tmp(ixOmin1:ixOmax1)=zero
      call divvector(tmpvec2,ixImin1,ixImax1,ixOmin1,ixOmax1,tmp)
      w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,e_)+tmp(ixOmin1:ixOmax1)*qdt
    end if

    if (fix_small_values)  call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,'add_source_hyperres')

  end subroutine add_source_hyperres

  subroutine add_source_glm(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-uawsom scheme or GLM-uawsom scheme
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision:: divb(ixImin1:ixImax1)
    integer          :: idim,idir
    double precision :: gradPsi(ixImin1:ixImax1)


    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (uawsom_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,psi_) = abs(uawsom_glm_alpha)*wCT(ixOmin1:ixOmax1,&
         psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixOmin1:ixOmax1,psi_) = dexp(&
           -qdt*cmax_global*uawsom_glm_alpha/minval(dxlevel(:)))*w(&
           ixOmin1:ixOmax1,psi_)
      else
        w(ixOmin1:ixOmax1,psi_) = dexp(&
           -qdt*cmax_global*uawsom_glm_alpha/minval(block%ds(ixOmin1:ixOmax1,&
           :),dim=ndim+1))*w(ixOmin1:ixOmax1,psi_)
      end if
    end if

    if(uawsom_glm_extended) then
      ! gradient of Psi
      if(total_energy) then
        do idim=1,ndim
          select case(typegrad)
          case("central")
            call gradient(wCT(ixImin1:ixImax1,psi_),ixImin1,ixImax1,ixOmin1,&
               ixOmax1,idim,gradPsi)
          case("limited")
            call gradientS(wCT(ixImin1:ixImax1,psi_),ixImin1,ixImax1,ixOmin1,&
               ixOmax1,idim,gradPsi)
          end select
          ! e  = e  -qdt (b . grad(Psi))
          w(ixOmin1:ixOmax1,e_) = w(ixOmin1:ixOmax1,&
             e_)-qdt*wCT(ixOmin1:ixOmax1,mag(idim))*gradPsi(ixOmin1:ixOmax1)
        end do
      end if

      ! We calculate now div B
      call get_divb(wCT,ixImin1,ixImax1,ixOmin1,ixOmax1,divb,&
          uawsom_divb_4thorder)

      ! m = m - qdt b div b
      do idir=1,ndir
        w(ixOmin1:ixOmax1,mom(idir))=w(ixOmin1:ixOmax1,&
           mom(idir))-qdt*uawsom_mag_i_all(w,ixImin1,ixImax1,ixOmin1,ixOmax1,&
           idir)*divb(ixOmin1:ixOmax1)
      end do
    end if

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,'add_source_glm')

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
    call get_divb(wCT,ixImin1,ixImax1,ixOmin1,ixOmax1,divb,&
        uawsom_divb_4thorder)

    ! calculate velocity
    call uawsom_get_v(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,v)

    if (total_energy) then
      ! e = e - qdt (v . b) * div b
      w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,e_)-qdt*sum(v(ixOmin1:ixOmax1,&
         :)*wCT(ixOmin1:ixOmax1,mag(:)),dim=ndim+1)*divb(ixOmin1:ixOmax1)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,&
         mag(idir))-qdt*v(ixOmin1:ixOmax1,idir)*divb(ixOmin1:ixOmax1)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,mom(idir))=w(ixOmin1:ixOmax1,&
         mom(idir))-qdt*uawsom_mag_i_all(w,ixImin1,ixImax1,ixOmin1,ixOmax1,&
         idir)*divb(ixOmin1:ixOmax1)
    end do

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision                :: divb(ixImin1:ixImax1),&
       vel(ixImin1:ixImax1,1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImax1,ixOmin1,ixOmax1,divb,&
        uawsom_divb_4thorder)

    call uawsom_get_v(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,vel)
    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,mag(idir))=w(ixOmin1:ixOmax1,&
         mag(idir))-qdt*vel(ixOmin1:ixOmax1,idir)*divb(ixOmin1:ixOmax1)
    end do

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,'add_source_janhunen')

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)
    ! Add Linde`s divB related sources to wnew within ixO
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
    call get_divb(wCT,ixImin1,ixImax1,ixpmin1,ixpmax1,divb,&
        uawsom_divb_4thorder)

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

    ! Add Linde`s diffusive terms
    do idim=1,ndim
       ! Calculate grad_idim(divb)
       select case(typegrad)
       case("central")
         call gradient(divb,ixImin1,ixImax1,ixpmin1,ixpmax1,idim,graddivb)
       case("limited")
         call gradientS(divb,ixImin1,ixImax1,ixpmin1,ixpmax1,idim,graddivb)
       end select

       ! Multiply by Linde`s eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab_uniform) then
          graddivb(ixpmin1:ixpmax1)=graddivb(ixpmin1:ixpmax1)*divbdiff/(&
             1.0d0/dxlevel(1)**2)
       else
          graddivb(ixpmin1:ixpmax1)=graddivb(ixpmin1:ixpmax1)*divbdiff &
             /(1.0d0/block%ds(ixpmin1:ixpmax1,1)**2)
       end if

       w(ixpmin1:ixpmax1,mag(idim))=w(ixpmin1:ixpmax1,&
          mag(idim))+graddivb(ixpmin1:ixpmax1)

       if (typedivbdiff=='all' .and. total_energy) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixpmin1:ixpmax1,e_)=w(ixpmin1:ixpmax1,e_)+wCT(ixpmin1:ixpmax1,&
            mag(idim))*graddivb(ixpmin1:ixpmax1)
       end if
    end do

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,'add_source_linde')

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
    invB(ixOmin1:ixOmax1)=sqrt(uawsom_mag_en_all(w,ixImin1,ixImax1,ixOmin1,&
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
  subroutine get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)  :: ixOmin1,ixOmax1, ixImin1,ixImax1
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw)
    integer, intent(out) :: idirmin
    integer :: idir, idirmin0

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for uawsom
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3)

    idirmin0 = 7-2*ndir

    call curlvector(w(ixImin1:ixImax1,mag(1:ndir)),ixImin1,ixImax1,ixOmin1,&
       ixOmax1,current,idirmin,idirmin0,ndir)

    if(B0field) current(ixOmin1:ixOmax1,idirmin0:3)=current(ixOmin1:ixOmax1,&
       idirmin0:3)+block%J0(ixOmin1:ixOmax1,idirmin0:3)
  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine uawsom_get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_cak_force, only: cak_get_dt

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
    if (uawsom_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/uawsom_eta
    else if (uawsom_eta<zero)then
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

    if(uawsom_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/uawsom_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixOmin1:ixOmax1,&
           1:ndim))**4/uawsom_eta_hyper,dtnew)
      end if
    end if

    if(uawsom_radiative_cooling) then
      call cooling_get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x,rc_fl)
    end if

    if(uawsom_viscosity) then
      call viscosity_get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
    end if

    if(uawsom_gravity) then
      call gravity_get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
    end if

    if(uawsom_ambipolar_exp) then
      dtnew=min(dtdiffpar*get_ambipolar_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,&
         dx1,x),dtnew)
    endif

    if (uawsom_cak_force) then
      call cak_get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
    end if

  end subroutine uawsom_get_dt

  ! Add geometrical source terms to w
  subroutine uawsom_add_source_geom(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,&
     x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,1:nw),&
        w(ixImin1:ixImax1,1:nw)

    integer          :: iw,idir, h1xmin1,h1xmax1
    double precision :: tmp(ixImin1:ixImax1),tmp1(ixImin1:ixImax1),&
       tmp2(ixImin1:ixImax1),invrho(ixOmin1:ixOmax1),invr(ixOmin1:ixOmax1)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    ! 1/rho
    invrho(ixOmin1:ixOmax1)=1.d0/wCT(ixOmin1:ixOmax1,rho_)
    invr(ixOmin1:ixOmax1)=1.d0/x(ixOmin1:ixOmax1,1)
    select case (coordinate)
    case (cylindrical)
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in uawsom")
      endif
      call uawsom_get_p_total(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,tmp)
      if(phi_>0) then
        w(ixOmin1:ixOmax1,mr_)=w(ixOmin1:ixOmax1,&
           mr_)+qdt*invr(ixOmin1:ixOmax1)*(tmp(ixOmin1:ixOmax1)-&
           wCT(ixOmin1:ixOmax1,bphi_)**2+wCT(ixOmin1:ixOmax1,&
           mphi_)**2*invrho(ixOmin1:ixOmax1))
        w(ixOmin1:ixOmax1,mphi_)=w(ixOmin1:ixOmax1,&
           mphi_)+qdt*invr(ixOmin1:ixOmax1)*(-wCT(ixOmin1:ixOmax1,&
           mphi_)*wCT(ixOmin1:ixOmax1,mr_)*invrho(ixOmin1:ixOmax1) &
           +wCT(ixOmin1:ixOmax1,bphi_)*wCT(ixOmin1:ixOmax1,br_))
        if(.not.stagger_grid) then
          w(ixOmin1:ixOmax1,bphi_)=w(ixOmin1:ixOmax1,&
             bphi_)+qdt*invr(ixOmin1:ixOmax1)*(wCT(ixOmin1:ixOmax1,&
             bphi_)*wCT(ixOmin1:ixOmax1,mr_) -wCT(ixOmin1:ixOmax1,&
             br_)*wCT(ixOmin1:ixOmax1,mphi_)) *invrho(ixOmin1:ixOmax1)
        end if
      else
        w(ixOmin1:ixOmax1,mr_)=w(ixOmin1:ixOmax1,&
           mr_)+qdt*invr(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1)
      end if
      if(uawsom_glm) w(ixOmin1:ixOmax1,br_)=w(ixOmin1:ixOmax1,&
         br_)+qdt*wCT(ixOmin1:ixOmax1,psi_)*invr(ixOmin1:ixOmax1)
    case (spherical)
       h1xmin1=ixOmin1-kr(1,1);h1xmax1=ixOmax1-kr(1,1); ;
       call uawsom_get_p_total(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,tmp1)
       ! m1
       tmp(ixOmin1:ixOmax1)=tmp1(ixOmin1:ixOmax1)*x(ixOmin1:ixOmax1,&
          1) *(block%surfaceC(ixOmin1:ixOmax1,&
          1)-block%surfaceC(h1xmin1:h1xmax1,1))/block%dvolume(ixOmin1:ixOmax1)
       do idir=2,ndir
         tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)+wCT(ixOmin1:ixOmax1,&
            mom(idir))**2*invrho(ixOmin1:ixOmax1)-wCT(ixOmin1:ixOmax1,&
            mag(idir))**2
       end do
       w(ixOmin1:ixOmax1,mom(1))=w(ixOmin1:ixOmax1,&
          mom(1))+qdt*tmp(ixOmin1:ixOmax1)*invr(ixOmin1:ixOmax1)
       ! b1
       if(uawsom_glm) then
         w(ixOmin1:ixOmax1,mag(1))=w(ixOmin1:ixOmax1,&
            mag(1))+qdt*invr(ixOmin1:ixOmax1)*2.0d0*wCT(ixOmin1:ixOmax1,psi_)
       end if

       

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixOmin1:ixOmax1)=-(wCT(ixOmin1:ixOmax1,&
              mom(3))*wCT(ixOmin1:ixOmax1,&
              mom(1))*invrho(ixOmin1:ixOmax1) -wCT(ixOmin1:ixOmax1,&
              mag(3))*wCT(ixOmin1:ixOmax1,mag(1))) 
           w(ixOmin1:ixOmax1,mom(3))=w(ixOmin1:ixOmax1,&
              mom(3))+qdt*tmp(ixOmin1:ixOmax1)*invr(ixOmin1:ixOmax1)
         else
           call mpistop("angmomfix not implemented yet in uawsom")
         end if
         ! b3
         if(.not.stagger_grid) then
           tmp(ixOmin1:ixOmax1)=(wCT(ixOmin1:ixOmax1,&
              mom(1))*wCT(ixOmin1:ixOmax1,mag(3)) -wCT(ixOmin1:ixOmax1,&
              mom(3))*wCT(ixOmin1:ixOmax1,mag(1)))*invrho(ixOmin1:ixOmax1) 
           w(ixOmin1:ixOmax1,mag(3))=w(ixOmin1:ixOmax1,&
              mag(3))+qdt*tmp(ixOmin1:ixOmax1)*invr(ixOmin1:ixOmax1)
         end if
       end if
    end select
  end subroutine uawsom_add_source_geom

  ! Add geometrical source terms to w
  subroutine uawsom_add_source_geom_split(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,&
     wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,1:nw),&
        w(ixImin1:ixImax1,1:nw)

    integer          :: iw,idir, h1xmin1,h1xmax1
    double precision :: tmp(ixImin1:ixImax1),tmp1(ixImin1:ixImax1),&
       tmp2(ixImin1:ixImax1),invrho(ixOmin1:ixOmax1),invr(ixOmin1:ixOmax1)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    if(has_equi_rho0) then
      invrho(ixOmin1:ixOmax1) = 1d0/(wCT(ixOmin1:ixOmax1,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,b0i))
    else
      invrho(ixOmin1:ixOmax1) = 1d0/wCT(ixOmin1:ixOmax1,rho_)
    end if
    invr(ixOmin1:ixOmax1)=1d0/x(ixOmin1:ixOmax1,1)

    select case (coordinate)
    case (cylindrical)
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in uawsom")
      endif
      call uawsom_get_p_total(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,tmp)
      if(phi_>0) then
        w(ixOmin1:ixOmax1,mr_)=w(ixOmin1:ixOmax1,&
           mr_)+qdt*invr(ixOmin1:ixOmax1)*(tmp(ixOmin1:ixOmax1)-&
           wCT(ixOmin1:ixOmax1,bphi_)**2+wCT(ixOmin1:ixOmax1,&
           mphi_)**2*invrho(ixOmin1:ixOmax1))
        w(ixOmin1:ixOmax1,mphi_)=w(ixOmin1:ixOmax1,&
           mphi_)+qdt*invr(ixOmin1:ixOmax1)*(-wCT(ixOmin1:ixOmax1,&
           mphi_)*wCT(ixOmin1:ixOmax1,mr_)*invrho(ixOmin1:ixOmax1) &
           +wCT(ixOmin1:ixOmax1,bphi_)*wCT(ixOmin1:ixOmax1,br_))
        if(.not.stagger_grid) then
          w(ixOmin1:ixOmax1,bphi_)=w(ixOmin1:ixOmax1,&
             bphi_)+qdt*invr(ixOmin1:ixOmax1)*(wCT(ixOmin1:ixOmax1,&
             bphi_)*wCT(ixOmin1:ixOmax1,mr_) -wCT(ixOmin1:ixOmax1,&
             br_)*wCT(ixOmin1:ixOmax1,mphi_)) *invrho(ixOmin1:ixOmax1)
        end if
      else
        w(ixOmin1:ixOmax1,mr_)=w(ixOmin1:ixOmax1,&
           mr_)+qdt*invr(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1)
      end if
      if(uawsom_glm) w(ixOmin1:ixOmax1,br_)=w(ixOmin1:ixOmax1,&
         br_)+qdt*wCT(ixOmin1:ixOmax1,psi_)*invr(ixOmin1:ixOmax1)
    case (spherical)
       h1xmin1=ixOmin1-kr(1,1);h1xmax1=ixOmax1-kr(1,1); ;
       call uawsom_get_p_total(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,tmp1)
       tmp(ixOmin1:ixOmax1)=tmp1(ixOmin1:ixOmax1)
       if(B0field) then
         tmp2(ixOmin1:ixOmax1)=sum(block%B0(ixOmin1:ixOmax1,:,&
            0)*wCT(ixOmin1:ixOmax1,mag(:)),dim=ndim+1)
         tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)+tmp2(ixOmin1:ixOmax1)
       end if
       ! m1
       tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)*x(ixOmin1:ixOmax1,&
          1) *(block%surfaceC(ixOmin1:ixOmax1,&
          1)-block%surfaceC(h1xmin1:h1xmax1,1))/block%dvolume(ixOmin1:ixOmax1)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)+wCT(ixOmin1:ixOmax1,&
              mom(idir))**2*invrho(ixOmin1:ixOmax1)-wCT(ixOmin1:ixOmax1,&
              mag(idir))**2
           if(B0field) tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)-&
              2.0d0*block%B0(ixOmin1:ixOmax1,idir,0)*wCT(ixOmin1:ixOmax1,&
              mag(idir))
         end do
       end if
       w(ixOmin1:ixOmax1,mom(1))=w(ixOmin1:ixOmax1,&
          mom(1))+qdt*tmp(ixOmin1:ixOmax1)*invr(ixOmin1:ixOmax1)
       ! b1
       if(uawsom_glm) then
         w(ixOmin1:ixOmax1,mag(1))=w(ixOmin1:ixOmax1,&
            mag(1))+qdt*invr(ixOmin1:ixOmax1)*2.0d0*wCT(ixOmin1:ixOmax1,psi_)
       end if

       

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixOmin1:ixOmax1)=-(wCT(ixOmin1:ixOmax1,&
              mom(3))*wCT(ixOmin1:ixOmax1,&
              mom(1))*invrho(ixOmin1:ixOmax1) -wCT(ixOmin1:ixOmax1,&
              mag(3))*wCT(ixOmin1:ixOmax1,mag(1))) 
           if (B0field) then
              tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)+&
                 block%B0(ixOmin1:ixOmax1,1,0)*wCT(ixOmin1:ixOmax1,&
                 mag(3)) +wCT(ixOmin1:ixOmax1,mag(1))*block%B0(ixOmin1:ixOmax1,&
                 3,0) 
           end if
           w(ixOmin1:ixOmax1,mom(3))=w(ixOmin1:ixOmax1,&
              mom(3))+qdt*tmp(ixOmin1:ixOmax1)*invr(ixOmin1:ixOmax1)
         else
           call mpistop("angmomfix not implemented yet in uawsom")
         end if
         ! b3
         if(.not.stagger_grid) then
           tmp(ixOmin1:ixOmax1)=(wCT(ixOmin1:ixOmax1,&
              mom(1))*wCT(ixOmin1:ixOmax1,mag(3)) -wCT(ixOmin1:ixOmax1,&
              mom(3))*wCT(ixOmin1:ixOmax1,mag(1)))*invrho(ixOmin1:ixOmax1) 
           if (B0field) then
              tmp(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)+(wCT(ixOmin1:ixOmax1,&
                 mom(1))*block%B0(ixOmin1:ixOmax1,3,0) -wCT(ixOmin1:ixOmax1,&
                 mom(3))*block%B0(ixOmin1:ixOmax1,1,&
                 0))*invrho(ixOmin1:ixOmax1)
           end if
           w(ixOmin1:ixOmax1,mag(3))=w(ixOmin1:ixOmax1,&
              mag(3))+qdt*tmp(ixOmin1:ixOmax1)*invr(ixOmin1:ixOmax1)
         end if
       end if
    end select
  end subroutine uawsom_add_source_geom_split

  !> Compute 2 times total magnetic energy
  function uawsom_mag_en_all(w, ixImin1,ixImax1, ixOmin1,ixOmax1) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: mge(ixOmin1:ixOmax1)

    if (B0field) then
      mge = sum((w(ixOmin1:ixOmax1, mag(:))+block%B0(ixOmin1:ixOmax1,:,&
         b0i))**2, dim=ndim+1)
    else
      mge = sum(w(ixOmin1:ixOmax1, mag(:))**2, dim=ndim+1)
    end if
  end function uawsom_mag_en_all

  !> Compute full magnetic field by direction
  function uawsom_mag_i_all(w, ixImin1,ixImax1, ixOmin1,ixOmax1,&
     idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1, idir
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: mgf(ixOmin1:ixOmax1)

    if (B0field) then
      mgf = w(ixOmin1:ixOmax1, mag(idir))+block%B0(ixOmin1:ixOmax1,idir,b0i)
    else
      mgf = w(ixOmin1:ixOmax1, mag(idir))
    end if
  end function uawsom_mag_i_all

  !> Compute evolving magnetic energy
  function uawsom_mag_en(w, ixImin1,ixImax1, ixOmin1,ixOmax1) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: mge(ixOmin1:ixOmax1)

    mge = 0.5d0 * sum(w(ixOmin1:ixOmax1, mag(:))**2, dim=ndim+1)
  end function uawsom_mag_en

  !> Compute wkminus and wkplus energy
  function uawsom_wk_en(w, ixImin1,ixImax1, ixOmin1,ixOmax1) result(we)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: we(ixOmin1:ixOmax1)

    we = w(ixOmin1:ixOmax1,wkminus_) + w(ixOmin1:ixOmax1,wkplus_)

  end function uawsom_wk_en

  !> Compute wAminus and wAplus energy
  function uawsom_wA_en(w, ixImin1,ixImax1, ixOmin1,ixOmax1) result(we)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: we(ixOmin1:ixOmax1)

    we = w(ixOmin1:ixOmax1,wAminus_) + w(ixOmin1:ixOmax1,wAplus_)

  end function uawsom_wA_en

  !> compute kinetic energy
  function uawsom_kin_en_origin(w, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: ke(ixOmin1:ixOmax1)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1)

    if (present(inv_rho)) then
      ke = 0.5d0 * sum(w(ixOmin1:ixOmax1, mom(:))**2, dim=ndim+1) * inv_rho
    else
      if(has_equi_rho0) then
        ke(ixOmin1:ixOmax1) = 0.5d0 * sum(w(ixOmin1:ixOmax1, mom(:))**2,&
            dim=ndim+1) / (w(ixOmin1:ixOmax1,&
            rho_) + block%equi_vars(ixOmin1:ixOmax1,equi_rho0_,0))
      else
        ke(ixOmin1:ixOmax1) = 0.5d0 * sum(w(ixOmin1:ixOmax1, mom(:))**2,&
            dim=ndim+1) / w(ixOmin1:ixOmax1, rho_)
      endif
    end if
  end function uawsom_kin_en_origin

  !> compute kinetic energy
  function uawsom_kin_en_boris(w, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      inv_rho) result(ke)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1, nw)
    double precision              :: ke(ixOmin1:ixOmax1)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1)

    if (present(inv_rho)) then
      ke=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,mag(:))**2,&
         dim=ndim+1)*inv_rho*inv_squared_c)
      ke=0.5d0*sum((w(ixOmin1:ixOmax1, mom(:)))**2,dim=ndim+1)*ke**2*inv_rho
    else
      ke=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,mag(:))**2,&
         dim=ndim+1)/w(ixOmin1:ixOmax1,rho_)*inv_squared_c)
      ke=0.5d0*sum(w(ixOmin1:ixOmax1, mom(:))**2,&
         dim=ndim+1)*ke**2/w(ixOmin1:ixOmax1, rho_)
    end if
  end function uawsom_kin_en_boris

  subroutine uawsom_getv_Hall(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: vHall(ixImin1:ixImax1,1:3)

    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3)
    double precision :: rho(ixImin1:ixImax1)

    call uawsom_get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,rho)
    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)
    vHall(ixOmin1:ixOmax1,1:3) = zero
    vHall(ixOmin1:ixOmax1,idirmin:3) = - uawsom_etah*current(ixOmin1:ixOmax1,&
       idirmin:3)
    do idir = idirmin, 3
       vHall(ixOmin1:ixOmax1,idir) = vHall(ixOmin1:ixOmax1,&
          idir)/rho(ixOmin1:ixOmax1)
    end do

  end subroutine uawsom_getv_Hall

  subroutine uawsom_get_Jambi(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,res)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, allocatable, intent(inout) :: res(:,:)


    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,7-2*ndir:3)

    res = 0d0

    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,current)
 
    res(ixOmin1:ixOmax1,idirmin:3)=-current(ixOmin1:ixOmax1,idirmin:3)
    do idir = idirmin, 3
      call multiplyAmbiCoef(ixImin1,ixImax1,ixOmin1,ixOmax1,&
         res(ixImin1:ixImax1,idir),w,x)
    enddo

  end subroutine uawsom_get_Jambi

    ! COMMENTED  because we have that in cmax now:
!  subroutine uawsom_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
!    use mod_global_parameters
!
!    integer, intent(in) :: ixI^L, ixO^L
!    double precision, intent(in)    :: dx^D
!    double precision, intent(in)    :: w(ixI^S,1:nw)
!    double precision, intent(in)    :: x(ixI^S,1:ndim)
!    double precision, intent(out)   :: dthall
!    !.. local ..
!    double precision :: dxarr(ndim)
!    double precision :: bmag(ixI^S)
!
!    dthall=bigdouble
!
!
!    ^D&dxarr(^D)=dx^D;
!
!    if (.not. B0field) then
!       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
!    else
!       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,b0i))**2))
!    end if
!
!    if(slab_uniform) then
!      dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(uawsom_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    else
!      dthall=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2.0d0/(uawsom_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    end if
!
!  end subroutine uawsom_getdt_Hall

  subroutine uawsom_modify_wLR(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,wLC,wRC,wLp,&
     wRp,s,idir)
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

    if(stagger_grid) then
      wLC(ixOmin1:ixOmax1,mag(idir))=s%ws(ixOmin1:ixOmax1,idir)
      wRC(ixOmin1:ixOmax1,mag(idir))=s%ws(ixOmin1:ixOmax1,idir)
      wLp(ixOmin1:ixOmax1,mag(idir))=s%ws(ixOmin1:ixOmax1,idir)
      wRp(ixOmin1:ixOmax1,mag(idir))=s%ws(ixOmin1:ixOmax1,idir)
    else
      ! Solve the Riemann problem for the linear 2x2 system for normal
      ! B-field and GLM_Psi according to Dedner 2002:
      ! This implements eq. (42) in Dedner et al. 2002 JcP 175
      ! Gives the Riemann solution on the interface
      ! for the normal B component and Psi in the GLM-uawsom system.
      ! 23/04/2013 Oliver Porth
      dB(ixOmin1:ixOmax1)   = wRp(ixOmin1:ixOmax1,&
         mag(idir)) - wLp(ixOmin1:ixOmax1,mag(idir))
      dPsi(ixOmin1:ixOmax1) = wRp(ixOmin1:ixOmax1,psi_) - wLp(ixOmin1:ixOmax1,&
         psi_)

      wLp(ixOmin1:ixOmax1,mag(idir))   = 0.5d0 * (wRp(ixOmin1:ixOmax1,&
         mag(idir)) + wLp(ixOmin1:ixOmax1,&
         mag(idir))) - 0.5d0/cmax_global * dPsi(ixOmin1:ixOmax1)
      wLp(ixOmin1:ixOmax1,psi_)       = 0.5d0 * (wRp(ixOmin1:ixOmax1,&
         psi_) + wLp(ixOmin1:ixOmax1,psi_)) - 0.5d0*cmax_global * &
         dB(ixOmin1:ixOmax1)

      wRp(ixOmin1:ixOmax1,mag(idir)) = wLp(ixOmin1:ixOmax1,mag(idir))
      wRp(ixOmin1:ixOmax1,psi_) = wLp(ixOmin1:ixOmax1,psi_)

      if(total_energy) then
        wRC(ixOmin1:ixOmax1,e_)=wRC(ixOmin1:ixOmax1,&
           e_)-half*wRC(ixOmin1:ixOmax1,mag(idir))**2
        wLC(ixOmin1:ixOmax1,e_)=wLC(ixOmin1:ixOmax1,&
           e_)-half*wLC(ixOmin1:ixOmax1,mag(idir))**2
      end if
      wRC(ixOmin1:ixOmax1,mag(idir)) = wLp(ixOmin1:ixOmax1,mag(idir))
      wRC(ixOmin1:ixOmax1,psi_) = wLp(ixOmin1:ixOmax1,psi_)
      wLC(ixOmin1:ixOmax1,mag(idir)) = wLp(ixOmin1:ixOmax1,mag(idir))
      wLC(ixOmin1:ixOmax1,psi_) = wLp(ixOmin1:ixOmax1,psi_)
      ! modify total energy according to the change of magnetic field
      if(total_energy) then
        wRC(ixOmin1:ixOmax1,e_)=wRC(ixOmin1:ixOmax1,&
           e_)+half*wRC(ixOmin1:ixOmax1,mag(idir))**2
        wLC(ixOmin1:ixOmax1,e_)=wLC(ixOmin1:ixOmax1,&
           e_)+half*wLC(ixOmin1:ixOmax1,mag(idir))**2
      end if
    end if

    if(associated(usr_set_wLR)) call usr_set_wLR(ixImin1,ixImax1,ixOmin1,&
       ixOmax1,qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine uawsom_modify_wLR

  subroutine uawsom_boundary_adjust(igrid,psb)
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

  end subroutine uawsom_boundary_adjust

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
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmax1,ixOmin1,&
          ixOmax1,w,x)
       
       
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmax1,ixOmin1,&
          ixOmax1,w,x)
     case(2)
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmax1,ixOmin1,&
          ixOmax1,w,x)
       
       
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmax1,ixOmin1,&
          ixOmax1,w,x)
     case(3)
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmax1,ixOmin1,&
          ixOmax1,w,x)
       
       
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmax1,ixOmin1,&
          ixOmax1,w,x)
     case(4)
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmax1,ixOmin1,&
          ixOmax1,w,x)
       
       
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmax1,ixOmin1,&
          ixOmax1,w,x)
     
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  

  subroutine uawsom_update_faces(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,qdt,wprim,&
     fC,fE,sCT,s,vcts)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)       :: qt,qdt
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

  end subroutine uawsom_update_faces

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
    ! non-ideal electric field on cell edges
    double precision, dimension(ixImin1:ixImax1,7-2*ndim:3) :: E_resi, E_ambi

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! if there is resistivity, get eta J
    if(uawsom_eta/=zero) call get_resistive_electric_field(ixImin1,ixImax1,&
       ixOmin1,ixOmax1,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(uawsom_ambipolar_exp) call get_ambipolar_electric_field(ixImin1,ixImax1,&
       ixOmin1,ixOmax1,sCT%w,x,E_ambi)

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

            ! add resistive electric field at cell edges E=-vxB+eta J
            if(uawsom_eta/=zero) fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
               idir)+E_resi(ixCmin1:ixCmax1,idir)
            ! add ambipolar electric field
            if(uawsom_ambipolar_exp) fE(ixCmin1:ixCmax1,&
               idir)=fE(ixCmin1:ixCmax1,idir)+E_ambi(ixCmin1:ixCmax1,idir)

            fE(ixCmin1:ixCmax1,idir)=qdt*s%dsC(ixCmin1:ixCmax1,&
               idir)*fE(ixCmin1:ixCmax1,idir)

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
    ! non-ideal electric field on cell edges
    double precision, dimension(ixImin1:ixImax1,7-2*ndim:3) :: E_resi, E_ambi
    ! total magnetic field at cell centers
    double precision                   :: Btot(ixImin1:ixImax1,1:ndim)
    integer                            :: hxCmin1,hxCmax1,ixCmin1,ixCmax1,&
       jxCmin1,jxCmax1,ixAmin1,ixAmax1,ixBmin1,ixBmax1
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm)

    if(B0field) then
      Btot(ixImin1:ixImax1,1:ndim)=wp(ixImin1:ixImax1,&
         mag(1:ndim))+block%B0(ixImin1:ixImax1,1:ndim,0)
    else
      Btot(ixImin1:ixImax1,1:ndim)=wp(ixImin1:ixImax1,mag(1:ndim))
    end if
    ECC=0.d0
    ! Calculate electric field at cell centers
    do idim1=1,ndim; do idim2=1,ndim; do idir=7-2*ndim,3
      if(lvc(idim1,idim2,idir)==1)then
         ECC(ixImin1:ixImax1,idir)=ECC(ixImin1:ixImax1,&
            idir)+Btot(ixImin1:ixImax1,idim1)*wp(ixImin1:ixImax1,mom(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixImin1:ixImax1,idir)=ECC(ixImin1:ixImax1,&
            idir)-Btot(ixImin1:ixImax1,idim1)*wp(ixImin1:ixImax1,mom(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(uawsom_eta/=zero) call get_resistive_electric_field(ixImin1,ixImax1,&
       ixOmin1,ixOmax1,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(uawsom_ambipolar_exp) call get_ambipolar_electric_field(ixImin1,ixImax1,&
       ixOmin1,ixOmax1,sCT%w,x,E_ambi)

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

            ! add resistive electric field at cell edges E=-vxB+eta J
            if(uawsom_eta/=zero) fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
               idir)+E_resi(ixCmin1:ixCmax1,idir)
            ! add ambipolar electric field
            if(uawsom_ambipolar_exp) fE(ixCmin1:ixCmax1,&
               idir)=fE(ixCmin1:ixCmax1,idir)+E_ambi(ixCmin1:ixCmax1,idir)

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
    double precision                   :: bfacetot(ixImin1:ixImax1,ndim)
    double precision                   :: btilL(ixImin1:ixImax1,ndim)
    double precision                   :: btilR(ixImin1:ixImax1,ndim)
    double precision                   :: cp(ixImin1:ixImax1,2)
    double precision                   :: cm(ixImin1:ixImax1,2)
    double precision                   :: circ(ixImin1:ixImax1,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixImin1:ixImax1,7-2*ndim:3) :: E_resi, E_ambi
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
    if(uawsom_eta/=zero) call get_resistive_electric_field(ixImin1,ixImax1,&
       ixOmin1,ixOmax1,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(uawsom_ambipolar_exp) call get_ambipolar_electric_field(ixImin1,ixImax1,&
       ixOmin1,ixOmax1,sCT%w,x,E_ambi)

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
      if(B0field) then
        bfacetot(ixImin1:ixImax1,idim1)=bfacesCT(ixImin1:ixImax1,&
           idim1)+block%B0(ixImin1:ixImax1,idim1,idim1)
        bfacetot(ixImin1:ixImax1,idim2)=bfacesCT(ixImin1:ixImax1,&
           idim2)+block%B0(ixImin1:ixImax1,idim2,idim2)
      else
        bfacetot(ixImin1:ixImax1,idim1)=bfacesCT(ixImin1:ixImax1,idim1)
        bfacetot(ixImin1:ixImax1,idim2)=bfacesCT(ixImin1:ixImax1,idim2)
      end if
      call reconstruct(ixImin1,ixImax1,ixCmin1,ixCmax1,idim2,&
         bfacetot(ixImin1:ixImax1,idim1),btilL(ixImin1:ixImax1,idim1),&
         btilR(ixImin1:ixImax1,idim1))

      call reconstruct(ixImin1,ixImax1,ixCmin1,ixCmax1,idim1,&
         bfacetot(ixImin1:ixImax1,idim2),btilL(ixImin1:ixImax1,idim2),&
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

      ! add resistive electric field at cell edges E=-vxB+eta J
      if(uawsom_eta/=zero) fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
         idir)+E_resi(ixCmin1:ixCmax1,idir)
      ! add ambipolar electric field
      if(uawsom_ambipolar_exp) fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
         idir)+E_ambi(ixCmin1:ixCmax1,idir)

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
      ! Time update cell-face magnetic field component
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
          ixCmax1=ixOmax1;
          ixCmin1=ixOmin1+kr(idir,1)-1;
          ixBmax1=ixCmax1-kr(idir,1)+1;
          ixBmin1=ixCmin1;
          ! current at transverse faces
          xs(ixBmin1:ixBmax1,:)=x(ixBmin1:ixBmax1,:)
          xs(ixBmin1:ixBmax1,idim2)=x(ixBmin1:ixBmax1,&
             idim2)+half*dx(ixBmin1:ixBmax1,idim2)
          call gradientx(wCTs(ixGslo1:ixGshi1,idim2),xs,ixGslo1,ixGshi1,&
             ixCmin1,ixCmax1,idim1,gradi,.true.)
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
    if(uawsom_eta>zero)then
      jce(ixImin1:ixImax1,:)=jce(ixImin1:ixImax1,:)*uawsom_eta
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

  !> get ambipolar electric field on cell edges
  subroutine get_ambipolar_electric_field(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,&
     fE)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)       :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out)      :: fE(ixImin1:ixImax1,7-2*ndim:3)

    double precision :: jxbxb(ixImin1:ixImax1,1:3)
    integer :: idir,ixAmin1,ixAmax1,ixCmin1,ixCmax1,ix1

    ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
    call uawsom_get_jxbxb(w,x,ixImin1,ixImax1,ixAmin1,ixAmax1,jxbxb)
    ! calcuate electric field on cell edges from cell centers
    do idir=7-2*ndim,3
      !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
      !jxbxb(ixA^S,i) = -(uawsom_eta_ambi/w(ixA^S, rho_)**2) * jxbxb(ixA^S,i)
      call multiplyAmbiCoef(ixImin1,ixImax1,ixAmin1,ixAmax1,&
         jxbxb(ixImin1:ixImax1,idir),w,x)
      ixCmax1=ixOmax1;
      ixCmin1=ixOmin1+kr(idir,1)-1;
      fE(ixCmin1:ixCmax1,idir)=0.d0
     do ix1=0,1
        if( ix1==1 .and. 1==idir ) cycle
        ixAmin1=ixCmin1+ix1;
        ixAmax1=ixCmax1+ix1;
        fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,&
           idir)+jxbxb(ixAmin1:ixAmax1,idir)
     end do
      fE(ixCmin1:ixCmax1,idir)=fE(ixCmin1:ixCmax1,idir)*0.25d0
    end do

  end subroutine get_ambipolar_electric_field

  !> calculate cell-center values from face-center values
  subroutine uawsom_face_to_center(ixOmin1,ixOmax1,s)
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

  end subroutine uawsom_face_to_center

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

end module mod_uawsom_phys
