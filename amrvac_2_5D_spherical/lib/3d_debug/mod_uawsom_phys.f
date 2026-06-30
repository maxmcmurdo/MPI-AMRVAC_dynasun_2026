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
  double precision, public, protected :: zeta0=5.0d0 
  
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
  logical, public, protected :: boundary_divbfix(2*3)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*3)=0

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

    subroutine mask_subroutine(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,res)
      use mod_global_parameters
      integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:nw)
      double precision, intent(inout) :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
    end subroutine mask_subroutine

    function fun_kin_en(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, inv_rho) result(ke)
      use mod_global_parameters, only: nw, ndim,block
      integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3, nw)
      double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
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
  
  public :: uawsom_clean_divb_multigrid
 

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

  subroutine uawsom_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim)
    use mod_global_parameters
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout)    :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux,1:ndim),  wnew(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    integer, intent(in)                :: idim
    integer                            :: hxOmin1,hxOmin2,hxOmin3,hxOmax1,&
       hxOmax2,hxOmax3, kxCmin1,kxCmin2,kxCmin3,kxCmax1,kxCmax2,kxCmax3, iw
    double precision                   :: inv_volume(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

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
    
    use mod_multigrid_coupling
   

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
    
    case ('multigrid')
       type_divb = divb_multigrid
       use_multigrid = .true.
       mg%operator_type = mg_laplacian
       phys_global_source_after => uawsom_clean_divb_multigrid
   
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

    
    ! clean initial divb
    if(clean_initial_divb) phys_clean_divb => uawsom_clean_divb_multigrid
   

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

    phys_te_images => uawsom_te_images

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


  subroutine uawsom_te_images
    use mod_global_parameters
    use mod_thermal_emission

    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_uawsom)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_uawsom)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_uawsom)
      case default
        call mpistop("Error in synthesize emission: Unknown convert_type")
      end select
  end subroutine uawsom_te_images


!!start th cond
  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  uawsom_sts_set_source_tc_mhd(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
     fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_mhd
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, igrid, nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_mhd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
       fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine uawsom_sts_set_source_tc_mhd

  function uawsom_get_tc_dt_mhd(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,&
     x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_mhd
 
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_mhd(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,x,tc_fl)
  end function uawsom_get_tc_dt_mhd

  subroutine uawsom_tc_handle_small_e(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, step)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    integer, intent(in)    :: step
    character(len=140) :: error_msg

    write(error_msg,"(a,i3)") "Thermal conduction step ", step
    call uawsom_handle_small_ei(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,e_,error_msg)
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
  subroutine set_equi_vars_grid_faces(igrid,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in) :: igrid, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    double precision :: delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: xC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim),xshift1,xshift2,xshift3
    integer :: idims, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ix, idims2

    if(slab_uniform)then
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1)=rnode(rpdx1_,&
         igrid)
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2)=rnode(rpdx2_,&
         igrid)
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3)=rnode(rpdx3_,&
         igrid)
    else
      ! for all non-cartesian and stretched cartesian coordinates
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)=ps(igrid)%dx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)
    endif

    do idims=1,ndim
      hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
      hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
      hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
      if(stagger_grid) then
        ! ct needs all transverse cells
        ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
        ixCmax2=ixOmax2+nghostcells-nghostcells*kr(idims,2)
        ixCmax3=ixOmax3+nghostcells-nghostcells*kr(idims,3)
        ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1)
        ixCmin2=hxOmin2-nghostcells+nghostcells*kr(idims,2)
        ixCmin3=hxOmin3-nghostcells+nghostcells*kr(idims,3);
      else
        ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
        ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=hxOmin1
        ixCmin2=hxOmin2;ixCmin3=hxOmin3;
      end if
      ! always xshift=0 or 1/2
      xshift1=half*(one-kr(1,idims));xshift2=half*(one-kr(2,idims))
      xshift3=half*(one-kr(3,idims));
      do idims2=1,ndim
        select case(idims2)
        case(1)
          do ix = ixCmin1,ixCmax1
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ix,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=x(ix,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,1)+(half-xshift1)*delx(ix,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,1)
          end do
        case(2)
          do ix = ixCmin2,ixCmax2
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ixCmin1:ixCmax1,ix,ixCmin3:ixCmax3,2)=x(ixCmin1:ixCmax1,ix,&
               ixCmin3:ixCmax3,2)+(half-xshift2)*delx(ixCmin1:ixCmax1,ix,&
               ixCmin3:ixCmax3,2)
          end do
        case(3)
          do ix = ixCmin3,ixCmax3
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ix,3)=x(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ix,3)+(half-xshift3)*delx(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ix,3)
          end do
        end select
      end do
      call usr_set_equi_vars(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,xC,&
         ps(igrid)%equi_vars(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:number_equi_vars,idims))
    end do

  end subroutine set_equi_vars_grid_faces

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid(igrid)
    use mod_global_parameters
    use mod_usr_methods
  
    integer, intent(in) :: igrid

    !values at the center
    call usr_set_equi_vars(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGlo1,&
       ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ps(igrid)%x,&
       ps(igrid)%equi_vars(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       1:number_equi_vars,0))

    !values at the interfaces
    call set_equi_vars_grid_faces(igrid,ps(igrid)%x,ixGlo1,ixGlo2,ixGlo3,&
       ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3)
 
  end subroutine set_equi_vars_grid

  ! w, wnew conserved, add splitted variables back to wnew
  function convert_vars_splitting(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
      nwc) result(wnew)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nwc
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision   :: wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        1:nwc)
    double precision   :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call  uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3))
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_) = rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(:)) =  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:))

    if (B0field) then
      ! add background magnetic field B0 to B
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0)
    else
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))
    end if

    if(uawsom_energy) then
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
      if(has_equi_pe0) then
        wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) = wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe0_,0)* inv_gamma_1
      end if
      if(B0field .and. .not. uawsom_internal_e) then
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)+0.5d0*sum(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,:,0)**2,dim=ndim+1) + sum(w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))*block%B0(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0),dim=ndim+1)
      end if
    end if

    !wnew(ixO^S,wminus_) = w(ixO^S,wminus_)
    !wnew(ixO^S,wplus_) = w(ixO^S,wplus_)

    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wkminus_)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wkplus_)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAminus_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAminus_)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAplus_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAplus_)


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
      miu0=4.d0*dpi ! G2,G2,G2 cm2,cm2,cm2 dyne-1,dyne-1,dyne-1
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

  if (mype==0 .and. mod(it,10000000)==0) then
    write(*,*) 'unit_magneticfield = ', unit_magneticfield
  end if 

  end subroutine uawsom_physical_units

  subroutine uawsom_check_w_origin(primitive,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    logical, intent(inout) :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    flag=.false.
    if(has_equi_rho0) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,0)
    else
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    endif
    where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_) = .true.

    if(uawsom_energy) then
      if(primitive) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_)
        if(has_equi_pe0) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe0_,0)
        endif
        where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) < small_pressure) flag(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_) = .true.
      else
        if(uawsom_internal_e) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_)
          if(has_equi_pe0) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)*inv_gamma_1
          endif
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_) = .true.
        else
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_)-uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
             ixOmax3)-uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
             ixOmax3)-uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
             ixOmax3)-uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) !Max: Added AW energy
          if(has_equi_pe0) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)*inv_gamma_1
          endif
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_) = .true.
        end if
      end if
    end if

  end subroutine uawsom_check_w_origin

  subroutine uawsom_check_w_hde(primitive,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    logical, intent(inout) :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    flag=.false.
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_) = .true.

    if(uawsom_energy) then
      if(primitive) then
        where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) = .true.
      else
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-uawsom_kin_en(w,ixImin1,ixImin2,&
           ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
           ixOmax2,ixOmax3)
        where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) = .true.
      end if
    end if

  end subroutine uawsom_check_w_hde

  !> Transform primitive variables into conservative ones
  subroutine uawsom_to_conserved_origin(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision :: inv_gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    !if (fix_small_values) then
    !  call uawsom_handle_small_values(.true., w, x, ixI^L, ixO^L, 'uawsom_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy and waves
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1)*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)+uawsom_mag_en(w, ixImin1,&
         ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
         ixOmax1,ixOmax2,ixOmax3)+uawsom_wk_en(w, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3)+uawsom_wA_en(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) !Max: Added AW energy 
      if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         eaux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         paux_)*inv_gamma_1
    end if
 
    if(uawsom_boris_simplification) then
      inv_gamma2=1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = inv_gamma2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(idir))
      end do
    else
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idir))
      end do
    end if
  end subroutine uawsom_to_conserved_origin

  !> Transform primitive variables into conservative ones
  subroutine uawsom_to_conserved_hde(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    integer                         :: idir

    !if (fix_small_values) then
    !  call uawsom_handle_small_values(.true., w, x, ixI^L, ixO^L, 'uawsom_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy 
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1)*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) !Max: hde so we do not consider contribution from magnetic waves AW and kink
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idir))
    end do
  end subroutine uawsom_to_conserved_hde

  !> Transform primitive variables into conservative ones
  subroutine uawsom_to_conserved_inte(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision :: inv_gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    !if (fix_small_values) then
    !  call uawsom_handle_small_values(.true., w, x, ixI^L, ixO^L, 'uawsom_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy ! Max: No kinetic energy or magnetic energy here?
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)*inv_gamma_1
    end if

    if(uawsom_boris_simplification) then
      inv_gamma2=1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = inv_gamma2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(idir))
      end do
    else
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idir))
      end do
    end if
  end subroutine uawsom_to_conserved_inte

  !> Transform primitive variables into conservative ones
  subroutine uawsom_to_conserved_split_rho(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    integer                         :: idir
    double precision                :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    !if (fix_small_values) then
    !  call uawsom_handle_small_values(.true., w, x, ixI^L, ixO^L, 'uawsom_to_conserved')
    !end if

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,b0i)
    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(uawsom_energy) then
      if(uawsom_internal_e) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           p_)*inv_gamma_1
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1)*rho(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)+uawsom_mag_en(w, ixImin1,ixImin2,&
           ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
           ixOmax2,ixOmax3)+uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3)+uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) !Max: Added AW energy
        if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,eaux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,paux_)*inv_gamma_1
      end if
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir)) = rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))
    end do
  end subroutine uawsom_to_conserved_split_rho

  !> Transform conservative variables into primitive ones
  subroutine uawsom_to_primitive_origin(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call uawsom_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3, 'uawsom_to_primitive_origin')
    end if

    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)

    ! Calculate pressure = (gamma-1) * (e-ek-eb) ! Max: In Norberts work it was (gamma-1) * (e-ek-eb - E_kink) since he included kink energy, 
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)-uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         inv_rho)-uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3)-uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3)-uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)) !Max: AW energy added
      if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         paux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         eaux_)*gamma_1
    end if

    if(uawsom_boris_simplification) then
      gamma2=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*inv_rho*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*inv_rho
      end do
    end if

  end subroutine uawsom_to_primitive_origin

  !> Transform conservative variables into primitive ones
  subroutine uawsom_to_primitive_hde(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call uawsom_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3, 'uawsom_to_primitive_hde')
    end if

    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)

    ! Calculate pressure = (gamma-1) * (e-ek)
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)-uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,inv_rho))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*inv_rho
    end do

  end subroutine uawsom_to_primitive_hde

  !> Transform conservative variables into primitive ones
  subroutine uawsom_to_primitive_inte(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call uawsom_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3, 'uawsom_to_primitive_inte')
    end if

    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)

    ! Calculate pressure = (gamma-1) * e_internal
    if(uawsom_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)*gamma_1
    end if

    if(uawsom_boris_simplification) then
      gamma2=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*inv_rho*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*inv_rho
      end do
    end if

  end subroutine uawsom_to_primitive_inte

  !> Transform conservative variables into primitive ones
  subroutine uawsom_to_primitive_split_rho(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision                :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call uawsom_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3, 'uawsom_to_primitive_split_rho')
    end if

    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,equi_rho0_,b0i))

    ! Calculate pressure = (gamma-1) * (e-ek-eb)
    if(uawsom_energy) then
      if(uawsom_internal_e) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           p_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)*gamma_1
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)-uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           inv_rho)-uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3)-uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3)-uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)) !Max: AW energy added
        if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,paux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,eaux_)*gamma_1
      end if
    end if
    
    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*inv_rho
    end do

  end subroutine uawsom_to_primitive_split_rho

  !> Transform internal energy to total energy
  subroutine uawsom_ei_to_e(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    ! Calculate total energy from internal, kinetic and magnetic energy
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)+uawsom_kin_en(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3)+uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3)+uawsom_wk_en(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3)+uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3) !Max: Added AW energy density

  end subroutine uawsom_ei_to_e

  !> Transform total energy to internal energy
  subroutine uawsom_e_to_ei(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    ! Calculate ei = e - ek - eb
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)-uawsom_kin_en(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3)-uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3)-uawsom_wk_en(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3)-uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3) !Max: Added AW energy

    if(fix_small_values) then
      call uawsom_handle_small_ei(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,e_,&
         'uawsom_e_to_ei')
    end if

  end subroutine uawsom_e_to_ei

  !> Transform hydrodynamic energy to internal energy
  subroutine uawsom_e_to_ei_hde(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    ! Calculate ei = e - ek
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)-uawsom_kin_en(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3)

    if(fix_small_values) then
      call uawsom_handle_small_ei(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,e_,&
         'uawsom_e_to_ei_hde')
    end if

  end subroutine uawsom_e_to_ei_hde

  !> Update eaux and transform internal energy to total energy
  subroutine uawsom_ei_to_e_aux(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
 
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,eaux_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)
    ! Calculate total energy from internal, kinetic and magnetic energy ! Max: Do we need Kink and AW energy here?
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)+uawsom_kin_en(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3)+uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3)

  end subroutine uawsom_ei_to_e_aux

  !> Transform total energy to internal energy via eaux as internal energy
  subroutine uawsom_e_to_ei_aux(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,eaux_)

  end subroutine uawsom_e_to_ei_aux

  subroutine uawsom_handle_small_values_origin(primitive, w, x, ixImin1,&
     ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
     ixOmax2,ixOmax3, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision :: tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    if(small_values_method == "ignore") return

    call phys_check_w(primitive, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_rho0) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_) = small_density-block%equi_vars(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,0)
        else
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_) = small_density
        endif
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom(idir)) = 0.0d0
          end if
        end do

        if(uawsom_energy) then
          if(primitive) then
           if(has_equi_pe0) then
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = small_pressure - &
               block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               equi_pe0_,0)
           else
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = small_pressure
           endif
           where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              e_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              p_) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          else
            ! conserved
            if(has_equi_pe0) then
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = small_e - block%equi_vars(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)*inv_gamma_1
            else
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = small_e
            endif
            if(uawsom_internal_e) then
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_))
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_)=tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
              end where
            else
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_))
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)+uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,&
                   ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                   ixOmax2,ixOmax3)+uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,&
                   ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                   ixOmax2,ixOmax3)+uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,&
                   ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                   ixOmax2,ixOmax3)+uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,&
                   ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                   ixOmax2,ixOmax3) !Max: Added AW energy
             end where
              if(uawsom_solve_eaux) then
                where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_))
                  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     eaux_)=tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)
                end where
              end if
            end if
          end if
        end if
      case ("average")
        if(primitive)then
          call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
              flag, rho_)
          if(uawsom_energy) then
            call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
                flag, p_)
          end if
        else
          ! do averaging of density
          call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
              flag, rho_)
          if(uawsom_energy) then
             ! do averaging of internal energy
            if(.not.uawsom_internal_e) then
              w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 e_)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 e_)-uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3)-uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3)-uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3)-uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3) !Max: Added AW energy 
            end if
            call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
                flag, e_)
             ! convert back
            if(.not.uawsom_internal_e) then
              w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 e_)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 e_)+uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3)+uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3)+uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3)+uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3) !Max: Added AW energy
            end if
            ! eaux
            if(uawsom_solve_eaux) then
              call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
                 ixOmax3, w, x, flag, paux_)
            end if
          end if
        endif
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek-eb)
          if(uawsom_energy) then
            if(uawsom_internal_e) then
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 p_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_)*gamma_1
            else
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,e_)-uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,&
                 ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                 ixOmax2,ixOmax3)-uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,&
                 ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                 ixOmax2,ixOmax3)-uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,&
                 ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                 ixOmax2,ixOmax3)-uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,&
                 ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                 ixOmax2,ixOmax3)) !Max: AW energy added
              if(uawsom_solve_eaux) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,paux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,eaux_)*gamma_1
            end if
          end if
          ! Convert momentum to velocity
          if(has_equi_rho0) then
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,0)
          else
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_)
          endif
          do idir = 1, ndir
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mom(idir))/tmp2(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do
        end if
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if
  end subroutine uawsom_handle_small_values_origin

  !> Calculate v vector
  subroutine uawsom_get_v_origin(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndir)

    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: idir

    call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1.d0/rho(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    ! Convert momentum to velocity
    do idir = 1, ndir
       v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

  end subroutine uawsom_get_v_origin

  !> Calculate v vector
  subroutine uawsom_get_v_boris(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndir)

    double precision              :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer :: idir

    call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1.d0/rho(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    gamma2=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mag(:))**2,dim=ndim+1)*rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*inv_squared_c)
    ! Convert momentum to velocity
    do idir = 1, ndir
       v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*gamma2
    end do

  end subroutine uawsom_get_v_boris

  !> Calculate v component
  subroutine uawsom_get_v_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision              :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)

    if(uawsom_boris_simplification) then
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idim)) / rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) /(1.d0+sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))**2,&
         dim=ndim+1)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*inv_squared_c)
    else
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idim)) / rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if

  end subroutine uawsom_get_v_idim

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine uawsom_get_cmax_origin(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call uawsom_get_csound(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,cmax)
    call uawsom_get_v_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,vel)

    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=abs(vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))+cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine uawsom_get_cmax_origin

  subroutine uawsom_get_a2max(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a2max)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(inout) :: a2max(ndim)
    double precision :: a2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ndim,nw)
    integer :: gxOmin1,gxOmin2,gxOmin3,gxOmax1,gxOmax2,gxOmax3,hxOmin1,hxOmin2,&
       hxOmin3,hxOmax1,hxOmax2,hxOmax3,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,&
       jxOmax3,kxOmin1,kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3,i,j

    a2=zero
    do i = 1,ndim
      !> 4th order
      hxOmin1=ixOmin1-kr(i,1);hxOmin2=ixOmin2-kr(i,2);hxOmin3=ixOmin3-kr(i,3)
      hxOmax1=ixOmax1-kr(i,1);hxOmax2=ixOmax2-kr(i,2);hxOmax3=ixOmax3-kr(i,3);
      gxOmin1=hxOmin1-kr(i,1);gxOmin2=hxOmin2-kr(i,2);gxOmin3=hxOmin3-kr(i,3)
      gxOmax1=hxOmax1-kr(i,1);gxOmax2=hxOmax2-kr(i,2);gxOmax3=hxOmax3-kr(i,3);
      jxOmin1=ixOmin1+kr(i,1);jxOmin2=ixOmin2+kr(i,2);jxOmin3=ixOmin3+kr(i,3)
      jxOmax1=ixOmax1+kr(i,1);jxOmax2=ixOmax2+kr(i,2);jxOmax3=ixOmax3+kr(i,3);
      kxOmin1=jxOmin1+kr(i,1);kxOmin2=jxOmin2+kr(i,2);kxOmin3=jxOmin3+kr(i,3)
      kxOmax1=jxOmax1+kr(i,1);kxOmax2=jxOmax2+kr(i,2);kxOmax3=jxOmax3+kr(i,3);
      a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,&
         1:nw)=abs(-w(kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3,&
         1:nw)+16.d0*w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
         1:nw)-30.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nw)+16.d0*w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
         1:nw)-w(gxOmin1:gxOmax1,gxOmin2:gxOmax2,gxOmin3:gxOmax3,1:nw))
      a2max(i)=maxval(a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,&
         1:nw))/12.d0/dxlevel(i)**2
    end do
  end subroutine uawsom_get_a2max

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine uawsom_get_tcutoff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,Tco_local,&
     Tmax_local)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(out) :: Tco_local,Tmax_local

    double precision, parameter :: trac_delta=0.25d0
    double precision :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       Te(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),lts(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir) :: bunitvec
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) :: gradT
    double precision :: Bdir(ndim)
    double precision :: ltrc,ltrp,altr(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer :: idims,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,ixA1,ixA2,ixA3,ixB1,ixB2,ixB3
    integer :: jxPmin1,jxPmin2,jxPmin3,jxPmax1,jxPmax2,jxPmax3,hxPmin1,hxPmin2,&
       hxPmin3,hxPmax1,hxPmax2,hxPmax3,ixPmin1,ixPmin2,ixPmin3,ixPmax1,ixPmax2,&
       ixPmax3
    logical :: lrlt(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    ! reuse lts as rhoc
    call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,lts)
    if(uawsom_internal_e) then
      tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,e_)*gamma_1
    else
      call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,tmp1)
    end if
    Te(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=tmp1(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)/lts(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    Tco_local=zero
    Tmax_local=maxval(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

    
    
    select case(uawsom_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      block%wextra(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         Tcoff_)=2.5d5/unit_temperature
    case(1,4,6)
      ! temperature gradient at cell centers
      do idims=1,ndim
        call gradient(Te,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,tmp1)
        gradT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do
      ! B vector
      if(B0field) then
        bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           :)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           :,0)
      else
        bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           :)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mag(:))
      end if
      if(uawsom_trac_type .gt. 1) then
        ! B direction at cell center
        Bdir=zero
        do ixA1=0,1
    do ixA2=0,1
    do ixA3=0,1
          ixB1=(ixOmin1+ixOmax1-1)/2+ixA1;ixB2=(ixOmin2+ixOmax2-1)/2+ixA2
          ixB3=(ixOmin3+ixOmax3-1)/2+ixA3;
          Bdir(1:ndim)=Bdir(1:ndim)+bunitvec(ixB1,ixB2,ixB3,1:ndim)
        end do
    end do
    end do
        if(sum(Bdir(:)**2) .gt. zero) then
          Bdir(1:ndim)=Bdir(1:ndim)/dsqrt(sum(Bdir(:)**2))
        end if
        block%special_values(3:ndim+2)=Bdir(1:ndim)
      end if
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=dsqrt(sum(bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)**2,dim=ndim+1))
      where(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/=0.d0)
        tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=1.d0/tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      elsewhere
        tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=bigdouble
      end where
      ! b unit vector: magnetic field direction vector
      do idims=1,ndim
        bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)=bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do
      ! temperature length scale inversed
      lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=abs(sum(gradT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndim)*bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndim),dim=ndim+1))/Te(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      ! fraction of cells size to temperature length scale
      if(slab_uniform) then
        lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=minval(dxlevel)*lts(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      else
        lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=minval(block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,:),dim=ndim+1)*lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
      lrlt=.false.
      where(lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) > trac_delta)
        lrlt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=.true.
      end where
      if(any(lrlt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))) then
        block%special_values(1)=maxval(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3), mask=lrlt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
      else
        block%special_values(1)=zero
      end if
      block%special_values(2)=Tmax_local
    case(2)
      !> iijima et al. 2021, LTRAC method
      ltrc=1.5d0
      ltrp=4.d0
      ixPmin1=ixOmin1-1;ixPmin2=ixOmin2-1;ixPmin3=ixOmin3-1;ixPmax1=ixOmax1+1
      ixPmax2=ixOmax2+1;ixPmax3=ixOmax3+1;
      ! temperature gradient at cell centers
      do idims=1,ndim
        call gradient(Te,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixPmin1,ixPmin2,ixPmin3,ixPmax1,ixPmax2,ixPmax3,idims,tmp1)
        gradT(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           idims)=tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)
      end do
      ! B vector
      if(B0field) then
        bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           :)=w(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           iw_mag(:))+block%B0(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           :,0)
      else
        bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           :)=w(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,iw_mag(:))
      end if
      tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3)=dsqrt(sum(bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3,:)**2,dim=ndim+1))
      where(tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)/=0.d0)
        tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)=1.d0/tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)
      elsewhere
        tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)=bigdouble
      end where
      ! b unit vector: magnetic field direction vector
      do idims=1,ndim
        bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           idims)=bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           idims)*tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)
      end do
      ! temperature length scale inversed
      lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3)=abs(sum(gradT(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3,1:ndim)*bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3,1:ndim),dim=ndim+1))/Te(ixPmin1:ixPmax1,&
         ixPmin2:ixPmax2,ixPmin3:ixPmax3)
      ! fraction of cells size to temperature length scale
      if(slab_uniform) then
        lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)=minval(dxlevel)*lts(ixPmin1:ixPmax1,&
           ixPmin2:ixPmax2,ixPmin3:ixPmax3)
      else
        lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)=minval(block%ds(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3,:),dim=ndim+1)*lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)
      end if
      lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)=max(one,&
          (exp(lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3))/ltrc)**ltrp)
  
      altr=zero
      do idims=1,ndim
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
        hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
        jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
        jxOmin3=ixOmin3+kr(idims,3);jxOmax1=ixOmax1+kr(idims,1)
        jxOmax2=ixOmax2+kr(idims,2);jxOmax3=ixOmax3+kr(idims,3);
        altr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=altr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+0.25d0*(lts(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3)+two*lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+lts(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3))*bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idims)**2
      end do
      block%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         Tcoff_)=Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*altr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**0.4d0
      ! need one ghost layer for thermal conductivity
      block%wextra(ixOmin1-1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
         Tcoff_)=block%wextra(ixOmin1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,Tcoff_) 
    block%wextra(ixPmin1:ixPmax1,ixOmin2-1,ixPmin3:ixPmax3,&
       Tcoff_)=block%wextra(ixPmin1:ixPmax1,ixOmin2,ixPmin3:ixPmax3,Tcoff_) 
    block%wextra(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixOmin3-1,&
       Tcoff_)=block%wextra(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixOmin3,Tcoff_) 
      block%wextra(ixOmax1+1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
         Tcoff_)=block%wextra(ixOmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,Tcoff_) 
    block%wextra(ixPmin1:ixPmax1,ixOmax2+1,ixPmin3:ixPmax3,&
       Tcoff_)=block%wextra(ixPmin1:ixPmax1,ixOmax2,ixPmin3:ixPmax3,Tcoff_) 
    block%wextra(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixOmax3+1,&
       Tcoff_)=block%wextra(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixOmax3,Tcoff_) 
    case(3,5)
      !> do nothing here
    case default
      call mpistop("unknown uawsom_trac_type")
    end select
   
  end subroutine uawsom_get_tcutoff

  !> get H speed for H-correction to fix the carbuncle problem at grid-aligned shock front
  subroutine uawsom_get_H_speed(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     Hspeed)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wprim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    double precision :: csound(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ndim),tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3, ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, ixAmin1,ixAmin2,ixAmin3,&
       ixAmax1,ixAmax2,ixAmax3, id, ix1,ix2,ix3

    Hspeed=0.d0
    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    do id=1,ndim
      call uawsom_get_csound_prim(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,id,&
         tmp)
      csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)=tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)
    end do
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
    ixCmin1=ixOmin1+kr(idim,1)-1;ixCmin2=ixOmin2+kr(idim,2)-1
    ixCmin3=ixOmin3+kr(idim,3)-1;
    jxCmax1=ixCmax1+kr(idim,1);jxCmax2=ixCmax2+kr(idim,2)
    jxCmax3=ixCmax3+kr(idim,3);
    jxCmin1=ixCmin1+kr(idim,1);jxCmin2=ixCmin2+kr(idim,2)
    jxCmin3=ixCmin3+kr(idim,3);
    Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       1)=0.5d0*abs(wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
       mom(idim))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
       idim)-wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       mom(idim))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       idim))

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=ixCmax1+kr(id,1);ixAmax2=ixCmax2+kr(id,2)
      ixAmax3=ixCmax3+kr(id,3);
      ixAmin1=ixCmin1+kr(id,1);ixAmin2=ixCmin2+kr(id,2)
      ixAmin3=ixCmin3+kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)-wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(id))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,id)))
      ixAmax1=ixCmax1-kr(id,1);ixAmax2=ixCmax2-kr(id,2)
      ixAmax3=ixCmax3-kr(id,3);
      ixAmin1=ixCmin1-kr(id,1);ixAmin2=ixCmin2-kr(id,2)
      ixAmin3=ixCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(id))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         id)-wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,id)))
    end do

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=jxCmax1+kr(id,1);ixAmax2=jxCmax2+kr(id,2)
      ixAmax3=jxCmax3+kr(id,3);
      ixAmin1=jxCmin1+kr(id,1);ixAmin2=jxCmin2+kr(id,2)
      ixAmin3=jxCmin3+kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)-wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         mom(id))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,id)))
      ixAmax1=jxCmax1-kr(id,1);ixAmax2=jxCmax2-kr(id,2)
      ixAmax3=jxCmax3-kr(id,3);
      ixAmin1=jxCmin1-kr(id,1);ixAmin2=jxCmin2-kr(id,2)
      ixAmin3=jxCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         mom(id))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         id)-wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,id)))
    end do

  end subroutine uawsom_get_H_speed

  !> Estimating bounds for the minimum and maximum signal velocities without split
  subroutine uawsom_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     idim,Hspeed,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    double precision :: wmean(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       nw)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix1,ix2,ix3

    select case (boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_))
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_))
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call uawsom_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundL)
      call uawsom_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundR)
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+0.5d0*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim))-wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim)))**2
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=abs(umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (2)
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux))
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim))/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)
      call uawsom_get_csound(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call uawsom_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundL)
      call uawsom_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundR)
      csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=max(csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))-csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    end select

  end subroutine uawsom_get_cbounds

  !> Estimating bounds for the minimum and maximum signal velocities with rho split
  subroutine uawsom_get_cbounds_split_rho(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3,idim,Hspeed,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    double precision :: wmean(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       nw)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix1,ix2,ix3
    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    select case (boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,b0i))
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,b0i))
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call uawsom_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundL)
      call uawsom_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundR)
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+0.5d0*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim))-wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim)))**2
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=abs(umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (2)
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux))
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim))/(wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,b0i))
      call uawsom_get_csound(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call uawsom_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundL)
      call uawsom_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundR)
      csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=max(csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))-csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    end select

  end subroutine uawsom_get_cbounds_split_rho

  !> prepare velocities for ct methods
  subroutine uawsom_get_ct_velocity(vcts,wLp,wRp,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(in), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    type(ct_velocity), intent(inout):: vcts

    integer                         :: idimE,idimN

    ! calculate velocities related to different UCT schemes
    select case(type_ct)
    case('average')
    case('uct_contact')
      if(.not.allocated(vcts%vnorm)) allocate(vcts%vnorm(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1:ndim))
      ! get average normal velocity at cell faces
      vcts%vnorm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)=0.5d0*(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim)))
    case('uct_hll')
      if(.not.allocated(vcts%vbarC)) then
        allocate(vcts%vbarC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1:ndir,2),vcts%vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndir,2),vcts%vbarRC(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,1:ndir,2))
        allocate(vcts%cbarmin(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1:ndim),vcts%cbarmax(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndim))
      end if
      ! Store magnitude of characteristics
      if(present(cmin)) then
        vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)=max(-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           zero)
        vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)=max( cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           zero)
      else
        vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)=max( cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           zero)
        vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)=vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)
      end if

      idimN=mod(idim,ndir)+1 ! 'Next' direction
      idimE=mod(idim+1,ndir)+1 ! Electric field direction
      ! Store velocities
      vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         1)=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idimN))
      vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         1)=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idimN))
      vcts%vbarC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         1)=(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,1) +vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,1))/(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)+vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim))

      vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         2)=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idimE))
      vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         2)=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idimE))
      vcts%vbarC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         2)=(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,2) +vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,1))/(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)+vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim))
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine uawsom_get_ct_velocity

  !> Calculate fast magnetosonic wave speed
  subroutine uawsom_get_csound(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), b2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        kmax
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if(has_equi_rho0) then
      inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,b0i))
    else
      inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)
    endif

    call uawsom_get_csound2(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound)

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = uawsom_mag_en_all(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)

    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2-4.0d0*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * uawsom_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       idim)**2 * inv_rho

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    if (.not. uawsom_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = sqrt(half*(cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),dxlevel(3),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3))),&
           uawsom_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))*inv_rho*kmax)
    end if

  end subroutine uawsom_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine uawsom_get_csound_prim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), b2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        kmax
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)


    if(uawsom_energy) then
      if(has_equi_pe0) then
        csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) + block%equi_vars(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,b0i)
      else
        csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_)
      endif
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=uawsom_gamma*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*inv_rho
    else
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=uawsom_gamma*uawsom_adiab*tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)**gamma_1
    end if

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)        = &
       uawsom_mag_en_all(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2-4.0d0*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * uawsom_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       idim)**2 * inv_rho

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    if (.not. uawsom_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = sqrt(half*(cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),dxlevel(3),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3))),&
           uawsom_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))*inv_rho*kmax)
    end if

  end subroutine uawsom_get_csound_prim

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine uawsom_get_pthermal_origin(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                      :: iw, ix1,ix2,ix3

    if(uawsom_energy) then
      if(uawsom_internal_e) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=gamma_1*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_)
      else
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_)- uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,&
           ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3)- uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3)- uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3)- uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)) !Max: AW energy added
      end if
      if(has_equi_pe0) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe0_,b0i)
      endif
    else
      call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=uawsom_adiab*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**uawsom_gamma
    end if

    if (fix_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
            pth(ix1,ix2,ix3)=small_pressure
         end if
      enddo
      enddo
      enddo
    end if

    if (check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix1,ix2,ix3),&
              " encountered when call uawsom_get_pthermal"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix1,ix2,ix3,:)
           write(*,*) "Cell number: ", ix1,ix2,ix3
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix1,ix2,ix3,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) sqrt(pth(ix1,ix2,ix3)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      enddo
      enddo
      enddo
    end if

  end subroutine uawsom_get_pthermal_origin

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine uawsom_get_pthermal_hde(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                      :: iw, ix1,ix2,ix3

    if(uawsom_energy) then
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,e_)-uawsom_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)) !Max: No wave energy contribution since we are considering only hydrodynamic 
    else
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=uawsom_adiab*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)**uawsom_gamma
    end if

    if (fix_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
            pth(ix1,ix2,ix3)=small_pressure
         end if
      enddo
      enddo
      enddo
    end if

    if (check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix1,ix2,ix3),&
              " encountered when call uawsom_get_pthermal_hde"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix1,ix2,ix3,:)
           write(*,*) "Cell number: ", ix1,ix2,ix3
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix1,ix2,ix3,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) sqrt(pth(ix1,ix2,ix3)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      enddo
      enddo
      enddo
    end if

  end subroutine uawsom_get_pthermal_hde

  !> Calculate temperature=p/rho when in e_ the internal energy is stored
  subroutine uawsom_get_temperature_from_eint(w, x, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = gamma_1 * &
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        e_) /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
  end subroutine uawsom_get_temperature_from_eint

  !> Calculate temperature=p/rho when in e_ the total energy is stored
  !> this does not check the values of uawsom_energy and uawsom_internal_e,
  !>  uawsom_energy = .true. and uawsom_internal_e = .false.
  !> also check small_values is avoided
  subroutine uawsom_get_temperature_from_etot(w, x, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(gamma_1*(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)- uawsom_kin_en(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)- uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3)- uawsom_wk_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3))- uawsom_wA_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) !Max: added AW energy here. What is res? Res seems to be Temperature? res = P/rho = T
  end subroutine uawsom_get_temperature_from_etot

  !> Calculate temperature from hydrodynamic energy
  subroutine uawsom_get_temperature_from_hde(w, x, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=gamma_1*(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)- uawsom_kin_en(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)
  end subroutine uawsom_get_temperature_from_hde

  subroutine uawsom_get_temperature_from_eint_with_equi(w, x, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = (gamma_1 * &
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        e_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_pe0_,b0i)) /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_) +block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_rho0_,b0i))
  end subroutine uawsom_get_temperature_from_eint_with_equi

  subroutine uawsom_get_temperature_equi(w,x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= &
       block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_pe0_,b0i)/block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,equi_rho0_,b0i)
  end subroutine uawsom_get_temperature_equi

  subroutine uawsom_get_rho_equi(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_rho0_,b0i)
  end subroutine uawsom_get_rho_equi

  subroutine uawsom_get_pe_equi(w,x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_pe0_,b0i)
  end subroutine uawsom_get_pe_equi

  subroutine uawsom_get_temperature_from_etot_with_equi(w, x, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(gamma_1*(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)- uawsom_kin_en(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)- uawsom_mag_en(w,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3)) +  block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,equi_pe0_,b0i))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_) +block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,equi_rho0_,b0i))
            
  end subroutine uawsom_get_temperature_from_etot_with_equi

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine uawsom_get_csound2(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision    :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    
    call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    if(uawsom_energy) then
      call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=uawsom_gamma*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    else
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=uawsom_gamma*uawsom_adiab*rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)**gamma_1
    end if
  end subroutine uawsom_get_csound2

  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine uawsom_get_p_total(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,p)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: p(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,p)

    p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = p(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) + 0.5d0 * sum(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(:))**2, dim=ndim+1)

  end subroutine uawsom_get_p_total

  !> Calculate fluxes within ixO^L without any splitting
  subroutine uawsom_get_flux(wC,w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nwflux)

    double precision             :: ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision             :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), zeta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision             :: vHall(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:,:,:,:) :: Jambi, btot
    double precision, allocatable, dimension(:,:,:) :: tmp2, tmp3

    ! Get flux of density
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

    call get_zeta(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,zeta)

    if(uawsom_energy) then 
      ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)+0.5d0*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)+(zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+1.d0)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,wkplus_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,wkminus_))/4.d0+(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,wAplus_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,wAminus_))/2.0d0 !Max: Added AW pressure
    else
      ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=uawsom_adiab*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)**uawsom_gamma+0.5d0*sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)
    end if

    if (uawsom_Hall) then
      call uawsom_getv_Hall(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vHall)
    end if

    ! Get flux of tracer
    do iw=1,uawsom_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))
      end if
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k) + Q-wkplus + Q+wkminus ! Max: + Z-wAplus + Z+wAminus for AWs TVD 2024 Eqn. 78
    if(uawsom_energy) then
      if (uawsom_internal_e) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
         ! Max: Adding AW to flux of energy. Note: mu (magnetic permeability) = 1 in cgs so does not appear here           
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))*(wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:)),&
           dim=ndim+1)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           wkminus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)*(zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+1.d0)/2.0d0)**0.5d0)+w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,wkplus_)*(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim)) - w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idim))/(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)*(zeta(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)+1.d0)/2.0d0)**0.5d0)+&
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           wAminus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)**0.5d0))+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           wAplus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)**0.5d0))


        if(uawsom_solve_eaux) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,eaux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,eaux_)
        if(uawsom_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim) * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(:))**2,dim=ndim+1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idim)) * sum(vHall(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)),dim=ndim+1)
        end if
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (uawsom_glm) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              psi_)
        else
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))
        if (uawsom_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir)) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir)) - vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idim)) + vHall(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)*w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))
        end if
      end if
    end do

    ! compute flux of wkminus Q+wkminus
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mag(idim))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)*(zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+1.d0)/2.d0)**0.5d0)

    ! compute flux of wkplus Q-wkplus 
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim)) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mag(idim))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)*(zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+1.d0)/2.d0)**0.5d0)
                              
                      
    ! Max: compute flux of wAminus  Z+wAminus
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAminus_)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAminus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mag(idim))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)**0.5d0))

    ! Max: compute flux of wAplus Z-wAplus
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAplus_)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAplus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim)) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mag(idim))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)**0.5d0))

    if (uawsom_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(uawsom_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * uawsom_eta_ambi/rho**2
      allocate(Jambi(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:3))
      call uawsom_get_Jambi(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Jambi)
      allocate(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3))
      btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1:3))
      allocate(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
         tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      !tmp2 = Btot^2
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)*btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
        case(2)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
        case(3)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
      endselect

      if(uawsom_energy .and. .not. uawsom_internal_e) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

  end subroutine uawsom_get_flux

  !> Calculate fluxes within ixO^L with possible splitting
  subroutine uawsom_get_flux_split(wC,w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nwflux)

    double precision             :: pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndir)
    double precision             :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), zeta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision             :: vHall(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:,:,:,:) :: Jambi, btot
    double precision, allocatable, dimension(:,:,:) :: tmp2, tmp3
    double precision :: tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)


    call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
    ! Get flux of density
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    ! pgas is time dependent only
    if(uawsom_energy) then
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
    else
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=uawsom_adiab*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**uawsom_gamma
      if(has_equi_pe0) then
        pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe0_,b0i)
      endif
    end if

    if (uawsom_Hall) then
      call uawsom_getv_Hall(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vHall)
    end if

    if(B0field) then
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,idim)
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:,idim),dim=ndim+1)
    else
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))
    end if

    call get_zeta(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,zeta)

    ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+0.5d0*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)+(zeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)+1.d0)*(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,wkplus_) + w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,wkminus_))/4.d0+(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAplus_) + w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAminus_))/2.d0 !Max: added AW pressure

    ! Get flux of tracer
    do iw=1,uawsom_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    if(B0field) then
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idim)-block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir,idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idim))
        else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)-block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir,idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idim))
        end if
      end do
    else
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idim)
        else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)
        end if
      end do
    end if

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k) ! Max: Eqn 78 TVD 2024, inside divergence ie div((rho*v^2/2 + p/gamma-1 + B^2/2mu)v - B*v.B/mu)
    !                                       ! Max: Norbert has also added wplus and wminus (Kink waves) as well to this equation-terms in Eqn 78.
    ! Max: This subtraction in the flux calculation a few lines down is because wkplus carries wave energy into the Sun
    if(uawsom_energy) then
      if (uawsom_internal_e) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idim))*(wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_)+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3))-B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,idim)*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(:)),dim=ndim+1)+w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,wkminus_)*B(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)/((w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)*(zeta(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3)+1.d0)/2.d0)**0.5d0) &
            -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            wkplus_)*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            idim)/((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_)*(zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+1.d0)/2.d0)**0.5d0) +w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAminus_)*B(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)/(w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)**0.5d0) -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAplus_)*B(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)/(w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)**0.5d0)

        if(uawsom_solve_eaux) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,eaux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,eaux_)
        if (uawsom_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
           if (uawsom_etah>zero) then
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idim) * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  mag(:))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:),&
                 dim=ndim+1) - B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,idim) * sum(vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)),dim=ndim+1)
           end if
        end if
      end if
      if(has_equi_pe0) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=  f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)) * block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe0_,idim) * inv_gamma_1
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (uawsom_glm) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              psi_)
        else
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)-B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))

        if (uawsom_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (uawsom_etah>zero) then
            f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(idir)) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(idir)) - vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,idir)*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,idim) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,idim)*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,idir)
          end if
        end if

      end if
    end do
     
    ! compute flux of wkminus
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim)) + B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idim)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)*(zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+1.d0)/2.d0)**0.5d0)

    ! compute flux of wkplus
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim)) - B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idim)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)*(zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+1.d0)/2.d0)**0.5d0)

    ! Max: compute flux of wAminus
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAminus_)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAminus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim)) + B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idim)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)**0.5d0))

    ! Max: compute flux of wAplus
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAplus_)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAplus_)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim)) - B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idim)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)**0.5d0))

    if (uawsom_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(uawsom_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * uawsom_eta_ambi/rho**2
      allocate(Jambi(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:3))
      call uawsom_get_Jambi(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Jambi)
      allocate(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3))
      if(B0field) then
        do idir=1,3
          btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir)) + block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir,idim)
        enddo
      else
        btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1:3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1:3))
      endif
      allocate(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
         tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      !tmp2 = Btot^2
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)*btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          if(B0field) tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
        case(2)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
          if(B0field) tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
        case(3)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          if(B0field) tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
      endselect

      if(uawsom_energy .and. .not. uawsom_internal_e) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) +  tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *  tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

  end subroutine uawsom_get_flux_split

  !> Source terms J.E in internal energy.
  !> For the ambipolar term E = ambiCoef * JxBxB=ambiCoef * B^2(-J_perpB)
  !=> the source term J.E = ambiCoef * B^2 * J_perpB^2 = ambiCoef * JxBxB^2/B^2
  !> ambiCoef is calculated as uawsom_ambi_coef/rho^2,  see also the subroutine uawsom_get_Jambi
  subroutine add_source_ambipolar_internal_energy(qdt,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     wCT,w,x,ie)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: jxbxb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)

    call uawsom_get_jxbxb(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,jxbxb)
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sum(jxbxb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3)**2,&
       dim=ndim+1) / uawsom_mag_en_all(wCT, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp,wCT,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie)+qdt * tmp

  end subroutine add_source_ambipolar_internal_energy

  subroutine uawsom_get_jxbxb(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: res(:,:,:,:)

    double precision  :: btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)
    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       b2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    res=0.d0
    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    !!!here we know that current has nonzero values only for components in the range idirmin, 3
 
    if(B0field) then
      do idir=1,3
        btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir)) + block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir,b0i)
      enddo
    else
      btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1:3))
    endif

    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3)*btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3),dim=ndim+1) !J.B
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sum(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3)**2,&
       dim=ndim+1) !B2,!B2,!B2
    do idir=1,idirmin-1
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) = btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    enddo
    do idir=idirmin,3
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) = btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) - current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idir) * b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    enddo
  end subroutine uawsom_get_jxbxb

  !> Sets the sources for the ambipolar
  !> this is used for the STS method
  ! The sources are added directly (instead of fluxes as in the explicit)
  !> at the corresponding indices
  !>  store_flux_var is explicitly called for each of the fluxes one by one
  subroutine sts_set_source_ambipolar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
     fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,igrid,nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3) :: tmp,ff
    double precision :: fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nflux,1:ndim)
    double precision :: fE(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       7-2*ndim:3)
    double precision  :: btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3),tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: i, ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3, ie_

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;

    fluxall=zero

    call uawsom_get_jxbxb(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmp)

    !set electric field in tmp: E=nuA * jxbxb, where nuA=-etaA/rho^2
    do i=1,3
      !tmp(ixA^S,i) = -(uawsom_eta_ambi/w(ixA^S, rho_)**2) * tmp(ixA^S,i)
      call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,i),w,x)
    enddo

    if(uawsom_energy .and. .not.uawsom_internal_e) then
      !btot should be only mag. pert.
      btot(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:3)=0.d0
      !if(B0field) then
      !  do i=1,ndir
      !    btot(ixA^S, i) = w(ixA^S,mag(i)) + block%B0(ixA^S,i,0)
      !  enddo
      !else
        btot(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           1:ndir) = w(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           mag(1:ndir))
      !endif
      call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmp,btot,ff)
      call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      !- sign comes from the fact that the flux divergence is a source now
      wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    endif

    if(stagger_grid) then
      if(ndir>ndim) then
        !!!Bz
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           1) = tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,2)
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           2) = -tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1)
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,3) = 0.d0
        call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1+ndir,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndim)
        wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(ndir))=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if
      fE=0.d0
      call update_faces_ambipolar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,tmp,fE,&
         btot)
      ixAmax1=ixOmax1;ixAmax2=ixOmax2;ixAmax3=ixOmax3;
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;
      wres(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mag(1:ndim))=-btot(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         1:ndim)
    else
      !write curl(ele) as the divergence
      !m1={0,ele[[3]],-ele[[2]]}
      !m2={-ele[[3]],0,ele[[1]]}
      !m3={ele[[2]],-ele[[1]],0}

      !!!Bx
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1) = 0.d0
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         2) = tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,3)
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         3) = -tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,2)
      call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,2,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      !flux divergence is a source now
      wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1))=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      !!!By
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         1) = -tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,3)
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,2) = 0.d0
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         3) = tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1)
      call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,3,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(2))=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

      if(ndir==3) then
        !!!Bz
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           1) = tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,2)
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           2) = -tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1)
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,3) = 0.d0
        call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1+ndir,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndim)
        wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(ndir))=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if

    end if

    if(fix_conserve_at_step) then
      fluxall=my_dt*fluxall
      call store_flux(igrid,fluxall,1,ndim,nflux)
      if(stagger_grid) then
        call store_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           my_dt*fE,1,ndim)
      end if
    end if

  end subroutine sts_set_source_ambipolar

  !> get ambipolar electric field and the integrals around cell faces
  subroutine update_faces_ambipolar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,ECC,fE,circ)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! amibipolar electric field at cell centers
    double precision, intent(in)       :: ECC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)
    double precision, intent(out)      :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3)
    double precision, intent(out)      :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixAmin1,&
       ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3
    integer                            :: idim1,idim2,idir,ix1,ix2,ix3

    fE=zero
    ! calcuate ambipolar electric field on cell edges from cell centers
    do idir=7-2*ndim,3
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
      ixCmin3=ixOmin3+kr(idir,3)-1;
     do ix3=0,1
     do ix2=0,1
     do ix1=0,1
        if( ix1==1 .and. 1==idir  .or. ix2==1 .and. 2==idir  .or. ix3==1 .and. &
           3==idir ) cycle
        ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;ixAmin3=ixCmin3+ix3;
        ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;ixAmax3=ixCmax3+ix3;
        fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)+ECC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,idir)
     end do
     end do
     end do
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)*0.25d0*block%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,idir)
    end do

    ! Calculate circulation on each face to get value of line integral of
    ! electric field in the positive idir direction.
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
    ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;

    circ=zero

    do idim1=1,ndim ! Coordinate perpendicular to face
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
          hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
          hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
          ! Add line integrals in direction idir
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
             hxCmin3:hxCmax3,idir))
        end do
      end do
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)/block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)
    end do

  end subroutine update_faces_ambipolar

  !> use cell-center flux to get cell-face flux
  !> and get the source term as the divergence of the flux
  subroutine get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,src)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(:,:,:,:), intent(inout) :: ff
    double precision, intent(out) :: src(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: ffc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: dxinv(ndim)
    integer :: idims, ix1,ix2,ix3, ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,&
       ixAmax3, ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3, ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    dxinv=1.d0/dxlevel
    ! cell corner flux in ffc
    ffc=0.d0
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=ixOmin1-1
    ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;
    do ix3=0,1
    do ix2=0,1
    do ix1=0,1
      ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;ixBmin3=ixCmin3+ix3;
      ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;ixBmax3=ixCmax3+ix3;
      ffc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:ndim)=ffc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:ndim)+ff(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1:ndim)
    end do
    end do
    end do
    ffc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       1:ndim)=0.5d0**ndim*ffc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       1:ndim)
    ! flux at cell face
    ff(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=0.d0
    do idims=1,ndim
      ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
      ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
      ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=ixBmin1
      ixCmin2=ixBmin2;ixCmin3=ixBmin3;
      do ix3=0,1 
      do ix2=0,1 
      do ix1=0,1 
         if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims  .or. ix3==0 &
            .and. 3==idims ) then
           ixBmin1=ixCmin1-ix1;ixBmin2=ixCmin2-ix2;ixBmin3=ixCmin3-ix3;
           ixBmax1=ixCmax1-ix1;ixBmax2=ixCmax2-ix2;ixBmax3=ixCmax3-ix3;
           ff(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims)=ff(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims)+ffc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
              idims)
         end if
      end do
      end do
      end do
      ff(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idims)=ff(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idims)*0.5d0**(ndim-1)
    end do
    src=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idims)=dxinv(idims)*ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+ff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)-ff(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,idims)
      end do
    else
      do idims=1,ndim
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idims)=ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idims)*block%surfaceC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+ff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)-ff(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,idims)
      end do
      src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=src(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%dvolume(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if
  end subroutine get_flux_on_cell_face

  !> Calculates the explicit dt for the ambipokar term
  !> This function is used by both explicit scheme and STS method
  function get_ambipolar_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,&
     x)  result(dtnew)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: dtnew

    double precision              :: coef
    double precision              :: dxarr(ndim)
    double precision              :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    dxarr(1)=dx1;dxarr(2)=dx2;dxarr(3)=dx3;
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = uawsom_mag_en_all(w,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp,w,x)
    coef = maxval(abs(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))
    if(coef/=0.d0) then
      coef=1.d0/coef
    else
      coef=bigdouble
    end if
    if(slab_uniform) then
      dtnew=minval(dxarr(1:ndim))**2.0d0*coef
    else
      dtnew=minval(block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndim))**2.0d0*coef
    end if

  end function get_ambipolar_dt

  !> multiply res by the ambipolar coefficient
  !> The ambipolar coefficient is calculated as -uawsom_eta_ambi/rho^2
  !> The user may mask its value in the user file
  !> by implemneting usr_mask_ambipolar subroutine
  subroutine multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    tmp=0.d0
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=-&
       uawsom_eta_ambi/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2
    if (associated(usr_mask_ambipolar)) then
      call usr_mask_ambipolar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,tmp)
    end if

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
  end subroutine multiplyAmbiCoef

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine uawsom_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
     qsourcesplit,active,wCTprim)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source
    use mod_cak_force, only: cak_add_source

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active
    double precision, intent(in), optional :: wCTprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    if (.not. qsourcesplit) then
      if(uawsom_internal_e) then
        ! Source for solving internal energy
        active = .true.
        call internal_energy_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
           w,x,e_)
      else
        if(uawsom_solve_eaux) then
          ! Source for auxiliary internal energy equation
          active = .true.
          call internal_energy_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
             ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
             wCT,w,x,eaux_)
        endif
        if(has_equi_pe0) then
          active = .true.
          call add_pe0_divv(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
        endif
      endif

      ! Source for B0 splitting
      if (B0field.or.B0field) then
        active = .true.
        call add_source_B0split(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      end if
    end if

      
    if(.not.source_split_divb .and. .not.qsourcesplit .and. istep==nstep) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm)
        active = .true.
        call add_source_glm(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pso(block%igrid)%w,&
           w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           pso(block%igrid)%w,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           pso(block%igrid)%w,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           pso(block%igrid)%w,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           pso(block%igrid)%w,w,x)
        call add_source_janhunen(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           pso(block%igrid)%w,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           pso(block%igrid)%w,w,x)
        call add_source_powel(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           pso(block%igrid)%w,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           pso(block%igrid)%w,w,x)
        call add_source_glm(dt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pso(block%igrid)%w,&
           w,x)
      case (divb_ct)
        continue ! Do nothing
      case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    else if(source_split_divb .and. qsourcesplit) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm)
        active = .true.
        call add_source_glm(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
        call add_source_janhunen(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
        call add_source_powel(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
        call add_source_glm(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_ct)
        continue ! Do nothing
      case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
   

    if(uawsom_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,&
         x,qsourcesplit,active, rc_fl)
    end if

    if(uawsom_viscosity) then
      call viscosity_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
         uawsom_energy,qsourcesplit,active)
    end if

    if(uawsom_gravity) then
      call gravity_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
         gravity_energy,qsourcesplit,active)
    end if

    call w_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,uawsom_energy,&
       qsourcesplit,active)

    if (uawsom_cak_force) then
      call cak_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,uawsom_energy,&
         qsourcesplit,active)
    end if

  end subroutine uawsom_add_source

  subroutine add_pe0_divv(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    double precision                :: divv(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call uawsom_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,v)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,&
           sixthorder=.true.)
      else
        call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,&
           fourthorder=.true.)
      end if
    else
     call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
        ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv)
    end if
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-qdt*block%equi_vars(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)*divv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine add_pe0_divv

  subroutine w_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,energy,&
     qsourcesplit,active)
    use mod_global_parameters
    use mod_geometry
    use, intrinsic :: ieee_arithmetic

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision                :: Tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
    double precision                :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    double precision                :: divv(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),Lperp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       Lperp_AW(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       Gamma_minus(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       Gamma_plus(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       R(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       radius(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       zeta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       pth(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),B(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,3)
    double precision :: Bsqr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       vA(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),gvA(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3),gradvA(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),refvA(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       Bfix(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       lnvA(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       lnrho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        gradlnvA(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) 
    integer :: half_window, ix1, ix2
    double precision :: smoothed_rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), sum_rho

    !> may use these later: ,B2(ixI^S),vA(ixI^S),gvA(ixI^S),gradvA(ixI^S),refvA(ixI^S)

    !if (mod(it,1000)==0 .and. mype==0) then
    !  write(*,*) 'qt = ', qt
    !end if

    if(B0field) then
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,1)
    else
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))
    end if
    
    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1.d8/unit_length * &
       ((Busr/unit_magneticfield)/B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1))**0.5d0 !Radius = R_0*(B0/B)0.5,R_0*(B0/B)0.5,R_0*(B0/B)0.5 TVD 2025 paper uses R_0 = 1Mm (1e8 cm)

    if(qsourcesplit .eqv. .false.) then
      active = .true.
    endif

    call uawsom_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,v)
    !call uawsom_get_pthermal_origin(w,x,ixI^L,ixO^L,pth)
    !pth(ixO^S) = pth(ixO^S)/w(ixO^S,rho_)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,&
           sixthorder=.true.)
      else
        call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,&
           fourthorder=.true.)
      end if
    else
     call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
        ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv)
    end if

    call get_zeta(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,zeta)

    !B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    !Bsqr(ixO^S) = B(ixO^S,1)**(2.d0)

    !do ix1 = ixImin1,ixImin1+1
    !  B(ix1,1) = B(ixOmin1,1) !> ixOmin1^%1ixO^S
    !end do

    !do ix1 = ixImax1-1,ixImax1  
    !  B(ix1,1) = B(ixOmax1,1)
    !end do

    ! output Alfven wave speed B/sqrt(rho)
    !vA(ixI^S) = B(ixI^S,1)/dsqrt(w(ixI^S,rho_))

    !lnvA(ixI^S) = log(vA(ixI^S))
    !lnrho(ixI^S) = log(w(ixI^S,rho_))

    !call gradient(lnvA,ixI^L,ixO^L,ndim,gradlnvA)
    !call gradient(vA,ixI^L,ixO^L,ndim,gradvA)
    !call gradient(lnrho,ixI^L,ixO^L,ndim,gradlnrho)
    
    !refvA(ixO^S) = (gradvA(ixO^S)**2)/(vA(ixO^S)**2) !> Old idea has wrong units  

    !Lperp(ixO^S) = dsqrt(10.d0*ff*dpi)*radius(ixO^S)*(zeta(ixO^S) + 1.d0 - ff)**(3.d0/2.d0)/((1.d0 - ff**(5.d0/2.d0))*(zeta(ixO^S) - 1.d0))

    Lperp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = (zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + 1.d0 - ff)**(3.d0/2.d0)/(1.d0 - &
       ff**(5.d0/2.d0))/(zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) - 1.d0)*3.1622776*(ff*dpi)**&
       0.5d0*radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) !/min(pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),1.d0) 
    !Lperp(ixO^S) = ((1/(radius(ixO^S)*(ff*dpi)**0.5d0))*((2.d0/5.d0)**0.5d0)*((zeta(ixO^S)-1)/2)*((1-ff**2.5d0)/(zeta(ixO^S)+1-ff)**1.5d0))**(-1.d0)

    ! Max: This is the last term in Eqn. 78 TVD 2024 contribution due to AWs is zero since mu*rho*alpha**2-1 = 0
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-qdt*(zeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)-1.d0)/(zeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)+1.d0)*(zeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)+1.d0)*(wCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,wkplus_) + wCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,wkminus_))/4.d0*divv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) 

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_)  - qdt*(divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_)/2.d0  + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkplus_)**(3.d0/2.d0)/(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)*(1+ff*zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-ff)**(-1.d0))**0.5d0/Lperp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_) - qdt*(divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_)/2.d0 + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_)**(3.d0/2.d0)/(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)*(1+ff*zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-ff)**(-1.d0))**0.5d0/Lperp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3))

    call uawsom_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,pth)
    Tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)

    Lperp_AW(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = (1/1.d8) * 1.5d5 * 1.0d2 * (Tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) / B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1))**0.5d0 

    !Lperp_AW(ixO^S) = 0.01d0*(20.0d0/B(ixO^S,1))**0.5d0 
    Gamma_plus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = (2.0d0 / Lperp_AW(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        wAminus_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_))**0.5d0
    Gamma_minus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = (2.0d0 / Lperp_AW(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        wAplus_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_))**0.5d0

    !do ix1 = ixOmin1,ixOmax1
      !refvA(ix1) = min(vA(ix1)*gradlnvA(ix1), max(Gamma_plus(ix1),Gamma_minus(ix1))) 
      !refvA(ix1^%1ixO^S) = vA(ix1^%1ixO^S)*gradlnvA(ix1^%1ixO^S)*0.0d0 !> Note vA*grad(ln(vA)) = gradvA in 1D
    !end do
  
    !> Single reflection ie waves cannot be re-reflected as per VDH2014 with moving frame taken into account (large upflows increase reflection)
    do ix1 = ixOmin1, ixOmax1 
      w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAminus_) = w(ix1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,wAminus_) - qdt*(divv(ix1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*wCT(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         wAminus_)/2.0d0 + wCT(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         wAminus_)*Gamma_minus(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) !+ ((w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1)) + vA(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))/vA(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*refvA(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*wCT(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAminus_))
      w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAplus_)  = w(ix1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,wAplus_)  - qdt*(divv(ix1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*wCT(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         wAplus_)/2.0d0 + wCT(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         wAplus_)*Gamma_plus(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) !- ((w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1)) + vA(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))/vA(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))*refvA(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*wCT(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAminus_))
    end do

  end subroutine w_add_source

  !> Compute the Lorentz force (JxB)
  subroutine get_Lorentz_force(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,JxB)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: JxB(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3)
    double precision                :: a(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3), b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       3)
    integer                         :: idir, idirmin
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for uawsom
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)

    b=0.0d0
    do idir = 1, ndir
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir) = uawsom_mag_i_all(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir)
    end do

    ! store J current in a
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)

    a=0.0d0
    do idir=7-2*ndir,3
      a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
    end do

    call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,JxB)
  end subroutine get_Lorentz_force

  subroutine internal_energy_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
     ie)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision                :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir),divv(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call uawsom_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,v)
    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,&
           sixthorder=.true.)
      else
        call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,&
           fourthorder=.true.)
      end if
    else
     call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
        ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv)
    end if
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie)-qdt*wCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie)*gamma_1*divv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    if(uawsom_ambipolar)then
       call add_source_ambipolar_internal_energy(qdt,ixImin1,ixImin2,ixImin3,&
          ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,wCT,w,x,ie)
    end if
    if(fix_small_values) then
      call uawsom_handle_small_ei(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ie,&
         'internal_energy_add_source')
    end if
  end subroutine internal_energy_add_source

  subroutine uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(has_equi_rho0) then
      rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,b0i)
    else
      rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    endif

  end subroutine uawsom_get_rho
  
  subroutine get_zeta(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,zeta)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: zeta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    
    !zeta(ixI^S) = zeta0*exp(-(x(ixI^S,1)-xprobmin1)/5.d0)
    zeta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = (zeta0-1.d0)*exp(-(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1)-xprobmin1)/5.d0)+1.d0
    !where(zeta(ixI^S) < 1.d0)
    !  zeta(ixI^S) = 1.d0
    !end where

  end subroutine get_zeta

  !> handle small or negative internal energy
  subroutine uawsom_handle_small_ei(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ie,&
      subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ie
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision              :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    flag=.false.
    if(has_equi_pe0) then
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         equi_pe0_,0)*inv_gamma_1<small_e)flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,ie)=.true.
    else
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)<small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)=.true.
    endif
    if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie))) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_pe0) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)=small_e - block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe0_,0)*inv_gamma_1
        else
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)=small_e
        endif
      case ("average")
        call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
            flag, ie)
      case default
        ! small values error shows primitive variables
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)*gamma_1
        call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
        do idir = 1, ndir
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mom(idir))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if

  end subroutine uawsom_handle_small_ei

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: a(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3),&
        b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3),&
        axb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3)
    integer :: idir

    a=0.d0
    b=0.d0
    ! for force-free field J0xB0 =0
    if(.not.B0field_forcefree) then
      ! store B0 magnetic field in b
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,0)

      ! store J0 current in a
      do idir=7-2*ndir,3
        a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)
      end do
      call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,axb)
      axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)=axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*qdt
      ! add J0xB0 source term in momentum equations
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1:ndir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1:ndir))+axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)
    end if

    if(total_energy) then
      a=0.d0
      ! for free-free field -(vxB0) dot J0 =0
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)=wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))
      ! store full magnetic field B0+B1 in b
      if(.not.B0field_forcefree) b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)=b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0)
      ! store velocity in a
      call uawsom_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1:ndir))
      call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,axb)
      axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)=axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*qdt
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)-axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)*block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)
      end do
      if(uawsom_ambipolar) then
        !reuse axb
        call uawsom_get_jxbxb(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,axb)
        ! source J0 * E
        do idir=7-2*ndim,3
          !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
          call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
             axb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idir),wCT,x)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)+axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)*block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)
        enddo
      endif
    end if

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_B0')

  end subroutine add_source_B0split

  !> Add resistive source to w within ixO Uses 3 point stencil (1 neighbour) in
  !> each direction, non-conservative. If the fourthorder precompiler flag is
  !> set, uses fourth order central difference for the laplacian. Then the
  !> stencil is 5 (2 neighbours).
  subroutine add_source_res1(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idir,jdir,kdir,&
       idirmin,idim,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,ix
    integer :: lxOmin1,lxOmin2,lxOmin3,lxOmax1,lxOmax2,lxOmax3, kxOmin1,&
       kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3

    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for uawsom
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: gradeta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), Bf(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)

    ! Calculating resistive sources involve one extra layer
    if (uawsom_4th_order) then
      ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmin3=ixOmin3-2;ixAmax1=ixOmax1+2
      ixAmax2=ixOmax2+2;ixAmax3=ixOmax3+2;
   else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
      ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    end if

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2.or.ixImin3>ixAmin3.or.ixImax3<ixAmax3) call &
       mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)

    if (uawsom_eta>zero)then
       eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)=uawsom_eta
       gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImin3,ixImax1,&
          ixImax2,ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,&
          idirmin,x,current,eta)
       ! assumes that eta is not function of current?
       do idim=1,ndim
          call gradient(eta,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,tmp)
          gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       end do
    end if

    if(B0field) then
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir,0)
    else
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndir))
    end if

    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (uawsom_4th_order) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3)=Bf(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,idir)
         do idim=1,ndim
            lxOmin1=ixOmin1+2*kr(idim,1);lxOmin2=ixOmin2+2*kr(idim,2)
            lxOmin3=ixOmin3+2*kr(idim,3);lxOmax1=ixOmax1+2*kr(idim,1)
            lxOmax2=ixOmax2+2*kr(idim,2);lxOmax3=ixOmax3+2*kr(idim,3);
            jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
            jxOmin3=ixOmin3+kr(idim,3);jxOmax1=ixOmax1+kr(idim,1)
            jxOmax2=ixOmax2+kr(idim,2);jxOmax3=ixOmax3+kr(idim,3);
            hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
            hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
            hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
            kxOmin1=ixOmin1-2*kr(idim,1);kxOmin2=ixOmin2-2*kr(idim,2)
            kxOmin3=ixOmin3-2*kr(idim,3);kxOmax1=ixOmax1-2*kr(idim,1)
            kxOmax2=ixOmax2-2*kr(idim,2);kxOmax3=ixOmax3-2*kr(idim,3);
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+(-tmp2(lxOmin1:lxOmax1,lxOmin2:lxOmax2,&
               lxOmin3:lxOmax3)+16.0d0*tmp2(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3)-30.0d0*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+16.0d0*tmp2(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3)-tmp2(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
               kxOmin3:kxOmax3)) /(12.0d0 * dxlevel(idim)**2)
         end do
       else
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3)=Bf(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,idir)
         do idim=1,ndim
            jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
            jxOmin3=ixOmin3+kr(idim,3);jxOmax1=ixOmax1+kr(idim,1)
            jxOmax2=ixOmax2+kr(idim,2);jxOmax3=ixOmax3+kr(idim,3);
            hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
            hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
            hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+(tmp2(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3)-2.0d0*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+tmp2(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3))/dxlevel(idim)**2
         end do
       end if

       ! Multiply by eta
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)

       ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
       if (uawsom_eta<zero)then
          do jdir=1,ndim; do kdir=idirmin,3
             if (lvc(idir,jdir,kdir)/=0)then
                if (lvc(idir,jdir,kdir)==1)then
                   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3)-gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3,jdir)*current(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
                else
                   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3)+gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3,jdir)*current(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
                end if
             end if
          end do; end do
       end if

       ! Add sources related to eta*laplB-grad(eta) x J to B and e
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mag(idir))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       if(total_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*Bf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir)
       end if
    end do ! idir

    if(uawsom_energy) then
      ! de/dt+=eta*J**2
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=qdt*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)**2,dim=ndim+1)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)+tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      if(uawsom_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           eaux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           eaux_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if
    end if

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_res1')

  end subroutine add_source_res1

  !> Add resistive source to w within ixO
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_res2(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for uawsom
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),curlj(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)
    double precision :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3),tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idir,idirmin,&
       idirmin1

    ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmin3=ixOmin3-2;ixAmax1=ixOmax1+2
    ixAmax2=ixOmax2+2;ixAmax3=ixOmax3+2;

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2.or.ixImin3>ixAmin3.or.ixImax3<ixAmax3) call &
       mpistop("Error in add_source_res2: Non-conforming input limits")

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,current)

    tmpvec=zero
    if(uawsom_eta>zero)then
      do idir=idirmin,3
        tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)=current(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)*uawsom_eta
      end do
    else
      call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,x,&
         current,eta)
      do idir=idirmin,3
        tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)=current(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)*eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)
      end do
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    call curlvector(tmpvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,curlj,idirmin1,1,3)
    if(stagger_grid) then
      if(ndim==2.and.ndir==3) then
        ! if 2.5D
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(ndir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(ndir))-qdt*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,ndir)
      end if
    else
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))-qdt*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir)
    end if

    if(uawsom_energy) then
      if(uawsom_eta>zero)then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=qdt*uawsom_eta*sum(current(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,dim=ndim+1)
      else
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=qdt*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,:)**2,dim=ndim+1)
      end if
      if(total_energy) then
        ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
        ! de1/dt= eta J^2 - B1 dot curl(eta J)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-qdt*sum(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mag(1:ndir))*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1:ndir),dim=ndim+1)
      else
        ! add eta*J**2 source term in the internal energy equation
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if
      if(uawsom_solve_eaux) then
        ! add eta*J**2 source term in the internal energy equation
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           eaux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           eaux_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if
    end if

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_res2')
  end subroutine add_source_res2

  !> Add Hyper-resistive source to w within ixO
  !> Uses 9 point stencil (4 neighbours) in each direction.
  subroutine add_source_hyperres(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    !.. local ..
    double precision                :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision                :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3),tmpvec2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3),tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),ehyper(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)
    integer                         :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,&
       ixAmax3,idir,jdir,kdir,idirmin,idirmin1

    ixAmin1=ixOmin1-3;ixAmin2=ixOmin2-3;ixAmin3=ixOmin3-3;ixAmax1=ixOmax1+3
    ixAmax2=ixOmax2+3;ixAmax3=ixOmax3+3;
    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2.or.ixImin3>ixAmin3.or.ixImax3<ixAmax3) call &
       mpistop("Error in add_source_hyperres: Non-conforming input limits")

    call get_current(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,current)
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
          jdir)=current(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,jdir)
    end do

    ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmin3=ixOmin3-2;ixAmax1=ixOmax1+2
    ixAmax2=ixOmax2+2;ixAmax3=ixOmax3+2;
    call curlvector(tmpvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmpvec2,idirmin1,1,3)

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndir)=zero
    call curlvector(tmpvec2,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmpvec,idirmin1,1,3)
    ehyper(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
       1:ndir) = - tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
       1:ndir)*uawsom_eta_hyper

    ixAmin1=ixOmin1;ixAmin2=ixOmin2;ixAmin3=ixOmin3;ixAmax1=ixOmax1
    ixAmax2=ixOmax2;ixAmax3=ixOmax3;
    tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndir)=zero
    call curlvector(ehyper,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmpvec2,idirmin1,1,3)

    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))-tmpvec2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)*qdt
    end do

    if(total_energy) then
      ! de/dt= +div(B x Ehyper)
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
      ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
      tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndir)=zero
      do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
        tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir) = tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)+ lvc(idir,jdir,kdir)*wCT(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,mag(jdir))*ehyper(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,kdir)
      end do; end do; end do
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
      call divvector(tmpvec2,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)+tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)*qdt
    end if

    if (fix_small_values)  call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_hyperres')

  end subroutine add_source_hyperres

  subroutine add_source_glm(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-uawsom scheme or GLM-uawsom scheme
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision:: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer          :: idim,idir
    double precision :: gradPsi(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)


    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (uawsom_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_) = abs(uawsom_glm_alpha)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           psi_) = dexp(-qdt*cmax_global*uawsom_glm_alpha/minval(dxlevel(&
           :)))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           psi_) = dexp(-qdt*cmax_global*uawsom_glm_alpha/minval(block%ds(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:),&
           dim=ndim+1))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           psi_)
      end if
    end if

    if(uawsom_glm_extended) then
      ! gradient of Psi
      if(total_energy) then
        do idim=1,ndim
          select case(typegrad)
          case("central")
            call gradient(wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               psi_),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
               ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,gradPsi)
          case("limited")
            call gradientS(wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               psi_),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
               ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,gradPsi)
          end select
          ! e  = e  -qdt (b . grad(Psi))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idim))*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        end do
      end if

      ! We calculate now div B
      call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb,&
          uawsom_divb_4thorder)

      ! m = m - qdt b div b
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))-qdt*uawsom_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do
    end if

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_glm')

  end subroutine add_source_glm

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb, uawsom_divb_4thorder)

    ! calculate velocity
    call uawsom_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)

    if (total_energy) then
      ! e = e - qdt (v . b) * div b
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-qdt*sum(v(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)),&
         dim=ndim+1)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))-qdt*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))-qdt*uawsom_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb, uawsom_divb_4thorder)

    call uawsom_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vel)
    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))-qdt*vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_janhunen')

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    ! Add Linde`s divB related sources to wnew within ixO
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer :: idim, idir, ixpmin1,ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3, i1,&
       i2,i3, iside
    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       graddivb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    logical, dimension(-1:1,-1:1,-1:1) :: leveljump

    ! Calculate div B
    ixpmin1=ixOmin1-1;ixpmin2=ixOmin2-1;ixpmin3=ixOmin3-1;ixpmax1=ixOmax1+1
    ixpmax2=ixOmax2+1;ixpmax3=ixOmax3+1;
    call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixpmin1,&
       ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3,divb, uawsom_divb_4thorder)

    ! for AMR stability, retreat one cell layer from the boarders of level jump
    do i3=-1,1
    do i2=-1,1
    do i1=-1,1
      if(i1==0.and.i2==0.and.i3==0) cycle
      if(neighbor_type(i1,i2,i3,block%igrid)==2 .or. neighbor_type(i1,i2,i3,&
         block%igrid)==4) then
        leveljump(i1,i2,i3)=.true.
      else
        leveljump(i1,i2,i3)=.false.
      end if
    end do
    end do
    end do

    ixpmin1=ixOmin1;ixpmin2=ixOmin2;ixpmin3=ixOmin3;ixpmax1=ixOmax1
    ixpmax2=ixOmax2;ixpmax3=ixOmax3;
    do idim=1,ndim
      select case(idim)
       case(1)
          do iside=1,2
            i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3)
            i3=kr(3,1)*(2*iside-3);
            if (leveljump(i1,i2,i3)) then
              if (iside==1) then
                ixpmin1=ixOmin1-i1
              else
                ixpmax1=ixOmax1-i1
              end if
            end if
          end do
       
       case(2)
          do iside=1,2
            i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3)
            i3=kr(3,2)*(2*iside-3);
            if (leveljump(i1,i2,i3)) then
              if (iside==1) then
                ixpmin2=ixOmin2-i2
              else
                ixpmax2=ixOmax2-i2
              end if
            end if
          end do
       
       case(3)
          do iside=1,2
            i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3)
            i3=kr(3,3)*(2*iside-3);
            if (leveljump(i1,i2,i3)) then
              if (iside==1) then
                ixpmin3=ixOmin3-i3
              else
                ixpmax3=ixOmax3-i3
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
         call gradient(divb,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixpmin1,ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3,idim,graddivb)
       case("limited")
         call gradientS(divb,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixpmin1,ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3,idim,graddivb)
       end select

       ! Multiply by Linde`s eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab_uniform) then
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)=graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)*divbdiff/(1.0d0/dxlevel(1)**2+&
             1.0d0/dxlevel(2)**2+1.0d0/dxlevel(3)**2)
       else
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)=graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)*divbdiff /(1.0d0/block%ds(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
             1)**2+1.0d0/block%ds(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3,2)**2+1.0d0/block%ds(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,ixpmin3:ixpmax3,3)**2)
       end if

       w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
          mag(idim))=w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
          mag(idim))+graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3)

       if (typedivbdiff=='all' .and. total_energy) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
            e_)=w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
            e_)+wCT(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
            mag(idim))*graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
            ixpmin3:ixpmax3)
       end if
    end do

    if (fix_small_values) call uawsom_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_linde')

  end subroutine add_source_linde

  !> Calculate div B within ixO
  subroutine get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb, fourthorder)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    logical, intent(in), optional   :: fourthorder

    integer                            :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
       ixCmax2,ixCmax3, idir

    if(stagger_grid) then
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=0.d0
      do idir=1,ndim
        ixCmin1=ixOmin1-kr(idir,1);ixCmin2=ixOmin2-kr(idir,2)
        ixCmin3=ixOmin3-kr(idir,3);ixCmax1=ixOmax1-kr(idir,1)
        ixCmax2=ixOmax2-kr(idir,2);ixCmax3=ixOmax3-kr(idir,3);
        divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+block%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir)*block%surfaceC(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)-block%ws(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)*block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idir)
      end do
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      select case(typediv)
      case("central")
        call divvector(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           mag(1:ndir)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb,fourthorder)
      case("limited")
        call divvectorS(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           mag(1:ndir)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
      end select
    end if

  end subroutine get_divb

  !> get dimensionless div B = |divB| * volume / area / |B|
  subroutine get_normalized_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)

    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                   :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), dsurface(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idims

    call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
    invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(uawsom_mag_en_all(w,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3))
    where(invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/=0.d0)
      invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    end where
    if(slab_uniform) then
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=0.5d0*abs(divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/sum(1.d0/dxlevel(:))
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;
      ixAmax1=ixOmax1-1;ixAmax2=ixOmax2-1;ixAmax3=ixOmax3-1;
      dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)= sum(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:),dim=ndim+1)
      do idims=1,ndim
        ixAmin1=ixOmin1-kr(idims,1);ixAmin2=ixOmin2-kr(idims,2)
        ixAmin3=ixOmin3-kr(idims,3);ixAmax1=ixOmax1-kr(idims,1)
        ixAmax2=ixOmax2-kr(idims,2);ixAmax3=ixOmax3-kr(idims,3);
        dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+block%surfaceC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)
      end do
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=abs(divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    end if

  end subroutine get_normalized_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)  :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer, intent(out) :: idirmin
    integer :: idir, idirmin0

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for uawsom
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)

    idirmin0 = 7-2*ndir

    call curlvector(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       mag(1:ndir)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,current,idirmin,idirmin0,ndir)

    if(B0field) current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin0:3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin0:3)+block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin0:3)
  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine uawsom_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_cak_force, only: cak_get_dt

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx1,dx2,dx3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    integer                       :: idirmin,idim
    double precision              :: dxarr(ndim)
    double precision              :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    dtnew = bigdouble

    dxarr(1)=dx1;dxarr(2)=dx2;dxarr(3)=dx3;
    if (uawsom_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/uawsom_eta
    else if (uawsom_eta<zero)then
       call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
       call usr_special_resistivity(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,x,&
          current,eta)
       dtnew=bigdouble
       do idim=1,ndim
         if(slab_uniform) then
           dtnew=min(dtnew,dtdiffpar/(smalldouble+maxval(eta(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)/dxarr(idim)**2)))
         else
           dtnew=min(dtnew,dtdiffpar/(smalldouble+maxval(eta(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%ds(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)**2)))
         end if
       end do
    end if

    if(uawsom_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/uawsom_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1:ndim))**4/uawsom_eta_hyper,dtnew)
      end if
    end if

    if(uawsom_radiative_cooling) then
      call cooling_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x,&
         rc_fl)
    end if

    if(uawsom_viscosity) then
      call viscosity_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    end if

    if(uawsom_gravity) then
      call gravity_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    end if

    if(uawsom_ambipolar_exp) then
      dtnew=min(dtdiffpar*get_ambipolar_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,&
         dx2,dx3,x),dtnew)
    endif

    if (uawsom_cak_force) then
      call cak_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    end if

  end subroutine uawsom_get_dt

  ! Add geometrical source terms to w
  subroutine uawsom_add_source_geom(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    integer          :: iw,idir, h1xmin1,h1xmin2,h1xmin3,h1xmax1,h1xmax2,&
       h1xmax3, h2xmin1,h2xmin2,h2xmin3,h2xmax1,h2xmax2,h2xmax3
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    ! 1/rho
    invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.d0/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)
    invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.d0/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)
    select case (coordinate)
    case (cylindrical)
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in uawsom")
      endif
      call uawsom_get_p_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
      if(phi_>0) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)**2+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mphi_)**2*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mr_)*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,br_))
        if(.not.stagger_grid) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mr_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mphi_)) *invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        end if
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
      if(uawsom_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)+qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    case (spherical)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmin3=ixOmin3-kr(1,3);h1xmax1=ixOmax1-kr(1,1)
       h1xmax2=ixOmax2-kr(1,2);h1xmax3=ixOmax3-kr(1,3)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmin3=ixOmin3-kr(2,3);h2xmax1=ixOmax1-kr(2,1)
       h2xmax2=ixOmax2-kr(2,2);h2xmax3=ixOmax3-kr(2,3);
       call uawsom_get_p_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp1)
       ! m1
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)-block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,h1xmin3:h1xmax3,&
          1))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       do idir=2,ndir
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(idir))**2*invrho(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3)-wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))**2
       end do
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       ! b1
       if(uawsom_glm) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*2.0d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,psi_)
       end if

       
       ! m2
       ! This will make hydrostatic p=const an exact solution
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))+qdt*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,2)-block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          h2xmin3:h2xmax3,2)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(2))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mag(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mag(2)))
       if(ndir==3) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(3))**2*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(3))**2)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(1)))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
         if(uawsom_glm) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) + dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,psi_)
         end if
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       end if
      

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(1))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(1)))  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3))) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mom(3))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)
         else
           call mpistop("angmomfix not implemented yet in uawsom")
         end if
         ! b3
         if(.not.stagger_grid) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(1)))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2)) /(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,rho_)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2)))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)
         end if
       end if
    end select
  end subroutine uawsom_add_source_geom

  ! Add geometrical source terms to w
  subroutine uawsom_add_source_geom_split(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    integer          :: iw,idir, h1xmin1,h1xmin2,h1xmin3,h1xmax1,h1xmax2,&
       h1xmax3, h2xmin1,h2xmin2,h2xmin3,h2xmax1,h2xmax2,h2xmax3
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    if(has_equi_rho0) then
      invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,b0i))
    else
      invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)
    end if
    invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1d0/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)

    select case (coordinate)
    case (cylindrical)
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in uawsom")
      endif
      call uawsom_get_p_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
      if(phi_>0) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)**2+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mphi_)**2*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mr_)*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,br_))
        if(.not.stagger_grid) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mr_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mphi_)) *invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        end if
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
      if(uawsom_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)+qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    case (spherical)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmin3=ixOmin3-kr(1,3);h1xmax1=ixOmax1-kr(1,1)
       h1xmax2=ixOmax2-kr(1,2);h1xmax3=ixOmax3-kr(1,3)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmin3=ixOmin3-kr(2,3);h2xmax1=ixOmax1-kr(2,1)
       h2xmax2=ixOmax2-kr(2,2);h2xmax3=ixOmax3-kr(2,3);
       call uawsom_get_p_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp1)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if(B0field) then
         tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=sum(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,:,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(:)),dim=ndim+1)
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       end if
       ! m1
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)-block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,h1xmin3:h1xmax3,&
          1))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(idir))**2*invrho(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)-wCT(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))**2
           if(B0field) tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)-2.0d0*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,idir,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(idir))
         end do
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       ! b1
       if(uawsom_glm) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*2.0d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,psi_)
       end if

       
       ! m2
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if(B0field) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       end if
       ! This will make hydrostatic p=const an exact solution
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,2)-block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          h2xmin3:h2xmax3,2)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(2))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mag(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mag(2)))
       if (B0field) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2,0)
       end if
       if(ndir==3) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(3))**2*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(3))**2)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)-2.0d0*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,3,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,mag(3))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2))
         end if
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(1)))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
         if(B0field) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2,0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,1,0))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)
         end if
         if(uawsom_glm) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) + dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,psi_)
         end if
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       end if
      

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(1))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(1)))  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3))) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))
           if (B0field) then
              tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,1,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mag(3)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mag(1))*block%B0(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,&
                 0)  +(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,2,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mag(3)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mag(2))*block%B0(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,0)) *dcos(x(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
           end if
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mom(3))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)
         else
           call mpistop("angmomfix not implemented yet in uawsom")
         end if
         ! b3
         if(.not.stagger_grid) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(1)))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2)) *invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))
           if (B0field) then
              tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mom(1))*block%B0(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,0) -wCT(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(3))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,1,0))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mom(3))*block%B0(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,0) -wCT(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,3,0))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,2)) *invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,2))
           end if
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)
         end if
       end if
    end select
  end subroutine uawsom_add_source_geom_split

  !> Compute 2 times total magnetic energy
  function uawsom_mag_en_all(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if (B0field) then
      mge = sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,&
         b0i))**2, dim=ndim+1)
    else
      mge = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(:))**2,&
          dim=ndim+1)
    end if
  end function uawsom_mag_en_all

  !> Compute full magnetic field by direction
  function uawsom_mag_i_all(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: mgf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if (B0field) then
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mag(idir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir,b0i)
    else
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(idir))
    end if
  end function uawsom_mag_i_all

  !> Compute evolving magnetic energy
  function uawsom_mag_en(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    mge = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        mag(:))**2, dim=ndim+1)
  end function uawsom_mag_en

  !> Compute wkminus and wkplus energy
  function uawsom_wk_en(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(we)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: we(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    we = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wkminus_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wkplus_)

  end function uawsom_wk_en

  !> Compute wAminus and wAplus energy
  function uawsom_wA_en(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(we)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: we(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    we = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       wAminus_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,wAplus_)

  end function uawsom_wA_en

  !> compute kinetic energy
  function uawsom_kin_en_origin(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    if (present(inv_rho)) then
      ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(:))**2, dim=ndim+1) * inv_rho
    else
      if(has_equi_rho0) then
        ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(:))**2, dim=ndim+1) / (w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_rho0_,0))
      else
        ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(:))**2, dim=ndim+1) / w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3, rho_)
      endif
    end if
  end function uawsom_kin_en_origin

  !> compute kinetic energy
  function uawsom_kin_en_boris(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      inv_rho) result(ke)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    if (present(inv_rho)) then
      ke=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ke=0.5d0*sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(:)))**2,dim=ndim+1)*ke**2*inv_rho
    else
      ke=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)*inv_squared_c)
      ke=0.5d0*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(:))**2,dim=ndim+1)*ke**2/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, rho_)
    end if
  end function uawsom_kin_en_boris

  subroutine uawsom_getv_Hall(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: vHall(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)

    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call uawsom_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3) = zero
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3) = - uawsom_etah*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,idirmin:3)
    do idir = idirmin, 3
       vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir) = vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

  end subroutine uawsom_getv_Hall

  subroutine uawsom_get_Jambi(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, allocatable, intent(inout) :: res(:,:,:,:)


    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)

    res = 0d0

    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
 
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3)=-current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3)
    do idir = idirmin, 3
      call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,idir),w,x)
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

  subroutine uawsom_modify_wLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,wLC,wRC,wLp,wRp,s,&
     idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    type(state)                     :: s
    double precision                :: dB(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), dPsi(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    if(stagger_grid) then
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=s%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=s%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=s%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=s%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
    else
      ! Solve the Riemann problem for the linear 2x2 system for normal
      ! B-field and GLM_Psi according to Dedner 2002:
      ! This implements eq. (42) in Dedner et al. 2002 JcP 175
      ! Gives the Riemann solution on the interface
      ! for the normal B component and Psi in the GLM-uawsom system.
      ! 23/04/2013 Oliver Porth
      dB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)   = wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idir)) - wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idir))
      dPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_) - wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_)

      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))   = 0.5d0 * (wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idir)) + wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idir))) - 0.5d0/cmax_global * &
         dPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)       = 0.5d0 * (wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_) + wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_)) - 0.5d0*cmax_global * dB(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)

      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir)) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))
      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)

      if(total_energy) then
        wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)-half*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
        wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)-half*wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
      end if
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir)) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir)) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)
      ! modify total energy according to the change of magnetic field
      if(total_energy) then
        wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+half*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
        wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+half*wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
      end if
    end if

    if(associated(usr_set_wLR)) call usr_set_wLR(ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine uawsom_modify_wLR

  subroutine uawsom_boundary_adjust(igrid,psb)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state), target :: psb(max_blocks)

    integer :: iB, idims, iside, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, i1,i2,i3

    block=>ps(igrid)
    dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
    dxlevel(3)=rnode(rpdx3_,igrid);
    do idims=1,ndim
       ! to avoid using as yet unknown corner info in more than 1D, we
       ! fill only interior mesh ranges of the ghost cell ranges at first,
       ! and progressively enlarge the ranges to include corners later
       do iside=1,2
          i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
          i3=kr(3,idims)*(2*iside-3);
          if (neighbor_type(i1,i2,i3,igrid)/=1) cycle
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
                  ixOmin1=ixGhi1+1-nghostcells+boundary_divbfix_skip(2*1)
                  ixOmin2=ixGlo2;ixOmin3=ixGlo3;
                  ixOmax1=ixGhi1;ixOmax2=ixGhi2;ixOmax3=ixGhi3;
               else
                  ! minimal boundary
                  ixOmin1=ixGlo1;ixOmin2=ixGlo2;ixOmin3=ixGlo3;
                  ixOmax1=ixGlo1-1+nghostcells-boundary_divbfix_skip(2*1-1)
                  ixOmax2=ixGhi2;ixOmax3=ixGhi3;
               end if 
            case (2)
               if (iside==2) then
                  ! maximal boundary
                  ixOmin1=ixGlo1
                  ixOmin2=ixGhi2+1-nghostcells+boundary_divbfix_skip(2*2)
                  ixOmin3=ixGlo3;
                  ixOmax1=ixGhi1;ixOmax2=ixGhi2;ixOmax3=ixGhi3;
               else
                  ! minimal boundary
                  ixOmin1=ixGlo1;ixOmin2=ixGlo2;ixOmin3=ixGlo3;
                  ixOmax1=ixGhi1
                  ixOmax2=ixGlo2-1+nghostcells-boundary_divbfix_skip(2*2-1)
                  ixOmax3=ixGhi3;
               end if 
            case (3)
               if (iside==2) then
                  ! maximal boundary
                  ixOmin1=ixGlo1;ixOmin2=ixGlo2
                  ixOmin3=ixGhi3+1-nghostcells+boundary_divbfix_skip(2*3);
                  ixOmax1=ixGhi1;ixOmax2=ixGhi2;ixOmax3=ixGhi3;
               else
                  ! minimal boundary
                  ixOmin1=ixGlo1;ixOmin2=ixGlo2;ixOmin3=ixGlo3;
                  ixOmax1=ixGhi1;ixOmax2=ixGhi2
                  ixOmax3=ixGlo3-1+nghostcells-boundary_divbfix_skip(2*3-1);
               end if 
            end select
            call fixdivB_boundary(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
               ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,psb(igrid)%w,&
               psb(igrid)%x,iB)
          end if
       end do
    end do

  end subroutine uawsom_boundary_adjust

  subroutine fixdivB_boundary(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,iB)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)

    double precision :: dx1x2,dx1x3,dx2x1,dx2x3,dx3x1,dx3x2
    integer :: ix1,ix2,ix3,ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3

    select case(iB)
     case(1)
       ! 2nd order CD for divB=0 to set normal B component better
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
       
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=w(ix1+1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) +dx1x2*(w(ix1,&
              ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) +dx1x3*(w(ix1,&
              ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-w(ix1,&
              ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=( (w(ix1+1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+w(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,mag(1)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,1)+(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,&
              mag(2))+w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              2)-(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+w(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*block%surfaceC(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,2)+(w(ix1,ixFmin2:ixFmax2,&
              ixFmin3+1:ixFmax3+1,mag(3))+w(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,mag(3)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,3)-(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              mag(3))+w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,&
              mag(3)))*block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,&
              3) )/block%surfaceC(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              1)-w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
      
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
     case(2)
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
       
       
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=w(ix1-1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) -dx1x2*(w(ix1,&
              ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) -dx1x3*(w(ix1,&
              ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-w(ix1,&
              ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=( (w(ix1-1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+w(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,mag(1)))*block%surfaceC(ix1-1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,1)-(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,&
              mag(2))+w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              2)+(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+w(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*block%surfaceC(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,2)-(w(ix1,ixFmin2:ixFmax2,&
              ixFmin3+1:ixFmax3+1,mag(3))+w(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,mag(3)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,3)+(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              mag(3))+w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,&
              mag(3)))*block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,&
              3) )/block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)-w(ix1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
      
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
     case(3)
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
       
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
              ix2+1,ixFmin3:ixFmax3,mag(2)) +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              ixFmin3:ixFmax3,mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,&
              ixFmin3:ixFmax3,mag(1))) +dx2x3*(w(ixFmin1:ixFmax1,ix2,&
              ixFmin3+1:ixFmax3+1,mag(3))-w(ixFmin1:ixFmax1,ix2,&
              ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,&
              mag(2))=( (w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              2)+(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,&
              mag(1))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              1)-(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,&
              1)+(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,&
              mag(3))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              3)-(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(3))+w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,&
              3) )/block%surfaceC(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,&
              2)-w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
      
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
     case(4)
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
       
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
              ix2-1,ixFmin3:ixFmax3,mag(2)) -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              ixFmin3:ixFmax3,mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,&
              ixFmin3:ixFmax3,mag(1))) -dx2x3*(w(ixFmin1:ixFmax1,ix2,&
              ixFmin3+1:ixFmax3+1,mag(3))-w(ixFmin1:ixFmax1,ix2,&
              ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,&
              mag(2))=( (w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,&
              2)-(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,&
              mag(1))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              1)+(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,&
              1)-(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,&
              mag(3))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              3)+(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(3))+w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,&
              3) )/block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              2)-w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
      
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
     
     case(5)
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3+1
       if(slab_uniform) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmax3,ixFmin3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))=w(ixFmin1:ixFmax1,&
              ixFmin2:ixFmax2,ix3+1,mag(3)) +dx3x1*(w(ixFmin1+1:ixFmax1+1,&
              ixFmin2:ixFmax2,ix3,mag(1))-w(ixFmin1-1:ixFmax1-1,&
              ixFmin2:ixFmax2,ix3,mag(1))) +dx3x2*(w(ixFmin1:ixFmax1,&
              ixFmin2+1:ixFmax2+1,ix3,mag(2))-w(ixFmin1:ixFmax1,&
              ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmax3,ixFmin3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,&
              mag(3))=( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,&
              mag(3))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              3)+(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,&
              mag(1))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              1)-(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(1))+w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,&
              1)+(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,&
              mag(2))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              2)-(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(2))+w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,&
              2) )/block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,&
              3)-w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
     case(6)
       if(total_energy) call uawsom_to_primitive(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3-1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))=w(ixFmin1:ixFmax1,&
              ixFmin2:ixFmax2,ix3-1,mag(3)) -dx3x1*(w(ixFmin1+1:ixFmax1+1,&
              ixFmin2:ixFmax2,ix3,mag(1))-w(ixFmin1-1:ixFmax1-1,&
              ixFmin2:ixFmax2,ix3,mag(1))) -dx3x2*(w(ixFmin1:ixFmax1,&
              ixFmin2+1:ixFmax2+1,ix3,mag(2))-w(ixFmin1:ixFmax1,&
              ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,&
              mag(3))=( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,&
              mag(3))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,&
              3)-(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,&
              mag(1))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              1)+(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(1))+w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,&
              1)-(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,&
              mag(2))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              2)+(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(2))+w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,&
              2) )/block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              3)-w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
       if(total_energy) call uawsom_to_conserved(ixGmin1,ixGmin2,ixGmin3,&
          ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,w,x)
    
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  
  subroutine uawsom_clean_divb_multigrid(qdt, qt, active)
    use mod_forest
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_geometry

    double precision, intent(in) :: qdt    !< Current time step
    double precision, intent(in) :: qt     !< Current time
    logical, intent(inout)       :: active !< Output if the source is active
    integer                      :: iigrid, igrid, id
    integer                      :: n, nc, lvl, ixmin1,ixmin2,ixmin3,ixmax1,&
       ixmax2,ixmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idim
    type(tree_node), pointer     :: pnode
    double precision             :: tmp(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3), grad(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3, ndim)
    double precision             :: res
    double precision, parameter  :: max_residual = 1d-3
    double precision, parameter  :: residual_reduction = 1d-10
    integer, parameter           :: max_its      = 50
    double precision             :: residual_it(max_its), max_divb

    mg%operator_type = mg_laplacian

    ! Set boundary conditions
    do n = 1, 2*ndim
       idim = (n+1)/2
       select case (typeboundary(mag(idim), n))
       case (bc_symm)
          ! d/dx B = 0, take phi = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_asymm)
          ! B = 0, so grad(phi) = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_cont)
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_special)
          ! Assume Dirichlet boundary conditions, derivative zero
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_periodic)
          ! Nothing to do here
       case default
          write(*,*)&
              "uawsom_clean_divb_multigrid warning: unknown boundary type"
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       end select
    end do

    ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmin3=ixMlo3-1;ixmax1=ixMhi1+1
    ixmax2=ixMhi2+1;ixmax3=ixMhi3+1;
    max_divb = 0.0d0

    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);

       call get_divb(ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
           1:nw), ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, ixMlo1,ixMlo2,&
          ixMlo3,ixMhi1,ixMhi2,ixMhi3, tmp, uawsom_divb_4thorder)
       mg%boxes(id)%cc(1:nc,1:nc,1:nc, mg_irhs) = tmp(ixMlo1:ixMhi1,&
          ixMlo2:ixMhi2,ixMlo3:ixMhi3)
       max_divb = max(max_divb, maxval(abs(tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
          ixMlo3:ixMhi3))))
    end do

    ! Solve laplacian(phi) = divB
    if(stagger_grid) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, max_divb, 1, MPI_DOUBLE_PRECISION,&
          MPI_MAX, icomm, ierrmpi)

      if (mype == 0) print *, "Performing multigrid divB cleaning"
      if (mype == 0) print *, "iteration vs residual"
      ! Solve laplacian(phi) = divB
      do n = 1, max_its
         call mg_fas_fmg(mg, n>1, max_res=residual_it(n))
         if (mype == 0) write(*, "(I4,E11.3)") n, residual_it(n)
         if (residual_it(n) < residual_reduction * max_divb) exit
      end do
      if (mype == 0 .and. n > max_its) then
         print *, "divb_multigrid warning: not fully converged"
         print *, "current amplitude of divb: ", residual_it(max_its)
         print *, "multigrid smallest grid: ", &
              mg%domain_size_lvl(:, mg%lowest_lvl)
         print *, "note: smallest grid ideally has <= 8 cells"
         print *, "multigrid dx/dy/dz ratio: ", mg%dr(:, 1)/mg%dr(1, 1)
         print *, "note: dx/dy/dz should be similar"
      end if
    else
      do n = 1, max_its
         call mg_fas_vcycle(mg, max_res=res)
         if (res < max_residual) exit
      end do
      if (res > max_residual) call mpistop("divb_multigrid: no convergence")
    end if


    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);

       ! Compute the gradient of phi
       tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = mg%boxes(id)%cc(:,:,:,&
           mg_iphi)

       if(stagger_grid) then
         do idim =1, ndim
           ixCmin1=ixMlo1-kr(idim,1);ixCmin2=ixMlo2-kr(idim,2)
           ixCmin3=ixMlo3-kr(idim,3);
           ixCmax1=ixMhi1;ixCmax2=ixMhi2;ixCmax3=ixMhi3;
           call gradientx(tmp,ps(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
              ixGhi3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim,&
              grad(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,idim),.false.)
           ! Apply the correction B* = B - gradient(phi)
           ps(igrid)%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idim)=ps(igrid)%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3,idim)-grad(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3,idim)
         end do
         ! store cell-center magnetic energy
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3) = sum(ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3, mag(1:ndim))**2, dim=ndim+1)
         ! change cell-center magnetic field
         call uawsom_face_to_center(ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
            ps(igrid))
       else
         do idim = 1, ndim
            call gradient(tmp,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
               ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,idim,grad(ixGlo1:ixGhi1,&
               ixGlo2:ixGhi2,ixGlo3:ixGhi3, idim))
         end do
         ! store cell-center magnetic energy
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3) = sum(ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3, mag(1:ndim))**2, dim=ndim+1)
         ! Apply the correction B* = B - gradient(phi)
         ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
             mag(1:ndim)) = ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3, mag(1:ndim)) - grad(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3, :)
       end if

       if(total_energy) then
         ! Determine magnetic energy difference
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3) = 0.5_dp * (sum(ps(igrid)%w(ixMlo1:ixMhi1,&
            ixMlo2:ixMhi2,ixMlo3:ixMhi3, mag(1:ndim))**2,&
             dim=ndim+1) - tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))
         ! Keep thermal pressure the same
         ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
             e_) = ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
             e_) + tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)
       end if
    end do

    active = .true.

  end subroutine uawsom_clean_divb_multigrid
 

  subroutine uawsom_update_faces(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,wprim,fC,&
     fE,sCT,s,vcts)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3)

    select case(type_ct)
    case('average')
      call update_faces_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fC,fE,&
         sCT,s)
    case('uct_contact')
      call update_faces_contact(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,wprim,&
         fC,fE,sCT,s,vcts)
    case('uct_hll')
      call update_faces_hll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fE,sCT,s,vcts)
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine uawsom_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fC,fE,sCT,&
     s)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3)

    integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,&
       jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3,ixCmmin1,ixCmmin2,ixCmmin3,&
       ixCmmax1,ixCmmax2,ixCmmax3
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2
    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3) :: E_resi, E_ambi

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! if there is resistivity, get eta J
    if(uawsom_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(uawsom_ambipolar_exp) call get_ambipolar_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT%w,x,E_ambi)

    do idim1=1,ndim
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=7-2*ndim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
            ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
            ixCmin3=ixOmin3+kr(idir,3)-1;
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
            jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
            jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmin3=ixCmin3+kr(idim2,3);hxCmax1=ixCmax1+kr(idim2,1)
            hxCmax2=ixCmax2+kr(idim2,2);hxCmax3=ixCmax3+kr(idim2,3);
            ! Interpolate to edges
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=quarter*(fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,iwdim1,idim2)+fC(jxCmin1:jxCmax1,&
               jxCmin2:jxCmax2,jxCmin3:jxCmax3,iwdim1,&
               idim2)-fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               iwdim2,idim1)-fC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
               hxCmin3:hxCmax3,iwdim2,idim1))

            ! add resistive electric field at cell edges E=-vxB+eta J
            if(uawsom_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_resi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)
            ! add ambipolar electric field
            if(uawsom_ambipolar_exp) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_ambi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)

            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=qdt*s%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)*fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)

            if (.not.slab) then
              where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                 r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
            hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
            hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
            ! Add line integrals in direction idir
            circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) > 1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)=bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
    end do

    end associate

  end subroutine update_faces_average

  !> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
  subroutine update_faces_contact(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,wp,fC,fE,&
     sCT,s,vcts)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3)

    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! electric field at cell centers
    double precision                   :: ECC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3)
    ! gradient of E at left and right side of a cell face
    double precision                   :: EL(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),ER(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    ! gradient of E at left and right side of a cell corner
    double precision                   :: ELC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),ERC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3) :: E_resi, E_ambi
    ! total magnetic field at cell centers
    double precision                   :: Btot(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,&
       jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,&
       ixAmax2,ixAmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm)

    if(B0field) then
      Btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)=wp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndim))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim,0)
    else
      Btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)=wp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndim))
    end if
    ECC=0.d0
    ! Calculate electric field at cell centers
    do idim1=1,ndim; do idim2=1,ndim; do idir=7-2*ndim,3
      if(lvc(idim1,idim2,idir)==1)then
         ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)=ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)+Btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idim1)*wp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            mom(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)=ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)-Btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idim1)*wp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            mom(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(uawsom_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(uawsom_ambipolar_exp) call get_ambipolar_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT%w,x,E_ambi)

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
            ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
            ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
            ixCmin3=ixOmin3+kr(idir,3)-1;
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
            jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
            jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmin3=ixCmin3+kr(idim2,3);hxCmax1=ixCmax1+kr(idim2,1)
            hxCmax2=ixCmax2+kr(idim2,2);hxCmax3=ixCmax3+kr(idim2,3);
            ! average cell-face electric field to cell edges
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=quarter*(fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,iwdim1,idim2)+fC(jxCmin1:jxCmax1,&
               jxCmin2:jxCmax2,jxCmin3:jxCmax3,iwdim1,&
               idim2)-fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               iwdim2,idim1)-fC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
               hxCmin3:hxCmax3,iwdim2,idim1))

            ! add slope in idim2 direction from equation (50)
            ixAmin1=ixCmin1;ixAmin2=ixCmin2;ixAmin3=ixCmin3;
            ixAmax1=ixCmax1+kr(idim1,1);ixAmax2=ixCmax2+kr(idim1,2)
            ixAmax3=ixCmax3+kr(idim1,3);
            EL(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=fC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,iwdim1,idim2)-ECC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,ixAmin3:ixAmax3,idir)
            hxCmin1=ixAmin1+kr(idim2,1);hxCmin2=ixAmin2+kr(idim2,2)
            hxCmin3=ixAmin3+kr(idim2,3);hxCmax1=ixAmax1+kr(idim2,1)
            hxCmax2=ixAmax2+kr(idim2,2);hxCmax3=ixAmax3+kr(idim2,3);
            ER(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=fC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,iwdim1,idim2)-ECC(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir)
            where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)>0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=EL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)
            else where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)<0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3)
            else where
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=0.5d0*(EL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)+EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3))
            end where
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmin3=ixCmin3+kr(idim2,3);hxCmax1=ixCmax1+kr(idim2,1)
            hxCmax2=ixCmax2+kr(idim2,2);hxCmax3=ixCmax3+kr(idim2,3);
            where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
               idim1)>0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=ER(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)
            else where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
               idim1)<0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3)
            else where
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=0.5d0*(ER(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)+ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3))
            end where
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)+0.25d0*(ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)+ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3))

            ! add slope in idim1 direction from equation (50)
            jxCmin1=ixCmin1+kr(idim2,1);jxCmin2=ixCmin2+kr(idim2,2)
            jxCmin3=ixCmin3+kr(idim2,3);jxCmax1=ixCmax1+kr(idim2,1)
            jxCmax2=ixCmax2+kr(idim2,2);jxCmax3=ixCmax3+kr(idim2,3);
            ixAmin1=ixCmin1;ixAmin2=ixCmin2;ixAmin3=ixCmin3;
            ixAmax1=ixCmax1+kr(idim2,1);ixAmax2=ixCmax2+kr(idim2,2)
            ixAmax3=ixCmax3+kr(idim2,3);
            EL(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=-fC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,iwdim2,idim1)-ECC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,ixAmin3:ixAmax3,idir)
            hxCmin1=ixAmin1+kr(idim1,1);hxCmin2=ixAmin2+kr(idim1,2)
            hxCmin3=ixAmin3+kr(idim1,3);hxCmax1=ixAmax1+kr(idim1,1)
            hxCmax2=ixAmax2+kr(idim1,2);hxCmax3=ixAmax3+kr(idim1,3);
            ER(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=-fC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,iwdim2,idim1)-ECC(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir)
            where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim2)>0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=EL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)
            else where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim2)<0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3)
            else where
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=0.5d0*(EL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)+EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3))
            end where
            hxCmin1=ixCmin1+kr(idim1,1);hxCmin2=ixCmin2+kr(idim1,2)
            hxCmin3=ixCmin3+kr(idim1,3);hxCmax1=ixCmax1+kr(idim1,1)
            hxCmax2=ixCmax2+kr(idim1,2);hxCmax3=ixCmax3+kr(idim1,3);
            where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
               idim2)>0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=ER(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)
            else where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
               idim2)<0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3)
            else where
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=0.5d0*(ER(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)+ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3))
            end where
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)+0.25d0*(ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)+ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3))

            ! add resistive electric field at cell edges E=-vxB+eta J
            if(uawsom_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_resi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)
            ! add ambipolar electric field
            if(uawsom_ambipolar_exp) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_ambi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)

            ! times time step and edge length
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)*qdt*s%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)
            if (.not.slab) then
              where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                 r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
            hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
            hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
            ! Add line integrals in direction idir
            circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) > 1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)=bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
    end do

    end associate

  end subroutine update_faces_contact

  !> update faces
  subroutine update_faces_hll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt, qdt
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts

    double precision                   :: vtilL(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,2)
    double precision                   :: vtilR(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,2)
    double precision                   :: bfacetot(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    double precision                   :: btilL(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    double precision                   :: btilR(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    double precision                   :: cp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,2)
    double precision                   :: cm(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,2)
    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3) :: E_resi, E_ambi
    integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
       ixCpmin1,ixCpmin2,ixCpmin3,ixCpmax1,ixCpmax2,ixCpmax3,jxCmin1,jxCmin2,&
       jxCmin3,jxCmax1,jxCmax2,jxCmax3,ixCmmin1,ixCmmin2,ixCmmin3,ixCmmax1,&
       ixCmmax2,ixCmmax3
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
    if(uawsom_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(uawsom_ambipolar_exp) call get_ambipolar_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT%w,x,E_ambi)

    do idir=7-2*ndim,3
      ! Indices
      ! idir: electric field component
      ! idim1: one surface
      ! idim2: the other surface
      ! cyclic permutation: idim1,idim2,idir=1,2,3
      ! Velocity components on the surface
      ! follow cyclic premutations:
      ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-1+kr(idir,1);ixCmin2=ixOmin2-1+kr(idir,2)
      ixCmin3=ixOmin3-1+kr(idir,3);

      ! Set indices and directions
      idim1=mod(idir,3)+1
      idim2=mod(idir+1,3)+1

      jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
      jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
      jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
      ixCpmin1=ixCmin1+kr(idim2,1);ixCpmin2=ixCmin2+kr(idim2,2)
      ixCpmin3=ixCmin3+kr(idim2,3);ixCpmax1=ixCmax1+kr(idim2,1)
      ixCpmax2=ixCmax2+kr(idim2,2);ixCpmax3=ixCmax3+kr(idim2,3);

      ! Reconstruct transverse transport velocities
      call reconstruct(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,vbarC(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,idim1,1),vtilL(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2),vtilR(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2))

      call reconstruct(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,vbarC(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,idim2,2),vtilL(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1),vtilR(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1))

      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      if(B0field) then
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1)+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1,idim1)
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2)+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2,idim2)
      else
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1)
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2)
      end if
      call reconstruct(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,&
         bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
         btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
         btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1))

      call reconstruct(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,&
         bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
         btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
         btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2))

      ! Take the maximum characteristic

      cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(cbarmin(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,&
         idim1),cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1))
      cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(cbarmax(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,&
         idim1),cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1))

      cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)=max(cbarmin(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
         cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2))
      cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)=max(cbarmax(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
         cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2))
     

      ! Calculate eletric field
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=-(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim2) + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim2) - cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*(btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim2)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim2)))/(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)+cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)) +(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) - cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*(btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)))/(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)+cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2))

      ! add resistive electric field at cell edges E=-vxB+eta J
      if(uawsom_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)+E_resi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)
      ! add ambipolar electric field
      if(uawsom_ambipolar_exp) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,idir)+E_ambi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,idir)

      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=qdt*s%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)*fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)

      if (.not.slab) then
        where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           r_)+half*dxlevel(r_)).lt.1.0d-9)
          fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=zero
        end where
      end if

    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

    ! Calculate circulation on each face: interal(fE dot dl)
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
            hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
            hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
            ! Add line integrals in direction idir
            circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) > 1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)=bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
    end do

    end associate
  end subroutine update_faces_hll

  !> calculate eta J at cell edges
  subroutine get_resistive_electric_field(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,sCT,s,&
     jce)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    type(state), intent(in)            :: sCT, s
    ! current on cell edges
    double precision :: jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       7-2*ndim:3)

    ! current on cell centers
    double precision :: jcc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       7-2*ndir:3)
    ! location at cell faces
    double precision :: xs(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3,&
       1:ndim)
    ! resistivity
    double precision :: eta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: gradi(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3)
    integer :: ix1,ix2,ix3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,ixBmin1,ixBmin2,ixBmin3,&
       ixBmax1,ixBmax2,ixBmax3,idir,idirmin,idim1,idim2

    associate(x=>s%x,dx=>s%dx,w=>s%w,wCT=>sCT%w,wCTs=>sCT%ws)
    ! calculate current density at cell edges
    jce=0.d0
    do idim1=1,ndim
      do idim2=1,ndim
        do idir=7-2*ndim,3
          if (lvc(idim1,idim2,idir)==0) cycle
          ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
          ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
          ixCmin3=ixOmin3+kr(idir,3)-1;
          ixBmax1=ixCmax1-kr(idir,1)+1;ixBmax2=ixCmax2-kr(idir,2)+1
          ixBmax3=ixCmax3-kr(idir,3)+1;
          ixBmin1=ixCmin1;ixBmin2=ixCmin2;ixBmin3=ixCmin3;
          ! current at transverse faces
          xs(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
             :)=x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,:)
          xs(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
             idim2)=x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
             idim2)+half*dx(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
             idim2)
          call gradientx(wCTs(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3,&
             idim2),xs,ixGslo1,ixGslo2,ixGslo3,ixGshi1,ixGshi2,ixGshi3,ixCmin1,&
             ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,gradi,.true.)
          if (lvc(idim1,idim2,idir)==1) then
            jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)+gradi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
          else
            jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)-gradi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
          end if
        end do
      end do
    end do
    ! get resistivity
    if(uawsom_eta>zero)then
      jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         :)=jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:)*uawsom_eta
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
      ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
      call get_current(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,jcc)
      call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,x,jcc,&
         eta)
      ! calcuate eta on cell edges
      do idir=7-2*ndim,3
        ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
        ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
        ixCmin3=ixOmin3+kr(idir,3)-1;
        jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=0.d0
       do ix3=0,1
       do ix2=0,1
       do ix1=0,1
          if( ix1==1 .and. 1==idir  .or. ix2==1 .and. 2==idir  .or. ix3==1 &
             .and. 3==idir ) cycle
          ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;ixAmin3=ixCmin3+ix3;
          ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;ixAmax3=ixCmax3+ix3;
          jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idir)=jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idir)+eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)
       end do
       end do
       end do
        jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)=jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)*0.25d0
        jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)=jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)*jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)
      enddo
    end if

    end associate
  end subroutine get_resistive_electric_field

  !> get ambipolar electric field on cell edges
  subroutine get_ambipolar_electric_field(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,fE)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)      :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndim:3)

    double precision :: jxbxb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)
    integer :: idir,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ix1,ix2,ix3

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    call uawsom_get_jxbxb(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,jxbxb)
    ! calcuate electric field on cell edges from cell centers
    do idir=7-2*ndim,3
      !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
      !jxbxb(ixA^S,i) = -(uawsom_eta_ambi/w(ixA^S, rho_)**2) * jxbxb(ixA^S,i)
      call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,jxbxb(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,idir),w,x)
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
      ixCmin3=ixOmin3+kr(idir,3)-1;
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=0.d0
     do ix3=0,1
     do ix2=0,1
     do ix1=0,1
        if( ix1==1 .and. 1==idir  .or. ix2==1 .and. 2==idir  .or. ix3==1 .and. &
           3==idir ) cycle
        ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;ixAmin3=ixCmin3+ix3;
        ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;ixAmax3=ixCmax3+ix3;
        fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)+jxbxb(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,idir)
     end do
     end do
     end do
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)*0.25d0
    end do

  end subroutine get_ambipolar_electric_field

  !> calculate cell-center values from face-center values
  subroutine uawsom_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    type(state)                        :: s

    integer                            :: fxOmin1,fxOmin2,fxOmin3,fxOmax1,&
       fxOmax2,fxOmax3, gxOmin1,gxOmin2,gxOmin3,gxOmax1,gxOmax2,gxOmax3,&
        hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, jxOmin1,jxOmin2,&
       jxOmin3,jxOmax1,jxOmax2,jxOmax3, kxOmin1,kxOmin2,kxOmin3,kxOmax1,&
       kxOmax2,kxOmax3, idim

    associate(w=>s%w, ws=>s%ws)

    ! calculate cell-center values from face-center values in 2nd order
    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the w arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
      hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
      hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idim))=half/s%surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)*(ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)*s%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)+ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         hxOmin3:hxOmax3,idim)*s%surfaceC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         hxOmin3:hxOmax3,idim))
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
  subroutine b_from_vector_potential(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
     ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x)
    use mod_global_parameters
    use mod_constrained_transport

    integer, intent(in)                :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
       ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout)    :: ws(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:nws)
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    double precision                   :: Adummy(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:3)

    call b_from_vector_potentialA(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,ixIsmax2,&
       ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x, Adummy)

  end subroutine b_from_vector_potential

end module mod_uawsom_phys
