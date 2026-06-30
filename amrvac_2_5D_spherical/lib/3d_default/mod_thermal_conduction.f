!> Thermal conduction for HD and MHD
!> Adaptation of mod_thermal_conduction for the mod_supertimestepping
!> In order to use it set use_mhd_tc=1 (for the mhd impl) or 2 (for the hd impl) in mhd_list  (for the mhd module both hd and mhd impl can be used)
!> or use_new_hd_tc in hd_list parameters to true
!> (for the hd module, hd implementation has to be used)
!> The TC is set by calling one
!> tc_init_hd_for_total_energy and tc_init_mhd_for_total_energy might
!> The second argument: ixArray has to be [rho_,e_,mag(1)] for mhd (Be aware that the other components of the mag field are assumed consecutive) and [rho_,e_] for hd
!> additionally when internal energy equation is solved, an additional element of this array is eaux_: the index of the internal energy variable.
!>
!> 10.07.2011 developed by Chun Xia and Rony Keppens
!> 01.09.2012 moved to modules folder by Oliver Porth
!> 13.10.2013 optimized further by Chun Xia
!> 12.03.2014 implemented RKL2 super timestepping scheme to reduce iterations
!> and improve stability and accuracy up to second order in time by Chun Xia.
!> 23.08.2014 implemented saturation and perpendicular TC by Chun Xia
!> 12.01.2017 modulized by Chun Xia
!>
!> PURPOSE:
!> IN MHD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(KAPPA_i,j . GRAD_j T)
!> where KAPPA_i,j = tc_k_para b_i b_j + tc_k_perp (I - b_i b_j)
!> b_i b_j = B_i B_j / B**2, I is the unit matrix, and i, j= 1, 2, 3 for 3D
!> IN HD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(tc_k_para . GRAD T)
!> USAGE:
!> 1. in mod_usr.t -> subroutine usr_init(), add
!>        unit_length=your length unit
!>        unit_numberdensity=your number density unit
!>        unit_velocity=your velocity unit
!>        unit_temperature=your temperature unit
!>    before call (m)hd_activate()
!> 2. to switch on thermal conduction in the (m)hd_list of amrvac.par add:
!>    (m)hd_thermal_conduction=.true.
!> 3. in the tc_list of amrvac.par :
!>    tc_perpendicular=.true.  ! (default .false.) turn on thermal conduction perpendicular to magnetic field
!>    tc_saturate=.false.  ! (default .true. ) turn off thermal conduction saturate effect
!>    tc_slope_limiter='MC' ! choose limiter for slope-limited anisotropic thermal conduction in MHD

module mod_thermal_conduction
  use mod_global_parameters, only: std_len
  use mod_geometry
  implicit none

    !> The adiabatic index
    double precision :: tc_gamma_1

  abstract interface
    subroutine get_var_subr(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,nw)
      double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
    end subroutine get_var_subr


  end interface

  type tc_fluid

    procedure (get_var_subr), pointer, nopass :: get_rho => null()
    procedure (get_var_subr), pointer, nopass :: get_rho_equi => null()
    procedure(get_var_subr), pointer,nopass :: get_temperature_from_eint => &
       null()
    procedure(get_var_subr), pointer,nopass :: get_temperature_from_conserved &
       => null()
    procedure(get_var_subr), pointer,nopass :: get_temperature_equi => null()
     !> Indices of the variables
    integer :: e_=-1
    !> Index of cut off temperature for TRAC
    integer :: Tcoff_
    ! if has_equi = .true. get_temperature_equi and get_rho_equi have to be set
    logical :: has_equi=.false.

    ! the following are read from param file or set in tc_read_hd_params or tc_read_mhd_params
    !> Coefficient of thermal conductivity (parallel to magnetic field)
    double precision :: tc_k_para

    !> Coefficient of thermal conductivity perpendicular to magnetic field
    double precision :: tc_k_perp

    !> Name of slope limiter for transverse component of thermal flux
    character(len=std_len)  :: tc_slope_limiter
    !> Logical switch for test constant conductivity
    logical :: tc_constant=.false.
    !> Calculate thermal conduction perpendicular to magnetic field (.true.) or not (.false.)
    logical :: tc_perpendicular=.false.

    !> Consider thermal conduction saturation effect (.true.) or not (.false.)
    logical :: tc_saturate=.true.
    ! END the following are read from param file or set in tc_read_hd_params or tc_read_mhd_params
  end type tc_fluid


  public :: tc_get_mhd_params
  public :: tc_get_hd_params
  public :: get_tc_dt_mhd
  public :: get_tc_dt_hd
  public :: sts_set_source_tc_mhd
  public :: sts_set_source_tc_hd

contains



  subroutine tc_init_params(phys_gamma)
    use mod_global_parameters
    double precision, intent(in) :: phys_gamma

    tc_gamma_1=phys_gamma-1d0
  end subroutine tc_init_params


  !> Init  TC coeffiecients: MHD case
  subroutine tc_get_mhd_params(fl,read_mhd_params)

    use mod_global_parameters

    interface
    subroutine read_mhd_params(fl)
      use mod_global_parameters, only: unitpar,par_files
      import tc_fluid
      type(tc_fluid), intent(inout) :: fl

    end subroutine read_mhd_params
    end interface

    type(tc_fluid), intent(inout) :: fl

    fl%tc_slope_limiter='MC'

    fl%tc_k_para=0.d0

    fl%tc_k_perp=0.d0

    call read_mhd_params(fl)

    if(fl%tc_k_para==0.d0 .and. fl%tc_k_perp==0.d0) then
      if(SI_unit) then
        ! Spitzer thermal conductivity with SI units
        fl%tc_k_para=8.d-12*unit_temperature**&
           3.5d0/unit_length/unit_density/unit_velocity**3
        ! thermal conductivity perpendicular to magnetic field
        fl%tc_k_perp=4.d-30*unit_numberdensity**2/unit_magneticfield**&
           2/unit_temperature**3*fl%tc_k_para
      else
        ! Spitzer thermal conductivity with cgs units
        fl%tc_k_para=8.d-7*unit_temperature**&
           3.5d0/unit_length/unit_density/unit_velocity**3
        ! thermal conductivity perpendicular to magnetic field
        fl%tc_k_perp=4.d-10*unit_numberdensity**2/unit_magneticfield**&
           2/unit_temperature**3*fl%tc_k_para
      end if
      if(mype .eq. 0) print*, "Spitzer MHD: par: ",fl%tc_k_para, " ,perp: ",&
         fl%tc_k_perp
    else
      fl%tc_constant=.true.
    end if

    contains

    !> Read tc module parameters from par file: MHD case

  end subroutine tc_get_mhd_params


  subroutine tc_get_hd_params(fl,read_hd_params)
    use mod_global_parameters

    interface
      subroutine read_hd_params(fl)
        use mod_global_parameters, only: unitpar,par_files
        import tc_fluid
        type(tc_fluid), intent(inout) :: fl

      end subroutine read_hd_params
    end interface
    type(tc_fluid), intent(inout) :: fl

    fl%tc_k_para=0.d0
    !> Read tc parameters from par file: HD case
    call read_hd_params(fl)
    if(fl%tc_k_para==0.d0 ) then
      if(SI_unit) then
        ! Spitzer thermal conductivity with SI units
        fl%tc_k_para=8.d-12*unit_temperature**&
           3.5d0/unit_length/unit_density/unit_velocity**3
      else
        ! Spitzer thermal conductivity with cgs units
        fl%tc_k_para=8.d-7*unit_temperature**&
           3.5d0/unit_length/unit_density/unit_velocity**3
      end if
      if(mype .eq. 0) print*, "Spitzer HD par: ",fl%tc_k_para
    end if
  end subroutine tc_get_hd_params

  !> Get the explicut timestep for the TC (mhd implementation)
  function get_tc_dt_mhd(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,x,&
     fl) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters

    type(tc_fluid), intent(in)  ::  fl
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: dtnew

    double precision :: dxinv(1:ndim),mf(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    double precision :: tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),Te(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3),B2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: dtdiff_tcond,maxtmp2
    integer          :: idim,ix1,ix2,ix3

    dxinv(1)=one/dx1;dxinv(2)=one/dx2;dxinv(3)=one/dx3;

    !temperature
    call fl%get_temperature_from_conserved(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Te)

    !tc_k_para_i
    if(fl%tc_constant) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=fl%tc_k_para
    else
      call fl%get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp2)
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=fl%tc_k_para*dsqrt(Te(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)**5)/tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if

    ! B
    if(B0field) then
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mag(:))+block%B0(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0)
    else
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mag(:))
    end if
    ! B^-2
    B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sum(mf(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,dim=ndim+1)
    ! B_i**2/B**2
    where(B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/=0.d0)
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)=mf(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2/B2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)=mf(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2/B2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)=mf(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)**2/B2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3);
    elsewhere
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)=1.d0
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)=1.d0
      mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)=1.d0;
    end where

    if(fl%tc_saturate) B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=22.d0*dsqrt(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    dtnew=bigdouble
    do idim=1,ndim
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*mf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)
      if(fl%tc_saturate) then
        where(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)>B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        end where
      end if
      maxtmp2=maxval(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      if(maxtmp2==0.d0) maxtmp2=smalldouble
      ! dt< dx_idim**2/((gamma-1)*tc_k_para_i/rho*B_i**2/B**2)
      dtdiff_tcond=1.d0/tc_gamma_1/(maxtmp2*dxinv(idim)**2)
      ! limit the time step
      dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)

  end  function get_tc_dt_mhd

  !> anisotropic thermal conduction with slope limited symmetric scheme
  !> Sharma 2007 Journal of Computational Physics 227, 123
  subroutine sts_set_source_tc_mhd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
     fix_conserve_at_step,my_dt,igrid,nflux,fl)
    use mod_global_parameters
    use mod_fix_conserve
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, igrid, nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in) ::   w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    type(tc_fluid), intent(in) :: fl

    !! qd store the heat conduction energy changing rate
    double precision :: qd(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       Te(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, allocatable, dimension(:,:,:,:) :: qvec_equi

    double precision, allocatable, dimension(:,:,:,:,:) :: fluxall
    double precision :: alpha,dxinv(ndim)
    integer :: idims,idir,ix1,ix2,ix3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixAmin1,ixAmin2,&
       ixAmin3,ixAmax1,ixAmax2,ixAmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,&
       ixBmax3

    ! coefficient of limiting on normal component
    if(ndim<3) then
      alpha=0.75d0
    else
      alpha=0.85d0
    end if
    ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
    ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;

    dxinv=1.d0/dxlevel

    call fl%get_temperature_from_eint(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, Te) !calculate Te in whole domain (+ghosts)
    call fl%get_rho(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, rho) !calculate rho in whole domain (+ghosts)
    call set_source_tc_mhd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,fl,qvec,rho,Te,&
       alpha)
    if(fl%has_equi) then
      allocate(qvec_equi(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim))
      call fl%get_temperature_equi(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, Te) !calculate Te in whole domain (+ghosts)
      call fl%get_rho_equi(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, rho) !calculate rho in whole domain (+ghosts)
      call set_source_tc_mhd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,fl,qvec_equi,rho,&
         Te,alpha)
      do idims=1,ndim
        qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)=qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims) - qvec_equi(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idims)
      end do
      deallocate(qvec_equi)
    endif

    qd=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)=dxinv(idims)*qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        qd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=qd(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)+qvec(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,idims)-qvec(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,ixBmin3:ixBmax3,idims)
      end do
    else
      do idims=1,ndim
        qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)=qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)*block%surfaceC(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        qd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=qd(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)+qvec(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,idims)-qvec(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,ixBmin3:ixBmax3,idims)
      end do
      qd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=qd(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%dvolume(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if

    if(fix_conserve_at_step) then
      allocate(fluxall(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1,&
         1:ndim))
      fluxall(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1,&
         1:ndim)=my_dt*qvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)
      call store_flux(igrid,fluxall,1,ndim,nflux)
      deallocate(fluxall)
    end if

    wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       fl%e_)=qd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    !do ix1 = ixImin1,ixImax1
    !  if (abs(wres(ix1,fl%e_)) > 10.d0) then
        !if (mype == 7) then      
        !  write(*,*) "it                 = ", it
        !  write(*,*) "TC turned off at x = ", x(ix1,1)
        !  write(*,*) "TC value was       = ", wres(ix1,fl%e_)
        !end if
    !    wres(ix1,fl%e_) = 0.d0
    !  end if 
    !end do

    !if (mod(it,5000) == 0 .and. mype == 0) then
    !  write(*,*) "max TC first block = ", maxval(wres(ixO^S,fl%e_))
    !  write(*,*) "min TC first block = ", minval(wres(ixO^S,fl%e_))
    !end if

    !if (mod(it,5000) == 0 .and. mype == 7) then
    !  write(*,*) "max TC last block = ", maxval(wres(ixO^S,fl%e_))
    !  write(*,*) "min TC last block = ", minval(wres(ixO^S,fl%e_))
    !end if

    !do ix1 = ixOmin1,ixOmax1
    !  if (abs(wres(ix1,fl%e_)) > 1.d0) then
    !    wres(ix1,fl%e_) = 0.d0
    !  end if
    !end do


  end subroutine sts_set_source_tc_mhd

  subroutine set_source_tc_mhd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,fl,qvec,rho,Te,alpha)
    use mod_global_parameters
    use mod_fix_conserve
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in) ::  w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    type(tc_fluid), intent(in) :: fl
    double precision, intent(in) :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),Te(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(in) :: alpha
    double precision, intent(out) :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    !! qd store the heat conduction energy changing rate
    double precision :: qd(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir) :: mf,Bc,Bcf
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) :: gradT
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: ka,kaf,ke,kef,qdd,qe,Binv,minq,maxq,Bnorm
    double precision, allocatable, dimension(:,:,:,:,:) :: fluxall
    integer, dimension(ndim) :: lowindex
    integer :: idims,idir,ix1,ix2,ix3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixAmin1,ixAmin2,&
       ixAmin3,ixAmax1,ixAmax2,ixAmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,&
       ixBmax3

    ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
    ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;

    ! T gradient at cell faces
    ! B vector
    if(B0field) then
      mf(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,iw_mag(:))+block%B0(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,:,0)
    else
      mf(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,iw_mag(:))
    end if
    ! |B|
    Binv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dsqrt(sum(mf(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,:)**2,dim=ndim+1))
    where(Binv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)/=0.d0)
      Binv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=1.d0/Binv(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)
    elsewhere
      Binv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=bigdouble
    end where
    ! b unit vector: magnetic field direction vector
    do idims=1,ndim
      mf(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idims)=mf(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3,idims)*Binv(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)
    end do
    ! ixC is cell-corner index
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=ixOmin1-1
    ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;
    ! b unit vector at cell corner
    Bc=0.d0
    do ix3=0,1
    do ix2=0,1
    do ix1=0,1
      ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;ixAmin3=ixCmin3+ix3;
      ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;ixAmax3=ixCmax3+ix3;
      Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:ndim)=Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:ndim)+mf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndim)
    end do
    end do
    end do
    Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       1:ndim)=Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       1:ndim)*0.5d0**ndim
    ! T gradient at cell faces
    gradT=0.d0
    do idims=1,ndim
      ixBmin1=ixmin1;ixBmin2=ixmin2;ixBmin3=ixmin3;
      ixBmax1=ixmax1-kr(idims,1);ixBmax2=ixmax2-kr(idims,2)
      ixBmax3=ixmax3-kr(idims,3);
      call gradientC(Te,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3,idims,minq)
      gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         idims)=minq(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3)
    end do
    if(fl%tc_constant) then
      if(fl%tc_perpendicular) then
        ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=fl%tc_k_para-fl%tc_k_perp
        ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=fl%tc_k_perp
      else
        ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=fl%tc_k_para
      end if
    else
      ! conductivity at cell center
      if(phys_trac) then
        minq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=Te(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)
        where(minq(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3) < block%wextra(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3,fl%Tcoff_))
          minq(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=block%wextra(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,fl%Tcoff_)
        end where
        minq(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)=fl%tc_k_para*sqrt(minq(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)**5)
      else
        minq(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)=fl%tc_k_para*sqrt(Te(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)**5)
      end if
      ka=0.d0
      do ix3=0,1
      do ix2=0,1
      do ix1=0,1
        ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;ixBmin3=ixCmin3+ix3;
        ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;ixBmax3=ixCmax3+ix3;
        ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=ka(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3)+minq(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,ixBmin3:ixBmax3)
      end do
      end do
      end do
      ! cell corner conductivity
      ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=0.5d0**ndim*ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      ! compensate with perpendicular conductivity
      if(fl%tc_perpendicular) then
        minq(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)=fl%tc_k_perp*rho(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)**2*Binv(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)**2/dsqrt(Te(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3))
        ke=0.d0
        do ix3=0,1
        do ix2=0,1
        do ix1=0,1
          ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;ixBmin3=ixCmin3+ix3;
          ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;ixBmax3=ixCmax3+ix3;
          ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)+minq(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
             ixBmin3:ixBmax3)
        end do
        end do
        end do
        ! cell corner conductivity: k_parallel-k_perpendicular
        ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=0.5d0**ndim*ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        where(ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)<ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3))
          ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)-ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)
        elsewhere
          ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=0.d0
          ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)
        end where
      end if
    end if
    if(fl%tc_slope_limiter=='no') then
      ! calculate thermal conduction flux with symmetric scheme
      do idims=1,ndim
        !qd corner values
        qd=0.d0
        do ix3=0,1 
        do ix2=0,1 
        do ix1=0,1 
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims  .or. ix3==0 &
              .and. 3==idims ) then
             ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;ixBmin3=ixCmin3+ix3;
             ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;ixBmax3=ixCmax3+ix3;
             qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                ixCmin3:ixCmax3)=qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                ixCmin3:ixCmax3)+gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                ixBmin3:ixBmax3,idims)
           end if
        end do
        end do
        end do
        ! temperature gradient at cell corner
        qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idims)=qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*0.5d0**(ndim-1)
      end do
      ! b grad T at cell corner
      qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=sum(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,1:ndim)*Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,1:ndim),dim=ndim+1)
      do idims=1,ndim
        ! TC flux at cell corner
        gradT(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idims)=ka(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*Bc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idims)*qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
        if(fl%tc_perpendicular) gradT(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idims)=gradT(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idims)+ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idims)
      end do
      ! TC flux at cell face
      qvec=0.d0
      do idims=1,ndim
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        ixAmax1=ixOmax1;ixAmax2=ixOmax2;ixAmax3=ixOmax3; ixAmin1=ixBmin1
        ixAmin2=ixBmin2;ixAmin3=ixBmin3;
        do ix3=0,1 
        do ix2=0,1 
        do ix1=0,1 
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims  .or. ix3==0 &
              .and. 3==idims ) then
             ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2;ixBmin3=ixAmin3-ix3;
             ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2;ixBmax3=ixAmax3-ix3;
             qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
                idims)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
                idims)+gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
                idims)
           end if
        end do
        end do
        end do
        qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idims)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idims)*0.5d0**(ndim-1)
        if(fl%tc_saturate) then
          ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
          ! unsigned saturated TC flux = 5 phi rho c**3, c is isothermal sound speed
          Bcf=0.d0
          do ix3=0,1 
          do ix2=0,1 
          do ix1=0,1 
             if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims  .or. &
                ix3==0 .and. 3==idims ) then
               ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2;ixBmin3=ixAmin3-ix3;
               ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2;ixBmax3=ixAmax3-ix3;
               Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
                  idims)=Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
                  idims)+Bc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
                  idims)
             end if
          end do
          end do
          end do
          ! averaged b at face centers
          Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
             idims)=Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
             idims)*0.5d0**(ndim-1)
          ixBmin1=ixAmin1+kr(idims,1);ixBmin2=ixAmin2+kr(idims,2)
          ixBmin3=ixAmin3+kr(idims,3);ixBmax1=ixAmax1+kr(idims,1)
          ixBmax2=ixAmax2+kr(idims,2);ixBmax3=ixAmax3+kr(idims,3);
          qd(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3)=2.75d0*(rho(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3)+rho(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
             ixBmin3:ixBmax3))*dsqrt(0.5d0*(Te(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3)+Te(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
             ixBmin3:ixBmax3)))**3*dabs(Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3,idims))
         do ix3=ixAmin3,ixAmax3
         do ix2=ixAmin2,ixAmax2
         do ix1=ixAmin1,ixAmax1
            if(dabs(qvec(ix1,ix2,ix3,idims))>qd(ix1,ix2,ix3)) then
              qvec(ix1,ix2,ix3,idims)=sign(1.d0,qvec(ix1,ix2,ix3,&
                 idims))*qd(ix1,ix2,ix3)
            end if
         end do
         end do
         end do
        end if
      end do
    else
      ! calculate thermal conduction flux with slope-limited symmetric scheme
      qvec=0.d0
      do idims=1,ndim
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        ixAmax1=ixOmax1;ixAmax2=ixOmax2;ixAmax3=ixOmax3; ixAmin1=ixBmin1
        ixAmin2=ixBmin2;ixAmin3=ixBmin3;
        ! calculate normal of magnetic field
        ixBmin1=ixAmin1+kr(idims,1);ixBmin2=ixAmin2+kr(idims,2)
        ixBmin3=ixAmin3+kr(idims,3);ixBmax1=ixAmax1+kr(idims,1)
        ixBmax2=ixAmax2+kr(idims,2);ixBmax3=ixAmax3+kr(idims,3);
        Bnorm(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)=0.5d0*(mf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)+mf(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
           ixBmin3:ixBmax3,idims))
        Bcf=0.d0
        kaf=0.d0
        kef=0.d0
        do ix3=0,1 
        do ix2=0,1 
        do ix1=0,1 
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims  .or. ix3==0 &
              .and. 3==idims ) then
             ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2;ixBmin3=ixAmin3-ix3;
             ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2;ixBmax3=ixAmax3-ix3;
             Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
                1:ndim)=Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
                1:ndim)+Bc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
                1:ndim)
             kaf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                ixAmin3:ixAmax3)=kaf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                ixAmin3:ixAmax3)+ka(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                ixBmin3:ixBmax3)
             if(fl%tc_perpendicular) kef(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                ixAmin3:ixAmax3)=kef(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                ixAmin3:ixAmax3)+ke(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                ixBmin3:ixBmax3)
           end if
        end do
        end do
        end do
        ! averaged b at face centers
        Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           1:ndim)=Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           1:ndim)*0.5d0**(ndim-1)
        ! averaged thermal conductivity at face centers
        kaf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)=kaf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)*0.5d0**(ndim-1)
        if(fl%tc_perpendicular) kef(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)=kef(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)*0.5d0**(ndim-1)
        ! limited normal component
        minq(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)=min(alpha*gradT(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims),gradT(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)/alpha)
        maxq(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)=max(alpha*gradT(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims),gradT(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)/alpha)
        ! eq (19)
        !corner values of gradT
        qdd=0.d0
        do ix3=0,1 
        do ix2=0,1 
        do ix1=0,1 
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims  .or. ix3==0 &
              .and. 3==idims ) then
             ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;ixBmin3=ixCmin3+ix3;
             ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;ixBmax3=ixCmax3+ix3;
             qdd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                ixCmin3:ixCmax3)=qdd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                ixCmin3:ixCmax3)+gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                ixBmin3:ixBmax3,idims)
           end if
        end do
        end do
        end do
        ! temperature gradient at cell corner
        qdd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=qdd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*0.5d0**(ndim-1)
        ! eq (21)
        qe=0.d0
        do ix3=0,1 
        do ix2=0,1 
        do ix1=0,1 
           qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3)=qdd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3)
           if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims  .or. ix3==0 &
              .and. 3==idims ) then
             ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2;ixBmin3=ixAmin3-ix3;
             ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2;ixBmax3=ixAmax3-ix3;
             where(qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                ixBmin3:ixBmax3)<=minq(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                ixAmin3:ixAmax3))
               qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                  ixBmin3:ixBmax3)=minq(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                  ixAmin3:ixAmax3)
             elsewhere(qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                ixBmin3:ixBmax3)>=maxq(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                ixAmin3:ixAmax3))
               qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                  ixBmin3:ixBmax3)=maxq(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                  ixAmin3:ixAmax3)
             end where
             qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
                idims)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
                idims)+Bc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
                idims)**2*qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3)
             if(fl%tc_perpendicular) qe(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                ixAmin3:ixAmax3)=qe(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
                ixAmin3:ixAmax3)+qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
                ixBmin3:ixBmax3)
           end if
        end do
        end do
        end do
        qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idims)=kaf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)*qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)*0.5d0**(ndim-1)
        ! add normal flux from perpendicular conduction
        if(fl%tc_perpendicular) qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)+kef(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)*qe(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3)*0.5d0**(ndim-1)
        ! limited transverse component, eq (17)
        ixBmin1=ixAmin1;ixBmin2=ixAmin2;ixBmin3=ixAmin3;
        ixBmax1=ixAmax1+kr(idims,1);ixBmax2=ixAmax2+kr(idims,2)
        ixBmax3=ixAmax3+kr(idims,3);
        do idir=1,ndim
          if(idir==idims) cycle
          qd(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)=slope_limiter(gradT(ixImin1:ixImax1,&
             ixImin2:ixImax2,ixImin3:ixImax3,idir),ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,&
             ixBmax3,idir,-1,fl%tc_slope_limiter)
          qd(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=slope_limiter(qd,&
             ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixAmin1,ixAmin2,&
             ixAmin3,ixAmax1,ixAmax2,ixAmax3,idims,1,fl%tc_slope_limiter)
          qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
             idims)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
             idims)+kaf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3)*Bnorm(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3)*Bcf(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3,idir)*qd(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3)
        end do

        ! consider magnetic null point
        !where(Binv(ixA^S)==0.d0)
        !  qvec(ixA^S,idims)=tc_k_para*(0.5d0*(Te(ixA^S)+Te(ixB^S)))**2.5d0*gradT(ixA^S,idims)
        !end where

        if(fl%tc_saturate) then
          ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
          ! unsigned saturated TC flux = 5 phi rho c**3, c is isothermal sound speed
          ixBmin1=ixAmin1+kr(idims,1);ixBmin2=ixAmin2+kr(idims,2)
          ixBmin3=ixAmin3+kr(idims,3);ixBmax1=ixAmax1+kr(idims,1)
          ixBmax2=ixAmax2+kr(idims,2);ixBmax3=ixAmax3+kr(idims,3);
          qd(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3)=2.75d0*(rho(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3)+rho(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
             ixBmin3:ixBmax3))*dsqrt(0.5d0*(Te(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3)+Te(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
             ixBmin3:ixBmax3)))**3*dabs(Bnorm(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
             ixAmin3:ixAmax3))
         do ix3=ixAmin3,ixAmax3
         do ix2=ixAmin2,ixAmax2
         do ix1=ixAmin1,ixAmax1
            if(dabs(qvec(ix1,ix2,ix3,idims))>qd(ix1,ix2,ix3)) then
        !      write(*,*) 'it',it,qvec(ix^D,idims),qd(ix^D),' TC saturated at ',&
        !      x(ix^D,:),' rho',rho(ix^D),' Te',Te(ix^D)
              qvec(ix1,ix2,ix3,idims)=sign(1.d0,qvec(ix1,ix2,ix3,&
                 idims))*qd(ix1,ix2,ix3)
            end if
         end do
         end do
         end do
        end if
      end do
    end if
  end subroutine set_source_tc_mhd

  function slope_limiter(f,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,pm,&
     tc_slope_limiter) result(lf)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idims, pm
    double precision, intent(in) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: lf(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    character(len=*), intent(in)  :: tc_slope_limiter

    double precision :: signf(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3

    ixBmin1=ixOmin1+pm*kr(idims,1);ixBmin2=ixOmin2+pm*kr(idims,2)
    ixBmin3=ixOmin3+pm*kr(idims,3);ixBmax1=ixOmax1+pm*kr(idims,1)
    ixBmax2=ixOmax2+pm*kr(idims,2);ixBmax3=ixOmax3+pm*kr(idims,3);
    signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sign(1.d0,&
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    select case(tc_slope_limiter)
     case('minmod')
       ! minmod limiter
       lf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*max(0.d0,min(abs(f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)),signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
          ixBmin3:ixBmax3)))
     case ('MC')
       ! montonized central limiter Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       lf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=two*signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)* max(zero,min(dabs(f(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)),signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
          ixBmin3:ixBmax3),signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*quarter*(f(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
          ixBmin3:ixBmax3)+f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))))
     case ('superbee')
       ! Roes superbee limiter (eq.3.51i)
       lf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)* max(zero,min(two*dabs(f(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)),signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
          ixBmin3:ixBmax3)),min(dabs(f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)),two*signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
          ixBmin3:ixBmax3)))
     case ('koren')
       ! Barry Koren Right variant
       lf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)* max(zero,min(two*dabs(f(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)),two*signf(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
          ixBmin3:ixBmax3),(two*f(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
          ixBmin3:ixBmax3)*signf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+dabs(f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))*third))
     case default
       call mpistop("Unknown slope limiter for thermal conduction")
    end select

  end function slope_limiter

  !> Calculate gradient of a scalar q at cell interfaces in direction idir
  subroutine gradientC(q,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
    double precision, intent(in)    :: q(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(inout) :: gradq(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                         :: jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,&
       jxOmax3

    associate(x=>block%x)

    jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
    jxOmin3=ixOmin3+kr(idir,3);jxOmax1=ixOmax1+kr(idir,1)
    jxOmax2=ixOmax2+kr(idir,2);jxOmax3=ixOmax3+kr(idir,3);
    select case(coordinate)
    case(Cartesian)
      gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
         jxOmin3:jxOmax3)-q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/dxlevel(idir)
    case(Cartesian_stretched)
      gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
         jxOmin3:jxOmax3)-q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/(x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
         idir)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))/(x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,1)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1))
        
      case(2)
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))/( (x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,2)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,2))*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1) )
       
        
      case(3)
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))/( (x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,3)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,3))*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,2)) )
       
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))/((x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,phi_)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,phi_))*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,r_))
      else
        gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3)-q(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))/(x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,idir)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

    end associate
  end subroutine gradientC

  function get_tc_dt_hd(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,x,&
     fl)  result(dtnew)
    ! Check diffusion time limit dt < dx_i**2 / ((gamma-1)*tc_k_para_i/rho)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    type(tc_fluid), intent(in) :: fl
    double precision :: dtnew

    double precision :: dxinv(1:ndim), tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), Te(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: dtdiff_tcond,dtdiff_tsat
    integer          :: idim,ix1,ix2,ix3

    dxinv(1)=one/dx1;dxinv(2)=one/dx2;dxinv(3)=one/dx3;

    call fl%get_temperature_from_conserved(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Te)
    call fl%get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)

    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)&
       =tc_gamma_1*fl%tc_k_para*dsqrt((Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))**5)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    dtnew = bigdouble

    do idim=1,ndim
       ! dt< dx_idim**2/((gamma-1)*tc_k_para_idim/rho)
       dtdiff_tcond=1d0/maxval(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*dxinv(idim)**2)
       if(fl%tc_saturate) then
         ! dt< dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
         dtdiff_tsat=1d0/maxval(tc_gamma_1*dsqrt(Te(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3))*5.d0*dxinv(idim)**2)
         ! choose the slower flux (bigger time scale) between classic and saturated
         dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
       end if
       ! limit the time step
       dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)

  end function get_tc_dt_hd

  subroutine sts_set_source_tc_hd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
     fix_conserve_at_step,my_dt,igrid,nflux,fl)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, igrid, nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in) ::  w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    type(tc_fluid), intent(in)    :: fl

    double precision :: Te(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim),qd(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, allocatable, dimension(:,:,:,:) :: qvec_equi
    double precision, allocatable, dimension(:,:,:,:,:) :: fluxall

    double precision :: dxinv(ndim)
    integer :: idims,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,ixBmin1,ixBmin2,&
       ixBmin3,ixBmax1,ixBmax2,ixBmax3

    ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
    ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;

    dxinv=1.d0/dxlevel

    !calculate Te in whole domain (+ghosts)
    call fl%get_temperature_from_eint(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, Te)
    call fl%get_rho(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, rho)
    call set_source_tc_hd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,fl,qvec,rho,Te)
    if(fl%has_equi) then
      allocate(qvec_equi(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim))
      call fl%get_temperature_equi(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, Te) !calculate Te in whole domain (+ghosts)
      call fl%get_rho_equi(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, rho) !calculate rho in whole domain (+ghosts)
      call set_source_tc_hd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,fl,qvec_equi,rho,&
         Te)
      do idims=1,ndim
        qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)=qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims) - qvec_equi(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idims)
      end do
      deallocate(qvec_equi)
    endif


    qd=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)=dxinv(idims)*qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        qd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=qd(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)+qvec(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,idims)-qvec(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,ixBmin3:ixBmax3,idims)
      end do
    else
      do idims=1,ndim
        qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)=qvec(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)*block%surfaceC(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        qd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=qd(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)+qvec(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,idims)-qvec(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,ixBmin3:ixBmax3,idims)
      end do
      qd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=qd(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%dvolume(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if

    if(fix_conserve_at_step) then
      allocate(fluxall(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1,&
         1:ndim))
      fluxall(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1,&
         1:ndim)=my_dt*qvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)
      call store_flux(igrid,fluxall,1,ndim,nflux)
      deallocate(fluxall)
    end if
    wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       fl%e_)=qd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine sts_set_source_tc_hd

  subroutine set_source_tc_hd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,fl,qvec,rho,Te)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in) ::  w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    type(tc_fluid), intent(in)    :: fl
    double precision, intent(in) :: Te(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(out) :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)


    double precision :: gradT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim),ke(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       qd(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    double precision :: dxinv(ndim)
    integer :: idims,ix1,ix2,ix3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
       ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixAmin1,ixAmin2,ixAmin3,&
       ixAmax1,ixAmax2,ixAmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3,&
       ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3

    ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
    ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;
    ! ixC is cell-corner index
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=ixOmin1-1
    ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;

    dxinv=1.d0/dxlevel

    !calculate Te in whole domain (+ghosts)
    ! cell corner temperature in ke
    ke=0.d0
    ixAmax1=ixmax1;ixAmax2=ixmax2;ixAmax3=ixmax3; ixAmin1=ixmin1-1
    ixAmin2=ixmin2-1;ixAmin3=ixmin3-1;
    do ix3=0,1
    do ix2=0,1
    do ix1=0,1
      ixBmin1=ixAmin1+ix1;ixBmin2=ixAmin2+ix2;ixBmin3=ixAmin3+ix3;
      ixBmax1=ixAmax1+ix1;ixBmax2=ixAmax2+ix2;ixBmax3=ixAmax3+ix3;
      ke(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)=ke(ixAmin1:ixAmax1,&
         ixAmin2:ixAmax2,ixAmin3:ixAmax3)+Te(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         ixBmin3:ixBmax3)
    end do
    end do
    end do
    ke(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)=0.5d0**&
       ndim*ke(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)
    ! T gradient (central difference) at cell corners
    gradT=0.d0
    do idims=1,ndim
      ixBmin1=ixmin1;ixBmin2=ixmin2;ixBmin3=ixmin3;
      ixBmax1=ixmax1-kr(idims,1);ixBmax2=ixmax2-kr(idims,2)
      ixBmax3=ixmax3-kr(idims,3);
      call gradient(ke,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixBmin1,&
         ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3,idims,qd)
      gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         idims)=qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3)
    end do
    ! transition region adaptive conduction
    if(phys_trac) then
      where(ke(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) < block%wextra(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,fl%Tcoff_))
        ke(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)=block%wextra(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,fl%Tcoff_)
      end where
    end if
    ! cell corner conduction flux
    do idims=1,ndim
      gradT(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idims)=gradT(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idims)*fl%tc_k_para*sqrt(ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)**5)
    end do

    if(fl%tc_saturate) then
      ! consider saturation with unsigned saturated TC flux = 5 phi rho c**3
      ! saturation flux at cell center
      qd(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=5.d0*rho(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)*dsqrt(Te(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)**3)
      !cell corner values of qd in ke
      ke=0.d0
      do ix3=0,1
      do ix2=0,1
      do ix1=0,1
        ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;ixBmin3=ixCmin3+ix3;
        ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;ixBmax3=ixCmax3+ix3;
        ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=ke(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3)+qd(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
           ixBmin3:ixBmax3)
      end do
      end do
      end do
      ! cell corner saturation flux
      ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=0.5d0**ndim*ke(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      ! magnitude of cell corner conduction flux
      qd(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=norm2(gradT(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,:),dim=ndim+1)
      do ix3=ixCmin3,ixCmax3
      do ix2=ixCmin2,ixCmax2
      do ix1=ixCmin1,ixCmax1
        if(qd(ix1,ix2,ix3)>ke(ix1,ix2,ix3)) then
          ke(ix1,ix2,ix3)=ke(ix1,ix2,ix3)/qd(ix1,ix2,ix3)
          do idims=1,ndim
            gradT(ix1,ix2,ix3,idims)=ke(ix1,ix2,ix3)*gradT(ix1,ix2,ix3,idims)
          end do
        end if
      end do
      end do
      end do
    end if

    ! conductionflux at cell face
    !face center values of gradT in qvec
    qvec=0.d0
    do idims=1,ndim
      ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
      ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
      ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
      ixAmax1=ixOmax1;ixAmax2=ixOmax2;ixAmax3=ixOmax3; ixAmin1=ixBmin1
      ixAmin2=ixBmin2;ixAmin3=ixBmin3;
      do ix3=0,1 
      do ix2=0,1 
      do ix1=0,1 
         if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims  .or. ix3==0 &
            .and. 3==idims ) then
           ixBmin1=ixAmin1-ix1;ixBmin2=ixAmin2-ix2;ixBmin3=ixAmin3-ix3;
           ixBmax1=ixAmax1-ix1;ixBmax2=ixAmax2-ix2;ixBmax3=ixAmax3-ix3;
           qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
              idims)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
              idims)+gradT(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
              idims)
         end if
      end do
      end do
      end do
      qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         idims)=qvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         idims)*0.5d0**(ndim-1)
    end do


  end subroutine set_source_tc_hd




end module mod_thermal_conduction
