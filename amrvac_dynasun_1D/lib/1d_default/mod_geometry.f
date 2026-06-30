!> Module with geometry-related routines (e.g., divergence, curl)
module mod_geometry
  implicit none
  public

  integer :: coordinate=-1
  integer, parameter :: Cartesian          = 0
  integer, parameter :: Cartesian_stretched= 1
  integer, parameter :: cylindrical        = 2
  integer, parameter :: spherical          = 3
  integer, parameter :: Cartesian_expansion= 4

  integer :: type_curl=0
  integer, parameter :: central=1
  integer, parameter :: Gaussbased=2
  integer, parameter :: Stokesbased=3

contains

  !> Set the coordinate system to be used
  subroutine set_coordinate_system(geom)
    use mod_global_parameters
    character(len=*), intent(in) :: geom !< Name of the coordinate system

    ! Store the geometry name
    geometry_name = geom

    select case (geom)
    case ("Cartesian","Cartesian_1D","Cartesian_2D","Cartesian_3D")
      ndir = ndim
      coordinate=Cartesian
    case ("Cartesian_1D_expansion")
      if (ndim /= 1) call mpistop&
         ("Geometry Cartesian_1D_expansion but ndim /= 1")
      ndir = ndim
      coordinate=Cartesian_expansion
    case ("Cartesian_1.5D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.5D but ndim /= 1")
      coordinate=Cartesian
      ndir = 2
    case ("Cartesian_1.75D")
      if (ndim /= 1) call mpistop("Geometry Cartesian_1.75D but ndim /= 1")
      coordinate=Cartesian
      ndir = 3
    case ("Cartesian_2.5D")
      if (ndim /= 2) call mpistop("Geometry Cartesian_2.5D but ndim /= 2")
      coordinate=Cartesian
      ndir = 3
    case ("cylindrical","cylindrical_2D","cylindrical_3D")
      ndir = ndim
      r_   = 1
      z_   = 2
      if(ndir==3) phi_ = 3
      coordinate=cylindrical
    case ("cylindrical_2.5D")
      if (ndim /= 2) call mpistop("Geometry cylindrical_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1
      z_   = 2
      phi_ = 3
      coordinate=cylindrical
    case ("polar","polar_2D","polar_3D")
      ndir = ndim
      r_   = 1
      phi_ = 2
      if(ndir==3) z_ = 3
      coordinate=cylindrical
    case ("polar_1.5D")
       if (ndim /= 1) call mpistop("Geometry polar_1.5D but ndim /= 1")
       ndir = 2
       r_   = 1
       phi_ = 2
       coordinate=cylindrical
    case ("polar_2.5D")
      if (ndim /= 2) call mpistop("Geometry polar_2.5D but ndim /= 2")
      ndir = 3
      r_   = 1
      phi_ = 2
      z_   = 3
      coordinate=cylindrical
    case ("spherical","spherical_2D","spherical_3D")
      ndir = ndim
      r_   = 1
      if(ndir==3) phi_ = 3
      z_   = -1
      coordinate=spherical
    case ("spherical_2.5D")
      if (ndim /= 2) call mpistop("Geometry spherical_2.5D requires ndim == 2")

      ndir = 3
      r_   = 1
      phi_ = 3
      z_   = -1
      coordinate=spherical
    case default
      call mpistop("Unknown geometry specified")
    end select
  end subroutine set_coordinate_system

  subroutine set_pole
    use mod_global_parameters

    select case (coordinate)
    case (spherical) 
    case (cylindrical)
      
      if (1 == phi_ .and. periodB(1)) then
        if(mod(ng1(1),2)/=0) then
          call mpistop("Number of meshes in phi-direction should be even!")
        end if

        if(abs(xprobmin1)<smalldouble) then
          if (mype==0) then
            write(unitterm,*) "Will apply pi-periodic conditions at r=0"
          end if
          poleB(1,1)=.true.
        else
          if (mype==0) then
            write(unitterm,*) "There is no cylindrical axis!"
          end if
        end if
      end if
    end select

  end subroutine set_pole

  !> Deallocate geometry-related variables
  subroutine putgridgeo(igrid)
    use mod_global_parameters
    integer, intent(in) :: igrid

    deallocate(ps(igrid)%surfaceC,ps(igrid)%surface,ps(igrid)%dvolume,&
       ps(igrid)%dx,psc(igrid)%dx,ps(igrid)%ds,psc(igrid)%dvolume)

  end subroutine putgridgeo

  !> calculate area of surfaces of cells
  subroutine get_surface_area(s,ixGmin1,ixGmax1)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_surface

    type(state) :: s
    integer, intent(in) :: ixGmin1,ixGmax1

    double precision :: x(ixGmin1:ixGmax1,ndim), xext(ixGmin1-1:ixGmax1,1),&
        drs(ixGmin1:ixGmax1), drs_ext(ixGmin1-1:ixGmax1), dx2(ixGmin1:ixGmax1),&
        dx3(ixGmin1:ixGmax1)
    double precision :: exp_factor_ext(ixGmin1-1:ixGmax1),&
       del_exp_factor_ext(ixGmin1-1:ixGmax1),&
       exp_factor_primitive_ext(ixGmin1-1:ixGmax1)
    double precision :: exp_factor(ixGmin1:ixGmax1),&
       del_exp_factor(ixGmin1:ixGmax1),exp_factor_primitive(ixGmin1:ixGmax1)

    select case (coordinate)

    case (Cartesian_expansion)
      drs(ixGmin1:ixGmax1)=s%dx(ixGmin1:ixGmax1,1)
      x(ixGmin1:ixGmax1,1)=s%x(ixGmin1:ixGmax1,1)
      
      if(associated(usr_set_surface))then
           call usr_set_surface(ixGmin1,ixGmax1,x,drs,exp_factor,&
              del_exp_factor,exp_factor_primitive)
           if (any(exp_factor <= zero)) call &
              mpistop("The area must always be strictly positive!")
      endif
      s%surface(ixGmin1:ixGmax1,1)=exp_factor(ixGmin1:ixGmax1)
      xext(ixGmin1-1,1)=x(1,1)-half*drs(1)
      xext(ixGmin1:ixGmax1,1)=x(ixGmin1:ixGmax1,1)+half*drs(ixGmin1:ixGmax1)
      drs_ext(ixGmin1-1)=drs(1)
      drs_ext(ixGmin1:ixGmax1)=drs(ixGmin1:ixGmax1)
      if(associated(usr_set_surface)) call usr_set_surface(ixGmin1-1,ixGmax1,&
         xext,drs_ext,exp_factor_ext,del_exp_factor_ext,&
         exp_factor_primitive_ext)
      s%surfaceC(ixGmin1-1:ixGmax1,1)=exp_factor_ext(ixGmin1-1:ixGmax1)
     

    case (Cartesian,Cartesian_stretched)
      drs(ixGmin1:ixGmax1)=s%dx(ixGmin1:ixGmax1,1)
      
      

      
      s%surfaceC(ixGmin1:ixGmax1,1)=1.d0
      s%surface(ixGmin1:ixGmax1,1) =1.d0
     
      
      
      s%surfaceC(0,1)=s%surfaceC(1,1);
    case (spherical)
      x(ixGmin1:ixGmax1,1)=s%x(ixGmin1:ixGmax1,1)
      
      drs(ixGmin1:ixGmax1)=s%dx(ixGmin1:ixGmax1,1)
      
      

      s%surfaceC(ixGmin1:ixGmax1,1)=(x(ixGmin1:ixGmax1,&
         1)+half*drs(ixGmin1:ixGmax1))**2 

      

      

      
      s%surfaceC(0,1)=dabs(x(1,1)-half*drs(1))**2
     
      
      

      s%surface(ixGmin1:ixGmax1,1)=x(ixGmin1:ixGmax1,1)**2 
      

      

    case (cylindrical)
      x(ixGmin1:ixGmax1,1)=s%x(ixGmin1:ixGmax1,1)
      drs(ixGmin1:ixGmax1)=s%dx(ixGmin1:ixGmax1,1)
      
      

      s%surfaceC(ixGmin1:ixGmax1,1)=dabs(x(ixGmin1:ixGmax1,&
         1)+half*drs(ixGmin1:ixGmax1))
      
      
      
      s%surfaceC(0,1)=dabs(x(1,1)-half*drs(1))
     
      
      
      

      s%surface(ixGmin1:ixGmax1,1)=dabs(x(ixGmin1:ixGmax1,1))
      
      

    case default
      call mpistop("Sorry, coordinate unknown")
    end select

  end subroutine get_surface_area

  !> Calculate gradient of a scalar q within ixL in direction idir
  subroutine gradient(q,ixImin1,ixImax1,ixOmin1,ixOmax1,idir,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idir
    double precision, intent(in)    :: q(ixImin1:ixImax1)
    double precision, intent(inout) :: gradq(ixImin1:ixImax1)
    double precision                :: x(ixImin1:ixImax1,1:ndim)
    integer                         :: jxOmin1,jxOmax1, hxOmin1,hxOmax1

    x(ixImin1:ixImax1,1:ndim)=block%x(ixImin1:ixImax1,1:ndim)

    hxOmin1=ixOmin1-kr(idir,1);hxOmax1=ixOmax1-kr(idir,1);
    jxOmin1=ixOmin1+kr(idir,1);jxOmax1=ixOmax1+kr(idir,1);
    select case(coordinate)
    case(Cartesian)
      gradq(ixOmin1:ixOmax1)=half*(q(jxOmin1:jxOmax1)-&
         q(hxOmin1:hxOmax1))/dxlevel(idir)
    case(Cartesian_stretched,Cartesian_expansion)
      gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
         q(hxOmin1:hxOmax1))/(x(jxOmin1:jxOmax1,idir)-x(hxOmin1:hxOmax1,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
           q(hxOmin1:hxOmax1))/((x(jxOmin1:jxOmax1,1)-x(hxOmin1:hxOmax1,1)))
        
        
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
           q(hxOmin1:hxOmax1))/((x(jxOmin1:jxOmax1,phi_)-x(hxOmin1:hxOmax1,&
           phi_))*x(ixOmin1:ixOmax1,r_))
      else
        gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
           q(hxOmin1:hxOmax1))/(x(jxOmin1:jxOmax1,idir)-x(hxOmin1:hxOmax1,&
           idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

  end subroutine gradient

  !> Calculate gradient of a scalar q in direction idir at cell interfaces
  subroutine gradientx(q,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idir,gradq,&
     fourth_order)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, idir
    double precision, intent(in)    :: q(ixImin1:ixImax1), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(inout) :: gradq(ixImin1:ixImax1)
    logical, intent(in)             :: fourth_order
    integer                         :: jxOmin1,jxOmax1, hxOmin1,hxOmax1,&
        kxOmin1,kxOmax1

    if(fourth_order) then
      ! Fourth order, stencil width is two
      kxOmin1=ixOmin1-2;kxOmax1=ixOmax1+2;
      if(ixImin1>kxOmin1.or.ixImax1<kxOmax1) call &
         mpistop("Error in gradientx: Non-conforming input limits")
      hxOmin1=ixOmin1-kr(idir,1);hxOmax1=ixOmax1-kr(idir,1);
      jxOmin1=ixOmin1+kr(idir,1);jxOmax1=ixOmax1+kr(idir,1);
      kxOmin1=ixOmin1+kr(idir,1)*2;kxOmax1=ixOmax1+kr(idir,1)*2;
    else
      hxOmin1=ixOmin1;hxOmax1=ixOmax1;
    end if
    jxOmin1=ixOmin1+kr(idir,1);jxOmax1=ixOmax1+kr(idir,1);
    select case(coordinate)
    case(Cartesian)
      if(fourth_order) then
        gradq(ixOmin1:ixOmax1)=(27.d0*(q(jxOmin1:jxOmax1)-q(ixOmin1:ixOmax1))-&
           q(kxOmin1:kxOmax1)+q(hxOmin1:hxOmax1))/24.d0/dxlevel(idir)
      else
        gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
           q(hxOmin1:hxOmax1))/dxlevel(idir)
      end if
    case(Cartesian_stretched,Cartesian_expansion)
      gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
         q(hxOmin1:hxOmax1))/(x(jxOmin1:jxOmax1,idir)-x(hxOmin1:hxOmax1,idir))
    case(spherical)
      select case(idir)
      case(1)
        gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
           q(hxOmin1:hxOmax1))/((x(jxOmin1:jxOmax1,1)-x(hxOmin1:hxOmax1,1)))
        
        
      end select
    case(cylindrical)
      if(idir==phi_) then
        gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
           q(hxOmin1:hxOmax1))/((x(jxOmin1:jxOmax1,phi_)-x(hxOmin1:hxOmax1,&
           phi_))*x(ixOmin1:ixOmax1,r_))
      else
        gradq(ixOmin1:ixOmax1)=(q(jxOmin1:jxOmax1)-&
           q(hxOmin1:hxOmax1))/(x(jxOmin1:jxOmax1,idir)-x(hxOmin1:hxOmax1,&
           idir))
      end if
    case default
      call mpistop('Unknown geometry')
    end select

  end subroutine gradientx

  !> Calculate gradient of a scalar q within ixL in direction idir
  !> first use limiter to go from cell center to edge
  subroutine gradientS(q,ixImin1,ixImax1,ixOmin1,ixOmax1,idir,gradq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
        idir
    double precision, intent(in)       :: q(ixImin1:ixImax1)
    double precision, intent(inout)    :: gradq(ixImin1:ixImax1)
    double precision ,dimension(ixImin1:ixImax1) :: qC,qL,qR,dqC,ldq,rdq

    double precision :: x(ixImin1:ixImax1,1:ndim)
    double precision :: invdx
    integer          :: hxOmin1,hxOmax1,ixCmin1,ixCmax1,jxCmin1,jxCmax1,&
       gxCmin1,gxCmax1,hxCmin1,hxCmax1

    x(ixImin1:ixImax1,1:ndim)=block%x(ixImin1:ixImax1,1:ndim)

    invdx=1.d0/dxlevel(idir)
    hxOmin1=ixOmin1-kr(idir,1);hxOmax1=ixOmax1-kr(idir,1);
    ixCmin1=hxOmin1;ixCmax1=ixOmax1;
    jxCmin1=ixCmin1+kr(idir,1);jxCmax1=ixCmax1+kr(idir,1);
    gxCmin1=ixCmin1-kr(idir,1);gxCmax1=jxCmax1;
    hxCmin1=gxCmin1+kr(idir,1);hxCmax1=gxCmax1+kr(idir,1);

    ! set the gradient limiter here
    qR(gxCmin1:gxCmax1) = q(hxCmin1:hxCmax1)
    qL(gxCmin1:gxCmax1) = q(gxCmin1:gxCmax1)
    if (type_gradient_limiter(block%level)/=limiter_ppm) then
      dqC(gxCmin1:gxCmax1)= qR(gxCmin1:gxCmax1)-qL(gxCmin1:gxCmax1)
      call dwlimiter2(dqC,ixImin1,ixImax1,gxCmin1,gxCmax1,idir,&
         type_gradient_limiter(block%level),ldw=ldq,rdw=rdq)
      qL(ixCmin1:ixCmax1) = qL(ixCmin1:ixCmax1) + half*ldq(ixCmin1:ixCmax1)
      qR(ixCmin1:ixCmax1) = qR(ixCmin1:ixCmax1) - half*rdq(jxCmin1:jxCmax1)
    else
      call PPMlimitervar(ixImin1,ixImax1,ixOmin1,ixOmax1,idir,q,q,qL,qR)
    endif

    select case(coordinate)
    case(Cartesian)
      gradq(ixOmin1:ixOmax1)=(qR(ixOmin1:ixOmax1)-qL(hxOmin1:hxOmax1))*invdx
    case(Cartesian_stretched,Cartesian_expansion)
      gradq(ixOmin1:ixOmax1)=(qR(ixOmin1:ixOmax1)-&
         qL(hxOmin1:hxOmax1))/block%dx(ixOmin1:ixOmax1,idir)
    case(spherical)
      gradq(ixOmin1:ixOmax1)=(qR(ixOmin1:ixOmax1)-&
         qL(hxOmin1:hxOmax1))/block%dx(ixOmin1:ixOmax1,idir)
      select case(idir)
      case(2)
        gradq(ixOmin1:ixOmax1)=gradq(ixOmin1:ixOmax1)/x(ixOmin1:ixOmax1,1)
        
      end select
    case(cylindrical)
      gradq(ixOmin1:ixOmax1)=(qR(ixOmin1:ixOmax1)-&
         qL(hxOmin1:hxOmax1))/block%dx(ixOmin1:ixOmax1,idir)
      if(idir==phi_) gradq(ixOmin1:ixOmax1)=gradq(ixOmin1:ixOmax1)/x(&
         ixOmin1:ixOmax1,1)
    end select

  end subroutine gradientS

  !> Calculate divergence of a vector qvec within ixL
  subroutine divvector(qvec,ixImin1,ixImax1,ixOmin1,ixOmax1,divq,fourthorder,&
     sixthorder)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,1:ndir)
    double precision, intent(inout) :: divq(ixImin1:ixImax1)
    logical, intent(in), optional   :: fourthorder !< Default: false
    logical, intent(in), optional   :: sixthorder !< Default: false
    logical                         :: use_4th_order
    logical                         :: use_6th_order
    double precision                :: qC(ixImin1:ixImax1), invdx(1:ndim)
    integer                         :: jxOmin1,jxOmax1, hxOmin1,hxOmax1,&
        ixCmin1,ixCmax1, jxCmin1,jxCmax1
    integer                         :: idims, ixmin1,ixmax1, gxOmin1,gxOmax1,&
        kxOmin1,kxOmax1
    integer                         :: lxOmin1,lxOmax1, fxOmin1,fxOmax1

    use_4th_order = .false.
    use_6th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder
    if (present(sixthorder))  use_6th_order = sixthorder
    if(use_4th_order .and. use_6th_order) call &
       mpistop("divvector: using 4th and 6th order at the same time")

    if(use_4th_order) then
      if (.not. slab_uniform) call mpistop(&
         "divvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ixmin1=ixOmin1-2;ixmax1=ixOmax1+2;
    else if(use_6th_order) then
      ! Sixth order, stencil width is three
      if (.not. slab_uniform) call mpistop(&
         "divvector: 6th order only supported for slab geometry")
      ixmin1=ixOmin1-3;ixmax1=ixOmax1+3;
    else
      ! Second order, stencil width is one
      ixmin1=ixOmin1-1;ixmax1=ixOmax1+1;
    end if

    if (ixImin1>ixmin1.or.ixImax1<ixmax1) call &
       mpistop("Error in divvector: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixOmin1:ixOmax1)=0.0d0

    if (slab_uniform) then
      do idims=1,ndim
        if(use_6th_order) then
          lxOmin1=ixOmin1+3*kr(idims,1);lxOmax1=ixOmax1+3*kr(idims,1);
          kxOmin1=ixOmin1+2*kr(idims,1);kxOmax1=ixOmax1+2*kr(idims,1);
          jxOmin1=ixOmin1+kr(idims,1);jxOmax1=ixOmax1+kr(idims,1);
          hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
          gxOmin1=ixOmin1-2*kr(idims,1);gxOmax1=ixOmax1-2*kr(idims,1);
          fxOmin1=ixOmin1-3*kr(idims,1);fxOmax1=ixOmax1-3*kr(idims,1);
          divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+(-qvec(fxOmin1:fxOmax1,&
             idims)+9.d0*qvec(gxOmin1:gxOmax1,&
             idims)-45.d0*qvec(hxOmin1:hxOmax1,idims)+qvec(lxOmin1:lxOmax1,&
             idims)-9.d0*qvec(kxOmin1:kxOmax1,&
             idims)+45.d0*qvec(jxOmin1:jxOmax1,idims))/(60.d0*dxlevel(idims))
        else if(use_4th_order) then
          ! Use fourth order scheme
          kxOmin1=ixOmin1+2*kr(idims,1);kxOmax1=ixOmax1+2*kr(idims,1);
          jxOmin1=ixOmin1+kr(idims,1);jxOmax1=ixOmax1+kr(idims,1);
          hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
          gxOmin1=ixOmin1-2*kr(idims,1);gxOmax1=ixOmax1-2*kr(idims,1);
          divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+(-qvec(kxOmin1:kxOmax1,&
             idims)+8.d0*qvec(jxOmin1:jxOmax1,idims)-8.d0*qvec(hxOmin1:hxOmax1,&
             idims)+qvec(gxOmin1:gxOmax1,idims))/(12.d0*dxlevel(idims))
        else
          ! Use second order scheme
          jxOmin1=ixOmin1+kr(idims,1);jxOmax1=ixOmax1+kr(idims,1);
          hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
          divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+&
             half*(qvec(jxOmin1:jxOmax1,idims) - qvec(hxOmin1:hxOmax1,&
             idims))*invdx(idims)
        end if
      end do
    else
      do idims=1,ndim
        hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
        ixCmin1=hxOmin1;ixCmax1=ixOmax1;
        jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);
        if(stretched_dim(idims) .and. stretch_uncentered) then
          ! linear interpolation at cell interface along stretched dimension
          qC(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
             idims)*(qvec(ixCmin1:ixCmax1,&
             idims)+0.5d0*block%dx(ixCmin1:ixCmax1,&
             idims)*(qvec(jxCmin1:jxCmax1,idims)-qvec(ixCmin1:ixCmax1,&
             idims))/(block%x(jxCmin1:jxCmax1,idims)-block%x(ixCmin1:ixCmax1,&
             idims)))
        else
          qC(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
             idims)*half*(qvec(ixCmin1:ixCmax1,idims)+qvec(jxCmin1:jxCmax1,&
             idims))
        end if
        divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+&
           qC(ixOmin1:ixOmax1)-qC(hxOmin1:hxOmax1)
      end do
      divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)/block%dvolume(&
         ixOmin1:ixOmax1)
    end if
  end subroutine divvector

  !> Calculate divergence of a vector qvec within ixL
  !> using limited extrapolation to cell edges
  subroutine divvectorS(qvec,ixImin1,ixImax1,ixOmin1,ixOmax1,divq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)       :: qvec(ixImin1:ixImax1,1:ndir)
    double precision, intent(inout)    :: divq(ixImin1:ixImax1)
    double precision, dimension(ixImin1:ixImax1) :: qL,qR,dqC,ldq,rdq

    double precision :: invdx(1:ndim)
    integer          :: hxOmin1,hxOmax1,ixCmin1,ixCmax1,jxCmin1,jxCmax1,idims,&
       ixmin1,ixmax1,gxCmin1,gxCmax1,hxCmin1,hxCmax1

    ixmin1=ixOmin1-2;ixmax1=ixOmax1+2;

    if (ixImin1>ixmin1.or.ixImax1<ixmax1) call &
       mpistop("Error in divvectorS: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixOmin1:ixOmax1)=zero
    do idims=1,ndim
      hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
      ixCmin1=hxOmin1;ixCmax1=ixOmax1;
      jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);
      gxCmin1=ixCmin1-kr(idims,1);gxCmax1=jxCmax1;
      hxCmin1=gxCmin1+kr(idims,1);hxCmax1=gxCmax1+kr(idims,1);
      qR(gxCmin1:gxCmax1) = qvec(hxCmin1:hxCmax1,idims)
      qL(gxCmin1:gxCmax1) = qvec(gxCmin1:gxCmax1,idims)
      if(type_gradient_limiter(block%level)/=limiter_ppm) then
        dqC(gxCmin1:gxCmax1)= qR(gxCmin1:gxCmax1)-qL(gxCmin1:gxCmax1)
        call dwlimiter2(dqC,ixImin1,ixImax1,gxCmin1,gxCmax1,idims,&
           type_gradient_limiter(block%level),ldw=ldq,rdw=rdq)
        qL(ixCmin1:ixCmax1) = qL(ixCmin1:ixCmax1) + half*ldq(ixCmin1:ixCmax1)
        qR(ixCmin1:ixCmax1) = qR(ixCmin1:ixCmax1) - half*rdq(jxCmin1:jxCmax1)
      else
        dqC(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,idims)
        call PPMlimitervar(ixImin1,ixImax1,ixOmin1,ixOmax1,idims,dqC,dqC,qL,&
           qR)
      endif

      if (slab_uniform) then
        divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+&
           half*(qR(ixOmin1:ixOmax1)-qL(hxOmin1:hxOmax1))*invdx(idims)
      else
        qR(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
           idims)*qR(ixCmin1:ixCmax1)
        qL(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
           idims)*qL(ixCmin1:ixCmax1)
        divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+&
           qR(ixOmin1:ixOmax1)-qL(hxOmin1:hxOmax1)
      end if
    end do
    if(.not.slab_uniform) divq(ixOmin1:ixOmax1)=divq(&
       ixOmin1:ixOmax1)/block%dvolume(ixOmin1:ixOmax1)

  end subroutine divvectorS

  !> Calculate curl of a vector qvec within ixL
  !> Options to
  !>        employ standard second order CD evaluations
  !>        use Gauss theorem for non-Cartesian grids
  !>        use Stokes theorem for non-Cartesian grids
  subroutine curlvector(qvec,ixImin1,ixImax1,ixOmin1,ixOmax1,curlvec,idirmin,&
     idirmin0,ndir0,fourthorder)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    integer, intent(in)             :: ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,1:ndir0)
    double precision, intent(inout) :: curlvec(ixImin1:ixImax1,idirmin0:3)
    logical, intent(in), optional   :: fourthorder !< Default: false

    integer          :: ixAmin1,ixAmax1,ixCmin1,ixCmax1,jxCmin1,jxCmax1,idir,&
       jdir,kdir,hxOmin1,hxOmax1,jxOmin1,jxOmax1,kxOmin1,kxOmax1,gxOmin1,&
       gxOmax1
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixImin1:ixImax1),tmp2(ixImin1:ixImax1),&
       xC(ixImin1:ixImax1),surface(ixImin1:ixImax1)
    logical          :: use_4th_order

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.

    use_4th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder

    if (use_4th_order) then
      if (.not. slab_uniform) call mpistop(&
         "curlvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ixAmin1=ixOmin1-2;ixAmax1=ixOmax1+2;
    else
      ! Second order, stencil width is one
      ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;
    end if

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1) call &
       mpistop("Error in curlvector: Non-conforming input limits")

    idirmin=4
    curlvec(ixOmin1:ixOmax1,idirmin0:3)=zero

    ! all non-Cartesian cases
    select case(coordinate)
      case(Cartesian) ! Cartesian grids
        invdx=1.d0/dxlevel
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            if (.not. use_4th_order) then
              ! Use second order scheme
              jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
              hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
              tmp(ixOmin1:ixOmax1)=half*(qvec(jxOmin1:jxOmax1,&
                 kdir) - qvec(hxOmin1:hxOmax1,kdir))*invdx(jdir)
            else
              ! Use fourth order scheme
              kxOmin1=ixOmin1+2*kr(jdir,1);kxOmax1=ixOmax1+2*kr(jdir,1);
              jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
              hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
              gxOmin1=ixOmin1-2*kr(jdir,1);gxOmax1=ixOmax1-2*kr(jdir,1);
              tmp(ixOmin1:ixOmax1)=(-qvec(kxOmin1:kxOmax1,&
                 kdir) + 8.0d0 * qvec(jxOmin1:jxOmax1,&
                 kdir) - 8.0d0 * qvec(hxOmin1:hxOmax1,&
                 kdir) + qvec(gxOmin1:gxOmax1,kdir))/(12.0d0 * dxlevel(jdir))
            end if
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                 idir)+tmp(ixOmin1:ixOmax1)
            else
              curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                 idir)-tmp(ixOmin1:ixOmax1)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(Cartesian_stretched) ! stretched Cartesian grids
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            select case(type_curl)
              case(central)
                tmp(ixAmin1:ixAmax1)=qvec(ixAmin1:ixAmax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
                ! second order centered differencing
                tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                   tmp(hxOmin1:hxOmax1))/(block%x(jxOmin1:jxOmax1,&
                   jdir)-block%x(hxOmin1:hxOmax1,jdir))
              case(Gaussbased)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                tmp(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
                   jdir)*(qvec(ixCmin1:ixCmax1,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                   kdir))/(block%x(jxCmin1:jxCmax1,&
                   jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                   tmp(hxOmin1:hxOmax1))/block%dvolume(ixOmin1:ixOmax1)
              case(Stokesbased)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                if(kdir<=ndim)then
                  tmp(ixCmin1:ixCmax1)=block%ds(ixOmin1:ixOmax1,&
                     kdir)*(qvec(ixCmin1:ixCmax1,&
                     kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                     kdir))/(block%x(jxCmin1:jxCmax1,&
                     jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                else
                  tmp(ixCmin1:ixCmax1)=(qvec(ixCmin1:ixCmax1,&
                     kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                     kdir))/(block%x(jxCmin1:jxCmax1,&
                     jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                endif
                if(idir<=ndim)then
                  tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))/block%surface(ixOmin1:ixOmax1,idir)
                else ! essentially for 2.5D case, idir=3 and jdir,kdir<=2
                  tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))/(block%ds(ixOmin1:ixOmax1,&
                     jdir)*block%ds(ixOmin1:ixOmax1,kdir))
                endif
              case default
                call mpistop('no such curl evaluator')
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                 idir)+tmp2(ixOmin1:ixOmax1)
            else
              curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                 idir)-tmp2(ixOmin1:ixOmax1)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(spherical) ! possibly stretched spherical grids
        select case(type_curl)
          case(central) ! ok for any dimensionality
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixAmin1:ixAmax1)=qvec(ixAmin1:ixAmax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
                select case(jdir)
                case(1)
                tmp(ixAmin1:ixAmax1)=tmp(ixAmin1:ixAmax1)*block%x(&
                   ixAmin1:ixAmax1,1)
                tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                   tmp(hxOmin1:hxOmax1))/((block%x(jxOmin1:jxOmax1,&
                   1)-block%x(hxOmin1:hxOmax1,1))*block%x(ixOmin1:ixOmax1,1))
                
                
                end select
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                endif
                if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased)
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                tmp(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
                   jdir)*(qvec(ixCmin1:ixCmax1,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                   kdir))/(block%x(jxCmin1:jxCmax1,&
                   jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                   tmp(hxOmin1:hxOmax1))/block%dvolume(ixOmin1:ixOmax1)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric terms
            if(idir==2.and.phi_>0) curlvec(ixOmin1:ixOmax1,&
               2)=curlvec(ixOmin1:ixOmax1,2)+qvec(ixOmin1:ixOmax1,&
               phi_)/block%x(ixOmin1:ixOmax1,r_)
            
            enddo;
          case(Stokesbased)
            !if(ndim<3) call mpistop("Stokesbased for 3D spherical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
             if(lvc(idir,jdir,kdir)/=0)then
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(dsin(xC(ixOmin1:ixOmax1))*tmp(ixOmin1:ixOmax1)-&
                     dsin(xC(hxOmin1:hxOmax1))*tmp(hxOmin1:hxOmax1))*block%dx(&
                     ixOmin1:ixOmax1,kdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,&
                     jdir))/block%surface(ixOmin1:ixOmax1,&
                     idir)*block%x(ixOmin1:ixOmax1,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))*block%dx(ixOmin1:ixOmax1,1)
                  !! integral along 3rd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(xC(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1)-&
                     xC(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1))*dsin(block%x(&
                     ixOmin1:ixOmax1,idir))*block%dx(ixOmin1:ixOmax1,&
                     kdir))/block%surface(ixOmin1:ixOmax1,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,jdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  if(ndim==3) then
                    surface(ixOmin1:ixOmax1)=block%surface(ixOmin1:ixOmax1,&
                       idir)
                  else
                    surface(ixOmin1:ixOmax1)=block%x(ixOmin1:ixOmax1,&
                       jdir)*block%dx(ixOmin1:ixOmax1,&
                       kdir)*block%dx(ixOmin1:ixOmax1,jdir)
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(xC(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1)-&
                     xC(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1))*block%dx(&
                     ixOmin1:ixOmax1,kdir))/surface(ixOmin1:ixOmax1)
                end if
              end select
              if(idir<idirmin)idirmin=idir
             endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
      case(cylindrical) ! possibly stretched cylindrical grids
        select case(type_curl)
          case(central)  ! works for any dimensionality, polar/cylindrical
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixAmin1:ixAmax1)=qvec(ixAmin1:ixAmax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
                if(z_==3.or.z_==-1) then
                  ! Case Polar_2D, Polar_2.5D or Polar_3D, i.e. R,phi,Z
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixAmin1:ixAmax1)=tmp(&
                     ixAmin1:ixAmax1)*block%x(ixAmin1:ixAmax1,1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                     tmp(hxOmin1:hxOmax1))/(block%x(jxOmin1:jxOmax1,&
                     1)-block%x(hxOmin1:hxOmax1,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1)=tmp2(&
                     ixOmin1:ixOmax1)/block%x(ixOmin1:ixOmax1,1) !(1/R)*d(R V_phi)/dR
                  
                  
                  end select
                end if
                if(phi_==3.or.phi_==-1) then
                  ! Case Cylindrical_2D, Cylindrical_2.5D or Cylindrical_3D, i.e. R,Z,phi
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixAmin1:ixAmax1)=tmp(&
                     ixAmin1:ixAmax1)*block%x(ixAmin1:ixAmax1,1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                     tmp(hxOmin1:hxOmax1))/(block%x(jxOmin1:jxOmax1,&
                     1)-block%x(hxOmin1:hxOmax1,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1)=tmp2(&
                     ixOmin1:ixOmax1)/block%x(ixOmin1:ixOmax1,1) !(1/R)*d(R V_phi)/dR
                  
                  
                  end select
                end if
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased) ! works for any dimensionality, polar/cylindrical
            if(ndim<2) call mpistop&
               ("Gaussbased for 2D, 2.5D or 3D polar or cylindrical only")
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                tmp(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
                   jdir)*(qvec(ixCmin1:ixCmax1,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                   kdir))/(block%x(jxCmin1:jxCmax1,&
                   jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                   tmp(hxOmin1:hxOmax1))/block%dvolume(ixOmin1:ixOmax1)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric term from d e_R/d phi= e_phi for unit vectors e_R, e_phi
            !       but minus sign appears due to R,Z,phi ordering (?)
            ! note that in cylindrical 2D (RZ), phi_ is -1
            ! note that in polar 2D     (Rphi), z_ is -1
            if((idir==phi_.or.(phi_==-1.and.idir==3)).and.z_>0) then
              ! cylindrical
              if(  z_==2) curlvec(ixOmin1:ixOmax1,&
                 idir)=curlvec(ixOmin1:ixOmax1,idir)-qvec(ixOmin1:ixOmax1,&
                 z_)/block%x(ixOmin1:ixOmax1,r_)
              ! polar
              if(phi_==2) curlvec(ixOmin1:ixOmax1,&
                 idir)=curlvec(ixOmin1:ixOmax1,idir)+qvec(ixOmin1:ixOmax1,&
                 z_)/block%x(ixOmin1:ixOmax1,r_)
            endif
            enddo;
          case(Stokesbased)
            if(ndim<3) call mpistop("Stokesbased for 3D cylindrical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
               if(idir==r_) then
                if(jdir==phi_) then
                  ! idir=r,jdir=phi,kdir=z
                  !! integral along z dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))*block%dx(ixOmin1:ixOmax1,kdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%x(ixOmin1:ixOmax1,&
                     idir)*block%dx(ixOmin1:ixOmax1,&
                     jdir))/block%surface(ixOmin1:ixOmax1,idir)
                end if
               else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))*block%dx(ixOmin1:ixOmax1,jdir)
                  !! integral along z dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,&
                     kdir))/block%surface(ixOmin1:ixOmax1,idir)
                end if
               else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,jdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(xC(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1)-&
                     xC(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1))*block%dx(&
                     ixOmin1:ixOmax1,kdir))/block%surface(ixOmin1:ixOmax1,&
                     idir)
                end if
               end if
               if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
        case default
          call mpistop('not possible to calculate curl')
    end select

  end subroutine curlvector

  !> Calculate idim transverse components of curl of a vector qvec within ixL
  !> Options to
  !>        employ standard second order CD evaluations
  !>        use Gauss theorem for non-Cartesian grids
  !>        use Stokes theorem for non-Cartesian grids
  subroutine curlvector_trans(qvec,qvecc,ixImin1,ixImax1,ixOmin1,ixOmax1,&
     curlvec,idim,idirmin,idirmin0,ndir0)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    integer, intent(in)             :: idim, ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixImin1:ixImax1,1:ndir0),&
       qvecc(ixImin1:ixImax1,1:ndir0)
    double precision, intent(inout) :: curlvec(ixImin1:ixImax1,idirmin0:3)

    integer          :: ixAmin1,ixAmax1,ixCmin1,ixCmax1,jxCmin1,jxCmax1,idir,&
       jdir,kdir,hxOmin1,hxOmax1,jxOmin1,jxOmax1
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixImin1:ixImax1),tmp2(ixImin1:ixImax1),&
       xC(ixImin1:ixImax1),surface(ixImin1:ixImax1)

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin:3)
    ! Determine exact value of idirmin while doing the loop.

    idirmin=4
    curlvec(ixOmin1:ixOmax1,idirmin0:3)=zero
    ! Second order, stencil width is one
    ixAmin1=ixOmin1-1;ixAmax1=ixOmax1+1;

    ! all non-Cartesian cases
    select case(coordinate)
      case(Cartesian) ! Cartesian grids
        invdx=1.d0/dxlevel
        do idir=idirmin0,3
          do jdir=1,ndim; do kdir=1,ndir0
            if(lvc(idir,jdir,kdir)/=0)then
              if(jdir/=idim) then
                tmp(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
              else
                ! use two cell-center values for gradient at interface of the two cells
                ! because left and right neighbor interface values is unavailable at block boundary faces
                tmp(ixImin1:ixImax1)=qvecc(ixImin1:ixImax1,kdir)
                hxOmin1=ixOmin1;hxOmax1=ixOmax1;
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
              end if
              ! second order centered differencing
              tmp(ixOmin1:ixOmax1)=half*(tmp(jxOmin1:jxOmax1)-&
                 tmp(hxOmin1:hxOmax1))*invdx(jdir)
              if(lvc(idir,jdir,kdir)==1)then
                curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                   idir)+tmp(ixOmin1:ixOmax1)
              else
                curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                   idir)-tmp(ixOmin1:ixOmax1)
              endif
              if(idir<idirmin)idirmin=idir
            endif
          enddo; enddo;
        end do
      case(Cartesian_stretched) ! stretched Cartesian grids
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            select case(type_curl)
              case(central)
                tmp(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
                ! second order centered differencing
                tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                   tmp(hxOmin1:hxOmax1))/(block%x(jxOmin1:jxOmax1,&
                   jdir)-block%x(hxOmin1:hxOmax1,jdir))
              case(Gaussbased)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                tmp(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
                   jdir)*(qvec(ixCmin1:ixCmax1,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                   kdir))/(block%x(jxCmin1:jxCmax1,&
                   jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                   tmp(hxOmin1:hxOmax1))/block%dvolume(ixOmin1:ixOmax1)
              case(Stokesbased)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                if(kdir<=ndim)then
                  tmp(ixCmin1:ixCmax1)=block%ds(ixOmin1:ixOmax1,&
                     kdir)*(qvec(ixCmin1:ixCmax1,&
                     kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                     kdir))/(block%x(jxCmin1:jxCmax1,&
                     jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                else
                  tmp(ixCmin1:ixCmax1)=(qvec(ixCmin1:ixCmax1,&
                     kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                     jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                     kdir))/(block%x(jxCmin1:jxCmax1,&
                     jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                endif
                if(idir<=ndim)then
                  tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))/block%surface(ixOmin1:ixOmax1,idir)
                else ! essentially for 2.5D case, idir=3 and jdir,kdir<=2
                  tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))/(block%ds(ixOmin1:ixOmax1,&
                     jdir)*block%ds(ixOmin1:ixOmax1,kdir))
                endif
              case default
                call mpistop('no such curl evaluator')
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                 idir)+tmp2(ixOmin1:ixOmax1)
            else
              curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                 idir)-tmp2(ixOmin1:ixOmax1)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      case(spherical) ! possibly stretched spherical grids
        select case(type_curl)
          case(central) ! ok for any dimensionality
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
                select case(jdir)
                case(1)
                tmp(ixAmin1:ixAmax1)=tmp(ixAmin1:ixAmax1)*block%x(&
                   ixAmin1:ixAmax1,1)
                tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                   tmp(hxOmin1:hxOmax1))/((block%x(jxOmin1:jxOmax1,&
                   1)-block%x(hxOmin1:hxOmax1,1))*block%x(ixOmin1:ixOmax1,1))
                
                
                end select
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                endif
                if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased)
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                tmp(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
                   jdir)*(qvec(ixCmin1:ixCmax1,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                   kdir))/(block%x(jxCmin1:jxCmax1,&
                   jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                   tmp(hxOmin1:hxOmax1))/block%dvolume(ixOmin1:ixOmax1)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric terms
            if(idir==2.and.phi_>0) curlvec(ixOmin1:ixOmax1,&
               2)=curlvec(ixOmin1:ixOmax1,2)+qvec(ixOmin1:ixOmax1,&
               phi_)/block%x(ixOmin1:ixOmax1,r_)
            
            enddo;
          case(Stokesbased)
            !if(ndim<3) call mpistop("Stokesbased for 3D spherical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
             if(lvc(idir,jdir,kdir)/=0)then
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(dsin(xC(ixOmin1:ixOmax1))*tmp(ixOmin1:ixOmax1)-&
                     dsin(xC(hxOmin1:hxOmax1))*tmp(hxOmin1:hxOmax1))*block%dx(&
                     ixOmin1:ixOmax1,kdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,&
                     jdir))/block%surface(ixOmin1:ixOmax1,&
                     idir)*block%x(ixOmin1:ixOmax1,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))*block%dx(ixOmin1:ixOmax1,1)
                  !! integral along 3rd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(xC(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1)-&
                     xC(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1))*dsin(block%x(&
                     ixOmin1:ixOmax1,idir))*block%dx(ixOmin1:ixOmax1,&
                     kdir))/block%surface(ixOmin1:ixOmax1,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,jdir)
                  !! integral along 2nd dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  if(ndim==3) then
                    surface(ixOmin1:ixOmax1)=block%surface(ixOmin1:ixOmax1,&
                       idir)
                  else
                    surface(ixOmin1:ixOmax1)=block%x(ixOmin1:ixOmax1,&
                       jdir)*block%dx(ixOmin1:ixOmax1,&
                       kdir)*block%dx(ixOmin1:ixOmax1,jdir)
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(xC(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1)-&
                     xC(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1))*block%dx(&
                     ixOmin1:ixOmax1,kdir))/surface(ixOmin1:ixOmax1)
                end if
              end select
              if(idir<idirmin)idirmin=idir
             endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
      case(cylindrical) ! possibly stretched cylindrical grids
        select case(type_curl)
          case(central)  ! works for any dimensionality, polar/cylindrical
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,kdir)
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
                if(z_==3.or.z_==-1) then
                  ! Case Polar_2D, Polar_2.5D or Polar_3D, i.e. R,phi,Z
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixAmin1:ixAmax1)=tmp(&
                     ixAmin1:ixAmax1)*block%x(ixAmin1:ixAmax1,1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                     tmp(hxOmin1:hxOmax1))/(block%x(jxOmin1:jxOmax1,&
                     1)-block%x(hxOmin1:hxOmax1,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1)=tmp2(&
                     ixOmin1:ixOmax1)/block%x(ixOmin1:ixOmax1,1) !(1/R)*d(R V_phi)/dR
                  
                  
                  end select
                end if
                if(phi_==3.or.phi_==-1) then
                  ! Case Cylindrical_2D, Cylindrical_2.5D or Cylindrical_3D, i.e. R,Z,phi
                  select case(jdir)
                  case(1)
                  if(idir==z_) tmp(ixAmin1:ixAmax1)=tmp(&
                     ixAmin1:ixAmax1)*block%x(ixAmin1:ixAmax1,1) !R V_phi
                  ! computes d(R V_phi)/dR or d V_Z/dR
                  tmp2(ixOmin1:ixOmax1)=(tmp(jxOmin1:jxOmax1)-&
                     tmp(hxOmin1:hxOmax1))/(block%x(jxOmin1:jxOmax1,&
                     1)-block%x(hxOmin1:hxOmax1,1))
                  if(idir==z_) tmp2(ixOmin1:ixOmax1)=tmp2(&
                     ixOmin1:ixOmax1)/block%x(ixOmin1:ixOmax1,1) !(1/R)*d(R V_phi)/dR
                  
                  
                  end select
                end if
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case(Gaussbased) ! works for any dimensionality, polar/cylindrical
            if(ndim<2) call mpistop&
               ("Gaussbased for 2D, 2.5D or 3D polar or cylindrical only")
            do idir=idirmin0,3;
            do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
                hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                tmp(ixCmin1:ixCmax1)=block%surfaceC(ixCmin1:ixCmax1,&
                   jdir)*(qvec(ixCmin1:ixCmax1,&
                   kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                   jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                   kdir))/(block%x(jxCmin1:jxCmax1,&
                   jdir)-block%x(ixCmin1:ixCmax1,jdir)))
                tmp2(ixOmin1:ixOmax1)=(tmp(ixOmin1:ixOmax1)-&
                   tmp(hxOmin1:hxOmax1))/block%dvolume(ixOmin1:ixOmax1)
                if(lvc(idir,jdir,kdir)==1)then
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)+tmp2(ixOmin1:ixOmax1)
                else
                  curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,&
                     idir)-tmp2(ixOmin1:ixOmax1)
                 endif
                 if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo;
            ! geometric term from d e_R/d phi= e_phi for unit vectors e_R, e_phi
            !       but minus sign appears due to R,Z,phi ordering (?)
            ! note that in cylindrical 2D (RZ), phi_ is -1
            ! note that in polar 2D     (Rphi), z_ is -1
            if((idir==phi_.or.(phi_==-1.and.idir==3)).and.z_>0) then
              ! cylindrical
              if(  z_==2) curlvec(ixOmin1:ixOmax1,&
                 idir)=curlvec(ixOmin1:ixOmax1,idir)-qvec(ixOmin1:ixOmax1,&
                 z_)/block%x(ixOmin1:ixOmax1,r_)
              ! polar
              if(phi_==2) curlvec(ixOmin1:ixOmax1,&
                 idir)=curlvec(ixOmin1:ixOmax1,idir)+qvec(ixOmin1:ixOmax1,&
                 z_)/block%x(ixOmin1:ixOmax1,r_)
            endif
            enddo;
          case(Stokesbased)
            if(ndim<3) call mpistop("Stokesbased for 3D cylindrical only")
            do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
              if(lvc(idir,jdir,kdir)/=0)then
               if(idir==r_) then
                if(jdir==phi_) then
                  ! idir=r,jdir=phi,kdir=z
                  !! integral along z dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))*block%dx(ixOmin1:ixOmax1,kdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%x(ixOmin1:ixOmax1,&
                     idir)*block%dx(ixOmin1:ixOmax1,&
                     jdir))/block%surface(ixOmin1:ixOmax1,idir)
                end if
               else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(ixOmin1:ixOmax1)-&
                     tmp(hxOmin1:hxOmax1))*block%dx(ixOmin1:ixOmax1,jdir)
                  !! integral along z dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,&
                     kdir))/block%surface(ixOmin1:ixOmax1,idir)
                end if
               else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxOmin1=ixOmin1-kr(kdir,1);hxOmax1=ixOmax1-kr(kdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(kdir,1);jxCmax1=ixCmax1+kr(kdir,1);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       kdir)*(qvec(jxCmin1:jxCmax1,jdir)-qvec(ixCmin1:ixCmax1,&
                       jdir))/(block%x(jxCmin1:jxCmax1,&
                       kdir)-block%x(ixCmin1:ixCmax1,kdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       jdir)+qvec(jxCmin1:jxCmax1,jdir))
                  end if
                  curlvec(ixOmin1:ixOmax1,&
                     idir)=(tmp(hxOmin1:hxOmax1)-&
                     tmp(ixOmin1:ixOmax1))*block%dx(ixOmin1:ixOmax1,jdir)
                  !! integral along phi dimension
                  hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
                  ixCmin1=hxOmin1;ixCmax1=ixOmax1;
                  jxCmin1=ixCmin1+kr(jdir,1);jxCmax1=ixCmax1+kr(jdir,1);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir) .and. stretch_uncentered) then
                    tmp(ixCmin1:ixCmax1)=qvec(ixCmin1:ixCmax1,&
                       kdir)+0.5d0*block%dx(ixCmin1:ixCmax1,&
                       jdir)*(qvec(jxCmin1:jxCmax1,kdir)-qvec(ixCmin1:ixCmax1,&
                       kdir))/(block%x(jxCmin1:jxCmax1,&
                       jdir)-block%x(ixCmin1:ixCmax1,jdir))
                  else
                    tmp(ixCmin1:ixCmax1)=0.5d0*(qvec(ixCmin1:ixCmax1,&
                       kdir)+qvec(jxCmin1:jxCmax1,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixCmin1:ixCmax1)=block%x(ixCmin1:ixCmax1,&
                     jdir)+0.5d0*block%dx(ixCmin1:ixCmax1,jdir)
                  curlvec(ixOmin1:ixOmax1,idir)=(curlvec(ixOmin1:ixOmax1,&
                     idir)+(xC(ixOmin1:ixOmax1)*tmp(ixOmin1:ixOmax1)-&
                     xC(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1))*block%dx(&
                     ixOmin1:ixOmax1,kdir))/block%surface(ixOmin1:ixOmax1,&
                     idir)
                end if
               end if
               if(idir<idirmin)idirmin=idir
              endif
            enddo; enddo; enddo;
          case default
            call mpistop('no such curl evaluator')
        end select
        case default
          call mpistop('not possible to calculate curl')
    end select

  end subroutine curlvector_trans

end module mod_geometry
