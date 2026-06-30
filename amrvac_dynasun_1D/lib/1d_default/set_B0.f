subroutine set_B0_grid(igrid)
  use mod_global_parameters
  use mod_mhd_phys, only: mhd_semirelativistic

  integer, intent(in) :: igrid

  integer :: ixCoGmin1,ixCoGmax1

  ixCoGmin1=1;
  ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells;

  call set_B0_cell(ps(igrid)%B0(ixGlo1:ixGhi1,:,0),ps(igrid)%x,ixGlo1,ixGhi1,&
     ixGlo1,ixGhi1)
  if(mhd_semirelativistic) call set_B0_cell(psc(igrid)%B0(ixCoGmin1:ixCoGmax1,&
     :,0),psc(igrid)%x,ixCoGmin1,ixCoGmax1,ixCoGmin1,ixCoGmax1)
  call set_J0_cell(igrid,ps(igrid)%J0,ixGlo1,ixGhi1,ixMlo1-1,ixMhi1+1)
  call set_B0_face(igrid,ps(igrid)%x,ixGlo1,ixGhi1,ixMlo1,ixMhi1)

end subroutine set_B0_grid

subroutine set_B0_cell(wB0,x,ixImin1,ixImax1,ixmin1,ixmax1)
  use mod_usr_methods, only: usr_set_B0
  use mod_global_parameters
  use mod_geometry
  
  integer, intent(in):: ixImin1,ixImax1,ixmin1,ixmax1
  double precision, intent(inout) :: wB0(ixImin1:ixImax1,1:ndir)
  double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)

  wB0(ixmin1:ixmax1,1:ndir)=zero

  ! approximate cell-averaged B0 as cell-centered B0
  select case (coordinate)
  case (spherical)
     
  end select
  
  if (associated(usr_set_B0)) call usr_set_B0(ixImin1,ixImax1,ixmin1,ixmax1,x,&
     wB0)

end subroutine set_B0_cell

subroutine set_J0_cell(igrid,wJ0,ixImin1,ixImax1,ixmin1,ixmax1)
  use mod_usr_methods, only: usr_set_J0
  use mod_global_parameters
  use mod_geometry

  integer, intent(in):: igrid,ixImin1,ixImax1,ixmin1,ixmax1
  double precision, intent(inout) :: wJ0(ixImin1:ixImax1,7-2*ndir:3)
  integer :: idirmin0, idirmin

  if(associated(usr_set_J0)) then
    call usr_set_J0(ixImin1,ixImax1,ixmin1,ixmax1,ps(igrid)%x,wJ0)
  else
    idirmin0 = 7-2*ndir
    call curlvector(ps(igrid)%B0(ixImin1:ixImax1,:,0),ixImin1,ixImax1,ixmin1,&
       ixmax1,wJ0,idirmin,idirmin0,ndir)
  end if

end subroutine set_J0_cell

subroutine set_B0_face(igrid,x,ixImin1,ixImax1,ixOmin1,ixOmax1)
  use mod_global_parameters

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
    call set_B0_cell(ps(igrid)%B0(ixImin1:ixImax1,:,idims),xC,ixImin1,ixImax1,&
       ixCmin1,ixCmax1)
  end do

end subroutine set_B0_face
