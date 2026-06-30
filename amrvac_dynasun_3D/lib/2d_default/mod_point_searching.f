! for Cartesian coordinate system
module mod_point_searching
  use mod_global_parameters
  use mod_physics
  use mod_particle_base
  implicit none

contains

  subroutine get_point_w_ingrid(igrid,xp,wp,variable_type)
    double precision :: xp(1:ndim),wp(1:nw)
    character(*) :: variable_type
    integer, intent(in) :: igrid

    double precision :: dxb1,dxb2,xd1,xd2,xbmin1,xbmin2,xbmax1,xbmax2,temp
    integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
    integer :: ixbl1,ixbl2,ix1,ix2,ixAmin1,ixAmin2,ixAmax1,ixAmax2,j,ingrid
    double precision :: factor(0:1,0:1)

    ixImin1=ixGlo1;ixImin2=ixGlo2;ixImax1=ixGhi1;ixImax2=ixGhi2;
    ixOmin1=ixMlo1;ixOmin2=ixMlo2;ixOmax1=ixMhi1;ixOmax2=ixMhi2;
    xbmin1=rnode(rpxmin1_,igrid)-rnode(rpdx1_,igrid)
    xbmin2=rnode(rpxmin2_,igrid)-rnode(rpdx2_,igrid);
    xbmax1=rnode(rpxmax1_,igrid)+rnode(rpdx1_,igrid)
    xbmax2=rnode(rpxmax2_,igrid)+rnode(rpdx2_,igrid);
    ingrid=0
    if (xp(2)>xbmin2 .and. xp(2)<xbmax2) ingrid=ingrid+1
    if (xp(1)>xbmin1 .and. xp(1)<xbmax1) ingrid=ingrid+1

    if (ingrid==ndim) then
      dxb1=rnode(rpdx1_,igrid)
      dxb2=rnode(rpdx2_,igrid)
      ixbl1=floor((xp(1)-ps(igrid)%x(ixOmin1,ixOmin2,1))/dxb1)+ixOmin1
      ixbl2=floor((xp(2)-ps(igrid)%x(ixOmin1,ixOmin2,2))/dxb2)+ixOmin2
      xd1=(xp(1)-ps(igrid)%x(ixbl1,ixbl2,1))/dxb1
      xd2=(xp(2)-ps(igrid)%x(ixbl1,ixbl2,2))/dxb2
      ixAmin1=ixbl1
      ixAmin2=ixbl2
      ixAmax1=ixbl1+1
      ixAmax2=ixbl2+1

      do ix1=0,1
      do ix2=0,1
        factor(ix1,ix2)=abs(1-ix1-xd1)*abs(1-ix2-xd2)
      enddo
      enddo

      select case(variable_type)
        case('primitive')
          call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,&
             ixAmin2,ixAmax1,ixAmax2,ps(igrid)%w,ps(igrid)%x)
        case('conserved')

        case default
          call mpistop("get_point_w_ingrid: Unknown variable type!")
      end select

      wp=0.d0
      do ix1=0,1
      do ix2=0,1
        do j=1,nw
          wp(j)=wp(j)+factor(ix1,ix2)*ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,j)
        enddo
      enddo
      enddo

      select case(variable_type)
        case('primitive')
        call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
           ixAmax1,ixAmax2,ps(igrid)%w,ps(igrid)%x)
        case('conserved')

        case default
          call mpistop("get_point_w_ingrid: Unknown variable type!")
      end select

      if(physics_type=='mhd') then
        wp(iw_mag(1):iw_mag(ndir))=0.d0
        do ix1=0,1
        do ix2=0,1
          do j=1,ndir
            if (b0field) then
              temp=ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,&
                 iw_mag(j))+ps(igrid)%B0(ixbl1+ix1,ixbl2+ix2,j,0)
              wp(iw_mag(j))=wp(iw_mag(j))+factor(ix1,ix2)*temp
            else
              temp=ps(igrid)%w(ixbl1+ix1,ixbl2+ix2,iw_mag(j))
              wp(iw_mag(j))=wp(iw_mag(j))+factor(ix1,ix2)*temp
            endif
          enddo
        enddo
        enddo
      endif
    else
      call mpistop("get_point_w_ingrid: The point is not in given grid!")
    endif

  end subroutine get_point_w_ingrid

  subroutine get_point_w(xp,wp,variable_type)
    ! for given point (xp), provide the plasma parameters (wp) at this point

    double precision :: xp(1:ndim),wp(1:nw)
    character(*) :: variable_type

    double precision :: x3d(3)
    integer :: indomain,ipe,igrid,j

    indomain=0
    if (xp(2)>=xprobmin2 .and. xp(2)<xprobmax2) indomain=indomain+1
    if (xp(1)>=xprobmin1 .and. xp(1)<xprobmax1) indomain=indomain+1
    if (indomain==ndim) then

      ! find pe and igrid
      x3d=0.d0
      do j=1,ndim
        x3d(j)=xp(j)
      enddo
      call find_particle_ipe(x3d,igrid,ipe)

      !do interpolation to get point values
      if (mype==ipe) then
        call get_point_w_ingrid(igrid,xp,wp,variable_type)
      endif

      call MPI_BCAST(wp,nw,MPI_DOUBLE_PRECISION,ipe,icomm,ierrmpi)
    endif

  end subroutine get_point_w

  subroutine get_cell_w(xp,wc,variable_type)
    ! for given point (xp), looking for corresponding cell and then provide
    ! the plasma parameters (wc) of this cell

    double precision :: xp(1:ndim),wc(1:nw)
    character(*) :: variable_type

    double precision :: x3d(3)
    double precision :: dxb1,dxb2,xbmin1,xbmin2,xbmax1,xbmax2
    integer :: indomain,ixOmin1,ixOmin2,ixOmax1,ixOmax2,ixb1,ixb2,ixAmin1,&
       ixAmin2,ixAmax1,ixAmax2,ixImin1,ixImin2,ixImax1,ixImax2
    integer :: ipe,igrid,j

    indomain=0
    if (xp(2)>=xprobmin2 .and. xp(2)<xprobmax2) indomain=indomain+1
    if (xp(1)>=xprobmin1 .and. xp(1)<xprobmax1) indomain=indomain+1
    if (indomain==ndim) then

      ! find pe and igrid
      x3d=0.d0
      do j=1,ndim
        x3d(j)=xp(j)
      enddo
      call find_particle_ipe(x3d,igrid,ipe)

      if (mype==ipe) then
        ! looking for the cell index
        xbmin1=rnode(rpxmin1_,igrid)
        xbmin2=rnode(rpxmin2_,igrid)
        xbmax1=rnode(rpxmax1_,igrid)
        xbmax2=rnode(rpxmax2_,igrid)
        dxb1=rnode(rpdx1_,igrid)
        dxb2=rnode(rpdx2_,igrid)
        ixOmin1=ixmlo1
        ixOmin2=ixmlo2
        ixImin1=ixglo1
        ixImin2=ixglo2
        ixImax1=ixghi1
        ixImax2=ixghi2
        ixb1=floor((xp(1)-xbmin1)/dxb1)+ixOmin1
        ixb2=floor((xp(2)-xbmin2)/dxb2)+ixOmin2
        ixAmin1=ixb1
        ixAmin2=ixb2
        ixAmax1=ixb1
        ixAmax2=ixb2

        wc=0.d0
        select case(variable_type)
          case('primitive')
            call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,&
               ixAmin2,ixAmax1,ixAmax2,ps(igrid)%w,ps(igrid)%x)
            do j=1,nw
              wc(j)=ps(igrid)%w(ixb1,ixb2,j)
            enddo
            call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,&
               ixAmin2,ixAmax1,ixAmax2,ps(igrid)%w,ps(igrid)%x)
          case('conserved')
            do j=1,nw
              wc(j)=ps(igrid)%w(ixb1,ixb2,j)
            enddo
          case default
            call mpistop("get_point_w: Unknown variable type!")
        end select

        ! for magnetic field
        if(physics_type=='mhd') then
          wc(iw_mag(1):iw_mag(ndir))=0.d0
          do j=1,ndir
            if (b0field) then
              wc(iw_mag(j))=ps(igrid)%w(ixb1,ixb2,iw_mag(j))+ps(igrid)%B0(ixb1,&
                 ixb2,j,0)
            else
              wc(iw_mag(j))=ps(igrid)%w(ixb1,ixb2,iw_mag(j))
            endif
          enddo
        endif
      endif

      call MPI_BCAST(wc,nw,MPI_DOUBLE_PRECISION,ipe,icomm,ierrmpi)
    endif

  end subroutine get_cell_w

  subroutine get_cell_index(xp,ipe,igrid,ixc1,ixc2)
    ! for given point (xp), provide the corrsponding processor (ipe), 
    ! grid number (igrid) and cell index (ixc^D)

    double precision :: xp(1:ndim)
    integer :: ipe,igrid,ixc1,ixc2

    double precision :: x3d(1:3)
    double precision :: dxb1,dxb2,xbmin1,xbmin2,xbmax1,xbmax2
    integer :: indomain,ixOmin1,ixOmin2,ixOmax1,ixOmax2,j
    integer :: datas(1:ndim+2)

    indomain=0
    if (xp(2)>=xprobmin2 .and. xp(2)<xprobmax2) indomain=indomain+1
    if (xp(1)>=xprobmin1 .and. xp(1)<xprobmax1) indomain=indomain+1

    if (indomain==ndim) then
      ! find pe and igrid
      x3d=0.d0
      do j=1,ndim
        x3d(j)=xp(j)
      enddo
      call find_particle_ipe(x3d,igrid,ipe)

      if (mype==ipe) then
        ! looking for the cell index
        xbmin1=rnode(rpxmin1_,igrid)
        xbmin2=rnode(rpxmin2_,igrid)
        xbmax1=rnode(rpxmax1_,igrid)
        xbmax2=rnode(rpxmax2_,igrid)
        dxb1=rnode(rpdx1_,igrid)
        dxb2=rnode(rpdx2_,igrid)
        ixOmin1=ixmlo1
        ixOmin2=ixmlo2
        ixc1=floor((xp(1)-xbmin1)/dxb1)+ixOmin1
        ixc2=floor((xp(2)-xbmin2)/dxb2)+ixOmin2
        datas(1)=ixc1
        datas(2)=ixc2
        datas(ndim+1)=ipe
        datas(ndim+2)=igrid
      endif

      call MPI_BCAST(datas,ndim+2,MPI_INTEGER,ipe,icomm,ierrmpi)
    endif

    ixc1=datas(1)
    ixc2=datas(2)
    ipe=datas(ndim+1)
    igrid=datas(ndim+2)

  end subroutine get_cell_index

end module mod_point_searching
