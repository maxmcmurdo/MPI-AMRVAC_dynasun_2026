module mod_constrained_transport
  implicit none
  public

contains
  !> re-calculate the magnetic field from the vector potential in a completely
  !> divergency free way
  subroutine recalculateB
    use mod_global_parameters
    use mod_fix_conserve
    use mod_physics
    use mod_ghostcells_update

    integer :: igrid,iigrid

    call init_comm_fix_conserve(1,ndim,1)

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! Make zero the magnetic fluxes
       ! Fake advance, storing electric fields at edges
       block=>ps(igrid)
       call fake_advance(igrid,1,1,ps(igrid))

    end do
    !$OMP END PARALLEL DO

    ! Do correction
    call recvflux(1,ndim)
    call sendflux(1,ndim)
    call fix_conserve(ps,1,ndim,1,1)

    call fix_edges(ps,1,1)

    ! Now we fill the centers for the staggered variables
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       call phys_to_primitive(ixGlo1,ixGhi1,ixMlo1,ixMhi1,ps(igrid)%w,&
          ps(igrid)%x)
       ! update cell center magnetic field
       call phys_face_to_center(ixMlo1,ixMhi1,ps(igrid))
       call phys_to_conserved(ixGlo1,ixGhi1,ixMlo1,ixMhi1,ps(igrid)%w,&
          ps(igrid)%x)
    end do
    !$OMP END PARALLEL DO

    call getbc(global_time,0.d0,ps,iwstart,nwgc)

  end subroutine recalculateB

  !> fake advance a step to calculate magnetic field
  subroutine fake_advance(igrid,idimmin,idimmax,s)
    use mod_global_parameters
    use mod_fix_conserve

    integer       :: igrid,idimmin,idimmax
    type(state)   :: s

    double precision             :: dx1
    double precision             :: fC(ixGlo1:ixGhi1,1:nwflux,1:ndim)
    double precision             :: fE(ixGlo1:ixGhi1,7-2*ndim:3)

    dx1=rnode(rpdx1_,igrid);

    call fake_update(ixGlo1,ixGhi1,s,fC,fE,dx1)

    call store_flux(igrid,fC,idimmin,idimmax,1)
    call store_edge(igrid,ixGlo1,ixGhi1,fE,idimmin,idimmax) 

  end subroutine fake_advance

  !>fake update magnetic field from vector potential
  subroutine fake_update(ixImin1,ixImax1,s,fC,fE,dx1)
    use mod_global_parameters

    integer       :: ixImin1,ixImax1
    type(state)   :: s
    double precision             :: fC(ixImin1:ixImax1,1:nwflux,1:ndim)
    double precision             :: fE(ixImin1:ixImax1,7-2*ndim:3)
    double precision             :: dx1

    integer                            :: ixIsmin1,ixIsmax1,ixOmin1,ixOmax1,&
       idir
    double precision                   :: A(s%ixGsmin1:s%ixGsmax1,1:3)

    associate(ws=>s%ws,x=>s%x)

    A=zero
    ws=zero

    ixIsmin1=s%ixGsmin1;ixIsmax1=s%ixGsmax1;
    ixOmin1=ixImin1+nghostcells;ixOmax1=ixImax1-nghostcells;

    fC=0.d0
    call b_from_vector_potentialA(ixIsmin1,ixIsmax1, ixImin1,ixImax1, ixOmin1,&
       ixOmax1, ws, x, A)

    ! This is important only in 3D
    do idir=7-2*ndim,3
       fE(ixImin1:ixImax1,idir) =-A(ixImin1:ixImax1,idir)
    end do

    end associate
  end subroutine fake_update

  !> calculate magnetic field from vector potential A at cell edges
  subroutine b_from_vector_potentialA(ixIsmin1,ixIsmax1, ixImin1,ixImax1,&
      ixOmin1,ixOmax1, ws, x, A)
    use mod_global_parameters
    use mod_usr_methods, only: usr_init_vector_potential

    integer, intent(in)                :: ixIsmin1,ixIsmax1, ixImin1,ixImax1,&
        ixOmin1,ixOmax1
    double precision, intent(inout)    :: ws(ixIsmin1:ixIsmax1,1:nws),&
       A(ixIsmin1:ixIsmax1,1:3)
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)

    integer                            :: ixCmin1,ixCmax1, hxCmin1,hxCmax1,&
        idim1, idim2, idir
    double precision                   :: xC(ixIsmin1:ixIsmax1,1:ndim),&
       xCC(ixIsmin1:ixIsmax1,1:ndim)
    double precision                   :: circ(ixIsmin1:ixIsmax1,1:ndim)

    A=zero
    ! extend one layer of cell center locations in xCC
    xC=0.d0
    xCC=0.d0
    xCC(ixImin1:ixImax1,1:ndim)=x(ixImin1:ixImax1,1:ndim)
    
    xCC(ixIsmin1,1:ndim)=x(ixImin1,1:ndim)
    xCC(ixIsmin1,1)=x(ixImin1,1)-block%dx(ixImin1,1)
    
    

    do idir=7-2*ndim,3
      ixCmax1=ixOmax1;
      ixCmin1=ixOmin1-1+kr(idir,1);
      do idim1=1,ndim
        ! Get edge coordinates
        if (idim1/=idir) then
          xC(ixCmin1:ixCmax1,idim1)=xCC(ixCmin1:ixCmax1,&
             idim1)+half*block%dx(ixCmin1:ixCmax1,idim1)
        else
          xC(ixCmin1:ixCmax1,idim1)=xCC(ixCmin1:ixCmax1,idim1)
        end if
      end do
      ! Initialise vector potential at the edge
      call usr_init_vector_potential(ixIsmin1,ixIsmax1, ixCmin1,ixCmax1, xC,&
          A(ixIsmin1:ixIsmax1,idir), idir)
      A(ixCmin1:ixCmax1,idir)=A(ixCmin1:ixCmax1,&
         idir)*block%dsC(ixCmin1:ixCmax1,idir)
    end do

    ! Take the curl of the vector potential 
    circ=zero
    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax1=ixOmax1;
      ixCmin1=ixOmin1-kr(idim1,1);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          if(lvc(idim1,idim2,idir)==0) cycle
          ! Assemble indices
          hxCmin1=ixCmin1-kr(idim2,1);hxCmax1=ixCmax1-kr(idim2,1);
          ! Add line integrals in direction idir
          circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,idim1)+lvc(idim1,&
             idim2,idir)* (A(ixCmin1:ixCmax1,idir)-A(hxCmin1:hxCmax1,idir))
        end do
      end do
      ! Divide by the area of the face to get B
      where(block%surfaceC(ixCmin1:ixCmax1,idim1)==0)
        circ(ixCmin1:ixCmax1,idim1)=zero
      elsewhere
        circ(ixCmin1:ixCmax1,idim1)=circ(ixCmin1:ixCmax1,&
           idim1)/block%surfaceC(ixCmin1:ixCmax1,idim1)
      end where
      ws(ixCmin1:ixCmax1,idim1) = circ(ixCmin1:ixCmax1,idim1)
    end do

  end subroutine b_from_vector_potentialA

  !> Reconstruct scalar q within ixO^L to 1/2 dx in direction idir
  !> Return both left and right reconstructed values 
  subroutine reconstruct(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,q,qL,qR)
    use mod_limiter
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImax1, ixCmin1,ixCmax1,&
        idir
    double precision, intent(in)       :: q(ixImin1:ixImax1)
    double precision, intent(out)      :: qL(ixImin1:ixImax1),&
        qR(ixImin1:ixImax1)

    double precision                   :: qC(ixImin1:ixImax1)
    double precision,dimension(ixImin1:ixImax1)  :: dqC,ldq,rdq
    integer                            :: ixOmin1,ixOmax1,jxCmin1,jxCmax1,&
       gxCmin1,gxCmax1,hxCmin1,hxCmax1

    jxCmin1=ixCmin1+kr(idir,1);jxCmax1=ixCmax1+kr(idir,1);
    gxCmin1=ixCmin1-kr(idir,1);gxCmax1=jxCmax1;
    hxCmin1=gxCmin1+kr(idir,1);hxCmax1=gxCmax1+kr(idir,1);

    qR(gxCmin1:gxCmax1) = q(hxCmin1:hxCmax1)
    qL(gxCmin1:gxCmax1) = q(gxCmin1:gxCmax1)

    select case (type_limiter(block%level))
    case (limiter_ppm)
       ! the ordinary grid-index:
       ixOmin1=ixCmin1+kr(idir,1);
       ixOmax1=ixCmax1;
       call PPMlimitervar(ixImin1,ixImax1,ixOmin1,ixOmax1,idir,q,q,qL,qR)
    case (limiter_mp5)
       call MP5limitervar(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,q,qL,qR)
    case (limiter_weno3)
       call WENO3limiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),q,&
          qL,qR,1)
    case (limiter_wenoyc3)
       call WENO3limiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),q,&
          qL,qR,2)
    case (limiter_weno5)
       call WENO5limiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),q,&
          qL,qR,1)
    case (limiter_weno5nm)
       call WENO5NMlimiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),&
          q,qL,qR,1)
    case (limiter_wenoz5)
       call WENO5limiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),q,&
          qL,qR,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),&
          q,qL,qR,2)
    case (limiter_wenozp5)
       call WENO5limiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),q,&
          qL,qR,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),&
          q,qL,qR,3)
    case (limiter_weno5cu6)
       call WENO5CU6limiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,q,qL,qR)
    case (limiter_teno5ad)
       call TENO5ADlimiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),&
          q,qL,qR)
    case (limiter_weno7)
       call WENO7limiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,q,qL,qR,1)
    case (limiter_mpweno7)
       call WENO7limiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,q,qL,qR,2)
    case (limiter_venk)
       call venklimiter(ixImin1,ixImax1,ixCmin1,ixCmax1,idir,dxlevel(idir),q,&
          qL,qR)
    case default
       dqC(gxCmin1:gxCmax1)= qR(gxCmin1:gxCmax1)-qL(gxCmin1:gxCmax1)
       call dwlimiter2(dqC,ixImin1,ixImax1,gxCmin1,gxCmax1,idir,&
          type_limiter(block%level),ldq,rdq)
       qL(ixCmin1:ixCmax1) = qL(ixCmin1:ixCmax1) + half*ldq(ixCmin1:ixCmax1)
       qR(ixCmin1:ixCmax1) = qR(ixCmin1:ixCmax1) - half*rdq(jxCmin1:jxCmax1)
    end select

  end subroutine reconstruct

end module mod_constrained_transport
