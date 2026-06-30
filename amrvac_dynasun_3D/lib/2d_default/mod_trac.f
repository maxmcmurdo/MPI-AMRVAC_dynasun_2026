module mod_trac
  use mod_global_parameters
  use mod_uawsom
  implicit none
  ! common
  integer :: numFL,numLP
  double precision :: dL,Tmax,trac_delta,T_bott
  double precision, allocatable :: xFi(:,:)
  ! TRACB
  integer :: numxT1,numxT2
  double precision :: dxT1,dxT2 
  double precision :: xTmin(ndim),xTmax(ndim)
  double precision, allocatable :: xT(:,:,:)
  ! TRACL mask
  integer, allocatable :: trac_grid(:),ground_grid(:)
  integer :: ngrid_trac,ngrid_ground
  logical, allocatable :: trac_pe(:)
contains

  subroutine init_trac_line(mask)
    logical, intent(in) :: mask 
    integer :: refine_factor,ix1,ix2,ix(ndim),j,iFL,numL(ndim),finegrid
    double precision :: lengthFL
    double precision :: xprobmin(ndim),xprobmax(ndim),domain_nx(ndim)
    integer :: memxFi

    trac_delta=0.25d0
    !----------------- init magnetic field -------------------!
    refine_factor=2**(refine_max_level-1)
    xprobmin(1)=xprobmin1
    xprobmin(2)=xprobmin2
    xprobmax(1)=xprobmax1
    xprobmax(2)=xprobmax2
    domain_nx(1)=domain_nx1
    domain_nx(2)=domain_nx2
    dL=(xprobmax(ndim)-xprobmin(ndim))/(domain_nx(ndim)*refine_factor)
    ! max length of a field line
    if (.not. mask) phys_trac_mask=xprobmax2
    
      lengthFL=(phys_trac_mask-xprobmin(2))*3.d0
   
    
    
    numLP=floor(lengthFL/dL)
    numL=1
    numFL=1
    do j=1,ndim-1
      ! number of field lines, every 4 finest grid, every direction
      finegrid=uawsom_trac_finegrid
      numL(j)=floor((xprobmax(j)-xprobmin(j))/dL/finegrid)
      numFL=numFL*numL(j)
    end do
    allocate(xFi(numFL,ndim))
    xFi(:,ndim)=xprobmin(ndim)+dL/50.d0
    do ix2=1,numL(2)
    do ix1=1,numL(1)
      ix(1)=ix1
      ix(2)=ix2
      iFL=0
      do j=ndim-1,1,-1
        iFL=iFL+(ix(j)-(ndim-1-j))*(numL(j))**(ndim-1-j)
      end do
      xFi(iFL,1:ndim-1)=xprobmin(1:ndim-1)+&
         finegrid*dL*ix(1:ndim-1)-finegrid*dL/2.d0
    end do
    end do
    
    if(mype .eq. 0) write(*,*)&
        'NOTE:  2D TRAC method take the y-dir == grav-dir'
   
    

    memxFi=floor(8*numFL*numLP*ndim/1e6)
    if (mype==0) write(*,*) 'Memory requirement for each processor in TRAC:'
    if (mype==0) write(*,*)  memxFi,' MB'

    allocate(trac_pe(0:npe-1),trac_grid(max_blocks),ground_grid(max_blocks))

  end subroutine init_trac_line

  subroutine init_trac_block(mask)
    logical, intent(in) :: mask 
    integer :: refine_factor,finegrid,iFL,j
    integer :: ix1,ix2,ixTmin1,ixTmin2,ixTmax1,ixTmax2
    integer :: numL(ndim),ix(ndim)
    double precision :: lengthFL
    double precision :: ration,a0
    double precision :: xprobmin(ndim),xprobmax(ndim),dxT(ndim)

    refine_factor=2**(refine_max_level-1)
    xprobmin(1)=xprobmin1
    xprobmin(2)=xprobmin2
    xprobmax(1)=xprobmax1
    xprobmax(2)=xprobmax2
    dxT1=(xprobmax1-xprobmin1)/(domain_nx1*refine_factor/block_nx1)
    dxT2=(xprobmax2-xprobmin2)/(domain_nx2*refine_factor/block_nx2)
    dxT(ndim)=dxT1
    dxT(ndim)=dxT2
    finegrid=uawsom_trac_finegrid
    
    
      dL=min(dxT1,dxT2)/finegrid
   
    ! table for interpolation
    xTmin(1)=xprobmin1
    xTmin(2)=xprobmin2
    xTmax(1)=xprobmax1
    xTmax(2)=xprobmax2
    if(mask) xTmax(ndim)=phys_trac_mask
    ! max length of a field line
    if(mask) then
      lengthFL=maxval(xTmax-xprobmin)*3.d0
    else
      lengthFL=maxval(xprobmax-xprobmin)*3.d0
    end if
    numLP=floor(lengthFL/dL)
    numxT1=ceiling((xTmax(1)-xTmin(1)-smalldouble)/dxT1)
    numxT2=ceiling((xTmax(2)-xTmin(2)-smalldouble)/dxT2)
    allocate(xT(numxT1,numxT2,ndim))
    ixTmin1=1
    ixTmin2=1
    ixTmax1=numxT1
    ixTmax2=numxT2
    do j=1,numxT1
      xT(j,ixTmin2:ixTmax2,1)=(j-0.5d0)*dxT1+xTmin(1)
    end do
    do j=1,numxT2
      xT(ixTmin1:ixTmax1,j,2)=(j-0.5d0)*dxT2+xTmin(2)
    end do
    if(mask) xTmax(ndim)=maxval(xT(:,:,ndim))+half*dxT(ndim)
    numL=1
    numFL=1
    do j=1,ndim-1
      ! number of field lines, every 4 finest grid, every direction
      numL(j)=floor((xprobmax(j)-xprobmin(j))/dL)
      numFL=numFL*numL(j)
    end do
    allocate(xFi(numFL,ndim))
    xFi(:,ndim)=xprobmin(ndim)+dL/50.d0
    do ix2=1,numL(2)
    do ix1=1,numL(1)
      ix(1)=ix1
      ix(2)=ix2
      iFL=0
      do j=ndim-1,1,-1
        iFL=iFL+(ix(j)-(ndim-1-j))*(numL(j))**(ndim-1-j)
      end do
      xFi(iFL,1:ndim-1)=xprobmin(1:ndim-1)+dL*ix(1:ndim-1)-dL/2.d0
    end do
    end do
    
    if(mype .eq. 0) write(*,*)&
        'NOTE:  2D TRAC method take the y-dir == grav-dir'
   
    
  end subroutine init_trac_block

  subroutine TRAC_simple(tco_global,trac_alfa,T_peak)
    double precision, intent(in) :: tco_global, trac_alfa,T_peak
    integer :: iigrid, igrid

    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      
      if(ps(igrid)%special_values(1)<trac_alfa*ps(igrid)%special_values(2)) &
         then
        ps(igrid)%special_values(1)=trac_alfa*ps(igrid)%special_values(2)
      end if
      if(ps(igrid)%special_values(1) .lt. T_bott) then
        ps(igrid)%special_values(1)=T_bott
      else if(ps(igrid)%special_values(1) .gt. 0.2d0*T_peak) then
        ps(igrid)%special_values(1)=0.2d0*T_peak
      end if
      ps(igrid)%wextra(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         iw_tcoff)=ps(igrid)%special_values(1)
      !> special values(2) to save old tcutoff
      ps(igrid)%special_values(2)=ps(igrid)%special_values(1)
    end do
  end subroutine TRAC_simple

  subroutine LTRAC(T_peak)
    double precision, intent(in) :: T_peak
    integer :: iigrid, igrid
    integer :: ixOmin1,ixOmin2,ixOmax1,ixOmax2,trac_tcoff

    ixOmin1=ixMlo1;ixOmin2=ixMlo2;ixOmax1=ixMhi1;ixOmax2=ixMhi2;
    trac_tcoff=iw_tcoff
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      where(ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         trac_tcoff) .lt. T_bott)
        ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,trac_tcoff)=T_bott
      else where(ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         trac_tcoff) .gt. 0.2d0*T_peak)
        ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           trac_tcoff)=0.2d0*T_peak
      end where
    end do
  end subroutine LTRAC

  subroutine TRACL(mask,T_peak)
    logical, intent(in) :: mask
    double precision, intent(in) :: T_peak
    integer :: ix1,ix2,j
    integer :: iigrid, igrid
    double precision :: xF(numFL,numLP,ndim)
    integer :: numR(numFL),ixmin1,ixmin2,ixmax1,ixmax2
    double precision :: Tlcoff(numFL)
    integer :: ipel(numFL,numLP),igridl(numFL,numLP)
    logical :: forwardl(numFL)

    Tmax=T_peak
    xF=zero
    
    do j=1,ndim
      xF(:,1,j)=xFi(:,j)
    end do

    call MPI_BARRIER(icomm,ierrmpi)
    ! record pe and grid for mask region
    call update_pegrid()
    ! get temperature for the calculation of Te gradient, 
    ! which will be stored in wextra(ixI^S,Tcoff_) temporarily
    call get_Te_grid()
    ! find out direction for tracing B field
    call get_Btracing_dir(ipel,igridl,forwardl)
    ! trace Bfield and get Tcoff for each field line
    call get_Tcoff_line(xF,numR,Tlcoff,ipel,igridl,forwardl,mask)
    ! init vairable Tcoff_ and Tweight_
    call init_trac_Tcoff()
    ! get cell Tcoff via interpolation
    call interp_Tcoff(xF,ipel,igridl,numR,Tlcoff)
    ! convert primitive data back to conserved data
    call MPI_BARRIER(icomm,ierrmpi)


  end subroutine TRACL

  subroutine TRACB(mask,T_peak)
    logical, intent(in) :: mask
    double precision, intent(in) :: T_peak
    integer :: peArr(numxT1,numxT2),gdArr(numxT1,numxT2),numR(numFL)
    double precision :: Tcoff(numxT1,numxT2),Tcmax(numxT1,numxT2),Bdir(numxT1,&
       numxT2,ndim)
    double precision :: xF(numFL,numLP,ndim),Tcoff_line(numFL)
    integer :: xpe(numFL,numLP,2**ndim)
    integer :: xgd(numFL,numLP,2**ndim)

    Tmax=T_peak
    Tcoff=zero
    Tcmax=zero
    Bdir=zero
    peArr=-1
    gdArr=-1
    call block_estable(mask,Tcoff,Tcmax,Bdir,peArr,gdArr)
    xF=zero
    numR=0
    Tcoff_line=zero
    call block_trace_mfl(mask,Tcoff,Tcoff_line,Tcmax,Bdir,peArr,gdArr,xF,numR,&
       xpe,xgd)
    call block_interp_grid(mask,xF,numR,xpe,xgd,Tcoff_line)
  end subroutine TRACB

  subroutine block_estable(mask,Tcoff,Tcmax,Bdir,peArr,gdArr)
    logical :: mask
    double precision :: Tcoff(numxT1,numxT2),Tcoff_recv(numxT1,numxT2)
    double precision :: Tcmax(numxT1,numxT2),Tcmax_recv(numxT1,numxT2)
    double precision :: Bdir(numxT1,numxT2,ndim),Bdir_recv(numxT1,numxT2,ndim)
    integer :: peArr(numxT1,numxT2),peArr_recv(numxT1,numxT2)
    integer :: gdArr(numxT1,numxT2),gdArr_recv(numxT1,numxT2)
    integer :: xcmin1,xcmin2,xcmax1,xcmax2,xdmin1,xdmin2,xdmax1,xdmax2,ix1,ix2
    integer :: iigrid,igrid,numxT,intab
    double precision :: xbmin1,xbmin2,xbmax1,xbmax2

    Tcoff_recv=zero
    Tcmax_recv=zero
    Bdir_recv=zero
    peArr_recv=-1
    gdArr_recv=-1
    !> combine table from different processors 
    xcmin1=nghostcells+1
    xcmin2=nghostcells+1
    xcmax1=block_nx1+nghostcells
    xcmax2=block_nx2+nghostcells
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ps(igrid)%wextra(:,:,Tweight_)=zero
      ps(igrid)%wextra(:,:,Tcoff_)=zero
      xbmin1=rnode(rpxmin1_,igrid)-xTmin(1)
      xbmin2=rnode(rpxmin2_,igrid)-xTmin(2)
      xbmax1=rnode(rpxmax1_,igrid)-xTmin(1)
      xbmax2=rnode(rpxmax2_,igrid)-xTmin(2)
      xdmin1=nint(xbmin1/dxT1)+1
      xdmin2=nint(xbmin2/dxT2)+1
      xdmax1=ceiling((xbmax1-smalldouble)/dxT1)
      xdmax2=ceiling((xbmax2-smalldouble)/dxT2)
      do ix1=xdmin1,xdmax1 
      do ix2=xdmin2,xdmax2 
        intab=0
        if (ix1 .le. numxT1) intab=intab+1 
        if (ix2 .le. numxT2) intab=intab+1 
        if(intab .eq. ndim) then
          !> in principle, no overlap will happen here
          Tcoff(ix1,ix2)=max(Tcoff(ix1,ix2),ps(igrid)%special_values(1))
          Tcmax(ix1,ix2)=ps(igrid)%special_values(2)
          !> abs(Bdir) <= 1, so that Bdir+2 should always be positive
          Bdir(ix1,ix2,1:ndim)=ps(igrid)%special_values(3:3+ndim-1)+2.d0
          peArr(ix1,ix2)=mype
          gdArr(ix1,ix2)=igrid
        end if
      end do
      end do
    end do 
    call MPI_BARRIER(icomm,ierrmpi)
    numxT=numxT1*numxT2
    call MPI_ALLREDUCE(peArr,peArr_recv,numxT,MPI_INTEGER,MPI_MAX,icomm,&
       ierrmpi)
    call MPI_ALLREDUCE(gdArr,gdArr_recv,numxT,MPI_INTEGER,MPI_MAX,icomm,&
       ierrmpi)
    call MPI_ALLREDUCE(Tcoff,Tcoff_recv,numxT,MPI_DOUBLE_PRECISION,MPI_MAX,&
       icomm,ierrmpi)
    call MPI_ALLREDUCE(Bdir,Bdir_recv,numxT*ndim,MPI_DOUBLE_PRECISION,MPI_MAX,&
       icomm,ierrmpi)
    if(.not. mask) then
      call MPI_ALLREDUCE(Tcmax,Tcmax_recv,numxT,MPI_DOUBLE_PRECISION,MPI_MAX,&
         icomm,ierrmpi)
    end if
    peArr=peArr_recv
    gdArr=gdArr_recv
    Tcoff=Tcoff_recv
    Bdir=Bdir_recv-2.d0
    if(.not. mask) Tcmax=Tcmax_recv
  end subroutine block_estable

  subroutine block_trace_mfl(mask,Tcoff,Tcoff_line,Tcmax,Bdir,peArr,gdArr,xF,&
     numR,xpe,xgd)
    integer :: i,j,k,k1,k2,ix_next1,ix_next2
    logical :: mask,flag,first
    double precision :: Tcoff(numxT1,numxT2),Tcoff_line(numFL)
    double precision :: Tcmax(numxT1,numxT2),Tcmax_line(numFL)
    double precision :: xF(numFL,numLP,ndim)
    integer :: ix_mod(ndim,2),numR(numFL)
    double precision :: alfa_mod(ndim,2)
    double precision :: nowpoint(ndim),nowgridc(ndim)
    double precision :: Bdir(numxT1,numxT2,ndim)
    double precision :: init_dir,now_dir1(ndim),now_dir2(ndim)
    integer :: peArr(numxT1,numxT2),xpe(numFL,numLP,2**ndim)
    integer :: gdArr(numxT1,numxT2),xgd(numFL,numLP,2**ndim)

    do i=1,numFL
      flag=.true.
      k1=ceiling((xFi(i,1)-xTmin(1)-smalldouble)/dxT1)
      k2=ceiling((xFi(i,2)-xTmin(2)-smalldouble)/dxT2)
      Tcoff_line(i)=Tcoff(k1,k2)
      if(.not. mask) Tcmax_line(i)=Tcmax(k1,k2)
      ix_next1=k1
      ix_next2=k2
      j=1
      xF(i,j,:)=xFi(i,:)
      do while(flag)
        nowpoint(:)=xF(i,j,:)
        nowgridc(:)=xT(ix_next1,ix_next2,:)
        first=.true.
        if(j .eq. 1) then 
          call RK_Bdir(nowgridc,nowpoint,ix_next1,ix_next2,now_dir1,Bdir,&
             ix_mod,first,init_dir)
        else
          call RK_Bdir(nowgridc,nowpoint,ix_next1,ix_next2,now_dir1,Bdir,&
             ix_mod,first)
        end if
        
          xgd(i,j,1)=gdArr(ix_mod(1,1),ix_mod(2,1))
          xgd(i,j,2)=gdArr(ix_mod(1,2),ix_mod(2,1))
          xgd(i,j,3)=gdArr(ix_mod(1,1),ix_mod(2,2))
          xgd(i,j,4)=gdArr(ix_mod(1,2),ix_mod(2,2))
          xpe(i,j,1)=peArr(ix_mod(1,1),ix_mod(2,1))
          xpe(i,j,2)=peArr(ix_mod(1,2),ix_mod(2,1))
          xpe(i,j,3)=peArr(ix_mod(1,1),ix_mod(2,2))
          xpe(i,j,4)=peArr(ix_mod(1,2),ix_mod(2,2))
       
        
        nowpoint(:)=nowpoint(:)+init_dir*now_dir1*dL
        if(nowpoint(1) .gt. xTmax(1) .or. nowpoint(1) .lt. xTmin(1)) then
          flag=.false.
        end if
        if(nowpoint(2) .gt. xTmax(2) .or. nowpoint(2) .lt. xTmin(2)) then
          flag=.false.
        end if
        if(mask .and. nowpoint(ndim) .gt. phys_trac_mask) then
          flag=.false.
        end if
        if(flag) then
          first=.false.
          ix_next1=ceiling((nowpoint(1)-xTmin(1)-smalldouble)/dxT1)
          ix_next2=ceiling((nowpoint(2)-xTmin(2)-smalldouble)/dxT2)
          nowgridc(:)=xT(ix_next1,ix_next2,:)
          call RK_Bdir(nowgridc,nowpoint,ix_next1,ix_next2,now_dir2,Bdir,&
             ix_mod,first)
          xF(i,j+1,:)=xF(i,j,:)+init_dir*dL*half*(now_dir1+now_dir2)
          if(xF(i,j+1,1) .gt. xTmax(1) .or. xF(i,j+1,1) .lt. xTmin(1)) then
            flag=.false.
          end if
          if(xF(i,j+1,2) .gt. xTmax(2) .or. xF(i,j+1,2) .lt. xTmin(2)) then
            flag=.false.
          end if
          if(mask .and. xF(i,j+1,ndim) .gt. phys_trac_mask) then
            flag=.false.
          end if
          if(flag) then
            ix_next1=ceiling((xF(i,j+1,1)-xTmin(1)-smalldouble)/dxT1)
            ix_next2=ceiling((xF(i,j+1,2)-xTmin(2)-smalldouble)/dxT2)
            j=j+1
            Tcoff_line(i)=max(Tcoff_line(i),Tcoff(ix_next1,ix_next2))
            if(.not.mask) Tcmax_line(i)=max(Tcmax_line(i),Tcmax(ix_next1,&
               ix_next2))
          end if
        end if
      end do
      numR(i)=j
      if(mask) then
        if(Tcoff_line(i) .gt. Tmax*0.2d0) then
          Tcoff_line(i)=Tmax*0.2d0
        end if
      else
        if(Tcoff_line(i) .gt. Tcmax_line(i)*0.2d0) then
          Tcoff_line(i)=Tcmax_line(i)*0.2d0
        end if
      end if
    end do
  end subroutine block_trace_mfl

  subroutine RK_Bdir(nowgridc,nowpoint,ix_next1,ix_next2,now_dir,Bdir,ix_mod,&
     first,init_dir)
    double precision :: nowpoint(ndim),nowgridc(ndim)
    integer :: ix_mod(ndim,2)
    double precision :: alfa_mod(ndim,2)
    integer :: ix_next1,ix_next2,k1,k2
    double precision :: now_dir(ndim)
    double precision :: Bdir(numxT1,numxT2,ndim)
    logical :: first
    double precision, optional :: init_dir

    if(nowpoint(1) .gt. xTmin(1)+half*dxT1 .and. nowpoint(1) .lt. &
       xTmax(1)-half*dxT1) then
      if(nowpoint(1) .le. nowgridc(1)) then
        ix_mod(1,1)=ix_next1-1
        ix_mod(1,2)=ix_next1
        alfa_mod(1,1)=abs(nowgridc(1)-nowpoint(1))/dxT1
        alfa_mod(1,2)=one-alfa_mod(1,1)
      else
        ix_mod(1,1)=ix_next1
        ix_mod(1,2)=ix_next1+1
        alfa_mod(1,2)=abs(nowgridc(1)-nowpoint(1))/dxT1
        alfa_mod(1,1)=one-alfa_mod(1,2)
      end if
    else
      ix_mod(1,:)=ix_next1
      alfa_mod(1,:)=half
    end if
    if(nowpoint(2) .gt. xTmin(2)+half*dxT2 .and. nowpoint(2) .lt. &
       xTmax(2)-half*dxT2) then
      if(nowpoint(2) .le. nowgridc(2)) then
        ix_mod(2,1)=ix_next2-1
        ix_mod(2,2)=ix_next2
        alfa_mod(2,1)=abs(nowgridc(2)-nowpoint(2))/dxT2
        alfa_mod(2,2)=one-alfa_mod(2,1)
      else
        ix_mod(2,1)=ix_next2
        ix_mod(2,2)=ix_next2+1
        alfa_mod(2,2)=abs(nowgridc(2)-nowpoint(2))/dxT2
        alfa_mod(2,1)=one-alfa_mod(2,2)
      end if
    else
      ix_mod(2,:)=ix_next2
      alfa_mod(2,:)=half
    end if
    now_dir=zero
    
    do k1=1,2
    do k2=1,2
      now_dir=now_dir + Bdir(ix_mod(1,k1),ix_mod(2,k2),:)*alfa_mod(1,&
         k1)*alfa_mod(2,k2)
    end do
    end do
    
    
    if(present(init_dir)) then
      init_dir=sign(one,now_dir(ndim))
    end if
  end subroutine RK_Bdir

  subroutine block_interp_grid(mask,xF,numR,xpe,xgd,Tcoff_line)
    logical :: mask
    double precision :: xF(numFL,numLP,ndim)
    integer :: numR(numFL)
    integer :: xpe(numFL,numLP,2**ndim)
    integer :: xgd(numFL,numLP,2**ndim)
    double precision :: Tcoff_line(numFL)
    double precision :: weightIndex,weight,ds
    integer :: i,j,k,igrid,iigrid,ixOmin1,ixOmin2,ixOmax1,ixOmax2,ixcmin1,&
       ixcmin2,ixcmax1,ixcmax2,ixc1,ixc2
    double precision :: dxMax1,dxMax2,dxb1,dxb2

    !> interpolate lines into grid cells
    weightIndex=2.d0
    dxMax1=2.d0*dL;dxMax2=2.d0*dL;
    ixOmin1=ixMlo1;ixOmin2=ixMlo2;ixOmax1=ixMhi1;ixOmax2=ixMhi2;
    do i=1,numFL
      do j=1,numR(i) 
        do k=1,2**ndim
          if(mype .eq. xpe(i,j,k)) then
            igrid=xgd(i,j,k)
            if(igrid .le. igrids(igridstail)) then
              dxb1=rnode(rpdx1_,igrid)
              dxb2=rnode(rpdx2_,igrid)
              ixcmin1=floor((xF(i,j,1)-dxMax1-ps(igrid)%x(ixOmin1,ixOmin2,&
                 1))/dxb1)+ixOmin1
              ixcmin2=floor((xF(i,j,2)-dxMax2-ps(igrid)%x(ixOmin1,ixOmin2,&
                 2))/dxb2)+ixOmin2
              ixcmax1=floor((xF(i,j,1)+dxMax1-ps(igrid)%x(ixOmin1,ixOmin2,&
                 1))/dxb1)+ixOmin1
              ixcmax2=floor((xF(i,j,2)+dxMax2-ps(igrid)%x(ixOmin1,ixOmin2,&
                 2))/dxb2)+ixOmin2
              if (ixcmin1<ixOmin1) ixcmin1=ixOmin1
              if (ixcmin2<ixOmin2) ixcmin2=ixOmin2
              if (ixcmax1>ixOmax1) ixcmax1=ixOmax1
              if (ixcmax2>ixOmax2) ixcmax2=ixOmax2
              do ixc1=ixcmin1,ixcmax1
              do ixc2=ixcmin2,ixcmax2
                ds=0.d0
                ds=ds+(xF(i,j,1)-ps(igrid)%x(ixc1,ixc2,1))**2
                ds=ds+(xF(i,j,2)-ps(igrid)%x(ixc1,ixc2,2))**2
                ds=sqrt(ds)
                if(ds .le. 0.099d0*dL) then
                  weight=(1/(0.099d0*dL))**weightIndex
                else
                  weight=(1/ds)**weightIndex
                endif
                ps(igrid)%wextra(ixc1,ixc2,Tweight_)=ps(igrid)%wextra(ixc1,&
                   ixc2,Tweight_)+weight
                ps(igrid)%wextra(ixc1,ixc2,Tcoff_)=ps(igrid)%wextra(ixc1,ixc2,&
                   Tcoff_)+weight*Tcoff_line(i)
              enddo
              enddo
            else
              call mpistop("we need to check here 366Line in mod_trac.t")
            end if
          end if
        end do
      end do
    end do
    ! finish interpolation
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      where (ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Tweight_)>0.d0)
        ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Tcoff_)=ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Tcoff_)/ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Tweight_)
      elsewhere
        ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Tcoff_)=0.2d0*Tmax
      endwhere
    enddo
  end subroutine block_interp_grid

   ! TODO remove, not used? 
!  subroutine get_Tmax_grid(x,w,ixI^L,ixO^L,Tmax_grid)
!    integer, intent(in) :: ixI^L,ixO^L
!    double precision, intent(in) :: x(ixI^S,1:ndim)
!    double precision, intent(out) :: w(ixI^S,1:nw)
!    double precision :: Tmax_grid
!    double precision :: pth(ixI^S),Te(ixI^S)
!
!    call uawsom_get_pthermal(w,x,ixI^L,ixO^L,pth)
!    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
!    Tmax_grid=maxval(Te(ixO^S))
!  end subroutine get_Tmax_grid  

  subroutine init_trac_Tcoff()
    integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       igrid,iigrid

    ixImin1=ixGlo1;ixImin2=ixGlo2;ixImax1=ixGhi1;ixImax2=ixGhi2;
    ixOmin1=ixMlo1;ixOmin2=ixMlo2;ixOmax1=ixMhi1;ixOmax2=ixMhi2;

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ps(igrid)%wextra(ixImin1:ixImax1,ixImin2:ixImax2,Tcoff_)=0.d0
      ps(igrid)%wextra(ixImin1:ixImax1,ixImin2:ixImax2,Tweight_)=0.d0
    end do

  end subroutine init_trac_Tcoff

  subroutine update_pegrid()

    ngrid_trac=0
    ngrid_ground=0
    trac_pe(:)=.false.
    trac_grid(:)=0
    ground_grid(:)=0

    ! traverse gridtable to information of grid inside mask region
    call traverse_gridtable()

  end subroutine update_pegrid

  subroutine traverse_gridtable()
    use mod_global_parameters

    double precision :: dxb1,dxb2,xbmin1,xbmin2,xbmax1,xbmax2
    integer :: iigrid,igrid,j
    logical, allocatable :: trac_pe_recv(:)
    double precision :: hcmax_bt

    allocate(trac_pe_recv(npe))

    hcmax_bt=xprobmin2+(xprobmax2-xprobmin2)/(domain_nx2*2**(refine_max_level-&
       1))

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      xbmin2=rnode(rpxmin2_,igrid)
      if (xbmin2<phys_trac_mask) then
        trac_pe(mype)=.true.
        ngrid_trac=ngrid_trac+1
        trac_grid(ngrid_trac)=igrid
        if (xbmin2<hcmax_bt) then
          ngrid_ground=ngrid_ground+1
          ground_grid(ngrid_ground)=igrid
        endif
      endif
    enddo

    call MPI_ALLREDUCE(trac_pe,trac_pe_recv,npe,MPI_LOGICAL,MPI_LOR,icomm,&
       ierrmpi)
    trac_pe=trac_pe_recv

    deallocate(trac_pe_recv)

  end subroutine traverse_gridtable

  subroutine get_Te_grid()
    integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       igrid,iigrid,j

    ixImin1=ixGlo1;ixImin2=ixGlo2;ixImax1=ixGhi1;ixImax2=ixGhi2;
    ixOmin1=ixMlo1;ixOmin2=ixMlo2;ixOmax1=ixMhi1;ixOmax2=ixMhi2;

    do j=1,ngrid_trac
      igrid=trac_grid(j)
      
      call uawsom_get_pthermal(ps(igrid)%w,ps(igrid)%x,ixImin1,ixImin2,ixImax1,&
         ixImax2,ixImin1,ixImin2,ixImax1,ixImax2,&
         ps(igrid)%wextra(ixImin1:ixImax1,ixImin2:ixImax2,Tcoff_))
      !TODO move  outside loop for optimziation? 
      if(has_equi_rho0) then
        ps(igrid)%wextra(ixImin1:ixImax1,ixImin2:ixImax2,&
           Tcoff_)=ps(igrid)%wextra(ixImin1:ixImax1,ixImin2:ixImax2,&
           Tcoff_)/(ps(igrid)%w(ixImin1:ixImax1,ixImin2:ixImax2,&
           rho_) + ps(igrid)%equi_vars(ixImin1:ixImax1,ixImin2:ixImax2,&
           equi_rho0_,0))
      else
        ps(igrid)%wextra(ixImin1:ixImax1,ixImin2:ixImax2,&
           Tcoff_)=ps(igrid)%wextra(ixImin1:ixImax1,ixImin2:ixImax2,&
           Tcoff_)/ps(igrid)%w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
      endif
    enddo

  end subroutine get_Te_grid

  subroutine get_Btracing_dir(ipel,igridl,forwardl)
    integer :: ipel(numFL,numLP),igridl(numFL,numLP)
    logical :: forwardl(numFL)

    integer :: igrid,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iFL,j,ix1,ix2,idir,ixb1,&
       ixb2,ixbb1,ixbb2
    double precision :: xbmin1,xbmin2,xbmax1,xbmax2,dxb1,dxb2,xd1,xd2,factor,&
       Bh
    integer :: numL(ndim),ixmin(ndim),ixmax(ndim),ix(ndim)
    logical :: forwardRC(numFL)

    ipel=-1
    igridl=0
    forwardl=.TRUE.

    ixOmin1=ixmlo1;ixOmin2=ixmlo2;
    ixOmax1=ixmhi1;ixOmax2=ixmhi2;

    do j=1,ngrid_ground
      igrid=ground_grid(j)
      dxb1=rnode(rpdx1_,igrid);dxb2=rnode(rpdx2_,igrid);
      xbmin1=rnode(rpxmin1_,igrid);xbmin2=rnode(rpxmin2_,igrid);
      xbmax1=rnode(rpxmax1_,igrid);xbmax2=rnode(rpxmax2_,igrid);
      ixmin(1)=floor((xbmin1-xprobmin1)/(uawsom_trac_finegrid*dL))+1
      ixmin(2)=floor((xbmin2-xprobmin2)/(uawsom_trac_finegrid*dL))+1;
      ixmax(1)=floor((xbmax1-xprobmin1)/(uawsom_trac_finegrid*dL))
      ixmax(2)=floor((xbmax2-xprobmin2)/(uawsom_trac_finegrid*dL));
      numL(1)=floor((xprobmax1-xprobmin1)/(uawsom_trac_finegrid*dL))
      numL(2)=floor((xprobmax2-xprobmin2)/(uawsom_trac_finegrid*dL));
      ixmin(2)=1
      ixmax(2)=1
      numL(2)=1

      do ix2=ixmin(2),ixmax(2)
      do ix1=ixmin(1),ixmax(1)
        ix(1)=ix1;ix(2)=ix2;
        iFL=0
        do idir=ndim-1,1,-1
          iFL=iFL+(ix(idir)-(ndim-1-idir))*(numL(idir))**(ndim-1-idir)
        enddo
        ipel(iFL,1)=mype
        igridl(iFL,1)=igrid

        ixb1=floor((xFi(iFL,1)-ps(igrid)%x(ixOmin1,ixOmin2,1))/dxb1)+ixOmin1
        ixb2=floor((xFi(iFL,2)-ps(igrid)%x(ixOmin1,ixOmin2,2))/dxb2)+ixOmin2;
        xd1=(xFi(iFL,1)-ps(igrid)%x(ixb1,ixb2,1))/dxb1
        xd2=(xFi(iFL,2)-ps(igrid)%x(ixb1,ixb2,2))/dxb2;
        Bh=0.d0
        do ixbb1=0,1
        do ixbb2=0,1
          factor=abs(1-ix1-xd1)*abs(1-ix2-xd2)
          if (B0field) then
            Bh=Bh+factor*(ps(igrid)%w(ixb1+ixbb1,ixb2+ixbb2,&
               mag(2))+ps(igrid)%B0(ixb1+ixbb1,ixb2+ixbb2,2,0))
          else
            Bh=Bh+factor*ps(igrid)%w(ixb1+ixbb1,ixb2+ixbb2,mag(2))
          endif
        enddo
        enddo
        if (Bh>0) then
          forwardl(iFL)=.TRUE.
        else
          forwardl(iFL)=.FALSE.
        endif
      enddo
      enddo
    enddo

    call MPI_ALLREDUCE(forwardl,forwardRC,numFL,MPI_LOGICAL,MPI_LAND,icomm,&
       ierrmpi)
    forwardl=forwardRC

  end subroutine get_Btracing_dir

  subroutine get_Tcoff_line(xFL,numR,TcoffFL,ipeFL,igridFL,forwardFL,mask)
    use mod_trace_field

    double precision :: xFL(numFL,numLP,ndim)
    integer :: numR(numFL)
    double precision :: TcoffFL(numFL),TmaxFL(numFL)
    integer :: ipeFL(numFL,numLP),igridFL(numFL,numLP)
    logical :: forwardFL(numFL)
    logical, intent(in) :: mask

    integer :: nwP,nwL,iFL,iLP
    double precision :: wPm(numFL,numLP,2),wLm(numFL,1+2)
    character(len=std_len) :: ftype,tcondi

    nwP=2
    nwL=2
    ftype='Bfield'
    tcondi='TRAC'    
    call trace_field_multi(xFL,wPm,wLm,dL,numFL,numLP,nwP,nwL,forwardFL,ftype,&
       tcondi)
    do iFL=1,numFL
      numR(iFL)=int(wLm(iFL,1))
      TcoffFL(iFL)=wLm(iFL,2)
      TmaxFL(iFL)=wLm(iFL,3)
      if(mask) then
        if(TcoffFL(iFL)>0.2d0*Tmax) TcoffFL(iFL)=0.2d0*Tmax
      else
        TmaxFL(iFL)=wLm(iFL,3)
        if(TcoffFL(iFL)>0.2d0*TmaxFL(iFL)) TcoffFL(iFL)=0.2d0*TmaxFL(iFL)
      end if

      if(TcoffFL(iFL)<T_bott) TcoffFL(iFL)=T_bott
    enddo

    do iFL=1,numFL
      if (numR(iFL)>0) then
        do iLP=1,numR(iFL)
          ipeFL(iFL,iLP)=int(wPm(iFL,iLP,1))
          igridFL(iFL,iLP)=int(wPm(iFL,iLP,2))
        enddo
      endif
    enddo


  end subroutine get_Tcoff_line

  subroutine interp_Tcoff(xF,ipel,igridl,numR,Tlcoff)
    double precision :: xF(numFL,numLP,ndim)
    integer :: numR(numFL),ipel(numFL,numLP),igridl(numFL,numLP)
    double precision :: Tlcoff(numFL)

    integer :: iFL,iLP,ixOmin1,ixOmin2,ixOmax1,ixOmax2,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixcmin1,ixcmin2,ixcmax1,ixcmax2,ixbmin1,ixbmin2,ixbmax1,ixbmax2,&
       ixc1,ixc2
    integer :: igrid,j,ipmin,ipmax,igrid_nb
    double precision :: dxb1,dxb2,dxMax1,dxMax2,xbmin1,xbmin2,xbmax1,xbmax2,&
       Tcnow
    double precision :: xFnow(ndim)
    integer :: weightIndex,idn1,idn2,ixmax2
    double precision :: ds,weight

    weightIndex=2

    ixImin1=ixglo1;ixImin2=ixglo2;
    ixImax1=ixghi1;ixImax2=ixghi2;
    ixOmin1=ixmlo1;ixOmin2=ixmlo2;
    ixOmax1=ixmhi1;ixOmax2=ixmhi2;

    do iFL=1,numFL
      iLP=1
      Tcnow=Tlcoff(iFL)
      do while (iLP<=numR(iFL))      
        ! find points in the same grid
        ipmin=iLP
        do while (ipel(iFL,ipmin)/=mype .and. ipmin<=numR(iFL))
          ipmin=ipmin+1
        enddo
        igrid=igridl(iFL,ipmin)
        ipmax=ipmin
        do while (ipel(iFL,ipmax)==mype .and. igridl(iFL,&
           ipmax+1)==igrid .and. ipmax<numR(iFL))
          ipmax=ipmax+1
        enddo

        ! parameters for the grid
        dxb1=rnode(rpdx1_,igrid);dxb2=rnode(rpdx2_,igrid);
        dxMax1=4*dxb1;dxMax2=4*dxb2;

        do iLP=ipmin,ipmax
          xFnow(:)=xF(iFL,iLP,:)
          ixbmin1=floor((xFnow(1)-dxMax1-ps(igrid)%x(ixOmin1,ixOmin2,&
             1))/dxb1)+ixOmin1
          ixbmin2=floor((xFnow(2)-dxMax2-ps(igrid)%x(ixOmin1,ixOmin2,&
             2))/dxb2)+ixOmin2;
          ixbmax1=floor((xFnow(1)+dxMax1-ps(igrid)%x(ixOmin1,ixOmin2,&
             1))/dxb1)+ixOmin1
          ixbmax2=floor((xFnow(2)+dxMax2-ps(igrid)%x(ixOmin1,ixOmin2,&
             2))/dxb2)+ixOmin2;

          ! do interpolation inside the grid
          ixcmin1=max(ixbmin1,ixOmin1)
          ixcmin2=max(ixbmin2,ixOmin2)
          ixcmax1=min(ixbmax1,ixOmax1)
          ixcmax2=min(ixbmax2,ixOmax2)
          xbmin2=rnode(rpxmin2_,igrid)
          xbmax2=rnode(rpxmax2_,igrid)
          ixmax2=floor((phys_trac_mask-xbmin2)/dxb2)+ixOmin2
          if (xbmax2>phys_trac_mask) ixcmax2=min(ixmax2,ixcmax2)
          do ixc1=ixcmin1,ixcmax1
          do ixc2=ixcmin2,ixcmax2
            ds=0.d0
            ds=ds+(xFnow(1)-ps(igrid)%x(ixc1,ixc2,1))**2
            ds=ds+(xFnow(2)-ps(igrid)%x(ixc1,ixc2,2))**2
            ds=sqrt(ds)
            if(ds<1.0d-2*dxb1) then
              weight=(1/(1.0d-2*dxb1))**weightIndex
            else
              weight=(1/ds)**weightIndex
            endif
            ps(igrid)%wextra(ixc1,ixc2,Tweight_)=ps(igrid)%wextra(ixc1,ixc2,&
               Tweight_)+weight
            ps(igrid)%wextra(ixc1,ixc2,Tcoff_)=ps(igrid)%wextra(ixc1,ixc2,&
               Tcoff_)+weight*Tcnow
          enddo
          enddo
          
          ! do interpolation in neighbor grids that have the same level and pe
          
            if (ixbmin1<ixOmin1) then
              idn1=0;idn2=0;
              idn1=-1
              if (neighbor(2,idn1,idn2,igrid)==mype .and. neighbor_type(idn1,&
                 idn2,igrid)==neighbor_sibling) then
                igrid_nb=neighbor(1,idn1,idn2,igrid)
                ixcmin1=max(ixbmin1,ixOmin1);ixcmin2=max(ixbmin2,ixOmin2);
                ixcmax1=min(ixbmax1,ixOmax1);ixcmax2=min(ixbmax2,ixOmax2);
                ixcmin1=ixOmax1+(ixbmin1-ixOmin1)
                ixcmax1=ixOmax1
                xbmin2=rnode(rpxmin2_,igrid_nb)
                xbmax2=rnode(rpxmax2_,igrid_nb)
                ixmax2=floor((phys_trac_mask-xbmin2)/dxb2)+ixOmin2
                if (xbmax2>phys_trac_mask) ixcmax2=min(ixmax2,ixcmax2)

                do ixc1=ixcmin1,ixcmax1;do ixc2=ixcmin2,ixcmax2;
                  ds=0.d0
                  ds=ds+(xFnow(1)-ps(igrid_nb)%x(ixc1,ixc1,1))**2
                  ds=ds+(xFnow(2)-ps(igrid_nb)%x(ixc2,ixc2,2))**2;
                  ds=sqrt(ds)
                  if(ds<1.0d-2*dxb1) then
                    weight=(1/(1.0d-2*dxb1))**weightIndex
                  else
                    weight=(1/ds)**weightIndex
                  endif
                  ps(igrid_nb)%wextra(ixc1,ixc2,&
                     Tweight_)=ps(igrid_nb)%wextra(ixc1,ixc2,Tweight_)+weight
                  ps(igrid_nb)%wextra(ixc1,ixc2,&
                     Tcoff_)=ps(igrid_nb)%wextra(ixc1,ixc2,&
                     Tcoff_)+weight*Tcnow
                enddo;enddo;
              endif
            endif

            if (ixbmax1>ixOmin1) then
              idn1=0;idn2=0;
              idn1=1
              if (neighbor(2,idn1,idn2,igrid)==mype .and. neighbor_type(idn1,&
                 idn2,igrid)==neighbor_sibling) then
                igrid_nb=neighbor(1,idn1,idn2,igrid)
                xbmin2=rnode(rpxmin2_,igrid_nb)
                if (xbmin2<phys_trac_mask) then
                  ixcmin1=max(ixbmin1,ixOmin1);ixcmin2=max(ixbmin2,ixOmin2);
                  ixcmax1=min(ixbmax1,ixOmax1);ixcmax2=min(ixbmax2,ixOmax2);
                  ixcmin1=ixOmin1
                  ixcmax1=ixOmin1+(ixbmax1-ixOmax1)
                  xbmax2=rnode(rpxmax2_,igrid_nb)
                  ixmax2=floor((phys_trac_mask-xbmin2)/dxb2)+ixOmin2
                  if (xbmax2>phys_trac_mask) ixcmax2=min(ixmax2,ixcmax2)

                  do ixc1=ixcmin1,ixcmax1;do ixc2=ixcmin2,ixcmax2;
                    ds=0.d0
                    ds=ds+(xFnow(1)-ps(igrid_nb)%x(ixc1,ixc1,1))**2
                    ds=ds+(xFnow(2)-ps(igrid_nb)%x(ixc2,ixc2,2))**2;
                    ds=sqrt(ds)
                    if(ds<1.0d-2*dxb1) then
                      weight=(1/(1.0d-2*dxb1))**weightIndex
                    else
                      weight=(1/ds)**weightIndex
                    endif
                    ps(igrid_nb)%wextra(ixc1,ixc2,&
                       Tweight_)=ps(igrid_nb)%wextra(ixc1,ixc2,&
                       Tweight_)+weight
                    ps(igrid_nb)%wextra(ixc1,ixc2,&
                       Tcoff_)=ps(igrid_nb)%wextra(ixc1,ixc2,&
                       Tcoff_)+weight*Tcnow
                  enddo;enddo;
                endif
              endif
            endif
          
          
            if (ixbmin2<ixOmin2) then
              idn1=0;idn2=0;
              idn2=-1
              if (neighbor(2,idn1,idn2,igrid)==mype .and. neighbor_type(idn1,&
                 idn2,igrid)==neighbor_sibling) then
                igrid_nb=neighbor(1,idn1,idn2,igrid)
                ixcmin1=max(ixbmin1,ixOmin1);ixcmin2=max(ixbmin2,ixOmin2);
                ixcmax1=min(ixbmax1,ixOmax1);ixcmax2=min(ixbmax2,ixOmax2);
                ixcmin2=ixOmax2+(ixbmin2-ixOmin2)
                ixcmax2=ixOmax2
                xbmin2=rnode(rpxmin2_,igrid_nb)
                xbmax2=rnode(rpxmax2_,igrid_nb)
                ixmax2=floor((phys_trac_mask-xbmin2)/dxb2)+ixOmin2
                if (xbmax2>phys_trac_mask) ixcmax2=min(ixmax2,ixcmax2)

                do ixc1=ixcmin1,ixcmax1;do ixc2=ixcmin2,ixcmax2;
                  ds=0.d0
                  ds=ds+(xFnow(1)-ps(igrid_nb)%x(ixc1,ixc1,1))**2
                  ds=ds+(xFnow(2)-ps(igrid_nb)%x(ixc2,ixc2,2))**2;
                  ds=sqrt(ds)
                  if(ds<1.0d-2*dxb1) then
                    weight=(1/(1.0d-2*dxb1))**weightIndex
                  else
                    weight=(1/ds)**weightIndex
                  endif
                  ps(igrid_nb)%wextra(ixc1,ixc2,&
                     Tweight_)=ps(igrid_nb)%wextra(ixc1,ixc2,Tweight_)+weight
                  ps(igrid_nb)%wextra(ixc1,ixc2,&
                     Tcoff_)=ps(igrid_nb)%wextra(ixc1,ixc2,&
                     Tcoff_)+weight*Tcnow
                enddo;enddo;
              endif
            endif

            if (ixbmax2>ixOmin2) then
              idn1=0;idn2=0;
              idn2=1
              if (neighbor(2,idn1,idn2,igrid)==mype .and. neighbor_type(idn1,&
                 idn2,igrid)==neighbor_sibling) then
                igrid_nb=neighbor(1,idn1,idn2,igrid)
                xbmin2=rnode(rpxmin2_,igrid_nb)
                if (xbmin2<phys_trac_mask) then
                  ixcmin1=max(ixbmin1,ixOmin1);ixcmin2=max(ixbmin2,ixOmin2);
                  ixcmax1=min(ixbmax1,ixOmax1);ixcmax2=min(ixbmax2,ixOmax2);
                  ixcmin2=ixOmin2
                  ixcmax2=ixOmin2+(ixbmax2-ixOmax2)
                  xbmax2=rnode(rpxmax2_,igrid_nb)
                  ixmax2=floor((phys_trac_mask-xbmin2)/dxb2)+ixOmin2
                  if (xbmax2>phys_trac_mask) ixcmax2=min(ixmax2,ixcmax2)

                  do ixc1=ixcmin1,ixcmax1;do ixc2=ixcmin2,ixcmax2;
                    ds=0.d0
                    ds=ds+(xFnow(1)-ps(igrid_nb)%x(ixc1,ixc1,1))**2
                    ds=ds+(xFnow(2)-ps(igrid_nb)%x(ixc2,ixc2,2))**2;
                    ds=sqrt(ds)
                    if(ds<1.0d-2*dxb1) then
                      weight=(1/(1.0d-2*dxb1))**weightIndex
                    else
                      weight=(1/ds)**weightIndex
                    endif
                    ps(igrid_nb)%wextra(ixc1,ixc2,&
                       Tweight_)=ps(igrid_nb)%wextra(ixc1,ixc2,&
                       Tweight_)+weight
                    ps(igrid_nb)%wextra(ixc1,ixc2,&
                       Tcoff_)=ps(igrid_nb)%wextra(ixc1,ixc2,&
                       Tcoff_)+weight*Tcnow
                  enddo;enddo;
                endif
              endif
            endif
          
        enddo

        iLP=ipmax+1
      enddo
    enddo


    do j=1,ngrid_trac
      igrid=trac_grid(j)
      where(ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Tweight_)>0.d0)
        ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Tcoff_)=ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           Tcoff_)/ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Tweight_)
      elsewhere
        ps(igrid)%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Tcoff_)=T_bott
      endwhere
    enddo

  end subroutine interp_Tcoff

end module mod_trac
