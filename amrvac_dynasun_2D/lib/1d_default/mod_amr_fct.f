module mod_amr_fct
  implicit none
  private

  type facealloc
    double precision, dimension(:), pointer :: face
  end type facealloc

  type fake_neighbors
    integer :: igrid
    integer :: ipe
  end type fake_neighbors

  type(facealloc), dimension(:,:,:), allocatable, public :: pface

  type(fake_neighbors), dimension(:,:,:), allocatable,&
      public :: fine_neighbors

  integer, dimension(:,:,:), allocatable, public :: old_neighbor

  integer :: itag, isend, irecv
  integer :: nrecv, nsend, ibuf_recv, ibuf_send, ibuf_send_next
  integer, dimension(1) :: isize
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable :: recvstatus, sendstatus
  double precision, allocatable :: recvbuffer(:), sendbuffer(:)

  public :: store_faces
  public :: comm_faces
  public :: end_comm_faces
  public :: deallocateBfaces
  public :: old_neighbors
  public :: prolong_2nd_stg
  public :: already_fine

contains
  !> This subroutine performs a 2nd order prolongation for a staggered field F,
  !> preserving the divergence of the coarse cell.
  !> This is useful for preserving DivF=0.
  !> If DivF=f(x), a different algorithm must be used.
  subroutine prolong_2nd_stg(sCo,sFi,ixComin1in,ixComax1in,ixFimin1in,&
     ixFimax1in,dxCo1,xComin1,dxFi1,xFimin1,ghost,fine_min1in,fine_max1in)
    use mod_global_parameters
    use mod_physics

    logical, intent(in)          :: ghost
    integer, intent(in)          :: ixComin1in,ixComax1in, ixFimin1in,&
       ixFimax1in
    double precision, intent(in) :: dxCo1, xComin1, dxFi1, xFimin1
    type(state), intent(in)      :: sCo
    type(state), intent(inout)   :: sFi

    logical, optional :: fine_min1in,fine_max1in
    logical           :: fine_min1,fine_max1

    double precision :: eta1, invdxCo1
    integer :: ixComin1,ixComax1,ixFimin1,ixFimax1
    integer :: idim1,idim2,idim3,ixFismin1,ixFismax1,ixGsmin1,ixGsmax1,&
       ixCosmin1,ixCosmax1,ixFisCmin1,ixFisCmax1
    integer :: ixCosVmin1(1:ndim),ixCosVmax1(1:ndim),ixFisVmin1(1:ndim),&
       ixFisVmax1(1:ndim)
    integer :: hxCosmin1,hxCosmax1,jxCosmin1,jxCosmax1,ixCosEmin1,ixCosEmax1,&
       ixFisEmin1,ixFisEmax1,hxFisCmin1,hxFisCmax1,jxFisCmin1,jxFisCmax1,&
       ipxFisCmin1,ipxFisCmax1,ixCosCmin1,ixCosCmax1,imxFisCmin1,imxFisCmax1,&
       jpxFisCmin1,jpxFisCmax1,jmxFisCmin1,jmxFisCmax1,hpxFisCmin1,hpxFisCmax1
    integer :: hxFimin1,hxFimax1,jxFimin1,jxFimax1,hijxFimin1,hijxFimax1,&
       hjixFimin1,hjixFimax1,hjjxFimin1,hjjxFimax1
    integer :: iihxFimin1,iihxFimax1,iijxFimin1,iijxFimax1,ijhxFimin1,&
       ijhxFimax1,ijjxFimin1,ijjxFimax1,ihixFimin1,ihixFimax1,ijixFimin1,&
       ijixFimax1,ihjxFimin1,ihjxFimax1
    integer :: jihxFimin1,jihxFimax1,jijxFimin1,jijxFimax1,jjhxFimin1,&
       jjhxFimax1,jjjxFimin1,jjjxFimax1,jhixFimin1,jhixFimax1,jjixFimin1,&
       jjixFimax1,jhjxFimin1,jhjxFimax1
    double precision :: bfluxCo(sCo%ixGsmin1:sCo%ixGsmax1,nws),&
       bfluxFi(sFi%ixGsmin1:sFi%ixGsmax1,nws)
    double precision :: slopes(sCo%ixGsmin1:sCo%ixGsmax1,ndim),&
       B_energy_change(ixGlo1:ixGhi1)
    

    
    call mpistop&
       ("CT prolongation not implemented in 1D. But CT is not needed.")
   

     ! END NOONED
  end subroutine prolong_2nd_stg

  !> To achive consistency and thus conservation of divergence,
  !> when refining a block we take into account the faces of the
  !> already fine neighbours, if any. This routine stores them.
  subroutine store_faces
    use mod_forest, only: refine 
    use mod_global_parameters
    integer :: igrid, iigrid, idims, iside, ineighbor, ipe_neighbor
    integer :: nx1, i1, ic1, inc1

    if (npe>1) then
      nsend_fc=0
      nrecv_fc=0
    end if

    ! Size of the block face
    nx1=ixMhi1-ixMlo1+1;

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ! Check whether it is necessary to store any block face, i.e.
      ! if any coarser neighbour is going to be refined 
     do iside=1,2
        i1=kr(1,1)*(2*iside-3);
        if (neighbor_pole(i1,igrid)/=0) cycle
        if (neighbor_type(i1,igrid)==neighbor_coarse) then
          ineighbor   =neighbor(1,i1,igrid)
          ipe_neighbor=neighbor(2,i1,igrid)
          if (refine(ineighbor,ipe_neighbor)) then
            allocate(pface(iside,1,igrid)%face(1))
            !! Store the faces
            if (iside==1) then !! left side
              pface(iside,1,igrid)%face(1)=ps(igrid)%ws(ixMlo1-1,1)
            else !! right side
              pface(iside,1,igrid)%face(1)=ps(igrid)%ws(ixMhi1,1)
            end if
            if (ipe_neighbor/=mype) nsend_fc(1)=nsend_fc(1)+1
          end if
        end if
      end do

      ! If a grid is going to be refined,
      ! remember what are its neighbours.
      if (refine(igrid,mype)) then
        ! Initialize neighbour array
        fine_neighbors(:,:,igrid)%igrid=-1
        fine_neighbors(:,:,igrid)%ipe=-1
        do idims=1,ndim
          do iside=1,2
            i1=kr(1,idims)*(2*iside-3);
            if (neighbor_pole(i1,igrid)/=0) cycle
            if (neighbor_type(i1,igrid)==neighbor_fine) then
             do ic1=1+int((1+i1)/2),2-int((1-i1)/2)
                inc1=ic1+i1
                ineighbor=neighbor_child(1,inc1,igrid)
                ipe_neighbor=neighbor_child(2,inc1,igrid)

                fine_neighbors(ic1,idims,igrid)%igrid= ineighbor
                fine_neighbors(ic1,idims,igrid)%ipe=ipe_neighbor

                if (ipe_neighbor/=mype) nrecv_fc(idims)=nrecv_fc(idims)+1
             end do
            end if
          end do
        end do
      end if

    end do

  end subroutine store_faces

  !> When refining a coarse block with fine neighbours, it is necessary
  !> prolong consistently with the already fine faces.
  !> This routine takes care of the communication of such faces.
  subroutine comm_faces
    use mod_forest, only: refine
    use mod_global_parameters

    integer                   :: iigrid,igrid,ineighbor,ipe_neighbor
    integer                   :: idims,iside,i1,ic1,inc1,nx1
    integer                   :: recvsize, sendsize

    ! Communicate the block faces to achieve consistency when refining
    ! Initialize communication structures

    nrecv=0
    nsend=0
    recvsize=0
    sendsize=0

    do idims=1,ndim
       select case (idims)
       case (1)
          nrecv=nrecv+nrecv_fc(1)
          nsend=nsend+nsend_fc(1)
          nx1=1;
          isize(1)=nx1
          recvsize=recvsize+nrecv_fc(1)*isize(1)
          sendsize=sendsize+nsend_fc(1)*isize(1)
       
       end select
    end do

    if (nrecv>0) then
    ! Allocate receive buffer
      allocate(recvbuffer(recvsize),recvstatus(MPI_STATUS_SIZE,nrecv),&
          recvrequest(nrecv))
      recvrequest=MPI_REQUEST_NULL
      ibuf_recv=1
      irecv=0

    ! Receive
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if (refine(igrid,mype)) then
         do ic1=1,2
            ! Only one of the sides will be necessary,
            ! so we do the loop only over dimensions, instead of
            ! over dimensions and sizes as in the routines
            ! old_neighbors and already_fine.
            do idims=1,ndim
              ipe_neighbor=fine_neighbors(ic1,idims,igrid)%ipe
              ineighbor   =fine_neighbors(ic1,idims,igrid)%igrid
              if (ineighbor>0.and.ipe_neighbor/=mype) then
               if (idims==1) iside=ic1
                !!! Check indices
                i1=kr(1,idims)*(2*iside-3);
                if (neighbor_pole(i1,igrid)/=0) cycle
                inc1=ic1+i1;
                irecv=irecv+1
                itag=4**1*(igrid-1)+inc1*4**(1-1)

                call MPI_IRECV(recvbuffer(ibuf_recv),isize(idims),&
                    MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                   recvrequest(irecv),ierrmpi)
                ibuf_recv=ibuf_recv+isize(idims)
              end if
            end do
         end do
        end if
      end do

    end if

    if (nsend>0) then
    ! Allocate send buffer
      allocate(sendbuffer(sendsize),sendstatus(MPI_STATUS_SIZE,nsend),&
         sendrequest(nsend))
      sendrequest=MPI_REQUEST_NULL
      isend=0
      ibuf_send=1
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ! Check whether it is necessary to store any block face, i.e.
        ! if any coarser neighbour is going to be refined 
       do iside=1,2
          i1=kr(1,1)*(2*iside-3);
 !When there is a pole, faces are always zero and this is not necessary
          if (neighbor_pole(i1,igrid)/=0) cycle
          if (neighbor_type(i1,igrid)==neighbor_coarse) then
            ineighbor   =neighbor(1,i1,igrid)
            ipe_neighbor=neighbor(2,i1,igrid)
            if (refine(ineighbor,ipe_neighbor)) then
              if (ipe_neighbor/=mype) then
                ic1=1+modulo(node(pig1_,igrid)-1,2);
                inc1=-2*i1+ic1;
                itag=4**1*(ineighbor-1)+inc1*4**(1-1)
                isend=isend+1
                ibuf_send_next=ibuf_send+isize(1)
                sendbuffer(ibuf_send:ibuf_send_next-1)=reshape(pface(iside,1,&
                   igrid)%face,(/isize(1)/))
                call MPI_ISEND(sendbuffer(ibuf_send),isize(1),&
                    MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                   sendrequest(isend),ierrmpi)
                ibuf_send=ibuf_send_next
              end if
            end if
          end if
        end do
      end do
    end if

    ! Waitalls
    if (nrecv>0) then
       call MPI_WAITALL(nrecv,recvrequest,recvstatus,ierrmpi)
       deallocate(recvstatus,recvrequest)
       ibuf_recv=1
    end if

    if (nsend>0) then
       call MPI_WAITALL(nsend,sendrequest,sendstatus,ierrmpi)
       deallocate(sendbuffer,sendstatus,sendrequest)
    end if

  end subroutine comm_faces

  subroutine end_comm_faces
    use mod_global_parameters
    ! Deallocate receive buffer
    if (nrecv>0) deallocate(recvbuffer)
  end subroutine end_comm_faces

  subroutine deallocateBfaces
    use mod_global_parameters
    integer :: igrid, iigrid, iside

    do iigrid=1,igridstail; igrid=igrids(iigrid);
     do iside=1,2
        if (associated(pface(iside,1,igrid)%face)) then
          deallocate(pface(iside,1,igrid)%face)
        end if
      end do
    end do

  end subroutine deallocateBfaces

  subroutine old_neighbors(child_igrid,child_ipe,igrid,ipe)
    use mod_global_parameters
    integer, dimension(2), intent(in) :: child_igrid, child_ipe
    integer, intent(in) :: igrid, ipe
    integer :: iside, i1, ic1

    do ic1=1,2
      old_neighbor(:,:,child_igrid(ic1))=-1
     do iside=1,2
        if (ic1==iside) then
          i1=kr(1,1)*(2*iside-3);
          old_neighbor(1,i1,child_igrid(ic1))=fine_neighbors(ic1,1,&
             igrid)%igrid
          old_neighbor(2,i1,child_igrid(ic1))=fine_neighbors(ic1,1,igrid)%ipe
        end if
      end do
    end do

  end subroutine old_neighbors

  !> This routine fills the fine faces before prolonging.
  !> It is the face equivalent of fix_conserve 
  subroutine already_fine(sFi,ichild,fine_min1,fine_max1)
    use mod_forest
    use mod_global_parameters
    type(tree_node_ptr) :: tree
    type(state) :: sFi
    integer, intent(in) :: ichild
    logical :: fine_min1,fine_max1

    integer :: ineighbor,ipe_neighbor,ibufnext
    integer :: iside,iotherside,i1,nx1

    ! Size of the block face
    nx1=ixMhi1-ixMlo1+1;

    ! Initialise everything to zero and false
    fine_min1=.false.;
    fine_max1=.false.;
    sFi%ws=zero

    
  end subroutine already_fine

end module mod_amr_fct
