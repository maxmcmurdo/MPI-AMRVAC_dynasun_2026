!> update ghost cells of all blocks including physical boundaries 
module mod_ghostcells_update

  implicit none
  ! Special buffer for pole copy
  type wbuffer
    double precision, dimension(:,:), allocatable :: w
  end type wbuffer

  ! A switch of update physical boundary or not
  logical, public :: bcphys=.true.

  integer :: ixMmin1,ixMmax1, ixCoGmin1,ixCoGmax1, ixCoMmin1,ixCoMmax1,&
      ixCoGsmin1,ixCoGsmax1

  ! The number of interleaving sending buffers for ghost cells
  integer, parameter :: npwbuf=2

  ! The first index goes from -1:2, where -1 is used when a block touches the
  ! lower boundary, 1 when a block touches an upper boundary, and 0 a situation
  ! away from boundary conditions, 2 when a block touched both lower and upper
  ! boundary

  ! index ranges to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(-1:2,-1:1) :: ixS_srl_min1,ixS_srl_max1, ixR_srl_min1,&
     ixR_srl_max1

  ! index ranges of staggered variables to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(1,-1:1) :: ixS_srl_stg_min1,ixS_srl_stg_max1,&
      ixR_srl_stg_min1,ixR_srl_stg_max1

  ! index ranges to send (S) restricted (r) ghost cells to coarser blocks 
  integer, dimension(-1:1,-1:1) :: ixS_r_min1,ixS_r_max1

  ! index ranges of staggered variables to send (S) restricted (r) ghost cells to coarser blocks 
  integer, dimension(1,-1:1) :: ixS_r_stg_min1,ixS_r_stg_max1

  ! index ranges to receive restriced ghost cells from finer blocks 
  integer, dimension(-1:1, 0:3) :: ixR_r_min1,ixR_r_max1

  ! index ranges of staggered variables to receive restriced ghost cells from finer blocks 
  integer, dimension(1,0:3)  :: ixR_r_stg_min1,ixR_r_stg_max1

  ! send prolongated (p) ghost cells to finer blocks, receive prolongated from coarser blocks
  integer, dimension(-1:1, 0:3) :: ixS_p_min1,ixS_p_max1, ixR_p_min1,&
     ixR_p_max1

  ! send prolongated (p) staggered ghost cells to finer blocks, receive prolongated from coarser blocks
  integer, dimension(1,0:3)  :: ixS_p_stg_min1,ixS_p_stg_max1, ixR_p_stg_min1,&
     ixR_p_stg_max1

  ! number of MPI receive-send pairs, srl: same refinement level; r: restrict; p: prolong
  integer :: nrecv_bc_srl, nsend_bc_srl, nrecv_bc_r, nsend_bc_r, nrecv_bc_p,&
      nsend_bc_p

  ! record index position of buffer arrays
  integer :: ibuf_send_srl, ibuf_recv_srl, ibuf_send_r, ibuf_recv_r,&
      ibuf_send_p, ibuf_recv_p

  ! count of times of send and receive
  integer :: isend_srl, irecv_srl, isend_r, irecv_r, isend_p, irecv_p

  ! count of times of send and receive for cell center ghost cells
  integer :: isend_c, irecv_c

  ! tag of MPI send and recv
  integer, private :: itag

  ! total sizes = cell-center normal flux + stagger-grid flux of send and receive
  integer, dimension(-1:1) :: sizes_srl_send_total, sizes_srl_recv_total

  ! sizes of buffer arrays for center-grid variable for siblings and restrict
  integer, dimension(:), allocatable :: recvrequest_c_sr, sendrequest_c_sr
  integer, dimension(:,:), allocatable :: recvstatus_c_sr, sendstatus_c_sr

  ! sizes of buffer arrays for center-grid variable for prolongation
  integer, dimension(:), allocatable :: recvrequest_c_p, sendrequest_c_p
  integer, dimension(:,:), allocatable :: recvstatus_c_p, sendstatus_c_p

  ! sizes of buffer arrays for stagger-grid variable
  integer, dimension(1,-1:1) :: sizes_srl_send_stg, sizes_srl_recv_stg

  integer, dimension(:), allocatable :: recvrequest_srl, sendrequest_srl
  integer, dimension(:,:), allocatable :: recvstatus_srl, sendstatus_srl
 
  ! buffer arrays for send and receive of siblings, allocate in build_connectivity
  double precision, dimension(:), allocatable :: recvbuffer_srl,&
      sendbuffer_srl

  integer, dimension(:), allocatable :: recvrequest_r, sendrequest_r
  integer, dimension(:,:), allocatable :: recvstatus_r, sendstatus_r

  ! buffer arrays for send and receive in restriction
  double precision, dimension(:), allocatable :: recvbuffer_r, sendbuffer_r

  integer, dimension(:), allocatable :: recvrequest_p, sendrequest_p
  integer, dimension(:,:), allocatable :: recvstatus_p, sendstatus_p

  ! buffer arrays for send and receive in prolongation
  double precision, dimension(:), allocatable :: recvbuffer_p, sendbuffer_p

  ! sizes to allocate buffer arrays for send and receive for restriction
  integer, dimension(-1:1)     :: sizes_r_send_total
  integer, dimension(0:3)      :: sizes_r_recv_total
  integer, dimension(1,-1:1) :: sizes_r_send_stg
  integer, dimension(1,0:3)  :: sizes_r_recv_stg

  ! sizes to allocate buffer arrays for send and receive for restriction
  integer, dimension(0:3)      :: sizes_p_send_total, sizes_p_recv_total
  integer, dimension(1,0:3)  :: sizes_p_send_stg, sizes_p_recv_stg

  ! There are two variants, _f indicates that all flux variables are filled,
  ! whereas _p means that part of the variables is filled 
  ! Furthermore _r_ stands for restrict, _p_ for prolongation.
  integer, dimension(-1:2,-1:1), target :: type_send_srl_f, type_recv_srl_f
  integer, dimension(-1:1,-1:1), target :: type_send_r_f
  integer, dimension(-1:1, 0:3), target :: type_recv_r_f, type_send_p_f,&
      type_recv_p_f
  integer, dimension(-1:2,-1:1), target :: type_send_srl_p1, type_recv_srl_p1
  integer, dimension(-1:1,-1:1), target :: type_send_r_p1
  integer, dimension(-1:1, 0:3), target :: type_recv_r_p1, type_send_p_p1,&
      type_recv_p_p1
  integer, dimension(-1:2,-1:1), target :: type_send_srl_p2, type_recv_srl_p2
  integer, dimension(-1:1,-1:1), target :: type_send_r_p2
  integer, dimension(-1:1, 0:3), target :: type_recv_r_p2, type_send_p_p2,&
      type_recv_p_p2
  integer, dimension(:,:), pointer :: type_send_srl, type_recv_srl,&
      type_send_r
  integer, dimension(:,:), pointer :: type_recv_r, type_send_p, type_recv_p

contains

  subroutine init_bc()
    use mod_global_parameters 
    use mod_physics, only: phys_req_diagonal, physics_type

    integer :: nghostcellsCo, interpolation_order
    integer :: nx1, nxCo1, ixGmin1,ixGmax1, i1, ic1, inc1, idir

    ixGmin1=ixGlo1;ixGmax1=ixGhi1;
    ixMmin1=ixGmin1+nghostcells;ixMmax1=ixGmax1-nghostcells;
    ixCoGmin1=1;
    ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells;
    ixCoGsmin1=0;
    ixCoGsmax1=ixCoGmax1;

    ixCoMmin1=ixCoGmin1+nghostcells;ixCoMmax1=ixCoGmax1-nghostcells;

    nx1=ixMmax1-ixMmin1+1;
    nxCo1=nx1/2;

    if(ghost_copy) then
       interpolation_order=1
    else
       interpolation_order=2
    end if
    nghostcellsCo=int((nghostcells+1)/2)
    
    if (nghostcellsCo+interpolation_order-1>nghostcells) then
       call mpistop("interpolation order for prolongation in getbc too high")
    end if

    ! (iib,i) index has following meanings: iib = 0 means it is not at any physical boundary
    ! iib=-1 means it is at the minimum side of a physical boundary  
    ! iib= 1 means it is at the maximum side of a physical boundary  
    ! i=-1 means subregion prepared for the neighbor at its minimum side 
    ! i= 1 means subregion prepared for the neighbor at its maximum side 
    
    ixS_srl_min1(:,-1)=ixMmin1
    ixS_srl_min1(:, 0)=ixMmin1
    ixS_srl_min1(:, 1)=ixMmax1+1-nghostcells
    ixS_srl_max1(:,-1)=ixMmin1-1+nghostcells
    ixS_srl_max1(:, 0)=ixMmax1
    ixS_srl_max1(:, 1)=ixMmax1
     
    ixR_srl_min1(:,-1)=1
    ixR_srl_min1(:, 0)=ixMmin1
    ixR_srl_min1(:, 1)=ixMmax1+1
    ixR_srl_max1(:,-1)=nghostcells
    ixR_srl_max1(:, 0)=ixMmax1
    ixR_srl_max1(:, 1)=ixGmax1
    
    ixS_r_min1(:,-1)=ixCoMmin1
    ixS_r_min1(:, 0)=ixCoMmin1
    ixS_r_min1(:, 1)=ixCoMmax1+1-nghostcells
    ixS_r_max1(:,-1)=ixCoMmin1-1+nghostcells
    ixS_r_max1(:, 0)=ixCoMmax1
    ixS_r_max1(:, 1)=ixCoMmax1
    
    ixR_r_min1(:, 0)=1
    ixR_r_min1(:, 1)=ixMmin1
    ixR_r_min1(:, 2)=ixMmin1+nxCo1
    ixR_r_min1(:, 3)=ixMmax1+1
    ixR_r_max1(:, 0)=nghostcells
    ixR_r_max1(:, 1)=ixMmin1-1+nxCo1
    ixR_r_max1(:, 2)=ixMmax1
    ixR_r_max1(:, 3)=ixGmax1

    ixS_p_min1(:, 0)=ixMmin1-(interpolation_order-1)
    ixS_p_min1(:, 1)=ixMmin1-(interpolation_order-1)
    ixS_p_min1(:, 2)=ixMmin1+nxCo1-nghostcellsCo-(interpolation_order-1)
    ixS_p_min1(:, 3)=ixMmax1+1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max1(:, 0)=ixMmin1-1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max1(:, 1)=ixMmin1-1+nxCo1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max1(:, 2)=ixMmax1+(interpolation_order-1)
    ixS_p_max1(:, 3)=ixMmax1+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
      ! exclude ghost-cell region when diagonal cells are unknown
      ixS_p_min1(:, 0)=ixMmin1
      ixS_p_max1(:, 3)=ixMmax1
      ixS_p_max1(:, 1)=ixMmin1-1+nxCo1+(interpolation_order-1)
      ixS_p_min1(:, 2)=ixMmin1+nxCo1-(interpolation_order-1)
    end if

    ixR_p_min1(:, 0)=ixCoMmin1-nghostcellsCo-(interpolation_order-1)
    ixR_p_min1(:, 1)=ixCoMmin1-(interpolation_order-1)
    ixR_p_min1(:, 2)=ixCoMmin1-nghostcellsCo-(interpolation_order-1)
    ixR_p_min1(:, 3)=ixCoMmax1+1-(interpolation_order-1)
    ixR_p_max1(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max1(:, 1)=ixCoMmax1+nghostcellsCo+(interpolation_order-1)
    ixR_p_max1(:, 2)=ixCoMmax1+(interpolation_order-1)
    ixR_p_max1(:, 3)=ixCoMmax1+nghostcellsCo+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
      ixR_p_max1(:, 0)=nghostcells
      ixR_p_min1(:, 3)=ixCoMmax1+1
      ixR_p_max1(:, 1)=ixCoMmax1+(interpolation_order-1)
      ixR_p_min1(:, 2)=ixCoMmin1-(interpolation_order-1)
    end if

    

    if (stagger_grid) then
      allocate(pole_buf%ws(ixGslo1:ixGshi1,nws))
      ! Staggered (face-allocated) variables
      do idir=1,ndim
       ixS_srl_stg_min1(idir,-1)=ixMmin1-kr(idir,1)
        ixS_srl_stg_max1(idir,-1)=ixMmin1-1+nghostcells
        ixS_srl_stg_min1(idir,0) =ixMmin1-kr(idir,1)
        ixS_srl_stg_max1(idir,0) =ixMmax1
        ixS_srl_stg_min1(idir,1) =ixMmax1-nghostcells+1-kr(idir,1)
        ixS_srl_stg_max1(idir,1) =ixMmax1
        
        ixR_srl_stg_min1(idir,-1)=1-kr(idir,1)
        ixR_srl_stg_max1(idir,-1)=nghostcells
        ixR_srl_stg_min1(idir,0) =ixMmin1-kr(idir,1)
        ixR_srl_stg_max1(idir,0) =ixMmax1
        ixR_srl_stg_min1(idir,1) =ixMmax1+1-kr(idir,1)
        ixR_srl_stg_max1(idir,1) =ixGmax1

        ixS_r_stg_min1(idir,-1)=ixCoMmin1-kr(idir,1)
        ixS_r_stg_max1(idir,-1)=ixCoMmin1-1+nghostcells
        ixS_r_stg_min1(idir,0) =ixCoMmin1-kr(idir,1)
        ixS_r_stg_max1(idir,0) =ixCoMmax1
        ixS_r_stg_min1(idir,1) =ixCoMmax1+1-nghostcells-kr(idir,1)
        ixS_r_stg_max1(idir,1) =ixCoMmax1
 
        ixR_r_stg_min1(idir,0)=1-kr(idir,1)
        ixR_r_stg_max1(idir,0)=nghostcells
        ixR_r_stg_min1(idir,1)=ixMmin1-kr(idir,1)
        ixR_r_stg_max1(idir,1)=ixMmin1-1+nxCo1
        ixR_r_stg_min1(idir,2)=ixMmin1+nxCo1-kr(idir,1)
        ixR_r_stg_max1(idir,2)=ixMmax1
        ixR_r_stg_min1(idir,3)=ixMmax1+1-kr(idir,1)
        ixR_r_stg_max1(idir,3)=ixGmax1
        
       if (idir==1) then
          ! Parallel components
          
          ixS_p_stg_min1(idir,0)=ixMmin1-1 ! -1 to make redundant 
          ixS_p_stg_max1(idir,0)=ixMmin1-1+nghostcellsCo
          ixS_p_stg_min1(idir,1)=ixMmin1-1 ! -1 to make redundant 
          ixS_p_stg_max1(idir,1)=ixMmin1-1+nxCo1+nghostcellsCo
          ixS_p_stg_min1(idir,2)=ixMmax1-nxCo1-nghostcellsCo
          ixS_p_stg_max1(idir,2)=ixMmax1
          ixS_p_stg_min1(idir,3)=ixMmax1-nghostcellsCo
          ixS_p_stg_max1(idir,3)=ixMmax1

          ixR_p_stg_min1(idir,0)=ixCoMmin1-1-nghostcellsCo
          ixR_p_stg_max1(idir,0)=ixCoMmin1-1
          ixR_p_stg_min1(idir,1)=ixCoMmin1-1 ! -1 to make redundant 
          ixR_p_stg_max1(idir,1)=ixCoMmax1+nghostcellsCo
          ixR_p_stg_min1(idir,2)=ixCoMmin1-1-nghostcellsCo
          ixR_p_stg_max1(idir,2)=ixCoMmax1
          ixR_p_stg_min1(idir,3)=ixCoMmax1+1-1 ! -1 to make redundant 
          ixR_p_stg_max1(idir,3)=ixCoMmax1+nghostcellsCo
          
        else
          
          ! Perpendicular component
          ixS_p_stg_min1(idir,0)=ixMmin1
          ixS_p_stg_max1(idir,0)=ixMmin1-1+nghostcellsCo+&
             (interpolation_order-1)
          ixS_p_stg_min1(idir,1)=ixMmin1
          ixS_p_stg_max1(idir,1)=ixMmin1-1+nxCo1+nghostcellsCo+&
             (interpolation_order-1)
          ixS_p_stg_min1(idir,2)=ixMmax1+&
             1-nxCo1-nghostcellsCo-(interpolation_order-1)
          ixS_p_stg_max1(idir,2)=ixMmax1
          ixS_p_stg_min1(idir,3)=ixMmax1+&
             1-nghostcellsCo-(interpolation_order-1)
          ixS_p_stg_max1(idir,3)=ixMmax1
 
          ixR_p_stg_min1(idir,0)=ixCoMmin1-nghostcellsCo-(interpolation_order-&
             1)
          ixR_p_stg_max1(idir,0)=ixCoMmin1-1
          ixR_p_stg_min1(idir,1)=ixCoMmin1
          ixR_p_stg_max1(idir,1)=ixCoMmax1+nghostcellsCo+&
             (interpolation_order-1)
          ixR_p_stg_min1(idir,2)=ixCoMmin1-nghostcellsCo-(interpolation_order-&
             1)
          ixR_p_stg_max1(idir,2)=ixCoMmax1
          ixR_p_stg_min1(idir,3)=ixCoMmax1+1
          ixR_p_stg_max1(idir,3)=ixCoMmax1+nghostcellsCo+&
             (interpolation_order-1)
          
        end if
       
      end do
      ! calculate sizes for buffer arrays for siblings
      do i1=-1,1
         ! Staggered (face-allocated) variables
         do idir=1,ndim
           sizes_srl_send_stg(idir,i1)=(ixS_srl_stg_max1(idir,&
              i1)-ixS_srl_stg_min1(idir,i1)+1)
           sizes_srl_recv_stg(idir,i1)=(ixR_srl_stg_max1(idir,&
              i1)-ixR_srl_stg_min1(idir,i1)+1)
           sizes_r_send_stg(idir,i1)=(ixS_r_stg_max1(idir,&
              i1)-ixS_r_stg_min1(idir,i1)+1)
         end do
         sizes_srl_send_total(i1)=sum(sizes_srl_send_stg(:,i1))
         sizes_srl_recv_total(i1)=sum(sizes_srl_recv_stg(:,i1))
         sizes_r_send_total(i1)=sum(sizes_r_send_stg(:,i1))
      end do

      do i1=0,3
         ! Staggered (face-allocated) variables
           do idir=1,ndim
             sizes_r_recv_stg(idir,i1)=(ixR_r_stg_max1(idir,&
                i1)-ixR_r_stg_min1(idir,i1)+1)
             sizes_p_send_stg(idir,i1)=(ixS_p_stg_max1(idir,&
                i1)-ixS_p_stg_min1(idir,i1)+1)
             sizes_p_recv_stg(idir,i1)=(ixR_p_stg_max1(idir,&
                i1)-ixR_p_stg_min1(idir,i1)+1)
           end do
           sizes_r_recv_total(i1)=sum(sizes_r_recv_stg(:,i1))
           sizes_p_send_total(i1)=sum(sizes_p_send_stg(:,i1))
           sizes_p_recv_total(i1)=sum(sizes_p_recv_stg(:,i1))
      end do
    end if
    if(.not.stagger_grid .or. physics_type=='mf') then
      ! extend index range to physical boundary
      
      ixS_srl_min1(-1,0)=1
      ixS_srl_min1( 1,0)=ixMmin1
      ixS_srl_min1( 2,0)=1
      ixS_srl_max1(-1,0)=ixMmax1
      ixS_srl_max1( 1,0)=ixGmax1
      ixS_srl_max1( 2,0)=ixGmax1
       
      ixR_srl_min1(-1,0)=1
      ixR_srl_min1( 1,0)=ixMmin1
      ixR_srl_min1( 2,0)=1
      ixR_srl_max1(-1,0)=ixMmax1
      ixR_srl_max1( 1,0)=ixGmax1
      ixR_srl_max1( 2,0)=ixGmax1
      
      ixS_r_min1(-1,0)=1
      ixS_r_min1( 1,0)=ixCoMmin1
      ixS_r_max1(-1,0)=ixCoMmax1
      ixS_r_max1( 1,0)=ixCoGmax1
      
      ixR_r_min1(-1,1)=1
      ixR_r_max1(-1,1)=ixMmin1-1+nxCo1
      ixR_r_min1( 1,2)=ixMmin1+nxCo1
      ixR_r_max1( 1,2)=ixGmax1

      ixS_p_min1(-1,1)=1
      ixS_p_max1( 1,2)=ixGmax1

      ixR_p_min1(-1,1)=1
      ixR_p_max1( 1,2)=ixCoGmax1
      
    end if

  end subroutine init_bc

  subroutine create_bc_mpi_datatype(nwstart,nwbc) 
    use mod_global_parameters 

    integer, intent(in) :: nwstart, nwbc
    integer :: i1, ic1, inc1, iib1

    do i1=-1,1
      if (i1==0) cycle
      do iib1=-1,2
         call get_bc_comm_type(type_send_srl(iib1,i1),ixS_srl_min1(iib1,i1),&
            ixS_srl_max1(iib1,i1),ixGlo1,ixGhi1,nwstart,nwbc)
         call get_bc_comm_type(type_recv_srl(iib1,i1),ixR_srl_min1(iib1,i1),&
            ixR_srl_max1(iib1,i1),ixGlo1,ixGhi1,nwstart,nwbc)
         if (iib1==2) cycle
         call get_bc_comm_type(type_send_r(iib1,i1),  ixS_r_min1(iib1,i1),&
            ixS_r_max1(iib1,i1),ixCoGmin1,ixCoGmax1,nwstart,nwbc)
         do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
            inc1=2*i1+ic1
            call get_bc_comm_type(type_recv_r(iib1,inc1),ixR_r_min1(iib1,inc1),&
               ixR_r_max1(iib1,inc1), ixGlo1,ixGhi1,nwstart,nwbc)
            call get_bc_comm_type(type_send_p(iib1,inc1),ixS_p_min1(iib1,inc1),&
               ixS_p_max1(iib1,inc1), ixGlo1,ixGhi1,nwstart,nwbc)
            call get_bc_comm_type(type_recv_p(iib1,inc1),ixR_p_min1(iib1,inc1),&
               ixR_p_max1(iib1,inc1),ixCoGmin1,ixCoGmax1,nwstart,nwbc)
         end do
      end do
    end do
  
  end subroutine create_bc_mpi_datatype

  subroutine get_bc_comm_type(comm_type,ixmin1,ixmax1,ixGmin1,ixGmax1,nwstart,&
     nwbc)
    use mod_global_parameters 
  
    integer, intent(inout) :: comm_type
    integer, intent(in) :: ixmin1,ixmax1, ixGmin1,ixGmax1, nwstart, nwbc
    
    integer, dimension(ndim+1) :: fullsize, subsize, start

    fullsize(1)=ixGmax1;
    fullsize(ndim+1)=nw
    subsize(1)=ixmax1-ixmin1+1;
    subsize(ndim+1)=nwbc
    start(1)=ixmin1-1;
    start(ndim+1)=nwstart-1
    
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,fullsize,subsize,start,&
       MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
    call MPI_TYPE_COMMIT(comm_type,ierrmpi)
    
  end subroutine get_bc_comm_type

  subroutine put_bc_comm_types()
    use mod_global_parameters 
 
    integer :: i1, ic1, inc1, iib1

    do i1=-1,1
       if (i1==0) cycle
       do iib1=-1,2
           call MPI_TYPE_FREE(type_send_srl(iib1,i1),ierrmpi)
           call MPI_TYPE_FREE(type_recv_srl(iib1,i1),ierrmpi)
           if (levmin==levmax) cycle
           if (iib1==2) cycle
           call MPI_TYPE_FREE(type_send_r(iib1,i1),ierrmpi)
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
              inc1=2*i1+ic1
              call MPI_TYPE_FREE(type_recv_r(iib1,inc1),ierrmpi)
              call MPI_TYPE_FREE(type_send_p(iib1,inc1),ierrmpi)
              call MPI_TYPE_FREE(type_recv_p(iib1,inc1),ierrmpi)
           end do
       end do
    end do
  
  end subroutine put_bc_comm_types

  !> do update ghost cells of all blocks including physical boundaries
  subroutine getbc(time,qdt,psb,nwstart,nwbc,req_diag)
    use mod_global_parameters
    use mod_physics

    double precision, intent(in)      :: time, qdt
    type(state), target               :: psb(max_blocks)
    integer, intent(in)               :: nwstart ! Fill from nwstart
    integer, intent(in)               :: nwbc    ! Number of variables to fill
    logical, intent(in), optional     :: req_diag !If false, skip diagonal ghost cells

    double precision :: time_bcin
    integer :: ipole, nwhead, nwtail
    integer :: iigrid, igrid, ineighbor, ipe_neighbor, isizes
    integer :: ixRmin1,ixRmax1, ixSmin1,ixSmax1
    integer :: i1, n_i1, ic1, inc1, n_inc1, iib1, idir
    ! store physical boundary indicating index
    integer :: idphyb(ndim,max_blocks)
    integer :: isend_buf(npwbuf), ipwbuf, nghostcellsco
    ! index pointer for buffer arrays as a start for a segment
    integer :: ibuf_start, ibuf_next
    ! shapes of reshape
    integer, dimension(1) :: shapes
    logical  :: req_diagonal
    type(wbuffer) :: pwbuf(npwbuf)

    time_bcin=MPI_WTIME()

    nwhead=nwstart
    nwtail=nwstart+nwbc-1

    req_diagonal = .true.
    if (present(req_diag)) req_diagonal = req_diag

    ! fill internal physical boundary
    if (internalboundary) then 
       call getintbc(time,ixGlo1,ixGhi1)
    end if

    ! fill physical-boundary ghost cells before internal ghost-cell values exchange
    if(bcphys.and. .not.stagger_grid) then
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        call fill_boundary_before_gc(igrid)
      end do
      !$OMP END PARALLEL DO
    end if

    ! prepare coarse values to send to coarser neighbors
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      if(any(neighbor_type(:,igrid)==neighbor_coarse)) then
        call coarsen_grid(psb(igrid),ixGlo1,ixGhi1,ixMmin1,ixMmax1,psc(igrid),&
           ixCoGmin1,ixCoGmax1,ixCoMmin1,ixCoMmax1)
       do i1=-1,1
          if(skip_direction([ i1 ])) cycle
          if(neighbor_type(i1,igrid)==neighbor_coarse) call &
             fill_coarse_boundary(igrid,i1)
       end do
      end if
    end do
    !$OMP END PARALLEL DO

    ! default : no singular axis
    ipole=0
    irecv_c=0
    isend_c=0
    isend_buf=0
    ipwbuf=1

    if(stagger_grid) then
      ibuf_recv_srl=1
      ibuf_recv_r=1
      ibuf_recv_p=1
      ibuf_send_srl=1
      ibuf_send_r=1
      ibuf_send_p=1
      irecv_srl=0
      irecv_r=0
      irecv_p=0
      isend_srl=0
      isend_r=0
      isend_p=0
    end if

    ! MPI receive ghost-cell values from sibling blocks and finer neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      call identifyphysbound(ps(igrid),iib1)
      idphyb(1,igrid)=iib1;
      do i1=-1,1
         if (skip_direction([ i1 ])) cycle
         select case (neighbor_type(i1,igrid))
         case (neighbor_sibling)
            call bc_recv_srl
         case (neighbor_fine)
            call bc_recv_restrict
         end select
      end do
    end do

    ! MPI send ghost-cell values to sibling blocks and coarser neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      iib1=idphyb(1,igrid);
      do i1=-1,1
         if(skip_direction([ i1 ])) cycle
         select case (neighbor_type(i1,igrid))
         case (neighbor_sibling)
            call bc_send_srl
         case (neighbor_coarse)
            call bc_send_restrict
         end select
      end do
    end do

    ! fill ghost-cell values of sibling blocks and coarser neighbors in the same processor
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid,iib1)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      iib1=idphyb(1,igrid);
      do i1=-1,1
        if(skip_direction([ i1 ])) cycle
         select case (neighbor_type(i1,igrid))
         case(neighbor_sibling)
           call bc_fill_srl(igrid,i1,iib1)
         case(neighbor_coarse)
           call bc_fill_restrict(igrid,i1,iib1)
         end select
      end do
    end do
    !$OMP END PARALLEL DO

    call MPI_WAITALL(irecv_c,recvrequest_c_sr,recvstatus_c_sr,ierrmpi)
    call MPI_WAITALL(isend_c,sendrequest_c_sr,sendstatus_c_sr,ierrmpi)

    if(stagger_grid) then
      call MPI_WAITALL(nrecv_bc_srl,recvrequest_srl,recvstatus_srl,ierrmpi)
      call MPI_WAITALL(nsend_bc_srl,sendrequest_srl,sendstatus_srl,ierrmpi)
      call MPI_WAITALL(nrecv_bc_r,recvrequest_r,recvstatus_r,ierrmpi)
      call MPI_WAITALL(nsend_bc_r,sendrequest_r,sendstatus_r,ierrmpi)
      ! unpack the received data from sibling blocks and finer neighbors to fill ghost-cell staggered values
      ibuf_recv_srl=1
      ibuf_recv_r=1
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        iib1=idphyb(1,igrid);
       do i1=-1,1
          if (skip_direction([ i1 ])) cycle
          select case (neighbor_type(i1,igrid))
          case (neighbor_sibling)
             call bc_fill_srl_stg
          case (neighbor_fine)
             call bc_fill_restrict_stg
          end select
       end do
      end do
    end if

    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
    end do

    irecv_c=0
    isend_c=0
    isend_buf=0
    ipwbuf=1

    ! MPI receive ghost-cell values from coarser neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      iib1=idphyb(1,igrid);
      do i1=-1,1
         if (skip_direction([ i1 ])) cycle
         if (neighbor_type(i1,igrid)==neighbor_coarse) call bc_recv_prolong
      end do
    end do
    ! MPI send ghost-cell values to finer neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      iib1=idphyb(1,igrid);
      do i1=-1,1
         if (skip_direction([ i1 ])) cycle
         if (neighbor_type(i1,igrid)==neighbor_fine) call bc_send_prolong
      end do
    end do

    ! fill coarse ghost-cell values of finer neighbors in the same processor
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid,iib1)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      iib1=idphyb(1,igrid);
      do i1=-1,1
         if (skip_direction([ i1 ])) cycle
         if (neighbor_type(i1,igrid)==neighbor_fine) call &
            bc_fill_prolong(igrid,i1,iib1)
      end do
    end do
    !$OMP END PARALLEL DO

    call MPI_WAITALL(irecv_c,recvrequest_c_p,recvstatus_c_p,ierrmpi)
    call MPI_WAITALL(isend_c,sendrequest_c_p,sendstatus_c_p,ierrmpi)

    if(stagger_grid) then
      call MPI_WAITALL(nrecv_bc_p,recvrequest_p,recvstatus_p,ierrmpi)
      call MPI_WAITALL(nsend_bc_p,sendrequest_p,sendstatus_p,ierrmpi)

      ! fill coarser representative ghost cells after receipt
      ibuf_recv_p=1
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        iib1=idphyb(1,igrid);
        do i1=-1,1
           if (skip_direction([ i1 ])) cycle
           if(neighbor_type(i1,igrid)==neighbor_coarse) call &
              bc_fill_prolong_stg
        end do
      end do
    end if
    ! do prolongation on the ghost-cell values based on the received coarse values from coarser neighbors
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      call gc_prolong(igrid)
    end do
    !$OMP END PARALLEL DO

    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
    end do

    ! fill physical boundary ghost cells after internal ghost-cell values exchange
    if(bcphys.and.stagger_grid) then
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        call fill_boundary_after_gc(igrid)
      end do
      !$OMP END PARALLEL DO
    end if

     ! modify normal component of magnetic field to fix divB=0 
    if(bcphys.and.associated(phys_boundary_adjust)) then
      !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        call phys_boundary_adjust(igrid,psb)
      end do
      !$OMP END PARALLEL DO
    end if

    time_bc=time_bc+(MPI_WTIME()-time_bcin)
    
    contains

      logical function skip_direction(dir)
        integer, intent(in) :: dir(1)

        if (all(dir == 0)) then
           skip_direction = .true.
        else if (.not. req_diagonal .and. count(dir /= 0) > 1) then
           skip_direction = .true.
        else
           skip_direction = .false.
        end if
      end function skip_direction

      !> Physical boundary conditions
      subroutine fill_boundary_before_gc(igrid)

        integer, intent(in) :: igrid

        integer :: idims,iside,i1,kmin1,kmax1,ixBmin1,ixBmax1

        block=>psb(igrid)
        dxlevel(1)=rnode(rpdx1_,igrid);
        do idims=1,ndim
          ! to avoid using as yet unknown corner info in more than 1D, we
          ! fill only interior mesh ranges of the ghost cell ranges at first,
          ! and progressively enlarge the ranges to include corners later
          
           kmin1=merge(0, 1, idims==1)
           kmax1=merge(0, 1, idims==1)
           ixBmin1=ixGlo1+kmin1*nghostcells
           ixBmax1=ixGhi1-kmax1*nghostcells
          
          
          
          do iside=1,2
            i1=kr(1,idims)*(2*iside-3);
            if (aperiodB(idims)) then
              if (neighbor_type(i1,igrid) /= neighbor_boundary .and. .not. &
                 psb(igrid)%is_physical_boundary(2*idims-2+iside)) cycle
            else
              if (neighbor_type(i1,igrid) /= neighbor_boundary) cycle
            end if
            call bc_phys(iside,idims,time,qdt,psb(igrid),ixGlo1,ixGhi1,ixBmin1,&
               ixBmax1)
          end do
        end do

      end subroutine fill_boundary_before_gc

      !> Physical boundary conditions
      subroutine fill_boundary_after_gc(igrid)

        integer, intent(in) :: igrid

        integer :: idims,iside,i1,kmin1,kmax1,ixBmin1,ixBmax1

        block=>psb(igrid)
        dxlevel(1)=rnode(rpdx1_,igrid);
        do idims=1,ndim
          ! to avoid using as yet unknown corner info in more than 1D, we
          ! fill only interior mesh ranges of the ghost cell ranges at first,
          ! and progressively enlarge the ranges to include corners later
          kmin1=0; kmax1=0;
          
          
          ixBmin1=ixGlo1+kmin1*nghostcells;
          ixBmax1=ixGhi1-kmax1*nghostcells;
          do iside=1,2
            i1=kr(1,idims)*(2*iside-3);
            if (aperiodB(idims)) then 
              if (neighbor_type(i1,igrid) /= neighbor_boundary .and. .not. &
                 psb(igrid)%is_physical_boundary(2*idims-2+iside)) cycle
            else 
              if (neighbor_type(i1,igrid) /= neighbor_boundary) cycle
            end if
            call bc_phys(iside,idims,time,qdt,psb(igrid),ixGlo1,ixGhi1,ixBmin1,&
               ixBmax1)
          end do
        end do

      end subroutine fill_boundary_after_gc

      !> Receive from sibling at same refinement level
      subroutine bc_recv_srl

        ipe_neighbor=neighbor(2,i1,igrid)
        if (ipe_neighbor/=mype) then
           irecv_c=irecv_c+1
           itag=(3**1+4**1)*(igrid-1)+(i1+1)*3**(1-1)
           call MPI_IRECV(psb(igrid)%w,1,type_recv_srl(iib1,i1), ipe_neighbor,&
              itag,icomm,recvrequest_c_sr(irecv_c),ierrmpi)
           if(stagger_grid) then
             irecv_srl=irecv_srl+1
             call MPI_IRECV(recvbuffer_srl(ibuf_recv_srl),&
                sizes_srl_recv_total(i1),MPI_DOUBLE_PRECISION, ipe_neighbor,&
                itag,icomm,recvrequest_srl(irecv_srl),ierrmpi)
             ibuf_recv_srl=ibuf_recv_srl+sizes_srl_recv_total(i1)
           end if
        end if

      end subroutine bc_recv_srl

      !> Receive from fine neighbor
      subroutine bc_recv_restrict

        do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
           inc1=2*i1+ic1
           ipe_neighbor=neighbor_child(2,inc1,igrid)
           if (ipe_neighbor/=mype) then
              irecv_c=irecv_c+1
              itag=(3**1+4**1)*(igrid-1)+3**1+inc1*4**(1-1)
              call MPI_IRECV(psb(igrid)%w,1,type_recv_r(iib1,inc1),&
                  ipe_neighbor,itag,icomm,recvrequest_c_sr(irecv_c),ierrmpi)
              if(stagger_grid) then
                irecv_r=irecv_r+1
                call MPI_IRECV(recvbuffer_r(ibuf_recv_r),&
                   sizes_r_recv_total(inc1), MPI_DOUBLE_PRECISION,ipe_neighbor,&
                   itag, icomm,recvrequest_r(irecv_r),ierrmpi)
                ibuf_recv_r=ibuf_recv_r+sizes_r_recv_total(inc1)
              end if
           end if
        end do

      end subroutine bc_recv_restrict

      !> Send to sibling at same refinement level
      subroutine bc_send_srl

        ipe_neighbor=neighbor(2,i1,igrid)

        if(ipe_neighbor/=mype) then
          ineighbor=neighbor(1,i1,igrid)
          ipole=neighbor_pole(i1,igrid)
          if(ipole==0) then
            n_i1=-i1;
            isend_c=isend_c+1
            itag=(3**1+4**1)*(ineighbor-1)+(n_i1+1)*3**(1-1)
            call MPI_ISEND(psb(igrid)%w,1,type_send_srl(iib1,i1), ipe_neighbor,&
               itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
            if(stagger_grid) then
              ibuf_start=ibuf_send_srl
              do idir=1,ndim
                ixSmin1=ixS_srl_stg_min1(idir,i1)
                ixSmax1=ixS_srl_stg_max1(idir,i1);
                ibuf_next=ibuf_start+sizes_srl_send_stg(idir,i1)
                shapes=(/sizes_srl_send_stg(idir,i1)/)
                sendbuffer_srl(ibuf_start:ibuf_next-&
                   1)=reshape(psb(igrid)%ws(ixSmin1:ixSmax1,idir),shapes)
                ibuf_start=ibuf_next
              end do
              isend_srl=isend_srl+1
              call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),&
                 sizes_srl_send_total(i1),MPI_DOUBLE_PRECISION, ipe_neighbor,&
                 itag,icomm,sendrequest_srl(isend_srl),ierrmpi)
              ibuf_send_srl=ibuf_next
            end if
          else
            ixSmin1=ixS_srl_min1(iib1,i1);ixSmax1=ixS_srl_max1(iib1,i1);
            select case (ipole)
            case (1)
               n_i1=i1;
            end select
            if (isend_buf(ipwbuf)/=0) then
               call MPI_WAIT(sendrequest_c_sr(isend_buf(ipwbuf)),&
                   sendstatus_c_sr(:,isend_buf(ipwbuf)),ierrmpi)
               deallocate(pwbuf(ipwbuf)%w)
            end if
            allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,nwhead:nwtail))
            call pole_buffer(pwbuf(ipwbuf)%w,ixSmin1,ixSmax1,ixSmin1,ixSmax1,&
               psb(igrid)%w,ixGlo1,ixGhi1,ixSmin1,ixSmax1)
            isend_c=isend_c+1
            isend_buf(ipwbuf)=isend_c
            itag=(3**1+4**1)*(ineighbor-1)+(n_i1+1)*3**(1-1)
            isizes=(ixSmax1-ixSmin1+1)*nwbc
            call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
                ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
            ipwbuf=1+modulo(ipwbuf,npwbuf)
            if(stagger_grid) then
              ibuf_start=ibuf_send_srl
              do idir=1,ndim
                ixSmin1=ixS_srl_stg_min1(idir,i1)
                ixSmax1=ixS_srl_stg_max1(idir,i1);
                ibuf_next=ibuf_start+sizes_srl_send_stg(idir,i1)
                shapes=(/sizes_srl_send_stg(idir,i1)/)
                sendbuffer_srl(ibuf_start:ibuf_next-&
                   1)=reshape(psb(igrid)%ws(ixSmin1:ixSmax1,idir),shapes)
                ibuf_start=ibuf_next
              end do
              isend_srl=isend_srl+1
              call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),&
                 sizes_srl_send_total(i1),MPI_DOUBLE_PRECISION, ipe_neighbor,&
                 itag,icomm,sendrequest_srl(isend_srl),ierrmpi)
              ibuf_send_srl=ibuf_next
            end if
          end if
        end if

      end subroutine bc_send_srl

      subroutine bc_fill_srl(igrid,i1,iib1)
        integer, intent(in) :: igrid,i1,iib1
        integer :: ineighbor,ipe_neighbor,ipole,ixSmin1,ixSmax1,ixRmin1,&
           ixRmax1,n_i1,idir

        ipe_neighbor=neighbor(2,i1,igrid)
        if(ipe_neighbor==mype) then
          ineighbor=neighbor(1,i1,igrid)
          ipole=neighbor_pole(i1,igrid)
          if(ipole==0) then
            n_i1=-i1;
            ixSmin1=ixS_srl_min1(iib1,i1);ixSmax1=ixS_srl_max1(iib1,i1);
            ixRmin1=ixR_srl_min1(iib1,n_i1);ixRmax1=ixR_srl_max1(iib1,n_i1);
            psb(ineighbor)%w(ixRmin1:ixRmax1,&
               nwhead:nwtail)=psb(igrid)%w(ixSmin1:ixSmax1,nwhead:nwtail)
            if(stagger_grid) then
              do idir=1,ndim
                ixSmin1=ixS_srl_stg_min1(idir,i1)
                ixSmax1=ixS_srl_stg_max1(idir,i1);
                ixRmin1=ixR_srl_stg_min1(idir,n_i1)
                ixRmax1=ixR_srl_stg_max1(idir,n_i1);
                psb(ineighbor)%ws(ixRmin1:ixRmax1,&
                   idir)=psb(igrid)%ws(ixSmin1:ixSmax1,idir)
              end do
            end if
          else
            ixSmin1=ixS_srl_min1(iib1,i1);ixSmax1=ixS_srl_max1(iib1,i1);
            select case (ipole)
            case (1)
              n_i1=i1;
            end select
            ixRmin1=ixR_srl_min1(iib1,n_i1);ixRmax1=ixR_srl_max1(iib1,n_i1);
            call pole_copy(psb(ineighbor)%w,ixGlo1,ixGhi1,ixRmin1,ixRmax1,&
               psb(igrid)%w,ixGlo1,ixGhi1,ixSmin1,ixSmax1,ipole)
            if(stagger_grid) then
              do idir=1,ndim
                ixSmin1=ixS_srl_stg_min1(idir,i1)
                ixSmax1=ixS_srl_stg_max1(idir,i1);
                ixRmin1=ixR_srl_stg_min1(idir,n_i1)
                ixRmax1=ixR_srl_stg_max1(idir,n_i1);
                call pole_copy_stg(psb(ineighbor)%ws,ixGslo1,ixGshi1,ixRmin1,&
                   ixRmax1,psb(igrid)%ws,ixGslo1,ixGshi1,ixSmin1,ixSmax1,idir,&
                   ipole)
              end do
            end if
          end if
        end if

      end subroutine bc_fill_srl

      subroutine fill_coarse_boundary(igrid,i1)
        integer, intent(in) :: igrid,i1

        integer :: idims,iside,kmin1,kmax1,ixBmin1,ixBmax1,ii1

        if(phyboundblock(igrid).and..not.stagger_grid.and.bcphys) then
          ! to use block in physical boundary setup for coarse representative
          block=>psc(igrid)
          ! filling physical boundary ghost cells of a coarser representative block for
          ! sending swap region with width of nghostcells to its coarser neighbor
          do idims=1,ndim
             ! to avoid using as yet unknown corner info in more than 1D, we
             ! fill only interior mesh ranges of the ghost cell ranges at first,
             ! and progressively enlarge the ranges to include corners later
             kmin1=merge(0, 1, idims==1)
             kmax1=merge(0, 1, idims==1)
             ixBmin1=ixCoGmin1+kmin1*nghostcells
             ixBmax1=ixCoGmax1-kmax1*nghostcells
             
             
             if(i1==-1) then
               ixBmin1=ixCoGmin1+nghostcells
               ixBmax1=ixCoGmin1+2*nghostcells-1
             else if(i1==1) then
               ixBmin1=ixCoGmax1-2*nghostcells+1
               ixBmax1=ixCoGmax1-nghostcells
             end if
             do iside=1,2
                ii1=kr(1,idims)*(2*iside-3);
                if (abs(i1)==1.and.abs(ii1)==1) cycle
                if (neighbor_type(ii1,igrid)/=neighbor_boundary) cycle
                call bc_phys(iside,idims,time,0.d0,psc(igrid),ixCoGmin1,&
                   ixCoGmax1,ixBmin1,ixBmax1)
             end do
          end do
        end if

      end subroutine fill_coarse_boundary

      !> Send to coarser neighbor
      subroutine bc_send_restrict

        ipe_neighbor=neighbor(2,i1,igrid)
        if(ipe_neighbor/=mype) then
          ic1=1+modulo(node(pig1_,igrid)-1,2);
          if(.not.(i1==0.or.i1==2*ic1-3)) return
          ineighbor=neighbor(1,i1,igrid)
          ipole=neighbor_pole(i1,igrid)
          if(ipole==0) then
            n_inc1=-2*i1+ic1;
            isend_c=isend_c+1
            itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
            call MPI_ISEND(psc(igrid)%w,1,type_send_r(iib1,i1), ipe_neighbor,&
               itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
            if(stagger_grid) then 
              ibuf_start=ibuf_send_r
              do idir=1,ndim
                ixSmin1=ixS_r_stg_min1(idir,i1)
                ixSmax1=ixS_r_stg_max1(idir,i1);
                ibuf_next=ibuf_start+sizes_r_send_stg(idir,i1)
                shapes=(/sizes_r_send_stg(idir,i1)/)
                sendbuffer_r(ibuf_start:ibuf_next-&
                   1)=reshape(psc(igrid)%ws(ixSmin1:ixSmax1,idir),shapes)
                ibuf_start=ibuf_next
              end do
              isend_r=isend_r+1
              call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i1),&
                 MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                 sendrequest_r(isend_r),ierrmpi)
              ibuf_send_r=ibuf_next
            end if
          else
            ixSmin1=ixS_r_min1(iib1,i1);ixSmax1=ixS_r_max1(iib1,i1);
            select case (ipole)
            case (1)
               n_inc1=2*i1+(3-ic1);
            end select
            if(isend_buf(ipwbuf)/=0) then
              call MPI_WAIT(sendrequest_c_sr(isend_buf(ipwbuf)),&
                  sendstatus_c_sr(:,isend_buf(ipwbuf)),ierrmpi)
              deallocate(pwbuf(ipwbuf)%w)
            end if
            allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,nwhead:nwtail))
            call pole_buffer(pwbuf(ipwbuf)%w,ixSmin1,ixSmax1,ixSmin1,ixSmax1,&
               psc(igrid)%w,ixCoGmin1,ixCoGmax1,ixSmin1,ixSmax1)
            isend_c=isend_c+1
            isend_buf(ipwbuf)=isend_c
            itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
            isizes=(ixSmax1-ixSmin1+1)*nwbc
            call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
                ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
            ipwbuf=1+modulo(ipwbuf,npwbuf)
            if(stagger_grid) then 
              ibuf_start=ibuf_send_r
              do idir=1,ndim
                ixSmin1=ixS_r_stg_min1(idir,i1)
                ixSmax1=ixS_r_stg_max1(idir,i1);
                ibuf_next=ibuf_start+sizes_r_send_stg(idir,i1)
                shapes=(/sizes_r_send_stg(idir,i1)/)
                sendbuffer_r(ibuf_start:ibuf_next-&
                   1)=reshape(psc(igrid)%ws(ixSmin1:ixSmax1,idir),shapes)
                ibuf_start=ibuf_next
              end do
              isend_r=isend_r+1
              call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i1),&
                 MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                 sendrequest_r(isend_r),ierrmpi)
              ibuf_send_r=ibuf_next
            end if
          end if
        end if

      end subroutine bc_send_restrict

      !> fill coarser neighbor's ghost cells
      subroutine bc_fill_restrict(igrid,i1,iib1)
        integer, intent(in) :: igrid,i1,iib1

        integer :: ic1,n_inc1,ixSmin1,ixSmax1,ixRmin1,ixRmax1,ipe_neighbor,&
           ineighbor,ipole,idir

        ipe_neighbor=neighbor(2,i1,igrid)
        if(ipe_neighbor==mype) then
          ic1=1+modulo(node(pig1_,igrid)-1,2);
          if(.not.(i1==0.or.i1==2*ic1-3)) return
          ineighbor=neighbor(1,i1,igrid)
          ipole=neighbor_pole(i1,igrid)
          if(ipole==0) then
            n_inc1=-2*i1+ic1;
            ixSmin1=ixS_r_min1(iib1,i1);ixSmax1=ixS_r_max1(iib1,i1);
            ixRmin1=ixR_r_min1(iib1,n_inc1);ixRmax1=ixR_r_max1(iib1,n_inc1);
            psb(ineighbor)%w(ixRmin1:ixRmax1,&
               nwhead:nwtail)=psc(igrid)%w(ixSmin1:ixSmax1,nwhead:nwtail)
            if(stagger_grid) then
              do idir=1,ndim
                 ixSmin1=ixS_r_stg_min1(idir,i1)
                 ixSmax1=ixS_r_stg_max1(idir,i1);
                 ixRmin1=ixR_r_stg_min1(idir,n_inc1)
                 ixRmax1=ixR_r_stg_max1(idir,n_inc1);
                 psb(ineighbor)%ws(ixRmin1:ixRmax1,&
                    idir)=psc(igrid)%ws(ixSmin1:ixSmax1,idir)
              end do
            end if
          else
            ixSmin1=ixS_r_min1(iib1,i1);ixSmax1=ixS_r_max1(iib1,i1);
            select case (ipole)
            case (1)
              n_inc1=2*i1+(3-ic1);
            end select
            ixRmin1=ixR_r_min1(iib1,n_inc1);ixRmax1=ixR_r_max1(iib1,n_inc1);
            call pole_copy(psb(ineighbor)%w,ixGlo1,ixGhi1,ixRmin1,ixRmax1,&
               psc(igrid)%w,ixCoGmin1,ixCoGmax1,ixSmin1,ixSmax1,ipole)
            if(stagger_grid) then
              do idir=1,ndim
                ixSmin1=ixS_r_stg_min1(idir,i1)
                ixSmax1=ixS_r_stg_max1(idir,i1);
                ixRmin1=ixR_r_stg_min1(idir,n_inc1)
                ixRmax1=ixR_r_stg_max1(idir,n_inc1);
                !! Fill ghost cells
                call pole_copy_stg(psb(ineighbor)%ws,ixGslo1,ixGshi1,ixRmin1,&
                   ixRmax1,psc(igrid)%ws,ixCoGsmin1,ixCoGsmax1,ixSmin1,ixSmax1,&
                   idir,ipole)
              end do
            end if
          end if
        end if

      end subroutine bc_fill_restrict

      !> fill siblings ghost cells with received data
      subroutine bc_fill_srl_stg
        double precision :: tmp(ixGslo1:ixGshi1)
        integer :: ixSmin1,ixSmax1,ixRmin1,ixRmax1,n_i1,ixSsyncmin1,&
           ixSsyncmax1,ixRsyncmin1,ixRsyncmax1
        integer :: idir

        ipe_neighbor=neighbor(2,i1,igrid)
        if(ipe_neighbor/=mype) then
          ineighbor=neighbor(1,i1,igrid)
          ipole=neighbor_pole(i1,igrid)

        !! Now the special treatment of the pole is done here, at the receive step
          if (ipole==0) then    
            ixRmin1=ixR_srl_min1(iib1,i1);ixRmax1=ixR_srl_max1(iib1,i1);
            !! Unpack the buffer and fill the ghost cells
            n_i1=-i1;
            do idir=1,ndim
              ixSmin1=ixS_srl_stg_min1(idir,n_i1)
              ixSmax1=ixS_srl_stg_max1(idir,n_i1);
              ixRmin1=ixR_srl_stg_min1(idir,i1)
              ixRmax1=ixR_srl_stg_max1(idir,i1);
              ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i1)
              tmp(ixSmin1:ixSmax1) = reshape(source=recvbuffer_srl(&
                 ibuf_recv_srl:ibuf_next-1),&
                 shape=shape(psb(igrid)%ws(ixSmin1:ixSmax1,idir)))       
              psb(igrid)%ws(ixRmin1:ixRmax1,idir) = tmp(ixSmin1:ixSmax1)
              ibuf_recv_srl=ibuf_next
            end do
          else ! There is a pole
            select case (ipole)
            case (1)
               n_i1=i1;
            end select
            pole_buf%ws=zero
            do idir=1,ndim
             ixRmin1=ixR_srl_stg_min1(idir,i1)
             ixRmax1=ixR_srl_stg_max1(idir,i1);
             ixSmin1=ixS_srl_stg_min1(idir,n_i1)
             ixSmax1=ixS_srl_stg_max1(idir,n_i1);
             ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i1)
             pole_buf%ws(ixSmin1:ixSmax1,&
                idir)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),&
                shape=shape(psb(igrid)%ws(ixSmin1:ixSmax1,idir)))
             ibuf_recv_srl=ibuf_next
             call pole_copy_stg(psb(igrid)%ws,ixGslo1,ixGshi1,ixRmin1,ixRmax1,&
                pole_buf%ws,ixGslo1,ixGshi1,ixSmin1,ixSmax1,idir,ipole)
            end do
          end if
        end if

      end subroutine bc_fill_srl_stg

      subroutine indices_for_syncing(idir,i1,ixRmin1,ixRmax1,ixSmin1,ixSmax1,&
         ixRsyncmin1,ixRsyncmax1,ixSsyncmin1,ixSsyncmax1)
        integer, intent(in)       :: i1,idir
        integer, intent(inout)    :: ixRmin1,ixRmax1,ixSmin1,ixSmax1
        integer, intent(out)      :: ixRsyncmin1,ixRsyncmax1,ixSsyncmin1,&
           ixSsyncmax1
      
        ixRsyncmin1=ixRmin1;ixRsyncmax1=ixRmax1;
        ixSsyncmin1=ixSmin1;ixSsyncmax1=ixSmax1;
        
        
        if (i1 == -1 .and. idir == 1) then
           ixRsyncmin1 = ixRmax1
           ixRsyncmax1 = ixRmax1
           ixSsyncmin1 = ixSmax1
           ixSsyncmax1 = ixSmax1
           ixRmax1 = ixRmax1 - 1
           ixSmax1 = ixSmax1 - 1
        else if (i1 == 1 .and. idir == 1) then
           ixRsyncmin1 = ixRmin1
           ixRsyncmax1 = ixRmin1
           ixSsyncmin1 = ixSmin1
           ixSsyncmax1 = ixSmin1
           ixRmin1 = ixRmin1 + 1
           ixSmin1 = ixSmin1 + 1
        end if
        

      end subroutine indices_for_syncing

      !> fill restricted ghost cells after receipt
      subroutine bc_fill_restrict_stg

        ipole=neighbor_pole(i1,igrid)
        if (ipole==0) then
          ! Loop over the children ic^D to and their neighbors inc^D
          do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
             inc1=2*i1+ic1
             ipe_neighbor=neighbor_child(2,inc1,igrid)
             if(ipe_neighbor/=mype) then
               ineighbor=neighbor_child(1,inc1,igrid)
               n_i1=-i1;
               !! Unpack the buffer and fill the ghost cells
               do idir=1,ndim
                 ixRmin1=ixR_r_stg_min1(idir,inc1)
                 ixRmax1=ixR_r_stg_max1(idir,inc1);
                 ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc1)
                 psb(igrid)%ws(ixRmin1:ixRmax1,&
                    idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                    shape=shape(psb(igrid)%ws(ixRmin1:ixRmax1,idir)))
                 ibuf_recv_r=ibuf_next
               end do
             end if
          end do
        else !! There is a pole
          do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
             inc1=2*i1+ic1
             ipe_neighbor=neighbor_child(2,inc1,igrid)
             if(ipe_neighbor/=mype) then
               ineighbor=neighbor_child(1,inc1,igrid)
               select case(ipole)
              case (1)
                 n_i1=i1;
               end select
               ixRmin1=ixR_r_min1(iib1,inc1);ixRmax1=ixR_r_max1(iib1,inc1);
               !! Unpack the buffer and fill an auxiliary array
               pole_buf%ws=zero
               do idir=1,ndim
                 ixSmin1=ixS_r_stg_min1(idir,n_i1)
                 ixSmax1=ixS_r_stg_max1(idir,n_i1);
                 ixRmin1=ixR_r_stg_min1(idir,inc1)
                 ixRmax1=ixR_r_stg_max1(idir,inc1);
                 ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc1)
                 pole_buf%ws(ixRmin1:ixRmax1,&
                    idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                    shape=shape(psb(igrid)%ws(ixRmin1:ixRmax1,idir)))
                 call pole_copy_stg(psb(igrid)%ws,ixGslo1,ixGshi1,ixRmin1,&
                    ixRmax1,pole_buf%ws,ixGslo1,ixGshi1,ixRmin1,ixRmax1,idir,&
                    ipole)
                 ibuf_recv_r=ibuf_next
               end do
             end if
          end do
        end if

      end subroutine bc_fill_restrict_stg

      !> Receive from coarse neighbor
      subroutine bc_recv_prolong

        ic1=1+modulo(node(pig1_,igrid)-1,2);
        if (.not.(i1==0.or.i1==2*ic1-3)) return

        ipe_neighbor=neighbor(2,i1,igrid)
        if (ipe_neighbor/=mype) then
           irecv_c=irecv_c+1
           inc1=ic1+i1;
           itag=(3**1+4**1)*(igrid-1)+3**1+inc1*4**(1-1)
           call MPI_IRECV(psc(igrid)%w,1,type_recv_p(iib1,inc1), ipe_neighbor,&
              itag,icomm,recvrequest_c_p(irecv_c),ierrmpi)  
           if(stagger_grid) then
             irecv_p=irecv_p+1
             call MPI_IRECV(recvbuffer_p(ibuf_recv_p),sizes_p_recv_total(inc1),&
                MPI_DOUBLE_PRECISION,ipe_neighbor,itag,icomm,&
                recvrequest_p(irecv_p),ierrmpi)
             ibuf_recv_p=ibuf_recv_p+sizes_p_recv_total(inc1)
           end if
        end if

      end subroutine bc_recv_prolong

      !> Send to finer neighbor
      subroutine bc_send_prolong
        integer :: ii1

        ipole=neighbor_pole(i1,igrid)

        do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
           inc1=2*i1+ic1
           ipe_neighbor=neighbor_child(2,inc1,igrid)
           if(ipe_neighbor/=mype) then
             ixSmin1=ixS_p_min1(iib1,inc1);ixSmax1=ixS_p_max1(iib1,inc1);
             ineighbor=neighbor_child(1,inc1,igrid)
             if(ipole==0) then
               n_i1=-i1;
               n_inc1=ic1+n_i1;
               isend_c=isend_c+1
               itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
               call MPI_ISEND(psb(igrid)%w,1,type_send_p(iib1,inc1),&
                   ipe_neighbor,itag,icomm,sendrequest_c_p(isend_c),ierrmpi)
               if(stagger_grid) then 
                 ibuf_start=ibuf_send_p
                 do idir=1,ndim
                   ixSmin1=ixS_p_stg_min1(idir,inc1)
                   ixSmax1=ixS_p_stg_max1(idir,inc1);
                   ibuf_next=ibuf_start+sizes_p_send_stg(idir,inc1)
                   shapes=(/sizes_p_send_stg(idir,inc1)/)
                   sendbuffer_p(ibuf_start:ibuf_next-&
                      1)=reshape(psb(igrid)%ws(ixSmin1:ixSmax1,idir),&
                      shapes)   
                   ibuf_start=ibuf_next
                 end do
                 isend_p=isend_p+1
                 call MPI_ISEND(sendbuffer_p(ibuf_send_p),&
                    sizes_p_send_total(inc1),MPI_DOUBLE_PRECISION,ipe_neighbor,&
                    itag, icomm,sendrequest_p(isend_p),ierrmpi)
                 ibuf_send_p=ibuf_next
               end if
             else
               select case (ipole)
               case (1)
                 n_inc1=inc1;
               end select
               if(isend_buf(ipwbuf)/=0) then
                 call MPI_WAIT(sendrequest_c_p(isend_buf(ipwbuf)),&
                     sendstatus_c_p(:,isend_buf(ipwbuf)),ierrmpi)
                 deallocate(pwbuf(ipwbuf)%w)
               end if
               allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,nwhead:nwtail))
               call pole_buffer(pwbuf(ipwbuf)%w,ixSmin1,ixSmax1,ixSmin1,&
                  ixSmax1,psb(igrid)%w,ixGlo1,ixGhi1,ixSmin1,ixSmax1)
               isend_c=isend_c+1
               isend_buf(ipwbuf)=isend_c
               itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
               isizes=(ixSmax1-ixSmin1+1)*nwbc
               call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
                   ipe_neighbor,itag,icomm,sendrequest_c_p(isend_c),ierrmpi)
               ipwbuf=1+modulo(ipwbuf,npwbuf)
               if(stagger_grid) then 
                 ibuf_start=ibuf_send_p
                 do idir=1,ndim
                   ixSmin1=ixS_p_stg_min1(idir,inc1)
                   ixSmax1=ixS_p_stg_max1(idir,inc1);
                   ibuf_next=ibuf_start+sizes_p_send_stg(idir,inc1)
                   shapes=(/sizes_p_send_stg(idir,inc1)/)
                   sendbuffer_p(ibuf_start:ibuf_next-&
                      1)=reshape(psb(igrid)%ws(ixSmin1:ixSmax1,idir),&
                      shapes)   
                   ibuf_start=ibuf_next
                 end do
                 isend_p=isend_p+1
                 call MPI_ISEND(sendbuffer_p(ibuf_send_p),&
                    sizes_p_send_total(inc1),MPI_DOUBLE_PRECISION,ipe_neighbor,&
                    itag, icomm,sendrequest_p(isend_p),ierrmpi)
                 ibuf_send_p=ibuf_next
               end if
             end if
           end if
        end do

      end subroutine bc_send_prolong

      !> Send to finer neighbor
      subroutine bc_fill_prolong(igrid,i1,iib1)
        integer, intent(in) :: igrid,i1,iib1

        integer :: ipe_neighbor,ineighbor,ixSmin1,ixSmax1,ixRmin1,ixRmax1,ic1,&
           inc1,ipole,idir

        ipole=neighbor_pole(i1,igrid)

        if(ipole==0) then
          do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
             inc1=2*i1+ic1
             ipe_neighbor=neighbor_child(2,inc1,igrid)
             if(ipe_neighbor==mype) then
               ixSmin1=ixS_p_min1(iib1,inc1);ixSmax1=ixS_p_max1(iib1,inc1);
               ineighbor=neighbor_child(1,inc1,igrid)
               ipole=neighbor_pole(i1,igrid)
               n_i1=-i1;
               n_inc1=ic1+n_i1;
               ixRmin1=ixR_p_min1(iib1,n_inc1)
               ixRmax1=ixR_p_max1(iib1,n_inc1);
               psc(ineighbor)%w(ixRmin1:ixRmax1,&
                  nwhead:nwtail) =psb(igrid)%w(ixSmin1:ixSmax1,nwhead:nwtail)
               if(stagger_grid) then
                 do idir=1,ndim
                   ixSmin1=ixS_p_stg_min1(idir,inc1)
                   ixSmax1=ixS_p_stg_max1(idir,inc1);
                   ixRmin1=ixR_p_stg_min1(idir,n_inc1)
                   ixRmax1=ixR_p_stg_max1(idir,n_inc1);
                   psc(ineighbor)%ws(ixRmin1:ixRmax1,&
                      idir)=psb(igrid)%ws(ixSmin1:ixSmax1,idir)
                 end do
               end if
             end if
          end do
        else
          do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
             inc1=2*i1+ic1
             ipe_neighbor=neighbor_child(2,inc1,igrid)
             if(ipe_neighbor==mype) then
               ixSmin1=ixS_p_min1(iib1,inc1);ixSmax1=ixS_p_max1(iib1,inc1);
               ineighbor=neighbor_child(1,inc1,igrid)
               ipole=neighbor_pole(i1,igrid)
               select case (ipole)
               case (1)
                  n_inc1=inc1;
               end select
               ixRmin1=ixR_p_min1(iib1,n_inc1)
               ixRmax1=ixR_p_max1(iib1,n_inc1);
               call pole_copy(psc(ineighbor)%w,ixCoGmin1,ixCoGmax1,ixRmin1,&
                  ixRmax1,psb(igrid)%w,ixGlo1,ixGhi1,ixSmin1,ixSmax1,ipole)
               if(stagger_grid) then
                 do idir=1,ndim
                   ixSmin1=ixS_p_stg_min1(idir,inc1)
                   ixSmax1=ixS_p_stg_max1(idir,inc1);
                   ixRmin1=ixR_p_stg_min1(idir,n_inc1)
                   ixRmax1=ixR_p_stg_max1(idir,n_inc1);
                   call pole_copy_stg(psc(ineighbor)%ws,ixCoGsmin1,ixCoGsmax1,&
                      ixRmin1,ixRmax1,psb(igrid)%ws,ixGslo1,ixGshi1,ixSmin1,&
                      ixSmax1,idir,ipole)
                 end do
               end if
             end if
          end do
        end if
      end subroutine bc_fill_prolong

      subroutine gc_prolong(igrid)
        integer, intent(in) :: igrid

        integer :: iib1,i1,idims,iside
        logical,dimension(-1:1) :: NeedProlong

        iib1=idphyb(1,igrid);
        NeedProlong=.false.
        do i1=-1,1
           if (skip_direction([ i1 ])) cycle
           if (neighbor_type(i1,igrid)==neighbor_coarse) then
             call bc_prolong(igrid,i1,iib1)
             NeedProlong(i1)=.true.
           end if
        end do
        if(stagger_grid) then
          ! Ghost cell prolongation for staggered variables
          ! must be done in a specific order.
          ! First the first neighbours, which have 2 indices=0 in 3D
          ! or one index=0 in 2D
          block=>psb(igrid)
          do idims=1,ndim
            i1=0;
            select case(idims)
           case(1)
              do i1=-1,1,2
                if (NeedProlong(i1)) call bc_prolong_stg(igrid,i1,iib1,&
                   NeedProlong)
              end do
            
            end select
          end do
          ! Then the second neighbours which have 1 index=0 in 3D
          ! (Only in 3D)
          
          ! Finally, the corners, that have no index=0
         do i1=-1,1,2
            if (NeedProlong(i1)) call bc_prolong_stg(igrid,i1,iib1,&
               NeedProlong)
         end do
        end if
      end subroutine gc_prolong

      !> fill coarser representative with data from coarser neighbors
      subroutine bc_fill_prolong_stg
        ic1=1+modulo(node(pig1_,igrid)-1,2);
        if (.not.(i1==0.or.i1==2*ic1-3)) return

        ipe_neighbor=neighbor(2,i1,igrid)
        if(ipe_neighbor/=mype) then
          ineighbor=neighbor(1,i1,igrid)
          ipole=neighbor_pole(i1,igrid)

          if (ipole==0) then   !! There is no pole 
            inc1=ic1+i1;
            ixRmin1=ixR_p_min1(iib1,inc1);ixRmax1=ixR_p_max1(iib1,inc1);
            do idir=1,ndim
              ixRmin1=ixR_p_stg_min1(idir,inc1)
              ixRmax1=ixR_p_stg_max1(idir,inc1);
              ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc1)
              psc(igrid)%ws(ixRmin1:ixRmax1,&
                 idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                 shape=shape(psc(igrid)%ws(ixRmin1:ixRmax1,idir)))
              ibuf_recv_p=ibuf_next
            end do
          else !! There is a pole
            inc1=ic1+i1;
            select case (ipole)
            case (1)
               n_inc1=2*i1+(3-ic1);
            end select
            !! Unpack the buffer and fill an auxiliary array
            pole_buf%ws=zero
            do idir=1,ndim
              ixRmin1=ixR_p_stg_min1(idir,inc1)
              ixRmax1=ixR_p_stg_max1(idir,inc1);
              ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc1)
              pole_buf%ws(ixRmin1:ixRmax1,&
                 idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                 shape=shape(psc(igrid)%ws(ixRmin1:ixRmax1,idir)))
              call pole_copy_stg(psc(igrid)%ws,ixCoGsmin1,ixCoGsmax1,ixRmin1,&
                 ixRmax1,pole_buf%ws,ixGslo1,ixGshi1,ixRmin1,ixRmax1,idir,&
                 ipole)
              ibuf_recv_p=ibuf_next
            end do
          end if
        end if

      end subroutine bc_fill_prolong_stg

      !> do prolongation for fine blocks after receipt data from coarse neighbors
      subroutine bc_prolong(igrid,i1,iib1)
        use mod_physics, only: phys_to_primitive, phys_to_conserved

        integer :: i1,iib1,igrid
        integer :: ixFimin1,ixFimax1,ixComin1,ixComax1,ii1, idims,iside,&
           ixBmin1,ixBmax1
        double precision :: dxFi1, dxCo1, xFimin1, xComin1, invdxCo1

        ixFimin1=ixR_srl_min1(iib1,i1);ixFimax1=ixR_srl_max1(iib1,i1);
        dxFi1=rnode(rpdx1_,igrid);
        dxCo1=two*dxFi1;
        invdxCo1=1.d0/dxCo1;

        ! compute the enlarged grid lower left corner coordinates
        ! these are true coordinates for an equidistant grid, 
        ! but we can temporarily also use them for getting indices 
        ! in stretched grids
        xFimin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxFi1;
        xComin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxCo1;

        if(stagger_grid.and.phyboundblock(igrid).and.bcphys) then
          block=>psc(igrid)
          do idims=1,ndim
            ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-&
               xComin1)*invdxCo1)+1-1;
            ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-&
               xComin1)*invdxCo1)+1+1;
            
            do iside=1,2
              ii1=kr(1,idims)*(2*iside-3);
              if(neighbor_type(ii1,igrid)/=neighbor_boundary) cycle
              if(( (iside==1.and.idims==1.and.ixComin1<ixCoGmin1+nghostcells) &
                 ) .or.( (iside==2.and.idims==&
                 1.and.ixComax1>ixCoGmax1-nghostcells))) then
                ixBmin1=merge(ixCoGmin1,ixComin1,idims==1);
                ixBmax1=merge(ixCoGmax1,ixComax1,idims==1);
                call bc_phys(iside,idims,time,0.d0,psc(igrid),ixCoGmin1,&
                   ixCoGmax1,ixBmin1,ixBmax1)
              end if
            end do
          end do
        end if

        if(prolongprimitive) then
          ! following line again assumes equidistant grid, but 
          ! just computes indices, so also ok for stretched case
          ! reason for +1-1 and +1+1: the coarse representation has 
          ! also nghostcells at each side. During
          ! prolongation, we need cells to left and right, hence -1/+1
          block=>psc(igrid)
          ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-xComin1)*invdxCo1)+&
             1-1;
          ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-xComin1)*invdxCo1)+&
             1+1;
          call phys_to_primitive(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,&
             psc(igrid)%w,psc(igrid)%x)
        end if

        if(ghost_copy) then
          call interpolation_copy(igrid,ixFimin1,ixFimax1,dxFi1,xFimin1,dxCo1,&
             invdxCo1,xComin1)
        else
          call interpolation_linear(igrid,ixFimin1,ixFimax1,dxFi1,xFimin1,&
             dxCo1,invdxCo1,xComin1)
        end if

        if(prolongprimitive) then
          block=>psc(igrid)
          call phys_to_conserved(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,&
             psc(igrid)%w,psc(igrid)%x)
        end if

      end subroutine bc_prolong

      subroutine bc_prolong_stg(igrid,i1,iib1,NeedProlong)
        use mod_amr_fct
        integer                    :: igrid,i1,iib1
        logical,dimension(-1:1) :: NeedProlong
        logical                    :: fine_min1in,fine_max1in
        integer                    :: ixFimin1,ixFimax1,ixComin1,ixComax1
        double precision           :: dxFi1,dxCo1,xFimin1,xComin1,invdxCo1
        ! Check what is already at the desired level
        fine_min1in=.false.;fine_max1in=.false.;
        
        if(i1>-1) fine_min1in=(.not.NeedProlong(i1-kr(1,&
           1)).and.neighbor_type(i1-kr(1,1),igrid)/=1)
        if(i1<1)  fine_max1in=(.not.NeedProlong(i1+kr(1,&
           1)).and.neighbor_type(i1+kr(1,1),igrid)/=1)
        

        ixFimin1=ixR_srl_min1(iib1,i1);ixFimax1=ixR_srl_max1(iib1,i1);

        dxFi1=rnode(rpdx1_,igrid);
        dxCo1=two*dxFi1;
        invdxCo1=1.d0/dxCo1;

        xFimin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxFi1;
        xComin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxCo1;

        ! moved the physical boundary filling here, to only fill the
        ! part needed

        ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-xComin1)*invdxCo1)+&
           1-1;
        ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-xComin1)*invdxCo1)+1+&
           1;

        if(prolongprimitive) call phys_to_primitive(ixGlo1,ixGhi1,ixFimin1,&
           ixFimax1,psb(igrid)%w,psb(igrid)%x)

        call prolong_2nd_stg(psc(igrid),psb(igrid),ixComin1,ixComax1,ixFimin1,&
           ixFimax1,dxCo1,xComin1,dxFi1,xFimin1,.true.,fine_min1in,&
           fine_max1in)

        if(prolongprimitive) call phys_to_conserved(ixGlo1,ixGhi1,ixFimin1,&
           ixFimax1,psb(igrid)%w,psb(igrid)%x)

        ! The current region has already been refined, so it does not need to be prolonged again
        NeedProlong(i1)=.false. 

      end subroutine bc_prolong_stg

      subroutine interpolation_linear(igrid,ixFimin1,ixFimax1,dxFi1,xFimin1,&
          dxCo1,invdxCo1,xComin1)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: igrid, ixFimin1,ixFimax1
        double precision, intent(in) :: dxFi1, xFimin1,dxCo1, invdxCo1,&
            xComin1

        integer :: ixCo1, jxCo1, hxCo1, ixFi1, ix1, iw, idims, nwmin,nwmax
        double precision :: xCo1, xFi1, eta1
        double precision :: slopeL, slopeR, slopeC, signC, signR
        double precision :: slope(1:nw,ndim)
        !!double precision :: local_invdxCo^D
        double precision :: signedfactorhalf1
        !integer :: ixshift^D, icase

        !icase=mod(nghostcells,2)

        if(prolongprimitive) then
          nwmin=1
          nwmax=nw
        else
          nwmin=nwhead
          nwmax=nwtail
        end if

        do ixFi1 = ixFimin1,ixFimax1
           ! cell-centered coordinates of fine grid point
           ! here we temporarily use an equidistant grid
           xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

           ! indices of coarse cell which contains the fine cell
           ! since we computed lower left corner earlier 
           ! in equidistant fashion: also ok for stretched case
           ixCo1=int((xFi1-xComin1)*invdxCo1)+1

           ! cell-centered coordinates of coarse grid point
           ! here we temporarily use an equidistant grid
           xCo1=xComin1+(dble(ixCo1)-half)*dxCo1 

           !if(.not.slab) then
           !   ^D&local_invdxCo^D=1.d0/psc(igrid)%dx({ixCo^DD},^D)\
           !endif

           if(slab_uniform) then
             ! actual cell-centered coordinates of fine grid point
             !!^D&xFi^D=block%x({ixFi^DD},^D)\
             ! actual cell-centered coordinates of coarse grid point
             !!^D&xCo^D=psc(igrid)%x({ixCo^DD},^D)\
             ! normalized distance between fine/coarse cell center
             ! in coarse cell: ranges from -0.5 to 0.5 in each direction
             ! (origin is coarse cell center)
             ! this is essentially +1/4 or -1/4 on cartesian mesh
             eta1=(xFi1-xCo1)*invdxCo1;
           else
             !select case(icase)
             ! case(0)
             !{! here we assume an even number of ghostcells!!!
             !ixshift^D=2*(mod(ixFi^D,2)-1)+1
             !if(ixshift^D>0.0d0)then
             !   ! oneven fine grid points
             !   eta^D=-0.5d0*(one-block%dvolume(ixFi^DD) &
             !     /sum(block%dvolume(ixFi^D:ixFi^D+1^D%ixFi^DD))) 
             !else
             !   ! even fine grid points
             !   eta^D=+0.5d0*(one-block%dvolume(ixFi^DD) &
             !     /sum(block%dvolume(ixFi^D-1:ixFi^D^D%ixFi^DD))) 
             !endif\}
             ! case(1)
             !{! here we assume an odd number of ghostcells!!!
             !ixshift^D=2*(mod(ixFi^D,2)-1)+1
             !if(ixshift^D>0.0d0)then
             !   ! oneven fine grid points
             !   eta^D=+0.5d0*(one-block%dvolume(ixFi^DD) &
             !     /sum(block%dvolume(ixFi^D-1:ixFi^D^D%ixFi^DD))) 
             !else
             !   ! even fine grid points
             !   eta^D=-0.5d0*(one-block%dvolume(ixFi^DD) &
             !     /sum(block%dvolume(ixFi^D:ixFi^D+1^D%ixFi^DD))) 
             !endif\}
             ! case default
             !  call mpistop("no such case")
             !end select
             ! the different cases for even/uneven number of ghost cells 
             ! are automatically handled using the relative index to ixMlo
             ! as well as the pseudo-coordinates xFi and xCo 
             ! these latter differ from actual cell centers when stretching is used
             ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1;
             if(xFi1>xCo1) then
                signedfactorhalf1=0.5d0
              else
                signedfactorhalf1=-0.5d0
              end if
              eta1=signedfactorhalf1*(one-psb(igrid)%dvolume(ixFi1) &
                 /sum(psb(igrid)%dvolume(ix1:ix1+1))) 
             !{eta^D=(xFi^D-xCo^D)*invdxCo^D &
             !      *two*(one-block%dvolume(ixFi^DD) &
             !      /sum(block%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
           end if

           do idims=1,ndim
              hxCo1=ixCo1-kr(1,idims)
              jxCo1=ixCo1+kr(1,idims)

              do iw=nwmin,nwmax
                 slopeL=psc(igrid)%w(ixCo1,iw)-psc(igrid)%w(hxCo1,iw)
                 slopeR=psc(igrid)%w(jxCo1,iw)-psc(igrid)%w(ixCo1,iw)
                 slopeC=half*(slopeR+slopeL)

                 ! get limited slope
                 signR=sign(one,slopeR)
                 signC=sign(one,slopeC)
                 !select case(prolong_limiter)
                 !case(1)
                 !  ! unlimit
                 !  slope(iw,idims)=slopeC
                 !case(2)
                 !  ! minmod
                 !  slope(iw,idims)=signR*max(zero,min(dabs(slopeR), &
                 !                                    signR*slopeL))
                 !case(3)
                 !  ! woodward
                 !  slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR), &
                 !                     signR*slopeL,signR*half*slopeC))
                 !case(4)
                 !  ! koren
                 !  slope(iw,idims)=signR*max(zero,min(two*signR*slopeL, &
                 !   (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
                 !case default
                   slope(iw,idims)=signC*max(zero,min(dabs(slopeC),&
                       signC*slopeL,signC*slopeR))
                 !end select
              end do
           end do

           ! Interpolate from coarse cell using limited slopes
           psb(igrid)%w(ixFi1,nwmin:nwmax)=psc(igrid)%w(ixCo1,&
              nwmin:nwmax)+(slope(nwmin:nwmax,1)*eta1)

        end do

        if(prolongprimitive) then
          block=>psb(igrid)
          call phys_to_conserved(ixGlo1,ixGhi1,ixFimin1,ixFimax1,psb(igrid)%w,&
             psb(igrid)%x)
        end if

      end subroutine interpolation_linear

      subroutine interpolation_copy(igrid, ixFimin1,ixFimax1,dxFi1,xFimin1,&
          dxCo1,invdxCo1,xComin1)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: igrid, ixFimin1,ixFimax1
        double precision, intent(in) :: dxFi1, xFimin1,dxCo1, invdxCo1,&
            xComin1

        integer :: ixCo1, ixFi1, nwmin,nwmax
        double precision :: xFi1

        if(prolongprimitive) then
          nwmin=1
          nwmax=nw
        else
          nwmin=nwhead
          nwmax=nwtail
        end if

        do ixFi1 = ixFimin1,ixFimax1
           ! cell-centered coordinates of fine grid point
           xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1
        
           ! indices of coarse cell which contains the fine cell
           ! note: this also works for stretched grids
           ixCo1=int((xFi1-xComin1)*invdxCo1)+1
        
           ! Copy from coarse cell
           psb(igrid)%w(ixFi1,nwmin:nwmax)=psc(igrid)%w(ixCo1,nwmin:nwmax)
        
        end do
        
        if(prolongprimitive) call phys_to_conserved(ixGlo1,ixGhi1,ixFimin1,&
           ixFimax1,psb(igrid)%w,psb(igrid)%x)
      
      end subroutine interpolation_copy

      subroutine pole_copy(wrecv,ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,wsend,&
         ixISmin1,ixISmax1,ixSmin1,ixSmax1,ipole)
      
        integer, intent(in) :: ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,ixISmin1,&
           ixISmax1,ixSmin1,ixSmax1,ipole
        double precision :: wrecv(ixIRmin1:ixIRmax1,1:nw),&
            wsend(ixISmin1:ixISmax1,1:nw)

        integer :: iw, iside, iB

        select case (ipole)
        case (1)
           iside=int((i1+3)/2)
           iB=2*(1-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case (bc_symm)
               wrecv(ixRmin1:ixRmax1,iw) = wsend(ixSmax1:ixSmin1:-1,iw)
             case (bc_asymm)
               wrecv(ixRmin1:ixRmax1,iw) =-wsend(ixSmax1:ixSmin1:-1,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do 
        end select
      
      end subroutine pole_copy

      subroutine pole_copy_stg(wrecv,ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,wsend,&
         ixISmin1,ixISmax1,ixSmin1,ixSmax1,idirs,ipole)
      
        integer, intent(in) :: ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,ixISmin1,&
           ixISmax1,ixSmin1,ixSmax1,idirs,ipole

        double precision :: wrecv(ixIRmin1:ixIRmax1,1:nws),&
            wsend(ixISmin1:ixISmax1,1:nws)
        integer :: iB, iside

        select case (ipole)
        case (1)
           iside=int((i1+3)/2)
           iB=2*(1-1)+iside
           select case (typeboundary(iw_mag(idirs),iB))
           case (bc_symm)
             wrecv(ixRmin1:ixRmax1,idirs) = wsend(ixSmax1:ixSmin1:-1,idirs)
           case (bc_asymm)
             wrecv(ixRmin1:ixRmax1,idirs) =-wsend(ixSmax1:ixSmin1:-1,idirs)
           case default
             call mpistop("Pole boundary condition should be symm or asymm")
           end select
         
        end select

      end subroutine pole_copy_stg

      subroutine pole_buffer(wrecv,ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,wsend,&
         ixISmin1,ixISmax1,ixSmin1,ixSmax1)
      
        integer, intent(in) :: ixIRmin1,ixIRmax1,ixRmin1,ixRmax1,ixISmin1,&
           ixISmax1,ixSmin1,ixSmax1
        double precision :: wrecv(ixIRmin1:ixIRmax1,nwhead:nwtail),&
            wsend(ixISmin1:ixISmax1,1:nw)

        integer :: iw, iside, iB

        select case (ipole)
        case (1)
           iside=int((i1+3)/2)
           iB=2*(1-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case (bc_symm)
               wrecv(ixRmin1:ixRmax1,iw) = wsend(ixSmax1:ixSmin1:-1,iw)
             case (bc_asymm)
               wrecv(ixRmin1:ixRmax1,iw) =-wsend(ixSmax1:ixSmin1:-1,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do 
        end select
      
      end subroutine pole_buffer

  end subroutine getbc

  subroutine identifyphysbound(s,iib1)
    use mod_global_parameters

    type(state)          :: s
    integer, intent(out) :: iib1

    
    if(s%is_physical_boundary(2*1) .and. s%is_physical_boundary(2*1-1)) then
      iib1=2
    else if(s%is_physical_boundary(2*1-1)) then
      iib1=-1
    else if(s%is_physical_boundary(2*1)) then
      iib1=1
    else
      iib1=0
    end if
    

  end subroutine identifyphysbound

end module mod_ghostcells_update
