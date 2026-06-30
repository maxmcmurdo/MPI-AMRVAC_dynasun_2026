!> Initialize the MPI environment
subroutine comm_start
  use mod_global_parameters

  integer(kind=MPI_ADDRESS_KIND) :: lb
  integer(kind=MPI_ADDRESS_KIND) :: sizes

  ! Initialize MPI
  call MPI_INIT(ierrmpi)

  ! Each process stores its rank, which ranges from 0 to N-1, where N is the
  ! number of processes.
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierrmpi)

  ! Store the number of processes
  call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierrmpi)

  ! Use the default communicator, which contains all the processes
  icomm = MPI_COMM_WORLD

  ! Get size of double/integer
  call MPI_TYPE_GET_EXTENT(MPI_REAL,lb,sizes,ierrmpi)
  if (sizes /= size_real) call mpistop("Incompatible real size")
  call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,sizes,ierrmpi)
  if (sizes /= size_double) call mpistop("Incompatible double size")
  call MPI_TYPE_GET_EXTENT(MPI_INTEGER,lb,sizes,ierrmpi)
  if (sizes /= size_int) call mpistop("Incompatible integer size")
  call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,sizes,ierrmpi)
  if (sizes /= size_logical) call mpistop("Incompatible logical size")

end subroutine comm_start

!> Finalize (or shutdown) the MPI environment
subroutine comm_finalize
  use mod_global_parameters
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
  call MPI_FINALIZE(ierrmpi)

end subroutine comm_finalize

!> Create and store the MPI types that will be used for parallel communication
subroutine init_comm_types
  use mod_global_parameters

  integer, dimension(ndim+1) :: sizes, subsizes, start
  integer :: i1, ic1, nx1, nxCo1, nxG1, idir

  nx1=ixMhi1-ixMlo1+1;
  nxG1=ixGhi1-ixGlo1+1;
  nxCo1=nx1/2;

  sizes(1)=ixGhi1;
  sizes(ndim+1)=nw
  subsizes(1)=nxG1;
  subsizes(ndim+1)=nw
  start(1)=ixGlo1-1;
  start(ndim+1)=0
  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
     MPI_DOUBLE_PRECISION, type_block,ierrmpi)
  call MPI_TYPE_COMMIT(type_block,ierrmpi)
  size_block=nxG1*nw*size_double

  sizes(1)=ixGhi1/2+nghostcells;
  sizes(ndim+1)=nw
  subsizes(1)=nxCo1;
  subsizes(ndim+1)=nw
  start(1)=ixMlo1-1;
  start(ndim+1)=0
  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
     MPI_DOUBLE_PRECISION, type_coarse_block,ierrmpi)
  call MPI_TYPE_COMMIT(type_coarse_block,ierrmpi)

  if(stagger_grid) then
    sizes(1)=ixGhi1+1;
    sizes(ndim+1)=nws
    subsizes(1)=nx1+1;
    subsizes(ndim+1)=nws
    start(1)=ixMlo1-1;
    start(ndim+1)=0
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
        MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_block_io_stg,ierrmpi)
    call MPI_TYPE_COMMIT(type_block_io_stg,ierrmpi)
    size_block_io_stg=(nx1+1)*nws*size_double

    sizes(1)=ixGhi1/2+nghostcells+1;
    sizes(ndim+1)=nws
   do ic1=1,2
      do idir=1,ndim
        subsizes(1)=nxCo1+kr(ic1,1)*kr(idir,1);
        subsizes(ndim+1)=1
        start(1)=ixMlo1-kr(ic1,1)*kr(idir,1);
        start(ndim+1)=idir-1

        call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_coarse_block_stg(idir,&
           ic1),ierrmpi)
        call MPI_TYPE_COMMIT(type_coarse_block_stg(idir,ic1),ierrmpi)
      end do
   end do

    sizes(1)=ixGhi1+1;
    sizes(ndim+1)=nws
   do ic1=1,2
      do idir=1,1
        subsizes(1)=nxCo1+kr(ic1,1)*kr(idir,1);
        subsizes(ndim+1)=1
        start(1)=ixMlo1-kr(ic1,1)*kr(idir,1)+(ic1-1)*nxCo1;
        start(ndim+1)=idir-1
        call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_sub_block_stg(idir,&
           ic1),ierrmpi)
        call MPI_TYPE_COMMIT(type_sub_block_stg(idir,ic1),ierrmpi)
      end do
   end do
  end if

  sizes(1)=ixGhi1;
  sizes(ndim+1)=nw
  do ic1=1,2
     subsizes(1)=nxCo1;
     subsizes(ndim+1)=nw
     start(1)=ixMlo1-1+(ic1-1)*nxCo1;
     start(ndim+1)=0
     call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION, type_sub_block(ic1),ierrmpi)
     call MPI_TYPE_COMMIT(type_sub_block(ic1),ierrmpi)
  end do

  sizes(1)=ixGhi1;
  sizes(ndim+1)=nw
  subsizes(1)=nx1;
  subsizes(ndim+1)=nw
  start(1)=ixMlo1-1;
  start(ndim+1)=0
  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
     MPI_DOUBLE_PRECISION, type_block_io,ierrmpi)
  call MPI_TYPE_COMMIT(type_block_io,ierrmpi)
  size_block_io=nx1*nw*size_double

  sizes(1)=ixMhi1-ixMlo1+1;
  sizes(ndim+1)=1
  subsizes(1)=sizes(1);
  subsizes(ndim+1)=1
  start(1)=0;
  start(ndim+1)=0
  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
     MPI_DOUBLE_PRECISION, type_block_xcc_io,ierrmpi)
  call MPI_TYPE_COMMIT(type_block_xcc_io,ierrmpi)

  sizes(1)=ixMhi1-ixMlo1+2;
  sizes(ndim+1)=1
  subsizes(1)=sizes(1);
  subsizes(ndim+1)=1
  start(1)=0;
  start(ndim+1)=0
  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
     MPI_DOUBLE_PRECISION, type_block_xc_io,ierrmpi)
  call MPI_TYPE_COMMIT(type_block_xc_io,ierrmpi)

  sizes(1)=ixMhi1-ixMlo1+1;
  sizes(ndim+1)=nw+nwauxio
  subsizes(1)=sizes(1);
  subsizes(ndim+1)=nw+nwauxio
  start(1)=0;
  start(ndim+1)=0
  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
     MPI_DOUBLE_PRECISION, type_block_wcc_io,ierrmpi)
  call MPI_TYPE_COMMIT(type_block_wcc_io,ierrmpi)

  sizes(1)=ixMhi1-ixMlo1+2;
  sizes(ndim+1)=nw+nwauxio
  subsizes(1)=sizes(1);
  subsizes(ndim+1)=nw+nwauxio
  start(1)=0;
  start(ndim+1)=0
  call MPI_TYPE_CREATE_SUBARRAY(ndim+1,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
     MPI_DOUBLE_PRECISION, type_block_wc_io,ierrmpi)
  call MPI_TYPE_COMMIT(type_block_wc_io,ierrmpi)

end subroutine init_comm_types

!> Exit MPI-AMRVAC with an error message
subroutine mpistop(message)
  use mod_global_parameters

  character(len=*), intent(in) :: message !< The error message
  integer                      :: ierrcode

  write(*, *) "ERROR for processor", mype, ":"
  write(*, *) trim(message)

  call MPI_ABORT(icomm, ierrcode, ierrmpi)

end subroutine mpistop
