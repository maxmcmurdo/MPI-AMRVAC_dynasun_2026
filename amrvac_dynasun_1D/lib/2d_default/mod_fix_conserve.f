!> Module for flux conservation near refinement boundaries
module mod_fix_conserve
  implicit none
  private

  type fluxalloc
     double precision, dimension(:,:,:), pointer:: flux => null()
     double precision, dimension(:,:,:), pointer:: edge => null()
  end type fluxalloc
  !> store flux to fix conservation
  type(fluxalloc), dimension(:,:,:), allocatable, public :: pflux

  integer, save                        :: nrecv, nsend
  double precision, allocatable, save  :: recvbuffer(:), sendbuffer(:)
  integer, dimension(:), allocatable   :: fc_recvreq, fc_sendreq
  integer, dimension(:,:), allocatable :: fc_recvstat, fc_sendstat
  integer, dimension(2), save        :: isize
  integer                              :: ibuf, ibuf_send
  ! ct for corner total
  integer, save                        :: nrecv_ct, nsend_ct
  ! buffer for corner coarse
  double precision, allocatable, save  :: recvbuffer_cc(:), sendbuffer_cc(:)
  integer, dimension(:), allocatable   :: cc_recvreq, cc_sendreq
  integer, dimension(:,:), allocatable :: cc_recvstat, cc_sendstat
  integer, dimension(2), save        :: isize_stg
  integer                              :: ibuf_cc, ibuf_cc_send
  integer                              :: itag, itag_cc, isend, isend_cc,&
      irecv, irecv_cc

  public :: init_comm_fix_conserve
  public :: allocateBflux
  public :: deallocateBflux
  public :: sendflux
  public :: recvflux
  public :: store_flux
  public :: store_edge
  public :: fix_conserve
  public :: fix_edges

 contains

   subroutine init_comm_fix_conserve(idimmin,idimmax,nwfluxin)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax,nwfluxin

     integer :: iigrid, igrid, idims, iside, i1,i2, nxCo1,nxCo2
     integer :: ic1,ic2, inc1,inc2, ipe_neighbor
     integer :: recvsize, sendsize
     integer :: recvsize_cc, sendsize_cc

     nsend    = 0
     nrecv    = 0
     recvsize = 0
     sendsize = 0
     if(stagger_grid) then
       ! Special communication for diagonal 'coarse corners'
       ! nrecv/send_cc (for 'coarse corners' is a dim=ndim-1 array which
       ! stores the faces that must be communicated in each direction.
       ! nrecv/send_ct (for 'corners total' is the total number of
       ! necessary communications. These special cases have their own
       ! send and receive buffers (send/recvbuffer_cc), their tags, etc.
       nsend_ct=0
       nrecv_ct=0
       recvsize_cc=0
       sendsize_cc=0
     end if

     do idims= idimmin,idimmax
       select case (idims)
         case (1)
         nrecv=nrecv+nrecv_fc(1)
         nsend=nsend+nsend_fc(1)
         nxCo1=1;nxCo2=ixGhi2/2-nghostcells;
         isize(1)=nxCo1*nxCo2*(nwfluxin)
         recvsize=recvsize+nrecv_fc(1)*isize(1)
         sendsize=sendsize+nsend_fc(1)*isize(1)
         if(stagger_grid) then
           ! This does not consider the 'coarse corner' case
           nxCo1=1;nxCo2=ixGhi2/2-nghostcells+1;
           isize_stg(1)=nxCo1*nxCo2*(2-1)
           ! the whole size is used (cell centered and staggered)
           isize(1)=isize(1)+isize_stg(1)      
           recvsize=recvsize+nrecv_fc(1)*isize_stg(1)
           sendsize=sendsize+nsend_fc(1)*isize_stg(1)
           ! Coarse corner case
           nrecv_ct=nrecv_ct+nrecv_cc(1)
           nsend_ct=nsend_ct+nsend_cc(1)
           recvsize_cc=recvsize_cc+nrecv_cc(1)*isize_stg(1)
           sendsize_cc=sendsize_cc+nsend_cc(1)*isize_stg(1)
         end if
         
         case (2)
         nrecv=nrecv+nrecv_fc(2)
         nsend=nsend+nsend_fc(2)
         nxCo1=ixGhi1/2-nghostcells;nxCo2=1;
         isize(2)=nxCo1*nxCo2*(nwfluxin)
         recvsize=recvsize+nrecv_fc(2)*isize(2)
         sendsize=sendsize+nsend_fc(2)*isize(2)
         if(stagger_grid) then
           ! This does not consider the 'coarse corner' case
           nxCo1=ixGhi1/2-nghostcells+1;nxCo2=1;
           isize_stg(2)=nxCo1*nxCo2*(2-1)
           ! the whole size is used (cell centered and staggered)
           isize(2)=isize(2)+isize_stg(2)      
           recvsize=recvsize+nrecv_fc(2)*isize_stg(2)
           sendsize=sendsize+nsend_fc(2)*isize_stg(2)
           ! Coarse corner case
           nrecv_ct=nrecv_ct+nrecv_cc(2)
           nsend_ct=nsend_ct+nsend_cc(2)
           recvsize_cc=recvsize_cc+nrecv_cc(2)*isize_stg(2)
           sendsize_cc=sendsize_cc+nsend_cc(2)*isize_stg(2)
         end if
         
       end select
     end do

     ! Reallocate buffers when size differs
     if (allocated(recvbuffer)) then
       if (recvsize /= size(recvbuffer)) then
         deallocate(recvbuffer)
         allocate(recvbuffer(recvsize))
       end if
     else
       allocate(recvbuffer(recvsize))
     end if

     if (allocated(fc_recvreq)) then
       if (nrecv /= size(fc_recvreq)) then
         deallocate(fc_recvreq, fc_recvstat)
         allocate(fc_recvstat(MPI_STATUS_SIZE,nrecv), fc_recvreq(nrecv))
       end if
     else
       allocate(fc_recvstat(MPI_STATUS_SIZE,nrecv), fc_recvreq(nrecv))
     end if

     if (allocated(fc_sendreq)) then
       if (nsend /= size(fc_sendreq)) then
         deallocate(fc_sendreq, fc_sendstat)
         allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
       end if
     else
       allocate(fc_sendstat(MPI_STATUS_SIZE,nsend), fc_sendreq(nsend))
     end if

     if(stagger_grid) then

       if (allocated(sendbuffer)) then
         if (sendsize /= size(sendbuffer)) then
           deallocate(sendbuffer)
           allocate(sendbuffer(sendsize))
         end if
       else
         allocate(sendbuffer(sendsize))
       end if

       if (allocated(recvbuffer_cc)) then
         if (recvsize_cc /= size(recvbuffer_cc)) then
           deallocate(recvbuffer_cc)
           allocate(recvbuffer_cc(recvsize_cc))
         end if
       else
         allocate(recvbuffer_cc(recvsize_cc))
       end if

       if (allocated(cc_recvreq)) then
         if (nrecv_ct /= size(cc_recvreq)) then
           deallocate(cc_recvreq, cc_recvstat)
           allocate(cc_recvstat(MPI_STATUS_SIZE,nrecv_ct),&
               cc_recvreq(nrecv_ct))
         end if
       else
         allocate(cc_recvstat(MPI_STATUS_SIZE,nrecv_ct), cc_recvreq(nrecv_ct))
       end if

       if (allocated(sendbuffer_cc)) then
         if (sendsize_cc /= size(sendbuffer_cc)) then
           deallocate(sendbuffer_cc)
           allocate(sendbuffer_cc(sendsize_cc))
         end if
       else
         allocate(sendbuffer_cc(sendsize_cc))
       end if

       if (allocated(cc_sendreq)) then
         if (nsend_ct /= size(cc_sendreq)) then
           deallocate(cc_sendreq, cc_sendstat)
           allocate(cc_sendstat(MPI_STATUS_SIZE,nsend_ct),&
               cc_sendreq(nsend_ct))
         end if
       else
         allocate(cc_sendstat(MPI_STATUS_SIZE,nsend_ct), cc_sendreq(nsend_ct))
       end if
     end if

   end subroutine init_comm_fix_conserve

   subroutine recvflux(idimmin,idimmax)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax

     integer :: iigrid, igrid, idims, iside, i1,i2, nxCo1,nxCo2
     integer :: ic1,ic2, inc1,inc2, ipe_neighbor
     integer :: pi1,pi2,mi1,mi2,ph1,ph2,mh1,mh2,idir

     if (nrecv>0) then
       fc_recvreq=MPI_REQUEST_NULL
       ibuf=1
       irecv=0

       do iigrid=1,igridstail; igrid=igrids(iigrid);
         do idims= idimmin,idimmax
           do iside=1,2
             i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3);

             if (neighbor_pole(i1,i2,igrid)/=0) cycle

             if (neighbor_type(i1,i2,igrid)/=4) cycle
             do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
             do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
               if (ipe_neighbor/=mype) then
                 irecv=irecv+1
                 itag=4**2*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)
                 call MPI_IRECV(recvbuffer(ibuf),isize(idims),&
                     MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                    fc_recvreq(irecv),ierrmpi)
                 ibuf=ibuf+isize(idims)
               end if
             end do
             end do
           end do
         end do
       end do
     end if

     if(stagger_grid) then
     ! receive corners
       if (nrecv_ct>0) then
         cc_recvreq=MPI_REQUEST_NULL
         ibuf_cc=1
         irecv_cc=0

         do iigrid=1,igridstail; igrid=igrids(iigrid);
           do idims= idimmin,idimmax
             do iside=1,2
               i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3);
               ! Check if there are special corners
               ! (Coarse block diagonal to a fine block)
               ! If there are, receive.
               ! Tags are calculated in the same way as for
               ! normal fluxes, but should not overlap because
               ! inc^D are different
               if (neighbor_type(i1,i2,igrid)==3) then
                 do idir=idims+1,ndim
                   pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);
                   mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);
                   ph1=pi1-kr(idims,1)*(2*iside-3)
                   ph2=pi2-kr(idims,2)*(2*iside-3);
                   mh1=mi1-kr(idims,1)*(2*iside-3)
                   mh2=mi2-kr(idims,2)*(2*iside-3);

                   if (neighbor_type(pi1,pi2,igrid)==4.and.neighbor_type(ph1,&
                      ph2,igrid)==3.and.neighbor_pole(pi1,pi2,igrid)==0) then
                      ! Loop on children (several in 3D)
                    do ic2=1+int((1-pi2)/2),2-int((1+pi2)/2)
                       inc2=2*pi2+ic2
                    do ic1=1+int((1-pi1)/2),2-int((1+pi1)/2)
                       inc1=2*pi1+ic1
                       ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
                       if (mype/=ipe_neighbor) then
                         irecv_cc=irecv_cc+1
                         itag_cc=4**2*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)
                         call MPI_IRECV(recvbuffer_cc(ibuf_cc),&
                            isize_stg(idims),MPI_DOUBLE_PRECISION,ipe_neighbor,&
                            itag_cc,icomm,cc_recvreq(irecv_cc),ierrmpi)
                         ibuf_cc=ibuf_cc+isize_stg(idims)
                       end if
                    end do
                    end do
                   end if

                   if (neighbor_type(mi1,mi2,igrid)==4.and.neighbor_type(mh1,&
                      mh2,igrid)==3.and.neighbor_pole(mi1,mi2,igrid)==0) then
                      ! Loop on children (several in 3D)
                    do ic2=1+int((1-mi2)/2),2-int((1+mi2)/2)
                        inc2=2*mi2+ic2
                    do ic1=1+int((1-mi1)/2),2-int((1+mi1)/2)
                        inc1=2*mi1+ic1
                       ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
                       if (mype/=ipe_neighbor) then
                         irecv_cc=irecv_cc+1
                         itag_cc=4**2*(igrid-1)+inc1*4**(1-1)+inc2*4**(2-1)
                         call MPI_IRECV(recvbuffer_cc(ibuf_cc),&
                            isize_stg(idims),MPI_DOUBLE_PRECISION,ipe_neighbor,&
                            itag_cc,icomm,cc_recvreq(irecv_cc),ierrmpi)
                         ibuf_cc=ibuf_cc+isize_stg(idims)
                       end if
                    end do
                    end do
                   end if
                 end do
               end if
             end do
           end do
         end do
       end if
     end if ! end if stagger grid

   end subroutine recvflux

   subroutine sendflux(idimmin,idimmax)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax

     integer :: idims, iside, i1,i2, ic1,ic2, inc1,inc2, ix1,ix2, ixCo1,ixCo2,&
         nxCo1,nxCo2, iw
     integer :: ineighbor, ipe_neighbor, igrid, iigrid, ibuf_send_next
     integer :: idir, ibuf_cc_send_next, pi1,pi2, ph1,ph2, mi1,mi2, mh1,mh2

     fc_sendreq = MPI_REQUEST_NULL
     isend      = 0
     if(stagger_grid) then
       ibuf_send  = 1
       cc_sendreq=MPI_REQUEST_NULL
       isend_cc=0
       ibuf_cc_send=1
     end if

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims = idimmin,idimmax
         select case (idims)
        case (1)
           do iside=1,2
             i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);

             if (neighbor_pole(i1,i2,igrid)/=0) cycle

             if (neighbor_type(i1,i2,igrid)==neighbor_coarse) then
               ! send flux to coarser neighbor
               ineighbor=neighbor(1,i1,i2,igrid)
               ipe_neighbor=neighbor(2,i1,i2,igrid)
               if (ipe_neighbor/=mype) then
                 ic1=1+modulo(node(pig1_,igrid)-1,2)
                 ic2=1+modulo(node(pig2_,igrid)-1,2);
                 inc1=-2*i1+ic1;inc2=ic2;
                 itag=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                 isend=isend+1

                 if(stagger_grid) then
                   ibuf_send_next=ibuf_send+isize(1)
                   sendbuffer(ibuf_send:ibuf_send_next-isize_stg(1)-&
                      1)=reshape(pflux(iside,1,igrid)%flux,&
                      (/isize(1)-isize_stg(1)/))

                   sendbuffer(ibuf_send_next-isize_stg(1):ibuf_send_next-&
                      1)=reshape(pflux(iside,1,igrid)%edge,(/isize_stg(1)/))
                   call MPI_ISEND(sendbuffer(ibuf_send),isize(1),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                   ibuf_send=ibuf_send_next
                 else
                   call MPI_ISEND(pflux(iside,1,igrid)%flux,isize(1),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                 end if
               end if

               if(stagger_grid) then
                 ! If we are in a fine block surrounded by coarse blocks
                 do idir=idims+1,ndim
                   pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);
                   mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);
                   ph1=pi1-kr(idims,1)*(2*iside-3)
                   ph2=pi2-kr(idims,2)*(2*iside-3);
                   mh1=mi1-kr(idims,1)*(2*iside-3)
                   mh2=mi2-kr(idims,2)*(2*iside-3);

                   if (neighbor_type(pi1,pi2,igrid)==2.and.neighbor_type(ph1,&
                      ph2,igrid)==2.and.mype/=neighbor(2,pi1,pi2,&
                      igrid).and.neighbor_pole(pi1,pi2,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,pi1,pi2,igrid)
                     ipe_neighbor=neighbor(2,pi1,pi2,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;
                     itag_cc=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(1)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,1,igrid)%edge,&
                        shape=(/isize_stg(1)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(1),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
       
                   if (neighbor_type(mi1,mi2,igrid)==2.and.neighbor_type(mh1,&
                      mh2,igrid)==2.and.mype/=neighbor(2,mi1,mi2,&
                      igrid).and.neighbor_pole(mi1,mi2,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,mi1,mi2,igrid)
                     ipe_neighbor=neighbor(2,mi1,mi2,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;
                     inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;
                     itag_cc=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(1)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,1,igrid)%edge,&
                        shape=(/isize_stg(1)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(1),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
                 end do
               end if ! end if stagger grid

             end if
           end do
        case (2)
           do iside=1,2
             i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);

             if (neighbor_pole(i1,i2,igrid)/=0) cycle

             if (neighbor_type(i1,i2,igrid)==neighbor_coarse) then
               ! send flux to coarser neighbor
               ineighbor=neighbor(1,i1,i2,igrid)
               ipe_neighbor=neighbor(2,i1,i2,igrid)
               if (ipe_neighbor/=mype) then
                 ic1=1+modulo(node(pig1_,igrid)-1,2)
                 ic2=1+modulo(node(pig2_,igrid)-1,2);
                 inc1=ic1;inc2=-2*i2+ic2;
                 itag=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                 isend=isend+1

                 if(stagger_grid) then
                   ibuf_send_next=ibuf_send+isize(2)
                   sendbuffer(ibuf_send:ibuf_send_next-isize_stg(2)-&
                      1)=reshape(pflux(iside,2,igrid)%flux,&
                      (/isize(2)-isize_stg(2)/))

                   sendbuffer(ibuf_send_next-isize_stg(2):ibuf_send_next-&
                      1)=reshape(pflux(iside,2,igrid)%edge,(/isize_stg(2)/))
                   call MPI_ISEND(sendbuffer(ibuf_send),isize(2),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                   ibuf_send=ibuf_send_next
                 else
                   call MPI_ISEND(pflux(iside,2,igrid)%flux,isize(2),&
                       MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                      fc_sendreq(isend),ierrmpi)
                 end if
               end if

               if(stagger_grid) then
                 ! If we are in a fine block surrounded by coarse blocks
                 do idir=idims+1,ndim
                   pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);
                   mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);
                   ph1=pi1-kr(idims,1)*(2*iside-3)
                   ph2=pi2-kr(idims,2)*(2*iside-3);
                   mh1=mi1-kr(idims,1)*(2*iside-3)
                   mh2=mi2-kr(idims,2)*(2*iside-3);

                   if (neighbor_type(pi1,pi2,igrid)==2.and.neighbor_type(ph1,&
                      ph2,igrid)==2.and.mype/=neighbor(2,pi1,pi2,&
                      igrid).and.neighbor_pole(pi1,pi2,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,pi1,pi2,igrid)
                     ipe_neighbor=neighbor(2,pi1,pi2,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;
                     itag_cc=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(2)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,2,igrid)%edge,&
                        shape=(/isize_stg(2)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(2),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
       
                   if (neighbor_type(mi1,mi2,igrid)==2.and.neighbor_type(mh1,&
                      mh2,igrid)==2.and.mype/=neighbor(2,mi1,mi2,&
                      igrid).and.neighbor_pole(mi1,mi2,igrid)==0) then
                     ! Get relative position in the grid for tags
                     ineighbor=neighbor(1,mi1,mi2,igrid)
                     ipe_neighbor=neighbor(2,mi1,mi2,igrid)
                     ic1=1+modulo(node(pig1_,igrid)-1,2)
                     ic2=1+modulo(node(pig2_,igrid)-1,2);
                     inc1=-2*pi1+ic1;inc2=-2*pi2+ic2;
                     inc1=-2*mi1+ic1;inc2=-2*mi2+ic2;
                     itag_cc=4**2*(ineighbor-1)+inc1*4**(1-1)+inc2*4**(2-1)
                     ! Reshape to buffer and send
                     isend_cc=isend_cc+1
                     ibuf_cc_send_next=ibuf_cc_send+isize_stg(2)
                     sendbuffer_cc(ibuf_cc_send:ibuf_cc_send_next-&
                        1)=reshape(pflux(iside,2,igrid)%edge,&
                        shape=(/isize_stg(2)/))
                     call MPI_ISEND(sendbuffer_cc(ibuf_cc_send),isize_stg(2),&
                        MPI_DOUBLE_PRECISION,ipe_neighbor,itag_cc,icomm,&
                        cc_sendreq(isend_cc),ierrmpi)
                     ibuf_cc_send=ibuf_cc_send_next
                   end if
                 end do
               end if ! end if stagger grid

             end if
           end do
         end select
       end do
     end do
   end subroutine sendflux

   subroutine allocateBflux
     use mod_global_parameters

     integer :: iigrid, igrid, iside, i1,i2, nx1,nx2, nxCo1,nxCo2
     integer :: idir,idim,pi1,pi2, mi1,mi2, ph1,ph2, mh1,mh2 !To detect corners

     nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
     nxCo1=nx1/2;nxCo2=nx2/2;

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! For every grid,
       ! arrays for the fluxes are allocated for every face direction(^D)
       ! and every side (1=left, 2=right)
       do iside=1,2
         i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);

         if (neighbor_pole(i1,i2,igrid)/=0) cycle

         select case (neighbor_type(i1,i2,igrid))
         case(neighbor_fine)
           allocate(pflux(iside,1,igrid)%flux(1,1:nx2,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,1,igrid)%edge(1,0:nx2,&
              1:ndim-1))
         case(neighbor_coarse)
           allocate(pflux(iside,1,igrid)%flux(1,1:nxCo2,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,1,igrid)%edge(1,0:nxCo2,&
              1:ndim-1))
         case(neighbor_sibling)
           if(stagger_grid) then
             idim=1
             do idir=idim+1,ndim
             !do idir=min(idim+1,ndim),ndim
               pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);
               mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);
               ph1=pi1-kr(1,1)*(2*iside-3);ph2=pi2-kr(1,2)*(2*iside-3);
               mh1=mi1-kr(1,1)*(2*iside-3);mh2=mi2-kr(1,2)*(2*iside-3);
               if ((neighbor_type(pi1,pi2,igrid)==4.and.neighbor_type(ph1,ph2,&
                  igrid)==3).or.(neighbor_type(mi1,mi2,&
                  igrid)==4.and.neighbor_type(mh1,mh2,igrid)==3)) then
                 allocate(pflux(iside,1,igrid)%edge(1,0:nx2,1:ndim-1))
                 exit
               end if
             end do
           end if
         end select
       end do
       do iside=1,2
         i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);

         if (neighbor_pole(i1,i2,igrid)/=0) cycle

         select case (neighbor_type(i1,i2,igrid))
         case(neighbor_fine)
           allocate(pflux(iside,2,igrid)%flux(1:nx1,1,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,2,igrid)%edge(0:nx1,1,&
              1:ndim-1))
         case(neighbor_coarse)
           allocate(pflux(iside,2,igrid)%flux(1:nxCo1,1,1:nwflux))
           if(stagger_grid) allocate(pflux(iside,2,igrid)%edge(0:nxCo1,1,&
              1:ndim-1))
         case(neighbor_sibling)
           if(stagger_grid) then
             idim=2
             do idir=idim+1,ndim
             !do idir=min(idim+1,ndim),ndim
               pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);
               mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);
               ph1=pi1-kr(2,1)*(2*iside-3);ph2=pi2-kr(2,2)*(2*iside-3);
               mh1=mi1-kr(2,1)*(2*iside-3);mh2=mi2-kr(2,2)*(2*iside-3);
               if ((neighbor_type(pi1,pi2,igrid)==4.and.neighbor_type(ph1,ph2,&
                  igrid)==3).or.(neighbor_type(mi1,mi2,&
                  igrid)==4.and.neighbor_type(mh1,mh2,igrid)==3)) then
                 allocate(pflux(iside,2,igrid)%edge(0:nx1,1,1:ndim-1))
                 exit
               end if
             end do
           end if
         end select
       end do
     end do

   end subroutine allocateBflux

   subroutine deallocateBflux
     use mod_global_parameters

     integer :: iigrid, igrid, iside

     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do iside=1,2
         if (associated(pflux(iside,1,igrid)%flux)) then
           deallocate(pflux(iside,1,igrid)%flux)
           nullify(pflux(iside,1,igrid)%flux)
         end if
         if (associated(pflux(iside,1,igrid)%edge)) then
           deallocate(pflux(iside,1,igrid)%edge)
           nullify(pflux(iside,1,igrid)%edge)
         end if
       end do
       do iside=1,2
         if (associated(pflux(iside,2,igrid)%flux)) then
           deallocate(pflux(iside,2,igrid)%flux)
           nullify(pflux(iside,2,igrid)%flux)
         end if
         if (associated(pflux(iside,2,igrid)%edge)) then
           deallocate(pflux(iside,2,igrid)%edge)
           nullify(pflux(iside,2,igrid)%edge)
         end if
       end do
     end do

   end subroutine deallocateBflux

   subroutine fix_conserve(psb,idimmin,idimmax,nw0,nwfluxin)
     use mod_global_parameters

     integer, intent(in) :: idimmin,idimmax, nw0, nwfluxin
     type(state) :: psb(max_blocks)

     integer :: iigrid, igrid, idims, iside, iotherside, i1,i2, ic1,ic2, inc1,&
        inc2, ixmin1,ixmin2,ixmax1,ixmax2
     integer :: nxCo1,nxCo2, iw, ix, ipe_neighbor, ineighbor, nbuf, ibufnext,&
         nw1
     double precision :: CoFiratio

     nw1=nw0-1+nwfluxin
     if (slab_uniform) then
       ! The flux is divided by volume of fine cell. We need, however,
       ! to divide by volume of coarse cell => muliply by volume ratio
       CoFiratio=one/dble(2**ndim)
     end if

     if (nrecv>0) then
       call MPI_WAITALL(nrecv,fc_recvreq,fc_recvstat,ierrmpi)
       ibuf=1
     end if

     nxCo1=(ixMhi1-ixMlo1+1)/2;nxCo2=(ixMhi2-ixMlo2+1)/2;

     ! for all grids: perform flux update at Coarse-Fine interfaces
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims= idimmin,idimmax
         select case (idims)
           case (1)
           do iside=1,2
             i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);

             if (neighbor_pole(i1,i2,igrid)/=0) cycle

             if (neighbor_type(i1,i2,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is 
             ! not computed, keep the buffer counter up to date:
             if (.not.neighbor_active(i1,i2,igrid).or..not.neighbor_active(0,0,&
                igrid) ) then
               do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
               if (ipe_neighbor/=mype) then
                 ibufnext=ibuf+isize(1)
                 ibuf=ibufnext
               end if
               end do
           end do
               cycle
             end if
             !

             select case (iside)
             case (1)
               ix=ixMlo1
             case (2)
               ix=ixMhi1
             end select

             ! remove coarse flux
             if (slab_uniform) then
               psb(igrid)%w(ix,ixMlo2:ixMhi2,nw0:nw1) = psb(igrid)%w(ix,&
                  ixMlo2:ixMhi2,nw0:nw1) -pflux(iside,1,igrid)%flux(1,:,&
                  1:nwfluxin)
             else
               do iw=nw0,nw1
                 psb(igrid)%w(ix,ixMlo2:ixMhi2,iw)=psb(igrid)%w(ix,&
                    ixMlo2:ixMhi2,iw)-pflux(iside,1,igrid)%flux(1,:,&
                    iw-nw0+1) /ps(igrid)%dvolume(ix,ixMlo2:ixMhi2)
               end do
             end if


             ! add fine flux
            do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
               ixmin1=ix;ixmin2=ixMlo2+(ic2-1)*nxCo2;
               ixmax1=ix;ixmax2=ixmin2-1+nxCo2;
               if (ipe_neighbor==mype) then
                 iotherside=3-iside
                 if (slab_uniform) then
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1) + pflux(iotherside,1,ineighbor)%flux(:,:,&
                      1:nwfluxin)* CoFiratio
                 else
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw) +pflux(iotherside,1,ineighbor)%flux(:,:,&
                        iw-nw0+1) /ps(igrid)%dvolume(ixmin1:ixmax1,&
                        ixmin2:ixmax2)
                   end do
                 end if
               else
                 if (slab_uniform) then
                   ibufnext=ibuf+isize(1)
                   if(stagger_grid) ibufnext=ibufnext-isize_stg(1)
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1)+CoFiratio &
                      *reshape(source=recvbuffer(ibuf:ibufnext-1),&
                       shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1)))
                   ibuf=ibuf+isize(1)
                 else
                   ibufnext=ibuf+isize(1)
                   if(stagger_grid) then
                     nbuf=(isize(1)-isize_stg(1))/nwfluxin
                   else
                     nbuf=isize(1)/nwfluxin
                   end if
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw) +reshape(source=recvbuffer(ibuf:ibufnext-1),&
                         shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw))) /ps(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2)
                     ibuf=ibuf+nbuf
                   end do
                   ibuf=ibufnext
                 end if
               end if
            end do
           end do
           end do
           case (2)
           do iside=1,2
             i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);

             if (neighbor_pole(i1,i2,igrid)/=0) cycle

             if (neighbor_type(i1,i2,igrid)/=4) cycle

 !opedit: skip over active/passive interface since flux for passive ones is 
             ! not computed, keep the buffer counter up to date:
             if (.not.neighbor_active(i1,i2,igrid).or..not.neighbor_active(0,0,&
                igrid) ) then
               do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
               if (ipe_neighbor/=mype) then
                 ibufnext=ibuf+isize(2)
                 ibuf=ibufnext
               end if
               end do
           end do
               cycle
             end if
             !

             select case (iside)
             case (1)
               ix=ixMlo2
             case (2)
               ix=ixMhi2
             end select

             ! remove coarse flux
             if (slab_uniform) then
               psb(igrid)%w(ixMlo1:ixMhi1,ix,&
                  nw0:nw1) = psb(igrid)%w(ixMlo1:ixMhi1,ix,&
                  nw0:nw1) -pflux(iside,2,igrid)%flux(:,1,1:nwfluxin)
             else
               do iw=nw0,nw1
                 psb(igrid)%w(ixMlo1:ixMhi1,ix,iw)=psb(igrid)%w(ixMlo1:ixMhi1,&
                    ix,iw)-pflux(iside,2,igrid)%flux(:,1,&
                    iw-nw0+1) /ps(igrid)%dvolume(ixMlo1:ixMhi1,ix)
               end do
             end if


             ! add fine flux
            do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
           do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
               ixmin1=ixMlo1+(ic1-1)*nxCo1;ixmin2=ix;
               ixmax1=ixmin1-1+nxCo1;ixmax2=ix;
               if (ipe_neighbor==mype) then
                 iotherside=3-iside
                 if (slab_uniform) then
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1) + pflux(iotherside,2,ineighbor)%flux(:,:,&
                      1:nwfluxin)* CoFiratio
                 else
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw) +pflux(iotherside,2,ineighbor)%flux(:,:,&
                        iw-nw0+1) /ps(igrid)%dvolume(ixmin1:ixmax1,&
                        ixmin2:ixmax2)
                   end do
                 end if
               else
                 if (slab_uniform) then
                   ibufnext=ibuf+isize(2)
                   if(stagger_grid) ibufnext=ibufnext-isize_stg(2)
                   psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1) = psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1)+CoFiratio &
                      *reshape(source=recvbuffer(ibuf:ibufnext-1),&
                       shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                      nw0:nw1)))
                   ibuf=ibuf+isize(2)
                 else
                   ibufnext=ibuf+isize(2)
                   if(stagger_grid) then
                     nbuf=(isize(2)-isize_stg(2))/nwfluxin
                   else
                     nbuf=isize(2)/nwfluxin
                   end if
                   do iw=nw0,nw1
                     psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw)=psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw) +reshape(source=recvbuffer(ibuf:ibufnext-1),&
                         shape=shape(psb(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                        iw))) /ps(igrid)%dvolume(ixmin1:ixmax1,ixmin2:ixmax2)
                     ibuf=ibuf+nbuf
                   end do
                   ibuf=ibufnext
                 end if
               end if
            end do
           end do
           end do
         end select
       end do
     end do

     if (nsend>0) then
       call MPI_WAITALL(nsend,fc_sendreq,fc_sendstat,ierrmpi)
     end if

   end subroutine fix_conserve

   subroutine store_flux(igrid,fC,idimmin,idimmax,nwfluxin)
     use mod_global_parameters

     integer, intent(in)          :: igrid, idimmin,idimmax, nwfluxin
     double precision, intent(in) :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwfluxin,&
        1:ndim)

     integer :: idims, iside, i1,i2, ic1,ic2, inc1,inc2, ix1,ix2, ixCo1,ixCo2,&
         nxCo1,nxCo2, iw

     do idims = idimmin,idimmax
       select case (idims)
         case (1)
         do iside=1,2
           i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);

           if (neighbor_pole(i1,i2,igrid)/=0) cycle

           select case (neighbor_type(i1,i2,igrid))
           case (neighbor_fine)
             select case (iside)
             case (1)
               pflux(iside,1,igrid)%flux(1,:,1:nwfluxin) = -fC(nghostcells,&
                  ixMlo2:ixMhi2,1:nwfluxin,1)
             case (2)
               pflux(iside,1,igrid)%flux(1,:,1:nwfluxin) = fC(ixMhi1,&
                  ixMlo2:ixMhi2,1:nwfluxin,1)
             end select
           case (neighbor_coarse)
             nxCo1=1;nxCo2=ixGhi2/2-nghostcells;
             select case (iside)
             case (1)
               do iw=1,nwfluxin
                do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=nghostcells;ix2=ixMlo2+2*(ixCo2-1);
                   pflux(iside,1,igrid)%flux(ixCo1,ixCo2,iw) = sum(fC(ix1,&
                      ix2:ix2+1,iw,1))
                end do
         end do
               end do
             case (2)
               do iw=1,nwfluxin
                do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=ixMhi1;ix2=ixMlo2+2*(ixCo2-1);
                   pflux(iside,1,igrid)%flux(ixCo1,ixCo2,iw) =-sum(fC(ix1,&
                      ix2:ix2+1,iw,1))
                end do
         end do
               end do
             end select
           end select
         end do
         case (2)
         do iside=1,2
           i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);

           if (neighbor_pole(i1,i2,igrid)/=0) cycle

           select case (neighbor_type(i1,i2,igrid))
           case (neighbor_fine)
             select case (iside)
             case (1)
               pflux(iside,2,igrid)%flux(:,1,1:nwfluxin) = -fC(ixMlo1:ixMhi1,&
                  nghostcells,1:nwfluxin,2)
             case (2)
               pflux(iside,2,igrid)%flux(:,1,1:nwfluxin) = fC(ixMlo1:ixMhi1,&
                  ixMhi2,1:nwfluxin,2)
             end select
           case (neighbor_coarse)
             nxCo1=ixGhi1/2-nghostcells;nxCo2=1;
             select case (iside)
             case (1)
               do iw=1,nwfluxin
                do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=ixMlo1+2*(ixCo1-1);ix2=nghostcells;
                   pflux(iside,2,igrid)%flux(ixCo1,ixCo2,&
                      iw) = sum(fC(ix1:ix1+1,ix2,iw,2))
                end do
         end do
               end do
             case (2)
               do iw=1,nwfluxin
                do ixCo2=1,nxCo2
         do ixCo1=1,nxCo1
                   ix1=ixMlo1+2*(ixCo1-1);ix2=ixMhi2;
                   pflux(iside,2,igrid)%flux(ixCo1,ixCo2,&
                      iw) =-sum(fC(ix1:ix1+1,ix2,iw,2))
                end do
         end do
               end do
             end select
           end select
         end do
       end select
     end do

   end subroutine store_flux

   subroutine store_edge(igrid,ixImin1,ixImin2,ixImax1,ixImax2,fE,idimmin,&
      idimmax)
     use mod_global_parameters
     
     integer, intent(in)          :: igrid, ixImin1,ixImin2,ixImax1,ixImax2,&
         idimmin,idimmax
     double precision, intent(in) :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
        7-2*ndim:3)
     
     integer :: idims, idir, iside, i1,i2
     integer :: pi1,pi2, mi1,mi2, ph1,ph2, mh1,mh2 ! To detect corners
     integer :: ixMcmin1,ixMcmin2,ixMcmax1,ixMcmax2

     do idims = idimmin,idimmax  !loop over face directions
       !! Loop over block faces
       do iside=1,2 
         i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3);
         if (neighbor_pole(i1,i2,igrid)/=0) cycle
         select case (neighbor_type(i1,i2,igrid))
         case (neighbor_fine)
           ! The neighbour is finer
           ! Face direction, side (left or right), restrict ==ired?, fE
           call flux_to_edge(igrid,ixImin1,ixImin2,ixImax1,ixImax2,idims,iside,&
              .false.,fE)
         case(neighbor_coarse)
           ! The neighbour is coarser
           call flux_to_edge(igrid,ixImin1,ixImin2,ixImax1,ixImax2,idims,iside,&
              .true.,fE)
         case(neighbor_sibling)
           ! If the neighbour is at the same level,
           ! check if there are corners
           ! If there is any corner, store the fluxes from that side
           do idir=idims+1,ndim
             pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);
             mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);
             ph1=pi1-kr(idims,1)*(2*iside-3);ph2=pi2-kr(idims,2)*(2*iside-3);
             mh1=mi1-kr(idims,1)*(2*iside-3);mh2=mi2-kr(idims,2)*(2*iside-3);
             if (neighbor_type(pi1,pi2,igrid)==4.and.neighbor_type(ph1,ph2,&
                igrid)==3) then
               call flux_to_edge(igrid,ixImin1,ixImin2,ixImax1,ixImax2,idims,&
                  iside,.false.,fE)
             end if
             if (neighbor_type(mi1,mi2,igrid)==4.and.neighbor_type(mh1,mh2,&
                igrid)==3) then
               call flux_to_edge(igrid,ixImin1,ixImin2,ixImax1,ixImax2,idims,&
                  iside,.false.,fE)
             end if
           end do
         end select
       end do
     end do

   end subroutine store_edge

   subroutine flux_to_edge(igrid,ixImin1,ixImin2,ixImax1,ixImax2,idims,iside,&
      restrict,fE)
     use mod_global_parameters

     integer                      :: igrid,ixImin1,ixImin2,ixImax1,ixImax2,&
        idims,iside
     logical                      :: restrict
     double precision, intent(in) :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
        7-2*ndim:3)

     integer                      :: idir1,idir2
     integer                      :: ixEmin1,ixEmin2,ixEmax1,ixEmax2,ixFmin1,&
        ixFmin2,ixFmax1,ixFmax2, nx1,nx2,nxCo1,nxCo2

     nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
     nxCo1=nx1/2;nxCo2=nx2/2;
     ! ixE are the indices on the 'edge' array.
     ! ixF are the indices on the 'fE' array
     ! jxF are indices advanced to perform the flux restriction (sum) in 3D
     ! A line integral of the electric field on the coarse side
     ! lies over two edges on the fine side. So, in 3D we restrict by summing
     ! over two cells on the fine side.

     do idir1=1,ndim-1
      
       ! Assign the only field component (3) to the only index (1)
       idir2=3

       if (restrict) then
         ! Set up indices for restriction
         ixFmin1=ixMlo1-1+kr(1,idir2);ixFmin2=ixMlo2-1+kr(2,idir2);
         ixFmax1=ixMhi1-kr(1,idir2);ixFmax2=ixMhi2-kr(2,idir2);
         ;

         ixEmin1=0+kr(1,idir2);ixEmin2=0+kr(2,idir2);
         ixEmax1=nxCo1;ixEmax2=nxCo2;
         select case(idims)
        case(1)
           ixEmin1=1;ixEmax1=1;
           select case(iside)
           case(1)
             ixFmax1=ixFmin1
             
           case(2)
             ixFmin1=ixFmax1
             
           end select
        
        case(2)
           ixEmin2=1;ixEmax2=1;
           select case(iside)
           case(1)
             ixFmax2=ixFmin2
             
           case(2)
             ixFmin2=ixFmax2
             
           end select
        
         end select

       pflux(iside,idims,igrid)%edge(ixEmin1:ixEmax1,ixEmin2:ixEmax2,&
          idir1)=fE(ixFmin1:ixFmax1:2,ixFmin2:ixFmax2:2,idir2);

       else
         ! Set up indices for copying 
         ixFmin1=ixMlo1-1+kr(1,idir2);ixFmin2=ixMlo2-1+kr(2,idir2);
         ixFmax1=ixMhi1;ixFmax2=ixMhi2;
         ixEmin1=0+kr(1,idir2);ixEmin2=0+kr(2,idir2);
         ixEmax1=nx1;ixEmax2=nx2;

         select case(idims)
        case(1)
           ixEmin1=1;ixEmax1=1;
           select case(iside)
           case(1)
             ixFmax1=ixFmin1
           case(2)
             ixFmin1=ixFmax1
           end select
        
        case(2)
           ixEmin2=1;ixEmax2=1;
           select case(iside)
           case(1)
             ixFmax2=ixFmin2
           case(2)
             ixFmin2=ixFmax2
           end select
        
         end select

         pflux(iside,idims,igrid)%edge(ixEmin1:ixEmax1,ixEmin2:ixEmax2,&
            idir1)=fE(ixFmin1:ixFmax1,ixFmin2:ixFmax2,idir2)

       end if

     end do

   end subroutine flux_to_edge

   subroutine fix_edges(psuse,idimmin,idimmax)
     use mod_global_parameters

     type(state) :: psuse(max_blocks)
     integer, intent(in) :: idimmin,idimmax

     integer :: iigrid, igrid, idims, iside, iotherside, i1,i2, ic1,ic2, inc1,&
        inc2, ixMcmin1,ixMcmin2,ixMcmax1,ixMcmax2
     integer :: nbuf, ibufnext
     integer :: ibufnext_cc
     integer :: pi1,pi2, mi1,mi2, ph1,ph2, mh1,mh2 ! To detect corners
     integer :: ixEmin1(1:3),ixEmin2(1:3),ixEmax1(1:3),ixEmax2(1:3), ixtEmin1,&
        ixtEmin2,ixtEmax1,ixtEmax2, ixFmin1(1:ndim),ixFmin2(1:ndim),&
        ixFmax1(1:ndim),ixFmax2(1:ndim), ixfEmin1(1:3),ixfEmin2(1:3),&
        ixfEmax1(1:3),ixfEmax2(1:3)
     integer :: nx1,nx2, idir, ix, ipe_neighbor, ineighbor
     logical :: pcorner(1:ndim),mcorner(1:ndim)

     if (nrecv_ct>0) then
        call MPI_WAITALL(nrecv_ct,cc_recvreq,cc_recvstat,ierrmpi)
     end if

     ! Initialise buffer counter again
     ibuf=1
     ibuf_cc=1
     do iigrid=1,igridstail; igrid=igrids(iigrid);
       do idims= idimmin,idimmax
         do iside=1,2
           i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3);
           if (neighbor_pole(i1,i2,igrid)/=0) cycle
           select case(neighbor_type(i1,i2,igrid))
           case(neighbor_fine)
             ! The first neighbour is finer
             if (.not.neighbor_active(i1,i2,igrid).or..not.neighbor_active(0,0,&
                igrid) ) then
               do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                  inc2=2*i2+ic2
               do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                  inc1=2*i1+ic1
                  ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
                  !! When the neighbour is in a different process
                  if (ipe_neighbor/=mype) then
                     ibufnext=ibuf+isize(idims)
                     ibuf=ibufnext
                     end if
               end do
               end do
                cycle
             end if

             ! Check if there are corners
             pcorner=.false.
             mcorner=.false.
             do idir=1,ndim
               pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);
               mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);
               ph1=pi1-kr(idims,1)*(2*iside-3)
               ph2=pi2-kr(idims,2)*(2*iside-3);
               mh1=mi1-kr(idims,1)*(2*iside-3)
               mh2=mi2-kr(idims,2)*(2*iside-3);
               if (neighbor_type(ph1,ph2,&
                  igrid)==neighbor_fine) pcorner(idir)=.true.
               if (neighbor_type(mh1,mh2,&
                  igrid)==neighbor_fine) mcorner(idir)=.true.
             end do
             ! Calculate indices range
             call set_ix_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,&
                ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,igrid,idims,iside,.false.,&
                .false.,0,0,pcorner,mcorner)
             ! Remove coarse part of circulation
             call add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,&
                ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,pflux(iside,idims,&
                igrid)%edge,idims,iside,.false.,psuse(igrid))
             ! Add fine part of the circulation
            do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
               inc2=2*i2+ic2
            do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
               inc1=2*i1+ic1
               ineighbor=neighbor_child(1,inc1,inc2,igrid)
               ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
               iotherside=3-iside
               nx1=(ixMhi1-ixMlo1+1)/2;nx2=(ixMhi2-ixMlo2+1)/2;
               call set_ix_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                  ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,&
                  ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,igrid,idims,iside,.true.,&
                  .false.,inc1,inc2,pcorner,mcorner)
               if (ipe_neighbor==mype) then
                 call add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                    ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,&
                    ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,pflux(iotherside,idims,&
                    ineighbor)%edge,idims,iside,.true.,psuse(igrid))
               else
                 ibufnext=ibuf+isize(idims)
                 call add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                    ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,&
                    ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,&
                    reshape(source=recvbuffer(ibufnext-&
                    isize_stg(idims):ibufnext-1),shape=(/ ixtEmax1-ixtEmin1+1,&
                    ixtEmax2-ixtEmin2+1 ,2-1 /)),idims,iside,.true.,&
                    psuse(igrid))
                 ibuf=ibufnext
               end if
            end do
            end do

           case(neighbor_sibling)
             ! The first neighbour is at the same level
             ! Check if there are corners
             do idir=idims+1,ndim
               pcorner=.false.
               mcorner=.false.
               pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);
               mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);
               ph1=pi1-kr(idims,1)*(2*iside-3)
               ph2=pi2-kr(idims,2)*(2*iside-3);
               mh1=mi1-kr(idims,1)*(2*iside-3)
               mh2=mi2-kr(idims,2)*(2*iside-3);
               if (neighbor_type(pi1,pi2,&
                  igrid)==neighbor_fine.and.neighbor_type(ph1,ph2,&
                  igrid)==neighbor_sibling.and.neighbor_pole(pi1,pi2,&
                  igrid)==0) then
                 pcorner(idir)=.true.
                 ! Remove coarse part
                 ! Set indices
                 call set_ix_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                    ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,&
                    ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,igrid,idims,iside,&
                    .false.,.true.,0,0,pcorner,mcorner)
                 ! Remove
                 call add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                    ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,&
                    ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,pflux(iside,idims,&
                    igrid)%edge,idims,iside,.false.,psuse(igrid))
                 ! Add fine part
                 ! Find relative position of finer grid
      
                 inc1=kr(idims,1)*3*(iside-1)+3*kr(idir,1)
                 inc2=kr(idims,2)*3*(iside-1)+3*kr(idir,2);
                 ineighbor=neighbor_child(1,inc1,inc2,igrid)
                 ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
                 iotherside=3-iside
                 ! Set indices
                 call set_ix_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                    ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,&
                    ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,igrid,idims,iside,&
                    .true.,.true.,inc1,inc2,pcorner,mcorner)
                 ! add
                 if (ipe_neighbor==mype) then
                   call add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                      ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,&
                      ixEmax2,ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,&
                      pflux(iotherside,idims,ineighbor)%edge,idims,iside,&
                      .true.,psuse(igrid))
                 else
                   ibufnext_cc=ibuf_cc+isize_stg(idims)
                   call add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                      ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,&
                      ixEmax2,ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,&
                      reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
                      shape=(/ ixtEmax1-ixtEmin1+1,ixtEmax2-ixtEmin2+1 ,&
                      2-1 /)),idims,iside,.true.,psuse(igrid))
                   ibuf_cc=ibufnext_cc
                 end if
      
               ! Set CoCorner to false again for next step
                 pcorner(idir)=.false.
               end if

               if (neighbor_type(mi1,mi2,&
                  igrid)==neighbor_fine.and.neighbor_type(mh1,mh2,&
                  igrid)==neighbor_sibling.and.neighbor_pole(mi1,mi2,&
                  igrid)==0) then
                   mcorner(idir)=.true.
                   ! Remove coarse part
                   ! Set indices
                   call set_ix_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                      ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,&
                      ixEmax2,ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,igrid,idims,&
                      iside,.false.,.true.,0,0,pcorner,mcorner)
                   ! Remove
                   call add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                      ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,&
                      ixEmax2,ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,pflux(iside,&
                      idims,igrid)%edge,idims,iside,.false.,psuse(igrid))
                   ! Add fine part
                   ! Find relative position of finer grid
        
                   inc1=kr(idims,1)*3*(iside-1);inc2=kr(idims,2)*3*(iside-1);
                   ineighbor=neighbor_child(1,inc1,inc2,igrid)
                   ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
                   iotherside=3-iside
                   ! Set indices
                   call set_ix_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,&
                      ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,&
                      ixEmax2,ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,igrid,idims,&
                      iside,.true.,.true.,inc1,inc2,pcorner,mcorner)
                   ! add
                   if (ipe_neighbor==mype) then
                     call add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,&
                        ixtEmin1,ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,&
                        ixEmax1,ixEmax2,ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,&
                        pflux(iotherside,idims,ineighbor)%edge,idims,iside,&
                        .true.,psuse(igrid))
                   else
                     ibufnext_cc=ibuf_cc+isize_stg(idims)
                     call add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,&
                        ixtEmin1,ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,&
                        ixEmax1,ixEmax2,ixfEmin1,ixfEmin2,ixfEmax1,ixfEmax2,&
                        reshape(source=recvbuffer_cc(ibuf_cc:ibufnext_cc-1),&
                        shape=(/ ixtEmax1-ixtEmin1+1,ixtEmax2-ixtEmin2+1 ,&
                        2-1 /)),idims,iside,.true.,psuse(igrid))
                     ibuf_cc=ibufnext_cc
                   end if
        
                 ! Set CoCorner to false again for next step
                  mcorner(idir)=.false.
               end if
             end do
           end select
         end do
       end do
     end do

     if (nsend_ct>0) call MPI_WAITALL(nsend_ct,cc_sendreq,cc_sendstat,ierrmpi)

   end subroutine fix_edges

   !> This routine sets the indexes for the correction
   !> of the circulation according to several different
   !> cases, as grids located in different cpus,
   !> presence of corners, and different relative locations
   !> of the fine grid respect to the coarse one
   subroutine set_ix_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,ixtEmin2,&
      ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,ixfEmin1,ixfEmin2,&
      ixfEmax1,ixfEmax2,igrid,idims,iside,add,CoCorner,inc1,inc2,pcorner,&
      mcorner)
     use mod_global_parameters
     
     integer,intent(in)    :: igrid,idims,iside,inc1,inc2
     logical,intent(in)    :: add,CoCorner
     logical,intent(inout) :: pcorner(1:ndim),mcorner(1:ndim)
     integer,intent(out)   :: ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmax1(1:ndim),&
        ixFmax2(1:ndim),ixtEmin1,ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1(1:3),&
        ixEmin2(1:3),ixEmax1(1:3),ixEmax2(1:3),ixfEmin1(1:3),ixfEmin2(1:3),&
        ixfEmax1(1:3),ixfEmax2(1:3) !Indices for faces and edges
     integer               :: icor1,icor2,idim1,idir,nx1,nx2,middle1,middle2
     integer               :: ixtfEmin1,ixtfEmin2,ixtfEmax1,ixtfEmax2

     ! ixF -> Indices for the _F_aces, and
     ! depends on the field component
     ! ixtE -> are the _t_otal range of the 'edge' array
     ! ixE -> are the ranges of the edge array,
     ! depending on the component
     ! ixfE -> are the ranges of the fE array (3D),
     ! and also depend on the component

     ! ... General ...
     ! Assign indices for the size of the E field array

     ixtfEmin1=ixMlo1-1;ixtfEmin2=ixMlo2-1;
     ixtfEmax1=ixMhi1;ixtfEmax2=ixMhi2;

     if(add) then
       nx1=(ixMhi1-ixMlo1+1)/2;nx2=(ixMhi2-ixMlo2+1)/2;
     else
       nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
     end if

     do idim1=1,ndim
       ixtEmin1=0;ixtEmin2=0;
       ixtEmax1=nx1;ixtEmax2=nx2;
       select case(idims)
       case(1)
         ixtEmin1=1;ixtEmax1=1;
         if (iside==1) ixtfEmax1=ixtfEmin1;
         if (iside==2) ixtfEmin1=ixtfEmax1;
       
       case(2)
         ixtEmin2=1;ixtEmax2=1;
         if (iside==1) ixtfEmax2=ixtfEmin2;
         if (iside==2) ixtfEmin2=ixtfEmax2;
       
       end select
     end do

     ! Assign indices, considering only the face
     ! (idims and iside)
     do idim1=1,ndim
       ixFmin1(idim1)=ixMlo1-kr(idim1,1);ixFmin2(idim1)=ixMlo2-kr(idim1,2);
       ixFmax1(idim1)=ixMhi1;ixFmax2(idim1)=ixMhi2;
       select case(idims)
       case(1)
          select case(iside)
          case(1)
          ixFmax1(idim1)=ixFmin1(idim1)
          case(2)
          ixFmin1(idim1)=ixFmax1(idim1)
          end select
       
       case(2)
          select case(iside)
          case(1)
          ixFmax2(idim1)=ixFmin2(idim1)
          case(2)
          ixFmin2(idim1)=ixFmax2(idim1)
          end select
       
       end select
     end do
     ! ... Relative position ...
     ! Restrict range using relative position
     if(add) then
       middle1=(ixMhi1+ixMlo1)/2;middle2=(ixMhi2+ixMlo2)/2;
       
       if(inc1==1) then
         ixFmax1(:)=middle1
         ixtfEmax1=middle1
       end if
       if(inc1==2) then
         ixFmin1(:)=middle1+1
         ixtfEmin1=middle1
       end if
       
       
       if(inc2==1) then
         ixFmax2(:)=middle2
         ixtfEmax2=middle2
       end if
       if(inc2==2) then
         ixFmin2(:)=middle2+1
         ixtfEmin2=middle2
       end if
       
     end if
     ! ... Adjust ranges of edges according to direction ...
     do idim1=1,3
       ixfEmax1(idim1)=ixtfEmax1;ixfEmax2(idim1)=ixtfEmax2;
       ixEmax1(idim1)=ixtEmax1;ixEmax2(idim1)=ixtEmax2;
       ixfEmin1(idim1)=ixtfEmin1+kr(idim1,1)
       ixfEmin2(idim1)=ixtfEmin2+kr(idim1,2);
       ixEmin1(idim1)=ixtEmin1+kr(idim1,1)
       ixEmin2(idim1)=ixtEmin2+kr(idim1,2);
     end do
     ! ... Corners ...
     ! 'Coarse' corners
     if (CoCorner) then
       do idim1=idims+1,ndim
         if (pcorner(idim1)) then
           do idir=1,3!Index arrays have size ndim
             if (idir==6-idim1-idims) then
              !!! Something here has to change
              !!! Array ixfE must have size 3, while
              !!! ixE must have size ndim
              if (1==idim1) then
                 ixfEmin1(idir)=ixfEmax1(idir)
                 if (add) then
                   ixEmax1(idir) =ixEmin1(idir)
                 else
                   ixEmin1(idir) =ixEmax1(idir)
                 end if
               end if
              if (2==idim1) then
                 ixfEmin2(idir)=ixfEmax2(idir)
                 if (add) then
                   ixEmax2(idir) =ixEmin2(idir)
                 else
                   ixEmin2(idir) =ixEmax2(idir)
                 end if
               end if
             else
               ixEmin1(idir)=1;ixEmin2(idir)=1;
               ixEmax1(idir)=0;ixEmax2(idir)=0;
               ixfEmin1(idir)=1;ixfEmin2(idir)=1;
               ixfEmax1(idir)=0;ixfEmax2(idir)=0;
             end if
           end do
         end if
         if (mcorner(idim1)) then
           do idir=1,3
             if (idir==6-idim1-idims) then
              if (1==idim1) then
                 ixfEmax1(idir)=ixfEmin1(idir)
                 if (add) then
                   ixEmin1(idir) =ixEmax1(idir)
                 else
                   ixEmax1(idir) =ixEmin1(idir)
                 end if
               end if
              if (2==idim1) then
                 ixfEmax2(idir)=ixfEmin2(idir)
                 if (add) then
                   ixEmin2(idir) =ixEmax2(idir)
                 else
                   ixEmax2(idir) =ixEmin2(idir)
                 end if
               end if
             else
               ixEmin1(idir)=1;ixEmin2(idir)=1;
               ixEmax1(idir)=0;ixEmax2(idir)=0;
               ixfEmin1(idir)=1;ixfEmin2(idir)=1;
               ixfEmax1(idir)=0;ixfEmax2(idir)=0;
             end if
           end do
         end if
       end do
     else
     ! Other kinds of corners
     ! Crop ranges to account for corners
     ! When the fine fluxes are added, we consider 
     ! whether they come from the same cpu or from
     ! a different one, in order to minimise the 
     ! amount of communication
     ! Case for different processors still not implemented!!!
      if((idims.gt.1).and.pcorner(1)) then
         if((.not.add).or.(inc1==2)) then
           !ixFmax1(:)=ixFmax1(:)-kr(1,1);!ixFmax2(:)=ixFmax2(:)-kr(1,2);
           do idir=1,3
             if ((idir==idims).or.(idir==1)) cycle
               ixfEmax1(idir)=ixfEmax1(idir)-1
               ixEmax1(idir)=ixEmax1(idir)-1
           end do
         end if
       end if
      if((idims.gt.2).and.pcorner(2)) then
         if((.not.add).or.(inc2==2)) then
           !ixFmax1(:)=ixFmax1(:)-kr(2,1);!ixFmax2(:)=ixFmax2(:)-kr(2,2);
           do idir=1,3
             if ((idir==idims).or.(idir==2)) cycle
               ixfEmax2(idir)=ixfEmax2(idir)-1
               ixEmax2(idir)=ixEmax2(idir)-1
           end do
         end if
       end if
      if((idims>1).and.mcorner(1)) then
         if((.not.add).or.(inc1==1)) then
           !ixFmin1(:)=ixFmin1(:)+kr(1,1);!ixFmin2(:)=ixFmin2(:)+kr(1,2);
           do idir=1,3
             if ((idir==idims).or.(idir==1)) cycle
               ixfEmin1(idir)=ixfEmin1(idir)+1
               ixEmin1(idir)=ixEmin1(idir)+1
           end do
         end if
       end if
      if((idims>2).and.mcorner(2)) then
         if((.not.add).or.(inc2==1)) then
           !ixFmin1(:)=ixFmin1(:)+kr(2,1);!ixFmin2(:)=ixFmin2(:)+kr(2,2);
           do idir=1,3
             if ((idir==idims).or.(idir==2)) cycle
               ixfEmin2(idir)=ixfEmin2(idir)+1
               ixEmin2(idir)=ixEmin2(idir)+1
           end do
         end if
       end if
     end if

   end subroutine set_ix_circ

   subroutine add_sub_circ(ixFmin1,ixFmin2,ixFmax1,ixFmax2,ixtEmin1,ixtEmin2,&
      ixtEmax1,ixtEmax2,ixEmin1,ixEmin2,ixEmax1,ixEmax2,ixfEmin1,ixfEmin2,&
      ixfEmax1,ixfEmax2,edge,idims,iside,add,s)
     use mod_global_parameters

     type(state)        :: s
     integer,intent(in) :: idims,iside
     integer            :: ixFmin1(1:ndim),ixFmin2(1:ndim),ixFmax1(1:ndim),&
        ixFmax2(1:ndim),ixtEmin1,ixtEmin2,ixtEmax1,ixtEmax2,ixEmin1(1:3),&
        ixEmin2(1:3),ixEmax1(1:3),ixEmax2(1:3),ixfEmin1(1:3),ixfEmin2(1:3),&
        ixfEmax1(1:3),ixfEmax2(1:3)
     double precision   :: edge(ixtEmin1:ixtEmax1,ixtEmin2:ixtEmax2,1:ndim-1)
     logical,intent(in) :: add

     integer            :: idim1,idim2,idir,middle1,middle2
     integer            :: ixfECmin1,ixfECmin2,ixfECmax1,ixfECmax2,ixECmin1,&
        ixECmin2,ixECmax1,ixECmax2
     double precision   :: fE(ixGlo1:ixGhi1,ixGlo2:ixGhi2,7-2*ndim:3) !!!!!!!!
     double precision   :: circ(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim) !!!!!!!!
     integer            :: ixmin1,ixmin2,ixmax1,ixmax2,hxmin1,hxmin2,hxmax1,&
        hxmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2 !Indices for edges

     ! ixF -> Indices for the faces, depends on the field component
     ! ixE -> Total range for the edges
     ! ixfE -> Edges in fE (3D) array
     ! ix,hx,ixC,hxC -> Auxiliary indices
     ! Assign quantities stored ad edges to make it as similar as 
     ! possible to the routine updatefaces.
     fE(:,:,:)=zero
     do idim1=1,ndim-1
       
          ! 2D: move E back to directon 3
       idir=3
       ixfECmin1=ixfEmin1(idir);ixfECmin2=ixfEmin2(idir)
       ixfECmax1=ixfEmax1(idir);ixfECmax2=ixfEmax2(idir);
       ixECmin1=ixEmin1(idir);ixECmin2=ixEmin2(idir);ixECmax1=ixEmax1(idir)
       ixECmax2=ixEmax2(idir);
       fE(ixfECmin1:ixfECmax1,ixfECmin2:ixfECmax2,idir)=edge(ixECmin1:ixECmax1,&
          ixECmin2:ixECmax2,idim1)
     end do

     ! Calculate part of circulation needed
     circ=zero
     do idim1=1,ndim
        do idim2=1,ndim
           do idir=7-2*ndim,3
             if (lvc(idim1,idim2,idir)==0) cycle
             ! Assemble indices
             ixCmin1=ixFmin1(idim1);ixCmin2=ixFmin2(idim1)
             ixCmax1=ixFmax1(idim1);ixCmax2=ixFmax2(idim1);
             hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
             hxCmax1=ixCmax1-kr(idim2,1);hxCmax2=ixCmax2-kr(idim2,2);
             if(idim1==idims) then
               circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)+lvc(idim1,&
                  idim2,idir)*(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  idir)-fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,idir))
             else
               select case(iside)
               case(2)
                 circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    idim1)+lvc(idim1,idim2,idir)*fE(ixCmin1:ixCmax1,&
                    ixCmin2:ixCmax2,idir)
               case(1)
                 circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                    idim1)-lvc(idim1,idim2,idir)*fE(hxCmin1:hxCmax1,&
                    hxCmin2:hxCmax2,idir)
               end select
             end if
           end do
        end do
     end do

     ! Divide circulation by surface and add
     do idim1=1,ndim
        ixCmin1=ixFmin1(idim1);ixCmin2=ixFmin2(idim1);ixCmax1=ixFmax1(idim1)
        ixCmax2=ixFmax2(idim1);
        where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           idim1)>1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=circ(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             idim1)
        elsewhere
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=zero
        end where
        ! Add/subtract to field at face
        if (add) then
          s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=s%ws(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             idim1)
        else
          s%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=s%ws(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idim1)+circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             idim1)
        end if
     end do

   end subroutine add_sub_circ

end module mod_fix_conserve
