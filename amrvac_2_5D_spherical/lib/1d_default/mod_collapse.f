!> Collapses D-dimensional output to D-1 view by line-of-sight integration
! for now: only for non-stretched cartesian cases, with LOS along coordinate
! writes out as either csv or vti file (collapse_type)
module mod_collapse
  use mod_global_parameters
  implicit none

contains
!=============================================================================
subroutine write_collapsed
use mod_global_parameters
! Writes a collapsed view of the data integrated over one grid-direction.  
! E.g. column density maps.  
! Uses flat interpolation throughout.
! by Oliver Porth
! 6.Nov 2013
integer :: idir
logical, save :: firstcollapse=.true.
!-----------------------------------------------------------------------------

if (firstcollapse) then
   firstcollapse=.false.
end if

do idir=1, ndim
   if (collapse(idir)) call put_collapse(idir)
end do

collapsenext=collapsenext+1
end subroutine write_collapsed
!=============================================================================
subroutine put_collapse(dir)
use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_slice, only:alloc_subnode, dealloc_subnode
use mod_global_parameters
integer, intent(in)                               :: dir
! .. local ..
integer                                           :: jgrid, igrid, Morton_no
double precision,dimension(0:nw+nwauxio)          :: normconv 
!-----------------------------------------------------------------------------
if(.not.slab) call mpistop("collapse only for slab cartesian cases")

call allocate_collapsed(dir)

jgrid=0
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   jgrid=jgrid+1
   call alloc_subnode(jgrid,dir,nwauxio)
   call collapse_subnode(igrid,jgrid,dir,normconv)
   call integrate_subnode(igrid,jgrid,dir)
end do

! Reduce to head-node:
if (mype==0) then
   call MPI_REDUCE(MPI_IN_PLACE,collapsedData,size(collapsedData),&
      MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
else
   call MPI_REDUCE(collapsedData,collapsedData,size(collapsedData),&
      MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
end if
call MPI_BARRIER(icomm, ierrmpi)


call mpistop&
("sorry, 1D collapsed output routine not yet implemented (should be easy)...")


! If we need the subnodes later, remove deallocation here:
do jgrid=1,Morton_stop(mype)-Morton_start(mype)+1
   call dealloc_subnode(jgrid)
end do
deallocate(collapsedData)

end subroutine put_collapse
!=============================================================================
subroutine output_collapsed_csv(dir,normconv)
use mod_global_parameters
use mod_calculate_xw
integer, intent(in)                               :: dir
double precision,dimension(0:nw+nwauxio),intent(in):: normconv 
character(len=1024) :: filename, outfilehead, line
logical             :: fileopen
integer                                           :: iw
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
integer, dimension(ndim)                          :: myshape

!-----------------------------------------------------------------------------

if (mype==0) then

! Get coordinates:



 inquire(unitcollapse,opened=fileopen)
 if(.not.fileopen)then
      ! generate filename: 
      write(filename,"(a,i1.1,a,i1.1,a,i4.4,a)") trim(base_filename) // '_d',&
          dir,'_l',collapseLevel,'_n',collapsenext,'.csv'
      open(unitcollapse,file=filename,status='unknown',form='formatted')
   end if
   ! get and write the header: 
   call getheadernames(wnamei,xandwnamei,outfilehead)
   line=''
   do iw=1,ndim+nw+nwauxio-1
      if (iw .eq. dir) cycle
      line = trim(line)//trim(xandwnamei(iw))//', '
   end do
   line = trim(line)//trim(xandwnamei(ndim+nw+nwauxio))
   write(unitcollapse,'(a)')trim(line)
   myshape = shape(collapsedData)

close(unitcollapse)

end if
end subroutine output_collapsed_csv
!=============================================================================
subroutine output_collapsed_vti(dir,normconv)
use mod_global_parameters
use mod_calculate_xw
integer, intent(in)                               :: dir
double precision,dimension(0:nw+nwauxio),intent(in):: normconv 
character(len=1024) :: filename, outfilehead, line
logical             :: fileopen
integer                                           :: iw
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
integer, dimension(ndim)                          :: myshape

double precision                                  :: origin(1:3), spacing(1:3)
integer                                           :: wholeExtent(1:6),&
    size_single, length, size_length
integer*8                                         :: offset
character                                         :: buf
!-----------------------------------------------------------------------------
if (mype==0) then
offset=0
size_single=4
size_length=4
! Get coordinates:



origin=0
spacing=zero
wholeExtent=0
myshape = shape(collapsedData)

length = size_single



 inquire(unitcollapse,opened=fileopen)
 if(.not.fileopen)then
      ! generate filename: 
    write(filename,"(a,i1.1,a,i1.1,a,i4.4,a)") trim(base_filename)//'_d',dir,&
       '_l',collapseLevel,'_n',collapsenext,'.vti'
      open(unitcollapse,file=filename,status='unknown',form='formatted')
 end if
! get the header: 
call getheadernames(wnamei,xandwnamei,outfilehead)

! generate xml header
write(unitcollapse,'(a)')'<?xml version="1.0"?>'
write(unitcollapse,'(a)',advance='no') '<VTKFile type="ImageData"'
if(type_endian==1)then
   write(unitcollapse,'(a)')' version="0.1" byte_order="LittleEndian">'
else
  write(unitcollapse,'(a)')' version="0.1" byte_order="BigEndian">'
endif
! following corresponds to uniform cartesian grid
write(unitcollapse,'(a,3(1pe14.6),a,6(i10),a,3(1pe14.6),a)')&
   '  <ImageData Origin="',origin,'" WholeExtent="',wholeExtent,'" Spacing="',&
   spacing,'">'
write(unitcollapse,'(a)')'<FieldData>'
write(unitcollapse,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
   'NumberOfTuples="1" format="ascii">'
write(unitcollapse,*) real(global_time*time_convert_factor)
write(unitcollapse,'(a)')'</DataArray>'
write(unitcollapse,'(a)')'</FieldData>'

! we write one VTK PIECE
write(unitcollapse,'(a,6(i10),a)') '<Piece Extent="',wholeExtent,'">'
write(unitcollapse,'(a)')'<CellData>'

do iw=1,nw+nwauxio
  if(iw<=nw) then 
    if(.not.w_write(iw)) cycle
  endif
  write(unitcollapse,'(a,a,a,i16,a)')'<DataArray type="Float32" Name="',&
     trim(wnamei(iw)),'" format="appended" offset="',offset,'"/>'
  offset = offset + length + size_length
enddo

write(unitcollapse,'(a)')'</CellData>'
write(unitcollapse,'(a)')'</Piece>'
write(unitcollapse,'(a)')'</ImageData>'
write(unitcollapse,'(a)')'<AppendedData encoding="raw">'
close(unitcollapse)

open(unitcollapse,file=filename,access='stream',form='unformatted',&
   position='append')
buf='_'
write(unitcollapse) TRIM(buf)

do iw=1,nw+nwauxio
  if(iw<=nw) then 
    if(.not.w_write(iw)) cycle
  endif
  write(unitcollapse) length
  
  write(unitcollapse) real(collapsedData(iw)*normconv(iw))
 
  
  
enddo

close(unitcollapse)
open(unitcollapse,file=filename,status='unknown',form='formatted',&
   position='append')
write(unitcollapse,'(a)')'</AppendedData>'
write(unitcollapse,'(a)')'</VTKFile>'
close(unitcollapse)
 
end if

end subroutine output_collapsed_vti
!=============================================================================
subroutine allocate_collapsed(dir)
use mod_global_parameters
integer, intent(in)                               :: dir
integer                                           :: dim1
integer                                           :: nx1
!-----------------------------------------------------------------------------
! allocate array for the collapsed data:
! number of cells per grid.
nx1=ixMhi1-ixMlo1+1;
dim1=ng1(1)*2**(collapselevel-1)*nx1



   allocate(collapsedData(nw+nwauxio))

collapsedData = zero

end subroutine allocate_collapsed
!=============================================================================
subroutine integrate_subnode(igrid,jgrid,dir)
use mod_forest, only: tree_node_ptr, igrid_to_node
use mod_global_parameters
integer, intent(in)                        :: igrid, jgrid, dir
! .. local ..
type(tree_node_ptr)                        :: tree
integer                                    :: nx1
integer                                    :: ig1, level, ig1targetmin,&
    ig1targetmax
integer                                    :: ix1targetmin,ix1targetmax,&
    idim1targetmin,idim1targetmax
integer                                    :: ixMdimlo1,ixMdimhi1, ix1orig,&
    ix1

!-----------------------------------------------------------------------------
tree%node => igrid_to_node(igrid, mype)%node
 ig1 = tree%node%ig1;
level = tree%node%level
! number of cells per grid.
nx1=ixMhi1-ixMlo1+1;

! assumes uniform cartesian grid with factor 2 refinement per dimension

if (level > collapseLevel) then
   ig1targetmin = int(dble(ig1-1)/2.0d0**(level-collapseLevel))+1
   ig1targetmax = ig1targetmin
else if (level < collapseLevel) then
   ig1targetmin = int(2.0d0**(collapseLevel-level))*(ig1-1)+1
   ig1targetmax = int(2.0d0**(collapseLevel-level))*ig1
else
   ig1targetmin = ig1
   ig1targetmax = ig1
end if
ix1targetmin = nx1*(ig1targetmin-1)+1
ix1targetmax = nx1*ig1targetmax





collapsedData(1:nw+nwauxio) = collapsedData(1:nw+nwauxio) + &
   ps_sub(jgrid)%w(1:nw+nwauxio)


end subroutine integrate_subnode
!=============================================================================
subroutine collapse_subnode(igrid,jgrid,dir,normconv)
  use mod_usr_methods, only: usr_aux_output
  use mod_global_parameters
  use mod_physics, only: phys_to_primitive
  use mod_calculate_xw

integer, intent(in) :: igrid, jgrid, dir
double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv 
! .. local ..
integer                :: ix, iw
double precision       :: dx1
double precision, dimension(ixMlo1-1:ixMhi1,ndim)       :: xC_TMP, xC
double precision, dimension(ixMlo1:ixMhi1,ndim)         :: xCC_TMP, xCC
double precision, dimension(ixMlo1-1:ixMhi1,nw+nwauxio) :: wC_TMP
double precision, dimension(ixMlo1:ixMhi1,nw+nwauxio)   :: wCC_TMP
integer                :: ixCmin1,ixCmax1, ixCCmin1,ixCCmax1
!-----------------------------------------------------------------------------
ps_sub(jgrid)%w=zero
dx1=rnode(rpdx1_,igrid);

call calc_x(igrid,xC,xCC)
call calc_grid(unitslice,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
   ixCmin1,ixCmax1,ixCCmin1,ixCCmax1,.true.)



   
if(slab_uniform) then
  do ix=ixMlo1,ixMhi1
     ps_sub(jgrid)%w(1:nw+nwauxio) = ps_sub(jgrid)%w(1:nw+nwauxio) + &
        wCC_TMP(ix,1:nw+nwauxio) * dx1
  end do
else
  do iw=1,nw+nwauxio
    do ix=ixMlo1,ixMhi1
      ps_sub(jgrid)%w(iw) = ps_sub(jgrid)%w(iw) + wCC_TMP(ix,&
         iw) * ps(igrid)%dx(ix,1)
    end do
  end do
end if
ps_sub(jgrid)%x(1:ndim) = ps(igrid)%x(ixMlo1,1:ndim)


end subroutine collapse_subnode
!=============================================================================
end module mod_collapse
