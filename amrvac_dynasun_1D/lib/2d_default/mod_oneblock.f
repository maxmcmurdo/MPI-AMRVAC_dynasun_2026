!=============================================================================
! Module oneblock to hold a datastructure as output by convert_type='oneblock'. 
! Currently only ASCII
! Can be handy to read in initial conditions.
! 2015-02-18 Oliver Porth
!=============================================================================
module mod_oneblock
implicit none
save

double precision, dimension(:,:,:), allocatable       :: woneblock
double precision, dimension(:,:,:), allocatable       :: xoneblock
integer                                                    :: nc1,nc2
integer                                                    :: unit=15 !file unit to read on
!-----------------------------------------------------------------------------


contains
!=============================================================================
subroutine read_oneblock(filename)

use mod_global_parameters

character(len=*), intent(in)             :: filename
! .. local ..
character(len=1024)                      :: outfilehead
integer                                  :: nctot, ix1,ix2, ixp1,ixp2
double precision                         :: time
integer                                  :: idim
integer,dimension(2)                   :: sendbuff
!-----------------------------------------------------------------------------

if (mype == 0) write(*,*) 'mod_oneblock: reading ',nw,' variables on unit:',&
    unit

!----------------------------------------
! Root does the reading:
!----------------------------------------
if (mype == 0) then 
   open(unit,file=filename,status='unknown')
   ! The header information:
   read(unit,'(A)') outfilehead
   read(unit,*) nctot,nc1,nc2
   read(unit,*) time
   
   ! Allocate and read the grid and variables:
   allocate(xoneblock(nc1,nc2,1:2))
   allocate(woneblock(nc1,nc2,1:nw))
   
   do ix2=1,nc2
   do ix1=1,nc1  
   
   read(unit,*) xoneblock(ix1,ix2,1:2), woneblock(ix1,ix2,1:nw)
   
   end do
   end do
   
! Close the file
   close(unit)
end if! mype==0

!----------------------------------------
! Broadcast what mype=0 read:
!----------------------------------------
if (npe>1) then 
   sendbuff(1)=nc1;sendbuff(2)=nc2;
   call MPI_BCAST(sendbuff,2,MPI_INTEGER,0,icomm,ierrmpi)
   if (mype .ne. 0) then 
      nc1=sendbuff(1);nc2=sendbuff(2);
      ! Allocate the grid and variables:
      allocate(xoneblock(nc1,nc2,1:2))
      allocate(woneblock(nc1,nc2,1:nw))
   end if
   call MPI_BCAST(xoneblock,nc1*nc2*2,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
   call MPI_BCAST(woneblock,nc1*nc2*nw,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
end if! npe>1



end subroutine read_oneblock
!=============================================================================
subroutine interpolate_oneblock(x,iw,out)

double precision, dimension(2),intent(in)             :: x 
integer, intent(in)                                     :: iw
double precision, intent(out)                           :: out
! .. local ..
double precision, dimension(2)                        :: xloc
integer                                                 :: ic1,ic2, ic11,ic12,&
    ic21,ic22
double precision                                        :: xd1,xd2

double precision                                        :: c00, c10


integer                                                 :: ipivot1,ipivot2,&
    idir
!-----------------------------------------------------------------------------

xloc=x

!--------------------------------------------
! Hunt for the index closest to the point
! This is a bit slow but allows for stretched grids
! (still need to be orthogonal for interpolation though)
!--------------------------------------------
ipivot1=1;ipivot2=1;


do idir = 1, 2
   select case (idir)
   case (1)
      ic1 = minloc(dabs(xloc(1)-xoneblock(:,ipivot2,idir)),1,mask=.true.)
   case (2)
      ic2 = minloc(dabs(xloc(2)-xoneblock(ipivot1,:,idir)),1,mask=.true.)
   case default
      call mpistop("error1 in interpolate_oneblock")
   end select
end do



! flat interpolation would simply be:
!out = woneblock(ic^D,iw)
!return

!-------------------------------------------
! Get the left and right indices
!-------------------------------------------

if (xoneblock(ic1,ic2,1) .lt. xloc(1)) then
   ic11 = ic1
else
   ic11 = ic1 -1
end if
ic21 = ic11 + 1


if (xoneblock(ic1,ic2,2) .lt. xloc(2)) then
   ic12 = ic2
else
   ic12 = ic2 -1
end if
ic22 = ic12 + 1


!--------------------------------------------
! apply flat interpolation if outside of range, 
! change point-location to make this easy!
!--------------------------------------------

if (ic11 .lt. 1) then
   ic11 = 1
   ic21 = ic11 + 1
   xloc(1) = xoneblock(ic1,ic2,1)
end if


if (ic12 .lt. 1) then
   ic12 = 1
   ic22 = ic12 + 1
   xloc(2) = xoneblock(ic1,ic2,2)
end if


if (ic21 .gt. nc1) then
   ic21 = nc1
   ic11 = ic21 - 1
   xloc(1) = xoneblock(ic1,ic2,1)
end if


if (ic22 .gt. nc2) then
   ic22 = nc2
   ic12 = ic22 - 1
   xloc(2) = xoneblock(ic1,ic2,2)
end if



!-------------------------------------------
! linear, bi- and tri- linear interpolations
!-------------------------------------------


xd1 = (xloc(1)-xoneblock(ic11,ic12,1)) / (xoneblock(ic21,ic12,&
   1) - xoneblock(ic11,ic12,1))      
xd2 = (xloc(2)-xoneblock(ic11,ic12,2)) / (xoneblock(ic11,ic22,&
   2) - xoneblock(ic11,ic12,2))
c00 = woneblock(ic11,ic12,iw) * (1.0d0 - xd1) + woneblock(ic21,ic12,iw) * xd1
c10 = woneblock(ic11,ic22,iw) * (1.0d0 - xd1) + woneblock(ic21,ic22,iw) * xd1
out = c00 * (1.0d0 - xd2) + c10 * xd2



end subroutine interpolate_oneblock
!=============================================================================
end module mod_oneblock
!=============================================================================
