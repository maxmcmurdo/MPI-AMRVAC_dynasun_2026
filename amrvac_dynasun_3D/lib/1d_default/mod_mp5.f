!> Module containing the MP5 (fifth order) flux scheme
module mod_mp5

  implicit none
  private

  public :: MP5limiter
  public :: MP5limitervar
  public :: MP5limiterL
  public :: MP5limiterR

contains

  !> MP5 limiter from Suresh & Huynh 1997 Following the convention of Mignone et
  !> al. 2010. Needs at least three ghost cells
  subroutine MP5limiter(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC,wRC)
    use mod_global_parameters
    use mod_weno, only: fix_limiter1

    integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)

    double precision, intent(inout) :: wRC(ixImin1:ixImax1,1:nw),&
       wLC(ixImin1:ixImax1,1:nw) 
    ! .. local ..
    integer                         :: iLmmin1,iLmmax1, iLmmmin1,iLmmmax1,&
        iLpmin1,iLpmax1, iLppmin1,iLppmax1, iLpppmin1,iLpppmax1
    integer                         :: idmin1,idmax1, idpmin1,idpmax1,&
        idppmin1,idppmax1, idmmin1,idmmax1, iemin1,iemax1, iemmin1,iemmax1,&
        iepmin1,iepmax1, ieppmin1,ieppmax1
    integer                         :: iw
    double precision, dimension(ixImin1:ixImax1,1:nw)  :: f, fmp, fmin, fmax,&
        ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixImin1:ixImax1,1:nw)  :: wRCtmp, wLCtmp
    double precision, dimension(ixImin1:ixImax1) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.d0, alpha=4.0d0
    !double precision                :: alpha
    !----------------------------------------------------------------------------

    ! Variable alpha:
    !alpha = float(nstep)/courantpar - one

    ! Left side:
    ! range to process:
    !iLmin^D=ixmin^D-kr(idims,^D);iLmax^D=ixmax^D;

    ! HALL
    ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
    ! also add one ghost zone!
    !   {iL^L=iL^L^LADD1;}

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

    iLmmin1=iLmin1-kr(idims,1);iLmmax1=iLmax1-kr(idims,1);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmax1=iLmmax1-kr(idims,1);
    iLpmin1=iLmin1+kr(idims,1);iLpmax1=iLmax1+kr(idims,1);
    iLppmin1=iLpmin1+kr(idims,1);iLppmax1=iLpmax1+kr(idims,1);

    f(iLmin1:iLmax1,1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLmmmin1:iLmmmax1,&
       1:nwflux) - 13.0d0* w(iLmmin1:iLmmax1,&
       1:nwflux) + 47.0d0* w(iLmin1:iLmax1,&
       1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,&
       1:nwflux) - 3.0d0*  w(iLppmin1:iLppmax1,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iLmin1:iLmax1) = w(iLpmin1:iLpmax1,iw)-w(iLmin1:iLmax1,iw)
       b(iLmin1:iLmax1) = alpha*(w(iLmin1:iLmax1,iw)-w(iLmmin1:iLmmax1,iw))
       call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
       fmp(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + tmp(iLmin1:iLmax1)
       ful(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + b(iLmin1:iLmax1)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1; idmin1=iLmin1-kr(idims,1);
    idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
    idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

    iemax1=idmax1+kr(idims,1); iemin1=idmin1;
    iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
    iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);

    d(iemin1:iemax1,1:nwflux) = w(iepmin1:iepmax1,&
       1:nwflux)-2.0d0*w(iemin1:iemax1,1:nwflux)+w(iemmin1:iemmax1,1:nwflux)

    do iw=1,nwflux
       a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idpmin1:idpmax1,iw)
       b(idmin1:idmax1) = 4.0d0*d(idpmin1:idpmax1,iw)-d(idmin1:idmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
       a(idmin1:idmax1) = d(idmin1:idmax1,iw)
       b(idmin1:idmax1) = d(idpmin1:idpmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,1:nwflux) = (w(iLmin1:iLmax1,1:nwflux)+w(iLpmin1:iLpmax1,&
       1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,1:nwflux)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,1:nwflux) = half*(3.0d0*w(iLmin1:iLmax1,&
       1:nwflux) - w(iLmmin1:iLmmax1,1:nwflux)) + &
       4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,1:nwflux)

    fmin(iLmin1:iLmax1,1:nwflux) = max(min(w(iLmin1:iLmax1,1:nwflux),&
       w(iLpmin1:iLpmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
       min(w(iLmin1:iLmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
       flc(iLmin1:iLmax1,1:nwflux)))

    fmax(iLmin1:iLmax1,1:nwflux) = min(max(w(iLmin1:iLmax1,1:nwflux),&
       w(iLpmin1:iLpmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
       max(w(iLmin1:iLmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
       flc(iLmin1:iLmax1,1:nwflux)))

    do iw=1,nwflux
       a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
       b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
       c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
       call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
       flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,1:nwflux)-w(iLmin1:iLmax1,&
       1:nwflux))*(f(iLmin1:iLmax1,1:nwflux)-fmp(iLmin1:iLmax1,&
       1:nwflux)) .le. eps)
       wLCtmp(iLmin1:iLmax1,1:nwflux) = f(iLmin1:iLmax1,1:nwflux)
    elsewhere
       wLCtmp(iLmin1:iLmax1,1:nwflux) = flim(iLmin1:iLmax1,1:nwflux)
    end where

    ! Right side:
    ! the interpolation from the right is obtained when the left-hand formula is applied to
    ! data mirrored about the interface.  
    ! thus substitute: 
    ! i-2 -> i+3
    ! i-1 -> i+2
    ! i   -> i+1
    ! i+1 -> i
    ! i+2 -> i-1

    iLpppmin1=iLppmin1+kr(idims,1);iLpppmax1=iLppmax1+kr(idims,1);

    f(iLmin1:iLmax1,1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLpppmin1:iLpppmax1,&
       1:nwflux) - 13.0d0* w(iLppmin1:iLppmax1,&
       1:nwflux) + 47.0d0* w(iLpmin1:iLpmax1,&
       1:nwflux) + 27.0d0* w(iLmin1:iLmax1,&
       1:nwflux) - 3.0d0*  w(iLmmin1:iLmmax1,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iLmin1:iLmax1) = w(iLmin1:iLmax1,iw)-w(iLpmin1:iLpmax1,iw)
       b(iLmin1:iLmax1) = alpha*(w(iLpmin1:iLpmax1,iw)-w(iLppmin1:iLppmax1,&
          iw))
       call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
       fmp(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + tmp(iLmin1:iLmax1)
       ful(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + b(iLmin1:iLmax1)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1+kr(idims,1); idmin1=iLmin1;
    idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
    idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

    iemax1=idmax1; iemin1=idmin1-kr(idims,1);
    iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
    iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);
    ieppmin1=iepmin1+kr(idims,1);ieppmax1=iepmax1+kr(idims,1);

    d(iemin1:iemax1,1:nwflux) = w(iemin1:iemax1,&
       1:nwflux)-2.0d0*w(iepmin1:iepmax1,1:nwflux)+w(ieppmin1:ieppmax1,&
       1:nwflux)

    do iw=1,nwflux
       a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idmmin1:idmmax1,iw)
       b(idmin1:idmax1) = 4.0d0*d(idmmin1:idmmax1,iw)-d(idmin1:idmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
       a(idmin1:idmax1) = d(idmin1:idmax1,iw)
       b(idmin1:idmax1) = d(idmmin1:idmmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,1:nwflux) = (w(iLmin1:iLmax1,1:nwflux)+w(iLpmin1:iLpmax1,&
       1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,1:nwflux)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,1:nwflux) = half*(3.0d0*w(iLpmin1:iLpmax1,&
       1:nwflux) - w(iLppmin1:iLppmax1,1:nwflux)) + &
       4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,1:nwflux)

    fmin(iLmin1:iLmax1,1:nwflux) = max(min(w(iLpmin1:iLpmax1,1:nwflux),&
       w(iLmin1:iLmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
       min(w(iLpmin1:iLpmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
       flc(iLmin1:iLmax1,1:nwflux)))

    fmax(iLmin1:iLmax1,1:nwflux) = min(max(w(iLpmin1:iLpmax1,1:nwflux),&
       w(iLmin1:iLmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
       max(w(iLpmin1:iLpmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
       flc(iLmin1:iLmax1,1:nwflux)))

    do iw=1,nwflux
       a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
       b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
       c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
       call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
       flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,1:nwflux)-w(iLpmin1:iLpmax1,&
       1:nwflux))*(f(iLmin1:iLmax1,1:nwflux)-fmp(iLmin1:iLmax1,&
       1:nwflux))  .le. eps)
       wRCtmp(iLmin1:iLmax1,1:nwflux) = f(iLmin1:iLmax1,1:nwflux)
    elsewhere
       wRCtmp(iLmin1:iLmax1,1:nwflux) = flim(iLmin1:iLmax1,1:nwflux)
    end where

    call fix_limiter1(ixImin1,ixImax1,iLmin1,iLmax1,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine MP5limiter

  subroutine MP5limiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC)
    use mod_global_parameters
    !use mod_weno, only: fix_onelimiter1

    integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)

    double precision, intent(inout) :: wLC(ixImin1:ixImax1,1:nw) 
    ! .. local ..
    integer                         :: iLmmin1,iLmmax1, iLmmmin1,iLmmmax1,&
        iLpmin1,iLpmax1, iLppmin1,iLppmax1
    integer                         :: idmin1,idmax1, idpmin1,idpmax1,&
        idppmin1,idppmax1, idmmin1,idmmax1, iemin1,iemax1, iemmin1,iemmax1,&
        iepmin1,iepmax1, ieppmin1,ieppmax1
    integer                         :: iw
    double precision, dimension(ixImin1:ixImax1,1:nw)  :: f, fmp, fmin, fmax,&
        ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixImin1:ixImax1) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.d0, alpha=4.0d0
    !double precision                :: alpha
    !double precision, dimension(ixI^S,1:nw)  :: wLCtmp

    ! Variable alpha:
    !alpha = float(nstep)/courantpar - one

    ! Left side:
    iLmmin1=iLmin1-kr(idims,1);iLmmax1=iLmax1-kr(idims,1);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmax1=iLmmax1-kr(idims,1);
    iLpmin1=iLmin1+kr(idims,1);iLpmax1=iLmax1+kr(idims,1);
    iLppmin1=iLpmin1+kr(idims,1);iLppmax1=iLpmax1+kr(idims,1);

    f(iLmin1:iLmax1,1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLmmmin1:iLmmmax1,&
       1:nwflux) - 13.0d0* w(iLmmin1:iLmmax1,&
       1:nwflux) + 47.0d0* w(iLmin1:iLmax1,&
       1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,&
       1:nwflux) - 3.0d0*  w(iLppmin1:iLppmax1,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iLmin1:iLmax1) = w(iLpmin1:iLpmax1,iw)-w(iLmin1:iLmax1,iw)
       b(iLmin1:iLmax1) = alpha*(w(iLmin1:iLmax1,iw)-w(iLmmin1:iLmmax1,iw))
       call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
       fmp(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + tmp(iLmin1:iLmax1)
       ful(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + b(iLmin1:iLmax1)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1; idmin1=iLmin1-kr(idims,1);
    idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
    idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

    iemax1=idmax1+kr(idims,1); iemin1=idmin1;
    iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
    iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);

    d(iemin1:iemax1,1:nwflux) = w(iepmin1:iepmax1,&
       1:nwflux)-2.0d0*w(iemin1:iemax1,1:nwflux)+w(iemmin1:iemmax1,1:nwflux)

    do iw=1,nwflux
       a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idpmin1:idpmax1,iw)
       b(idmin1:idmax1) = 4.0d0*d(idpmin1:idpmax1,iw)-d(idmin1:idmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
       a(idmin1:idmax1) = d(idmin1:idmax1,iw)
       b(idmin1:idmax1) = d(idpmin1:idpmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,1:nwflux) = (w(iLmin1:iLmax1,1:nwflux)+w(iLpmin1:iLpmax1,&
       1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,1:nwflux)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,1:nwflux) = half*(3.0d0*w(iLmin1:iLmax1,&
       1:nwflux) - w(iLmmin1:iLmmax1,1:nwflux)) + &
       4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,1:nwflux)

    fmin(iLmin1:iLmax1,1:nwflux) = max(min(w(iLmin1:iLmax1,1:nwflux),&
       w(iLpmin1:iLpmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
       min(w(iLmin1:iLmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
       flc(iLmin1:iLmax1,1:nwflux)))

    fmax(iLmin1:iLmax1,1:nwflux) = min(max(w(iLmin1:iLmax1,1:nwflux),&
       w(iLpmin1:iLpmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
       max(w(iLmin1:iLmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
       flc(iLmin1:iLmax1,1:nwflux)))

    do iw=1,nwflux
       a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
       b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
       c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
       call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
       flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,1:nwflux)-w(iLmin1:iLmax1,&
       1:nwflux))*(f(iLmin1:iLmax1,1:nwflux)-fmp(iLmin1:iLmax1,&
       1:nwflux)) .le. eps)
       wLC(iLmin1:iLmax1,1:nwflux) = f(iLmin1:iLmax1,1:nwflux)
    elsewhere
       wLC(iLmin1:iLmax1,1:nwflux) = flim(iLmin1:iLmax1,1:nwflux)
    end where

    !call fix_onelimiter1(ixI^L,iL^L,wLCtmp,wLC)

  end subroutine MP5limiterL

  subroutine MP5limiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC)
    use mod_global_parameters
    !use mod_weno, only: fix_onelimiter1

    integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)

    double precision, intent(inout) :: wRC(ixImin1:ixImax1,1:nw)
    ! .. local ..
    integer                         :: iLmmin1,iLmmax1, iLpmin1,iLpmax1,&
        iLppmin1,iLppmax1, iLpppmin1,iLpppmax1
    integer                         :: idmin1,idmax1, idpmin1,idpmax1,&
        idppmin1,idppmax1, idmmin1,idmmax1, iemin1,iemax1, iemmin1,iemmax1,&
        iepmin1,iepmax1, ieppmin1,ieppmax1
    integer                         :: iw
    double precision, dimension(ixImin1:ixImax1,1:nw)  :: f, fmp, fmin, fmax,&
        ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixImin1:ixImax1) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.d0, alpha=4.0d0
    !double precision                :: alpha
    !double precision, dimension(ixI^S,1:nw)  :: wRCtmp

    ! Right side:
    ! the interpolation from the right is obtained when the left-hand formula is applied to
    ! data mirrored about the interface.  
    ! thus substitute: 
    ! i-2 -> i+3
    ! i-1 -> i+2
    ! i   -> i+1
    ! i+1 -> i
    ! i+2 -> i-1

    iLmmin1=iLmin1-kr(idims,1);iLmmax1=iLmax1-kr(idims,1);
    iLpmin1=iLmin1+kr(idims,1);iLpmax1=iLmax1+kr(idims,1);
    iLppmin1=iLpmin1+kr(idims,1);iLppmax1=iLpmax1+kr(idims,1);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmax1=iLppmax1+kr(idims,1);

    f(iLmin1:iLmax1,1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLpppmin1:iLpppmax1,&
       1:nwflux) - 13.0d0* w(iLppmin1:iLppmax1,&
       1:nwflux) + 47.0d0* w(iLpmin1:iLpmax1,&
       1:nwflux) + 27.0d0* w(iLmin1:iLmax1,&
       1:nwflux) - 3.0d0*  w(iLmmin1:iLmmax1,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iLmin1:iLmax1) = w(iLmin1:iLmax1,iw)-w(iLpmin1:iLpmax1,iw)
       b(iLmin1:iLmax1) = alpha*(w(iLpmin1:iLpmax1,iw)-w(iLppmin1:iLppmax1,&
          iw))
       call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
       fmp(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + tmp(iLmin1:iLmax1)
       ful(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + b(iLmin1:iLmax1)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1+kr(idims,1); idmin1=iLmin1;
    idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
    idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

    iemax1=idmax1; iemin1=idmin1-kr(idims,1);
    iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
    iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);
    ieppmin1=iepmin1+kr(idims,1);ieppmax1=iepmax1+kr(idims,1);

    d(iemin1:iemax1,1:nwflux) = w(iemin1:iemax1,&
       1:nwflux)-2.0d0*w(iepmin1:iepmax1,1:nwflux)+w(ieppmin1:ieppmax1,&
       1:nwflux)

    do iw=1,nwflux
       a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idmmin1:idmmax1,iw)
       b(idmin1:idmax1) = 4.0d0*d(idmmin1:idmmax1,iw)-d(idmin1:idmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
       a(idmin1:idmax1) = d(idmin1:idmax1,iw)
       b(idmin1:idmax1) = d(idmmin1:idmmax1,iw)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
       call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,1:nwflux) = (w(iLmin1:iLmax1,1:nwflux)+w(iLpmin1:iLpmax1,&
       1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,1:nwflux)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,1:nwflux) = half*(3.0d0*w(iLpmin1:iLpmax1,&
       1:nwflux) - w(iLppmin1:iLppmax1,1:nwflux)) + &
       4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,1:nwflux)

    fmin(iLmin1:iLmax1,1:nwflux) = max(min(w(iLpmin1:iLpmax1,1:nwflux),&
       w(iLmin1:iLmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
       min(w(iLpmin1:iLpmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
       flc(iLmin1:iLmax1,1:nwflux)))

    fmax(iLmin1:iLmax1,1:nwflux) = min(max(w(iLpmin1:iLpmax1,1:nwflux),&
       w(iLmin1:iLmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
       max(w(iLpmin1:iLpmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
       flc(iLmin1:iLmax1,1:nwflux)))

    do iw=1,nwflux
       a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
       b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
       c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
       call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
       flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,1:nwflux)-w(iLpmin1:iLpmax1,&
       1:nwflux))*(f(iLmin1:iLmax1,1:nwflux)-fmp(iLmin1:iLmax1,&
       1:nwflux))  .le. eps)
       wRC(iLmin1:iLmax1,1:nwflux) = f(iLmin1:iLmax1,1:nwflux)
    elsewhere
       wRC(iLmin1:iLmax1,1:nwflux) = flim(iLmin1:iLmax1,1:nwflux)
    end where
  
    !call fix_onelimiter1(ixI^L,iL^L,wRCtmp,wRC)

  end subroutine MP5limiterR

  subroutine minmod(ixImin1,ixImax1,ixOmin1,ixOmax1,a,b,minm)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: a(ixImin1:ixImax1), b(ixImin1:ixImax1)
    double precision, intent(out):: minm(ixImin1:ixImax1)

    minm(ixOmin1:ixOmax1) = (sign(one,a(ixOmin1:ixOmax1))+sign(one,&
       b(ixOmin1:ixOmax1)))/2.0d0 * min(abs(a(ixOmin1:ixOmax1)),&
       abs(b(ixOmin1:ixOmax1)))

  end subroutine minmod

  subroutine median(ixImin1,ixImax1,ixOmin1,ixOmax1,a,b,c,med)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: a(ixImin1:ixImax1), b(ixImin1:ixImax1),&
        c(ixImin1:ixImax1)
    double precision, intent(out):: med(ixImin1:ixImax1)

    double precision             :: tmp1(ixImin1:ixImax1),&
       tmp2(ixImin1:ixImax1)

    tmp1(ixOmin1:ixOmax1) = b(ixOmin1:ixOmax1) - a(ixOmin1:ixOmax1)
    tmp2(ixOmin1:ixOmax1) = c(ixOmin1:ixOmax1) - a(ixOmin1:ixOmax1)

    med(ixOmin1:ixOmax1) = a(ixOmin1:ixOmax1) + (sign(one,&
       tmp1(ixOmin1:ixOmax1))+sign(one,tmp2(ixOmin1:ixOmax1)))/2.0d0 * &
       min(abs(tmp1(ixOmin1:ixOmax1)),abs(tmp2(ixOmin1:ixOmax1)))

  end subroutine median

  !> MP5 limiter from Suresh & Huynh 1997
  !> Following the convention of Mignone et al. 2010.
  !> Needs at least three ghost cells.  Set nghostcells=3.
  ! for one variable only: no fixing applied
  subroutine MP5limitervar(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC,wRC)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1),&
       wLC(ixImin1:ixImax1) 

    integer                         :: iLmmin1,iLmmax1, iLmmmin1,iLmmmax1,&
        iLpmin1,iLpmax1, iLppmin1,iLppmax1, iLpppmin1,iLpppmax1
    integer                         :: idmin1,idmax1, idpmin1,idpmax1,&
        idppmin1,idppmax1, idmmin1,idmmax1, iemin1,iemax1, iemmin1,iemmax1,&
        iepmin1,iepmax1, ieppmin1,ieppmax1
    double precision, dimension(ixImin1:ixImax1)  :: f, fmp, fmin, fmax, ful,&
        dm4, d, fmd, flc, flim
    double precision, dimension(ixImin1:ixImax1)  :: wRCtmp, wLCtmp
    double precision, dimension(ixImin1:ixImax1) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.0d0, alpha=4.0d0
    !double precision                :: alpha

    ! Variable alpha:
    !alpha = float(nstep)/courantpar - one

    ! Left side:
    ! range to process:
    !iLmin^D=ixmin^D-kr(idims,^D);iLmax^D=ixmax^D;

    ! HALL
    ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
    ! also add one ghost zone!
    !   {iL^L=iL^L^LADD1;}

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

    iLmmin1=iLmin1-kr(idims,1);iLmmax1=iLmax1-kr(idims,1);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmax1=iLmmax1-kr(idims,1);
    iLpmin1=iLmin1+kr(idims,1);iLpmax1=iLmax1+kr(idims,1);
    iLppmin1=iLpmin1+kr(idims,1);iLppmax1=iLpmax1+kr(idims,1);

    f(iLmin1:iLmax1) = 1.0d0/60.0d0 * (2.0d0* w(iLmmmin1:iLmmmax1) - 13.0d0* &
       w(iLmmin1:iLmmax1) + 47.0d0* w(iLmin1:iLmax1) + 27.0d0* &
       w(iLpmin1:iLpmax1) - 3.0d0*  w(iLppmin1:iLppmax1))

    ! get fmp and ful:
    a(iLmin1:iLmax1) = w(iLpmin1:iLpmax1)-w(iLmin1:iLmax1)
    b(iLmin1:iLmax1) = alpha*(w(iLmin1:iLmax1)-w(iLmmin1:iLmmax1))
    call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
    fmp(iLmin1:iLmax1) = w(iLmin1:iLmax1) + tmp(iLmin1:iLmax1)
    ful(iLmin1:iLmax1) = w(iLmin1:iLmax1) + b(iLmin1:iLmax1)

    ! get dm4:
    idmax1=iLmax1; idmin1=iLmin1-kr(idims,1);
    idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
    idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

    iemax1=idmax1+kr(idims,1); iemin1=idmin1;
    iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
    iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);

    d(iemin1:iemax1) = w(iepmin1:iepmax1)-2.0d0*w(iemin1:iemax1)+&
       w(iemmin1:iemmax1)

    a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1)-d(idpmin1:idpmax1)
    b(idmin1:idmax1) = 4.0d0*d(idpmin1:idpmax1)-d(idmin1:idmax1)
    call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
    a(idmin1:idmax1) = d(idmin1:idmax1)
    b(idmin1:idmax1) = d(idpmin1:idpmax1)
    call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
    call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
    dm4(idmin1:idmax1) = tmp3(idmin1:idmax1)

    ! get fmd:
    fmd(iLmin1:iLmax1) = (w(iLmin1:iLmax1)+&
       w(iLpmin1:iLpmax1))/2.0d0-dm4(iLmin1:iLmax1)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1) = half*(3.0d0*w(iLmin1:iLmax1) - w(iLmmin1:iLmmax1)) + &
       4.0d0/3.0d0*dm4(iLmmin1:iLmmax1)

    fmin(iLmin1:iLmax1) = max(min(w(iLmin1:iLmax1),w(iLpmin1:iLpmax1),&
       fmd(iLmin1:iLmax1)),min(w(iLmin1:iLmax1),ful(iLmin1:iLmax1),&
       flc(iLmin1:iLmax1)))

    fmax(iLmin1:iLmax1) = min(max(w(iLmin1:iLmax1),w(iLpmin1:iLpmax1),&
       fmd(iLmin1:iLmax1)),max(w(iLmin1:iLmax1),ful(iLmin1:iLmax1),&
       flc(iLmin1:iLmax1)))

    call median(ixImin1,ixImax1,iLmin1,iLmax1,fmin,f,fmax,tmp)
    flim(iLmin1:iLmax1) = tmp(iLmin1:iLmax1)

    ! check case
    where ((f(iLmin1:iLmax1)-w(iLmin1:iLmax1))*(f(iLmin1:iLmax1)-&
       fmp(iLmin1:iLmax1)) .le. eps)
       wLC(iLmin1:iLmax1) = f(iLmin1:iLmax1)
    elsewhere
       wLC(iLmin1:iLmax1) = flim(iLmin1:iLmax1)
    end where

    ! Right side:
    ! the interpolation from the right is obtained when the left-hand formula is applied to
    ! data mirrored about the interface.  
    ! thus substitute: 
    ! i-2 -> i+3
    ! i-1 -> i+2
    ! i   -> i+1
    ! i+1 -> i
    ! i+2 -> i-1

    iLpppmin1=iLppmin1+kr(idims,1);iLpppmax1=iLppmax1+kr(idims,1);

    f(iLmin1:iLmax1) = 1.0d0/60.0d0 * (2.0d0* w(iLpppmin1:iLpppmax1) - 13.0d0* &
       w(iLppmin1:iLppmax1) + 47.0d0* w(iLpmin1:iLpmax1) + 27.0d0* &
       w(iLmin1:iLmax1) - 3.0d0*  w(iLmmin1:iLmmax1))

    ! get fmp and ful:
    a(iLmin1:iLmax1) = w(iLmin1:iLmax1)-w(iLpmin1:iLpmax1)
    b(iLmin1:iLmax1) = alpha*(w(iLpmin1:iLpmax1)-w(iLppmin1:iLppmax1))
    call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
    fmp(iLmin1:iLmax1) = w(iLpmin1:iLpmax1) + tmp(iLmin1:iLmax1)
    ful(iLmin1:iLmax1) = w(iLpmin1:iLpmax1) + b(iLmin1:iLmax1)

    ! get dm4:
    idmax1=iLmax1+kr(idims,1); idmin1=iLmin1;
    idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
    idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

    iemax1=idmax1; iemin1=idmin1-kr(idims,1);
    iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
    iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);
    ieppmin1=iepmin1+kr(idims,1);ieppmax1=iepmax1+kr(idims,1);

    d(iemin1:iemax1) = w(iemin1:iemax1)-2.0d0*w(iepmin1:iepmax1)+&
       w(ieppmin1:ieppmax1)

    a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1)-d(idmmin1:idmmax1)
    b(idmin1:idmax1) = 4.0d0*d(idmmin1:idmmax1)-d(idmin1:idmax1)
    call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
    a(idmin1:idmax1) = d(idmin1:idmax1)
    b(idmin1:idmax1) = d(idmmin1:idmmax1)
    call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
    call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
    dm4(idmin1:idmax1) = tmp3(idmin1:idmax1)

    ! get fmd:
    fmd(iLmin1:iLmax1) = (w(iLmin1:iLmax1)+&
       w(iLpmin1:iLpmax1))/2.0d0-dm4(iLmin1:iLmax1)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1) = half*(3.0d0*w(iLpmin1:iLpmax1) - &
       w(iLppmin1:iLppmax1)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1)

    fmin(iLmin1:iLmax1) = max(min(w(iLpmin1:iLpmax1),w(iLmin1:iLmax1),&
       fmd(iLmin1:iLmax1)),min(w(iLpmin1:iLpmax1),ful(iLmin1:iLmax1),&
       flc(iLmin1:iLmax1)))

    fmax(iLmin1:iLmax1) = min(max(w(iLpmin1:iLpmax1),w(iLmin1:iLmax1),&
       fmd(iLmin1:iLmax1)),max(w(iLpmin1:iLpmax1),ful(iLmin1:iLmax1),&
       flc(iLmin1:iLmax1)))

    call median(ixImin1,ixImax1,iLmin1,iLmax1,fmin,f,fmax,flim)

    ! check case
    where ((f(iLmin1:iLmax1)-w(iLpmin1:iLpmax1))*(f(iLmin1:iLmax1)-&
       fmp(iLmin1:iLmax1))  .le. eps)
       wRC(iLmin1:iLmax1) = f(iLmin1:iLmax1)
    elsewhere
       wRC(iLmin1:iLmax1) = flim(iLmin1:iLmax1)
    end where

  end subroutine MP5limitervar

end module mod_mp5
