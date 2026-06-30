module mod_ppm

  implicit none
  private

  public :: PPMlimiter
  public :: PPMlimitervar

contains

  subroutine PPMlimitervar(ixImin1,ixImax1,ixmin1,ixmax1,idims,q,qCT,qLC,qRC)

    ! references:
    ! Mignone et al 2005, ApJS 160, 199,
    ! Miller and Colella 2002, JCP 183, 26
    ! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
    ! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
    ! version : april 2009
    ! author: zakaria.meliani@wis.kuleuven.be

    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixmin1,ixmax1, idims
    double precision, intent(in)    :: q(ixImin1:ixImax1),qCT(ixImin1:ixImax1)

    double precision, intent(inout) :: qRC(ixGlo1:ixGhi1),qLC(ixGlo1:ixGhi1)

    double precision,dimension(ixGlo1:ixGhi1)  :: dqC,d2qC,ldq
    double precision,dimension(ixGlo1:ixGhi1)  :: qMin,qMax,tmp

    integer   :: lxCmin1,lxCmax1,lxRmin1,lxRmax1
    integer   :: ixLLmin1,ixLLmax1,ixLmin1,ixLmax1,ixOmin1,ixOmax1,ixRmin1,&
       ixRmax1,ixRRmin1,ixRRmax1
    integer   :: hxLmin1,hxLmax1,hxCmin1,hxCmax1,hxRmin1,hxRmax1
    integer   :: kxLLmin1,kxLLmax1,kxLmin1,kxLmax1,kxCmin1,kxCmax1,kxRmin1,&
       kxRmax1,kxRRmin1,kxRRmax1

    ixOmin1=ixmin1-kr(idims,1);ixOmax1=ixmax1+kr(idims,1); !ixO[ixMmin1-1,ixMmax1+1]
    ixLmin1=ixOmin1-kr(idims,1);ixLmax1=ixOmax1-kr(idims,1); !ixL[ixMmin1-2,ixMmax1]
    ixLLmin1=ixLmin1-kr(idims,1);ixLLmax1=ixLmax1-kr(idims,1); !ixLL[ixMmin1-3,ixMmax1-1]
    ixRmin1=ixOmin1+kr(idims,1);ixRmax1=ixOmax1+kr(idims,1); !ixR=[iMmin1,ixMmax+2]
    ixRRmin1=ixRmin1+kr(idims,1);ixRRmax1=ixRmax1+kr(idims,1); !ixRR=[iMmin1+1,ixMmax+3]

    hxCmin1=ixOmin1;hxCmax1=ixmax1; ! hxC = [ixMmin-1,ixMmax]
    hxLmin1=hxCmin1-kr(idims,1);hxLmax1=hxCmax1-kr(idims,1); !hxL = [ixMmin-2,ixMmax-1]
    hxRmin1=hxCmin1+kr(idims,1);hxRmax1=hxCmax1+kr(idims,1); !hxR = [ixMmin,ixMmax+1]

    kxCmin1=ixLLmin1; kxCmax1=ixRmax1; ! kxC=[iMmin1-3,ixMmax1+2]
    kxLmin1=kxCmin1-kr(idims,1);kxLmax1=kxCmax1-kr(idims,1); !kxL=[iMmin1-4,ixMmax1+1]
    kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1); !kxR=[iMmin1-2,ixMmax1+3]

    lxCmin1=ixLLmin1-kr(idims,1);lxCmax1=ixRRmax1;! ixC=[iMmin1-4,ixMmax1+3]
    lxRmin1=lxCmin1+kr(idims,1);lxRmax1=lxCmax1+kr(idims,1); !lxR=[iMmin1-3,ixMmax1+4]


    dqC(lxCmin1:lxCmax1)=q(lxRmin1:lxRmax1)-q(lxCmin1:lxCmax1)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26
    d2qC(kxCmin1:kxCmax1)=half*(q(kxRmin1:kxRmax1)-q(kxLmin1:kxLmax1))
    where(dqC(kxCmin1:kxCmax1)*dqC(kxLmin1:kxLmax1)>zero)
       ! Store the sign of d2qC in qMin
       qMin(kxCmin1:kxCmax1)= sign(one,d2qC(kxCmin1:kxCmax1))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26
       ldq(kxCmin1:kxCmax1)= qMin(kxCmin1:kxCmax1)*min(dabs(d2qC(&
          kxCmin1:kxCmax1)),2.0d0*dabs(dqC(kxLmin1:kxLmax1)),&
          2.0d0*dabs(dqC(kxCmin1:kxCmax1)))
    elsewhere
       ldq(kxCmin1:kxCmax1)=zero
    end where

    ! Eq. 66, Miller and Colella 2002, JCP 183, 26
    qLC(ixOmin1:ixOmax1)=qLC(ixOmin1:ixOmax1)+half*dqC(ixOmin1:ixOmax1)+&
       (ldq(ixOmin1:ixOmax1)-ldq(ixRmin1:ixRmax1))/6.0d0
    qRC(ixLmin1:ixLmax1)=qRC(ixLmin1:ixLmax1) &
       -(half*dqC(ixLmin1:ixLmax1)-(ldq(ixLmin1:ixLmax1)-&
       ldq(ixOmin1:ixOmax1))/6.0d0)

    ! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
    call extremaq(ixImin1,ixImax1,kxCmin1,kxCmax1,qCT,1,qMax,qMin)

    ! Eq. B8, page 217, Mignone et al 2005, ApJS
    qRC(ixLmin1:ixLmax1)=max(qMin(ixOmin1:ixOmax1),min(qMax(ixOmin1:ixOmax1),&
       qRC(ixLmin1:ixLmax1)))
    qLC(ixOmin1:ixOmax1)=max(qMin(ixOmin1:ixOmax1),min(qMax(ixOmin1:ixOmax1),&
       qLC(ixOmin1:ixOmax1)))

    ! Eq. B9, page 217, Mignone et al 2005, ApJS
    where((qRC(ixLmin1:ixLmax1)-qCT(ixOmin1:ixOmax1))*(qCT(ixOmin1:ixOmax1)-&
       qLC(ixOmin1:ixOmax1))<=zero)
       qRC(ixLmin1:ixLmax1)=qCT(ixOmin1:ixOmax1)
       qLC(ixOmin1:ixOmax1)=qCT(ixOmin1:ixOmax1)
    end where

    qMax(ixOmin1:ixOmax1)=(qLC(ixOmin1:ixOmax1)-&
       qRC(ixLmin1:ixLmax1))*(qCT(ixOmin1:ixOmax1)-(qLC(ixOmin1:ixOmax1)+&
       qRC(ixLmin1:ixLmax1))/2.0d0)
    qMin(ixOmin1:ixOmax1)=(qLC(ixOmin1:ixOmax1)-&
       qRC(ixLmin1:ixLmax1))**2.0d0/6.0d0
    tmp(ixLmin1:ixLmax1)=qRC(ixLmin1:ixLmax1)

    ! Eq. B10, page 218, Mignone et al 2005, ApJS
    where(qMax(hxRmin1:hxRmax1)>qMin(hxRmin1:hxRmax1))
       qRC(hxCmin1:hxCmax1)= 3.0d0*qCT(hxRmin1:hxRmax1)-&
          2.0d0*qLC(hxRmin1:hxRmax1)
    end where
    ! Eq. B11, page 218, Mignone et al 2005, ApJS
    where(qMax(hxCmin1:hxCmax1)<-qMin(hxCmin1:hxCmax1))
       qLC(hxCmin1:hxCmax1)= 3.0d0*qCT(hxCmin1:hxCmax1)-&
          2.0d0*tmp(hxLmin1:hxLmax1)
    end where

  end subroutine PPMlimitervar

  subroutine PPMlimiter(ixImin1,ixImax1,ixmin1,ixmax1,idims,w,wCT,wLC,wRC)

    ! references:
    ! Mignone et al 2005, ApJS 160, 199, 
    ! Miller and Colella 2002, JCP 183, 26 
    ! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
    ! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
    ! version : april 2009
    ! author: zakaria.meliani@wis.kuleuven.be

    use mod_global_parameters
    use mod_physics, only: phys_ppm_flatcd, phys_ppm_flatsh

    integer, intent(in)             :: ixImin1,ixImax1, ixmin1,ixmax1, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw),&
       wCT(ixImin1:ixImax1,1:nw)

    double precision, intent(inout) :: wRC(ixGlo1:ixGhi1,1:nw),&
       wLC(ixGlo1:ixGhi1,1:nw) 

    double precision,dimension(ixGlo1:ixGhi1,1:nwflux)  :: dwC,d2wC,ldw
    double precision,dimension(ixGlo1:ixGhi1,1:nwflux)  :: wMin,wMax,tmp
    double precision,dimension(ixGlo1:ixGhi1) :: aa, ab, ac, dv
    double precision,dimension(ixGlo1:ixGhi1,1:ndim) ::  exi

    integer   :: lxCmin1,lxCmax1,lxRmin1,lxRmax1
    integer   :: ixLLmin1,ixLLmax1,ixLmin1,ixLmax1,ixOmin1,ixOmax1,ixRmin1,&
       ixRmax1,ixRRmin1,ixRRmax1
    integer   :: hxLmin1,hxLmax1,hxCmin1,hxCmax1,hxRmin1,hxRmax1
    integer   :: kxLLmin1,kxLLmax1,kxLmin1,kxLmax1,kxCmin1,kxCmax1,kxRmin1,&
       kxRmax1,kxRRmin1,kxRRmax1
    integer   :: iw, idimss

    double precision, parameter :: betamin=0.75d0, betamax=0.85d0,Zmin=0.25d0,&
        Zmax=0.75d0,eta1=20.0d0,eta2=0.05d0,eps=0.01d0,kappa=0.1d0

    ixOmin1=ixmin1-kr(idims,1);ixOmax1=ixmax1+kr(idims,1); !ixO[ixMmin1-1,ixMmax1+1]
    ixLmin1=ixOmin1-kr(idims,1);ixLmax1=ixOmax1-kr(idims,1); !ixL[ixMmin1-2,ixMmax1]
    ixLLmin1=ixLmin1-kr(idims,1);ixLLmax1=ixLmax1-kr(idims,1); !ixLL[ixMmin1-3,ixMmax1-1]
    ixRmin1=ixOmin1+kr(idims,1);ixRmax1=ixOmax1+kr(idims,1); !ixR=[iMmin1,ixMmax+2]
    ixRRmin1=ixRmin1+kr(idims,1);ixRRmax1=ixRmax1+kr(idims,1); !ixRR=[iMmin1+1,ixMmax+3]

    hxCmin1=ixOmin1;hxCmax1=ixmax1; ! hxC = [ixMmin-1,ixMmax]
    hxLmin1=hxCmin1-kr(idims,1);hxLmax1=hxCmax1-kr(idims,1); !hxL = [ixMmin-2,ixMmax-1]
    hxRmin1=hxCmin1+kr(idims,1);hxRmax1=hxCmax1+kr(idims,1); !hxR = [ixMmin,ixMmax+1]

    kxCmin1=ixLLmin1; kxCmax1=ixRmax1; ! kxC=[iMmin1-3,ixMmax1+2]
    kxLmin1=kxCmin1-kr(idims,1);kxLmax1=kxCmax1-kr(idims,1); !kxL=[iMmin1-4,ixMmax1+1]
    kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1); !kxR=[iMmin1-2,ixMmax1+3]

    lxCmin1=ixLLmin1-kr(idims,1);lxCmax1=ixRRmax1;! ixC=[iMmin1-4,ixMmax1+3]
    lxRmin1=lxCmin1+kr(idims,1);lxRmax1=lxCmax1+kr(idims,1); !lxR=[iMmin1-3,ixMmax1+4]

    dwC(lxCmin1:lxCmax1,1:nwflux)=w(lxRmin1:lxRmax1,&
       1:nwflux)-w(lxCmin1:lxCmax1,1:nwflux)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
    d2wC(kxCmin1:kxCmax1,1:nwflux)=half*(w(kxRmin1:kxRmax1,&
       1:nwflux)-w(kxLmin1:kxLmax1,1:nwflux))
    where(dwC(kxCmin1:kxCmax1,1:nwflux)*dwC(kxLmin1:kxLmax1,1:nwflux)>zero)
       ! Store the sign of dwC in wMin
       wMin(kxCmin1:kxCmax1,1:nwflux)= sign(one,d2wC(kxCmin1:kxCmax1,&
          1:nwflux))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
       ldw(kxCmin1:kxCmax1,1:nwflux)= wMin(kxCmin1:kxCmax1,&
          1:nwflux)*min(dabs(d2wC(kxCmin1:kxCmax1,1:nwflux)),&
          2.0d0*dabs(dwC(kxLmin1:kxLmax1,1:nwflux)),&
          2.0d0*dabs(dwC(kxCmin1:kxCmax1,1:nwflux)))
    elsewhere
       ldw(kxCmin1:kxCmax1,1:nwflux)=zero
    endwhere

    ! Eq. 66,  Miller and Colella 2002, JCP 183, 26 
    wLC(ixOmin1:ixOmax1,1:nwflux)=wLC(ixOmin1:ixOmax1,&
       1:nwflux)+half*dwC(ixOmin1:ixOmax1,1:nwflux)+(ldw(ixOmin1:ixOmax1,&
       1:nwflux)-ldw(ixRmin1:ixRmax1,1:nwflux))/6.0d0

    wRC(ixLmin1:ixLmax1,1:nwflux)=wRC(ixLmin1:ixLmax1,&
       1:nwflux)-(half*dwC(ixLmin1:ixLmax1,1:nwflux)-(ldw(ixLmin1:ixLmax1,&
       1:nwflux)-ldw(ixOmin1:ixOmax1,1:nwflux))/6.0d0)

    ! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
    call extremaw(ixImin1,ixImax1,kxCmin1,kxCmax1,wCT,1,wMax,wMin)

    ! Eq. B8, page 217, Mignone et al 2005, ApJS
    wRC(ixLmin1:ixLmax1,1:nwflux)=max(wMin(ixOmin1:ixOmax1,1:nwflux),&
       min(wMax(ixOmin1:ixOmax1,1:nwflux),wRC(ixLmin1:ixLmax1,1:nwflux))) 
    wLC(ixOmin1:ixOmax1,1:nwflux)=max(wMin(ixOmin1:ixOmax1,1:nwflux),&
       min(wMax(ixOmin1:ixOmax1,1:nwflux),wLC(ixOmin1:ixOmax1,1:nwflux)))


    ! Eq. B9, page 217, Mignone et al 2005, ApJS
    where((wRC(ixLmin1:ixLmax1,1:nwflux)-wCT(ixOmin1:ixOmax1,&
       1:nwflux))*(wCT(ixOmin1:ixOmax1,1:nwflux)-wLC(ixOmin1:ixOmax1,&
       1:nwflux))<=zero)
       wRC(ixLmin1:ixLmax1,1:nwflux)=wCT(ixOmin1:ixOmax1,1:nwflux)
       wLC(ixOmin1:ixOmax1,1:nwflux)=wCT(ixOmin1:ixOmax1,1:nwflux)
    end where

    wMax(ixOmin1:ixOmax1,1:nwflux)=(wLC(ixOmin1:ixOmax1,&
       1:nwflux)-wRC(ixLmin1:ixLmax1,1:nwflux))*(wCT(ixOmin1:ixOmax1,&
       1:nwflux)-(wLC(ixOmin1:ixOmax1,1:nwflux)+wRC(ixLmin1:ixLmax1,&
       1:nwflux))/2.0d0)
    wMin(ixOmin1:ixOmax1,1:nwflux)=(wLC(ixOmin1:ixOmax1,&
       1:nwflux)-wRC(ixLmin1:ixLmax1,1:nwflux))**2.0d0/6.0d0
    tmp(ixLmin1:ixLmax1,1:nwflux)=wRC(ixLmin1:ixLmax1,1:nwflux)
    ! Eq. B10, page 218, Mignone et al 2005, ApJS
    where(wMax(hxRmin1:hxRmax1,1:nwflux)>wMin(hxRmin1:hxRmax1,1:nwflux))
       wRC(hxCmin1:hxCmax1,1:nwflux)= 3.0d0*wCT(hxRmin1:hxRmax1,&
          1:nwflux)-2.0d0*wLC(hxRmin1:hxRmax1,1:nwflux)
    endwhere
    ! Eq. B11, page 218, Mignone et al 2005, ApJS
    where(wMax(hxCmin1:hxCmax1,1:nwflux)<-wMin(hxCmin1:hxCmax1,1:nwflux))
       wLC(hxCmin1:hxCmax1,1:nwflux)= 3.0d0*wCT(hxCmin1:hxCmax1,&
          1:nwflux)-2.0d0*tmp(hxLmin1:hxLmax1,1:nwflux)
    endwhere

    ! flattening at the contact discontinuities
    if(flatcd)then
       call phys_ppm_flatcd(ixImin1,ixImax1,kxCmin1,kxCmax1,kxLmin1,kxLmax1,&
          kxRmin1,kxRmax1,wCT,d2wC,aa,ab)
       if(any(kappa*aa(kxCmin1:kxCmax1)>=ab(kxCmin1:kxCmax1)))then
          do iw=1,nwflux
             where(kappa*aa(kxCmin1:kxCmax1)>=ab(kxCmin1:kxCmax1).and. &
                dabs(dwC(kxCmin1:kxCmax1,iw))>smalldouble)
                wMax(kxCmin1:kxCmax1,iw) = wCT(kxRmin1:kxRmax1,&
                   iw)-2.0d0*wCT(kxCmin1:kxCmax1,iw)+wCT(kxLmin1:kxLmax1,iw)
             end where

             where(wMax(ixRmin1:ixRmax1,iw)*wMax(ixLmin1:ixLmax1,&
                iw)<zero .and.dabs(wCT(ixRmin1:ixRmax1,iw)-wCT(ixLmin1:ixLmax1,&
                iw))-eps*min(dabs(wCT(ixRmin1:ixRmax1,iw)),&
                dabs(wCT(ixLmin1:ixLmax1,&
                iw)))>zero .and. kappa*aa(ixOmin1:ixOmax1)>=&
                ab(ixOmin1:ixOmax1).and. dabs(dwC(ixOmin1:ixOmax1,&
                iw))>smalldouble)

                ac(ixOmin1:ixOmax1)=(wCT(ixLLmin1:ixLLmax1,&
                   iw)-wCT(ixRRmin1:ixRRmax1,iw)+4.0d0*dwC(ixOmin1:ixOmax1,&
                   iw))/(12.0d0*dwC(ixOmin1:ixOmax1,iw))
                wMin(ixOmin1:ixOmax1,iw)=max(zero,&
                   min(eta1*(ac(ixOmin1:ixOmax1)-eta2),one))
             elsewhere
                wMin(ixOmin1:ixOmax1,iw)=zero
             end where

             where(wMin(hxCmin1:hxCmax1,iw)>zero)
                wLC(hxCmin1:hxCmax1,iw) = wLC(hxCmin1:hxCmax1,&
                   iw)*(one-wMin(hxCmin1:hxCmax1,iw))+(wCT(hxCmin1:hxCmax1,&
                   iw)+half*ldw(hxCmin1:hxCmax1,iw))*wMin(hxCmin1:hxCmax1,iw)
             end where
             where(wMin(hxRmin1:hxRmax1,iw)>zero)
                wRC(hxCmin1:hxCmax1,iw) = wRC(hxCmin1:hxCmax1,&
                   iw)*(one-wMin(hxRmin1:hxRmax1,iw))+(wCT(hxRmin1:hxRmax1,&
                   iw)-half*ldw(hxRmin1:hxRmax1,iw))*wMin(hxRmin1:hxRmax1,iw)
             end where
          end do
       endif
    endif

    ! flattening at the shocks
    if(flatsh)then
       ! following MILLER and COLELLA 2002 JCP 183, 26
       kxCmin1=ixmin1-2; kxCmax1=ixmax1+2; ! kxC=[ixMmin1-2,ixMmax1+2]
       do idimss=1,ndim
          kxLmin1=kxCmin1-kr(idimss,1);kxLmax1=kxCmax1-kr(idimss,1); !kxL=[ixMmin1-3,ixMmax1+1]
          kxRmin1=kxCmin1+kr(idimss,1);kxRmax1=kxCmax1+kr(idimss,1); !kxR=[ixMmin1-1,ixMmax1+3]
          kxLLmin1=kxLmin1-kr(idimss,1);kxLLmax1=kxLmax1-kr(idimss,1); !kxLL=[ixMmin-4,ixMmax]
          kxRRmin1=kxRmin1+kr(idimss,1);kxRRmax1=kxRmax1+kr(idimss,1); !kxRR=[ixMmin,ixMmax+4]

          call phys_ppm_flatsh(ixImin1,ixImax1,kxCmin1,kxCmax1,kxLLmin1,&
             kxLLmax1,kxLmin1,kxLmax1,kxRmin1,kxRmax1,kxRRmin1,kxRRmax1,idimss,&
             wCT,aa,ab,dv)

          ! eq. B17, page 218, Mignone et al 2005, ApJS (had(Xi1))
          ac(kxCmin1:kxCmax1) = max(zero,min(one,&
             (betamax-aa(kxCmin1:kxCmax1))/(betamax-betamin)))
          ! eq. B18, page 218, Mignone et al 2005, ApJS (had(Xi1))
          ! recycling aa(ixL^S)
          where (dabs(dv(kxCmin1:kxCmax1))<smalldouble)
             aa(kxCmin1:kxCmax1) = max(ac(kxCmin1:kxCmax1), min(one,&
                (Zmax-ab(kxCmin1:kxCmax1))/(Zmax-Zmin)))
          elsewhere
             aa(kxCmin1:kxCmax1) = one
          endwhere

           call extremaa(ixImin1,ixImax1,ixOmin1,ixOmax1,aa,1,ab)
          
       enddo
       
       ! recycling wMax
       do iw=1,nwflux
          where(dabs(ab(ixOmin1:ixOmax1)-one)>smalldouble)
             wMax(ixOmin1:ixOmax1,iw) = (one-&
                ab(ixOmin1:ixOmax1))*wCT(ixOmin1:ixOmax1,iw)
          endwhere

          where(dabs(ab(hxCmin1:hxCmax1)-one)>smalldouble)
             wLC(hxCmin1:hxCmax1,iw) = ab(hxCmin1:hxCmax1)*wLC(hxCmin1:hxCmax1,&
                iw)+wMax(hxCmin1:hxCmax1,iw)
          endwhere

          where(dabs(ab(hxRmin1:hxRmax1)-one)>smalldouble)
             wRC(hxCmin1:hxCmax1,iw) = ab(hxRmin1:hxRmax1)*wRC(hxCmin1:hxCmax1,&
                iw)+wMax(hxRmin1:hxRmax1,iw)
          endwhere
       enddo
    endif

  end subroutine PPMlimiter

  subroutine extremaq(ixImin1,ixImax1,ixOmin1,ixOmax1,q,nshift,qMax,qMin)

    use mod_global_parameters

    integer,intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in) :: q(ixImin1:ixImax1)
    integer,intent(in)           :: nshift

    double precision, intent(out) :: qMax(ixImin1:ixImax1),&
       qMin(ixImin1:ixImax1)

    integer           :: ixsmin1,ixsmax1,ixsRmin1,ixsRmax1,ixsLmin1,ixsLmax1,&
       idims,jdims,kdims,ishift,i,j

    do ishift=1,nshift
      idims=1
      ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmax1=ixOmax1+ishift*kr(idims,1);
      ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmax1=ixOmax1-ishift*kr(idims,1);
      if (ishift==1) then
        qMax(ixOmin1:ixOmax1)=max(q(ixOmin1:ixOmax1),q(ixsRmin1:ixsRmax1),&
           q(ixsLmin1:ixsLmax1))
        qMin(ixOmin1:ixOmax1)=min(q(ixOmin1:ixOmax1),q(ixsRmin1:ixsRmax1),&
           q(ixsLmin1:ixsLmax1))
      else
        qMax(ixOmin1:ixOmax1)=max(qMax(ixOmin1:ixOmax1),q(ixsRmin1:ixsRmax1),&
           q(ixsLmin1:ixsLmax1))
        qMin(ixOmin1:ixOmax1)=min(qMin(ixOmin1:ixOmax1),q(ixsRmin1:ixsRmax1),&
           q(ixsLmin1:ixsLmax1))
      end if
      
      
    enddo

  end subroutine  extremaq

  subroutine extremaa(ixImin1,ixImax1,ixOmin1,ixOmax1,a,nshift,aMin)
    use mod_global_parameters

    integer,intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in) :: a(ixImin1:ixImax1)
    integer,intent(in)           :: nshift

    double precision, intent(out) :: aMin(ixImin1:ixImax1)

    integer          :: ixsmin1,ixsmax1,ixsRmin1,ixsRmax1,ixsLmin1,ixsLmax1,&
       idims,jdims,kdims,ishift,i,j

    do ishift=1,nshift
      idims=1
      ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmax1=ixOmax1+ishift*kr(idims,1);
      ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmax1=ixOmax1-ishift*kr(idims,1);
      aMin(ixOmin1:ixOmax1)=min(a(ixsRmin1:ixsRmax1),a(ixOmin1:ixOmax1),&
         a(ixsLmin1:ixsLmax1))
      
      
    end do

  end subroutine extremaa

  subroutine extremaw(ixImin1,ixImax1,ixOmin1,ixOmax1,w,nshift,wMax,wMin)
    use mod_global_parameters

    integer,intent(in)            :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw)
    integer,intent(in)            :: nshift

    double precision, intent(out) :: wMax(ixImin1:ixImax1,1:nwflux),&
       wMin(ixImin1:ixImax1,1:nwflux)

    integer          :: ixsmin1,ixsmax1,ixsRmin1,ixsRmax1,ixsLmin1,ixsLmax1,&
       idims,jdims,kdims,ishift,i,j

    do ishift=1,nshift
      idims=1
      ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmax1=ixOmax1+ishift*kr(idims,1);
      ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmax1=ixOmax1-ishift*kr(idims,1);
      if (ishift==1) then
        wMax(ixOmin1:ixOmax1,1:nwflux)= max(w(ixOmin1:ixOmax1,1:nwflux),&
           w(ixsRmin1:ixsRmax1,1:nwflux),w(ixsLmin1:ixsLmax1,1:nwflux))
        wMin(ixOmin1:ixOmax1,1:nwflux)= min(w(ixOmin1:ixOmax1,1:nwflux),&
           w(ixsRmin1:ixsRmax1,1:nwflux),w(ixsLmin1:ixsLmax1,1:nwflux))
      else
        wMax(ixOmin1:ixOmax1,1:nwflux)= max(wMax(ixOmin1:ixOmax1,1:nwflux),&
           w(ixsRmin1:ixsRmax1,1:nwflux),w(ixsLmin1:ixsLmax1,1:nwflux))
        wMin(ixOmin1:ixOmax1,1:nwflux)= min(wMin(ixOmin1:ixOmax1,1:nwflux),&
           w(ixsRmin1:ixsRmax1,1:nwflux),w(ixsLmin1:ixsLmax1,1:nwflux))
      end if
      
      
    enddo

  end subroutine  extremaw

end module mod_ppm
