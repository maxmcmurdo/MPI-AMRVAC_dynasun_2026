!> Subroutines for TVD-MUSCL schemes
module mod_tvd

  implicit none
  private

  public :: tvdlimit
  public :: tvdlimit2
  public :: entropyfix

contains

  subroutine tvdlimit(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idimmin,idimmax,s,qt,snew,fC,dxs,x)
    use mod_global_parameters

    integer, intent(in) :: method
    double precision, intent(in) :: qdt, qt, dxs(ndim)
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idimmin,idimmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,nw) :: w, wnew
    type(state) :: s, snew
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,1:ndim)

    integer :: idims, ixICmin1,ixICmin2,ixICmax1,ixICmax2, jxICmin1,jxICmin2,&
       jxICmax1,jxICmax2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,nw) :: wR, wL

    associate(w=>s%w,wnew=>snew%w)
    do idims= idimmin,idimmax
       ixICmax1=ixOmax1+kr(idims,1);ixICmax2=ixOmax2+kr(idims,2)
       ixICmin1=ixOmin1-2*kr(idims,1);ixICmin2=ixOmin2-2*kr(idims,2);
       wL(ixICmin1:ixICmax1,ixICmin2:ixICmax2,1:nw)=w(ixICmin1:ixICmax1,&
          ixICmin2:ixICmax2,1:nw)
       jxICmin1=ixICmin1+kr(idims,1);jxICmin2=ixICmin2+kr(idims,2)
       jxICmax1=ixICmax1+kr(idims,1);jxICmax2=ixICmax2+kr(idims,2);
       wR(ixICmin1:ixICmax1,ixICmin2:ixICmax2,1:nw)=w(jxICmin1:jxICmax1,&
          jxICmin2:jxICmax2,1:nw)
       call tvdlimit2(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixICmin1,&
          ixICmin2,ixICmax1,ixICmax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims,wL,&
          wR,wnew,x,fC,dxs)
    end do
    end associate
  end subroutine tvdlimit

  subroutine tvdlimit2(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixICmin1,&
     ixICmin2,ixICmax1,ixICmax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims,wL,wR,&
     wnew,x,fC,dxs)

    ! Limit the flow variables in wnew according to typetvd. 
    ! wroeC is based on wL and wR.
    ! If method=fs_tvd an extra adtdx**2*jumpC is added to phiC for 2nd order
    ! accuracy in time.

    use mod_global_parameters
    use mod_physics_roe

    integer, intent(in) :: method
    double precision, intent(in) :: qdt, dxs(ndim)
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixICmin1,ixICmin2,&
       ixICmax1,ixICmax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idims
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw) :: wL, wR
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision :: wnew(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,1:ndim)

    double precision:: workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nworkroe)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw) :: wroeC
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: phiC, rphiC,&
        jumpC, adtdxC, smallaC
    double precision :: dxinv(1:ndim)
    integer :: hxOmin1,hxOmin2,hxOmax1,hxOmax2, ixCmin1,ixCmin2,ixCmax1,&
       ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, jxICmin1,jxICmin2,jxICmax1,&
       jxICmax2, iw, il
    !-----------------------------------------------------------------------------

    hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
    hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
    ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2; 

    jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
    jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);
    jxICmin1=ixICmin1+kr(idims,1);jxICmin2=ixICmin2+kr(idims,2)
    jxICmax1=ixICmax1+kr(idims,1);jxICmax2=ixICmax2+kr(idims,2);

    call phys_average(wL,wR,x,ixICmin1,ixICmin2,ixICmax1,ixICmax2,idims,wroeC,&
       workroe)

    dxinv=qdt/dxs

    ! A loop on characteristic variables to calculate the dissipative flux phiC.
    do il=1,nwflux
       !Calculate the jump in the il-th characteristic variable: L(wroe)*dw
       call phys_get_eigenjump(wL,wR,wroeC,x,ixICmin1,ixICmin2,ixICmax1,&
          ixICmax2,il,idims,smallaC,adtdxC,jumpC,workroe)

       ! Normalize the eigenvalue "a" (and its limit "smalla" if needed):
       if (slab_uniform) then
          adtdxC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=adtdxC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2)*dxinv(idims)
          if (typeentropy(il)=='harten' .or. &
             typeentropy(il)=='powell')smallaC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2)=smallaC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2)*dxinv(idims)
       else
          adtdxC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=adtdxC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2)*qdt*block%surfaceC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2,idims)*2.0d0/(block%dvolume(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2)+block%dvolume(jxICmin1:jxICmax1,&
             jxICmin2:jxICmax2))
          if (typeentropy(il)=='harten' .or. &
             typeentropy(il)=='powell')smallaC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2)=smallaC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2)*qdt*block%surfaceC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2,idims)*2.0d0/(block%dvolume(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2)+block%dvolume(jxICmin1:jxICmax1,&
             jxICmin2:jxICmax2))
       endif

       ! Calculate the flux limiter function phi
       call getphi(method,jumpC,adtdxC,smallaC,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixICmin1,ixICmin2,ixICmax1,ixICmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
          il,idims,phiC)

       !Add R(iw,il)*phiC(il) to each variable iw in wnew
       do iw=1,nwflux
          call phys_rtimes(phiC,wroeC,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,il,&
             idims,rphiC,workroe)

          if (slab_uniform) then
             rphiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=rphiC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)*half
             fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2,iw,idims)+rphiC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+rphiC(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2)-rphiC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)
          else
             rphiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=rphiC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)*quarter* (block%dvolume(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)+block%dvolume(jxCmin1:jxCmax1,&
                jxCmin2:jxCmax2))
             fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2,iw,idims)+rphiC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+(rphiC(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2)-rphiC(hxOmin1:hxOmax1,&
                hxOmin2:hxOmax2)) /block%dvolume(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2)
          endif
       end do  !iw
    end do     !il

  end subroutine tvdlimit2

  subroutine getphi(method,jumpC,adtdxC,smallaC,ixImin1,ixImin2,ixImax1,&
     ixImax2,ixICmin1,ixICmin2,ixICmax1,ixICmax2,ixCmin1,ixCmin2,ixCmax1,&
     ixCmax2,il,idims,phiC)

    ! Calculate the dissipative flux from jumpC=L*dw and adtdx=eigenvalue*dt/dx.
    ! Add Lax-Wendroff type correction if method=fs_tvd.
    ! Limit according to method and typetvd.
    use mod_limiter
    use mod_global_parameters

    integer, intent(in) :: method
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixICmin1,ixICmin2,&
       ixICmax1,ixICmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, il, idims
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: jumpC, adtdxC,&
        smallaC, phiC

    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: ljumpC, tmp
    integer :: jxCmin1,jxCmin2,jxCmax1,jxCmax2, ixmin1,ixmin2,ixmax1,ixmax2,&
        hxmin1,hxmin2,hxmax1,hxmax2, typelimiter
    !-----------------------------------------------------------------------------

    typelimiter=type_limiter(block%level)
    if(method==fs_tvdmu)then
       ! In the MUSCL scheme phi=|a|*jump, apply entropy fix to it
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=abs(adtdxC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
       else
          where(abs(adtdxC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))>=smallaC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
             phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=abs(adtdxC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
          elsewhere
             phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=half*(smallaC(&
                ixCmin1:ixCmax1,ixCmin2:ixCmax2)+adtdxC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)**2/smallaC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
          endwhere
       endif
       ! That's all for the MUSCL scheme
       return
    endif

    if(method==fs_tvd)then
       !Entropy fix to |a|-a**2
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=abs(adtdxC(&
             ixICmin1:ixICmax1,ixICmin2:ixICmax2))-adtdxC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2)**2
       else
          where(abs(adtdxC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2))>=smallaC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2))
             phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=abs(adtdxC(&
                ixICmin1:ixICmax1,ixICmin2:ixICmax2))-adtdxC(ixICmin1:ixICmax1,&
                ixICmin2:ixICmax2)**2
          elsewhere
             phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=half*smallaC(&
                ixICmin1:ixICmax1,ixICmin2:ixICmax2)+&
                (half/smallaC(ixICmin1:ixICmax1,&
                ixICmin2:ixICmax2)-one)*adtdxC(ixICmin1:ixICmax1,&
                ixICmin2:ixICmax2)**2
          endwhere
       endif
    else
       !Entropy fix to |a|
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=abs(adtdxC(&
             ixICmin1:ixICmax1,ixICmin2:ixICmax2))
       else
          where(abs(adtdxC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2))>=smallaC(ixICmin1:ixICmax1,&
             ixICmin2:ixICmax2))
             phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=abs(adtdxC(&
                ixICmin1:ixICmax1,ixICmin2:ixICmax2))
          elsewhere
             phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=half*smallaC(&
                ixICmin1:ixICmax1,ixICmin2:ixICmax2)+&
                half/smallaC(ixICmin1:ixICmax1,&
                ixICmin2:ixICmax2)*adtdxC(ixICmin1:ixICmax1,&
                ixICmin2:ixICmax2)**2
          endwhere
       endif
    endif

    jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
    jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);
    hxmin1=ixICmin1;hxmin2=ixICmin2; hxmax1=ixICmax1-kr(idims,1)
    hxmax2=ixICmax2-kr(idims,2);
    ixmin1=hxmin1+kr(idims,1);ixmin2=hxmin2+kr(idims,2)
    ixmax1=hxmax1+kr(idims,1);ixmax2=hxmax2+kr(idims,2);

    if (.not. limiter_symmetric(typelimiter)) then
       call mpistop("TVD only supports symmetric limiters")
    end if

    select case(typetvd)
    case('roe')
       call dwlimiter2(jumpC,ixImin1,ixImin2,ixImax1,ixImax2,ixICmin1,ixICmin2,&
          ixICmax1,ixICmax2,idims,typelimiter,ldw=ljumpC)
       where(adtdxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)<=0)
          phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(jumpC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-ljumpC(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
       elsewhere
          phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(jumpC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
       end where
       !extra (a*lambda)**2*delta
       if(method==fs_tvd)phiC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2)+adtdxC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2)**2*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
    case('sweby')
       !Sweby eqs.4.11-4.15, but no 0.5 ?!
       phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=phiC(ixICmin1:ixICmax1,&
          ixICmin2:ixICmax2)*jumpC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)
       call dwlimiter2(phiC,ixImin1,ixImin2,ixImax1,ixImax2,ixICmin1,ixICmin2,&
          ixICmax1,ixICmax2,idims,typelimiter,ldw=ljumpC)
       where(adtdxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)<=0)
          phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-ljumpC(jxCmin1:jxCmax1,jxCmin2:jxCmax2)
       elsewhere
          phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-ljumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
       end where
       !extra (a*lambda)**2*delta
       if(method==fs_tvd)phiC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2)=phiC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2)+adtdxC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2)**2*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
    case('yee')
       !eq.3.51 with correction
       call dwlimiter2(jumpC,ixImin1,ixImin2,ixImax1,ixImax2,ixICmin1,ixICmin2,&
          ixICmax1,ixICmax2,idims,typelimiter,ldw=ljumpC)

       !Use phiC as 0.5*(|nu|-nu**2) eq.3.45e for tvd otherwise 0.5*|nu|
       phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=half*phiC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2)
       !gamma*lambda eq.3.51d, use tmp to store agdtdxC
       where(abs(jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>smalldouble)
          tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=adtdxC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)+phiC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(ljumpC(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2)-ljumpC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))/jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
       elsewhere
          tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=adtdxC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)
       end where

       !eq.3.51a
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-phiC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(ljumpC(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2)+ljumpC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))+abs(tmp(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
       else
          where(abs(tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>=&
             smallaC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
             phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-phiC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)*(ljumpC(jxCmin1:jxCmax1,&
                jxCmin2:jxCmax2)+ljumpC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2))+abs(tmp(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2))*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
          elsewhere
             phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-phiC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)*(ljumpC(jxCmin1:jxCmax1,&
                jxCmin2:jxCmax2)+ljumpC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2))+(half*smallaC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)+half/smallaC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)*tmp(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)**2)*jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
          endwhere
       endif
    case('harten')
       !See Ryu, section 2.3
       !Use phiC as 0.5*(|nu|-nu**2)*jumpC eq.3.45b,e
       phiC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)=half*phiC(ixICmin1:ixICmax1,&
          ixICmin2:ixICmax2)*jumpC(ixICmin1:ixICmax1,ixICmin2:ixICmax2)
       call dwlimiter2(phiC,ixImin1,ixImin2,ixImax1,ixImax2,ixICmin1,ixICmin2,&
          ixICmax1,ixICmax2,idims,typelimiter,ldw=ljumpC)

       !gamma*lambda eq.3.45d, use tmp as agdtdxC
       where(abs(jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>smalldouble)
          tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=adtdxC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)+(ljumpC(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2)-ljumpC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))/jumpC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
       elsewhere
          tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=adtdxC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)
       end where
       !eq.3.45a with correction
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-ljumpC(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2)-ljumpC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)+jumpC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*abs(tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
       else
          where(abs(tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>=&
             smallaC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
             phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-ljumpC(jxCmin1:jxCmax1,&
                jxCmin2:jxCmax2)-ljumpC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)+jumpC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)*abs(tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
          elsewhere
             phiC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-ljumpC(jxCmin1:jxCmax1,&
                jxCmin2:jxCmax2)-ljumpC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)+jumpC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)*(half*smallaC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)+half/smallaC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2)*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2)
          endwhere
       endif
       !extra -(a*lambda)**2*delta
    case default
       call mpistop("Error in TVDLimit: Unknown TVD type")
    end select

  end subroutine getphi

  subroutine entropyfix(ixmin1,ixmin2,ixmax1,ixmax2,il,aL,aR,a,smalla)

    ! Apply entropyfix based on typeentropy(il),aL,aR, and a
    ! Calculate "smalla" (Harten,Powell) or modify "a" (ratio)

    use mod_global_parameters

    integer, intent(in) :: ixmin1,ixmin2,ixmax1,ixmax2, il
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: aL, aR, a,&
        smalla
    !-----------------------------------------------------------------------------

    select case(typeentropy(il))
    case('harten')
       smalla(ixmin1:ixmax1,ixmin2:ixmax2)=max(zero,a(ixmin1:ixmax1,&
          ixmin2:ixmax2)-aL(ixmin1:ixmax1,ixmin2:ixmax2),aR(ixmin1:ixmax1,&
          ixmin2:ixmax2)-a(ixmin1:ixmax1,ixmin2:ixmax2))
    case('powell')
       smalla(ixmin1:ixmax1,ixmin2:ixmax2)=max(zero,two*(aR(ixmin1:ixmax1,&
          ixmin2:ixmax2)-aL(ixmin1:ixmax1,ixmin2:ixmax2)))
       !!case('ratio')
       !!   where(aL(ix^S)<zero .and. aR(ix^S)>zero)&
       !!      a(ix^S)=a(ix^S)-2*aR(ix^S)*aL(ix^S)/(aR(ix^S)-aL(ix^S))
    case('yee')
       ! This has been done in geteigenjump already
    case('nul')
       ! No entropyfix is applied
    case default
       call mpistop("No such type of entropy fix")
    end select

  end subroutine entropyfix

end module mod_tvd
