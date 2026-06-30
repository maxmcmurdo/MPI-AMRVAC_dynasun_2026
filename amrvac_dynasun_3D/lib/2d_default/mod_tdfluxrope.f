module mod_tdfluxrope
  implicit none
  double precision :: d_TD99,L_TD99,R_TD99,a_TD99,q_TD99,Izero_TD99,Li_TD99
  double precision :: p_Bt_ratio

contains


  ! numercial recipe's subroutines
        function elle(phi,ak)
        implicit none
        double precision, intent(in) :: phi,ak
        double precision :: elle
        double precision :: cc,q,s
        integer, parameter :: dp = kind(1.d0)
        s=sin(phi)
        cc=cos(phi)**2
        q=(1.-s*ak)*(1.+s*ak)
        elle=s*(rf(cc,q,1.0_dp)-((s*ak)**2)*rd(cc,q,1.0_dp)/3.0)
        return
        end function elle
  
        function ellf(phi,ak)
        double precision, intent(in) :: phi,ak
        double precision :: ellf
        double precision :: s
        integer, parameter :: dp = kind(1.d0)
        s=sin(phi)
        ellf=s*rf(cos(phi)**2,(1.-s*ak)*(1.+s*ak),1.0_dp)
        return
        end function ellf
  
        function ellpi(phi,en,ak)
        double precision, intent(in) :: phi,en,ak
        double precision :: ellpi
        double precision :: cc,enss,q,s
        integer, parameter :: dp = kind(1.d0)
        s=sin(phi)
        enss=en*s*s
        cc=cos(phi)**2
        q=(1.-s*ak)*(1.+s*ak)
        ellpi=s*(rf(cc,q,1.0_dp)-enss*rj(cc,q,1.0_dp,1.0_dp+enss)/3.)
        return
        end function ellpi
  
        FUNCTION rc(x,y)
        double precision, intent(in) :: x,y
        double precision :: rc
        double precision :: ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD, C1,&
           C2,C3,C4
        PARAMETER (ERRTOL=.04,TINY=1.69e-38,SQRTNY=1.3e-19,BIG=3.E37,&
            TNBG=TINY*BIG,COMP1=2.236/SQRTNY,COMP2=TNBG*TNBG/25.,THIRD=1./3.,&
            C1=.3,C2=1./7.,C3=.375,C4=9./22.)
        double precision :: alamb,ave,s,w,xt,yt
        if(x.lt.0..or.y.eq.0..or.(x+abs(y)).lt.TINY.or.(x+ &
           abs(y)).gt.BIG.or.(y.lt.-COMP1.and.x.gt.0..and.x.lt.COMP2)) then 
        write(*,*)'invalid arguments in rc'
        read(*,*)
        endif
        if(y.gt.0.)then
          xt=x
          yt=y
          w=1.
        else
          xt=x-y
          yt=-y
          w=sqrt(x)/sqrt(xt)
        endif
  1     continue
          alamb=2.*sqrt(xt)*sqrt(yt)+yt
          xt=.25*(xt+alamb)
          yt=.25*(yt+alamb)
          ave=THIRD*(xt+yt+yt)
          s=(yt-ave)/ave
        if(abs(s).gt.ERRTOL)goto 1
        rc=w*(1.+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
        return
        end function rc
  
        FUNCTION rd(x,y,z)
        double precision, intent(in) :: x,y,z
        double precision :: rd
        double precision :: ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
        PARAMETER (ERRTOL=.05,TINY=1.e-25,BIG=4.5E21,C1=3./14.,C2=1./6.,&
            C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)
        double precision :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,&
           sqrty, sqrtz,sum,xt,yt,zt
        if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y, z).gt.BIG) then 
        write(*,*)'invalid arguments in rd'
        read(*,*)
        endif
        xt=x
        yt=y
        zt=z
        sum=0.
        fac=1.
  1     continue
          sqrtx=sqrt(xt)
          sqrty=sqrt(yt)
          sqrtz=sqrt(zt)
          alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          sum=sum+fac/(sqrtz*(zt+alamb))
          fac=.25*fac
          xt=.25*(xt+alamb)
          yt=.25*(yt+alamb)
          zt=.25*(zt+alamb)
          ave=.2*(xt+yt+3.*zt)
          delx=(ave-xt)/ave
          dely=(ave-yt)/ave
          delz=(ave-zt)/ave
        if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.*eb
        ee=ed+ec+ec
        rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3* &
           ec+delz*C4*ea)))/(ave*sqrt(ave))
        return
        end function rd
  
        FUNCTION rf(x,y,z)
        double precision, intent(in) :: x,y,z
        double precision :: rf
        double precision :: ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
        PARAMETER (ERRTOL=.08,TINY=1.5e-38,BIG=3.E37,THIRD=1./3., C1=1./24.,&
           C2=.1,C3=3./44.,C4=1./14.)
        double precision :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,&
           xt,yt,zt
        if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,&
            z).gt.BIG) then 
        write(*,*)'invalid arguments in rf'
        read(*,*)
        endif
        xt=x
        yt=y
        zt=z
  1     continue
          sqrtx=sqrt(xt)
          sqrty=sqrt(yt)
          sqrtz=sqrt(zt)
          alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          xt=.25*(xt+alamb)
          yt=.25*(yt+alamb)
          zt=.25*(zt+alamb)
          ave=THIRD*(xt+yt+zt)
          delx=(ave-xt)/ave
          dely=(ave-yt)/ave
          delz=(ave-zt)/ave
        if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        return
        end function rf
  
        FUNCTION rj(x,y,z,p)
        double precision, intent(in) :: x,y,z,p
        double precision :: rj
        double precision :: ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8
        PARAMETER (ERRTOL=.05,TINY=2.5e-13,BIG=9.E11,C1=3./14.,C2=1./3.,&
            C3=3./22.,C4=3./26.,C5=.75*C3,C6=1.5*C4,C7=.5*C2,C8=C3+C3)
        double precision :: a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,&
           ec,ed,ee, fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt
        if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY.or.max(x,y, z,&
           abs(p)).gt.BIG) then
        write(*,*)'invalid arguments in rj'
        read(*,*)
        endif
        sum=0.
        fac=1.
        if(p.gt.0.)then
          xt=x
          yt=y
          zt=z
          pt=p
        else
          xt=min(x,y,z)
          zt=max(x,y,z)
          yt=x+y+z-xt-zt
          a=1./(yt-p)
          b=a*(zt-yt)*(yt-xt)
          pt=yt+b
          rho=xt*zt/yt
          tau=p*pt/yt
          rcx=rc(rho,tau)
        endif
  1     continue
          sqrtx=sqrt(xt)
          sqrty=sqrt(yt)
          sqrtz=sqrt(zt)
          alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
          beta=pt*(pt+alamb)**2
          sum=sum+fac*rc(alpha,beta)
          fac=.25*fac
          xt=.25*(xt+alamb)
          yt=.25*(yt+alamb)
          zt=.25*(zt+alamb)
          pt=.25*(pt+alamb)
          ave=.2*(xt+yt+zt+pt+pt)
          delx=(ave-xt)/ave
          dely=(ave-yt)/ave
          delz=(ave-zt)/ave
          delp=(ave-pt)/ave
        if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1
        ea=delx*(dely+delz)+dely*delz
        eb=delx*dely*delz
        ec=delp**2
        ed=ea-3.*ec
        ee=eb+2.*delp*(ea-ec)
        rj=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))+ &
           delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
        if (p.le.0.) rj=a*(b*rj+3.*(rcx-rf(xt,yt,zt)))
        return
        end function rj
end module mod_tdfluxrope
