!> Hydrodynamics HLLC module
module mod_hd_hllc
  use mod_hd_phys

  implicit none
  private

  public :: hd_hllc_init

contains

  subroutine hd_hllc_init()
    use mod_physics_hllc

    phys_diffuse_hllcd => hd_diffuse_hllcd
    phys_get_lCD => hd_get_lCD
    phys_get_wCD => hd_get_wCD
  end subroutine hd_hllc_init

  subroutine hd_diffuse_hllcd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idim,wLC,wRC,fLC,fRC,patchf)

    ! when method is hllcd or hllcd1 then: 

    ! this subroutine is to enforce regions where we AVOID HLLC
    ! and use TVDLF instead: this is achieved by setting patchf to 4 in
    ! certain regions. An additional input parameter is nxdiffusehllc
    ! which sets the size of the fallback region.

    use mod_global_parameters

    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        intent(in)      :: wRC,wLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
       intent(in) :: fLC, fRC
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(inout)        :: patchf

    integer                                           :: ixOO1,ixOO2,TxOOmin1,&
       TxOOmin2,TxOOmax1,TxOOmax2

    ! In a user-controlled region around any point with flux sign change between
    ! left and right, ensure fallback to TVDLF
    do ixOO1= ixOmin1,ixOmax1
    do ixOO2= ixOmin2,ixOmax2
    
    TxOOmin1= max(ixOO1 - nxdiffusehllc*kr(idim,1), ixOmin1);
    TxOOmax1= min(ixOO1 + nxdiffusehllc*kr(idim,1), ixOmax1);
    
    
    TxOOmin2= max(ixOO2 - nxdiffusehllc*kr(idim,2), ixOmin2);
    TxOOmax2= min(ixOO2 + nxdiffusehllc*kr(idim,2), ixOmax2);
    
    if(abs(patchf(ixOO1,ixOO2)) == 1 .or. abs(patchf(ixOO1,ixOO2)) == 4)Then
       if(any(fRC(ixOO1,ixOO2,1:nwflux)*fLC(ixOO1,ixOO2,&
          1:nwflux)<-smalldouble))Then
          where(Abs(patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2))==1)
             patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2) = 4
          endwhere
       endif
    endif
    enddo
    enddo

  end subroutine hd_diffuse_hllcd

  subroutine hd_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idim,ixImin1,ixImin2,ixImax1,&
     ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2, whll,Fhll,lambdaCD,patchf)
    ! Calculate lambda at CD and set the patchf to know the orientation
    ! of the riemann fan and decide on the flux choice
    ! We also compute here the HLL flux and w value, for fallback strategy

    ! In this nul version, we simply compute nothing and ensure TVDLF fallback
    use mod_global_parameters

    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        intent(in)      :: wLC,wRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
        intent(in)  :: fLC,fRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(in)           :: cmax,cmin
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(inout)        :: patchf
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
        intent(out) :: Fhll,whll
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(out)          :: lambdaCD

    logical         , dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)     :: Cond_patchf
    double precision                       :: Epsilon
    integer                                :: iw,ix1,ix2

    !--------------------------------------------
    ! on entry, patch is preset to contain values from -2,1,2,4
    !      -2: take left flux, no computation here
    !      +2: take right flux, no computation here
    !      +4: take TVDLF flux, no computation here
    !       1: compute the characteristic speed for the CD

    Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(abs(patchf(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))==1)
    lambdaCD=0.d0

    do iw=1,nwflux
      where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
        !============= compute HLL flux ==============!
        Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)= (cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iw)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*fRC(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,iw) + cmin(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)*cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)*(wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iw)-wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)))/(cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
        !======== compute intermediate HLL state =======!
        whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = (cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iw)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*wLC(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,iw)+fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iw)-fRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw))/(cmax(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      end where
    end do

    ! deduce the characteristic speed at the CD
    where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
       lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=whll(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(idim))/whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          rho_)
    end where

    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
       if(Cond_patchf(ix1,ix2)) then
         ! double check whether obtained speed is in between min and max speeds given
         ! and identify in which part of the Riemann fan the time-axis is
         if(cmin(ix1,ix2) < zero .and. lambdaCD(ix1,ix2)>zero.and.lambdaCD(ix1,&
            ix2)<cmax(ix1,ix2)) then
            patchf(ix1,ix2) = -1
         else if(cmax(ix1,ix2) > zero .and. lambdaCD(ix1,&
            ix2) < zero.and.lambdaCD(ix1,ix2)>cmin(ix1,ix2)) then
            patchf(ix1,ix2) =  1
         else if(lambdaCD(ix1,ix2) >= cmax(ix1,ix2) .or. lambdaCD(ix1,&
            ix2) <= cmin(ix1,ix2)) then
            lambdaCD(ix1,ix2) = zero
            ! we will fall back to HLL flux case in this degeneracy
            patchf(ix1,ix2) =  3
         end if
       end if
    end do
    end do

    where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)== 3)
      Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=.false.
    end where

    ! handle the specific case where the time axis is exactly on the CD 
    if(any(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2).and.&
       lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==zero))then
      ! determine which sector (forward or backward) of the Riemann fan is smallest
      ! and select left or right flux accordingly
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        if(lambdaCD(ix1,ix2)==zero.and.Cond_patchf(ix1,ix2)) then
          if(-cmin(ix1,ix2)>=cmax(ix1,ix2)) then
            patchf(ix1,ix2) =  1
          else
            patchf(ix1,ix2) = -1
          end if
        end if
      end do
      end do
    end if

    ! eigenvalue lambda for contact is near zero: decrease noise by this trick
    if(flathllc)then
      Epsilon=1.0d-6
      where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2).and. &
         dabs(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))/max(cmax(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
         Epsilon)< Epsilon  .and. dabs(lambdaCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))/max(dabs(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
         Epsilon)< Epsilon)
        lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  zero
      end where
    end if

  end subroutine hd_get_lCD

  subroutine hd_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
     ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,f)
    ! compute the intermediate state U*
    ! only needed where patchf=-1/1

    ! reference Li S., JCP, 203, 2005, 344-357
    ! reference T. Miyoski, Kusano JCP, 2008, 2005.

    use mod_global_parameters
    use mod_dust

    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        intent(in)      :: wRC,wLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
        intent(in):: whll, Fhll
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(in)           :: lambdaCD
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(in)           :: cmax,cmin
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
        intent(in):: fRC,fLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
       intent(out):: f
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)      :: wCD,wSub
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux)  :: fSub
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)           :: vSub,cspeed,pCD
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
        intent(in)           :: patchf

    double precision :: csmls
    integer                                      :: n, iw, ix1,ix2

    !-------------- auxiliary Speed and array-------------!
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
       if(patchf(ix1,ix2)==1) then
         cspeed(ix1,ix2)=cmax(ix1,ix2)
         vSub(ix1,ix2)=wRC(ix1,ix2,mom(idim))/wRC(ix1,ix2,rho_)
         wSub(ix1,ix2,:)=wRC(ix1,ix2,:)
         fSub(ix1,ix2,:)=fRC(ix1,ix2,:)
       else if(patchf(ix1,ix2)==-1) then
         cspeed(ix1,ix2)=cmin(ix1,ix2)
         vSub(ix1,ix2)=wLC(ix1,ix2,mom(idim))/wLC(ix1,ix2,rho_)
         wSub(ix1,ix2,:)=wLC(ix1,ix2,:)
         fSub(ix1,ix2,:)=fLC(ix1,ix2,:)
       end if
    end do
    end do

    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
      if(abs(patchf(ix1,ix2))==1) then
        csmls=one/(cspeed(ix1,ix2)-lambdaCD(ix1,ix2))
        wCD(ix1,ix2,rho_) = wSub(ix1,ix2,rho_)*(cspeed(ix1,ix2)-vSub(ix1,&
           ix2))*csmls
        do n=1,hd_n_tracer
          iw = tracer(n)
          wCD(ix1,ix2,iw) = wSub(ix1,ix2,iw)*(cspeed(ix1,ix2)-vSub(ix1,&
             ix2))*csmls
        end do

        !------- Momentum ------!
        do iw=1, ndir
          if(iw /= idim)then
            ! eq. 21 22
            wCD(ix1,ix2,mom(iw))=(cspeed(ix1,ix2)*wSub(ix1,ix2,&
               mom(iw))-fSub(ix1,ix2,mom(iw)))*csmls
          else
            ! eq. 20
            wCD(ix1,ix2,mom(iw)) =  wCD(ix1,ix2,rho_) * lambdaCD(ix1,ix2)
          endif
        enddo
        if(hd_energy) then
          pCD(ix1,ix2) = wSub(ix1,ix2,rho_)*(cspeed(ix1,ix2)-vSub(ix1,&
             ix2))*(lambdaCD(ix1,ix2)-vSub(ix1,ix2))+fSub(ix1,ix2,&
             mom(idim))-wSub(ix1,ix2,mom(idim))*vSub(ix1,ix2)
          ! Eq 31
          wCD(ix1,ix2,e_) = (cspeed(ix1,ix2)*wSub(ix1,ix2,e_) -fSub(ix1,ix2,&
             e_)+lambdaCD(ix1,ix2)*pCD(ix1,ix2))*csmls
        end if
        !if(hd_dust) then
        !  do n=1,dust_n_species
        !    wCD(ix^D,dust_rho(n)) = wSub(ix^D,dust_rho(n))*(cspeed(ix^D)-vSub(ix^D))*csmls
        !    do iw=1,ndir
        !      if(iw /= idim)then
        !        ! eq. 21 22
        !        wCD(ix^D,dust_mom(iw,n))=wSub(ix^D,dust_mom(iw,n))*(cspeed(ix^D)-vSub(ix^D))*csmls
        !      else
        !        ! eq. 20
        !        wCD(ix^D,dust_mom(iw,n)) =  wCD(ix^D,dust_rho(n)) * lambdaCD(ix^D)
        !      endif
        !    end do
        !  end do
        !end if

        if(hd_dust) then
          do iw=1,nwflux
            if(iw>=dust_rho(1).and.iw<=dust_mom(ndir,dust_n_species)) then
              ! use HLL flux for dust
              f(ix1,ix2,iw)=Fhll(ix1,ix2,iw)
            else
              ! f_i=fsub+lambda (wCD-wSub)
              f(ix1,ix2,iw)=fsub(ix1,ix2,iw)+cspeed(ix1,ix2)*(wCD(ix1,ix2,&
                 iw)-wSub(ix1,ix2,iw))
            end if
          end do
        else
          do iw=1,nwflux
            ! f_i=fsub+lambda (wCD-wSub)
            f(ix1,ix2,iw)=fsub(ix1,ix2,iw)+cspeed(ix1,ix2)*(wCD(ix1,ix2,&
               iw)-wSub(ix1,ix2,iw))
          end do
        end if
      end if
    end do
    end do

  end subroutine hd_get_wCD

end module mod_hd_hllc
