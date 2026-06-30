module mod_mhd_hllc
  use mod_mhd_phys

  implicit none
  private

  public :: mhd_hllc_init

contains

  subroutine mhd_hllc_init()
    use mod_physics_hllc

    phys_diffuse_hllcd => mhd_diffuse_hllcd
    phys_get_lCD => mhd_get_lCD
    phys_get_wCD => mhd_get_wCD

  end subroutine mhd_hllc_init

  subroutine mhd_diffuse_hllcd(ixImin1,ixImax1,ixOmin1,ixOmax1,idim,wLC,wRC,&
     fLC,fRC,patchf)
  ! when method is hllcd or hllcd1 then: 
  ! this subroutine is to enforce regions where we AVOID HLLC
  ! and use TVDLF instead: this is achieved by setting patchf to 4 in
  ! certain regions. An additional input parameter is nxdiffusehllc
  ! which sets the size of the fallback region.
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixImin1,&
       ixImax1,ixOmin1,ixOmax1,idim
    double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)      :: wRC,&
       wLC
    double precision, dimension(ixImin1:ixImax1,1:nwflux),intent(in) :: fLC,&
        fRC
    integer, dimension(ixImin1:ixImax1), intent(inout) :: patchf
    
    integer :: ixOO1,TxOOmin1,TxOOmax1

    
    ! In a user-controlled region around any point with flux sign change between
    ! left and right, ensure fallback to TVDLF
    do ixOO1= ixOmin1,ixOmax1
      
      TxOOmin1= max(ixOO1 - nxdiffusehllc*kr(idim,1), ixOmin1);
      TxOOmax1= min(ixOO1 + nxdiffusehllc*kr(idim,1), ixOmax1);
      
      if(abs(patchf(ixOO1)) == 1 .or. abs(patchf(ixOO1)) == 4)Then
         if(any(fRC(ixOO1,1:nwflux)*fLC(ixOO1,1:nwflux)<-smalldouble))Then
           where(Abs(patchf(TxOOmin1:TxOOmax1))==1)
             patchf(TxOOmin1:TxOOmax1) = 4
           endwhere
         endif
      endif
    enddo
  
  end subroutine mhd_diffuse_hllcd

  subroutine mhd_get_lCD(wLC,wRC,fLC,fRC,cmin,cmax,idim,ixImin1,ixImax1,&
     ixOmin1,ixOmax1, whll,Fhll,lambdaCD,patchf)
  ! Calculate lambda at CD and set the patchf to know the orientation
  ! of the riemann fan and decide on the flux choice
  ! We also compute here the HLL flux and w value, for fallback strategy
    use mod_global_parameters

    integer, intent(in)                                      :: ixImin1,&
       ixImax1,ixOmin1,ixOmax1,idim
    double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)      :: wLC,&
       wRC
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(in)  :: fLC,&
       fRC
    double precision, dimension(ixImin1:ixImax1), intent(in)           :: cmax,&
       cmin
    integer         , dimension(ixImin1:ixImax1),&
        intent(inout)        :: patchf
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(out) :: Fhll,&
       whll
    double precision, dimension(ixImin1:ixImax1),&
        intent(out)          :: lambdaCD

    logical         , dimension(ixOmin1:ixOmax1)     :: Cond_patchf
    double precision                       :: Epsilon
    integer                                :: iw,ix1

    ! on entry, patch is preset to contain values from -2,1,2,4
    !      -2: take left flux, no computation here
    !      +2: take right flux, no computation here
    !      +4: take TVDLF flux, no computation here
    !       1: compute the characteristic speed for the CD

    Cond_patchf(ixOmin1:ixOmax1)=(abs(patchf(ixOmin1:ixOmax1))==1)
    lambdaCD=0.d0

    do iw=1,nwflux
      where(Cond_patchf(ixOmin1:ixOmax1))
      !============= compute HLL flux ==============!
      Fhll(ixOmin1:ixOmax1,iw)= (cmax(ixOmin1:ixOmax1)*fLC(ixOmin1:ixOmax1,&
         iw)-cmin(ixOmin1:ixOmax1)*fRC(ixOmin1:ixOmax1,&
         iw) + cmin(ixOmin1:ixOmax1)*cmax(ixOmin1:ixOmax1)*(wRC(&
         ixOmin1:ixOmax1,iw)-wLC(ixOmin1:ixOmax1,&
         iw)))/(cmax(ixOmin1:ixOmax1)-cmin(ixOmin1:ixOmax1))
      !======== compute intermediate HLL state =======!
      whll(ixOmin1:ixOmax1,iw) = (cmax(ixOmin1:ixOmax1)*wRC(ixOmin1:ixOmax1,&
         iw)-cmin(ixOmin1:ixOmax1)*wLC(ixOmin1:ixOmax1,iw)+fLC(ixOmin1:ixOmax1,&
         iw)-fRC(ixOmin1:ixOmax1,iw))/(cmax(ixOmin1:ixOmax1)-&
         cmin(ixOmin1:ixOmax1))
      end where
    end do

    ! deduce the characteristic speed at the CD
    where(Cond_patchf(ixOmin1:ixOmax1))
      lambdaCD(ixOmin1:ixOmax1)=whll(ixOmin1:ixOmax1,&
         mom(idim))/whll(ixOmin1:ixOmax1,rho_)
    end where

    do ix1=ixOmin1,ixOmax1
      if(Cond_patchf(ix1)) then
        ! double check whether obtained speed is in between min and max speeds given
        ! and identify in which part of the Riemann fan the time-axis is
        if(cmin(ix1)<zero.and.lambdaCD(ix1)>zero.and.lambdaCD(ix1)<cmax(ix1)) &
           then
          patchf(ix1) = -1
        else if(cmax(ix1)>zero.and.lambdaCD(ix1)<zero.and.&
           lambdaCD(ix1)>cmin(ix1)) then
          patchf(ix1) =  1
        else if(lambdaCD(ix1)>=cmax(ix1).or.lambdaCD(ix1) <= cmin(ix1)) then
          lambdaCD(ix1) = zero
          ! we will fall back to HLL flux case in this degeneracy
          patchf(ix1) =  3
        end if
      end if
    end do

    where(patchf(ixOmin1:ixOmax1)== 3)
      Cond_patchf(ixOmin1:ixOmax1)=.false.
    end where

    ! handle the specific case where the time axis is exactly on the CD 
    if(any(lambdaCD(ixOmin1:ixOmax1)==zero.and.&
       Cond_patchf(ixOmin1:ixOmax1)))then
      ! determine which sector (forward or backward) of the Riemann fan is smallest
      ! and select left or right flux accordingly
      do ix1=ixOmin1,ixOmax1
        if(lambdaCD(ix1)==zero.and.Cond_patchf(ix1)) then
          if(-cmin(ix1)>=cmax(ix1)) then
            patchf(ix1) =  1
          else
            patchf(ix1) = -1
          end if
        end if
      end do
    end if

    ! eigenvalue lambda for contact is near zero: decrease noise by this trick
    if(flathllc)then
      Epsilon=1.0d-6
      where(Cond_patchf(ixOmin1:ixOmax1).and. &
         dabs(lambdaCD(ixOmin1:ixOmax1))/max(cmax(ixOmin1:ixOmax1),&
         Epsilon)< Epsilon  .and. dabs(lambdaCD(ixOmin1:ixOmax1))/max(dabs(&
         cmin(ixOmin1:ixOmax1)),Epsilon)< Epsilon)
        lambdaCD(ixOmin1:ixOmax1) =  zero
      end where
    end if 

  end subroutine mhd_get_lCD

  subroutine mhd_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
     ixImin1,ixImax1,ixOmin1,ixOmax1,idim,f)
  ! compute the intermediate state U*
  ! only needed where patchf=-1/1
  
  ! reference Li S., JCP, 203, 2005, 344-357
  ! reference T. Miyoski, Kusano JCP, 2008, 2005.
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixImin1,&
       ixImax1,ixOmin1,ixOmax1,idim
    double precision, dimension(ixImin1:ixImax1,1:nw), intent(in)      :: wRC,&
       wLC
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(in):: whll,&
        Fhll
    double precision, dimension(ixImin1:ixImax1),&
        intent(in)           :: lambdaCD
    double precision, dimension(ixImin1:ixImax1), intent(in)           :: cmax,&
       cmin
    double precision, dimension(ixImin1:ixImax1,1:nwflux), intent(in):: fRC,&
       fLC
    double precision, dimension(ixImin1:ixImax1,1:nwflux),intent(out):: f
    double precision, dimension(ixImin1:ixImax1,1:nw)        :: wCD,wSub
    double precision, dimension(ixImin1:ixImax1,1:nwflux)    :: fSub
    double precision, dimension(ixImin1:ixImax1)             :: vSub,cspeed,&
       pCD,VdotBCD
    integer         , dimension(ixImin1:ixImax1),&
        intent(in)           :: patchf

    integer                                        :: n, iw, idir,ix1

    !-------------- auxiliary Speed and array-------------!
    do ix1=ixOmin1,ixOmax1
       if(patchf(ix1)==1) then
         cspeed(ix1)=cmax(ix1)
         vSub(ix1)=wRC(ix1,mom(idim))/wRC(ix1,rho_)
         wSub(ix1,:)=wRC(ix1,:)
         fSub(ix1,:)=fRC(ix1,:)
       else if(patchf(ix1)==-1) then
         cspeed(ix1)=cmin(ix1)
         vSub(ix1)=wLC(ix1,mom(idim))/wLC(ix1,rho_)
         wSub(ix1,:)=wLC(ix1,:)
         fSub(ix1,:)=fLC(ix1,:)
       end if
    end do
    
    do ix1=ixOmin1,ixOmax1
      if(abs(patchf(ix1))==1) then
        wCD(ix1,rho_) = wSub(ix1,rho_)*(cspeed(ix1)-vSub(ix1))/(cspeed(ix1)-&
           lambdaCD(ix1))
        do n=1,mhd_n_tracer
          iw = tracer(n)
          wCD(ix1,iw) = wSub(ix1,iw)*(cspeed(ix1)-vSub(ix1))/(cspeed(ix1)-&
             lambdaCD(ix1))
        end do
        !==== Magnetic field ====!
        do idir=1,ndir
          ! case from eq 31
          wCD(ix1,mag(idir)) = whll(ix1,mag(idir))
        end do
        !------- Momentum ------!
        do iw=1, ndir
          if(iw /= idim)then
            ! eq. 21 22
            wCD(ix1,mom(iw))=(cspeed(ix1)*wSub(ix1,mom(iw))-fSub(ix1,&
               mom(iw))-wCD(ix1,mag(idim))*wCD(ix1,&
               mag(iw)))/(cspeed(ix1)-lambdaCD(ix1))
          else
            ! eq. 20
            wCD(ix1,mom(iw)) =  wCD(ix1,rho_) * lambdaCD(ix1)
          endif
        enddo
        if(mhd_energy) then
          VdotBCD(ix1) = sum(whll(ix1,mom(:))*whll(ix1,mag(:)))/whll(ix1,rho_)
          ! Eq 17
          pCD(ix1)  = wsub(ix1,rho_)*(cspeed(ix1)-vSub(ix1))*(lambdaCD(ix1)-&
             vSub(ix1))+fSub(ix1,mom(idim))-wsub(ix1,&
             mom(idim))*vSub(ix1)+ wCD(ix1,mag(idim))**2
          ! Eq 31
          wCD(ix1,e_) = (cspeed(ix1)*wSub(ix1,e_) -fSub(ix1,&
             e_)+lambdaCD(ix1)*pCD(ix1)-VdotBCD(ix1)*wCD(ix1,&
             mag(idim)))/(cspeed(ix1)-lambdaCD(ix1))
        end if
      end if
    end do

    do iw=1,nwflux
     if(iw == mag(idim)) then
       f(ixOmin1:ixOmax1,iw)=zero
     else if(mhd_glm .and. iw == psi_) then
       f(ixOmin1:ixOmax1,iw)=zero
     else
       where(abs(patchf(ixOmin1:ixOmax1))==1)
         ! f_i=fsub+lambda (wCD-wSub)
         f(ixOmin1:ixOmax1,iw)=fsub(ixOmin1:ixOmax1,&
            iw)+cspeed(ixOmin1:ixOmax1)*(wCD(ixOmin1:ixOmax1,&
            iw)-wsub(ixOmin1:ixOmax1,iw))
       endwhere
     end if
    end do

  end subroutine mhd_get_wCD

end module mod_mhd_hllc
