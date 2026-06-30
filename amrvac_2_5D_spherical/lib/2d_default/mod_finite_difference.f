!> Module with finite difference methods for fluxes
module mod_finite_difference

  implicit none
  private

  public :: fd
  public :: centdiff

contains

  subroutine fd(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,idimsmin,idimsmax,qtC,sCT,qt,snew,fC,fE,dxs,x)
    use mod_physics
    use mod_source, only: addsource2
    use mod_finite_volume, only: reconstruct_LR
    use mod_global_parameters
    use mod_usr_methods

    double precision, intent(in)                                     :: qdt,&
        qtC, qt, dxs(ndim)
    integer, intent(in)                                              :: &
       ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
        idimsmin,idimsmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        intent(in)            :: x

    type(state)                                                      :: sCT,&
        snew
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
       1:ndim), intent(out)  :: fC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)                    :: fE

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux)                      :: fCT
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)                          :: fm, fp, fmR, fpL, wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC,&
        wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cmaxC,&
        cminC
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)      :: Hspeed
    double precision, dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv
    logical :: transport, active
    integer :: idims, iw, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixmin1,ixmin2,&
       ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,hxOmax2, kxCmin1,kxCmin2,kxCmax1,&
       kxCmax2, kxRmin1,kxRmin2,kxRmax1,kxRmax2
    type(ct_velocity) :: vcts

    associate(wCT=>sCT%w,wnew=>snew%w)

    fC=0.d0
    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)

    b0i = 0 ! fd uses centered values in phys_get_flux  
    do idims= idimsmin,idimsmax

       !b0i=idims

       ! Get fluxes for the whole grid (mesh+nghostcells)
        ixmin1 = ixOmin1 - nghostcells * kr(idims,1)
         ixmin2 = ixOmin2 - nghostcells * kr(idims,2)
        ixmax1 = ixOmax1 + nghostcells * kr(idims,1)
         ixmax2 = ixOmax2 + nghostcells * kr(idims,2)

       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
         ixCmax2=ixOmax2+nghostcells-nghostcells*kr(idims,2)
         ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1)
         ixCmin2=hxOmin2-nghostcells+nghostcells*kr(idims,2);
         ixmax1=ixmax1+nghostcells-nghostcells*kr(idims,1)
         ixmax2=ixmax2+nghostcells-nghostcells*kr(idims,2)
         ixmin1=ixmin1-nghostcells+nghostcells*kr(idims,1)
         ixmin2=ixmin2-nghostcells+nghostcells*kr(idims,2);
         kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
         kxCmax2=ixImax2-kr(idims,2);
         kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
         kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2);
         ! wRp and wLp are defined at the same locations, and will correspond to
         ! the left and right reconstructed values at a cell face. Their indexing
         ! is similar to cell-centered values, but in direction idims they are
         ! shifted half a cell towards the 'lower' direction.
         wRp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nw)=wprim(kxRmin1:kxRmax1,&
            kxRmin2:kxRmax2,1:nw)
         wLp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nw)=wprim(kxCmin1:kxCmax1,&
            kxCmin2:kxCmax2,1:nw)
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2;
       end if

       call phys_get_flux(wCT,wprim,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,&
          ixmin2,ixmax1,ixmax2,idims,fCT)

       do iw=iwstart,nwflux
          ! Lax-Friedrich splitting:
          fp(ixmin1:ixmax1,ixmin2:ixmax2,iw) = half * (fCT(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw) + tvdlfeps * cmax_global * wCT(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw))
          fm(ixmin1:ixmax1,ixmin2:ixmax2,iw) = half * (fCT(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw) - tvdlfeps * cmax_global * wCT(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw))
       end do ! iw loop

       ! now do the reconstruction of fp and fm:
       call reconstructL(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
          ixCmax1,ixCmax2,idims,fp,fpL)
       call reconstructR(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
          ixCmax1,ixCmax2,idims,fm,fmR)

       fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iwstart:nwflux,&
          idims) = fpL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
          iwstart:nwflux) + fmR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
          iwstart:nwflux)

       if(stagger_grid) then
         ! apply limited reconstruction for left and right status at cell interfaces
         call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
            ixCmax1,ixCmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wprim,wLC,&
            wRC,wLp,wRp,x,dxs(idims))
         if(H_correction) then
           call phys_get_H_speed(wprim,x,ixImin1,ixImin2,ixImax1,ixImax2,&
              ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims,Hspeed)
         end if
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,Hspeed,cmaxC,cminC)
         call phys_get_ct_velocity(vcts,wLp,wRp,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,cminC)
       end if

    end do !idims loop
    !b0i=0

    if(stagger_grid) call phys_update_faces(ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,qdt,wprim,fC,fE,sCT,snew,vcts)

    if(slab_uniform) then
      dxinv=-qdt/dxs
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
        do iw=iwstart,nwflux
          fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
             idims) = dxinv(idims) * fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
             idims)
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw)+(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
        end do ! iw loop
      end do ! Next idims
    else
      inv_volume=1.d0/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
        do iw=iwstart,nwflux
          fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)=-qdt*fC(ixImin1:ixImax1,&
             ixImin2:ixImax2,iw,idims)*block%surfaceC(ixImin1:ixImax1,&
             ixImin2:ixImax2,idims)
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw)+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idims))*inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        end do ! iw loop
      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)

    if(stagger_grid) call phys_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       snew)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImax1,&
          ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,'fd')
    endif

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
       .false.,active,wprim)

    if(phys_solve_eaux.and.levmin==levmax) then
      ! synchronize internal energy for uniform grid
      call phys_energy_synchro(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,wnew,x)
    endif
    end associate

  end subroutine fd

  subroutine reconstructL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
     iLmax2,idims,w,wLC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision, intent(out)   :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw) 

    double precision                :: ldw(ixImin1:ixImax1,ixImin2:ixImax2),&
        dwC(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: jxRmin1,jxRmin2,jxRmax1,jxRmax2,&
        ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2,&
        kxCmin1,kxCmin2,kxCmax1,kxCmax2, iw
    double precision                :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wLC)
    case (limiter_weno5)
       call WENO5limiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wLC,1)
    case (limiter_weno5nm)
       call WENO5NMlimiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
          iLmax1,iLmax2,idims,w,wLC,1)
    case (limiter_wenoz5)
       call WENO5limiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wLC,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
          iLmax1,iLmax2,idims,w,wLC,2)
    case (limiter_wenozp5)
       call WENO5limiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wLC,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
          iLmax1,iLmax2,idims,w,wLC,3)
    case default 

       kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
       kxCmax2=ixImax2-kr(idims,2);

       wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iwstart:nwflux) = w(kxCmin1:kxCmax1,&
          kxCmin2:kxCmax2,iwstart:nwflux)

       jxRmin1=iLmin1+kr(idims,1);jxRmin2=iLmin2+kr(idims,2)
       jxRmax1=iLmax1+kr(idims,1);jxRmax2=iLmax2+kr(idims,2);

       ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=iLmin1-kr(idims,1)
       ixCmin2=iLmin2-kr(idims,2);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

       do iw=iwstart,nwflux
          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)

           if(need_global_a2max) then 
             a2max=a2max_global(idims)
           else
             select case(idims)
             case(1)
               a2max=schmid_rad1
             
             case(2)
               a2max=schmid_rad2
             
             case default
               call mpistop("idims is wrong in mod_limiter")
             end select
           end if

          call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
             ixCmax1,ixCmax2,idims,type_limiter(block%level),ldw=ldw,&
             a2max=a2max)

          wLC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wLC(iLmin1:iLmax1,iLmin2:iLmax2,&
             iw)+half*ldw(iLmin1:iLmax1,iLmin2:iLmax2)
       end do
    end select

  end subroutine reconstructL

  subroutine reconstructR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
     iLmax2,idims,w,wRC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision, intent(out)   :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw) 

    double precision                :: rdw(ixImin1:ixImax1,ixImin2:ixImax2),&
        dwC(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: jxRmin1,jxRmin2,jxRmax1,jxRmax2,&
        ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2,&
        kxCmin1,kxCmin2,kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,kxRmax2, iw
    double precision                :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wRC)
    case (limiter_weno5)
       call WENO5limiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wRC,1)
    case (limiter_weno5nm)
       call WENO5NMlimiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
          iLmax1,iLmax2,idims,w,wRC,1)
    case (limiter_wenoz5)
       call WENO5limiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wRC,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
          iLmax1,iLmax2,idims,w,wRC,2)
    case (limiter_wenozp5)
       call WENO5limiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
          iLmax2,idims,w,wRC,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
          iLmax1,iLmax2,idims,w,wRC,3)
    case default 

       kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
       kxCmax2=ixImax2-kr(idims,2);
       kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
       kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2);

       wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iwstart:nwflux)=w(kxRmin1:kxRmax1,&
          kxRmin2:kxRmax2,iwstart:nwflux)

       jxRmin1=iLmin1+kr(idims,1);jxRmin2=iLmin2+kr(idims,2)
       jxRmax1=iLmax1+kr(idims,1);jxRmax2=iLmax2+kr(idims,2);
       ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=iLmin1-kr(idims,1)
       ixCmin2=iLmin2-kr(idims,2);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

       do iw=iwstart,nwflux
          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)

           if(need_global_a2max) then
             a2max=a2max_global(idims)
           else
             select case(idims)
             case(1)
               a2max=schmid_rad1
             
             case(2)
               a2max=schmid_rad2
             
             case default
               call mpistop("idims is wrong in mod_limiter")
             end select
           end if

          call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
             ixCmax1,ixCmax2,idims,type_limiter(block%level),rdw=rdw,&
             a2max=a2max)

          wRC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wRC(iLmin1:iLmax1,iLmin2:iLmax2,&
             iw)-half*rdw(jxRmin1:jxRmax1,jxRmin2:jxRmax2)
       end do
    end select

  end subroutine reconstructR

  subroutine centdiff(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idimsmin,idimsmax,qtC,sCT,qt,s,fC,fE,dxs,x)

    ! Advance the flow variables from global_time to global_time+qdt within ixO^L by
    ! fourth order centered differencing in space 
    ! for the dw/dt+dF_i(w)/dx_i=S type equation.
    ! wCT contains the time centered variables at time qtC for flux and source.
    ! w is the old value at qt on input and the new value at qt+qdt on output.
    use mod_physics
    use mod_finite_volume, only: reconstruct_LR
    use mod_global_parameters
    use mod_source, only: addsource2
    use mod_usr_methods
    use mod_variables

    integer, intent(in) :: method
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idimsmin,idimsmax
    double precision, intent(in) :: qdt, qtC, qt, dxs(ndim)
    type(state)      :: sCT, s
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,1:ndim)
    double precision :: fE(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndim:3)

    double precision :: v(ixImin1:ixImax1,ixImin2:ixImax2,ndim),&
        f(ixImin1:ixImax1,ixImin2:ixImax2, nwflux)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wprim,&
        wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: vLC,&
        cmaxLC, cmaxRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:number_species)      :: cmaxC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:number_species)      :: cminC
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)      :: Hspeed
    double precision, dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)      :: inv_volume

    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,&
       jxCmax2, hxCmin1,hxCmin2,hxCmax1,hxCmax2, kxCmin1,kxCmin2,kxCmax1,&
       kxCmax2, kkxCmin1,kkxCmin2,kkxCmax1,kkxCmax2, kkxRmin1,kkxRmin2,&
       kkxRmax1,kkxRmax2
    type(ct_velocity) :: vcts
    logical :: transport, new_cmax, patchw(ixImin1:ixImax1,ixImin2:ixImax2),&
        active

    associate(wCT=>sCT%w,w=>s%w)
    ! two extra layers are needed in each direction for which fluxes are added.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idims= idimsmin,idimsmax
       ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
       ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
    end do

    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       then
       call mpistop("Error in centdiff: Non-conforming input limits")
    end if

    fC=0.d0
    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)

    ! get fluxes
    do idims= idimsmin,idimsmax
       b0i=idims

       ixmin1=ixOmin1-2*kr(idims,1);ixmin2=ixOmin2-2*kr(idims,2)
       ixmax1=ixOmax1+2*kr(idims,1);ixmax2=ixOmax2+2*kr(idims,2); 
       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
         ixCmax2=ixOmax2+nghostcells-nghostcells*kr(idims,2);
         ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1)
         ixCmin2=hxOmin2-nghostcells+nghostcells*kr(idims,2);
         ixmax1=ixmax1+nghostcells-nghostcells*kr(idims,1)
         ixmax2=ixmax2+nghostcells-nghostcells*kr(idims,2);
         ixmin1=ixmin1-nghostcells+nghostcells*kr(idims,1)
         ixmin2=ixmin2-nghostcells+nghostcells*kr(idims,2);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2;
       end if
       hxCmin1=ixCmin1-kr(idims,1);hxCmin2=ixCmin2-kr(idims,2)
       hxCmax1=ixCmax1-kr(idims,1);hxCmax2=ixCmax2-kr(idims,2); 
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2); 
       kxCmin1=ixCmin1+2*kr(idims,1);kxCmin2=ixCmin2+2*kr(idims,2)
       kxCmax1=ixCmax1+2*kr(idims,1);kxCmax2=ixCmax2+2*kr(idims,2); 

       kkxCmin1=ixImin1;kkxCmin2=ixImin2; kkxCmax1=ixImax1-kr(idims,1)
       kkxCmax2=ixImax2-kr(idims,2);
       kkxRmin1=kkxCmin1+kr(idims,1);kkxRmin2=kkxCmin2+kr(idims,2)
       kkxRmax1=kkxCmax1+kr(idims,1);kkxRmax2=kkxCmax2+kr(idims,2);
       wRp(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,&
          1:nwflux)=wprim(kkxRmin1:kkxRmax1,kkxRmin2:kkxRmax2,1:nwflux)
       wLp(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,&
          1:nwflux)=wprim(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,1:nwflux)

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
          ixCmax1,ixCmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wprim,wLC,wRC,&
          wLp,wRp,x,dxs(idims))

       if(stagger_grid) then
         if(H_correction) then
           call phys_get_H_speed(wprim,x,ixImin1,ixImin2,ixImax1,ixImax2,&
              ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims,Hspeed)
         end if
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,Hspeed,cmaxC,cminC)
         call phys_get_ct_velocity(vcts,wLp,wRp,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,cminC)
       end if

       ! Calculate velocities from upwinded values
       call phys_get_cmax(wLC,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
          ixCmax1,ixCmax2,idims,cmaxLC)
       call phys_get_cmax(wRC,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
          ixCmax1,ixCmax2,idims,cmaxRC)
       ! now take the maximum of left and right states
       vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=max(cmaxRC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2),cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

       call phys_get_flux(wCT,wprim,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,&
          ixmin2,ixmax1,ixmax2,idims,f)

       ! Center flux to interface
       if(method==fs_cd) then
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iwstart:nwflux,&
             idims)=half*(f(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iwstart:nwflux)+f(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             iwstart:nwflux))
       else
          ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iwstart:nwflux,&
             idims)=(-f(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
             iwstart:nwflux)+7.0d0*(f(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             iwstart:nwflux) + f(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iwstart:nwflux))-f(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
             iwstart:nwflux))/12.0d0

          do iw=iwstart,nwflux
             ! add rempel dissipative flux, only second order version for now
             fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2,iw,idims)-tvdlfeps*half*vLC(ixCmin1:ixCmax1,&
                ixCmin2:ixCmax2) *(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
          end do
       end if

    end do       !next idims
    b0i=0

    if(stagger_grid) call phys_update_faces(ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,qdt,wprim,fC,fE,sCT,s,vcts)

    if(slab_uniform) then
      dxinv=-qdt/dxs
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
        do iw=iwstart,nwflux
          fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
             idims)=dxinv(idims)*fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)
          ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)-f_(i-1))+f_(i-2)]/12
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw)+(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
        end do    !next iw
      end do ! Next idims
    else
      inv_volume=1.d0/block%dvolume
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
        do iw=iwstart,nwflux
          fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
             idims)=-qdt*block%surfaceC(ixImin1:ixImax1,ixImin2:ixImax2,&
             idims)*fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw)+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idims))*inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        end do    !next iw
      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)

    if(stagger_grid) call phys_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       s)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,w,x,ixImin1,ixImin2,ixImax1,&
          ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,'centdiff')
    endif

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,w,x,&
       .false.,active,wprim)

    if(phys_solve_eaux.and.levmin==levmax) then
      ! synchronize internal energy for uniform grid
      call phys_energy_synchro(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)
    endif
    end associate
  end subroutine centdiff

end module mod_finite_difference
