!> Module with finite difference methods for fluxes
module mod_finite_difference

  implicit none
  private

  public :: fd
  public :: centdiff

contains

  subroutine fd(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimsmin,idimsmax,qtC,sCT,&
     qt,snew,fC,fE,dxs,x)
    use mod_physics
    use mod_source, only: addsource2
    use mod_finite_volume, only: reconstruct_LR
    use mod_global_parameters
    use mod_usr_methods

    double precision, intent(in)                                     :: qdt,&
        qtC, qt, dxs(ndim)
    integer, intent(in)                                              :: &
       ixImin1,ixImax1, ixOmin1,ixOmax1, idimsmin,idimsmax
    double precision, dimension(ixImin1:ixImax1,1:ndim),&
        intent(in)            :: x

    type(state)                                                      :: sCT,&
        snew
    double precision, dimension(ixImin1:ixImax1,1:nwflux,1:ndim),&
        intent(out)  :: fC
    double precision, dimension(ixImin1:ixImax1,&
       7-2*ndim:3)                    :: fE

    double precision, dimension(ixImin1:ixImax1,&
       1:nwflux)                      :: fCT
    double precision, dimension(ixImin1:ixImax1,&
       1:nw)                          :: fm, fp, fmR, fpL, wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,1:nw) :: wLp, wRp
    double precision, dimension(ixImin1:ixImax1)      :: cmaxC, cminC
    double precision, dimension(ixImin1:ixImax1)      :: Hspeed
    double precision, dimension(ixOmin1:ixOmax1)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv
    logical :: transport, active
    integer :: idims, iw, ixCmin1,ixCmax1, ixmin1,ixmax1, hxOmin1,hxOmax1,&
        kxCmin1,kxCmax1, kxRmin1,kxRmax1
    type(ct_velocity) :: vcts

    associate(wCT=>sCT%w,wnew=>snew%w)

    fC=0.d0
    wprim=wCT
    call phys_to_primitive(ixImin1,ixImax1,ixImin1,ixImax1,wprim,x)

    b0i = 0 ! fd uses centered values in phys_get_flux  
    do idims= idimsmin,idimsmax

       !b0i=idims

       ! Get fluxes for the whole grid (mesh+nghostcells)
        ixmin1 = ixOmin1 - nghostcells * kr(idims,1)
        ixmax1 = ixOmax1 + nghostcells * kr(idims,1)

       hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
         ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1);
         ixmax1=ixmax1+nghostcells-nghostcells*kr(idims,1)
         ixmin1=ixmin1-nghostcells+nghostcells*kr(idims,1);
         kxCmin1=ixImin1; kxCmax1=ixImax1-kr(idims,1);
         kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1);
         ! wRp and wLp are defined at the same locations, and will correspond to
         ! the left and right reconstructed values at a cell face. Their indexing
         ! is similar to cell-centered values, but in direction idims they are
         ! shifted half a cell towards the 'lower' direction.
         wRp(kxCmin1:kxCmax1,1:nw)=wprim(kxRmin1:kxRmax1,1:nw)
         wLp(kxCmin1:kxCmax1,1:nw)=wprim(kxCmin1:kxCmax1,1:nw)
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax1=ixOmax1; ixCmin1=hxOmin1;
       end if

       call phys_get_flux(wCT,wprim,x,ixGlo1,ixGhi1,ixmin1,ixmax1,idims,fCT)

       do iw=iwstart,nwflux
          ! Lax-Friedrich splitting:
          fp(ixmin1:ixmax1,iw) = half * (fCT(ixmin1:ixmax1,&
             iw) + tvdlfeps * cmax_global * wCT(ixmin1:ixmax1,iw))
          fm(ixmin1:ixmax1,iw) = half * (fCT(ixmin1:ixmax1,&
             iw) - tvdlfeps * cmax_global * wCT(ixmin1:ixmax1,iw))
       end do ! iw loop

       ! now do the reconstruction of fp and fm:
       call reconstructL(ixImin1,ixImax1,ixCmin1,ixCmax1,idims,fp,fpL)
       call reconstructR(ixImin1,ixImax1,ixCmin1,ixCmax1,idims,fm,fmR)

       fC(ixCmin1:ixCmax1,iwstart:nwflux,idims) = fpL(ixCmin1:ixCmax1,&
          iwstart:nwflux) + fmR(ixCmin1:ixCmax1,iwstart:nwflux)

       if(stagger_grid) then
         ! apply limited reconstruction for left and right status at cell interfaces
         call reconstruct_LR(ixImin1,ixImax1,ixCmin1,ixCmax1,ixCmin1,ixCmax1,&
            idims,wprim,wLC,wRC,wLp,wRp,x,dxs(idims))
         if(H_correction) then
           call phys_get_H_speed(wprim,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,&
              Hspeed)
         end if
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImax1,ixCmin1,&
            ixCmax1,idims,Hspeed,cmaxC,cminC)
         call phys_get_ct_velocity(vcts,wLp,wRp,ixImin1,ixImax1,ixCmin1,&
            ixCmax1,idims,cmaxC,cminC)
       end if

    end do !idims loop
    !b0i=0

    if(stagger_grid) call phys_update_faces(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,&
       qdt,wprim,fC,fE,sCT,snew,vcts)

    if(slab_uniform) then
      dxinv=-qdt/dxs
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
        do iw=iwstart,nwflux
          fC(ixImin1:ixImax1,iw,idims) = dxinv(idims) * fC(ixImin1:ixImax1,iw,&
             idims)
          wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,&
             iw)+(fC(ixOmin1:ixOmax1,iw,idims)-fC(hxOmin1:hxOmax1,iw,idims))
        end do ! iw loop
      end do ! Next idims
    else
      inv_volume=1.d0/block%dvolume(ixOmin1:ixOmax1)
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
        do iw=iwstart,nwflux
          fC(ixImin1:ixImax1,iw,idims)=-qdt*fC(ixImin1:ixImax1,iw,&
             idims)*block%surfaceC(ixImin1:ixImax1,idims)
          wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,&
             iw)+ (fC(ixOmin1:ixOmax1,iw,idims)-fC(hxOmin1:hxOmax1,iw,&
             idims))*inv_volume(ixOmin1:ixOmax1)
        end do ! iw loop
      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,wCT,wnew,x)

    if(stagger_grid) call phys_face_to_center(ixOmin1,ixOmax1,snew)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImax1,ixOmin1,&
          ixOmax1,'fd')
    endif

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), ixImin1,ixImax1,&
       ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,wnew,x,.false.,active,wprim)

    if(phys_solve_eaux.and.levmin==levmax) then
      ! synchronize internal energy for uniform grid
      call phys_energy_synchro(ixImin1,ixImax1,ixOmin1,ixOmax1,wnew,x)
    endif
    end associate

  end subroutine fd

  subroutine reconstructL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)

    double precision, intent(out)   :: wLC(ixImin1:ixImax1,1:nw) 

    double precision                :: ldw(ixImin1:ixImax1),&
        dwC(ixImin1:ixImax1)
    integer                         :: jxRmin1,jxRmax1, ixCmin1,ixCmax1,&
        jxCmin1,jxCmax1, kxCmin1,kxCmax1, iw
    double precision                :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC)
    case (limiter_weno5)
       call WENO5limiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC,1)
    case (limiter_weno5nm)
       call WENO5NMlimiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC,1)
    case (limiter_wenoz5)
       call WENO5limiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC,2)
    case (limiter_wenozp5)
       call WENO5limiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC,3)
    case default 

       kxCmin1=ixImin1; kxCmax1=ixImax1-kr(idims,1);

       wLC(kxCmin1:kxCmax1,iwstart:nwflux) = w(kxCmin1:kxCmax1,iwstart:nwflux)

       jxRmin1=iLmin1+kr(idims,1);jxRmax1=iLmax1+kr(idims,1);

       ixCmax1=jxRmax1; ixCmin1=iLmin1-kr(idims,1);
       jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);

       do iw=iwstart,nwflux
          dwC(ixCmin1:ixCmax1)=w(jxCmin1:jxCmax1,iw)-w(ixCmin1:ixCmax1,iw)

           if(need_global_a2max) then 
             a2max=a2max_global(idims)
           else
             select case(idims)
             case(1)
               a2max=schmid_rad1
             
             
             case default
               call mpistop("idims is wrong in mod_limiter")
             end select
           end if

          call dwlimiter2(dwC,ixImin1,ixImax1,ixCmin1,ixCmax1,idims,&
             type_limiter(block%level),ldw=ldw,a2max=a2max)

          wLC(iLmin1:iLmax1,iw)=wLC(iLmin1:iLmax1,iw)+half*ldw(iLmin1:iLmax1)
       end do
    end select

  end subroutine reconstructL

  subroutine reconstructR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)

    double precision, intent(out)   :: wRC(ixImin1:ixImax1,1:nw) 

    double precision                :: rdw(ixImin1:ixImax1),&
        dwC(ixImin1:ixImax1)
    integer                         :: jxRmin1,jxRmax1, ixCmin1,ixCmax1,&
        jxCmin1,jxCmax1, kxCmin1,kxCmax1, kxRmin1,kxRmax1, iw
    double precision                :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC)
    case (limiter_weno5)
       call WENO5limiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC,1)
    case (limiter_weno5nm)
       call WENO5NMlimiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC,1)
    case (limiter_wenoz5)
       call WENO5limiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC,2)
    case (limiter_wenozp5)
       call WENO5limiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC,3)
    case default 

       kxCmin1=ixImin1; kxCmax1=ixImax1-kr(idims,1);
       kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1);

       wRC(kxCmin1:kxCmax1,iwstart:nwflux)=w(kxRmin1:kxRmax1,iwstart:nwflux)

       jxRmin1=iLmin1+kr(idims,1);jxRmax1=iLmax1+kr(idims,1);
       ixCmax1=jxRmax1; ixCmin1=iLmin1-kr(idims,1);
       jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);

       do iw=iwstart,nwflux
          dwC(ixCmin1:ixCmax1)=w(jxCmin1:jxCmax1,iw)-w(ixCmin1:ixCmax1,iw)

           if(need_global_a2max) then
             a2max=a2max_global(idims)
           else
             select case(idims)
             case(1)
               a2max=schmid_rad1
             
             
             case default
               call mpistop("idims is wrong in mod_limiter")
             end select
           end if

          call dwlimiter2(dwC,ixImin1,ixImax1,ixCmin1,ixCmax1,idims,&
             type_limiter(block%level),rdw=rdw,a2max=a2max)

          wRC(iLmin1:iLmax1,iw)=wRC(iLmin1:iLmax1,&
             iw)-half*rdw(jxRmin1:jxRmax1)
       end do
    end select

  end subroutine reconstructR

  subroutine centdiff(method,qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimsmin,&
     idimsmax,qtC,sCT,qt,s,fC,fE,dxs,x)

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
    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, idimsmin,idimsmax
    double precision, intent(in) :: qdt, qtC, qt, dxs(ndim)
    type(state)      :: sCT, s
    double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
    double precision :: fC(ixImin1:ixImax1,1:nwflux,1:ndim)
    double precision :: fE(ixImin1:ixImax1,7-2*ndim:3)

    double precision :: v(ixImin1:ixImax1,ndim), f(ixImin1:ixImax1, nwflux)

    double precision, dimension(ixImin1:ixImax1,1:nw) :: wprim, wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,1:nw) :: wLp, wRp
    double precision, dimension(ixImin1:ixImax1)      :: vLC, cmaxLC, cmaxRC
    double precision, dimension(ixImin1:ixImax1,&
       1:number_species)      :: cmaxC
    double precision, dimension(ixImin1:ixImax1,&
       1:number_species)      :: cminC
    double precision, dimension(ixImin1:ixImax1)      :: Hspeed
    double precision, dimension(ixOmin1:ixOmax1)      :: inv_volume

    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ixmin1,ixmax1, hxOmin1,hxOmax1, ixCmin1,ixCmax1,&
        jxCmin1,jxCmax1, hxCmin1,hxCmax1, kxCmin1,kxCmax1, kkxCmin1,kkxCmax1,&
        kkxRmin1,kkxRmax1
    type(ct_velocity) :: vcts
    logical :: transport, new_cmax, patchw(ixImin1:ixImax1), active

    associate(wCT=>sCT%w,w=>s%w)
    ! two extra layers are needed in each direction for which fluxes are added.
    ixmin1=ixOmin1;ixmax1=ixOmax1;
    do idims= idimsmin,idimsmax
       ixmin1=ixmin1-2*kr(idims,1);ixmax1=ixmax1+2*kr(idims,1);
    end do

    if (ixImin1>ixmin1.or.ixImax1<ixmax1) then
       call mpistop("Error in centdiff: Non-conforming input limits")
    end if

    fC=0.d0
    wprim=wCT
    call phys_to_primitive(ixImin1,ixImax1,ixImin1,ixImax1,wprim,x)

    ! get fluxes
    do idims= idimsmin,idimsmax
       b0i=idims

       ixmin1=ixOmin1-2*kr(idims,1);ixmax1=ixOmax1+2*kr(idims,1); 
       hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1);
         ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1);
         ixmax1=ixmax1+nghostcells-nghostcells*kr(idims,1);
         ixmin1=ixmin1-nghostcells+nghostcells*kr(idims,1);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax1=ixOmax1; ixCmin1=hxOmin1;
       end if
       hxCmin1=ixCmin1-kr(idims,1);hxCmax1=ixCmax1-kr(idims,1); 
       jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1); 
       kxCmin1=ixCmin1+2*kr(idims,1);kxCmax1=ixCmax1+2*kr(idims,1); 

       kkxCmin1=ixImin1; kkxCmax1=ixImax1-kr(idims,1);
       kkxRmin1=kkxCmin1+kr(idims,1);kkxRmax1=kkxCmax1+kr(idims,1);
       wRp(kkxCmin1:kkxCmax1,1:nwflux)=wprim(kkxRmin1:kkxRmax1,1:nwflux)
       wLp(kkxCmin1:kkxCmax1,1:nwflux)=wprim(kkxCmin1:kkxCmax1,1:nwflux)

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(ixImin1,ixImax1,ixCmin1,ixCmax1,ixCmin1,ixCmax1,&
          idims,wprim,wLC,wRC,wLp,wRp,x,dxs(idims))

       if(stagger_grid) then
         if(H_correction) then
           call phys_get_H_speed(wprim,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,&
              Hspeed)
         end if
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImax1,ixCmin1,&
            ixCmax1,idims,Hspeed,cmaxC,cminC)
         call phys_get_ct_velocity(vcts,wLp,wRp,ixImin1,ixImax1,ixCmin1,&
            ixCmax1,idims,cmaxC,cminC)
       end if

       ! Calculate velocities from upwinded values
       call phys_get_cmax(wLC,x,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,cmaxLC)
       call phys_get_cmax(wRC,x,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,cmaxRC)
       ! now take the maximum of left and right states
       vLC(ixCmin1:ixCmax1)=max(cmaxRC(ixCmin1:ixCmax1),&
          cmaxLC(ixCmin1:ixCmax1))

       call phys_get_flux(wCT,wprim,x,ixImin1,ixImax1,ixmin1,ixmax1,idims,f)

       ! Center flux to interface
       if(method==fs_cd) then
          fC(ixCmin1:ixCmax1,iwstart:nwflux,idims)=half*(f(ixCmin1:ixCmax1,&
             iwstart:nwflux)+f(jxCmin1:jxCmax1,iwstart:nwflux))
       else
          ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
          fC(ixCmin1:ixCmax1,iwstart:nwflux,idims)=(-f(kxCmin1:kxCmax1,&
             iwstart:nwflux)+7.0d0*(f(jxCmin1:jxCmax1,&
             iwstart:nwflux) + f(ixCmin1:ixCmax1,&
             iwstart:nwflux))-f(hxCmin1:hxCmax1,iwstart:nwflux))/12.0d0

          do iw=iwstart,nwflux
             ! add rempel dissipative flux, only second order version for now
             fC(ixCmin1:ixCmax1,iw,idims)=fC(ixCmin1:ixCmax1,iw,&
                idims)-tvdlfeps*half*vLC(ixCmin1:ixCmax1) &
                *(wRC(ixCmin1:ixCmax1,iw)-wLC(ixCmin1:ixCmax1,iw))
          end do
       end if

    end do       !next idims
    b0i=0

    if(stagger_grid) call phys_update_faces(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,&
       qdt,wprim,fC,fE,sCT,s,vcts)

    if(slab_uniform) then
      dxinv=-qdt/dxs
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
        do iw=iwstart,nwflux
          fC(ixImin1:ixImax1,iw,idims)=dxinv(idims)*fC(ixImin1:ixImax1,iw,&
             idims)
          ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)-f_(i-1))+f_(i-2)]/12
          w(ixOmin1:ixOmax1,iw)=w(ixOmin1:ixOmax1,iw)+(fC(ixOmin1:ixOmax1,iw,&
             idims)-fC(hxOmin1:hxOmax1,iw,idims))
        end do    !next iw
      end do ! Next idims
    else
      inv_volume=1.d0/block%dvolume
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
        do iw=iwstart,nwflux
          fC(ixImin1:ixImax1,iw,idims)=-qdt*block%surfaceC(ixImin1:ixImax1,&
             idims)*fC(ixImin1:ixImax1,iw,idims)
          w(ixOmin1:ixOmax1,iw)=w(ixOmin1:ixOmax1,iw)+ (fC(ixOmin1:ixOmax1,iw,&
             idims)-fC(hxOmin1:hxOmax1,iw,idims))*inv_volume(ixOmin1:ixOmax1)
        end do    !next iw
      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImax1,ixOmin1,ixOmax1,wCT,w,x)

    if(stagger_grid) call phys_face_to_center(ixOmin1,ixOmax1,s)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,w,x,ixImin1,ixImax1,ixOmin1,&
          ixOmax1,'centdiff')
    endif

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), ixImin1,ixImax1,&
       ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,w,x,.false.,active,wprim)

    if(phys_solve_eaux.and.levmin==levmax) then
      ! synchronize internal energy for uniform grid
      call phys_energy_synchro(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    endif
    end associate
  end subroutine centdiff

end module mod_finite_difference
