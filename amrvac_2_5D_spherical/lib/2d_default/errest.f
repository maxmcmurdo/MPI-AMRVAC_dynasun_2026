!> Do all local error estimation which determines (de)refinement
subroutine errest
  use mod_forest, only: refine, buffer
  use mod_global_parameters

  integer :: igrid, iigrid, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2
  double precision :: factor
  logical, dimension(:,:), allocatable :: refine2

  if (igridstail==0) return

  select case (refine_criterion)
  case (0) 
     ! all refinement solely based on user routine usr_refine_grid
  case (1) 
     ! simply compare w_n-1 with w_n and trigger refinement on relative
     ! differences
  !$OMP PARALLEL DO PRIVATE(igrid)
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        call compare1_grid(igrid,pso(igrid)%w,ps(igrid)%w)
     end do
  !$OMP END PARALLEL DO
  case (2)
     ! Error estimation is based on Lohner's original scheme
  !$OMP PARALLEL DO PRIVATE(igrid)
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        call lohner_orig_grid(igrid)
     end do
  !$OMP END PARALLEL DO
  case (3)
     ! Error estimation is based on Lohner's scheme
  !$OMP PARALLEL DO PRIVATE(igrid)
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        call lohner_grid(igrid)
     end do
  !$OMP END PARALLEL DO
  case default
     call mpistop("Unknown error estimator")
  end select

  ! enforce additional refinement on e.g. coordinate and/or time info here
  if (nbufferx1/=0.or.nbufferx2/=0) then
     allocate(refine2(max_blocks,npe))
     call MPI_ALLREDUCE(refine,refine2,max_blocks*npe,MPI_LOGICAL,MPI_LOR,&
         icomm,ierrmpi)
     refine=refine2
  end if
  !$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     block=>ps(igrid)
     call forcedrefine_grid(igrid,ps(igrid)%w)
  end do
  !$OMP END PARALLEL DO

  if (nbufferx1/=0.or.nbufferx2/=0) buffer=.false.

end subroutine errest

subroutine lohner_grid(igrid)
  use mod_usr_methods, only: usr_var_for_errest, usr_refine_threshold
  use mod_forest, only: coarsen, refine
  use mod_physics, only: phys_energy
  use mod_global_parameters

  integer, intent(in) :: igrid

  integer                            :: iflag, idims, idims2, level
  integer                            :: ixmin1,ixmin2,ixmax1,ixmax2, hxmin1,&
     hxmin2,hxmax1,hxmax2, jxmin1,jxmin2,jxmax1,jxmax2, h2xmin1,h2xmin2,&
     h2xmax1,h2xmax2, j2xmin1,j2xmin2,j2xmax1,j2xmax2, ix1,ix2
  double precision                   :: epsilon, threshold, wtol(1:nw),&
      xtol(1:ndim)
  double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2) :: numerator,&
      denominator, error
  double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: tmp, tmp1, tmp2
  double precision                   :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)
  logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)          :: refineflag,&
      coarsenflag

  epsilon = 1.0d-6
  level   = node(plevel_,igrid)
  ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmax1=ixMhi1+1;ixmax2=ixMhi2+1;

  error=zero

  w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)=ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
     1:nw)

  if(B0field) then
    if(phys_energy) w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iw_e)=w(ixGlo1:ixGhi1,&
       ixGlo2:ixGhi2,iw_e)+0.5d0*sum(ps(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       :,0)**2,dim=ndim+1) + sum(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       iw_mag(:))*ps(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,:,0),dim=ndim+1)
    w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iw_mag(:))=w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       iw_mag(:))+ps(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,:,0)
  end if

  do iflag=1,nw+1

     if(w_refine_weight(iflag)==0.d0) cycle
     numerator=zero

     if (iflag > nw) then
        if (.not. associated(usr_var_for_errest)) then
           call mpistop("usr_var_for_errest not defined")
        else
           call usr_var_for_errest(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,ixGlo2,&
              ixGhi1,ixGhi2,iflag,ps(igrid)%w,ps(igrid)%x,tmp1)
        end if
     end if

     do idims=1,ndim
        hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
        hxmax1=ixmax1-kr(1,idims);hxmax2=ixmax2-kr(2,idims);
        jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
        jxmax1=ixmax1+kr(1,idims);jxmax2=ixmax2+kr(2,idims);
        if (iflag<=nw) then
          if (logflag(iflag)) then
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(w(jxmin1:jxmax1,&
               jxmin2:jxmax2,iflag))-dlog10(w(hxmin1:hxmax1,hxmin2:hxmax2,&
               iflag))
          else
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=w(jxmin1:jxmax1,jxmin2:jxmax2,&
               iflag)-w(hxmin1:hxmax1,hxmin2:hxmax2,iflag)
          end if
        else
          if (logflag(iflag)) then
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(tmp1(jxmin1:jxmax1,&
               jxmin2:jxmax2))-dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2))
          else
            tmp(ixmin1:ixmax1,ixmin2:ixmax2)=tmp1(jxmin1:jxmax1,&
               jxmin2:jxmax2)-tmp1(hxmin1:hxmax1,hxmin2:hxmax2)
          end if
        end if
        do idims2=1,ndim
           h2xmin1=ixMlo1-kr(1,idims2);h2xmin2=ixMlo2-kr(2,idims2)
           h2xmax1=ixMhi1-kr(1,idims2);h2xmax2=ixMhi2-kr(2,idims2);
           j2xmin1=ixMlo1+kr(1,idims2);j2xmin2=ixMlo2+kr(2,idims2)
           j2xmax1=ixMhi1+kr(1,idims2);j2xmax2=ixMhi2+kr(2,idims2);
           numerator=numerator+(tmp(j2xmin1:j2xmax1,&
              j2xmin2:j2xmax2)-tmp(h2xmin1:h2xmax1,h2xmin2:h2xmax2))**2.0d0
        end do
     end do
     denominator=zero
     do idims=1,ndim
        if (iflag<=nw) then
           if (logflag(iflag)) then
            tmp=dabs(dlog10(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iflag)))
           else
            tmp=dabs(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iflag))
           end if
        else
           if (logflag(iflag)) then
            tmp=dabs(dlog10(tmp1(ixGlo1:ixGhi1,ixGlo2:ixGhi2)))
           else
            tmp=dabs(tmp1(ixGlo1:ixGhi1,ixGlo2:ixGhi2))
           end if
        end if
        hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
        hxmax1=ixmax1-kr(1,idims);hxmax2=ixmax2-kr(2,idims);
        jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
        jxmax1=ixmax1+kr(1,idims);jxmax2=ixmax2+kr(2,idims);
        tmp2(ixmin1:ixmax1,ixmin2:ixmax2)=tmp(jxmin1:jxmax1,&
           jxmin2:jxmax2)+tmp(hxmin1:hxmax1,hxmin2:hxmax2)
        hxmin1=ixMlo1-2*kr(1,idims);hxmin2=ixMlo2-2*kr(2,idims)
        hxmax1=ixMhi1-2*kr(1,idims);hxmax2=ixMhi2-2*kr(2,idims);
        jxmin1=ixMlo1+2*kr(1,idims);jxmin2=ixMlo2+2*kr(2,idims)
        jxmax1=ixMhi1+2*kr(1,idims);jxmax2=ixMhi2+2*kr(2,idims);
        if (iflag<=nw) then
          if (logflag(iflag)) then
            tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(dlog10(w(jxmin1:jxmax1,&
               jxmin2:jxmax2,iflag))-dlog10(w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
               iflag))) +dabs(dlog10(w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
               iflag))-dlog10(w(hxmin1:hxmax1,hxmin2:hxmax2,iflag)))
          else
             tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(w(jxmin1:jxmax1,&
                jxmin2:jxmax2,iflag)-w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                iflag)) +dabs(w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                iflag)-w(hxmin1:hxmax1,hxmin2:hxmax2,iflag))
          end if
        else
          if (logflag(iflag)) then
            tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(dlog10(tmp1(jxmin1:jxmax1,&
               jxmin2:jxmax2))-dlog10(tmp1(ixMlo1:ixMhi1,&
               ixMlo2:ixMhi2))) +dabs(dlog10(tmp1(ixMlo1:ixMhi1,&
               ixMlo2:ixMhi2))-dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2)))
          else
             tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(tmp1(jxmin1:jxmax1,&
                jxmin2:jxmax2)-tmp1(ixMlo1:ixMhi1,&
                ixMlo2:ixMhi2)) +dabs(tmp1(ixMlo1:ixMhi1,&
                ixMlo2:ixMhi2)-tmp1(hxmin1:hxmax1,hxmin2:hxmax2))
          end if
        end if
        do idims2=1,ndim
           h2xmin1=ixMlo1-kr(1,idims2);h2xmin2=ixMlo2-kr(2,idims2)
           h2xmax1=ixMhi1-kr(1,idims2);h2xmax2=ixMhi2-kr(2,idims2);
           j2xmin1=ixMlo1+kr(1,idims2);j2xmin2=ixMlo2+kr(2,idims2)
           j2xmax1=ixMhi1+kr(1,idims2);j2xmax2=ixMhi2+kr(2,idims2);
           denominator=denominator +(tmp(ixMlo1:ixMhi1,&
              ixMlo2:ixMhi2)+amr_wavefilter(level)*(tmp2(j2xmin1:j2xmax1,&
              j2xmin2:j2xmax2)+tmp2(h2xmin1:h2xmax1,h2xmin2:h2xmax2)))**2
        end do
     end do
     error=error+w_refine_weight(iflag)*dsqrt(numerator/max(denominator,&
        epsilon))
  end do
  
  refineflag=.false.
  coarsenflag=.false.
  threshold=refine_threshold(level)
  do ix2=ixMlo2,ixMhi2
  do ix1=ixMlo1,ixMhi1
  
     if (associated(usr_refine_threshold)) then
        wtol(1:nw)   = w(ix1,ix2,1:nw)
        xtol(1:ndim) = ps(igrid)%x(ix1,ix2,1:ndim)
        call usr_refine_threshold(wtol, xtol, threshold, global_time, level)
     end if
  
     if (error(ix1,ix2) >= threshold) then
        refineflag(ix1,ix2) = .true.
     else if (error(ix1,ix2) <= derefine_ratio(level)*threshold) then
        coarsenflag(ix1,ix2) = .true.
     end if
  end do
  end do
  
  if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level<refine_max_level) &
     refine(igrid,mype)=.true.
  if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level>1) coarsen(igrid,&
     mype)=.true.

end subroutine lohner_grid

subroutine lohner_orig_grid(igrid)
  use mod_usr_methods, only: usr_var_for_errest, usr_refine_threshold
  use mod_forest, only: coarsen, refine
  use mod_global_parameters

  integer, intent(in) :: igrid

  integer :: iflag, idims, level
  integer :: ixmin1,ixmin2,ixmax1,ixmax2, hxmin1,hxmin2,hxmax1,hxmax2, jxmin1,&
     jxmin2,jxmax1,jxmax2, ix1,ix2
  double precision :: epsilon, threshold, wtol(1:nw), xtol(1:ndim)
  double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2) :: numerator,&
      denominator, error
  double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: dp, dm, dref,&
      tmp1
  logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag, coarsenflag

  epsilon=1.0d-6
  level=node(plevel_,igrid)
  ixmin1=ixMlo1;ixmin2=ixMlo2;ixmax1=ixMhi1;ixmax2=ixMhi2;

  error=zero
  do iflag=1,nw+1
     if(w_refine_weight(iflag)==0.d0) cycle
     numerator=zero
     denominator=zero

     if (iflag > nw) then
        if (.not. associated(usr_var_for_errest)) then
           call mpistop("usr_var_for_errest not defined")
        else
           call usr_var_for_errest(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,ixGlo2,&
              ixGhi1,ixGhi2,iflag,ps(igrid)%w,ps(igrid)%x,tmp1)
        end if
     end if

     do idims=1,ndim
        hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
        hxmax1=ixmax1-kr(1,idims);hxmax2=ixmax2-kr(2,idims);
        jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
        jxmax1=ixmax1+kr(1,idims);jxmax2=ixmax2+kr(2,idims);
        if (iflag<=nw) then
          if (logflag(iflag)) then
            dp(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(ps(igrid)%w(jxmin1:jxmax1,&
               jxmin2:jxmax2,iflag))-dlog10(ps(igrid)%w(ixmin1:ixmax1,&
               ixmin2:ixmax2,iflag))
            dm(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(ps(igrid)%w(ixmin1:ixmax1,&
               ixmin2:ixmax2,iflag))-dlog10(ps(igrid)%w(hxmin1:hxmax1,&
               hxmin2:hxmax2,iflag))
            dref(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(dlog10(ps(igrid)%w(&
               jxmin1:jxmax1,jxmin2:jxmax2,&
               iflag)))+ 2.0d0 * dabs(dlog10(ps(igrid)%w(ixMlo1:ixMhi1,&
               ixMlo2:ixMhi2,iflag))) + dabs(dlog10(ps(igrid)%w(hxmin1:hxmax1,&
               hxmin2:hxmax2,iflag)))
          else
            dp(ixmin1:ixmax1,ixmin2:ixmax2)=ps(igrid)%w(jxmin1:jxmax1,&
               jxmin2:jxmax2,iflag)-ps(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
               iflag)
            dm(ixmin1:ixmax1,ixmin2:ixmax2)=ps(igrid)%w(ixmin1:ixmax1,&
               ixmin2:ixmax2,iflag)-ps(igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,&
               iflag)
            dref(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(ps(igrid)%w(jxmin1:jxmax1,&
               jxmin2:jxmax2,iflag))+2.0d0*dabs(ps(igrid)%w(ixMlo1:ixMhi1,&
               ixMlo2:ixMhi2,iflag)) +dabs(ps(igrid)%w(hxmin1:hxmax1,&
               hxmin2:hxmax2,iflag))
          end if
        else
          if (logflag(iflag)) then
            dp(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(tmp1(jxmin1:jxmax1,&
               jxmin2:jxmax2))-dlog10(tmp1(ixmin1:ixmax1,ixmin2:ixmax2))
            dm(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(tmp1(ixmin1:ixmax1,&
               ixmin2:ixmax2))-dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2))
            dref(ixmin1:ixmax1,ixmin2:ixmax2)=dabs(dlog10(tmp1(jxmin1:jxmax1,&
               jxmin2:jxmax2)))+ 2.0d0 * dabs(dlog10(tmp1(ixmin1:ixmax1,&
               ixmin2:ixmax2))) + dabs(dlog10(tmp1(hxmin1:hxmax1,&
               hxmin2:hxmax2)))
          else
            dp(ixmin1:ixmax1,ixmin2:ixmax2)=tmp1(jxmin1:jxmax1,&
               jxmin2:jxmax2)-tmp1(ixmin1:ixmax1,ixmin2:ixmax2)
            dm(ixmin1:ixmax1,ixmin2:ixmax2)=tmp1(ixmin1:ixmax1,&
               ixmin2:ixmax2)-tmp1(hxmin1:hxmax1,hxmin2:hxmax2)
            dref(ixmin1:ixmax1,ixmin2:ixmax2)=dabs(tmp1(jxmin1:jxmax1,&
               jxmin2:jxmax2))+2.0d0*dabs(tmp1(ixmin1:ixmax1,&
               ixmin2:ixmax2)) +dabs(tmp1(hxmin1:hxmax1,hxmin2:hxmax2))
          end if
        end if

        numerator(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=numerator+(dp(ixMlo1:ixMhi1,&
           ixMlo2:ixMhi2)-dm(ixMlo1:ixMhi1,ixMlo2:ixMhi2))**2.0d0

        denominator(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=denominator + &
           (dabs(dp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)) + dabs(dm(ixMlo1:ixMhi1,&
           ixMlo2:ixMhi2)) + amr_wavefilter(level)*dref(ixMlo1:ixMhi1,&
           ixMlo2:ixMhi2))**2.0d0

     end do
     error=error+w_refine_weight(iflag)*dsqrt(numerator/max(denominator,&
        epsilon))
  end do

  refineflag=.false.
  coarsenflag=.false.

  threshold=refine_threshold(level)
  do ix2=ixMlo2,ixMhi2
  do ix1=ixMlo1,ixMhi1

     if (associated(usr_refine_threshold)) then
        wtol(1:nw)   = ps(igrid)%w(ix1,ix2,1:nw)
        xtol(1:ndim) = ps(igrid)%x(ix1,ix2,1:ndim)
        call usr_refine_threshold(wtol, xtol, threshold, global_time, level)
     end if

     if (error(ix1,ix2) >= threshold) then
        refineflag(ix1,ix2) = .true.
     else if (error(ix1,ix2) <= derefine_ratio(level)*threshold) then
        coarsenflag(ix1,ix2) = .true.
     end if
  end do
  end do

  if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level<refine_max_level) &
     refine(igrid,mype)=.true.
  if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level>1) coarsen(igrid,&
     mype)=.true.

end subroutine lohner_orig_grid

subroutine compare1_grid(igrid,wold,w)
  use mod_usr_methods, only: usr_refine_threshold
  use mod_forest, only: coarsen, refine
  use mod_global_parameters

  integer, intent(in) :: igrid
  double precision, intent(in) :: wold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
      w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)

  integer :: ix1,ix2, iflag, level
  double precision :: epsilon, threshold, wtol(1:nw), xtol(1:ndim)
  double precision :: average, error
  double precision :: averages(nw)
  logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag, coarsenflag

  ! identify the points to be flagged in two steps:
  !  step I: compare w_n-1 with w_n solution, store w_for_refine in auxiliary
  !  step II: transfer w_for_refine from auxiliary to refine and coarsen

  epsilon=1.0d-6

  refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = .false.
  coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = .false.
  level=node(plevel_,igrid)
  threshold=refine_threshold(level)
  do ix2=ixMlo2,ixMhi2 
  do ix1=ixMlo1,ixMhi1 
     average=zero
     error=zero
     do iflag=1,nw+1
        if(w_refine_weight(iflag)==0) cycle
        averages(iflag) = w(ix1,ix2,iflag)-wold(ix1,ix2,iflag)
        average=average+w_refine_weight(iflag)*abs(averages(iflag))
        if (abs(wold(ix1,ix2,iflag))<smalldouble)then
           error=error+w_refine_weight(iflag)* &
              abs(averages(iflag))/(abs(wold(ix1,ix2,iflag))+epsilon)
        else
           error=error+w_refine_weight(iflag)* &
              abs(averages(iflag))/(abs(wold(ix1,ix2,iflag)))
        end if
     end do

     if (associated(usr_refine_threshold)) then
        wtol(1:nw)   = ps(igrid)%w(ix1,ix2,1:nw)
        xtol(1:ndim) = ps(igrid)%x(ix1,ix2,1:ndim)
        call usr_refine_threshold(wtol, xtol, threshold, global_time, level)
     end if

     if (error >= threshold) then
        refineflag(ix1,ix2) = .true.
     else if (error <= derefine_ratio(level)*threshold) then
        coarsenflag(ix1,ix2) = .true.
     end if
  end do
  end do

  if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2))) then
     if (level<refine_max_level) refine(igrid,mype)=.true.
  end if
  if (time_advance) then
     if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level>1) &
        coarsen(igrid,mype)=.true.
  end if

end subroutine compare1_grid

subroutine forcedrefine_grid(igrid,w)
  use mod_usr_methods, only: usr_refine_grid
  use mod_forest, only: coarsen, refine, buffer
  use mod_global_parameters

  integer, intent(in) :: igrid
  double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

  integer :: level
  integer :: my_refine, my_coarsen
  double precision :: qt
  logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag

  level=node(plevel_,igrid)

  ! initialize to 0
  my_refine   = 0
  my_coarsen  = 0

  if (time_advance) then
     qt=global_time+dt
  else
     if (refine_criterion==1) then
        qt=global_time+dt
     else
        qt=global_time
     end if
  end if

  if (associated(usr_refine_grid)) then
     call usr_refine_grid(igrid,level,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,&
        ixMlo2,ixMhi1,ixMhi2,qt,w,ps(igrid)%x, my_refine,my_coarsen)
  end if

  if (my_coarsen==1) then
     if (level>1) then
        refine(igrid,mype)=.false.
        coarsen(igrid,mype)=.true.
     else
        refine(igrid,mype)=.false.
        coarsen(igrid,mype)=.false.
     end if
  endif

  if (my_coarsen==-1)then
     coarsen(igrid,mype)=.false.
  end if

  if (my_refine==1) then
     if (level<refine_max_level) then
        refine(igrid,mype)=.true.
        coarsen(igrid,mype)=.false.
     else
        refine(igrid,mype)=.false.
        coarsen(igrid,mype)=.false.
     end if
  end if

  if (my_refine==-1) then
    refine(igrid,mype)=.false.
  end if

  if (nbufferx1/=0.or.nbufferx2/=0) then
     if (refine(igrid,mype) .and. .not.buffer(igrid,mype)) then
        refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=.true.
        call refinebuffer(igrid,refineflag)
     end if
  end if

end subroutine forcedrefine_grid

subroutine forcedrefine_grid_io(igrid,w)
  use mod_forest, only: coarsen, refine
  use mod_global_parameters

  integer, intent(in)          :: igrid
  double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

  integer                   :: level, my_levmin, my_levmax
  logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag

  level=node(plevel_,igrid)

  if (level_io > 0) then
     my_levmin = level_io
     my_levmax = level_io
  else
     my_levmin = max(1,level_io_min)
     my_levmax = min(refine_max_level,level_io_max)
  end if


  if (level>my_levmax) then
        refine(igrid,mype)=.false.
        coarsen(igrid,mype)=.true.
  elseif (level<my_levmin) then
        refine(igrid,mype)=.true.
        coarsen(igrid,mype)=.false.
  end if

  if (level==my_levmin .or. level==my_levmax) then
    refine(igrid,mype)=.false.
    coarsen(igrid,mype)=.false.
  end if


  if(refine(igrid,mype).and.level>=refine_max_level)refine(igrid,mype)=.false.
  if(coarsen(igrid,mype).and.level<=1)coarsen(igrid,mype)=.false.

end subroutine forcedrefine_grid_io

subroutine refinebuffer(igrid,refineflag)
  use mod_forest, only: refine, buffer
  use mod_global_parameters

  integer, intent(in) :: igrid
  logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2), intent(in) :: refineflag

  integer :: ishiftbuf1,ishiftbuf2, i1,i2, ixmin1,ixmin2,ixmax1,ixmax2,&
      ineighbor, ipe_neighbor, level

  ishiftbuf1=ixMhi1-ixMlo1-nbufferx1+1;ishiftbuf2=ixMhi2-ixMlo2-nbufferx2+1;
  do i2=-1,1
  do i1=-1,1
     ixmin1=max(ixMlo1,ixMlo1+i1*ishiftbuf1)
     ixmin2=max(ixMlo2,ixMlo2+i2*ishiftbuf2);
     ixmax1=min(ixMhi1,ixMhi1+i1*ishiftbuf1)
     ixmax2=min(ixMhi2,ixMhi2+i2*ishiftbuf2);
     if (ixmax1<ixmin1.or.ixmax2<ixmin2) cycle
     if (any(refineflag(ixmin1:ixmax1,ixmin2:ixmax2))) then
        select case (neighbor_type(i1,i2,igrid))
        case (neighbor_coarse)
           ineighbor=neighbor(1,i1,i2,igrid)
           ipe_neighbor=neighbor(2,i1,i2,igrid)
           if (.not.refine(ineighbor,ipe_neighbor)) then
              buffer(ineighbor,ipe_neighbor)=.true.
              refine(ineighbor,ipe_neighbor)=.true.
           end if
        case (neighbor_sibling)
           level=node(plevel_,igrid)
           if (level<refine_max_level) then
              ineighbor=neighbor(1,i1,i2,igrid)
              ipe_neighbor=neighbor(2,i1,i2,igrid)
              if (.not.refine(ineighbor,ipe_neighbor)) then
                 buffer(ineighbor,ipe_neighbor)=.true.
                 refine(ineighbor,ipe_neighbor)=.true.
              end if
           end if
        end select
     end if
  end do
  end do

end subroutine refinebuffer
