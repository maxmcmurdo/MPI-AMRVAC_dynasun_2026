!> fill ghost cells at a physical boundary
subroutine bc_phys(iside,idims,time,qdt,s,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
   ixBmin1,ixBmin2,ixBmax1,ixBmax2)
  use mod_usr_methods, only: usr_special_bc
  use mod_bc_data, only: bc_data_set
  use mod_global_parameters
  use mod_physics
  use mod_mhd_phys, only: get_divb

  integer, intent(in) :: iside, idims, ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2
  double precision, intent(in) :: time,qdt
  type(state), intent(inout) :: s
  double precision :: wtmp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nwflux)

  integer :: idir, is
  integer :: ixOsmin1,ixOsmin2,ixOsmax1,ixOsmax2,hxOmin1,hxOmin2,hxOmax1,&
     hxOmax2,jxOmin1,jxOmin2,jxOmax1,jxOmax2
  double precision :: Q(ixGmin1:ixGmax1,ixGmin2:ixGmax2),Qp(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2) 
  integer :: iw, iB, ix1,ix2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixMmin1,ixMmin2,&
     ixMmax1,ixMmax2, nghostcellsi,iib1,iib2
  logical  :: isphysbound

  associate(x=>s%x,w=>s%w,ws=>s%ws)
  select case (idims)
  case (1)
     if (iside==2) then
        ! maximal boundary
        iB=2*1
        ixOmin1=ixBmax1+1-nghostcells;ixOmin2=ixBmin2;
        ixOmax1=ixBmax1;ixOmax2=ixBmax2;
        ! cont/symm/asymm types
        do iw=1,nwflux+nwaux
           select case (typeboundary(iw,iB))
           case (bc_special)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_cont)
              do ix1=ixOmin1,ixOmax1
                 w(ix1,ixOmin2:ixOmax2,iw) = w(ixOmin1-1,ixOmin2:ixOmax2,iw)
              end do
           case (bc_symm)
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 iw) = w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,iw)
           case (bc_asymm)
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 iw) =-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,iw)
           case (bc_periodic)
              ! skip it here, periodic bc info should come from neighbors
           case(bc_noinflow)
              if (iw==1+1)then
                do ix1=ixOmin1,ixOmax1
                    w(ix1,ixOmin2:ixOmax2,iw) = max(w(ixOmin1-1,&
                       ixOmin2:ixOmax2,iw),zero)
                end do
              else
                do ix1=ixOmin1,ixOmax1
                    w(ix1,ixOmin2:ixOmax2,iw) = w(ixOmin1-1,ixOmin2:ixOmax2,&
                       iw)
                end do
              end if
           case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = - w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw)
           case (bc_data)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_character)
              ! skip it here, do AFTER all normal type boundaries are set
           case default
              write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                  "for variable iw=",iw," and side iB=",iB
           end select
        end do
        if(stagger_grid) then
          do idir=1,nws
 !At this stage, extrapolation is applied only to the tangential components
            if(idir==1) cycle 
            ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
            ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir);
            select case(typeboundary(iw_mag(idir),iB))
            case (bc_special)
               ! skip it here, do AFTER all normal type boundaries are set
            case (bc_cont)
              do ix1=ixOsmin1,ixOsmax1
                 ws(ix1,ixOsmin2:ixOsmax2,idir) = ws(ixOsmin1-1,&
                    ixOsmin2:ixOsmax2,idir)
              end do
            case (bc_symm)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir) = ws(ixOsmin1-1:ixOsmin1-nghostcells:-1,&
                 ixOsmin2:ixOsmax2,idir)
            case (bc_asymm)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir) =-ws(ixOsmin1-1:ixOsmin1-nghostcells:-1,&
                 ixOsmin2:ixOsmax2,idir)
            case (bc_periodic)
            case default
              write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                  "for variable iw=",iw," and side iB=",iB
            end select
          end do
          ! Now that the tangential components are set,
 !fill the normal components using a prescription for the divergence.
          ! This prescription is given by the typeB for the normal component.
          do idir=1,nws
            ! Consider only normal direction
            if (idir/=1) cycle
            ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmax1=ixOmax1
            ixOsmax2=ixOmax2;
            select case(typeboundary(iw_mag(idir),iB))
            case(bc_cont)
              hxOmin1=ixOmin1-nghostcells*kr(1,1)
              hxOmin2=ixOmin2-nghostcells*kr(2,1)
              hxOmax1=ixOmax1-nghostcells*kr(1,1)
              hxOmax2=ixOmax2-nghostcells*kr(2,1);
              ! Calculate divergence and partial divergence
              call get_divb(s%w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,hxOmin1,&
                 hxOmin2,hxOmax1,hxOmax2,Q)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,idir)=zero
              do ix1=0,nghostcells-1
                call get_divb(s%w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
                   ixOmin2,ixOmax1,ixOmax2,Qp)
                ws(ixOsmin1+ix1,ixOsmin2:ixOsmax2,idir)=(Q(hxOmax1,&
                   hxOmin2:hxOmax2)*s%dvolume(hxOmax1,&
                   hxOmin2:hxOmax2)-Qp(ixOmin1+ix1,&
                   ixOmin2:ixOmax2)*s%dvolume(ixOmin1+ix1,&
                   ixOmin2:ixOmax2))/s%surfaceC(ixOsmin1+ix1,ixOsmin2:ixOsmax2,&
                   1)
              end do
            case(bc_symm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir)= ws(ixOsmin1-2:ixOsmin1-nghostcells-1:-1,&
                 ixOsmin2:ixOsmax2,idir)
            case(bc_asymm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir)=-ws(ixOsmin1-2:ixOsmin1-nghostcells-1:-1,&
                 ixOsmin2:ixOsmax2,idir)
            case(bc_periodic)
            end select
          end do
          ! Fill cell averages
          call phys_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,s)
        end if
     else
        ! minimal boundary
        iB=2*1-1
        ixOmin1=ixBmin1;ixOmin2=ixBmin2;
        ixOmax1=ixBmin1-1+nghostcells;ixOmax2=ixBmax2;
        ! cont/symm/asymm types
        do iw=1,nwflux+nwaux
           select case (typeboundary(iw,iB))
           case (bc_special)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_cont)
              do ix1=ixOmin1,ixOmax1
                 w(ix1,ixOmin2:ixOmax2,iw) = w(ixOmax1+1,ixOmin2:ixOmax2,iw)
              end do
           case (bc_symm)
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 iw) = w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,iw)
           case (bc_asymm)
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 iw) =-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,iw)
           case (bc_periodic)
              ! skip it here, periodic bc info should come from neighbors
           case(bc_noinflow)
              if (iw==1+1)then
                 do ix1=ixOmin1,ixOmax1
                   w(ix1,ixOmin2:ixOmax2,iw) = min(w(ixOmax1+1,ixOmin2:ixOmax2,&
                      iw),zero)
                 end do
              else
                 do ix1=ixOmin1,ixOmax1
                   w(ix1,ixOmin2:ixOmax2,iw) = w(ixOmax1+1,ixOmin2:ixOmax2,iw)
                 end do
              end if
           case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = - w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw)
           case (bc_data)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_character)
              ! skip it here, do AFTER all normal type boundaries are set
           case default
              write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                  "for variable iw=",iw," and side iB=",iB
           end select
        end do
        if(stagger_grid) then
          do idir=1,nws
 !At this stage, extrapolation is applied only to the tangential components
            if(idir==1) cycle 
            ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
            ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir);
            select case(typeboundary(iw_mag(idir),iB))
            case (bc_special)
              ! skip it here, periodic bc info should come from neighbors
            case (bc_cont)
              do ix1=ixOsmin1,ixOsmax1
                 ws(ix1,ixOsmin2:ixOsmax2,idir) = ws(ixOsmax1+1,&
                    ixOsmin2:ixOsmax2,idir)
              end do
            case (bc_symm)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir) = ws(ixOsmax1+nghostcells:ixOsmax1+1:-1,&
                 ixOsmin2:ixOsmax2,idir)
            case (bc_asymm)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir) =-ws(ixOsmax1+nghostcells:ixOsmax1+1:-1,&
                 ixOsmin2:ixOsmax2,idir)
            case (bc_periodic)
            case default
              write (unitterm,*) "Undefined boundarytype in bc_phys",&
                  "for variable iw=",iw," and side iB=",iB
            end select
          end do
          ! Now that the tangential components are set,
 !fill the normal components using a prescription for the divergence.
          ! This prescription is given by the typeB for the normal component.
          do idir=1,nws
            ! Consider only normal direction
            if (idir/=1) cycle
            ixOsmin1=ixOmin1-kr(1,1);ixOsmin2=ixOmin2-kr(2,1)
            ixOsmax1=ixOmax1-kr(1,1);ixOsmax2=ixOmax2-kr(2,1);
            select case(typeboundary(iw_mag(idir),iB))
            case(bc_cont)
              jxOmin1=ixOmin1+nghostcells*kr(1,1)
              jxOmin2=ixOmin2+nghostcells*kr(2,1)
              jxOmax1=ixOmax1+nghostcells*kr(1,1)
              jxOmax2=ixOmax2+nghostcells*kr(2,1);
              ! Calculate divergence and partial divergence
              call get_divb(s%w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,jxOmin1,&
                 jxOmin2,jxOmax1,jxOmax2,Q)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,idir)=zero
              do ix1=0,nghostcells-1
                call get_divb(s%w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
                   ixOmin2,ixOmax1,ixOmax2,Qp)
                ws(ixOsmax1-ix1,ixOsmin2:ixOsmax2,idir)=-(Q(jxOmin1,&
                   jxOmin2:jxOmax2)*s%dvolume(jxOmin1,&
                   jxOmin2:jxOmax2)-Qp(ixOmax1-ix1,&
                   ixOmin2:ixOmax2)*s%dvolume(ixOmax1-ix1,&
                   ixOmin2:ixOmax2))/s%surfaceC(ixOsmax1-ix1,ixOsmin2:ixOsmax2,&
                   1)
              end do
            case(bc_symm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir)= ws(ixOsmax1+nghostcells+1:ixOsmax1+2:-1,&
                 ixOsmin2:ixOsmax2,idir)
            case(bc_asymm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir)=-ws(ixOsmax1+nghostcells+1:ixOsmax1+2:-1,&
                 ixOsmin2:ixOsmax2,idir)
            case(bc_periodic)
            end select
          end do
          ! Fill cell averages
          call phys_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,s)
        end if
     end if 
  case (2)
     if (iside==2) then
        ! maximal boundary
        iB=2*2
        ixOmin1=ixBmin1;ixOmin2=ixBmax2+1-nghostcells;
        ixOmax1=ixBmax1;ixOmax2=ixBmax2;
        ! cont/symm/asymm types
        do iw=1,nwflux+nwaux
           select case (typeboundary(iw,iB))
           case (bc_special)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_cont)
              do ix2=ixOmin2,ixOmax2
                 w(ixOmin1:ixOmax1,ix2,iw) = w(ixOmin1:ixOmax1,ixOmin2-1,iw)
              end do
           case (bc_symm)
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = w(ixOmin1:ixOmax1,&
                 ixOmin2-1:ixOmin2-nghostcells:-1,iw)
           case (bc_asymm)
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =-w(ixOmin1:ixOmax1,&
                 ixOmin2-1:ixOmin2-nghostcells:-1,iw)
           case (bc_periodic)
              ! skip it here, periodic bc info should come from neighbors
           case(bc_noinflow)
              if (iw==1+2)then
                do ix2=ixOmin2,ixOmax2
                    w(ixOmin1:ixOmax1,ix2,iw) = max(w(ixOmin1:ixOmax1,&
                       ixOmin2-1,iw),zero)
                end do
              else
                do ix2=ixOmin2,ixOmax2
                    w(ixOmin1:ixOmax1,ix2,iw) = w(ixOmin1:ixOmax1,ixOmin2-1,&
                       iw)
                end do
              end if
           case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = - w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw)
           case (bc_data)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_character)
              ! skip it here, do AFTER all normal type boundaries are set
           case default
              write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                  "for variable iw=",iw," and side iB=",iB
           end select
        end do
        if(stagger_grid) then
          do idir=1,nws
 !At this stage, extrapolation is applied only to the tangential components
            if(idir==2) cycle 
            ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
            ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir);
            select case(typeboundary(iw_mag(idir),iB))
            case (bc_special)
               ! skip it here, do AFTER all normal type boundaries are set
            case (bc_cont)
              do ix2=ixOsmin2,ixOsmax2
                 ws(ixOsmin1:ixOsmax1,ix2,idir) = ws(ixOsmin1:ixOsmax1,&
                    ixOsmin2-1,idir)
              end do
            case (bc_symm)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir) = ws(ixOsmin1:ixOsmax1,&
                 ixOsmin2-1:ixOsmin2-nghostcells:-1,idir)
            case (bc_asymm)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir) =-ws(ixOsmin1:ixOsmax1,&
                 ixOsmin2-1:ixOsmin2-nghostcells:-1,idir)
            case (bc_periodic)
            case default
              write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                  "for variable iw=",iw," and side iB=",iB
            end select
          end do
          ! Now that the tangential components are set,
 !fill the normal components using a prescription for the divergence.
          ! This prescription is given by the typeB for the normal component.
          do idir=1,nws
            ! Consider only normal direction
            if (idir/=2) cycle
            ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmax1=ixOmax1
            ixOsmax2=ixOmax2;
            select case(typeboundary(iw_mag(idir),iB))
            case(bc_cont)
              hxOmin1=ixOmin1-nghostcells*kr(1,2)
              hxOmin2=ixOmin2-nghostcells*kr(2,2)
              hxOmax1=ixOmax1-nghostcells*kr(1,2)
              hxOmax2=ixOmax2-nghostcells*kr(2,2);
              ! Calculate divergence and partial divergence
              call get_divb(s%w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,hxOmin1,&
                 hxOmin2,hxOmax1,hxOmax2,Q)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,idir)=zero
              do ix2=0,nghostcells-1
                call get_divb(s%w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
                   ixOmin2,ixOmax1,ixOmax2,Qp)
                ws(ixOsmin1:ixOsmax1,ixOsmin2+ix2,idir)=(Q(hxOmin1:hxOmax1,&
                   hxOmax2)*s%dvolume(hxOmin1:hxOmax1,&
                   hxOmax2)-Qp(ixOmin1:ixOmax1,&
                   ixOmin2+ix2)*s%dvolume(ixOmin1:ixOmax1,&
                   ixOmin2+ix2))/s%surfaceC(ixOsmin1:ixOsmax1,ixOsmin2+ix2,2)
              end do
            case(bc_symm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir)= ws(ixOsmin1:ixOsmax1,&
                 ixOsmin2-2:ixOsmin2-nghostcells-1:-1,idir)
            case(bc_asymm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir)=-ws(ixOsmin1:ixOsmax1,&
                 ixOsmin2-2:ixOsmin2-nghostcells-1:-1,idir)
            case(bc_periodic)
            end select
          end do
          ! Fill cell averages
          call phys_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,s)
        end if
     else
        ! minimal boundary
        iB=2*2-1
        ixOmin1=ixBmin1;ixOmin2=ixBmin2;
        ixOmax1=ixBmax1;ixOmax2=ixBmin2-1+nghostcells;
        ! cont/symm/asymm types
        do iw=1,nwflux+nwaux
           select case (typeboundary(iw,iB))
           case (bc_special)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_cont)
              do ix2=ixOmin2,ixOmax2
                 w(ixOmin1:ixOmax1,ix2,iw) = w(ixOmin1:ixOmax1,ixOmax2+1,iw)
              end do
           case (bc_symm)
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = w(ixOmin1:ixOmax1,&
                 ixOmax2+nghostcells:ixOmax2+1:-1,iw)
           case (bc_asymm)
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =-w(ixOmin1:ixOmax1,&
                 ixOmax2+nghostcells:ixOmax2+1:-1,iw)
           case (bc_periodic)
              ! skip it here, periodic bc info should come from neighbors
           case(bc_noinflow)
              if (iw==1+2)then
                 do ix2=ixOmin2,ixOmax2
                   w(ixOmin1:ixOmax1,ix2,iw) = min(w(ixOmin1:ixOmax1,ixOmax2+1,&
                      iw),zero)
                 end do
              else
                 do ix2=ixOmin2,ixOmax2
                   w(ixOmin1:ixOmax1,ix2,iw) = w(ixOmin1:ixOmax1,ixOmax2+1,iw)
                 end do
              end if
           case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = - w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw)
           case (bc_data)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_character)
              ! skip it here, do AFTER all normal type boundaries are set
           case default
              write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                  "for variable iw=",iw," and side iB=",iB
           end select
        end do
        if(stagger_grid) then
          do idir=1,nws
 !At this stage, extrapolation is applied only to the tangential components
            if(idir==2) cycle 
            ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
            ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir);
            select case(typeboundary(iw_mag(idir),iB))
            case (bc_special)
              ! skip it here, periodic bc info should come from neighbors
            case (bc_cont)
              do ix2=ixOsmin2,ixOsmax2
                 ws(ixOsmin1:ixOsmax1,ix2,idir) = ws(ixOsmin1:ixOsmax1,&
                    ixOsmax2+1,idir)
              end do
            case (bc_symm)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir) = ws(ixOsmin1:ixOsmax1,&
                 ixOsmax2+nghostcells:ixOsmax2+1:-1,idir)
            case (bc_asymm)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir) =-ws(ixOsmin1:ixOsmax1,&
                 ixOsmax2+nghostcells:ixOsmax2+1:-1,idir)
            case (bc_periodic)
            case default
              write (unitterm,*) "Undefined boundarytype in bc_phys",&
                  "for variable iw=",iw," and side iB=",iB
            end select
          end do
          ! Now that the tangential components are set,
 !fill the normal components using a prescription for the divergence.
          ! This prescription is given by the typeB for the normal component.
          do idir=1,nws
            ! Consider only normal direction
            if (idir/=2) cycle
            ixOsmin1=ixOmin1-kr(1,2);ixOsmin2=ixOmin2-kr(2,2)
            ixOsmax1=ixOmax1-kr(1,2);ixOsmax2=ixOmax2-kr(2,2);
            select case(typeboundary(iw_mag(idir),iB))
            case(bc_cont)
              jxOmin1=ixOmin1+nghostcells*kr(1,2)
              jxOmin2=ixOmin2+nghostcells*kr(2,2)
              jxOmax1=ixOmax1+nghostcells*kr(1,2)
              jxOmax2=ixOmax2+nghostcells*kr(2,2);
              ! Calculate divergence and partial divergence
              call get_divb(s%w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,jxOmin1,&
                 jxOmin2,jxOmax1,jxOmax2,Q)
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,idir)=zero
              do ix2=0,nghostcells-1
                call get_divb(s%w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
                   ixOmin2,ixOmax1,ixOmax2,Qp)
                ws(ixOsmin1:ixOsmax1,ixOsmax2-ix2,idir)=-(Q(jxOmin1:jxOmax1,&
                   jxOmin2)*s%dvolume(jxOmin1:jxOmax1,&
                   jxOmin2)-Qp(ixOmin1:ixOmax1,&
                   ixOmax2-ix2)*s%dvolume(ixOmin1:ixOmax1,&
                   ixOmax2-ix2))/s%surfaceC(ixOsmin1:ixOsmax1,ixOsmax2-ix2,2)
              end do
            case(bc_symm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir)= ws(ixOsmin1:ixOsmax1,&
                 ixOsmax2+nghostcells+1:ixOsmax2+2:-1,idir)
            case(bc_asymm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,&
                 idir)=-ws(ixOsmin1:ixOsmax1,&
                 ixOsmax2+nghostcells+1:ixOsmax2+2:-1,idir)
            case(bc_periodic)
            end select
          end do
          ! Fill cell averages
          call phys_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,s)
        end if
     end if 
  end select

  ! do user defined special boundary conditions
  if (any(typeboundary(1:nwflux+nwaux,iB)==bc_special)) then
     if (.not. associated(usr_special_bc)) call &
        mpistop("usr_special_bc not defined")
     call usr_special_bc(time,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,iB,w,x)
  end if

  ! fill boundary conditions from external data vtk files
  if (any(typeboundary(1:nwflux+nwaux,iB)==bc_data)) then
     call bc_data_set(time,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
        ixOmax1,ixOmax2,iB,w,x)
  end if

  
  !end do

  end associate
end subroutine bc_phys

!> fill inner boundary values
subroutine getintbc(time,ixGmin1,ixGmin2,ixGmax1,ixGmax2)
  use mod_usr_methods, only: usr_internal_bc
  use mod_global_parameters

  double precision, intent(in)   :: time
  integer, intent(in)            :: ixGmin1,ixGmin2,ixGmax1,ixGmax2

  integer :: iigrid, igrid, ixOmin1,ixOmin2,ixOmax1,ixOmax2

  ixOmin1=ixGmin1+nghostcells;ixOmin2=ixGmin2+nghostcells
  ixOmax1=ixGmax1-nghostcells;ixOmax2=ixGmax2-nghostcells;

  !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
  do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
     block=>ps(igrid)
     dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);

     if (associated(usr_internal_bc)) then
        call usr_internal_bc(node(plevel_,igrid),time,ixGmin1,ixGmin2,ixGmax1,&
           ixGmax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,ps(igrid)%w,ps(igrid)%x)
     end if
  end do
  !$OMP END PARALLEL DO

end subroutine getintbc
