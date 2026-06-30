!> coarsen one grid to its coarser representative
subroutine coarsen_grid(sFi,ixFiGmin1,ixFiGmin2,ixFiGmax1,ixFiGmax2,ixFimin1,&
   ixFimin2,ixFimax1,ixFimax2,sCo,ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
   ixComin1,ixComin2,ixComax1,ixComax2)
  use mod_global_parameters
  use mod_physics

  type(state), intent(inout)      :: sFi, sCo
  integer, intent(in) :: ixFiGmin1,ixFiGmin2,ixFiGmax1,ixFiGmax2, ixFimin1,&
     ixFimin2,ixFimax1,ixFimax2, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
      ixComin1,ixComin2,ixComax1,ixComax2

  integer :: ixCo1,ixCo2, ixFi1,ixFi2, iw
  double precision :: CoFiratio
  double precision :: B_energy_change(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2)

  associate(wFi=>sFi%w(ixFiGmin1:ixFiGmax1,ixFiGmin2:ixFiGmax2,1:nw),&
      wCo=>sCo%w(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nw))
  staggered: associate(wFis=>sFi%ws,wCos=>sCo%ws)
  ! coarsen by 2 in every direction - conservatively

  if(coarsenprimitive) call phys_to_primitive(ixFiGmin1,ixFiGmin2,ixFiGmax1,&
     ixFiGmax2,ixFimin1,ixFimin2,ixFimax1,ixFimax2,wFi,sFi%x)

  if(slab_uniform) then
    CoFiratio=one/dble(2**ndim)
    do iw=1,nw
       do ixCo2 = ixComin2,ixComax2
          ixFi2=2*(ixCo2-ixComin2)+ixFimin2
       do ixCo1 = ixComin1,ixComax1
          ixFi1=2*(ixCo1-ixComin1)+ixFimin1
          wCo(ixCo1,ixCo2,iw)=sum(wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
             iw))*CoFiratio
       end do
       end do
    end do
  else
    do iw=1,nw
      do ixCo2 = ixComin2,ixComax2
         ixFi2=2*(ixCo2-ixComin2)+ixFimin2
      do ixCo1 = ixComin1,ixComax1
         ixFi1=2*(ixCo1-ixComin1)+ixFimin1
         wCo(ixCo1,ixCo2,iw)= sum(sFi%dvolume(ixFi1:ixFi1+1,&
            ixFi2:ixFi2+1)*wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,&
            iw)) /sCo%dvolume(ixCo1,ixCo2)
      end do
      end do
    end do
  end if

  if(stagger_grid) then
    do iw=1,nws
      ! Start one layer before
      do ixCo2 = ixComin2-kr(2,iw),ixComax2
         ixFi2=2*(ixCo2-ixComin2+kr(2,iw))+ixFimin2-kr(2,iw)
      do ixCo1 = ixComin1-kr(1,iw),ixComax1
         ixFi1=2*(ixCo1-ixComin1+kr(1,iw))+ixFimin1-kr(1,iw)
         ! This if statement catches the axis where surface is zero:
         if (sCo%surfaceC(ixCo1,ixCo2,iw)>1.0d-9*sCo%dvolume(ixCo1,&
            ixCo2)) then !Normal case
           wCos(ixCo1,ixCo2,iw)=sum(sFi%surfaceC(ixFi1:ixFi1+1-kr(iw,1),&
              ixFi2:ixFi2+1-kr(iw,2),iw)*wFis(ixFi1:ixFi1+1-kr(iw,1),&
              ixFi2:ixFi2+1-kr(iw,2),iw)) /sCo%surfaceC(ixCo1,ixCo2,iw)
         else ! On axis
           wCos(ixCo1,ixCo2,iw)=zero
         end if
      end do
      end do
    end do
    if(phys_total_energy.and. .not.coarsenprimitive) then
      B_energy_change(ixComin1:ixComax1,ixComin2:ixComax2)=0.5d0*sum(wCo(&
         ixComin1:ixComax1,ixComin2:ixComax2,iw_mag(:))**2,dim=ndim+1)
    end if
    ! average to fill cell-centred values
    call phys_face_to_center(ixComin1,ixComin2,ixComax1,ixComax2,sCo)
    if(phys_total_energy.and. .not.coarsenprimitive) then
      wCo(ixComin1:ixComax1,ixComin2:ixComax2,iw_e)=wCo(ixComin1:ixComax1,&
         ixComin2:ixComax2,iw_e)-B_energy_change(ixComin1:ixComax1,&
         ixComin2:ixComax2)+0.5d0*sum(wCo(ixComin1:ixComax1,ixComin2:ixComax2,&
         iw_mag(:))**2,dim=ndim+1)
    end if
  end if

  if(coarsenprimitive) then
    call phys_to_conserved(ixFiGmin1,ixFiGmin2,ixFiGmax1,ixFiGmax2,ixFimin1,&
       ixFimin2,ixFimax1,ixFimax2,wFi,sFi%x)
    call phys_to_conserved(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,ixComin1,&
       ixComin2,ixComax1,ixComax2,wCo,sCo%x)
  end if
  end associate staggered
  end associate
end subroutine coarsen_grid
