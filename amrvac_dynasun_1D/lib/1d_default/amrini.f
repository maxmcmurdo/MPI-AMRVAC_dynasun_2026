!> Generate and initialize all grids at the coarsest level (level one)
subroutine initlevelone
  use mod_global_parameters
  use mod_ghostcells_update

  integer :: iigrid, igrid

  levmin=1
  levmax=1

  call init_forest_root

  call getigrids
  call build_connectivity

  ! fill solution space of all root grids
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     call alloc_node(igrid)
     ! in case gradient routine used in initial condition, ensure geometry known
     call initial_condition(igrid)
  end do
  

  ! update ghost cells
  call getbc(global_time,0.d0,ps,iwstart,nwgc)

end subroutine initlevelone

!> fill in initial condition
subroutine initial_condition(igrid)
  ! Need only to set the mesh values (can leave ghost cells untouched)
  use mod_usr_methods, only: usr_init_one_grid
  use mod_global_parameters

  integer, intent(in) :: igrid

  ! in case gradient routine used in initial condition, ensure geometry known
  block=>ps(igrid)
  dxlevel(1)=rnode(rpdx1_,igrid);

  if (.not. associated(usr_init_one_grid)) then
     call mpistop("usr_init_one_grid not defined")
  else
     call usr_init_one_grid(ixGlo1,ixGhi1,ixMlo1,ixMhi1,ps(igrid)%w,&
        ps(igrid)%x)
  end if

end subroutine initial_condition

!> modify initial condition
subroutine modify_IC
  use mod_usr_methods, only: usr_init_one_grid
  use mod_global_parameters

  integer :: iigrid, igrid

  do iigrid=1,igridstail; igrid=igrids(iigrid);
     block=>ps(igrid)
     dxlevel(1)=rnode(rpdx1_,igrid);

     if (.not. associated(usr_init_one_grid)) then
        call mpistop("usr_init_one_grid not defined")
     else
        call usr_init_one_grid(ixGlo1,ixGhi1,ixMlo1,ixMhi1,ps(igrid)%w,&
           ps(igrid)%x)
     end if
  end do

end subroutine modify_IC


