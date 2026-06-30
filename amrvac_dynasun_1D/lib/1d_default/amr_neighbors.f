!> find neighors of level-one root blocks
subroutine find_root_neighbor(tree_neighbor,tree,i1)
  use mod_forest
  use mod_global_parameters
  use mod_geometry

  type(tree_node_ptr) :: tree_neighbor, tree
  integer, intent(in) :: i1

  integer :: jg1

  jg1=tree%node%ig1+i1;

  ! find the periodic grid indices, modulo(-1,10)=9
  if (periodB(1)) jg1=1+modulo(jg1-1,ng1(1))

  ! pi-periodicity at pole
  select case (coordinate)
  case (spherical) 
  case (cylindrical)
     if (poleB(1,1).and.jg1==0) then ! cylindrical axis
        jg1=1
        if (1==phi_) jg1=1+modulo(jg1+ng1(1)/2-1,ng1(1))
     end if
  end select

  if (jg1>=1.and.jg1<=ng1(1)) then
     tree_neighbor%node => tree_root(jg1)%node
  else
     nullify(tree_neighbor%node)
  end if

end subroutine find_root_neighbor

!> find neighors of all blocks
subroutine find_neighbor(my_neighbor,my_neighbor_type,tree,i1,pole)
  use mod_forest
  use mod_global_parameters

  type(tree_node_ptr) :: tree, my_neighbor
  integer, intent(in) :: i1
  integer, intent(out) :: my_neighbor_type
  logical, dimension(ndim), intent(out) :: pole

  integer :: level, ig1, ic1, n_ic1, inp1

  pole=.false.
  level=tree%node%level
  if (level==1) then
     call find_root_neighbor(my_neighbor,tree,i1)
     if (associated(my_neighbor%node)) then
       if (phi_ > 0) then
         ig1=tree%node%ig1;
         if ((poleB(2,1).and.ig1==ng1(1).and.i1==1) .or. (poleB(1,&
            1).and.ig1==1.and.i1==-1)) pole(1)=.true.
       end if
        if (my_neighbor%node%leaf) then
           my_neighbor_type=3
        else
           my_neighbor_type=4
        end if
     else
        my_neighbor_type=1
        return
     end if
  else
     ig1=tree%node%ig1;

     if (phi_ > 0) then
       if ((poleB(2,1).and.ig1==ng1(level).and.i1==1) .or. (poleB(1,&
          1).and.ig1==1.and.i1==-1)) pole(1)=.true.
     end if

     ! ic^D is 1 when ig^D is odd, is 2 when ig^D is even
     ic1=1+modulo(ig1-1,2);
     inp1=int((ic1+i1+1)/2)-1;
     my_neighbor%node => tree%node%parent%node
     if (inp1/=0) then
        my_neighbor%node => my_neighbor%node%neighbor(ic1,1)%node
        if (.not.associated(my_neighbor%node)) then
           my_neighbor_type=1
           return
        end if
     end if
     if (my_neighbor%node%leaf) then
        my_neighbor_type=2
     else
        if (i1==0 .or. pole(1)) then
           n_ic1=ic1
        else
           n_ic1=3-ic1 ! switch 1 <--> 2
        end if
        my_neighbor%node => my_neighbor%node%child(n_ic1)%node
        if (associated(my_neighbor%node)) then
           if (my_neighbor%node%leaf) then
              my_neighbor_type=3
           else
              my_neighbor_type=4
           end if
        else
           my_neighbor_type=0
        end if
     end if
  end if

end subroutine find_neighbor

!> asign tree node neighor
subroutine asign_tree_neighbor(tree)
  use mod_forest
  use mod_global_parameters

  type(tree_node_ptr) :: tree

  logical, dimension(ndim) :: pole
  integer :: my_neighbor_type, i1, iside
  type(tree_node_ptr) :: my_neighbor

  do iside=1,2
     i1=kr(1,1)*(2*iside-3);
     call find_neighbor(my_neighbor,my_neighbor_type,tree,i1,pole)
     select case (my_neighbor_type)
     case (neighbor_sibling, neighbor_fine)
        tree%node%neighbor(iside,1)%node => my_neighbor%node
        if (associated(my_neighbor%node)) then
           if (pole(1)) then
              my_neighbor%node%neighbor(iside,1)%node => tree%node
           else
              my_neighbor%node%neighbor(3-iside,1)%node => tree%node
           end if
        end if
     case default
        nullify(tree%node%neighbor(iside,1)%node)
     end select
  end do

end subroutine asign_tree_neighbor
