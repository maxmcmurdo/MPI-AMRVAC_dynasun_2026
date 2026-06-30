module mod_physics_roe

  implicit none
  public

  procedure(sub_average), pointer         :: phys_average => null()
  procedure(sub_get_eigenjump), pointer   :: phys_get_eigenjump => null()
  procedure(sub_rtimes), pointer          :: phys_rtimes => null()

  integer :: nworkroe = -1

  abstract interface
     subroutine sub_average(wL, wR, x, ixmin1,ixmax1, idim, wroe, workroe)
       use mod_global_parameters
       import
       integer, intent(in)             :: ixmin1,ixmax1, idim
       double precision, intent(in)    :: wL(ixGlo1:ixGhi1, nw),&
           wR(ixGlo1:ixGhi1, nw)
       double precision, intent(inout) :: wroe(ixGlo1:ixGhi1, nw)
       double precision, intent(inout) :: workroe(ixGlo1:ixGhi1, nworkroe)
       double precision, intent(in)    :: x(ixGlo1:ixGhi1, 1:1)
     end subroutine sub_average

     subroutine sub_get_eigenjump(wL, wR, wC, x, ixmin1,ixmax1, il, idim,&
         smalla, a, jump, workroe)
       use mod_global_parameters
       import
       integer, intent(in)                          :: ixmin1,ixmax1, il, idim
       double precision, dimension(ixGlo1:ixGhi1, nw)       :: wL, wR, wC
       double precision, dimension(ixGlo1:ixGhi1)           :: smalla, a, jump
       double precision, dimension(ixGlo1:ixGhi1, nworkroe) :: workroe
       double precision, intent(in)                 :: x(ixGlo1:ixGhi1, 1:1)
     end subroutine sub_get_eigenjump

     subroutine sub_rtimes(q, w, ixmin1,ixmax1, iw, il, idim, rq, workroe)
       use mod_global_parameters
       import
       integer, intent(in)             :: ixmin1,ixmax1, iw, il, idim
       double precision, intent(in)    :: w(ixGlo1:ixGhi1, nw),&
           q(ixGlo1:ixGhi1)
       double precision, intent(inout) :: rq(ixGlo1:ixGhi1)
       double precision, intent(inout) :: workroe(ixGlo1:ixGhi1, nworkroe)
     end subroutine sub_rtimes
  end interface

contains

  subroutine phys_roe_check()
    if (.not. associated(phys_average)) phys_average => dummy_roe_average

    if (.not. associated(phys_get_eigenjump)) phys_get_eigenjump => &
       dummy_roe_get_eigenjump

    if (.not. associated(phys_rtimes)) phys_rtimes => dummy_roe_rtimes

  end subroutine phys_roe_check

  subroutine dummy_roe_average(wL, wR, x, ixmin1,ixmax1, idim, wroe, workroe)
    use mod_global_parameters
    integer, intent(in)             :: ixmin1,ixmax1, idim
    double precision, intent(in)    :: wL(ixGlo1:ixGhi1, nw), wR(ixGlo1:ixGhi1,&
        nw)
    double precision, intent(inout) :: wroe(ixGlo1:ixGhi1, nw)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1, nworkroe)
    double precision, intent(in)    :: x(ixGlo1:ixGhi1, 1:1)

    call mpistop("dummy_roe_average should not be called")
  end subroutine dummy_roe_average

  subroutine dummy_roe_get_eigenjump(wL, wR, wC, x, ixmin1,ixmax1, il, idim,&
      smalla, a, jump, workroe)
    use mod_global_parameters
    integer, intent(in)                          :: ixmin1,ixmax1, il, idim
    double precision, dimension(ixGlo1:ixGhi1, nw)       :: wL, wR, wC
    double precision, dimension(ixGlo1:ixGhi1)           :: smalla, a, jump
    double precision, dimension(ixGlo1:ixGhi1, nworkroe) :: workroe
    double precision, intent(in)                 :: x(ixGlo1:ixGhi1, 1:1)

    call mpistop("dummy_roe_get_eigenjump should not be called")
  end subroutine dummy_roe_get_eigenjump

  subroutine dummy_roe_rtimes(q, w, ixmin1,ixmax1, iw, il, idim, rq, workroe)
    use mod_global_parameters
    integer, intent(in)             :: ixmin1,ixmax1, iw, il, idim
    double precision, intent(in)    :: w(ixGlo1:ixGhi1, nw), q(ixGlo1:ixGhi1)
    double precision, intent(inout) :: rq(ixGlo1:ixGhi1)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1, nworkroe)

    call mpistop("dummy_roe_rtimes should not be called")
  end subroutine dummy_roe_rtimes

end module mod_physics_roe
