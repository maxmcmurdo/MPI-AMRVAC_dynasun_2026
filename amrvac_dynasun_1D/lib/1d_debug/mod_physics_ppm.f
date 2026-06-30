module mod_physics_ppm

  implicit none
  public

  procedure(sub_ppm_flatcd), pointer      :: phys_ppm_flatcd => null()
  procedure(sub_ppm_flatsh), pointer      :: phys_ppm_flatsh => null()

  abstract interface
     subroutine sub_ppm_flatcd(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,ixLmax1,&
        ixRmin1,ixRmax1,w,d2w,drho,dp)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,&
          ixLmin1,ixLmax1,ixRmin1,ixRmax1
       double precision, intent(in)    :: w(ixImin1:ixImax1,nw),&
          d2w(ixImin1:ixImax1,1:nwflux)
       double precision, intent(inout) :: drho(ixImin1:ixImax1),&
          dp(ixImin1:ixImax1)
     end subroutine sub_ppm_flatcd

     subroutine sub_ppm_flatsh(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLLmin1,&
        ixLLmax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1,idims,w,&
        drho,dp,dv)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,&
          ixLLmin1,ixLLmax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1
       integer, intent(in)             :: idims
       double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
       double precision, intent(inout) :: drho(ixImin1:ixImax1),&
          dp(ixImin1:ixImax1),dv(ixImin1:ixImax1)
     end subroutine sub_ppm_flatsh
  end interface

contains

  subroutine phys_ppm_check
    if (.not. associated(phys_ppm_flatcd)) phys_ppm_flatcd => dummy_ppm_flatcd

    if (.not. associated(phys_ppm_flatsh)) phys_ppm_flatsh => dummy_ppm_flatsh
  end subroutine phys_ppm_check

  subroutine dummy_ppm_flatcd(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,ixLmax1,&
     ixRmin1,ixRmax1,w,d2w,drho,dp)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,&
       ixLmax1,ixRmin1,ixRmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw),&
       d2w(ixImin1:ixImax1,1:nwflux)
    double precision, intent(inout) :: drho(ixImin1:ixImax1),&
       dp(ixImin1:ixImax1)
    drho(ixOmin1:ixOmax1)=zero
    dp(ixOmin1:ixOmax1)=zero
  end subroutine dummy_ppm_flatcd

  subroutine dummy_ppm_flatsh(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLLmin1,&
     ixLLmax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1,idims,w,drho,&
     dp,dv)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,&
       ixLLmin1,ixLLmax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
    double precision, intent(inout) :: drho(ixImin1:ixImax1),&
       dp(ixImin1:ixImax1),dv(ixImin1:ixImax1)
    drho(ixOmin1:ixOmax1)=zero
    dp(ixOmin1:ixOmax1)=zero
    dv(ixOmin1:ixOmax1)=zero
  end subroutine dummy_ppm_flatsh

end module mod_physics_ppm
