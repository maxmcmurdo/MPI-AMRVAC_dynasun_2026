module mod_physics_ppm

  implicit none
  public

  procedure(sub_ppm_flatcd), pointer      :: phys_ppm_flatcd => null()
  procedure(sub_ppm_flatsh), pointer      :: phys_ppm_flatsh => null()

  abstract interface
     subroutine sub_ppm_flatcd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLmin1,ixLmin2,&
        ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
        ixRmax2,ixRmax3,w,d2w,drho,dp)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,&
          ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,&
          ixRmin3,ixRmax1,ixRmax2,ixRmax3
       double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
          ixImin3:ixImax3,nw),d2w(ixImin1:ixImax1,ixImin2:ixImax2,&
          ixImin3:ixImax3,1:nwflux)
       double precision, intent(inout) :: drho(ixImin1:ixImax1,ixImin2:ixImax2,&
          ixImin3:ixImax3),dp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
     end subroutine sub_ppm_flatcd

     subroutine sub_ppm_flatsh(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLLmin1,ixLLmin2,&
        ixLLmin3,ixLLmax1,ixLLmax2,ixLLmax3,ixLmin1,ixLmin2,ixLmin3,ixLmax1,&
        ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
        ixRRmin1,ixRRmin2,ixRRmin3,ixRRmax1,ixRRmax2,ixRRmax3,idims,w,drho,dp,&
        dv)
       use mod_global_parameters
       integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,&
          ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
          ixLLmin1,ixLLmin2,ixLLmin3,ixLLmax1,ixLLmax2,ixLLmax3,ixLmin1,&
          ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,&
          ixRmax1,ixRmax2,ixRmax3,ixRRmin1,ixRRmin2,ixRRmin3,ixRRmax1,ixRRmax2,&
          ixRRmax3
       integer, intent(in)             :: idims
       double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
          ixImin3:ixImax3,nw)
       double precision, intent(inout) :: drho(ixImin1:ixImax1,ixImin2:ixImax2,&
          ixImin3:ixImax3),dp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
          dv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
     end subroutine sub_ppm_flatsh
  end interface

contains

  subroutine phys_ppm_check
    if (.not. associated(phys_ppm_flatcd)) phys_ppm_flatcd => dummy_ppm_flatcd

    if (.not. associated(phys_ppm_flatsh)) phys_ppm_flatsh => dummy_ppm_flatsh
  end subroutine phys_ppm_check

  subroutine dummy_ppm_flatcd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLmin1,ixLmin2,ixLmin3,&
     ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,w,&
     d2w,drho,dp)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLmin1,ixLmin2,&
       ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
       ixRmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw),d2w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:nwflux)
    double precision, intent(inout) :: drho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),dp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
  end subroutine dummy_ppm_flatcd

  subroutine dummy_ppm_flatsh(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLLmin1,ixLLmin2,&
     ixLLmin3,ixLLmax1,ixLLmax2,ixLLmax3,ixLmin1,ixLmin2,ixLmin3,ixLmax1,&
     ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixRRmin1,&
     ixRRmin2,ixRRmin3,ixRRmax1,ixRRmax2,ixRRmax3,idims,w,drho,dp,dv)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLLmin1,&
       ixLLmin2,ixLLmin3,ixLLmax1,ixLLmax2,ixLLmax3,ixLmin1,ixLmin2,ixLmin3,&
       ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
       ixRRmin1,ixRRmin2,ixRRmin3,ixRRmax1,ixRRmax2,ixRRmax3
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(inout) :: drho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),dp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       dv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    dv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
  end subroutine dummy_ppm_flatsh

end module mod_physics_ppm
