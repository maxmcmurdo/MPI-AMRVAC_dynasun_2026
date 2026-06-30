module mod_twofl_ppm

#include "amrvac.h"
  use mod_twofl_phys

  implicit none
  private

  public :: twofl_ppm_init

contains

  subroutine twofl_ppm_init()
    use mod_physics_ppm
    phys_ppm_flatcd => twofl_ppm_flatcd
    phys_ppm_flatsh => twofl_ppm_flatsh
  end subroutine twofl_ppm_init

  subroutine twofl_ppm_flatcd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,&
     ixRmax2,w,d2w,drho,dp)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,&
       ixRmax1,ixRmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw),&
       d2w(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux)
    double precision, intent(inout) :: drho(ixImin1:ixImax1,ixImin2:ixImax2),&
       dp(ixImin1:ixImax1,ixImin2:ixImax2)

    if(twofl_eq_energy/=0) then
      drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
         =twofl_gamma*dabs(d2w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_c_))/min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,rho_c_),&
         w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,rho_c_))
      dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dabs(d2w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,e_c_))/min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,e_c_),&
         w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,e_c_))
    else
      call mpistop&
         ("PPM with flatcd=.true. can not be used without energy equation!")
    end if
  end subroutine twofl_ppm_flatcd

  ! based on Mignone and Miller and Collela 2002
  ! PPM flattening at shocks: we use total pressure and not thermal pressure
  subroutine twofl_ppm_flatsh(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,ixLmin1,ixLmin2,&
     ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,ixRRmin1,ixRRmin2,&
     ixRRmax1,ixRRmax2,idims,w,drho,dp,dv)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,ixLmin1,&
       ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,ixRRmin1,&
       ixRRmin2,ixRRmax1,ixRRmax2
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(inout) :: drho(ixImin1:ixImax1,ixImin2:ixImax2),&
       dp(ixImin1:ixImax1,ixImin2:ixImax2),dv(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                :: ptot(ixImin1:ixImax1,ixImin2:ixImax2)

    if(twofl_eq_energy/=0) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_c_)+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))**2,&
         dim=ndim+1)
      where (dabs(ptot(ixRRmin1:ixRRmax1,&
         ixRRmin2:ixRRmax2)-ptot(ixLLmin1:ixLLmax1,&
         ixLLmin2:ixLLmax2))>smalldouble)
         drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dabs((ptot(ixRmin1:ixRmax1,&
            ixRmin2:ixRmax2)-ptot(ixLmin1:ixLmax1,&
            ixLmin2:ixLmax2))/(ptot(ixRRmin1:ixRRmax1,&
            ixRRmin2:ixRRmax2)-ptot(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2)))
      elsewhere
         drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dp" to save squared sound speed, assume primitive in w
      dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(twofl_gamma*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,e_c_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_c_))

      dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = dabs(ptot(ixRmin1:ixRmax1,&
         ixRmin2:ixRmax2)-ptot(ixLmin1:ixLmax1,&
         ixLmin2:ixLmax2))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_c_)*dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      ! recycle ptot to store v
      ptot(ixImin1:ixImax1,ixImin2:ixImax2)= w(ixImin1:ixImax1,ixImin2:ixImax2,&
         mom_c(idims))
      call gradient(ptot,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,idims,dv)
    else
      call mpistop&
         ("PPM with flatsh=.true. can not be used without energy equation!")
    end if

  end subroutine twofl_ppm_flatsh

end module mod_twofl_ppm
