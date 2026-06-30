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

  subroutine twofl_ppm_flatcd(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,ixLmax1,&
     ixRmin1,ixRmax1,w,d2w,drho,dp)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,&
       ixLmax1,ixRmin1,ixRmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw),&
       d2w(ixImin1:ixImax1,1:nwflux)
    double precision, intent(inout) :: drho(ixImin1:ixImax1),&
       dp(ixImin1:ixImax1)

    if(twofl_eq_energy/=0) then
      drho(ixOmin1:ixOmax1) =twofl_gamma*dabs(d2w(ixOmin1:ixOmax1,&
         rho_c_))/min(w(ixLmin1:ixLmax1,rho_c_),w(ixRmin1:ixRmax1,rho_c_))
      dp(ixOmin1:ixOmax1) = dabs(d2w(ixOmin1:ixOmax1,&
         e_c_))/min(w(ixLmin1:ixLmax1,e_c_),w(ixRmin1:ixRmax1,e_c_))
    else
      call mpistop&
         ("PPM with flatcd=.true. can not be used without energy equation!")
    end if
  end subroutine twofl_ppm_flatcd

  ! based on Mignone and Miller and Collela 2002
  ! PPM flattening at shocks: we use total pressure and not thermal pressure
  subroutine twofl_ppm_flatsh(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLLmin1,&
     ixLLmax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1,idims,w,drho,&
     dp,dv)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,&
       ixLLmin1,ixLLmax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,nw)
    double precision, intent(inout) :: drho(ixImin1:ixImax1),&
       dp(ixImin1:ixImax1),dv(ixImin1:ixImax1)
    double precision                :: ptot(ixImin1:ixImax1)

    if(twofl_eq_energy/=0) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      ptot(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,e_c_)+half*sum(w(ixOmin1:ixOmax1,&
         mag(:))**2,dim=ndim+1)
      where (dabs(ptot(ixRRmin1:ixRRmax1)-&
         ptot(ixLLmin1:ixLLmax1))>smalldouble)
         drho(ixOmin1:ixOmax1) = dabs((ptot(ixRmin1:ixRmax1)-&
            ptot(ixLmin1:ixLmax1))/(ptot(ixRRmin1:ixRRmax1)-&
            ptot(ixLLmin1:ixLLmax1)))
      elsewhere
         drho(ixOmin1:ixOmax1) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dp" to save squared sound speed, assume primitive in w
      dp(ixOmin1:ixOmax1)=(twofl_gamma*w(ixOmin1:ixOmax1,&
         e_c_)/w(ixOmin1:ixOmax1,rho_c_))

      dp(ixOmin1:ixOmax1)  = dabs(ptot(ixRmin1:ixRmax1)-&
         ptot(ixLmin1:ixLmax1))/(w(ixOmin1:ixOmax1,&
         rho_c_)*dp(ixOmin1:ixOmax1))
      ! recycle ptot to store v
      ptot(ixImin1:ixImax1)= w(ixImin1:ixImax1,mom_c(idims))
      call gradient(ptot,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,dv)
    else
      call mpistop&
         ("PPM with flatsh=.true. can not be used without energy equation!")
    end if

  end subroutine twofl_ppm_flatsh

end module mod_twofl_ppm
