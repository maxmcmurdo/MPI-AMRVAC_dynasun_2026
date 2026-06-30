module mod_mhd_ppm
  use mod_mhd_phys

  implicit none
  private

  public :: mhd_ppm_init

contains

  subroutine mhd_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => mhd_ppm_flatcd
    phys_ppm_flatsh => mhd_ppm_flatsh
  end subroutine mhd_ppm_init

  subroutine mhd_ppm_flatcd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
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

    if(mhd_energy) then
      drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) =mhd_gamma*dabs(d2w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_))/min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
         ixLmin3:ixLmax3,rho_),w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
         ixRmin3:ixRmax3,rho_))
      dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = dabs(d2w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,p_))/min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
         ixLmin3:ixLmax3,p_),w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
         p_))
    else
      call mpistop&
         ("PPM with flatcd=.true. can not be used without energy equation!")
    end if
  end subroutine mhd_ppm_flatcd

  ! based on Mignone and Miller and Collela 2002
  ! PPM flattening at shocks: we use total pressure and not thermal pressure
  subroutine mhd_ppm_flatsh(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ixLLmin1,ixLLmin2,&
     ixLLmin3,ixLLmax1,ixLLmax2,ixLLmax3,ixLmin1,ixLmin2,ixLmin3,ixLmax1,&
     ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixRRmin1,&
     ixRRmin2,ixRRmin3,ixRRmax1,ixRRmax2,ixRRmax3,idims,w,drho,dp,dv)
    use mod_global_parameters
    use mod_geometry

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
    double precision                :: ptot(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(mhd_energy) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)+half*sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)
      where (dabs(ptot(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,&
         ixRRmin3:ixRRmax3)-ptot(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2,&
         ixLLmin3:ixLLmax3))>smalldouble)
         drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) = dabs((ptot(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
            ixRmin3:ixRmax3)-ptot(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
            ixLmin3:ixLmax3))/(ptot(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,&
            ixRRmin3:ixRRmax3)-ptot(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2,&
            ixLLmin3:ixLLmax3)))
      elsewhere
         drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dp" to save squared sound speed, assume primitive in w
      dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(mhd_gamma*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,p_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_))

      dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)  = dabs(ptot(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
         ixRmin3:ixRmax3)-ptot(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
         ixLmin3:ixLmax3))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)*dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      ! recycle ptot to store v
      ptot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)= w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,mom(idims))
      call gradient(ptot,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,dv)
    else
      call mpistop&
         ("PPM with flatsh=.true. can not be used without energy equation!")
    end if

  end subroutine mhd_ppm_flatsh

end module mod_mhd_ppm
