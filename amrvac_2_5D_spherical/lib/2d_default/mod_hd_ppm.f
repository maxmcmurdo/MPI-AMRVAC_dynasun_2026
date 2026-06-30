!> Hydrodynamics PPM module
module mod_hd_ppm
  use mod_hd_phys

  implicit none
  private

  public :: hd_ppm_init

contains

  subroutine hd_ppm_init()
    use mod_physics_ppm

    phys_ppm_flatcd => hd_ppm_flatcd
    phys_ppm_flatsh => hd_ppm_flatsh
  end subroutine hd_ppm_init

  subroutine hd_ppm_flatcd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,&
     ixRmax2,w,d2w,drho,dp)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixLmin1,ixLmin2,ixLmax1,ixLmax2,&
        ixRmin1,ixRmin2,ixRmax1,ixRmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        d2w(ixGlo1:ixGhi1,ixGlo2:ixGhi2, 1:nwflux)
    double precision, intent(inout) :: drho(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
        dp(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

    if(hd_energy) then
      drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = hd_gamma*abs(d2w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, rho_))/min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2, rho_),&
          w(ixRmin1:ixRmax1,ixRmin2:ixRmax2, rho_))
      dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs(d2w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, e_))/min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2, e_),&
          w(ixRmin1:ixRmax1,ixRmin2:ixRmax2, e_))
    else
      call mpistop("PPM with flatcd=.true. can not be used with -eos = iso !")
    end if
  end subroutine hd_ppm_flatcd

  subroutine hd_ppm_flatsh(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,ixLmin1,ixLmin2,&
     ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,ixRRmin1,ixRRmin2,&
     ixRRmax1,ixRRmax2,idims,w,drho,dp,dv)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,&
        ixLmin1,ixLmin2,ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
        ixRRmin1,ixRRmin2,ixRRmax1,ixRRmax2
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(inout) :: drho(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
        dp(ixGlo1:ixGhi1,ixGlo2:ixGhi2), dv(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
    double precision                :: v(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

    if(hd_energy) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      where (abs(w(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,&
          e_)-w(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2, e_))>smalldouble)
         drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs((w(ixRmin1:ixRmax1,&
            ixRmin2:ixRmax2, e_)-w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
             e_))/(w(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,&
             e_)-w(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2, e_)))
      else where
         drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dp" to save squared sound speed, assuming primitives
      dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =(hd_gamma*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, e_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_))

      dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs(w(ixRmin1:ixRmax1,&
         ixRmin2:ixRmax2, e_)-w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
          e_))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)*dp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))
      v(ixImin1:ixImax1,ixImin2:ixImax2)  = w(ixImin1:ixImax1,ixImin2:ixImax2,&
          mom(idims))
      call gradient(v, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, idims, dv)
    else
      call mpistop("PPM with flatsh=.true. can not be used with -eos = iso !")
    end if
  end subroutine hd_ppm_flatsh

end module mod_hd_ppm
