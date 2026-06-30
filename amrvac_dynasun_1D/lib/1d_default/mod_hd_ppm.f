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

  subroutine hd_ppm_flatcd(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,ixLmax1,&
     ixRmin1,ixRmax1,w,d2w,drho,dp)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
        ixLmin1,ixLmax1, ixRmin1,ixRmax1
    double precision, intent(in)    :: w(ixImin1:ixImax1, nw),&
        d2w(ixGlo1:ixGhi1, 1:nwflux)
    double precision, intent(inout) :: drho(ixGlo1:ixGhi1), dp(ixGlo1:ixGhi1)

    if(hd_energy) then
      drho(ixOmin1:ixOmax1) = hd_gamma*abs(d2w(ixOmin1:ixOmax1,&
          rho_))/min(w(ixLmin1:ixLmax1, rho_), w(ixRmin1:ixRmax1, rho_))
      dp(ixOmin1:ixOmax1) = abs(d2w(ixOmin1:ixOmax1,&
          e_))/min(w(ixLmin1:ixLmax1, e_), w(ixRmin1:ixRmax1, e_))
    else
      call mpistop("PPM with flatcd=.true. can not be used with -eos = iso !")
    end if
  end subroutine hd_ppm_flatcd

  subroutine hd_ppm_flatsh(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLLmin1,ixLLmax1,&
     ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1,idims,w,drho,dp,dv)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
        ixLLmin1,ixLLmax1, ixLmin1,ixLmax1, ixRmin1,ixRmax1, ixRRmin1,ixRRmax1
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixImin1:ixImax1, nw)
    double precision, intent(inout) :: drho(ixGlo1:ixGhi1), dp(ixGlo1:ixGhi1),&
        dv(ixGlo1:ixGhi1)
    double precision                :: v(ixGlo1:ixGhi1)

    if(hd_energy) then
      ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
      where (abs(w(ixRRmin1:ixRRmax1, e_)-w(ixLLmin1:ixLLmax1,&
          e_))>smalldouble)
         drho(ixOmin1:ixOmax1) = abs((w(ixRmin1:ixRmax1, e_)-w(ixLmin1:ixLmax1,&
             e_))/(w(ixRRmin1:ixRRmax1, e_)-w(ixLLmin1:ixLLmax1, e_)))
      else where
         drho(ixOmin1:ixOmax1) = zero
      end where

      !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
      !  use "dp" to save squared sound speed, assuming primitives
      dp(ixOmin1:ixOmax1) =(hd_gamma*w(ixOmin1:ixOmax1, e_)/w(ixOmin1:ixOmax1,&
          rho_))

      dp(ixOmin1:ixOmax1) = abs(w(ixRmin1:ixRmax1, e_)-w(ixLmin1:ixLmax1,&
          e_))/(w(ixOmin1:ixOmax1, rho_)*dp(ixOmin1:ixOmax1))
      v(ixImin1:ixImax1)  = w(ixImin1:ixImax1, mom(idims))
      call gradient(v, ixImin1,ixImax1, ixOmin1,ixOmax1, idims, dv)
    else
      call mpistop("PPM with flatsh=.true. can not be used with -eos = iso !")
    end if
  end subroutine hd_ppm_flatsh

end module mod_hd_ppm
