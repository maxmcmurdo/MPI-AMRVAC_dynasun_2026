!> Subroutines for Roe-type Riemann solver for HD
module mod_hyperdiffusivity


  implicit none

  !all public
  !private
  !public :: hyperdiffusivity_init


contains

  subroutine hyperdiffusivity_init()
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: phys_req_diagonal

    !print*, slab_uniform !!THIS IS FALSE, why??

    !if(coordinate .ne. Cartesian .or. .not. slab_uniform) then
    if(coordinate .ne. Cartesian) then
      call mpistop&
         ("Hyperdiffusivity only implemented for Cartesian uniform grid")
    endif

    nghostcells = max(nghostcells,3)
    phys_req_diagonal = .true.

  end subroutine hyperdiffusivity_init


  subroutine hyp_coeff(ixImin1,ixImax1, ixOmin1,ixOmax1, var, idimm, nu_hyp)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImax1, idimm
    integer, intent(out)                :: ixOmin1,ixOmax1
    double precision, intent(in)       :: var(ixImin1:ixImax1)
    double precision, intent(out)       :: nu_hyp(ixImin1:ixImax1)

    double precision        :: tmp(ixImin1:ixImax1), tmp2(ixImin1:ixImax1)
    integer :: hxmin1,hxmax1, hx1fmin1,hx1fmax1, hx1bmin1,hx1bmax1, hx2bmin1,&
       hx2bmax1

    double precision, parameter :: eps=1e-8


   hxmin1=ixImin1+2*kr(idimm,1); 
   hxmax1=ixImax1-kr(idimm,1); 

   hx1fmin1=hxmin1+kr(idimm,1);hx1fmax1=hxmax1+kr(idimm,1);  
   hx1bmin1=hxmin1-kr(idimm,1);hx1bmax1=hxmax1-kr(idimm,1);  
   hx2bmin1=hxmin1-2*kr(idimm,1);hx2bmax1=hxmax1-2*kr(idimm,1);  

  !store d3 in tmp
  tmp(hxmin1:hxmax1)= abs(3d0 * (var(hxmin1:hxmax1) - var(hx1bmin1:hx1bmax1)) &
     - (var(hx1fmin1:hx1fmax1)-var(hx2bmin1:hx2bmax1)))
  !store d1 in tmp2
  tmp2(hxmin1:hxmax1)= abs((var(hxmin1:hxmax1) - var(hx1bmin1:hx1bmax1)))

  ixOmin1=hxmin1+kr(idimm,1);
  ixOmax1=hxmax1-kr(idimm,1);
  hx1fmin1=ixOmin1+kr(idimm,1);hx1fmax1=ixOmax1+kr(idimm,1);  
  hx1bmin1=ixOmin1-kr(idimm,1);hx1bmax1=ixOmax1-kr(idimm,1);  

  nu_hyp(ixOmin1:ixOmax1) = max(tmp(hx1bmin1:hx1bmax1),tmp(ixOmin1:ixOmax1),&
     tmp(hx1fmin1:hx1fmax1))/max(tmp2(hx1bmin1:hx1bmax1),tmp2(ixOmin1:ixOmax1),&
     tmp2(hx1fmin1:hx1fmax1),eps)  

  !print*, "HYP IXO ", ixO^L 

  end subroutine hyp_coeff



  subroutine div_vel_coeff(ixImin1,ixImax1, ixOmin1,ixOmax1, vel, idimm,&
      nu_vel)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImax1, idimm
    integer, intent(out)                :: ixOmin1,ixOmax1
    double precision, intent(in)       :: vel(ixImin1:ixImax1,1:ndir)
    double precision, intent(out)       :: nu_vel(ixImin1:ixImax1)
    integer :: hx1fmin1,hx1fmax1, hx1bmin1,hx1bmax1

    integer :: ii

   ixOmin1=ixImin1+1; 
   ixOmax1=ixImax1-1; 

   ixOmin1=ixOmin1+kr(idimm,1); 

   nu_vel(ixOmin1:ixOmax1) = 0d0 

  ! ndim should always be smaller or equal to ndir !
   do ii = 1, ndim
    hx1bmin1=ixOmin1-kr(ii,1);hx1bmax1=ixOmax1-kr(ii,1);
    if(ii==idimm) then
      nu_vel(ixOmin1:ixOmax1) = nu_vel(ixOmin1:ixOmax1) + (vel(ixOmin1:ixOmax1,&
         ii) - vel(hx1bmin1:hx1bmax1,ii))/dxlevel(ii) 
    else
      hx1fmin1=ixOmin1+kr(ii,1);hx1fmax1=ixOmax1+kr(ii,1);  
      nu_vel(ixOmin1:ixOmax1) = nu_vel(ixOmin1:ixOmax1) + &
         (vel(hx1fmin1:hx1fmax1,ii) - vel(hx1bmin1:hx1bmax1,&
         ii))/(2d0*dxlevel(ii)) 
    endif

    where(nu_vel(ixOmin1:ixOmax1) < 0d0)
      nu_vel(ixOmin1:ixOmax1) = -nu_vel(ixOmin1:ixOmax1)
    elsewhere
      nu_vel(ixOmin1:ixOmax1) = 0d0
    endwhere

   enddo 

  !print*, "DIV_VEL IXO ", ixO^L 

  end subroutine div_vel_coeff



  !var has cell centered values
  !nu_hyper is defined at the interfaces
  subroutine second_same_deriv(ixImin1,ixImax1, ixOmin1,ixOmax1, nu_hyper, var,&
      idimm, res)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImax1, idimm
    integer, intent(out)                :: ixOmin1,ixOmax1
    double precision, intent(in)       :: nu_hyper(ixImin1:ixImax1),&
       var(ixImin1:ixImax1)
    double precision, intent(out)       :: res(ixImin1:ixImax1)

    integer :: hxfmin1,hxfmax1, hxbmin1,hxbmax1

    ixOmin1=ixImin1+3;
    ixOmax1=ixImax1-3;

    hxfmin1=ixOmin1+kr(idimm,1);hxfmax1=ixOmax1+kr(idimm,1);
    hxbmin1=ixOmin1-kr(idimm,1);hxbmax1=ixOmax1-kr(idimm,1);

    res(ixOmin1:ixOmax1) = 1d0/(dxlevel(idimm)**2)*(nu_hyper(hxfmin1:hxfmax1) &
       * (var(hxfmin1:hxfmax1)-var(ixOmin1:ixOmax1))-nu_hyper(ixOmin1:ixOmax1) &
       * (var(ixOmin1:ixOmax1)-var(hxbmin1:hxbmax1)))
   !print*, "SECOND SAME DERIV IXO ", ixO^L 

  end subroutine second_same_deriv

  !var has cell centered values
  !var2 has cell centered values
  !nu_hyper is defined at the interfaces
  subroutine second_same_deriv2(ixImin1,ixImax1, ixOmin1,ixOmax1, nu_hyper,&
      var2, var, idimm, res)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImax1, idimm
    integer, intent(out)                :: ixOmin1,ixOmax1
    double precision, intent(in)       :: nu_hyper(ixImin1:ixImax1),&
       var(ixImin1:ixImax1), var2(ixImin1:ixImax1)
    double precision, intent(out)       :: res(ixImin1:ixImax1)

    integer :: hxfmin1,hxfmax1, hxbmin1,hxbmax1

    ixOmin1=ixImin1+3;
    ixOmax1=ixImax1-3;

    hxfmin1=ixOmin1+kr(idimm,1);hxfmax1=ixOmax1+kr(idimm,1);
    hxbmin1=ixOmin1-kr(idimm,1);hxbmax1=ixOmax1-kr(idimm,1);

    res(ixOmin1:ixOmax1) = 1d0/(2d0*dxlevel(idimm)**&
       2)*(nu_hyper(hxfmin1:hxfmax1) * (var2(hxfmin1:hxfmax1)+&
       var2(ixOmin1:ixOmax1)) *  (var(hxfmin1:hxfmax1)-var(ixOmin1:ixOmax1))-&
       nu_hyper(ixOmin1:ixOmax1) * (var2(hxbmin1:hxbmax1)+&
       var2(ixOmin1:ixOmax1)) * (var(ixOmin1:ixOmax1)-var(hxbmin1:hxbmax1)))
    !print*, "SECOND SAME DERIV2 IXO ", ixO^L 

  end subroutine second_same_deriv2

  !idimm inner derivative, idimm2 outer
  ! deriv_idimm2(nu * deriv_idimm (u) )
  !var has cell centered values
  !nu_hyper is defined at the interfaces
  subroutine second_cross_deriv(ixImin1,ixImax1, ixOmin1,ixOmax1, nu_hyper,&
      var, idimm, idimm2, res)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImax1, idimm,idimm2
    integer, intent(out)                :: ixOmin1,ixOmax1
    double precision, intent(in)       :: nu_hyper(ixImin1:ixImax1),&
       var(ixImin1:ixImax1)
    double precision, intent(out)       :: res(ixImin1:ixImax1)

    integer :: hxfimin1,hxfimax1, hxbimin1,hxbimax1, hxfifjmin1,hxfifjmax1,&
        hxbibjmin1,hxbibjmax1, hxfibjmin1,hxfibjmax1, hxbifjmin1,hxbifjmax1

    ixOmin1=ixImin1+3;
    ixOmax1=ixImax1-3;

    hxfimin1=ixOmin1+kr(idimm2,1);hxfimax1=ixOmax1+kr(idimm2,1);
    hxbimin1=ixOmin1-kr(idimm2,1);hxbimax1=ixOmax1-kr(idimm2,1);

    hxfifjmin1=hxfimin1+kr(idimm,1);hxfifjmax1=hxfimax1+kr(idimm,1);
    hxfibjmin1=hxfimin1-kr(idimm,1);hxfibjmax1=hxfimax1-kr(idimm,1);

    hxbifjmin1=hxbimin1+kr(idimm,1);hxbifjmax1=hxbimax1+kr(idimm,1);
    hxbibjmin1=hxbimin1-kr(idimm,1);hxbibjmax1=hxbimax1-kr(idimm,1);

    res(ixOmin1:ixOmax1) = 1d0/(8d0*dxlevel(idimm) * &
       dxlevel(idimm2))*((nu_hyper(hxfifjmin1:hxfifjmax1) + &
       nu_hyper(hxfimin1:hxfimax1)) * (var(hxfifjmin1:hxfifjmax1)-&
       var(hxfibjmin1:hxfibjmax1))-(nu_hyper(hxbifjmin1:hxbifjmax1) + &
       nu_hyper(hxbimin1:hxbimax1)) * (var(hxbifjmin1:hxbifjmax1)-&
       var(hxbibjmin1:hxbibjmax1)))

    !print*, "SECOND CROSS DERIV IXO ", ixO^L 

  end subroutine second_cross_deriv


  !idimm inner derivative, idimm2 outer
  ! deriv_idimm2(nu * var2 * deriv_idimm (u) )
  !var has cell centered values
  !var2 has cell centered values
  !nu_hyper is defined at the interfaces (with numbering index left center): center_i-1 interface_i center_i
  subroutine second_cross_deriv2(ixImin1,ixImax1, ixOmin1,ixOmax1, nu_hyper,&
      var2, var, idimm, idimm2, res)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImax1, idimm,idimm2
    integer, intent(out)                :: ixOmin1,ixOmax1
    double precision, intent(in)       :: nu_hyper(ixImin1:ixImax1),&
       var(ixImin1:ixImax1),var2(ixImin1:ixImax1)
    double precision, intent(out)       :: res(ixImin1:ixImax1)

    integer :: hxfimin1,hxfimax1, hxbimin1,hxbimax1, hxfifjmin1,hxfifjmax1,&
        hxbibjmin1,hxbibjmax1, hxfibjmin1,hxfibjmax1, hxbifjmin1,hxbifjmax1

    ixOmin1=ixImin1+3;
    ixOmax1=ixImax1-3;

    hxfimin1=ixOmin1+kr(idimm2,1);hxfimax1=ixOmax1+kr(idimm2,1);
    hxbimin1=ixOmin1-kr(idimm2,1);hxbimax1=ixOmax1-kr(idimm2,1);

    hxfifjmin1=hxfimin1+kr(idimm,1);hxfifjmax1=hxfimax1+kr(idimm,1);
    hxfibjmin1=hxfimin1-kr(idimm,1);hxfibjmax1=hxfimax1-kr(idimm,1);

    hxbifjmin1=hxbimin1+kr(idimm,1);hxbifjmax1=hxbimax1+kr(idimm,1);
    hxbibjmin1=hxbimin1-kr(idimm,1);hxbibjmax1=hxbimax1-kr(idimm,1);

    res(ixOmin1:ixOmax1) = 1d0/(8d0*dxlevel(idimm) * &
       dxlevel(idimm2))*((nu_hyper(hxfifjmin1:hxfifjmax1) + &
       nu_hyper(hxfimin1:hxfimax1)) * var2(hxfimin1:hxfimax1) * &
       (var(hxfifjmin1:hxfifjmax1)-var(hxfibjmin1:hxfibjmax1))-&
       (nu_hyper(hxbifjmin1:hxbifjmax1) + nu_hyper(hxbimin1:hxbimax1)) * &
       var2(hxbimin1:hxbimax1) * (var(hxbifjmin1:hxbifjmax1)-&
       var(hxbibjmin1:hxbibjmax1)))

    !print*, "SECOND CROSS DERIV2 IXO ", ixO^L 

  end subroutine second_cross_deriv2

end module mod_hyperdiffusivity
