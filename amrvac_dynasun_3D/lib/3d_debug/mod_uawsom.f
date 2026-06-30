module mod_uawsom
  use mod_uawsom_phys
  use mod_amrvac

  implicit none
  public

contains

  subroutine uawsom_activate()
    call uawsom_phys_init()
  end subroutine uawsom_activate

end module mod_uawsom
