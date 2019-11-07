
module soot_model_module

  use iso_c_binding
  implicit none

  interface
     subroutine fi_init_soot_vars(moms) bind(c)
       import
       implicit none
       real(c_double), intent(inout), dimension(*) :: moms
     end subroutine fi_init_soot_vars
  end interface

end module soot_model_module
