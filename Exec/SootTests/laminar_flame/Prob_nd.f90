subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use amrex_error_module
  use flame_data_module
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(3), probhi(3), L
  character(len=72) :: flame_datafile = ''

  integer untin,i

  namelist /fortin/ flame_datafile, p_ref, T_ref, u_ref
    
  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)
  if (namlen .gt. maxlen) then
     call amrex_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  p_ref = 1.01325d6
  u_ref = 0.
  T_ref = 300.d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  L = probhi(1) - problo(1)
  call initialize_flame_data(flame_datafile, L)
end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
  subroutine pc_initdata(level,time,lo,hi,nvar, &
       state,state_lo,state_hi, &
       delta,xlo,xhi) bind(C, name="pc_initdata")

  use probdata_module
  use network, only: nspecies
  use chemistry_module, only : get_species_index
  use eos_type_module
  use meth_params_module, only : URHO, UMX, UMY, UMZ, &
       UEDEN, UEINT, UFS, UTEMP, nsoot, UFSOOT
  use amrex_constants_module, only: ZERO, HALF, M_PI
  use prob_params_module, only: problo, probhi
  use soot_model_module
  use eos_module
  
  implicit none

  integer :: level, nvar
  integer :: lo(3), hi(3)
  integer :: state_lo(3),state_hi(3)
  double precision :: xlo(3), xhi(3), time, delta(3)
  double precision :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),nvar)

  integer :: i,j,k,isoot,sootcomp
  double precision :: x,y,z,rho,u,v,w,eint
  double precision :: L, p0
  double precision :: moms(7)
  integer, parameter :: out_unit=20
  double precision, allocatable :: flame_Y(:)
  double precision :: flame_T, flame_U, flame_Rho
  
  type(eos_t) :: eos_state

  call build(eos_state)
  allocate(flame_Y(nspecies))

  ! Initialize soot moments to zero
  call fi_init_soot_vars(moms)
  do isoot = 1, nsoot
     sootcomp = UFSOOT + isoot - 1
     state(:,:,:,sootcomp) = moms(isoot)
  enddo

  do k = lo(3), hi(3)
     z = (k+HALF)*delta(3)

     do j = lo(2), hi(2)
        y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

        do i = lo(1), hi(1)
           x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

           call flame(x, flame_Rho, flame_T, flame_U, flame_Y)
           eos_state % massfrac(1:nspecies) = flame_Y(1:nspecies)
           u = flame_U
           v = 0.d0
           w = 0.d0

           eos_state % rho =  flame_rho
           eos_state % T = flame_T
           ! Call EOS by specifying the temperature and density
           call eos_rt(eos_state)
           eint = eos_state % e
           rho = eos_state % rho

           ! Fill the states
           state(i,j,k,URHO)            = rho
           state(i,j,k,UFS:UFS+nspecies-1) = rho * eos_state % massfrac(1:nspecies)
           state(i,j,k,UMX)             = rho * u
           state(i,j,k,UMY)             = rho * v
           state(i,j,k,UMZ)             = rho * w
           state(i,j,k,UEINT)           = rho * eint
           state(i,j,k,UEDEN)           = rho * (eint + HALF * (u**2 + v**2 + w**2))
           state(i,j,k,UTEMP)           = flame_T
        enddo
     enddo
  enddo

end subroutine pc_initdata


subroutine pc_prob_close() &
     bind(C, name="pc_prob_close")

end subroutine pc_prob_close
