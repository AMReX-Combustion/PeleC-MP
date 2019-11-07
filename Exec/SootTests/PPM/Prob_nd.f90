subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use amrex_error_module
  use flamelet_data_module
  use hit_turb_data_module
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(3), probhi(3)
  character(len=72) :: flamelet_datafile = ''
  character(len=72) :: turbulence_datafile = ''

  integer untin,i

  namelist /fortin/ hz, delz, urms0, inres, flamelet_datafile, &
       turbulence_datafile, p_ref, T_ref, u_ref
    
  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  hz = 1.d0
  delz = 0.2d0
  p_ref = 1.01325d6
  u_ref = 0.
  T_ref = 300.d0
  urms0 = 1.d0
  inres = 0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  call initialize_flamelet(flamelet_datafile)
  call initialize_hit_turb(turbulence_datafile)

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
  use hit_turb_data_module
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
  integer :: n, np1, m, mp1
  double precision :: x,y,z,rho,u,v,w,eint
  double precision :: sval, rval
  double precision :: L, p0, zmix
  double precision :: moms(7)
  integer, parameter :: out_unit=20
  double precision, allocatable :: flamelet_Y(:)
  double precision :: flamelet_T
  
  type(eos_t) :: eos_state
  
  integer :: iN2, iO2, iFuel

  iN2 = get_species_index("N2")
  iO2 = get_species_index("O2")
  iFuel = get_species_index("N-C7H16")

  call build(eos_state)
  allocate(flamelet_Y(nspecies))

  ! Define the length
  L = probhi(2) - problo(2)

  ! Initial pressure and temperature
  p0 = p_ref ! [erg cm^-3]

  ! Set the equation of state variables
  eos_state % p = p0
  eos_state % T = 300.d0
  eos_state % massfrac    = 0.d0
  eos_state % massfrac(iN2) = 0.767
  eos_state % massfrac(iO2) = one - 0.767

  call eos_tp(eos_state)

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
        call locate(y, n, np1, sval)

        do i = lo(1), hi(1)
           x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)
           call locate(x, m, mp1, rval)

           zmix = 1.d0/4.d0*(1.d0 + tanh((L + hz + delz - 2.d0*y)/delz))* &
                (1.d0 + tanh((-L + hz + delz + 2.d0*y)/delz))
           call flamelet(zmix,flamelet_Y,flamelet_T)
           eos_state % massfrac(1:nspecies) = flamelet_Y(1:nspecies)
           call interpvel(rval, sval, m, mp1, n, np1, u, v)
           w = 0.d0

           eos_state % p =  p0
           eos_state % T = flamelet_T
           ! Call EOS by specifying the temperature and pressure
           call eos_tp(eos_state)
           rho  = eos_state % rho
           eint = eos_state % e

           ! Fill the states
           state(i,j,k,URHO)            = rho
           state(i,j,k,UFS:UFS+nspecies-1) = rho * eos_state % massfrac(1:nspecies)
           state(i,j,k,UMX)             = rho * u
           state(i,j,k,UMY)             = rho * v
           state(i,j,k,UMZ)             = rho * w
           state(i,j,k,UEINT)           = rho * eint
           state(i,j,k,UEDEN)           = rho * (eint + HALF * (u**2 + v**2 + w**2))
           state(i,j,k,UTEMP)           = flamelet_T
        enddo
     enddo
  enddo

end subroutine pc_initdata


subroutine pc_prob_close() &
     bind(C, name="pc_prob_close")

end subroutine pc_prob_close
