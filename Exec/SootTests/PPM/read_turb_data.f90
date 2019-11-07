module hit_turb_data_module
  implicit none
  character(len=72) :: turbulence_filename = ''
  integer :: turb_init = 0
  integer(kind=8) :: nx, ny
  double precision, dimension(:,:), allocatable :: uinput, vinput
  double precision, dimension(:), allocatable :: xarray, xdiff
  double precision :: Linput

contains

  subroutine read_turb()

    use amrex_constants_module, only: zero, one, half
    use prob_params_module, only: dim
    use probdata_module
    integer :: i
    double precision, dimension(:), allocatable :: data
    double precision, dimension(:,:), allocatable :: xinput

    nx = int8(inres)
    ny = int8(1)
    if (dim .ge. 2) then
       ny = int8(inres)
       if (dim .ge. 3) then
          call bl_abort('must be 2D')
       endif
    endif
    allocate(data(0:nx*ny*4-1))
    allocate(xinput(0:nx-1,0:ny-1))
    allocate(uinput(0:nx-1,0:ny-1))
    allocate(vinput(0:nx-1,0:ny-1))
    call read_csv(turbulence_filename, nx, ny, data)
    uinput = urms0 * reshape(data(2::4), (/nx, ny/))
    vinput = urms0 * reshape(data(3::4), (/nx, ny/))
    xinput = reshape(data(0::4), (/nx, ny/))
    allocate(xarray(0:nx-1))
    allocate(xdiff(0:nx-1))
    xarray(0:nx-1) = xinput(:,0)
    xdiff(:nx-2) = xarray(1:) - xarray(:nx-2)
    xdiff(nx-1) = xarray(nx-1) - xarray(nx-2)
    Linput = maxval(xinput(:,0)) + HALF*xdiff(nx-1)

    ! Deallocate some stuff
    deallocate(data)
    deallocate(xinput)
  end subroutine read_turb

  ! ::: -----------------------------------------------------------
  ! ::: Read a csv file
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: iname => filename
  ! ::: nx    => input resolution
  ! ::: ny    => input resolution
  ! ::: nz    => input resolution
  ! ::: data  <= output data
  ! ::: -----------------------------------------------------------
  subroutine read_csv(iname,nx,ny,data)

    implicit none

    character(len=72), intent(in) :: iname
    integer(kind=8), intent(in) :: nx, ny
    double precision, intent(out) :: data(0:nx*ny*4-1)

    integer :: i
    integer, parameter :: in_unit=1
    integer :: nlines = 0, ios = 0

    ! Get number of lines in file
    open(in_unit,file=trim(iname), access='sequential', form='formatted', status='old', action='read')
    read(in_unit,*) ! skip header
    nlines = 0
    do
       read(in_unit,*,iostat=ios)
       if (ios .ne. 0) exit
       nlines = nlines + 1
    enddo
    ios = 0

    ! Quick sanity check
    if (nlines .ne. nx*ny) then
       write(*,'("Number of lines in the input file (=",I0,") does not ")')nlines
       write(*,'("  match the input resolution (n=",I0,") in the probin file")')nx
       stop 99
    endif

    ! Read the data from the file
    rewind(in_unit)
    read(in_unit,*) ! skip header
    do i = 0, nlines-1
       read(in_unit, *, iostat=ios)data(i*4:(i+1)*4-1)
       if (ios .ne. 0) then
          write(*,*)'Error in CSV input file read. Exiting with read error', ios
          stop 99
       endif
    enddo
    close(in_unit)

  end subroutine read_csv

  ! ::: -----------------------------------------------------------
  ! ::: Search for the closest index in an array to a given value
  ! ::: using the bisection technique.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: x             => x location
  ! ::: idxlo        <=> output st. xarray(idxlo) <= x < xarray(idxlo+1)
  ! ::: idxho        <=> output mod(idxlo+1,inres)
  ! ::: rr           <=> output (x - xarray(idxlo))/xdiff(idxlo)
  ! ::: -----------------------------------------------------------
  subroutine locate(xin, idxlo, idxhi,  rr)

    use probdata_module
    implicit none
    double precision, intent(out) :: rr
    double precision, intent(in) :: xin
    integer, intent(out) :: idxlo, idxhi

    ! Local variables
    integer :: idxmid
    logical :: notdone
    double precision :: x

    ! If x is out of bounds, return boundary index
    x = mod(xin, Linput)
    if (x >= xarray(inres-1)) then
       idxlo=inres-1
       return
    elseif (x <= xarray(0)) then
       idxlo=0
       return
    endif

    ! Make sure the search array is increasing
    if (xarray(0) > xarray(inres-1)) then
       write(*,'("Error in locate: non ascending input search array.")')
       stop 99
    endif

    ! Do the bisection
    idxlo = 0
    idxhi = inres-1
    notdone = .true.
    do while (notdone)
       if (idxhi-idxlo <= 1) then
          notdone = .false.
       else
          idxmid = (idxhi+idxlo)/2
          if (x >= xarray(idxmid)) then
             idxlo = idxmid
          else
             idxhi = idxmid
          endif
       endif
    enddo
    rr = (x - xarray(idxlo))/xdiff(idxlo)
    idxhi = mod(idxlo, inres)
    return
  end subroutine locate

  subroutine interpvel(rval, sval, m, mp1, n, np1, u, v)
    use amrex_constants_module, only: one
    implicit none

    double precision, intent(in) :: rval, sval
    double precision, intent(out) :: u, v
    integer, intent(in) :: n, np1, m, mp1

    double precision :: f0, f1, f2, f3
    double precision :: u0, u1, u2, u3
    double precision :: v0, v1, v2, v3

    f0 = (one - rval)*(one - sval)
    f1 = rval*(one - sval)
    f2 = (one - rval)*sval
    f3 = rval*sval
    u0 = uinput(m,n)
    u1 = uinput(mp1,n)
    u2 = uinput(m,np1)
    u3 = uinput(mp1,np1)
    u = u0*f0 + u1*f1 + u2*f2 + u3*f3
    v0 = vinput(m,n)
    v1 = vinput(mp1,n)
    v2 = vinput(m,np1)
    v3 = vinput(mp1,np1)
    v = v0*f0 + v1*f1 + v2*f2 + v3*f3
    return
  end subroutine interpvel
end module hit_turb_data_module

subroutine initialize_hit_turb(filename)
  use hit_turb_data_module
  character(len=*) :: filename
  integer :: n
  turbulence_filename = filename
  call read_turb()

end subroutine initialize_hit_turb
