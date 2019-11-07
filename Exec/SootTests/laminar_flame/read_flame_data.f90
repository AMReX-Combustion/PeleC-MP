module flame_data_module
  implicit none
  character(len=72) :: flame_filename = ''
  integer :: flame_init = 0, flame_M, flame_N
  double precision, dimension(:), allocatable :: flame_X, flame_T, flame_U, &
       flame_Rho
  double precision, allocatable ::  flame_Y(:,:)
  character (len=20), allocatable :: flame_names(:)

contains

  subroutine read_flame()

    use amrex_constants_module, only: zero, one
    character (len=20) :: ctmp1, ctmp2, fmt
    character (len=1) :: ctmp
    integer index, found, ist, reason, i, j, lsize, pos1, pos2, n

    integer, parameter :: maxlinelen = 1024
    character(len=maxlinelen) :: line
    character(len=maxlinelen) :: tline, name

    integer :: idbg
    double precision :: ynsum

    !     Read 2 header lines, first looks like VARIABLES = NAME1 NAME2 NAME3..., we dont care about second

    open(unit=32,file=flame_filename,status='old')

    reason = 0

    !  Get length of first line
    lsize = 1
    do while (.true.)
       read(32,'(A)',EOR=5,END=5,ADVANCE='NO') ctmp
       lsize = lsize+1
    enddo
5   if (lsize .gt. maxlinelen) then
       print *,'maxlinelen in read_flame_data not large enough for first line of data file'
       stop
    endif
    rewind(32)

    write(fmt,'("(A", I0, ")")') lsize
    read(32,trim(fmt)) line

    ! Get number of variables
    flame_M = 0

    tline = line
    DO
       do while (tline(1:1) .eq. ' ' .and. len(tline).gt.0)
          tline = tline(2:)
       enddo
       pos1 = INDEX(tline, '"')

       IF (len(tline).eq.0 .or. pos1 == 0) THEN
          EXIT
       END IF
       flame_M = flame_M + 1

       tline = tline(pos1+1:)

       pos1 = INDEX(tline, '"')

       if (pos1.eq.0) then
          print *,'variable name in data file missing a closing quote'
       endif

       tline = tline(pos1+1:)

       if (len(trim(tline)).eq.0) EXIT
    END DO

    ! Allocate space for names, and do all this all over again
    ! (fortran text parsing is torture!)
    allocate(flame_names(flame_M))
    flame_M = flame_M - 4 ! remove the X, U, T, and rho

    n = 0
    tline = line
    DO
       do while (tline(1:1) .eq. ' ' .and. len(tline).gt.0)
          tline = tline(2:)
       enddo
       pos1 = INDEX(tline, '"')

       IF (len(tline).eq.0 .or. pos1 == 0) THEN
          EXIT
       END IF

       tline = tline(pos1+1:)

       pos1 = INDEX(tline, '"')

       n = n + 1
       flame_names(n) = tline(:pos1-1)

       tline = tline(pos1+1:)

       if (len(trim(tline)).eq.0) EXIT
    END DO

    !  Count the lines of data
    flame_N = 0
    do while (.true.)
       read(32,*,END=2) line
       flame_N = flame_N+1
    enddo
2   continue

    allocate(flame_X(flame_N))
    allocate(flame_U(flame_N))
    allocate(flame_T(flame_N))
    allocate(flame_Rho(flame_N))
    allocate(flame_Y(flame_N,flame_M))
    
    !  Now read data
    rewind(32)
    read(32,'(A)',advance='yes') line ! Skip line with variable names
    do i = 1,flame_N
       read(32,*) flame_X(i),flame_U(i),flame_T(i),flame_Rho(i),&
            (flame_Y(i,j),j=1,flame_M)
    enddo
    !  Re-normalize mass fraction data to ensure they sum to 1
    do i = 1,flame_N
       ynsum = zero
       do j = 1,flame_M
          flame_Y(i,j) = max(zero, min(one, flame_Y(i,j)))
          ynsum = ynsum + flame_Y(i,j)
       enddo
       ynsum = one/ynsum
       do j = 1, flame_M
          flame_Y(i,j) = flame_Y(i,j)*ynsum
       enddo
    enddo
    close(32)
    !  Now mark that we have read the data
    flame_init = 1
  end subroutine read_flame

  function flame_ncomp() result(ncomp)
    integer :: ncomp
    ncomp = flame_M
  end function flame_ncomp

  function flame_npts() result(npts)
    integer :: npts
    npts = flame_N
  end function flame_npts

  ! Interpolate values from current grid
  subroutine interp_flame(xval,rho,T,U,y_vector)
    double precision, intent(in) :: xval
    double precision, intent(out) :: rho, T, U, y_vector(*)
    integer i,j, loside,hiside
    double precision :: x1,x2,ddx,val1,val2

    if (flame_init .eq. 0) then
       if (trim(flame_filename) .eq. '') then
          call bl_abort('must set flame_filename prior to calling interp_flame')
       endif
       call read_flame()
    endif
    loside = 0
    hiside = 0
    if (xval .le. flame_X(1)) then
       loside = 1
       hiside = 1
    end if
    if (xval .ge. flame_X(flame_N)) then
       loside = flame_N
       hiside = flame_N
    end if
    if (loside.eq.0) then
       do i = 1, flame_N-1                           
          if ( (xval .ge. flame_X(i)) .and.  &
               (xval .le. flame_X(i+1)) ) then
             loside  = i
             hiside  = i+1
          end if
       end do
    end if
    x1 = flame_X(loside)
    x2 = flame_X(hiside)
    if (loside.eq.hiside) then
       ddx = 0.d0
    else
       ddx = (xval-x1)/(x2-x1)
    end if
    ! Interpolate rho
    val1 = flame_Rho(loside)
    val2 = flame_Rho(hiside)
    rho = val1 + (val2-val1)*ddx
    ! Interpolate temperature
    val1 = flame_T(loside)
    val2 = flame_T(hiside)
    T = val1 + (val2-val1)*ddx
    ! Interpolate velocity
    val1 = flame_U(loside)
    val2 = flame_U(hiside)
    U = val1 + (val2-val1)*ddx
    ! Interpolate mass fractions
    do j = 1, flame_M
       val1 = flame_Y(loside,j)
       val2 = flame_Y(hiside,j)
       y_vector(j) = val1 + (val2-val1)*ddx
    end do
  end subroutine interp_flame

end module flame_data_module

subroutine initialize_flame_data(filename, Lx)
  use flame_data_module
  use chemistry_module, only : nspecies, get_species_index
  use amrex_constants_module, only : half, zero
  double precision, intent(in) :: Lx ! Length of computational domain
  character (len=*) :: filename
  integer :: n
  double precision :: Linput
  flame_filename = filename
  call read_flame()
  if (flame_ncomp() .ne. nspecies) then
     print *,'Number of dependent variables in file:',flame_ncomp()
     print *,'Number expected:',nspecies
     stop 'flame data file not compatible with current chemistry model, wrong number of species'
  endif

  if (flame_names(1) .ne. "X")  stop 'flame data file not compatible with flame data reader, X must be first variable'
  if (flame_names(2) .ne. "U" ) stop 'flame data file not compatible with flame data reader, U must be second variable'
  if (flame_names(3) .ne. "T")  stop 'flame data file not compatible with flame data reader, T must be third variable'
  if (flame_names(4) .ne. "rho")  stop 'flame data file not compatible with flame data reader, rho must be fourth variable'
  do n=1,nspecies
     if (n .ne. get_species_index(flame_names(4+n))) then
        stop 'flame data file not compatible with current chemistry model, wrong species'
     endif
  enddo
  Linput = flame_X(flame_N)
  if (abs(Linput - Lx) > 1.d-8) then
     call bl_abort('Lengths of domain and initial flame data domain are not the same')
  endif
end subroutine initialize_flame_data

subroutine flame(xval,rho,T,U,y_vector)
  use flame_data_module
  double precision :: xval, rho, T, U, y_vector(*)
  call interp_flame(xval,rho,T,U,y_vector)
end subroutine flame

subroutine flame_inflow(rho,T,U,y_vector)
  use flame_data_module
  use network, only: nspecies
  double precision :: rho, T, U, y_vector(*)
  rho = flame_Rho(1)
  T = flame_T(1)
  U = flame_U(1)
  y_vector(1:nspecies) = flame_Y(1,1:nspecies)
end subroutine flame_inflow

subroutine flame_outflow(rho,T,U,y_vector)
  use flame_data_module
  use network, only: nspecies
  double precision :: rho, T, U, y_vector(*)
  rho = flame_Rho(flame_N)
  T = flame_T(flame_N)
  U = flame_U(flame_N)
  y_vector(1:nspecies) = flame_Y(flame_N,1:nspecies)
end subroutine flame_outflow
