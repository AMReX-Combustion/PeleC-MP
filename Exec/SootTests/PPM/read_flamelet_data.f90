module flamelet_data_module
  implicit none
  character(len=72) :: flamelet_filename = ''
  integer :: flamelet_init = 0, flamelet_M, flamelet_N
  double precision, allocatable :: flamelet_Z(:)
  double precision, allocatable :: flamelet_T(:)
  double precision, allocatable ::  flamelet_Y(:,:)
  character (len=20), allocatable :: flamelet_names(:)

contains

  subroutine read_flamelet()

    use amrex_constants_module, only: zero, one
    character (len=20) :: ctmp1, ctmp2, fmt
    character (len=1) :: ctmp
    integer index, found, ist, reason, i, j, lsize, pos1, pos2, n

    integer, parameter :: maxlinelen = 512
    character(len=maxlinelen) :: line
    character(len=maxlinelen) :: tline, name

    integer :: idbg
    double precision :: ynsum

    !     Read 2 header lines, first looks like VARIABLES = NAME1 NAME2 NAME3..., we dont care about second
    open(unit=32,file=flamelet_filename,status='old')

    reason = 0

    !  Get length of first line
    lsize = 1
    do while (.true.)
       read(32,'(A)',EOR=5,END=5,ADVANCE='NO') ctmp
       lsize = lsize+1
    enddo
5   if (lsize .gt. maxlinelen) then
       print *,'maxlinelen in read_flamelet_data not large enough for first line of data file'
       stop
    endif
    rewind(32)

    write(fmt,'("(A", I0, ")")') lsize
    read(32,trim(fmt)) line

    ! Get number of variables
    ist = INDEX(line, "=")
    line = trim(line(ist+1:))
    flamelet_M = 0

    tline = line
    DO
       do while (tline(1:1) .eq. ' ' .and. len(tline).gt.0)
          tline = tline(2:)
       enddo
       pos1 = INDEX(tline, '"')

       IF (len(tline).eq.0 .or. pos1 == 0) THEN
          EXIT
       END IF
       flamelet_M = flamelet_M + 1

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
    allocate(flamelet_names(flamelet_M))
    flamelet_M = flamelet_M - 2 ! remove the Z and T

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
       flamelet_names(n) = tline(:pos1-1)

       tline = tline(pos1+1:)

       if (len(trim(tline)).eq.0) EXIT
    END DO

    !  Count the lines of data
    flamelet_N = 0
    do while (.true.)
       read(32,*,END=2) line
       flamelet_N = flamelet_N+1
    enddo
2   continue

    allocate(flamelet_Z(flamelet_N))
    allocate(flamelet_T(flamelet_N))
    allocate(flamelet_Y(flamelet_N,flamelet_M))
    
    !  Now read data
    rewind(32)
    read(32,'(A)',advance='yes') line ! Skip line with variable names
    do i = 1,flamelet_N
       read(32,*) flamelet_Z(i),flamelet_T(i),(flamelet_Y(i,j),j=1,flamelet_M)
    enddo
    !  Re-normalize mass fraction data to ensure they sum to 1
    do i = 1,flamelet_N
       ynsum = zero
       do j = 1,flamelet_M
          flamelet_Y(i,j) = max(zero, min(one, flamelet_Y(i,j)))
          ynsum = ynsum + flamelet_Y(i,j)
       enddo
       ynsum = one/ynsum
       do j = 1, flamelet_M
          flamelet_Y(i,j) = flamelet_Y(i,j)*ynsum
       enddo
    enddo
    close(32)
    !  Now mark that we have read the data
    flamelet_init = 1
  end subroutine read_flamelet

  function flamelet_ncomp() result(ncomp)
    integer :: ncomp
    ncomp = flamelet_M
  end function flamelet_ncomp

  function flamelet_npts() result(npts)
    integer :: npts
    npts = flamelet_N
  end function flamelet_npts

  subroutine interp_flamelet(zval,y_vector,Tval)
    double precision, intent(in) :: zval
    double precision, intent(out) :: y_vector(*), Tval
    integer i,j, loside,hiside
    double precision ylo,yhi,z1,y1,z2,y2,dydx,T1,T2

    if (flamelet_init .eq. 0) then
       if (trim(flamelet_filename) .eq. '') then
          call bl_abort('must set flamelet_filename prior to calling interp_flamelet')
       endif
       call read_flamelet()
    endif
    loside = 0
    hiside = 0
    if (zval .le. flamelet_Z(1)) then
       loside = 1
       hiside = 1
    end if
    if (zval .ge. flamelet_Z(flamelet_N)) then
       loside = flamelet_N
       hiside = flamelet_N
    end if
    if (loside.eq.0) then
       do i = 1, flamelet_N-1                           
          if ( (zval .ge. flamelet_Z(i)) .and.  &
               (zval .le. flamelet_Z(i+1)) ) then
             loside  = i
             hiside  = i+1
          end if
       end do
    end if
    z1 = flamelet_Z(loside)
    z2 = flamelet_Z(hiside)
    T1 = flamelet_T(loside)
    T2 = flamelet_T(hiside)
    if (loside.eq.hiside) then
       dydx = 0.d0
    else
       dydx = (T2-T1)/(z2-z1)
    end if
    Tval = T1 + dydx*(zval-z1)
    do j = 1, flamelet_M
       y1 = flamelet_Y(loside,j)
       y2 = flamelet_Y(hiside,j)

       if (loside.eq.hiside) then
          dydx = 0.d0
       else
          dydx = (y2-y1)/(z2-z1)
       end if

       y_vector(j) = y1 + dydx*(zval - z1)
    end do
  end subroutine interp_flamelet

end module flamelet_data_module

subroutine initialize_flamelet(filename)
  use flamelet_data_module
  use chemistry_module, only : nspecies, get_species_index
  character (len=*) :: filename
  integer :: n
  flamelet_filename = filename
  call read_flamelet()
  if (flamelet_ncomp() .ne. nspecies) then
     print *,'Number of dependent variables in file:',flamelet_ncomp()
     print *,'Number expected:',nspecies
     stop 'flamelet data file not compatible with current chemistry model, wrong number of species'
  endif

  if (flamelet_names(1) .ne. "Z")  stop 'flamelet data file not compatible with flamelet data reader, Z must be first variable'
  if (flamelet_names(2) .ne. "T")  stop 'flamelet data file not compatible with flamelet data reader, T must be second variable'
  do n=1,nspecies
     if (n .ne. get_species_index(flamelet_names(2+n))) then
        stop 'flamelet data file not compatible with current chemistry model, wrong species'
     endif
  enddo
end subroutine initialize_flamelet

subroutine flamelet(zval,y_vector,Tval)
  use flamelet_data_module
  double precision :: zval,y_vector(*),Tval
  call interp_flamelet(zval,y_vector,Tval)
end subroutine flamelet

