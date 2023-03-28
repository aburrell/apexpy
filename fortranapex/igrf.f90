! Uncomment below if using an Intel compiler
! use ifport

module igrf
  ! Read in relevant data from IGRF coefficient file.
  !
  ! INPUTS:
  ! filename_in = Filename for IGRF coefficient files
  !
  ! RETURNS:
  ! None (all output saved as module variables)


  implicit none

  real(8)                  :: datel
  integer(4)               :: nepo, nght
  real(8), allocatable     :: epoch(:), nmxe(:)
  real(8), allocatable     :: gyr(:,:,:), hyr(:,:,:)
  real(8), allocatable     :: gt(:,:), ht(:,:)

contains

  subroutine read_igrf(filename_in)

    implicit none

    character(len=*), intent(in) :: filename_in
    character(len=10000)         :: s
    integer(8)                   :: offset
    integer(4)                   :: state, i, pos
    integer(4)                   :: num_sh, L_max
    integer(4)                   :: num_epochs
    integer(4)                   :: l, m, e
    real(8), allocatable         :: g(:,:)
    integer(4), allocatable      :: nm(:,:)

    ! Open IGRF file
    open(unit=100, file=filename_in, status='old', iostat=state)
    if (state /= 0) then
       stop "File open error"
    end if

    ! Skip comment lines
    do
       read(unit = 100, fmt = '(A)') s
       if (s(1:1) /= '#') exit
    end do

    ! Read epochs
    num_epochs = 0
    do i=1, len_trim(s)
      if ((s(i:i+3) == 'IGRF').or.(s(i:i+3) == 'DGRF')) then
        num_epochs = num_epochs + 1
      end if
    end do
    allocate(epoch(1:num_epochs))
    allocate(nmxe(1:num_epochs))

    ! Read epochs
    read(unit = 100, fmt = '(A)', iostat = state) s
    read(s(8:), *) epoch


    ! Number of coefficients
    do i = 1, num_epochs
       if (epoch(i) .ge. 2000.0d0) then
          nmxe(i) = 13
       elseif (epoch(i) .ge. 1900.0d0) then
          nmxe(i) = 10
       else
          write(*,*) 'ERROR: Epoch unavailable!'
          exit
       end if
    end do

    ! Save file position
    offset = ftell(100)

    ! Get the number of lines (coefficients)
    num_sh = 0
    do
       read(unit = 100, fmt =*, iostat = state)
       if (state < 0) exit
       num_sh = num_sh + 1
    end do
    L_max = int(sqrt(num_sh + 1.0d0)) - 1
    close(100)

    ! Restore file position, must re-open after reaching EOF
    open(unit=100, file=filename_in, status='old', iostat=state)
    call fseek(100, offset, 0, state)

    ! Assign the variables for Gauss coefficients
    allocate(g(1:num_sh, 1:num_epochs+1))
    allocate(nm(1:num_sh, 2))

    ! Read coefficients
    do i = 1, num_sh
       read(100,*, iostat = state) s, nm(i, 1), nm(i, 2), g(i,:)
       if (state < 0) exit
    end do
    close(100)

    ! Assign the return values
    nght = (int(sqrt(num_sh + 1.0)) + 1) ** 2
    nepo = num_epochs

    ! Assign the Gauss coefficients
    allocate(gyr(nght, nght, nepo), gt(nght, nght))
    allocate(hyr(nght, nght, nepo), ht(nght, nght))
    gyr = 0.0d0
    gt = 0.0d0
    hyr = 0.0d0
    ht = 0.0d0
    do e = 1, nepo + 1
       do l = 1, L_max
          if (e .le. nepo) then
             gyr(l + 1, 1, e) = g(l ** 2, e)
          else
             gt(l + 1, 1)    = g(l ** 2, e)
          end if
          do m = 1, l
             if (m .le. l) then
                pos = 1 + m * (L_max + 2) + l
                if (e .le. nepo) then
                   gyr(l + 1, m + 1, e) = g(l ** 2 + 2 * m - 1, e)
                   hyr(l + 1, m + 1, e) = g(l ** 2 + 2 * m  , e)
                else
                   gt(l + 1, m + 1) = g(l ** 2 + 2 * m - 1, e)
                   ht(l + 1, m + 1) = g(l ** 2 + 2 * m  , e)
                end if
             end if
          end do
       end do
    end do
    return
  end subroutine read_igrf

end module igrf
