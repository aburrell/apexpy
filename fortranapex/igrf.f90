! Uncomment below if using an Intel compiler
! use ifport

module igrf
      
  implicit none

contains

  subroutine read_igrf(filename_in,GYR,HYR,GT,HT,NEPO,NGHT,EPOCH,NMXE)

    implicit none

    real*8,allocatable,intent(inout) :: GYR(:,:,:),HYR(:,:,:)
    real*8,allocatable,intent(inout) :: GT(:,:),HT(:,:)
    real*4,allocatable,intent(inout) :: EPOCH(:),NMXE(:)
    integer, intent(out) :: NEPO,NGHT
    character(len=*),intent(in) :: filename_in

    character(len=10000) :: s
    character(len=100)     :: junk
    integer   :: state, i, offset, pos,o
    integer   :: num_sh, L_max
    integer   :: num_epochs
    integer   :: l,m,e
    real*8,allocatable    :: g(:,:)
    integer,allocatable   :: nm(:,:)

    ! Get number of Gauss coefficients
    !num_sh = NGHT-2*sqrt(real(NGHT))
    !write(*,*) num_sh

    ! Open IGRF file
    open(unit=100, file=filename_in,status='old',iostat=state)
    if (state /= 0) then
       stop "File open error"
    end if

    ! Skip comment lines
    do
       read(unit=100,fmt='(A)') s
       if (s(1:1) .ne. '#') exit
    enddo

    ! Read epochs
    num_epochs=count([(s(i:i+3),i=1,len_trim(s))].eq.'IGRF')
    num_epochs=count([(s(i:i+3),i=1,len_trim(s))].eq.'DGRF')+num_epochs
    allocate(EPOCH(1:num_epochs))
    allocate(NMXE(1:num_epochs))

    ! Read epochs
    read(100,*,iostat=state) s
    do i=1,num_epochs
       EPOCH(i) = 1900+(i-1)*5.0d0
    enddo

    ! Number of coefficients
    do i=1,num_epochs
       if (EPOCH(i) .ge. 2000.0d0) then
          NMXE(i) = 13
       elseif (EPOCH(i) .ge. 1900.0d0) then
          NMXE(i) = 10
       else
          write(*,*) 'ERROR: Epoch unavailable!'
          exit
       endif
    enddo

    ! Save file position
    offset = ftell(100)

    ! Get the number of lines (coefficients)
    num_sh = 0
    do
       read(unit=100,fmt=*,iostat=state) 
       if (state < 0) exit
       num_sh = num_sh+1 
    enddo
    L_max = sqrt(num_sh+1.0d0)-1
    close(100)

    ! Restore file position, must re-open after reaching EOF
    open(unit=100, file=filename_in,status='old',iostat=state)
    call fseek(100,offset,0,state)

    ! Assign the variables for Gauss coefficients
    allocate(g(1:num_sh,1:num_epochs))
    allocate(nm(1:num_sh,2))

    ! Read coefficients
    do i=1,num_sh
       read(100,*,iostat=state) s,nm(i,1),nm(i,2),g(i,:)
       if (state < 0) exit
    enddo
    close(100)

    ! Assign the return values
    NGHT = (sqrt(num_sh+1.0)+1)**2
    NEPO = num_epochs

    ! Assign the Gauss coefficients
    allocate(GYR(NGHT,NGHT,NEPO),GT(NGHT,NGHT))
    allocate(HYR(NGHT,NGHT,NEPO),HT(NGHT,NGHT))
    GYR=0.0d0
    GT =0.0d0
    HYR=0.0d0
    HT =0.0d0
    do e=1,NEPO+1
       do l=1,L_max
          if (e .le. NEPO) then
             GYR(l+1,1,e) = g(l**2,e)
          else
             GT(l+1,1)    = g(l**2,e)
          endif
          do m=1,l
             if (m .le. l) then
                pos = 1 + m*(L_max+2) + l
                if (e .le. NEPO) then
                   GYR(l+1,m+1,e)=g(l**2+2*m-1,e)
                   HYR(l+1,m+1,e)=g(l**2+2*m  ,e)
                else
                   GT(l+1,m+1)=g(l**2+2*m-1,e)
                   HT(l+1,m+1)=g(l**2+2*m  ,e)
                endif
             endif
          enddo
       enddo
    enddo
    return
  end subroutine read_igrf

end module igrf
