program PlotWaveform
    use iso_c_binding
    implicit none

    type, bind(c) :: Header
        character(kind=c_char), dimension(4) :: Code        ! 4-character array
        real(c_double)                       :: orgintime   ! Origin time (8 bytes)
        integer(c_int16_t)                   :: ncom        ! Number of components (2 bytes)
        integer(c_int32_t)                   :: ndata       ! Number of data points (4 bytes)
        real(c_float)                        :: dt          ! Sampling interval (4 bytes)
    end type Header

    ! Declare variables
    type(Header) :: wh                      ! Waveform header
    real(c_float), allocatable :: wd(:,:)   ! Waveform data array
    integer :: i, j
    real(c_float), allocatable :: T(:)      ! Time array
    character(len=4) :: CodeString

    ! Variables for PGPLOT
    integer :: n
    integer :: pgid
    real(c_float), allocatable :: E1(:), E2(:), E3(:)
    real :: min_T, max_T, min_E, max_E

    ! Open the binary file for unformatted stream access with appropriate endianness
    open(8, file='../data/seismicdata.bin', status='old', &
        form='unformatted', access='stream', convert='little_endian')

    ! Read the header
    read(8) wh

    ! Convert the Code array to a string
    CodeString = ''
    do i = 1, 4
        CodeString(i:i) = wh%Code(i)
    end do

    ! Print the header values for debugging
    print *, 'Code:', trim(CodeString)
    print *, 'Origin Time:', wh%orgintime
    print *, 'Number of Components:', wh%ncom
    print *, 'Number of Data Points:', wh%ndata
    print *, 'Sampling Interval:', wh%dt

    ! Allocate the waveform data array based on header information
    allocate(wd(wh%ncom, wh%ndata))

    ! Read the waveform data
    read(8) wd

    ! Close the file
    close(8)

    ! Create time array
    allocate(T(wh%ndata))
    do i = 1, wh%ndata
        T(i) = (i - 1) * wh%dt
    end do

    ! Assign n
    n = wh%ndata

    ! Allocate arrays for components
    allocate(E1(n), E2(n), E3(n))
    E1 = wd(1, :)
    E2 = wd(2, :)
    E3 = wd(3, :)

    ! Initialize PGPLOT
    call pgopen('../plot/waveform.ps/VCPS')

    call pgsubp(1, 3)            ! Divide the page into 1 row and 3 columns

    call pgsci(1)                ! Set color index to 1 (default foreground)
    call pgslw(2)                ! Set line width
    call pgsch(1.0)              ! Set character height
    call pgscf(2)                ! Set character font

    ! Plot E1
    call pgsci(1)
    min_E = minval(E1)
    max_E = maxval(E1)
    call pgenv(minval(T), maxval(T), min_E, max_E, 0, 0)
    call pglab('Time (sec)', 'Amplitude', 'Waveform E1')
    call pgsci(2)
    call pgline(n, T, E1)

    ! Plot E2
    call pgsci(1)
    min_E = minval(E2)
    max_E = maxval(E2)
    call pgenv(minval(T), maxval(T), min_E, max_E, 0, 0)
    call pglab('Time (sec)', 'Amplitude', 'Waveform E2')
    call pgsci(3)
    call pgline(n, T, E2)

    ! Plot E3
    call pgsci(1)
    min_E = minval(E3)
    max_E = maxval(E3)
    call pgenv(minval(T), maxval(T), min_E, max_E, 0, 0)
    call pglab('Time (sec)', 'Amplitude', 'Waveform E3')
    call pgsci(4)
    call pgline(n, T, E3)

    call pgclos()

    deallocate(E1, E2, E3)

end program PlotWaveform
