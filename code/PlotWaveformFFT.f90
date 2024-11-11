program PlotWaveformFFT
    use iso_c_binding
    implicit none

    type, bind(c) :: Header
        character(kind=c_char), dimension(4) :: Code
        real(c_double)                       :: orgintime
        integer(c_int16_t)                   :: ncom
        integer(c_int32_t)                   :: ndata
        real(c_float)                        :: dt
    end type Header

    ! Declare variables
    type(Header) :: wh
    real(c_float), allocatable :: wd(:,:)
    integer :: i, j, k
    real(c_float), allocatable :: T(:)
    character(len=4) :: CodeString
    
    ! Variables for FFT
    integer :: n_fft
    complex, allocatable :: fft_data(:,:)
    real, allocatable :: freq(:)
    real, allocatable :: amplitude(:,:)
    
    ! Variables for PGPLOT
    integer :: pgid
    real(c_float), allocatable :: E1(:), E2(:), E3(:)
    real :: min_T, max_T, min_E, max_E
    real :: min_F, max_F, min_A, max_A

    ! Open and read binary file
    open(8, file='../data/seismicdata.bin', status='old', &
        form='unformatted', access='stream', convert='little_endian')
    read(8) wh
    
    ! Convert Code array to string
    CodeString = ''
    do i = 1, 4
        CodeString(i:i) = wh%Code(i)
    end do
    
    ! Print header info
    print *, 'Code:', trim(CodeString)
    print *, 'Origin Time:', wh%orgintime
    print *, 'Number of Components:', wh%ncom
    print *, 'Number of Data Points:', wh%ndata
    print *, 'Sampling Interval:', wh%dt
    
    ! Allocate and read waveform data
    allocate(wd(wh%ncom, wh%ndata))
    read(8) wd
    close(8)
    
    ! Create time array
    allocate(T(wh%ndata))
    do i = 1, wh%ndata
        T(i) = (i - 1) * wh%dt
    end do
    
    ! Prepare for FFT
    n_fft = wh%ndata
    allocate(fft_data(wh%ncom, n_fft))
    allocate(freq(n_fft/2))
    allocate(amplitude(wh%ncom, n_fft/2))
    
    ! Calculate frequency array
    do i = 1, n_fft/2
        freq(i) = (i-1) * (1.0/(wh%dt * n_fft))
    end do
    
    ! Perform FFT for each component
    do i = 1, wh%ncom
        ! Copy data to complex array
        do j = 1, n_fft
            fft_data(i,j) = cmplx(wd(i,j), 0.0)
        end do
        
        ! Perform FFT
        call fft(fft_data(i,:), n_fft, 1)
        
        ! Calculate amplitude spectrum
        do j = 1, n_fft/2
            amplitude(i,j) = sqrt(real(fft_data(i,j))**2 + aimag(fft_data(i,j))**2)
            ! Convert to dB
            if (amplitude(i,j) > 0.0) then
                amplitude(i,j) = 20.0 * log10(amplitude(i,j))
            else
                amplitude(i,j) = -200.0  ! Set minimum dB value
            end if
        end do
    end do
    
    ! Initialize PGPLOT
    call pgopen('../plot/waveform_fft.ps/VCPS')
    call pgsubp(2, 3)  ! 2 rows, 3 columns
    
    ! Plot time domain waveforms
    do i = 1, wh%ncom
        call pgsci(1)
        min_E = minval(wd(i,:))
        max_E = maxval(wd(i,:))
        call pgenv(minval(T), maxval(T), min_E, max_E, 0, 0)
        call pglab('Time (sec)', 'Amplitude', 'Waveform E'//char(i+48))
        call pgsci(i+1)
        call pgline(n_fft, T, wd(i,:))
    end do
    
    ! Plot frequency domain spectra
    do i = 1, wh%ncom
        call pgsci(1)
        min_A = minval(amplitude(i,:))
        max_A = maxval(amplitude(i,:))
        call pgenv(0.0, maxval(freq), min_A, max_A, 0, 0)
        call pglab('Frequency (Hz)', 'Amplitude (dB)', 'Spectrum E'//char(i+48))
        call pgsci(i+1)
        call pgline(n_fft/2, freq, amplitude(i,:))
    end do
    
    call pgclos()
    
    ! Cleanup
    deallocate(wd, T, fft_data, freq, amplitude)

contains
    ! FFT subroutine (Cooley-Tukey algorithm)
    subroutine fft(x, n, direction)
        complex, intent(inout) :: x(:)
        integer, intent(in) :: n, direction
        complex :: temp
        integer :: i, j, k, m, lstep
        real :: angle
        real, parameter :: pi = 3.14159265359
        
        ! Bit-reverse ordering
        j = 1
        do i = 1, n-1
            if (i < j) then
                temp = x(j)
                x(j) = x(i)
                x(i) = temp
            end if
            m = n/2
            do while (j > m .and. m > 0)
                j = j - m
                m = m/2
            end do
            j = j + m
        end do
        
        ! Compute FFT
        do lstep = 1, log2(real(n))
            m = 2**lstep
            do k = 1, n, m
                do i = k, k+m/2-1
                    angle = -2.0 * pi * direction * (i-k) / m
                    temp = x(i+m/2) * cmplx(cos(angle), sin(angle))
                    x(i+m/2) = x(i) - temp
                    x(i) = x(i) + temp
                end do
            end do
        end do
        
        ! Scale if inverse FFT
        if (direction == -1) then
            x = x / n
        end if
    end subroutine fft
    
    ! Function to compute log2
    function log2(x)
        real, intent(in) :: x
        integer :: log2
        log2 = int(log(x)/log(2.0))
    end function log2
    
end program PlotWaveformFFT
