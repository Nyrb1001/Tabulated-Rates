module fermi_rates_module
  implicit none
contains

  ! Subroutine to get data from a file
  subroutine get_data(n, z, Tmin, Tmax, mumin, mumax, rate)
    integer, intent(in) :: n, z, Tmin, Tmax, mumin, mumax
    real, dimension(4,4), intent(out):: rate
    character(len=100) :: root, path_rate
    integer :: iunit
    character(len=100) :: filename
    real, dimension(29,89) :: temp_rate

    ! Generate the file path for the rate data
    write(root, '("Rates/z", I3.3)') z
    write(filename, '("z", I3.3, "n", I3.3, ".dat")') z, n
    path_rate = trim(root) // '/' // trim(filename)
    iunit = 10
    ! Open and read the rate data
    open(unit=iunit, file=path_rate, status="old")
    read(iunit, *) temp_rate
    close(iunit)

    ! Extract the relevant subset of the data
    rate = temp_rate(Tmin:Tmax, mumin:mumax)
    !!In the case of small inputs, repeat the lowest value to ensure the second derivative is defined.
    select case(mumin)
        case (1:86)
        case (0)
            rate(1:4, 1) = rate(1:4, 2)
        case(87)
            rate(1:4, 4) = rate(1:4, 3)
    end select
    
    select case(Tmin)
        case (1:26)
        case (0)
            rate(1, 1:4) = rate(2, 1:4)
        case(27)
            rate(4, 1:4) = rate(3,1:4)
    end select
end subroutine get_data

  subroutine get_indices(T_in, mu_in, Tmin, Tmax, mumin, mumax)
  ! Routine to return T and chemical potnentialindicies that enclose a 4x4 subset of the entire array of rates
  ! The central 2x2 is used for interpolation, while the bordering slices 1: and 4: are used to calculate the second derivatives
    implicit none
    real, intent(in) :: T_in, mu_in
    integer :: i
    !real, parameter :: mu(101) = [(i, i = 0, 100)] * 0.04
    real(8), parameter :: mu(89) = (/ &
     (0.0d0 + 0.09d0 * (i - 1) / 9.0d0, i = 1, 10), &
     (0.10d0 + (4.0d0 - 0.10d0) * (i - 1) / 78.0d0, i = 1, 79) &
  /)
    real, parameter :: T_Kelvin(29) = &
    (/ 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, &
    0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, &
    7.0, 8.0, 9.0, 10.0 /)
    integer, intent(out) :: Tmin, Tmax, mumin, mumax
    write(*,*) mu
    if (T_in > 10.0 .or. mu_in > 4.0) then
        write(*,*) 'Out of bounds input'
        stop
    end if
    do i = 2, size(T_Kelvin)
      if (T_in <= T_Kelvin(i)) then
        Tmin = i - 2
        Tmax = i + 1
        exit
      end if
    end do

    do i = 2, size(mu)
      if (mu_in <= mu(i)) then
        mumin = i - 2
        mumax = i + 1
        exit
      end if
    end do
  end subroutine get_indices

  ! Function for bilinear interpolation
  function bilinear_interpolation(T_in, mu_in, rate, T_sub, mu_sub) result(f)
    real, intent(in) :: mu_in, T_in
    real, dimension(:,:), intent(in) :: rate
    real, dimension(4), intent(in) :: mu_sub, T_sub
    real :: f
    real :: f_1, f_2
    real :: x1, x2, y1, y2

    x1 = mu_sub(2)
    x2 = mu_sub(3)
    y1 = log10(T_sub(2))
    y2 = log10(T_sub(3))
    ! Perform the interpolation along chem axis
    f_1 = rate(2,2)*(x2 - mu_in) / (x2 - x1) + rate(2,3)*(mu_in - x1)/(x2 - x1)
    f_2 = rate(3,2)*(x2 - mu_in) / (x2 - x1) + rate(3,3)*(mu_in - x1)/(x2 - x1)
    ! (Log) Inteprolate along T
    f = 10**(log10(f_1)*(y2 - log10(T_in)) / (y2 - y1) + log10(f_2)*(log10(T_in) - y1)/(y2 - y1))
    
  end function bilinear_interpolation

  ! Subroutine for calculating the second derivative for error estimation
  real function second_diff(arr, x, ind)
    real, dimension(:), intent(in) :: arr, x
    integer, intent(in) :: ind
    real :: hs, hd
    hs = x(ind) - x(ind - 1)
    hd = x(ind + 1) - x(ind)
    second_diff = 2.0*(hs*arr(ind + 1) - (hd + hs)*arr(ind) + hd*arr(ind - 1)) / (hs * hd * (hd + hs))
  end function second_diff


  ! Subroutine to estimate the error
  real function estimate_error(rate, T, mu)
    real, dimension(4,4), intent(in) :: rate
    real, dimension(:), intent(in) :: mu, T
    real :: ftt_1, ftt_2, fmm_1, fmm_2, err_T, err_mu

    ftt_1 = second_diff(rate(:,2), T, 2)
    ftt_2 = second_diff(rate(:,2), T, 3)
    
    fmm_1 = second_diff(rate(2,:), mu, 2)
    fmm_2 = second_diff(rate(2,:), mu, 3)

    err_T = (1.0 / 8.0) * max(abs(ftt_1), abs(ftt_2)) * (T(3) - T(2))**2
    err_mu = (1.0 / 8.0) * max(abs(fmm_1), abs(fmm_2)) * (mu(3) - mu(2))**2
    estimate_error = sqrt(err_T**2 + err_mu**2)

  end function estimate_error

  ! Main subroutine to get rate and error
  subroutine get_rate(n, z, T_in, mu_in, fermi_rate, error)
    integer, intent(in) :: n, z
    real, intent(in) :: T_in, mu_in
    real, intent(out) :: fermi_rate, error
    integer :: Tmin, Tmax, mumin, mumax, i
    
    real(8), parameter :: mu(89) = (/ &
     (0.0d0 + (0.09d0 - 0.0d0) * (i - 1) / 9.0d0, i = 1, 10), &
     (0.10d0 + (4.0d0 - 0.10d0) * (i - 1) / 78.0d0, i = 1, 79) &
  /)
    real, parameter :: T_Kelvin(29) = &
    (/ 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, &
    0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, &
    7.0, 8.0, 9.0, 10.0 /) 
    
    real, dimension(4) :: mu_sub, T_sub
    real, dimension(4,4) :: rate

    ! Get the indices for temperature and chemical potential
    call get_indices(T_in, mu_in, Tmin, Tmax, mumin, mumax)

    ! Get the subarrays for temperature and chemical potential
    
    select case(mumin)
        case (1:86)
            mu_sub = mu(mumin:mumax)
        case (0)
            mu_sub = (/-0.0001d0, mu(1), mu(2), mu(3)/)
        case(87)
            mu_sub = (/mu(87), mu(88), mu(89), 4.1d0/)
    end select
    
    select case(Tmin)
        case (1:26)
            T_sub = T_Kelvin(Tmin:Tmax)
        case (0)
            T_sub = (/-0.0001, T_Kelvin(1), T_Kelvin(2), T_Kelvin(3)/)
        case(27)
            T_sub = (/T_Kelvin(27), T_Kelvin(28), T_Kelvin(29), 15.0/)
    end select

    ! Get the rate data from file
    call get_data(n, z, Tmin, Tmax, mumin, mumax, rate)

    ! Perform bilinear interpolation and estimate error
    fermi_rate = bilinear_interpolation(T_in, mu_in, rate, T_sub, mu_sub)
    error = estimate_error(rate, T_sub, mu_sub)

  end subroutine get_rate

end module fermi_rates_module