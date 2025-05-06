program main
    use fermi_rates_module
    implicit none
    real :: val, error
    !! Neutron number, proton number
    integer :: n, z 
    !!T is in units of GK, mu is units of MeV
	real :: T, mu 
    
	
	T = 0.5
	mu = 2.2
    
	z = 73
	n = 149
    
    call get_rate(n, z, T, mu, val, error)
    write(*,*) val, error


    
end program