program QSHO
    implicit none
    integer, parameter :: N = 30000
    real(8), parameter :: x_min = -5.0, x_max = 5.0
    real(8), parameter :: h = (x_max - x_min) / N
    real(8) :: x(N+1), psi(N+1), K(N+1), E
    real(8), parameter :: h2 = h**2
    real(8), parameter :: b = h2 / 12.0
    
    integer :: i
    integer :: unit_wavefunction
    real(8) :: norm_factor 
    integer :: last_digit

    E = 2.5
    
    unit_wavefunction = 10
    open(unit=unit_wavefunction, file='sho1.txt', status='replace')
  
    do i = 1, N+1
      x(i) = x_min + (i - 1) * h
      K(i) = 2.0 * E - x(i)**2 
    end do
    
    psi(1) = 0.0d0
    psi(2) = 0.001
  
    do i = 2, N
     ! psi(i+1) = (2.0 * (1.0 - 5.0 * b * K(i)) * psi(i) - (1.0 + b * K(i-1)) * psi(i-1)) / (1.0 + b * K(i+1))
        psi(i+1) = (2 * (1 - 5 * b * K(i)) * psi(i) - (1 + b * K(i-1)) * psi(i-1)) / (1 + b * K(i+1))
    end do
  
    norm_factor = 0.0
  
    do i = 1, N+1
        norm_factor = norm_factor + abs(psi(i)**2)
    end do
    norm_factor = sqrt(norm_factor * h)  
      
    psi = psi / norm_factor
    
    last_digit = mod(int(E), 10)

    if (mod(last_digit, 2) /= 0) then
        psi = -psi
    end if    
    do i = 1, N+1
        write(unit_wavefunction, '(F10.5, 2x, F10.5)') x(i), psi(i)
    end do
    close(unit_wavefunction)
  
    print *, 'Wavefunction data has been saved to sho1.txt'
  
  
  
  end program QSHO