module asymptote
    use kind_type
    use global 
    implicit none
    complex(dp), allocatable, private, protected :: f(:)
contains 
    

! ==================================================
! MATRIX
! ==================================================
! matrix f -----------------------------------------
subroutine SUB_f
    use math_const, only: i => math_i
    use boundary, only: S 
    real   (dp) :: k, tmp1, tmp2 
    integer(i4) :: j 
    k = (2.d0*Scatt)**0.50
    do j = 0, L 
        tmp1 = aimag(S(j))/2.d0 
        tmp2 = (1.d0 -real(S(j)))/2.d0 
        f(j) = (2.d0*dble(j) +1.d0)/k*(tmp1 +i*tmp2)
    end do 
end subroutine SUB_f 
! end matrix ---------------------------------------









! ==================================================
! PROCESS
! ==================================================
! cs plot ------------------------------------------
subroutine PROC_CS_plot 
    use math_const, only: pi => math_pi, degree => math_degree
    use unit_const, only: au_bohr
    use hamiltonian, only: theta_grid
    use fgsl, only: fgsl_sf_legendre_Pl
    integer  (i1), parameter :: file_dcs = 101, file_tcs = 102 
    character(30), parameter :: form_cs  = '(30ES25.10)'
    character(30), parameter :: form_out = '(1A15, 5X, 1ES25.10)'
    real     (dp), parameter :: radian_to_degree = 1.d0/degree 
    real     (dp), parameter :: au_to_AA = au_bohr*10.d0**10.d0 
    complex  (dp) :: tmp1 
    real     (dp) :: unit_theta, unit_cs, k, tmp2 
    complex  (qp) :: sum
    integer  (i4) :: i, j 

    unit_theta = 1.d0 
    if(op_degree == 1) unit_theta = radian_to_degree
    unit_cs    = 1.d0 
    if(op_aa == 1) unit_cs = (au_to_AA)**2.d0

    if(allocated(f)) deallocate(f)
    allocate(f(0:L))
    call SUB_f

    open(file_tcs, file = "output/total_cs.d")
    sum = 0.d0 
    k   = (2.d0*Scatt)**0.50
    do i = 0, L 
        tmp1 = 4.d0*pi/k*f(i)
        sum  = sum +tmp1 
        write(file_tcs, form_cs) dble(i), aimag(tmp1)
    end do 
    tmp1 = sum 
    write(file_log, form_out) "total sigma: ", aimag(tmp1)*unit_cs
    close(file_tcs)

    open(file_dcs, file = "output/diff_cs.d")
    do j = 0, p3d 
        sum = 0.d0 
        do i = 0, L 
            tmp2 = cos(theta_grid(j))
            sum  = sum +f(i)*fgsl_sf_legendre_Pl(i, tmp2)
        end do 
        write(file_dcs, form_cs) theta_grid(j)*unit_theta, abs(sum)**2.d0*unit_cs
    end do 
    close(file_dcs)
    if(allocated(f)) deallocate(f)
end subroutine PROC_CS_plot
! end cs plot --------------------------------------
end module asymptote
