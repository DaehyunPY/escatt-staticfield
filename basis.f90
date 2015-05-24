module basis
    use kind_type 
    use global 
    implicit none
    real(dp), save, allocatable, protected :: EL(:) 
    real(dp), save, pointer, protected :: HL(:, :, :)
contains


! ==================================================
! HARMILTONIAN TERMS 
! ==================================================
! kinetic term -------------------------------------
function term_kinet(i, j)
    use hamiltonian, only: Delta_grid, dr_pdrho
    integer(i4), intent(in) :: i, j 
    real(dp) :: term_kinet
    term_kinet = -Delta_grid(i, j)/dr_pdrho**2.d0/2.d0 
end function term_kinet 
! potential term -----------------------------------
function term_poten(i)
    use hamiltonian, only: r_grid, poten_r 
    integer(i4), intent(in) :: i
    real(dp) :: term_poten
    term_poten = poten_r(r_grid(i))
end function term_poten 
! angular term -------------------------------------
function term_angular(i, l)  
    use hamiltonian, only: r_grid
    integer(i4), intent(in) :: i, l 
    real(dp) :: term_angular
    term_angular = dble(l)*(dble(l) +1.d0)/r_grid(i)**2.d0/2.d0 
end function term_angular
! dipole l +1 term ---------------------------------
function term_dipole(i, l)
    use hamiltonian, only: r_grid
    integer(i4), intent(in) :: i, l 
    real(dp) :: term_dipole, filter, alpha = 0.003d0, cut
!     term_dipole = dble(l +1)/(dble(2*l +1)*dble(2*l +3))**0.5d0*Amp*r_grid(i)
    cut = Ba*0.75d0 
    if(r_grid(i) <= cut) then 
!         filter = 1.d0 
        filter = (1.d0 -exp(-alpha*(r_grid(i) -cut)**2.d0))
        term_dipole = dble(l +1)/(dble(2*l +1)*dble(2*l +3))**0.5d0*Amp &
                        *r_grid(i)*filter
    else 
        term_dipole = 0.d0 
    end if 
end function term_dipole
! end funciton -------------------------------------










! ==================================================
! SUB-CALCULATE
! ==================================================
! full size hamiltonian ----------------------------
subroutine SUB_full_angular
    use linear, only: diag_sym
    use plot, only: plot_mat, plot_vec 
    real(dp), pointer :: HL_p0(:), HL_p1(:, :, :, :), HL_p2(:, :)
    integer(i4) :: i, j, k 
    nullify(HL_p0)
    nullify(HL_p1)
    nullify(HL_p2)
    allocate(HL_p0(1:N*(L +1)*N*(L +1)))
    allocate(HL_p1(1:N, 0:L, 1:N, 0:L))
    allocate(HL_p2(1:N*(L +1), 1:N*(L +1)))
    HL(1:N, 0:L, 1:N*(L +1)) => HL_p0(1:N*(L +1)*N*(L +1))
    HL_p1(1:N, 0:L, 1:N, 0:L) => HL_p0(1:N*(L +1)*N*(L +1))
    HL_p2(1:N*(L +1), 1:N*(L +1)) => HL_p0(1:N*(L +1)*N*(L +1))

    ! note: x1, angular block, x2, angular block 
    HL_p1(:, :, :, :) = 0.d0 
    do k = 0, L 
        do j = 1, N 
            do i = 1, j  
                HL_p1(i, k, j, k) = term_kinet(i, j)
            end do 
            HL_p1(j, k, j, k) = HL_p1(j, k, j, k) +term_poten(j) +term_angular(j, k) 
            if(.not. k == L) then 
                HL_p1(j, k, j, k +1) = term_dipole(j, k) 
            end if 
        end do 
    end do 
    call plot_mat(50, HL_p2)
    call diag_sym(HL_p2, EL)
    call plot_vec(51, EL)
    call plot_mat(52, HL_p2)
    nullify(HL_p0)
    nullify(HL_p1)
    nullify(HL_p2)
end subroutine SUB_full_angular
! descale full size wave function ------------------
subroutine SUB_descale_HL
    use hamiltonian, only: weight_grid, dr_pdrho
    integer(i4) :: i, j, k 
    real(dp) :: sum 

!     do j = 1, N*(L +1) 
!         do i = 0, L 
!             sum = 0.d0 
!             do k = 1, N 
!                 sum = sum +abs(HL(k, i, j))**2.d0 
!             end do 
!             if(1.d0 +dble(sum) == 1.d0) then 
!                 HL(1:N, i, j) = 0.d0 
!             else 
!                 HL(1:N, i, j) = HL(1:N, i, j)/sum**0.5d0
!             end if 
!         end do 
!     end do 

    do j = 1, N*(L +1) 
        do i = 0, L 
            HL(1:N, i, j) = HL(1:N, i, j)/(weight_grid(1:N)*dr_pdrho)**0.5d0 
        end do 
    end do 
end subroutine SUB_descale_HL
! end sub-calculate --------------------------------










! ==================================================
! PROCESS
! ==================================================
! hamiltonian --------------------------------------
subroutine PROC_H
    character(30), parameter :: form_out1 = '(1A15, 5F9.3)', form_out2 = '(1A15, 1ES15.3, 1ES15.3)'
    nullify(HL)
    if(allocated(EL)) deallocate(EL)
    allocate(HL(1:N, 0:L, 1:N*(L +1)))
    allocate(EL(1:N*(L +1)))
    call SUB_full_angular
    call SUB_descale_HL
    write(file_log, form_out2) "E: ", minval(EL(:)), maxval(EL(:)) 
end subroutine PROC_H
! break --------------------------------------------
subroutine PROC_basis_break
    nullify(HL)
    if(allocated(EL)) deallocate(EL)
end subroutine PROC_basis_break
! out ----------------------------------------------
subroutine PROC_basis_out
    nullify(HL)
    if(allocated(EL)) deallocate(EL)
end subroutine PROC_basis_out
! end process --------------------------------------
end module basis
