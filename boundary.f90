module boundary
    use kind_type
    use global
    implicit none
    real(dp), save, allocatable, protected :: R(:), K(:)
    complex(dp), save, allocatable, protected :: S(:) 
contains


! ==================================================
! FUNCTIONS
! ==================================================
! full size R matrix -------------------------------
function matching_RL(l1, l2, E)
    use basis, only: HL, EL 
    integer(i4) :: i
    real(dp) :: matching_RL
    real(dp), intent(in) :: E
    integer(i4), intent(in) :: l1, l2 
    real(qp) :: sum
    sum = 0.d0
    do i = 1, N*(L +1)
        sum = sum +HL(N, l1, i)*HL(N, l2, i)/(EL(i) -E)
    end do
    matching_RL = sum/(2.d0*Ba)
end function matching_RL
! end functions (single) ---------------------------










! ==================================================
! SUB-CALCULATION
! ==================================================
! matrix R -----------------------------------------
subroutine SUB_R
    integer(i4) :: i, j 

!     do i = 0, L 
!         do j = 0, L 
!             print *, i, j, matching_RL(i, j, Scatt)
!         end do 
!     end do 
!     stop 

    do j = 0, L 
        R(j) = matching_RL(j, j, Scatt)
    end do 
end subroutine SUB_R
! matrix K -----------------------------------------
subroutine SUB_K
    use fgsl, only: fgsl_sf_bessel_jsl, fgsl_sf_bessel_ysl
    real(dp) :: ka, sb_j, sb_y, diff_j, diff_y 
    real(dp) :: agamma, tmp1, tmp2 
    integer(i4) :: j 
    ka = (2.d0*Scatt)**0.5d0*Ba
    do j = 0, L
        print *, "Here #2-1", j 
        sb_j = fgsl_sf_bessel_jsl(j, ka)
        sb_y = fgsl_sf_bessel_ysl(j, ka)
        if(j /= 0) then 
            diff_j = fgsl_sf_bessel_jsl(j -1_i4, ka) -dble(j +1)/ka*fgsl_sf_bessel_jsl(j, ka)
            diff_y = fgsl_sf_bessel_ysl(j -1_i4, ka) -dble(j +1)/ka*fgsl_sf_bessel_ysl(j, ka)
        else if(j == 0) then 
            diff_j = -fgsl_sf_bessel_jsl(1_i4, ka) -dble(1)/ka*fgsl_sf_bessel_jsl(0_i4, ka)
            diff_y = -fgsl_sf_bessel_ysl(1_i4, ka) -dble(1)/ka*fgsl_sf_bessel_ysl(0_i4, ka)
        end if 
        print *, "Here #2-2", R(j)
        agamma = 1.d0/R(j) -1.d0
        tmp1   = ka*diff_j -agamma*sb_j 
        tmp2   = ka*diff_y -agamma*sb_y 
        K(j)   = tmp1/tmp2
    end do 
end subroutine SUB_K
! matrix S -----------------------------------------
subroutine SUB_S
    use math_const, only: i => math_i 
    integer(i4) :: j 
    do j = 0, L 
        S(j) = (1.d0 +i*K(j))/(1.d0 -i*K(j))
    end do 
end subroutine SUB_S
! end sub-calculation ------------------------------










! ==================================================
! PROCESS
! ==================================================
! boundary matrix ----------------------------------
subroutine PROC_matching
    if(allocated(R)) deallocate(R)
    if(allocated(K)) deallocate(K)
    if(allocated(S)) deallocate(S)
    allocate(R(0:L))
    allocate(K(0:L))
    allocate(S(0:L))
    print *, "Here #1"
    call SUB_R
    print *, "Here #2"
    call SUB_K
    print *, "Here #3"
    call SUB_S 
end subroutine PROC_matching
! out ----------------------------------------------
subroutine PROC_boundary_out
    if(allocated(R)) deallocate(R)
    if(allocated(K)) deallocate(K)
    if(allocated(S)) deallocate(S)
end subroutine PROC_boundary_out
! end out ------------------------------------------
end module boundary
