module blas_interface
    use base_parameters
    use modelparms
    !Not actually F90, but LAPACK95 interface is so much easier than calling functions for work mtx size
    !Be aware that there are type conversions used here to go from qp -> dp -> qp, as we assume we have 
    !qp operands and blas/lapack are built w/double precision
    use lapack95, only : geev
    implicit none

    private :: stationary_distn_2ev

contains

!Pretty self-explanatory - for qp though, must be very careful
subroutine left_mult_by_UT(uppertriangular,rightoperand)

	real(kind=qp),dimension(:,:),INTENT(IN) :: uppertriangular
    real(kind=qp),dimension(:,:),INTENT(INOUT) :: rightoperand
    real(kind=dp),DIMENSION(size(uppertriangular,1),size(uppertriangular,2)) :: uppertriangular_dbl
    real(kind=dp),DIMENSION(size(rightoperand,1),size(rightoperand,2)) :: rightoperand_dbl
	integer :: rowsUT,colsUT,rowsRO,colsRO !the dimensions of both matrices

	rowsUT = size(uppertriangular_dbl,1)
	colsUT = size(uppertriangular_dbl,2)
	rowsRO = size(rightoperand_dbl,1)
	colsRO = size(rightoperand_dbl,2)
    
    uppertriangular_dbl = real(uppertriangular,dp)
    rightoperand_dbl = real(rightoperand,dp)
 
	!just make sure dims match up
	if( colsUT .ne. rowsRO) then
		STOP "DIMS DO NOT MATCH UP IN left_mult_by_UT!"
	else
		call dtrmm('L','U','N','N',rowsRO,colsRO,1.0D0,uppertriangular_dbl,&
			rowsUT,rightoperand_dbl,rowsRO)
	endif

    rightoperand = real(rightoperand_dbl,qp)

end subroutine left_mult_by_UT

subroutine gen_mult_by_mtx(leftoperand,rightoperand,res)

	real(kind=qp),dimension(:,:),INTENT(IN) :: leftoperand
	real(kind=qp),dimension(:,:),INTENT(IN) :: rightoperand
    real(kind=qp),dimension(:,:),INTENT(OUT) :: res
    
    real(kind=dp),dimension(size(leftoperand,1),size(leftoperand,2)) :: leftoperand_dbl
	real(kind=dp),dimension(size(rightoperand,1),size(rightoperand,2)) :: rightoperand_dbl
	real(kind=dp),dimension(size(res,1),size(res,2)) :: res_dbl

	integer :: rowsLO,colsLO,rowsRO,colsRO,rowsRes,colsRes !the dimensions of both matrices

	rowsLO = size(leftoperand_dbl,1)
	colsLO = size(leftoperand_dbl,2)
	rowsRO = size(rightoperand_dbl,1)
	colsRO = size(rightoperand_dbl,2)
	rowsRes = size(res_dbl,1)
	colsRes = size(res_dbl,2)
    
    res = 0.0_qp
    leftoperand_dbl = real(leftoperand,dp)
    rightoperand_dbl = real(rightoperand,dp)
    res_dbl = real(res,dp)
    if (any(res_dbl .ne. 0.0)) then
        stop "Conversion quad to double not yielding 0.0D0, may not have proper quad support!"
    endif
    !just make sure dims match up
	if( colsLO .ne. rowsRO) then
		STOP "COLSLO != ROWSRO gen_mult_by_mtx!"
	else if (rowsRes .ne. rowsLO) then
		STOP "ROWSRES != ROWSLO gen_mult_by_mtx!"
	else if (colsRes .ne. colsRO) then
		STOP "COLSRES != COLSRO gen_mult_by_mtx!"
	else
		call dgemm('N','N',rowsRes,colsRes,colsLO,1.0D0,leftoperand_dbl,&
			rowsLO,rightoperand_dbl,rowsRO,0.0D0,res_dbl,rowsRes)
    endif
    res = real(res_dbl,qp)

end subroutine gen_mult_by_mtx

!Make an identity mtx
subroutine identity_mtx(sz,M)
	integer,intent(in) :: sz
	real(kind=qp),dimension(0:,0:),intent(out) :: M
	integer :: rows,cols,i
	rows = size(M,1)
    cols = size(M,2)
    if (sz > min(rows,cols)) then
        stop "Identity_mtx size exceeds row/col dims!"
    endif
	M=0.0_qp
    do i =0,sz
        M(i,i) = 1.0_qp
    enddo
end subroutine identity_mtx



subroutine stationary_distn(mat,ev,res)
    !integer,INTENT(IN) :: nMessp1
    real(kind=qp),dimension(:,:),INTENT(IN) :: mat
    real(kind=qp),dimension(size(mat,1)),INTENT(OUT) ::res
    real(kind=qp),INTENT(OUT) :: ev
    real(kind=qp),DIMENSION(2) :: evals
    real(kind=qp),DIMENSION(size(mat,1),2) :: evecs
    logical :: e_one_eq_one, e_two_eq_one,same_signs_one,same_signs_two
    call stationary_distn_2ev(mat,evals,evecs)
    e_one_eq_one = .false.
    e_two_eq_one = .false.
    same_signs_one = .false.
    same_signs_two = .false.

    if (abs(evals(1) - 1.0_qp) .lt. 1e-10_qp) then
        e_one_eq_one = .true.
    endif
    if (abs(evals(2) - 1.0_qp) .lt. 1e-10_qp) then
        e_two_eq_one = .true.
    endif
    if (all(evecs(:,1)/sum(evecs(:,1)).gt. -1e-10_qp)) then
        same_signs_one = .true.
    endif
    if (all(evecs(:,2)/sum(evecs(:,2)) .gt. -1e-10_qp)) then
        same_signs_two = .true.
    endif

    ev = 0.0_qp
    res = 0.0_qp
    if (e_one_eq_one .and. same_signs_one) then
        if (e_two_eq_one .and. same_signs_two) then
            ev = evals(2)
            res = evecs(:,2)
        else
            ev = evals(1)
            res = evecs(:,1)
        endif
    else if (e_two_eq_one .and. same_signs_two) then
        ev = evals(2)
        res = evecs(:,2)
    endif
    if (ev .eq. 0.0_qp) then
        print*, "stationary distn could not be found"
    endif
    res = res / sum(res)
end subroutine stationary_distn

subroutine stationary_distn_2ev(mat,eval12,evec12)
    !integer,INTENT(IN) :: nMessp1
    real(kind=qp),dimension(:,:),INTENT(IN) :: mat
    real(kind=qp),dimension(size(mat,1),2),INTENT(OUT) ::evec12
    real(kind=qp),dimension(2),INTENT(OUT) :: eval12
    real(kind=dp),dimension(size(mat,1),size(mat,2)) :: matcopy
    
    real(kind=dp),DIMENSION(size(mat,1)) :: s,si

    real(kind=dp),DIMENSION(size(mat,2),size(mat,2)) :: vt
    integer :: status

    matcopy = real(mat,dp)
    !print*,"SUMMAT",sum(matcopy,1)
    call geev(matcopy,wr=s,wi=si,vr=vt,info=status)
    if (status .ne. 0) then
        evec12 = 0.0_qp
        eval12 = 0.0_qp
    else
        eval12= real(s(1:2),qp)
        evec12= real(vt(:,1:2),qp)
    endif

end subroutine stationary_distn_2ev


end module blas_interface
