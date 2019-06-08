module operatorsolutions
    use base_parameters
use blas_interface
use modelparms
use combinatorialsubs
use zetas
implicit none
	
PRIVATE
PUBLIC :: solve_system
contains



!generate an UT matrix corresponding to the action of a single miRNA on
!the system. 
!Write_mtx_r specifies whether to print out the R matrix for testing
!nMess - number of mRNAs
!R_mat - output matrix
!RNAStrength - a parameter used in the operation
subroutine generate_operator_R(nMess,R_mat,RNAStrength)
	integer,intent(in) :: nMess
	real(kind=qp),intent(in) :: RNAStrength
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(out) :: R_mat
	integer :: i
	R_mat = 0.0_qp

	do i=0,nMess
		call binomialDistribution(i,RNAStrength,R_mat(0:i,i))
	end do
end subroutine generate_operator_R


!k is number of active mRNAs
!pset is parameters
!V1 is a column of C_mat_mRNA
!V2 is a column of C_mat_micro
subroutine operator_C_for_k(k,nMess,nMicro,tf_deg_low,tf_deg_high,tf_p_nz,tf_pmrProb,V1,V2,hyptensor)
!Operator C is a mapping from 
!(# nonzero elements of mRNA) -> (#nonzero elements of mRNA) \otimes (#nonzero elems of miRNA)
!Given the input k = #nonzero mRNAs, it generates the columns of transitions from k
!which are given by C_ij;k = P(Binom(nMess,pp)=i|k nzelems)*(Binom(nMicro,pp)=j|k nzelems)
!where pp= the return of pPlusForModel(k,pp). This is expressible as a tensor product, so
!we leave it so we can cross it with L and lose the nMicro-dim part
		integer,intent(in) :: k
		integer,intent(in) :: nMess,nMicro,tf_deg_low,tf_deg_high
        real(kind=qp),intent(in) :: tf_p_nz,tf_pmrProb
        real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),INTENT(IN) :: hyptensor
		real(kind=qp),dimension(0:MESSMAX),intent(out) :: V1 !dim = nMess
		real(kind=qp),dimension(0:MICROMAX),intent(out) :: V2 !dim = nMicro
		real(kind=qp) :: pp

		!function pPlusForModel(N1,nMess,tf_deg_low,tf_deg_high,tf_p_nz,tf_pmrProb)
		pp= pPlusForModel(k,tf_deg_low,tf_deg_high,tf_p_nz,tf_pmrProb,hyptensor)


		V1=0.0_qp
		V2=0.0_qp
        !print*,"pp(",k,") = ",pp
		call binomialDistribution(nMess,pp,V1(0:nMess))
		call binomialDistribution(nMicro,pp,V2(0:nMicro))
end subroutine operator_C_for_k


!Loops operator_C_for_k to generate the C matrices
subroutine generate_operator_C(nMess,nMicro,tf_deg_low,tf_deg_high,tf_p_nz,tf_pmrProb,C_mat_mRNA,C_mat_micro,hyptensor)
	integer,intent(in) :: nMess,nMicro,tf_deg_low,tf_deg_high
    real(kind=qp),intent(in) :: tf_p_nz,tf_pmrProb
    real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),INTENT(IN) :: hyptensor
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(out) :: C_mat_mRNA
	real(kind=qp),dimension(0:MICROMAX,0:MESSMAX),intent(out) :: C_mat_micro
	integer :: i
	C_mat_mRNA = 0.0_qp
    C_mat_micro = 0.0_qp
   
	do i=0,nMess
		call operator_C_for_k(i,nMess,nMicro,tf_deg_low,tf_deg_high,tf_p_nz,tf_pmrProb,C_mat_mRNA(:,i),C_mat_micro(:,i),hyptensor)
	enddo
	!print*,sum(C_mat_micro,dim=1)
	!call dbgmatsum(pset%nMess,pset%nMess,size(C_mat_mRNA,1),size(C_mat_mRNA,2),C_mat_mRNA)



end subroutine generate_operator_C


!This is the subroutine denoting the effect of a single miRNA on the nMess elts
!nmirtarg_probs is the probabilities of inhibition w.r.t. the number of targets (see zetas.f90)
!M_mat is the single operator
!M_pow_mat is M_mat iteratively multiplied, it is set to the identity here
subroutine generate_operator_M(nMess,m_d_low,m_d_high,nmirtarg_probs,M_mat,M_pow_mat,hyps,pkills)
	integer,intent(in) :: nMess,m_d_low,m_d_high
	real(kind=qp),DIMENSION(0:MESSMAX),intent(in) :: nmirtarg_probs
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: hyps,pkills
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(out) :: M_mat
	real(kind=qp),dimension(MESSMAX+1,0:MESSMAX),intent(out) :: M_pow_mat
	integer :: intrans, outtrans
	real(kind=qp) :: ptrans
	M_mat = 0.0_qp
    M_pow_mat = 0.0_qp

	do intrans=0,nMess
		do outtrans=0,nMess

			ptrans= transprob_zeta(m_d_low,m_d_high,nmirtarg_probs,intrans,outtrans,hyps,pkills)
			M_mat(outtrans,intrans) = ptrans
		enddo
		!call dbgprint(intrans)
	enddo
	!Init M_pow to identity
	call identity_mtx(nMess,M_pow_mat)


end subroutine generate_operator_M

!TxL operations - uses M_matrices and R_matrix
subroutine applyOperator_L(nMess,nMicro,R_mtx,C_mtx_mRNA,C_mtx_micro,M_mtx,M_pow,total_res)
	!We take our matrix for C_mtx_RNA. Each column arises as a linear combination of
	!applications of M to it. So we will have "M_new." It starts at zero, and then
	!we first multiply each column by the 1st row and add to stash
	integer,intent(in) :: nMess,nMicro
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(in) :: R_mtx
	real(kind=qp),dimension(0:MICROMAX,0:MESSMAX),intent(in) :: C_mtx_micro
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(in) :: M_mtx
	real(kind=qp),dimension(MESSMAX+1,0:MESSMAX),intent(inout) :: M_pow
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(inout) :: C_mtx_mRNA
	real(kind=qp),dimension(:,:),intent(out) :: total_res
	integer :: i,j
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX) :: C_scratch
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX) :: combined_operators

	C_scratch = 0.0_qp
	combined_operators = 0.0_qp
	total_res = 0.0_qp
	!Apply operator_R - this happens regardless of miRNAs
	call left_mult_by_UT(R_mtx,C_mtx_mRNA)
	if(nMicro .gt. 0) then
		!On the first iteration (0 miRNAs active), we don't want to update M_pow, it should be
		!the identity going into the coming loop
			do j=0,nMess
				C_scratch(:,j)=C_mtx_micro(0,j)*C_mtx_mRNA(:,j)
			enddo
		!Now we have multiplied the columns of the micro & mRNA vectors and put in C_scratch
		!M_pow is identity so this does nothing, kept for consistency
		call left_mult_by_UT(M_pow,C_scratch)
		
		!Left multiply by C_scratch and add to combined_operators
		combined_operators = C_scratch
		
		!Now do for active miRNAs
		do i=1,nMicro
			!This comes at the beginning since nActiveMicro = 0 already done
			call left_mult_by_UT(M_mtx,M_pow)

		!Multiply each column of C_mtx_mRNA by ith row of C_mtx_micro and stash in C_scratch
		!The loop is faster than a matmul here, but equivalent
			do j=0,nMess
				C_scratch(:,j)=C_mtx_micro(i,j)*C_mtx_mRNA(:,j)
			enddo
		!Now we have multiplied the columns and put in C_scratch

		!Apply dtrmm with M_pow on C_scratch and add it to combined_operators
			call left_mult_by_UT(M_pow,C_scratch)
			combined_operators=combined_operators+C_scratch
		enddo
	!If no miRNAs, jump to here
	else
		combined_operators = C_mtx_mRNA
	endif
	total_res = 0.0_qp
	!Just setting this to make sure it is sized correctly
	total_res(1:size(combined_operators,1),1:size(combined_operators,2)) = combined_operators
	
end subroutine applyOperator_L
	




!total solver
subroutine solve_system(write_matrices,pset,total_res,nmirtarg_probs,hyps,pkills,noise_mtx_mess,noise_mtx_micro)

	logical,intent(in) :: write_matrices
	type(modelparams),intent(in) :: pset
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),INTENT(IN):: noise_mtx_mess
    real(kind=qp),dimension(0:MICROMAX,0:MICROMAX),INTENT(IN) :: noise_mtx_micro
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: hyps,pkills
	real(kind=qp) :: tf_pmrProb
	real(kind=qp),dimension(:,:),intent(inout) :: total_res
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX):: M_mat
	real(kind=qp),dimension(MESSMAX+1,0:MESSMAX):: M_pow_mat
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX) :: R_mat
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX) :: C_mat_mRNA
    real(kind=qp),dimension(0:MICROMAX,0:MESSMAX) :: C_mat_micro
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX) :: noised_C_op_mess
    real(kind=qp),dimension(0:MICROMAX,0:MESSMAX) :: noised_C_op_micro
	real(kind=qp),DIMENSION(0:MESSMAX),INTENT(IN):: nmirtarg_probs
	integer :: nMess,nMicro,tf_deg_low,tf_deg_high,m_deg_low,m_deg_high
	real(kind=qp) :: tf_p_nz,p_zero,p_one,tmpr,mirstrength,rho
	if (pset%nMess .gt. MESSMAX) then
		STOP "NMESS TOO BIG!"
	endif
	if (pset%nMicro .gt. MICROMAX) then
		STOP "NMICRO TOO BIG!"
	endif
	!call zetas_init(pset,nmirtarg_probs)
	M_mat = 0.0_qp
	M_pow_mat = 0.0_qp
	R_mat = 0.0_qp
	C_mat_micro =0.0_qp
    C_mat_mRNA = 0.0_qp
    noised_C_op_mess = 0.0_qp
    noised_C_op_micro = 0.0_qp
	nMess = pset%nMess
	nMicro = pset%nMicro
	tf_deg_low = pset%tf_deg_low
    tf_deg_high = pset%tf_deg_high
    tf_p_nz = pset%tf_p_nz
	m_deg_low = pset%m_deg_low
	m_deg_high = pset%m_deg_high
	p_zero=pset%noise_pzero
	p_one = pset%noise_pone
	tmpr = pset%tmpr
	tf_pmrProb = tmpr / (1.0_qp + tmpr)
	mirstrength = 1.0_qp - pset%defect
	rho  =pset%rho

	call generate_operator_C(nMess,nMicro,tf_deg_low,tf_deg_high, & 
		tf_p_nz,tf_pmrProb,C_mat_mRNA,C_mat_micro,hyps)
	!call dbgprint(1)
	!call dbgmatsum(pset%nMess,pset%nMess,size(C_mat_mRNA,1),size(C_mat_mRNA,2),C_mat_mRNA)

	if (write_matrices .eqv. .TRUE.) then
		call write_mtx(C_mat_mRNA,"testoutput/C_mtx_mRNA.mtx.txt")
		call write_mtx(C_mat_micro,"testoutput/C_mtx_micro.mtx.txt")
	endif

    !Add noise here
    call gen_mult_by_mtx(noise_mtx_mess,C_mat_mRNA,noised_C_op_mess)
    call gen_mult_by_mtx(noise_mtx_micro,C_mat_micro,noised_C_op_micro)
    C_mat_mRNA = noised_C_op_mess
    C_mat_micro = noised_C_op_micro

    if (write_matrices .eqv. .TRUE.) then
		call write_mtx(C_mat_mRNA,"testoutput/C_mtx_mRNA_noised.mtx.txt")
		call write_mtx(C_mat_micro,"testoutput/C_mtx_micro_noised.mtx.txt")
	endif
	
	!call dbgprint(2)

	call generate_operator_M(nMess,m_deg_low,m_deg_high,nmirtarg_probs, &
		M_mat, M_pow_mat,hyps,pkills)

	if (write_matrices .eqv. .TRUE.) then
		call write_mtx(M_mat,"testoutput/operator_M.mtx.txt")
		call write_mtx(M_pow_mat,"testoutput/M_pow_ident.mtx.txt")
	endif

	!call dbgprint(3)
	call generate_operator_R(nMess,R_mat,rho)
	if (write_matrices .eqv. .TRUE.) then
		!print*,shape(R_mat)
		call write_mtx(R_mat,"testoutput/operator_R.mtx.txt")
	endif
	call applyOperator_L(nMess,nMicro,R_mat,C_mat_mRNA, &
		C_mat_micro,M_mat,M_pow_mat,total_res)
	if (write_matrices .eqv. .TRUE.) then
		call write_mtx(total_res,"testoutput/sys_sln.mtx.txt")
	endif

end subroutine solve_system


end module operatorsolutions
