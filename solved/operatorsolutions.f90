module operatorsolutions
	use blas_interface
	use modelparms
	use combinatorialsubs
	use zetas
	implicit none
	
	

contains



!generate an UT matrix corresponding to the action of a single miRNA on
!the system. 
!Write_mtx_r specifies whether to print out the R matrix for testing
!nMess - number of mRNAs
!R_mat - output matrix
!RNAStrength - a parameter used in the operation
subroutine generate_operator_R(write_mtx_R,nMess,R_mat,RNAStrength)
	logical,intent(in) :: write_mtx_R
	integer,intent(in) :: nMess
	real(kind=qp),intent(in) :: RNAStrength
	real(kind=qp),dimension(0:,0:),intent(out) :: R_mat
	integer :: i
	R_mat = 0.0_qp

	do i=0,nMess
		call binomialDistribution(i,RNAStrength,R_mat(0:i,i))
	end do

	if (write_mtx_R .eqv. .TRUE.) then
		call write_mtx(R_mat,"operator_R")
	endif

end subroutine generate_operator_R


!k is number of active mRNAs
!pset is parameters
!V1 is a column of C_mat_mRNA
!V2 is a column of C_mat_micro
subroutine operator_C_for_k(k,pset,V1,V2)
!Operator C is a mapping from 
!(# nonzero elements of mRNA) -> (#nonzero elements of mRNA) \otimes (#nonzero elems of miRNA)
!Given the input k = #nonzero mRNAs, it generates the columns of transitions from k
!which are given by C_ij;k = P(Binom(nMess,pp)=i|k nzelems)*(Binom(nMicro,pp)=j|k nzelems)
!where pp= the return of pPlusForModel(k,pp). This is expressible as a tensor product, so
!we leave it so we can cross it with L and lose the nMicro-dim part
		integer,intent(in) :: k
		type(parmset),intent(in) :: pset
		real(kind=qp),dimension(0:MESSMAX),intent(out) :: V1 !dim = nMess
		real(kind=qp),dimension(0:MICROMAX),intent(out) :: V2 !dim = nMicro
		real(kind=qp) :: pp

		call pPlusForModel(k,pset,pp)


		V1=0.0_qp
		V2=0.0_qp

		call binomialDistribution(pset%nMess,pp,V1(0:pset%nMess))
		call binomialDistribution(pset%nMicro,pp,V2(0:pset%nMicro))
	end subroutine operator_C_for_k


!Loops operator_C_for_k to generate the C matrices
subroutine generate_operator_C(write_mtx_C,pset,C_mat_mRNA,C_mat_micro)
	logical,intent(in) :: write_mtx_C
	type(parmset),intent(in) :: pset

	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(out) :: C_mat_mRNA
	real(kind=qp),dimension(0:MICROMAX,0:MESSMAX),intent(out) :: C_mat_micro
	integer :: i
	C_mat_mRNA = 0.0_qp
	C_mat_micro = 0.0_qp
	do i=0,pset%nMess
		call operator_C_for_k(i,pset,C_mat_mRNA(:,i),C_mat_micro(:,i))
	enddo

	!call dbgmatsum(pset%nMess,pset%nMess,size(C_mat_mRNA,1),size(C_mat_mRNA,2),C_mat_mRNA)

	if (write_mtx_C .eqv. .TRUE.) then
		call write_mtx(C_mat_mRNA,"C_mtx_mRNA")
		call write_mtx(C_mat_micro,"C_mtx_micro")
	endif

end subroutine generate_operator_C


!This is the subroutine denoting the effect of a single miRNA on the nMess elts
!nmirtarg_probs is the probabilities of inhibition w.r.t. the number of targets (see zetas.f90)
!M_mat is the single operator
!M_pow_mat is M_mat iteratively multiplied, it is set to the identity here
subroutine generate_operator_M(write_mtx_M,pset,nmirtarg_probs,M_mat,M_pow_mat)

	logical,intent(in) :: write_mtx_M
	type(parmset),intent(in) :: pset
	real(kind=qp),DIMENSION(0:MESSMAX),intent(in) :: nmirtarg_probs
	real(kind=qp),dimension(0:pset%nMess,0:pset%nMess),intent(out) :: M_mat
	real(kind=qp),dimension(pset%nMess+1,0:pset%nMess),intent(out) :: M_pow_mat
	integer :: intrans, outtrans
	real(kind=qp) :: ptrans,rsum
	M_mat = 0.0_qp
	M_pow_mat = 0.0_qp
	do intrans=0,pset%nMess
		do outtrans=0,pset%nMess
			call transprob_zeta(pset,nmirtarg_probs,intrans,outtrans,ptrans)
			M_mat(outtrans,intrans) = ptrans
		enddo
		!call dbgprint(intrans)
	enddo
	!Init M_pow to identity
	call identity_mtx(pset%nMess,M_pow_mat)

	if (write_mtx_M .eqv. .TRUE.) then
		call write_mtx(M_mat,"operator_M")
	endif
end subroutine generate_operator_M


subroutine applyOperator_L(pset,R_mtx,C_mtx_mRNA,C_mtx_micro,M_mtx,M_pow,total_res)
	!We take our matrix for C_mtx_RNA. Each column arises as a linear combination of
	!applications of M to it. So we will have "M_new." It starts at zero, and then
	!we first multiply each column by the 1st row and add to stash
	type(parmset),intent(in) :: pset
	real(kind=qp),dimension(0:pset%nMess,0:pset%nMess),intent(in) :: R_mtx
	real(kind=qp),dimension(0:pset%nMicro,0:pset%nMess),intent(in) :: C_mtx_micro
	real(kind=qp),dimension(0:pset%nMess,0:pset%nMess),intent(in) :: M_mtx
	real(kind=qp),dimension(pset%nMess+1,0:pset%nMess),intent(inout) :: M_pow
	real(kind=qp),dimension(0:pset%nMess,0:pset%nMess),intent(inout) :: C_mtx_mRNA
	real(kind=qp),dimension(:,:),intent(out) :: total_res
	integer :: i,j
	real(kind=qp),dimension(0:pset%nMess,0:pset%nMess) :: C_scratch
	real(kind=qp),dimension(0:pset%nMess,0:pset%nMess) :: combined_operators

	C_scratch = 0.0_qp

	!Apply operator_R - this happens regardless of miRNAs
	call left_mult_by_UT(R_mtx,C_mtx_mRNA)
	if(pset%nMicro .gt. 0) then
		!On the first iteration (0 miRNAs active), we don't want to update M_pow, it should be
		!the identity going into the coming loop
			do j=0,pset%nMess
				C_scratch(:,j)=C_mtx_micro(0,j)*C_mtx_mRNA(:,j)
			enddo
		!Now we have multiplied the columns of the micro & mRNA vectors and put in C_scratch
		!M_pow is identity so this does nothing, kept for consistency
		call left_mult_by_UT(M_pow,C_scratch)
		
		!Left multiply by C_scratch and add to combined_operators
		combined_operators = C_scratch
		
		!Now do for active miRNAs
		do i=1,pset%nMicro
			!This comes at the beginning since nActiveMicro = 0 already done
			call left_mult_by_UT(M_mtx,M_pow)

		!Multiply each column of C_mtx_mRNA by ith row of C_mtx_micro and stash in C_scratch
		!The loop is faster than a matmul here, but equivalent
			do j=0,pset%nMess
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
	

	!The "noise matrix operation"
	!transcriptional noise is applied in the following way
	!1) Operator_C is calculated
	!2) Prior to applying operator R (and later M^q, powers of M), we insert transcription_noise_mtx
	!That is, make a mRNA-noise and a miRNA-noise matrix and left multiply the elements of C
	!3) Proceed as usual

subroutine add_transcriptional_noise(write_noise_mtx,pset,C_mtx_mRNA,C_mtx_micro)
	logical,intent(in) :: write_noise_mtx
	type(parmset),intent(in) :: pset
	real(kind=qp),dimension(:,:),intent(inout) :: C_mtx_mRNA,C_mtx_micro
	real(kind=qp),dimension(0:pset%nMess,0:pset%nMess) :: noised_C_op_mess
	real(kind=qp),dimension(0:pset%nMess,0:pset%nMess) :: noise_mtx_mess
	real(kind=qp),dimension(0:pset%nMicro,0:pset%nMess) :: noised_C_op_micro
	real(kind=qp),dimension(0:pset%nMicro,0:pset%nMicro) :: noised_mtx_micro


		!First, generate noise mtx
		call generate_transcriptional_noise_operator(pset,noise_mtx_mess)
		!Apply it to C_mtx_mRNA
		call gen_mult_by_mtx(noise_mtx_mess,C_mtx_mRNA,noised_C_op_mess)
		!Overwrite C_mtx_mRNA
		C_mtx_mRNA = noised_C_op_mess
		if (write_noise_mtx .eqv. .TRUE.) then
			call write_mtx(noise_mtx_mess,"transcriptional_noise_mtx_mRNA")
		endif
		!Do the same for C_micro
		call generate_transcriptional_noise_operator(pset,noised_mtx_micro)
		
		call gen_mult_by_mtx(noised_mtx_micro,C_mtx_micro,noised_C_op_micro)
		!overwrite C_mtx_micro
		C_mtx_micro = noised_C_op_micro
		if (write_noise_mtx .eqv. .TRUE.) then
			!call write_mtx(noised_mtx_micro,"transcriptional_noise_mtx_micro")
		endif

end subroutine add_transcriptional_noise


!Gets transition probabilities for the "noise step"
!Dim will either be nMess or nMicro depending on application
!This is cumbersome and expensive
subroutine transprob_noise(dim,s1,s2,p_zero,p_one,ptrans)
	integer :: s1,s2,dim
	real(kind=qp) :: ptrans,ptemp,p_zero,p_one

	integer :: nkill,n10,n01

	!number of 1's to kill off
	nkill = s1-s2

	ptrans = 0.0_qp

	if (nkill .ge. 0) then
		!This is how many ones turn to zeros
		do n10=nkill,s1
			!so then n10-n01=nkill
			n01=n10-nkill
			!turn this many of the 1's into zeros
			ptemp = binomialProbability(s1,n10,p_zero)
			!turn this many of the 0's into ones and multiply
			ptemp = ptemp * binomialProbability(dim-s1,n01,p_one)
				if(s1 .eq. dim) then
		
            	endif
			ptrans = ptrans+ptemp
		enddo
	else
		!more zeros than ones, so we need to turn on some elems
		!the min we need to turn on is -1*nkill; e.g.  14-16 = -2
		!so n01 starts at 2 and goes up to the total # of zeros that can
		!be changed, i.e. (dim-s1)
		do n01=(-1*nkill),(dim-s1)
			n10 = n01 + nkill
			ptemp = binomialProbability(s1,n10,p_zero)
			ptemp = ptemp * binomialProbability(dim-s1,n01,p_one)
			ptrans=ptrans+ptemp
		enddo
	endif

end subroutine transprob_noise

!Pretty basic method of generating transcriptional noise, expensive and 
!doesn't need to be rerun for each solution
subroutine generate_transcriptional_noise_operator(pset,noise_mtx)
	type(parmset),intent(in) :: pset
	integer :: dim,i,j
	real(kind=qp) :: ptrans
	real(kind=qp),dimension(0:,0:) :: noise_mtx

	dim= size(noise_mtx,1)-1

	if (pset%transcription_noise_present .eqv. .TRUE.) then
		do i=0,dim
			do j=0,dim
				call transprob_noise(dim,i,j,pset%transcription_p0,pset%transcription_p1,ptrans)
				noise_mtx(j,i) = ptrans
			enddo
		enddo
		
	endif
end subroutine generate_transcriptional_noise_operator



!total solver
subroutine solve_system(write_matrices,pset,total_res)

	logical,intent(in) :: write_matrices
	type(parmset),intent(in) :: pset

	real(kind=qp),dimension(:,:),intent(inout) :: total_res
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX):: M_mat
	real(kind=qp),dimension(MESSMAX+1,0:MESSMAX):: M_pow_mat
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX) :: R_mat
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX) :: C_mat_mRNA
	real(kind=qp),dimension(0:MICROMAX,0:MESSMAX) :: C_mat_micro
	real(kind=qp),DIMENSION(0:MESSMAX):: nmirtarg_probs
	if (pset%nMess .gt. MESSMAX) then
	!STOP "NMESS TOO BIG!"
	endif
	if (pset%nMicro .gt. MICROMAX) then
	!STOP "NMICRO TOO BIG!"
	endif
	call zetas_init(pset,nmirtarg_probs)
	M_mat = 0.0_qp
	M_pow_mat = 0.0_qp
	R_mat = 0.0_qp
	C_mat_micro =0.0_qp
	C_mat_mRNA = 0.0_qp
	call generate_operator_C(write_matrices,pset,C_mat_mRNA,C_mat_micro)
	!call dbgprint(1)
	!call dbgmatsum(pset%nMess,pset%nMess,size(C_mat_mRNA,1),size(C_mat_mRNA,2),C_mat_mRNA)
	call add_transcriptional_noise(write_matrices,pset,C_mat_mRNA(0:pset%nMess,0:pset%nMess),C_mat_micro(0:pset%nMicro,0:pset%nMess))
	!call dbgprint(2)
	call generate_operator_M(write_matrices,pset,nmirtarg_probs,M_mat(0:pset%nMess,0:pset%nMess),M_pow_mat(1:pset%nMess+1,0:pset%nMess))
	!call dbgprint(3)
	call generate_operator_R(write_matrices,pset%nMess,R_mat(0:pset%nMess,0:pset%nMess),pset%RNAStrength)
	
	call applyOperator_L(pset,R_mat(0:pset%nMess,0:pset%nMess),C_mat_mRNA(0:pset%nMess,0:pset%nMess),C_mat_micro(0:pset%nMicro,0:pset%nMess),M_mat(0:pset%nMess,0:pset%nMess),M_pow_mat(1:pset%nMess+1,0:pset%nMess),total_res)


	!call write_mtx(total_res(0:nMess,0:nMess),"total_sln")
end subroutine solve_system


end module operatorsolutions
