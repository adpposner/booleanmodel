module generateoperators
    use modelparms
    use combinatorialsubs
    implicit none
    
    contains


!Gets transition probabilities for the "noise step"
!Dim will either be nMess or nMicro depending on application
!This is cumbersome and expensive
subroutine transprob_noise_ext(dim,s1,s2,p_zero,p_one,ptrans)
    integer,intent(in) :: s1,s2,dim
    real(kind=qp),intent(in) :: p_zero,p_one
    real(kind=qp),intent(out) :: ptrans
    real(kind=qp) :: ptemp

    integer :: nkill,n10,n01

    !number of 1's to kill off
    nkill = s1-s2

    ptrans = 0.0D0

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
        !more zeros than ones, so we need to turn on some stuff
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

end subroutine transprob_noise_ext

!Pretty basic method of generating transcriptional noise, expensive and 
!doesn't need to be rerun for each solution
subroutine generate_transcriptional_noise_mess(nMess,pzero,pone,noise_mtx_mess)
	integer,intent(in) :: nMess
    real(kind=qp),intent(in) ::pzero,pone
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(out) :: noise_mtx_mess
	integer :: i,j
	real(kind=qp) :: ptrans
    
    noise_mtx_mess = 0.0d0
		do i=0,nMess
			do j=0,nMess
				call transprob_noise_ext(nMess,i,j,pzero,pone,ptrans)
				noise_mtx_mess(j,i) = ptrans
			enddo
		enddo
		
end subroutine generate_transcriptional_noise_mess

!Pretty basic method of generating transcriptional noise, expensive and 
!doesn't need to be rerun for each solution
subroutine generate_transcriptional_noise_micro(nMicro,pzero,pone,noise_mtx_micro)
	integer,intent(in) :: nMicro
    real(kind=qp),intent(in) ::pzero,pone
    real(kind=qp),dimension(0:MICROMAX,0:MICROMAX),intent(out) :: noise_mtx_micro
	integer :: i,j
	real(kind=qp) :: ptrans
    
    noise_mtx_micro = 0.0d0
		do i=0,nMicro
			do j=0,nMicro
				call transprob_noise_ext(nMicro,i,j,pzero,pone,ptrans)
				noise_mtx_micro(j,i) = ptrans
			enddo
		enddo
		
end subroutine generate_transcriptional_noise_micro



!generate an UT matrix corresponding to the action of a single miRNA on
!the system. 
!Write_mtx_r specifies whether to print out the R matrix for testing
!nMess - number of mRNAs
!R_mat - output matrix
!RNAStrength - a parameter used in the operation
subroutine generate_operator_R_ext(nMess,R_mat,RNAStrength)
	logical,intent(in) :: write_mtx_R
	integer,intent(in) :: nMess
	real(kind=qp),intent(in) :: RNAStrength
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(out) :: R_mat
	integer :: i
	R_mat = 0.0d0

	do i=0,nMess
		call binomialDistribution(i,RNAStrength,R_mat(0:i,i))
	end do

end subroutine generate_operator_R_ext




!Loops operator_C_for_k_ext to generate the C matrices
subroutine generate_operator_C_ext(nMess,nMicro,deg_low,deg_high,p_nz,pmrProb,C_mat_mRNA,C_mat_micro)
    integer,intent(in) :: nMess,nMicro,deg_low,deg_high
    real(kind=qp),intent(in) :: p_nz,pmrProb
	real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(out) :: C_mat_mRNA
	real(kind=qp),dimension(0:MICROMAX,0:MESSMAX),intent(out) :: C_mat_micro
	real(kind=qp) :: pp
    integer :: i
	C_mat_mRNA = 0.0d0
	C_mat_micro = 0.0d0
    

	do i=0,nMess
        call pPlusForModel_ext(i,nMess,deg_low,deg_high,p_nz,pmrProb,pp)
    	call binomialDistribution(nMess,pp,C_mat_mRNA(0:nMess,i))
		call binomialDistribution(nMicro,pp,C_mat_micro(0:nMicro,i))
		!call operator_C_for_k(i,nMess,nMicro,C_mat_mRNA(:,i),C_mat_micro(:,i))
	enddo

end subroutine generate_operator_C_ext

subroutine generate_operators_external(write_mtx,nMess,nMicro,deg_low,deg_high,p_nz,pmrProb,pzero,pone,RNAStrength,noise_mtx_mess, &
        noise_mtx_micro,R_mat,C_mat_mRNA,C_mat_micro)
           

    integer,intent(in) :: nMess,nMicro,deg_low,deg_high
    real(kind=qp),intent(in) :: pzero,pone,RNAStrength,p_nz,pmrProb
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),intent(out) :: noise_mtx_mess,R_mat,C_mat_mRNA
    real(kind=qp),dimension(0:MICROMAX,0:MESSMAX),intent(out) :: C_mat_micro
    real(kind=qp),dimension(0:MICROMAX,0:MICROMAX),intent(out) :: noise_mtx_micro

    call generate_operator_R_ext(nMess,R_mat,RNAStrength)
    call generate_operator_C_ext(nMess,nMicro,deg_low,deg_high,p_nz,pmrProb,C_mat_mRNA,C_mat_micro)
    call generate_transcriptional_noise_mess(nMess,pzero,pone,noise_mtx_mess)
    call generate_transcriptional_noise_micro(nMicro,pzero,pone,noise_mtx_micro)


end subroutine generate_operators_external

! !k is number of active mRNAs
! !pset is parameters
! !V1 is a column of C_mat_mRNA
! !V2 is a column of C_mat_micro
! subroutine operator_C_for_k_ext(k,nMess,nMicro,V1,V2)
! !Operator C is a mapping from 
! !(# nonzero elements of mRNA) -> (#nonzero elements of mRNA) \otimes (#nonzero elems of miRNA)
! !Given the input k = #nonzero mRNAs, it generates the columns of transitions from k
! !which are given by C_ij;k = P(Binom(nMess,pp)=i|k nzelems)*(Binom(nMicro,pp)=j|k nzelems)
! !where pp= the return of pPlusForModel(k,pp). This is expressible as a tensor product, so
! !we leave it so we can cross it with L and lose the nMicro-dim part
! 		integer,intent(in) :: k,nMess,nMicro
! 		real(kind=qp),dimension(0:MESSMAX),intent(out) :: V1 !dim = nMess
! 		real(kind=qp),dimension(0:MICROMAX),intent(out) :: V2 !dim = nMicro
! 		real(kind=qp) :: pp
!         V1 = 0.0d0
!         V2 = 0.0d0
! 		call pPlusForModel(k,pset,pp)
! 		call binomialDistribution(nMess,pp,V1(0:nMess))
! 		call binomialDistribution(nMicro,pp,V2(0:nMicro))
! end subroutine operator_C_for_k_ext



end module generateoperators
