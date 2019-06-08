module combinatorialsubs
use base_parameters
implicit none

PRIVATE
public :: binomialDistribution
public :: make_hyptensor
public :: check_hyptensor
public :: make_mirtarg_probs
PUBLIC :: make_pkill_tensor
PUBLIC :: check_mirtarg_probs
PUBLIC :: check_pkill_tensor
public :: pPlusForModel
public :: generate_mess_noise
public :: generate_micro_noise

contains



!Multiplication method for binomial dist'n
function binomialProbability(n,x,p)
    integer :: x,n,xcopy
    real(kind=qp) :: p,pcopy,oneminusp
    real(kind=qp) :: binomialProbability
    integer :: a,b,c

    if ((x .gt. n) .or. (x .lt. 0)) then
        binomialProbability = 0.0_qp
    else
        !Make sure we do as few moves as possible
        if ((2*x) .gt. n) then
            xcopy = n-x
            pcopy = 1.0_qp-p
            oneminusp = p
        else
            xcopy = x
            pcopy = p
            oneminusp = 1.0_qp - p
        endif

        a = 0
        b = 0
        c = 0
        binomialProbability = 1.0_qp

        do while ( (max(a,b) .lt. xcopy) .or. (c .lt. n-xcopy))

        if ( (a .lt. xcopy) .and. (binomialProbability .lt. 1.0_qp)) then
            a=a+1
            binomialProbability = binomialProbability*DBLE(n-xcopy+a)/DBLE(a)
        else if (b < xcopy) then
            b = b+1
            binomialProbability = binomialProbability * pcopy
        else
            c = c + 1
            binomialProbability = binomialProbability * oneminusp
        endif

        enddo

    endif
end function binomialProbability

!Loop of the above
subroutine binomialDistribution(n,p,distn)
    integer,INTENT(IN) :: n
    real(kind=qp),INTENT(IN) :: p
    real(kind=qp),dimension(0:n),INTENT(OUT) :: distn

    !loopvar
    integer :: i

    do i=0,n
    distn(i) = binomialProbability(n,i,p)
    enddo
end subroutine binomialDistribution


!hypergeometric function
function hyp2(Nc,N1c,N2c,kc)
    integer,INTENT(IN) :: Nc,N1c,N2c,kc
    integer :: N,N1,N2,k,NN1,N2k

    real(kind=qp) :: hyp2

    integer :: arrIndex,i,j,numerIdx,denomIdx
    integer,dimension(2*Nc) :: numerator,denominator

    hyp2 = 1.0_qp
    !If we don't have max(0,N1+N2-N) <= k <= min(N1,N2), coeff is = 0
    N=Nc
    if((max(0,N1c+N2c-Nc) .gt. kc) .or. (kc .gt. (min(N1c,N2c)))) then
        hyp2 = 0.0_qp
    else
    !debug, check init'l params
    !print*,"hyp2: Nc= ",Nc," N1c = ",N1c," N2c= ",N2c," kc = ",kc

        N=Nc
        N1=N1c
        NN1 = Nc-N1c
        N2k = N2c-kc

    !for denom. factorial
        if(N2c .gt. (Nc/2)) then
            N2 = Nc-N2c
        else
            N2=N2c
        endif

        !numerator factorials
        if (kc .gt. (N1c/2)) then
            k =  N1c-kc
        else
            k=kc
        endif
        if ((N2c-kc) .gt. ((Nc-N1c)/2)) then
            N2k = Nc-N1c - N2c +kc
        else
            N2k = N2c-kc
        endif

        !debug, check new params
        !print*,"hyp2: N= ",N," N1 = ",N1," N2= ",N2," k = ",k

        numerator = 1
        denominator = 1

        arrIndex = 0
        !fill in numerator and denominator
        do i=1,k
            arrIndex = arrIndex+1
            numerator(arrIndex)=max(N1-k+i,1)
            denominator(arrIndex)=i
        enddo
        do i=1,N2k
            arrIndex = arrIndex+1
            numerator(arrIndex)=max(NN1-N2k+i,1)
            denominator(arrIndex)=i
        enddo
        do i=1,N2
            arrIndex = arrIndex+1
            numerator(arrIndex) = i
            denominator(arrIndex) = max(N-N2+i,1)
        enddo
            !print*,"hyp2: calculating (",N1,",",k,")*(",NN1,",",N2k,") / (",N,",",N2,")"
            !print*,"hyp2: numerator = ",numerator
            !print*,"hyp2: denominator = ",denominator
            !print*,"hyp2: arrIndex = ",arrIndex
        !Do a cancellation
        do i=1,arrIndex
            do j=1,arrIndex
            if (numerator(i) .eq. denominator(j)) then
                numerator(i)=1
                denominator(j)=1
            endif
            enddo
        enddo

        numerIdx=1
        denomIdx=1
        do while (min(numerIdx,denomIdx) .le. arrIndex)
            !print*,hyp2
            if ((hyp2 .le. 1.0_qp) .and. (numerIdx .le. arrIndex) ) then
                !print*,DBLE(numerator(numerIdx))
                hyp2 = hyp2 * numerator(numerIdx)
                numerIdx = numerIdx+1
            else
                !print*,DBLE(denominator(numerIdx))
                hyp2 = hyp2 / denominator(denomIdx)
                denomIdx = denomIdx+1
            endif
        enddo
    endif
end function hyp2



!probability of N2 nonzero elements in row A given deg_lo,deg_hi,p_nz,p_nzbar
function pNZN2(N2,dl,dh,pnz)
    integer,INTENT(IN) :: N2,dl,dh
    real(kind=qp),INTENT(IN) :: pnz
    real(kind=qp) ::pNZN2
    !zero conditions - assume dl & dh both <= N
    !test if N2 between dl dh and pnz != 0
    if ((N2 .lt. dl) .or. (N2 .gt.dh)) then
        pNZN2 = 0.0_qp
    else if ((min(abs(pnz),abs(1.0_qp-pnz))) .lt. epsilon(pnz)) then
        pNZN2 = 0.0_qp
    else
        pNZN2 = binomialProbability(dh-dl,N2-dl,pnz)
    endif
    !print*,"N2 = ",N2," dl = ",dl, " dh = ",dh, " pnz = ",pnz
  end function pNZN2

!probability function evaluates to >0, i.e. P+
!given k, p, pbar only.
function pPlusgivenk(k,p)
    integer,INTENT(IN) :: k
    real(kind=qp),INTENT(IN) ::p
    integer :: l,lmin
    real(kind=qp) :: pPlusgivenk

    pPlusgivenk = 0
    !since we only want l > k/2, we actually don't
    !need to test for even/oddness
    lmin = k/2 + 1
    do l=lmin,k
        pPlusgivenk = pPlusgivenk + binomialProbability(k,l,p)
    enddo
end function pPlusgivenk

function pPlusgivenN2(N1,N2,p,hyptensor)
    integer,INTENT(IN) :: N1,N2
    real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: hyptensor
    real(kind=qp),INTENT(IN) :: p
    integer :: k,kmax
    real(kind=qp) :: hypNN1N2k,ppgk
    real(kind=qp) :: pPlusgivenN2

    
    kmax=min(N1,N2)
    pPlusgivenN2 = 0.0_qp

    do k=1,kmax
        !htest = hyp2(N,N1,N2,k)
        hypNN1N2k = hyptensor(N1,N2,k)
        ! if (abs(hypNN1N2k - htest) .gt. epsilon(1.0_sp)) then
        !     print*,"ppgivenn2 N = ",N," N1 = ",N1," N2= ",N2," k = ",k
        !     print*,"direct = ",hypNN1N2k," indirect = ",htest
        !     stop
        ! endif
        ppgk = pPlusgivenk(k,p)
        pPlusgivenN2 = pPlusgivenN2 + hypNN1N2k*ppgk
    end do

end function pPlusgivenN2


function pPlusgivendegs(N1,dl,dh,pnz,pmr,hyptensor)
    integer,INTENT(IN) :: N1,dl,dh
    real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: hyptensor
    real(kind=qp),INTENT(IN) :: pnz,pmr
    real(kind=qp) :: pNZN2val,pOnN2
    real(kind=qp) :: pPlusgivendegs
    integer  :: N2
    pPlusgivendegs = 0.0_qp

    do N2=dl,dh
        pNZN2val = pNZN2(N2,dl,dh,pnz)
        pOnN2 = pPlusgivenN2(N1,N2,pmr,hyptensor)
        
        !print*,"pnzn2=",pNZN2val," pOnN2=",pOnN2
        pPlusgivendegs = pPlusgivendegs + pNZN2val*pOnN2
    end do
end function pPlusgivendegs


!This generates our probability of transition for a given number of nonzero inputs N1
function pPlusForModel(N1,tf_deg_low,tf_deg_high,tf_p_nz,tf_pmrProb,hyptensor)
    integer,intent(in) :: N1,tf_deg_low,tf_deg_high
    real(kind=qp),INTENT(IN) :: tf_p_nz,tf_pmrProb
    real(kind=qp):: pPlusForModel
    real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: hyptensor
    pPlusForModel = pPlusgivendegs(N1,tf_deg_low,tf_deg_high,tf_p_nz,tf_pmrProb,hyptensor)

end function pPlusForModel


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
subroutine generate_transcriptional_noise_operator(dim,pzero,pone,noise_mtx)
    real(kind=qp),intent(in) :: pzero,pone
    integer,INTENT(IN) :: dim
	integer :: i,j
	real(kind=qp) :: ptrans
	real(kind=qp),dimension(0:,0:),INTENT(OUT) :: noise_mtx

    if (dim .gt. ubound(noise_mtx, dim=1)) then
        print*,"dim = ",dim," noise_mtx ubounds = ",ubound(noise_mtx,dim=1),ubound(noise_mtx,dim=2)
        stop
    endif
	noise_mtx = 0.0_qp
	! if ( ubound(noise_mtx,1) .ne. ubound(noise_mtx,2) ) then
	! 	print*,"NOISE MATRIX NOT SQUARE - generate_transcriptional_noise_operator:263"
	! end if
		do i=0,dim
			do j=0,dim
				!transprob_noise(dim,s1,s2,p_zero,p_one,ptrans)
				call transprob_noise(dim,i,j,pzero,pone,ptrans)
				noise_mtx(j,i) = ptrans
			enddo
		enddo
		
end subroutine generate_transcriptional_noise_operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The following are used to generate the persistent (workspace) variables, w/memory malloc'd
!on C side.
!Check functions are used for checking PERSISTENCE, not validity
!!!!!!!!!!!

subroutine generate_mess_noise(nMess,pzero,pone,noise_mtx_mess)
    INTEGER,INTENT(IN) :: nMess
    real(kind=qp),INTENT(IN) :: pzero,pone
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX),INTENT(OUT) :: noise_mtx_mess

    call generate_transcriptional_noise_operator(nMess,pzero,pone,noise_mtx_mess)
end subroutine generate_mess_noise

subroutine generate_micro_noise(nMicro,pzero,pone,noise_mtx_micro)
    INTEGER,INTENT(IN) :: nMicro
    real(kind=qp),INTENT(IN) :: pzero,pone
    real(kind=qp),dimension(0:MICROMAX,0:MICROMAX),INTENT(OUT) :: noise_mtx_micro

    call generate_transcriptional_noise_operator(nMicro,pzero,pone,noise_mtx_micro)
end subroutine generate_micro_noise


subroutine make_hyptensor(nMe,hyptensor)
    integer,intent(in) :: nMe
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX), intent(out) :: hyptensor
    integer :: i,j,k
    !print*,nMe,m_d_lo,m_d_hi
    hyptensor = 0.0_qp
    !print*,"makehyp"
    do i=0,nMe+1
        do j=0,nMe+1
            do k=0,nMe+1
                hyptensor(i,j,k) = hyp2(nMe,i,j,k)
            enddo
        enddo
    enddo
    !print*,"size of hyptensor element = ", sizeof(hyptensor(0,0,0)), " size of hyptensor = ",sizeof(hyptensor)
    !print*,"maxval of hyptensor = ",maxval(hyptensor,mask=hyptensor .lt. 0.99_qp)
    !print*,"minval of hyptensor = ",minval(hyptensor,mask=hyptensor .gt. 0.01_qp)
end subroutine make_hyptensor

subroutine check_hyptensor(nMe,htorig)
    integer,intent(in) :: nMe
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX), intent(in) :: htorig
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX) :: hyptensor
    integer :: i,j,k
    real(kind=qp) :: maxdiff
    print*,"Checking tensor"
    hyptensor = 0.0_qp
    do i=0,nMe+1
        do j=0,nMe+1
            do k=0,nMe+1
                hyptensor(i,j,k) = hyp2(nMe,i,j,k)
            enddo
        enddo
    enddo
    maxdiff = maxval(abs(hyptensor-htorig))
    
    if (maxdiff .gt. epsilon(1.0_sp)) then
        print*,"maxdiff exceeds tolerance for hyptensor = ",maxdiff
        stop
    endif
!print*,"size of hyptensor element = ", sizeof(hyptensor(0,0,0)), " size of hyptensor = ",sizeof(hyptensor)
end subroutine check_hyptensor


subroutine make_mirtarg_probs(mdl,mdh,effprob,nmirtarg_probs)
    integer,intent(in) :: mdl,mdh
    real(kind=qp),intent(in) :: effprob
    real(kind=qp),dimension(0:MESSMAX),intent(out) :: nmirtarg_probs
    integer :: miRNA_deg_range
    nmirtarg_probs = 0.0_qp
    miRNA_deg_range = mdh - mdl

    call binomialDistribution(miRNA_deg_range,effprob,nmirtarg_probs(mdl:mdh))
end subroutine make_mirtarg_probs

subroutine check_mirtarg_probs(mdl,mdh,effprob,nmirtarg_orig)
    integer,intent(in) :: mdl,mdh
    real(kind=qp),intent(in) :: effprob
    real(kind=qp),dimension(0:MESSMAX),intent(in) :: nmirtarg_orig
    real(kind=qp),dimension(0:MESSMAX) :: nmirtarg_check
    integer :: miRNA_deg_range
    nmirtarg_check = 0.0_qp
    miRNA_deg_range = mdh - mdl

    call binomialDistribution(miRNA_deg_range,effprob,nmirtarg_check(mdl:mdh))
    nmirtarg_check = abs(nmirtarg_check - nmirtarg_orig)
    if (maxval(nmirtarg_check) .gt. epsilon(1.0_sp)) then
        print*,"mirtarg_probs, maxdiff exceeds tol = ",maxval(nmirtarg_check)
        stop
    endif


end subroutine check_mirtarg_probs



subroutine make_pkill_tensor(mdl,mdh,mirstrength,pkillten)
    INTEGER,INTENT(IN) :: mdl,mdh
    REAL(kind=qp),INTENT(IN) :: mirstrength
    real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),INTENT(OUT) :: pkillten
    integer :: i,j
    real(kind=qp),DIMENSION(0:MESSMAX) :: elimProbs

    pkillten=0.0_qp
    !print*,"pkilltensorparams: mdl = ",mdl," mdh= ",mdh," mirstrength= ",mirstrength

    elimProbs = 0.0_qp
    do i=mdl,mdh
        elimProbs(i) = mirstrength / i
    end do

    do i=0,MESSMAX
        do j=0,MESSMAX
            ! print*,i,j
            ! print*,elimProbs(j)
            ! print*,pkillten(i,0:1,j)
            call binomialDistribution(i,elimProbs(j),pkillten(i,0:i,j))
            ! if ((i .eq. 0) .and. (j .eq. 0)) then
            !     print*,"pkillzerovals = "
            !     print*,pkillten(i,:,j)
            ! endif
            if (abs(sum(pkillten(i,0:i,j)) -1.0_qp) .gt. epsilon(1.0_sp)) then
    !            print*,"i=",i,"j=",j,"elimProbs(j)=",elimProbs(j)
    !            print*,pkillten(i,:,j)
                stop
            endif
        end do
    end do
!    print*,"maxval of pkill = ",maxval(pkillten,mask=pkillten .lt. 0.99_qp)
!    print*,"minval of pkill = ",minval(pkillten,mask=pkillten .gt. 0.01_qp)
    !print*,pkillten(20,:,10)
end subroutine make_pkill_tensor

subroutine check_pkill_tensor(mdl,mdh,mirstrength,pkorig)
    integer,intent(in) ::mdl,mdh
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX), intent(in) :: pkorig
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX) :: cmparison
    integer :: i,j
    real(kind=qp),DIMENSION(0:MESSMAX) :: elimProbs
    real(kind=qp),INTENT(IN) :: mirstrength
    cmparison=0.0_qp
   
    elimProbs = 0.0_qp
    do i=mdl,mdh
        elimProbs(i) = mirstrength / i
    end do
 !   print*,"generating comparison kill tensor...."
    do i=0,MESSMAX
        do j=0,MESSMAX
            call binomialDistribution(i,elimProbs(j),cmparison(i,0:i,j))
        end do
    end do
    cmparison = abs(cmparison - pkorig)
    if (maxval(cmparison) > epsilon(1.0_sp)) then
        print*,"pkill comparison maxdiff exceeds tol = ",maxval(cmparison)
        print*,"maxval of pkorig = ",maxval(pkorig,mask=pkorig .lt. 0.99_qp)
        print*,"minval of pkorig = ",minval(pkorig,mask=pkorig .gt. 0.01_qp)
        stop
    endif
end subroutine check_pkill_tensor



end module combinatorialsubs