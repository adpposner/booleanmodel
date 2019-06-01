module combinatorialsubs
  use modelparms
implicit none



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
function pNZN2(N,N2,dl,dh,pnz)
    integer,INTENT(IN) :: N,N2,dl,dh
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

function pPlusgivenN2(N,N1,N2,p)
    integer,INTENT(IN) :: N,N1,N2
    real(kind=qp),INTENT(IN) :: p
    integer :: k,kmax
    real(kind=qp) :: hypNN1N2k,ppgk
    real(kind=qp) :: pPlusgivenN2

    
    kmax=min(N1,N2)
    pPlusgivenN2 = 0.0_qp

    do k=1,kmax
        hypNN1N2k = hyp2(N,N1,N2,k)
        ppgk = pPlusgivenk(k,p)
        pPlusgivenN2 = pPlusgivenN2 + hypNN1N2k*ppgk
    end do

end function pPlusgivenN2


function pPlusgivendegs(N,N1,dl,dh,pnz,pmr)
    integer,INTENT(IN) :: N,N1,dl,dh
    real(kind=qp),INTENT(IN) :: pnz,pmr
    real(kind=qp) :: pNZN2val,pOnN2
    real(kind=qp) :: pPlusgivendegs
    integer  :: N2
    pPlusgivendegs = 0.0_qp

    do N2=dl,dh
        pNZN2val = pNZN2(N,N2,dl,dh,pnz)
        pOnN2 = pPlusgivenN2(N,N1,N2,pmr)
        pPlusgivendegs = pPlusgivendegs + pNZN2val*pOnN2
    end do
end function pPlusgivendegs




!This generates our probability of transition for a given number of nonzero inputs N1
subroutine pPlusForModel(N1,pset,pPlus)
    integer,intent(in) :: N1
    type(parmset),intent(in) :: pset
    real(kind=qp),intent(out) :: pPlus

    pPlus = pPlusgivendegs(pset%nMess,N1,pset%deg_low,pset%deg_high,pset%p_nz,pset%pmrProb)

end subroutine pPlusForModel


end module combinatorialsubs