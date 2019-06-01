module zetas
!This somewhat cryptic name refers to the variables zeta_ij in the manuscript. These correspond
!to the probabilities that a given miRNA will knock down its target in a given step. We have variables
!miRNAeffectprob, miRNA_deg_low, and miRNA_deg_high
    use modelparms
    
    use combinatorialsubs
    implicit none


    !This is a 2D array corresponding to the distribution of zeta
    !This ia binomially dist'd r.v. where zeta= (1-defect)/k with 
    !Probability \binom(miRNA_deg_high-miRNA_deg_low,k-miRNA_deg_low)*(miRNAeffectprob)^(k-miRNA_deg_low)*(1-miRNAeffectprob)^(miRNA_deg_high-k)
    !and this will use exact arithmetic
    real(kind=qp),dimension(0:MESSMAX) :: elimProbs
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX) :: probofkillLGivenNTargsAndMatches
save
contains

    subroutine zetas_init(pset,nmirtarg_probs)
        type(parmset),intent(in) :: pset
        real(kind=qp),DIMENSION(0:MESSMAX),intent(out) :: nmirtarg_probs
        real(kind=qp) :: eff_prob,mir_strength
		integer :: miRNA_deg_range,d_lo,d_hi,i,j

        
        d_lo = pset%miRNA_deg_low
        d_hi = pset%miRNA_deg_low
        eff_prob = pset%miRNAeffectprob
        mir_strength = pset%miRStrength
        nmirtarg_probs = 0.0d0
    
        !Important to recognize that the actual range is [miRNA_deg_low,miRNA_deg_high] inclusive
        miRNA_deg_range = d_hi - d_lo

        elimProbs = 0.0_qp
        do i=d_lo,d_hi
            elimProbs(i) = mir_strength / i
        end do
      
        !Can shorten subroutine pkillLGivenNTargsAndMatches
        probofkillLGivenNTargsAndMatches = 0.0_qp
        do i=0,MESSMAX
            do j=0,MESSMAX
                call binomialDistribution(i,elimProbs(j),probofkillLGivenNTargsAndMatches(i,0:i,j))
            end do
        end do
              

        !Fill probs for # targets
        call binomialDistribution(miRNA_deg_range,eff_prob,nmirtarg_probs(d_lo:d_hi))
        
    end subroutine zetas_init

    !here for completeness - prefer to always keep an end with an init
    subroutine zetas_end()
    
    end subroutine zetas_end

    !If N mRNAs are active, and miRNA targets k elements, what is probability of targeting j active elements
    !This is also pretty slow and can be sped up a lot by initializing table of hypergeometrics
    function PKtargetsmatchJactive(Nme,NactiveMess,Nmirtargs,jmatches)
        integer,intent(in) :: Nme,NactiveMess,Nmirtargs,jmatches
        real(kind=qp) :: PKtargetsmatchJactive
        PKtargetsmatchJactive = hyp2(Nme,NactiveMess,Nmirtargs,jmatches)
    end function PKtargetsmatchJactive

    function pkillLGivenNTargsAndMatches(miRStrength,Ntargs,jmatches,lToKill)
        integer,intent(in) :: Ntargs,jmatches,lToKill
        real(kind=qp),intent(in) :: miRStrength
        real(kind=qp) :: pkillLGivenNTargsAndMatches,diff

        !new method, just look up
        pkillLGivenNTargsAndMatches = probofkillLGivenNTargsAndMatches(jmatches,lToKill,Ntargs)
        !old method - much slower
        !pkillLGivenNTargsAndMatches = binomialProbability(jmatches,lToKill,elimProbs(Ntargs))
        
        !Check
        ! if (abs(diff) .gt. epsilon(1.0_sp)) then
        !     print*,"diff too big",diff,probofkillLGivenNTargsAndMatches(jmatches,lToKill,Ntargs),pkillLGivenNTargsAndMatches
        !     print*,jmatches,lToKill,Ntargs,elimProbs(Ntargs)
        !     !print*,probofkillLGivenNTargsAndMatches(jmatches,:,Ntargs)
        !     stop
        ! end if
    end function pkillLGivenNTargsAndMatches

    function PKillLGivenNTargsAndNumActivemRNA(nMe,miRStrength,NactiveMess,Nmirtargs,lToKill)
        integer,intent(in) :: nMe,NactiveMess,Nmirtargs,lToKill
        real(kind=qp),intent(in) :: miRStrength
        real(kind=qp) :: PKillLGivenNTargsAndNumActivemRNA
        integer :: nmatches,maxpossiblematches

        maxpossiblematches = MIN(NactiveMess,Nmirtargs)
        PKillLGivenNTargsAndNumActivemRNA = 0.0d0
        if (lToKill .gt. maxpossiblematches) then
            return
        endif
            do nmatches=lToKill,maxpossiblematches
                PKillLGivenNTargsAndNumActivemRNA = PKillLGivenNTargsAndNumActivemRNA + &
                PKtargetsmatchJactive(nMe,NactiveMess,Nmirtargs,nmatches)* &
                pkillLGivenNTargsAndMatches(miRStrength,Nmirtargs,nmatches,lToKill)
            enddo
        if (PKillLGivenNTargsAndNumActivemRNA .gt. 1.0d0) then
            PKillLGivenNTargsAndNumActivemRNA = 1.0d0
        endif
    end function PKillLGivenNTargsAndNumActivemRNA

    function PKillLGivenActiveMessAndParms(pset,nmirtarg_probs,NactiveMess,lToKill)
        type(parmset),intent(in) :: pset
        integer,intent(in) :: NactiveMess,lToKill
        real(kind=qp),DIMENSION(0:MESSMAX),intent(in) :: nmirtarg_probs
        real(kind=qp) :: PKillLGivenActiveMessAndParms
        integer :: ntarg,kmin,kmax


        kmin = MIN(lToKill,pset%miRNA_deg_low)
        kmax = MAX(NactiveMess,pset%miRNA_deg_high)
        PKillLGivenActiveMessAndParms = 0.0d0
        do ntarg=kmin,kmax
            PKillLGivenActiveMessAndParms = PKillLGivenActiveMessAndParms + &
            PKillLGivenNTargsAndNumActivemRNA(pset%nMess,pset%miRStrength,NactiveMess,ntarg,lToKill)*nmirtarg_probs(ntarg)
        enddo
        !if (PKillLGivenActiveMessAndParms .gt. 1.0d0) then
        !print*, "PKillMessParms",nmirtarg_probs,NactiveMess,lToKill,PKillLGivenActiveMessAndParms
        !stop "PKillgivenparms > 1.0"
        !endif
    end function PKillLGivenActiveMessAndParms


    !This is the major function which determines the M operator
    subroutine transprob_zeta(pset,nmirtarg_probs,s1,s2,p12)
        type(parmset),intent(in) :: pset
        real(kind=qp),DIMENSION(0:MESSMAX),intent(in) :: nmirtarg_probs
        integer,intent(in) :: s1,s2
        real(kind=qp),intent(out) :: p12
        integer :: ntarg,lToKill
        lToKill = s1-s2
        p12=0.0d0
        if (lToKill .lt. 0) then
            return
        endif


        p12 = PKillLGivenActiveMessAndParms(pset,nmirtarg_probs,s1,lToKill)

    end subroutine transprob_zeta
 
end module zetas
