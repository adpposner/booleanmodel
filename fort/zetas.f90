module zetas
!This somewhat cryptic name refers to the variables zeta_ij in the manuscript. These correspond
!to the probabilities that a given miRNA will knock down its target in a given step. We have variables
!miRNAeffectprob, miRNA_deg_low, and miRNA_deg_high
    use base_parameters
    use modelparms
    use combinatorialsubs
    implicit none

    !used for testing purposes only
    real(kind=qp),DIMENSION(0:MESSMAX) :: elp
save
contains


    subroutine zetas_init(mdl,mdh,mstrength)
        integer,intent(in) :: mdl,mdh
        real(kind=qp),intent(in) :: mstrength
        integer :: i
        elp = 0.0_qp
        do i=mdl,mdh
            elp(i) = mstrength / i
        end do
    end subroutine zetas_init


    !here for completeness - prefer to always keep an end with an init
    subroutine zetas_end()
    
    end subroutine zetas_end

    !If N mRNAs are active, and miRNA targets k elements, what is probability of targeting j active elements
    !This is also pretty slow and can be sped up a lot by initializing table of hypergeometrics
    function PKtargetsmatchJactive(NactiveMess,Nmirtargs,jmatches,hyps)
        integer,intent(in) :: NactiveMess,Nmirtargs,jmatches
        real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: hyps
        real(kind=qp) :: PKtargetsmatchJactive
        real(kind=qp) :: PKhyptest
        PKhyptest = hyps(NactiveMess,Nmirtargs,jmatches)
        !PKtargetsmatchJactive = hyp2(Nme,NactiveMess,Nmirtargs,jmatches)
        ! if (abs(PKhyptest - PKtargetsmatchJactive) > epsilon(1.0_sp)) then
        !     print*, "nonequal, hyptensor = ",PKhyptest," func = ",PKtargetsmatchJactive
        !     print*, "nactive = ", NactiveMess," nmirtarg = ", Nmirtargs, " jmatches = ",jmatches
        !     stop
        ! end if
        PKtargetsmatchJactive=PKhyptest
    end function PKtargetsmatchJactive

    function pkillLGivenNTargsAndMatches(Ntargs,jmatches,lToKill,pkills)
        integer,intent(in) :: Ntargs,jmatches,lToKill

        real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: pkills
        real(kind=qp) :: pkillLGivenNTargsAndMatches
        

        pkillLGivenNTargsAndMatches = pkills(jmatches,lToKill,Ntargs)
        !old method - much slower
        !diff = binomialProbability(jmatches,lToKill,elp(Ntargs)) - pkillLGivenNTargsAndMatches
        
        
        !Check
        ! if (abs(diff) .gt. epsilon(1.0_sp)) then
        !     print*,"diff too big",diff,pkills(jmatches,lToKill,Ntargs),binomialProbability(jmatches,lToKill,elp(Ntargs))
        !     print*,jmatches,lToKill,Ntargs,elp(Ntargs)
        !     print*,"tensor value = "
        !     !print*,pkills(jmatches,:,Ntargs)
        !     print*,pkills(0,0,1)
        !     !print*,probofkillLGivenNTargsAndMatches(jmatches,:,Ntargs)
        !     stop
        ! end if
    
    end function pkillLGivenNTargsAndMatches

    function PKillLGivenNTargsAndNumActivemRNA(NactiveMess,Nmirtargs,lToKill,hyps,pkills)
        integer,intent(in) :: NactiveMess,Nmirtargs,lToKill
        real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: pkills,hyps
        real(kind=qp) :: PKillLGivenNTargsAndNumActivemRNA
        integer :: nmatches,maxpossiblematches

        maxpossiblematches = MIN(NactiveMess,Nmirtargs)
        PKillLGivenNTargsAndNumActivemRNA = 0.0_qp
        if (lToKill .gt. maxpossiblematches) then
            return
        endif
            do nmatches=lToKill,maxpossiblematches
                PKillLGivenNTargsAndNumActivemRNA = PKillLGivenNTargsAndNumActivemRNA + &
                PKtargetsmatchJactive(NactiveMess,Nmirtargs,nmatches,hyps)* &
                pkillLGivenNTargsAndMatches(Nmirtargs,nmatches,lToKill,pkills)
            enddo
        if (PKillLGivenNTargsAndNumActivemRNA .gt. 1.0d0) then
            PKillLGivenNTargsAndNumActivemRNA = 1.0d0
        endif
    end function PKillLGivenNTargsAndNumActivemRNA

    function PKillLGivenActiveMessAndParms(m_d_low,m_d_high,nmirtarg_probs,NactiveMess,lToKill,hyps,pkills)
        integer,intent(in) :: NactiveMess,lToKill,m_d_low,m_d_high
        real(kind=qp),DIMENSION(0:MESSMAX),intent(in) :: nmirtarg_probs
        real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: pkills,hyps
        real(kind=qp) :: PKillLGivenActiveMessAndParms
        integer :: ntarg,kmin,kmax


        kmin = MIN(lToKill,m_d_low)
        kmax = MAX(NactiveMess,m_d_high)
        PKillLGivenActiveMessAndParms = 0.0_qp
        do ntarg=kmin,kmax
            PKillLGivenActiveMessAndParms = PKillLGivenActiveMessAndParms + &
            PKillLGivenNTargsAndNumActivemRNA(NactiveMess,ntarg,lToKill,hyps,pkills)*nmirtarg_probs(ntarg)
        enddo
        !if (PKillLGivenActiveMessAndParms .gt. 1.0d0) then
        !print*, "PKillMessParms",nmirtarg_probs,NactiveMess,lToKill,PKillLGivenActiveMessAndParms
        !stop "PKillgivenparms > 1.0"
        !endif
    end function PKillLGivenActiveMessAndParms


    !This is the major function which determines the M operator
    function transprob_zeta(m_d_low,m_d_high,nmirtarg_probs,s1,s2,hyps,pkills)
        integer,intent(in) :: m_d_low,m_d_high
        real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: pkills,hyps
        real(kind=qp),DIMENSION(0:MESSMAX),intent(in) :: nmirtarg_probs
        integer,intent(in) :: s1,s2
        real(kind=qp) :: transprob_zeta
        integer :: lToKill
        lToKill = s1-s2
        transprob_zeta=0.0_qp
        if (lToKill .lt. 0) then
            return
        endif
       

        transprob_zeta = PKillLGivenActiveMessAndParms(m_d_low,m_d_high,nmirtarg_probs,s1,lToKill,hyps,pkills)

    end function transprob_zeta
 
end module zetas
