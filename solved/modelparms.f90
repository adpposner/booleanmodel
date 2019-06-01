module modelparms



integer,parameter :: MESSMAX = 150 
integer,parameter :: MICROMAX = 50
integer,parameter :: sp = kind(1.0), &
                     dp = selected_real_kind(2*precision(1.0_sp)), &
                     qp = selected_real_kind(2*precision(1.0_dp))

! integer :: error
! real(kind=qp) :: pmrProb

! integer :: integer_temp

! integer :: deg_low , deg_high, ncycles, n_sim, n_samples, errcode
type parmset
integer :: nMess,nMicro,deg_low,deg_high
 real(kind=qp) :: miRStrength

real(kind=qp) :: zero_one_ratio

real(kind=qp) :: pmrProb


real(kind=qp) :: miRNAeffectprob
integer :: miRNA_deg_high,miRNA_deg_low


real(kind=qp) :: RNAStrength


! !transcriptional noise
 logical :: transcription_noise_present
 real(kind=qp) :: transcription_p0
 real(kind=qp) :: transcription_p1
real(kind=qp) :: p_nz
end type parmset


contains

subroutine dbgprint(val)
  integer,intent(in) :: val
  print*,"Debug item ",val
end subroutine dbgprint

function allcloseto(val,mat)
  logical :: allcloseto
  real(kind=qp),INTENT(IN) :: val
  real(kind=qp),DIMENSION(:,:),INTENT(IN) :: mat
  
  if (maxval(abs(mat - val))>0.1_qp) then
    allcloseto = .false.
  ! print*,maxval(abs(mat - val))
   !print*,size(mat,1),size(mat,2)
   ! print*,maxloc(abs(mat - val))
   ! print*,mat(maxloc(mat,1),maxloc(mat,2))
  else
    allcloseto = .true.
  endif 
  end function allcloseto


  !Check if matrices are stochastic and fit according to nMess/nMicro
subroutine dbgmatsum(nzrow,nzcol,dimrow,dimcol,mat)
  integer,INTENT(IN) :: nzrow,nzcol,dimrow,dimcol
  real(kind=qp),dimension(0:,0:),INTENT(IN) :: mat
  real(kind=qp),dimension(size(mat,2)) :: sumcol,testcols
  real(kind=qp) :: maxdiff
  logical :: satisfies
  INTEGER :: i
  satisfies = .true.
  if (size(mat,2) .ne. dimcol) then
    print*,"mat invalid dimension 2 = ",size(mat,2), " should be = ", dimcol
    satisfies = .false.
  else if (size(mat,1) .ne. dimrow) then
    print*,"mat invalid dimension 1 = ",ubound(mat,1), " should be = ",dimrow
    satisfies = .false.
  else if (allcloseto(0.0_qp,mat(:,nzcol+1:ubound(mat,2))) .eqv. .false.) then
  !check if there are any nonzer
    print*,"mat has nonzero elements in columns from ",nzcol
    ! do i=nzcol+1,MESSMAX
    ! print*, i,pack(mat(:,i),mat(:,i) .eq. 0.0_qp)
    ! end do
    satisfies = .false.
  else if (allcloseto(0.0_qp,mat(nzrow+1:ubound(mat,1),:)) .eqv. .false.) then
    print*,"mat has nonzero elements in rows from ", nzrow
    !print*,mat(nzrow+1:size(mat,1),:)
    satisfies = .false.
  else
    sumcol = sum(mat,1)
    testcols(1:nzcol+1) = 1.0_qp
    testcols(nzcol+2:size(testcols)) = 0.0_qp
    maxdiff = maxval(abs(sumcol -testcols))
    print*,maxdiff
    if ( maxdiff > epsilon(1.0_sp)) then
      print*,"sumcols are invalid "
      !print("(151E36.25)"),sumcol
      !print*,"testcols ", testcols, size(mat,2)
      satisfies = .false.
    endif
  endif
  if  (satisfies .eqv. .false.) then
    print*,"unsatisfied"
    stop
  endif

end subroutine dbgmatsum


subroutine modelparameters_init_from_external(nMe, nMi, d_low , d_high,def,tmpr,miReffP,rnastr,pnz, &
mir_deg_low,mir_deg_high,Cnoisepresent,cp0,cp1,pset)
integer,intent(in)::nMe, nMi, d_low , d_high
real(kind=qp),intent(in) :: def,tmpr,miReffP,rnastr,pnz
integer,intent(in)::mir_deg_low,mir_deg_high
logical,intent(in) :: Cnoisepresent
real(kind=qp),intent(in) :: cp0,cp1

type(parmset),INTENT(OUT) :: pset
pset%nMess = nMe
pset%nMicro = nMi
pset%deg_low = d_low
pset%deg_high = d_high
pset%miRStrength = 1.0_qp - def

pset%miRNAeffectprob = miReffP
pset%miRNA_deg_high =mir_deg_high 
pset%miRNA_deg_low = mir_deg_low
pset%RNAStrength = rnastr

pset%transcription_noise_present = Cnoisepresent
pset%transcription_p0 = cp0
pset%transcription_p1 = cp1

pset%p_nz = pnz
pset%pmrProb = tmpr /(1.0d0+tmpr)

end subroutine modelparameters_init_from_external

subroutine modelparameters_end()
end subroutine modelparameters_end


subroutine write_mtx(mat,fname)
character(len=*) :: fname
real(kind=qp),dimension(0:,0:) :: mat
CHARACTER(len=35) :: fmtstring
write(fmtString,'( "(",I4,"(E36.25))")') size(mat,1)
open(unit=107,FILE=fname)
write(107,fmtstring) ((mat(i,j),i=0,ubound(mat,1)),j=0,ubound(mat,2))
close(107)
end subroutine write_mtx


end module modelparms
