module modelparms
use iso_c_binding
use base_parameters
use combinatorialsubs
implicit none



type,bind(c) :: modelparams
integer(c_int) :: nMess
integer(c_int) :: nMicro
integer(c_int) :: tf_deg_low
integer(c_int) :: tf_deg_high
real(C_QUAD128) :: tf_p_nz
real(C_QUAD128) :: defect
real(C_QUAD128) :: tmpr

integer(c_int) :: m_deg_high
integer(c_int)  :: m_deg_low
real(C_QUAD128) :: m_p_nz
real(C_QUAD128) :: rho

! !transcriptional noise
real(C_QUAD128) :: noise_pzero
real(C_QUAD128) :: noise_pone

end type modelparams


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


subroutine modelparameters_init_from_external(nMe, nMi, tf_d_low , tf_d_high,tf_p_nz,defect,tmpr,m_d_low,m_d_high, &
        m_p_nz,rho,noise_pzero,noise_pone,pset)
integer,intent(in) :: nMe, nMi, tf_d_low , tf_d_high,m_d_low,m_d_high
real(kind=qp),intent(in) :: tf_p_nz,defect,tmpr,m_p_nz,rho,noise_pzero,noise_pone
type(modelparams),INTENT(OUT) :: pset

pset%nMess = nMe
pset%nMicro = nMi
pset%tf_deg_low = tf_d_low
pset%tf_deg_high = tf_d_high
pset%defect = defect

pset%m_p_nz = m_p_nz
pset%m_deg_high =m_d_high 
pset%m_deg_low = m_d_low
pset%rho = rho

pset%noise_pzero = noise_pzero
pset%noise_pone = noise_pone

pset%tf_p_nz = tf_p_nz
pset%tmpr = tmpr

end subroutine modelparameters_init_from_external

subroutine modelparameters_end()
end subroutine modelparameters_end

subroutine print_modelparams(pset)
  type(modelparams),INTENT(IN) :: pset
  print*,"pset%nMess = ",pset%nMess
  print*,"pset%nMicro = ",pset%nMicro
  print*,"pset%tf_deg_low = ",pset%tf_deg_low
  print*,"pset%tf_deg_high = ",pset%tf_deg_high
  print*,"pset%defect = ",pset%defect
  print*,"pset%m_p_nz = ",pset%m_p_nz
  print*,"pset%m_deg_high =",pset%m_deg_high
  print*,"pset%m_deg_low = ",pset%m_deg_low
  print*,"pset%rho = ",pset%rho
  print*,"pset%noise_pzero = ",pset%noise_pzero
  print*,"pset%noise_pone = ",pset%noise_pone
  print*,"pset%tf_p_nz = ",pset%tf_p_nz
  print*,"pset%tmpr = ",pset%tmpr

end subroutine print_modelparams


subroutine write_mtx(mat,fname)
character(len=*) :: fname
real(kind=qp),dimension(0:,0:) :: mat
CHARACTER(len=35) :: fmtstring
integer :: i,j
write(fmtString,'( "(",I4,"(E48.25E3))")') size(mat,2)
open(unit=107,FILE=fname)
write(107,fmtstring) ((mat(i,j),j=0,ubound(mat,2)),i=0,ubound(mat,1))
close(107)
end subroutine write_mtx


subroutine globals_init(paramset,nmirtarg_probs,pkillLvalues,hypsvalues,messnoise,micronoise)
  type(modelparams),intent(in) :: paramset
  real(kind=qp),DIMENSION(0:MESSMAX),intent(out) :: nmirtarg_probs
  real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(out) :: pkillLvalues
  real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(out) :: hypsvalues
  real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX),INTENT(OUT) :: messnoise
  real(kind=qp),DIMENSION(0:MICROMAX,0:MICROMAX),INTENT(OUT) :: micronoise
  real(kind=qp) :: eff_prob,mir_strength,pzero,pone
  integer :: m_d_lo,m_d_hi,nMess,nMicro

  !real(kind=qp),dimension(0:MESSMAX) :: elimProbs

  nMess = paramset%nMess  
  m_d_lo = paramset%m_deg_low
  m_d_hi = paramset%m_deg_high
  eff_prob = paramset%m_p_nz
  mir_strength = 1.0_qp - paramset%defect
  pzero = paramset%noise_pzero
  pone = paramset%noise_pone
  nMicro = paramset%nMicro
  nmirtarg_probs = 0.0_qp
  hypsvalues = 0.0_qp
  pkillLvalues(0:MESSMAX,0:MESSMAX,0:MESSMAX) = 0.0_qp
  !print*,"making noise operators..."
    call generate_mess_noise(nMess,pzero,pone,messnoise)
    call generate_micro_noise(nMicro,pzero,pone,micronoise)
  !print*,"making noise operators complete"
  !print*,"making nmirtargprobs..."
  call make_mirtarg_probs(m_d_lo, m_d_hi, eff_prob, nmirtarg_probs)
  !print*,"mirtargprobs making complete"
  !print*,"making hyptensor..."
  call make_hyptensor(nMess,hypsvalues)
  !print*, "init hyptensor complete"
  !print*,"making pkill tensor..."
  call make_pkill_tensor(m_d_lo, m_d_hi, mir_strength, pkillLvalues)
  !print*, "making pkill complete"
  
end subroutine globals_init

end module modelparms
