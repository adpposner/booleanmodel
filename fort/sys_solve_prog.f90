program sys_solve_prog
use blas_interface
use combinatorialsubs
use zetas
use operatorsolutions
use modelparms
use sys_solve_mod

implicit none
REAL(kind=qp), DIMENSION(0:MESSMAX) :: nmirtargsten
 REAL(kind=qp), DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX) :: pkillten
 REAL(kind=qp), DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX) :: hypsten
 REAL(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX) :: noise_mess
 REAL(kind=qp),DIMENSION(0:MICROMAX,0:MICROMAX) :: noise_micro

type(modelparams) :: pset
real(kind=qp),dimension(0:MESSMAX,0:MESSMAX):: sys_sln
integer :: nMess,nMicro,tf_d_low,tf_d_high,m_d_lo,m_d_hi
real(kind=qp) :: tf_p_nz,m_p_nz,rho,defect,tmpr,noise_pzero,noise_pone
    nMess = 130
    nMicro = 20
    tf_d_low = 1
    tf_d_high = 14
    m_d_lo = 3
    m_d_hi = 26
    tf_p_nz = 0.201_qp
    m_p_nz = 0.666_qp
    rho = 0.85_qp
    defect = 0.3_qp
    tmpr = 3.0_qp
    noise_pzero = 0.01_qp
    noise_pone = 0.01_qp
    call zetas_init(m_d_lo,m_d_hi,1.0_qp - defect)
call modelparameters_init_from_external(nMess, nMicro, tf_d_low, tf_d_high, tf_p_nz, defect, &
         tmpr, m_d_lo, m_d_hi, m_p_nz, rho, noise_pzero, noise_pone, pset)
!  call modelparameters_init_from_external(135,10,1,14,0.3_qp,1.5_qp,0.666_qp,0.75_qp, &
!     0.201_qp,3,26,.TRUE.,0.01_qp,0.01_qp,pset)
print*,"init"
call globals_init(pset, nmirtargsten, pkillten, hypsten,noise_mess,noise_micro)
!call check_pkill_tensor(m_d_lo,m_d_hi,1.0_qp - defect,pkillten)
!print*,"of 0,0,1",pkillten(0,0,1)
call solve_system(.TRUE.,pset,sys_sln,nmirtargsten,hypsten,pkillten,noise_mess,noise_micro)
!call write_mtx(sys_sln,"sys_sln.mtx.txt")

end program sys_solve_prog