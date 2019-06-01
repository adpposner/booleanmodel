program sys_solve_prog
use blas_interface
use combinatorialsubs
use zetas
use operatorsolutions
use modelparms
use sys_solve_mod

implicit none
real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX) :: bindist
integer :: i,k,j
character(len=35) :: fmtname
type(parmset) :: pset
real(kind=qp),dimension(0:MESSMAX,0:MESSMAX):: sys_sln
 call modelparameters_init_from_external(135,10,1,14,0.3_qp,1.5_qp,0.666_qp,0.75_qp, &
    0.201_qp,3,26,.TRUE.,0.01_qp,0.01_qp,pset)
print*,"init"
call solve_system(.TRUE.,pset,sys_sln)


end program sys_solve_prog