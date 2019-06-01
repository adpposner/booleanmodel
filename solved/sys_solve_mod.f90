module sys_solve_mod
    use blas_interface
	use combinatorialsubs
	use zetas
	use operatorsolutions
    implicit none


    contains

    subroutine condentropy(nMess,mat,res)
    integer,INTENT(IN) :: nMess
    real(kind=qp),dimension(:,:),INTENT(IN) :: mat
    real(kind=qp),dimension(:), INTENT(INOUT) :: res
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX) :: tmp
    
    tmp = 0.0d0

    where (mat .gt. EPSILON(0.0d0))
        tmp = mat*log(mat)
    elsewhere
        tmp = 0.0d0
    endwhere
    res=0.0d0
    res = sum(tmp,1)
    res = -1.0d0 * res
    !print*,"res=",res
    end subroutine condentropy

    subroutine sys_solve(nMe, nMi, d_low , d_high,def,tmpr,miReffP,rnastr,pnz, &
    mir_deg_low,mir_deg_high,Cnoisepresent,cp0,cp1,sys_sln,condent)
    integer::nMe, nMi, d_low , d_high
    real(kind=qp) :: def,tmpr,miReffP,rnastr,pnz
    integer::mir_deg_low,mir_deg_high
    logical :: Lnoisepresent,Cnoisepresent
    real(kind=qp) :: cp0,cp1
    real(kind=qp) :: lp0,lp1
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX), intent(inout) :: sys_sln
    real(kind=qp),dimension(0:MESSMAX),INTENT(INOUT) :: condent
    type(parmset) :: pset
    sys_sln = 0.0d0
    call modelparameters_init_from_external(nMe,nMi,d_low,d_high,def,tmpr,miReffP,rnastr, &
    pnz,mir_deg_low,mir_deg_high,Cnoisepresent,cp0,cp1,pset)

    call solve_system(.FALSE.,pset,sys_sln)
    !print*,sum(sys_sln,1)
   ! call stationary_distn(nMe+1,sys_sln(0:nMe,0:nMe),condent)
   call condentropy(nMe,sys_sln,condent)
    sys_sln = TRANSPOSE(sys_sln)
    end subroutine sys_solve

    subroutine sys_solve_stationary_eig2(nMe, nMi, d_low , d_high,def,tmpr,miReffP,rnastr,pnz, &
    mir_deg_low,mir_deg_high,Cnoisepresent,cp0,cp1,eigv12,statdist12,condent)
    integer::nMe, nMi, d_low , d_high
    real(kind=qp) :: def,tmpr,miReffP,rnastr,pnz
    integer::mir_deg_low,mir_deg_high
    logical :: Cnoisepresent
    real(kind=qp) :: cp0,cp1

    real(kind=qp),dimension(0:MESSMAX,2), intent(inout) :: statdist12
    real(kind=qp),dimension(2),intent(inout) :: eigv12
    real(kind=qp),dimension(0:MESSMAX),INTENT(INOUT) :: condent
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX) :: sys_sln
    type(parmset) :: pset
    sys_sln = 0.0d0
    call modelparameters_init_from_external(nMe,nMi,d_low,d_high,def,tmpr,miReffP,rnastr, &
    pnz,mir_deg_low,mir_deg_high,Cnoisepresent,cp0,cp1,pset)

    call solve_system(.FALSE.,pset,sys_sln)
   !print*,"SYSSLN:",sum(sys_sln,1)
   call stationary_distn_2ev(nMe+1,sys_sln,eigv12,statdist12(0:nMe+1,:))
   call condentropy(nMe,sys_sln,condent)
    end subroutine sys_solve_stationary_eig2

    subroutine sys_solve_stationary_no_eig(nMe, nMi, d_low , d_high,def,tmpr,miReffP,rnastr,pnz, &
    mir_deg_low,mir_deg_high,Cnoisepresent,cp0,cp1,statdist,condent)
    integer,INTENT(IN)::nMe, nMi, d_low , d_high
    real(kind=qp),INTENT(IN) :: def,tmpr,miReffP,rnastr,pnz
    integer,INTENT(IN)::mir_deg_low,mir_deg_high
    logical,INTENT(IN) :: Cnoisepresent
    real(kind=qp),INTENT(IN) :: cp0,cp1
    real(kind=qp),dimension(0:MESSMAX),INTENT(INOUT) :: statdist,condent
    real(kind=qp),dimension(0:MESSMAX,2) :: stddisttmp
    real(kind=qp),dimension(2) :: evaltmp

    call sys_solve_stationary_eig2(nMe, nMi, d_low , d_high,def,tmpr,miReffP,rnastr,pnz, &
    mir_deg_low,mir_deg_high,Cnoisepresent,cp0,cp1,evaltmp,stddisttmp,condent)    
    if (Cnoisepresent .eqv. .TRUE.) then
    statdist = stddisttmp(:,1)
    else
    statdist = stddisttmp(:,2)
    endif

    end subroutine sys_solve_stationary_no_eig


end module sys_solve_mod