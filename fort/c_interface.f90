module c_interface
    use,intrinsic :: iso_c_binding
    use modelparms
    use combinatorialsubs
    use operatorsolutions
    use sys_solve_mod,only : condentropy
    use blas_interface,only : stationary_distn
    use zetas,only : zetas_init
implicit none

contains


subroutine check_parms_struct_values(pset) bind(c)
    type(modelparams) :: pset
    call print_modelparams(pset)
  
end subroutine check_parms_struct_values
  
  
subroutine get_sizes_bytes(mirtargprobsize,hyptensorsize,pkilltensorsize,messnoisesize,micronoisesize, &
    condentsize, statdistsize, totalressize) bind(c)
    INTEGER(c_int),INTENT(OUT) :: mirtargprobsize,hyptensorsize,pkilltensorsize,condentsize
    INTEGER(c_int),INTENT(OUT) :: statdistsize,totalressize
    integer(c_int),INTENT(OUT) :: messnoisesize,micronoisesize
    real(c_double),DIMENSION(0:MESSMAX,0:MESSMAX) :: dummyslnmatrix
    real(c_double),DIMENSION(0:MESSMAX) :: dummystatdistcondentvec
    real(kind=qp),dimension(0:MESSMAX,0:MESSMAX,0:MESSMAX) :: dummypkilltensor,dummyhyptensor
    real(kind=qp),DIMENSION(0:MESSMAX) :: dummymirtarg
    real(kind=qp),DIMENSION(0:MESSMAX,0:MESSMAX) :: dummymessnoise
    real(kind=qp),DIMENSION(0:MICROMAX,0:MICROMAX) :: dummymicronoise
    
    mirtargprobsize = sizeof(dummymirtarg)
    hyptensorsize=sizeof(dummyhyptensor)
    pkilltensorsize = sizeof(dummypkilltensor)
    messnoisesize = sizeof(dummymessnoise)
    micronoisesize = sizeof(dummymicronoise)
    statdistsize = sizeof(dummystatdistcondentvec)
    condentsize = sizeof(dummystatdistcondentvec)
    totalressize = sizeof(dummyslnmatrix)

end subroutine get_sizes_bytes

subroutine update_global_vars(params_old,params_new,nmirtargsData,pkillData, &
    hypsData, noisemess,noisemicro) bind(C)
    type(modelparams),INTENT(IN) :: params_old
    type(modelparams),INTENT(IN) :: params_new
    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),INTENT(INOUT) :: pkillData,hypsData
    real(C_QUAD128),DIMENSION(0:MESSMAX),INTENT(INOUT) :: nmirtargsData
    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX),INTENT(INOUT) :: noisemess
    real(C_QUAD128),DIMENSION(0:MICROMAX,0:MICROMAX),INTENT(INOUT) :: noisemicro
    logical :: new_mess_noise,new_micro_noise,new_mirtargs,new_hyps,new_pkill
    new_mess_noise = .false.
    new_micro_noise = .false.
    new_mirtargs = .false.
    new_hyps = .false.
    new_pkill = .false.
    if (params_old%nMess .ne. params_new%nMess) then
        new_mess_noise = .true.
        new_hyps = .true.
    endif

    if (params_old%nMicro .ne. params_new%nMicro) then
        new_micro_noise = .true.
    endif
    !THESE ARE REGENERATED EVERY TIME
    if (params_old%tf_deg_low .ne. params_new%tf_deg_low) then
    endif
    if (params_old%tf_deg_high .ne. params_new%tf_deg_high) then
    endif
    if (params_old%tf_p_nz .ne. params_new%tf_p_nz) then
    endif
    if (params_old%tmpr .ne. params_new%tmpr) then
    endif
    if (params_old%rho .ne. params_new%rho) then
    endif
    !END REGENERATED EVERY TIME
    if (params_old%defect .ne. params_new%defect) then
        new_pkill = .true.
    endif

    if (params_old%m_deg_high .ne. params_new%m_deg_high) then
        new_pkill = .true.
        new_mirtargs = .true.
        !not necessary despite params
        !new_hyps = .true.
    endif
    if (params_old%m_deg_low .ne. params_new%m_deg_low) then
        new_pkill = .true.
        new_mirtargs = .true.
        !not necessary despite params
        !new_hyps = .true.

    endif
    if (params_old%m_p_nz .ne. params_new%m_p_nz) then
        new_mirtargs = .true.
        !not necessary despite params
        !new_hyps = .true.
    endif

    if (params_old%noise_pzero .ne. params_new%noise_pzero) then
        new_mess_noise = .true.
        new_micro_noise = .true.
    endif
    if (params_old%noise_pone .ne. params_new%noise_pone) then
        new_mess_noise = .true.
        new_micro_noise = .true.
    endif

    if (new_mess_noise) then
    call generate_mess_noise(params_new%nMess,params_new%noise_pzero,params_new%noise_pone,noisemess)
    endif
    if (new_micro_noise) then
    call generate_micro_noise(params_new%nMicro,params_new%noise_pzero,params_new%noise_pone,noisemicro)
    endif
    if (new_mirtargs) then
    call make_mirtarg_probs(params_new%m_deg_low, params_new%m_deg_high,params_new%m_p_nz, &
        nmirtargsData)
    endif
    if (new_hyps) then
        call make_hyptensor(params_new%nMess,hypsData)
    endif
    if (new_pkill) then
        call make_pkill_tensor(params_new%m_deg_low, params_new%m_deg_high, &
            1.0_qp - params_new%defect, pkillData)
    endif
end subroutine update_global_vars

subroutine generate_global_vars(pset,nmirtargsData,pkillData,hypsData,noisemess,noisemicro) bind(c)
    type(modelparams),INTENT(IN) :: pset
    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(out) :: pkillData,hypsData
    real(C_QUAD128),DIMENSION(0:MESSMAX),INTENT(OUT) :: nmirtargsData
    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX),INTENT(OUT) :: noisemess
    real(C_QUAD128),DIMENSION(0:MICROMAX,0:MICROMAX),INTENT(OUT) :: noisemicro
    call globals_init(pset,nmirtargsData,pkillData,hypsData,noisemess,noisemicro)

end subroutine generate_global_vars

subroutine check_global_vars(pset,nmtorig,pkillorig,hypsorig) bind(c)
type(modelparams),INTENT(IN) :: pset
real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: pkillorig,hypsorig
real(C_QUAD128),DIMENSION(0:MESSMAX),INTENT(in) :: nmtorig
    !TODO_add_body

call check_mirtarg_probs(pset%m_deg_low,pset%m_deg_high,pset%m_p_nz,nmtorig)
!print*,"checking pkill values"
call check_pkill_tensor(pset%m_deg_low,pset%m_deg_high,1.0_qp - pset%defect,pkillorig)
!print*,"pkill values consistent"
!print*,"checking hyptensor values"
call check_hyptensor(pset%nMess,hypsorig)
!print*,"hyptensor values consistent"
end subroutine check_global_vars

subroutine solve_sys_base(write_mats,pset,nmirtarg_probs,pkills,hyps, &
    noise_mess,noise_micro,condent_out,meanexp_out,stddev_out,entrate_out) bind(c)
    type(modelparams),INTENT(IN) :: pset
    INTEGER(c_int),INTENT(IN) :: write_mats
    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: pkills,hyps
    real(C_QUAD128),DIMENSION(0:MESSMAX),INTENT(in) :: nmirtarg_probs
    real(C_QUAD128),DIMENSION(0:MESSMAX) :: condent_internal
    real(C_QUAD128) :: meanexp_internal,stddev_internal,entrate_internal
    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX),INTENT(IN) :: noise_mess
    real(C_QUAD128),DIMENSION(0:MICROMAX,0:MICROMAX),INTENT(IN) :: noise_micro

    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX):: total_res_internal
    real(C_QUAD128),DIMENSION(0:MESSMAX) :: statdist_internal    
    real(C_QUAD128) :: eigval_internal,twom_internal
    real(C_QUAD128),dimension(0:MESSMAX) :: mnvals
    real(C_QUAD128),dimension(0:MESSMAX) :: mnvals2 
    real(c_double),DIMENSION(0:MESSMAX),INTENT(OUT) :: condent_out
    real(c_double),INTENT(OUT) :: meanexp_out,stddev_out,entrate_out
 
    integer :: i
    LOGICAL :: write_matrices
    if (write_mats .eq. 0) then
        write_matrices = .false.
    else
        write_matrices = .true.
    endif
    call zetas_init(pset%m_deg_low,pset%m_deg_high,1.0_qp - pset%defect)
    mnvals =  (/ (i,i=0,MESSMAX) /)
    mnvals2 = mnvals * mnvals
    condent_internal = 0.0_qp
    meanexp_internal = 0.0_qp
    stddev_internal = 0.0_qp
    entrate_internal = 0.0_qp
    total_res_internal = 0.0_qp
    statdist_internal = 0.0_qp
    eigval_internal = 0.0_qp

    condent_out = 0.0_dp
    meanexp_out = 0.0_dp
    stddev_out = 0.0_dp
    entrate_out = 0.0_dp
    

    call solve_system(write_matrices,pset,total_res_internal,nmirtarg_probs, &
    hyps,pkills,noise_mess,noise_micro)
    call condentropy(total_res_internal,condent_internal)
    call stationary_distn(total_res_internal,eigval_internal,statdist_internal)
    !print*,"stationary_distn",statdist_internal
    !print*,"condent",condent_internal
    entrate_internal = dot_product(condent_internal,statdist_internal)
    meanexp_internal = dot_product(statdist_internal,mnvals)
    twom_internal = dot_product(statdist_internal,mnvals2)
    stddev_internal = sqrt(twom_internal - (meanexp_internal * meanexp_internal))
    !print*,"entrate = ",entrate_internal, " meanexp = ", meanexp_internal, &
     !   " twom = ",twom_internal, " stddev = ",stddev_internal
    condent_out = real(condent_internal,dp)
    meanexp_out = real(meanexp_internal,dp)
    stddev_out = real(stddev_internal,dp)
    entrate_out = real(entrate_internal,dp)
end subroutine solve_sys_base

subroutine solve_sys_cplte(write_mats,pset,nmirtarg_probs,pkills,hyps, &
    noise_mess,noise_micro, total_res_out, condent_out, &
    statdist_out,eigval_out,entrate_out) bind(c)
    integer(c_int),INTENT(IN) :: write_mats
    type(modelparams),INTENT(IN) :: pset
    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX),INTENT(IN) :: noise_mess
    real(C_QUAD128),DIMENSION(0:MICROMAX,0:MICROMAX),INTENT(IN) :: noise_micro
    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX,0:MESSMAX),intent(in) :: pkills,hyps
    real(C_QUAD128),DIMENSION(0:MESSMAX),INTENT(in) :: nmirtarg_probs
    real(c_double),DIMENSION(0:MESSMAX,0:MESSMAX),INTENT(OUT) :: total_res_out
    real(C_QUAD128),DIMENSION(0:MESSMAX,0:MESSMAX) :: total_res_internal
    real(kind=qp),DIMENSION(0:MESSMAX) :: condent_internal
    real(kind=qp),DIMENSION(0:MESSMAX) :: statdist_internal
    real(c_double),DIMENSION(0:MESSMAX),intent(out) :: condent_out
    real(c_double),DIMENSION(0:MESSMAX),intent(out) :: statdist_out
    real(kind=qp) :: eigval_internal,entrate_internal
    real(c_double),INTENT(OUT) :: eigval_out,entrate_out
   
    LOGICAL :: write_matrices

    condent_internal = 0.0_qp
    statdist_internal = 0.0_qp
    eigval_internal = 0.0_qp
    entrate_internal = 0.0_qp
    condent_out = 0.0_dp
    statdist_out = 0.0_dp
    eigval_out = 0.0_dp
    entrate_out = 0.0_dp
    call zetas_init(pset%m_deg_low,pset%m_deg_high,1.0_qp - pset%defect)
    if (write_mats .eq. 0) then
        write_matrices = .false.
    else
        write_matrices = .true.
    endif
    call solve_system(write_matrices,pset,total_res_internal,nmirtarg_probs, &
    hyps,pkills,noise_mess,noise_micro)
    call condentropy(total_res_internal,condent_internal)
    call stationary_distn(total_res_internal,eigval_internal,statdist_internal)
    entrate_internal = dot_product(condent_internal,statdist_internal)
    total_res_internal = TRANSPOSE(total_res_internal)
    condent_out = real(condent_internal,dp)
    entrate_out = real(entrate_internal,dp)
    total_res_out = real(total_res_internal,dp)
    statdist_out = real(statdist_internal,dp)
    eigval_out = real(eigval_internal,dp)
end subroutine solve_sys_cplte


end module c_interface
