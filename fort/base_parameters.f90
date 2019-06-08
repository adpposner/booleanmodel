module base_parameters
use,intrinsic :: iso_c_binding
    implicit none

integer,parameter :: MESSMAX = 150 
integer,parameter :: MICROMAX = 80
integer,parameter :: sp = kind(1.0), &
                     dp = selected_real_kind(2*precision(1.0_sp)), &
                     qp = selected_real_kind(2*precision(1.0_dp))

    
integer(kind=4),PARAMETER :: C_QUAD128 = 16



contains
end module base_parameters