                                                                  
!     .. module file of global variables 

module globals

    use precision_definition
    use mathconst
  
    implicit none

    ! .. variables  
    
    integer  :: nsize         ! size lattice, numer of layers
  
    !  .. variables of input nucloesome/ncp conformation

    integer  :: nnucl         ! number of nucleosome
    integer  :: nseg          ! number of segment 
    integer  :: nsegAA        ! number of AA segments
    integer  :: nsegtypes     ! number of segment types
    integer  :: nsegtypesAA   ! number of AA segment types 

    integer  :: cuantas       ! number of configurations
    integer  :: cuantas_no_overlap ! number of configurations with no overlap
    integer  :: max_confor    ! maximum number of configurations
    integer  :: set_confor    ! number that indicates which set of configuration to use in 
    integer  :: local_conf    ! conformation number used from set

    !  .. variables of output nucloesome/ncp conformation

    integer :: nseg_sph      ! number of spheriod/nucleosme 
    integer :: nnucl_sph 
    real(dp) :: lseg_sph     
    integer :: max_confor_sph


    ! .. development variables
    integer :: nsegsource     ! number of segment in source traj files
    integer :: s_begin        ! first segment to actually use from source traj file
    integer :: s_end          ! last segement to actually use from source traj file 
 
    integer(8) :: neq         ! number of non-linear equations
    integer(8) :: neqmax      ! maximum number of non-linear equations
    integer  :: neqint        ! number of non-linear equations, for mpi fnc bindings
    
    character(len=15) :: systype   ! systype selects fcn    
    character(len=15) :: runtype   ! runtype
    character(len=2)  :: bctype(2) ! bctype selects electrostatic boundary condition of surface  LEFT and RIGHT

    integer, parameter :: LEFT = 1. ! z=0 bottom plane 
    integer, parameter :: RIGHT = 2 ! z=nz delta  top plane

    logical, parameter :: DEBUG = .false. ! switch for debug information

end module globals

