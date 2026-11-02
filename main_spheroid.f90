
! ---------------------------------------------------------------|
! Solves the SCMFT eqs for nucleosomes chains and                |
! weak  polyelectrolytes polymers                                |
! input/output: see myio.f90                                     |
! ---------------------------------------------------------------|


program main

    !  .. variable and constant declaractions
    ! use mpivars
    use globals       ! parameters definitions
    use physconst
    use mathconst
    use volume
    use random
    use field
    use parameters
    use matrices
    use energy
    use conform_entropy
    use chains
    use chaingenerator
    use VdW
    use listfcn
    use fcnpointer
    use initxvector
    use surface
    use myio
    use myutils
    use dielectric_const
    use modfcnMg
    use modfcnMgexpl
    use modfcnFeexpl
    use chains_spheroid, only : init_chain_NCP
    use chaingenerator_spheroid, only : init_spheroid, make_chain_spheroid_mc, read_SREV, build_NCP_chain
    
    implicit none

    real(dp),  dimension(:), allocatable :: x         ! iteration vector
    real(dp),  dimension(:), allocatable :: xguess    ! guess iteration vector
    real(dp),  dimension(:), allocatable :: xstored   ! stored iteration vector
    real(dp),  dimension(:), allocatable :: fvec

    integer :: i,c,num
    logical :: use_xstored
    logical :: isfirstguess
    logical :: issolution
    integer :: info
    character(len=lenText) :: text, istr
    character(len=20) :: fname
    type (looplist), pointer :: loop
    real(dp) :: loopbegin,  loopstepsizebegin
    real(dp), parameter :: loopeps = 1.0e-10_dp 
    real(dp), parameter :: listeps = 1.0e-7_dp   
    real(dp), dimension(:),  pointer :: list
    real(dp), pointer :: list_val
    real(dp) :: list_step
    integer  :: nlist_elem, maxlist_elem, nlist_step
    integer :: phoscutoff
    logical :: isVolfracLargerOne

    ! .. executable statements

    ! .. logfile
    
    fname='status.log'
    call open_logfile(logUnit,fname)

    
    text='program begins'
    call print_to_log(LogUnit,text)
    write(istr,'(A40)')VERSION
    text='program version     = '//trim(adjustl(istr))
    call print_to_log(LogUnit,text)
    print*,text

    ! .. init
    call read_inputfile(info)
    call error_handler(info,"read_inputfile")


    call init_constants()
    call make_geometry()                        ! generate volume elements lattice
    call allocate_chain_parameters()            !  
    call init_matrices()                        ! init matrices for MC chain generation
    call allocate_chains(cuantas,nnucl,nseg,nsegAA,nsegtypes,nsegtypesAA,maxnchainsrotations,maxnchainsrotationsxy)
    call init_chain_parameters()                                   ! reads pKa, pKaions zpol and vpol from file
    call make_sequence_chain(chainperiod,chaintype)                ! set isAmonomer, type_of_monomer and  type_of_monomer_char
    call make_charge_table(ismonomer_chargeable,zpol,nsegtypes)    ! set ismonomer_chargeable,
    call make_type_of_charge_table(type_of_charge,zpol,nsegtypes)  ! set type_of_charge
    call make_segcom(segcm,nnucl,segcmfname)
    call set_properties_chain(chainperiod,chaintype) 
    call set_mapping_num_to_char(mapping_num_to_char)
   

    ! init distributed volume  
    if(systype=="nucl_ionbin_sv".or.systype=="nucl_neutral_sv".or.systype=="nucl_ionbin_Mg".or. &
        systype=="nucl_ionbin_MgA".or.systype=="nucl_ionbin_Fe") then 
        call init_vnucl_type(info) ! ismonomer_chargeable etc needs to be set
        call error_handler(info,"init_vnucl_type")
    endif 
 
    call make_VdWeps(info) 
    call error_handler(info,"make_VdWeps")
    ! call set_value_isVdW_on_values(nsegtypes, VdWeps, isVdW) 
    call set_value_isVdW(systype,isVdW)
    write(istr,'(L2)')isVdW
    text='VdW interaction: isVdW = '//trim(adjustl(istr))
    call print_to_log(LogUnit,text) 
    print*,text

    call init_chain_NCP(info)                                   ! needs  to be called before make_chains
    call error_handler(info,"init_chain_NCP")

    call make_chains(chainmethod,systype)   
    call allocate_field(nx,ny,nz,nsegtypes)
   
    if(systype=="nucl_ionbin_Mg") then ! auxiliary array index_phos
        call find_phosphate_location(index_phos,inverse_index_phos,len_index_phos) 
        call allocate_index_neighbors_phos(maxneigh,len_index_phos)
        phoscutoff=int(distphoscutoff/delta)+2
        call make_table_index_neighbors_phos(phoscutoff,len_index_phos,index_phos) 
        call allocate_field_pairs(nx,ny,nz,maxneigh,7,len_index_phos) ! internal systype switch !
    endif

    if(systype=="nucl_ionbin_MgA") then 
        phoscutoff=int(distphoscutoff/delta)+2
        call allocate_field_pairs(nx,ny,nz,maxneigh,7,len_index_phos) ! internal systype switch !
    endif   
    
    if(systype=="nucl_ionbin_Fe") then 
        phoscutoff=int(distphoscutoff/delta)+2
        call allocate_field_pairs(nx,ny,nz,maxneigh,7,len_index_phos) ! internal systype switch !
        call allocate_field_triplets(11)
        call init_var_compute_fdisPPP
    endif   
    
    if(systype=="nonuclcp") then 
        no_overlapchain(1)=.true. ! otherwise loop does not start !!
    endif    
    
    call init_field()
    call init_surface(bctype,nsurf)
    call make_isrhoselfconsistent(info)
    call set_size_neq()             ! number of non-linear equation neq
    call set_fcn()
    call set_dielect_fcn(dielect_env)
    call write_chain_config()
    call write_chain_struct(write_struct,info) 
   ! if(write_struct) call error_handler(1,"stop after write_chains_struct")

    ! call test_index_histone(info)  

    ! .. computation starts


    print*,"test-spheriod" 
    call init_spheroid(info)
    call error_handler(info,"init_spheriod")
   
    ! call make_chain_spheroid_mc()

    call build_NCP_chain(info)

  
    ! call write_chain_config()
   ! call write_chain_struct(write_struct,info) 
   ! if(write_struct) call error_handler(1,"stop after write_chains_struct")


    text="program end"

    call print_to_log(LogUnit,text)
    call close_logfile(LogUnit)

end program main
