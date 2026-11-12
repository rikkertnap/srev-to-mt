
! ---------------------------------------------------------------|
! Solves the SCMFT eqs for nucleosomes chains and                |
! weak  polyelectrolytes polymers                                |
! input/output: see myio.f90                                     |
! ---------------------------------------------------------------|


program main

    !  .. variable and constant declaractions
  
    use globals, only : cuantas, nnucl, nseg, nsegAA, nsegtypes, nsegtypesAA, systype, bctype
   ! use field, only : allocate_field_pairs, allocate_field, allocate_field_pairs, allocate_field_triplets, init_field
    use parameters, only : allocate_chain_parameters, init_chain_parameters, init_constants, &
        init_vnucl_type, make_isrhoselfconsistent, set_size_neq, chainperiod, chaintype, &
        chainmethod, maxnchainsrotations, maxnchainsrotationsxy, zpol, segcmfname, &
        isVdW, dielect_env, write_struct
    use matrices, only : init_matrices
    use chains, only : allocate_chains, ismonomer_chargeable, type_of_charge, segcm, mapping_num_to_char, &
        index_phos, inverse_index_phos, len_index_phos, maxneigh, no_overlapchain, distphoscutoff
    use chaingenerator, only : make_sequence_chain, make_charge_table, make_type_of_charge_table, &
        make_segcom, set_properties_chain, set_mapping_num_to_char, make_chains, &
        find_phosphate_location, write_chain_struct
    use VdW, only : make_VdWeps
    ! use listfcn, only : set_fcn
    ! use surface, only : init_surface
    use myio, only : read_inputfile, set_value_isVdW, write_chain_config
    use myutils, only : lenText, LogUnit, open_logfile, close_logfile, print_to_log, error_handler
    !use dielectric_const, only : set_dielect_fcn
    use modfcnFeexpl, only : init_var_compute_fdisPPP
    use chains_spheroid, only : init_chain_NCP
    use chaingenerator_spheroid, only : init_spheroid, build_NCP_chain
    use volume, only : make_geometry, nx, ny, nz, nsurf, delta, &
        allocate_index_neighbors_phos, make_table_index_neighbors_phos
    
    implicit none

    integer :: info
    character(len=lenText) :: text, istr
    character(len=20) :: fname
    integer :: phoscutoff

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

   
    if(systype=="nucl_ionbin_Mg") then ! auxiliary array index_phos
        call find_phosphate_location(index_phos,inverse_index_phos,len_index_phos) 
        call allocate_index_neighbors_phos(maxneigh,len_index_phos)
        phoscutoff=int(distphoscutoff/delta)+2
        call make_table_index_neighbors_phos(phoscutoff,len_index_phos,index_phos) 
    endif

    if(systype=="nucl_ionbin_MgA") then 
        phoscutoff=int(distphoscutoff/delta)+2
    endif   
    
    if(systype=="nucl_ionbin_Fe") then 
        phoscutoff=int(distphoscutoff/delta)+2
   
    endif   
    
    if(systype=="nonuclcp") then 
        no_overlapchain(1)=.true. ! otherwise loop does not start !!
    endif    
    
   
    call make_isrhoselfconsistent(info)
    call set_size_neq()             ! number of non-linear equation neq
   
    call write_chain_config()
    call write_chain_struct(write_struct,info) 
   ! if(write_struct) call error_handler(1,"stop after write_chains_struct")

   

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
