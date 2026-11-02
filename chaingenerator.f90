!---------------------------------------------------------------| 
! chainsgenerator.f90:                                          |       
! Generator chains conformations either using                   | 
! 1. RIS algorithm                                              | 
! 2. Reads in conforamtion from file                            |
! --------------------------------------------------------------|


module chaingenerator 

    use precision_definition   
    implicit none

    ! integer, parameter :: lenfname=40
    integer :: conf_write
    integer :: conf_write_chain
    real(dp) :: xgraftloop(3,2)

    real(dp),          parameter :: eps_equilat=1.0e-8_dp
    character(len=80), parameter :: fmt3xyz = "(3ES15.5E2)"
    ! 

    logical, parameter :: COMOLD =.FALSE.

    private 

    public :: make_chains, make_chains_mc,read_chains_xyz
    public :: make_charge_table, make_segcom, make_sequence_chain, make_type_of_charge_table
    public :: set_mapping_num_to_char, set_properties_chain, write_chain_struct
    public :: find_phosphate_location
    public :: write_indexchain_histone, test_index_histone

    ! made public to be used in chaingenerator_spheroid
    public :: open_chain_xyz, write_chain_xyz
    public :: conf_write_chain, open_chain_energy, open_chain_lammps_trj, write_chain_lammps_trj

    private :: eps_equilat,fmt3xyz, COMOLD

contains

! Main chain generator routine selects type of method to use 
! based on value 
! input character(len=15) chainmethod
!       character(len=15) systype

subroutine make_chains(chainmethod,systype)

    use myutils, only : print_to_log, LogUnit, lenText, error_handler
    use myio, only : myio_err_chainmethod

    character(len=15), intent(in) :: chainmethod
    character(len=15), intent(in) :: systype 

    integer :: i, info
    character(len=lenText) :: text, istr

    if(systype=="nonuclcp") return ! alternate return

    info=0
    select case (chainmethod)
    case ("MC")
        call make_chains_mc()
    case ("FILE_XYZ")      
        call read_chains_xyz(systype,info)       
    case default
        text="chainmethod not equal to MC or FILE_XYZ"            
        call print_to_log(LogUnit,text)
        print*,text
        info = myio_err_chainmethod
    end select

    call error_handler(info,"make_chains")
    
    
end subroutine make_chains


!  Uses Monte Carlo RIS algorithme to generate of 'cuantas' 
!  polymer configurations of polymer chain anchored onto a flat surface 

subroutine make_chains_mc()
  
    use globals
    use chains
    use random, only : seed
    use parameters, only : lseg, write_mc_chains, isVdW, isVdWintEne
    use parameters, only : maxnchainsrotations, maxnchainsrotationsxy
    use parameters, only : pbc_chains
    use volume, only : nx, ny, nz, delta, geometry
    use volume, only : coordinateFromLinearIndex, linearIndexFromCoordinate
    use volume, only : ut, vt
    use myutils, only : print_to_log, lenText, LogUnit
    use cadenas_linear
    use cadenas_sequence

    !  use chain_spheroid, only : orientation_spheriod_chain

    !     .. variable and constant declaractions      

    integer :: i,j,k,s,g,gn      ! dummy indices
    integer :: idx               ! index label
    integer :: ix,iy,idxtmp,ntheta
    integer :: nchains           ! number of rotations
    integer :: maxnchains        ! number of rotations
    integer :: maxntheta         ! maximum number of rotation in xy-plane
    integer :: conf              ! counts number of conformations
    integer :: allowedconf
    real(dp) :: chain(3,nseg,200) ! chain(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp) :: chain_rot(3,nseg)
    real(dp) :: x(nseg), y(nseg), z(nseg) ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates
    real(dp) :: xpp(nseg), ypp(nseg), zpp(nseg)  
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm         ! sizes box
    real(dp) :: xpt,ypt          ! coordinates
    real(dp) :: theta, theta_angle
    character(len=lenText) :: text, istr
    integer  :: xi,yi,zi ,un_trj, un_ene, segcenter
    real(dp) :: energy   
    logical :: saw 
    integer :: info    
   
    !  .. executable statements
    !  .. initializations of variables     
       
    conf = 1                 ! counter for conformations
    !seed = 435672*(rank+1)   ! seed for random number generator  different on each node
    seed = 435672
    maxnchains = maxnchainsrotations
    maxntheta = maxnchainsrotationsxy         ! maximum number of rotation in xy-plane  
    !theta_angle = 2.0_dp*pi/maxntheta
    
    Lz = nz*delta            ! maximum height box 
    Lx = nx*delta            ! maximum width box 
    Ly = ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp           ! center box
    ycm= Ly/2.0_dp
    zcm= LZ/2.0_dp
   
    energy=0.0_dp

    if(write_mc_chains) then 
        conf_write_chain=0
        un_trj= open_chain_lammps_trj(info)
        un_ene= open_chain_energy()
    endif   

    if(isHomopolymer.eqv..FALSE.) then 
        allocate(lsegseq(nseg))
        call make_lsegseq(lsegseq,nseg)
    endif    
  
    do while (conf<=max_confor)

        nchains= 0      ! init zero 
        if(isHomopolymer) then 
            call make_linear_chains(chain,nchains,maxnchains,nseg,lseg) ! chain generator f90
        else
            call make_linear_seq_chains(chain,nchains,maxnchains,nseg) 
        endif  

        if(write_mc_chains) then 
            energy=0.0_dp
            write(un_ene,*)energy  
            call write_chain_lammps_trj(un_trj,chain,nchains)  
        endif    
   
            
        if(geometry=="cubic") then

            do j=1,nchains   
            
                do s=1,nseg                          !  transforming form real- to lattice coordinates
                    z(s) = chain(1,s,j)+zcm
                    x(s) = chain(2,s,j)+xcm
                    y(s) = chain(3,s,j)+ycm
                enddo
                
                do s=1,nseg 

                    if(pbc_chains) then 
                        ! .. periodic boundary conditions in x-direction and y-direction 
                        chain_rot(2,s) = pbc(x(s),Lx)
                        chain_rot(3,s) = pbc(y(s),Ly)
                        chain_rot(1,s) = pbc(z(s),Lz)        
                    else 
                        chain_rot(2,s) = x(s)
                        chain_rot(3,s) = y(s)
                        chain_rot(1,s) = z(s)   
                    endif 

                    ! .. transforming form real- to lattice coordinates                 
                    xi = int(chain_rot(2,s)/delta)+1
                    yi = int(chain_rot(3,s)/delta)+1
                    zi = int(chain_rot(1,s)/delta)+1
                    
                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                    indexchain(s,conf) = idx
                    if(idx<=0) then
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                    endif
                enddo   

               ! call orientation_spheriod_chain(chain_rot,nseg)            
        
                conf = conf +1 

            
            enddo             ! end loop over nchains

        else if(geometry=="prism") then
            
            do j=1,nchains   
            
                do s=1,nseg                     !  transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)+zcm
                    xp(s) = chain(2,s,j)+xcm
                    yp(s) = chain(3,s,j)+ycm
                enddo

                do s=1,nseg
                        
                    ! .. prism coordinates

                    x(s) = ut(xp(s),yp(s))
                    y(s) = vt(xp(s),yp(s))
                    
                    ! .. periodic boundary conditions in x-direction and y-direction an z-direction 
                    chain_rot(2,s) = pbc(x(s),Lx)
                    chain_rot(3,s) = pbc(y(s),Ly)
                    chain_rot(1,s) = pbc(zp(s),Lz)       

                    ! .. transforming form real- to lattice coordinates                 
                    xi = int(chain_rot(2,s)/delta)+1
                    yi = int(chain_rot(3,s)/delta)+1
                    zi = int(chain_rot(1,s)/delta)+1

                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                    indexchain(s,conf) = idx
                    if(idx<=0) then
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                    endif
                    
                enddo         ! end loop over graft points

                conf = conf +1 
             
            enddo                 ! end loop over nchains

        else 
           print*,"Error: in make_chains_mc: geometry not cubic or prism"
           print*,"stopping program"
           stop
        endif
     
    enddo                     ! end while loop
    
    !  .. end chains generation 
      
   
    text='AB Chains generated'
    call print_to_log(LogUnit,text)

  
    if(isHomopolymer.eqv..FALSE.) deallocate(lsegseq)

    if(write_mc_chains) then
        close(un_trj)
        close(un_ene)
    endif   

    energychain=0.0_dp ! no internal energy 
    energychainLJ= 0.0_dp
    
end subroutine make_chains_mc


subroutine read_chains_xyz(systype,info)

    ! .. argument
    character(len=15), intent(in) :: systype
    integer, intent(out) :: info

    if(systype=="nucl_neutral_sv".or.systype=="nucl_ionbin_sv".or.systype=="nucl_ionbin_Mg" &
       .or. systype=="nucl_ionbin_MgA" .or. systype=="nucl_ionbin_Fe") then 
       print*,"read_chains_xyz_nucl_volume"
        call read_chains_xyz_nucl_volume(info)
    else
        call read_chains_xyz_nucl(info)
    endif    

end subroutine


!elemental subroutine sp_to_dp( lhs, rhs )
!      real(kind=dp), intent(out) :: lhs
!      real(kind=sp), intent(in)  :: rhs
!      character(len=*), parameter :: fmtdp = "(1pg15.7)"
!      character(len=22) :: s 
!      write( s, fmtdp ) rhs
!      read( s, fmt=*) lhs
!end subroutine


! Reads conformations from a file called traj.<rank>.xyz
! Reads energy of the conformation from a file called energy.<rank>.ene
! Reading of energy file only if isChainEnergyFile==.true.
! Format conformation file : lammps trajectory file xyz file
! number of ATOMS much equal nseg  
! 
! input/output  : integer  info : info = 0 chain read in correct inf0/=0 an error occured
! returns : indexchain and weightchain
!         : structural quantities:  
!           Rgsqr   = radius of gyration,
!           Rendsqr = end_to_end distance,
!           bond_angle,
!           dihedral_angle,
!           nucl_spacing = nucleosomal_spacing
!           no_overalchain = logical indicating overlap 

subroutine read_chains_xyz_nucl(info)

    !     .. variable and constant declaractions  
                                                                                  
    use globals, only : nsize, nseg, nsegsource, nsegtypes, s_begin, s_end, nsegAA 
    use globals, only : nnucl, cuantas, cuantas_no_overlap, max_confor , runtype
    use chains, only : indexchain, logweightchain, no_overlapchain, segcm, sgraftpts
    use chains, only : energychain, energychainLJ, energychainLJ0, unitvector_triplets
    use chains, only : Rgsqr, Rendsqr, bond_angle, dihedral_angle, nucl_spacing
    use chains, only :  Asphparam
    use eigenvalues, only : asphericty_parameter
    use parameters, only : VdWscale, isChainEnergyFile, isVdW , isvdwintene 
    use parameters, only : pbc_chains, readinchains, unit_conv, write_rotations
    use volume, only :  nx, ny,nz, delta, geometry, ut, vt, linearindexfromcoordinate
    use chain_rotation, only : rotate_nucl_chain, rotate_nucl_chain_test
    use myio, only : myio_err_chainsfile, myio_err_energyfile, myio_err_index
    use myio, only : myio_err_conf, myio_err_nseg, myio_err_geometry, myio_err_equilat
    use myutils,  only :  print_to_log, LogUnit, lenText, newunit
    use Eqtriangle, only : equilateralness_vectors
    use VdW_potential, only : GBenergyeffective, init_GBenergyeffective
    use VdW_potential, only : make_com_nucl_rotation

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    integer :: i,j,s,rot,g,gn      ! dummy indices
    integer :: idx                 ! index label
    integer :: ix,iy,iz,idxtmp,ntheta
    integer :: nchains              ! number of rotations
    integer :: maxnchains           ! number of rotations
    integer :: maxntheta            ! maximum number of rotation in xy-plane
    integer :: conf,conffile        ! counts number of conformations  
    integer :: nsegfile             ! nseg in chain file      
    integer :: cuantasfile          ! cuantas in chain file                                              
    real(dp) :: chain(3,nseg),chain_rot(3,nseg),chain_pbc(3,nseg)  ! chains(x,i)= coordinate x of segement i
    real(dp) :: xseg(3,nseg)
    real(dp) :: x(nseg), y(nseg), z(nseg)    ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates 
    real(dp) :: xpp(nseg),ypp(nseg)
    integer  :: xi,yi,zi
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm ! sizes box and center of mass box
    real(dp) :: xpt,ypt              ! coordinates
    real(dp) :: xc,yc,zc               
    real(dp) :: energy, energyLJ                                            
    character(len=25) :: fname
    integer :: ios, rankfile, iosene
    character(len=30) :: str
    real(dp) :: scalefactor
    integer :: un,unw,un_ene ! unit number
    logical :: exist
    character(len=lenText) :: text,istr
    real(dp) :: d_type_num, d_atom_num
    integer :: i_type_num, i_atom_num
    logical :: isReadGood
    real(dp) :: equilat, equilat_rot
    integer :: ii
    integer :: s_local
    integer :: segnumAAstart(nnucl), segnumAAend(nnucl) ! segment numbers first/last AAs 
    logical :: no_overlap
    real(dp) :: rcom(3,nnucl)       ! hold mean or com of nucleosome coordinate
    real(dp) :: gyr_tensor(3,3) 

    ! .. executable statements   

    info=0

    ! .. open file   

    rankfile=0                                                                                    
    
    write(istr,'(I4)')rankfile
    fname='traj.'//trim(adjustl(istr))//'.xyz'
    
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,' traj file :',fname,' does not exit'
        info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    endif

    if(isChainEnergyFile) then
        write(istr,'(I4)')rankfile
        fname='energy.'//trim(adjustl(istr))//'.ene'
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un_ene),file=fname,status='old',iostat=ios)
        else
            text='energy.'//trim(adjustl(istr))//'.ene file does not exit'
            print*,text
            info = myio_err_energyfile
            return
        endif
        if(ios >0 ) then
            print*, 'Error opening file : iostat =', ios
            info = myio_err_energyfile
            return
        endif
    endif    

    conf=1                    ! counter for conformations                                                           
    conffile=0                ! counter for conformations in file    
    ios=0
    scalefactor=unit_conv

    energy=0.0_dp
    energyLJ= 0.0_dp
    no_overlap = .true.
   
    ios=0

    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp          ! center box
    ycm= Ly/2.0_dp
    zcm= Lz/2.0_dp

    
    if(isVdW) then 
                
        call compute_segnumAAstart(nseg,nsegtypes,nnucl,segnumAAstart)
        call compute_segnumAAend(nnucl,nsegAA,segnumAAstart,segnumAAend)

        call init_GBenergyeffective(segcm,nnucl,segnumAAstart,segnumAAend)
    endif
        
    isReadGood=.true. 

    do while ((conf<=max_confor).and.isReadGood)
    

        if(conf.ne.1) then ! skip lines
            read(un,*,iostat=ios)
            read(un,*,iostat=ios) 
        else               ! read preamble
            read(un,*,iostat=ios)nsegfile
            if(ios/=0) isReadGood=.false.
            read(un,*,iostat=ios) ! skip line
            if(ios/=0) isReadGood=.false.   
        
            if(nsegfile.ne.nsegsource) then
                text="nseg chain file not equal internal nseg : stop program"
                print*,"nsegfile=",nsegfile," nsegsource=",nsegsource
                call print_to_log(LogUnit,text)
                info=myio_err_nseg
                return
            endif
        endif

        s_local=0

        do s=1,nsegsource              ! .. read form  trajectory file

            read(un,fmt3xyz,iostat=ios)xc,yc,zc
            if(ios/=0) isReadGood=.false.

            if(s_begin<=s.and.s<=s_end) then ! filter 

                s_local=s_local+1
                xseg(1,s_local) = xc*scalefactor
                xseg(2,s_local) = yc*scalefactor
                xseg(3,s_local) = zc*scalefactor

            endif
        enddo


        if(isChainEnergyFile) read(un_ene,*,iostat=ios)energy

        !if(isVdW) EnergyLJ = VdWpotentialenergy(xseg)

        if(isReadGood) then ! read was succesfull  

            conffile=conffile +1 
           
            do s=1,nseg        
                chain(1,s) = xseg(1,s)-xseg(1,sgraftpts(1)) 
                chain(2,s) = xseg(2,s)-xseg(2,sgraftpts(1)) 
                chain(3,s) = xseg(3,s)-xseg(3,sgraftpts(1)) 
            enddo

            ! rotate chain 
            
            ! call rotate_nucl_chain(chain,chain_rot,sgraftpts,nseg)
            call rotate_nucl_chain_test(chain,chain_rot,sgraftpts,nseg,write_rotations)                

            ! extra test rotation
            equilat=equilateralness_vectors(chain(:,sgraftpts(1)),chain(:,sgraftpts(2)),chain(:,sgraftpts(3)))
            equilat_rot=equilateralness_vectors(chain_rot(:,sgraftpts(1)),chain_rot(:,sgraftpts(2)),chain_rot(:,sgraftpts(3)))
 
            if( abs(equilat-equilat_rot)>=eps_equilat )  then
                text="Equilateralness test failed : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_equilat 
                return
            endif    
            

            select case (geometry)
            case ("cubic")

                do s=1,nseg

                    if(pbc_chains) then 
                        chain_pbc(1,s) = pbc(chain_rot(1,s)+xcm,Lx) ! periodic boundary conditions in x and y and z direcxtion 
                        chain_pbc(2,s) = pbc(chain_rot(2,s)+ycm,Ly)
                        chain_pbc(3,s) = pbc(chain_rot(3,s)+zcm,Lz)  
                    else
                        chain_pbc(1,s) = chain_rot(1,s)+xcm 
                        chain_pbc(2,s) = chain_rot(2,s)+ycm
                        chain_pbc(3,s) = chain_rot(3,s)+zcm
                    endif                    

                    ! transforming form real- to lattice coordinates                 
                    xi = int(chain_pbc(1,s)/delta)+1
                    yi = int(chain_pbc(2,s)/delta)+1
                    zi = int(chain_pbc(3,s)/delta)+1
                        
                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                        
                    indexchain(s,conf) = idx

                    if(idx<=0.or.idx>nsize) then   

                        text="Conformation outside box:"
                        call print_to_log(LogUnit,text)  
                        print*,text                          
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        info= myio_err_index
                        return
                    endif
                
                enddo
               
                if(isVdW) EnergyLJ = GBenergyeffective(chain_pbc,nnucl,no_overlap)

                call make_com_nucl_rotation(chain_pbc,nnucl,unitvector_triplets,rcom)

                energychainLJ(conf)    = energyLJ
                energychain(conf)      = energy
                no_overlapchain(conf)  = no_overlap 
                Rgsqr(conf)            = radius_gyration_com_rotation(rcom,nnucl)
                Rendsqr(conf)          = end_to_end_distance_com_rotation(rcom,nnucl)
                bond_angle(:,conf)     = bond_angles_com_rotation(rcom,nnucl)
                dihedral_angle(:,conf) = dihedral_angles_com_rotation(rcom,nnucl)
                nucl_spacing(:,conf)   = nucleosomal_spacing_com_rotation(rcom,nnucl)
                gyr_tensor             = gyr_tensor_com_rotation(rcom,nnucl)
                Asphparam(conf)        = Asphericty_parameter(Rgsqr(conf),gyr_tensor)
         
                if(COMOLD) then 
                    Rgsqr(conf)            = radius_gyration_com(chain_pbc,nnucl,segcm)
                    Rendsqr(conf)          = end_to_end_distance_com(chain_pbc,nnucl,segcm)
                    bond_angle(:,conf)     = bond_angles_com(chain_pbc,nnucl,segcm)
                    dihedral_angle(:,conf) = dihedral_angles_com(chain_pbc,nnucl,segcm)
                    nucl_spacing(:,conf)   = nucleosomal_spacing_com(chain_pbc,nnucl,segcm)
                    gyr_tensor             = gyr_tensor_com(chain_pbc,nnucl,segcm)
                    Asphparam(conf)        = Asphericty_parameter(Rgsqr(conf),gyr_tensor)
                endif    

                conf=conf+1   
                                    
            case("prism") 
                    
                do s=1,nseg
                    
                    xp(s) = chain_rot(1,s)+xcm
                    yp(s) = chain_rot(2,s)+ycm
                    zp(s) = chain_rot(3,s)+zcm

                    ! .. transformation to prism coordinates 
                    xpp(s) = ut(xp(s),yp(s))
                    ypp(s) = vt(xp(s)+ycm,yp(s))

                    if(pbc_chains) then 
                         ! .. periodic boundary conditions in u and v ands z direction
                         chain_pbc(1,s) = pbc(xpp(s),Lx) 
                         chain_pbc(2,s) = pbc(ypp(s),Ly)
                         chain_pbc(3,s) = pbc(zp(s),Lz)        
                    else
                         chain_pbc(1,s) = xpp(s)  
                         chain_pbc(2,s) = ypp(s)
                         chain_pbc(3,s) = zp(s)
                    endif

                    ! .. transforming form real- to lattice coordinates                 
                    xi = int(chain_pbc(1,s)/delta)+1
                    yi = int(chain_pbc(2,s)/delta)+1
                    zi = int(chain_pbc(3,s)/delta)+1
                        
                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                        
                    indexchain(s,conf) = idx

                    if(idx<=0.or.idx>nsize) then    
                        text="Conformation outside box:"
                        call print_to_log(LogUnit,text)  
                        print*,text                          
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        info= myio_err_index
                        return
                    endif

                enddo
                
                if(isVdW) EnergyLJ = GBenergyeffective(chain_pbc,nnucl,no_overlap)

                call make_com_nucl_rotation(chain_pbc,nnucl,unitvector_triplets,rcom)

                energychainLJ(conf)    = energyLJ
                energychain(conf)      = energy
                no_overlapchain(conf)  = no_overlap 
                Rgsqr(conf)            = radius_gyration_com_rotation(rcom,nnucl)
                Rendsqr(conf)          = end_to_end_distance_com_rotation(rcom,nnucl)
                bond_angle(:,conf)     = bond_angles_com_rotation(rcom,nnucl)
                dihedral_angle(:,conf) = dihedral_angles_com_rotation(rcom,nnucl)
                nucl_spacing(:,conf)   = nucleosomal_spacing_com_rotation(rcom,nnucl)
                gyr_tensor             = gyr_tensor_com_rotation(rcom,nnucl)
                Asphparam(conf)        = Asphericty_parameter(Rgsqr(conf),gyr_tensor)
      
                
                if(COMOLD) then 
                    Rgsqr(conf)            = radius_gyration_com(chain_pbc,nnucl,segcm)
                    Rendsqr(conf)          = end_to_end_distance_com(chain_pbc,nnucl,segcm)
                    bond_angle(:,conf)     = bond_angles_com(chain_pbc,nnucl,segcm)
                    dihedral_angle(:,conf) = dihedral_angles_com(chain_pbc,nnucl,segcm)
                    nucl_spacing(:,conf)   = nucleosomal_spacing_com(chain_pbc,nnucl,segcm)
                    gyr_tensor             = gyr_tensor_com(chain_pbc,nnucl,segcm)
                    Asphparam(conf)        = Asphericty_parameter(Rgsqr(conf),gyr_tensor)
      
                endif    
                    
                conf = conf + 1   

            case default
                text="Error: in make_chains_XYZ_nucl geometry not cubic or prism: stopping program"
                call print_to_log(LogUnit,text)
                info = myio_err_geometry
                return 
                    
            end select

        endif   ! read 

    enddo       ! end while loop                                                                                                          
    ! end chains generation    
    
    conf=conf-1  ! lower by one  

    if(conf<max_confor) then
        print*,"subroutine make_chains_XYZ_nucl :" 
        print*,"conf     = ",conf," less then imposed max cuantas     = ",max_confor
        print*,"conffile = ",conffile
        cuantas=conf   
        info=myio_err_conf        
    else
        text="Chains generated: subroutine make_chains_XYZ_nucl"
        call print_to_log(LogUnit,text)
        readinchains=conffile
        info=0
    endif

    write(istr,'(L2)')isVdWintEne
    text="isVdWintEne = "//trim(adjustl(istr))
    call print_to_log(LogUnit,text)

    if(.not.(isChainEnergyFile)) energychain=0.0_dp

    close(un) 
    if(isChainEnergyFile) close(un_ene)

    cuantas_no_overlap = compute_cuantas_no_overlap(cuantas,no_overlapchain)

    !call normed_weightchains() 
    text="No reweighting chains."
    call print_to_log(LogUnit,text)
    logweightchain=0.0_dp    

    deallocate(energychain) ! free unused variables 

    if(runtype=="rangeVdWeps") then 
        allocate(energychainLJ0(cuantas))
        do i=1,cuantas
            energychainLJ0(i)=energychainLJ(i)
        enddo 
    else 
        ! energy of conformation fixed: VdWscale%val set the energy scale relative to epsilonS
        ! for runtype=="rangVdWeps" loop over values of VdWscale%val     
        do i=1,cuantas
            energychainLJ(i)=VdWscale%val*energychainLJ(i)
        enddo
    endif           
    
end subroutine read_chains_xyz_nucl


! Reads base conformations from a file called traj.<rank>.xyz.
! Reads AA side chains of  histone from a file ????
! Reads energy of the conformation from a file called energy.<rank>.ene.
! Reading of energy file only if isChainEnergyFile==.true.
! Format conformation file : lammps trajectory file xyz file
! number of ATOMS much equal nseg  
! 
! input/output  : integer  info : info = 0 chain read in correct info/=0 an error occured
! returns : indexconf and weightchain for fcnnucl_neutral_sv 
!         : structural quantities:  
!           Rgsqr   = radius of gyration,
!           Rendsqr = end_to_end distance,
!           bond_angle =
!           dihedral_angle =
!           nucl_spacing = nucleosomal_spacing
!           no_overalchain = logical indicating overlap 

subroutine read_chains_xyz_nucl_volume(info)

    !     .. variable and constant declaractions                                                                                    
    use globals, only : nsize,nseg, nsegsource, nsegtypes, s_begin, s_end , nsegAA, DEBUG 
    use globals, only : nnucl, cuantas, cuantas_no_overlap, max_confor, set_confor, runtype, systype
    use chains, only : var_darray
    use chains, only : indexconf, nelem, nelemAA, typeAA, elem_charge, nucl_elem_type, distphoscutoff
    use chains, only : logweightchain, no_overlapchain, segcm, sgraftpts
    use chains, only : energychain, energychainLJ, energychainLJ0, unitvector_triplets, orientation_triplets 
    use chains, only : Rgsqr, Rendsqr, bond_angle, dihedral_angle, nucl_spacing, Asphparam
    use chains, only : allocate_indexconf, allocate_indexconfpair, allocate_nneighbor
    use chains, only : allocate_indexconftriplet, allocate_ntriplet
    use eigenvalues, only : asphericty_parameter
    use parameters, only : VdWscale, isChainEnergyFile, isVdW, isvdwintene, mtpdbfname, orientfname
    use parameters, only : pbc_chains, readinchains, unit_conv, vnucl, write_rotations
    use volume, only   :  nx, ny,nz, delta, geometry, ut, vt, linearindexfromcoordinate, coordtoindex
    use chain_rotation, only : rotate_nucl_chain, rotate_nucl_chain_test, orientation_coordinates
    use chain_rotation, only : orientation_vector_ref, orientation_vector, rotate_chain_elem
    use chain_rotation, only : check_chain_elem_index_and_chain
    use myio, only     : myio_err_chainsfile, myio_err_energyfile, myio_err_index
    use myio, only     : myio_err_conf, myio_err_nseg, myio_err_geometry, myio_err_equilat
    use myutils, only  : print_to_log, LogUnit, lenText, newunit
    use Eqtriangle, only : equilateralness_vectors
    use myutils, only  : error_handler
    use VdW_potential, only : GBenergyeffective, init_GBenergyeffective
    use VdW_potential, only : make_com_nucl_rotation
    use listfcn, only : maximum_xnucl
    !  new 
    use chains_spheroid, only : chain_NCP

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    integer :: i,j,s,sprime,rot,g,gn,k,sAA ,tPhos ! dummy indices
    integer :: idx, idx_tmp                   ! index label
    integer :: ix,iy,iz 
    integer :: conf,conffile        ! counts number of conformations  
    integer :: nsegfile             ! nseg in chain file      
    integer :: cuantasfile          ! cuantas in chain file                                              
    real(dp) :: chain(3,nseg),chain_rot(3,nseg),chain_pbc(3,nseg),chain_pbc_tmp(3),chain_tmp(3)    ! chains(x,i)= coordinate x of segement i
    real(dp) :: xseg(3,nseg)
    real(dp) :: x(nseg), y(nseg), z(nseg)    ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates 
    real(dp) :: xpp(nseg),ypp(nseg)
    real(dp) :: rtranslate(3)
    integer  :: xi,yi,zi, ri(3)
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm, Lr(3), rcm(3) ! sizes box and center of mass box
    real(dp) :: xpt,ypt              ! coordinates
    real(dp) :: xc,yc,zc             ! coordinates            
    real(dp) :: energy, energyLJ, energyLJ_comb                                            
    character(len=25) :: fname
    character(lenText):: fname2
    integer :: ios, rankfile, iosene
    character(len=30) :: str
    real(dp) :: scalefactor
    integer :: un,unw,un_ene ! unit number
    logical :: exist
    character(len=lenText) :: text,istr
    real(dp) :: d_type_num, d_atom_num
    integer :: i_type_num, i_atom_num
    logical :: isReadGood, isVolfracLargerOne
    integer :: nrotpts
    real(dp) :: equilat,equilat_rot
    integer :: nelem2(3),nsegAA2 
    integer  :: segnumAAstart(nnucl), segnumAAend(nnucl) ! segment numbers first/last AAs 
    integer  :: orient_triplet_ref(3)
    real(dp) :: orient_vector_ref(3)
    real(dp) :: orient_vectors(3,nnucl)
    real(dp) :: triangle_orient_vectors(3,3,nnucl)

    ! chain_elem(k,s)%elem(j) the k=(1,2,3)=(x,yz) cartesian coordinate of jth chain element of sth AA of reference histone protein 
    ! chain_elem_rot(k,s,n)%elem(j) the k=(1,2,3)=(x,yz) cartesian coordinate of jth chain element of sth AA of nth nucleosome 
    ! chain_elem_index(k,s)%elem(j) the k=(1,2,3)=(x,yz) cartesian coordinate of jth chain element of sth segement (AA and DNA) 
   
    type(var_darray), dimension(:,:), allocatable   :: chain_elem  
    type(var_darray), dimension(:,:,:), allocatable :: chain_elem_rot
    type(var_darray), dimension(:,:), allocatable   :: chain_elem_index
    real(dp) :: rcom(3,nnucl)
    
    integer :: un_traj, info_traj
    real(dp) :: chain_lammps(3,nseg,1)
    real(dp) :: sqrdist, sqrDphoscutoff ! square distance and square cutoff for pair distances of phosphates
    ! integer, dimension(:,:), allocatable   :: list_of_pairs
    integer :: max_range_nneigh 
    integer :: s_local 
    logical :: no_overlap, isCheck
    real(dp) :: gyr_tensor(3,3)
    ! .. executable statements   


    real(dp) :: chain_NCP_lammps(3,nseg,1)

    info=0

    ! .. open file   

    rankfile=set_confor                                                                                     
    
    write(istr,'(I4)')rankfile
    fname='traj.'//trim(adjustl(istr))//'.xyz'
    
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,' traj file :',fname,' does not exit'
        info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    endif

    if(isChainEnergyFile) then
        write(istr,'(I4)')rankfile
        fname='energy.'//trim(adjustl(istr))//'.ene'
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un_ene),file=fname,status='old',iostat=ios)
        else
            text='energy.'//trim(adjustl(istr))//'.ene file does not exit'
            print*,text
            info = myio_err_energyfile
            return
        endif
        if(ios >0 ) then
            print*, 'Error opening file : iostat =', ios
            info = myio_err_energyfile
            return
        endif
    endif    

    conf=1                    ! counter for conformations                                                           
    conffile=0                ! counter for conformations in file 
    conf_write=0  
    conf_write_chain=0 
    
    ios=0
    scalefactor=unit_conv

    energy = 0.0_dp
    energyLJ = 0.0_dp
    no_overlap = .true.
  
    ios=0

    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp          ! center box
    ycm= Ly/2.0_dp
    zcm= Lz/2.0_dp

    rcm=(/xcm,ycm,zcm/)
    Lr=(/Lx,Ly,Lz/)

    sqrDphoscutoff=distphoscutoff**2
    !print*,"sqrDphoscutoff=",sqrDphoscutoff
    max_range_nneigh = 10                             ! assume at maximum there  at 10 ! 
    ! allocate(list_of_pairs(nseg,max_range_nneigh))   ! temporaly storage of neighbor index

    nrotpts=sgraftpts(1)  ! nucleosome id /segment number around which to rotate whole conformation
    !nrotpts=rotation_triplets(1)
    
    ! return position (chain_elem) and number (nelem) of elements of every AA segment
                          

    call read_nucl_elements(mtpdbfname,nsegAA,nelemAA,chain_elem,typeAA,vnucl,nucl_elem_type,elem_charge,info)
    if(info/=0) return

    if(DEBUG) call print_nucl_elements(nsegAA,nelemAA,chain_elem)
    
    call read_nucl_orient_triplets(orientfname,nnucl,orientation_triplets,info)

    call compute_segnumAAstart(nseg,nsegtypes,nnucl,segnumAAstart)
    call compute_segnumAAend(nnucl,nsegAA,segnumAAstart,segnumAAend)

    if(DEBUG) then 
        do i=1,nnucl
            print*," sAAstart ",segnumAAstart(i)," sAAend ",segnumAAend(i)
        enddo    
    endif    

    ! compare ORDER of sgraftpts and orientation triplets !!!!!!!
    ! check indices

    do i=1,3
        orient_triplet_ref(i)=orientation_triplets(1,i)-segnumAAstart(1)+1
    enddo 

    if(DEBUG) then
        print*,"Module : chaingenerator : read_chains_xyz_nucl_volume "
        print*,"orient_triplets"
        do i=1,3
            print*,"orient_triplet_ref(",i,")=",orient_triplet_ref(i)," original ",orientation_triplets(1,i)
            print*,"orientation_triplets=",orientation_triplets
        enddo
    endif    

    ! make chain_elem relative to CA of AA of CM rotation i.e., sgraftpts
    call shift_nucl_elements(nsegAA,nelemAA,orient_triplet_ref(1),chain_elem)
   
    if(DEBUG) call print_nucl_elements(nsegAA,nelemAA,chain_elem)
   
    call make_nelem(nseg,nsegAA,nnucl,segnumAAstart,nelemAA,nelem)
    call allocate_indexconf(cuantas,nseg,nelem)
    call allocate_nucl_chain_elements(nnucl,nsegAA,nelemAA,chain_elem_rot) 
    call allocate_chain_elements(nseg,nelem,chain_elem_index)
    call orientation_vector_ref(chain_elem,orient_triplet_ref,orient_vector_ref)
    
    ! pairs variable 
    if(systype=="nucl_ionbin_Mg".or. systype=="nucl_ionbin_MgA".or. systype=="nucl_ionbin_Fe") then 
        call allocate_indexconfpair(cuantas,nseg)
        call allocate_nneighbor(cuantas,nseg)
        tPhos = find_type_phosphate()
    endif
    ! triplet  variable 
    if(systype=="nucl_ionbin_Fe") then 
        print*,"allocate_indexconftriplet"
        call allocate_indexconftriplet(cuantas,nseg)
        call allocate_ntriplet(cuantas,nseg)
    endif

    if(isVdW) call init_GBenergyeffective(segcm,nnucl,segnumAAstart,segnumAAend) 
    
    isReadGood=.true. 
   

    do while ((conf<=max_confor).and.isReadGood)
    
        if(conf.ne.1) then ! skip lines
            read(un,*,iostat=ios)
            read(un,*,iostat=ios) 
        else               ! read preamble
            read(un,*,iostat=ios)nsegfile
            if(ios/=0) isReadGood=.false.
            read(un,*,iostat=ios) ! skip line
            if(ios/=0) isReadGood=.false.   
        
            if(nsegfile.ne.nsegsource) then 
                text="nseg chain file not equal internal nseg : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_nseg 
                return
            endif    
        endif
        
        s_local=0

        do s=1,nsegsource              ! .. read form  trajectory file
    
            read(un,fmt3xyz,iostat=ios)xc,yc,zc
            if(ios/=0) isReadGood=.false. 
                          
            if(s_begin<=s.and.s<=s_end) then ! filter 

                 s_local=s_local+1
                 xseg(1,s_local) = xc*scalefactor 
                 xseg(2,s_local) = yc*scalefactor  
                 xseg(3,s_local) = zc*scalefactor 

            endif
        enddo

        if(isChainEnergyFile) read(un_ene,*,iostat=ios)energy

        !if(isVdW)  energyLJ    = GBenergyeffective(xseg,nnucl,segcm,segunitvector) 
        
        if(isReadGood) then ! read was succesfull  

            chain_lammps(:,:,1)=xseg

            if(DEBUG)then
                un_traj=open_chain_lammps_trj(info_traj,"input")
                call write_chain_lammps_trj(un_traj,chain_lammps,1)
                close(un_traj)
            endif      

            conffile=conffile +1 
           
            ! -1. translation such that chain(:,nrotpts) in origin

            do s=1,nseg        
                chain(1,s) = xseg(1,s)-xseg(1,nrotpts) 
                chain(2,s) = xseg(2,s)-xseg(2,nrotpts) 
                chain(3,s) = xseg(3,s)-xseg(3,nrotpts) 
            enddo
          

            ! 0. rotate chain 
            
            call rotate_nucl_chain_test(chain,chain_rot,sgraftpts,nseg,write_rotations) 
            
            ! 0.a. extra test rotation
            equilat=equilateralness_vectors(chain(:,sgraftpts(1)),chain(:,sgraftpts(2)),chain(:,sgraftpts(3)))
            equilat_rot=equilateralness_vectors(chain_rot(:,sgraftpts(1)),chain_rot(:,sgraftpts(2)),chain_rot(:,sgraftpts(3)))
 
            if( abs(equilat-equilat_rot)>=eps_equilat )  then
                text="Equilateralness test failed : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_equilat 
                return
            endif  
              
            ! 1. determine orientation chain_rot : n 
            ! 2. get orientation vector of all nnucl nucleosomes in chain_rot and of reference in chain_elem
                
            call orientation_vector(chain_rot,orientation_triplets,orient_vectors)

            ! 2.a. get coordinates of for all orientation triplets of all nucleosomes   

            call orientation_coordinates(chain_rot,orientation_triplets,triangle_orient_vectors)

            ! 3. make rotation axis/angle for all nnuclnucleosome
            ! 4. apply rotation for all nucleosome to chain_elem
                
            call rotate_chain_elem(orient_vector_ref,orient_vectors,nelemAA,chain_elem,chain_elem_rot,&
                orientation_triplets,triangle_orient_vectors,segnumAAstart)


            if(DEBUG)then
                print*,"Module : chaingenerator : read_chains_xyz_nucl_volume "
                print*,"orientation triangle: "
                print*," sAA      k          chain_elem_rot(k,sAA,1)%elem(1)  chain_rot(k,sAA+3)"
                ! indices 
                do sAA=1,nsegAA
                    do k=1,3
                        print*," ",sAA," " ,k," ",chain_elem_rot(k,sAA,1)%elem(1),chain_rot(k,sAA+3)
                    enddo 
                enddo
            endif      

            ! 5. add chain+chain_elem_rot together
            
            call add_chain_rot_and_chain_elem_rot(nseg,nsegAA,nnucl,segnumAAstart,segnumAAend,nelemAA,&
                chain_rot,chain_elem_rot,chain_elem_index) 

            chain=chain_rot ! this alligns variables chain_elem_index and chain : contain same backbone conformation

            call check_chain_elem_index_and_chain(nseg,nelem,chain_elem_index,chain,isCheck)
            if(.not.isCheck) then 
                print*,"chain=chain_rot => chain_elem_index==chain",isCheck
            endif    

            ! 5.b export chain to chain_NCP used make nucleosome build
            chain_NCP=chain_rot
            chain_NCP_lammps(:,:,1)=chain_NCP

!           if(DEBUG)then
                un_traj=open_chain_elem_index_lammps_trj(info_traj)
                call write_chain_elem_index_lammps_trj(un_traj,chain_elem_index)
                close(un_traj)
            !endif     


            ! 6.a translate conformation 

            rtranslate(1) =  xcm
            rtranslate(2) =  ycm
            rtranslate(3) =  zcm
        
            call translate_chain_elem_index(nseg,nelem,chain_elem_index,rtranslate) 
            call translate_chain(nseg,chain,rtranslate)

            call check_chain_elem_index_and_chain(nseg,nelem,chain_elem_index,chain,isCheck)
            if(.not.isCheck) then 
                print*,"post translation => chain_elem_index==chain",isCheck
            endif    


            if(DEBUG)then
                un_traj=open_chain_elem_index_lammps_trj(info_traj)
                call write_chain_elem_index_lammps_trj(un_traj,chain_elem_index)
                close(un_traj)
            endif     

            !un_traj = open_chain_lammps_trj(info,"NCP")
            !call write_chain_lammps_trj(un_traj,chain_NCP_lammps,1)
            !close(un_traj)

            !un_traj = open_chain_xyz(info,"NCP")
            !call write_chain_xyz(un_traj,chain_NCP_lammps,1)
            ! close(un_traj)


            ! 6. make indexconfig i.e. place conformation on lattice
              

            select case (geometry)
            case ("cubic")

                do s=1,nseg

                    ! transforming form real- to lattice coordinates  
                    if(pbc_chains) then 
                        chain_pbc(1,s) = pbc(chain(1,s),Lx) ! periodic boundary conditions  
                        chain_pbc(2,s) = pbc(chain(2,s),Ly) 
                        chain_pbc(3,s) = pbc(chain(3,s),Lz) 
                    else
                        chain_pbc(1,s) = chain(1,s)  
                        chain_pbc(2,s) = chain(2,s)
                        chain_pbc(3,s) = chain(3,s)
                    endif                    

                    xi  = int(chain_pbc(1,s)/delta)+1
                    yi  = int(chain_pbc(2,s)/delta)+1
                    zi  = int(chain_pbc(3,s)/delta)+1

                    if(isOutsideLattice(xi,yi,zi,nx,ny,nz)) then   

                        text="Conformation outside box:"
                        call print_to_log(LogUnit,text)  
                        print*,text    
                        print*,"chain_pbc=",chain_pbc(:,s)                      
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        info= myio_err_index
                        return

                    endif

                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                    
                    idx_tmp=coordtoindex(xi,yi,zi)

                    ! test 
                    if(idx/=idx_tmp) then
                        text="idx/=idx_tmp"
                        call print_to_log(LogUnit,text)  
                        print*,text 
                        info= myio_err_index
                        return
                    endif
                        
                    indexconf(s,conf)%elem(1) = idx ! CA element

                    ! apply elements other than CA
                    do j=2,nelem(s) 

                        if(pbc_chains) then 
                            chain_pbc_tmp(1) = pbc(chain_elem_index(1,s)%elem(j),Lx) ! periodic boundary conditions 
                            chain_pbc_tmp(2) = pbc(chain_elem_index(2,s)%elem(j),Ly)
                            chain_pbc_tmp(3) = pbc(chain_elem_index(3,s)%elem(j),Lz)  
                        else
                            chain_pbc_tmp(1) = chain_elem_index(1,s)%elem(j)   
                            chain_pbc_tmp(2) = chain_elem_index(2,s)%elem(j)
                            chain_pbc_tmp(3) = chain_elem_index(3,s)%elem(j)
                        endif 

                        xi = int(chain_pbc_tmp(1)/delta)+1
                        yi = int(chain_pbc_tmp(2)/delta)+1
                        zi = int(chain_pbc_tmp(3)/delta)+1

                        if(isOutsideLattice(xi,yi,zi,nx,ny,nz)) then   

                            text="Conformation outside box:"
                            call print_to_log(LogUnit,text)  
                            print*,text  
                            print*,"chain_elem_index= ",(chain_elem_index(k,s)%elem(j),k=1,3)
                            print*,"chain_pbc= ",chain_pbc_tmp(:),"s= ",s," j= ",j                          
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                            info= myio_err_index
                            return
                        endif

                        call linearIndexFromCoordinate(xi,yi,zi,idx)

                        ! idx=coordtoindex(xi,yi,zi)
                        
                        indexconf(s,conf)%elem(j) = idx ! all other element

                    enddo   
                
                enddo ! end s loop

                if(systype=="nucl_ionbin_Mg".or.systype=="nucl_ionbin_MgA")  then
    
                    call find_phosphate_pairs(nseg,conf,tPhos,sqrDphoscutoff,chain,Lx,Ly,Lz)
                    call write_indexconfpair(nseg,conf,tPhos,sqrDphoscutoff)
                    call write_phosphate_pairs(.true.,conf,info)
                
                endif    

                if(systype=="nucl_ionbin_Fe") then
                
                    call find_phosphate_triplets(nseg,conf,tPhos,sqrDphoscutoff,chain,Lx,Ly,Lz)
                    call write_phosphate_triplets(.true.,conf,info)
                    call find_phosphate_pairs_exclude_triplets(nseg,conf,tPhos,sqrDphoscutoff,chain,Lx,Ly,Lz)
                    call write_indexconfpair(nseg,conf,tPhos,sqrDphoscutoff)
                    call write_phosphate_pairs(.true.,conf,info)
                    
                endif  

                ! if(isVdW) energyLJ = GBenergyeffective(chain_rot,nnucl,no_overlap)

                energyLJ = GBenergyeffective(chain_rot,nnucl,no_overlap)

                call maximum_xnucl(conf,isVolfracLargerOne)
                print*,"conf=",conf,"isVolfracLargerOne=",isVolfracLargerOne

                if(isVolfracLargerOne) no_overlap=.false. 
                
                call make_com_nucl_rotation(chain_pbc,nnucl,unitvector_triplets,rcom)

                energychainLJ(conf)    = energyLJ
                energychain(conf)      = energy
                no_overlapchain(conf)  = no_overlap 

                Rgsqr(conf)            = radius_gyration_com_rotation(rcom,nnucl)
                Rendsqr(conf)          = end_to_end_distance_com_rotation(rcom,nnucl)
                bond_angle(:,conf)     = bond_angles_com_rotation(rcom,nnucl)
                dihedral_angle(:,conf) = dihedral_angles_com_rotation(rcom,nnucl)
                nucl_spacing(:,conf)   = nucleosomal_spacing_com_rotation(rcom,nnucl)
                gyr_tensor             = gyr_tensor_com_rotation(rcom,nnucl)
    
                Asphparam(conf)        = Asphericty_parameter(Rgsqr(conf),gyr_tensor)
               
                if(COMOLD) then     
                    Rgsqr(conf)            = radius_gyration_com(chain_pbc,nnucl,segcm)
                    Rendsqr(conf)          = end_to_end_distance_com(chain_pbc,nnucl,segcm)
                    bond_angle(:,conf)     = bond_angles_com(chain_pbc,nnucl,segcm)
                    dihedral_angle(:,conf) = dihedral_angles_com(chain_pbc,nnucl,segcm)
                    nucl_spacing(:,conf)   = nucleosomal_spacing_com(chain_pbc,nnucl,segcm)
                    gyr_tensor             = gyr_tensor_com(chain_pbc,nnucl,segcm)
                    Asphparam(conf)        = Asphericty_parameter(Rgsqr(conf),gyr_tensor)
                endif 

                conf=conf+1       
                                    
            case("prism") 
                    
                do s=1,nseg
                    
                    xp(s) = chain_rot(1,s)+xcm
                    yp(s) = chain_rot(2,s)+ycm
                    zp(s) = chain_rot(3,s)+zcm

                    ! .. transformation to prism coordinates 
                    xpp(s) = ut(xp(s),yp(s))
                    ypp(s) = vt(xp(s),yp(s))

                    ! .. periodic boundary conditions in u and v ands z direction
                    if(pbc_chains) then 
                        chain_pbc(1,s) = pbc(xpp(s),Lx) 
                        chain_pbc(2,s) = pbc(ypp(s),Ly)
                        chain_pbc(3,s) = pbc(zp(s),Lz)
                    else
                        chain_pbc(1,s) = xpp(s) 
                        chain_pbc(2,s) = ypp(s)
                        chain_pbc(3,s) = zp(s)
                    endif        

                    ! .. transforming form real- to lattice coordinates                 
                    xi = int(chain_pbc(1,s)/delta)+1
                    yi = int(chain_pbc(2,s)/delta)+1
                    zi = int(chain_pbc(3,s)/delta)+1
                        
                    if(isOutsideLattice(xi,yi,zi,nx,ny,nz)) then   
                        text="Conformation outside box:"
                        call print_to_log(LogUnit,text)  
                        print*,text                          
                        print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        info= myio_err_index
                        return
                    endif 


                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                    
                    indexconf(s,conf)%elem(1) = idx ! CA element    


                    ! apply elements other than CA
                    do j=2,nelem(s) 

                        chain_tmp(1) = chain_elem_index(1,s)%elem(j)+xcm  
                        chain_tmp(2) = chain_elem_index(2,s)%elem(j)+ycm
                        chain_tmp(3) = chain_elem_index(3,s)%elem(j)+zcm  

                        if(pbc_chains) then 
                            chain_pbc_tmp(1) = pbc(ut(chain_tmp(1),chain_tmp(2)),Lx) 
                            chain_pbc_tmp(2) = pbc(vt(chain_tmp(1),chain_tmp(2)),Ly)
                            chain_pbc_tmp(3) = pbc(chain_tmp(3) ,Lz)
                        else
                            chain_pbc_tmp(1) = ut(chain_tmp(1),chain_tmp(2))
                            chain_pbc_tmp(2) = vt(chain_tmp(1),chain_tmp(2))
                            chain_pbc_tmp(3) = chain_tmp(3)
                        endif

                        xi = int(chain_pbc_tmp(1)/delta)+1
                        yi = int(chain_pbc_tmp(2)/delta)+1
                        zi = int(chain_pbc_tmp(3)/delta)+1
                        print*,isOutsideLattice(xi,yi,zi,nx,ny,nz)

                        if(isOutsideLattice(xi,yi,zi,nx,ny,nz)) then   

                            text="Conformation outside box:"
                            call print_to_log(LogUnit,text)  
                            print*,text  
                            print*,"chain_elem_index= ",(chain_elem_index(k,s)%elem(j),k=1,3)
                            print*,"chain_pbc= ",chain_pbc_tmp(:),"s= ",s," j= ",j                          
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                            info= myio_err_index
                            return
                        endif

                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                        
                        indexconf(s,conf)%elem(j) = idx ! all other element

                    enddo 
                    
                enddo

                ! energyLJ    = GBenergyeffective(chain_pbc,nnucl,no_overlap)  

                call make_com_nucl_rotation(chain_rot,nnucl,unitvector_triplets,rcom)
                
                call maximum_xnucl(conf,isVolfracLargerOne)

                print*,"conf=",conf,"isVolfracLargerOne=",isVolfracLargerOne
                
                if(isVolfracLargerOne) no_overlap=.false. 

                energychainLJ(conf)    = energyLJ
                energychain(conf)      = energy
                no_overlapchain(conf)  = no_overlap 
                Rgsqr(conf)            = radius_gyration_com_rotation(rcom,nnucl)
                Rendsqr(conf)          = end_to_end_distance_com_rotation(rcom,nnucl)
                bond_angle(:,conf)     = bond_angles_com_rotation(rcom,nnucl)
                dihedral_angle(:,conf) = dihedral_angles_com_rotation(rcom,nnucl)
                nucl_spacing(:,conf)   = nucleosomal_spacing_com_rotation(rcom,nnucl)
                gyr_tensor             = gyr_tensor_com_rotation(rcom,nnucl)
                Asphparam(conf)        = Asphericty_parameter(Rgsqr(conf),gyr_tensor)
          
                ! export make_spheriod

                if(COMOLD) then 
                    Rgsqr(conf)            = radius_gyration_com(chain_pbc,nnucl,segcm)
                    Rendsqr(conf)          = end_to_end_distance_com(chain_pbc,nnucl,segcm)
                    bond_angle(:,conf)     = bond_angles_com(chain_pbc,nnucl,segcm)
                    dihedral_angle(:,conf) = dihedral_angles_com(chain_pbc,nnucl,segcm)
                    nucl_spacing(:,conf)   = nucleosomal_spacing_com(chain_pbc,nnucl,segcm)
                    gyr_tensor             = gyr_tensor_com(chain_pbc,nnucl,segcm)
                    Asphparam(conf)        = Asphericty_parameter(Rgsqr(conf),gyr_tensor)
                endif    

                conf=conf+1   

            case default

                text="Error: in make_chains_XYZ_nucl geometry not cubic or prism: stopping program"
                call print_to_log(LogUnit,text)
                info = myio_err_geometry
                return 
                    
            end select

        endif   ! read 

    enddo       ! end while loop                                                                                                          
    ! end chains generation    
    
    conf=conf-1  ! lower by one  

    if(conf<max_confor) then
        print*,"subroutine read_chains_xyz_nucl_volume :" 
        print*,"conf     = ",conf," less then imposed max cuantas     = ",max_confor
        print*,"conffile = ",conffile
        cuantas=conf   
        info = myio_err_conf        
    else
        text="Chains generated: subroutine read_chains_xyz_nucl_volume"
        call print_to_log(LogUnit,text)
        readinchains=conffile
        info = 0
    endif


    write(istr,'(L2)')isVdWintEne
    text="isVdWintEne = "//trim(adjustl(istr))
    call print_to_log(LogUnit,text)

    if(.not.(isChainEnergyFile)) energychain=0.0_dp

    close(un) 
    if(isChainEnergyFile) close(un_ene)

    cuantas_no_overlap = compute_cuantas_no_overlap(cuantas,no_overlapchain)

    call normed_weightchains()
    text="No reweighting chains."
    call print_to_log(LogUnit,text)
    logweightchain=0.0_dp        

    deallocate(energychain) ! free unused variables 

    if(DEBUG) call write_indexconf_lammps_trj(info_traj)

    if(runtype=="rangeVdWeps") then 
        allocate(energychainLJ0(cuantas))
        do i=1,cuantas
            energychainLJ0(i)=energychainLJ(i)
        enddo 
    else 
        ! energy of conformation fixed: VdWscale%val set the energy scale relative to epsilonS
        ! for runtype=="rangVdWeps" loop over values of VdWscale%val     
        do i=1,cuantas
            energychainLJ(i)=VdWscale%val*energychainLJ(i)
        enddo
    endif            

end subroutine read_chains_xyz_nucl_volume

subroutine read_graftpts_xyz_nucl(info)

    !use mpivars, only : rank
    use parameters, only : unit_conv
    use myio, only : myio_err_chainsfile, myio_err_graft
    use myutils,  only : newunit
    use chains, only : sgraftpts

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    character(len=25) :: fname
    integer :: ios 
    real(dp) :: xc,yc,zc          
    integer :: ix,iy,iz,un,s, i,t
    integer :: rankfile
    integer :: item,moltype,nsegfile,idatom
    character(len=30) :: istr,str
    !real(dp) :: xbox0,xbox1
    real(dp) :: scalefactor
    logical :: exist, isGraftItem

    ! .. executable statements 
    info=0

    !. . open file    
    rankfile=0                                                                                           
    write(istr,'(I4)')rankfile
    fname='traj-graft.'//trim(adjustl(istr))//'.xyz'
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj-graft.rank.xyz file does not exit'
        info = myio_err_chainsfile
        return
    end if
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    end if

    ! read preamble/header
   
    scalefactor=unit_conv
    isGraftItem=.false.
    do s=1,2
        read(un,*,iostat=ios)item,xc,yc,zc
        if(item==sgraftpts(1)) then
            t=1
            isGraftItem=.true.
        else
            t=2
        end if    
        xgraftloop(1,t)=xc*scalefactor
        xgraftloop(2,t)=yc*scalefactor
        xgraftloop(3,t)=zc*scalefactor    
    end do
    
    close(un)

    if(.not.isGraftItem) then
        print*,'Error read_graftpts_xyz_loop: sgraft not in traj file.'
        info = myio_err_graft
    end if
        
end subroutine


! post:  make isAmonomer, type_of_monomer and  type_of_monomer_char 
! pre: nseg, nsegtypes, freq=chainperiod and chaintype are set

subroutine make_sequence_chain(freq,chaintype)
  
    use globals, only : nseg, nsegtypes
    use chains, only :  isAmonomer,type_of_monomer_char, type_of_monomer,ismonomer_of_type
    use parameters, only : typesfname

    integer, intent(in) :: freq
    character(len=8),intent(in)  :: chaintype

    !     .. local variables
    integer :: info
    integer :: s

    select case (chaintype ) 
    case ('altA')
        do s=1,nseg
            if(mod(s,freq).ne.0) then ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo
    case('altB') 
        do s=1,nseg
            if(mod(s,freq).eq.0) then ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo
    case('diblockA')
        do s=1,nseg
            if(s<=freq) then   ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo
    case('diblockB') 
        do s=1,nseg
            if(s>=freq) then   ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo  
    case('copolyAB')

        call read_sequence_copoly_from_file(info)

         do s=1,nseg
           if(isAmonomer(s)) then 
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            end if
        enddo

    case('multi')
    
        call read_type_of_monomer(type_of_monomer,type_of_monomer_char,typesfname, nseg) 
        call make_type_table(ismonomer_of_type,type_of_monomer,nseg,nsegtypes)    
        
        do s=1,nseg
            isAmonomer(s)=(type_of_monomer_char(s)=="AA")
        end do

    case default
        print*,"Wrong chaintype: aborting program"
        stop
    end select

end subroutine make_sequence_chain


subroutine read_sequence_copoly_from_file(info)

    use globals, only : nseg
    use chains,  only : isAmonomer
    use myutils, only : newunit

    integer, intent(out),optional :: info

    character(len=11) :: fname
    integer :: ios, un_seq, s
    character :: char

    if (present(info)) info = 0

    write(fname,'(A11)')'sequence.in'
    open(unit=newunit(un_seq),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening file sequence.in : iostat =', ios
        stop
    end if
        
    s=0
    ios=0

    do while (s<nseg.and.ios==0)

        read(un_seq,*,iostat=ios)char

        if (ios==0) then
            s=s+1 
            if(char=="A") then   ! A segment
                isAmonomer(s)=.TRUE.
            else
                isAmonomer(s)=.FALSE.
            end if
        end if    
    end do 
    if(s/=nseg) then 
        print*,"reached end of file before all elements read"
        info = 1
        stop "read sequence file failed"
    end if

end subroutine read_sequence_copoly_from_file


! Compute weight chain w=e^E/Tre^E and normalize
! layout conformations : there are  nset of conformations thus total nodes = nset 
! Each set hold different conformations

subroutine normed_weightchains()

    !use mpivars
    use globals, only : cuantas
    use chains, only : energychain, logweightchain
   
    integer :: un, c, k
    real(dp) :: localsum, totalsum, logtotalsum

        
    localsum=0.0_dp    

    do c=1,cuantas
        localsum=localsum+exp(energychain(c))   
    enddo    
  
    ! call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
    ! call MPI_ALLREDUCE(localsum, totalsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr) 
    totalsum=localsum
    ! normalize
    logtotalsum=log(totalsum)
    do c=1,cuantas
        logweightchain(c)=energychain(c)-logtotalsum
    enddo    

end subroutine



function minimum_chainenergy() result(min_chainenergy)

   ! use  mpivars
    use  globals, only :cuantas
    use  chains, only : energychain

    real(dp) :: min_chainenergy

    integer :: conf
    
    min_chainenergy =0.0_dp
    do conf=1,cuantas      
        if(min_chainenergy > energychain(conf)) min_chainenergy=energychain(conf)
    enddo
   
end function


subroutine global_minimum_chainenergy()

    !use  mpivars
    use  globals, only : cuantas
    use  chains, only : energychain, energychain_min
    use  parameters, only: isEnergyShift

    real(dp) :: localmin(2), globalmin(2)
    integer :: i

    localmin(1)=minimum_chainenergy()
    localmin(2)=0 !rank   
    print*,"rank ",0 ," has local lowest value of", localmin(1)
    
!    call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
!    call MPI_ALLREDUCE(localmin, globalmin, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
    
    !if (rank == 0) then  
        print*,"0 ",globalmin(2)," has lowest value of", globalmin(1)  
   ! endif
    
    energychain_min =globalmin(1)

    if(isEnergyShift) then        
        call shift_energy_chain(globalmin(1))
        print*,"Shift minimum by :",energychain_min 
    endif
end subroutine



! post: substract minimum internal energy from all internal chain energies
subroutine shift_energy_chain(min_chainenergy)
    
    use  globals, only :cuantas
    use  chains, only : energychain
    
    real(dp), intent(in) :: min_chainenergy

    ! local
    integer :: conf

    do conf=1,cuantas        
        energychain(conf)=energychain(conf)-min_chainenergy
    enddo

end subroutine shift_energy_chain


subroutine set_properties_chain(freq,chaintype)
    
    use globals
    use chains
    use parameters, only : lseg

    integer, intent(in) :: freq
    character(len=8), intent(in)   :: chaintype

    call set_isHomopolymer(freq,chaintype)
    call set_lsegAA()

    if(isHomopolymer) lseg=set_lseg_Homopolymer()

end subroutine set_properties_chain

! set logical isHomoolymer 
! pre : type_of_monomer need to be initialized

subroutine set_isHomopolymer(freq,chaintype)

    use globals, only : nseg
    use chains, only : type_of_monomer, isHomopolymer

    integer, intent(in) :: freq
    character(len=8), intent(in)   :: chaintype
    
    integer :: s,  type_number
    logical :: flag

    if(chaintype=="multi") then
        type_number=type_of_monomer(1)
        flag=.true.
        do s=2,nseg
            if((type_of_monomer(s)/=type_number)) flag=.false.
        enddo
        if(flag) then 
          isHomopolymer=.True.
        else
          isHomopolymer=.False.
        endif 
                  
    else if(freq>nseg) then 
        isHomopolymer=.TRUE.
    else 
        if(freq==0.and.(chaintype.eq."diblockA".or. chaintype=="diblockB")) then 
            isHomopolymer=.TRUE.
        else
            isHomopolymer=.FALSE.
        endif      
    endif    


end subroutine set_isHomopolymer


!  Assigns specific values to lsegAA initial set in init_lseg in parameters

subroutine set_lsegAA

    use globals
    use chains
    use parameters, only : lseg, lsegAA,lsegPAA, lsegPAMPS, lsegPEG
    use parameters, only : chainmethod 

    if(chainmethod=='MC') then  ! chain are not read in from file 
        
        select case (systype)
        case ("elect","electVdWAB","electdouble") ! diblock copolymer lseg determined in cadenas_sequence      
            ! assume A is present and type_of_monomer =1 
            lsegAA(1) = lsegPAA
            lsegAA(2) = lsegPAMPS        
        case ("electA")                 ! homopolymer weak polyacid VdW     
            lsegAA = lsegPAA
        case ("neutral","neutralnoVdW")                ! homopolymer neutral
            !lsegAA = lsegPEG
        case ("brush_mul","brush_mulnoVdW","brush","brush_neq","brushvarelec","brushborn","brushdna")
            ! lsegAA = lsegPAA !0.36287_dp
        case default
            print*,"Error: in set_lsegAA, systype=",systype
            print*,"stopping program"
            stop
        end select  
    endif   

end subroutine set_lsegAA

     
function set_lseg_Homopolymer()result(lengthseg)

    use parameters, only : lsegAA

    real(dp) :: lengthseg
    
    lengthseg=lsegAA(1)
    
end function set_lseg_Homopolymer



! .. compute periodic boundary condition in z direction
! .. maps z onto interval [0,deltaz]

function pbc(z,deltaz) result(zpbc)
    implicit none 
    real(dp), intent(in) :: deltaz
    real(dp), intent(in) :: z
    real(dp) :: zpbc 
    if(z>0) then 
        zpbc=z-int(z/deltaz)*deltaz
    else
        zpbc=z-(int(z/deltaz)-1)*deltaz
    endif   
end function


! function ipbc(ival,imax) result(intpbc)
!     implicit none 
!     integer, intent(in) :: ival
!     integer, intent(in) :: imax
!     integer :: intpbc

!     if(ival>0) then
!         intpbc=ival-int((ival-1)/imax)*imax
!     else
!         intpbc=ival-(int((ival-1)/imax)-1)*imax
!     endif

! end function


!     .. returns number of A monomer
!     .. this works only for binary AB polymers

function number_Amonomers(nseg)result(npolA)
  
    use chains, only : isAmonomer

    integer , intent(in) :: nseg
    integer :: npolA

    !     .. local variables
    integer :: s

    npolA=0
    do s=1,nseg
        if(isAmonomer(s).eqv..true.) npolA=npolA+1
    enddo

end function number_Amonomers



!  .. assign type_of_monomer and type_of_monomer_char from values in file named filename

subroutine read_type_of_monomer(type_of_monomer, type_of_monomer_char,filename, nseg)

    !use mpivars
    use myutils, only : newunit
    use parameters, only : lenfname
    
    !     .. arguments 
    integer, intent(inout) :: type_of_monomer(:)
    character(len=3), intent(inout) :: type_of_monomer_char(:)
    character(lenfname), intent(in) :: filename
    integer,  intent(in) :: nseg

    !      .. local variables
    integer :: ios, un  ! un = unit number
    integer :: s
    character(80) :: istr,str,letter

    !     .. reading in of variables from file
    open(unit=newunit(un),file=filename,iostat=ios,status='old')
    if(ios/=0 ) then
        write(istr,'(I2)')ios
        str='Error opening file '//trim(adjustl(filename))//' : iostat = '//istr
        print*,str
!        call MPI_FINALIZE(ierr)
        stop
    endif
    
    s=0
    ios=0

    do while (s<nseg.and.ios==0)
        s=s+1
        read(un,*,iostat=ios)type_of_monomer(s),type_of_monomer_char(s)
    enddo

    if(s/=nseg.or.ios/=0) then 
        str="reached end of file before all elements read or ios error"
        print*,str
        str="read file "//trim(adjustl(filename))//" failed"
        print*,str
!        call MPI_FINALIZE(ierr)
        stop
    endif

    close(un)

end subroutine read_type_of_monomer



! ismonomer_of_type is a table which row index is the segment  and column index correspond to the segment type  
! pre: type_of_monomer needs to be initialized
! post: table list of logicals indication is segement s is of type t

subroutine make_type_table(ismonomer_of_type,type_of_monomer,nseg,nsegtypes)
    
    logical, intent(out) :: ismonomer_of_type(:,:)
    integer, intent(in) :: type_of_monomer(:) 
    integer, intent(in) :: nseg
    integer, intent(in) :: nsegtypes

    ! local variable 
    integer :: s,t 

    do s=1,nseg 
        do t=1,nsegtypes   
            ismonomer_of_type(s,t)=.false.
        enddo
        ismonomer_of_type(s,type_of_monomer(s))=.true.    
    enddo

end subroutine make_type_table


! routine determines is segment type t is chargeable
! pre: zpol needs to be initialized
! post: ismonomer_chargeable list of logicals

subroutine make_charge_table(ismonomer_chargeable,zpol,nsegtypes)
 
    logical, intent(out) :: ismonomer_chargeable(:)
    integer, intent(in) :: zpol(:,:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t 

    do t=1,nsegtypes
        if(zpol(t,1)==0.and.zpol(t,2)==0) then 
            ismonomer_chargeable(t)=.false.
        else
            ismonomer_chargeable(t)=.true.
        endif
    enddo

end subroutine make_charge_table


! routine determines is type of charge of segment type t 
! pre: zpol needs to be initialized
! post: type_of_charge list of characters

subroutine make_type_of_charge_table(type_of_charge,zpol,nsegtypes)
 
    character(len=1), intent(out) :: type_of_charge(:)
    integer, intent(in) :: zpol(:,:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t 

    do t=1,nsegtypes
        if(zpol(t,1)==0.and.zpol(t,2)==0) then 
            type_of_charge(t)="N"
        else if(zpol(t,1)==0.and.zpol(t,2)==-1) then
             type_of_charge(t)="A"
        else if(zpol(t,1)==1.and.zpol(t,2)==0) then
             type_of_charge(t)="B"  
        endif
    enddo

end subroutine make_type_of_charge_table

! Checks if (xi,yi,zi) is inside lattice
! 0<=xi<=nx  and 0<= yi <= ny and o
! returns: logical = .true. if outside 

function isOutsideLattice(xi,yi,zi,nx,ny,nz)result(isOutside)
 
    integer, intent(in) :: xi,yi,zi
    integer, intent(in) :: nx,ny,nz

    logical :: isOutside

    integer :: intxi,intyi,intzi
    logical :: isInside
    
    intxi=int(xi/nx)
    intyi=int(yi/ny)
    intzi=int(zi/nz)

    isInside=(0<=intxi).and.(intxi<=1).and.(0<=intyi).and.(intyi<=1).and.(0<=intzi).and.(intzi<=1)
    isOutside=.not.isInside
    
end function isOutsideLattice


logical function is_polymer_neutral(ismonomer_chargeable, nsegtypes)
 
    logical, intent(in) :: ismonomer_chargeable(:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t,flag 

    flag=0
    do t=1,nsegtypes
        if(.not.ismonomer_chargeable(t)) flag=flag+1
    enddo

    if(flag==nsegtypes) then 
        is_polymer_neutral=.true. 
    else
        is_polymer_neutral=.false. 
    endif

end function is_polymer_neutral


! Make a lammps trajectory file based on indexchain
! The lammps trajectory uses a molecular ATOM Style
! Unit length in Angstrom

subroutine write_indexchain_lammps_trj(info)

    !use mpivars, only : rank
    use globals, only : nseg, cuantas
    use myutils, only : newunit, lenText
    use volume, only : delta, nz, nx, ny, indextocoord
    use chains, only : indexchain, type_of_monomer
    use myio, only   : myio_err_chainsfile

    integer, optional, intent(inout) :: info

    character(len=lenText) :: text, istr
    character(len=25) :: fname
    integer :: ios, un_trj 
    real(dp):: x, y, z
    integer :: ix, iy, iz, ic(3), idx
    real(dp) :: xbox0, xbox1(3)
    integer :: idatom, item, moltype, conf

    if (present(info)) info=0

    ! open file 
    write(istr,'(I4)')0
    fname='traj_indexchain_'//trim(adjustl(istr))//'.lammpstrj'
    open(unit=newunit(un_trj),file=fname,status='new',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

    xbox0=0.0_dp
    xbox1(1)=nx*delta*10.0_dp
    xbox1(2)=ny*delta*10.0_dp
    xbox1(3)=nz*delta*10.0_dp
    
    do conf=1,cuantas
        ! write preamble 
        write(un_trj,*)'ITEM: TIMESTEP' 
        write(un_trj,*)conf 
        write(un_trj,*)'ITEM: NUMBER OF ATOMS' 
        write(un_trj,*)nseg 
        write(un_trj,*)'ITEM: BOX BOUNDS ff ff ff'
        write(un_trj,*)xbox0,xbox1(1)
        write(un_trj,*)xbox0,xbox1(2)
        write(un_trj,*)xbox0,xbox1(3)
        write(un_trj,*)'ITEM: ATOMS id mol type x y z '
        !  determine x, y, z,
        do item=1,nseg
            idx=indexchain(item,conf)
            ic=indextocoord(idx,:)
            x=(ic(1)-0.5_dp)*delta*10.0_dp
            y=(ic(2)-0.5_dp)*delta*10.0_dp
            z=(ic(3)-0.5_dp)*delta*10.0_dp
            idatom=1
            moltype=type_of_monomer(item) 
            write(un_trj,*)item,idatom,moltype,x,y,z
        enddo    
    enddo
        
    close(un_trj)    


end subroutine write_indexchain_lammps_trj



! Make a lammps trajectory file based on indexchain
! the lammps trajectory uses a molecular ATOM Style
! Unit length in Angstrom

subroutine write_indexconf_lammps_trj(info)

    !use mpivars, only : rank
    use globals, only : nseg, cuantas
    use myutils, only : newunit, lenText
    use volume, only : delta, nx, ny, nz, coordinateFromLinearIndex
    use chains, only : indexconf, type_of_monomer, nucl_elem_type, nelem
    use myio, only   : myio_err_chainsfile

    integer, optional, intent(inout) :: info

    character(len=lenText) :: text, istr
    character(len=25) :: fname
    integer :: ios, un_trj 
    real(dp):: x, y, z
    integer :: ix, iy, iz, i, j, k, idx, s, em
    real(dp) :: xbox0,xbox1(3)
    integer :: idatom, item, conf
    character(len=3) :: moltype
    integer :: nsegtot

    if (present(info)) info=0

    ! open file 
    write(istr,'(I4)')0
    fname='traj_indexconf'//trim(adjustl(istr))//'.lammpstrj'
    open(unit=newunit(un_trj),file=fname,status='new',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

    xbox0=0.0_dp

    xbox1(1)=nx*delta*10.0_dp ! conversion to Angstrom
    xbox1(2)=ny*delta*10.0_dp 
    xbox1(3)=nz*delta*10.0_dp 

    nsegtot=0
    do s=1,nseg
        nsegtot=nsegtot+nelem(s)
    enddo   


    do conf=1,cuantas
        ! write preamble 
        write(un_trj,'(14A)')'ITEM: TIMESTEP' 
        write(un_trj,*)conf 
        write(un_trj,'(21A)')'ITEM: NUMBER OF ATOMS' 
        write(un_trj,*)nsegtot 
        write(un_trj,'(25A)')'ITEM: BOX BOUNDS ff ff ff'
        write(un_trj,*)xbox0,xbox1(1)
        write(un_trj,*)xbox0,xbox1(2)
        write(un_trj,*)xbox0,xbox1(3)
        write(un_trj,'(29A)')'ITEM: ATOMS id mol type x y z '
        !  determine x, y, z,

        do s=1,nseg
            do em=1,nelem(s)
                idx = indexconf(s,conf)%elem(em)
                
                call coordinateFromLinearIndex(idx, i, j, k)
                x=(i-0.5_dp)*delta*10.0_dp
                y=(j-0.5_dp)*delta*10.0_dp
                z=(k-0.5_dp)*delta*10.0_dp
                
                item=s
                idatom=type_of_monomer(s)
                moltype=nucl_elem_type(em,type_of_monomer(s))  
                write(un_trj,*)item,idatom,moltype,x,y,z

            enddo 
        enddo       
    enddo
        
    close(un_trj)    

end subroutine write_indexconf_lammps_trj

! Creates and opens file for lammps trajectory film
! Writing do by write_chain_lammps_trj
! input:        integer : info : value 0 no errors
! output/return integer : un_trj : unit nubmer of file

function open_chain_lammps_trj(info,text)result(un_trj)

    !use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    integer, intent(inout) :: info
    character(len=*), intent(in),optional :: text ! extra filename descriptor 

    integer :: un_trj
    ! local
    character(len=lenText) :: istr
    character(len=30) :: fname
    integer :: ios
    logical :: exist
    
    info = 0    

    if(present(text)) then
        fname='trajchain_'//trim(adjustl(text))//'_'
    else
        fname='trajchain_'
    endif

    write(istr,'(I4)')0
    fname=trim(fname)//trim(adjustl(istr))//'.lammpstrj'
    inquire(file=fname,exist=exist)
    if(.not.exist) then
        open(unit=newunit(un_trj),file=fname,status='new',action="write",iostat=ios)
        if(ios > 0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            info = myio_err_chainsfile
            return
        endif
    else
        open(unit=newunit(un_trj),file=fname,position="append",status='old',action="write",iostat=ios)
        if(ios > 0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            info = myio_err_chainsfile
            return
        endif
    endif   

end function



! Creates and opens file for lammps trajectory film
! Writing done by write_chain_elem_index_lammps_trj
! input:        integer : info : value 0 no errors
! output/return integer : un_trj : unit nubmer of file

function open_chain_elem_index_lammps_trj(info,text)result(un_trj)

    !use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    integer, intent(inout) :: info
    character(len=*), intent(in),optional :: text ! extra filename descriptor 

    integer :: un_trj

    ! local
    character(len=lenText) :: istr
    character(len=40) :: fname
    integer :: ios
    logical :: exist
    
    info = 0    

    if(present(text)) then
        fname='trajchainelemindex_'//trim(adjustl(text))//'_'
    else
        fname='trajchainelemindex_'
    endif

    write(istr,'(I4)')0
    fname=trim(fname)//trim(adjustl(istr))//'.lammpstrj'
    inquire(file=fname,exist=exist)
    if(.not.exist) then
        open(unit=newunit(un_trj),file=fname,status='new',action="write",iostat=ios)
        if(ios > 0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            info = myio_err_chainsfile
            return
        endif
    else
        open(unit=newunit(un_trj),file=fname,position="append",status='old',action="write", iostat=ios)
        if(ios > 0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            info = myio_err_chainsfile
            return
        endif
    endif   

end function

! Creates and opens file for lammps xyz trajectory film
! Writing do by write_chain_xyz
! input:        integer : info : value 0 no errors
! output/return integer : un_trj : unit nubmer of file

function open_chain_xyz(info,text)result(un_trj)

    !use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    integer, intent(inout) :: info
    character(len=*), intent(in),optional :: text ! extra filename descriptor 

    integer :: un_trj
    ! local
    character(len=lenText) :: istr
    character(len=30) :: fname
    integer :: ios
    logical :: exist
    
    info = 0    

    if(present(text)) then
        fname='trajchain_'//trim(adjustl(text))//'_'
    else
        fname='trajchain_'
    endif

    write(istr,'(I4)')0
    fname=trim(fname)//trim(adjustl(istr))//'.xyz'
    inquire(file=fname,exist=exist)
    if(.not.exist) then
        open(unit=newunit(un_trj),file=fname,status='new',action="write",iostat=ios)
        if(ios > 0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            info = myio_err_chainsfile
            return
        endif
    else
        open(unit=newunit(un_trj),file=fname,position="append",status='old',action="write",iostat=ios)
        if(ios > 0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            info = myio_err_chainsfile
            return
        endif
    endif   

end function


! Writes the coordinates of chain object to lammps a trajectory file.
! Unit length in Angstrom
! pre : file opened with open_chain_lammps_trj
! input: integer un_trj : unit number of file to write to.
!        real(dp) chain(:,:,:) : index 1= coordinate 2= segment/atom numnber 3 number of dublicates conformations 
!        integer nchains : number of confomations
! return: a lammpstrj file

subroutine write_chain_xyz(un_trj,chain,nchains)

    use myutils, only : newunit, lenText
    use volume, only : delta, nx, ny, nz
    use chains, only : type_of_monomer

    real(dp), intent(in) :: chain(:,:,:)
    integer, intent(in) :: nchains
    integer , intent(in) :: un_trj
    
   
    real(dp) :: x, y, z
    integer ::  i, j, item
    integer :: ndim ! local  nseg

    ndim = size(chain,2) 

    do j=1,nchains   
        
        conf_write_chain=conf_write_chain+1 
        ! write preamble 
        write(un_trj,*)ndim
        write(un_trj,'(5A)')'ATOMS'
        do item=1,ndim
            x = chain(1,item,j)
            y = chain(2,item,j)
            z = chain(3,item,j)
            write(un_trj,*)x,y,z
        enddo    
    enddo
        
end subroutine write_chain_xyz

! Writes the coordinates of chain object to lammps a trajectory file.
! Unit length in Angstrom
! pre : file opened with open_chain_lammps_trj
! input: integer un_trj : unit number of file to write to.
!        real(dp) chain(:,:,:) : index 1= coordinate 2= segment/atom numnber 3 number of dublicates conformations 
!        integer nchains : number of confomations
! return: a lammpstrj file

subroutine write_chain_lammps_trj(un_trj,chain,nchains, type_of_monomer_loc)

    use myutils, only : newunit, lenText
    use volume, only : delta, nx, ny, nz
    use chains, only : type_of_monomer_default => type_of_monomer, mapping_num_to_char

    real(dp), intent(in) :: chain(:,:,:)
    integer, intent(in) :: nchains
    integer, intent(in) :: un_trj
    integer, intent(in), optional :: type_of_monomer_loc(:)
    
    
    character(len=lenText) :: istr
    real(dp) :: x, y, z
    integer ::  i, j, k
    real(dp) :: xbox0, xbox1(3)
    integer :: idatom, item!, conf_write
    integer :: nseg ! local  nseg
    integer, allocatable :: type_of_monomer(:)
    character(len=3) :: moltype

    nseg = size(chain,2) 

   ! allocate(type_of_monomer)

    if(present(type_of_monomer_loc)) then
        type_of_monomer = type_of_monomer_loc
    else
        type_of_monomer = type_of_monomer_default 
    endif

    ! Check consistency
    if (size(type_of_monomer) /= nseg) then
         error stop "write_chain_lammps_trj: inconsistent monomer type array size"
    endif
    if (size(chain,1) /= 3) error stop "Invalid chain coordinate dimension"
    if (size(chain,3) < nchains) error stop "Invalid number of chains"


    xbox0=0.0_dp
    xbox1(1)=nx*delta*10.0_dp ! conversion to Angstrom 
    xbox1(2)=ny*delta*10.0_dp
    xbox1(3)=nz*delta*10.0_dp

    do j=1,nchains   
        
        conf_write_chain=conf_write_chain+1 
        ! write preamble 
        write(un_trj,'(14A)')'ITEM: TIMESTEP' 
        write(un_trj,*)conf_write 
        write(un_trj,'(21A)')'ITEM: NUMBER OF ATOMS' 
        write(un_trj,*)nseg
        write(un_trj,'(25A)')'ITEM: BOX BOUNDS ff ff ff'
        write(un_trj,*)xbox0,xbox1(1)
        write(un_trj,*)xbox0,xbox1(2)
        write(un_trj,*)xbox0,xbox1(3)
        write(un_trj,'(29A)')'ITEM: ATOMS id mol type x y z'
        do item=1,nseg
            x = chain(1,item,j)*10.0_dp
            y = chain(2,item,j)*10.0_dp
            z = chain(3,item,j)*10.0_dp
            idatom=type_of_monomer(item)
            moltype= mapping_num_to_char(type_of_monomer(item) )
            write(un_trj,*)item,idatom,moltype,x,y,z
        enddo   

        ! do item=1,nseg
        !     x = chain(1,item,j)*10.0_dp
        !     y = chain(2,item,j)*10.0_dp
        !     z = chain(3,item,j)*10.0_dp
        !     idatom=1
        !     moltype=type_of_monomer(item) 
        !     write(un_trj,*)item,idatom,moltype,x,y,z
        ! enddo    

    enddo
        
end subroutine write_chain_lammps_trj

! Writes the coordinates of chain_elem_index object to a lammps trajectory file.
! chain_elem_index(k,s)%elem(j) the k=(1,2,3)=(x,yz) cartesian coordinate of jth chain element of sth segement (AA and DNA) 
!    type(var_darray), dimension(:,:), allocatable :: chain_elem_index
! Unit length in Angstrom
! pre : file opened with open_chain_lammps_trj
! input: integer un_trj : unit number of file to write to.
!        chain_elem_index (:,:,:) : index 1= coordinate 2= segment/atom number 3 number of dublicates conformations 
!    
! return: a lammpstrj file

subroutine write_chain_elem_index_lammps_trj(un_trj,chain_elem_index)

    use globals, only : nseg
    use myutils, only : newunit, lenText
    use volume, only : delta, nx, ny, nz, coordinateFromLinearIndex
    use chains, only : type_of_monomer, nelem,  nucl_elem_type, type_of_monomer_char
    use chains, only : var_darray ! type def  

    type(var_darray), allocatable, dimension(:,:), intent(in) :: chain_elem_index
    integer, intent(in) :: un_trj
    
    character(len=lenText) :: istr
    real(dp) :: x, y, z, xcm, ycm, zcm
    integer ::  j, s
    real(dp) :: xbox0, xbox1(3)
    integer :: idatom, item, conf, nsegtot
    character(len=3) :: moltype,idatom_char

    xbox0=0.0_dp 
    xbox1(1)=nx*delta*10.0_dp
    xbox1(2)=ny*delta*10.0_dp
    xbox1(3)=nz*delta*10.0_dp
    
    xcm= nx*delta*10.0_dp/2.0_dp           
    ycm= ny*delta*10.0_dp/2.0_dp          ! center box
    zcm= nz*delta*10.0_dp/2.0_dp

    item=0
    conf_write=conf_write+1

    nsegtot=0
    do s=1,nseg
        nsegtot=nsegtot+nelem(s)
    enddo    

    ! write preamble 
    write(un_trj,'(14A)')'ITEM: TIMESTEP' 
    write(un_trj,*)conf_write 
    write(un_trj,'(21A)')'ITEM: NUMBER OF ATOMS' 
    write(un_trj,*)nsegtot 
    write(un_trj,'(25A)')'ITEM: BOX BOUNDS ff ff ff'
    write(un_trj,*)xbox0,xbox1(1)
    write(un_trj,*)xbox0,xbox1(2)
    write(un_trj,*)xbox0,xbox1(3)
    write(un_trj,'(29A)')'ITEM: ATOMS id mol type x y z'
    do s=1,nseg
        do j=1,nelem(s)
            x = chain_elem_index(1,s)%elem(j)*10.0_dp !+xcm
            y = chain_elem_index(2,s)%elem(j)*10.0_dp !+ycm
            z = chain_elem_index(3,s)%elem(j)*10.0_dp !+zcm
        
            idatom=type_of_monomer(s)
            moltype=nucl_elem_type(j,type_of_monomer(s))  
            
            item=item+1
            write(un_trj,*)item,idatom," ",moltype,x,y,z
        enddo    
    enddo    
    
end subroutine write_chain_elem_index_lammps_trj


function open_chain_energy(info)result(un_ene)

   ! use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    integer, optional, intent(inout) :: info

    integer :: un_ene
    ! local
    character(len=lenText) :: istr
    character(len=25) :: fname
    integer :: ios

    write(istr,'(I4)')0
    fname='traj.'//trim(adjustl(istr))//'.ene'
    open(unit=newunit(un_ene),file=fname,status='replace',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

end function

! .. Compute radius of gyration of a sub chain conformation
! .. sub chain confomation defined by sequence segcm = segment number 
! .. denoting center mass of  the nmer histone of the nucleosome nmer-array
! .. use rcom estalished with rotation alogrithme

function radius_gyration_com_rotation(rcom,nmer) result(Rgsqr)

    real(dp), intent(in) :: rcom(:,:) ! dimension : (3,nnucl)
    integer, intent(in) :: nmer

    real(dp) :: Rgsqr
    integer :: i,j,k

    Rgsqr=0.0_dp
    
    do i=1,nmer
        do j=1,nmer
            do k=1,3
                Rgsqr = Rgsqr+(rcom(k,i)-rcom(k,j))**2
            enddo
        enddo
    enddo     

    Rgsqr=Rgsqr/(2.0_dp*nmer*nmer)

end function radius_gyration_com_rotation


! .. Compute radius of gyration of a sub chain conformation
! .. sub chain confomation defined by sequence segcm = segment number 
! .. denoting center mass of  the nmer histone of the nucleosome nmer-array

function radius_gyration_com(chain,nmer,segcm) result(Rgsqr)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcm(:) 

    real(dp) :: Rgsqr
    integer :: i,j,k,isegcm,jsegcm

    Rgsqr=0.0_dp
    
    do i=1,nmer
        isegcm=segcm(i)
        do j=1,nmer
            jsegcm=segcm(j)
            do k=1,3
                Rgsqr=Rgsqr+(chain(k,isegcm)-chain(k,jsegcm))**2
            enddo
        enddo
    enddo     

    Rgsqr=Rgsqr/(2.0_dp*nmer*nmer)

end function radius_gyration_com


function radius_gyration(chain,nseg) result(Rgsqr)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nseg
    real(dp) :: Rgsqr
    integer :: i,j,k

    Rgsqr=0.0_dp
    
    do i=1,nseg
        do j=1,nseg
            do k=1,3
                Rgsqr=Rgsqr+(chain(k,i)-chain(k,j))**2
            enddo
        enddo
    enddo     

    Rgsqr=Rgsqr/(2.0_dp*nseg*nseg)

end function radius_gyration

function end_to_end_distance_com_rotation(rcom,nmer) result(Rendsqr)

    real(dp), intent(in) :: rcom(:,:)! dimension : (3,nnucl)
    integer, intent(in) :: nmer

    real(dp) :: Rendsqr
    integer :: k

    Rendsqr=0.0_dp
       
    do k=1,3
        Rendsqr=Rendsqr+(rcom(k,1)-rcom(k,nmer))**2
    enddo     

end function end_to_end_distance_com_rotation

function end_to_end_distance_com(chain,nmer,segcm) result(Rendsqr)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcm(:)  

    real(dp) :: Rendsqr
    integer :: isegcm, jsegcm, k

    Rendsqr=0.0_dp
       
    isegcm=segcm(1)
    jsegcm=segcm(nmer)

    do k=1,3
        Rendsqr=Rendsqr+(chain(k,isegcm)-chain(k,jsegcm))**2
    enddo     

end function end_to_end_distance_com

function nucleosomal_spacing_com_rotation(rcom,nmer) result(spacing)

    real(dp), intent(in) :: rcom(:,:) ! dimension : (3,nnucl)
    integer, intent(in) :: nmer

    real(dp) :: spacing(nmer-1)
    integer :: i,k
    real(dp) :: spacingsqr
        
    spacing=0.0_dp
    if(nmer>=2) then  ! need at least 2 unit/nucleosomes
        do i=1,nmer-1
            spacingsqr=0.0_dp
            do k=1,3
                spacingsqr=spacingsqr+(rcom(k,i+1)-rcom(k,i))**2
            enddo
            spacing(i)=sqrt(spacingsqr)
        enddo  
    endif           

end function nucleosomal_spacing_com_rotation

function nucleosomal_spacing_com(chain,nmer,segcom) result(spacing)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcom(:)  

    real(dp) :: spacing(nmer-1)
    integer :: i,k
    real(dp) :: spacingsqr
        
    spacing=0.0_dp
    if(nmer>=2) then  ! need at least 2 unit/nucleosomes
        do i=1,nmer-1
            spacingsqr=0.0_dp
            do k=1,3
                spacingsqr=spacingsqr+(chain(k,segcom(i+1))-chain(k,segcom(i)))**2
            enddo
            spacing(i)=sqrt(spacingsqr)
        enddo  
    endif           

end function nucleosomal_spacing_com


function gyr_tensor_com(chain,nmer,segcom) result(mat)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcom(:)  
    
    real(dp) :: mat(3,3)
    integer :: m, n, i, j, ii, jj
        
    mat=0.0_dp

    do m=1,3
        do n=m,3
            do i=1,nmer
                ii=segcom(i)
                do j=1,nmer
                    jj=segcom(j)
                    mat(m,n)=mat(m,n)+(chain(m,ii)-chain(m,jj))*(chain(n,ii)-chain(n,jj))
                enddo 
            enddo 
            mat(m,n)=mat(m,n)/(2.0_dp*nmer*nmer)
            mat(n,m)=mat(m,n) ! symmetric 
        enddo    
    enddo      

end function gyr_tensor_com

function gyr_tensor_com_rotation(rcom,nmer) result(mat)

    real(dp), intent(inout) :: rcom(:,:) ! dimension : (3,nnucl)
    integer, intent(in) :: nmer  
    
    real(dp) :: mat(3,3)
    integer :: m, n, i, j, ii, jj
        
    mat=0.0_dp

    do m=1,3
        do n=1,3
            mat(m,n)=0.0_dp
            do i=1,nmer
                do j=1,nmer
                    mat(m,n)=mat(m,n)+(rcom(m,i)-rcom(m,j))*(rcom(n,i)-rcom(n,j ))
                enddo 
            enddo 
            mat(m,n)=mat(m,n)/(2.0_dp*nmer*nmer)
            !mat(n,m)=mat(m,n) ! symmetric 
        enddo    
    enddo     
   

end function gyr_tensor_com_rotation

function bond_angles_com_rotation(rcom,nmer) result(bondangle)

    real(dp), intent(in) :: rcom(:,:) ! dimension : (3,nnucl)
    integer, intent(in) :: nmer

    real(dp) :: bondangle(nmer-2)

    ! .. local variables
    real(dp) :: u1(3), u2(3), absu1, absu2
    integer :: i, j, k ,isegcom, ipls1segcom


    bondangle=0.0_dp    

    if(nmer>=3) then  ! need at least 3 unit/nucleosomes
        isegcom=1
        ipls1segcom=2
        do k=1,3
            u1(k)=(rcom(k,ipls1segcom)-rcom(k,isegcom))
        enddo

        absu1=sqrt(dot_product(u1,u1))

        do i=3,nmer
            isegcom=i-1
            ipls1segcom=i
            do k=1,3
                u2(k)=(rcom(k,ipls1segcom)-rcom(k,isegcom))
            enddo
            absu2=sqrt(dot_product(u2,u2))
            bondangle(i-2)=acos(dot_product(u1,u2)/(absu2*absu1))
            u1=u2
            absu1=absu2
        enddo
    endif         

end function bond_angles_com_rotation


function bond_angles_com(chain,nmer,segcom) result(bondangle)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcom(:)
    real(dp) :: bondangle(nmer-2)

    ! .. local variables
    real(dp) :: u1(3), u2(3), absu1, absu2
    integer :: i,j,k ,isegcom, ipls1segcom


    bondangle=0.0_dp    

    if(nmer>=3) then  ! need at least 3 unit/nucleosomes
        isegcom=segcom(1)
        ipls1segcom=segcom(2)
        do k=1,3
            u1(k)=(chain(k,ipls1segcom)-chain(k,isegcom))
        enddo

        absu1=sqrt(dot_product(u1,u1))

        do i=3,nmer
            isegcom=segcom(i-1)
            ipls1segcom=segcom(i)
            do k=1,3
                u2(k)=(chain(k,ipls1segcom)-chain(k,isegcom))
            enddo
            absu2=sqrt(dot_product(u2,u2))
            bondangle(i-2)=acos(dot_product(u1,u2)/(absu2*absu1))
            u1=u2
            absu1=absu2
        enddo
    endif         

end function bond_angles_com


function dihedral_angles_com(chain,nmer,segcm) result(dihedral)

    real(dp), intent(in) :: chain(:,:)
    integer, intent(in) :: nmer
    integer, intent(in) :: segcm(:) 
    
    real(dp) :: dihedral(nmer-3)
    
    ! .. local variables
    integer  :: i,j,k 
    real(dp) :: u1(3),u2(3),u3(3),n123(3),n234(3)
    real(dp) :: absn123, absn234, absu2
    integer  :: isegcm,ipls1segcm,ipls2segcm,ipls3segcm
    real(dp) :: theta(nmer-3), costheta, x, y, sintheta, sintheta2 ,theta2,sintheta3, sintheta4
    


    if(nmer>=4) then  ! need at least 4 unit/nucleosomes

        isegcm    =segcm(1)
        ipls1segcm=segcm(2)
        ipls2segcm=segcm(3)
        ipls3segcm=segcm(4)

        do k=1,3
            u1(k)=(chain(k,ipls1segcm)-chain(k,isegcm    ))
            u2(k)=(chain(k,ipls2segcm)-chain(k,ipls1segcm))
            u3(k)=(chain(k,ipls3segcm)-chain(k,ipls2segcm))
        enddo

        n123=crossproduct(u1,u2)
        n234=crossproduct(u2,u3)

        absn123=sqrt(dot_product(n123,n123))  ! or sqrt(sum(n123*n123)) 
        absn234=sqrt(dot_product(n234,n234))
        absu2=sqrt(dot_product(u2,u2))
   
        y = dot_product(u2,crossproduct(n123,n234)) 
        x = absu2*dot_product(n123,n234)   
            
        theta(1)=atan2(y,x)  !formula See e.g. Bondel and Karplus j. comp. chem. 17 1132

        costheta= dot_product(n123,n234)/(absn123*absn234) ! cos first dihedral angle
        sintheta= dot_product(u2,crossproduct(n123,n234))/(absu2*absn123*absn234)

       ! following dihedral angles 

        do i=5,nmer 
        
            ! advance one unit
            n123=n234 
            absn123=absn234 
            u1=u2
            u2=u3 

            ipls2segcm=ipls3segcm
            ipls3segcm=segcm(i) 

            do k=1,3
                u3(k)=(chain(k,ipls3segcm)-chain(k,ipls2segcm))
            enddo

            n234=crossproduct(u2,u3)
            absn234=sqrt(dot_product(n234,n234))
            absu2=sqrt(dot_product(u2,u2))
           
            y = dot_product(u2,crossproduct(n123,n234)) 
            x = absu2*dot_product(n123,n234)
            
            theta(i-3)=atan2(y,x) ! i-3 dihedral angle

        enddo    
   
    endif   

    dihedral=theta    

end function dihedral_angles_com


function dihedral_angles_com_rotation(rcom,nmer) result(dihedral)

    real(dp), intent(in) :: rcom(:,:)
    integer, intent(in) :: nmer
    
    real(dp) :: dihedral(nmer-3)
    
    ! .. local variables
    integer  :: i,j,k 
    real(dp) :: u1(3),u2(3),u3(3),n123(3),n234(3)
    real(dp) :: absn123, absn234, absu2
    integer  :: isegcm,ipls1segcm,ipls2segcm,ipls3segcm
    real(dp) :: theta(nmer-3), costheta, x, y, sintheta, sintheta2 ,theta2,sintheta3, sintheta4
    


    if(nmer>=4) then  ! need at least 4 unit/nucleosomes

        isegcm    =1
        ipls1segcm=2
        ipls2segcm=3
        ipls3segcm=4

        do k=1,3
            u1(k)=(rcom(k,ipls1segcm)-rcom(k,isegcm    ))
            u2(k)=(rcom(k,ipls2segcm)-rcom(k,ipls1segcm))
            u3(k)=(rcom(k,ipls3segcm)-rcom(k,ipls2segcm))
        enddo

        n123=crossproduct(u1,u2)
        n234=crossproduct(u2,u3)

        absn123=sqrt(dot_product(n123,n123))  ! or sqrt(sum(n123*n123)) 
        absn234=sqrt(dot_product(n234,n234))
        absu2=sqrt(dot_product(u2,u2))
   
        y = dot_product(u2,crossproduct(n123,n234)) 
        x = absu2*dot_product(n123,n234)   
            
        theta(1)=atan2(y,x)  !formula See e.g. Bondel and Karplus j. comp. chem. 17 1132

        costheta= dot_product(n123,n234)/(absn123*absn234) ! cos first dihedral angle
        sintheta= dot_product(u2,crossproduct(n123,n234))/(absu2*absn123*absn234)

       ! following dihedral angles 

        do i=5,nmer 
        
            ! advance one unit
            n123=n234 
            absn123=absn234 
            u1=u2
            u2=u3 

            ipls2segcm=ipls3segcm
            ipls3segcm= i 

            do k=1,3
                u3(k)=(rcom(k,ipls3segcm)-rcom(k,ipls2segcm))
            enddo

            n234=crossproduct(u2,u3)
            absn234=sqrt(dot_product(n234,n234))
            absu2=sqrt(dot_product(u2,u2))
           
            y = dot_product(u2,crossproduct(n123,n234)) 
            x = absu2*dot_product(n123,n234)
            
            theta(i-3)=atan2(y,x) ! i-3 dihedral angle

        enddo    
   
    endif   

    dihedral=theta    

end function dihedral_angles_com_rotation



! inner product of vectors a and b with dimenstion 3
! similar to intrisic function dotproduct

function dotproduct(a,b) result(dotprod)
   
    real(dp), dimension(3) :: a, b 
    real(dp) :: dotprod
  
    dotprod = sum(a*b)

end function dotproduct

function crossproduct(a,b)result(crossprod)
  
  real(dp), intent(in), dimension(3) :: a, b
  real(dp), dimension(3) :: crossprod
  
  crossprod(1) = a(2)*b(3) - a(3)*b(2)
  crossprod(2) =-a(1)*b(3) + a(3)*b(1)
  crossprod(3) = a(1)*b(2) - a(2)*b(1)

end function crossproduct



!  Read segment id number to which we assign center of mass of one histone from file
!  segment id numbers stored in array segcm. nucleosome has nnucl histones 
!  input integer nnucl : number of nucleosome 
!        character(len=lenfname) filename : file that contains segcom 
!  input/output integer  segcom(n=nnucl) 


!subroutine read_segcom_id(segcom,nnucl,filename)

subroutine make_segcom(segcom,nnucl,filename)

    !use mpivars
    use myutils
    use parameters, only : lenfname
    
    !     .. arguments 
    integer, intent(inout) :: segcom(:)
    character(lenfname), intent(in) :: filename
    integer,  intent(in) :: nnucl

    !      .. local variables
    integer :: ios, un  ! un = unit number
    integer :: s
    character(80) :: istr,str,letter

    !     .. reading in of variables from file
    open(unit=newunit(un),file=filename,iostat=ios,status='old')
    if(ios/=0 ) then
        write(istr,'(I2)')ios
        str='Error opening file '//trim(adjustl(filename))//' : iostat = '//istr
        print*,str
        stop
    endif
    
    s=0
    ios=0

    do while (s<nnucl.and.ios==0)
        s=s+1
        read(un,*,iostat=ios)segcom(s)
    enddo

    if(s/=nnucl.or.ios/=0) then 
        str="reached end of file before all elements read or ios error"
        print*,str
        str="read file "//trim(adjustl(filename))//" failed"
        print*,str
!        call MPI_FINALIZE(ierr)
        stop
    endif

    close(un)

end subroutine make_segcom

subroutine write_chain_struct(write_struct,info)

    use globals, only : cuantas,nnucl
    use myutils, only : lenText
    use chains, only : Rgsqr,Rendsqr,bond_angle,dihedral_angle,nucl_spacing,energychainLJ
    use chains, only : no_overlapchain, Asphparam
    implicit none 

    logical, intent(in) :: write_struct
    integer, intent(inout) :: info
 
    ! .. local
    character(len=lenText) :: filename
    integer :: un_dihedral,un_bond,un_Rg,un_Rend, un_dist, un_VdW, un_Asph
    integer :: c,s,nbonds,nangles,ndihedrals

    info=0

    if(write_struct) then
    
        filename="dihedral."
        un_dihedral=open_chain_struct_file(filename,info)
        filename="bond."
        un_bond=open_chain_struct_file(filename,info)
        filename="Rgalpha."
        un_Rg=open_chain_struct_file(filename,info)
        filename="Rend."
        un_Rend=open_chain_struct_file(filename,info)
        filename="spacing."
        un_dist=open_chain_struct_file(filename,info)
        filename="VdWenergy."
        un_VdW=open_chain_struct_file(filename,info)
        filename="Aspherparam."
        un_Asph=open_chain_struct_file(filename,info)
        
        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3
    
        do c=1,cuantas
            write(un_Rg,*)Rgsqr(c)
            write(un_Rend,*)Rendsqr(c)
            write(un_bond,*)(bond_angle(s,c),s=1,nangles)
            write(un_dihedral,*)(dihedral_angle(s,c),s=1,ndihedrals)
            write(un_dist,*)(nucl_spacing(s,c),s=1,nbonds)
            write(un_VdW,*)energychainLJ(c),no_overlapchain(c)
            write(un_Asph,*)Asphparam(c)
        enddo 

        close(un_dihedral)
        close(un_bond)
        close(un_Rg)
        close(un_Rend)  
        close(un_dist)
        close(un_VdW)
        close(un_Asph)
      
    endif

        
end subroutine write_chain_struct


function open_chain_struct_file(filename,info)result(un)

    !use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    character(len=lenText) :: filename
    integer, optional, intent(inout) :: info

    integer :: un

    ! local
    character(len=lenText) :: istr
    character(len=25) :: fname
    integer :: ios
    logical :: exist
    
    if (present(info))  info=0 ! init
    
    write(istr,'(I4)')0
    fname=trim(adjustl(filename))//trim(adjustl(istr))//'.dat'

    inquire(file=fname,exist=exist)
    if(.not.exist) then
        open(unit=newunit(un),file=fname,status='new',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    else
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    endif    
end function open_chain_struct_file

function open_output_file(filename,info)result(un)

    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    character(len=lenText) :: filename
    integer, optional, intent(inout) :: info

    integer :: un

    ! local
    character(len=lenText) :: istr
    character(len=25) :: fname
    integer :: ios
    logical :: exist
    
    if (present(info))  info=0 ! init
    
    fname=trim(adjustl(filename))

    inquire(file=fname,exist=exist)
    if(.not.exist) then
        open(unit=newunit(un),file=fname,status='new',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    else
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    endif    
end function open_output_file


subroutine write_phosphate_pairs(write_pairs,conf,info)

    use globals, only : cuantas,nseg
    use myutils, only : lenText, newunit
    use chains, only : type_of_monomer, nneigh, indexconfpair, distphoscutoff
    use parameters, only : ta


    implicit none 

    logical, intent(in) :: write_pairs
    integer, intent(in) :: conf
    integer, intent(inout) :: info

 
    ! .. local
    character(len=lenText) :: filename
    integer :: un_pp
    integer :: s, j, tPhos
    character(len=10) ::istr

    info=0

    if(write_pairs) then

        tPhos = find_type_phosphate()
        
        write(istr,'(I4)')conf
        filename='phosphate_pairs.'//trim(adjustl(istr))//'.log'
        
        !     .. opening file
        open(unit=newunit(un_pp),file=filename)

        !un_pp=open_chain_struct_file(filename,info)

        write(un_pp,*)"phosphate pairs"
        write(un_pp,*)"value ta=",ta, "value tPhos=",tPhos
        write(un_pp,*)"distphoscutoff=",distphoscutoff
        
        write(un_pp,*)"conf=",conf
        do s=1,nseg
            write(un_pp,*)s,type_of_monomer(s),nneigh(s,conf),&
                (indexconfpair(s,conf)%elem(j),j=1,nneigh(s,conf))
        enddo
            
        close(un_pp)
          
    endif
        
end subroutine write_phosphate_pairs


subroutine write_phosphate_triplets(write_triplets,conf,info)

    use globals, only : cuantas,nseg
    use myutils, only : lenText, newunit
    use chains, only : type_of_monomer, ntriplet, indexconftriplet, distphoscutoff
    use parameters, only : ta

    implicit none 

    logical, intent(in) :: write_triplets
    integer, intent(in) :: conf
    integer, intent(inout) :: info
 
    ! .. local
    character(len=lenText) :: filename
    integer :: un_pp, num_triplets
    integer :: s, j, tPhos
    character(len=10) ::istr

    info=0

    num_triplets=0
    do s=1,nseg
        if(ntriplet(s,conf)>0) num_triplets=num_triplets+ntriplet(s,conf)
    enddo
    num_triplets = num_triplets/3


    print*,"write_triplets=",write_triplets
    
    if(write_triplets) then

        tPhos = find_type_phosphate()
        
        write(istr,'(I4)')conf
        filename='phosphate_triplets.'//trim(adjustl(istr))//'.log'
        
        !     .. opening file
        open(unit=newunit(un_pp),file=filename)

        write(un_pp,*)"phosphate triplets"
        write(un_pp,*)"value ta=",ta, "value tPhos=",tPhos
        write(un_pp,*)"distphoscutoff=",distphoscutoff
        write(un_pp,*)"conf=",conf
        write(Un_pp,*)"number of triplets = ", num_triplets

        do s=1,nseg
            if(ntriplet(s,conf)>0) then 
                write(un_pp,*)s,type_of_monomer(s),ntriplet(s,conf),&
                    (indexconftriplet(1,s,conf)%elem(j),indexconftriplet(2,s,conf)%elem(j),j=1,ntriplet(s,conf))
            endif
        enddo
            
        close(un_pp)
          
    endif
        
end subroutine write_phosphate_triplets


subroutine set_mapping_num_to_char(type_of_monomer_num_to_char)

    use globals, only : nsegtypes,nseg
    use chains, only : type_of_monomer, type_of_monomer_char

    character(len=3), intent(inout) :: type_of_monomer_num_to_char(:)  
        
    integer :: t, k
    logical :: notfound

    do t=1,nsegtypes
        k=1
        notfound=.true.
        do while (k<=nseg .and. notfound)
            if(type_of_monomer(k)==t) then
                type_of_monomer_num_to_char(t)=type_of_monomer_char(k)
                notfound=.false. ! stops
                k=1 ! reset
            endif 
            k=k+1
        enddo        
    enddo

end subroutine set_mapping_num_to_char

! reads in MT pdb file to assign AA chain elements. 
! input  : fname   = char(*)  filename of MTpdb file ) 
!        : nsegAA  = integer = number of AA residues
!        : typesAA = array of integer:  the s element gives number of AA type
! output : nelemAA = array of integer:  the s element gives number of element the sth AA consits of
!        : chain_elem(3,s)%elem(j) = position of jthe element of AA number s. Enumarate the AA in order
!        : vnucl(j,t) = real(dp) volume of DNA or AA  element j of DNA or DNA type number t
!        : nucl_elem_type(j,t) = character(len=3) type of DNA or AA element j of DNA or AA type number t
!        : elem_charge(s) = array of integers: the element of AA element of AA s that is chargeable
!        : info     = integer val=0 assignment  succesfull val/=0 error              

subroutine read_nucl_elements(fname,nsegAA,nelemAA,chain_elem,typeAA,vnucl,nucl_elem_type,elem_charge,info)

    use globals, only    : DEBUG, nsegtypes,nsegtypesAA
    use myutils, only    : newunit, lenText, error_handler
    use myio, only       : myio_err_chainsfile
    use parameters, only : lenfname, vpol, vsol
    use chains, only     : var_darray, ismonomer_chargeable, mapping_num_to_char

    character(lenfname),   intent(in)                            :: fname
    integer,               intent(in)                            :: nsegAA
    integer, dimension(:), intent(in)                            :: typeAA
    integer, dimension(:), intent(inout)                         :: nelemAA
    type(var_darray), dimension(:,:), allocatable, intent(inout) :: chain_elem
    real(dp), dimension(:,:) , intent(inout)                     :: vnucl
    character(len=3), dimension(:,:), allocatable  , intent(inout) :: nucl_elem_type
    integer, dimension(:), allocatable, intent(inout)            :: elem_charge
    integer, intent(out),optional                                :: info 

    ! local variables

    integer          :: ios, un, s, sAA, j, k, ttype ,t
    character(len=3) :: elem_type(nsegtypesAA), voltype
    real(dp)         :: x(3)
    logical          :: isReadGood
    integer          :: num_elem_CA
    character(len=3) :: list_type_char_DNA(6)  ! list of DNA elements in char
    integer          :: list_type_int_DNA(6)   ! list of DNA elements in integer 
    integer          :: AAid, DNAid, nelem_types
    character(len=3) :: type_char  

    if (present(info)) info = 0


    open(unit=newunit(un),file=fname,status='old',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

    ! allocate objects chain_elem and nucl_elem_type

    allocate(chain_elem(3,nsegAA))
    nelem_types=size(vnucl,1)
    allocate(nucl_elem_type(nelem_types,nsegtypes))
    nucl_elem_type=" "

    do sAA=1,nsegAA

        read(un,*,iostat=ios) ! comment
        if(ios/=0) isReadGood=.false.
        read(un,*,iostat=ios) AAid ! type number
        if(ios/=0) isReadGood=.false.
        read(un,*,iostat=ios) nelemAA(sAA)
        if(ios/=0) isReadGood=.false.
        
        do k=1,3
            allocate(chain_elem(k,sAA)%elem(nelemAA(sAA)))
        enddo     
        
        elem_type=""
        do j=1,nelemAA(sAA)    
            
            read(un,*,iostat=ios)elem_type(j),x(1),x(2),x(3) 
            !write(100,*)elem_type(j),x(1),x(2),x(3)            
            do k=1,3
                chain_elem(k,sAA)%elem(j)=x(k)
            enddo  
            
            vnucl(j,AAid)=find_vol_elem(elem_type(j))  
            nucl_elem_type(j,AAid)=elem_type(j) 

            ! assignment doubles since multiple sAA  corresponds to same type 

            if(ios/=0) isReadGood=.false.

        enddo 

        num_elem_CA=find_CA_elem(elem_type,nelemAA(sAA))

        ! assignment doubles since multiple sAA  corresponds to same type 

        if(ismonomer_chargeable(AAid)) then ! ismonomer_chargeable has type number as input !! 
            elem_charge(AAid)=find_chargeable_elem(elem_type,nelemAA(sAA))
            if(elem_charge(AAid)==-1) then 
                print*,"=> sAA=",sAA, "AAid=",AAid,"char=",mapping_num_to_char(AAid)
                print*,"=> elem_type = ", elem_type
                call error_handler(1,"read_nucl_elements")
            endif    
        else 
            elem_charge(AAid)=-1
        endif      
       
        
        if(num_elem_CA/=1) then 
            ! swapping
            call swap_char3(elem_type(num_elem_CA),elem_type(1))
            call swap_real(vnucl(num_elem_CA,AAid), vnucl(1,AAid))

            do k=1,3
                call swap_real(chain_elem(k,sAA)%elem(num_elem_CA),chain_elem(k,sAA)%elem(1))
            enddo
            call swap_char3(nucl_elem_type(num_elem_CA,AAid),nucl_elem_type(1,AAid)) 
        endif  

        num_elem_CA=find_CA_elem(elem_type,nelemAA(sAA))

        if(DEBUG) then 
            print*,"check swapping "
            print*,"elem_type   = ",elem_type
            print*,"num_elem_CA = ",num_elem_CA
        endif

    enddo
    

    ! .. assign volume to vnucl that is from  DNA elements 
    ! .. above loops assign vnucl only for the AAs 

    list_type_char_DNA = (/"P  ","S  ","A  ","C  ","G  ","T  "/)
    list_type_int_DNA  = 0

    ! .. make list_type_int_DNA"

    do k=1,6                                ! 6 = number DNA elements 
        type_char= list_type_char_DNA(k)
        do t=1,nsegtypes
            if( mapping_num_to_char(t) == type_char)  then 
                list_type_int_DNA(k) = t
            endif
        enddo 
    enddo

    do k=1,6
        DNAid=list_type_int_DNA(k)
        if(DNAid/=0) then  ! in case the conformation is incomplete and does not 
            vnucl(1,DNAid) = vpol(DNAid)*vsol
            nucl_elem_type(1,DNAid)=list_type_char_DNA(k)  
        endif
    enddo 

    ! .. end assignment DNA volume 

    ! .. assign index to elem_charge  that is from  DNA elements 
    ! .. above loops assign elem_charge only for the AAs

    do k=1,6
        DNAid=list_type_int_DNA(k)
        if(DNAid/=0) then  ! in case the conformation is incomplete and does not 
            if(list_type_char_DNA(k)=="P") then 
                elem_charge(DNAid) = 1
            else 
                elem_charge(DNAid) = -1
            endif
        endif
    enddo 



    if(DEBUG) then 
        print*,"DEBUG: output of read_nucl_elements : chain_elem"
        print*," list DNAid"
        do k=1,6
            print*,"DNAid=",list_type_int_DNA(k)
        enddo        
        print*,"DEBUG: output of read_nucl_elements "
        print*,"s              k                  j        chain_elem(k,s)%elem(j)" 
        do s=1,nsegAA
            do k=1,3
                do j=1,nelemAA(s)
                    print*,"s=",s,"k=",k,"j=",j," ",chain_elem(k,s)%elem(j)  
                enddo 
            enddo 
            print*," "
        enddo
        print*," "
        print*,"DEBUG: output of read_nucl_elements : vnucl  "
        print*,"j               t               vnucl(j,t)"
        do t=1,nsegtypes
            do j=1,nsegtypesAA
                print*,"j=",j," t=",t," ",vnucl(j,t)
            enddo 
        enddo

        print*," "

        print*,"DEBUG: output of read_nucl_elements : nucl_elem_type  "
        print*,"j               t              nucl_elem_type  (j,t)"
        do t=1,nsegtypes
            do j=1,nsegtypesAA
                print*,"j=",j," t=",t," ",nucl_elem_type(j,t)
            enddo 
        enddo
    endif       

end subroutine read_nucl_elements


! Shift chain_elem(k,s)%elem(j) by position of the CA of AA that is the origin of the "CM" rotation triplet
! translation vector is chain_CA(k)=chain_elem(k,segnumcm)%elem(1)
! input  : integer nsegAA : number of AA 
!          integer nelemAA : array hold
!          integer segnumcm  : segment number of cm or ration origin           
! input/output : type(var_darray) : chain_elem

subroutine shift_nucl_elements(nsegAA,nelemAA,segnumcm,chain_elem)
    
    use chains, only : var_darray ! type def  

    integer, intent(in)                  :: nsegAA
    integer, dimension(:), intent(in)    :: nelemAA
    integer, intent(in)                  :: segnumcm
    type(var_darray), dimension(:,:), allocatable, intent(inout) :: chain_elem

    integer :: s,k,j,i
    real(dp) :: chain_CA(3)
 
    ! translation vector
    do i=1,3 
        chain_CA(i)=chain_elem(i,segnumcm)%elem(1) ! 1= CA
    enddo     
    ! substract 
    do s=1,nsegAA
        do j=1,nelemAA(s)
            do k=1,3
                chain_elem(k,s)%elem(j)=chain_elem(k,s)%elem(j) - chain_CA(k)  
            enddo
        enddo        
   enddo     

end subroutine shift_nucl_elements

! Function find array index of  elem_type that matches CA
! input  : character(3) elem_type : array of characters denoting type
! return : integer elem : index 

function find_CA_elem(elem_type,nelem)result(elem)

    character(len=3), intent(in) :: elem_type(:)
    integer, intent(in) :: nelem

    integer :: elem, j 
    logical :: IsNotFound

    IsNotFound=.true.
    j=0
    do while (IsNotFound.and.j<=nelem)
        j=j+1
        if(elem_type(j)=="CA") IsNotFound=.false.
    enddo 

    elem=j 
    if(j==0) print*,"error No CA found "

end function find_CA_elem



! Function assign/finds volume that matches with elem_type
! input  : character(3) elem_type : array of characters denoting type
! return : real(dp) vol_elem      : if no match found vol_elem is set to 0   

function find_vol_elem(elem_type)result(vol_elem)

    use parameters, only : vnucl_type, vnucl_type_char

    character(len=3), intent(in) :: elem_type

    real(dp):: vol_elem

    integer :: elem, j , nelem_types
    logical :: IsNotFound

    nelem_types=size(vnucl_type)
    IsNotFound=.true.
    j=0
    do while (IsNotFound.and.j<=nelem_types)
        j=j+1
        if(vnucl_type_char(j)==elem_type) IsNotFound=.false.
    enddo 

    if(j/=0) then
        vol_elem=vnucl_type(j)
    else 
        vol_elem=0.0_dp
        print*,"error vol nucl not found "
    endif    

end function find_vol_elem

! Function locates the array index of the elem_type element of AA that is chargeable 
! input  : character(3) elem_type(*) : sequence or array of characters denoting type 
!        : integer nelem : number of elements
! return : integer elem : array index          

function find_chargeable_elem(elem_type,nelem)result(elem)

    use parameters, only : vnucl_type, vnucl_type_char,vnucl_type_isChargeable

    character(len=3), intent(in) :: elem_type(:)
    integer, intent(in) :: nelem

    integer :: elem, j, k, ncharge, num_elem_charge, nelem_types 
    logical :: IsNotFound

    nelem_types = size(vnucl_type)
    ncharge = 0
    do k=1,nelem
        IsNotFound=.true.
        j=0
        do while (IsNotFound.and.j<=nelem_types)
            j=j+1
            if(elem_type(k)==vnucl_type_char(j)) IsNotFound=.false.
        enddo 
        if(vnucl_type_isChargeable(j)) then  ! kth element of elem_type is charged 
            ncharge = ncharge+1  ! raise
            num_elem_charge = k  ! this the index of elem_type
       endif
    enddo

    if(ncharge==1) then !  Should be only one chargeable element  in sequence
        elem = num_elem_charge
    else if (ncharge >1 ) then 
        elem =  0
        print*,"Error found more then one element AA charged"
    else 
        elem = -1
        print*,"Error found no element AA charged"
    endif    

end function find_chargeable_elem

function find_type_phosphate()result(tPhos)

    use myutils, only : lenText, error_handler
    use globals, only : nseg
    use chains, only : type_of_monomer,type_of_monomer_char
    integer :: tPhos
    integer :: s, ios
    character(len=lenText) :: text

    ios=0
    tPhos=0
    do s=1, nseg
        if(type_of_monomer_char(s)=="P")  tPhos =type_of_monomer(s)
    enddo   
    
    if(tPhos==0) ios=1  ! phosphate monomer "P" not found in type_of_monomer_char 

    if(ios>0) then
        text='find_type_phosphate : phosphate monomer is not defined'
        call error_handler(ios,text)
    endif

end function find_type_phosphate

subroutine swap_int(a,b)
    
    integer, intent(inout) :: a, b
    integer :: tmp 

    tmp = a
    a = b 
    b = tmp 

end subroutine swap_int   


subroutine swap_real(a,b)
    
    real(dp), intent(inout) :: a, b
    real(dp) :: tmp 

    tmp = a
    a = b 
    b = tmp 
    
end subroutine 

subroutine swap_char3(a,b)
    
    character(len=3), intent(inout) :: a, b
    character(len=3) :: tmp 

    tmp = a
    a = b 
    b = tmp 

end subroutine swap_char3



subroutine read_nucl_orient_triplets(fname,nnucl,nuc_orient_triplet,info)

    use myutils, only : newunit, lenText
    use parameters, only : lenfname
    use myio, only : myio_err_chainsfile, myio_err_readfile 

    character(lenfname),    intent(in)      :: fname
    integer, intent(in)                    :: nnucl
    integer, dimension(:,:), intent(inout) :: nuc_orient_triplet
    integer, intent(out), optional         :: info

    ! local arguments

    integer  :: ios, un, s,  i
    integer :: ix(3)
    logical :: isReadGood

    if (present(info)) info = 0
    isReadGood=.true.

    open(unit=newunit(un),file=fname,status='old',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_readfile
        return
    endif

    do i=1,nnucl
        read(un,*,iostat=ios)ix(1),ix(2),ix(3) 
        if(ios/=0) isReadGood=.false.
        nuc_orient_triplet(i,:)=ix     
    enddo
 
    if(isReadGood.eqv..False. ) then
        print*, 'Error reading : ',fname,' file : iostat =', ios
        if (present(info)) info= myio_err_readfile 
    endif    

    close(un)
    
end subroutine read_nucl_orient_triplets

! Allocates object chain_elem
! input : integer : nsegAA :  number of AAs, has to been assigned
!         integer : nelem(s) : number of elems for AA number s, has been assigned 
! output: type(var_darray) chain_elem(k,s)%elem(s) allocated

subroutine allocate_chain_elements(nsegAA,nelemAA,chain_elem)

    use chains, only : var_darray

    integer, intent(in)               :: nsegAA
    integer, dimension(:), intent(in) :: nelemAA
    type(var_darray), dimension(:,:), allocatable, intent(inout) :: chain_elem

    ! local arguments

    integer   :: s, k, j
    
    allocate(chain_elem(3,nsegAA))

    do s=1,nsegAA
        do k=1,3
            allocate(chain_elem(k,s)%elem(nelemAA(s)))
        enddo     
    enddo 

    !call print_nucl_elements(nsegAA,nelemAA,chain_elem)

end subroutine allocate_chain_elements


subroutine allocate_nucl_chain_elements(nnucl,nsegAA,nelemAA,chain_elem_rot)

    use chains, only : var_darray

    integer, intent(in)               :: nnucl
    integer, intent(in)               :: nsegAA
    integer, dimension(:), intent(in) :: nelemAA
    type(var_darray), dimension(:,:,:), allocatable, intent(inout) :: chain_elem_rot

    ! local arguments

    integer   :: s, k, n
    
    allocate(chain_elem_rot(3,nsegAA,nnucl))

    do n=1,nnucl
        do s=1,nsegAA
            do k=1,3
                allocate(chain_elem_rot(k,s,n)%elem(nelemAA(s)))
            enddo     
        enddo
    enddo     

end subroutine allocate_nucl_chain_elements

subroutine print_nucl_elements(nsegAA,nelemAA,chain_elem)

    use chains, only : var_darray

    integer, intent(in)               :: nsegAA
    integer, dimension(:), intent(in) :: nelemAA
    type(var_darray), dimension(:,:), allocatable, intent(inout) :: chain_elem

    integer   :: s, k, j

    do s=1,nsegAA
        do k=1,3
            do j=1,nelemAA(s)
                print*,"s=",s,"k=",k,"j=",j," ",chain_elem(k,s)%elem(j)
            enddo     
        enddo
    enddo  

end subroutine print_nucl_elements

! Computes segment numbers that correspond to first AA of nth nucleosome.
! inputs : nseg      = number of segments
!          nsegtypes = number of segment types
!          nnucl     = number of nuclesomes
! output : segAAstart = array of integer = nth val= segment number that is first AA of nth nucleosome

subroutine compute_segnumAAstart(nseg,nsegtypes,nnucl,segnumAAstart)

    use chains, only : type_of_monomer, type_of_monomer_char, mapping_num_to_char  

    integer, intent(in) :: nseg 
    integer, intent(in) :: nsegtypes
    integer, intent(in) :: nnucl 
    integer, dimension(:), intent(inout) :: segnumAAstart

    ! .. local arguments

    character(len=3) :: list_type_char_DNA(6)  ! list of DNA elements in char
    integer          :: list_type_int_DNA(6)   ! list of DNA elements in integer 
    integer          :: i, s, t, k, nnucl_counter, type_int
    character(len=3) :: type_char  
    logical          :: isMonomerDNA, isPreviousMonomerDNA

    ! .. init list DNA segments  
    
    list_type_char_DNA = (/"P  ","S  ","A  ","C  ","G  ","T  "/)
    list_type_int_DNA  = 0

    ! .. make list_type_init_DNA"

    do k=1,6                                ! 6 = number DNA elements 
        type_char= list_type_char_DNA(k)
        do t=1,nsegtypes
            if( mapping_num_to_char(t) == type_char)  then 
                list_type_int_DNA(k)=t
            endif
        enddo 
    enddo

    ! locate segment numbers that is correspond to first AA of nth nucl 

    nnucl_counter=1
    isPreviousMonomerDNA=.true.
    
    do s=1,nseg 
        type_int=type_of_monomer(s)
        isMonomerDNA=.false.
        do k=1,6
            if(type_int == list_type_int_DNA(k)) isMonomerDNA=.true.        
        enddo     
        if(isPreviousMonomerDNA .and. .not. isMonomerDNA) then 
            segnumAAstart(nnucl_counter)=s 
            isPreviousMonomerDNA =.false.
            nnucl_counter=nnucl_counter+1
        endif     
        if(isMonomerDNA) isPreviousMonomerDNA =.true.
    enddo 
                
end subroutine compute_segnumAAstart


! Computes segment numbers that correspond to last AA of nth nucleosome.
! inputs : nnucl     = number of nuclesomes
!          nsegAA    = number of AA segments per nucleosome
!          segnumAAstart = segment numbers that correspond to first AA of nth nucleosome.
! output : segnumAAsend = array of integer = nth val= segment number that is last AA of nth nucleosome

subroutine compute_segnumAAend(nnucl,nsegAA,segnumAAstart,segnumAAend)

    integer, intent(in) :: nnucl
    integer, intent(in) :: nsegAA
    integer, dimension(:), intent(in) :: segnumAAstart
    integer, dimension(:), intent(inout) :: segnumAAend 

    integer :: n 

    do n=1,nnucl 
        segnumAAend(n)=segnumAAstart(n)+nsegAA-1
    enddo 

end subroutine compute_segnumAAend 

subroutine make_nelem(nseg,nsegAA,nnucl,segnumAAstart,nelemAA,nelem)

    integer,               intent(in) :: nseg
    integer,               intent(in) :: nsegAA
    integer,               intent(in) :: nnucl
    integer, dimension(:), intent(in) :: segnumAAstart
    integer, dimension(:), intent(in) :: nelemAA
    integer, dimension(:), intent(inout) :: nelem

    integer :: n , sn, s

    nelem=1

    do n=1,nnucl 
        do sn=1,nsegAA
            s=sn+segnumAAstart(n)-1
            nelem(s)=nelemAA(sn)
        enddo 
    enddo 
    
end subroutine make_nelem


! Computes chain_index(k,s)%elem(j) coordinate (component k) of segment s and element j
! by adding chain_elem_rot to chain_rot(k,s)
! Note chain_elem_rot distance relative to cm such be relative to location CA of every AA  
! inputs : integer nseg      = number of segments
!          integer nsegAA    = number of AA segments per nucleosome
!          integer nnucl     = number of nucleosomes
!          integer segnumAAstart = segment numbers that corresponds to first AA of nth nucleosome
!          integer segnumAAsend  = segment numbers that corresponds to last AA of nth nucleosome
!          integer nelemAA    = array list number elements for each segments
!          type(var_darray) chain_rot = array of coordinates of chain for AA only the CA elements
!          type(var_darray) chain_elem_rot => chain_elem_rot(k,sn,n)%elem(j) oriented AA elements for nth nucleosome 
! output : type(var_darray) chain_elem_index => chain_index(k,s)%elem(j) = coordinate (component k) of segment s and element j


subroutine add_chain_rot_and_chain_elem_rot(nseg,nsegAA,nnucl,segnumAAstart,segnumAAend,nelemAA,&
    chain_rot,chain_elem_rot,chain_elem_index) 

    use chains, only : var_darray

    integer, intent(in) :: nseg
    integer, intent(in) :: nsegAA
    integer, intent(in) :: nnucl
    integer, dimension(:), intent(in) :: segnumAAstart
    integer, dimension(:), intent(in) :: segnumAAend 
    integer, dimension(:), intent(in) :: nelemAA
    real(dp), dimension(:,:), intent(in) :: chain_rot 
    type(var_darray), dimension(:,:,:), allocatable, intent(in) :: chain_elem_rot
    type(var_darray), dimension(:,:), allocatable, intent(inout) ::chain_elem_index

    integer :: n, s, j, k, sn

    ! chain_elem_rot distance relative to cm such be relative to location CA of every AA  

    ! AA coordinates
    do n=1,nnucl                    ! loop over number of nucleosomes
        do sn=1,nsegAA              ! loop over number of AA 
            s=sn+segnumAAstart(n)-1 ! segment number 
            do j=1,nelemAA(sn)      ! loop over number of chain elements for segment s 
                do k=1,3
                    chain_elem_index(k,s)%elem(j)= chain_rot(k,s)+(chain_elem_rot(k,sn,n)%elem(j)-chain_elem_rot(k,sn,n)%elem(1))
                enddo
            enddo 
        enddo 
    enddo  

    ! DNA coordinates

    do s=1,segnumAAstart(1)-1 
        do k=1,3
            chain_elem_index(k,s)%elem(1)= chain_rot(k,s)
        enddo  
    enddo        
    do n=1,nnucl-1 
        do s=segnumAAend(n)+1,segnumAAstart(n+1)-1 
            do k=1,3
                chain_elem_index(k,s)%elem(1)= chain_rot(k,s)
            enddo 
        enddo 
    enddo 
    do s=segnumAAend(nnucl)+1,nseg 
        do k=1,3
            chain_elem_index(k,s)%elem(1)= chain_rot(k,s)
        enddo  
    enddo           

end subroutine add_chain_rot_and_chain_elem_rot

! Translate chain(s,k)  coordinate (component k) of segment s  
! translation_vector(k)

subroutine translate_chain(nseg,chain,translation_vector) 

    integer, intent(in) :: nseg
    real(dp), intent(inout) :: chain(:,:)
    real(dp), dimension(3), intent(in) :: translation_vector

    integer :: s, k

    do s=1,nseg              ! loop over segment 
        do k=1,3
            chain(k,s) = chain(k,s)+translation_vector(k)
        enddo     
    enddo     

end subroutine translate_chain


! Translate chain_index(k,s)%elem(j) coordinate (component k) of segment s and element j by 
! translation_vector(k)

subroutine translate_chain_elem_index(nseg,nelem,chain_elem_index,translation_vector) 

    use chains, only : var_darray

    integer, intent(in) :: nseg
    integer, dimension(:), intent(in) :: nelem
    type(var_darray), dimension(:,:), allocatable, intent(inout) ::chain_elem_index
    real(dp), dimension(3), intent(in) :: translation_vector

    integer :: s, j, k

    do s=1,nseg              ! loop over segment
        do j=1,nelem(s)      ! loop over number of chain elements for segment s 
            do k=1,3
                chain_elem_index(k,s)%elem(j)= chain_elem_index(k,s)%elem(j)+translation_vector(k)
            enddo
        enddo 
    enddo     

end subroutine translate_chain_elem_index


! Finds phosphate pairs for given conformation number conf 
! Conformation is stored in chain
! Assigns  nneigh(s,conf) and indexconfpair(s,conf)%elem(j) with 0<=j<=neigh(s,conf)
! input integer :: nseg : number atoms/segment
!       integer :: conf : conformation number
!       integer :: tPphos : number associated with type of phosphates
!       real(dp) :: sqrDphoscutoff : squared distance of cutoffdistance citeria for pair
!       real(dp) :: chain(3,nseg) : hold coordiante of backbone confomation for all atom/segment
!       real(dp) :: LX,Ly,Lz : dimension lattice/box in nm 
! output assigment  
!       integer :: nneight(nseg,conf) : number of 'neighbors' that a phosphate s has for conformation conf
!                  if type of s in not a phospate then value zero
!       type(var_iarray) ::  indexconfpair(s,conf)%elem(j) : layer number of neigbor j of conf alpha and segment number s

subroutine find_phosphate_pairs(nseg,conf,tPhos,sqrDphoscutoff,chain,Lx,Ly,Lz)

    use chains, only :  type_of_monomer,indexconfpair
    use chains, only : nneigh, indexconfpair, distphoscutoff
    use parameters, only : tA 
    use parameters, only : pbc_chains
    use volume, only : delta, linearIndexFromCoordinate
    use myutils, only : newunit
    
    integer, intent(in) :: nseg
    integer, intent(in) :: conf
    integer, intent(in) :: tPhos
    real(dp), intent(in) :: sqrDphoscutoff
    real(dp), intent(in) :: chain(3,nseg)
    real(dp), intent(in) :: Lx,Ly,Lz

    integer, parameter :: maxnneigh = 10

    integer :: s, sprime, i, j
    integer :: xi, yi, zi , idx
    integer, dimension(:,:), allocatable :: list_of_pairs, index_of_pairs
    real(dp) :: sqrdist
    
    character(len=100) :: fname
    integer :: un_pp
    character(len=10) ::istr
    
    allocate(list_of_pairs(nseg,maxnneigh))
    allocate(index_of_pairs(nseg,maxnneigh))

    do s=1,nseg 
        nneigh(s,conf)=0
        if(type_of_monomer(s)==tPhos) then ! tPhos equiv to ta which is not set yet 
            do sprime=1,nseg
                if(type_of_monomer(sprime)==tPhos ) then
                    if(s/=sprime) then ! prevent s=sprime being counted as a pair
                       
                        sqrdist=0.0_dp
                        do i=1,3
                            sqrdist=sqrdist+(chain(i,s)-chain(i,sprime))**2
                        enddo
    
                        if(sqrdist<=sqrDphoscutoff) then ! comparing square of distance to square of cutoff  
                            ! accept s and sprime are a pair
                            nneigh(s,conf)=nneigh(s,conf)+1
                            list_of_pairs(s,nneigh(s,conf))=sprime ! temporarily storage of  segment number of neighbor to (s,conf)

                            ! transforming form real- to lattice coordinates                 
                            
                            if(pbc_chains) then  
                                xi = int(pbc(chain(1,sprime),Lx)/delta)+1
                                yi = int(pbc(chain(2,sprime),Ly)/delta)+1
                                zi = int(pbc(chain(3,sprime),Lz)/delta)+1
                            else 
                                xi = int(chain(1,sprime)/delta)+1
                                yi = int(chain(2,sprime)/delta)+1
                                zi = int(chain(3,sprime)/delta)+1
                            endif    

                            call linearIndexFromCoordinate(xi,yi,zi,idx)
                            !idx = coordtoindex(xi,yi,zi) ! hash-table look up
                            
                            index_of_pairs(s,nneigh(s,conf))=idx  ! temporarily storage of index of neighbor to (s, conf)
                        endif
                    endif    
                endif
           enddo             
        endif    
    enddo 

    ! print 
    if(.true.) then
        
        fname='phosphate_pairs.log'
        fname=trim(adjustl(fname))
        !     .. opening file
        open(unit=newunit(un_pp),file=fname)

        write(un_pp,*)"phosphate pairs"
        write(un_pp,*)"value ta=",ta, "value tPhos=",tPhos
        write(un_pp,*)"distphoscutoff=",distphoscutoff
        write(un_pp,*)"conf=",conf
        do s=1,nseg
             write(un_pp,*)s,type_of_monomer(s),nneigh(s,conf),(list_of_pairs(s,j),j=1,nneigh(s,conf))
        enddo
        close(un_pp)

    endif        

    ! allocate indexconfpair
    do s=1,nseg
        allocate(indexconfpair(s,conf)%elem(nneigh(s,conf)))
    enddo
    
    ! ..assign indexconfpair
    do s=1,nseg
        do j=1,nneigh(s,conf)
            indexconfpair(s,conf)%elem(j)=index_of_pairs(s,j)
        enddo 
    enddo 

    deallocate(list_of_pairs)
    deallocate(index_of_pairs)


end subroutine find_phosphate_pairs


subroutine write_indexconfpair(nseg,conf,tPhos,sqrDphoscutoff)
    
    use chains, only :  type_of_monomer,indexconfpair, indexconf
    use chains, only : nneigh, distphoscutoff
    use parameters, only : tA 
    use volume, only : delta, linearIndexFromCoordinate,coordinateFromLinearIndex
    use myutils, only : newunit
    
    integer, intent(in) :: nseg
    integer, intent(in) :: conf
    integer, intent(in) :: tPhos
    real(dp), intent(in) :: sqrDphoscutoff


    integer :: s, sprime, i, j
    integer :: xi, yi, zi , idx    
    character(len=100) :: fname
    integer :: un_pp
    character(len=10) ::istr
    
    write(istr,'(I4)')conf
    fname='phosphate_indexconfpair.'//trim(adjustl(istr))//'.log'
    !     .. opening file
    open(unit=newunit(un_pp),file=fname)

    write(un_pp,*)"phosphate pairs"
    write(un_pp,*)"value ta=",ta, "value tPhos=",tPhos
    write(un_pp,*)"distphoscutoff=",distphoscutoff
    write(un_pp,*)"conf=",conf
    
    do s=1,nseg
        write(un_pp,*)s,type_of_monomer(s),nneigh(s,conf),indexconf(s,conf)%elem(1),&
        (indexconfpair(s,conf)%elem(j),j=1,nneigh(s,conf)) 
    enddo
    
    write(un_pp,*)" "

    do s=1,nseg
        write(un_pp,'(I5,I5,I5,I8)',advance="no")s,type_of_monomer(s),nneigh(s,conf),indexconf(s,conf)%elem(1)
        do j=1,nneigh(s,conf)
            call coordinateFromLinearIndex(indexconfpair(s,conf)%elem(j),xi,yi,zi)
            write(un_pp,'(A, 3I3, A)',advance="no")" ",xi,yi,zi," "
        enddo 
        write(un_pp,*)   
    enddo

    close(un_pp)    

end subroutine write_indexconfpair


! Like subroutine find_phosphate_pairs
! Here phosphate monomers belonging to tripelt are exluded from being pairs
! pre  find_phosphate_triplets need to establihsed first : use ntriplet

subroutine find_phosphate_pairs_exclude_triplets(nseg,conf,tPhos,sqrDphoscutoff,chain,Lx,Ly,Lz)

    use chains, only :  type_of_monomer,indexconfpair
    use chains, only : nneigh, indexconfpair, distphoscutoff, ntriplet
    use parameters, only : tA 
    use parameters, only : pbc_chains
    use volume, only : delta, linearIndexFromCoordinate
    use myutils, only : newunit,LogUnit, error_handler, print_to_log
    
    integer, intent(in) :: nseg
    integer, intent(in) :: conf
    integer, intent(in) :: tPhos
    real(dp), intent(in) :: sqrDphoscutoff
    real(dp), intent(in) :: chain(3,nseg)
    real(dp), intent(in) :: Lx,Ly,Lz

    integer, parameter :: maxnneigh = 10

    integer :: s, sprime, i, j, xi, yi, zi , idx
    integer, dimension(:,:), allocatable :: list_of_pairs, index_of_pairs
    real(dp) :: sqrdist

    character(len=100) :: fname, text
    integer :: un_pp, info
    character(len=10) ::istr
    logical ::  pair_triplet_disjoined
    
    allocate(list_of_pairs(nseg,maxnneigh))
    allocate(index_of_pairs(nseg,maxnneigh))

    info=0

    do s=1,nseg 
        nneigh(s,conf)=0
        if(type_of_monomer(s)==tPhos) then ! tPhos equiv to ta which is not set yet 
            if(ntriplet(s,conf)==0) then 
                do sprime=1,nseg
                    if(type_of_monomer(sprime)==tPhos ) then
                        if(s/=sprime) then ! prevent s=sprime being counted as a pair
                            if(ntriplet(sprime,conf)==0) then  

                                sqrdist=0.0_dp
                                do i=1,3
                                    sqrdist=sqrdist+(chain(i,s)-chain(i,sprime))**2
                                enddo
    
                                if(sqrdist<=sqrDphoscutoff) then ! comparing square of distance to square of cutoff  
                                    ! accept s and sprime are a pair
                                    nneigh(s,conf)=nneigh(s,conf)+1
                                    list_of_pairs(s,nneigh(s,conf))=sprime ! temporarily storage of  segment number of neighbor to (s,conf)

                                    ! transforming form real- to lattice coordinates                 
                                    
                                    if(pbc_chains) then  
                                        xi = int(pbc(chain(1,sprime),Lx)/delta)+1
                                        yi = int(pbc(chain(2,sprime),Ly)/delta)+1
                                        zi = int(pbc(chain(3,sprime),Lz)/delta)+1
                                    else 
                                        xi = int(chain(1,sprime)/delta)+1
                                        yi = int(chain(2,sprime)/delta)+1
                                        zi = int(chain(3,sprime)/delta)+1
                                    endif    

                                    call linearIndexFromCoordinate(xi,yi,zi,idx)
                                    !idx = coordtoindex(xi,yi,zi) ! hash-table look up
                                    
                                    index_of_pairs(s,nneigh(s,conf))=idx  ! temporarily storage of index of neighbor to (s, conf)
                                endif
                            endif
                        endif 
                    endif
                enddo
            endif                    
        endif    
    enddo 

    ! print 
    if(.true.) then
        
        fname='phosphate_pairs.log'
        fname=trim(adjustl(fname))
        !     .. opening file
        open(unit=newunit(un_pp),file=fname)

        write(un_pp,*)"phosphate pairs"
        write(un_pp,*)"value ta=",ta, "value tPhos=",tPhos
        write(un_pp,*)"distphoscutoff=",distphoscutoff
        write(un_pp,*)"conf=",conf
        do s=1,nseg
             write(un_pp,*)s,type_of_monomer(s),nneigh(s,conf),(list_of_pairs(s,j),j=1,nneigh(s,conf))
        enddo
        close(un_pp)

    endif        

    ! allocate indexconfpair
    do s=1,nseg
        allocate(indexconfpair(s,conf)%elem(nneigh(s,conf)))
    enddo
    
    ! ..assign indexconfpair
    do s=1,nseg
        do j=1,nneigh(s,conf)
            indexconfpair(s,conf)%elem(j)=index_of_pairs(s,j)
        enddo 
    enddo 

    deallocate(list_of_pairs)
    deallocate(index_of_pairs)

    ! check that union of sets of pairs and triplet is empty 
    ! ntriplet(s,conf) nonozero then nneigh(s,conf) zero and reverse

    pair_triplet_disjoined=.true.
    do s=1,nseg
        if(ntriplet(s,conf)>0)then
            if(nneigh(s,conf)>0) pair_triplet_disjoined=.false.
        endif    
    enddo

    if(.not.pair_triplet_disjoined) then
        info=1
        text="Warning: union of sets of pairs and triplets is NOT empty"
        print*,text
        call print_to_log(LogUnit,text) 
        call error_handler(info,"find_phosphate_pairs_exclude_triplets")
    endif

end subroutine find_phosphate_pairs_exclude_triplets


! Finds triplets of phosphate pairs for given conformation number conf 
! Conformation is stored in chain
! input integer :: nseg : number atoms/segment
!       integer :: conf : conformation number
!       integer :: tPphos : number associated with type of phosphates
!       real(dp) :: sqrDphoscutoff : squared distance of cutoffdistance citeria for pair
!       real(dp) :: chain(3,nseg) : hold coordiante of backbone confomation for all atom/segment
!       real(dp) :: LX,Ly,Lz : dimension lattice/box in nm 


subroutine find_phosphate_triplets(nseg,conf,tPhos,sqrDphoscutoff,chain,Lx,Ly,Lz)

    use chains, only :  type_of_monomer,indexconfpair
    use chains, only : nneigh, indexconfpair, distphoscutoff
    use chains, only : ntriplet,indexconftriplet
    use parameters, only : tA 
    use parameters, only : pbc_chains
    use volume, only : delta, linearIndexFromCoordinate
    use myutils, only : newunit
    
    integer, intent(in) :: nseg
    integer, intent(in) :: conf
    integer, intent(in) :: tPhos
    real(dp), intent(in) :: sqrDphoscutoff
    real(dp), intent(in) :: chain(3,nseg)
    real(dp), intent(in) :: Lx,Ly,Lz

    ! local variables 
     integer, parameter :: maxntriplet = 1


    integer :: s, sprime, sdblprime, i, j, idx, idxdble
    real(dp) :: sqrdist, sqrdists, sqrdistsprime
    integer :: num_triplets, countseg
    integer :: xi, yi, zi 

    integer, dimension(:,:,:), allocatable ::  list_of_triplets
    integer, dimension(:,:,:), allocatable :: index_of_triplets

    allocate(list_of_triplets(nseg,maxntriplet,2))
    allocate(index_of_triplets(2,nseg,maxntriplet))
    
    ntriplet=0
    
    do s=1,nseg 
        if(type_of_monomer(s)==tPhos) then ! tPhos equiv to ta which is not set yet 
            do sprime=1,nseg
                if(type_of_monomer(sprime)==tPhos ) then
                    if(s/=sprime) then ! prevent s=sprime being counted as a pair
                       
                        sqrdist=0.0_dp
                        do i=1,3
                            sqrdist=sqrdist+(chain(i,s)-chain(i,sprime))**2
                        enddo
    
                        if(sqrdist<=sqrDphoscutoff) then ! comparing square of distance to square of cutoff  
                            ! accept s and sprime are a pair

                            do sdblprime=sprime+1,nseg ! sdbprime=1,nseg
                                if(type_of_monomer(sdblprime)==tPhos ) then
                                    if((s/=sdblprime).and.(sprime/=sdblprime )) then ! prevent sdblprime =s or s=sprime to count counted as a pair
                                        sqrdists=0.0_dp
                                        sqrdistsprime=0.0_dp
                                        do i=1,3
                                            sqrdists=sqrdists+(chain(i,sdblprime)-chain(i,s))**2
                                            sqrdistsprime=sqrdistsprime+(chain(i,sdblprime)-chain(i,sprime))**2
                                        enddo
                                        if((sqrdists<=sqrDphoscutoff) .and. (sqrdistsprime<=sqrDphoscutoff) ) then 
                                            ! accept sdlbprime as a pair with s and sprime 
                                            ! print*,"triplet =",s," ",sprime," ",sdblprime
                                            ntriplet(s,conf)= ntriplet(s,conf)+1

                                            list_of_triplets(s,ntriplet(s,conf),1) = sprime ! temporarily storage 
                                            list_of_triplets(s,ntriplet(s,conf),2) = sdblprime
                                            
                                            ! transforming form real- to lattice coordinates                 
                            
                                            if(pbc_chains) then  
                                                xi = int(pbc(chain(1,sprime),Lx)/delta)+1
                                                yi = int(pbc(chain(2,sprime),Ly)/delta)+1
                                                zi = int(pbc(chain(3,sprime),Lz)/delta)+1
                                            else 
                                                xi = int(chain(1,sprime)/delta)+1
                                                yi = int(chain(2,sprime)/delta)+1
                                                zi = int(chain(3,sprime)/delta)+1
                                            endif    

                                            call linearIndexFromCoordinate(xi,yi,zi,idx)
                                            !idx = coordtoindex(xi,yi,zi) ! hash-table look up
                                        
                                             if(pbc_chains) then  
                                                xi = int(pbc(chain(1,sdblprime),Lx)/delta)+1
                                                yi = int(pbc(chain(2,sdblprime),Ly)/delta)+1
                                                zi = int(pbc(chain(3,sdblprime),Lz)/delta)+1
                                            else 
                                                xi = int(chain(1,sdblprime)/delta)+1
                                                yi = int(chain(2,sdblprime)/delta)+1
                                                zi = int(chain(3,sdblprime)/delta)+1
                                            endif    

                                            call linearIndexFromCoordinate(xi,yi,zi,idxdble)
                                            !idx dble= coordtoindex(xi,yi,zi) ! hash-table look up
                                            
                                            index_of_triplets(1,s,ntriplet(s,conf))=idx     ! temporarily storage 
                                            index_of_triplets(2,s,ntriplet(s,conf))=idxdble  ! temporarily storage 
                                        endif
                                    endif
                                endif
                            enddo

                        endif

                    endif
                endif                        
            enddo
        endif
    enddo 

    ! print 
    if(.false.) then 
        num_triplets=0
        do s=1,nseg
            if(ntriplet(s,conf)>0) then
                num_triplets=num_triplets+ntriplet(s,conf)
                print*,"ntriplet(",s,conf,")=",ntriplet(s,conf)
            endif     
        enddo

        num_triplets = num_triplets/3
        print*,"number of triplets = ", num_triplets
    endif


    ! allocate indexconftriplet
    do s=1,nseg
        allocate(indexconftriplet(1,s,conf)%elem(ntriplet(s,conf)))
        allocate(indexconftriplet(2,s,conf)%elem(ntriplet(s,conf)))
    enddo
    
    ! ..assign indexconftriplet
    do s=1,nseg
        do j=1,ntriplet(s,conf)
            indexconftriplet(1,s,conf)%elem(j)=index_of_triplets(1,s,j)
            indexconftriplet(2,s,conf)%elem(j)=index_of_triplets(2,s,j)
        enddo 
    enddo 

    deallocate(list_of_triplets)
    deallocate(index_of_triplets)
    
end subroutine find_phosphate_triplets


! Compute index_phos and len_phos:
! index_phos is a list representing all latice element that contain phosphates
! len_phos : length of index_phos

subroutine find_phosphate_loc(index_phos,len_index_phos)
    
    use chains, only :  type_of_monomer,indexconf
    use globals, only : nseg, nsize, cuantas
    
    integer, allocatable, dimension(:), intent(inout) :: index_phos
    integer, intent(inout) :: len_index_Phos

    integer :: location_list(nsize) ! location_list(i) if > 0 then there is a or multiple phosphate located at latiice element i
    integer :: num_phos ! number of phophates
    integer :: cell_num, conf, s, i, k, tPhos
    integer :: index_phos_tmp(nsize)
    
    ! find location 

    tPhos=find_type_phosphate() ! tA is not set until set call to  init_dna !!!!
    
    location_list=0
    do conf=1, cuantas
        do s=1,nseg
            if(tPhos==type_of_monomer(s)) then
            
                cell_num = indexconf(s,conf)%elem(1)
                location_list(cell_num)=location_list(cell_num) +1
            endif
        enddo
    enddo
    
    k=0
    do i=1, nsize
        if(location_list(i)>0) then
            k=k+1
            index_phos_tmp(k)=i
        endif
    enddo               
    
    len_index_phos=k 

    allocate(index_phos(len_index_phos))
    
    do i=1, len_index_phos
        index_phos(i)= index_phos_tmp(i)
        ! write(222,*)index_phos(i)
    enddo

    do s=1,nseg
        if(tPhos==type_of_monomer(s)) then
            i = indexconf(s,1)%elem(1)
            ! write(111,*),i
        endif
    enddo
       
end subroutine find_phosphate_loc

! Compute index_phos and len_phos:
! index_phos is a list representing all latice element that contain phosphates
! inverse_index_phos : inverse 
! len_phos : length of index_phos

subroutine find_phosphate_location(index_phos,inverse_index_phos,len_index_phos)
    
    use chains, only :  type_of_monomer,indexconf
    use globals, only : nseg, nsize, cuantas
    
    integer, allocatable, dimension(:), intent(inout) :: index_phos
    integer, allocatable, dimension(:), intent(inout) :: inverse_index_phos
    integer, intent(inout) :: len_index_Phos

    integer :: location_list(nsize) ! location_list(i) if > 0 then there is a or multiple phosphate located at latiice element i
   ! integer :: num_phos ! number of phophates
    integer :: cell_num
    integer ::  conf, s, i, k, tPhos
    integer :: index_phos_tmp(nsize)
    
    allocate(inverse_index_phos(nsize))
    inverse_index_phos=0

    ! find location 

    tPhos=find_type_phosphate() ! tA is not set until by set call to  init_dna !!!!
    
    location_list=0
    do conf=1, cuantas
        do s=1,nseg
            if(tPhos==type_of_monomer(s)) then           
                cell_num = indexconf(s,conf)%elem(1)
                location_list(cell_num)=location_list(cell_num) +1
            endif
        enddo
    enddo
    
    k=0
    do i=1, nsize
        if(location_list(i)>0) then
            k=k+1
            index_phos_tmp(k)=i
            inverse_index_phos(i)=k
        endif
    enddo               
    
    len_index_phos=k 

    allocate(index_phos(len_index_phos))
    
    do i=1, len_index_phos
        index_phos(i)= index_phos_tmp(i)
    enddo

       
end subroutine find_phosphate_location


! print index coordinates of conf between sbegin and send ( i..e histone ) 
! for conformations conf_begin through conf_end

subroutine write_indexchain_histone(sbegin,send,conf_begin,conf_end)

    use globals, only : cuantas
    use chains,  only : indexchain
    use myutils, only : lenText, newunit, error_handler
    use volume, only : indextocoord

    integer, intent(in) :: sbegin,send, conf_begin
    integer, intent(inout) :: conf_end

    character(len=lenText) :: fname, istr, text
    integer :: s, conf, un, ios, info, ind
    integer :: ivec(3)

    info=0
    if(cuantas<conf_end) then
        print*,"write_indexchain: cuantas < conf_end"
        conf_end=cuantas
    endif
    if(conf_begin>cuantas)then 
        info=1
        text="conf_begin>cuantas: stopping"
        call error_handler(info,"write_indexchain")
    endif
    do conf=conf_begin,conf_end
       
        fname='index-chain.cuantas'
        write(istr,'(I4)')conf
        fname=trim(fname)//trim(adjustl(istr))//'.dat'
        open(unit=newunit(un),file=fname,status='new',iostat=ios)
        do s=sbegin,send
           ind=indexchain(s,conf)
           ivec=indextocoord(ind,:)
           write(un,*,iostat=ios)s,ind,ivec
        enddo
        close(un)
    enddo


end subroutine 

! Check indexchain of histone ( CA of AA) atoms,  which are between sbegin and send, that 
! are  translated and rotate in the same position and orientation for all  conformantion
 
subroutine compare_indexchain_histone(sbegin,send,info)

    use globals, only : cuantas
    use chains,  only : indexchain
   ! use mpivars, only : rank

    integer, intent(in) :: sbegin,send
    integer, intent(out) :: info

    integer :: s, conf,ind, indref

    info=0
    do conf=1,cuantas
        do s=sbegin,send
           indref= indexchain(s,1) ! assume conf 1 is oke  
           ind = indexchain(s,conf)
           if(ind/=indref) then 
              info=info+1
              print*,"rank=",0," conf=",conf,"s=",s,"ind=",ind," indref=",indref
           endif
        enddo
    enddo
    print*,"rank=",0," info=",info


end subroutine



! Check indexconf of histone ( CA of AA) atoms,  which are between sbegin and send, that 
! are  translated and rotate in the same position and orientation for all  conformantion

subroutine compare_indexconf_histone(sbegin,send,info)

    use globals, only : cuantas
    use chains,  only : indexconf, nelem
    !use mpivars, only : rank

    integer, intent(in) :: sbegin,send
    integer, intent(out) :: info

    integer :: s, conf, j, ind, indref

    info=0
    do conf=1,cuantas
        do s=sbegin,send
           do j=1,nelem(s)
              indref= indexconf(s,1)%elem(j)  ! assume conf 1 is oke  
              ind   = indexconf(s,conf)%elem(j)
              if(ind/=indref) then
                 info=info+1
                 print*,"rank=",0," conf=",conf," s=",s," j= ",j," ind=",ind," indref=",indref
              endif
           enddo
        enddo
    enddo

    print*,"rank=",0," info=",info


end subroutine

subroutine test_nmer_indexchain_histone(s0,s1,info)

    use globals, only : systype, nnucl

    integer, intent(inout) :: info     
    integer, intent(in) :: s0,s1


    call compare_indexchain_histone(s0,s1,info)
    print*,"info=",info
    info=-1 ! Warning

    if(systype=="nucl_ionbin_sv".or.systype=="nucl_neutral_sv".or.systype=="nucl_ionbin_Mg".or. &
       systype=="nucl_ionbin_MgA") then
       call compare_indexconf_histone(s0,s1,info)
       print*,"info=",info
       info=-1 ! Warning
    endif

end subroutine

subroutine test_index_histone(info) 
    
    use globals, only : systype, nnucl 
 
    integer, intent(inout) :: info

    integer :: s0,s1
    logical :: flag

    flag=.true.

    if(nnucl==3) then
       s0=2663
       s1=3494
    else if(nnucl==5) then 
       s0=4509
       s1=5340
    else if(nnucl==8) then
       s0=6355
       s1=7188  
    else  
        flag=.false.
        print*,"nnucl not equal to 3,5 ore 8: test_index_histone not applicable"
        info=-1 ! Warning  
    endif
    if(flag) call test_nmer_indexchain_histone(s0,s1,info)

end subroutine


function compute_cuantas_no_overlap(cuantas,no_overlapchain)result(num_no_overlap)
    
    integer, intent(in) :: cuantas 
    logical, intent(in) :: no_overlapchain(:)
    integer :: num_no_overlap


    integer :: num, conf

    num=0
    do conf=1,cuantas 
        if(no_overlapchain(conf)) num=num+1
    enddo    

    num_no_overlap = num

end function compute_cuantas_no_overlap


end module chaingenerator
