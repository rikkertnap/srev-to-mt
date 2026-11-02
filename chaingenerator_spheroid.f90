! ----------------------------------------------------------|

module chaingenerator_spheroid

    use precision_definition

    implicit none

    character(len=80), parameter :: fmt3xyz = "(3ES15.5E2)"

    private
    public   :: init_spheroid, make_chain_spheroid_mc, read_SREV, build_NCP_chain

contains  

    subroutine init_spheroid(info)

        use globals, only : nseg_sph, nnucl_sph, lseg_sph, max_confor_sph, nseg
        use GB_potential, only : init_GB_const_defaults, init_GB_const, print_GB_const, init_GB_general_const
        use VdW_potential, only : read_GBinputfile !, EGB_threshold
       
        integer, intent(out) :: info 

        info = 0 
        
        ! init 
        lseg_sph = 10.6_dp
        max_confor_sph = 100

        nseg_sph = nseg * nnucl_sph 

        ! set default for GB potential
        call init_GB_const_defaults()
        ! overide default 
        call read_GBinputfile(info)
        ! compute dependent auxilary functions
        call init_GB_general_const() ! check which one ???
        ! call init_GB_const()

        call print_GB_const()
    
        call make_seed_random_unit_vec

    end subroutine init_spheroid


    function unit_vec_3d(theta,phi) result(unitvector)

        real(dp), intent(in)  :: theta, phi 
        real(dp)  :: unitvector(3) 

        unitvector(1) = sin(theta) * cos(phi) 
        unitvector(2) = sin(theta) * sin(phi)
        unitvector(3) = cos(theta) 

    end function  unit_vec_3d


    subroutine make_seed_random_unit_vec

        integer :: state_size
        integer, allocatable, dimension(:) :: state
 
        ! init random number generator
        call random_seed( size=state_size )
        allocate(state(state_size))
        state = 20180815
        call random_seed( put=state )

    end subroutine make_seed_random_unit_vec

    ! generate s sequence of random oriented unit vectors
    ! input unvector(3,n)
    subroutine make_random_unit_vec(uvector)

        use mathconst 

        real(dp), intent(inout)  :: uvector(:,:) 

        integer :: i, nmer
        real(dp) :: theta, phi, ran
        real(dp) :: unvec(3)

        nmer =size(uvector,2)

        do i = 1, nmer
            call random_number(ran) ! random number between 0<= ran < 1
            theta = acos(1.0_dp-2.0_dp * ran) 
            call random_number(ran) 
            phi = ran * 2.0_dp * pi
            unvec = unit_vec_3d(theta,phi)
            uvector(:,i) = unvec 
        enddo    

    end subroutine make_random_unit_vec


    subroutine orientation_spheriod_chain(chains,nmer)

        use VdW_potential, only  : write_chain_com_unitvec_lammps_trj, open_chain_com_unitvec_lammps_trj
        use VdW_potential, only  : GBenergyeff

        real(dp), intent(in) :: chains(:,:)
        integer, intent(in) :: nmer 

        ! .. local variables

        real(dp) :: uvector(3,nmer)
        integer :: info , un_traj_com
        real(dp) :: GBEnergy
        logical :: no_overlap

        call make_random_unit_vec(uvector )
        GBEnergy = GBenergyeff(nmer,chains,uvector,no_overlap)

        info = 0
        
        ! output conformation
    
        un_traj_com = open_chain_com_unitvec_lammps_trj(info,"sph")
        
        call write_chain_com_unitvec_lammps_trj(un_traj_com,chains,uvector)
        close(un_traj_com)
        
        ! print*,"GBEnergy=",GBEnergy," nooverlap=",no_overlap

    end subroutine orientation_spheriod_chain



    subroutine make_chain_spheroid_mc

        use globals, only : nnucl_sph, lseg_sph, max_confor_sph
        use random, only : seed
        use parameters, only : write_mc_chains, pbc_chains
        use volume, only : nx, ny, nz, delta, pbc
        use myutils, only : print_to_log, lenText, LogUnit
        use cadenas_linear, only : make_linear_chains
        use chaingenerator, only : conf_write_chain, open_chain_energy, open_chain_lammps_trj
        use chaingenerator, only : write_chain_lammps_trj


        !  .. variable and constant declaractions      

        integer :: j,s              ! dummy indices
        integer :: nchains           ! number of rotations
        integer :: maxnchains        ! number of rotations
        integer :: conf              ! counts number of conformations
        real(dp) :: chain(3,nnucl_sph,200) ! chain(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
        real(dp) :: chain_rot(3,nnucl_sph), x(nnucl_sph), y(nnucl_sph), z(nnucl_sph) ! coordinates
        real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm  ! sizes box
        character(len=lenText) :: text
        integer  :: un_trj, un_ene
        integer :: info    
        real(dp) :: energy 

        !  .. executable statements
        !  .. initializations of variables     
        
        conf = 1                 ! counter for conformations
        seed = 435672            ! seed for random number generator  
      
        Lz = nz*delta            ! maximum height box 
        Lx = nx*delta            ! maximum width box 
        Ly = ny*delta            ! maximum depth box 
        xcm= Lx/2.0_dp           ! center box
        ycm= Ly/2.0_dp
        zcm= Lz/2.0_dp
    
        energy=0.0_dp

        if(write_mc_chains) then 
            conf_write_chain=0
            un_trj = open_chain_lammps_trj(info,"sph")
            un_ene = open_chain_energy()
        endif   

        do while (conf<=max_confor_sph)

            nchains= 0      ! init zero 
            call make_linear_chains(chain,nchains,maxnchains,nnucl_sph,lseg_sph) ! chain generator f90
            
            if(write_mc_chains) then 
                energy=0.0_dp
                write(un_ene,*)energy  
                call write_chain_lammps_trj(un_trj,chain,nchains)  
            endif    
    
            
            do j=1,nchains   
                
                do s=1,nnucl_sph                          !  transforming form real- to lattice coordinates
                    z(s) = chain(1,s,j)+zcm
                    x(s) = chain(2,s,j)+xcm
                    y(s) = chain(3,s,j)+ycm
                enddo
                
                do s=1,nnucl_sph

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

                enddo   

                call orientation_spheriod_chain(chain_rot,nnucl_sph)            
        
                conf = conf +1 

            enddo             ! end loop over nchains

        enddo                     ! end while loop
        
        !  .. end chains generation 
    
        text='Spheroid Chains generated'
        call print_to_log(LogUnit,text)

        if(write_mc_chains) then
            close(un_trj)
            close(un_ene)
        endif   
    
    end subroutine


    ! reed file containing center of mass/center of volume from a SREV chain
    ! file named trajSREV.0.xyz
    ! format file : first line <nnucl_sph>
    !               second line a comment
    !               remnaining lines: x y z      

    subroutine read_SREV(rcom,info)

        use globals, only : nnucl_sph
        use myio, only : myio_err_chainsfile, myio_err_nnucl_sph
        use myutils, only : lenText,LogUnit, print_to_log, newunit
    
        real(dp), intent(inout) :: rcom(:,:)
        integer, intent(inout) :: info

        ! .. local varaibles

        character(len=12) :: fname
        integer :: ios, un
        logical :: exist, isReadGood
        character(len=lenText) text
        integer :: nnucl_sph_file, s
        real(dp) :: xc, yc, zc 

        ! .. executable statements   

        info=0

        ! .. open file   

        fname='trajSREV.xyz'
            
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un),file=fname,status='old',iostat=ios)
        else
            print*,'traj file :',fname,' does not exit'
            info = myio_err_chainsfile
            return
        endif
        if(ios >0 ) then
            print*, 'Error opening file : iostat =', ios
            info = myio_err_chainsfile
            return
        endif

        ! .. read preamble
        read(un,*,iostat=ios)nnucl_sph_file
        if(ios/=0) isReadGood=.false.
        read(un,*,iostat=ios) ! skip line
        if(ios/=0) isReadGood=.false.   
        
        if(nnucl_sph_file.ne.nnucl_sph) then
            text="nnucl_sph chain file not equal internal nnucl_sph : stop program"
            print*,"nnucl_sph_file=",nnucl_sph_file," nnucl_sph=",nnucl_sph
            call print_to_log(LogUnit,text)
            info=myio_err_nnucl_sph
            return
        endif

        do s=1,nnucl_sph  ! .. read trajectory coordinates

            !read(un,fmt3xyz,iostat=ios)xc,yc,zc
            read(un,*,iostat=ios)xc,yc,zc ! read unformatted file
            if(ios/=0) isReadGood=.false.
        
            rcom(1,s) = xc 
            rcom(2,s) = yc 
            rcom(3,s) = zc 
        !    print*,xc,yc,zc
        enddo

        close(un) 

    end subroutine read_SREV


    ! Make a nucleosome chain of nnucl_sph unconnected NCPs.
    ! pre: Assumes one NCP equals  a 1-mer chain read by subroutine read_chains_xyz_nucl_volume in module chaingenerator
    ! read_chains_xyz_nucl_volume generates 
    ! chain_NCP(:,:) and  unitvector_triplets(:)  needed to build the NPC chains
    ! chain_NCP(:,:) holds the 1-mer chain conformation 
    ! unitvector_triplets(:) =atom_id used to make unit vector and com of NCP must be set 
    ! nnucl for read n-mer must be n=1 !!
    ! post: 1. xyx confomration file holding  nnucl_sph unconnected NCPs. 
    !       2. types squence file
    !       3. orient element 
    !       4. atom id of 'com' of of each NCP       


    subroutine build_NCP_chain(info)

        use globals, only : nseg, nnucl_sph, nseg_sph, nnucl
        use chains_spheroid, only : chain_NCP, init_chain_sph, chain_sph
        use myutils, only : error_handler
        use chains, only : unitvector_triplets 
        ! hold the unit numbers of the triplet of AA that is used to build the com 
        use VdW_potential, only : make_com_unit_vector_nucl_rotation, GBenergyeff
        use VdW_potential, only  : write_chain_com_unitvec_lammps_trj, open_chain_com_unitvec_lammps_trj
        use chains, only : type_of_monomer
        use chain_rotation, only : rotate_chain
        use chaingenerator, only : open_chain_xyz, write_chain_xyz, write_chain_lammps_trj
    

        integer, intent(out) :: info

        integer :: s, n, un_traj                ! dummy variable
        integer ::  nseg_shift                  ! shift in index for assignmen tof chain_spj
        real(dp) :: rcom_NCP(3,1)               ! com of NCP chain
        real(dp) :: uvector_NCP(3,1)            ! oreintation vector of NCP chain
        real(dp) :: rcom_sph(3,nnucl_sph)       ! com of each nucloesome of chain
        real(dp) :: uvector_sph(3,nnucl_sph)    ! oreintation vector of each nucloesome of chain
        real(dp) :: chain_NCP_rot(3,nseg)       ! rotated chain_NCP
        real(dp) :: chain_sph_xyz(3,nseg_sph,1) ! outputed NCP chains 
        logical :: DEBUG_local, no_overlap
        integer :: type_of_monomer_loc(nseg_sph) 
        real(dp) :: GBEnergy

        DEBUG_local = .true. 
 
        info = 0
        print*,"nseg_sph=",nseg_sph

        ! 0. 
        if(nnucl/=1) then 
            info=1
            print*,"Read n-mer not a NCP: nnucl=",nnucl
            return 
        endif 

        ! 1.  make data variable to hold new nucleosome chain of nnucl_sph 
        
       !  nseg_sph = nseg * nnucl_sph 

        call init_chain_sph(info)  
        call error_handler(info,"init_chain_sph")

        ! 2. get com and orientation of chain_NCP
        !
        call make_com_unit_vector_nucl_rotation(chain_NCP,nnucl,unitvector_triplets,rcom_NCP,uvector_NCP)
       
        print*,"rcom_NCP=", rcom_NCP
        print*,"uvector_NCP=", uvector_NCP

        ! call make_com_unit_vector_nucl_simpleCOM(chain_NCP,nnucl,segcm,segunitvector,rcom_NCP,uvector_NCP)
        ! call make_com_unit_vector_nucl_simple(chain_NCP,nnucl,segcm,segunitvector,rcom_NCP,uvector_NCP)

        ! 2.a. translate to make com to origin
        do s=1,nseg
            chain_NCP(:,s) = chain_NCP(:,s) - rcom_NCP(:,1)
        enddo

        ! 3. get SREV COM over every nucleosome 

        call read_SREV(rcom_sph, info)
        call error_handler(info,"read_SREV")

        ! 4. make orientation of every nucleosome (rotation)
        
        call make_random_unit_vec(uvector_sph)
  
        ! 5. loop over all nnucl_sph  

        do n=1,nnucl_sph 
            
            ! 5.a. rotate each nucleosome NCP 
            call rotate_chain(chain_NCP,chain_NCP_rot,uvector_NCP,uvector_sph(:,n))

            ! 5.b. translate each nucleoseme into com of SREV

            nseg_shift = nseg*(n-1)
            do s=1,nseg
                chain_sph(:,s+nseg_shift) = chain_NCP_rot(:,s)+rcom_sph(:,n)
            enddo
        
        enddo  

        ! 6. 
        
        GBEnergy = GBenergyeff(nnucl_sph,rcom_sph,uvector_sph,no_overlap)
        print*,"GBEnergy=",GBEnergy," nooverlap=",no_overlap
 
        ! 7. output xyz trajectory 
        
        chain_sph_xyz(:,:,1) = chain_sph
        un_traj = open_chain_xyz(info,"sph")
        call write_chain_xyz(un_traj,chain_sph_xyz,1)
        close(un_traj)

    
        ! 7. output of auxilary file types file and orient files 

        call make_types_file_NCP_chain(info)
        call make_segcm_file_NCP_chain(info)
        call make_nuc_orient_file_NCP_chain(info)

        ! 9. output lammps trajectory

         do n=1,nnucl_sph
            do s=1,nseg
                type_of_monomer_loc(s+nseg*(n-1))=type_of_monomer(s)
            enddo
        enddo 

        un_traj=open_outfile(info,"trajchainelem_sph.lammpstrj")
        call write_chain_lammps_trj(un_traj,chain_sph_xyz,1, type_of_monomer_loc)
        close(un_traj)

        ! 10. output com conformation
    
        un_traj= open_chain_com_unitvec_lammps_trj(info,"sph")
        call write_chain_com_unitvec_lammps_trj(un_traj, rcom_sph ,uvector_sph)
        close(un_traj)

    end subroutine build_NCP_chain 

    subroutine make_types_file_NCP_chain(info)

        use chains, only : type_of_monomer, type_of_monomer_char
        use globals, only : nseg, nnucl_sph

        integer, intent(out) :: info

        character(len=12) ::  filename
        integer :: un, s, n, ios

        info = 0 
        filename="types_NCP.in"
        un = open_outfile(info,filename)

        do n=1,nnucl_sph
            do s=1,nseg
                write(un,*,iostat=ios)type_of_monomer(s),type_of_monomer_char(s)
            enddo
        enddo

        if(ios/=0) info=ios
        close(un)

    end subroutine make_types_file_NCP_chain

    subroutine make_segcm_file_NCP_chain(info)

        use globals, only : nnucl_sph, nseg
        use chains, only : segcm

        integer, intent(out) :: info

        character(len=12) :: filename
        integer :: un, n, ios
        integer :: segcm_sph 

        info = 0 
        filename="segcm_NCP.in"
        un = open_outfile(info,filename)

        do n=1,nnucl_sph
            segcm_sph = segcm(1)+nseg*(n-1)
            write(un,*,iostat=ios) segcm_sph
        enddo

        if(ios/=0) info=ios
        close(un)

    end subroutine make_segcm_file_NCP_chain


    subroutine make_nuc_orient_file_NCP_chain(info)

        use globals, only : nnucl_sph, nseg
        use chains, only : orientation_triplets

        integer, intent(out) :: info

        character(len=17) :: filename
        integer :: un, n, ios
        real(dp) :: triplet_NCP(3), triplet_shift(3), triplet_sph(3)

        info = 0 
        filename="nuc_orient_NCP.in"
        un = open_outfile(info,filename)

        ! orientation triplet of NCP , only one
        triplet_NCP = orientation_triplets(1,:) ! shape(orientation_triplets)=(nnucl,3)   
        
        triplet_sph = triplet_NCP ! first 
        triplet_shift(1) = nseg   ! stride in atom ID
        triplet_shift(2) = nseg
        triplet_shift(3) = nseg

        do n=1,nnucl_sph
            write(un,*,iostat=ios)triplet_sph
            triplet_sph = triplet_sph + triplet_shift
        enddo

        if(ios/=0) info=ios
        close(un)

    end subroutine make_nuc_orient_file_NCP_chain




    ! Creates and opens file 
    ! input/outpout:     integer : info : value 0 no errors
    !                    character(lne=*): filename 
    ! output/return:     integer : un  : unit nubmer of file

    function open_outfile(info,filename)result(un)

        use myutils, only : newunit, lenText
        use myio, only : myio_err_file_exist, myio_err_outputfile 

        integer, intent(inout) :: info
        character(len=*), intent(in) :: filename  

        integer :: un

        ! local
        character(len=30) :: fname
        integer :: ios
        logical :: exist
        
        info = 0    

        fname=trim(filename)
        inquire(file=fname,exist=exist)
        if(.not.exist) then
            open(unit=newunit(un),file=fname,status='new',action="write",iostat=ios)
            if(ios > 0 ) then
                print*, 'Error opening : ',fname,' file : iostat =', ios
                info = myio_err_file_exist
                return
            endif
        else
            open(unit=newunit(un),file=fname,position="append",status='old',action="write",iostat=ios)
            if(ios > 0 ) then
                print*, 'Error opening : ',fname,' file : iostat =', ios
                info = myio_err_outputfile
                return
            endif
        endif   

    end function



     
end module  chaingenerator_spheroid