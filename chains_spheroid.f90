module chains_spheroid
  
    use precision_definition

    implicit none

    real(dp), dimension(:,:), allocatable :: chain_NCP                      ! cooridnate for elements of NCP particle    
    real(dp), dimension(:,:), allocatable :: chain_sph                      ! cooridnate for elements of complete nucleosome chain   


    ! chain stuctural quantities

    integer,  dimension(:), allocatable   :: segcm_sph                    ! monomer or units of chain closest to cm of histone  
    real(dp)                              :: Rgsqr_sph                    ! radius of gyration 
    real(dp)                              :: Rendsqr_sph                  ! end-to-end distance
    real(dp), dimension(:), allocatable   :: bond_angle_sph               ! bond angle
    real(dp), dimension(:), allocatable   :: dihedral_angle_sph           ! dihedralangle
    real(dp), dimension(:), allocatable   :: nucl_spacing_sph             ! spacing or distance between Nuclesome 
    real(dp), dimension(:,:), allocatable :: gyr_tensor_sph               ! asphericity matrix/gyration tensor 
    real(dp)                              :: Asphparam_sph                ! asphericity parameter  invariant of gyration tensor


    private

    public :: chain_NCP, init_chain_NCP, chain_sph, init_chain_sph  

    contains 

        subroutine init_chain_NCP(info)
            
            use globals, only : nseg

            integer, intent(out) :: info 

            integer :: ier
            
            info=0

            if (.not. allocated(chain_NCP)) then
                allocate(chain_NCP(3, nseg),stat=ier)
                if(ier/=0) info=ier
            endif
        
        end subroutine init_chain_NCP

        subroutine init_chain_sph(info)

              use globals, only : nseg_sph

            integer, intent(out) :: info 

            integer :: ier
            
            info=0

            if (.not. allocated(chain_sph)) then
                allocate(chain_sph(3, nseg_sph),stat=ier)
                if(ier/=0) info=ier
            endif
        
        end subroutine init_chain_sph

        subroutine allocate_chains_sph(nnucl) 

            integer, intent(in) :: nnucl

            ! chain_sph stuctural quantities
            allocate(segcm_sph(nnucl))
            allocate(bond_angle_sph(nnucl-2))
            allocate(dihedral_angle_sph(nnucl-3))
            allocate(nucl_spacing_sph(nnucl-1)) 
            allocate(gyr_tensor_sph(3,3)) 
  
        end subroutine allocate_chains_sph

end module chains_spheroid