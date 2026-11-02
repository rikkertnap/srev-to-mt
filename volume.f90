!     makes volume elements 

module volume      

    use precision_definition
!    use mpivars

    implicit none
  
    !     .. variables
   
    real(dp) :: delta               ! delta  spacing of lattice site in x-, y- and z-direction  
    integer  :: nz                  ! nz number of lattice sites in z-direction 
    integer  :: nx                  ! nx number of lattice sites in x-direction 
    integer  :: ny                  ! ny number of lattice sites in y-direction 
    integer  :: nsurf               ! nsurf number of lattice site  at z=0 ore z=nz*delta
    real(dp) :: volcell             ! volcell=delta**3
    real(dp) :: areacell            ! areacell=delta**2
    real(dp) :: areasurf            ! area surfaces spanned in x- and y- direction
    real(dp) :: gamma               ! angle between oblique basis vectors u and v: cubic = gamma=90=pi/2 hexagonl gamma=60=2pi/3                          
    real(dp) :: beta                ! related beta = (pi/2- gamma)/2, angle between basis vector u and x and v and y
    real(dp) :: cos_two_beta        ! sqrt(cos(beta)**2 - sin(beta)**2)=cos(2beta) scales u and v coordinates 
    real(dp) :: sin_two_beta        ! sin(2beta)  
    character(len=11) :: geometry

    ! variable for grafting position

    ! integer :: sgraftpts(3)         ! triplet of units number of histone that is rotated into fixed orientation
    ! integer, dimension(:), allocatable :: sRg ! unit numbers of AA in histone close to cm of : number of number = nnucl

    ! hash table
    integer, dimension(:,:,:), allocatable :: coordtoindex 
    integer, dimension(:,:),  allocatable  :: indextocoord

   integer, dimension(:,:), allocatable   :: indexneighbor
!    integer, dimension(:,:), allocatable   :: inverse_indexneighbor
    integer, dimension(:,:), allocatable   :: inverse_indexneighbor_phos

    private
    
    public :: delta,nx,ny,nz,volcell,areacell,geometry, nsurf
    public :: gamma,cos_two_beta, sin_two_beta
    public :: indexneighbor !, inverse_indexneighbor, 
    public :: inverse_indexneighbor_phos
    public :: coordtoindex,indextocoord
    public :: coordinateFromLinearIndex, linearIndexFromCoordinate
    public :: xt, yt, ut, vt, ipbc, pbc
    public :: make_geometry

    public :: allocate_index_neighbors_phos, make_table_index_neighbors_phos

contains

    subroutine init_lattice

        use globals, only : nsize, systype
        use mathconst
        use chains, only : distphoscutoff ,maxneigh

        implicit none 

        integer :: rangecutoff


        if(geometry=="cubic") then 
            gamma = pi/2.0_dp
        else
            gamma = gamma * pi /180.0_dp ! convert from degree to radians
        endif        

        beta = (pi/2.0_dp - gamma) / 2.0_dp       ! used by ut, vt, xt and yt functions
        cos_two_beta=cos(2.0_dp*beta) ! sqrt(cos(beta)**2 - sin(beta)**2)  ! scaling of u and v coordinates
        sin_two_beta=sin(2.0_dp*beta)

        ! cubic lattice  or prism surface in x-y direction at z=0 and z=nz  
        ! nz=nzmax
        nsize = nx*ny*nz                             ! total number of cells or layers
        volcell = delta * delta * delta * 1.0_dp     ! volume of one latice volume 
        areacell = delta * delta
        nsurf = nx * ny  

        areasurf=nsurf*delta*delta
      
        ! use of hash table 
        call allocate_hashtable(nx,ny,nz)
        call make_hashtable()

        ! this should go somewhere else
 
        if(systype=="nucl_ionbin_Mg".or.systype=="nucl_ionbin_MgA".or.systype=="nucl_ionbin_Fe") then
            rangecutoff=int(distphoscutoff/delta)+2 
            maxneigh = (2*rangecutoff+1)**3
            !call allocate_index_neighbors(maxneigh)
            !call make_table_index_neighbors(rangecutoff)
        endif    

    end subroutine init_lattice
         

    subroutine  make_geometry()
    
        use globals
        use myutils

        implicit none
    
        character(len=lenText) :: text, str 

        call init_lattice

        write(str,'(A20)')geometry
        text="geometry="//trim(adjustl(str)) 
        call print_to_log(LogUnit,text) 


    end subroutine make_geometry

    subroutine linearIndexFromCoordinate(x,y,z,idx)
     
      implicit none 
      
        integer, intent(in)  :: x,y,z
        integer, intent(out) :: idx
      
        integer :: a,b,c,d

        a = 1
        b = nx 
        c = nx * ny 
        d = 1-a-b-c
        idx = a*x + b*y + c*z + d

    end subroutine linearIndexFromCoordinate


    subroutine coordinateFromLinearIndex(idx, x,y, z)

        implicit none

        integer, intent(out)  :: x,y,z
        integer, intent(in)   :: idx
        integer :: idxtmp
   
        idxtmp=idx
        x =  mod(idxtmp-1,nx)+1
        idxtmp =int((idxtmp-1)/nx)+1
        y = mod(idxtmp-1,ny)+1
        idxtmp = int((idxtmp-1)/ny)+1
        z = idxtmp
    
    end subroutine coordinateFromLinearIndex

     
    function mirror_index(idx,nz) result(idx_mirror)

        implicit none 
      
        integer, intent(in)  :: idx, nz
        integer :: idx_mirror

        integer :: x,y,z

        call coordinateFromLinearIndex(idx,x,y,z)
        call LinearIndexFromCoordinate(x,y,nz+1-z,idx_mirror)

    end function mirror_index

    ! coordinate transformation from cartesian x,y to prism or oblique u,v coordinates
    ! u = s(cos(gamma)*x -sin(gamma)*y); s=1/srqt(cos(gamma)^2-sin(gamma)^2)
    ! v = s(-sin(gamma)*x +sin(gamma)*y); 
    ! metric ds^2= s^2 du^2 + s^2 dv^2 + ss^2 sin(2 gamma) du dv +dz^2

    function ut(x, y) result(ut_val)
       
        real(dp), intent(in) :: x, y
        real(dp) :: ut_val
        
        ut_val = cos(beta)*x - sin(beta)*y
        ut_val = ut_val / cos_two_beta
    
    end function ut


    ! coordinate transformation from cartesian x,y to prism u,v  coordinates

    function vt(x, y) result(vt_val)
    
        real(dp), intent(in) :: x, y
        real(dp) :: vt_val
       
        vt_val = -sin(beta)*x + cos(beta)*y
        vt_val = vt_val / cos_two_beta
    
    end function vt



    ! inverse coordinate transformation from prism u,v to  cartesian x,y  coordinates

    function xt(u, v) result(xt_val)
        
        real(dp), intent(in) :: u, v
        real(dp) :: xt_val
    
        xt_val = cos(beta)*u + sin(beta)*v
        xt_val = xt_val / cos_two_beta
    
    end function xt

    ! inverse coordinate transformation from prism u,v to  cartesian x,y  coordinates
    
    function yt(u, v) result(yt_val)

        real(dp), intent(in) :: u, v
        real(dp) :: yt_val
    
        yt_val = sin(beta)*u + cos(beta)*v
        yt_val = yt_val / cos_two_beta
    
    end function yt
 
    function isEven(number) result(val) 
        integer, intent(in) :: number
        logical :: val

        if ( mod(number, 2) ==0) then
            val= .True. 
        else
            val= .False.
        endif 
   
    end function

    subroutine allocate_hashtable(nx,ny,nz)

        use globals, only : nsize
        integer, intent(in) :: nx,ny,nz

        allocate(coordtoindex(nx,ny,nz))
        allocate(indextocoord(nsize,3))
        
    end subroutine allocate_hashtable

    subroutine make_hashtable

        use globals, only : nsize

        integer :: idx, ix, iy, iz
            
        do idx=1,nsize
            
            call coordinateFromLinearIndex(idx, ix, iy,iz)

            coordtoindex(ix,iy,iz)=idx
            
            indextocoord(idx,1)=ix
            indextocoord(idx,2)=iy
            indextocoord(idx,3)=iz

        enddo
        
    end subroutine make_hashtable

    ! compute periodic boundary condition in integer units 
    ! real(dp) version  : pbbc is located in chaingenerator  

    function ipbc(ival,imax) result(intpbc)
        implicit none 
        integer, intent(in) :: ival
        integer, intent(in) :: imax
        integer :: intpbc

        if(ival>0) then
            intpbc=ival-int((ival-1)/imax)*imax
        else
            intpbc=ival-(int((ival-1)/imax)-1)*imax
        endif

    end function
         
   
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


!   subroutine allocate_index_neighbors(maxneigh)
!
!        use globals, only : nsize
!
!       integer, intent(in) :: maxneigh
!
!        allocate(indexneighbor(nsize,maxneigh))
!        allocate(inverse_indexneighbor(nsize,nsize))
!        
!    end subroutine allocate_index_neighbors

    
    ! Calculate indexneighbor(idx,k) 
    ! indexneighbor(idx,k) = return lattice cell index of neighbor number k of lattice cell index idx
    ! range : 1<=idx<=nsize 
    !       : 1<=k<=(2*rangecutoff+1)**3 
    !       : k=1 is same lattice cell         

    subroutine make_table_index_neighbors(rangecutoff)

        use globals, only : nsize

        integer, intent(in) :: rangecutoff
        
        integer :: idx, i, j, k, ix, iy, iz, deltaix,deltaiy, deltaiz
        integer :: neighbornumber, idxneigh
        integer :: ip, jp, kp


        do idx=1,nsize

            ! use hash table to get coordiantes 
            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            neighbornumber=0 ! 

            do deltaix=-rangecutoff,rangecutoff
                i=ix+deltaix
                do deltaiy=-rangecutoff,rangecutoff
                    j=iy+deltaiy
                    do deltaiz=-rangecutoff,rangecutoff
                        k=iz+deltaiz
                        ! apply pbc 
                        ip=ipbc(i,nx)
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)
                        neighbornumber=neighbornumber+1
                        idxneigh=coordtoindex(ip,jp,kp) ! index of neighbour
                        indexneighbor(idx,neighbornumber) = idxneigh
                    !    inverse_indexneighbor(idx,idxneigh) = neighbornumber
                    enddo
                enddo
            enddo

        enddo    

    end subroutine make_table_index_neighbors


    ! Stores in a memory efficient way inverse_index_neighbors
    !  
    ! inverse_indexneighbor(ind,k) evaluated for all lattice positions  
    ! inverse_indexneighbor_phos only evaluate for those position that contain a phosphate   
    ! inverse_indexneighbor_phos(ind,k): argument ind is different. 
    ! Here "ind" is ind th element of array that enumarated 
    ! all different phosphate postions 
    ! The assocaited lattice number = index_phos(ind)  ! give the lattice location 

    subroutine allocate_index_neighbors_phos(maxneigh,len_index_phos)

        use globals, only : nsize

        integer, intent(in) :: maxneigh
        integer, intent(in) :: len_index_phos

        allocate(indexneighbor(nsize,maxneigh))
        allocate(inverse_indexneighbor_phos(len_index_phos,nsize))

    end subroutine allocate_index_neighbors_phos



    ! Calculate indexneighbor(idx,k) and inverse_indexneighbor(idx 
    ! indexneighbor(idx,k) = return lattice cell index of neighbor number k of lattice cell index idx
    ! range : 1<=idx<=len_phos 
    !       : 1<=k<=(2*rangecutoff+1)**3 
    !       : k=1 is same lattice cell         

    subroutine make_table_index_neighbors_phos(rangecutoff,len_index_phos,index_phos)

        integer, intent(in) :: rangecutoff
        integer, intent(in) :: len_index_phos
        integer, intent(in) :: index_phos(len_index_phos)

        integer :: idx, i, j, k, ix, iy, iz, deltaix,deltaiy, deltaiz
        integer :: neighbornumber, idxneigh
        integer :: ip, jp, kp, n

        do n=1,len_index_phos    ! loop over element of unique phospate lattice positions 
            
            idx = index_phos(n)  ! give the lattice number of list in unique phosphate 

            ! use hash table to get coordinates 
            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            neighbornumber=0 ! 

            do deltaix=-rangecutoff,rangecutoff
                i=ix+deltaix
                do deltaiy=-rangecutoff,rangecutoff
                    j=iy+deltaiy
                    do deltaiz=-rangecutoff,rangecutoff
                        k=iz+deltaiz
                        ! apply pbc 
                        ip=ipbc(i,nx)
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)
                        neighbornumber=neighbornumber+1
                        idxneigh=coordtoindex(ip,jp,kp) ! index of neighbor
                        indexneighbor(idx,neighbornumber) = idxneigh
                        inverse_indexneighbor_phos(n,idxneigh) = neighbornumber

                    enddo
                enddo
            enddo

        enddo

    end subroutine make_table_index_neighbors_phos

end module volume
  
