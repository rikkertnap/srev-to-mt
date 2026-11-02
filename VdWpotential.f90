! Module to compute VdW/LJ/GB potential energy for nucleosome conformation
! 
! Modudle has following function and subroutines
! 
! function VdWpotentialenergy_MC(chain,nchains)result(Energy)
! function VdWpotentialenergy_SAW(chain)result(Energy)
! function VdWpotentialenergy(chain)result(Energy)

! function LJpotentialenergy(chain) result(Energy)
! function LJenergyeffective(chain,nmer)result(Energy)

! subroutine init_GBenergyeffective(segcm,nmer,segnumAAstartGBcom, segnumAAendGBcom) 
! function  GBenergyeffective_comb(chain,nmer)result(Energy)
! function  GBenergyeff(nmer,rcom,uvector)result(Energy)

! subroutine make_com_unit_vector_nucl_simple(chain,nmer,segcm,segunitvector,rcom,uvector)
! subroutine make_com_unit_vector_nucl_simpleCOM(chain,nmer,segcm,segunitvector,rcom,uvector)
! function com_chains(chains,nmer)
! subroutine parameter_com_ref(vec,s,t,vec_origin)
! subroutine parameter_com(vec,s_ref,t_ref,vec_origin)
! function normal_vector(vec)result(norm_vector)
! subroutine make_com_unit_vector_nucl_rotation(chain,nmer,unitvector_triplets,rcom,uvector)
! subroutine write_chain_com_unitvec_lammps_trj(un_trj,segcom,segunitvector)
! function open_chain_com_unitvec_lammps_trj(info)result(un_trj)

module VdW_potential     

    use precision_definition
  
    implicit none

    integer :: conf_write_com = 0         ! counter for write_chain_com_unitvec_lammps_trj
    !integer :: un_traj_com                ! unit number for file 
    
    ! auxilary variable used in determine COM and unitvectors 

    integer, allocatable :: segnumAAstartGBcom(:), segnumAAendGBcom(:) ! local variable 
    real(dp) :: vec_ref(3,3)              ! reference orientation vectors spanning plane through origin
                                          ! first index coordiantes second index number 
    real(dp) :: unit_vec_ref(3)           ! init unitvector from python program pca  
    integer  :: atom_id_unit_relative(3)  ! relative to atomid/segnumber of first AA

    integer :: deltaseg
    integer, allocatable :: segunitvector(:) 

    ! 
    real(dp) :: EGB_threshold = 0.5_dp 

    private :: dotproduct, pbc
    !private :: make_com_unit_vector_nucl_simple !,make_com_unit_vector_nucl_simpleCOM 
    ! private :: make_com_unit_vector_nucl_rotation
    private :: parameter_com_ref, parameter_com, normal_vector
    !private :: un_traj_com
    private :: segnumAAstartGBcom, segnumAAendGBcom
    private :: vec_ref, unit_vec_ref, atom_id_unit_relative
    private :: EGB_threshold 

contains 

    
    function VdWpotentialenergy_MC(chain,nchains)result(Energy)

        use globals, only : nseg
        use myutils, only : newunit, lenText
        use volume, only : delta, nx, ny, nz, coordinateFromLinearIndex
        use chains, only : type_of_monomer
        use parameters, only :  lsegAA,VdWeps

        real(dp), intent(in) :: chain(:,:,:)
        integer, intent(in) :: nchains

        real(dp) :: Energy(nchains)
        real(dp) :: Ene,sqrlseg,sqrdist
        integer :: k,i,j,s,t
        real(dp) :: xi,xj,yi,yj,zi,zj
        real(dp) :: Lz,Ly,Lx

        Lz= nz*delta            ! maximum height box s
        Lx= nx*delta            ! maximum width box 
        Ly= ny*delta            ! maximum depth box 

        do k=1,nchains
            Ene=0.0_dp      
            do i=1,nseg
                do j=i+1,nseg

                    s=type_of_monomer(i)
                    t=type_of_monomer(j)
                   
                    zi = pbc(chain(1,i,k),Lz)
                    xi = pbc(chain(2,i,k),Lx)
                    yi = pbc(chain(3,i,k),Ly) 
                    zj = pbc(chain(1,j,k),Lz)
                    xj = pbc(chain(2,j,k),Lx)
                    yj = pbc(chain(3,j,k),Ly) 

                    sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
                    sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

                    Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3

                enddo
            enddo 
            Energy(k) = Ene            
        enddo

    end function VdWpotentialenergy_MC


    ! .. compute periodic boundary condition in z direction
    ! .. maps z onto interval [0,deltaz]
    ! .. see also in chaingeneration.f90

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


    ! pre : chain conformation 
    ! post: VdW of conformation and if conformation is saw or not

    function VdWpotentialenergy_SAW(chain)result(Energy)

        use globals, only : nseg
        use chains, only : type_of_monomer
        use parameters, only :  lsegAA,VdWeps

        real(dp), intent(in)  :: chain(:,:)
        real(dp) :: Energy
       
        logical :: saw

        real(dp) :: Ene,sqrlseg,sqrdist
        integer :: i,j,s,t
        real(dp) :: xi,xj,yi,yj,zi,zj
       
        Ene=0.0_dp      
        saw=.true.

        do i=1,nseg
            do j=i+1,nseg

                s=type_of_monomer(i)
                t=type_of_monomer(j)

                zi = chain(1,i)
                xi = chain(2,i)
                yi = chain(3,i)

                zj = chain(1,j)
                xj = chain(2,j)
                yj = chain(3,j) 

                sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
                sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

                Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3
                
                if(j/=(i+1).and.sqrdist<sqrlseg) then 
                    saw=.false.
                    print*,"overlap occured for i= ",i," and j= ",j 
                endif    
            enddo
        enddo 

        !print*,"total ",Ene," saw ",saw,' ',(lsegAA(t),t=1,nsegtypes),VdWeps 
        Energy = Ene            

    end function VdWpotentialenergy_SAW
   

    ! pre : chain conformation 
    ! post: VdW of conformation 
    ! V_VdW(r) = -epsilon * (sigma/r)^6 : r >= 2*(1/6) sigma
    !          = -epsilon               : r  < 2*(1/6) sigma

    function VdWpotentialenergy(chain)result(Energy)

        use globals, only : nseg
        use chains, only : type_of_monomer
        use parameters, only :  lsegAA,VdWeps

        real(dp), intent(in)  :: chain(:,:)
        real(dp) :: Energy
        
        real(dp) :: Ene,sqrlseg,sqrdist !,maxdist
        integer ::  i,j,s,t
        real(dp) :: xi,xj,yi,yj,zi,zj
        
        Ene=0.0_dp 

        !maxdist=VdWcutoff*delta 
       
        do i=1,nseg

                s=type_of_monomer(i)
                
                zi = chain(1,i)
                xi = chain(2,i)
                yi = chain(3,i)


            do j=i+1,nseg

                t=type_of_monomer(j)

                zj = chain(1,j)
                xj = chain(2,j)
                yj = chain(3,j) 

                sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
                sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

                if(sqrdist<sqrlseg) sqrdist=sqrlseg

                Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3
                
            enddo
        enddo 

        Energy = Ene            

    end function  VdWpotentialenergy

     
    ! pre : chain conformation 
    ! post: VLJ energy of conformation 
    ! VLJ(r) = 4 epsilon * [ (sigma/r)^12-(sigma/r)^6 

    function  LJpotentialenergy(chain) result(Energy)

        use globals, only : nseg
        use chains, only : type_of_monomer
        use parameters, only :  lsegAA,VdWeps

        real(dp), intent(in)  :: chain(:,:)
        real(dp) :: Energy
        
        real(dp) :: Ene,sqrlseg,sqrdist
        integer ::  i,j,s,t
        real(dp) :: xi,xj,yi,yj,zi,zj
        
        Ene=0.0_dp 
       
        do i=1,nseg
            s=type_of_monomer(i)
            
            zi = chain(1,i)
            xi = chain(2,i)
            yi = chain(3,i)
            
            do j=i+1,nseg

                t=type_of_monomer(j)

                zj = chain(1,j)
                xj = chain(2,j)
                yj = chain(3,j) 

                sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
                sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2
        
                Ene=Ene + 4.0_dp*VdWeps(s,t)*((sqrlseg/sqrdist)**6-(sqrlseg/sqrdist)**3)
            
            enddo
        enddo 

        Energy = Ene            

    end function LJpotentialenergy

    ! pre : chain conformation 
    ! post: effective VLJ energy conformation using location COM of nucleosome
    ! VLJ(r) = 4 epsilon * [ (sigma/r)^12-(sigma/r)^6 

    function  LJenergyeffective(chain,nmer)result(Energy)

        use chains, only : segcm
        !use parameters, only :  omegaLJ,epsLJ

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer

        real(dp) :: Energy
        
        real(dp) :: Ene,sqrdist,sqromega
        integer ::  i,j
        real(dp) :: xi,xj,yi,yj,zi,zj
        integer :: isegcm, jsegcm

        real(dp) :: omegaLJ  ! omega of LJ should be in parameters or GBpotential module
        real(dp) :: epsLJ     

        Ene=0.0_dp 
        sqromega=omegaLJ**2
        
        print*,"Warning : LJenergyeffective : omegaLJ and epsLJ uninitialized"

       
        do i=1,nmer
            isegcm=segcm(i)
            
            xi = chain(1,isegcm)
            yi = chain(2,isegcm)
            zi = chain(3,isegcm)

            do j=i+1,nmer
                jsegcm=segcm(j)

                xj = chain(1,jsegcm)
                yj = chain(2,jsegcm)
                zj = chain(3,jsegcm) 

                sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2

                Ene=Ene + 4.0_dp*epsLJ*((sqromega/sqrdist)**6-(sqromega/sqrdist)**3)
            
            enddo
        enddo     

        Energy = Ene            

    end function LJenergyeffective


    ! Inner product of vectors a and b with dimenstion 3
    ! similar to intrisic function dot_product

    function dotproduct(a,b) result(dotprod)
       
        real(dp), dimension(3) :: a, b 
        real(dp) :: dotprod
      
        dotprod = sum(a*b)

    end function dotproduct


    ! Inits constant variable used VGBenergyeffective: makes array segunitvector 
    ! Inits constant parameter of module GB_potential
    ! input : integer segcm
    !         integer nmer
    !         integer segnumAAstart(nmer)
    !         integer segnumAAend(nmer)
    ! output: integer segunitvector(nmer) 
    !         and all variable of Gay-Bern Potential in module GB_potential

    subroutine init_GBenergyeffective(segcm,nmer,segnumAAstart,segnumAAend)
        
        use GB_potential, only : init_GB_const_defaults, init_GB_const
        use chains, only : unitvector_triplets

        integer, intent(in) :: segcm(:)
        integer, intent(in) :: nmer
        integer, intent(in)  :: segnumAAstart(nmer), segnumAAend(nmer) ! segment numbers first/last AAs 

        integer :: i, info

        info = 0

        allocate(segunitvector(nmer)) 

        ! set default for GB potential

        call init_GB_const_defaults()

        ! set default auxilary variable for determine COM and unit vector 

        deltaseg= (1577-1069) ! using only AA number from pca python program

        ! init unitvector from python program pca 
        
        unit_vec_ref= (/0.99857452_dp, 0.02966983_dp, 0.0443692_dp/) 

        ! set default  reference orientation vectors spanning plane through origin

        vec_ref(:,1)=(/1.08398777_dp, -5.05775356_dp,  3.57485842_dp/)
        vec_ref(:,2)=(/-1.95081223_dp,  0.48525144_dp, -2.62704158_dp/)
        vec_ref(:,3)=(/1.43088777_dp, -1.73315356_dp,  2.53525842_dp/)
   
        ! relative to atomid/segnumber of first AA

        atom_id_unit_relative=(/0, 233, 720/)  

        ! overide default 
        call read_GBinputfile(info)

        ! compute dependent auxilary functions
        
        call init_GB_const()

        do i=1,nmer
            segunitvector(i)=segcm(i)+deltaseg
        enddo
        
        allocate(unitvector_triplets(nmer,3))
        allocate(segnumAAstartGBcom(nmer))
        allocate(segnumAAendGBcom(nmer))

        segnumAAstartGBcom = segnumAAstart
        segnumAAendGBcom   = segnumAAend

        print*,"EGB_threshold=",EGB_threshold     
 
    end subroutine init_GBenergyeffective

    
    ! pre : chain conformation and segcm and segunitvector
    ! post: effective VGB energy conformation using location COM of nucleosome
    ! VLJ(u1,u2,r) = Gay-Berne potential 
    ! Assumes that the directionof the vector connecting AAs number segcom and segunitvector 
    ! is an unitvector. 

    function  GBenergyeffective_comb(chain,nmer)result(Energy)

        use chains, only : segcm 
        use GB_potential, only : GBpotential_general

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer
        real(dp) :: Energy
        
        real(dp) :: Ene
        integer ::  i,j
        integer :: isegcm, isegunitvector, jsegcm, jsegunitvector
        real(dp) :: xi, xj, yi, yj, zi, zj
        real(dp) :: xiunit, xjunit, yiunit, yjunit, ziunit,zjunit
        real(dp) :: u1(3), u2(3), norm_u1, norm_u2 ! unitvector and its norm
        real(dp) :: rvec(3)  ! distance between com

        Ene=0.0_dp 
       
        do i=1,nmer

            isegcm=segcm(i) ! atom_id of AA that is (closest to) COM
            xi = chain(1,isegcm)
            yi = chain(2,isegcm)
            zi = chain(3,isegcm)

            isegunitvector=segunitvector(i) ! atom_id of AA with isegcm to make unit vector
            xiunit = chain(1,isegunitvector)
            yiunit = chain(2,isegunitvector)
            ziunit = chain(3,isegunitvector)

            u1=(/xiunit-xi,yiunit-yi,ziunit-zi/) ! unit vector

            norm_u1=sqrt(dotproduct(u1,u1)) 
            u1=u1/norm_u1

            do j=i+1,nmer

                jsegcm=segcm(j)
                xj = chain(1,jsegcm)
                yj = chain(2,jsegcm)
                zj = chain(3,jsegcm) 

                jsegunitvector=segunitvector(j)

                xjunit = chain(1,jsegunitvector)
                yjunit = chain(2,jsegunitvector)
                zjunit = chain(3,jsegunitvector)

                u2=(/xjunit-xj,yjunit-yj,zjunit-zj/)
                norm_u2=sqrt(dotproduct(u2,u2)) 
                u2=u2/norm_u2

                rvec=(/xj-xi,yj-yi,zj-zi/)
               
                Ene=Ene + GBpotential_general(u1,u2,rvec)
            
            enddo
        enddo     

        Energy = Ene          

    end function GBenergyeffective_comb

    ! Calculate sum value Gay-Berne potential for nucleosome conformation 
    ! input :  real(dp :: chainchain(:,:) conformation 
    !          integer :: nmer : nnucl 
    ! output:  logical :: no_overlap
    ! output/return:  real(dp) :: energy 

    function  GBenergyeffective(chain,nmer,no_overlap)result(Energy)
 
        use chains, only : segcm, unitvector_triplets
        use GB_potential, only : GBCOMtype

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer
        logical, intent(inout) :: no_overlap
        real(dp) :: Energy

        real(dp) :: rcom(3,nmer)
        real(dp) :: uvector(3,nmer)
        integer :: info,  un_traj_com

        select case (GBCOMtype) ! GBCOMtype selects algorithm to make com and unitvectors
        case ("simple") 

            call make_com_unit_vector_nucl_simple(chain,nmer,segcm,segunitvector,rcom,uvector)

        case ("simpleCOM") 

            call make_com_unit_vector_nucl_simpleCOM(chain,nmer,segcm,segunitvector,rcom,uvector)

        case("rotation")

            call make_com_unit_vector_nucl_rotation(chain,nmer,unitvector_triplets,rcom,uvector)
        
        case default
        
            print*,"Error in GBenergyeffective"
            print*,"Wrong value GBCOMtype : ", GBCOMtype
        
        end select

        Energy = GBenergyeff(nmer,rcom,uvector,no_overlap)
  
        print*,"rcom   =", rcom
        print*,"uvector=", uvector
        ! 
        info = 1
        un_traj_com= open_chain_com_unitvec_lammps_trj(info)
        call write_chain_com_unitvec_lammps_trj(un_traj_com,rcom,uvector)
        close(un_traj_com)

    end function GBenergyeffective


    !  Calculate sum of value Gay-Berne potential for given set of com and unitvectors
    !  input : real(dp) : rcom(:,:) : com of nucleosome 
    !          real(dp) : uvector(:,:) : orientation vector ( normalized) of each nucleosome 
    !  output:  logical : no_overlap
    !  output/return : real(dp)  energy : effective GB energy = VLJ(u1,u2,r)  

    function  GBenergyeff(nmer,rcom,uvector,no_overlap)result(Energy)

        use GB_potential, only : GBpotential_general

        integer, intent(in)  :: nmer
        real(dp), intent(in) :: rcom(:,:)
        real(dp), intent(in)  :: uvector(:,:) 
        logical,  intent(inout) :: no_overlap

        real(dp) :: Energy
        
        real(dp) :: Ene, EneGB
        integer ::  i,j
        real(dp) :: uveci(3), uvecj(3) ! unitvector
        real(dp) :: rvec(3), rveci(3), rvecj(3) ! distance between com's

        no_overlap=.true.

        Ene=0.0_dp 
       
        do i=1,nmer

            rveci = rcom(:,i)
            uveci = uvector(:,i)

            do j=i+1,nmer

                rvecj =rcom(:,j)
                uvecj =uvector(:,j) 
                rvec = rvecj-rveci
                EneGB = GBpotential_general(uveci,uvecj,rvec)
                Ene = Ene + EneGB
                if( EneGB > EGB_threshold) no_overlap=.false.
            enddo
        enddo     

        Energy = Ene            

    end function GBenergyeff


    ! Determine com and unit vector for a given nucloesome conformation
    ! pre:  subroutine init_GBenergyeffective
    ! input :  real(dp :: chainchain(:,:) conformation 
    !          integer :: nmer : nnucl
    !          integer :: segcm(:) : atom_id closest to com 
    !          integer :: segunitvector(:) atom_id used to make unit vector
    ! output:  real(dp) :: rcom(:,:)    : com of each nucleosome 
    !          real(dp) :: uvector(:,:) : unit vector of each nucleosome

    subroutine make_com_unit_vector_nucl_simple(chain,nmer,segcm,segunitvector,rcom,uvector)

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer
        integer, intent(in) :: segcm(:)
        integer, intent(in) :: segunitvector(:) 
        real(dp), intent(inout) :: rcom(:,:)
        real(dp), intent(inout) :: uvector(:,:)
      

        integer :: i, isegcm, isegunitvec
        real(dp) :: rvec(3), rvecunitvec(3), uvec(3), norm_uvec

        do i=1,nmer

            isegcm=segcm(i) ! atom_id of AA that is (closest to) COM
            rvec = chain(:,isegcm)
        
            isegunitvec=segunitvector(i) ! atom_id of AA with isegcm to make unit vector
            rvecunitvec = chain(:,isegunitvec)

            uvec=rvecunitvec-rvec   ! unit vector
            norm_uvec=sqrt(dot_product(uvec,uvec))  
            uvec=uvec/norm_uvec  ! normalize

            rcom(:,i) = rvec
            uvector(:,i) = uvec

        enddo     

    end subroutine make_com_unit_vector_nucl_simple


    ! Determine com and unit vector for a given nucloesome conformation
    ! pre:  call to   init_GBenergyeffective
    ! input :  real(dp :: chainchain(:,:) conformation 
    !          integer :: nmer : nnucl
    !          integer :: segcm(:) : atom_id closest to com 
    !          integer :: segunitvector(:) atom_id used to make unit vector
    ! output:  real(dp) :: rcom(:,:)    : com of each nucleosome 
    !          real(dp) :: uvector(:,:) : unit vector of each nucleosome

    subroutine make_com_unit_vector_nucl_simpleCOM(chain,nmer,segcm,segunitvector,rcom,uvector)

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer
        integer, intent(in) :: segcm(:)
        integer, intent(in) :: segunitvector(:) 
        real(dp), intent(inout) :: rcom(:,:)
        real(dp), intent(inout) :: uvector(:,:)
      
        integer :: i, isegcm, isegunitvec
        real(dp) :: rvec(3), rvecunitvec(3), uvec(3), norm_uvec

        rcom = com_chains(chain,nmer)
       
        do i=1,nmer

            isegcm=segcm(i) ! atom_id of AA that is (closest to) COM
            rvec = chain(:,isegcm)
        
            isegunitvec=segunitvector(i) ! atom_id of AA with isegcm to make unit vector
            rvecunitvec = chain(:,isegunitvec)

            uvec=rvecunitvec-rvec   ! unit vector
            norm_uvec=sqrt(dot_product(uvec,uvec))  
            uvec=uvec/norm_uvec     ! normalize

            uvector(:,i) = uvec
        enddo     

    end subroutine make_com_unit_vector_nucl_simpleCOM



    ! Determine com  for a given nucloesome conformation
    ! input :  real(dp :: chainchain(:,:) conformation 
    !          integer :: nmer : nnucl      
    ! return output:  real(dp) :: rcom(:,:)    : com of each nucleosome 

    function com_chains(chain,nmer)result(rcom)

        real(dp), intent(in) :: chain(:,:)
        integer, intent(in) :: nmer
       
        real(dp) :: rcom(3,nmer)
       
        integer :: n, i, numelements

        do n=1,nmer
            
            rcom(:,n)=0.0_dp
            numelements=segnumAAendGBcom(n)-segnumAAstartGBcom(n)+1
           
            do i=segnumAAstartGBcom(n),segnumAAendGBcom(n)
                rcom(:,n) = rcom(:,n)+chain(:,i)
            enddo    

            rcom(:,n) = rcom(:,n)/(numelements*1.0_dp)

        enddo     

    end function com_chains



    ! Computes value of parameteric function s and t of plane that is spanned by 
    ! v=vec[2]-vec[1] and w=vec[3]-vec[1] at origin: x= s*v + t*w +vec[1]=0  
    ! input : real(dp) :: vec(3,3) three xyz-coordinate triplets/vectors 
    ! output: real(dp) ::  sval, tval:  value of parameteric function plane as origin
    !        real(dp) :: vec_origin(3)

    subroutine parameter_com_ref(vec,sval,tval,vec_origin)

        real(dp), intent(in)   :: vec(3,3) 
        real(dp), intent(inout) :: sval, tval
        real(dp), intent(inout) :: vec_origin(3) 

        real(dp) :: v1(3),v2(3), u1(3), u2(3),  sqr_norm_u1,  sqr_norm_u2

        v1=vec(:,2)-vec(:,1)
        v2=vec(:,3)-vec(:,1)
        
        ! orthogonal vectors
        u1=v1
        u2=v2-(dot_product(v2,u1)/dot_product(u1,u1))*u1
        

        sqr_norm_u1=dot_product(u1,u1)
        sqr_norm_u2=dot_product(u2,u2)
    
        ! value of 

        sval = -dot_product(u1,vec(:,1))/sqr_norm_u1
        tval = -dot_product(u2,vec(:,1))/sqr_norm_u2
       
        ! print*,"parameter_com_ref"
        ! print*,"u1=",u1 
        ! print*,"u2=",u2 
        ! print*,"(u1,vec(:,1)=",dot_product(u1,vec(:,1)),"(u1,u1)=",sqr_norm_u1
        ! print*,"s=",sval," t=",tval

        ! compute (relative) origin/check
   
        vec_origin= sval*u1 + tval*u2 +vec(:,1)

    end subroutine parameter_com_ref



    ! Computes value of parameteric function s and t of plane that is spanned by 
    ! v=vec[2]-vec[1] and w=vec[3]-vec[1] at origin: x= s*v + t*w +vec[1]=0  
    ! input : real(dp) :: vec(3,3) three xyz-coordinate triplets/vectors 
    ! output: integer ::  s_ref, t_ref:  value of parameteric function plane as origin
    !        real(dp) :: vec_origin(3)

    subroutine parameter_com(vec,s_ref,t_ref,vec_origin)

        real(dp), intent(in)   :: vec(3,3) 
        real(dp), intent(in) ::  s_ref, t_ref
        real(dp), intent(inout) :: vec_origin(3) 

        real(dp) :: v1(3),v2(3), u1(3), u2(3)

        v1=vec(:,2)-vec(:,1)
        v2=vec(:,3)-vec(:,1)
        
        ! orthogonal vectors
        u1=v1
        u2=v2-(dot_product(v2,u1)/dot_product(u1,u1))*u1
    
        ! compute (relative) origin
   
        vec_origin= s_ref*u1 + t_ref*u2 + vec(:,1)
        
        ! print*,"parameter_com"
        ! print*,"u1=",u1 
        ! print*,"u2=",u2 
        ! print*,"(u1,vec(:,1)=",dot_product(u1,vec(:,1)),"(u1,u1)=",dot_product(u1,u1)
        ! print*,"s=",s_ref," t=",t_ref

    end subroutine parameter_com     
    
    ! Calculation the normal vector to plane spanned by three points
    ! input: real(dp)  vec(3,3)       : three xyz-coordinate triplets/vectors
    ! output: real(dp) norm_vector(3) : normal vector

    function normal_vector(vec)result(norm_vector)

        use quaternions, only : crossproduct

        real(dp), intent(in) :: vec(3,3)  
        real(dp) :: norm_vector(3)

        real(dp) :: v1(3), v2(3)

        v1=vec(:,2)-vec(:,1)
        v2=vec(:,3)-vec(:,1)
        norm_vector=crossproduct(v1,v2)
        
    end function normal_vector
   

    ! Determine com and unit vector for a given nucleosome conformation
    ! 1. com : using set of 3 points==(unitvector_triplets) as auxilary points to span a plane that include the
    !    origin plane or com
    ! 2. orientation vector: using 3 points and 3 point of a reference to compute rotation and apply rotation 
    !    one referece orientation vector 
    ! pre:  subroutine  init_GBenergyeffective
    ! input :  real(dp :: chain(:,:) conformation 
    !          integer :: nmer : nnucl
    !          integer :: unitvector_triplets(:) atom_id used to make unit vector
    ! output:  real(dp) :: rcom(:,:)    : com of each nucleosome 
    !          real(dp) :: uvector(:,:) : unit vector of each nucleosome

    subroutine make_com_unit_vector_nucl_rotation(chain,nmer,unitvector_triplets,rcom,uvector)

        use quaternions, only : rot_axis_angle, rot_axis_angle_to_quat, rotation_matrix_from_quat, vec_norm

        real(dp), intent(in) :: chain(:,:)               ! dimension : (3,s) 
        integer, intent(in) :: nmer
        integer, intent(inout) :: unitvector_triplets(:,:)  ! dimension : (nnucl,3)=(nmer,3)
        real(dp), intent(inout) :: rcom(:,:)             ! dimension : (nnucl,3)
        real(dp), intent(inout) :: uvector(:,:)          ! dimension : (nnucl,3)
      
        real(dp) :: s_ref, t_ref
        real(dp) :: unit_vec(3), vec(3,3)
        real(dp) :: vec_origin_ref(3), vec_origin(3)
        real(dp) :: triangle_vec_ref(3,3), triangle_vec(3,3,nmer)
        real(dp) :: norm_vec_ref(3), norm_vec(3,nmer), norm_tri(3,nmer)
        real(dp) :: a(3), a_rot(3), b(3), u(3)
        real(dp) :: rcom_ref(3), unorm, angle
        real(dp) :: qu(4)
        real(dp) :: Rmat1(3,3), Rmat2(3,3), Rmat_comb(3,3)
        integer  :: k, n
         
        real(dp), parameter  :: epsAngle = 0.0000000000001_dp ! init

        ! determine com using unitvectortriplets
        ! determine s_ref and t_ref value of parameter function of plane of reference vectors. 
        ! use vec_ref and unit_vec_ref and atom_id_unit_relative in init_GBenergyeffective

        ! atom id AA of all nucleosome
        do n=1,nmer 
            do k=1,3
                unitvector_triplets(n,k)= atom_id_unit_relative(k)+segnumAAstartGBcom(n)
            enddo
        enddo    

        ! end init 

        call parameter_com_ref(vec_ref,s_ref,t_ref,vec_origin_ref)
        rcom_ref=vec_origin_ref
        norm_vec_ref=normal_vector(vec_ref)

        do k=1,3
            triangle_vec_ref(:,k) = vec_ref(:,k)-rcom_ref
        enddo

        ! print*,"s_ref=",s_ref," t_ref=",t_ref
        ! print*,"rcom_ref        =",rcom_ref
        ! print*,"vec_origin_ref  =",vec_origin_ref

        ! do k=1,3
        !     print*,"vec_ref         =",vec_ref(:,k)
        !     print*,"triangle_vec_ref=",triangle_vec_ref(:,k)
        ! enddo    
        
        ! compute com 

        do n=1,nmer 
            do k=1,3 ! get vector spanning plane
                ! print*,'n=',n,'k=',k,'unitvector_triplets(n,k)=',unitvector_triplets(n,k)
                vec(:,k)=chain(:,unitvector_triplets(n,k))
            enddo

            call parameter_com(vec,s_ref,t_ref,vec_origin)
            rcom(:,n) = vec_origin
            norm_vec(:,n)=normal_vector(vec)
            
            do k=1,3 ! get vectors from 'origin' plane
               triangle_vec(:,k,n) =vec(:,k)-rcom(:,n)
            enddo   
            
            norm_tri(:,n)=normal_vector(triangle_vec(:,:,n))

            ! print*,"norm_vec(:,n)=",norm_vec(:,n)
            ! print*,"norm_tri(:,n)=",norm_tri(:,n)
            ! print*,""
            ! print*,"rcom(:,n)=",rcom(:,n)


        enddo

        ! compute unitvectors

        do n=1,nmer                      ! loop over nucleosomes 

            ! rotation of norm vector  

            call rot_axis_angle(norm_vec_ref, norm_vec(:,n), u, angle)
            unorm = vec_norm(u)
            u(1:3) = u(1:3)/unorm         ! rotation axis 

            if(abs(angle)<epsAngle) u=a   ! angle ==0 no rotation u=>a  

            qu = rot_axis_angle_to_quat(u, angle)   ! rotation quaternion
            call rotation_matrix_from_quat(qu,Rmat1)  

            ! apply second in plane rotation to get coordinate to coincide

            a =triangle_vec_ref(:,2)  

            a_rot = matmul(Rmat1,a) ! apply rotation
            b = triangle_vec(:,2,n)
            call rot_axis_angle(a_rot, b, u, angle)
            unorm = vec_norm(u)
            u(1:3) = u(1:3)/unorm         ! rotation axis 

            if(abs(angle)<epsAngle) u=a   ! angle ==0 no rotation u=>a  

            qu = rot_axis_angle_to_quat(u, angle)   ! rotation quaternion
            call rotation_matrix_from_quat(qu,Rmat2)  

            Rmat_comb=matmul(Rmat2,Rmat1) ! Rmat_comb= Rmat2 x Rmat1
     
            ! apply rotation to unit_vec_ref

            unit_vec=matmul(Rmat_comb,unit_vec_ref)

            ! assign vector 
            uvector(:,n)=unit_vec
            print*,"unit_vec_ref=",unit_vec_ref    

        enddo   
  
    end subroutine make_com_unit_vector_nucl_rotation

    ! Determine com a for a given nucleosome conformation
    ! 1. com : using set of 3 points==(unitvector_triplets) as auxilary points to span a plane that include the
    !    origin plane or com
    !  
    ! pre:  subroutine  init_GBenergyeffective
    ! input :  real(dp :: chain(:,:) conformation 
    !          integer :: nmer : nnucl
    !          integer :: unitvector_triplets(:) atom_id used to make unit vector
    ! output:  real(dp) :: rcom(:,:)    : com of each nucleosome 

    subroutine make_com_nucl_rotation(chain,nmer,unitvector_triplets,rcom)

        use quaternions, only : rot_axis_angle, rot_axis_angle_to_quat, rotation_matrix_from_quat, vec_norm

        real(dp), intent(in) :: chain(:,:)               ! dimension : (3,s) 
        integer, intent(in) :: nmer
        integer, intent(inout) :: unitvector_triplets(:,:)  ! dimension : (nnucl,3)=(nmer,3)
        real(dp), intent(inout) :: rcom(:,:)             ! dimension : (nnucl,3)
      
        real(dp) :: s_ref, t_ref
        real(dp) :: vec(3,3) !, unit_vec(3)
        real(dp) :: vec_origin_ref(3), vec_origin(3)
        real(dp) :: triangle_vec_ref(3,3) !, triangle_vec(3,3,nmer)
        real(dp) :: norm_vec_ref(3) !, norm_vec(3,nmer)
        real(dp) :: rcom_ref(3)
        integer  :: k, n
         
        ! determine com using unitvectortriplets
        ! determine s_ref and t_ref value of parameter function of plane of reference vectors. 
        ! use vec_ref and unit_vec_ref and atom_id_unit_relative in init_GBenergyeffective

        ! atom id AA of all nucleosome
        do n=1,nmer 
            do k=1,3
                unitvector_triplets(n,k)= atom_id_unit_relative(k)+segnumAAstartGBcom(n)
            enddo
        enddo    

        ! end init 

        call parameter_com_ref(vec_ref,s_ref,t_ref,vec_origin_ref)
        rcom_ref=vec_origin_ref
        norm_vec_ref=normal_vector(vec_ref)

        do k=1,3
            triangle_vec_ref(:,k) = vec_ref(:,k)-rcom_ref
        enddo

        ! print*,"s_ref=",s_ref," t_ref=",t_ref
        ! print*,"rcom_ref        =",rcom_ref
        ! print*,"vec_origin_ref  =",vec_origin_ref

        ! do k=1,3
        !     print*,"vec_ref         =",vec_ref(:,k)
        !     print*,"triangle_vec_ref=",triangle_vec_ref(:,k)
        ! enddo    
        
        ! compute com 

        do n=1,nmer 
            do k=1,3 ! get vector spanning plane
                vec(:,k)=chain(:,unitvector_triplets(n,k))
            enddo

            call parameter_com(vec,s_ref,t_ref,vec_origin)
            rcom(:,n) = vec_origin
           
            ! print*,""
            ! print*,"rcom(:,n)=",rcom(:,n)
        enddo

  
    end subroutine make_com_nucl_rotation


    ! Writes the coordinates of com of nucleosomes orientation vector associated with com to a lammps trajectory file.
    ! Unit length in Angstrom
    ! pre : file opened with open_com_unitvec_lammps_trj
    ! input: integer un_trj : unit number of file to write to.
    !        real(dp) rcom(:) : coordinates of com of nucleosomes
    !        real(dp) unitvector : orientation vector associated with com 
    !    
    ! return: a lammpstrj file

    subroutine write_chain_com_unitvec_lammps_trj(un_trj,rveccom,unitvector)

        ! use globals, only : nnucl
        use volume, only : delta, nx, ny, nz 
        use quaternions, only :  rot_axis_angle, rot_axis_angle_to_quat
        use GB_potential, only : sigmaS,sigmaE  
        use quaternions, only : quat_norm, vec_norm

        integer, intent(in) :: un_trj
        real(dp), dimension(:,:), intent(in) :: rveccom
        real(dp), dimension(:,:), intent(in) :: unitvector
        
        integer ::  n, numnucl
        real(dp) :: xbox0, xbox1(3)
        integer :: idatom, item
        real(dp) :: shapex,shapey,shapez
        real(dp) :: rcom(3), uvec(3), uref(3), qu(4), rotaxis(3), theta

        numnucl = size(rveccom, 2) ! asummes dimension rveccom represent correct 

        xbox0 = 0.0_dp 
        xbox1(1) = nx*delta*10.0_dp ! conversion form nm -> Angstrom
        xbox1(2) = ny*delta*10.0_dp
        xbox1(3) = nz*delta*10.0_dp

        item = 0
        idatom = 1
        conf_write_com = conf_write_com+1
        
        ! axis length of Gay-Berne potential
        shapex = sigmaS*10.0_dp ! conversion form nm -> Angstrom
        shapey = sigmaS*10.0_dp 
        shapez = sigmaE*10.0_dp 

        ! reference unit vector
        uref=(/0.0_dp,0.0_dp,1.0_dp/) 

        ! write preamble 
        write(un_trj,'(14A)')'ITEM: TIMESTEP' 
        write(un_trj,*)conf_write_com 
        write(un_trj,'(21A)')'ITEM: NUMBER OF ATOMS' 
        write(un_trj,*)numnucl
        write(un_trj,'(25A)')'ITEM: BOX BOUNDS ff ff ff'
        write(un_trj,*)xbox0,xbox1(1)
        write(un_trj,*)xbox0,xbox1(2)
        write(un_trj,*)xbox0,xbox1(3)
        write(un_trj,'(70A)')'ITEM: ATOMS id type x y z quatw quati quatj quatk shapex shapey shapez'
        
        do n=1, numnucl 
            
            rcom = rveccom(:,n) * 10.0_dp 
            uvec = unitvector(:,n) ! unit less

            call rot_axis_angle(uref, uvec, rotaxis, theta) 
            rotaxis=rotaxis/vec_norm(rotaxis)
            qu = rot_axis_angle_to_quat(rotaxis, theta) 

            ! print*,"norm uvec=",vec_norm(uvec)," norm  qu=",quat_norm(qu),"rotaxis ",vec_norm(rotaxis)
            ! print*,"qu=",qu," norm  qu=",quat_norm(qu)
            item=item+1  
            write(un_trj,*)item,idatom,rcom(1),rcom(2),rcom(3),qu(1),qu(2),qu(3),qu(4),shapex,shapey,shapez

        enddo  

    end subroutine write_chain_com_unitvec_lammps_trj


    ! Creates and opens file for lammps trajectory film
    ! Writing done by write_chain_com_unitvec_lammps_trj
    ! input:        integer : info : value 0 no errors
    ! output/return integer : un_trj : unit nubmer of file

    function open_chain_com_unitvec_lammps_trj(info,text)result(un_trj)

    
        use myutils, only : newunit, lenText
        use myio, only : myio_err_chainsfile
        use GB_potential, only : GBtype, GBCOMtype

        integer, intent(inout) :: info
        character(len=*), intent(in),optional :: text ! extra filename descriptor 

        integer :: un_trj

        ! local
        character(len=lenText) :: fname
        integer :: ios
        logical :: exist
        
        info = 0    

        if(present(text)) then
            fname='trajcomunit_'//trim(adjustl(text))//'_'
        else
            fname='trajcomunit_'
        endif

        fname=trim(fname)//trim(adjustl(GBtype))//trim(adjustl(GBCOMtype))
        fname=trim(adjustl(fname))//'.lammpstrj'
       
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

    end function open_chain_com_unitvec_lammps_trj


    subroutine read_GBinputfile(info)

        use myutils, only : newunit
        use GB_potential, only : epsilonE, epsilonS 
        use myio, only : myio_err_GBinputfile, myio_err_GBinputlabel

        integer, intent(out), optional :: info

        ! .. local arguments

        character(len=10) :: fname
        integer :: ios,un_input  ! un = unit number
        character(len=100) :: buffer, label
        integer :: pos
        integer :: line
        logical :: file_exist

        if (present(info)) info = 0

        !     .. reading in of variables from file
        write(fname,'(A10)')'GBinput.in'

        inquire(file=fname,exist=file_exist)
        print*,"file_exist=",file_exist

        if(file_exist) then
            open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
        else
            print*,'GBinput file does not exit: use default value for GB potential'
            info = 0
            return
        endif

        if(ios >0 ) then
            print*, 'Error opening GBinput.in file : iostat =', ios
            if (present(info)) info = myio_err_GBinputfile
            return
        endif

        ios = 0
        line = 0

        ! ios<0 : if an end of record condition is encountered or if an end of file condition was detected.
        ! ios>0 : if an error occured
        ! ios=0 : otherwise.

        do while (ios == 0)

            read(un_input, '(A)', iostat=ios) buffer

            if (ios == 0) then

                line = line + 1

                !  Split label and data based on first occurence of a whitespace
                pos = scan(buffer, '     ')
                label = buffer(1:pos)
                buffer = buffer(pos+1:)

                select case (label) !list-directed The character variable is treated as an 'internal file'
                case ('epsilonE')
                    read(buffer, *,iostat=ios) epsilonE
                case ('epsilonS')
                    read(buffer, *,iostat=ios) epsilonS
                case ('deltaseg')
                    read(buffer, *,iostat=ios) deltaseg  
                case ('EGB_threshold')
                    read(buffer, *,iostat=ios) EGB_threshold
                case default
                    if(pos>1) then
                        print *, 'Invalid label at line', line  ! empty lines are skipped
                        if (present(info)) info = myio_err_GBinputlabel
                        return
                    endif
                end select
            endif
        enddo

        if(ios >0 ) then
            print*, 'Error parsing file : iostat =', ios
            print*, 'Read error at line =', line
            if (present(info)) info = myio_err_GBinputfile
            return
        endif

        close(un_input)

        ! .. check error flag

    end subroutine read_GBinputfile

end module VdW_potential    