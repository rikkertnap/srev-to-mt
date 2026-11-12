 module parameters

    use volume
    use molecules
    use loopvar

    implicit none

    !  .. list of parameters

    type(moleclist) :: xbulk
    type(moleclist) :: expmu
    type(moleclist) :: ion_excess  
    type(bornmoleclist) :: bornrad,bornbulk 
    type(moleclist) :: beta_ion_excess

    ! .. index for different chemical states of phosphate used in vPP and qPP 
    ! .. in for systype nucl_ionbin_Mg and systype nucl_ionbin_Fe
    ! .. warning numbering different  from avfdisA 

    integer, parameter :: Phos=1 
    integer, parameter :: PhosH=2
    integer, parameter :: PhosK=3
    integer, parameter :: PhosNa=4
    integer, parameter :: PhosMg=5
    integer, parameter :: PhosFe2=6
    integer, parameter :: PhosFe3=7
    integer, parameter :: Phos2Mg=8   ! warning renumbered this from 6 to 7 to 8  !!!
    integer, parameter :: Phos2Fe2=9
    integer, parameter :: Phos2Fe3=10
    integer, parameter :: Phos3Fe3=11


    !  .. volume 
  
    real(dp) :: vsol                 ! volume of solvent  in nm^3       
    real(dp) :: vpolA(5),deltavA(4)  ! volume of one polymer segment, vpol  in units of vsol
    real(dp) :: vpolB(5),deltavB(4)
    real(dp) :: vpolAA(13),deltavAA(13)        
    real(dp), dimension(:), allocatable         :: vpol   ! volume of polymer segment of given type, vpol in units of vsol
    real(dp), dimension(:,:), allocatable       :: vnucl  ! volume of segment type t element j  
    character(len=3), dimension(:), allocatable :: vnucl_type_char
    real(dp), dimension(:),   allocatable       :: vnucl_type
    logical , dimension(:),   allocatable       :: vnucl_type_isChargeable
    real(dp), dimension(10) :: vPP
    real(dp), dimension(11) :: vPPP ! volume of different chemical states of phosphate 
    
    real(dp) :: vNa                ! volume Na+ ion in units of vsol
    real(dp) :: vK                 ! volume K+  ion in units of vsol
    real(dp) :: vCl                ! volume Cl_ion in units of vsol   
    real(dp) :: vCa                ! volume positive Ca2+ ion in units of vsol
    real(dp) :: vMg                ! volume positive Mg2+ ion in units of vsol
    real(dp) :: vFe2               ! volume positive Fe2+ ion in units of vsol  
    real(dp) :: vFe3               ! volume positive Fe3+ ion in units of vsol
    real(dp) :: vNaCl              ! volume NaCl ionpair in units of vsol 
    real(dp) :: vKCl               ! volume KCl ionpair in units of vsol 
    real(dp) :: vO2                ! volume dissolved O2 in units of vsol
  
    !  .. radii
  
    real(dp) :: RNa
    real(dp) :: RK
    real(dp) :: RCl
    real(dp) :: RCa
    real(dp) :: RMg
    real(dp) :: RFe2 
    real(dp) :: RFe3
    real(dp) :: RO2

    ! .. charges 

    integer, dimension(:,:), allocatable :: zpol ! valence charge polymer
    integer :: zpolAA(8)
    integer :: zpolA(5)          ! valence charge polymer
    integer :: zpolB(5)          ! valence charge polymer
    integer, dimension(10) :: qPP  ! charge of different phosphate chemical state in pairs
    real, dimension(11)    :: qPPP ! charge of different phophatec chemical state in triplet

    integer :: zNa               ! valence charge Na+ ion 
    integer :: zK                ! valence charge K+ ion 
    integer :: zCa               ! valence charge Ca++ ion 
    integer :: zMg               ! valence charge Mg++ ion 
    integer :: zCl               ! valence charge Cl- ion 
    integer :: zFe2              ! valence charge Fe++ ion    
    integer :: zFe3              ! valence charge Fe+++ ion 
  
    !  .. VdW variables

    real(dp), dimension(:,:), allocatable :: VdWeps, VdWepsin    ! strenght VdW interaction in units of kT
    real(dp) :: VdWepsAA, VdWepsBB,VdWepsAB            ! strenght VdW interaction in units of kT
    logical :: isVdW              ! if true VdW energy 
    logical :: isVdWintEne        ! if true VdWpotentialenergy is used to compute internal VdW energy chain
    logical :: isChainEnergyFile
    real(dp) :: VdWcutoff         ! cutoff VdW interaction in units of lseg 	
    real(dp) :: VdWcutoffdelta    ! cutoff VdW interaction in units of delta
    real(dp), parameter :: Vdwepsilon=1.0e-5_dp ! thresholds below which VdWeps is assumed to be zero
    logical, dimension(:), allocatable :: isrhoselfconsistent
    type(looplist), target :: VdWscale     ! scale factor in VdW interaction
    type(looplist), target :: dielectscale ! scale factor in dielectric constant: used for loop of eps to descrease elect interactions

     !  .. input filenames 
    integer, parameter :: lenfname=40
    character(len=lenfname) :: chainfname,vpolfname,pKafname,pKaionfname,typesfname,lsegfname,segcmfname
    character(len=lenfname) :: mtpdbfname,vnuclfname,orientfname

    ! .. other physical quanties
  
    real(dp) :: Tref             ! temperature in K
    real(dp) :: dielectW         ! dielectric constant of water
    real(dp) :: dielectP         ! dielectric constant of hydrocarbons/PA
    character(len=15) :: dielect_env ! selects dielectric fun 
    real(dp) :: lb,lb0           ! Bjerrum lengtin water and vacuum   
    real(dp) :: constqW          ! constant in Poisson eq dielectric constant of water  
    real(dp) :: constqWin        ! constant in Poisson eq dielectric constant of water  stored for loop dielect   
    real(dp) :: constq0          ! constant in Poisson eq dielectric constant of vacuum 
    real(dp) :: constqE          ! electrostatic pre-factor in pdf 
    real(dp) :: beta_times_e     ! conversion factor  e/ kBT 
  
    !  .. solver variables

    integer  :: itmax            ! maximum number of iterations
    real(dp) :: tol_conv         ! error imposed accuaracy
    real(dp) :: fnorm            ! L2 norm of residual vector function fcn  
    integer  :: infile           ! infile=1 read input files infile!=1 no input files 
    integer  :: iter             ! counts number of iterations
    integer(8) :: maxniter       ! maximum of nonlinear iterations  
    integer  :: maxfkfunevals    ! maximum number of fkun function evaluations
    character(len=8) :: method   ! method="kinsol"  
    logical ::  precondition     ! controls use of precondtioner 

    !  .. output control

    logical :: write_localcharge     ! if .true. rhoq and rhoqpol also outputted
    logical :: write_iondensities    ! if .true. ion density also outputted
    logical :: write_rotations       ! if .true. extra information by chain rotation test_rotate_nucl_chain written 
    logical :: write_frac            ! if .true. densityfrac, densityfracP and densityfacion also outputted
    logical :: write_Palpha          ! if .true. Palpha is outputted
    logical :: write_sys_only        ! if .true. only systemfile and energy file outputted
    logical :: write_oxygen          ! if .true. xO2 is outputted
    
    ! .. chain variables 
    real(dp) :: lseg              ! segment length of A polymer in nm
    real(dp), dimension(:), allocatable :: lsegAA 
    real(dp) :: lsegA             ! segment length of A polymer in nm
    real(dp) :: lsegB             ! segment length of B polymer in nm
    real(dp) :: lsegPAA   
    real(dp) :: lsegPEG  
    real(dp) :: lsegPAMPS
    real(dp) :: lsegPS  

    character(len=15) :: chainmethod  ! method of generating chains ="MC" or "FILE" 
    character(len=8)  :: chaintype    ! type of chain: diblock,alt
    integer :: readinchains           ! nunmber of used/readin chains
    integer :: chainperiod            ! peridociy of repeat of A or B block 
    integer :: maxnchainsrotations    ! number of rotations read in from input.in, assigned to maxnchain in chaingenerator default 12  
    integer :: maxnchainsrotationsxy  ! number of rotations read in from input.in, assigned to maxnctheta in chaingenerator default 1  
    integer :: tA                     ! segment number type of monomer type A  
    logical :: write_mc_chains        ! if true MC chain write of file
    logical :: write_struct           ! if true structural quantities writting to file
    logical :: isEnergyShift          ! if true energychain is shifted by energychain_min see chaingenerator
    logical :: pbc_chains             ! if true apply pbc to chain conformation
 
    ! ..average structural properties of layer

    real(dp) :: avRsqr             ! average Radius of Gyration Nucleosome chain
    real(dp) :: qpolA              ! charge poly A of layer 
    real(dp) :: qpolB              ! charge poly B of layer 
    real(dp) :: qpol_tot           ! charge poly A+B of layer 
  
    real(dp), dimension(:), allocatable :: qpol                ! charge poly of layer 
    real(dp), dimension(:), allocatable :: avfdis              ! average degree of dissociation of monomer of type t
    real(dp), dimension(:,:), allocatable :: avgdisA,avgdisB   ! average fraction of Acidic and Basic AA in state A,AH,ANa etc 
    
    real(dp) :: avfdisP2Mg, avfdisP2Fe2, avfdisP2Fe3, avfdisPP(7,7)   ! average fraction monomer of phopsphate pairs in chemical state PP,PPH, , etc  
    real(dp) :: avfdisPPP(11,11,11)
    real(dp) :: avfdisA_pairs(13)
    real(dp) :: avfdisA_triplets(13)
    real(dp) :: avfdisA(13)        
    ! average fraction of monomer ta=phosphate in state 
    ! A, AH, ANa, ACa, A2Ca, AMg, A2Mg, AK, AFe2, A2Fe2, AFe3, A2Fe3 
    !  index of avfdispairs <=> index avfdisA 
    !  1 == P^- 2 == PH  , 3  == PNa,   4  == PCa,  5  == P2Ca,  6  == PMg, 7 == P2Mg 
    !  8 == PK, 9 == PFe2, 10 == P2Fe2, 11 == PFe3, 12 == P2Fe3, 13 == P3Fe3  
    real(dp) :: avfdisB(5)         ! average fraction of monomer state
    real(dp) :: sum_ion_excess     ! sum of ion_excess of all ions weighted with valence of ion
    real(dp) :: max_psi(6)         ! maximum electrostatic potential on each face of lattice 

    !  .. weak polyelectrolyte variables 
    !  .. equibrium constant
    real(dp), dimension(:), allocatable :: K0a               ! intrinsic equilibruim constant
    real(dp), dimension(:), allocatable :: Ka                ! experimemtal equilibruim constant 
    real(dp), dimension(:), allocatable :: pKa               ! experimental equilibruim constant pKa= -log[Ka]

    real(dp), dimension(:,:), allocatable :: K0aion          ! intrinsic equilibruim constant
    real(dp), dimension(:,:), allocatable :: Kaion           ! experimemtal equilibruim constant 
    real(dp), dimension(:,:), allocatable :: pKaion          ! experimental equilibruim constant pKaion= -log[Kaion]

    real(dp) :: KaA(4),K0aA(4),pKaA(4)     ! constant for  acid 
    real(dp) :: KaB(4),K0aB(4),pKaB(4)   
    real(dp) :: KaAA(12),K0aAA(12), pKaAA(12) 
    type (looplist), target :: pKd         ! binding constants 
    type (looplist), target :: deltaGd 
      
     ! water equilibruim constant pKw= -log[Kw] ,Kw=[H+][OH-]   
    real(dp) :: pKw                 
  
    real(dp) :: K0ionNa             ! intrinsic equilibruim constant
    real(dp) :: KionNa              ! experimemtal equilibruim constant 
    real(dp) :: pKionNa             ! experimental equilibruim constant pKion= -log[Kion]	 
  
    real(dp) :: K0ionK              ! intrinsic equilibruim constant
    real(dp) :: KionK               ! experimemtal equilibruim constant 
    real(dp) :: pKionK              ! experimental equilibruim constant pKion= -log[Kion]	 

    !     .. bulk concentrations 
   
    real(dp) :: cHplus             ! concentration of H+ in bulk in mol/liter
    real(dp) :: cOHmin             ! concentration of OH- in bulk in mol/liter
    real(dp),target :: cNaCl       ! concentration of salt in bulk in mol/liter
    real(dp),target :: cKCl        ! concentration of salt in bulk in mol/liter
    real(dp) :: cCaCl2             ! concentration of CaCl2 in bulk in mol/liter
    real(dp),target :: cMgCl2      ! concentration of MgCl2 in bulk in mol/liter
    type (looplist), target :: pH
    real(dp) :: pHbulk             ! pH of bulk pH = -log([H+])
    real(dp) :: pOHbulk            ! p0H of bulk p0H = -log([0H-])
    real(dp), target :: cFeCl2     ! concentration of FeCl2 in bulk in mol/liter 
                                   ! FEeCl2 does not truely exit, used here to get number of counter ion correct
    real(dp), target :: cFeCl3     ! concentration of FeCl2 in bulk in mol/liter 
                                   ! FEeCl3 does not truely exit, used here to get number of counter ion correct
    real(dp) :: IS                 ! ion strength solution in mol/liter  

    ! .. bulk variable for oxygen     
    real(dp) :: pO2                ! partial pressure oxygen in Pa
    real(dp) :: Hcp_s              ! Henry constant of oxygen Hcp_s = ca /p unit : mol/(m^3 Pa)
    real(dp) :: Hxp_s              ! Henry constant of oxygen Hxp_s = x /p  unit : 1/Pa
    real(dp) :: rhoidO2            ! ideal number desnity of oxygen 
    logical :: isO2present         ! logical switch to determine if O2 is consider or not            
   
    ! surface charge 
    type (looplist), target :: psiS
    real(dp) :: psiSL    ! surface potenital in reduced units for bc=cp at z = 0 
    real(dp) :: psiSR    ! surface potenital in reduced units for bc=cp at z = nz delta

    ! unit conversion : converts unit of input conformation to nm unit!
    real(dp) :: unit_conv

    !  return error of subroutine read_pKds 
    integer, parameter ::  err_pKdfile_noexist = 1
    integer, parameter ::  err_pKdfile         = 2 
    integer, parameter ::  err_pKderror        = 3

    !  return error  
    integer, parameter ::  err_file_noexist = 1
    integer, parameter ::  err_file         = 2 
    integer, parameter ::  err_error        = 3

    ! eps value for sum xbulk and sum charge in bulk 
    
    real(dp), parameter  :: eps_val =1.0e-14_dp 
    
    private :: err_pKdfile_noexist,err_pKdfile,err_pKderror
    private :: err_file_noexist,err_file,err_error    
    private :: eps_val

contains

    ! determine total number of non linear equations

    subroutine set_size_neq()

        use globals, only: systype,nsegtypes, nsize,LEFT,RIGHT, neq, neqint
        use volume, only : nx, ny, nz
        use myutils, only  : error_handler

        integer :: numeq, t

        nsize= nx*ny*nz

        select case (systype)
            case ("brush_mul") 
                neq = (2+nsegtypes) * nsize 
            case ("brush_mulnoVdW") 
                neq = 2 * nsize     
            case ("nucl_neutral_sv")
                neq =  nsize 
            case ("nucl_ionbin_sv","nucl_ionbin_Mg","nucl_ionbin_MgA","nucl_ionbin_Fe","nonuclcp")
                neq = 2 * nsize 
            case ("brushdna","nucl_ionbin")
                numeq=0 
                do t=1,nsegtypes
                    if(isrhoselfconsistent(t)) numeq=numeq+1
                enddo    
                neq = (2+numeq) * nsize 
            case ("brushborn")
                numeq=0 
                do t=1,nsegtypes
                    if(isrhoselfconsistent(t)) numeq=numeq+1
                enddo    
                neq = nsize*(6+numeq)    
            case ("elect") 
                neq = 4 * nsize 
            case ("neutral") 
                numeq=0
                do t=1,nsegtypes
                    if(isrhoselfconsistent(t)) numeq=numeq+1
                enddo   
                neq = (1+numeq) * nsize
            case ("neutralnoVdW") 
                neq = nsize  
            case ("bulk_water") 
                neq = 5 
            case ("bulk_water_ox") 
                neq = 1
            case default
                print*,"Wrong value systype:  ",systype
                call error_handler(1,"set_size_neq")
        end select  

        neqint =int(neq,kind(neqint)) ! used for MPI func binding, MPI has no integer(8)
         
    end subroutine set_size_neq

    
    function BjerrumLenght(T)result(lb)

        use mathconst
        use physconst
        
        implicit none
        
        real(dp), intent(in) :: T     
        real(dp) :: lb

        lb=(elemcharge**2)/(4.0_dp*pi*dielectW*dielect0*kBoltzmann*T) ! bjerrum length in water=solvent in m
        lb=lb/1.0e-9_dp              ! bjerrum length in water in nm

    end function BjerrumLenght
        

    !     purpose: initialize all constants parameter 
    !     pre: first read_inputfile has to be called   
    
    subroutine init_constants()

        use globals
        use volume
        use random
        use physconst
        
        implicit none      
        
        real(dp) :: vA,vB, vAA, vAMPS, vPEG
        
        !  .. initializations of variables
 
        pi=acos(-1.0_dp)          ! pi = arccos(-1)
        itmax=2000                ! maximum number of iterations      
      
        !  .. charges  
        zNa   = 1                 ! valence positive charged ion
        zK    = 1                 ! valence positive charged ion
        zCa   = 2                 ! valence divalent positive charged ion
        zMg   = 2                 ! valence divalent positive charged ion
        zCl   =-1                 ! valence negative charged ion
        zFe2  = 2                ! valence divalent positive charged ion
        zFe3  = 3                ! valence divalent positive charged ion

        zpolA(1)=-1 ! A-
        zpolA(2)= 0 ! AH
        zpolA(3)= 0 ! ANa
        zpolA(4)= 1 ! ACa+
        zpolA(5)= 0 ! A2Ca
        
        zpolB(1)=-1 ! B-
        zpolB(2)= 0 ! BH
        zpolB(3)= 0 ! BNa
        zpolB(4)= 1 ! BCa+
        zpolB(5)= 0 ! B2Ca
        
        zpolAA(1)=-1 ! A-
        zpolAA(2)= 0 ! AH
        zpolAA(3)= 0 ! ANa
        zpolAA(4)= 1 ! ACa+
        zpolAA(5)= 0 ! A2Ca
        zpolAA(6)= 1 ! AMg+
        zpolAA(7)= 0 ! A2Mg
        zpolAA(8)= 0 ! AK

        !  .. radii
        !  .. ionic radii
        !  .. https://www.chemguide.co.uk/atoms/properties/atradius.html and http://abulafia.mt.ic.ac.uk/shannon/ptable.php

        RNa = 0.102_dp             ! radius of Na+ in nm
        RK  = 0.138_dp             ! radius of K+ in nm
        RCl = 0.181_dp             ! radius of Cl- in nm
        RCa = 0.106_dp             ! radius of Ca2+ in nm
        RMg = 0.072_dp             ! radius of Mg2+ in nm 
        RFe2= 0.078_dp             ! radius of Fe2+ in nm 
        RFe3= 0.064_dp             ! radius of Fe3+ in nm
        RO2 = 0.3_dp               ! radius of  O2 in nm

        print*,"Warning: init_contant:  value of volume radius of O2 ~ 0.3 nm!!"        
        
        ! .. volume
        
        vsol = 0.030_dp              ! volume water solvent molecule in (nm)^3

        vNa  = ((4.0_dp/3.0_dp)*pi*(RNa)**3)/vsol 
        vK   = ((4.0_dp/3.0_dp)*pi*(RK)**3)/vsol 
        vCl  = ((4.0_dp/3.0_dp)*pi*(RCl)**3)/vsol 
        vCa  = ((4.0_dp/3.0_dp)*pi*(RCa)**3)/vsol 
        vMg  = ((4.0_dp/3.0_dp)*pi*(RMg)**3)/vsol 
        vFe2 = ((4.0_dp/3.0_dp)*pi*(RFe2)**3)/vsol 
        vFe3 = ((4.0_dp/3.0_dp)*pi*(RFe3)**3)/vsol 
        vO2  = ((4.0_dp/3.0_dp)*pi*(RO2)**3)/vsol 

        vNaCl= (vNa+vCl)          ! contact ion pair
        vKCl = (vK+vCl)           ! contact ion pair
        
        ! .. volume polymer segments
        ! .. all volume scaled by vsol
        
        vAA  =  0.07448_dp/vsol ! volume based on VdW radii 
        vAMPS = 0.2134_dp/vsol
        vPEG  = 0.065_dp/vsol

        ! .. volume AA and vAMPS

        vA = vAA
        vB = vAMPS

        vpolA(1) = vA              ! vA-
        vpolA(2) = vA              ! vAH
        vpolA(3) = vA+vNa          ! vANa
        vpolA(4) = vA+vCa          ! vACa
        vpolA(5) = 2.0_dp*vA+vCa   ! vA2Ca
        
        vpolB(1) = vB              ! vB-
        vpolB(2) = vB              ! vBH
        vpolB(3) = vB+vNa          ! vBNa
        vpolB(4) = vB+vCa          ! vBCa
        vpolB(5) = 2.0_dp*vB+vCa   ! vB2Ca
        
        deltavA(1) = vpolA(1)+1.0_dp-vpolA(2) ! vA-+vH+-vAH
        deltavA(2) = vpolA(1)+vNa-vpolA(3)    ! vA-+vNa+-vANa+
        deltavA(3) = vpolA(1)+vCa-vpolA(4)    ! vA- + vCa2+ -vACa+
        deltavA(4) = 2.0_dp*vpolA(1)+vCa-vpolA(5) ! 2vA- + vCa2+ -vA2Ca
        
        deltavB(1) = vpolB(1)+1.0_dp-vpolB(2) ! vB-+vH+-vBH
        deltavB(2) = vpolB(1)+vNa-vpolB(3)    ! vB-+vNa+-vBNa+
        deltavB(3) = vpolB(1)+vCa-vpolB(4)    ! vB- +vCa2+ -vBCa+
        deltavB(4) = 2.0_dp*vpolB(1)+vCa-vpolB(5) ! 2vB- + vCa2+ -vB2Ca+

        ! dissociation constant AA and AMPS

        pKaA(1) =  5.0_dp
        pKaA(2) = -0.4_dp
        pKaA(3) =  1.0_dp
        pKaA(4) =  4.0_dp

        pKaB(1) = -2.0_dp
        pKaB(2) = -0.42_dp
        pKaB(3) = -0.72243_dp
        pKaB(4) = -10.0_dp

        ! .. other physical varaibles
        lsegPAA  = 0.36287_dp       ! segment length in nm
        lsegPAMPS = 0.545_dp        ! segment length in nm
        lsegPEG =  0.3_dp           ! segment length in nm 
              
        lsegA = lsegPAA            
        lsegB = lsegPAMPS 

        ! .. see also subroutine set_chain_properties 

        pKw = 14.0_dp                 ! water equilibruim constant
        Tref = 310.0_dp               ! temperature in Kelvin
        dielectW = 78.54_dp         ! dielectric constant water
        dielectP = 2.0_dp 

        seed = 435672              ! seed for random number generator

        call init_elect_constants(Tref)  

        if(systype=="brushborn") then 
            ! bornrad%pol, bonrrad%polCa and bornrad%polMg can not assign born radius of charge monomers yet, see init_dna

            bornrad%Na  = RNa
            bornrad%Cl  = RCl
            bornrad%K   = RK 
            bornrad%Ca  = RCa
            bornrad%Mg  = RMg
            bornrad%Hplus = radiussphere(vsol)
            bornrad%OHmin = radiussphere(vsol)
            bornrad%Fe2 = RFe2
            bornrad%Fe3 = RFe3
            
        endif    

        cuantas=max_confor

    end subroutine init_constants


    function volumesphere(radius)result(volume)

        use mathconst

        real(dp), intent(in) :: radius
        real(dp) :: volume

        volume=(4.0_dp/3.0_dp)*pi*(radius**3)
    
    end function
      

    function radiussphere(volume)result(radius)

        use mathconst

        real(dp), intent(in) :: volume
        real(dp) :: radius

        radius=(volume*3.0_dp/(4.0_dp*pi))**(1.0_dp/3.0_dp)
    
    end function
      

    ! Init variables specific for DNA used in systype=brush, brushborn etc 
    ! variable are  constants, deltaG and K for sytype=brush,bruhborn etc.
    ! pre : nsegtype, vsol,vpol, vNa etc and ismonomer_chargable need to be set 
    ! post : equlibrium constant and volume set for charge state of 
    !        carboxylic group of systype ==brushborn are initliazed 

    subroutine init_dna  

        use globals, only : nsegtypes,nseg,systype
        use chains, only : type_of_monomer_char,type_of_monomer,ismonomer_chargeable
        use physconst, only : Na
        use myutils, only : error_handler

        real(dp) :: vA    
        integer  :: i,tt,s,flag_one
        logical  :: isApresent
        integer  :: info

        ! determine segment type number of phosphate constaining segments

        do s=1, nseg
            if(type_of_monomer_char(s)=="P")  tA =type_of_monomer(s)
        enddo

        isApresent=(tA/=0) ! check if phosphate acid monomer is defined in list of typesfname

        if(.not.isApresent) then
            print*,"Warning in init_dna:"
            print*,"A momomer is not defined in typesfname"
            print*,"tA= ",tA
        endif    

        call read_pKds(pKaAA,info)
        if(info/=0) then 
            if(info==err_pKdfile_noexist) then 
                ! set equilibrium constants for acrylic acid 
                pKaAA(1)=  5.0_dp   ! AH    <=> A^- +H^+
                pKaAA(2)= -0.4_dp   ! ANa   <=> A^- +Na^+
                pKaAA(3)=  1.0_dp   ! ACa+  <=> A^- +Ca^++
                pKaAA(4)=  4.0_dp   ! A_2Ca <=> 2A^- + Ca^++
                pKaAA(5)=  1.0_dp   ! AMg^+ <=> A^- + Mg^++
                pKaAA(6)=  4.0_dp   ! A_2Mg <=> 2A^- + Mg^++
                pKaAA(7)= -0.4_dp   ! AK    <=> A^- + K^+
                pKaAA(8)= -0.4_dp   ! AFe^+  <=> A^- + Fe^++ pK for binding of Fe++ not determined for acrylic acid 
                pKaAA(9)=  4.0_dp   ! A2Fe   <=> 2A^- + Fe^++  
                pKaAA(10)= -0.4_dp  ! AFe^2+ <=> A^-  + Fe^+++  
                pKaAA(11)= 4.0_dp   ! A2Fe^+ <=> 2A^- + Fe^+++  
                pKaAA(12)= 4.0_dp   ! A3Fe   <=> 3A^- + Fe^+++  

            else ! errro
                print*,"Error in init_dna:"
                print*,"Failure to init pKaAA: info=",info
                stop
            endif    
        endif
        
        do i=1,size(K0aAA)
            KaAA(i) = 10.0_dp**(-pKaAA(i))  
            K0aAA(i) = KaAA(i)*(vsol*Na/1.0e24_dp)
        enddo

        if(systype/="nucl_ionbin_Mg".and. systype/="nucl_ionbin_MgA".and. systype/="nucl_ionbin_Fe") then ! old binding model
            K0aAA(4) = K0aAA(4)*(vsol*Na/1.0e24_dp) ! A2Ca
            K0aAA(6) = K0aAA(6)*(vsol*Na/1.0e24_dp) ! A2Mg
        endif 

        ! set volumes 
         
        vA=vpol(tA) 
        vpolAA(1) = vA              ! vA-
        vpolAA(2) = vA              ! vAH
        vpolAA(3) = vA+vNa          ! vANa
        vpolAA(4) = vA+vCa          ! vACa
        vpolAA(5) = 2.0_dp*vA+vCa   ! vA2Ca 
        vpolAA(6) = vA+vMg          ! vAMg       
        vpolAA(7) = 2.0_dp*vA+vMg   ! vA2Mg 
        vpolAA(8) = vA+vK           ! vAK 
        vpolAA(9) = vA+vFe2         ! vAFe2
        vpolAA(10) = 2.0_dp*vA+vFe2 ! vA2Fe2 
        vpolAA(11) = vA+vFe3        ! vAFe3
        vpolAA(12) = 2.0_dp*vA+vFe3 ! vA2Fe3 
        vpolAA(13) = 3.0_dp*vA+vFe3 ! vA3Fe3 
        
        
        deltavAA(1) = vpolAA(1)+1.0_dp-vpolAA(2)        ! vA- + vH+ - vAH
        deltavAA(2) = vpolAA(1)+vNa-vpolAA(3)           ! vA- + vNa+ - vANa
        deltavAA(3) = vpolAA(1)+vCa-vpolAA(4)           ! vA- + vCa2+ - vACa+
        deltavAA(4) = 2.0_dp*vpolAA(1)+vCa-vpolAA(5)    ! 2vA- + vCa2+ -vA2Ca 
        deltavAA(5) = vpolAA(1)+vMg-vpolAA(6)           ! vA- + vMg2+ - vAMg+
        deltavAA(6) = 2.0_dp*vpolAA(1)+vMg-vpolAA(7)    ! 2vA- + vMg2+ -vA2Mg
        deltavAA(7) = vpolAA(1)+vK-vpolAA(8)            ! vA- + vK+ - vAK
        deltavAA(8) = vpolAA(1)+vFe2-vpolAA(9)          ! vA- + vFe2+ - vAFe2
        deltavAA(9) = 2.0_dp*vpolAA(1)+vFe2-vpolAA(10)  ! 2vA- + vFe2+ -vA2Fe2
        deltavAA(10) = vpolAA(1)+vFe3-vpolAA(11)        ! vA- + vFe3+ - vAFe3(2+)
        deltavAA(11) = 2.0_dp*vpolAA(1)+vFe3-vpolAA(12) ! 2vA- + vFe3+ -vA2Fe3(1+) 
        deltavAA(12) = 3.0_dp*vpolAA(1)+vFe3-vpolAA(13) ! 3vA- + vFe3+ -vA2Fe3


        if(systype=="nucl_ionbin_Mg" .or. systype=="nucl_ionbin_MgA".or. systype=="nucl_ionbin_Fe") then
            call init_vPP(info)
            call error_handler(info,"init_vPP")
            call init_qpp()
        endif
            

        if(systype=="nucl_ionbin_Fe") then
            call init_vPPP(info)
            call error_handler(info,"init_vPPP")
            call init_qppp()
        endif
            

        ! determine if there is only one seg type is chargeable
        flag_one=0
        do tt=1,nsegtypes
            if(ismonomer_chargeable(tt)) flag_one=flag_one+1    
        enddo 

        if(flag_one==0) then
            print*,"Warning: init_dna: zero chargeable acid monomers types"
        else if(flag_one>1) then
            print*,"Warning: init_dna: more then one chargeable acid monomer types"
        endif    

        if(systype=="brushborn") then

            ! make zpolAA zero if monomor t=tA has zpol(1)=zpol(2)=0
            ! see routine read_pKas_and_zpol
            if(.not.ismonomer_chargeable(tA)) zpolAA=0
        
            ! assign born radius charged states of acrylic acid monomer  
            ! see also init_constants
       
            bornrad%pol   = radiussphere(vpolAA(1)*vsol)  ! A^-
            bornrad%polCa = radiussphere(vpolAA(4)*vsol)  ! ACa^+ 
            bornrad%polMg = radiussphere(vpolAA(6)*vsol)  ! AMg^+
            
        endif

    end subroutine init_dna
     

    function dielectric_constant_water(Temp) result(eps)
        implicit none

        real(dp), intent(in) :: Temp ! temperature in Kelvin
        real(dp) :: eps, Tc

        Tc=Temp-273.15_dp
        eps=(3.70886e4_dp - 8.2168e1_dp*Tc)/( 4.21854e2_dp + Tc)        !   Meissner and  Wentz 
        !eps= 87.7410_dp - 0.4008_dp*Tc + 9.398e-4_dp*Tc**2 - 1.410e-6_dp*Tc**3 ! Malmber and Maryott 

    end function

    subroutine solubility_oxygen_water(Temp, Hcp_s, Hxp_s)
       
        use physconst, only : Na

        real(dp), intent(in) :: Temp ! temperature in Kelvin
        real(dp), intent(inout) :: Hcp_s, Hxp_s

        real(dp) :: Tc, Tstd, B , Hcps_std

        Tstd = 298.15_dp      ! standard temperature
        Tc = Temp 
        Hcps_std = 1.3e-5_dp  ! unit mol/(m^3 Pa) at standard temperature
        B = 1500.0_dp         ! unit Kelvin  B = d ln(H_s)/d(1/T)
        ! source Sander Atmos. Phys. 23 10901 , 2023

        Hcp_s = Hcps_std * exp( B *(1.0_dp/Tc -1.0/Tstd)) ! mol/(m^3 Pa) 
        
        Hxp_s = Hcp_s * (Na/1.0e27_dp) * vO2 * vsol         ! 1/ Pa

        print*,"Hcp_s = ",Hcp_s," Hxp_s= ",Hxp_s

    end subroutine  solubility_oxygen_water



    subroutine solve_xbulk_with_oxygen(xbulk)

        use globals, only : neq, systype
        use physconst, only : Na
        use myutils, only : lenText, print_to_log, LogUnit

        type(moleclist), intent(inout) :: xbulk

        ! local variables
        real(dp) ::  xbulkidO2, xbulksalt, conc, tolerance
        real(dp),  dimension(:), allocatable ::  x, xguess
        character(len=15) :: systype_old
        logical :: issolution
        character(len=lenText) :: text

        
        call solubility_oxygen_water(Tref, Hcp_s, Hxp_s) 

        ! convert mol/(m^3 Pa)
        conc= pO2 * (Hcp_s/1000.0_dp)
        print*,"concentration in Mol/l=", conc, " pressure in Pa = " ,pO2  
        rhoidO2 = pO2 * (Hcp_s * Na / 1.0e27_dp ) ! density in 1/nm^3 :concentration of ideal solution in M=mol/l
      
        xbulkidO2 = rhoidO2 * vO2 * vsol ! xbulk = rho * Hcp_s also !!
        print*,"x volume fraction =", xbulkidO2, " pressure again=", Hxp_s* xbulkidO2 

        tolerance = 1.0e-10_dp
        systype_old = systype 
        systype = "bulk_water_ox"        ! set solver to fcnbulk
        call set_size_neq()              ! number of nonlinear equations
      
       !  call set_fcn()                 ! set fcnptr to correct fcn   
        allocate(x(neq))
        allocate(xguess(neq))
         
        ! guess
        xbulksalt=xbulk%Hplus +xbulk%OHmin +xbulk%Cl +xbulk%Na +xbulk%K+xbulk%NaCl+xbulk%KCl & 
            +xbulk%Ca +xbulk%Fe2 +xbulk%Fe3 +xbulk%Mg 
        x(1)= 2.0 ! 1.0_dp -(xbulksalt +xbulkidO2)           
        xguess(1)=x(1)
        
        ! call solver(x, xguess, tolerance, fnorm, issolution)
            
        !  .. return solution
            
        xbulk%sol = x(1)
        xbulk%O2 = rhoidO2 * (vO2 * vsol) * ( xbulk%sol ** vO2 )


         if(xbulk%sol<0) then
            text="xsol%bulk negative : stop program."
            call print_to_log(LogUnit,text) 
            print*,text
            stop
        endif   

        ! reset systype
        systype = systype_old  
        call set_size_neq()         ! number of nonlinear equations
        ! call set_fcn() 

        deallocate(x)
        deallocate(xguess)
        
    end subroutine solve_xbulk_with_oxygen


    subroutine solve_xbulk_with_ionbinding(xbulk)

        use globals, only : neq, systype
        use physconst, only : Na
        use myutils, only : lenText, print_to_log, LogUnit

        type(moleclist), intent(inout) :: xbulk

        real(dp),  dimension(:), allocatable ::  x, xguess
        character(len=15) :: systype_old
        logical :: issolution
        character(len=lenText) :: text

        systype_old=systype 
        systype="bulk_water"        ! set solver to fcnbulk
        call set_size_neq()         ! number of nonlinear equations
        allocate(x(neq))
        allocate(xguess(neq))
        
        x(1)=xbulk%Na
        x(2)=xbulk%Cl
        x(3)=xbulk%NaCl
        x(4)=xbulk%K
        x(5)=xbulk%KCl
        
        xguess(1)=x(1)
        xguess(2)=x(2)
        xguess(3)=x(3)
        xguess(4)=x(4)
        xguess(5)=x(5)
        
        ! call solver(x, xguess, tol_conv, fnorm, issolution)
        
        !     .. return solution
        
        xbulk%Na  =x(1)
        xbulk%Cl  =x(2)
        xbulk%NaCl=x(3)
        xbulk%K   =x(4)
        xbulk%KCl =x(5)

        ! .. reset of flags
        iter=0
        systype=systype_old         ! switch solver back
        call set_size_neq()         ! set number of non-linear equation  
        
        xbulk%sol=1.0_dp-xbulk%Hplus-xbulk%OHmin - xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca 

        if(xbulk%sol<0) then
            text="xsol%bulk negative : wrong pH and or salt concentration,stop program."
            call print_to_log(LogUnit,text) 
            print*,text
            !call MPI_FINALIZE(ierr)
            stop
        endif   

        ! reset systype
        systype = systype_old  
        call set_size_neq()         ! number of nonlinear equations
        ! call set_fcn() 

        deallocate(x)
        deallocate(xguess)

    end subroutine solve_xbulk_with_ionbinding

    subroutine init_elect_constants(Temp)
        
        use globals
        use volume, only : delta
        use physconst

        real(dp), intent(in) :: Temp ! temperature in Kelvin

        dielectW=dielectric_constant_water(Temp)

        lb = (elemcharge**2)/(4.0_dp*pi*dielectW*dielect0*kBoltzmann*Temp) ! bjerrum length in water=solvent in m
        lb = lb/1.0e-9_dp                           ! bjerrum length in water in nm
        constqW = delta*delta*(4.0_dp*pi*lb)/vsol   ! multiplicative constant Poisson Eq. 

        lb0 =(elemcharge**2)/(4.0_dp*pi*dielect0*kBoltzmann*Temp) ! bjerrum length in vacum in m
        lb0 = lb0/1.0e-9_dp                         ! bjerrum length in vacum in nm
        constq0 = delta*delta*(4.0_dp*pi*lb0)/vsol  ! multiplicative constant Poisson Eq. 
        constqE = 1.0_dp /( 8.0_dp *constqW)        ! factor in PDF
        constqWin = constqW                         ! assignment for  loop of dielect

        beta_times_e = elemcharge/(kBoltzmann*Temp) !  e/ kBT conversion factor  between potential in Volt-J/C to dimensionless potential 
    
    !    print*,"beta_times_e=",beta_times_e, "Temp=", Temp

    end subroutine init_elect_constants
   
    ! compute number density density polymer

    function init_denspol()result(denspol)

        use globals, only : nseg, nsize
        use volume, only : delta

        real(dp) :: denspol

        ! local 
        real(dp) :: vol

        vol = nsize*(delta**3)
        denspol = nseg*1.0_dp/vol   
 
    end function

    ! compute global volume fraction polymer-chain

    function init_xvolpol(rhopol)result(xvolpol)

        use globals, only : nsize, nsegtypes
        use volume, only : delta
    
        real(dp), intent(in) :: rhopol(:,:)
        real(dp) :: xvolpol


        ! local 
        real(dp) :: vol
        integer :: i,t


        vol=nsize*(delta**3) ! volume lattice
        
        xvolpol=0.0_dp
        do t=1,nsegtypes
            do i=1,nsize
                xvolpol=xvolpol+rhopol(i,t)*vpol(t)
            enddo
        enddo

        xvolpol = xvolpol/vol   
 
    end function


    ! computes ion strenth solution 
    ! input type(moleclist) :: xbulk
    ! return read(dp) :: IS 

    function ion_strength(xbulk)result(IS)

        use physconst, only : Na
        use myutils, only : lenText, print_to_log, LogUnit

        type(moleclist), intent(in) :: xbulk
        real(dp) :: IS

        type(moleclist) :: conc  ! concentratio ions 

        IS = 0.0_dp

        call init_zero_moleclist(conc)
        
        ! concentration in mol/l
        conc%Na = xbulk%Na * (1.0e24_dp/Na)/(vNa*vsol) 
        conc%K  = xbulk%K * (1.0e24_dp/Na)/(vK*vsol)
        conc%Mg = xbulk%Mg * (1.0e24_dp/Na)/(vMg*vsol)
        conc%Ca = xbulk%Ca * (1.0e24_dp/Na)/(vCa*vsol)
        conc%Fe2 = xbulk%Fe2 * (1.0e24_dp/Na)/(vFe2*vsol) 
        conc%Fe3 = xbulk%Fe3 * (1.0e24_dp/Na)/(vFe3*vsol)
        conc%Cl = xbulk%Cl * (1.0e24_dp/Na)/(vCl*vsol) 
        conc%Hplus = xbulk%Hplus * (1.0e24_dp/Na)/(vsol)
        conc%OHmin = xbulk%OHmin * (1.0e24_dp/Na)/(vsol)

        ! ion strength 
        IS = IS + conc%Na + conc%K + conc%Cl+  conc%Hplus +conc%OHmin  ! monovalent ions
        IS = IS + 4.0_dp * ( conc%Ca + conc%Mg + conc%Fe2)             ! divalent ions
        IS = IS + 9.0_dp *  conc%Fe3                                   ! trivalent ions
        IS = IS/ 2.0_dp 
        
    end function ion_strength

   
    !     purpose: initialize expmu needed by fcn 
    !     pre: first read_inputfile has to be called

    subroutine init_expmu_elect()
 
        use globals
        use physconst, only : Na
        use dielectric_const
        use myutils, only : print_to_log,LogUnit,lenText
        
        !   .. local variable
        
        integer :: i
        character(len=lenText) :: text
        
        real(dp) :: xNaClsalt          ! volume fraction of NaCl salt in bulk
        real(dp) :: xKClsalt           ! volume fraction of KCl salt in bulk
        real(dp) :: xCaCl2salt         ! volume fraction of CaCl2 salt in bulk
        real(dp) :: xMgCl2salt         ! volume fraction of MgCl2 salt in bulk
        real(dp) :: xFeCl2salt         ! volume fraction of "FeCl2" salt in bulk 
        real(dp) :: xFeCl3salt         ! volume fraction of "FeCl3" salt in bulk

        real(dp) :: xtmp, rhoqtmp
        real(dp) :: KaAA6              ! auxilary varialbe
        
        !     .. initializations of input dependent variables, electrostatic part 
        
        pHbulk=pH%val ! transfer pH value 

        cHplus = (10.0_dp)**(-pHbulk) ! concentration H+ in bulk
        pOHbulk = pKw -pHbulk       
        cOHmin  = (10.0_dp)**(-pOHbulk) ! concentration OH- in bulk
        
        xbulk%Hplus = (cHplus*Na/(1.0e24_dp))*(vsol) ! volume fraction H+ in bulk vH+=vsol
        xbulk%OHmin = (cOHmin*Na/(1.0e24_dp))*(vsol) ! volume fraction OH- in bulk vOH-=vsol
        
        xNaClsalt = (cNaCl*Na/(1.0d24))*((vNa+vCl)*vsol) ! volume fraction NaCl salt in mol/l
        
        if(pHbulk.le.7) then      ! pH<= 7
            xbulk%Na=xNaClsalt*vNa/(vNa+vCl)  
            xbulk%Cl=xNaClsalt*vCl/(vNa+vCl) +(xbulk%Hplus -xbulk%OHmin)*vCl  ! NaCl+ HCl
        else                      ! pH >7
            xbulk%Na=xNaClsalt*vNa/(vNa+vCl) +(xbulk%OHmin -xbulk%Hplus)*vNa ! NaCl+ NaOH  
            xbulk%Cl=xNaClsalt*vCl/(vNa+vCl)  
        endif
        

        xKClsalt = (cKCl*Na/(1.0e24_dp))*((vK+vCl)*vsol) ! volume fraction KCl salt 
        xbulk%K = xKClsalt*vK/(vK+vCl)  
        xbulk%Cl = xbulk%Cl+xKClsalt*vCl/(vK+vCl) 

        xFeCl2salt = (cFeCl2 * Na/(1.0e24_dp))*((vFe2+2.0_dp*vCl)*vsol) ! volume fraction FeCl2 
        xbulk%Fe2 = xFeCl2salt * vFe2/(vFe2+2.0_dp*vCl)
        xbulk%Cl = xbulk%Cl+ xFeCl2salt*2.0_dp*vCl/(vFe2+2.0_dp*vCl)

        xFeCl3salt = (cFeCl3 * Na/(1.0e24_dp))*((vFe3+3.0_dp*vCl)*vsol) ! volume fraction FeCl3 
        xbulk%Fe3 = xFeCl3salt * vFe3/(vFe3+3.0_dp*vCl)
        xbulk%Cl = xbulk%Cl+ xFeCl3salt*3.0_dp*vCl/(vFe3+3.0_dp*vCl)
        
        xCaCl2salt = (cCaCl2*Na/(1.0e24_dp))*((vCa+2.0_dp*vCl)*vsol) ! volume fraction CaCl2 
        xbulk%Ca = xCaCl2salt*vCa/(vCa+2.0_dp*vCl)
        xbulk%Cl = xbulk%Cl+ xCaCl2salt*2.0_dp*vCl/(vCa+2.0_dp*vCl)
        
        xMgCl2salt = (cMgCl2*Na/(1.0e24_dp))*((vMg+2.0_dp*vCl)*vsol) ! volume fraction MgCl2
        xbulk%Mg = xMgCl2salt*vMg/(vMg+2.0_dp*vCl)
        xbulk%Cl = xbulk%Cl+ xMgCl2salt*2.0_dp*vCl/(vMg+2.0_dp*vCl)

        xbulk%NaCl = 0.0_dp    ! no ion pairing
        xbulk%KCl = 0.0_dp     ! no ion pairing

        ! established if xO2 in solution  
        if(pO2>0.0_dp) then 
            
            call solve_xbulk_with_oxygen(xbulk)

            xtmp=1.0_dp -xbulk%Hplus -xbulk%OHmin -xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl & 
                -xbulk%Ca -xbulk%Fe2 -xbulk%Fe3 -xbulk%Mg -xbulk%O2 -xbulk%sol
            !print*,"xtmp=",xtmp
        else 
          
            xbulk%O2 = 0.0_dp        
            xbulk%sol=1.0_dp -xbulk%Hplus -xbulk%OHmin -xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl & 
                -xbulk%Ca -xbulk%Fe2 -xbulk%Fe3 -xbulk%Mg -xbulk%O2
            
            if(xbulk%sol<0) then
                text="xsol%bulk negative : wrong pH and or salt concentration,stop program."
                call print_to_log(LogUnit,text)
                print*,text
                !call MPI_FINALIZE(ierr)
                stop
            endif   
        endif    

        !  .. if Kion == 0 ion pairing !
        !  .. intrinstic equilibruim constant acid        
        !  Kion  = 0.246_dp ! unit 1/M= liter per mol !!!
        K0ionK  = KionK /(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        K0ionNa = KionNa/(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        
        if((KionNa/=0.0_dp).or.(KionK/=0.0_dp)) then 
            call solve_xbulk_with_ionbinding(xbulk)
        endif
         
        !  .. intrinstic equilibruim constants      
        do i=1,4
             KaA(i)  = 10.0_dp**(-pKaA(i)) ! experimental equilibruim constant acid 
             KaB(i)  = 10.0_dp**(-pKaB(i)) ! experimental equilibruim constant acid
             K0aA(i) = (KaA(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
             K0aB(i) = (KaB(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo
        !     .. rescale for i=4 2A- Ca <=> A2Ca
          
        K0aA(4) = (K0aA(4)*vsol)*(Na/1.0e24_dp)
        K0aB(4) = (K0aB(4)*vsol)*(Na/1.0e24_dp)
         
        ! pKa constant from file assigned from read_pKas_and_zpol

        Ka  = 10.0_dp**(-pKa)                       ! experimental equilibruim constant acid 
        K0a = (Ka*vsol)*(Na/1.0e24_dp)              ! intrinstic equilibruim constant 
 
        if(systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv".or.systype=="nucl_ionbin_Mg".or.& 
           systype=="nucl_ionbin_MgA".or.systype=="nucl_ionbin_Fe") then 
            Kaion  = 10.0_dp**(-pKaion)             ! experimental equilibruim ionbinding 
            K0aion = (Kaion*vsol)*(Na/1.0e24_dp)    ! intrinstic equilibruim 
        endif    

        if(systype=="brushborn") then 

            bornbulk%pol   = born(lb,bornrad%pol,-1)
            bornbulk%polCa = born(lb,bornrad%polCa,1)
            bornbulk%polMg = born(lb,bornrad%polMg,1)
            
            bornbulk%Hplus = born(lb,bornrad%Hplus,1)
            bornbulk%Na    = born(lb,bornrad%Na,zNa)
            bornbulk%K     = born(lb,bornrad%K,zK)
            bornbulk%Ca    = born(lb,bornrad%Ca,zCa)
            bornbulk%Mg    = born(lb,bornrad%Mg,zMg)
            bornbulk%Cl    = born(lb,bornrad%Cl,zCl)
            bornbulk%Fe2   = born(lb,bornrad%Fe2,zFe2)
            bornbulk%Fe3   = born(lb,bornrad%Fe3,zFe3)
            bornbulk%OHmin = born(lb,bornrad%OHmin,-1)

            expmu%Na    = (xbulk%Na   /(xbulk%sol**vNa))*exp(bornbulk%Na) 
            expmu%Cl    = (xbulk%Cl   /(xbulk%sol**vCl))*exp(bornbulk%Cl) 
            expmu%K     = (xbulk%K    /(xbulk%sol**vK) )*exp(bornbulk%K) 
            expmu%Ca    = (xbulk%Ca   /(xbulk%sol**vCa))*exp(bornbulk%Ca) 
            expmu%Mg    = (xbulk%Mg   /(xbulk%sol**vMg))*exp(bornbulk%Mg) 
            expmu%Fe2   = (xbulk%Fe2  /(xbulk%sol**vFe2))*exp(bornbulk%Fe2) 
            expmu%Fe3   = (xbulk%Fe3  /(xbulk%sol**vFe3))*exp(bornbulk%Fe3) 
            expmu%Hplus = (xbulk%Hplus/xbulk%sol) *      exp(bornbulk%Hplus)  
            expmu%OHmin = (xbulk%OHmin/xbulk%sol) *      exp(bornbulk%OHmin)  
            expmu%O2    = 0.0_dp ! no oxygen 

        else

            ! exp(beta mu_i) = (rhobulk_i v_i) / exp(- beta pibulk v_i) 

            expmu%Na    = xbulk%Na   /(xbulk%sol**vNa) 
            expmu%K     = xbulk%K    /(xbulk%sol**vK)
            expmu%Fe2   = xbulk%Fe2  /(xbulk%sol**vFe2)
            expmu%Fe3   = xbulk%Fe3  /(xbulk%sol**vFe3)
            expmu%Ca    = xbulk%Ca   /(xbulk%sol**vCa) 
            expmu%Mg    = xbulk%Mg   /(xbulk%sol**vMg) 
            expmu%Cl    = xbulk%Cl   /(xbulk%sol**vCl)
            expmu%NaCl  = xbulk%NaCl /(xbulk%sol**vNaCl)
            expmu%KCl   = xbulk%KCl  /(xbulk%sol**vKCl)
            expmu%Hplus = xbulk%Hplus/xbulk%sol ! vsol = vHplus 
            expmu%OHmin = xbulk%OHmin/xbulk%sol ! vsol = vOHmin 
            expmu%O2    = xbulk%O2   /(xbulk%sol**vO2)      

        endif    
              
        if(runtype=="rangepKd") then 
            
            pKaAA(6)=pKd%val ! this override value read in.

            KaAA6=10.0_dp**(-pKaAA(6))  
            K0aAA(6) = KaAA6*(vsol*Na/1.0e24_dp)
            if(systype/="nucl_ionbin_Mg".and.systype/="nucl_ionbin_sv_Mg") then  
                K0aAA(6) = K0aAA(6)*(vsol*Na/1.0e24_dp) ! A2Mg
            endif    
        endif   
              
        ! ionic strength 
        IS= ion_strength(xbulk)

        ! test

        xtmp = sum_value_moleclist(xbulk)

        if(abs(xtmp-1.0_dp)>eps_val) then 
            print*,"Error: in init_expmu_elect : sum xbulk not one sum xbulk =", xtmp
            stop
        endif    

        rhoqtmp= xbulk%Hplus -xbulk%OHmin +zCl*xbulk%Cl/vCl +xbulk%Na/vNa +xbulk%K/vK & 
                +zCa*xbulk%Ca/vCa +zFe2*xbulk%Fe2/vFe2 +zFe3*xbulk%Fe3/vFe3 +zMg*xbulk%Mg/vMg
        
        if(abs(rhoqtmp)>eps_val) then 
            print*,"Error: in init_expmu_elect : bulk not charge neutral : residual charge =", rhoqtmp
        endif

        
        ! .. end init electrostatic part 
        
    end subroutine init_expmu_elect


    subroutine init_expmu_neutral

        use precision_definition
        
        xbulk%sol=1.0_dp                  ! volume fraction solvent 

    end subroutine init_expmu_neutral



    ! inits chem potential 

    subroutine init_vars_input()

        use globals, only : systype
    
        select case (systype)
        case ("elect")
            call init_expmu_elect()
            call set_VdWepsAAandBB() ! special assigment of VdWepsAA etc  
            call set_VdWeps_scale(VdWscale)
            call set_energychainLJ_scale(VdWscale)
            call set_dielect_scale(dielectscale)
        case ("neutral","neutralnoVdW")
            call init_expmu_neutral()   
            call set_VdWeps_scale(VdWscale)
            call set_energychainLJ_scale(VdWscale)
        case ("brush_mul","brush_mulnoVdW") 
            call init_expmu_elect() 
            call set_VdWeps_scale(VdWscale)
            call set_energychainLJ_scale(VdWscale) 
            call set_dielect_scale(dielectscale)    
        case ("brushdna","nucl_ionbin","nucl_ionbin_sv","nucl_ionbin_Mg","nucl_ionbin_MgA","nucl_ionbin_Fe") 
            call init_dna() 
            call init_expmu_elect()
            call set_VdWeps_scale(VdWscale)
            call set_energychainLJ_scale(VdWscale)
            call set_dielect_scale(dielectscale)
        case ("nucl_neutral_sv") 
            ! call init_dna() 
            call init_expmu_neutral()
            call set_VdWeps_scale(VdWscale)
            call set_energychainLJ_scale(VdWscale)
        case("brushborn") 
            call init_dna()
            call init_expmu_elect()  
            call set_VdWeps_scale(VdWscale)
            call set_energychainLJ_scale(VdWscale)
        case ("nonuclcp") 
            call init_expmu_elect()
        case default   
            print*,"Error: systype incorrect at init_vars_input" 
            print*,"Wrong value systype : ", systype
            print*,"stopping program"
            stop
         end select
             
    end subroutine init_vars_input

   
    subroutine allocate_chain_parameters
        
        use globals, only : nsegtypes
        
        !  allocate array depending on nsegtypes

        allocate(vpol(nsegtypes))       !  volume polymer segments, all volume scaled by vsol
        allocate(pKa(nsegtypes))        !  equilibrium constants
        allocate(pKaion(nsegtypes,4))   !  equilibrium constants  ion binding
        allocate(Ka(nsegtypes)) 
        allocate(Kaion(nsegtypes,4))  
        allocate(K0a(nsegtypes))  
        allocate(K0aion(nsegtypes,4))  
        allocate(zpol(nsegtypes,2))     !  charge of segment of given type       
        allocate(qpol(nsegtypes))       !  total charge of polymer type 
        allocate(avfdis(nsegtypes))     !  average fraction of charge of polymer type 
        allocate(avgdisA(nsegtypes,4))  !  fraction of acidic AA in state A, AH, ANa or AK
        allocate(avgdisB(nsegtypes,3))  !  fraction of basic AA in state BH, B, BCl
        allocate(lsegAA(nsegtypes))

    end subroutine allocate_chain_parameters

    !   call to init_chain_parameters 
    !   post/after  allocate_chain_parameters and init_constants    

    subroutine init_chain_parameters
        
        call init_volume_pol 
        call init_pKas_and_zpol
        call init_lseg  ! init segment length
        call init_pKaions

    end subroutine init_chain_parameters


    subroutine init_volume_pol

        use globals, only : nsegtypes

        call read_volume_pol(vpol,vsol,vpolfname, nsegtypes)
        
    end subroutine init_volume_pol
        

    subroutine allocate_vnucl(nelemtypes)

        use globals, only : nsegtypes,systype

        integer, intent(in) :: nelemtypes

        if(systype=="nucl_neutral_sv".or.systype=="nucl_ionbin_sv".or.systype=="nucl_ionbin_Mg".or.&
           systype=="nucl_ionbin_MgA" .or. systype=="nucl_ionbin_Fe") then 
            allocate(vnucl(nelemtypes,nsegtypes))
        endif    

    end subroutine allocate_vnucl  


    subroutine allocate_vnucl_type(nelemtypes)

        use globals, only : systype

        integer, intent(in) :: nelemtypes

        if(systype=="nucl_neutral_sv".or.systype=="nucl_ionbin_sv".or.systype=="nucl_ionbin_Mg".or.&
           systype=="nucl_ionbin_MgA".or. systype=="nucl_ionbin_Fe") then 
            allocate(vnucl_type(nelemtypes))
            allocate(vnucl_type_char(nelemtypes))
            allocate(vnucl_type_isChargeable(nelemtypes))
        endif    

    end subroutine allocate_vnucl_type  

    subroutine init_vnucl_type(info) 

        integer, intent(inout) :: info

        call read_vnucl_type(vnucl_type,vnucl_type_char,vnucl_type_isChargeable,vnuclfname,info)

    end subroutine init_vnucl_type

    subroutine init_vnucl

        use globals, only : systype

        if(systype=="nucl_neutral_sv") vnucl=0.0_dp !vnucl=0.01_dp ! init
        if(systype=="nucl_ionbin_sv") vnucl=0.0_dp 
        if(systype=="nucl_ionbin_Mg") vnucl=0.0_dp
        if(systype=="nucl_ionbin_MgA") vnucl=0.0_dp
        if(systype=="nucl_ionbin_Fe") vnucl=0.0_dp

    end subroutine init_vnucl
       

    subroutine init_pKas_and_zpol

        use globals, only : nsegtypes

        pKa  = 0.0_dp
        zpol = 0.0_dp
    
        call read_pKas_and_zpol(pKa,zpol,pKafname, nsegtypes) 
       
    end subroutine init_pKas_and_zpol


    subroutine init_pKaions

        use globals, only : nsegtypes, systype

        pKaion = 0.0_dp
        
        if(systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv".or.systype=="nucl_ionbin_Mg".or.&
           systype=="nucl_ionbin_MgA".or. systype=="nucl_ionbin_Fe") then 
            call read_pKaions(pKaion,zpol,pKaionfname, nsegtypes) 
        endif    
       
    end subroutine init_pKaions


    ! Warning !!this routine is semi redundant because lseg not used in VdW yet !!!

    subroutine init_lseg 

       use globals, only :nsegtypes
       
       call read_lseg(lsegAA,lsegfname, nsegtypes)

    end subroutine init_lseg


    ! vPP is equal to volume of indivual phosphate of pair when NOT in a briged state 
    ! vPP in a bridged state i.e. state= Phos2Mg, Phos2Fe2 or Phos2Fe3 , equal to totat volume of pairs
    ! Meaning e.g vPP(Phos2Mg)/2 is volume per phosphate 
    ! vPP uses vpol(tphos) in units of terms of vsol 
    ! pre:  vpol and vsol and tA need to be set 
    ! post: vPP 

    subroutine init_vPP(info)

        integer, intent(inout) :: info
        integer :: tPhos
        real(dp), parameter :: eps_vpol=1.0e-5_dp

        tPhos = tA
        info=0

        if(tPhos==0) then
            info=err_error
            return
        endif

        if(abs(vpol(tPhos))<eps_vpol) then 
            info=err_error
            return
        endif    

        vPP(Phos) = vpol(tPhos) * vsol
        vPP(PhosH) = vpol(tPhos) * vsol
        vPP(PhosK) = (vpol(tPhos)+vK ) * vsol
        vPP(PhosNa) = (vpol(tPhos)+vNa ) * vsol 
        vPP(PhosMg) = (vpol(tPhos)+vMg ) * vsol
        vPP(Phos2Mg) = (2.0_dp*vpol(tPhos)+vMg ) * vsol
        vPP(PhosFe2) = (vpol(tPhos)+vFe2 ) * vsol
        vPP(Phos2Fe2) = (2.0_dp*vpol(tPhos)+vFe2 ) * vsol
        vPP(PhosFe3) = (vpol(tPhos)+vFe3 ) * vsol
        vPP(Phos2Fe3) = (2.0_dp*vpol(tPhos)+vFe3 ) * vsol  
       !  vPP(Phos3Fe3) = (3.0_dp*vpol(tPhos)+vFe3 ) * vsol
 
    end subroutine   init_vPP 

    ! vPPP equal vPP for nonovalent condensation reaction : singlet
    ! vPPP is not vPP for pair and triplet
    ! vPPP is volume of one phosphate that is part of a triplet in a given chamical state
    ! vPPP indivual volumer 
    ! for doublet 1/2 of the total volume of total volume doublet state 
    ! for triplet 1?3 of the total volume of total volume triplet state  
    ! vPPP uses vpol(tphos) in units of terms of vsol 
    ! pre:  vpol and vsol and tA need to be set 
    ! post: vPPP 

    subroutine init_vPPP(info)

        integer, intent(inout) :: info
        integer :: tPhos
        real(dp), parameter :: eps_vpol=1.0e-5_dp

        tPhos = tA
        info = 0

        if(tPhos==0) then
            info=err_error
            return
        endif

        if(abs(vpol(tPhos))<eps_vpol) then 
            info=err_error
            return
        endif    

        ! singlet state of triplet

        vPPP(Phos) = vpol(tPhos) * vsol
        vPPP(PhosH) = vpol(tPhos) * vsol
        vPPP(PhosK) = (vpol(tPhos)+vK ) * vsol
        vPPP(PhosNa) = (vpol(tPhos)+vNa ) * vsol 
        vPPP(PhosMg) = (vpol(tPhos)+vMg ) * vsol
        vPPP(PhosFe2) = (vpol(tPhos)+vFe2 ) * vsol
        vPPP(PhosFe3) = (vpol(tPhos)+vFe3 ) * vsol

        ! doublet state of triplet 
        vPPP(Phos2Mg) = (2.0_dp*vpol(tPhos)+vMg ) * vsol /2.0_dp
        vPPP(Phos2Fe2) = (2.0_dp*vpol(tPhos)+vFe2 ) * vsol /2.0_dp
        vPPP(Phos2Fe3) = (2.0_dp*vpol(tPhos)+vFe3 ) * vsol  /2.0_dp

        ! triplet state of triplet
        vPPP(Phos3Fe3) = (3.0_dp*vpol(tPhos)+vFe3 ) * vsol/3.0_dp
 
    end subroutine   init_vPPP

    subroutine init_qPP()
        
        ! charge of phophate as part of pair
        ! singlet charge of pair
        qPP(Phos) = -1
        qPP(PhosH) = 0
        qPP(PhosK) = 0
        qPP(PhosNa) = 0
        qPP(PhosMg) = 1
        qPP(PhosFe2) = 1
        qPP(PhosFe3) = 2

        ! doublet charge of pair
        qPP(Phos2Mg) = 0
        qPP(Phos2Fe2) = 0
        qPP(Phos2Fe3) = 1

    end subroutine   init_qPP
    
    subroutine init_qPPP()

        ! charge of phospahte as part of triplet state 
        ! singlet charge of triplet  
        
        qPPP(Phos) = -1.0_dp
        qPPP(PhosH) = 0.0_dp
        qPPP(PhosK) = 0.0_dp
        qPPP(PhosNa) = 0.0_dp
        qPPP(PhosMg) = 1.0_dp
        qPPP(PhosFe2) = 1.0_dp 
        qPPP(PhosFe3) = 2.0_dp

        ! doublet charges of triplet 
 
        qPPP(Phos2Mg)  = 0.0_dp
        qPPP(Phos2Fe2) = 0.0_dp
        qPPP(Phos2Fe3) = 0.5_dp ! half of qPP(Phos2Fe3) distribute over pair

        ! charge of phospahte that for a trivalent triplet state 
        qPPP(Phos3Fe3) = 0.0_dp


    end subroutine   init_qPPP


    !  .. assign vpol from values in file named filename
    !  .. values vpol are normalized by vsol
    !  .. checked if file exists- content and length not checked     

    subroutine read_lseg(lsegAA,filename, ntypes)

        use  myutils

        implicit none 
        
        !     .. arguments 
        real(dp), intent(inout)  :: lsegAA(:) 
        character(40), intent(in) :: filename  
        integer,  intent(in) :: ntypes 

        integer :: ios, un  ! un = unit number
        integer :: t 
        character(80) :: istr,str

        !     .. reading in of variables from file
        open(unit=newunit(un),file=filename,iostat=ios,status='old')
        if(ios/=0 ) then
            write(istr,'(I2)')ios
            str='Error read_lseg: opening file '//trim(adjustl(filename))//' : iostat = '//istr
            print*,str
            stop
        endif
    
        do t=1,ntypes   
            read(un,*)lsegAA(t)
        enddo    

        close(un)
   
    end subroutine read_lseg


    !  .. assign vpol from values in file named filename
    !  .. values vpol are normalized by vsol
    !  .. checked if file exists- content and length not checked     

    subroutine read_volume_pol(vpol,vsol,filename, ntypes)

        use  myutils

        implicit none 
        
        ! .. arguments 
        real(dp), intent(inout)  :: vpol(:) 
        real(dp), intent(in) :: vsol 
        character(40), intent(in) :: filename  
        integer,  intent(in) :: ntypes 

        ! .. local variables
        integer :: ios, un  ! un = unit number
        integer :: t 
        character(80) :: istr,str

        ! .. reading in of variables from file
        open(unit=newunit(un),file=filename,iostat=ios,status='old')
        if(ios/=0 ) then
            write(istr,'(I2)')ios
            str='Error read_volume_pol: opening file '//trim(adjustl(filename))//' : iostat = '//istr
            print*,str
        !  .. local variables
            stop
        endif
    
        do t=1,ntypes   
            read(un,*)vpol(t)
            vpol(t)=vpol(t)/vsol
        enddo    

        close(un)
   
    end subroutine read_volume_pol

    
    !  .. assign  pKa and zpol from values in file named filename
   
    subroutine read_pKas_and_zpol(pKa,zpol,filename, ntypes)
   
        use  myutils
        implicit none 
        
        !     .. arguments 
        real(dp), intent(inout) :: pKa(:)
        integer, intent(inout) ::  zpol(:,:)
        character(lenfname), intent(in) :: filename
        integer,  intent(in) :: ntypes

        !      .. local variables
        integer :: ios, un  ! un = unit number
        integer :: t 
        character(80) :: istr,str

        !     .. reading in of variables from file
        open(unit=newunit(un),file=filename,iostat=ios,status='old')
        if(ios/=0 ) then
            write(istr,'(I2)')ios
            str='Error read_pKas_and_zpol opening file '//trim(adjustl(filename))//' : iostat = '//istr
            print*,str
            stop
        endif
        
        do t=1,ntypes   
            read(un,*)pKa(t),zpol(t,1),zpol(t,2)
        enddo    
        
        close(un)


    end subroutine read_pKas_and_zpol

    !  .. assign  pKaion from values in file named filename
   
    subroutine read_pKaions(pKaion,zpol,filename, ntypes)
   
        use  myutils
        implicit none 
        
        !     .. arguments 
        real(dp), intent(inout) :: pKaion(:,:)
        integer, intent(in) ::  zpol(:,:)
        character(lenfname), intent(in) :: filename
        integer,  intent(in) :: ntypes

        !      .. local variables
        integer :: ios, un  ! un = unit number
        integer :: t, k, info ,zpol1,zpol2
        character(len=80) :: istr,str
        character(len=2)  :: istr2
        logical :: exist

        info = 0

        inquire(file=filename,exist=exist)
        if(exist) then
            !     .. reading in of variables from file
            open(unit=newunit(un),file=filename,iostat=ios,status='old')
            if(ios/=0 ) then
                write(istr,'(I2)')ios
                str='Error read_pKaions opening file '//trim(adjustl(filename))//' : iostat = '//istr
                print*,str
                info =  err_pKdfile
            endif
            
        else
            str='pKaions file does not exist: '//trim(adjustl(filename))
            print*,str
            info= err_pKdfile_noexist
        endif
        
        if(info==0) then        
            do t=1,ntypes 
                if((zpol(t,1)==0).and.(zpol(t,2)==-1)) then !acid 
                    read(un,*,iostat=ios)pKaion(t,1),zpol1,zpol2,(pKaion(t,k),k=2,3)
                !    print*,pKaion(t,1),zpol1,zpol2,(pKaion(t,k),k=2,3)
                else if((zpol(t,1)==1).and.(zpol(t,2)==0)) then ! base
                    read(un,*,iostat=ios)pKaion(t,1),zpol1,zpol2,pKaion(t,2)
                !    print*,pKaion(t,1),zpol1,zpol2,pKaion(t,2)
                else
                    read(un,*,iostat=ios)
                endif   
                if(ios>0) then 
                    write(istr2,'(I2)')t
                    str='Error occur reading line '//istr2//' of pKaion inputfile'
                    info =  err_pKderror
                endif    
            enddo    
                
            close(un)
 
        else  ! something went wrong
            stop
        endif    
            

    end subroutine read_pKaions

    ! Read pKd for acid group tA including acid-base equilbrium Na condensation etc
    ! Four return values 
    ! info=0 , correct, 
    ! info= err_pKdfile,err_pKdfile_noexist or err_pKderror : failure

    subroutine read_pKds(pKd,info)
   
        use  myutils
        implicit none 
        
        !     .. arguments 
        real(dp), intent(inout) :: pKd(:)
        integer, intent(out) :: info 
       
        !      .. local variables
        integer :: ios, un  ! un = unit number
        integer :: line, maxline
        character(len=80) :: str
        logical :: exist
        character(len=10) :: fname

        info = 0
        !     .. reading in of variables from file
        write(fname,'(A10)')'pKdacid.in'
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un),file=fname,status='old',iostat=ios)
            if(ios >0 ) then
                print*, 'Error opening file : iostat =', ios
                info = err_pKdfile
                return
            endif
        else
            info= err_pKdfile_noexist 
            return
        endif
        
        line=0
        ios=0
        maxline=size(pKd) ! 12
        
        do while (line<maxline.and.ios==0)
            line=line+1
            read(un,*,iostat=ios)pKd(line)
            !print*,line,pKd(line)
        enddo
        
        close(un)

    
        if(line/=maxline.or.ios/=0) then 
            str="reached end of file before all elements read or ios error"
            print*,str
            str="read file "//trim(adjustl(fname))//" failed"
            print*,str
            info = err_pKderror
            return
        endif
        

    end subroutine read_pKds

    ! Processes input file containing the type, volume and chargeable of volume elements of AA
    ! input: character(lenfname) fname = filename of input file    
    ! input/output: real(dp)         vnucl_type(:)
    !               character(len=3) vnucl_type_char(:)
    !               logical vnucl_type_isChargeable(:)
    ! output : integer info
    !          return value  info=0 , correct, 
    !                        info/=0, failure info= err_error or err_file_noexist

    subroutine read_vnucl_type(vnucl_type,vnucl_type_char,vnucl_type_isChargeable,fname,info)

        use myutils
        use globals, only : DEBUG

        ! .. arguments 
        real(dp), intent(inout), allocatable          :: vnucl_type(:)
        character(len=3), intent(inout),allocatable   :: vnucl_type_char(:)
        logical, intent(inout),allocatable            :: vnucl_type_isChargeable(:) 
        character(lenfname), intent(in)               :: fname  
        integer, intent(out)                          :: info   
        
        ! .. local variables
        integer :: ios, un, line, maxline, i
        integer :: nelemtypes
        character(len=90) :: istr,str
        logical :: exist
        character(len=3) :: vol_char
        integer :: vol_int
        real(dp) :: vol
        logical  :: isChargeable

        info = 0 ! return 

        ! .. reading in of variables from file
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un),file=fname,iostat=ios,status='old')
            if(ios/=0 ) then
                write(istr,'(I2)')ios
                str='Error read_vnucl_type: opening file '//trim(adjustl(fname))//' : iostat = '//istr
                print*,str
                info = err_error
            endif
        else 
            str='Error read_vnucl_type: file '//trim(adjustl(fname))//' : does not exist'
            print*,str
            info= err_file_noexist 
            return
        endif    
    
        ios  = 0
        read(un,*,iostat=ios)nelemtypes                    ! read first line
        call allocate_vnucl_type(nelemtypes) ! allocate vnucl_type,vnucl_type_char
        
        call allocate_vnucl(nelemtypes) ! allocate vnucl this is wrong palace to allocate it th is is a side effect 
        call init_vnucl()

        line = 0
        ios  = 0
        maxline = nelemtypes
        do while (line<maxline.and.ios==0)
            line=line+1
            read(un,*,iostat=ios)vol_char,vol_int,vol,isChargeable
            vnucl_type(vol_int)=vol 
            vnucl_type_char(vol_int)=vol_char
            vnucl_type_isChargeable(vol_int)=isChargeable
        enddo  

        if(line/=maxline.or.ios/=0) then 
            str="reached end of file before all elements read or ios error"
            print*,str
            str="read file "//trim(adjustl(fname))//" failed"
            print*,str
            info = err_error
            return
        endif

        close(un)
    
        if(DEBUG) then
            print*,"Module : parameters : read_vnucl_type"
            print*,"i          vnucl_type        vnucl_type_char" 
            do i=1,nelemtypes
                print*,i," ",vnucl_type(i), " ",vnucl_type_char(i)
            enddo
        endif        

    end subroutine read_vnucl_type



    subroutine allocate_isrhoselfconsistent(info)
    
        use globals, only : nsegtypes

        integer,  intent(out), optional :: info

        integer :: ier

        ier = 0 
        if (present(info)) info = 0

        if (.not. allocated(isrhoselfconsistent))  then 
            allocate(isrhoselfconsistent(nsegtypes),stat=ier)
        endif        

        if(ier/=0) then 
            print*,'Allocation error: allocate_isrhoselfconsistent failed'
            if (present(info)) info = ier

        endif    
       
    end subroutine allocate_isrhoselfconsistent

    ! Determines for every segement is need to selfconsistently solved
    ! Internally it also determines the segment type number of A constaining segments
    ! this is also done in routine init_dna 
    ! isrhoselfconsistent needs to be know before neq can be  determined !!!
    ! pre: Vdweps and isVdW is allready intialized
    ! post isrhoselfconsistent is set  

    subroutine make_isrhoselfconsistent(info)

        use globals, only : nsegtypes,nseg,systype
        use chains, only : type_of_monomer_char,type_of_monomer
        use myutils, only : print_to_log, LogUnit, lenText
 
        !     .. arguments 
    
        integer,  intent(out), optional :: info

        integer :: info_alloc
        integer :: i, s
        integer :: ttAA, ttP ! local location of A and P segment
        character(len=lenText) :: text
        
        if(present(info)) info=0

        call allocate_isrhoselfconsistent(info_alloc)
        if(info_alloc/=0) then 
            print*,"Error: in allocate_isrhoselfconsistent"
            if(present(info)) info=info_alloc
            return
        endif    

        ! init all element of isrhoselfconsistent to false
        do i=1,nsegtypes
            isrhoselfconsistent(i)=.false.
        enddo

        ! Mg/Ca self consistent equation

        ! determine segment type number of P = phosphate dsDNA or ssDNA  
        ttP=0
        ttAA=0
        do s=1, nseg
            if(type_of_monomer_char(s)=="AA") ttAA=type_of_monomer(s)
            if(type_of_monomer_char(s)=="P")  ttP=type_of_monomer(s)
        enddo
           
        if(systype=='brush_born'.or. systype=='brush_dna') then 
           if(ttAA>0) isrhoselfconsistent(ttAA)=.true.  ! check condition ttA==0
           if(ttP>0) isrhoselfconsistent(ttP)=.true.    ! check condition ttP==0
        endif    
        
        if(systype=='nucl_ionbin'.or.systype=='nucl_ionbin_sv') then

           ! if [Mg]/=0 for segment type==ttp is selfc nsistent, but the equations are not correct in 3D.  
           ! Use systype=="nucl_ionbin instead, Here a warning is issued 
           text="Warning: in make_isrhoselfconsistent: "
           call print_to_log(LogUnit,text)
           print*,text
           text="Mg concentration must be zero for systypes: "
           text=trim(adjustl(text))//trim(adjustl(systype))
           call print_to_log(LogUnit,text)
           print*,text
        endif
 
        ! check that we do not have simultenoeus A=Acrylic acid and P=phosphate
        if((ttAA>0).and.(ttP>0)) then
            text="Error: in make_isrhoselfconsistent: both A and P segments defined"
           call print_to_log(LogUnit,text)
           print*,text
        endif   

        !print*,"ttAA=",ttAA,"ttP=",ttP,ttAA>0

    end subroutine make_isrhoselfconsistent


    ! special assignment for runtype==rangeVdWeps
    ! pre VdWeps and VdWepsin allocated 

    subroutine set_VdWepsin

        VdWepsin= VdWeps

    end subroutine set_VdWepsin


    ! special assignment for certain systype values
    ! use VdWeps to assigns specific values VdWepsAA etc
    ! pre VdWeps and VdWepsin allocated 
    subroutine set_VdWepsAAandBB

        use globals, only : systype
            
        select case (systype)
        case ("elect")  ! diblock copolymer lseg determined in cadenas_sequence      
            VdWepsAA = VdWeps(1,1) 
            VdWepsAB = VdWeps(1,2) 
            VdWepsBB = VdWeps(2,1) 
        case ("neutral","neutralnoVdW","brush_mul","brush_mulnoVdW","brushvarelec","brushborn","brushdna",&
                "nucl_ionbin","nucl_ionbin_sv","nucl_neutral_sv","nucl_ionbin_Mg","nucl_ionbin_MgA","nucl_ionbin_Fe",&
                "nonuclcp")
        case default
            print*,"Error: in set_VdWepsAAandBB, systype=",systype
            print*,"stopping program"
            stop
        end select  

    end subroutine set_VdWepsAAandBB

    ! special assignment for runtype==rangeVdWeps
    
    subroutine set_energychainLJ_scale(VdWscale)

        use globals, only : cuantas
        use chains, only : energychainLJ, energychainLJ0
        use globals, only : runtype 

        TYPE(looplist), intent(in) :: VdWscale

        integer :: c

        if(runtype=="rangeVdWeps") then 
            do c=1,cuantas
                energychainLJ(c)=VdWscale%val*energychainLJ0(c)
            enddo
        endif    


    end subroutine set_energychainLJ_scale

    ! special assignment for runtype==rangeVdWeps
    
    subroutine set_VdWeps_scale(VdWscale)

        TYPE(looplist), intent(in) :: VdWscale

        VdWeps=VdWscale%val*VdWepsin

        call set_VdWepsAAandBB

    end subroutine set_VdWeps_scale

    
    ! special assignment of constqw for runtype==rangedielect
    
    subroutine set_dielect_scale(dielectscale)

        use globals, only : runtype

        type(looplist), intent(in) :: dielectscale

        if(runtype=="rangedielect") then 
            constqW=constqWin/dielectscale%val
        endif    

    end subroutine set_dielect_scale

 end module parameters
