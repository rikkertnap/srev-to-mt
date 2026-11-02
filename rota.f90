! ------------------------------------------------------------|
! rotates a given chains conformation                         |  
! pre: xend = input chain                                     |
! post: xendr= rotated chain                                  |
!       return= is_postive                                    |
! ------------------------------------------------------------|

module chain_rotation

    use precision_definition
  
    implicit none

    real(dp), parameter  :: epsAngle = 0.0000000000001_dp

    private :: epsAngle 

contains  

!Euler rotation Z1X2Z3

function rotation(xend,xendr,nseg)result(is_positive_z)
  
    use random
    use mathconst
    
    implicit none

    integer,  intent(in)    :: nseg  
    real(dp), intent(in)    :: xend(:,:)
    real(dp), intent(inout) :: xendr(:,:)  
    logical                 :: is_positive_z  ! output 
  
    ! ..local argument 

    integer :: i
    real(dp)  :: fac,fac1,fac2,sbe,cbe,sal,cal,sga,dist
    real(dp)  :: alfa,gama,cga,a,b,c
    
    fac=rands(seed)
    fac1=rands(seed)
    fac2=rands(seed)
    alfa=fac*2.0_dp*pi
    cbe=fac1*2.0_dp-1.0_dp
    gama=fac2*2.0_dp*pi
  
    sbe=(1.0_dp-cbe**2)**0.5_dp
    cal=cos(alfa)
    sal=sin(alfa)
    cga=cos(gama)
    sga=sin(gama)
  
    do i=1,nseg+1     ! rotation of xend result stored in xendr
        a=xend(1,i)
        b=xend(2,i)
        c=xend(3,i)
        xendr(1,i)=a*(-cbe*sal*sga+cal*cga) -b*(cbe*sal*cga+cal*sga)+c*sbe*sal
        xendr(2,i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
        xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe 
    enddo
  
    is_positive_z=.true.

    i=2 ! i=1 skipped because (0,0,0)

    do while((i.le.(nseg+1)).and.(is_positive_z)) 
        if(xendr(1,i)<=0.0_dp) is_positive_z=.false.
        i=i+1
    enddo

end function rotation


! rotation only around x-axis
! Euler rotatation X1Z2X3 

function rotationXaxis(xend,xendr,nseg)result(is_positive_z)
  
    use random
    use mathconst
    
    implicit none

    integer,  intent(in)    :: nseg  
    real(dp), intent(in)    :: xend(:,:)
    real(dp), intent(inout) :: xendr(:,:)  
    logical                 :: is_positive_z  ! output 
  
    ! ..local argument 

    integer :: i
    real(dp)  :: fac,fac1,fac2,sbe,cbe,sal,cal,sga
    real(dp)  :: alfa,gama,cga,a,b,c
    
    fac = rands(seed)
    fac1= rands(seed)
    fac2= rands(seed)

    alfa=fac*2.0_dp*pi
    cbe=1.0_dp ! fac1*2.0_dp-1.0_dp
    gama=fac2*2.0_dp*pi
  
    sbe=(1.0_dp-cbe**2)**0.5_dp
    cal=cos(alfa)
    sal=sin(alfa)
    cga=cos(gama)
    sga=sin(gama)
  
    do i=1,nseg+1     ! rotation of xend result stored in xendr
        a=xend(1,i)
        b=xend(2,i)
        c=xend(3,i)

        xendr(1,i)=a ! a*cbe + b*(-cga*sbe)+c*sbe*sga
        xendr(2,i)=a*cal*sbe +b*(cal*cbe*cga-sal*sga) +c*(-cga*sal-cal*cbe*sga)
        xendr(3,i)=a*sbe*sga +b*(cal*sga+cbe*cga*sal) +c*( cal*cga-cbe*sal*sga)
    enddo
  
    is_positive_z=.true.
    i=2 
    do while((i.le.(nseg+1)).and.(is_positive_z)) 
        if(xendr(3,i)<=0.0_dp) is_positive_z=.false.
        i=i+1
    enddo

end function rotationXaxis


! returns Rmat=Rz(alpha)Rx(beta)Rz(gamma)
! exterior rotation over axis z-x-z over (gamma, beta, alpha) 
! equivalently to interior Euler rotatation over z-x'-z'' by (alpha, beta, gamma)

function rotation_matrix(alpha,beta,gamma)result(Rmat)

    real(dp), intent(in)    :: alpha,beta,gamma
    real(dp) :: Rmat(3,3)

    ! local variables
    real(dp) :: cal, sal, cbe, sbe, cga, sga

    cal=cos(alpha)
    sal=sin(alpha)
    cbe=cos(beta)    
    sbe=sin(beta)
    cga=cos(gamma)
    sga=sin(gamma)

    ! Rmat=Rz(alpha)Rx(beta)Rz(gamma)


    Rmat(1,1)= cal*cga-sal*sga*cbe
    Rmat(2,1)= sal*cga+cal*sga*cbe
    Rmat(3,1)= sbe*sga

    Rmat(1,2)= -cal*sga-sal*cbe*cga
    Rmat(2,2)= cal*cbe*cga -sal*sga
    Rmat(3,2)= sbe*cga

    Rmat(1,3)= sal*sbe
    Rmat(2,3)=-cal*sbe
    Rmat(3,3)= cbe


end function rotation_matrix

! Pre: chain(s1,:)={0,0,0)
! Post: Rotate nuclesome conformation chain onto chain_rot 
! such that the three segments (sgraftpts=(s1,s2,s3)
! are located at postion chain_rot(:,s1)=(0,0,0) and 
! chains_rot(:,s2) =(0,0,zval) ( on z-axis) and chain_rot(s3)=(0,yval,uval) ( in z-y  plane) 


subroutine rotate_nucl_chain(chain,chain_rot,sgraftpts,nseg)

    use vectornorm

    real(dp), intent(in)    :: chain(:,:)
    real(dp), intent(inout) :: chain_rot(:,:)  
    integer , intent(in)    :: sgraftpts(3) 
    integer , intent(in)    :: nseg   

    ! .. local variable
    real(dp) :: Rmat(3,3) ,vec2(3), vec3(3),vec3prime(3)
    real(dp) :: alpha, beta, gamma
    integer :: s

    ! .. determine rotation angles 
    
    vec2(:) = chain(:,sgraftpts(2))
    vec3(:) = chain(:,sgraftpts(3))

    ! atan2(Y, X) angle of X+iY between  -pi< angle < pi

    gamma = atan2(vec2(1),vec2(2)) 
    beta  = atan2(sqrt(vec2(1)**2+vec2(2)**2),vec2(3))

    Rmat = rotation_matrix(0.0_dp,beta,gamma) 

    vec3prime = Matmul(Rmat,vec3) ! new coordinate third vector

    alpha=atan2(vec3prime(1),vec3prime(2)) ! rotation such that new third vector is in z-y plane
    
    ! .. make rotation matrix

    Rmat = rotation_matrix(alpha,beta,gamma)

    ! .. apply Rmat onto chain

    do s=1,nseg
        chain_rot(:,s)=Matmul(Rmat,chain(:,s))
    enddo  


end subroutine rotate_nucl_chain

subroutine unit_test_rotation_matrix(info)

    integer, intent(out) :: info

    real(dp) :: Rmat(3,3),Rmattest(3,3),Rdiffsqr(3,3)
    real(dp) :: alpha, beta, gamma
    integer :: i,j
    real(dp) :: chain(3,4),chain_rot(3,4)
    real(dp) :: x(3), Rdiffmax
    real(dp), parameter :: epsRdiffmax= 0.00001_dp


    !  Euler matrix from Mathematica for alpha,beta, gamma
    alpha=1.20_dp
    beta=0.714_dp
    gamma=0.5_dp

    Rmat=transpose(reshape((/-0.0197023_dp, -0.791881_dp, 0.610357_dp,&
        0.949233_dp, -0.206516_dp, -0.237294_dp,&
        0.313958_dp, 0.574695_dp, 0.755748_dp /),shape(Rmat)))

    ! Rotatation
    Rmattest=rotation_matrix(alpha,beta,gamma)
    
    ! difference
    Rdiffsqr=(Rmattest-Rmat)**2
    Rdiffmax=sqrt(sum(Rdiffsqr))


    if(Rdiffmax<epsRdiffmax) then
        info=0
    else
        print*,"Test rotation matrix failed: epsRdiffmax> Rdiffmax: Rdiffmax=",Rdiffmax
        info=1
    endif       
   
end subroutine unit_test_rotation_matrix
 
subroutine unit_test_rotate_nucl_chain(info)

    use vectornorm

    integer, intent(out) :: info

    ! .. local variables

    real(dp) :: chain(3,5),chain_rot(3,5)
    integer :: nseg,sgraftpts(3),i
    real(dp) :: radius
    real(dp), parameter :: epsradius = 0.000001_dp

    ! example of a 'chain conformation'
    chain=reshape((/0.0_dp,0.0_dp,0.0_dp,&
        -1.0_dp, 0.1230_dp,-1.0_dp,&
        -0.0197023_dp, -0.791881_dp, 0.610357_dp,&
        0.949233_dp, -0.206516_dp, -0.237294_dp,&
        0.313958_dp, 0.574695_dp, 0.755748_dp /),shape(chain))

    nseg=5
    sgraftpts=(/1,2,3/)

    call rotate_nucl_chain(chain,chain_rot,sgraftpts,nseg)

        
    ! test
    info=0

    radius=L2norm(chain_rot(:,1),3)
    if(radius> epsradius) then
        print*,"First chain conformation not in origin" 
        info=1
    endif  
    
    radius=sqrt(chain_rot(1,2)**2+chain_rot(2,2)**2)
    if( radius>epsradius) then 
        print*,"Second chain conformation not on z-axis" 
        info=1
    endif

    radius=sqrt(chain_rot(1,3)**2)
    if(radius> epsradius) then
        print*,"Third chain conformation not in z-y plane" 
        info=1
    endif
           
end subroutine unit_test_rotate_nucl_chain

! identical to rotate_nucl_chain plus test of rotated conformation

subroutine rotate_nucl_chain_test(chain,chain_rot,sgraftpts,nseg,write_rotations)

    use vectornorm

    real(dp), intent(in)    :: chain(:,:)
    real(dp), intent(inout) :: chain_rot(:,:)  
    integer , intent(in)    :: sgraftpts(3) 
    integer , intent(in)    :: nseg  
    logical , intent(in)    :: write_rotations
    
    ! .. local variable
    real(dp) :: Rmat(3,3) ,vec2(3), vec3(3),vec3prime(3)
    real(dp) :: alpha, beta, gamma
    integer :: s,i

    integer   :: info 
    real(dp) :: radius
    real(dp), parameter :: epsradius = 0.000001_dp

    ! .. determine rotation angles 
    
    vec2(:) = chain(:,sgraftpts(2))
    vec3(:) = chain(:,sgraftpts(3))

    ! atan2(Y, X) angle of X+iY between  -pi< angle < pi

    gamma = atan2(vec2(1),vec2(2)) 
    beta  = atan2(sqrt(vec2(1)**2+vec2(2)**2),vec2(3))

    Rmat = rotation_matrix(0.0_dp,beta,gamma) 

    vec3prime = Matmul(Rmat,vec3) ! new coordinate third vector

    alpha=atan2(vec3prime(1),vec3prime(2)) ! rotation such that new third vector is in z-y plane
    
    ! .. make rotation matrix

    Rmat = rotation_matrix(alpha,beta,gamma)

    ! .. apply Rmat onto chain

    do s=1,nseg
        chain_rot(:,s)=Matmul(Rmat,chain(:,s))
    enddo   


    ! test
    info=0

    radius=L2norm(chain_rot(:,sgraftpts(1)),3)
    if(radius> epsradius) then
        print*,"First chain conformation not in origin" 
        info=1
    endif  
    
    radius=sqrt(chain_rot(1,sgraftpts(1))**2+chain_rot(2,sgraftpts(2))**2)
    if( radius>epsradius) then 
        print*,"Second chain conformation not on z-axis" 
        info=2
    endif

    radius=sqrt(chain_rot(1,sgraftpts(3))**2)
    if(radius> epsradius) then
        print*,"Third chain conformation not in z-y plane" 
        info=3
    endif
    
    if(info/=0.or.write_rotations) then
        print*,"Write Euler rotations:"
        print*,"alpha=",alpha
        print*,"beta =",beta
        print*,"gamma=",gamma
        do s=1,3
            do i=1,3 
                print*,s,i,chain(i,sgraftpts(s)),chain_rot(i,sgraftpts(s))
            enddo  
            print*,""
        enddo    
    endif
             
end subroutine rotate_nucl_chain_test

! Calculation the orientation vector: 
! Orientation vector defined as  the normal vector to the plane spanned three points
! input: real(dp)  chain(3,nseg)           : coordiantes of all nseg segments 
!        integer   orient_triplet(nnucl,3) : segment id (s) of triplet of point ( plane ) 1=origin=CA
! output: real(dp) orient_vectors(3,nnucl) : normal vector; first index cartesian coordinate second index nucleosome  

subroutine orientation_vector(chain,orient_triplets,orient_vectors)

    use globals, only    : nnucl, nseg
    use quaternions, only : crossproduct

    real(dp), intent(in) :: chain(3,nseg)  
    integer , intent(in) :: orient_triplets(nnucl,3) 
    real(dp), intent(inout) :: orient_vectors(3,nnucl)

    integer :: triplet(3), n
    real(dp) :: u1(3),u2(3),u3(3)

    do n=1,nnucl     
        triplet=orient_triplets(n,1:3) 
        u1(1:3)=chain(1:3,triplet(1))
        u2(1:3)=chain(1:3,triplet(2))
        u3(1:3)=chain(1:3,triplet(3))
        orient_vectors(:,n)=crossproduct(u2-u1, u3-u1)
    enddo

end subroutine orientation_vector


! Calculation the orientation vector of reference histone  
! Orientation vector defined as the normal vector to the plane spanned three points
! input: real(dp)  chain_elem(3,nseg)%elem()  : coordiantes of AA chain elems segments need on CA=> elem(1) 
!        integer   orient_triplet_ref(3)      : AA segment id (sAA) of triplet of point ( plane ) 1=origin=CA 
! output: real(dp) orient_vector(3)           : normal vector  

subroutine orientation_vector_ref(chain_elem,orient_triplet_ref,orient_vector)

    use quaternions, only : crossproduct
    use chains, only : var_darray

    type(var_darray), dimension(:,:), intent(in) :: chain_elem
    integer , intent(in) :: orient_triplet_ref(3) 
    real(dp), intent(inout) :: orient_vector(3)

    real(dp) :: u1(3),u2(3),u3(3)
    integer :: k

    do k=1,3
        u1(k)=chain_elem(k,orient_triplet_ref(1))%elem(1)
        u2(k)=chain_elem(k,orient_triplet_ref(2))%elem(1)
        u3(k)=chain_elem(k,orient_triplet_ref(3))%elem(1)
    enddo    
    orient_vector=crossproduct(u2-u1, u3-u1)

end subroutine orientation_vector_ref

! Get coordinates of for all orientation triplets of all nucleosomes,
! All triplets coordinate translated such that first coordiant ( CA) over every triplet is (0,0,0)
! input:  real(dp)  chain_rot(3,nseg)       : coordiantes of all nseg segments 
!         integer   orient_triplet(nnucl,3) : segment id of triplet of point ( plane ) 1=origin=CA
! output: real(dp)  orient_vectors(3,nnucl) : normal vector; first index cartesian coordinate second index nucleosome

subroutine orientation_coordinates(chain_rot,orientation_triplets,triangle_orient_vectors)
    
    use globals, only : nnucl

    real(dp), dimension(:,:),  intent(in)    :: chain_rot               ! (k,s) 
    integer, dimension(:,:) ,  intent(in)    :: orientation_triplets    ! (nnucl,3)
    real(dp),dimension(:,:,:), intent(inout) :: triangle_orient_vectors ! (3,3,nnucl)

    ! local arguments

    integer :: n, tri, s, s0

    do n=1,nnucl 
        do tri=1,3
            s  = orientation_triplets(n,tri)
            s0 = orientation_triplets(n,1) 
            triangle_orient_vectors(:,tri,n) = chain_rot(:,s)-chain_rot(:,s0)
            !print*,"triangle ="
            !print*,triangle_orient_vectors(:,tri,n)  
        enddo
    enddo

end subroutine 



! Make rotation axis/angle for all nnucl nucleosome
! Apply rotation for all nucleosome to chain_elem retruned in chain_elem_rot


subroutine rotate_chain_elem(orient_vector_ref,orient_vectors,nelemAA,chain_elem,chain_elem_rot,&
    orientation_triplets,triangle_orient_vectors,segnumAAstart)
    
    use mathconst
    use globals, only     : nnucl, nseg, nsegAA, DEBUG
    use quaternions, only : rot_axis_angle, rot_axis_angle_to_quat, rotation_matrix_from_quat, vec_norm 
    use quaternions, only : print_rotation_matrix
    use chains, only      : var_darray

    real(dp), intent(in) :: orient_vector_ref(3) 
    real(dp), intent(in) :: orient_vectors(3,nnucl)
    integer, intent(in)  :: nelemAA(nsegAA)
    type(var_darray), dimension(:,:), intent(in) :: chain_elem
    type(var_darray), dimension(:,:,:), intent(inout) :: chain_elem_rot 
    integer, dimension(:,:) ,  intent(in) :: orientation_triplets       ! (nnucl,3)
    real(dp),dimension(:,:,:), intent(in) :: triangle_orient_vectors    ! (3,3,nnucl)
    integer, dimension(:), intent(in) :: segnumAAstart

    ! local arguments 

    integer  :: sAA, j, n, k
    real(dp) :: a(3), b(3), u(3), unorm, angle, a_rot(3)
    real(dp) :: q(4)
    real(dp) :: Rmat1(3,3),Rmat2(3,3),Rmat_comb(3,3)
    real(dp) :: xvec(3), wvec(3)

    do n=1,nnucl                      ! loop over nucleosomes 

        ! rotation of orientation vector  

        a = orient_vector_ref         ! vector to rotate
        b = orient_vectors(:,n)       ! target vector

        call rot_axis_angle(a, b, u, angle)
        unorm = vec_norm(u)
        u(1:3) = u(1:3)/unorm         ! rotation axis 

        if(abs(angle)<epsAngle) u=a   ! angle ==0 no rotation u=>a  

        q = rot_axis_angle_to_quat(u, angle)   ! rotation quaternion
        call rotation_matrix_from_quat(q,Rmat1)  

        ! apply second in plane rotation to get coordinate to coincede

        sAA=orientation_triplets(n,2) -segnumAAstart(n)+1 ! second element of triplet 
        ! print*,"rotate_chain_elem: sAA=",sAA," orientation_triplets(n,2)=",orientation_triplets(n,2)
        do k=1,3 
            a(k)=chain_elem(k,sAA)%elem(1) ! vector to rotate
        enddo 

        a_rot = matmul(Rmat1,a) ! apply rotation
        b = triangle_orient_vectors(:,2,n)
        
        call rot_axis_angle(a_rot, b, u, angle)
        unorm = vec_norm(u)
        u(1:3) = u(1:3)/unorm         ! rotation axis 

        if(abs(angle)<epsAngle) u=a   ! angle ==0 no rotation u=>a  

        q = rot_axis_angle_to_quat(u, angle)   ! rotation quaternion
        call rotation_matrix_from_quat(q,Rmat2)  

        Rmat_comb=matmul(Rmat2,Rmat1) ! Rmat_comb= Rmat2 x Rmat1
     

        ! apply rotation to chain_elem

        do sAA=1, nsegAA
            do j = 1, nelemAA(sAA)  
                do k=1,3 
                    xvec(k)=chain_elem(k,sAA)%elem(j)
                enddo       
        
                wvec=matmul(Rmat_comb,xvec)
        
                do k=1,3
                    chain_elem_rot(k,sAA,n)%elem(j)=wvec(k)  ! assign rotation
                enddo    
            enddo 
        enddo   

        if(DEBUG) then 
            print*,"Module : chain_rotation : rotate_chain_elem "
            print*,"print chain_elem(k,s) and rotated chain_elem_rot(k,s,n)"

            do sAA=1,nsegAA
                do j=1,nelemAA(sAA)
                    do k=1,3
                        print*,"sAA=",sAA,"j=",j,"k=",k," ",chain_elem(k,sAA)%elem(j),chain_elem_rot(k,sAA,n)%elem(j)
                    enddo 
                    print*,""    
                enddo
            enddo  
        endif    
        
    enddo  

end subroutine rotate_chain_elem 


subroutine unit_test_rotation(info)

    use quaternions, only : rot_axis_angle, rot_axis_angle_to_quat, rotation_matrix_from_quat, vec_norm
    use quaternions, only : print_rotation_matrix

    integer, intent(out) :: info

    ! .. local variables

    real(dp), parameter :: eps = 0.000001_dp

    real(dp) :: a(3), b(3), anorm, bnorm, u(3), unorm, angle, a_rot(3)
    real(dp) :: q(4)
    real(dp) :: Rmat(3,3)
    real(dp) :: wvec(3)
    integer :: i

    info=0

    a = (/1.0_dp,1.0_dp,0.0_dp/)
    anorm= vec_norm(a)
    a(1:3)=a(1:3)/anorm
    b = (/1.0_dp,-1.0_dp,0.0_dp/)   ! target vector
    bnorm = vec_norm(b)
    b(1:3)=b(1:3)/bnorm

    call rot_axis_angle(a, b, u, angle)
    unorm = vec_norm(u)
    u(1:3) = u(1:3)/unorm         ! rotation axis 

    if(abs(angle)<epsAngle) u=a   ! angle ==0 no rotation u=>a  

    print*,"a=",a
    print*,"b=",b
    print*,"angle=",angle
    print*,"nvec=",u

    q = rot_axis_angle_to_quat(u, angle)   ! rotation quaternion
    call rotation_matrix_from_quat(q,Rmat)
    a_rot = matmul(Rmat,a) ! apply rotation
    wvec = a_rot-b
    print*,"arot-b =",wvec

    do i=1,3
        if(abs( wvec(i))>=eps) info=1
    enddo


end subroutine unit_test_rotation


! Check if chain_index(k,s)%elem(1) == nseg(k,s) 
! input :  integer :: nseg       
!          integer :: nelem(nseg)
!          type(var_darray):: chain_elem_index
!          real(dp) :: chain(3,nseg)
! output : logical ::isSame  : isSame =.true chain_index(k,s)%elem(1) == nseg(k,s) .false. otherwise

subroutine check_chain_elem_index_and_chain(nseg,nelem,chain_elem_index,chain,isSame)
     
    use chains, only : var_darray

    integer, intent(in) :: nseg       
    integer, dimension(:), intent(in) :: nelem
    type(var_darray), dimension(:,:), allocatable, intent(in) ::chain_elem_index
    real(dp),dimension(:,:), intent(in)    :: chain
    logical, intent(inout) :: isSame

    integer :: s, i
    logical :: flag
    real(dp) :: sqrdiff 
    real(dp), parameter  :: epssqrdiff = 0.00000001_dp 

    s=1              ! loop over segments
    flag=.true.

    do while(s<=nseg.and.flag)    
       
        sqrdiff=0.0_dp
        do i=1,3
            sqrdiff=sqrdiff+(chain_elem_index(i,s)%elem(1)-chain(i,s))**2
        enddo
        if(sqrdiff>=epssqrdiff) then 
            flag=.false.
            ! print*,"sqrdiff=",sqrdiff,"s= ",s
        endif    
        s=s+1    
    enddo    

    isSame=flag

end subroutine check_chain_elem_index_and_chain


! Rotates chain(3,n) into  chain(3,nelem) for rotation that roates vector a into vector b 
! input: real(dp)  chain(3,n) : cartesian coordiantes chain of  n elememnts  
!        real(dp )  a         : vector to rotate 
!        real(dp )  b         : target vector 
! output: real(dp) chain_rot(3,n)  : rotated chain

subroutine rotate_chain(chain,chain_rot,a,b)
    
    use globals, only     : DEBUG
    use quaternions, only : rot_axis_angle, rot_axis_angle_to_quat, rotation_matrix_from_quat, vec_norm 

    real(dp), intent(in)    :: chain(:,:)
    real(dp), intent(inout) :: chain_rot(:,:)
    real(dp), intent(in)    :: a(3) ! vector to rotate
    real(dp), intent(in)    :: b(3) ! target vector

    ! local arguments 

    integer  :: s, k, ndim
    real(dp) :: u(3), unorm, angle
    real(dp) :: q(4)
    real(dp) :: Rmat(3,3)
    real(dp) :: xvec(3), wvec(3)

    ndim = size(chain,2)

    ! rotation of orientation vector  

    call rot_axis_angle(a, b, u, angle)
    unorm = vec_norm(u)
    u(1:3) = u(1:3)/unorm         ! rotation axis 

    if(abs(angle)<epsAngle) u=a   ! angle ==0 no rotation u=>a  

    q = rot_axis_angle_to_quat(u, angle)   ! rotation quaternion
    call rotation_matrix_from_quat(q,Rmat)  

    ! apply rotation to chain

    do s=1,ndim 
        xvec =chain(:,s)       
        wvec = matmul(Rmat,xvec)
        chain_rot(:,s) = wvec  ! assign rotation
    enddo 
    
    if(DEBUG) then 
        print*,"Module : chain_rotation : rotate_chain"
        do s=1,ndim
            print*,"s=",s," ",(chain(k,s),k=1,3)
        enddo    
    endif

end subroutine rotate_chain
    


end module chain_rotation
    
