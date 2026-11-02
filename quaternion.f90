module quaternions

    use precision_definition
  
    implicit none

contains  

!  Computes conjugates of quaternion q = a + bi + cj + dk 
!  The conjugate of q is q* =  = a - bi - cj - dk
!  input real(dp) q(4) 
!  output real(dp) q_conj(4)

    function  quat_conjugate(q) result(q_conj)

        real(dp), intent(in) :: q(4)
        real(dp) :: q_conj(4)

        q_conj(1)   =  q(1)
        q_conj(2:4) = -q(2:4)

    end function quat_conjugate
   
!   Computes exponential of quaternion q = a + bi + cj + dk = a + v ; v=bi +cj  + dk 
!   Exponential is given by 
!      e^q = e^a * ( cos ( ||v|| ) + v/||v|| sin ||v|| )
!      with v = sqrt ( b^2 + c^2 + d^2 )

    function  quat_exponential(q) result(exp_q)
        
        real(dp), intent(in) :: q(4)
        real(dp) :: exp_q(4)

        real(dp) :: v(3),vnorm ,a

        a = q(1)
        v = q(2:4)
        vnorm = sqrt(sum( q(2:4) ** 2 ))

        exp_q(1) = exp(a)*cos(vnorm)
        
        if(vnorm /= 0.0_dp) then
            exp_q(2:4) = exp(a)*sin(vnorm) * v / vnorm
        else
            exp_q(2:4) = 0.0_dp
        endif

    end function quat_exponential

!   Computes logarithm of quaternion q = a + bi + cj + dk = a + v ; v=bi +cj  + dk 
!   Logarithm is given by 
!      ln(q) = ln(||q||) + v/||v|| arccos(a/||q||)
!      with v = sqrt ( b^2 + c^2 + d^2 )

    function  quat_logarithm(q) result(ln_q)
        
        real(dp), intent(in) :: q(4)
        real(dp) :: ln_q(4)

        real(dp) :: v(3),vnorm,a,qnorm

        a = q(1)
        v = q(2:4)
        vnorm = sqrt(sum(q(2:4)**2))
        qnorm = sqrt(sum(q(1:4)**2))

        ln_q(1) = log(qnorm)
        
        if(vnorm /= 0.0_dp) then
            ln_q(2:4) = acos(a/qnorm) * v / vnorm
        else
            ln_q(2:4) = 0.0_dp
        endif

    end function quat_logarithm




!   Computes the inverse of quaternion q is q = a + bi + cj + dk.
!      q_inv= inverse( q ) = conjugate ( q ) / ( norm ( q) )^2.
!   input : real(dp) q(4), output real(dp) q(4)
!   assumes quaternion has non-zero norm 

    function  quat_inverse(q) result(q_inv)

        real(dp), intent(in) :: q(4)
        real(dp) :: q_inv(4)
        real(dp) :: qnorm 

        qnorm=sum(q(1:4) ** 2)  ! norm

        q_inv(1:4) = q(1:4) / qnorm  ! norm
        q_inv(2:4) = - q_inv(2:4)    ! conjugation 

    end function quat_inverse


!   Computes the Hamilton product of two quaternions p q i.e., the multiplciation p q 
!   q_prod= p q = (a1 + b1i + c1j + d1k)(a2 + b2i + c2j + d2k) 
!     with        i j = -j i = k
!       j k = -k j = i
!       k i = -i k = j
!       i i =  j j = k * k = -1
!   input : real(dp) p(4), q(4) output real(dp) q_prod(4)

    function quat_mult( p, q )result(q_prod)

        real(dp), intent(in) :: p(4), q(4)
        real(dp) :: q_prod(4)

        q_prod(1) = p(1) * q(1) - p(2) * q(2) - p(3) * q(3) - p(4) * q(4)
        q_prod(2) = p(1) * q(2) + p(2) * q(1) + p(3) * q(4) - p(4) * q(3)
        q_prod(3) = p(1) * q(3) - p(2) * q(4) + p(3) * q(1) + p(4) * q(2)
        q_prod(4) = p(1) * q(4) + p(2) * q(3) - p(3) * q(2) + p(4) * q(1)

    end function quat_mult

!  Computes norm of a quaternion q = a + bi + cj + dk.
!  norm = sqrt( a^2 +b^2+c^2+d ^2)

    function quat_norm(q)result(qnorm)

        real(dp), intent(in) :: q(4)
        real(dp) :: qnorm

        qnorm = sqrt(sum(q(1:4)**2))

    end function quat_norm

!  Computes 3D norm of a vector u=bi + cj + dk.
!  norm = sqrt( b^2+c^2+d ^2)

    function vec_norm(u)result(unorm)

        real(dp), intent(in) :: u(3)
        real(dp) :: unorm

        unorm = sqrt(sum(u(1:3)**2))

    end function vec_norm

!  Computes a quaternion q = a + bi + cj + dk
!  associated with rotation axis u  and rotation angle  theta
!  q = cos(theta/2) + (u1 i + u2 j + u3 k ) sin(theta/2)
!  input  real(dp) u(3) and real(dp) theta 
!  vector u assumed to have unit length||u||=1) and theta in radians

    function rot_axis_angle_to_quat(u, theta) result(q)

        real(dp), intent(in) :: u(3), theta
        real(dp) :: q(4)

        q(1)   = cos( 0.5_dp * theta )
        q(2:4) = sin( 0.5_dp * theta ) * u(1:3) 

    end function rot_axis_angle_to_quat

!  Computes 3D rotatation matrix associated with quaternion q = s + xi + yj + zk
!  with quaternion multiplication p'=q p q^-1 = R p   and p = xp i +yp j+zp k scalar(p) = sp =0  
!  assumes ||q||=1 

    subroutine rotation_matrix_from_quat(q,R)
 
        real(dp), intent(in) :: q(4)
        real(dp), intent(inout) :: R(3,3)

        real(dp) :: s,x,y,z

        s = q(1)
        x = q(2)
        y = q(3)
        z = q(4)

        R(1,1) = 1.0_dp - 2.0_dp * ( y * y + z * z)
        R(1,2) = 2.0_dp * ( x * y - s * z )
        R(1,3) = 2.0_dp * ( x * z + s * y )
        R(2,1) = 2.0_dp * ( x * y + s * z )
        R(2,2) = 1.0_dp - 2.0_dp * ( x * x + z * z )
        R(2,3) = 2.0_dp * ( y * z - s * x )
        R(3,1) = 2.0_dp * ( x * z - s * y )
        R(3,2) = 2.0_dp * ( y * z + s * x )
        R(3,3) = 1.0_dp - 2.0_dp * ( x * x + y * y )

    end subroutine rotation_matrix_from_quat


    
    !  Computes rotation axis and angle between vector u and v
    !  tan(theta)=sin(theta)/cos(theta)= ||u x v||/ (u,v)
    !  Here atan2(y, x) computes the principal value of the argument
    !  function of the complex number x + i y

    subroutine rot_axis_angle(u, v, rotaxis, theta)

        real(dp), dimension(3), intent(in) :: u, v 
        real(dp), dimension(3),intent(inout) :: rotaxis 
        real(dp), intent(inout) :: theta

        real(dp) :: norm_rotaxis, dot

        rotaxis      = crossproduct(u, v)
        dot          = dotproduct(u, v)
        norm_rotaxis = vec_norm(rotaxis)
        
        theta  = atan2(norm_rotaxis, dot)

        ! rotaxis  = rotaxis  /norm_rotaxis

    end subroutine rot_axis_angle


    ! Inner product of vectors a and b with dimenstion 3
    ! similar to intrisic function dotproduct
    ! dotprod = a.b

    function dotproduct(a,b) result(dotprod)
       
        real(dp), dimension(3) :: a, b 
        real(dp) :: dotprod
      
        dotprod = sum(a*b)

    end function dotproduct

    ! Cross product of vectors a and b with dimenstion 3
    ! crossprod = a x b 
    
    function crossproduct(a,b)result(crossprod)
      
      real(dp), intent(in), dimension(3) :: a, b
      real(dp), dimension(3) :: crossprod
      
      crossprod(1) = a(2)*b(3) - a(3)*b(2)
      crossprod(2) =-a(1)*b(3) + a(3)*b(1)
      crossprod(3) = a(1)*b(2) - a(2)*b(1)

    end function crossproduct

    ! print rotation matrix RMat 3x3
     
    subroutine print_rotation_matrix(Rmat)

        real(dp), intent(in) :: Rmat(3,3)
          
        integer :: i,j 

        do i=1,3
            print*,(Rmat(i,j),j=1,3)
        enddo  
    end subroutine print_rotation_matrix


end module quaternions

