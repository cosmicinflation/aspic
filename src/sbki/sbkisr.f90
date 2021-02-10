!slow-roll functions for the symmetry breaking Kahler inflation potential
!
!V(phi) = M**4 * x**2 exp(alpha x**2 + alpha**2/6 x**4)
!
!x = phi/Mp
!
!alpha is a small parameter that can be negative or positive

module sbkisr
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use specialinf, only : log_otherbranchcut
  implicit none

  private

  public  sbki_norm_potential, sbki_epsilon_one, sbki_epsilon_two, sbki_epsilon_three
  public  sbki_x_endinf, sbki_efold_primitive, sbki_x_trajectory
  public  sbki_epsilon_one_min, sbki_efoldmax
  public  sbki_norm_deriv_potential, sbki_norm_deriv_second_potential
  public  sbki_xinimax, sbki_check_params, sbki_alphamin, sbki_alphamax

contains

! return true if there is more than efoldnum of inflation
  function sbki_check_params(efoldNum, alpha)
    implicit none
    logical :: sbki_check_params
    real(kp), intent(in) :: efoldNum, alpha    
    
    sbki_check_params = sbki_efoldmax(alpha) .ge. efoldNum

  end function sbki_check_params


!ensures that epsonemin < 1 for alpha > 0
  function sbki_alphamax()
    implicit none
    real(kp) :: sbki_alphamax
    sbki_alphamax = 9._kp/4._kp/(11._kp + 5._kp*sqrt(5._kp))
    
  end function sbki_alphamax


!ensures that epsonemin < 1 for alpha < 0
  function sbki_alphamin()
    real(kp) :: sbki_alphamin

    sbki_alphamin = 9._kp/4._kp/(11._kp - 5._kp*sqrt(5._kp))

  end function sbki_alphamin


!returns V/M**4
  function sbki_norm_potential(x, alpha)
    implicit none
    real(kp) :: sbki_norm_potential
    real(kp), intent(in) :: x, alpha

    sbki_norm_potential = x**2*exp(alpha*x**2+alpha**2/6._kp*x**4)

  end function sbki_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function sbki_norm_deriv_potential(x, alpha)
    implicit none
    real(kp) :: sbki_norm_deriv_potential
    real(kp), intent(in) :: x, alpha

   sbki_norm_deriv_potential = (2._kp*exp((alpha*x**2*(6._kp+alpha*x**2))/6._kp)* &
                                x*(3._kp+3._kp*alpha*x**2+alpha**2*x**4))/3._kp

  end function sbki_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function sbki_norm_deriv_second_potential(x, alpha)
    implicit none
    real(kp) :: sbki_norm_deriv_second_potential
    real(kp), intent(in) :: x, alpha

    sbki_norm_deriv_second_potential = (2._kp*exp((alpha*x**2*(6._kp+alpha*x**2))/6.)* &
                                    (9._kp+45._kp*alpha*x**2+39._kp*alpha**2*x**4+12._kp* &
                                    alpha**3*x**6+2._kp*alpha**4*x**8))/9._kp

  end function sbki_norm_deriv_second_potential



!epsilon_one(x)
  function sbki_epsilon_one(x, alpha)    
    implicit none
    real(kp) :: sbki_epsilon_one
    real(kp), intent(in) :: x, alpha
    
    sbki_epsilon_one = (2._kp*(3._kp+alpha*x**2*(3._kp+alpha*x**2))**2)/(9._kp*x**2)
    
  end function sbki_epsilon_one


!epsilon_two(x)
  function sbki_epsilon_two(x, alpha)    
    implicit none
    real(kp) :: sbki_epsilon_two
    real(kp), intent(in) :: x, alpha
    
    sbki_epsilon_two = 4._kp/x**2-4._kp*alpha*(1._kp+alpha*x**2)
    
  end function sbki_epsilon_two


!epsilon_three(x)
  function sbki_epsilon_three(x, alpha)    
    implicit none
    real(kp) :: sbki_epsilon_three
    real(kp), intent(in) :: x, alpha
    
    sbki_epsilon_three = (-4._kp*(1._kp+alpha**2*x**4)*(3._kp+alpha*x**2* &
                        (3._kp+alpha*x**2)))/(3._kp*x**2*(-1._kp+alpha*x**2+alpha**2*x**4))
    
  end function sbki_epsilon_three

! This is the minimum value taken by epsilon_one
  function sbki_epsilon_one_min(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: sbki_epsilon_one_min

    if (alpha > 0._kp) then
        sbki_epsilon_one_min = 4._kp/9._kp*(11._kp+5._kp*sqrt(5._kp))*alpha
    else
        sbki_epsilon_one_min = 4._kp/9._kp*(11._kp-5._kp*sqrt(5._kp))*alpha
    endif

  end function sbki_epsilon_one_min



!returns the 2 solutions of eps1=1
  function sbki_x_epsoneunity(alpha)
    implicit none
    real(kp), dimension(2) :: sbki_x_epsoneunity
    real(kp), intent(in) :: alpha

    if (alpha.ge.sbki_alphamax()) then
       stop 'sbki_x_epsoneunity: alpha >= alphamax!'
    elseif (alpha.lt.sbki_alphamin()) then
       stop 'sbki_x_epsoneunity: alpha < alphamin!'
    endif
    
    sbki_x_epsoneunity(1) = sqrt((-4._kp*alpha**3+(10._kp*2**(2._kp/3._kp)*alpha**6)/ &
         (9._kp*alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp* &
         alpha+64._kp*alpha**2)))**(1._kp/3._kp)+2._kp**(1._kp/3._kp)*(9._kp* &
         alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp*alpha+ &
         64._kp*alpha**2)))**(1._kp/3._kp))/alpha**4)/(2._kp*sqrt(2._kp))- &
         1._kp/2._kp*sqrt(-(4._kp/alpha)-(5._kp*2**(2._kp/3._kp)*alpha**2)/ &
         (9._kp*alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp* &
         alpha+64._kp*alpha**2)))**(1._kp/3._kp)-(9._kp*alpha**8-44._kp* &
         alpha**9+sqrt(-alpha**16*(-81._kp+792._kp*alpha+64._kp*alpha**2)))** &
         (1._kp/3._kp)/(2._kp**(2._kp/3._kp)*alpha**4)+6._kp/(alpha**2* &
         sqrt((-4._kp*alpha**3+(10._kp*2._kp**(2._kp/3._kp)*alpha**6)/(9._kp* &
         alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp*alpha+64._kp* &
         alpha**2)))**(1._kp/3._kp)+2._kp**(1._kp/3._kp)*(9._kp* &
         alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp*alpha+64._kp* &
         alpha**2)))**(1._kp/3._kp))/alpha**4)))



    
    sbki_x_epsoneunity(2) = sqrt((-4._kp*alpha**3+(10._kp*2**(2._kp/3._kp)*alpha**6)/ &
         (9._kp*alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp* &
         alpha+64._kp*alpha**2)))**(1._kp/3._kp)+2._kp**(1._kp/3._kp)*(9._kp* &
         alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp*alpha+ &
         64._kp*alpha**2)))**(1._kp/3._kp))/alpha**4)/(2._kp*sqrt(2._kp))+ &
         1._kp/2._kp*sqrt(-(4._kp/alpha)-(5._kp*2._kp**(2._kp/3._kp)*alpha**2)/ &
         (9._kp*alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp* &
         alpha+64._kp*alpha**2)))**(1._kp/3._kp)-(9._kp*alpha**8-44._kp* &
         alpha**9+sqrt(-alpha**16*(-81._kp+792._kp*alpha+64._kp*alpha**2)))** &
         (1._kp/3._kp)/(2._kp**(2._kp/3._kp)*alpha**4)+6._kp/(alpha**2* &
         sqrt((-4._kp*alpha**3+(10._kp*2._kp**(2._kp/3._kp)*alpha**6)/(9._kp* &
         alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp*alpha+64._kp* &
         alpha**2)))**(1._kp/3._kp)+2._kp**(1._kp/3._kp)*(9._kp* &
         alpha**8-44._kp*alpha**9+sqrt(-alpha**16*(-81._kp+792._kp*alpha+64._kp* &
         alpha**2)))**(1._kp/3._kp))/alpha**4)))


  end function sbki_x_epsoneunity

  

!returns x at the end of inflation defined as epsilon1=1
  function sbki_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: sbki_x_endinf

    real(kp), dimension(2) :: xepsoneone

    xepsoneone = sbki_x_epsoneunity(alpha)

    sbki_x_endinf = xepsoneone(1)
   
  end function sbki_x_endinf


  
!returns the other value x at which epsilon1=1 which places an upper bound on xstar
  function sbki_xinimax(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: sbki_xinimax

    real(kp), dimension(2) :: xepsoneone

    xepsoneone = sbki_x_epsoneunity(alpha)
    
    sbki_xinimax = xepsoneone(2)

    
  end function sbki_xinimax



!this is integral(V(phi)/V'(phi) dphi)
  function sbki_efold_primitive(x, alpha)
    implicit none
    real(kp), intent(in) :: x, alpha
    real(kp) :: sbki_efold_primitive
    complex(kp) :: a, astar

    a = (-3._kp+sqrt(3._kp)*cmplx(0._kp,1._kp,kp))/(2._kp*alpha)
    astar = conjg(a)

    sbki_efold_primitive = -sqrt(3._kp)/(4._kp*alpha)*real(cmplx(0._kp,1._kp,kp)* &
                           log_otherbranchcut((x**2-a)/(x**2-astar)),kp)

  end function sbki_efold_primitive


!returns x at bfold=-efolds before the end of inflation
  function sbki_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha
    real(kp) :: sbki_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,xplus
    type(transfert) :: sbkiData

    ! Analytical Inversion
    complex(kp) :: a, astar, c
    a = (-3._kp+sqrt(3._kp)*cmplx(0._kp,1._kp,kp))/(2._kp*alpha)
    astar = conjg(a)
    c = (xend**2-astar)/(xend**2-a)*exp(4._kp*cmplx(0._kp,1._kp,kp)*alpha/sqrt(3._kp)*bfold)
    sbki_x_trajectory = real(sqrt((astar-c*a)/(1._kp-c)),kp)

  end function sbki_x_trajectory


! This is the maximum number of efolds one can realise
  function sbki_efoldmax(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: sbki_efoldmax, xend, xplus

    xend = sbki_x_endinf(alpha)
    xplus = sbki_xinimax(alpha)

    print *,'test',alpha,xend,xplus
    
    sbki_efoldmax = -sbki_efold_primitive(xend,alpha) + sbki_efold_primitive(xplus,alpha)

    
  end function sbki_efoldmax



  
end module sbkisr
