!Common slow-roll functions for the spontaneous symmetry breaking Inflation 1,2,3,4,5,6 potential
!
!
!V(phi) = M^4 [ 1 + alpha x^2 + beta x^4 ]
!
!1: alpha>0, beta>0
!2: alpha<0, beta<0
!3: alpha>0, beta<0, inflation proceeds from the right to the left
!4: alpha>0, beta<0, inflation proceeds from the left to the right
!5: alpha<0, beta>0, inflation proceeds from the left to the right
!6: alpha<0, beta>0, inflation proceeds from the right to the left
!
!x = phi/Mp


module ssbicommon
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public ssbi_norm_potential
  public ssbi_norm_deriv_potential, ssbi_norm_deriv_second_potential
  public ssbi_epsilon_one, ssbi_epsilon_two, ssbi_epsilon_three
  public ssbi_efold_primitive, find_ssbi_x_trajectory
  public ssbi136_x_epstwozero, ssbi245_x_potzero, ssbi3456_x_derivpotzero

contains


!returns V/M^4
  function ssbi_norm_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssbi_norm_potential = 1._kp+alpha*x**2+beta*x**4

  end function ssbi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function ssbi_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta

    ssbi_norm_deriv_potential = 2._kp*x*(alpha+2._kp*beta*x**2)

  end function ssbi_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M^4
  function ssbi_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta

   
    ssbi_norm_deriv_second_potential = 2._kp*(alpha+6._kp*beta*x**2)
    

  end function ssbi_norm_deriv_second_potential



!epsilon_one(x)
  function ssbi_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi_epsilon_one
    real(kp), intent(in) :: x,alpha,beta


    ssbi_epsilon_one = 2._kp*(alpha*x+2._kp*beta*x**3)**2/(1._kp+alpha*x**2+beta*x**4)**2
    
  end function ssbi_epsilon_one


!epsilon_two(x)
  function ssbi_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    
    ssbi_epsilon_two = 4._kp*(-alpha+(alpha**2-6._kp*beta)*x**2+alpha*beta*x**4+ & 
         2._kp*beta**2*x**6)/(1._kp+alpha*x**2+beta*x**4)**2
    
  end function ssbi_epsilon_two


!epsilon_three(x)
  function ssbi_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi_epsilon_three
    real(kp), intent(in) :: x,alpha,beta

    
    ssbi_epsilon_three = (4._kp*x**2*(alpha+2._kp*beta*x**2)*(-3._kp*alpha**2+6._kp*beta+ & 
         alpha*(alpha**2-12._kp*beta)*x**2+3._kp*(alpha**2-8._kp*beta)*beta* &
         x**4+2._kp*beta**3*x**8))/((1._kp+alpha*x**2+beta*x**4)**2* &
         (-alpha+(alpha**2-6._kp*beta)*x**2+alpha*beta*x**4+2._kp*beta**2*x**6))
    
  end function ssbi_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function ssbi_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssbi_efold_primitive


    if (alpha*beta .eq. 0._kp) stop 'ssbi_efold_primitive: alpha*beta=0!'

    ssbi_efold_primitive = 1._kp/(2._kp*alpha)*log(x)+x**2/8._kp &
         +(alpha**2-4._kp*beta)/(16._kp*alpha*beta) &
         *log(abs(1._kp+2._kp*beta/alpha*x**2)) 

  end function ssbi_efold_primitive



  function find_ssbi_x_trajectory(x,ssbiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssbiData
    real(kp) :: find_ssbi_x_trajectory
    real(kp) :: alpha,beta,NplusNuend

    alpha = ssbiData%real1
    beta = ssbiData%real2
    NplusNuend = ssbiData%real3

    find_ssbi_x_trajectory = ssbi_efold_primitive(x,alpha,beta) - NplusNuend
   
  end function find_ssbi_x_trajectory

! Return the position x at which epsilon2 vanishes for ssbi1, ssbi3, ssbi6
  function ssbi136_x_epstwozero(alpha,beta)
  implicit none
    real(kp) :: ssbi136_x_epstwozero
    real(kp), intent(in) :: alpha,beta

    ssbi136_x_epstwozero = sqrt(-(alpha/(6._kp*beta))+(16._kp*alpha**3*beta**3+2._kp* &
         sqrt(cmplx((64._kp*alpha**6+(5._kp*alpha**2-36._kp*beta)**3)* &
         beta**6,0._kp)))**(1._kp/3._kp)/(6._kp*2._kp**(1._kp/3._kp)*beta**2)+ &
         (-5._kp*alpha**2+36._kp*beta)/(6._kp*(8._kp*alpha**3*beta**3+ & 
         3._kp*sqrt(3._kp)*sqrt(cmplx((alpha**2-4._kp*beta)*beta**6* &
         (7._kp*alpha**4-72._kp*alpha**2*beta+432._kp*beta**2),0._kp)))**(1._kp/3._kp)))

   end function ssbi136_x_epstwozero



! Return the position x at which the potential vanishes for ssbi2, ssbi4 and ssbi5
  function ssbi245_x_potzero(alpha,beta)
  implicit none
    real(kp) :: ssbi245_x_potzero
    real(kp), intent(in) :: alpha,beta

    ssbi245_x_potzero = sqrt(-(alpha+sqrt(alpha**2-4._kp*beta))/(2._kp*beta))

   end function ssbi245_x_potzero


! Return the position x at which the first derivative of the potential with respect to x vanishes for ssbi3&4 and ssbi5&6
  function ssbi3456_x_derivpotzero(alpha,beta)
  implicit none
    real(kp) :: ssbi3456_x_derivpotzero
    real(kp), intent(in) :: alpha,beta

    ssbi3456_x_derivpotzero = sqrt(-alpha/(2._kp*beta))

   end function ssbi3456_x_derivpotzero


end module ssbicommon
