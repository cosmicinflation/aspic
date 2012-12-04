!Common slow-roll functions for the sneutrino supersymmetric 1,2,3,4,5,6 potential
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


module ssicommon
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public ssi_norm_potential
  public ssi_norm_deriv_potential, ssi_norm_deriv_second_potential
  public ssi_epsilon_one, ssi_epsilon_two, ssi_epsilon_three
  public ssi_efold_primitive, find_ssitraj
  public ssi136_x_epsilon2_Equals_0, ssi245_x_V_Equals_0, ssi3456_x_Vprime_Equals_0

contains


!returns V/M^4
  function ssi_norm_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssi_norm_potential = 1._kp+alpha*x**2+beta*x**4

  end function ssi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function ssi_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta

    ssi_norm_deriv_potential = 2._kp*x*(alpha+2._kp*beta*x**2)

  end function ssi_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M^4
  function ssi_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta

   
    ssi_norm_deriv_second_potential = 2._kp*(alpha+6._kp*beta*x**2)
    

  end function ssi_norm_deriv_second_potential



!epsilon_one(x)
  function ssi_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssi_epsilon_one
    real(kp), intent(in) :: x,alpha,beta


    ssi_epsilon_one = 2._kp*(alpha*x+2._kp*beta*x**3)**2/(1._kp+alpha*x**2+beta*x**4)**2
    
  end function ssi_epsilon_one


!epsilon_two(x)
  function ssi_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssi_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    
    ssi_epsilon_two = 4._kp*(-alpha+(alpha**2-6._kp*beta)*x**2+alpha*beta*x**4+ & 
                      2._kp*beta**2*x**6)/(1._kp+alpha*x**2+beta*x**4)**2
    
  end function ssi_epsilon_two


!epsilon_three(x)
  function ssi_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssi_epsilon_three
    real(kp), intent(in) :: x,alpha,beta

    
    ssi_epsilon_three = (4._kp*x**2*(alpha+2._kp*beta*x**2)*(-3._kp*alpha**2+6._kp*beta+ & 
                        alpha*(alpha**2-12._kp*beta)*x**2+3._kp*(alpha**2-8._kp*beta)*beta* &
                        x**4+2._kp*beta**3*x**8))/((1._kp+alpha*x**2+beta*x**4)**2* &
                        (-alpha+(alpha**2-6._kp*beta)*x**2+alpha*beta*x**4+2._kp*beta**2*x**6))
    
  end function ssi_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function ssi_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssi_efold_primitive


    if (alpha*beta .eq. 0._kp) stop 'ssi_efold_primitive: alpha*beta=0!'

    ssi_efold_primitive = 1._kp/(2._kp*alpha)*log(x)+x**2/8._kp &
                          +(alpha**2-4._kp*beta)/(16._kp*alpha*beta) &
                          *log(abs(1._kp+2._kp*beta/alpha*x**2)) 

  end function ssi_efold_primitive



  function find_ssitraj(x,ssiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssiData
    real(kp) :: find_ssitraj
    real(kp) :: alpha,beta,NplusNuend

    alpha = ssiData%real1
    beta = ssiData%real2
    NplusNuend = ssiData%real3

    find_ssitraj = ssi_efold_primitive(x,alpha,beta) - NplusNuend
   
  end function find_ssitraj

! Return the position x at which epsilon2 vanishes for SSI1, SSI3, SSI6
  function ssi136_x_epsilon2_Equals_0(alpha,beta)
  implicit none
    real(kp) :: ssi136_x_epsilon2_Equals_0
    real(kp), intent(in) :: alpha,beta

    ssi136_x_epsilon2_Equals_0 = sqrt(-(alpha/(6._kp*beta))+(16._kp*alpha**3*beta**3+2._kp* &
                               sqrt(complex((64._kp*alpha**6+(5._kp*alpha**2-36._kp*beta)**3)* &
                               beta**6,0._kp)))**(1._kp/3._kp)/(6._kp*2._kp**(1._kp/3._kp)*beta**2)+ &
                               (-5._kp*alpha**2+36._kp*beta)/(6._kp*(8._kp*alpha**3*beta**3+ & 
                               3._kp*sqrt(3._kp)*sqrt(complex((alpha**2-4._kp*beta)*beta**6* &
                               (7._kp*alpha**4-72._kp*alpha**2*beta+432._kp*beta**2),0._kp)))**(1._kp/3._kp)))

   end function ssi136_x_epsilon2_Equals_0


! Return the position x at which the potential vanishes for SSI2, SSI4 and SSI5
  function ssi245_x_V_Equals_0(alpha,beta)
  implicit none
    real(kp) :: ssi245_x_V_Equals_0
    real(kp), intent(in) :: alpha,beta

    ssi245_x_V_Equals_0 = sqrt(-(alpha+sqrt(alpha**2-4._kp*beta))/(2._kp*beta))

   end function ssi245_x_V_Equals_0


! Return the position x at which the first derivative of the potential with respect to x vanishes for SSI3&4 and SSI5&6
  function ssi3456_x_Vprime_Equals_0(alpha,beta)
  implicit none
    real(kp) :: ssi3456_x_Vprime_Equals_0
    real(kp), intent(in) :: alpha,beta

    ssi3456_x_Vprime_Equals_0 = sqrt(-alpha/(2._kp*beta))

   end function ssi3456_x_Vprime_Equals_0


end module ssicommon
