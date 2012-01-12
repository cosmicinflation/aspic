!slow-roll functions for the Kähler moduli inflation II potential
!
!V(phi) = M^4 [1-alpha x^(4/3) exp(-beta x^(4/3) ) ]
!
!x = phi/Mp

module kmiiisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : lambert
  implicit none

  private

  public  kmiii_norm_potential, kmiii_epsilon_one, kmiii_epsilon_two, kmiii_epsilon_three
  public  kmiii_x_endinf, kmiii_efold_primitive, kmiii_x_trajectory
  public  kmiii_norm_deriv_potential, kmiii_norm_deriv_second_potential
  public  kmiii_alphamin
 
contains
!returns V/M^4
  function kmiii_norm_potential(x,alpha,beta)
    implicit none
    real(kp) :: kmiii_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    kmiii_norm_potential = 1._kp-alpha*x**(4._kp/3._kp)*exp(-beta*x**(4._kp/3._kp))

  end function kmiii_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function kmiii_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: kmiii_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta

   kmiii_norm_deriv_potential = -4._kp/3._kp*alpha*exp(-beta*x**(4._kp/3._kp))* &
                   x**(1./3.)*(1._kp-beta*x**(4._kp/3._kp))

  end function kmiii_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function kmiii_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: kmiii_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta

    kmiii_norm_deriv_second_potential = 4._kp*alpha/9._kp*exp(-beta*x**(4._kp/3._kp)) &
                    *(9._kp*beta*x**(2._kp/3._kp)-4._kp*beta**2*x**2-x**(-2._kp/3._kp))

  end function kmiii_norm_deriv_second_potential



!epsilon_one(x)
  function kmiii_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: kmiii_epsilon_one
    real(kp), intent(in) :: x,alpha,beta
    
    kmiii_epsilon_one = 8._kp*alpha**2*x**(2._kp/3._kp)*(-1._kp+beta*x**(4._kp/3._kp))**2 &
                        /(9._kp*(exp(beta*x**(4._kp/3._kp))-alpha*x**(4._kp/3._kp))**2)
    
  end function kmiii_epsilon_one


!epsilon_two(x)
  function kmiii_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: kmiii_epsilon_two
    real(kp), intent(in) :: x,alpha,beta
    
    kmiii_epsilon_two =(8._kp*alpha*(alpha*x**(4._kp/3._kp)*(3._kp+beta*x**(4._kp/3._kp))+ & 
                       exp(beta*x**(4._kp/3._kp))*(1._kp-9._kp*beta*x**(4._kp/3._kp)+4._kp &
                       *beta**2*x**(8._kp/3._kp))))/(9._kp*x**(2._kp/3._kp)*(exp(beta* &
                       x**(4._kp/3._kp))-alpha*x**(4._kp/3._kp))**2)
    
  end function kmiii_epsilon_two


!epsilon_three(x)
  function kmiii_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: kmiii_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
    
    kmiii_epsilon_three =(8._kp*alpha*(-1._kp+beta*x**(4._kp/3._kp))*(-alpha**2*x**(8._kp/3._kp) &
                         *(9._kp+beta*x**(4._kp/3._kp))+2._kp*alpha*exp(beta*x**(4._kp/3._kp))* &
                         x**(4/3)*(-4._kp+19._kp*beta*x**(4._kp/3._kp)-9._kp*beta**2* &
                         x**(8._kp/3._kp)+4._kp*beta**3*x**4)+exp(2._kp*beta*x**(4._kp/3._kp))* &
                         (1._kp+11._kp*beta*x**(4._kp/3._kp)-30._kp*beta**2*x**(8._kp/3._kp)+ &
                         8._kp*beta**3*x**4)))/(9._kp*x**(2._kp/3._kp)*(exp(beta*x**(4._kp/3._kp))- &
                         alpha*x**(4._kp/3._kp))**2*(alpha*x**(4._kp/3._kp)*(3._kp+beta* &
                         x**(4._kp/3._kp))+exp(beta*x**(4._kp/3._kp))*(1._kp+beta*x**(4._kp/3._kp)* &
                         (-9._kp+4._kp*beta*x**(4._kp/3._kp)))))
  end function kmiii_epsilon_three


!returns x at the end of inflation defined as epsilon1=1, if inflation proceeds from the right to the left
  function kmiii_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: kmiii_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,xminus,xplus
    type(transfert) :: kmiiiData

    kmiiiData%real1 = alpha
    kmiiiData%real2 = beta

    !Calculates the position of the two extrema of epsilon1
    IF (beta/alpha.lt.exp(-1._kp)) then
      xminus=(-1._kp/beta*lambert(-beta/alpha,0))**(3._kp/4._kp)
      xplus=(-1._kp/beta*lambert(-beta/alpha,-1))**(3._kp/4._kp)
    ELSE
      xminus=kmiii_xminus_eps1max(alpha,beta)
      xplus=kmiii_xplus_eps1max(alpha,beta)
    ENDIF

    IF(kmiii_epsilon_one(xplus,alpha,beta).gt.1._kp) THEN
          mini = xplus
          maxi = 1._kp/epsilon(1._kp)
          kmiii_x_endinf = zbrent(find_kmiiiendinf,mini,maxi,tolFind,kmiiiData)

    ELSE IF(kmiii_epsilon_one(xminus,alpha,beta).gt.1._kp) THEN
          mini = xminus
          maxi = xplus
          kmiii_x_endinf = zbrent(find_kmiiiendinf,mini,maxi,tolFind,kmiiiData)

    ELSE

	STOP 'module kmiiisr: function kmiii_x_endinf: Inflation cannot stop by violation of slow-roll in KMIII for this set of parameters'

    ENDIF

  
  end function kmiii_x_endinf


  function find_kmiiiendinf(x,kmiiiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: kmiiiData
    real(kp) :: find_kmiiiendinf
    real(kp) :: alpha,beta
    
    alpha = kmiiiData%real1
    beta = kmiiiData%real2
    
    find_kmiiiendinf = kmiii_epsilon_one(x,alpha,beta)-1._kp
    
  end function find_kmiiiendinf


!Returns the position of the smallest maximum of epsilon1
  function kmiii_xminus_eps1max(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: kmiii_xminus_eps1max
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kmiiiData

    mini = epsilon(1._kp)
    maxi = (1._kp/beta)**(3._kp/4._kp) !exact position of the minimum of epsilon1

    kmiiiData%real1 = alpha
    kmiiiData%real2 = beta

    kmiii_xminus_eps1max = zbrent(find_kmiiieps1max,mini,maxi,tolFind,kmiiiData)

  end function kmiii_xminus_eps1max

!Returns the position of the highest maximum of epsilon1
  function kmiii_xplus_eps1max(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: kmiii_xplus_eps1max
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kmiiiData

    mini = (1._kp/beta)**(3._kp/4._kp) !exact position of the minimum of epsilon1
    maxi = mini*10._kp

    kmiiiData%real1 = alpha
    kmiiiData%real2 = beta

    kmiii_xplus_eps1max = zbrent(find_kmiiieps1max,mini,maxi,tolFind,kmiiiData)

  end function kmiii_xplus_eps1max

  function find_kmiiieps1max(x,kmiiiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: kmiiiData
    real(kp) :: find_kmiiieps1max
    real(kp) :: alpha,beta
    
    alpha = kmiiiData%real1
    beta = kmiiiData%real2
    
    find_kmiiieps1max = (alpha*x**(4._kp/3._kp)*(3._kp+beta*x**(4._kp/3._kp))+ &
                        exp(beta*x**(4._kp/3._kp))*(1._kp-9._kp*beta*x**(4._kp/3._kp) &
                        +4._kp*beta**2*x**(8._kp/3._kp)))/(27._kp*x**(1._kp/3._kp) &
                        *(-exp(beta*x**(4._kp/3._kp))+alpha*x**(4._kp/3._kp))**3)
    
  end function find_kmiiieps1max


!this is integral[V(phi)/V'(phi) dphi], valid in the beta x^(4/3) >> 1 and alpha x^(4/3) exp(-beta x^(4/3) ) << 1 region
  function kmiii_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: kmiii_efold_primitive

    if (alpha.eq.0._kp.and.beta.eq.0._kp) stop 'kmiii_efold_primitive: alpha=0 or beta=0!'
!1/x^2 -> -2 log(x) to prevent overflowing the exponential
    kmiii_efold_primitive = 9._kp/(16._kp*alpha*beta**2)*exp(beta*x**(4._kp/3._kp) - 2._kp*log(x))

  end function kmiii_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend, if inflation proceeds from the right to the left
  function kmiii_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha,beta, xend
    real(kp) :: kmiii_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kmiiiData

  
    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)
  

    kmiiiData%real1 = alpha
    kmiiiData%real2 = beta
    kmiiiData%real3 = -bfold + kmiii_efold_primitive(xend,alpha,beta)
    
    kmiii_x_trajectory = zbrent(find_kmiiitraj,mini,maxi,tolFind,kmiiiData)
       
  end function kmiii_x_trajectory

  function find_kmiiitraj(x,kmiiiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: kmiiiData
    real(kp) :: find_kmiiitraj
    real(kp) :: alpha,beta,NplusNuend

    alpha = kmiiiData%real1
    beta = kmiiiData%real2
    NplusNuend = kmiiiData%real3

    find_kmiiitraj = kmiii_efold_primitive(x,alpha,beta) - NplusNuend
   
  end function find_kmiiitraj

  function kmiii_alphamin(beta) !Given beta, returns the minimum value of alpha such that inflation can be stopped by violation of the slow roll conditions
    implicit none
    real(kp), intent(in) :: beta
    real(kp) ::kmiii_alphamin,alphamini,alphamaxi
    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: kmiiiData

    alphamini=epsilon(1._kp)
    alphamaxi=beta*exp(1._kp) - epsilon(1._kp) !maximum allowed value such that the potential is positive everywhere (with a numerical safety)
          

    if(kmiii_epsilon_one(kmiii_xplus_eps1max(alphamaxi,beta),alphamaxi,beta) .lt. 1._kp) then !in that case the prior space in empty
         kmiii_alphamin=beta*exp(1._kp)
    else
         kmiiiData%real1 = beta
         kmiii_alphamin=zbrent(find_kmiii_alphamin,alphamini,alphamaxi,tolFind,kmiiiData)

    endif

  end function kmiii_alphamin

  function find_kmiii_alphamin(x,kmiiiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: kmiiiData
    real(kp) :: find_kmiii_alphamin
    real(kp) :: beta

    beta = kmiiiData%real1

    find_kmiii_alphamin = kmiii_epsilon_one(kmiii_xplus_eps1max(x,beta),x,beta)-1._kp
   
  end function find_kmiii_alphamin


  
end module kmiiisr
