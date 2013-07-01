!slow-roll functions for the Kähler moduli inflation II potential
!
!V(phi) = M^4 { 1 - alpha x^(4/3) exp[-beta x^(4/3)] }
!
!x = phi/Mp

module kmiiisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent, easydverk
  use specialinf, only : lambert
  implicit none

  private

  public  kmiii_norm_potential, kmiii_epsilon_one, kmiii_epsilon_two, kmiii_epsilon_three
  public  kmiii_x_endinf, kmiii_efold_primitive, kmiii_x_trajectory
  public  kmiii_norm_deriv_potential, kmiii_norm_deriv_second_potential
  public  kmiii_x_endinf_appr


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

! Returns the location x where the potential vanishes (in its increasing branch) 
  function kmiii_x_potzero(alpha,beta) 
    implicit none
    real(kp) :: kmiii_x_potzero
    real(kp), intent(in) :: alpha,beta

    kmiii_x_potzero=(-1._kp/beta*lambert(-beta/alpha,-1))**(0.75_kp)

  end function kmiii_x_potzero


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

   mini = kmiii_x_potzero(alpha,beta)
   maxi = 1._kp/epsilon(1._kp)
   kmiii_x_endinf = zbrent(find_kmiii_x_endinf,mini,maxi,tolFind,kmiiiData)

  end function kmiii_x_endinf


  function find_kmiii_x_endinf(x,kmiiiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: kmiiiData
    real(kp) :: find_kmiii_x_endinf
    real(kp) :: alpha,beta
    
    alpha = kmiiiData%real1
    beta = kmiiiData%real2
    
    find_kmiii_x_endinf = kmiii_epsilon_one(x,alpha,beta)-1._kp
    
  end function find_kmiii_x_endinf

!Approximated analytical formula for xend, in the large field limit
  function kmiii_x_endinf_appr(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: kmiii_x_endinf_appr

    kmiii_x_endinf_appr=(-1.25_kp/beta*lambert(-4._kp* &
           9._kp**(2._kp/5._kp)/(5._kp*8._kp**(2._kp/5._kp))* &
           alpha**(-4._kp/5._kp)*beta**(1._kp/5._kp),-1))**(3._kp/4._kp)

  end function kmiii_x_endinf_appr



  !this is integral[V(phi)/V'(phi) dphi]
  function kmiii_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: kmiii_efold_primitive
    type(transfert) :: kmiiiData

    real(kp), parameter :: tolInt = tolkp
    integer, parameter :: neq = 1

    real(kp) :: xvar
    real(kp), dimension(neq) :: yvar
    real(kp) :: primapprox

    real(kp), parameter :: switchToExactAt = 1._kp

!approx in the beta x^(4/3) >> 1
!and alpha x^(4/3) exp(-beta x^(4/3) ) << 1 region

    if (x**(4._kp/3._kp)*beta.gt.switchToExactAt) then

       kmiii_efold_primitive = 9._kp/(16._kp*alpha*beta**2) &
            *exp(beta*x**(4._kp/3._kp) - 2._kp*log(x))   

       return
    endif

    write(*,*)'kmiii_efold_primitive:'
    write(*,*)'switching to exact slow-roll trajectory...'

    xvar = beta**(-3._kp/4._kp)
    yvar(1) = 9._kp/(16._kp*alpha*beta**2) &
         *exp(beta*xvar**(4._kp/3._kp) - 2._kp*log(xvar))
   
    kmiiiData%real1 = alpha
    kmiiiData%real2 = beta
  
    call easydverk(neq,find_kmiii_efold_primitive,xvar,yvar,x,tolInt,kmiiiData)

    kmiii_efold_primitive = yvar(1)

!    print *,'test',kmiii_efold_primitive,9._kp/(16._kp*alpha*beta**2) &
!            *exp(beta*x**(4._kp/3._kp) - 2._kp*log(x))
    
  end function kmiii_efold_primitive


  subroutine find_kmiii_efold_primitive(n,x,y,yprime,kmiiiData)
    implicit none          
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: kmiiiData
    real(kp) :: alpha,beta

    alpha = kmiiiData%real1
    beta = kmiiiData%real2
!regularized to avoid eps1=0
    yprime(1) = 1._kp/sqrt(epsilon(1._kp) + 2._kp*kmiii_epsilon_one(x,alpha,beta))


  end subroutine find_kmiii_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend, if
!inflation proceeds from the right to the left
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
    
    kmiii_x_trajectory = zbrent(find_kmiii_x_trajectory,mini,maxi,tolFind,kmiiiData)
       
  end function kmiii_x_trajectory

  function find_kmiii_x_trajectory(x,kmiiiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: kmiiiData
    real(kp) :: find_kmiii_x_trajectory
    real(kp) :: alpha,beta,NplusNuend

    alpha = kmiiiData%real1
    beta = kmiiiData%real2
    NplusNuend = kmiiiData%real3

    find_kmiii_x_trajectory = kmiii_efold_primitive(x,alpha,beta) - NplusNuend
   
  end function find_kmiii_x_trajectory
  
end module kmiiisr
