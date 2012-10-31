!slow-roll functions for the MSSMI and GMSSMI potential
!("GMSSMI" means that the condition alpha^2/beta is relaxed)
!
!V(phi) = M^4 [ x^2 - alpha x^6 + beta x^10 ]
!
!x = phi/Mp


module mssmicommon
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public gmssmi_gen_norm_potential, gmssmi_gen_norm_deriv_potential, gmssmi_gen_norm_deriv_second_potential
  public gmssmi_gen_epsilon_one, gmssmi_gen_epsilon_two, gmssmi_gen_epsilon_three
  public gmssmi_gen_x_epsilon1_min, gmssmi_gen_x_endinf


contains

 
!returns V/M**4
  function gmssmi_gen_norm_potential(x,alpha,beta)
    implicit none
    real(kp) :: gmssmi_gen_norm_potential
    real(kp), intent(in) :: x,alpha,beta  
    
    gmssmi_gen_norm_potential = x**2-alpha*x**6+beta*x**10

  end function gmssmi_gen_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function gmssmi_gen_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: gmssmi_gen_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
  
   gmssmi_gen_norm_deriv_potential = 2._kp*(x-3._kp*alpha*x**5+5._kp*beta*x**9)

  end function gmssmi_gen_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function gmssmi_gen_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: gmssmi_gen_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta
   
    gmssmi_gen_norm_deriv_second_potential = 2._kp*(1._kp-15._kp*alpha*x**4+45._kp*beta*x**8)

  end function gmssmi_gen_norm_deriv_second_potential



!epsilon_one(x)
  function gmssmi_gen_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: gmssmi_gen_epsilon_one
    real(kp), intent(in) :: x,alpha,beta
      
    gmssmi_gen_epsilon_one =(2._kp*(1._kp-3._kp*alpha*x**4+5._kp*beta*x**8)**2) &
         /(x-alpha*x**5+beta*x**9)**2
    
  end function gmssmi_gen_epsilon_one


!epsilon_two(x)
  function gmssmi_gen_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: gmssmi_gen_epsilon_two
    real(kp), intent(in) :: x,alpha,beta
     
    gmssmi_gen_epsilon_two =(4._kp*(1._kp+4._kp*alpha*x**4+3._kp*alpha**2*x**8 + & 
                       beta*x**8*(-26._kp+5._kp*beta*x**8)))/(x-alpha*x**5+beta*x**9)**2
    
  end function gmssmi_gen_epsilon_two


!epsilon_three(x)
  function gmssmi_gen_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: gmssmi_gen_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
       
    gmssmi_gen_epsilon_three = (4._kp*(-1._kp+3._kp*alpha*x**4-5._kp*beta*x**8)* &
                          (-1._kp+3._kp*alpha**3*x**12+3._kp*alpha**2*x**8* &
                          (7._kp-5._kp*beta*x**8)+beta*x**8*(-87._kp-5._kp*beta* &
                          x**8*(-33._kp+beta*x**8))-3._kp*alpha*x**4*(-3._kp+ &
                          beta*x**8*(18._kp+5._kp*beta*x**8))))/((x-alpha*x**5+ &
                          beta*x**9)**2*(1._kp+4._kp*alpha*x**4+3._kp*alpha**2*x**8 + & 
                          beta*x**8*(-26._kp+5._kp*beta*x**8)))
    
  end function gmssmi_gen_epsilon_three

!Returns the position of the first local minimum of epsilon1
  function gmssmi_gen_x_epsilon1_min(alpha,beta)   
    implicit none
    real(kp) :: gmssmi_gen_x_epsilon1_min
    real(kp), intent(in) :: alpha,beta
  
    complex(kp) :: delta,BigDelta,sigma,BigSigma,x_eps2NUL

    
    if (alpha**2/beta<20._kp/9._kp) then

       delta=9._kp*alpha**4-156._kp*alpha**2*beta+736._kp*beta**2
       BigDelta=27._kp*alpha**8-11808._kp*alpha**4*beta**2+ &
            153088._kp*alpha**2*beta**3-430336._kp*beta**4
       sigma=27._kp*alpha**6-702._kp*alpha**4*beta+6624._kp*alpha**2* &
            beta**2-12896._kp*beta**3+6._kp*sqrt(15._kp)*beta*sqrt(BigDelta)
       BigSigma=-6._kp*alpha**2+52._kp*beta+delta/(sigma**(1._kp/3._kp))+sigma**(1._kp/3._kp)

       x_eps2NUL=(1._kp/(2._kp*sqrt(15._kp)*beta)*(sqrt(BigSigma)-&
            sqrt(156._kp*beta-18._kp*alpha**2-BigSigma &
            -24._kp*sqrt(15._kp)*alpha*beta/(sqrt(BigSigma)))))**(0.25_kp)

       gmssmi_gen_x_epsilon1_min = real(x_eps2NUL,kp)

    else

       gmssmi_gen_x_epsilon1_min = (3._kp*alpha/(10._kp*beta) &
            *(1._kp-sqrt(1._kp-20._kp*beta/(9._kp*alpha**2))))**(0.25_kp)

    endif
    
  end function gmssmi_gen_x_epsilon1_min


!returns x at the end of inflation defined as epsilon1=1
  function gmssmi_gen_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    
    real(kp) :: gmssmi_gen_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmssmiData

   
    mini = 0._kp
    maxi = gmssmi_gen_x_epsilon1_min(alpha,beta)*(1._kp-epsilon(1._kp)) !Position of the first local minimum of epsilon1
  
    gmssmiData%real1 = alpha
    gmssmiData%real2 = beta
    
    gmssmi_gen_x_endinf = zbrent(find_gmssmi_gen_x_endinf,mini,maxi,tolFind,gmssmiData)
   

  end function gmssmi_gen_x_endinf

  function find_gmssmi_gen_x_endinf(x,gmssmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gmssmiData
    real(kp) :: find_gmssmi_gen_x_endinf
    real(kp) :: alpha,beta

    alpha = gmssmiData%real1
    beta = gmssmiData%real2

    find_gmssmi_gen_x_endinf = gmssmi_gen_epsilon_one(x,alpha,beta)-1._kp
   
  end function find_gmssmi_gen_x_endinf




end module mssmicommon
