!slow-roll functions for the MSSMI and GMSSMI potential
!("GMSSMI" means that the condition alpha=1 is relaxed)
!
!V(phi) = M**4 [ x**2 - 2/3 alpha x**6 + alpha/5 x**10 ]
!
!x = phi/phi0
!phi0=phi0/Mp


module gmssmicommon
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public gmssmi_norm_potential, gmssmi_norm_deriv_potential, gmssmi_norm_deriv_second_potential
  public gmssmi_epsilon_one, gmssmi_epsilon_two, gmssmi_epsilon_three
  public gmssmi_x_epsonemin, gmssmi_x_endinf
  public gmssmi_x_epsonezero, gmssmi_x_epstwozero


contains

 
!returns V/M**4
  function gmssmi_norm_potential(x,alpha,phi0)
    implicit none
    real(kp) :: gmssmi_norm_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: phi0
    
    gmssmi_norm_potential = x**2-2._kp/3._kp*alpha*x**6+alpha/5._kp*x**10

  end function gmssmi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function gmssmi_norm_deriv_potential(x,alpha,phi0)
    implicit none
    real(kp) :: gmssmi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: phi0
  
   gmssmi_norm_deriv_potential = 2._kp*(x-2._kp*alpha*x**5+alpha*x**9)

  end function gmssmi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function gmssmi_norm_deriv_second_potential(x,alpha,phi0)
    implicit none
    real(kp) :: gmssmi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: phi0
   
    gmssmi_norm_deriv_second_potential = 2._kp+2._kp*alpha*x**4*(-10._kp+9._kp*x**4)

  end function gmssmi_norm_deriv_second_potential



!epsilon_one(x)
  function gmssmi_epsilon_one(x,alpha,phi0)    
    implicit none
    real(kp) :: gmssmi_epsilon_one
    real(kp), intent(in) :: x,alpha,phi0
      
    gmssmi_epsilon_one =(450._kp*(1._kp+alpha*x**4*(-2._kp+x**4))**2)/ &
         (phi0**2*x**2*(15._kp+alpha*x**4*(-10._kp+3._kp*x**4))**2)
    
  end function gmssmi_epsilon_one


!epsilon_two(x)
  function gmssmi_epsilon_two(x,alpha,phi0)    
    implicit none
    real(kp) :: gmssmi_epsilon_two
    real(kp), intent(in) :: x,alpha,phi0
     
    gmssmi_epsilon_two =(60._kp*(15._kp+alpha*x**4*(40._kp+x**4*(-78._kp+alpha* &
         (20._kp+3._kp*x**8)))))/(phi0**2*x**2*(15._kp+ &
         alpha*x**4*(-10._kp+3._kp*x**4))**2)
    
  end function gmssmi_epsilon_two


!epsilon_three(x)
  function gmssmi_epsilon_three(x,alpha,phi0)    
    implicit none
    real(kp) :: gmssmi_epsilon_three
    real(kp), intent(in) :: x,alpha,phi0
       
    gmssmi_epsilon_three = (60._kp*(1._kp+alpha*x**4*(-2._kp+x**4))*(225._kp+ &
         alpha*x**4*(-1350._kp+x**4*(3915._kp+alpha*(-2100._kp+ &
         20._kp*(81._kp-10._kp*alpha)*x**4+15._kp*(-99._kp+20._kp*alpha)* &
         x**8+90._kp*alpha*x**(12)+9._kp*alpha*x**(16))))))/(x**2*(15._kp+ &
         alpha*x**4*(-10._kp+3._kp*x**4))**2*(15._kp+ &
         alpha*x**4*(40._kp+x**4*(-78._kp+alpha*(20._kp+3._kp*x**8)))))/phi0**2

  end function gmssmi_epsilon_three


 function gmssmi_x_epstwozero(alpha) 
    real(kp), dimension(2) :: gmssmi_x_epstwozero
    real(kp), intent(in) :: alpha
    complex(kp) :: delta,BigDelta,sigma,BigSigma,x_eps2NULMinus, x_eps2NULPlus

    if (alpha .gt. 9._kp/5._kp) then

       gmssmi_x_epstwozero = -1._kp !error value

    else

       delta=(736._kp*alpha**2)/25._kp-(208._kp*alpha**3)/15._kp+(16._kp*alpha**4)/9._kp
       BigDelta=-((430336._kp*alpha**4)/625._kp)+(612352._kp*alpha**5)/1125._kp- &
            (20992._kp*alpha**6)/225._kp+(256._kp*alpha**8)/243._kp
       sigma=-((12896._kp*alpha**3)/125._kp)+(2944._kp*alpha**4)/25._kp- &
            (416._kp*alpha**5)/15._kp+(64._kp*alpha**6)/27._kp+6._kp* &
            sqrt(15._kp)*alpha/5._kp*sqrt(BigDelta)
       BigSigma=(52._kp*alpha)/5._kp-(8._kp*alpha**2)/3._kp+delta/ &
            (sigma**(1._kp/3._kp))+sigma**(1._kp/3._kp)

       x_eps2NULMinus=(1._kp/(2._kp*alpha)*sqrt(5._kp/3._kp)*(sqrt(BigSigma)-2._kp* &
            sqrt(39._kp*alpha/5._kp-2._kp*alpha**2-BigSigma/4._kp- &
            12._kp*alpha**2/sqrt(15._kp*BigSigma))))**(0.25_kp)
       x_eps2NULPlus=(1._kp/(2._kp*alpha)*sqrt(5._kp/3._kp)*(sqrt(BigSigma)+2._kp* &
            sqrt(39._kp*alpha/5._kp-2._kp*alpha**2-BigSigma/4._kp- &
            12._kp*alpha**2/sqrt(15._kp*BigSigma))))**(0.25_kp)

       gmssmi_x_epstwozero(1) = real(x_eps2NULMinus,kp)
       gmssmi_x_epstwozero(2) = real(x_eps2NULPlus,kp)

    endif


  end function gmssmi_x_epstwozero


   function gmssmi_x_epsonezero(alpha)
    real(kp), dimension(2) :: gmssmi_x_epsonezero
    real(kp), intent(in) :: alpha 

    if (alpha .lt. 1._kp) then

       gmssmi_x_epsonezero=-1._kp !error value

    else

       gmssmi_x_epsonezero(1) = (1._kp-sqrt(1._kp-1._kp/alpha))**(0.25_kp)
       gmssmi_x_epsonezero(2) = (1._kp+sqrt(1._kp-1._kp/alpha))**(0.25_kp)

    endif

  end function gmssmi_x_epsonezero



!Returns the position of the first local minimum of epsilon1
  function gmssmi_x_epsonemin(alpha)   
    implicit none
    real(kp) :: gmssmi_x_epsonemin
    real(kp), intent(in) :: alpha
    real(kp), dimension(2) :: xEpsTwoZero, xEpsOneZero
    
    if (alpha .lt. 1._kp) then

       xEpsTwoZero = gmssmi_x_epstwozero(alpha)

       gmssmi_x_epsonemin = xEpsTwoZero(1)

    else if (alpha .eq. 1._kp) then

       gmssmi_x_epsonemin =  1._kp

    else

       xEpsOneZero = gmssmi_x_epsonezero(alpha)
       
       gmssmi_x_epsonemin = xEpsOneZero(1)

    endif

  end function gmssmi_x_epsonemin
  
  
  !returns x at the end of inflation defined as epsilon1=1
  function gmssmi_x_endinf(alpha,phi0)
    implicit none
    real(kp), intent(in) :: alpha,phi0

    real(kp) :: gmssmi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmssmiData


    mini = epsilon(1._kp)
    maxi = gmssmi_x_epsonemin(alpha)*(1._kp-epsilon(1._kp)) !Position of the first local minimum of epsilon1

    gmssmiData%real1 = alpha
    gmssmiData%real2 = phi0

    gmssmi_x_endinf = zbrent(find_gmssmi_x_endinf,mini,maxi,tolFind,gmssmiData)


  end function gmssmi_x_endinf

  function find_gmssmi_x_endinf(x,gmssmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gmssmiData
    real(kp) :: find_gmssmi_x_endinf
    real(kp) :: alpha,phi0

    alpha = gmssmiData%real1
    phi0 = gmssmiData%real2

    find_gmssmi_x_endinf = gmssmi_epsilon_one(x,alpha,phi0)-1._kp

  end function find_gmssmi_x_endinf


end module gmssmicommon


