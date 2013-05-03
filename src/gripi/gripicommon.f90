!slow-roll functions for the GRIP inflation potential
!
!V(phi) = M**4 [ x**2 - 4/3 alpha x**3 + alpha/2 x**4 ]
!
!x = phi/phi0

module gripicommon
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent

  implicit none

  private

  public gripi_norm_potential, gripi_norm_deriv_potential
  public gripi_norm_deriv_second_potential
  public gripi_epsilon_one, gripi_epsilon_two, gripi_epsilon_three  
  public gripi_x_epsonemin, gripi_x_endinf
  public gripi_x_epstwozero,gripi_x_epsonezero


contains

!returns V/M**4
  function gripi_norm_potential(x,alpha,phi0)
    implicit none
    real(kp) :: gripi_norm_potential
    real(kp), intent(in) :: x, alpha
    real(kp), intent(in), optional :: phi0

    gripi_norm_potential = x**2-4._kp/3._kp*alpha*x**3+0.5_kp*alpha*x**4

  end function gripi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function gripi_norm_deriv_potential(x,alpha,phi0)
    implicit none
    real(kp) :: gripi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in), optional :: phi0

    gripi_norm_deriv_potential = 2._kp*x*(1._kp+alpha*(-2._kp+x)*x)

  end function gripi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function gripi_norm_deriv_second_potential(x,alpha,phi0)
    implicit none
    real(kp) :: gripi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in), optional :: phi0

    gripi_norm_deriv_second_potential = 2._kp+2._kp*alpha*x*(-4._kp+3._kp*x)

  end function gripi_norm_deriv_second_potential


!epsilon_one(x)
  function gripi_epsilon_one(x,alpha,phi0)    
    implicit none
    real(kp) :: gripi_epsilon_one
    real(kp), intent(in) :: x,alpha,phi0

    gripi_epsilon_one =(72._kp*(1._kp+alpha*(-2._kp+x)*x)**2)/ &
         (x**2*(6._kp+alpha*x*(-8._kp+3._kp*x))**2)/phi0**2

  end function gripi_epsilon_one


!epsilon_two(x)
  function gripi_epsilon_two(x,alpha,phi0)    
    implicit none
    real(kp) :: gripi_epsilon_two
    real(kp), intent(in) :: x,alpha,phi0

    gripi_epsilon_two =(24._kp*(6._kp+alpha*x*(-16._kp+x*(3._kp+alpha* &
         (16._kp+3._kp*(-4._kp+x)*x)))))/(x**2* &
         (6._kp+alpha*x*(-8._kp+3._kp*x))**2)/phi0**2

  end function gripi_epsilon_two


!epsilon_three(x)
  function gripi_epsilon_three(x,alpha,phi0)    
    implicit none
    real(kp) :: gripi_epsilon_three
    real(kp), intent(in) :: x,alpha,phi0

    gripi_epsilon_three = (24._kp*(1._kp+alpha*(-2._kp+x)*x)*(36._kp+alpha*x* & 
         (-144._kp+x*(54._kp+alpha*(192._kp+x*(-108._kp+alpha* &
         (-128._kp+9._kp*x*(16._kp+(-6._kp+x)*x))))))))/(x**2* &
         (6._kp+alpha*x*(-8._kp+3._kp*x))**2*(6._kp+alpha*x* &
         (-16._kp+x*(3._kp+alpha*(16._kp+3._kp*(-4._kp+x)*x)))))/ &
         phi0**2

  end function gripi_epsilon_three

  function gripi_x_epstwozero(alpha) 
    real(kp), dimension(2) :: gripi_x_epstwozero
    real(kp), intent(in) :: alpha

    if (alpha .gt. 9._kp/8._kp) then

       gripi_x_epstwozero = -1._kp !error value

    else

! Lower solution of eps2=0
       gripi_x_epstwozero(1)=1._kp/6._kp*(6._kp+sqrt(4._kp-6._kp/alpha+1._kp/alpha* &
            (-1917._kp+9504._kp*alpha-11520._kp*alpha**2+4096._kp* &
            alpha**3+18._kp*sqrt(2._kp)*sqrt(-11907._kp+56268._kp* &
            alpha-92448._kp*alpha**2+64512._kp*alpha**3- &
            16384._kp*alpha**4))**(1._kp/3._kp)+(15._kp-16._kp* &
            alpha)**2/(alpha**3*(-1917._kp+9504._kp*alpha- &
            11520._kp*alpha**2+4096._kp*alpha**3+18._kp*sqrt(2._kp)* &
            sqrt(-11907._kp+56268._kp*alpha-92448._kp*alpha**2+ &
            64512._kp*alpha**3-16384._kp*alpha**4)))**(1._kp/ &
            3._kp))-sqrt(8._kp-12._kp/alpha-(15._kp- &
            16._kp*alpha)**2/(alpha*(-1917._kp+9504._kp* &
            alpha-11520._kp*alpha**2+4096._kp*alpha**3+18._kp* &
            sqrt(2._kp)*sqrt(-11907._kp+56268._kp*alpha- &
            92448._kp*alpha**2+64512._kp*alpha**3- &
            16384._kp*alpha**4))**(1._kp/3._kp))-1._kp/alpha* &
            (-1917._kp+9504._kp*alpha-11520._kp*alpha**2+4096._kp* &
            alpha**3+18._kp*sqrt(2._kp)*sqrt(-11907._kp+56268._kp* &
            alpha-92448._kp*alpha**2+64512._kp*alpha**3-16384._kp* &
            alpha**4))**(1._kp/3._kp)+(36._kp*(-4._kp+5._kp/ &
            alpha))/(sqrt(4._kp-6._kp/alpha+1._kp/alpha*(-1917._kp &
            +9504._kp*alpha-11520._kp*alpha**2+4096._kp*alpha**3 &
            +18._kp*sqrt(2._kp)*sqrt(-11907._kp+56268._kp*alpha- &
            92448._kp*alpha**2+64512._kp*alpha**3-16384._kp* &
            alpha**4))**(1._kp/3._kp)+(15._kp-16._kp*alpha)**2/ &
            (alpha**3*(-1917._kp+9504._kp*alpha-11520._kp*alpha**2 &
            +4096._kp*alpha**3+18*sqrt(2._kp)*sqrt(-11907._kp+ &
            56268._kp* alpha-92448._kp*alpha**2+64512._kp*alpha**3 &
            - 16384._kp*alpha**4)))**(1._kp/3._kp)))))

! Upper solution of eps2=0 
       gripi_x_epstwozero(2)=1._kp/6._kp*(6._kp+sqrt(4._kp-6._kp/alpha+1._kp/alpha* &
            (-1917._kp+9504._kp*alpha-11520._kp*alpha**2+4096._kp* &
            alpha**3+18._kp*sqrt(2._kp)*sqrt(-11907._kp+56268._kp* &
            alpha-92448._kp*alpha**2+64512._kp*alpha**3- &
            16384._kp*alpha**4))**(1._kp/3._kp)+(15._kp-16._kp* &
            alpha)**2/(alpha**3*(-1917._kp+9504._kp*alpha- &
            11520._kp*alpha**2+4096._kp*alpha**3+18._kp*sqrt(2._kp)* &
            sqrt(-11907._kp+56268._kp*alpha-92448._kp*alpha**2+ &
            64512._kp*alpha**3-16384._kp*alpha**4)))**(1._kp/ &
            3._kp))+sqrt(8._kp-12._kp/alpha-(15._kp- &
            16._kp*alpha)**2/(alpha*(-1917._kp+9504._kp* &
            alpha-11520._kp*alpha**2+4096._kp*alpha**3+18._kp* &
            sqrt(2._kp)*sqrt(-11907._kp+56268._kp*alpha- &
            92448._kp*alpha**2+64512._kp*alpha**3- &
            16384._kp*alpha**4))**(1._kp/3._kp))-1._kp/alpha* &
            (-1917._kp+9504._kp*alpha-11520._kp*alpha**2+4096._kp* &
            alpha**3+18._kp*sqrt(2._kp)*sqrt(-11907._kp+56268._kp* &
            alpha-92448._kp*alpha**2+64512._kp*alpha**3-16384._kp* &
            alpha**4))**(1._kp/3._kp)+(36._kp*(-4._kp+5._kp/ &
            alpha))/(sqrt(4._kp-6._kp/alpha+1._kp/alpha*(-1917._kp &
            +9504._kp*alpha-11520._kp*alpha**2+4096._kp*alpha**3 &
            +18._kp*sqrt(2._kp)*sqrt(-11907._kp+56268._kp*alpha- &
            92448._kp*alpha**2+64512._kp*alpha**3-16384._kp* &
            alpha**4))**(1._kp/3._kp)+(15._kp-16._kp*alpha)**2/ &
            (alpha**3*(-1917._kp+9504._kp*alpha-11520._kp*alpha**2 &
            +4096._kp*alpha**3+18*sqrt(2._kp)*sqrt(-11907._kp+ &
            56268._kp* alpha-92448._kp*alpha**2+64512._kp*alpha**3 &
            - 16384._kp*alpha**4)))**(1._kp/3._kp)))))

    endif


  end function gripi_x_epstwozero

  function gripi_x_epsonezero(alpha)
    real(kp), dimension(2) :: gripi_x_epsonezero
    real(kp), intent(in) :: alpha 

    if (alpha .lt. 1._kp) then

       gripi_x_epsonezero=-1._kp !error value

    elseif (alpha.eq.0._kp) then

       stop 'gripi_x_epsonezero: alpha=0 is singular!'

    else

       gripi_x_epsonezero(1) = 1._kp-sqrt(1._kp-1._kp/alpha)
       gripi_x_epsonezero(2) = 1._kp+sqrt(1._kp-1._kp/alpha)

    endif

  end function gripi_x_epsonezero

!Returns the position of the first local minimum of epsilon1
  function gripi_x_epsonemin(alpha)   
    implicit none
    real(kp) :: gripi_x_epsonemin
    real(kp), intent(in) :: alpha
    real(kp), dimension(2) :: xEpsOneZero,xEpsTwoZero

    if (alpha .lt. 1._kp) then

       xEpsTwoZero=gripi_x_epstwozero(alpha)
       gripi_x_epsonemin = xEpsTwoZero(1)

    else if (alpha .eq. 1._kp) then

       gripi_x_epsonemin =  1._kp

    else

       xEpsOneZero = gripi_x_epsonezero(alpha)

       gripi_x_epsonemin = xEpsOneZero(1)

    endif

  end function gripi_x_epsonemin

!returns x at the end of inflation defined as epsilon1=1
  function gripi_x_endinf(alpha,phi0)
    implicit none
    real(kp), intent(in) :: alpha,phi0

    real(kp) :: gripi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gripiData


    mini = epsilon(1._kp)
!Position of the first local minimum of epsilon1
    maxi = gripi_x_epsonemin(alpha)*(1._kp-epsilon(1._kp)) 

    gripiData%real1 = alpha
    gripiData%real2 = phi0

    gripi_x_endinf = zbrent(find_gripi_x_endinf,mini,maxi,tolFind,gripiData)


  end function gripi_x_endinf

  function find_gripi_x_endinf(x,gripiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gripiData
    real(kp) :: find_gripi_x_endinf
    real(kp) :: alpha,phi0

    alpha = gripiData%real1
    phi0 = gripiData%real2

    find_gripi_x_endinf = gripi_epsilon_one(x,alpha,phi0)-1._kp

  end function find_gripi_x_endinf


end module gripicommon
