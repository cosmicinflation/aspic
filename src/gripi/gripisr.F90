!slow-roll functions for the GRIP inflation potential
!
!V(phi) = M**4 [ x**2 - 4/3 alpha x**3 + alpha/2 x**4 ]
!
!x = phi/phi0

module gripisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
#ifdef NOF08
  use specialinf, only : atan
#endif
  implicit none

  private

  public  gripi_norm_potential, gripi_epsilon_one, gripi_epsilon_two
  public  gripi_epsilon_three,gripi_x_endinf
  public  gripi_efold_primitive, gripi_x_trajectory
  public  gripi_norm_deriv_potential, gripi_norm_deriv_second_potential
  public  gripi_alphamin, gripi_alphamax, gripi_x_epsonemin
  public  gripi_x_epstwozero,gripi_x_epsonezero


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
         (phi0**2*x**2*(6._kp+alpha*x*(-8._kp+3._kp*x))**2)

  end function gripi_epsilon_one


  !epsilon_two(x)
  function gripi_epsilon_two(x,alpha,phi0)    
    implicit none
    real(kp) :: gripi_epsilon_two
    real(kp), intent(in) :: x,alpha,phi0

    gripi_epsilon_two =(24._kp*(6._kp+alpha*x*(-16._kp+x*(3._kp+alpha* &
         (16._kp+3._kp*(-4._kp+x)*x)))))/(phi0**2*x**2* &
         (6._kp+alpha*x*(-8._kp+3._kp*x))**2)

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
    maxi = gripi_x_epsonemin(alpha)*(1._kp-epsilon(1._kp)) !Position of the first local minimum of epsilon1

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


  !this is integral[V(phi)/V'(phi) dphi]
  function gripi_efold_primitive(x,alpha,phi0)
    implicit none
    real(kp), intent(in) :: x,alpha,phi0
    real(kp) :: gripi_efold_primitive
    complex(kp) :: carg

    if (alpha.eq.0._kp) then 
       gripi_efold_primitive = x**2/4._kp

    elseif (alpha.eq.1._kp) then

       stop 'gripi_efold_primitive: you must use ripi for alpha=1!'

    else

       carg = (x-1._kp)/(sqrt(cmplx(1._kp/alpha-1._kp,0._kp,kp)))

       gripi_efold_primitive = phi0**2*real((5._kp-4._kp)/ &
            (12._kp*sqrt(cmplx(alpha*(1._kp-alpha),0._kp,kp)))* &
            atan(carg)+0.5_kp*x*(0.25_kp*x-1._kp/3._kp)+ &
            (1._kp/(8._kp*alpha)-1._kp/6._kp)*log(1._kp+alpha*x &
            *(x-2._kp)),kp)

    endif

  end function gripi_efold_primitive



  !returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function gripi_x_trajectory(bfold,xend,alpha,phi0)
    implicit none
    real(kp), intent(in) :: bfold, alpha, phi0, xend
    real(kp) :: gripi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gripiData

    mini = xend

    if (alpha .lt. 1._kp) then

       maxi = gripi_x_epsonemin(alpha)*100._kp 

    elseif (alpha.eq.1._kp) then

       stop 'gripi_x_trajectory: you must use ripi for alpha=1!'

    else

       maxi = gripi_x_epsonemin(alpha) !local maximum of the potential

    endif

    gripiData%real1 = alpha
    gripiData%real2 = phi0
    gripiData%real3 = -bfold + gripi_efold_primitive(xend,alpha,phi0)

    gripi_x_trajectory = zbrent(find_gripi_x_trajectory,mini,maxi,tolFind,gripiData)

  end function gripi_x_trajectory

  function find_gripi_x_trajectory(x,gripiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gripiData
    real(kp) :: find_gripi_x_trajectory
    real(kp) :: alpha,phi0,NplusNuend

    alpha= gripiData%real1
    phi0 = gripiData%real2
    NplusNuend = gripiData%real3

    find_gripi_x_trajectory = gripi_efold_primitive(x,alpha,phi0) - NplusNuend

  end function find_gripi_x_trajectory


  !Returns the prior alphamin(phi0) such that, when alpha<1, the minimum of epsilon1 is less than one
  function gripi_alphamin(phi0)    
    implicit none
    real(kp), intent(in) :: phi0   
    real(kp) :: gripi_alphamin
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gripiData

    mini = epsilon(1._kp)
    mini = 0.9_kp
    maxi = 1._kp*(1._kp-epsilon(1._kp))

    gripiData%real1 = phi0

    gripi_alphamin = zbrent(find_gripi_alphamin,mini,maxi,tolFind,gripiData)


  end function gripi_alphamin

  function find_gripi_alphamin(alpha,gripiData)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: gripiData
    real(kp) :: find_gripi_alphamin
    real(kp) :: phi0

    phi0= gripiData%real1

    find_gripi_alphamin = gripi_epsilon_one(gripi_x_epsonemin(alpha),alpha,phi0)-1._kp

  end function find_gripi_alphamin

  !Returns the prior alphamax(phi0,efold) such that, when alpha>1,
  !at least efold e-folds can be realized
  function gripi_alphamax(phi0,efold)    
    implicit none
    real(kp), intent(in) :: phi0,efold 
    real(kp) :: gripi_alphamax
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gripiData

    mini = 1._kp*(1._kp+epsilon(1._kp))
    mini = 1._kp*(1._kp+5._kp*epsilon(1._kp))
    maxi = 2._kp

    gripiData%real1 = phi0
    gripiData%real2 = efold

    gripi_alphamax = zbrent(find_gripi_alphamax,mini,maxi,tolFind,gripiData)


  end function gripi_alphamax

  function find_gripi_alphamax(alpha,gripiData)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: gripiData
    real(kp) :: find_gripi_alphamax
    real(kp) :: phi0,efold

    phi0= gripiData%real1
    efold= gripiData%real2

    find_gripi_alphamax = gripi_efold_primitive(gripi_x_epsonemin(alpha),alpha,phi0) &
         -gripi_efold_primitive(gripi_x_endinf(alpha,phi0),alpha,phi0) &
         -efold

  end function find_gripi_alphamax


end module gripisr
