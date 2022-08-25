!common slow-roll functions for string axion inflation
!
!V(phi) = M^4 {1 - cos(x) + alpha [ x sin(x) + (1/2) beta x^2 ]}
!
!x=phi/mu
!
!with no assumptions on alpha, beta and mu, but one must have: alpha beta > 0
!
module saiiicommon
  use infprec, only : kp, tolkp, toldp, transfert, pi
  use inftools, only : zbrent, easydverk

  implicit none

  
  private

  public saiii_norm_potential, saiii_norm_deriv_potential
  public saiii_norm_deriv_second_potential
  public saiii_epsilon_one, saiii_epsilon_two
  public saiii_epsilon_three, saiii_x_epsoneunity
  public saiii_efold_primitive, find_saiii_x_trajectory
  public saiii_check_params, saiii_check_minima, saiii_check_zero_negpot
  public saiii_alpha_potneg, saiii_x_derivpotzero, saiii_x_potmax
  public saiii_x_potzero, saiii_alpha_one, saiii_alpha_two, saiii_alpha_three

!the solution of sqrt(beta^2-1)[ pi + arccos(1/beta) ] = 1  
  real(kp), parameter :: beta0 = 1.0417370933450080562297020810858930632202316681725_kp

!the solution of sqrt(beta^2-1)[ 2pi + arccos(-1/beta) ] = 1  
  real(kp), parameter :: beta1 = -1.0119940554981193879678434388031569608462871828152_kp

!-2 sinc(x) at its first minimum
  real(kp), parameter :: beta2 = 0.4344672564224433148165586511249414684460898308712_kp

!-2 sinc(x) at its first non-vanishing maximum
  real(kp), parameter :: beta3 = -0.256749107051798273386061556467739004321436872134_kp

  real(kp), parameter :: saiiiXBig = log(epsilon(1._8)*huge(1._8))
  real(kp), parameter :: saiiiAlphaBig = 1._kp/toldp
  
  public beta0, beta1, beta2, beta3, saiiiXBig
  
contains


!returns true for sane parameters alpha beta > 0
  function saiii_check_params(alpha,beta,mu)
    implicit none
    logical :: saiii_check_params
    real(kp), intent(in) :: alpha,beta,mu

    saiii_check_params = (alpha*beta.ge.0._kp)

  end function saiii_check_params

  
  
  function saiii_norm_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii_norm_potential
    real(kp), intent(in) :: x,alpha,beta,mu

    saiii_norm_potential = 1._kp - cos(x) &
         + alpha*(x*sin(x) + beta*x*x/2._kp)
    
  end function saiii_norm_potential


  
!derivative with respect to x (not phi!)  
  function saiii_norm_deriv_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta,mu

    saiii_norm_deriv_potential = (1._kp + alpha)*sin(x) &
         + alpha*x*(cos(x) + beta)


  end function saiii_norm_deriv_potential

  

!second derivative with respect to x (not phi!)    
  function saiii_norm_deriv_second_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii_norm_deriv_second_potential = (1._kp + 2._kp*alpha)*cos(x) &
         + alpha*(beta - x*sin(x))

  end function saiii_norm_deriv_second_potential



  
  function saiii_epsilon_one(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii_epsilon_one
    real(kp), intent(in) :: x,alpha,beta,mu
    real(kp) :: sinx,cosx
    
    sinx = sin(x)
    cosx = cos(x)

    saiii_epsilon_one = ((1._kp + alpha)*sinx &
         + alpha*x*(cosx + beta))**2/(1._kp - cosx &
         + alpha*(x*sinx + beta*x*x/2._kp))**2/2._kp/mu/mu
    
   
    
  end function saiii_epsilon_one
 
  
  
  function saiii_epsilon_two(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii_epsilon_two
    real(kp), intent(in) :: x,alpha,beta,mu
    real(kp) :: sinx,cosx
    
    sinx = sin(x)
    cosx = cos(x)

    saiii_epsilon_two= (4._kp*(2._kp + alpha*(4._kp - 2._kp*beta &
         + alpha*(1._kp + x**2*(2._kp + beta**2))) &
         + (-2._kp - 4._kp*alpha + alpha*(2._kp + x**2*(-1._kp + 2._kp*alpha))*beta) &
         * cosx - alpha**2*cos(2*x) + x*alpha*(2._kp + (4._kp + (2._kp + x**2)*alpha)*beta) &
         * sinx)) / ((2._kp + x**2*alpha*beta - 2._kp*cosx + 2._kp*x*alpha*sinx)**2)/mu/mu

  end function saiii_epsilon_two


  
  function saiii_epsilon_three(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii_epsilon_three
    real(kp), intent(in) :: x,alpha,beta,mu
    real(kp) :: sinx, cosx, cos2x, sin2x

    cosx = cos(x)
    sinx = sin(x)
    sin2x = 2._kp*sinx*cosx
    cos2x = cosx*cosx - sinx*sinx

    saiii_epsilon_three = -(((2._kp*x*alpha*(beta + cosx) + 2._kp*(1._kp + alpha)*sinx) &
         *(-4._kp*(x*alpha*beta + x*alpha*cosx + (1._kp + alpha)*sinx)*(2._kp + alpha &
         *(4._kp - 2._kp*beta + alpha*(1._kp + x**2*(2._kp + beta**2))) + (-2._kp - 4._kp*alpha &
         + alpha*(2._kp + x**2*(-1._kp + 2._kp*alpha))*beta)*cosx - alpha**2*cos2x &
         +  x*alpha*(2._kp + (4._kp + (2._kp + x**2)*alpha)*beta)*sinx) &
         + (2._kp + x**2*alpha*beta - 2._kp*cosx + 2._kp*x*alpha*sinx)*(x*alpha &
         * (2._kp + (2._kp + (6._kp + x**2)*alpha)*beta)*cosx + (2._kp + (2._kp + x**2) &
         *alpha**2*beta + alpha*(6._kp + (2._kp + x**2)*beta))*sinx + 2._kp*alpha**2 &
         *(x*(2._kp + beta**2) + sin2x)))) &
         /((2._kp + x**2*alpha*beta - 2._kp*cosx + 2._kp*x*alpha*sinx)**2 &
         *(2._kp + alpha*(4._kp - 2._kp*beta + alpha*(1._kp + x**2*(2._kp + beta**2))) &
         +(-2._kp - 4._kp*alpha + alpha*(2._kp + x**2*(-1._kp + 2*alpha))*beta)*cosx &
         - alpha**2*cos2x + x*alpha*(2._kp + (4._kp + (2._kp + x**2)*alpha)*beta)*sinx)))/mu/mu
   

  end function saiii_epsilon_three




  function saiii_alpha_one(beta)
    implicit none
    real(kp) :: saiii_alpha_one
    real(kp), intent(in) :: beta

    if (abs(beta).lt.1._kp) then
       stop 'saiii_alpha_one: exist only for beta^2 >= 1'
    endif
    
    saiii_alpha_one = -1._kp/(sqrt(beta*beta-1._kp)*(pi+acos(1._kp/beta)) + 1._kp)
    
  end function saiii_alpha_one

  

  
  function saiii_alpha_two(beta)
    implicit none
    real(kp) :: saiii_alpha_two
    real(kp), intent(in) :: beta

    if (abs(beta).lt.1._kp) then
       stop 'saiii_alpha_two: exist only for beta^2 >= 1'
    endif
    
    saiii_alpha_two = 1._kp/(sqrt(beta*beta-1._kp)*(pi+acos(1._kp/beta)) - 1._kp)
    
  end function saiii_alpha_two


  

  function saiii_alpha_three(beta)
    implicit none
    real(kp) :: saiii_alpha_three
    real(kp), intent(in) :: beta

    if (abs(beta).lt.1._kp) then
       stop 'saiii_alpha_three: exist only for beta^2 >= 1'
    endif

    saiii_alpha_three = 1._kp/(sqrt(beta*beta-1._kp)*(2._kp*pi+acos(-1._kp/beta)) - 1._kp)

  end function saiii_alpha_three
  

  

!returns true is the potential has minima, false otherwise
  function saiii_check_minima(alpha,beta,mu)
    implicit none
    logical :: saiii_check_minima
    real(kp), intent(in) :: alpha, beta, mu

    if (.not.saiii_check_params(alpha,beta,mu)) then
       stop 'saiii_check_minima: alpha.beta < 0!'
    endif

    
    if ((beta.ge.-1._kp).and.(beta.le.beta0)) then

       saiii_check_minima = .true.

    elseif (beta.gt.beta0) then

       saiii_check_minima = .not.(alpha.gt.saiii_alpha_two(beta))

    elseif ((beta.lt.-1._kp).and.(beta.gt.beta1)) then

       saiii_check_minima = .not. &
            ((alpha.lt.saiii_alpha_one(beta)).and.(alpha.gt.saiii_alpha_three(beta)))

    elseif (beta.le.beta1) then

       saiii_check_minima = .not. (alpha.lt.saiii_alpha_one(beta))

    else

       stop 'saiii_check_minima: this should not have happened!'

    endif
       
  end function saiii_check_minima



  
!returns true is the potential is negative in zero
  function saiii_check_zero_negpot(alpha,beta,mu)
    implicit none
    logical :: saiii_check_zero_negpot
    real(kp), intent(in) :: alpha,beta,mu

    if (.not.saiii_check_params(alpha,beta,mu)) then
       stop 'saiii_check_negpot: alpha.beta < 0!'
    endif
    
    if (beta.gt.-2._kp) then
       saiii_check_zero_negpot = (alpha.lt.-1._kp/(2._kp+beta))
    else
       saiii_check_zero_negpot = .false.
    endif
       
  end function saiii_check_zero_negpot
  
  

!return the minimal value of alpha > 0 and maximal value of
!alpha < 0 such that the potential becomes negative in a domain
!non-connected to x=0. Such a domain exists only for beta3 < beta < beta2
  function saiii_alpha_potneg(beta)
    implicit none
    real(kp) :: saiii_alpha_potneg
    real(kp), intent(in) :: beta

    real(kp), save :: stobeta = huge(1._kp)
    real(kp), save :: stoalpha = huge(1._kp)
!$omp threadprivate(stobeta,stoalpha)
    
    real(kp) :: mini, maxi
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: saiiiData
    
    if ((beta.gt.beta2).or.(beta.lt.beta3)) then
       stop 'saiii_alpha_potneg: beta out of [beta3,beta2], second minimum of V positive!'
    elseif ((beta.eq.beta2).or.(beta.eq.beta3)) then
       stop 'saiii_alpha_potneg: beta=beta2/3, second minimum of V is exactly null!'
    endif
    
    if (beta.eq.stobeta) then
       saiii_alpha_potneg = stoalpha
       return
    else
       stobeta = beta
    endif

    
    if (beta.gt.0._kp) then

       mini = tolFind
       maxi = saiiiAlphaBig
       
       saiiiData%real2 = beta
       stoalpha = zbrent(find_saiii_alphapos_potneg,mini,maxi,tolFind,saiiiData)

    elseif (beta.lt.0._kp) then

       mini = -saiiiAlphaBig
       maxi = -tolFind

       saiiiData%real2 = beta
       stoalpha = zbrent(find_saiii_alphaneg_potneg,mini,maxi,tolFind,saiiiData)

    else

       stop 'saiii_alpha_potneg: beta = 0, use SAII module'

    endif

       
    saiii_alpha_potneg = stoalpha
    
  end function saiii_alpha_potneg
  
     
  function find_saiii_alphapos_potneg(alpha,saiiiData)
    implicit none
    real(kp) :: find_saiii_alphapos_potneg
    real(kp), intent(in) :: alpha
    type(transfert), optional, intent(inout) :: saiiiData
    real(kp), parameter :: tolFind = tolkp
    real(kp) :: beta,xmin,xmax,xVmin

    saiiiData%real1 = alpha
    beta = saiiiData%real2
    xmin = pi + acos(beta)
    xmax = 2._kp*pi
    xVmin = zbrent(find_saiii_x_derivpotzero,xmin,xmax,tolFind,saiiiData)
        
    find_saiii_alphapos_potneg = saiii_norm_potential(xVmin,alpha,beta,1._kp)

  end function find_saiii_alphapos_potneg
       

  function find_saiii_alphaneg_potneg(alpha,saiiiData)
    implicit none
    real(kp) :: find_saiii_alphaneg_potneg
    real(kp), intent(in) :: alpha
    type(transfert), optional, intent(inout) :: saiiiData
    real(kp), parameter :: tolFind = tolkp
    real(kp) :: beta,xmin,xmax,xVmin

    beta = saiiiData%real2
    saiiiData%real1 = alpha
    
    if (alpha.ge.-1._kp) then
       xmin = 2._kp*pi 
       xmax = 2._kp*pi + acos(-beta) 
    else
       xmin = 2._kp*pi + acos(-beta)
       xmax = 3._kp*pi 
    endif       
    
    xVmin = zbrent(find_saiii_x_derivpotzero,xmin,xmax,tolFind,saiiiData)

    find_saiii_alphaneg_potneg = saiii_norm_potential(xVmin,alpha,beta,1._kp)

    
  end function find_saiii_alphaneg_potneg
  


!non vanishing field values at which the potential is positive and extremal
  function saiii_x_derivpotzero(alpha,beta,mu)
    implicit none
    real(kp), dimension(2) :: saiii_x_derivpotzero
    real(kp), intent(in) :: alpha,beta,mu

    real(kp) :: mini, maxi
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: saiiiData

    real(kp), save :: stoalpha = huge(1._kp)
    real(kp), save :: stobeta = huge(1._kp)
    real(kp), dimension(2), save :: stox = huge(1._kp)
!$omp threadprivate(stoalpha,stobeta,stox)    


    if (.not.saiii_check_minima(alpha,beta,mu)) then
       stop 'saiii_x_derivpotzero: V has no minima for these values of alpha/beta!'
    endif
    

    if ((alpha.eq.stoalpha).and.(beta.eq.stobeta)) then
       saiii_x_derivpotzero = stox
       return
    else
       stoalpha = alpha
       stobeta = beta
    endif

    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    
    if ((beta.ge.1._kp).or.((beta.le.-1._kp).and.(alpha.ge.-1._kp))) then

       mini = pi
       maxi = pi + acos(1._kp/beta)
       stox(1) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)       

       mini = pi + acos(1._kp/beta)
       maxi = 2._kp*pi
       stox(2) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)  
              
    elseif (abs(beta).lt.1._kp) then

       if (alpha.ge.0._kp) then

          mini = acos(-beta)
          maxi = pi
          stox(1) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)

          mini = pi + acos(beta)
          maxi = 2._kp*pi
          stox(2) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)

       elseif ((alpha.lt.-1._kp/(beta+2._kp)).and.(alpha.lt.-1._kp)) then

!the first maximum, positive, we keep it          
          mini = pi + acos(-beta)
          maxi = 2._kp*pi
          stox(1) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)

!the second minimum, positive, we keep it
          mini = 2._kp*pi + acos(-beta)
          maxi = 3._kp*pi
          stox(2) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)

       elseif ((alpha.lt.-1._kp/(beta+2._kp)).and.(alpha.ge.-1._kp)) then
          
!the first maximum, positive, we keep it          
          mini = pi
          maxi = pi + acos(-beta)
          stox(1) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)

!the second minimum, positive, we keep it
          mini = 2._kp*pi
          maxi = 2._kp*pi + acos(-beta)
          stox(2) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)

          
       else

          mini = pi
          maxi = pi + acos(beta)
          stox(1) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)

          mini = 2._kp*pi
          maxi = 2._kp*pi + acos(-beta)
          stox(2) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)

       endif
                    
    elseif ((beta.gt.beta1).and.(alpha.le.saiii_alpha_three(beta))) then

       mini = 2._kp*pi
       maxi = 2._kp*pi + acos(-1._kp/beta)
       stox(1) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)

       mini = 2._kp*pi + acos(-1._kp/beta)
       maxi = 3._kp*pi
       stox(2) = zbrent(find_saiii_x_derivpotzero,mini,maxi,tolFind,saiiiData)
       
    else

       write(*,*)'alpha= beta= ',alpha,beta
       stop ' saiii_x_derivpotzero: this should not happen!'

    endif
                  
    saiii_x_derivpotzero = stox
    
  end function saiii_x_derivpotzero

  
  function find_saiii_x_derivpotzero(x,saiiiData)
    implicit none
    real(kp) :: find_saiii_x_derivpotzero
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiiData
    real(kp) :: alpha,beta

    alpha = saiiiData%real1
    beta = saiiiData%real2
        
    find_saiii_x_derivpotzero = saiii_norm_deriv_potential(x,alpha,beta,1._kp)
    
  end function find_saiii_x_derivpotzero


    
  function saiii_x_potmax(alpha,beta,mu)
    implicit none
    real(kp), intent(in) :: alpha,beta,mu
    real(kp) :: saiii_x_potmax
    real(kp), dimension(2) :: xdVzero

    xdVzero = saiii_x_derivpotzero(alpha,beta,mu)

    saiii_x_potmax = xdVzero(1)    
    
  end function saiii_x_potmax
  


  function saiii_x_potzero(alpha,beta,mu)
    implicit none
    real(kp), dimension(3) :: saiii_x_potzero
    real(kp), intent(in) :: alpha,beta,mu

    real(kp) :: mini, maxi
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: saiiiData

    real(kp), dimension(2) :: xdVzero
    real(kp) :: alpthresh
    
    logical :: negInZero, posInZero
    logical :: hasMinima, hasNotMinimum
    
    negInZero = saiii_check_zero_negpot(alpha,beta,mu)
    posInZero = .not.negInZero
    hasMinima = saiii_check_minima(alpha,beta,mu)
    hasNotMinimum = .not.hasMinima

    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    
    if (hasNotMinimum) then
    
!monotonous increasing function of x vanishing in x=0 only    
       if (posInZero) then
          saiii_x_potzero(1:3) = 0._kp

!only one zero, and monotonous increasing function of x at higher x values           
       else

          mini = tolFind
          maxi = 2._kp*pi
          saiii_x_potzero(1:3) = zbrent(find_saiii_x_potzero,mini,maxi,tolFind,saiiiData)

       endif


!has minima, that may or maynot be negative       
    else

       xdVzero = saiii_x_derivpotzero(alpha,beta,mu)

       if (negInZero) then

          mini = 0._kp + epsilon(1._kp)
          maxi = xdVzero(1)
          saiii_x_potzero(1) = zbrent(find_saiii_x_potzero,mini,maxi,tolFind,saiiiData)
          
!positive in zero
       else

          saiii_x_potzero(1) = 0._kp

       endif
          
       if ((beta.le.beta2).and.(beta.ge.beta3)) then

!to test saiii_alpha_potneg
!          alpthresh = saiii_alpha_potneg(beta)
!          if (((alpha.gt.0._kp).and.(alpha.ge.alpthresh)).or. &
!               ((alpha.lt.0._kp).and.(alpha.le.alpthresh))) then             
!             if (saiii_norm_potential(xdVzero(2),alpha,beta,mu).gt.0._kp) then
!                print *,'alpha,beta= ',alpha,beta
!                print *,'V=',saiii_norm_potential(xdVzero(2),alpha,beta,mu)
!                print *,'x=',xdVzero(2)
!                print *,'dV=',saiii_norm_deriv_potential(xdVzero(2),alpha,beta,mu)
!                stop 'saiii_x_potzero: mistake inside!'
!             endif

          
          if (saiii_norm_potential(xdVzero(2),alpha,beta,mu).lt.0._kp) then
          
             mini = xdVzero(1)
             maxi = xdVzero(2)
             saiii_x_potzero(2) = zbrent(find_saiii_x_potzero,mini,maxi,tolFind,saiiiData)

             mini = xdVzero(2)
             maxi = pi + xdVzero(2)
             saiii_x_potzero(3) = zbrent(find_saiii_x_potzero,mini,maxi,tolFind,saiiiData)

          elseif (saiii_norm_potential(xdVzero(2),alpha,beta,mu).lt.epsilon(1._kp)) then

             saiii_x_potzero(2:3) = xdVzero(2) - tolFind

          else
             
             saiii_x_potzero(2:3) = saiii_x_potzero(1)

          endif

       else

          saiii_x_potzero(2:3) = saiii_x_potzero(1)

       endif


    endif
    
  end function saiii_x_potzero



  function find_saiii_x_potzero(x,saiiiData)
    implicit none
    real(kp) :: find_saiii_x_potzero
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiiData
    real(kp) :: alpha,beta

    alpha = saiiiData%real1
    beta = saiiiData%real2
        
    find_saiii_x_potzero = saiii_norm_potential(x,alpha,beta,1._kp)

  end function find_saiii_x_potzero




  
!returns the two roots of epsilon1(x)=1, in the domain for which V>0
!and around its first maximum. If V has not any extrema, return the
!lowest root
  function saiii_x_epsoneunity(alpha,beta,mu)
    implicit none
    real(kp), dimension(2) :: saiii_x_epsoneunity
    real(kp), intent(in) :: alpha,beta,mu

    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: saiiiData
    
    real(kp) :: mini, maxi, xpotmax, xguess
    logical :: negInZero, hasMinima, newtSuccess

    real(kp), dimension(3) :: xpotzero    
    
    real(kp), save, dimension(2) :: xeps
    real(kp), save :: stoalpha = huge(1._kp)
    real(kp), save :: stobeta = huge(1._kp)
    real(kp), save :: stomu = huge(1._kp)
!$omp threadprivate(stoalpha,stobeta,stomu,xeps)

    if ((alpha.eq.stoalpha).and.(beta.eq.stobeta).and.(mu.eq.stomu)) then

       saiii_x_epsoneunity = xeps
       return

    else

       stoalpha = alpha
       stobeta = beta
       stomu = mu

    endif

    
    negInZero = saiii_check_zero_negpot(alpha,beta,mu)
    hasMinima = saiii_check_minima(alpha,beta,mu)
    
    xpotzero = saiii_x_potzero(alpha,beta,mu)
    
    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    saiiiData%real3 = mu


    if (.not.hasMinima) then
!brent should give the lowest root, but oscillating non-extremal
!potential may have a lot of eps-1 roots at higher value. Inflation
!would actually ends and starts a few times before reaching the
!computed root here and slow-roll violation should be expected before.
       mini = xpotzero(1) + tolkp
       maxi = saiiiXBig
       saiiiData%real4 = +1._kp
       xeps(1:2)=zbrent(find_saiii_x_epsoneunity,mini,maxi,tolFind,saiiiData)

    else

       xpotmax = saiii_x_potmax(alpha,beta,mu)
       
       mini = xpotzero(1) + epsilon(1._kp)
       maxi = xpotmax
       saiiiData%real4 = +1._kp
       xeps(1) = zbrent(find_saiii_x_epsoneunity,mini,maxi,tolFind,saiiiData)

       if (xpotzero(2).ne.xpotzero(1)) then

          mini = xpotmax
          maxi = xpotzero(2) - epsilon(1._kp)
          saiiiData%real4 = -1._kp
          xeps(2) = zbrent(find_saiii_x_epsoneunity,mini,maxi,tolFind,saiiiData)

       else

          xeps(2) = xeps(1)

       endif


    endif
    
    saiii_x_epsoneunity = xeps

  end function saiii_x_epsoneunity



  function eval_saiii_x_epsoneunity(x,saiiiData)
    implicit none
    real(kp) :: eval_saiii_x_epsoneunity
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiiData

    real(kp) :: alpha,beta,mu

    alpha = saiiiData%real1
    beta = saiiiData%real2
    mu = saiiiData%real3

    eval_saiii_x_epsoneunity = log(saiii_epsilon_one(x,alpha,beta,mu))
    
  end function eval_saiii_x_epsoneunity

  
  function deriv_saiii_x_epsoneunity(x,saiiiData)
    implicit none
    real(kp) :: deriv_saiii_x_epsoneunity
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiiData

    real(kp) :: alpha,beta,mu,V,dV,ddV

    alpha = saiiiData%real1
    beta = saiiiData%real2
    mu = saiiiData%real3

    V = saiii_norm_potential(x,alpha,beta,mu)
    dV = saiii_norm_deriv_potential(x,alpha,beta,mu)
    ddV = saiii_norm_deriv_second_potential(x,alpha,beta,mu)

    if (V*dV.eq.0._kp) stop 'deriv_saiii_x_epsoneunity: VdV = 0!'
    
    deriv_saiii_x_epsoneunity = 2._kp*(V*ddV - dV*dV)/(V*dV)
    
    
  end function deriv_saiii_x_epsoneunity


  
  function find_saiii_x_epsoneunity(x,saiiiData)
    implicit none
    real(kp) :: find_saiii_x_epsoneunity
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiiData
    real(kp) :: alpha,beta,mu
    real(kp), parameter :: sqr2 = sqrt(2._kp)
    real(kp) :: sinx, cosx, sqrmu, pm

    
    alpha = saiiiData%real1
    beta = saiiiData%real2
    mu = saiiiData%real3
    pm = saiiiData%real4
    
    sinx = sin(x)
    cosx = cos(x)
    sqrmu = sqrt(mu)
    
    find_saiii_x_epsoneunity = pm*sqrmu*sqr2*(1._kp - cosx + alpha*x*sinx +0.5_kp*alpha*beta*x*x) &
         - ((1._kp+alpha)*sinx + alpha*x*(beta+cosx))/sqrmu
    
  end function find_saiii_x_epsoneunity



  
!this is integral[V(phi)/V'(phi) dphi]
  function saiii_efold_primitive(x,alpha,beta,mu)
    implicit none
    real(kp), intent(in) :: x, alpha,beta,mu
    real(kp) :: saiii_efold_primitive

    type(transfert) :: saiiiData
!too long to integrate in QUADPREC
    real(kp), parameter :: tolInt = max(tolkp,toldp)
    integer, parameter :: neq = 1

    real(kp), dimension(3) :: xpotzero
    real(kp) :: xpotmax
    real(kp) :: xvar
    real(kp), dimension(neq) :: yvar

    logical :: hasMinima
        
    hasMinima = saiii_check_minima(alpha,beta,mu)

    xpotzero = saiii_x_potzero(alpha,beta,mu)
    
    yvar(1) = 0._kp
    
    if (.not.hasMinima) then

       xvar = xpotzero(1)                     

    else

       xpotmax = saiii_x_potmax(alpha,beta,mu)

       if (x.lt.xpotmax) then
          xvar = xpotzero(1)
       else
          xvar = xpotzero(2)
          if (xpotzero(2).eq.xpotzero(1)) then
             write(*,*)'x= xpotzero(1:3)= '
             write(*,*)'saiii_efold_primitive: x > xVmax with Vmin>0!'
          endif
       endif
          
       
    endif
       
    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    
    call easydverk(neq,find_saiii_efold_primitive,xvar,yvar,x,tolInt,saiiiData)

    saiii_efold_primitive = yvar(1) * mu * mu
    
  end function saiii_efold_primitive

  
  subroutine find_saiii_efold_primitive(n,x,y,yprime,saiiiData)
    implicit none          
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: saiiiData
    real(kp) :: alpha,beta,cosx,sinx
!the expansion is accurate up to order 3, and we take a factor of ten
!margin
    real(kp), parameter :: xtaylor = 10._kp*epsilon(1._kp)**(1._kp/3._kp)

    alpha = saiiiData%real1
    beta = saiiiData%real2
    
    cosx = cos(x)
    sinx = sin(x)
    
    if (abs(x).lt.xtaylor) then
       if (alpha.eq.-1._kp/(beta+2._kp)) then
          if (beta.eq.2._kp) then
             yprime(1) = x/6._kp
          else
             yprime(1) = 0.25_kp*x
          endif
       else
          yprime(1) = 0.5_kp*x
       endif
    else
       yprime(1) = (1._kp - cosx + alpha * x * sinx + 0.5_kp*alpha*beta*x*x) &
            / ( (1._kp+alpha)*sinx + alpha * x *(cosx + beta) )
    endif

  end subroutine find_saiii_efold_primitive



  function find_saiii_x_trajectory(x,saiiiData)
    implicit none
    real(kp) :: find_saiii_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiiData

    real(kp) :: alpha, beta, mu, NplusNuend

    alpha = saiiiData%real1
    beta = saiiiData%real2
    mu = saiiiData%real3    
    NplusNuend = saiiiData%real4

    find_saiii_x_trajectory = saiii_efold_primitive(x,alpha,beta,mu) - NplusNuend
    
  end function find_saiii_x_trajectory

  
end module saiiicommon
