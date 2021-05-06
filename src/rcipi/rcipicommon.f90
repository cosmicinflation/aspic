!common slow-roll function for radiatively corrected plateau inflation
!
!V(phi) = M^4 x^p [1 + alpha ln(x) + beta ln(x)^2]
!
!x=phi/Mp
!
!with no assumptions on p, alpha, beta
!
module rcipicommon
  use infprec, only : kp, tolkp, transfert
  use inftools, only : selectsort, zbrent

  implicit none

  
  private

  public rcipi_norm_potential, rcipi_norm_deriv_potential
  public rcipi_norm_deriv_second_potential
  public rcipi_epsilon_one, rcipi_epsilon_two
  public rcipi_epsilon_three, rcipi_x_epstwozero, rcipi_x_epsoneunity
  public rcipi_efold_primitive, find_rcipi_x_trajectory
  public rcipi_check_derivpotzero, rcipi_check_potzero
  public rcipi_x_potzero, rcipi_x_derivpotzero
  public rcipi_numacc_x_potbig, rcipi_alpha_zero

  real(kp), parameter :: lnzero = -huge(1._kp)

  
contains

  
  function rcipi_norm_potential(x,p,alpha,beta)
    implicit none
    real(kp) :: rcipi_norm_potential
    real(kp), intent(in) :: x,alpha,beta,p
    real(kp) :: lnx

    lnx = 0.5_kp*log(x*x)
    
    rcipi_norm_potential = x**p*(1._kp + alpha*lnx + beta *lnx**2)
    
  end function rcipi_norm_potential


  
  function rcipi_norm_deriv_potential(x,p,alpha,beta)
    implicit none
    real(kp) :: rcipi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta,p
    real(kp) :: lnx

    lnx = 0.5_kp*log(x*x)
    
    rcipi_norm_deriv_potential = x**(p-1._kp)* &
         (p + alpha + (alpha*p + 2*beta)*lnx + p*beta*lnx**2)

  end function rcipi_norm_deriv_potential


  
  function rcipi_norm_deriv_second_potential(x,p,alpha,beta)
    implicit none
    real(kp) :: rcipi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p,alpha,beta
    real(kp) :: lnx

    lnx = 0.5_kp*log(x*x)
    
    rcipi_norm_deriv_second_potential = x**(p-2._kp) &
         * (p*p - alpha + p*(2*alpha - 1._kp) + 2*beta &
         + (-p*alpha + p*p*alpha - 2*beta + 4*p*beta)*lnx &
         + (p-1._kp)*p*beta*lnx**2)

  end function rcipi_norm_deriv_second_potential


  
  function rcipi_epsilon_one(x,p,alpha,beta)
    implicit none
    real(kp) :: rcipi_epsilon_one
    real(kp), intent(in) :: x,p,alpha,beta
    real(kp) :: lnx, denom

    lnx = 0.5_kp*log(x*x)

    denom = 1 + lnx*(alpha + beta*lnx)

!NaN regularizator for extreme cases
    if (denom.eq.0._kp) denom = tiny(1._kp)
        
    rcipi_epsilon_one = (p + (alpha + 2*beta*lnx) &
         / denom)**2 /(2._kp*x**2)
    
  end function rcipi_epsilon_one
 
  
  
  function rcipi_epsilon_two(x,p,alpha,beta)
    implicit none
    real(kp) :: rcipi_epsilon_two
    real(kp), intent(in) :: x,p,alpha,beta
    real(kp) :: lnx

    lnx = 0.5_kp*log(x*x)

    rcipi_epsilon_two = (2*(p + (alpha**2 - 4*beta) &
         /(1._kp + lnx*(alpha + beta*lnx))**2 &
         + (alpha + 2*beta + 2*beta*lnx) &
         /(1._kp + lnx*(alpha + beta*lnx))))/x**2
    
  end function rcipi_epsilon_two


  
  function rcipi_epsilon_three(x,p,alpha,beta)
    implicit none
    real(kp) :: rcipi_epsilon_three
    real(kp), intent(in) :: x,p,alpha,beta
    real(kp) :: lnx

    lnx = 0.5_kp*log(x*x)
    
    rcipi_epsilon_three = -(((p + (alpha + 2*beta*lnx) &
         /(1._kp + lnx*(alpha + beta*lnx)))*(-2*p &
         - (2*(alpha**2 - 4*beta)*(alpha + 2*beta*lnx)) &
         /(1._kp + lnx*(alpha + beta*lnx))**3 &
         - (3*alpha**2 + 2*(-6 + alpha)*beta + 4*beta**2*lnx) &
         / (1 + lnx*(alpha + beta*lnx))**2 - (2*(alpha + 3*beta &
         + 2*beta*lnx))/(1._kp + lnx*(alpha + beta*lnx)))) &
         /(x**2*(p + (alpha**2 - 4*beta)/(1._kp + lnx*(alpha &
         + beta*lnx))**2 + (alpha + 2*beta + 2*beta*lnx) &
         /(1._kp + lnx*(alpha + beta*lnx)))))

  end function rcipi_epsilon_three

  
!true if the potential vanishes (and becomes negative in some regions)  
  function rcipi_check_potzero(alpha,beta)
    implicit none
    logical :: rcipi_check_potzero
    real(kp), intent(in) :: alpha, beta
    
    rcipi_check_potzero = ((alpha**2 - 4*beta).ge.0._kp)
    
  end function rcipi_check_potzero

  
  
!non-vanishing field values at which the potential vanishes
  function rcipi_x_potzero(alpha,beta)
    implicit none
    real(kp), dimension(2) :: rcipi_x_potzero
    real(kp), intent(in) :: alpha,beta

    rcipi_x_potzero = exp(rcipi_lnx_potzero(alpha,beta))
    
  end function rcipi_x_potzero


  function rcipi_lnx_potzero(alpha,beta)
    implicit none
    real(kp), dimension(2) :: rcipi_lnx_potzero
    real(kp), intent(in) :: alpha,beta

    if (.not.rcipi_check_potzero(alpha,beta)) then
       stop 'rcipi_x_potzero: V > 0'
    endif

    if ((beta.eq.0._kp).and.(alpha.ne.0._kp)) then
       rcipi_lnx_potzero(1) = 1._kp/alpha
       rcipi_lnx_potzero(2) = 1._kp/alpha
       return
    endif
           
    if (beta.gt.0._kp) then
       rcipi_lnx_potzero(1) = 0.5_kp*(-alpha - sqrt(alpha*alpha - 4*beta))/beta
       rcipi_lnx_potzero(2) = 0.5_kp*(-alpha + sqrt(alpha*alpha - 4*beta))/beta
    else
       rcipi_lnx_potzero(2) = 0.5_kp*(-alpha - sqrt(alpha*alpha - 4*beta))/beta
       rcipi_lnx_potzero(1) = 0.5_kp*(-alpha + sqrt(alpha*alpha - 4*beta))/beta
    endif
       
       
  end function rcipi_lnx_potzero


!For beta <=p^2, the value of alpha for which the inflection point
!appears
  function rcipi_alpha_zero(p,beta)
    implicit none
    real(kp) :: rcipi_alpha_zero
    real(kp), intent(in) :: p,beta

    if (beta.gt.p*p) then
       stop 'rcipi_alpha_zero: beta > p^2!'
    endif
    
    rcipi_alpha_zero = 2._kp*sqrt(beta)*sqrt(1._kp - beta/p/p)
    
  end function rcipi_alpha_zero


  
!true if the potential admits local extrema
  function rcipi_check_derivpotzero(p,alpha,beta)
    implicit none
    logical :: rcipi_check_derivpotzero
    real(kp), intent(in) :: p,alpha, beta

    if (p.eq.0._kp) then
       rcipi_check_derivpotzero = .true.
       return
    endif
    
    rcipi_check_derivpotzero = (4._kp*beta*beta + p*p*(alpha*alpha-4._kp*beta)).ge.0._kp
    
  end function rcipi_check_derivpotzero
  

!non vanishing field values at which the potential is extremal
  function rcipi_x_derivpotzero(p,alpha,beta)
    implicit none
    real(kp), dimension(2) :: rcipi_x_derivpotzero
    real(kp), intent(in) :: p,alpha,beta

    rcipi_x_derivpotzero = exp(rcipi_lnx_derivpotzero(p,alpha,beta))
    
  end function rcipi_x_derivpotzero
  

  
  function rcipi_lnx_derivpotzero(p,alpha,beta)
    implicit none
    real(kp), dimension(2) :: rcipi_lnx_derivpotzero
    real(kp), intent(in) :: p,alpha,beta

    real(kp) :: delta
    
    if (.not.rcipi_check_derivpotzero(p,alpha,beta)) then
       stop 'rcipi_x_potmax: dV/dx > 0'
    endif

    if (p.eq.0._kp) then
       rcipi_lnx_derivpotzero = -0.5_kp * alpha/beta
       return
    endif

    delta = (2*beta/p+alpha)**2 - 4*beta*(1._kp+alpha/p)
    
    rcipi_lnx_derivpotzero(1) = 0.5_kp*(-2*beta/p-alpha &
         - sqrt(delta))/beta
    rcipi_lnx_derivpotzero(2) = 0.5_kp*(-2*beta/p-alpha &
         + sqrt(delta))/beta

    
    call selectsort(rcipi_lnx_derivpotzero)


  end function rcipi_lnx_derivpotzero

  

!return the non-vanishing field values at which epsilon1 is extremal
!(epsilon2 vanishes). At most two, zero values should be ignored.
  function rcipi_x_epstwozero(p,alpha,beta)
    use inftools, only : quarticroots
    implicit none
    real(kp), dimension(2) :: rcipi_x_epstwozero
    real(kp), intent(in) :: p,alpha,beta

    integer :: i
    real(kp), dimension(2) :: lnxeps
    
    lnxeps = rcipi_lnx_epstwozero(p,alpha,beta)

    do i=1,2
       if (lnxeps(i).eq.lnzero) then
          rcipi_x_epstwozero(i) = 0._kp
       else
          rcipi_x_epstwozero(i) = exp(lnxeps(i))
       endif
    enddo
    
  end function rcipi_x_epstwozero


  function rcipi_lnx_epstwozero(p,alpha,beta)
    use inftools, only : quarticroots
    implicit none
    real(kp), dimension(2) :: rcipi_lnx_epstwozero
    real(kp), intent(in) :: p,alpha,beta

    real(kp), dimension(5) :: coefficients
    complex(kp), dimension(4) :: croots

    integer :: i,j

    rcipi_lnx_epstwozero = lnzero
    
    coefficients(1) = p + alpha + alpha*alpha - 2*beta
    coefficients(2) = alpha*alpha + 2*(1._kp + alpha)*beta + 2*alpha*p
    coefficients(3) = beta*(3*alpha + 2*beta) + (alpha*alpha + 2*beta)*p
    coefficients(4) = 2*beta*(beta + alpha*p)
    coefficients(5) = beta*beta*p

    call quarticroots(coefficients,croots)

    j=0
    do i=1,4
       if (aimag(croots(i)).ne.0._kp) cycle
       j=j+1
       if (j.gt.2) stop 'rcipi_x_epstwozero: internal error!'
       rcipi_lnx_epstwozero(j) = real(croots(i))
    enddo
    
    call selectsort(rcipi_lnx_epstwozero)
    
  end function rcipi_lnx_epstwozero

  
!numerical accuracy limitation for large field values  
  function rcipi_numacc_x_potbig(p)
    implicit none
    real(kp) :: rcipi_numacc_x_potbig
    real(kp), intent(in) :: p
    real(kp), parameter :: big = epsilon(1._kp)*huge(1._kp)

    if (p.ne.0._kp) then    
       rcipi_numacc_x_potbig = big**(1._kp/p)
    else
       rcipi_numacc_x_potbig = big
    end if
    
  end function rcipi_numacc_x_potbig
  

!returns the solution of epsilon1(x)=1, in the domain in which V>0. At
!most five, zero values should be ignored.
  function rcipi_x_epsoneunity(p,alpha,beta)
    implicit none
    integer, parameter :: nepsmax = 5
    real(kp), dimension(nepsmax) :: rcipi_x_epsoneunity
    real(kp), intent(in) :: p,alpha,beta

    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: rcipiData

    real(kp), dimension(5) :: xeps
    real(kp), dimension(2) :: xepstwozero, epsones
    real(kp), dimension(2) :: xdpotzero, xpotzero

    real(kp) :: mini, maxi
    integer :: neps2

    rcipiData%real1 = p
    rcipiData%real2 = alpha
    rcipiData%real3 = beta

!the potential is piecewise positive    
    if (rcipi_check_potzero(alpha,beta)) then

       xpotzero = rcipi_x_potzero(alpha,beta)
       xdpotzero = rcipi_x_derivpotzero(p,alpha,beta)
       
       if (beta.ge.0._kp) then
!positive around x=0 and at infinity
          mini = epsilon(1._kp)
          maxi = xdpotzero(1)
          xeps(1) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

          mini = xdpotzero(1)
          maxi = xpotzero(1)
          xeps(2) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

          xeps(3) = 0._kp
          xeps(4) = 0._kp
          
          mini = xpotzero(2)
          maxi = rcipi_numacc_x_potbig(p)
          xeps(5) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

       else
!positive only between the potential zeros
          xeps(1) = 0._kp
          xeps(2) = 0._kp
          
          mini = xpotzero(1)
          maxi = xdpotzero(2)
          xeps(3) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

          mini = xdpotzero(2)
          maxi = xpotzero(2)
          xeps(4) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

          xeps(5) = 0._kp
          
       endif
          
    else
!the potential is never negative
       xepstwozero = rcipi_x_epstwozero(p,alpha,beta)
       neps2 = count(xepstwozero.ne.0._kp)
       
       if (neps2.eq.2) then
          epsones(1) = rcipi_epsilon_one(xepstwozero(1),p,alpha,beta)
          epsones(2) = rcipi_epsilon_one(xepstwozero(2),p,alpha,beta)
       else
          mini = epsilon(1._kp)
          maxi = rcipi_numacc_x_potbig(p)
          xeps(1:4) = 0._kp
          xeps(5) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)
          rcipi_x_epsoneunity = xeps
          return
       endif
       
       if (rcipi_check_derivpotzero(p,alpha,beta)) then
!at most five roots
          xdpotzero = rcipi_x_derivpotzero(p,alpha,beta)
          if (epsones(2).lt.1._kp) then
             xeps(5) = 0._kp
             xeps(4) = 0._kp
          else
             mini = xepstwozero(2)
             maxi = rcipi_numacc_x_potbig(p)
             xeps(5) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)
             mini = xdpotzero(2)
             maxi = xepstwozero(2)
             xeps(4) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)
          endif

          if (epsones(1).lt.1._kp) then
             xeps(3) = 0._kp
             xeps(2) = 0._kp
          else
             mini = xepstwozero(1)
             maxi = xdpotzero(2)
             xeps(3) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)
             mini = xdpotzero(1)
             maxi = xepstwozero(1)
             xeps(2) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)
          endif

          mini = epsilon(1._kp)
          maxi = xdpotzero(1)
          xeps(1) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

          !at most three roots       
       else

          xeps(3) = 0._kp
          xeps(2) = 0._kp
                    
          if (epsones(2).lt.1._kp) then
             xeps(5) = 0._kp
             xeps(4) = 0._kp
             mini = epsilon(1._kp)
             maxi = xepstwozero(1)
             xeps(1) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

          elseif (epsones(1).lt.1._kp) then
             mini = xepstwozero(2)
             maxi = rcipi_numacc_x_potbig(p)
             xeps(5) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

             mini = xepstwozero(1)
             maxi = xepstwozero(2)
             xeps(4) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

             mini = epsilon(1._kp)
             maxi = xepstwozero(1)
             xeps(1) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)

          elseif (epsones(1).ge.1._kp) then
             xeps(4) = 0._kp
             xeps(1) = 0._kp
             mini = xepstwozero(2)
             maxi = rcipi_numacc_x_potbig(p)
             xeps(5) = zbrent(find_rcipi_x_epsoneunity,mini,maxi,tolFind,rcipiData)
          endif

       endif

    endif

    rcipi_x_epsoneunity = xeps
    
    
  end function rcipi_x_epsoneunity

  function find_rcipi_x_epsoneunity(x,rcipiData)
    implicit none
    real(kp) :: find_rcipi_x_epsoneunity
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcipiData
    real(kp) :: p,alpha,beta

    p = rcipiData%real1
    alpha = rcipiData%real2
    beta = rcipiData%real3

    find_rcipi_x_epsoneunity = rcipi_epsilon_one(x,p,alpha,beta)-1._kp

  end function find_rcipi_x_epsoneunity


  
  function rcipi_efold_primitive(x,p,alpha,beta)
    use specialinf, only : cei, ei
    implicit none
    real(kp) :: rcipi_efold_primitive
    real(kp), intent(in) :: x,p,alpha,beta
    complex(kp) :: sqrtdelta, za, zb
    complex(kp) :: expinta, expintb,primitive
    real(kp) :: lnx, sqrt4bma2

    lnx = 0.5_kp*log(x*x)

!the bloody special cases
    if ((p.eq.0._kp).and.(alpha.eq.0._kp).and.(beta.eq.0._kp)) then
       stop 'rcipi_efold_primitive: de-Sitter model not supported!'
    endif
    
    if ((p.eq.0._kp).and.(beta.ne.0._kp)) then
       rcipi_efold_primitive = (-(((alpha**2 - 4*beta) &
            *ei(alpha/beta + 2*lnx))/exp(alpha/beta)) &
            + beta*x**2*(alpha - beta + 2*beta*lnx))/(8*beta**2)
       return
    endif

    if ((beta.eq.0._kp).and.(p.ne.0._kp).and.(alpha.ne.0._kp)) then
       rcipi_efold_primitive = x**2/(2._kp*p) - ei(2*(1._kp/alpha + 1._kp/p + lnx)) &
            /(exp(2*(1._kp/alpha + 1._kp/p))*p**2)
       return
    endif

    if ((beta.eq.0._kp).and.(p.eq.0._kp).and.alpha.ne.0._kp) then
       rcipi_efold_primitive = (x**2*(2 - alpha + 2*alpha*lnx))/(4._kp*alpha)
       return
    endif

    sqrtdelta = sqrt(cmplx(p*p*(alpha*alpha-4*beta) + 4*beta*beta))

    if (sqrtdelta.eq.0._kp) then
       sqrt4bma2 = sqrt(4._kp*beta - alpha*alpha)
       rcipi_efold_primitive = 0.25_kp*(-((exp((-(alpha*beta) &
            + sqrt(-((alpha**2 - 4*beta)* beta**2)))/beta**2)*sqrt4bma2 &
            *(alpha**2*beta - beta* sqrt(beta**2*(-alpha**2 + 4*beta)) + beta**2*(-4._kp + sqrt4bma2) &
            + sqrt(-((alpha**2 - 4*beta)*beta**2))*sqrt4bma2) &
            * ei((alpha*beta - sqrt(beta**2*(-alpha**2 + 4*beta)))/beta**2 + 2*lnx)&
            )/beta**4) - (sqrt4bma2*(alpha**2*beta + beta*sqrt(beta**2*(-alpha**2 + 4*beta)) &
            + beta**2*(-4._kp + sqrt4bma2) - sqrt(-((alpha**2 - 4*beta)*beta**2))*sqrt4bma2) &
            * ei((alpha*beta + sqrt(beta**2*(-alpha**2 + 4*beta)))/beta**2 + 2*lnx) &
            )/(beta**4*exp((alpha*beta + sqrt(beta**2*(-alpha**2 + 4*beta)))/beta**2)) &
            + (x**2*(alpha**4 + alpha**2*beta*(-8._kp + sqrt4bma2) - 2*beta**2*(-8._kp + sqrt4bma2) &
            + alpha**3*sqrt4bma2 - 4*alpha*beta*sqrt4bma2 + 2*beta*(alpha**2 - 4*beta + alpha*beta) &
            * sqrt4bma2*lnx + 2*beta**3*sqrt4bma2*lnx**2))/(beta**2*(alpha**2 - 2*beta &
            + 2*alpha*beta*lnx + 2*beta**2*lnx**2)))
       return
    endif

!the generic case    
    za = (-sqrtdelta + p*alpha + 2*beta)/p/beta
    zb = (sqrtdelta + p*alpha + 2*beta)/p/beta

    if (aimag(za).ne.0._kp) then
       expinta = cei(za+2._kp*lnx)
    else
       expinta = ei(real(za+2._kp*lnx,kp))
    endif
    if (aimag(zb).ne.0._kp) then
       expintb = cei(zb+2._kp*lnx)
    else
       expintb = ei(real(zb+2._kp*lnx,kp))
    endif
       
    primitive = (exp(zb)*p*sqrtdelta*x**2 - 2*exp((2*sqrtdelta)/(beta*p)) &
         *(-2*beta + sqrtdelta)*expinta - 2*(2*beta + sqrtdelta)*expintb) &
         /(2._kp*exp(zb)*p**2*sqrtdelta)
    
    rcipi_efold_primitive = real(primitive,kp)

    if (isnan(rcipi_efold_primitive)) then
       print *,'za,zb= ', za+2._kp*lnx, zb+2._kp*lnx
       print *,'sqrdelta= ',sqrtdelta
       print *,'cei',expinta, expintb
       stop 'rcipi_efold_primitive: Congrats, you got a NaN!'
    endif

  end function rcipi_efold_primitive


  function find_rcipi_x_trajectory(x,rcipiData)
    implicit none
    real(kp) :: find_rcipi_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcipiData

    real(kp) :: p, alpha, beta, NplusNuend

    p = rcipiData%real1
    alpha = rcipiData%real2
    beta = rcipiData%real3
    NplusNuend = rcipiData%real4

    find_rcipi_x_trajectory = rcipi_efold_primitive(x,p,alpha,beta) - NplusNuend
    
  end function find_rcipi_x_trajectory

  
end module rcipicommon
