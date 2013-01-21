!slow-roll functions for the effective one field valley hybrid inflation potential
!
!V(phi) = M^4 [1 + x^p]
!
!x = phi/mu
!mu = mu/Mp


module vhisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use inftools, only : zbrent
  implicit none

  private

  public vhi_norm_potential, vhi_epsilon_one, vhi_epsilon_two, vhi_epsilon_three
  public vhi_efold_primitive, vhi_x_trajectory
  public vhi_norm_deriv_potential, vhi_norm_deriv_second_potential
  public vhi_epsilon_one_max, vhi_x_epsonemax, vhi_x_epsoneunity
  public vhi_xinimax, vhi_xendmax,vhi_xendmin



  real(kp), parameter :: phihuge = 10000._kp


contains



!returns V/M**4
  function vhi_norm_potential(x,p,mu)
    implicit none
    real(kp) :: vhi_norm_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

    vhi_norm_potential = 1._kp+x**p


  end function vhi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function vhi_norm_deriv_potential(x,p,mu)
    implicit none
    real(kp) :: vhi_norm_deriv_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

   vhi_norm_deriv_potential = p*x**(p-1._kp)

  end function vhi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function vhi_norm_deriv_second_potential(x,p,mu)
    implicit none
    real(kp) :: vhi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

    vhi_norm_deriv_second_potential = p*(p-1._kp)*x**(p-2._kp)

  end function vhi_norm_deriv_second_potential



!epsilon_one(x)
  function vhi_epsilon_one(x,p,mu)    
    implicit none
    real(kp) :: vhi_epsilon_one
    real(kp), intent(in) :: x,p,mu

  
    vhi_epsilon_one = p**2/(2._kp*mu**2)*x**(2._kp*p-2._kp)/&
                      ((1._kp+x**p)**2)

    
  end function vhi_epsilon_one


!epsilon_two(x)
  function vhi_epsilon_two(x,p,mu)    
    implicit none
    real(kp) :: vhi_epsilon_two
    real(kp), intent(in) :: x,p,mu
    
    vhi_epsilon_two = 2._kp*p/(mu**2)*x**(p-2._kp)*(x**p-p+1._kp)/&
                      ((1._kp+x**p)**2)
    
  end function vhi_epsilon_two


!epsilon_three(x)
  function vhi_epsilon_three(x,p,mu)    
    implicit none
    real(kp) :: vhi_epsilon_three
    real(kp), intent(in) :: x,p,mu
    
    vhi_epsilon_three = p/(mu**2)*x**(p-2._kp)/((1._kp+x**p)**2*(x**p-p+1._kp))*&
                        (2._kp*x**(2._kp*p)-(p-1._kp)*(p+4._kp)*x**p &
                        +(p-1._kp)*(p-2._kp))
    
  end function vhi_epsilon_three



!the maximal value of eps1 
  function vhi_epsilon_one_max(p,mu)
    implicit none
    real(kp) :: vhi_epsilon_one_max
    real(kp), intent(in) :: p,mu


    if (p.lt.1._kp) stop 'vhi_epsilon_one_max: no local maximum for eps1 with p<1!'

    if (p.eq.1._kp) then
       vhi_epsilon_one_max = 0.5_kp/mu/mu
       return
    endif

    vhi_epsilon_one_max= 0.5_kp*((p-1._kp)*(p-1._kp))**(1._kp-1._kp/p)/mu/mu


  end function vhi_epsilon_one_max


!the field value at which epsilon is maximal
 function vhi_x_epsonemax(p,mu)
    implicit none
    real(kp) :: vhi_x_epsonemax
    real(kp), intent(in) :: p
    real(kp), intent(in), optional :: mu    

    if (p.lt.1._kp) stop 'vhi_x_epsonemax: no local maximum for eps1 with p<1!'

    if (p.eq.1._kp) then
       vhi_x_epsonemax = 0._kp
    else
       vhi_x_epsonemax = (p-1._kp)**(1._kp/p)
    endif

  end function vhi_x_epsonemax



!returns a vector containing the two roots of eps1=1
  function vhi_x_epsoneunity(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp), dimension(2) :: vhi_x_epsoneunity
    real(kp) :: xeps1max
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: vhiData


    vhiData%real1 = p
    vhiData%real2 = mu

    if (p.lt.1._kp) then
!only one root as eps1 blows up for x->0
       
       mini=epsilon(1._kp)
       maxi = phihuge*p/mu
       
       vhi_x_epsoneunity(2) &
            = zbrent(find_vhi_x_epsoneunity,mini,maxi,tolFind,vhiData)
       vhi_x_epsoneunity(1) = vhi_x_epsoneunity(2)

    elseif (p.eq.1_kp) then
!eps1 is maximal in x=0, only one root
       if (vhi_epsilon_one_max(p,mu).lt.1._kp) then
          stop 'vhi_x_epsoneunity: no solution for eps1=1'
       endif

       vhi_x_epsoneunity(2) = 1._kp/(mu*sqrt(2._kp)) - 1._kp
       vhi_x_epsoneunity(1) = vhi_x_epsoneunity(2)
       

    elseif (p.gt.1_kp) then

       if (vhi_epsilon_one_max(p,mu).lt.1._kp) then
          stop 'vhi_x_epsoneunity: no solution for eps1=1'
       endif
       xeps1max = vhi_x_epsonemax(p,mu)    
!smallest root
       mini = epsilon(1._kp)
       maxi = xeps1max + epsilon(1._kp)

       vhi_x_epsoneunity(1) &
            = zbrent(find_vhi_x_epsoneunity,mini,maxi,tolFind,vhiData)

!largest root
       mini = xeps1max - epsilon(1._kp)
       maxi = phihuge* p/mu
           
       vhi_x_epsoneunity(2) &
            = zbrent(find_vhi_x_epsoneunity,mini,maxi,tolFind,vhiData)

    else

       stop 'vhi_x_epsoneunity: internal error!'
        
    endif


  end function vhi_x_epsoneunity
  
  function find_vhi_x_epsoneunity(x,vhiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: vhiData
    real(kp) :: find_vhi_x_epsoneunity
    real(kp) :: p,mu

    p = vhiData%real1
    mu = vhiData%real2
    
    find_vhi_x_epsoneunity = vhi_epsilon_one(x,p,mu) - 1._kp
  
  end function find_vhi_x_epsoneunity


!returns the minimum value for xend as a function of p and mu such
!that inflation stops by instability.
  function vhi_xendmin(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: vhi_xendmin
    real(kp), parameter :: tolFind=tolkp
    real(kp), dimension(2) :: xepsones   
    
    vhi_xendmin=epsilon(1._kp)

    if (p.eq.1._kp) then

       if (vhi_epsilon_one_max(p,mu).gt.1._kp) then
          xepsones = vhi_x_epsoneunity(p,mu)
          vhi_xendmin = xepsones(2)
       endif
    elseif (p.lt.1._kp) then
       
       xepsones = vhi_x_epsoneunity(p,mu)
       vhi_xendmin = xepsones(2)

    endif

  end function vhi_xendmin



!returns the maximum value for xini as a function of p and mu
  function vhi_xinimax(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: vhi_xinimax
    real(kp), parameter :: tolFind=tolkp
    real(kp), dimension(2) :: xepsones 
   
    vhi_xinimax = phihuge*p/mu

    if ((p.gt.1._kp).and.(vhi_epsilon_one_max(p,mu).gt.1._kp)) then
       xepsones = vhi_x_epsoneunity(p,mu)
       vhi_xinimax = xepsones(1)
    endif

  end function vhi_xinimax


!returns the maximal value for xend such that there are efold number
!of inflation from xinimax
  function vhi_xendmax(efold,p,mu)
    implicit none
    real(kp), intent(in) :: efold,p,mu
    real(kp) :: vhi_xendmax, xiniMax, xendMin
    real(kp), parameter :: tolFind=tolkp

    real(kp) :: efoldMax
    
    xiniMax = vhi_xinimax(p,mu)
    xendMin = vhi_xendmin(p,mu)

    efoldMax = -vhi_efold_primitive(xendMin,p,mu) &
         + vhi_efold_primitive(xiniMax,p,mu)

    if (efold.gt.efoldMax) then
       write(*,*)'vhi_xendmax: not enough efolds!'
       write(*,*)'efold requested=   efold maxi= ',efold,efoldMax
       stop
    endif

    vhi_xendmax = vhi_x_trajectory(efold,xiniMax,p,mu)
       

  end function vhi_xendmax

  

!this is integral[V(phi)/V'(phi) dphi]
  function vhi_efold_primitive(x,p,mu)
    implicit none
    real(kp), intent(in) :: x,p,mu
    real(kp) :: vhi_efold_primitive

    if (p .eq. 2._kp) then

      vhi_efold_primitive = mu**2/4._kp*(x**2+2._kp*log(x))

    else

      vhi_efold_primitive = mu**2/(2._kp*p)*(x**2-2._kp/(p-2._kp)*x**(2._kp-p))

    endif

  end function vhi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function vhi_x_trajectory(bfold,xend,p,mu)
    implicit none
    real(kp), intent(in) :: bfold, p,mu, xend
    real(kp) :: vhi_x_trajectory
    real(kp), parameter :: tolFind = tolkp
    real(kp) :: mini,maxi
    logical :: approx = .false. !set to true if one wants to use the approximate phi/mu << 1 formula
    type(transfert) :: vhiData

!trick land here: if bfold>0 is inputed, allows to return xend by
!input of xini (ah ah ah)

    if (p .eq. 1._kp) then
       
       vhi_x_trajectory=-1._kp+sqrt(1._kp-2._kp/(mu**2)*bfold+xend**2+2._kp*xend)
      
    else if (p .eq. 2._kp) then 

       vhi_x_trajectory=sqrt(lambert(xend**2*exp(xend**2-4._kp*bfold/(mu**2)),0))
       
    else if (approx) then

       vhi_x_trajectory=(xend**(2._kp-p)+p*(p-2._kp)*bfold/(mu**2))**(1._kp/(2._kp-p))

    else

 
       if (bfold.le.0._kp) then
          mini = xend*(1._kp+epsilon(1._kp))
          maxi= vhi_xinimax(p,mu)*(1._kp-epsilon(1._kp))
       else
          maxi = xend*(1._kp+epsilon(1._kp))
          mini = vhi_xendmin(p,mu)*(1._kp-epsilon(1._kp))
       endif

       vhiData%real1 = p
       vhiData%real2 = mu
       vhiData%real3 = -bfold + vhi_efold_primitive(xend,p,mu)

       vhi_x_trajectory = zbrent(find_vhi_x_trajectory,mini,maxi,tolFind,vhiData)

    end if
       
  end function vhi_x_trajectory


  function find_vhi_x_trajectory(x,vhiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: vhiData
    real(kp) :: find_vhi_x_trajectory
    real(kp) :: p,mu,NplusNuend

    p= vhiData%real1
    mu= vhiData%real2
    NplusNuend = vhiData%real3

    find_vhi_x_trajectory = vhi_efold_primitive(x,p,mu) - NplusNuend
   
  end function find_vhi_x_trajectory



end module vhisr
