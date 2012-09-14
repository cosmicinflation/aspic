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

  public  vhi_norm_potential, vhi_epsilon_one, vhi_epsilon_two, vhi_epsilon_three
  public  vhi_efold_primitive, vhi_x_trajectory
  public  vhi_norm_deriv_potential, vhi_norm_deriv_second_potential
  public  vhi_xend_max,vhi_xend_min


contains



!returns V/M**4
  function vhi_norm_potential(x,p)
    implicit none
    real(kp) :: vhi_norm_potential
    real(kp), intent(in) :: x,p

    vhi_norm_potential = 1._kp+x**p


  end function vhi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function vhi_norm_deriv_potential(x,p)
    implicit none
    real(kp) :: vhi_norm_deriv_potential
    real(kp), intent(in) :: x,p

   vhi_norm_deriv_potential = p*x**(p-1._kp)

  end function vhi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function vhi_norm_deriv_second_potential(x,p)
    implicit none
    real(kp) :: vhi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p

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



!returns the minimum value for xend as a function of p and mu, so that xend should vary between vhi_xend_min and vhi_xend_max if inflation stops by critical instability.
  function vhi_xend_min(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: vhi_xend_min
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: vhiData
    
    if (p .gt. 1._kp) then

      vhi_xend_min=epsilon(1._kp)
  
    else if (p .eq. 1._kp) then

      if (sqrt(2._kp)*mu .lt. 1._kp) then !eps1>1 at x=0. Return x_eps.

        vhi_xend_min = 1._kp/(sqrt(2._kp)*mu)-1._kp

      else !eps<1 everywhere

        vhi_xend_min = epsilon(1._kp)

      endif

    else if (p .lt. 1._kp) then

        maxi= 10000._kp*p/mu
        mini=epsilon(1._kp)
        vhiData%real1 = p
        vhiData%real2 = mu
        vhi_xend_min = zbrent(find_vhi_xend_bound,mini,maxi,tolFind,vhiData)

    endif

  end function vhi_xend_min



!returns the maximum value for xend as a function of p and mu, so that xend should vary between vhi_xend_min and vhi_xend_max if inflation stops by critical instability.
  function vhi_xend_max(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: vhi_xend_max
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: vhiData
    
    if (p .gt. 1._kp) then

      if (mu*sqrt(2._kp) .lt. (p-1._kp)**((p-1._kp)/p)) then !eps1>1 for some values of x. Return x_eps^-.

        if (p .eq. 2._kp) then
          vhi_xend_max=1._kp/(sqrt(2._kp)*mu)*(1._kp-sqrt(1._kp-2._kp*mu**2))
        else
      
        maxi= (p-1._kp)**(1._kp/p) !value of x such that epsilon1 is maximum
        mini=0._kp
        vhiData%real1 = p
        vhiData%real2 = mu
        vhi_xend_max = zbrent(find_vhi_xend_bound,mini,maxi,tolFind,vhiData)

        endif

      else !eps1<1 everywhere
       
        vhi_xend_max = 10000._kp*p/mu
      
      endif


    else if (p .le. 1._kp) then

        vhi_xend_max = 10000._kp*p/mu

    endif

  end function vhi_xend_max

  function find_vhi_xend_bound(x,vhiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: vhiData
    real(kp) :: find_vhi_xend_bound
    real(kp) :: p,mu

    p = vhiData%real1
    mu = vhiData%real2
    
    find_vhi_xend_bound = vhi_epsilon_one(x,p,mu) - 1._kp
  
  end function find_vhi_xend_bound



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

    if (p .eq. 1._kp) then

      vhi_x_trajectory=-1._kp+sqrt(1._kp-2._kp/(mu**2)*bfold+xend**2+2._kp*xend)

    else if (p .eq. 2._kp) then 

      vhi_x_trajectory=sqrt(lambert(xend**2*exp(xend**2-4._kp*bfold/(mu**2)),0))

    else if (approx .eqv. .true.) then

      vhi_x_trajectory=(xend**(2._kp-p)+p*(p-2._kp)*bfold/(mu**2))**(1._kp/(2._kp-p))

    else

    mini = xend*(1._kp+epsilon(1._kp))
    maxi= vhi_xend_max(p,mu)*(1._kp-epsilon(1._kp))

    vhiData%real1 = p
    vhiData%real2 = mu
    vhiData%real3 = -bfold + vhi_efold_primitive(xend,p,mu)
    
    vhi_x_trajectory = zbrent(find_vhitraj,mini,maxi,tolFind,vhiData)

    end if
       
  end function vhi_x_trajectory


  function find_vhitraj(x,vhiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: vhiData
    real(kp) :: find_vhitraj
    real(kp) :: p,mu,NplusNuend

    p= vhiData%real1
    mu= vhiData%real2
    NplusNuend = vhiData%real3

    find_vhitraj = vhi_efold_primitive(x,p,mu) - NplusNuend
   
  end function find_vhitraj



end module vhisr
