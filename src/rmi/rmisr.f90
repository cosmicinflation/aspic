!slow-roll functions for the running mass potential
!
!V(phi) = M^4 { 1 + nu [1/p - ln(x)](mu*x)^p }
!
!x = phi/mu
!
!These models are sorted according to the sign of nu and x-1.
!RM1 is nu>0 x<1 (x decreases during inflation)
!RM2 is nu>0 x>1 (x increases during inflation)
!RM3 is nu<0 x<1 (x increases during inflation)
!RM4 is nu<0 x>1 (x decreases during inflation)

module rmsrevol
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : hypergeom_2F1, dei
  implicit none

  private

  public rm_norm_potential, rm_epsilon_one, rm_epsilon_two
  public rm_x_endinf, rm_nufunc, rm_x_trajectory
 
contains
!returns V/M^4
  function rm_norm_potential(x,p,mu,nu)
    implicit none
    real(kp) :: rm_norm_potential
    real(kp), intent(in) :: x,p,mu,nu

    if (p.lt.2._kp) stop 'rm_norm_potential: p<2 not valid!'

    rm_norm_potential = 1._kp + nu*(1._kp/p &
         - log(x))*(mu*x)**p 

 end function rm_norm_potential


!epsilon1(x)
  function rm_epsilon_one(x,p,mu,nu)
    implicit none
    real(kp) :: rm_epsilon_one
    real(kp), intent(in) :: x,p,mu,nu
    real(kp) :: y, lnx
    
    y = mu*x
    lnx = log(x)
    
    rm_epsilon_one = (nu**2*p**2 * y**(-2._kp + 2*p) * lnx**2) &
         /(2._kp*(1._kp + nu*y**p * (1._kp/p - lnx))**2)
    
  end function rm_epsilon_one


!epsilon2(x)
  function rm_epsilon_two(x,p,mu,nu)
    implicit none
    real(kp) :: rm_epsilon_two
    real(kp), intent(in) :: x,p,mu,nu
    real(kp) :: y,lnx

    y = x*mu
    lnx = log(x)

    rm_epsilon_two = (2._kp*nu*p**2*y**(-2._kp + p) &
         * (p + nu*y**p + lnx &
         * ((-1._kp + p)*p + nu*y**p* (-1 + p*lnx)))) &
         /(p - nu*y**p * (-1 + p*lnx))**2
    
  end function rm_epsilon_two



!this is nu(x)=integral[V(phi)/V'(phi) dphi]
  function rm_nufunc(x,p,mu,nu)
    implicit none
    real(kp), intent(in) :: x,p,mu,nu
    real(kp) :: rm_nufunc

    real(kp) :: y,lnx,l
    real(kp) :: argei1, argei2

    l=nu*mu**p
    y = x*mu
    lnx = log(x)

    if (p == 2._kp) then
       argei1 = 2._kp*lnx
       rm_nufunc = x**2 - 2._kp*log(abs(lnx))/l &
            - dei(argei1)
    else
       argei1 = (2._kp-p)*lnx
       argei2 = 2._kp*lnx
       rm_nufunc = x**2 - (2._kp/l)*dei(argei1) - (2._kp/p)*dei(argei2)
    endif

    rm_nufunc = (0.5_kp*mu*mu/p) * rm_nufunc
   
  end function rm_nufunc


!return xstop or x such as eps1(x)=1
  function rm_x_endinf(p,mu,nu,xstop)
    implicit none
    real(kp) :: rm_x_endinf
    real(kp), intent(in) :: p,mu,nu,xstop
    real(kp) :: xend

!check if we are in RM2
    if ((nu.gt.0._kp).and.(xstop.gt.1._kp)) then
       xEnd = rm_x_epsoneisone(p,mu,nu)              
    else
       xEnd = xstop
    endif
!even in RM2, one may want to stop inflation by instability before
!reaching eps1=1
    if (xEnd.lt.xstop) then
       write(*,*)'rm_x_endinf: slow-roll violations stop inflation before'
       write(*,*)'xend= xstop= ',xend, xstop
    endif
    
    rm_x_endinf = xend

  end function rm_x_endinf


!returns x when eps1=1. However this equation has one solution only in
![mu,xZero] where V(xZero) = 0 and only for RM2. So in general, this
!is not the end of inflation
  function rm_x_epsoneisone(p,mu,nu)
    implicit none
    real(kp), intent(in) :: p,mu,nu
    real(kp) :: rm_x_epsoneisone
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,xzero
    type(transfert) :: rmData

    if (nu.lt.0._kp) stop 'rm_x_endinf: nu<0, xend=xstop'

!find xzero first (for x>1 RM2)   
    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    rmData%real1 = p
    rmData%real2 = mu
    rmData%real3 = nu              

    xzero = zbrent(find_rmzeropot,mini,maxi,tolFind,rmData)

!we know were to look for xend
    mini = 1._kp + epsilon(1._kp)
    maxi = xzero

    rmData%real1 = p
    rmData%real2 = mu
    rmData%real3 = nu

    rm_x_epsoneisone = zbrent(find_rmepsoneisone,mini,maxi,tolFind,rmData)

  end function rm_x_epsoneisone

  function find_rmzeropot(x,rmData)
!vanished for V=0
    implicit none
    real(kp) :: find_rmzeropot
    real(kp), intent(in) :: x 
    type(transfert), optional, intent(inout) :: rmData
    real(kp) :: p,mu,nu
    
    p=rmData%real1
    mu=rmData%real2
    nu=rmData%real3
    
    find_rmzeropot = rm_norm_potential(x,p,mu,nu)

  end function find_rmzeropot
  
  function find_rmepsoneisone(x,rmData)
!vanishes for eps1=1
    implicit none
    real(kp), intent(in) :: x 
    type(transfert), optional, intent(inout) :: rmData
    real(kp) :: find_rmepsoneisone
    real(kp) :: p,mu,nu,lnx,y

    p=rmData%real1
    mu=rmData%real2
    nu=rmData%real3
       
    y = mu*x
    lnx = log(x)

    find_rmepsoneisone = 1._kp + nu*(1._kp/p - lnx)*y**p &
         - (nu*p/sqrt(2._kp))*lnx*y**(p-1._kp)

  end function find_rmepsoneisone
 
       
  
!returns x at bfold=-efolds before the end of inflation
  function rm_x_trajectory(bfold,xend,p,mu,nu)
    implicit none
    real(kp), intent(in) :: bfold, p, mu,nu, xend
    real(kp) :: rm_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rmData
  
    if (xend.eq.1._kp) stop'rm_x_trajectory: xend = 1'


    if (xend.lt.1._kp) then
       mini = epsilon(1._kp)
       maxi = 1._kp - epsilon(1._kp)
    elseif (xend.gt.1._kp) then
       mini = 1._kp + epsilon(1._kp)
       if (nu.gt.0._kp) then
          maxi = 1._kp
       else
          maxi = 1._kp/epsilon(1._kp)
       endif
    else
       stop 'rm_x_trajectory: internal error'
    endif    

    rmData%real1 = p
    rmData%real2 = mu
    rmData%real3 = nu
    rmData%real4 = -bfold + rm_nufunc(xend,p,mu,nu)
    
    rm_x_trajectory = zbrent(find_rmtraj,mini,maxi,tolFind,rmData)
       
  end function rm_x_trajectory

  function find_rmtraj(x,rmData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: rmData
    real(kp) :: find_rmtraj
    real(kp) :: p,mu,nu,NplusNuend

    p=rmData%real1
    mu = rmData%real2
    nu = rmData%real3
    NplusNuend = rmData%real4

    find_rmtraj = rm_nufunc(x,p,mu,nu) - NplusNuend
   
  end function find_rmtraj

  
end module rmsrevol
