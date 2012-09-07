!slow-roll functions for the one field effective hybrid potential
!
!V(phi) = M^4 [1 + x^p]
!
!x = phi/mu


module hysrevol
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public hy_norm_potential, hy_epsilon_one, hy_epsilon_two
  public hy_x_endinf, hy_nufunc, hy_x_trajectory, hy_check_slowroll
  public hy_x_epsoneismax, hy_mu_connex, hy_x_epsoneisone
 

contains
!returns V/M^4
  function hy_norm_potential(x,p)
    implicit none
    real(kp) :: hy_norm_potential
    real(kp), intent(in) :: x,p
    
    hy_norm_potential = 1._kp + x**p

 end function hy_norm_potential


!epsilon1(x)
  function hy_epsilon_one(x,p,mu)
    implicit none
    real(kp) :: hy_epsilon_one
    real(kp), intent(in) :: x,p,mu
    
    hy_epsilon_one = 0.5_kp*((p/mu)*x**(p-1._kp)/(1._kp+x**p))**2
    
  end function hy_epsilon_one


!epsilon2(x)
  function hy_epsilon_two(x,p,mu)
    implicit none
    real(kp) :: hy_epsilon_two
    real(kp), intent(in) :: x,p,mu
       
    hy_epsilon_two = 2._kp*p*x**(p-2._kp)*(1._kp-p+x**p) &
         /(mu*(1._kp+x**p))**2
    
  end function hy_epsilon_two



!this is nu(x)=integral[V(phi)/V'(phi) dphi]
  function hy_nufunc(x,p,mu)
    implicit none
    real(kp), intent(in) :: x,p,mu
    real(kp) :: hy_nufunc

    if (p == 2._kp) then
       hy_nufunc = x**2 + 2._kp * log(x)
    else
       hy_nufunc = x**2 + 2._kp/(2._kp-p) * x**(2._kp-p)
    endif

    if (p.eq.0._kp) stop 'hy_nufunc: p=0 is singular'

    hy_nufunc = 0.5_kp*hy_nufunc*mu*mu/p
      
  end function hy_nufunc

  

!return xstop or x such as eps1(x)=1
  function hy_x_endinf(p,mu,xstop)
    implicit none
    real(kp) :: hy_x_endinf
    real(kp), intent(in) :: p,mu,xstop
    real(kp) :: xtrans, xend
    real(kp), dimension(2) :: xeps

!no solution for eps1=1
    if (mu.gt.hy_mu_connex(p)) then
       hy_x_endinf = xstop
       return
    endif

!otherwise, we could be in the large field or small field regime
    xeps = hy_x_epsoneisone(p,mu)
   
    if ((xstop.gt.xeps(1)).and.(xstop.lt.xeps(2))) then
       write(*,*)'hy_x_endinf: xstop violates slow-roll!'
       write(*,*)'xstop= xeps= ',xstop,xeps
       write(*,*)'xend set at ',xeps(2)
       xend = xeps(2)
    else
       xend = xstop
    endif
       
    hy_x_endinf = xend

  end function hy_x_endinf



!returns x when eps1 is maximal, i.e. eps2=0
  function hy_x_epsoneismax(p)
    implicit none
    real(kp) :: hy_x_epsoneismax
    real(kp), intent(in) :: p

    hy_x_epsoneismax = (p - 1._kp)**(1._kp/p)

  end function hy_x_epsoneismax


!returns x when eps1=1, if it exists. This is not the end of inflation
!in general and there are two roots: large field and small field 
  function hy_x_epsoneisone(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp), dimension(2) :: hy_x_epsoneisone
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi, xtrans
    type(transfert) :: hyData

    if (mu.gt.hy_mu_connex(p)) stop 'hy_x_epsoneisone: mu > mu_connex'
    
    xtrans = hy_x_epsoneismax(p)

    hyData%real1 = p
    hyData%real2 = mu

    mini = 0._kp
    maxi = xtrans
    hy_x_epsoneisone(1) = zbrent(find_hyepsoneisone,mini,maxi,tolFind,hyData)

    mini = xtrans
    maxi = 1._kp/epsilon(1._kp)
    hy_x_epsoneisone(2) = zbrent(find_hyepsoneisone,mini,maxi,tolFind,hyData)

  end function hy_x_epsoneisone

   

  function hy_mu_connex(p)
    implicit none
    real(kp) :: hy_mu_connex
    real(kp), intent(in) :: p

    hy_mu_connex = p/sqrt(8._kp)

  end function hy_mu_connex



  function find_hyepsoneisone(x,hyData)
!vanishes for eps1=1
    implicit none
    real(kp), intent(in) :: x 
    type(transfert), optional, intent(inout) :: hyData
    real(kp) :: find_hyepsoneisone
    real(kp) :: p,mu

    p=hyData%real1
    mu=hyData%real2    
    
    find_hyepsoneisone = x**(p-1._kp) - sqrt(2._kp)*mu/p * (x**p + 1._kp)
        
  end function find_hyepsoneisone


!informs about some potential problem 
  subroutine hy_check_slowroll(x,xend,p,mu)
    implicit none
    real(kp), intent(in) :: x,xend,p,mu
    real(kp), dimension(2) :: xeps

    if (mu.gt.hy_mu_connex(p)) return

    xeps = hy_x_epsoneisone(p,mu)

    if ((xend.gt.xeps(1)).and.(xend.lt.xeps(2))) then
       stop 'hy_check_slowroll: eps(xend) > 1!'
    endif
    
    if ((xend.le.xeps(1)).and.(x.gt.xeps(1))) then
       write(*,*)'hy_check_slowroll: slow-roll will be violated:'
       write(*,*)'x= xend= xeps= ',x,xend,xeps
    endif

  end subroutine hy_check_slowroll
  

!returns x at bfold=-efolds before the end of inflation
  function hy_x_trajectory(bfold,xend,p,mu)
    implicit none
    real(kp), intent(in) :: bfold, p, mu, xend
    real(kp) :: hy_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi, x
    real(kp), dimension(2) :: xeps
    type(transfert) :: hyData
    
    hyData%real1 = p
    hyData%real2 = mu   
    hyData%real4 = -bfold + hy_nufunc(xend,p,mu)
    
    x = zbrent(find_hytraj,mini,maxi,tolFind,hyData)
    
    call hy_check_slowroll(x,xend,p,mu)
    
    hy_x_trajectory = x

  end function hy_x_trajectory


  function find_hytraj(x,hyData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: hyData
    real(kp) :: find_hytraj
    real(kp) :: p,mu,NplusNuend

    p=hyData%real1
    mu = hyData%real2
    NplusNuend = hyData%real4

    find_hytraj = hy_nufunc(x,p,mu) - NplusNuend
   
  end function find_hytraj

  
end module hysrevol
