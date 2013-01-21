!slow-roll functions for the KKLT potential in the small field
!potential limit: valid for phi>>mu only. These functions are
!for testing purpose
!
!V(phi) = M^4/[1 + (mu/phi)^p] ~ M^4 [1 - (phi/mu)^(-p)]
!
!x = phi/mu

module kksfsrevol
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use sfsrevol, only : sf_norm_potential, sf_epsilon_one, sf_epsilon_two
  use sfsrevol, only : sf_nufunc, sf_x_endinf, sf_x_trajectory
  use sfsrevol, only : sf_x_epstwounity
  use sfsrevol, only : find_sfbi_x_epstwounity, find_sfbi_x_trajectory, find_sfbi_x_endinf

  implicit none

  private

  public  kksf_norm_potential, kksf_epsilon_one, kksf_epsilon_two
  public  kksf_x_endinf, kksf_nufunc, kksf_x_trajectory, kksf_x_epstwounity
 
contains
!returns V/M^4
  function kksf_norm_potential(x,p)
    implicit none
    real(kp) :: kksf_norm_potential
    real(kp), intent(in) :: x,p

    kksf_norm_potential = sf_norm_potential(x,-p)
  end function kksf_norm_potential


!epsilon1(x)
  function kksf_epsilon_one(x,p,mu)    
    implicit none
    real(kp) :: kksf_epsilon_one
    real(kp), intent(in) :: x,p,mu
    
    kksf_epsilon_one = sf_epsilon_one(x,-p,mu)
    
  end function kksf_epsilon_one


!epsilon2(x)
  function kksf_epsilon_two(x,p,mu)    
    implicit none
    real(kp) :: kksf_epsilon_two
    real(kp), intent(in) :: x,p,mu
    
    kksf_epsilon_two = sf_epsilon_two(x,-p,mu)
    
  end function kksf_epsilon_two



!this is nu(x)=integral[V(phi)/V'(phi) dphi]
  function kksf_nufunc(x,p,mu)
    implicit none
    real(kp), intent(in) :: x,p,mu
    real(kp) :: kksf_nufunc
        
    if (p.eq.0._kp) stop 'kksf_nufunc: p=0 is singular'

    kksf_nufunc = sf_nufunc(x,-p,mu)

  end function kksf_nufunc


  

!returns x at the end of inflation defined as epsilon1=1
  function kksf_x_endinf(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: kksf_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kksfData

    mini = epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    kksfData%real1 = -p
    kksfData%real2 = mu

    kksf_x_endinf = zbrent(find_sfbi_x_endinf,mini,maxi,tolFind,kksfData)
   
  end function kksf_x_endinf
  



!returns x at bfold=-efolds before the end of inflation
  function kksf_x_trajectory(bfold,xend,p,mu)
    implicit none
    real(kp), intent(in) :: bfold, p, mu, xend
    real(kp) :: kksf_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kksfData

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    kksfData%real1 = -p
    kksfData%real2 = mu
    kksfData%real3 = -bfold + kksf_nufunc(xend,p,mu)
    
    kksf_x_trajectory = zbrent(find_sfbi_x_trajectory,mini,maxi,tolFind,kksfData)
    
   
  end function kksf_x_trajectory

    
!returns x given epsilon2. If x~1, the small field approx may be not
!good enough, use kklt instead
  function kksf_x_epstwounity(eps2,p,mu)   
    implicit none
    real(kp), intent(in) :: p,mu,eps2
    real(kp) :: kksf_x_epstwounity
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kksfData

    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    kksfData%real1 = -p
    kksfData%real2 = mu
    kksfData%real3 = eps2

    kksf_x_epstwounity = zbrent(find_sfbi_x_epstwounity,mini,maxi,tolFind,kksfData)
   
 end function kksf_x_epstwounity

    
end module kksfsrevol
