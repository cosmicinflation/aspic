!Mixed large field reheating functions in the slow-roll

module mlfireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use mlfisr, only : mlfi_epsilon_one, mlfi_epsilon_two
  use mlfisr, only : mlfi_norm_potential
  use mlfisr, only : mlfi_x_endinf, mlfi_efold_primitive
  implicit none

  private

  public mlfi_x_star, mlfi_lnrhoend 

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function mlfi_x_star(p,q,alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: mlfi_x_star
    real(kp), intent(in) :: p,q,alpha,lnRhoReh,w,Pstar
    real(kp), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: mlfiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = mlfi_x_endinf(p,q,alpha)
    epsOneEnd = mlfi_epsilon_one(xEnd,p,q,alpha)
    potEnd = mlfi_norm_potential(xEnd,p,q,alpha)
    primEnd = mlfi_efold_primitive(xEnd,p,q,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    mlfiData%real1 = p
    mlfiData%real2 = q
    mlfiData%real3 = alpha
    mlfiData%real4 = w
    mlfiData%real5 = calF + primEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_mlfi_x_star,mini,maxi,tolzbrent,mlfiData)
    mlfi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (mlfi_efold_primitive(x,p,q,alpha) - primEnd)
    endif
    
  end function mlfi_x_star

  function find_mlfi_x_star(x,mlfiData)   
    implicit none
    real(kp) :: find_mlfi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mlfiData

    real(kp) :: primStar,alpha,p,q,w,CalFplusprimEnd,potStar,epsOneStar

    p=mlfiData%real1
    q = mlfiData%real2
    alpha = mlfiData%real3
    w = mlfiData%real4
    CalFplusprimEnd = mlfiData%real5

    primStar = mlfi_efold_primitive(x,p,q,alpha)
    epsOneStar = mlfi_epsilon_one(x,p,q,alpha)
    potStar = mlfi_norm_potential(x,p,q,alpha)

    find_mlfi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)

  end function find_mlfi_x_star



  function mlfi_lnrhoend(p,q,alpha,Pstar) 
    implicit none
    real(kp) :: mlfi_lnrhoend
    real(kp), intent(in) :: p,q,alpha,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk = 0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = mlfi_x_endinf(p,q,alpha)      
    potEnd  = mlfi_norm_potential(xEnd,p,q,alpha)
    epsEnd = mlfi_epsilon_one(xEnd,p,q,alpha)

!   Trick to return x such that rho_reh=rho_end
       
    x = mlfi_x_star(p,q,alpha,wrad,junk,Pstar)    
    potStar = mlfi_norm_potential(x,p,q,alpha)
    epsOneStar = mlfi_epsilon_one(x,p,q,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'mlfi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsEnd,potEnd/potStar)

    mlfi_lnrhoend = lnRhoEnd

  end function mlfi_lnrhoend

  

    
end module mlfireheat
