!Natural Inflation with the Plus Sign reheating functions in the slow-roll

module pnireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use pnisr, only : pni_epsilon_one, pni_epsilon_two
  use pnisr, only : pni_norm_potential
  use pnisr, only : pni_x_endinf, pni_efold_primitive
  implicit none

  private

  public pni_x_star, pni_lnrhoend 

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function pni_x_star(f,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: pni_x_star
    real(kp), intent(in) :: f,lnRhoReh,w,Pstar
    real(kp), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: pniData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = pni_x_endinf(f)
    epsOneEnd = pni_epsilon_one(xEnd,f)
    potEnd = pni_norm_potential(xEnd,f)
    primEnd = pni_efold_primitive(xEnd,f)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    pniData%real1 = f
    pniData%real2 = w
    pniData%real3 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd

    x = zbrent(find_pni_x_star,mini,maxi,tolzbrent,pniData)
    pni_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (pni_efold_primitive(x,f) - primEnd)
    endif
    
  end function pni_x_star

  function find_pni_x_star(x,pniData)   
    implicit none
    real(kp) :: find_pni_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: pniData

    real(kp) :: primStar,f,w,CalFplusprimEnd,potStar,epsOneStar

    f = pniData%real1
    w = pniData%real2
    CalFplusprimEnd = pniData%real3

    primStar = pni_efold_primitive(x,f)
    epsOneStar = pni_epsilon_one(x,f)
    potStar = pni_norm_potential(x,f)

    find_pni_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)

  end function find_pni_x_star



  function pni_lnrhoend(f,Pstar) 
    implicit none
    real(kp) :: pni_lnrhoend
    real(kp), intent(in) :: f,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk = 0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = pni_x_endinf(f)      
    potEnd  = pni_norm_potential(xEnd,f)
    epsEnd = pni_epsilon_one(xEnd,f)

!   Trick to return x such that rho_reh=rho_end
       
    x = pni_x_star(f,wrad,junk,Pstar)    
    potStar = pni_norm_potential(x,f)
    epsOneStar = pni_epsilon_one(x,f)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'pni_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsEnd,potEnd/potStar)

    pni_lnrhoend = lnRhoEnd

  end function pni_lnrhoend

  

    
end module pnireheat
