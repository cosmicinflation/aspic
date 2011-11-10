!Natural Inflation with the Plus Sign reheating functions in the slow-roll

module mnireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use mnisr, only : mni_epsilon_one, mni_epsilon_two
  use mnisr, only : mni_norm_potential
  use mnisr, only : mni_x_endinf, mni_efold_primitive
  implicit none

  private

  public mni_x_star, mni_lnrhoend 

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function mni_x_star(f,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: mni_x_star
    real(kp), intent(in) :: f,lnRhoReh,w,Pstar
    real(kp), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: mniData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = mni_x_endinf(f)
    epsOneEnd = mni_epsilon_one(xEnd,f)
    potEnd = mni_norm_potential(xEnd,f)
    primEnd = mni_efold_primitive(xEnd,f)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    mniData%real1 = f
    mniData%real2 = w
    mniData%real3 = calF + primEnd

    mini = xEnd
    maxi = f*acos(-1._kp)*(1._kp-epsilon(1._kp))

    x = zbrent(find_mni_x_star,mini,maxi,tolzbrent,mniData)
    mni_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (mni_efold_primitive(x,f) - primEnd)
    endif
    
  end function mni_x_star

  function find_mni_x_star(x,mniData)   
    implicit none
    real(kp) :: find_mni_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mniData

    real(kp) :: primStar,f,w,CalFplusprimEnd,potStar,epsOneStar

    f = mniData%real1
    w = mniData%real2
    CalFplusprimEnd = mniData%real3

    primStar = mni_efold_primitive(x,f)
    epsOneStar = mni_epsilon_one(x,f)
    potStar = mni_norm_potential(x,f)

    find_mni_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)

  end function find_mni_x_star



  function mni_lnrhoend(f,Pstar) 
    implicit none
    real(kp) :: mni_lnrhoend
    real(kp), intent(in) :: f,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk = 0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = mni_x_endinf(f)      
    potEnd  = mni_norm_potential(xEnd,f)
    epsEnd = mni_epsilon_one(xEnd,f)

!   Trick to return x such that rho_reh=rho_end
       
    x = mni_x_star(f,wrad,junk,Pstar)    
    potStar = mni_norm_potential(x,f)
    epsOneStar = mni_epsilon_one(x,f)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'mni_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsEnd,potEnd/potStar)

    mni_lnrhoend = lnRhoEnd

  end function mni_lnrhoend

  

    
end module mnireheat
