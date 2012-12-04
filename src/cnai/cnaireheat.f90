!constant ns A inflation reheating functions in the slow-roll approximations

module cnaireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use cnaisr, only : cnai_epsilon_one, cnai_epsilon_two, cnai_epsilon_three
  use cnaisr, only : cnai_norm_potential
  use cnaisr, only : cnai_x_endinf, cnai_efold_primitive
  implicit none

  private

  public cnai_x_star, cnai_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function cnai_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: cnai_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: cnaiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = cnai_x_endinf(alpha)
    epsOneEnd = cnai_epsilon_one(xEnd,alpha)
    potEnd = cnai_norm_potential(xEnd,alpha)
    primEnd = cnai_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    cnaiData%real1 = alpha    
    cnaiData%real2 = w
    cnaiData%real3 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd*(1._kp-epsilon(1._kp))

    x = zbrent(find_cnai_x_star,mini,maxi,tolzbrent,cnaiData)
    cnai_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (cnai_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnai_x_star

  function find_cnai_x_star(x,cnaiData)   
    implicit none
    real(kp) :: find_cnai_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnaiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=cnaiData%real1
    w = cnaiData%real2
    CalFplusprimEnd = cnaiData%real3

    primStar = cnai_efold_primitive(x,alpha)
    epsOneStar = cnai_epsilon_one(x,alpha)
    potStar = cnai_norm_potential(x,alpha)

    find_cnai_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_cnai_x_star



  function cnai_lnrhoend(alpha,Pstar) 
    implicit none
    real(kp) :: cnai_lnrhoend
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = cnai_x_endinf(alpha)
    potEnd  = cnai_norm_potential(xEnd,alpha)
    epsOneEnd = cnai_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = cnai_x_star(alpha,wrad,junk,Pstar)    
    potStar = cnai_norm_potential(x,alpha)
    epsOneStar = cnai_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'cnai_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    cnai_lnrhoend = lnRhoEnd

  end function cnai_lnrhoend

  
end module cnaireheat
