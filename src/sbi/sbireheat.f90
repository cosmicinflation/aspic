!supergravity brane inflation reheating functions in the slow-roll approximations

module sbireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use sbisr, only : sbi_epsilon_one, sbi_epsilon_two, sbi_epsilon_three
  use sbisr, only : sbi_norm_potential
  use sbisr, only : sbi_x_endinf, sbi_efold_primitive
  implicit none

  private

  public sbi_x_star, sbi_lnrhoreh_max
  public sbi_x_rrad, sbi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function sbi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: sbi_x_star
    real(kp), intent(in) :: alpha,beta,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sbiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = sbi_epsilon_one(xEnd,alpha,beta)
    potEnd = sbi_norm_potential(xEnd,alpha,beta)
    primEnd = sbi_efold_primitive(xEnd,alpha,beta)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    sbiData%real1 = alpha 
    sbiData%real2 = beta
    sbiData%real3 = w
    sbiData%real4 = calF + primEnd

    maxi = xEnd*(1._kp-epsilon(1._kp))
    mini = epsilon(1._kp)

    x = zbrent(find_sbi_x_star,mini,maxi,tolzbrent,sbiData)
    sbi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (sbi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function sbi_x_star

  function find_sbi_x_star(x,sbiData)   
    implicit none
    real(kp) :: find_sbi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sbiData

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=sbiData%real1
    beta=sbiData%real2
    w = sbiData%real3
    CalFplusprimEnd = sbiData%real4

    primStar = sbi_efold_primitive(x,alpha,beta)
    epsOneStar = sbi_epsilon_one(x,alpha,beta)
    potStar = sbi_norm_potential(x,alpha,beta)

    find_sbi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_sbi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function sbi_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: sbi_x_rrad
    real(kp), intent(in) :: alpha,beta,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sbiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = sbi_epsilon_one(xEnd,alpha,beta)
    potEnd = sbi_norm_potential(xEnd,alpha,beta)
    primEnd = sbi_efold_primitive(xEnd,alpha,beta)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    sbiData%real1 = alpha 
    sbiData%real2 = beta
    sbiData%real3 = calF + primEnd

    maxi = xEnd*(1._kp-epsilon(1._kp))
    mini = epsilon(1._kp)

    x = zbrent(find_sbi_x_rrad,mini,maxi,tolzbrent,sbiData)
    sbi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (sbi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function sbi_x_rrad

  function find_sbi_x_rrad(x,sbiData)   
    implicit none
    real(kp) :: find_sbi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sbiData

    real(kp) :: primStar,alpha,beta,CalFplusprimEnd,potStar,epsOneStar

    alpha=sbiData%real1
    beta=sbiData%real2
    CalFplusprimEnd = sbiData%real3

    primStar = sbi_efold_primitive(x,alpha,beta)
    epsOneStar = sbi_epsilon_one(x,alpha,beta)
    potStar = sbi_norm_potential(x,alpha,beta)

    find_sbi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_sbi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function sbi_x_rreh(alpha,beta,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: sbi_x_rreh
    real(kp), intent(in) :: alpha,beta,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sbiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = sbi_epsilon_one(xEnd,alpha,beta)
    potEnd = sbi_norm_potential(xEnd,alpha,beta)
    primEnd = sbi_efold_primitive(xEnd,alpha,beta)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    sbiData%real1 = alpha 
    sbiData%real2 = beta
    sbiData%real3 = calF + primEnd

    maxi = xEnd*(1._kp-epsilon(1._kp))
    mini = epsilon(1._kp)

    x = zbrent(find_sbi_x_rreh,mini,maxi,tolzbrent,sbiData)
    sbi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (sbi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function sbi_x_rreh

  function find_sbi_x_rreh(x,sbiData)   
    implicit none
    real(kp) :: find_sbi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sbiData

    real(kp) :: primStar,alpha,beta,CalFplusprimEnd,potStar

    alpha=sbiData%real1
    beta=sbiData%real2
    CalFplusprimEnd = sbiData%real3

    primStar = sbi_efold_primitive(x,alpha,beta)   
    potStar = sbi_norm_potential(x,alpha,beta)

    find_sbi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_sbi_x_rreh



  function sbi_lnrhoreh_max(alpha,beta,xend,Pstar) 
    implicit none
    real(kp) :: sbi_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = sbi_norm_potential(xEnd,alpha,beta)
    epsOneEnd = sbi_epsilon_one(xEnd,alpha,beta)

!   Trick to return x such that rho_reh=rho_end

    x = sbi_x_star(alpha,beta,xend,wrad,junk,Pstar)    
    potStar = sbi_norm_potential(x,alpha,beta)
    epsOneStar = sbi_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'sbi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    sbi_lnrhoreh_max = lnRhoEnd

  end function sbi_lnrhoreh_max

  
end module sbireheat
