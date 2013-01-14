!supergravity brane inflation reheating functions in the slow-roll approximations

module sbireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use sbisr, only : sbi_epsilon_one, sbi_epsilon_two, sbi_epsilon_three
  use sbisr, only : sbi_norm_potential
  use sbisr, only : sbi_x_endinf, sbi_efold_primitive
  implicit none

  private

  public sbi_x_star, sbi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function sbi_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: sbi_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: sbiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = sbi_x_endinf(alpha,beta)
    epsOneEnd = sbi_epsilon_one(xEnd,alpha,beta)
    potEnd = sbi_norm_potential(xEnd,alpha,beta)
    primEnd = sbi_efold_primitive(xEnd,alpha,beta)


    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    sbiData%real1 = alpha 
    sbiData%real2 = beta
    sbiData%real3 = xEnd
    sbiData%real4 = w
    sbiData%real5 = calF + primEnd

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

    real(kp) :: primStar,alpha,beta,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=sbiData%real1
    beta=sbiData%real2
    xEnd=sbiData%real3
    w = sbiData%real4
    CalFplusprimEnd = sbiData%real5

    primStar = sbi_efold_primitive(x,alpha,beta)
    epsOneStar = sbi_epsilon_one(x,alpha,beta)
    potStar = sbi_norm_potential(x,alpha,beta)

    find_sbi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_sbi_x_star



  function sbi_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: sbi_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = sbi_x_endinf(alpha,beta)
    potEnd  = sbi_norm_potential(xEnd,alpha,beta)
    epsOneEnd = sbi_epsilon_one(xEnd,alpha,beta)

!   Trick to return x such that rho_reh=rho_end

    x = sbi_x_star(alpha,beta,wrad,junk,Pstar)    
    potStar = sbi_norm_potential(x,alpha,beta)
    epsOneStar = sbi_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'sbi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    sbi_lnrhoend = lnRhoEnd

  end function sbi_lnrhoend

  
end module sbireheat
