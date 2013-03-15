!radiatively corrected quartic inflation reheating functions in the slow-roll approximations

module rchireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rchisr, only : rchi_epsilon_one, rchi_epsilon_two, rchi_epsilon_three
  use rchisr, only : rchi_norm_potential
  use rchisr, only : rchi_x_endinf, rchi_efold_primitive
  implicit none

  private

  public rchi_x_star, rchi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rchi_x_star(AI,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rchi_x_star
    real(kp), intent(in) :: AI,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: rchiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = rchi_x_endinf(AI)
    epsOneEnd = rchi_epsilon_one(xEnd,AI)
    potEnd = rchi_norm_potential(xEnd,AI)
    primEnd = rchi_efold_primitive(xEnd,AI)

!   print*,'xEnd=',xEnd,'epsOneEnd=',epsOneEnd,'potEnd=',potEnd,'primEnd=',primEnd,'AI=',AI
!   pause
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rchiData%real1 = AI    
    rchiData%real2 = w
    rchiData%real3 = calF + primEnd

    mini=xend*(1._kp+epsilon(1._kp))
    if (AI .lt. 0._kp .and. abs(AI) .lt. 64._kp*acos(-1._kp)**2) then !maximum of the potential
      maxi=min(abs(sqrt(3._kp/2._kp)*log(abs(AI)/(64._kp*acos(-1._kp)**2))),xend+100._kp)
    else
      maxi = max(abs(AI)*1000._kp/(48._kp*acos(-1._kp)**2)*sqrt(3._kp/2._kp)+xend,xend+10._kp,100._kp)
    endif

    x = zbrent(find_rchi_x_star,mini,maxi,tolzbrent,rchiData)
    rchi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rchi_efold_primitive(x,AI) - primEnd)
    endif

!    print*,'xstar=',x,'epsOneStar=',rchi_epsilon_one(x,AI),'primStar=',rchi_efold_primitive(x,AI),'AI=',AI
!    pause

  end function rchi_x_star

  function find_rchi_x_star(x,rchiData)   
    implicit none
    real(kp) :: find_rchi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rchiData

    real(kp) :: primStar,AI,w,CalFplusprimEnd,potStar,epsOneStar

    AI=rchiData%real1
    w = rchiData%real2
    CalFplusprimEnd = rchiData%real3

    primStar = rchi_efold_primitive(x,AI)
    epsOneStar = rchi_epsilon_one(x,AI)
    potStar = rchi_norm_potential(x,AI)

    find_rchi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rchi_x_star



  function rchi_lnrhoend(AI,Pstar) 
    implicit none
    real(kp) :: rchi_lnrhoend
    real(kp), intent(in) :: AI,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = rchi_x_endinf(AI)
    potEnd  = rchi_norm_potential(xEnd,AI)
    epsOneEnd = rchi_epsilon_one(xEnd,AI)

!   Trick to return x such that rho_reh=rho_end

    x = rchi_x_star(AI,wrad,junk,Pstar)    
    potStar = rchi_norm_potential(x,AI)
    epsOneStar = rchi_epsilon_one(x,AI)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'rchi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rchi_lnrhoend = lnRhoEnd

  end function rchi_lnrhoend

  
end module rchireheat
