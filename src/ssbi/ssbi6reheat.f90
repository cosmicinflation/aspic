!spontaneous symmetry breaking 6 reheating functions in the slow-roll approximations

module ssbi6reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssbi6sr, only : ssbi6_epsilon_one, ssbi6_epsilon_two, ssbi6_epsilon_three
  use ssbi6sr, only : ssbi6_norm_potential
  use ssbi6sr, only : ssbi6_x_endinf, ssbi6_efold_primitive
  use cosmopar, only : QrmsOverT

  implicit none

  private

  public ssbi6_x_star, ssbi6_lnrhoend

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi6_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssbi6_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssbi6Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssbi6_x_endinf(alpha,beta)
    epsOneEnd = ssbi6_epsilon_one(xEnd,alpha,beta)
    potEnd = ssbi6_norm_potential(xEnd,alpha,beta)

    primEnd = ssbi6_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssbi6Data%real1 = alpha
    ssbi6Data%real2 = beta
    ssbi6Data%real3 = w
    ssbi6Data%real4 = calF + primEnd

    mini = ssbi6_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = 10._kp**(6._kp)*mini

    x = zbrent(find_ssbi6_x_star,mini,maxi,tolzbrent,ssbi6Data)
    ssbi6_x_star = x

!   print*,'ssbi6_x_star:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd, &
!       '   primEnd=',primEnd,'   xstar=',x
!    pause

    if (present(bfoldstar)) then
       bfoldstar = - (ssbi6_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssbi6_x_star

  function find_ssbi6_x_star(x,ssbi6Data)   
    implicit none
    real(kp) :: find_ssbi6_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssbi6Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssbi6Data%real1
    beta=ssbi6Data%real2
    w = ssbi6Data%real3
    CalFplusprimEnd = ssbi6Data%real4

    primStar = ssbi6_efold_primitive(x,alpha,beta)
    epsOneStar = ssbi6_epsilon_one(x,alpha,beta)
    potStar = ssbi6_norm_potential(x,alpha,beta)

    find_ssbi6_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssbi6_x_star



  function ssbi6_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssbi6_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssbi6_x_endinf(alpha,beta)

    potEnd  = ssbi6_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssbi6_epsilon_one(xEnd,alpha,beta)

!   print*,'ssbi6_lnrhoend:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssbi6_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssbi6_norm_potential(x,alpha,beta)
    epsOneStar = ssbi6_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssbi6_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssbi6_lnrhoend = lnRhoEnd

  end function ssbi6_lnrhoend

  
end module ssbi6reheat
