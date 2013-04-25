!spontaneous symmetry breaking 3 reheating functions in the slow-roll approximations

module ssbi3reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssbi3sr, only : ssbi3_epsilon_one, ssbi3_epsilon_two, ssbi3_epsilon_three
  use ssbi3sr, only : ssbi3_norm_potential, ssbi3_x_potmax
  use ssbi3sr, only : ssbi3_x_endinf, ssbi3_efold_primitive  
  use cosmopar, only : QrmsOverT

  implicit none

  private

  public ssbi3_x_star, ssbi3_lnrhoreh_max

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi3_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssbi3_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssbi3Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssbi3_x_endinf(alpha,beta)
    epsOneEnd = ssbi3_epsilon_one(xEnd,alpha,beta)
    potEnd = ssbi3_norm_potential(xEnd,alpha,beta)

    primEnd = ssbi3_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssbi3Data%real1 = alpha
    ssbi3Data%real2 = beta
    ssbi3Data%real3 = w
    ssbi3Data%real4 = calF + primEnd

    mini = ssbi3_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssbi3_x_potmax(alpha,beta)*(1._kp-epsilon(1._kp))

    x = zbrent(find_ssbi3_x_star,mini,maxi,tolzbrent,ssbi3Data)
    ssbi3_x_star = x

!   print*,'ssbi3_x_star:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd, &
!       '   primEnd=',primEnd,'   xstar=',x
!    pause

    if (present(bfoldstar)) then
       bfoldstar = - (ssbi3_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssbi3_x_star

  function find_ssbi3_x_star(x,ssbi3Data)   
    implicit none
    real(kp) :: find_ssbi3_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssbi3Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssbi3Data%real1
    beta=ssbi3Data%real2
    w = ssbi3Data%real3
    CalFplusprimEnd = ssbi3Data%real4

    primStar = ssbi3_efold_primitive(x,alpha,beta)
    epsOneStar = ssbi3_epsilon_one(x,alpha,beta)
    potStar = ssbi3_norm_potential(x,alpha,beta)

    find_ssbi3_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssbi3_x_star



  function ssbi3_lnrhoreh_max(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssbi3_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssbi3_x_endinf(alpha,beta)

    potEnd  = ssbi3_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssbi3_epsilon_one(xEnd,alpha,beta)

!   print*,'ssbi3_lnrhoreh_max:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssbi3_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssbi3_norm_potential(x,alpha,beta)
    epsOneStar = ssbi3_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssbi3_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssbi3_lnrhoreh_max = lnRhoEnd

  end function ssbi3_lnrhoreh_max

  
end module ssbi3reheat
