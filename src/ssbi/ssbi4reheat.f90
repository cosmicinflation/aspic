!spontaneous symmetry breaking 4 reheating functions in the slow-roll approximations

module ssbi4reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssbi4sr, only : ssbi4_epsilon_one, ssbi4_epsilon_two, ssbi4_epsilon_three
  use ssbi4sr, only : ssbi4_norm_potential, ssbi4_x_potmax
  use ssbi4sr, only : ssbi4_x_endinf, ssbi4_efold_primitive
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public ssbi4_x_star, ssbi4_lnrhoreh_max

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi4_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssbi4_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssbi4Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssbi4_x_endinf(alpha,beta)
    epsOneEnd = ssbi4_epsilon_one(xEnd,alpha,beta)
    potEnd = ssbi4_norm_potential(xEnd,alpha,beta)

    primEnd = ssbi4_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssbi4Data%real1 = alpha
    ssbi4Data%real2 = beta
    ssbi4Data%real3 = w
    ssbi4Data%real4 = calF + primEnd


    mini = ssbi4_x_potmax(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssbi4_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))

    x = zbrent(find_ssbi4_x_star,mini,maxi,tolzbrent,ssbi4Data)
    ssbi4_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ssbi4_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssbi4_x_star

  function find_ssbi4_x_star(x,ssbi4Data)   
    implicit none
    real(kp) :: find_ssbi4_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssbi4Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssbi4Data%real1
    beta=ssbi4Data%real2
    w = ssbi4Data%real3
    CalFplusprimEnd = ssbi4Data%real4

    primStar = ssbi4_efold_primitive(x,alpha,beta)
    epsOneStar = ssbi4_epsilon_one(x,alpha,beta)
    potStar = ssbi4_norm_potential(x,alpha,beta)

    find_ssbi4_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssbi4_x_star



  function ssbi4_lnrhoreh_max(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssbi4_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssbi4_x_endinf(alpha,beta)

    potEnd  = ssbi4_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssbi4_epsilon_one(xEnd,alpha,beta)


!   Trick to return x such that rho_reh=rho_end

    x = ssbi4_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssbi4_norm_potential(x,alpha,beta)
    epsOneStar = ssbi4_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssbi4_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssbi4_lnrhoreh_max = lnRhoEnd

  end function ssbi4_lnrhoreh_max

  
end module ssbi4reheat
