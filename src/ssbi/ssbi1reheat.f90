!spontaneous symmetry breaking 1 reheating functions in the slow-roll approximations

module ssbi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssbi1sr, only : ssbi1_epsilon_one, ssbi1_epsilon_two, ssbi1_epsilon_three
  use ssbi1sr, only : ssbi1_norm_potential
  use ssbi1sr, only : ssbi1_x_endinf, ssbi1_efold_primitive
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public ssbi1_x_star, ssbi1_lnrhoreh_max

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi1_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssbi1_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssbi1Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssbi1_x_endinf(alpha,beta)
    epsOneEnd = ssbi1_epsilon_one(xEnd,alpha,beta)
    potEnd = ssbi1_norm_potential(xEnd,alpha,beta)

    primEnd = ssbi1_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssbi1Data%real1 = alpha
    ssbi1Data%real2 = beta
    ssbi1Data%real3 = w
    ssbi1Data%real4 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)

    x = zbrent(find_ssbi1_x_star,mini,maxi,tolzbrent,ssbi1Data)
    ssbi1_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ssbi1_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssbi1_x_star

  function find_ssbi1_x_star(x,ssbi1Data)   
    implicit none
    real(kp) :: find_ssbi1_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssbi1Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssbi1Data%real1
    beta=ssbi1Data%real2
    w = ssbi1Data%real3
    CalFplusprimEnd = ssbi1Data%real4

    primStar = ssbi1_efold_primitive(x,alpha,beta)
    epsOneStar = ssbi1_epsilon_one(x,alpha,beta)
    potStar = ssbi1_norm_potential(x,alpha,beta)

    find_ssbi1_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssbi1_x_star



  function ssbi1_lnrhoreh_max(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssbi1_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssbi1_x_endinf(alpha,beta)

    potEnd  = ssbi1_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssbi1_epsilon_one(xEnd,alpha,beta)

!    print*,'ssbi1_lnrhoreh_max:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssbi1_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssbi1_norm_potential(x,alpha,beta)
    epsOneStar = ssbi1_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssbi1_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssbi1_lnrhoreh_max = lnRhoEnd

  end function ssbi1_lnrhoreh_max

  
end module ssbi1reheat
