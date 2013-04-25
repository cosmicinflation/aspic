!spontaneous symmetry breaking 2 reheating functions in the slow-roll approximations

module ssbi2reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssbi2sr, only : ssbi2_epsilon_one, ssbi2_epsilon_two, ssbi2_epsilon_three
  use ssbi2sr, only : ssbi2_norm_potential
  use ssbi2sr, only : ssbi2_x_endinf, ssbi2_efold_primitive
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public ssbi2_x_star, ssbi2_lnrhoreh_max

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi2_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssbi2_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssbi2Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssbi2_x_endinf(alpha,beta)
    epsOneEnd = ssbi2_epsilon_one(xEnd,alpha,beta)
    potEnd = ssbi2_norm_potential(xEnd,alpha,beta)

    primEnd = ssbi2_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssbi2Data%real1 = alpha
    ssbi2Data%real2 = beta
    ssbi2Data%real3 = w
    ssbi2Data%real4 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = ssbi2_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))

    x = zbrent(find_ssbi2_x_star,mini,maxi,tolzbrent,ssbi2Data)
    ssbi2_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ssbi2_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssbi2_x_star

  function find_ssbi2_x_star(x,ssbi2Data)   
    implicit none
    real(kp) :: find_ssbi2_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssbi2Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssbi2Data%real1
    beta=ssbi2Data%real2
    w = ssbi2Data%real3
    CalFplusprimEnd = ssbi2Data%real4

    primStar = ssbi2_efold_primitive(x,alpha,beta)
    epsOneStar = ssbi2_epsilon_one(x,alpha,beta)
    potStar = ssbi2_norm_potential(x,alpha,beta)

    find_ssbi2_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssbi2_x_star



  function ssbi2_lnrhoreh_max(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssbi2_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssbi2_x_endinf(alpha,beta)

    potEnd  = ssbi2_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssbi2_epsilon_one(xEnd,alpha,beta)

!    print*,'ssbi2_lnrhoreh_max:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssbi2_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssbi2_norm_potential(x,alpha,beta)
    epsOneStar = ssbi2_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssbi2_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssbi2_lnrhoreh_max = lnRhoEnd

  end function ssbi2_lnrhoreh_max

  
end module ssbi2reheat
