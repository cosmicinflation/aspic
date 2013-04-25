!spontaneous symmetry breaking 5 reheating functions in the slow-roll approximations

module ssbi5reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssbi5sr, only : ssbi5_epsilon_one, ssbi5_epsilon_two, ssbi5_epsilon_three
  use ssbi5sr, only : ssbi5_norm_potential
  use ssbi5sr, only : ssbi5_x_endinf, ssbi5_efold_primitive
  use cosmopar, only : QrmsOverT

  implicit none

  private

  public ssbi5_x_star, ssbi5_lnrhoreh_max

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi5_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssbi5_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssbi5Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssbi5_x_endinf(alpha,beta)
    epsOneEnd = ssbi5_epsilon_one(xEnd,alpha,beta)
    potEnd = ssbi5_norm_potential(xEnd,alpha,beta)

    primEnd = ssbi5_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssbi5Data%real1 = alpha
    ssbi5Data%real2 = beta
    ssbi5Data%real3 = w
    ssbi5Data%real4 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = ssbi5_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))

    x = zbrent(find_ssbi5_x_star,mini,maxi,tolzbrent,ssbi5Data)
    ssbi5_x_star = x

!   print*,'ssbi5_x_star:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd, &
!       '   primEnd=',primEnd,'   xstar=',x
!    pause

    if (present(bfoldstar)) then
       bfoldstar = - (ssbi5_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssbi5_x_star

  function find_ssbi5_x_star(x,ssbi5Data)   
    implicit none
    real(kp) :: find_ssbi5_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssbi5Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssbi5Data%real1
    beta=ssbi5Data%real2
    w = ssbi5Data%real3
    CalFplusprimEnd = ssbi5Data%real4

    primStar = ssbi5_efold_primitive(x,alpha,beta)
    epsOneStar = ssbi5_epsilon_one(x,alpha,beta)
    potStar = ssbi5_norm_potential(x,alpha,beta)

    find_ssbi5_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssbi5_x_star



  function ssbi5_lnrhoreh_max(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssbi5_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssbi5_x_endinf(alpha,beta)

    potEnd  = ssbi5_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssbi5_epsilon_one(xEnd,alpha,beta)

!   print*,'ssbi5_lnrhoreh_max:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssbi5_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssbi5_norm_potential(x,alpha,beta)
    epsOneStar = ssbi5_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssbi5_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssbi5_lnrhoreh_max = lnRhoEnd

  end function ssbi5_lnrhoreh_max

  
end module ssbi5reheat
