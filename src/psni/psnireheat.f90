!pseudo natural inflation reheating functions in the slow-roll approximations

module psnireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use psnisr, only : psni_epsilon_one, psni_epsilon_two, psni_epsilon_three
  use psnisr, only : psni_norm_potential
  use psnisr, only : psni_x_endinf, psni_efold_primitive
  implicit none

  private

  public psni_x_star, psni_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function psni_x_star(alpha,mu,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: psni_x_star
    real(kp), intent(in) :: alpha,mu,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: psniData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = psni_x_endinf(alpha,mu)
    epsOneEnd = psni_epsilon_one(xEnd,alpha,mu)
    potEnd = psni_norm_potential(xEnd,alpha)
    primEnd = psni_efold_primitive(xEnd,alpha,mu)

!    print*,'psni_x_star:    xEnd=',xEnd,'  epsOneEnd=',epsOneEnd,'  potEnd=',potEnd,'  primEnd=',primEnd
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    psniData%real1 = alpha 
    psniData%real2 = mu
    psniData%real3 = xEnd
    psniData%real4 = w
    psniData%real5 = calF + primEnd

    mini=epsilon(1._kp)
    maxi = xEnd


    x = zbrent(find_psni_x_star,mini,maxi,tolzbrent,psniData)
    psni_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (psni_efold_primitive(x,alpha,mu) - primEnd)
    endif

  end function psni_x_star

  function find_psni_x_star(x,psniData)   
    implicit none
    real(kp) :: find_psni_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: psniData

    real(kp) :: primStar,alpha,mu,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=psniData%real1
    mu=psniData%real2
    xEnd=psniData%real3
    w = psniData%real4
    CalFplusprimEnd = psniData%real5

    primStar = psni_efold_primitive(x,alpha,mu)
    epsOneStar = psni_epsilon_one(x,alpha,mu)
    potStar = psni_norm_potential(x,alpha)

    find_psni_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_psni_x_star



  function psni_lnrhoend(alpha,mu,Pstar) 
    implicit none
    real(kp) :: psni_lnrhoend
    real(kp), intent(in) :: alpha,mu,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = psni_x_endinf(alpha,mu)
    potEnd  = psni_norm_potential(xEnd,alpha)
    epsOneEnd = psni_epsilon_one(xEnd,alpha,mu)

!   Trick to return x such that rho_reh=rho_end

    x = psni_x_star(alpha,mu,wrad,junk,Pstar)    
    potStar = psni_norm_potential(x,alpha)
    epsOneStar = psni_epsilon_one(x,alpha,mu)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'psni_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    psni_lnrhoend = lnRhoEnd

  end function psni_lnrhoend

  
end module psnireheat
