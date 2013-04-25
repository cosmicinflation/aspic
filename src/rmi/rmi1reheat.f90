!running mass 1 reheating functions in the slow-roll approximations

module rmi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rmi1sr, only : rmi1_epsilon_one, rmi1_epsilon_two, rmi1_epsilon_three
  use rmi1sr, only : rmi1_norm_potential, rmi1_efold_primitive
  use cosmopar, only : QrmsOverT

  implicit none

  private

  public rmi1_x_star, rmi1_lnrhoreh_max

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rmi1_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rmi1_x_star
    real(kp), intent(in) :: c,phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: rmi1Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = rmi1_epsilon_one(xEnd,c,phi0)
    potEnd = rmi1_norm_potential(xEnd,c,phi0)

    primEnd = rmi1_efold_primitive(xEnd,c,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rmi1Data%real1 = c
    rmi1Data%real2 = phi0
    rmi1Data%real3 = w
    rmi1Data%real4 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = 1._kp*(1._kp-epsilon(1._kp))

!    print*,'rmi1_x_star:   mini=',mini,'maxi=',maxi,'f(mini)=', &
!           find_reheat(rmi1_efold_primitive(mini,c,phi0),rmi1Data%real4,w, &
!           rmi1_epsilon_one(mini,c,phi0),rmi1_norm_potential(mini,c,phi0)), &
!           'f(maxi)=', find_reheat(rmi1_efold_primitive(maxi,c,phi0),rmi1Data%real4,w, &
!           rmi1_epsilon_one(maxi,c,phi0),rmi1_norm_potential(maxi,c,phi0))
!    print*,'rmi1_x_star:   efoldprimitive(mini)=',rmi1_efold_primitive(mini,c,phi0), &
!           'epsilonOne(mini)=',rmi1_epsilon_one(mini,c,phi0), &
!           'pot(mini)=',rmi1_norm_potential(mini,c,phi0)
!    print*,'rmi1_x_star:   efoldprimitive(maxi)=',rmi1_efold_primitive(maxi,c,phi0), &
!           'epsilonOne(maxi)=',rmi1_epsilon_one(maxi,c,phi0), &
!           'pot(maxi)=',rmi1_norm_potential(maxi,c,phi0)
!    print*,'primEnd=',primEnd,'calF+primEnd=',calF + primEnd

!    pause

    x = zbrent(find_rmi1_x_star,mini,maxi,tolzbrent,rmi1Data)
    rmi1_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rmi1_efold_primitive(x,c,phi0) - primEnd)
    endif

  end function rmi1_x_star

  function find_rmi1_x_star(x,rmi1Data)   
    implicit none
    real(kp) :: find_rmi1_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rmi1Data

    real(kp) :: primStar,c,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    c=rmi1Data%real1
    phi0=rmi1Data%real2
    w = rmi1Data%real3
    CalFplusprimEnd = rmi1Data%real4

    primStar = rmi1_efold_primitive(x,c,phi0)
    epsOneStar = rmi1_epsilon_one(x,c,phi0)
    potStar = rmi1_norm_potential(x,c,phi0)

    find_rmi1_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rmi1_x_star



  function rmi1_lnrhoreh_max(c,phi0,xend,Pstar) 
    implicit none
    real(kp) :: rmi1_lnrhoreh_max
    real(kp), intent(in) :: c, phi0, xend, Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd

    potEnd  = rmi1_norm_potential(xEnd,c,phi0)

    epsOneEnd = rmi1_epsilon_one(xEnd,c,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = rmi1_x_star(c,phi0,xend,wrad,junk,Pstar)  


    potStar = rmi1_norm_potential(x,c,phi0)
    epsOneStar = rmi1_epsilon_one(x,c,phi0)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'rmi1_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rmi1_lnrhoreh_max = lnRhoEnd

  end function rmi1_lnrhoreh_max

  
end module rmi1reheat
