!Logamediate inflation 1 reheating functions in the slow-roll approximations

module lmi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lmi1sr, only : lmi1_epsilon_one, lmi1_epsilon_two, lmi1_epsilon_three
  use lmi1sr, only : lmi1_norm_potential
  use lmi1sr, only : lmi1_x_endinf, lmi1_efold_primitive
  use lmi1sr, only : lmi1_alpha, lmi1_beta
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public lmi1_x_star, lmi1_lnrhoend,lmi1_M

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lmi1_x_star(gamma_lmi,M,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lmi1_x_star
    real(kp), intent(in) :: gamma_lmi,M,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: lmi1Data

    real(kp) ::alpha,beta
    alpha=lmi1_alpha(gamma_lmi,M)
    beta=lmi1_beta(gamma_lmi,M)

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lmi1_x_endinf(gamma_lmi,M)

    epsOneEnd = lmi1_epsilon_one(xEnd,gamma_lmi,M)
    potEnd = lmi1_norm_potential(xEnd,gamma_lmi,M)

    primEnd = lmi1_efold_primitive(xEnd,gamma_lmi,M)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    lmi1Data%real1 = gamma_lmi
    lmi1Data%real2 = M
    lmi1Data%real3 = w
    lmi1Data%real4 = calF + primEnd

    maxi = (alpha/(beta*gamma_lmi))**(1._kp/gamma_lmi)*(1._kp-100._kp*epsilon(1._kp))
    mini = lmi1_x_endinf(gamma_lmi,M)*(1._kp+100._kp*epsilon(1._kp))

    x = zbrent(find_lmi1_x_star,mini,maxi,tolzbrent,lmi1Data)
    lmi1_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi1_efold_primitive(x,gamma_lmi,M) - primEnd)
    endif

  end function lmi1_x_star

  function find_lmi1_x_star(x,lmi1Data)   
    implicit none
    real(kp) :: find_lmi1_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi1Data

    real(kp) :: primStar,gamma_lmi,M,w,CalFplusprimEnd,potStar,epsOneStar

    gamma_lmi=lmi1Data%real1
    M=lmi1Data%real2
    w = lmi1Data%real3
    CalFplusprimEnd = lmi1Data%real4

    primStar = lmi1_efold_primitive(x,gamma_lmi,M)
    epsOneStar = lmi1_epsilon_one(x,gamma_lmi,M)
    potStar = lmi1_norm_potential(x,gamma_lmi,M)

    find_lmi1_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lmi1_x_star



  function lmi1_lnrhoend(gamma_lmi,M,Pstar) 
    implicit none
    real(kp) :: lmi1_lnrhoend
    real(kp), intent(in) :: gamma_lmi,M,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = lmi1_x_endinf(gamma_lmi,M)

    potEnd  = lmi1_norm_potential(xEnd,gamma_lmi,M)

    epsOneEnd = lmi1_epsilon_one(xEnd,gamma_lmi,M)


!   Trick to return x such that rho_reh=rho_end

    x = lmi1_x_star(gamma_lmi,M,wrad,junk,Pstar)  


    potStar = lmi1_norm_potential(x,gamma_lmi,M)
    epsOneStar = lmi1_epsilon_one(x,gamma_lmi,M)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'lmi1_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lmi1_lnrhoend = lnRhoEnd

  end function lmi1_lnrhoend



!Finds the value of M compatible with the COBE normalization, given gamma,lnRhoReh,w and Pstar
    function lmi1_M(gamma_lmi,lnRhoReh,w,Pstar)
    implicit none
    real(kp), intent(in) :: gamma_lmi,lnRhoReh,w,Pstar
    real(kp) :: lmi1_M
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi1Data

    maxi = 10000.
    mini = 0.0001

    lmi1Data%real1 = gamma_lmi
    lmi1Data%real2 = lnRhoReh
    lmi1Data%real3 = w
    lmi1Data%real4 = Pstar
    
    lmi1_M = zbrent(find_lmi1M,mini,maxi,tolFind,lmi1Data)


    end function lmi1_M

  function find_lmi1M(x,lmi1Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi1Data
    real(kp) :: find_lmi1M
    real(kp) :: gamma_lmi,lnRhoReh,w,Pstar,xstar,alpha,beta

    gamma_lmi = lmi1Data%real1
    lnRhoReh = lmi1Data%real2
    w = lmi1Data%real3
    Pstar = lmi1Data%real4
    xstar = lmi1_x_star(gamma_lmi,x,w,lnRhoReh,Pstar)
    alpha=lmi1_alpha(gamma_lmi,x)
    beta=lmi1_beta(gamma_lmi,x)
    
    find_lmi1M = 720._kp*acos(-1._kp)**2*(alpha-beta*gamma_lmi*xstar**gamma_lmi)**2* &
                exp(beta*xstar**gamma_lmi)*xstar**(-alpha-2._kp)*QrmsOverT**2-x**4
   
  end function find_lmi1M

  
end module lmi1reheat
