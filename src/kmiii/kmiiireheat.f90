!Kähler moduli inflation II reheating functions in the slow-roll approximations

module kmiiireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use kmiiisr, only : kmiii_epsilon_one, kmiii_epsilon_two, kmiii_epsilon_three
  use kmiiisr, only : kmiii_norm_potential
  use kmiiisr, only : kmiii_x_endinf, kmiii_efold_primitive
  implicit none

  private

  public kmiii_x_star, kmiii_lnrhoreh_max 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function kmiii_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: kmiii_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: kmiiiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = kmiii_x_endinf(alpha,beta)
    epsOneEnd = kmiii_epsilon_one(xEnd,alpha,beta)
    potEnd = kmiii_norm_potential(xEnd,alpha,beta)
    primEnd = kmiii_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    kmiiiData%real1 = alpha    
    kmiiiData%real2 = beta  
    kmiiiData%real3 = w
    kmiiiData%real4 = calF + primEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_kmiii_x_star,mini,maxi,tolzbrent,kmiiiData)
    kmiii_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (kmiii_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function kmiii_x_star

  function find_kmiii_x_star(x,kmiiiData)   
    implicit none
    real(kp) :: find_kmiii_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kmiiiData

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=kmiiiData%real1
    beta=kmiiiData%real2
    w = kmiiiData%real3
    CalFplusprimEnd = kmiiiData%real4

    primStar = kmiii_efold_primitive(x,alpha,beta)
    epsOneStar = kmiii_epsilon_one(x,alpha,beta)
    potStar = kmiii_norm_potential(x,alpha,beta)

    find_kmiii_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_kmiii_x_star



  function kmiii_lnrhoreh_max(alpha,beta,Pstar) 
    implicit none
    real(kp) :: kmiii_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = kmiii_x_endinf(alpha,beta)
    potEnd  = kmiii_norm_potential(xEnd,alpha,beta)
    epsOneEnd = kmiii_epsilon_one(xEnd,alpha,beta)


!   Trick to return x such that rho_reh=rho_end

    x = kmiii_x_star(alpha,beta,wrad,junk,Pstar)    
    potStar = kmiii_norm_potential(x,alpha,beta)
    epsOneStar = kmiii_epsilon_one(x,alpha,beta)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'kmiii_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    kmiii_lnrhoreh_max = lnRhoEnd

  end function kmiii_lnrhoreh_max

  
end module kmiiireheat
