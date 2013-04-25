!kälher moduli inflation reheating functions in the slow-roll approximations

module kmiireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use kmiisr, only : kmii_epsilon_one, kmii_epsilon_two, kmii_epsilon_three
  use kmiisr, only : kmii_norm_potential
  use kmiisr, only : kmii_x_endinf, kmii_efold_primitive
  implicit none

  private

  public kmii_x_star, kmii_lnrhoreh_max 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function kmii_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: kmii_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: kmiiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = kmii_x_endinf(alpha)
    epsOneEnd = kmii_epsilon_one(xEnd,alpha)
    potEnd = kmii_norm_potential(xEnd,alpha)
    primEnd = kmii_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    kmiiData%real1 = alpha    
    kmiiData%real2 = w
    kmiiData%real3 = calF + primEnd

    !Assuming that inflation proceeds from the right to the left,
    !from initial high values of the field compared with the Planck mass
    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)
  
    x = zbrent(find_kmii_x_star,mini,maxi,tolzbrent,kmiiData)
    kmii_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (kmii_efold_primitive(x,alpha) - primEnd)
    endif

  end function kmii_x_star

  function find_kmii_x_star(x,kmiiData)   
    implicit none
    real(kp) :: find_kmii_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kmiiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=kmiiData%real1
    w = kmiiData%real2
    CalFplusprimEnd = kmiiData%real3

    primStar = kmii_efold_primitive(x,alpha)
    epsOneStar = kmii_epsilon_one(x,alpha)
    potStar = kmii_norm_potential(x,alpha)

    find_kmii_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_kmii_x_star



  function kmii_lnrhoreh_max(alpha,Pstar) 
    implicit none
    real(kp) :: kmii_lnrhoreh_max
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = kmii_x_endinf(alpha)
    potEnd  = kmii_norm_potential(xEnd,alpha)
    epsOneEnd = kmii_epsilon_one(xEnd,alpha)

   print*, 'xEnd=',xEnd,'potEnd=',potEnd,'epsOneEnd=',epsOneEnd

!   Trick to return x such that rho_reh=rho_end

    x = kmii_x_star(alpha,wrad,junk,Pstar)    
    potStar = kmii_norm_potential(x,alpha)
    epsOneStar = kmii_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'kmii_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    kmii_lnrhoreh_max = lnRhoEnd

  end function kmii_lnrhoreh_max

  
end module kmiireheat
