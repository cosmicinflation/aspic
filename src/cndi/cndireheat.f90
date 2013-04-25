!non canonical Kahler inflation reheating functions in the slow-roll approximations

module cndireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use cndisr, only : cndi_epsilon_one, cndi_epsilon_two, cndi_epsilon_three
  use cndisr, only : cndi_norm_potential, cndi_efold_primitive, cndi_xinimax
  implicit none

  private

  public cndi_x_star, cndi_lnrhoreh_max 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function cndi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: cndi_x_star
    real(kp), intent(in) :: alpha,beta,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: cndiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = cndi_epsilon_one(xEnd,alpha,beta)
    potEnd = cndi_norm_potential(xEnd,alpha,beta)
    primEnd = cndi_efold_primitive(xEnd,alpha,beta)


    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    cndiData%real1 = alpha 
    cndiData%real2 = beta
    cndiData%real3 = xEnd
    cndiData%real4 = w
    cndiData%real5 = calF + primEnd

    maxi = cndi_xinimax(alpha,beta)
    mini = xend

    x = zbrent(find_cndi_x_star,mini,maxi,tolzbrent,cndiData)
    cndi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (cndi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function cndi_x_star

  function find_cndi_x_star(x,cndiData)   
    implicit none
    real(kp) :: find_cndi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cndiData

    real(kp) :: primStar,alpha,beta,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=cndiData%real1
    beta=cndiData%real2
    xEnd=cndiData%real3
    w = cndiData%real4
    CalFplusprimEnd = cndiData%real5

    primStar = cndi_efold_primitive(x,alpha,beta)
    epsOneStar = cndi_epsilon_one(x,alpha,beta)
    potStar = cndi_norm_potential(x,alpha,beta)

    find_cndi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_cndi_x_star



  function cndi_lnrhoreh_max(alpha,beta,xend,Pstar) 
    implicit none
    real(kp) :: cndi_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    potEnd  = cndi_norm_potential(xEnd,alpha,beta)
    epsOneEnd = cndi_epsilon_one(xEnd,alpha,beta)
  

!   Trick to return x such that rho_reh=rho_end

    x = cndi_x_star(alpha,beta,xend,wrad,junk,Pstar)    
    potStar = cndi_norm_potential(x,alpha,beta)
    epsOneStar = cndi_epsilon_one(x,alpha,beta)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'cndi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    cndi_lnrhoreh_max = lnRhoEnd

  end function cndi_lnrhoreh_max

  
end module cndireheat
