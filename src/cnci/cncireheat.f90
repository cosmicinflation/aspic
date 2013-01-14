!constant ns C model reheating functions in the slow-roll approximations

module cncireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use cncisr, only : cnci_epsilon_one, cnci_epsilon_two, cnci_epsilon_three
  use cncisr, only : cnci_norm_potential,cnci_efold_primitive,cnci_x_epsOne_equals_one
  implicit none

  private

  public cnci_x_star, cnci_lnrhoend,find_cnci_x_star

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function cnci_x_star(alpha,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: cnci_x_star
    real(kp), intent(in) :: alpha,xEnd,w,lnRhoReh,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: cnciData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = cnci_epsilon_one(xEnd,alpha)
    potEnd = cnci_norm_potential(xEnd,alpha)
    primEnd = cnci_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    cnciData%real1 = alpha
    cnciData%real2 = w
    cnciData%real3 = calF + primEnd

    maxi = xend*(1._kp-epsilon(1._kp))
    mini = cnci_x_epsOne_equals_one(alpha)    !Value of x such that epsilon1=1 below which inflation cannot proceed

    x = zbrent(find_cnci_x_star,mini,maxi,tolFind,cnciData)
    cnci_x_star = x

    if (present(bfold)) then
       bfold = -(cnci_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnci_x_star

  function find_cnci_x_star(x,cnciData)   
    implicit none
    real(kp) :: find_cnci_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnciData

    real(kp) :: primStar,alpha,xend,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=cnciData%real1
    w = cnciData%real2
    CalFplusPrimEnd = cnciData%real3

    primStar = cnci_efold_primitive(x,alpha)
    epsOneStar = cnci_epsilon_one(x,alpha)
    potStar = cnci_norm_potential(x,alpha)


    find_cnci_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)


  end function find_cnci_x_star


  function cnci_lnrhoend(alpha,xEnd,Pstar) 
    implicit none
    real(kp) :: cnci_lnrhoend
    real(kp), intent(in) :: alpha,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
         
    potEnd  = cnci_norm_potential(xEnd,alpha)
    epsOneEnd = cnci_epsilon_one(xEnd,alpha)
       
    x = cnci_x_star(alpha,xEnd,wrad,junk,Pstar)    
    potStar = cnci_norm_potential(x,alpha)
    epsOneStar = cnci_epsilon_one(x,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'cnci_lnrhoend: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    cnci_lnrhoend = lnRhoEnd

  end function cnci_lnrhoend

  
  
end module cncireheat
