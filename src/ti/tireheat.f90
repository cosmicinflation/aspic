!tip inflation reheating functions in the slow-roll approximations

module tireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use tisr, only : ti_epsilon_one, ti_epsilon_two, ti_epsilon_three
  use tisr, only : ti_norm_potential,ti_efold_primitive,ti_x_endinf
  implicit none

  private

  public ti_x_star, ti_lnrhoend,find_ti_x_star

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function ti_x_star(alpha,mu,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: ti_x_star
    real(kp), intent(in) :: alpha,mu,w,lnRhoReh,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd

    type(transfert) :: tiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ti_x_endinf(alpha,mu)
    
    epsOneEnd = ti_epsilon_one(xEnd,alpha,mu)
    potEnd = ti_norm_potential(xEnd,alpha)
    primEnd = ti_efold_primitive(xEnd,alpha,mu)
   

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    tiData%real1 = alpha
    tiData%real2 = mu
    tiData%real3 = w
    tiData%real4 = calF + primEnd

    maxi = xEnd*(1._kp-epsilon(1._kp))
    mini=maxi/10000._kp !To avoid numerical infinity


    x = zbrent(find_ti_x_star,mini,maxi,tolFind,tiData)
    ti_x_star = x

    if (present(bfold)) then
       bfold = -(ti_efold_primitive(x,alpha,mu) - primEnd)
    endif


  end function ti_x_star

  function find_ti_x_star(x,tiData)   
    implicit none
    real(kp) :: find_ti_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: tiData

    real(kp) :: primStar,alpha,mu,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=tiData%real1
    mu=tiData%real2
    w = tiData%real3
    CalFplusPrimEnd = tiData%real4

    primStar = ti_efold_primitive(x,alpha,mu)
    epsOneStar = ti_epsilon_one(x,alpha,mu)
    potStar = ti_norm_potential(x,alpha)


    find_ti_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)


  end function find_ti_x_star


  function ti_lnrhoend(alpha,mu,Pstar) 
    implicit none
    real(kp) :: ti_lnrhoend
    real(kp), intent(in) :: alpha,mu,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    xEnd=ti_x_endinf(alpha,mu) 
    potEnd  = ti_norm_potential(xEnd,alpha)
    epsOneEnd = ti_epsilon_one(xEnd,alpha,mu)

       
    x = ti_x_star(alpha,mu,wrad,junk,Pstar)    
    potStar = ti_norm_potential(x,alpha)
    epsOneStar = ti_epsilon_one(x,alpha,mu)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'ti_lnrhoend: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ti_lnrhoend = lnRhoEnd

  end function ti_lnrhoend

  
  
end module tireheat
