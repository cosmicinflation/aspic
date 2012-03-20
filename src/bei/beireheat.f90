!beta inflation reheating functions in the slow-roll approximabeions

module beireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use beisr, only : bei_epsilon_one, bei_epsilon_two, bei_epsilon_three
  use beisr, only : bei_norm_potential,bei_efold_primitive,bei_x_endinf
  implicit none

  private

  public bei_x_star, bei_lnrhoend,find_bei_x_star

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function bei_x_star(lambda,beta,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: bei_x_star
    real(kp), intent(in) :: lambda,beta,w,lnRhoReh,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd

    type(transfert) :: beiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=bei_x_endinf(lambda,beta)
    
    epsOneEnd = bei_epsilon_one(xEnd,lambda,beta)
    potEnd = bei_norm_potential(xEnd,lambda,beta)
    primEnd = bei_efold_primitive(xEnd,lambda,beta)

!    print*,'bei_x_star:   xEnd=',xEnd,'  epsOneEnd=',epsOneEnd,'  potEnd=',potEnd,'  primEnd=',primEnd
   

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    beiData%real1 = lambda
    beiData%real2 = beta
    beiData%real3 = w
    beiData%real4 = calF + primEnd

    maxi = xEnd*(1._kp-epsilon(1._kp))
    mini=maxi-1000._kp


!    print*,'bei_x_star:   mini=',mini,'  maxi=',maxi,'  f(mini)=',find_bei_x_star(mini,beiData), &
!               '  f(maxi)=',find_bei_x_star(maxi,beiData)


    x = zbrent(find_bei_x_star,mini,maxi,tolFind,beiData)
    bei_x_star = x

    if (present(bfold)) then
       bfold = -(bei_efold_primitive(x,lambda,beta) - primEnd)
    endif


  end function bei_x_star

  function find_bei_x_star(x,beiData)   
    implicit none
    real(kp) :: find_bei_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: beiData

    real(kp) :: primStar,lambda,beta,w,CalFplusPrimEnd,potStar,epsOneStar

    lambda=beiData%real1
    beta=beiData%real2
    w = beiData%real3
    CalFplusPrimEnd = beiData%real4

    primStar = bei_efold_primitive(x,lambda,beta)
    epsOneStar = bei_epsilon_one(x,lambda,beta)
    potStar = bei_norm_potential(x,lambda,beta)


    find_bei_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)


  end function find_bei_x_star


  function bei_lnrhoend(lambda,beta,Pstar) 
    implicit none
    real(kp) :: bei_lnrhoend
    real(kp), intent(in) :: lambda,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    xEnd=bei_x_endinf(lambda,beta) 
    potEnd  = bei_norm_potential(xEnd,lambda,beta)
    epsOneEnd = bei_epsilon_one(xEnd,lambda,beta)

!    print*,'bei_lnrhoend:   xEnd=',xEnd,'  potEnd=',potEnd,'  epsOneEnd=',epsOneEnd

       
    x = bei_x_star(lambda,beta,wrad,junk,Pstar)    
    potStar = bei_norm_potential(x,lambda,beta)
    epsOneStar = bei_epsilon_one(x,lambda,beta)

!    print*,'bei_lnrhoend:   xstar=',x,'  potStar=',potStar,'  epsOneStar=',epsOneStar
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'bei_lnrhoend: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    bei_lnrhoend = lnRhoEnd

  end function bei_lnrhoend

  
  
end module beireheat
