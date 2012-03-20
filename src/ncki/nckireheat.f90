!non canonical Kahler inflation reheating functions in the slow-roll approximations

module nckireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use nckisr, only : ncki_epsilon_one, ncki_epsilon_two, ncki_epsilon_three
  use nckisr, only : ncki_norm_potential
  use nckisr, only : ncki_x_endinf, ncki_efold_primitive
  implicit none

  private

  public ncki_x_star, ncki_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ncki_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ncki_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: nckiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = ncki_x_endinf(alpha,beta)
    epsOneEnd = ncki_epsilon_one(xEnd,alpha,beta)
    potEnd = ncki_norm_potential(xEnd,alpha,beta)
    primEnd = ncki_efold_primitive(xEnd,alpha,beta)

!    print*,'ncki_x_star:    xEnd=',xEnd,'  epsOneEnd=',epsOneEnd,'  potEnd=',potEnd,'  primEnd=',primEnd

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    nckiData%real1 = alpha 
    nckiData%real2 = beta
    nckiData%real3 = xEnd
    nckiData%real4 = w
    nckiData%real5 = calF + primEnd

    mini=xEnd
    if (beta.lt.0._kp) then
      maxi = sqrt(alpha/(2._kp*abs(beta)))*(1._kp-epsilon(1._kp)) !position of the maximum of the potential if beta<0
    else
      maxi = 1000._kp* sqrt(alpha/(2._kp*abs(beta))) !several times the position of the inflexion point if beta>0
    endif

!    print*,'ncki_x_star:    mini=',mini,'  maxi=',maxi,'  f(mini)=',find_ncki_x_star(mini,nckiData), &
!            '  f(maxi)=',find_ncki_x_star(maxi,nckiData),'  prim(mini)=', & 
!            ncki_efold_primitive(mini,alpha,beta),'  prim(maxi)=',ncki_efold_primitive(maxi,alpha,beta)



        x = zbrent(find_ncki_x_star,mini,maxi,tolzbrent,nckiData)
        ncki_x_star = x








    if (present(bfoldstar)) then
       bfoldstar = - (ncki_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ncki_x_star

  function find_ncki_x_star(x,nckiData)   
    implicit none
    real(kp) :: find_ncki_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: nckiData

    real(kp) :: primStar,alpha,beta,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=nckiData%real1
    beta=nckiData%real2
    xEnd=nckiData%real3
    w = nckiData%real4
    CalFplusprimEnd = nckiData%real5

    primStar = ncki_efold_primitive(x,alpha,beta)
    epsOneStar = ncki_epsilon_one(x,alpha,beta)
    potStar = ncki_norm_potential(x,alpha,beta)

    find_ncki_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ncki_x_star



  function ncki_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ncki_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ncki_x_endinf(alpha,beta)
    potEnd  = ncki_norm_potential(xEnd,alpha,beta)
    epsOneEnd = ncki_epsilon_one(xEnd,alpha,beta)

!   Trick to return x such that rho_reh=rho_end

    x = ncki_x_star(alpha,beta,wrad,junk,Pstar)    
    potStar = ncki_norm_potential(x,alpha,beta)
    epsOneStar = ncki_epsilon_one(x,alpha,beta)

    !print*,'ncki_lnrhoend:   xstar=',x,'  potStar=',potStar,'  epsOneStar=',epsOneStar


    
    if (.not.slowroll_validity(epsOneStar)) stop 'ncki_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ncki_lnrhoend = lnRhoEnd

  end function ncki_lnrhoend

  
end module nckireheat
