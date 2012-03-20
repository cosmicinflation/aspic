!Logamediate inflation 1 reheating functions in the slow-roll approximations

module lmi3reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lmi3sr, only : lmi3_epsilon_one, lmi3_epsilon_two, lmi3_epsilon_three
  use lmi3sr, only : lmi3_norm_potential
  use lmi3sr, only : lmi3_efold_primitive
  use lmi3sr, only : lmi3_alpha, lmi3_beta
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public lmi3_x_star, lmi3_lnrhoend,lmi3_M

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lmi3_x_star(gamma_lmi,M,xEnd,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lmi3_x_star
    real(kp), intent(in) :: gamma_lmi,M,xEnd,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: lmi3Data

    real(kp) ::alpha,beta
    alpha=lmi3_alpha(gamma_lmi,M)
    beta=lmi3_beta(gamma_lmi,M)

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = lmi3_epsilon_one(xEnd,gamma_lmi,M)
    potEnd = lmi3_norm_potential(xEnd,gamma_lmi,M)

    primEnd = lmi3_efold_primitive(xEnd,gamma_lmi,M)
   

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    lmi3Data%real1 = gamma_lmi
    lmi3Data%real2 = M
    lmi3Data%real3 = w
    lmi3Data%real4 = calF + primEnd

    maxi = xEnd*(1._kp-epsilon(1._kp))
    mini = 10._kp*((alpha/(beta*gamma_lmi))**(1._kp/gamma_lmi))*(1._kp+epsilon(1._kp))

!    print*,'lmi3_x_star:   ','mini=',mini,'maxi=',maxi,'xEnd=',xEnd, &
!          'f(mini)=',find_lmi3_x_star(mini,lmi3Data), &
!                 'f(maxi)=',find_lmi3_x_star(maxi,lmi3Data)


!    print*,'primStar(mini)=',lmi3_efold_primitive(mini,gamma_lmi,M), &
!           'epsOneStar(mini)=',lmi3_epsilon_one(mini,gamma_lmi,M), &
!           'potStar(mini)=',lmi3_norm_potential(mini,gamma_lmi,M)

    x = zbrent(find_lmi3_x_star,mini,maxi,tolzbrent,lmi3Data)
    lmi3_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi3_efold_primitive(x,gamma_lmi,M) - primEnd)
    endif

  end function lmi3_x_star

  function find_lmi3_x_star(x,lmi3Data)   
    implicit none
    real(kp) :: find_lmi3_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi3Data

    real(kp) :: primStar,gamma_lmi,M,w,CalFplusprimEnd,potStar,epsOneStar

    gamma_lmi=lmi3Data%real1
    M=lmi3Data%real2
    w = lmi3Data%real3
    CalFplusprimEnd = lmi3Data%real4

    primStar = lmi3_efold_primitive(x,gamma_lmi,M)
    epsOneStar = lmi3_epsilon_one(x,gamma_lmi,M)
    potStar = lmi3_norm_potential(x,gamma_lmi,M)

    find_lmi3_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lmi3_x_star



  function lmi3_lnrhoend(gamma_lmi,M,xEnd,Pstar) 
    implicit none
    real(kp) :: lmi3_lnrhoend
    real(kp), intent(in) :: gamma_lmi,M,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd

    potEnd  = lmi3_norm_potential(xEnd,gamma_lmi,M)

    epsOneEnd = lmi3_epsilon_one(xEnd,gamma_lmi,M)


!   Trick to return x such that rho_reh=rho_end

    x = lmi3_x_star(gamma_lmi,M,xEnd,wrad,junk,Pstar)  


    potStar = lmi3_norm_potential(x,gamma_lmi,M)
    epsOneStar = lmi3_epsilon_one(x,gamma_lmi,M)

    print*,'xstar=',x,'potStar=',potStar,'epsOneStar=',epsOneStar

    print*,'x0=',(4._kp*(1._kp-gamma_lmi/2._kp)/sqrt(2._kp*sqrt(3._kp)) &
                *M)**(1._kp/gamma_lmi)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'lmi3_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lmi3_lnrhoend = lnRhoEnd

  end function lmi3_lnrhoend



!Finds the value of M compatible with the COBE normalization, given gamma,xEnd,lnRhoReh,w and Pstar
    function lmi3_M(gamma_lmi,xEnd,lnRhoReh,w,Pstar)
    implicit none
    real(kp), intent(in) :: gamma_lmi,xEnd,lnRhoReh,w,Pstar
    real(kp) :: lmi3_M,alpha,junk
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi3Data
    integer ::i

    junk=0._kp
    alpha=lmi3_alpha(gamma_lmi,junk)

    maxi = (xEnd**gamma_lmi*sqrt(2._kp*sqrt(3._kp))/alpha)
    mini = maxi/100.


    lmi3Data%real1 = gamma_lmi
    lmi3Data%real2 = xEnd
    lmi3Data%real3 = lnRhoReh
    lmi3Data%real4 = w
    lmi3Data%real5 = Pstar

 !   print*,'lmi3_M:        mini=',mini,'maxi=',maxi
 !   print*,'f(mini)=',find_lmi3M(1._kp,lmi3Data)
 !   print*,'f(maxi)=',find_lmi3M(maxi/10.,lmi3Data) 


      
    lmi3_M = zbrent(find_lmi3M,mini,maxi,tolFind,lmi3Data)


    end function lmi3_M

  function find_lmi3M(x,lmi3Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi3Data
    real(kp) :: find_lmi3M
    real(kp) :: gamma_lmi,xEnd,lnRhoReh,w,Pstar,xstar,alpha,beta

    gamma_lmi = lmi3Data%real1
    xEnd = lmi3Data%real2
    lnRhoReh = lmi3Data%real3
    w = lmi3Data%real4
    Pstar = lmi3Data%real5
    xstar = lmi3_x_star(gamma_lmi,x,xEnd,w,lnRhoReh,Pstar)
    alpha=lmi3_alpha(gamma_lmi,x)
    beta=lmi3_beta(gamma_lmi,x)
    
    find_lmi3M = 720._kp*acos(-1._kp)**2*(alpha-beta*gamma_lmi*xstar**gamma_lmi)**2* &
                exp(beta*xstar**gamma_lmi)*xstar**(-alpha-2._kp)*QrmsOverT**2-x**4
   
  end function find_lmi3M

  
end module lmi3reheat
