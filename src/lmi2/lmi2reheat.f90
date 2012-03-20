!Logamediate inflation 1 reheating functions in the slow-roll approximations

module lmi2reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lmi2sr, only : lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three
  use lmi2sr, only : lmi2_norm_potential
  use lmi2sr, only : lmi2_x_endinf, lmi2_efold_primitive
  use lmi2sr, only : lmi2_alpha, lmi2_beta
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public lmi2_x_star, lmi2_lnrhoend,lmi2_M

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lmi2_x_star(gamma_lmi,M,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lmi2_x_star
    real(kp), intent(in) :: gamma_lmi,M,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: lmi2Data

    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lmi2_x_endinf(gamma_lmi,M)

    epsOneEnd = lmi2_epsilon_one(xEnd,gamma_lmi,M)
    potEnd = lmi2_norm_potential(xEnd,gamma_lmi,M)

    primEnd = lmi2_efold_primitive(xEnd,gamma_lmi,M)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    lmi2Data%real1 = gamma_lmi
    lmi2Data%real2 = M
    lmi2Data%real3 = w
    lmi2Data%real4 = calF + primEnd

    maxi = lmi2_x_endinf(gamma_lmi,M)*(1._kp-epsilon(1._kp))
    mini = (alpha/(beta*gamma_lmi))**(1._kp/gamma_lmi)*(1._kp+epsilon(1._kp))

    !print*,'lmi2_x_star','mini=',mini,'maxi=',maxi,'find...(mini)=',find_lmi2_x_star(mini,lmi2Data) &
    !  ,'find...(maxi)=',find_lmi2_x_star(maxi,lmi2Data)

    x = zbrent(find_lmi2_x_star,mini,maxi,tolzbrent,lmi2Data)
    lmi2_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi2_efold_primitive(x,gamma_lmi,M) - primEnd)
    endif

  end function lmi2_x_star

  function find_lmi2_x_star(x,lmi2Data)   
    implicit none
    real(kp) :: find_lmi2_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi2Data

    real(kp) :: primStar,gamma_lmi,M,w,CalFplusprimEnd,potStar,epsOneStar

    gamma_lmi=lmi2Data%real1
    M=lmi2Data%real2
    w = lmi2Data%real3
    CalFplusprimEnd = lmi2Data%real4

    primStar = lmi2_efold_primitive(x,gamma_lmi,M)
    epsOneStar = lmi2_epsilon_one(x,gamma_lmi,M)
    potStar = lmi2_norm_potential(x,gamma_lmi,M)

    find_lmi2_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lmi2_x_star



  function lmi2_lnrhoend(gamma_lmi,M,Pstar) 
    implicit none
    real(kp) :: lmi2_lnrhoend
    real(kp), intent(in) :: gamma_lmi,M,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = lmi2_x_endinf(gamma_lmi,M)

    potEnd  = lmi2_norm_potential(xEnd,gamma_lmi,M)

    epsOneEnd = lmi2_epsilon_one(xEnd,gamma_lmi,M)

    print*,'xEnd=',xEnd,'potEnd=',potEnd,'epsOneEnd=',epsOneEnd


!   Trick to return x such that rho_reh=rho_end

    x = lmi2_x_star(gamma_lmi,M,wrad,junk,Pstar)  


    potStar = lmi2_norm_potential(x,gamma_lmi,M)
    epsOneStar = lmi2_epsilon_one(x,gamma_lmi,M)

    print*,'xstar=',x,'potStar=',potStar,'epsOneStar=',epsOneStar

    print*,'x0=',(4._kp*(1._kp-gamma_lmi/2._kp)/sqrt(2._kp*sqrt(3._kp)) &
                *M)**(1._kp/gamma_lmi)

   
    !print*,'x**gamma=',x**gamma_lmi,'beta*gamma*x**gamma',sqrt(2._kp*sqrt(3._kp))/(M*gamma_lmi) &
    !          *gamma_lmi*x**gamma_lmi,'alpha=',4._kp*(1._kp-gamma_lmi/2._kp), &
    !          'alpha-beta*gamma*x**gamma=',4._kp*(1._kp-gamma_lmi/2._kp)- &
    !      sqrt(2._kp*sqrt(3._kp))/(M*gamma_lmi)*gamma_lmi*x**gamma_lmi
    !print*,'num=',(4._kp*(1._kp-gamma_lmi/2._kp)- &
    !      sqrt(2._kp*sqrt(3._kp))/(M*gamma_lmi)*gamma_lmi*x**gamma_lmi)**2
    !print*,'denum=',2._kp*x**2
    !print*,'eps1star=',(4._kp*(1._kp-gamma_lmi/2._kp)- &
    !      sqrt(2._kp*sqrt(3._kp))/(M*gamma_lmi)*gamma_lmi*x**gamma_lmi)**2/ &
    !        (2._kp*x**2) 
    
    if (.not.slowroll_validity(epsOneStar)) stop 'lmi2_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lmi2_lnrhoend = lnRhoEnd

  end function lmi2_lnrhoend



!Finds the value of M compatible with the COBE normalization, given gamma,lnRhoReh,w and Pstar
    function lmi2_M(gamma_lmi,lnRhoReh,w,Pstar)
    implicit none
    real(kp), intent(in) :: gamma_lmi,lnRhoReh,w,Pstar
    real(kp) :: lmi2_M,alpha,junk
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi2Data
    integer ::i

    junk=0._kp
    alpha=lmi2_alpha(gamma_lmi,junk)

    maxi = 2._kp/(alpha**2)*(alpha/((1._kp-gamma_lmi)*sqrt(2._kp*sqrt(3._kp)))) &
            **(2._kp/gamma_lmi)*((1._kp-gamma_lmi)/gamma_lmi)**2
    mini = maxi/1000._kp

    lmi2Data%real1 = gamma_lmi
    lmi2Data%real2 = lnRhoReh
    lmi2Data%real3 = w
    lmi2Data%real4 = Pstar

    print*,'lmi2_M:     mini=',mini,'maxi=',maxi,'f(mini)=',find_lmi2M(mini,lmi2Data), &
                 'f(maxi)=',find_lmi2M(maxi,lmi2Data)

   do i=1,100
    print*,'M=',mini*(maxi/mini)**(real(i,kp)/real(100,kp)), &
            'f(M)=',find_lmi2M(mini*(maxi/mini)**(real(i,kp)/real(100,kp)),lmi2Data)
   end do
    
    lmi2_M = zbrent(find_lmi2M,mini,maxi,tolFind,lmi2Data)


    end function lmi2_M

  function find_lmi2M(x,lmi2Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2M
    real(kp) :: gamma_lmi,lnRhoReh,w,Pstar,xstar,alpha,beta

    gamma_lmi = lmi2Data%real1
    lnRhoReh = lmi2Data%real2
    w = lmi2Data%real3
    Pstar = lmi2Data%real4
    xstar = lmi2_x_star(gamma_lmi,x,w,lnRhoReh,Pstar)
    alpha=lmi2_alpha(gamma_lmi,x)
    beta=lmi2_beta(gamma_lmi,x)
    
    find_lmi2M = 720._kp*acos(-1._kp)**2*(alpha-beta*gamma_lmi*xstar**gamma_lmi)**2* &
                exp(beta*xstar**gamma_lmi)*xstar**(-alpha-2._kp)*QrmsOverT**2-x**4
   
  end function find_lmi2M

  
end module lmi2reheat
