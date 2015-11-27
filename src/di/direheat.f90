!Reheating functions for Dual Inflation.
!
!These are non-standard as M^4, the potential normalization, is
!completely determined by the two model parameters f and Lambda. As a
!result, for a given f, Lambda is completely determined by the CMB
!amplitude knowing the reheating parameter, or phistar. However, at
!given reheating parameter, phistar can only be solved knowing
!Lambda. So the reheating equations and potential normalisation are
!two coupled non-linear algebraic equations to be solved
!simultaneously given the reheating parameter and f. The solutions of
!which are phistar and Lambda. Fortunately, the equation for the CMB
!to M^4 normalisation can analytically be solved for Lambda and reads:
!
! Lambda^6 = 24 pi^4 Pstar * [pepsOne(k2star)/(v(k2star)*f^2]
!
!Here pepsOne = epsilon1/Lambda^2 (which is not a function of Lambda),
!is the parametric epsilon1 function which depends on f and k2
!only. However, on top of that, the k2end value at which inflation
!ends depends on Lambda thereby rendering the numerical determination
!of the solution challenging. The reheating equations are accordingly
!modified by using the above expression, taking care of k2end
!consistently such that "k2star" is solved given the reheating and
!f. Lambda is then determined uniquely using the above equation.

module direheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use dicommon, only : di_parametric_epsilon_one, di_norm_parametric_potential
  use dicommon, only : di_parametric_efold_primitive, di_k2_potmin
  use disr, only : di_k2_epsoneunity
  
  implicit none

  real(kp), parameter :: sixth = 1._kp/6._kp
  real(kp), parameter :: third = 1._kp/3._kp


  private
  
  public di_lambda_star, di_lnrhoreh_max
  public di_k2_star, di_k2_rrad, di_k2_rreh
!  public di_x_star, di_x_rrad, di_x_rreh

contains

!return lambda from k2, f and Pstar(CMB normalised only for k2=k2star)
  function di_lambda_star(k2,f,Pstar)
    implicit none
    real(kp) :: di_lambda_star
    real(kp), intent(in) :: k2,f
    real(kp) :: peps, ppot, Pstar

    peps = di_parametric_epsilon_one(k2,f)
    ppot = di_norm_parametric_potential(k2,f)

    di_lambda_star = (24._kp*pi**4*Pstar*peps/ppot/f/f)**sixth

  end function di_lambda_star


!returns k2star given potential parameters, scalar power, wreh and
!lnrhoreh
  function di_k2_star(f,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: di_k2_star
    real(kp), intent(in) :: f,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF, pi4Pstar24
    real(kp) :: epsOneEnd,ppotEnd
    real(kp) :: k2min, k2max, k2potmin, k2star, k2end, lambda
    type(transfert) :: diData
  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = 1._kp
!trick to use the std calfconst function
    ppotEnd = 1._kp

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,ppotEnd)

    pi4Pstar24 = 24._kp*pi**4*Pstar
         
!calF(Vend=1) + DI corrections coming from replacing Lambda by the
!CMB normalization
    diData%real1 = f
    diData%real2 = calF - log(pi4Pstar24/f/f)*(1._kp+3._kp*w)/(18._kp*(1._kp+w))
    diData%real3 = pi4Pstar24
    diData%real4 = w

    k2potmin = di_k2_potmin(f)

    k2min= epsilon(1._kp)
    k2max = (1._kp-tolkp)*k2potmin


    k2star = zbrent(find_di_k2_star,k2min,k2max,tolzbrent,diData)    

    di_k2_star = k2star

    if (present(bfoldstar)) then
       lambda = di_lambda_star(k2star,f,Pstar)
       k2end = di_k2_epsoneunity(f,lambda)

       bfoldstar = -( di_parametric_efold_primitive(k2star,f) &
            - di_parametric_efold_primitive(k2end,f) )*lambda**2
    end if

  end function di_k2_star


  function find_di_k2_star(k2,diData)   
    implicit none
    real(kp) :: find_di_k2_star
    real(kp), intent(in) :: k2
    type(transfert), optional, intent(inout) :: diData
    real(kp) :: pi4Pstar24
    real(kp) :: f, pprimStar,ppotStar,pepsOneStar, lambdaStar
    real(kp) :: k2end, w
    real(kp) :: effCalFend, pprimEnd, ppotEnd, effprimStar

    f=diData%real1
    effCalFend = diData%real2
    pi4Pstar24 = diData%real3  
    w = diData%real4

    pepsOneStar = di_parametric_epsilon_one(k2,f)
    if (pepsOneStar.eq.0._kp) pepsOneStar = epsilon(1._kp)

    ppotStar = di_norm_parametric_potential(k2,f)
!numacc towards minimum
    if (ppotStar.eq.0._kp) ppotStar = epsilon(1._kp)

    lambdaStar = (pi4Pstar24*pepsOneStar/ppotStar/f/f)**sixth    
    
    k2end = di_k2_epsoneunity(f,lambdaStar)

    pprimStar = di_parametric_efold_primitive(k2,f)

    pprimEnd = di_parametric_efold_primitive(k2end,f)

    ppotEnd = di_norm_parametric_potential(k2end,f)
    
    effprimStar = (pprimStar - pprimEnd)*lambdaStar**2 &
         + log(pepsOneStar/ppotStar)*(1._kp+3._kp*w)/(18._kp*(1._kp+w)) &
         + log(ppotEnd)/(3._kp+3._kp*w)

    find_di_k2_star = find_reheat(effprimStar,effcalFend,w,pepsOneStar,ppotStar)
  
  end function find_di_k2_star




!returns k2 given potential parameters, scalar power, and lnRrad
  function di_k2_rrad(f,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: di_k2_rrad
    real(kp), intent(in) :: f,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF, pi4Pstar24
    real(kp) :: epsOneEnd,ppotEnd
    real(kp) :: k2min, k2max, k2potmin, k2star, k2end, lambda
    type(transfert) :: diData

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = 1._kp
!trick to use the std calfconst function
    ppotEnd = 1._kp
    pi4Pstar24 = 24._kp*pi**4*Pstar
    
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,ppotEnd)
         
!calF standard + DI corrections coming from replacing Lambda by the
!CMB normalization
    diData%real1 = f
    diData%real2 = calF - log(pi4Pstar24/f/f)/12._kp
    diData%real3 = pi4Pstar24
  
    k2potmin = di_k2_potmin(f)

    k2min= epsilon(1._kp)
    k2max = (1._kp-tolkp)*k2potmin


    k2star = zbrent(find_di_k2_rrad,k2min,k2max,tolzbrent,diData)

    di_k2_rrad = k2star 


    if (present(bfoldstar)) then
       lambda = di_lambda_star(k2star,f,Pstar)
       k2end = di_k2_epsoneunity(f,lambda)

       bfoldstar = -( di_parametric_efold_primitive(k2star,f) &
            - di_parametric_efold_primitive(k2end,f) )*lambda**2
    end if

  end function di_k2_rrad

  function find_di_k2_rrad(k2,diData)   
    implicit none
    real(kp) :: find_di_k2_rrad
    real(kp), intent(in) :: k2
    type(transfert), optional, intent(inout) :: diData
    real(kp) :: pi4Pstar24, k2end
    real(kp) :: f, pprimStar,ppotStar,pepsOneStar, lambdaStar
    real(kp) :: effCalFend, pprimEnd, ppotEnd, effprimStar

    f=diData%real1
    effCalFend = diData%real2
    pi4Pstar24 = diData%real3

    pepsOneStar = di_parametric_epsilon_one(k2,f)

    ppotStar = di_norm_parametric_potential(k2,f)
    if (ppotStar.eq.0._kp) ppotStar = epsilon(1._kp)

    lambdaStar = (pi4Pstar24*pepsOneStar/ppotStar/f/f)**sixth

    k2end = di_k2_epsoneunity(f,lambdaStar)

    pprimStar = di_parametric_efold_primitive(k2,f)

    pprimEnd = di_parametric_efold_primitive(k2end,f)

    ppotEnd = di_norm_parametric_potential(k2end,f)
    
    effprimStar = (pprimStar - pprimEnd)*lambdaStar**2 &
         + log(pepsOneStar/ppotStar)/12._kp  + 0.25_kp * log(ppotEnd)   

    find_di_k2_rrad = find_reheat_rrad(effprimStar,effCalFend,pepsOneStar,ppotStar)
  
  end function find_di_k2_rrad
  


!returns k2 given potential parameters, scalar power, and lnRreh
  function di_k2_rreh(f,lnRreh,Pstar,bfoldstar)
    implicit none
    real(kp) :: di_k2_rreh
    real(kp), intent(in) :: f,lnRreh,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF, pi4Pstar24
    real(kp) :: epsOneEnd,ppotEnd
    real(kp) :: k2min, k2max, k2potmin, k2star, k2end, lambda
    type(transfert) :: diData

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = 1._kp
!trick to use the std calfconst function
    ppotEnd = 1._kp
    pi4Pstar24 = 24._kp*pi**4*Pstar
    
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,ppotEnd)
         
    diData%real1 = f
    diData%real2 = calF
    diData%real3 = pi4Pstar24
  
    k2potmin = di_k2_potmin(f)

    k2min= epsilon(1._kp)
    k2max = (1._kp-tolkp)*k2potmin

    k2star = zbrent(find_di_k2_rreh,k2min,k2max,tolzbrent,diData)

    di_k2_rreh = k2star

    if (present(bfoldstar)) then
       lambda = di_lambda_star(k2star,f,Pstar)
       k2end = di_k2_epsoneunity(f,lambda)
       bfoldstar = -( di_parametric_efold_primitive(k2star,f) &
            - di_parametric_efold_primitive(k2end,f) )*lambda**2
    end if

  end function di_k2_rreh

  function find_di_k2_rreh(k2,diData)   
    implicit none
    real(kp) :: find_di_k2_rreh
    real(kp), intent(in) :: k2
    type(transfert), optional, intent(inout) :: diData
    real(kp) :: pi4Pstar24, k2end
    real(kp) :: f, pprimStar,ppotStar,pepsOneStar, lambdaStar
    real(kp) :: effCalFend, pprimEnd, ppotEnd, effprimStar

    f=diData%real1
    effCalFend = diData%real2
    pi4Pstar24 = diData%real3

    pepsOneStar = di_parametric_epsilon_one(k2,f)

    ppotStar = di_norm_parametric_potential(k2,f)
    if (ppotStar.eq.0._kp) ppotStar = epsilon(1._kp)

    lambdaStar = (pi4Pstar24*pepsOneStar/ppotStar/f/f)**sixth
    k2end = di_k2_epsoneunity(f,lambdaStar)

    pprimStar = di_parametric_efold_primitive(k2,f)

    pprimEnd = di_parametric_efold_primitive(k2end,f)

    ppotEnd = di_norm_parametric_potential(k2end,f)
    
    effprimStar =  (pprimStar - pprimEnd)*lambdaStar**2 + 0.5_kp*log(ppotEnd)

    find_di_k2_rreh = find_reheat_rreh(effprimStar,effCalFend,ppotStar)  

  end function find_di_k2_rreh
  

  function di_lnrhoreh_max(f,Pstar)
    implicit none
    real(kp) :: di_lnrhoreh_max
    real(kp), intent(in) :: f, Pstar

    real(kp), parameter :: lnRreh=0._kp

    real(kp) :: pi4Pstar24
    real(kp) :: k2star, k2end, lambdaStar
    real(kp) :: ppotStar, pepsOneStar, ppotEnd
    real(kp) :: epsOneStar, lnRhoEnd, epsOneEnd

    pi4Pstar24 = 24._kp*pi**4*Pstar

    k2star = di_k2_rreh(f,lnRreh,Pstar)

    pepsOneStar = di_parametric_epsilon_one(k2star,f)
    ppotStar = di_norm_parametric_potential(k2star,f)
    !numacc towards minimum
    if (ppotStar.eq.0._kp) ppotStar = epsilon(1._kp)

    lambdaStar = (pi4Pstar24*pepsOneStar/ppotStar/f/f)**sixth

    k2end = di_k2_epsoneunity(f,lambdaStar)
    ppotEnd = di_norm_parametric_potential(k2end,f)

    epsOneStar = pepsOneStar/lambdaStar/lambdaStar

    epsOneEnd = di_parametric_epsilon_one(k2end,f)/lambdaStar/lambdaStar

    if (.not.slowroll_validity(epsOneStar)) stop 'di_lnrhoreh_max: slow-roll violated!'

    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,ppotEnd/ppotStar)

    di_lnrhoreh_max = lnRhoEnd

  end function di_lnrhoreh_max


end module direheat
