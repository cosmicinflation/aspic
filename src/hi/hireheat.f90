!Reheating functions for Higgs Inflation.
!
!These are non-standard as M^4, the potential normalization, is
!completely determined by the only model parameters xi and xstar (or
!hbarstar). So the reheating equations and potential normalisation are
!two coupled non-linear algebraic equations to be solved
!simultaneously given the reheating parameter and xi. The solutions of
!which are hbarstar (->xstar) and xi. The equation for the CMB to M^4
!normalisation cannot analytically be solved for xi, but assuming the
!Higgs vev to be zero, there is an analytical solution. This solution
!is used to start a recursion allowing to determinte the exact value
!of xi given hstar.
!
!This value of xi is everytime recomputed when zbrenting the reheating
!equation. Let us also notice that the reheating equation are
!non-standard due to the fact that the reheating parameters to be
!plugged in are in the Jordan Frame (see Encyclopedia
!Inflationaris). As such, expressed in terms of rhoend (Einstein
!Frame), there is an extra conformal factor in (1+hbarend^2)^2.

module hireheat
  use infprec, only :  kp, tolkp, toldp, pidp, transfert
  use cosmopar, only : HiggsCoupling
  use inftools, only : zbrent
  use srreheat, only : display
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : pi, ln_rho_enhinf, ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh

  use hicommon, only : vev2
  use hicommon, only : hi_parametric_epsilon_one, hi_norm_parametric_potential
  use hicommon, only : hi_parametric_efold_primitive, hi_parametric_hbar_endinf
  use hisr, only : hi_x
  
  implicit none

  real(kp), parameter :: lambda = HiggsCoupling
  real(kp), parameter :: pstarconst = 56._kp*pidp*pidp/lambda

  real(kp), parameter :: onethird = 1._kp/3._kp
  real(kp), parameter :: twothird = 2._kp/3._kp

  private
  
  public hi_lambda_star, hi_lnrhoreh_max
  public hi_k2_star, hi_k2_rrad, hi_k2_rreh
  public hi_x_star, hi_x_rrad, hi_x_rreh

contains


!xistar is an output
  function hi_x_star(w,lnRhoReh,Pstar,bfoldstar,xistar)
    implicit none
    real(kp) :: hi_x_star
    real(kp), intent(in) :: w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfoldstar
    real(kp), intent(out),optional :: xistar
    
    real(kp) :: hbarstar

    hbarstar = hi_hbar_star(xistar,w,lnRhoReh,Pstar,bfoldstar)

    hi_x_star = hi_x(hbarstar,xistar)
    
  end function hi_x_star



  function hi_x_rrad(lnRrad,Pstar,bfoldstar,xistar)    
    implicit none
    real(kp) :: hi_x_rrad
    real(kp), intent(in) :: lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar
    real(kp), intent(out), optional :: xistar
    
    real(kp) :: hbarstar

    hbarstar = hi_hbar_rrad(xistar,lnRrad,Pstar,bfoldstar)

    hi_x_rrad = hi_x(hbarstar,xistar)

  end function hi_x_rrad

    
!calling arguments are non-standard and require Pstar (due to Lambda)
!as opposed to the usual *_x_rreh functions
  function hi_x_rreh(lnRreh,Pstar,bfoldstar,xistar)
    implicit none
    real(kp) :: hi_x_rreh
    real(kp), intent(in) :: lnRreh,Pstar
    real(kp), intent(out), optional :: bfoldstar
    real(kp), intent(out), optional :: xistar
    
    real(kp) :: hbarstar

    hbarstar = hi_hbar_rreh(xistar,lnRreh,Pstar,bfoldstar)

    hi_x_rreh = hi_x(hbarstar,xistar)
    
  end function hi_x_rreh



!return xistar, matching the CMB normalisation, from any input values
!of hbar (not necessarily solution of the reheating equation that we
!refer to as hbarstar). We guess first with the solution obtained for
!vev=0 and then recurse to get the right CMB normalisation at
!non-vanishing vev.
  recursive function hi_xi_star(hbar,Pstar) result (xistar)
    implicit none
    real(kp) :: xistar
    real(kp), intent(in) :: hbar,Pstar
    real(kp) :: hbar2,C

    real(kp), parameter :: tolfind = epsilon(1._kp)
    logical, parameter :: display = .true.
    
!the only posituve solution defined for all hbar and obtained by assuming
!HiggsVev = 0

    hbar2 = hbar*hbar
    C = pstarconst * Pstar
    
    xistar = real( ((2._kp + cmplx(0._kp,2._kp)*sqrt(3._kp))*hbar2**4._kp &
         + (2._kp**onethird*(1._kp - cmplx(0._kp,1._kp)*sqrt(3._kp)) &
         * (-(C**2._kp*hbar2**3._kp*(1._kp + hbar2)**5._kp) &
         + sqrt(C**3._kp*hbar2**6._kp*(1._kp + hbar2)**6._kp &
         * (-4._kp*hbar2**6._kp + C*(1._kp + hbar2)**4._kp)))**twothird) &
         / (C*(1._kp + hbar2)**2._kp))/(4._kp*2._kp**twothird &
         * (- (C**2._kp*hbar2**3._kp*(1._kp + hbar2)**5._kp) &
         + sqrt(C**3._kp*hbar2**6._kp*(1._kp + hbar2)**6._kp &
         * (-4._kp*hbar2**6._kp + C*(1._kp + hbar2)**4._kp)))**onethird) &
         ,kp)

!W(hbar)/eps1(hbar)/xi^2
    potstar = hi_norm_parametric_potential(hbar,xistar)
    epsstar = hi_parametric_epsilon_one(hbar,xistar)

    Cexact = potstar/(xistar*xistar*epsstar)
    
    if ((2._kp*abs(Cexact-C)/(Cexact+C)).gt.tolfind) then
       if (display) then
          write(*,*)'hi_xi_star:'
          write(*,*)'C= Cexact= ',C,Cexact
          write(*,*)'recursing on xistar = ',xistar
          write(*,*)
       endif
       xistar = hi_xi_star(hbar,Pstar*Cexact/C)
    endif
           
  end function hi_xi_star


!returns hbarstar, and xistar, given the scalar power, wreh and
!lnrhoreh
  function hi_hbar_star(w,lnRhoReh,Pstar,bfoldstar,xistar)
    implicit none
    real(kp) :: hi_hbar_star
    real(kp), intent(in) :: lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar, xistar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF, xi
    real(kp) :: epsOneEnd,potEnd

    real(kp) :: hbarmin,hbarmax

    type(transfert) :: hiData
  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = 1._kp
!trick to use the std calfconst function not including the term in potend
    potEnd = 1._kp

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)
         
!calF(Vend=1). Missing part and HI correction added in find()
    hiData%real1 = calF
    hiData%real2 = Pstar
    hiData%real3 = w
    
    hbarmin= epsilon(1._kp)
    hbarmax = huge(1._kp)

    hbarstar = zbrent(find_hi_hbar_star,hbarmin,hbarmax,tolzbrent,hiData)    

    hi_hbar_star = hbarstar        
    xi = hi_xi_star(hbarstar,Pstar)
    
    if (present(bfoldstar)) then

       hbarend = hi_parametric_hbar_endinf(xi)

       bfoldstar = -( hi_parametric_efold_primitive(hbarstar,xi) &
            - hi_parametric_efold_primitive(hbarend,xi) )
       
    end if

    if (present(xistar)) then
       xistar = xi
    endif
    
  end function hi_hbar_star


  function find_hi_hbar_star(hbar,hiData)   
    implicit none
    real(kp) :: find_hi_hbar_star
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: hiData
    
    real(kp) :: xi, primStar,potStar,epsOneStar
    real(kp) :: hbarend, w
    real(kp) :: effCalFend, primEnd, potEnd

    calF=hiData%real1
    Pstar = hiData%real2
    w = hiData%real3

    xi = hi_xi_star(hbar,Pstar)
    
    epsOneStar = hi_parametric_epsilon_one(hbar,xi)
!    if (epsOneStar.eq.0._kp) epsOneStar = epsilon(1._kp)

    potStar = hi_norm_parametric_potential(hbar,xi)
!numacc towards minimum
!    if (potStar.eq.0._kp) potStar = epsilon(1._kp)

    primStar = hi_parametric_efold_primitive(hbar,xi)
    
    hbarend = hi_hbar_endinf(xi)

    potEnd = hi_norm_parametric_potential(hbarend,xi)
    
    primEnd = hi_parametric_efold_primitive(hbarend,xi)

    
!the last term is coming from the convertion of JF rhoend -> EF rhoend
    effcalFplusNuEnd = calF - log(potEnd)/(3._kp+3._kp*w) &
         - (1._kp - 3._kp*w)/(6._kp+6._kp*w)*log(1._kp + hbarend*hbarend) &
         + primEnd
    
    find_hi_hbar_star = find_reheat(primStar,effcalFend,w,epsOneStar,potStar)
  
  end function find_hi_hbar_star




!returns k2 given potential parameters, scalar power, and lnRrad
  function hi_hbar_rrad(lnRrad,Pstar,bfoldstar,xistar)    
    implicit none
    real(kp) :: hi_hbar_rrad
    real(kp), intent(in) :: lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar,xistar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF, xi
    real(kp) :: epsOneEnd,potEnd
    real(kp) :: hbarmin, hbarmax, hbarstar
    type(transfert) :: hiData

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = 1._kp
!trick to use the std calfconst function
    potEnd = 1._kp
        
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)
         
!calF(Vend=1). Missing parts added in find()
    hiData%real1 = calF
    hiData%real2 = Pstar
    
    hbarmin= epsilon(1._kp)
    hbarmax = huge(1._kp)
     
    hbarstar = zbrent(find_hi_hbar_rrad,hbarmin,hbarmax,tolzbrent,hiData)

    hi_hbar_rrad = hbarstar
    xi = hi_xi_star(hbarstar,Pstar)
    
    if (present(bfoldstar)) then

       hbarend = hi_parametric_hbar_endinf(xi)

       bfoldstar = -( hi_parametric_efold_primitive(hbarstar,xi) &
            - hi_parametric_efold_primitive(hbarend,xi) )
       
    end if

    if (present(xistar)) then
       xistar = xi
    endif

  end function hi_hbar_rrad

  
  function find_hi_hbar_rrad(hbar,hiData)   
    implicit none
    real(kp) :: find_hi_hbar_rrad
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: hiData

    real(kp) :: xi, primStar,potStar,epsOneStar
    real(kp) :: hbarend
    real(kp) :: effCalFend, primEnd, potEnd
    
    calF=hiData%real1
    Pstar = hiData%real2

    xi = hi_xi_star(hbar,Pstar)
    
    epsOneStar = hi_parametric_epsilon_one(hbar,xi)

    potStar = hi_norm_parametric_potential(hbar,xi)

    primStar = hi_parametric_efold_primitive(hbar,xi)
    
    hbarend = hi_hbar_endinf(xi)

    potEnd = hi_norm_parametric_potential(hbarend,xi)
    
    primEnd = hi_parametric_efold_primitive(hbarend,xi)

    
!no last term as the JF rhoend is included in lnRrad
    effcalFplusNuEnd = calF - 0.25_kp*log(potEnd) + primEnd    
    
    find_hi_hbar_rrad = find_reheat_rrad(primStar,effCalFplusNuEnd,epsOneStar,potStar)
  
  end function find_hi_hbar_rrad
  


!returns hbar, and xistar, given calar power and lnRreh (non-standard,
!we need Pstar)
  function hi_hbar_rreh(lnRreh,Pstar,bfoldstar,xistar)
    implicit none
    real(kp) :: hi_hbar_rreh
    real(kp), intent(in) :: f,lnRreh,Pstar
    real(kp), intent(out), optional :: bfoldstar, xistar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF, xi
    real(kp) :: epsOneEnd,potEnd
    real(kp) :: hbarmin, hbarmax, hbarstar
    type(transfert) :: hiData

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = 1._kp
!trick to use the std calfconst function
    potEnd = 1._kp
    
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

!calF(Vend=1). Missing parts added in find()
    hiData%real1 = calF
    hiData%real2 = Pstar
    
    hbarmin= epsilon(1._kp)
    hbarmax = huge(1._kp)
         
    hbarstar = zbrent(find_hi_hbar_rreh,hbarmin,hbarmax,tolzbrent,hiData)

    hi_hbar_rreh = hbarstar        
    xi = hi_xi_star(hbarstar,Pstar)

    if (present(bfoldstar)) then

       hbarend = hi_parametric_hbar_endinf(xi)

       bfoldstar = -( hi_parametric_efold_primitive(hbarstar,xi) &
            - hi_parametric_efold_primitive(hbarend,xi) )

    endif

    
    if (present(xistar)) then
       xistar = xi
    endif

  end function hi_hbar_rreh

  
  function find_hi_hbar_rreh(hbar,hiData)   
    implicit none
    real(kp) :: find_hi_hbar_rreh
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: hiData

    real(kp) :: xi, primStar,potStar,epsOneStar
    real(kp) :: hbarend
    real(kp) :: effCalFend, primEnd, potEnd
    
    calF=hiData%real1
    Pstar = hiData%real2

    xi = hi_xi_star(hbar,Pstar)
    
    epsOneStar = hi_parametric_epsilon_one(hbar,xi)

    potStar = hi_norm_parametric_potential(hbar,xi)

    primStar = hi_parametric_efold_primitive(hbar,xi)
    
    hbarend = hi_hbar_endinf(xi)

    potEnd = hi_norm_parametric_potential(hbarend,xi)
    
    primEnd = hi_parametric_efold_primitive(hbarend,xi)
    
!The last term is present because we define lnRreh = lnRrad + 1/4
!ln(rhoendbarJF)
    effcalFplusNuEnd = calF - 0.5_kp*log(potEnd) + primEnd &
         -0.5_kp*log(1._kp + hbarend*hbarend)
  
    find_hi_hbar_rreh = find_reheat_rreh(primStar,effCalFplusNuEnd,potStar)

    
  end function find_hi_hbar_rreh
  

  function hi_lnrhoreh_max(Pstar)
    implicit none
    real(kp) :: hi_lnrhoreh_max
    real(kp), intent(in) :: Pstar

    real(kp), parameter :: lnRreh=0._kp

    real(kp) :: xistar
    real(kp) :: hbarstar, hbarend
    real(kp) :: potStar, potEnd
    real(kp) :: epsOneStar, lnRhoEnd, epsOneEnd

    hbarstar = hi_hbar_rreh(lnRreh,Pstar,xistar=xistar)

    epsOneStar = hi_parametric_epsilon_one(hbarstar,xistar)
    potStar = hi_norm_parametric_potential(hbarstar,xistar)

    hbarend = hi_hbar_endinf(xistar)
    potEnd = hi_norm_parametric_potential(hbarend,xistar)
!should be one
    epsOneEnd = hi_parametric_epsilon_one(hbarend,xistar)
    print *,'test',epsoneend
    
    if (.not.slowroll_validity(epsOneStar)) stop 'hi_lnrhoreh_max: slow-roll violated!'

    lnRhoEnd = ln_rho_enhinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    hi_lnrhoreh_max = lnRhoEnd

  end function hi_lnrhoreh_max


end module hireheat
