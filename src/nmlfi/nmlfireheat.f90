!Reheating functions for Higgs Inflation.
!
!These are non-standard as M^4, the potential normalization, is
!completely determined by the only model parameters xi and xstar (or
!hbarstar). So the reheating equations and potential normalisation are
!two coupled non-linear algebraic equations to be solved
!simultaneously given the reheating parameter and xi. The solutions of
!wnmlfich are hbarstar (->xstar) and xistar. The equation for the CMB to
!M^4 normalisation cannot analytically be solved for xi, but assuming
!the Higgs vev to be zero, there is an analytical solution. Tnmlfis
!solution is used to start a recursion allowing to determinte the
!exact value of xi given hbarstar at non-vanisnmlfing Higgs vev. Tnmlfis
!value of xi is everytime recomputed when zbrenting the reheating
!equation.
!
!Let us also notice that the reheating equation are non-standard due
!to the fact that lnRrad has to be evaluated in the Jordan Frame (see
!Encyclopedia Inflationaris). As such, expressed in terms of rhoend
!(Einstein Frame), there is an extra conformal factor in
!Omega4End=(1+hbarend^2)^2. Therefore, energies are returned in the
!Jordan Frame and tnmlfis matters for determining the maximal possible
!value of rhoreh.
!
!Mind the new definition:
!
! lnRreh = lnRrad + (1/4) ln(rhoendJF)
!
! to be consistent with lnRrad defined with JF energy densities and to
! ensure than all reheating-related energies are Jordan Frame-defined.
!
module nmlfireheat
  use infprec, only :  kp, tolkp, toldp, pi, transfert
  use cosmopar, only : HiggsCoupling
  use inftools, only : zbrent
  use srreheat, only : display
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : pi, ln_rho_endinf, ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh

  use nmlficommon, only : vev2, nmlfi_x, nmlfi_hbar, hbarBig, hbarSmall
  use nmlficommon, only : nmlfi_parametric_epsilon_one, nmlfi_norm_parametric_potential
  use nmlficommon, only : nmlfi_parametric_efold_primitive, nmlfi_parametric_hbar_endinf
  
  implicit none

!see cosmopar.f90 for Higgs Lagrangian convention  
  real(kp), parameter :: lambda = HiggsCoupling
  real(kp), parameter :: pstarconst = 96._kp*pi*pi/lambda

  real(kp), parameter :: onetnmlfird = 1._kp/3._kp
  real(kp), parameter :: twotnmlfird = 2._kp/3._kp

  private
  
  public nmlfi_lnrhoreh_max, HiggsCoupling
  public nmlfi_hbar_star, nmlfi_hbar_rrad, nmlfi_hbar_rreh
  public nmlfi_x_star, nmlfi_x_rrad, nmlfi_x_rreh

contains


!xistar is an output
  function nmlfi_x_star(w,lnRhoReh,Pstar,bfoldstar,xistar)
    implicit none
    real(kp) :: nmlfi_x_star
    real(kp), intent(in) :: w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfoldstar,xistar
    
    real(kp) :: hbarstar

    hbarstar = nmlfi_hbar_star(w,lnRhoReh,Pstar,bfoldstar,xistar)

    nmlfi_x_star = nmlfi_x(hbarstar,xistar)
    
  end function nmlfi_x_star



  function nmlfi_x_rrad(lnRrad,Pstar,bfoldstar,xistar)    
    implicit none
    real(kp) :: nmlfi_x_rrad
    real(kp), intent(in) :: lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar,xistar
    
    real(kp) :: hbarstar

    hbarstar = nmlfi_hbar_rrad(lnRrad,Pstar,bfoldstar,xistar)

    nmlfi_x_rrad = nmlfi_x(hbarstar,xistar)

  end function nmlfi_x_rrad

    
!calling arguments are non-standard and require Pstar
  function nmlfi_x_rreh(lnRreh,Pstar,bfoldstar,xistar)
    implicit none
    real(kp) :: nmlfi_x_rreh
    real(kp), intent(in) :: lnRreh,Pstar
    real(kp), intent(out), optional :: bfoldstar,xistar
    
    real(kp) :: hbarstar

    hbarstar = nmlfi_hbar_rreh(lnRreh,Pstar,bfoldstar,xistar)

    nmlfi_x_rreh = nmlfi_x(hbarstar,xistar)
    
  end function nmlfi_x_rreh



!return xistar, matcnmlfing the CMB normalisation, from any input values
!of hbar (not necessarily solution of the reheating equation that we
!refer to as hbarstar). We guess first with the solution by using the
!analytical one for vev=0 and then recurse to get the right CMB
!normalisation at non-vanisnmlfing vev.
  recursive function nmlfi_xi_star(hbar,Pstar,xiprev) result (xistar)
    implicit none
    real(kp) :: xistar
    real(kp), intent(in) :: hbar,Pstar
    real(kp), intent(in), optional :: xiprev
    real(kp) :: epsstar, potstar

    complex(kp) :: hbar2
    real(kp) :: C,Cstar,xistarzero
    
    real(kp), parameter :: tolfind = 10._kp*epsilon(1._kp)
    logical, parameter :: debug = .false.

    C = pstarconst * Pstar

    hbar2 = cmplx(hbar*hbar,0._kp,kp)
    
    if (.not.present(xiprev)) then
!initial guess, the only positive solution defined for all hbar and
!obtained by assuming HiggsVev = 0

       xistar = real (&
            (6*2**onetnmlfird*hbar2**4)/(1728*C**2*hbar2**3 + 8640*C**2*hbar2**4 + 17280*C**2*hbar2**5 &
            + 17280*C**2*hbar2**6 + 8640*C**2*hbar2**7 + 1728*C**2*hbar2**8 &
            + sqrt(-11943936*hbar2**12*(C + 2*C*hbar2 + C*hbar2**2)**3 + (1728*C**2*hbar2**3 &
            + 8640*C**2*hbar2**4 + 17280*C**2*hbar2**5 + 17280*C**2*hbar2**6 + 8640*C**2*hbar2**7 &
            + 1728*C**2*hbar2**8)**2))**onetnmlfird + (1728*C**2*hbar2**3 + 8640*C**2*hbar2**4 &
            + 17280*C**2*hbar2**5 + 17280*C**2*hbar2**6 + 8640*C**2*hbar2**7 + 1728*C**2*hbar2**8 &
            + sqrt(-11943936*hbar2**12*(C + 2*C*hbar2 + C*hbar2**2)**3 + (1728*C**2*hbar2**3 &
            + 8640*C**2*hbar2**4 + 17280*C**2*hbar2**5 + 17280*C**2*hbar2**6 + 8640*C**2*hbar2**7 &
            + 1728*C**2*hbar2**8)**2))**onetnmlfird/(24._kp*2**onetnmlfird*(C + 2*C*hbar2 + C*hbar2**2)) &
            ,kp )
       
    else

       xistar = xiprev

    endif
       
        
!Cstar = W(hbar)/eps1(hbar)/xi^2
    potstar = nmlfi_norm_parametric_potential(hbar,xistar)
    epsstar = nmlfi_parametric_epsilon_one(hbar,xistar)    
    Cstar = potstar/(xistar*xistar*epsstar)

!new guess for recursion
    xistar = sqrt(potstar/(C*epsstar))
    
!stop recursion when Cstar matches the actual COBE normalisation C    
    if ((2._kp*abs(Cstar-C)/(Cstar+C)).gt.tolfind) then
       if (debug) then
          write(*,*)'nmlfi_xi_star:'
          write(*,*)'tol(C)= ',2._kp*abs(Cstar-C)/(Cstar+C)
          write(*,*)'recursing on xistar = ',xistar
          write(*,*)
       endif
       xistar = nmlfi_xi_star(hbar,Pstar,xistar)
    endif
    
  end function nmlfi_xi_star


!returns hbarstar, and xistar, given the scalar power, wreh and
!lnrhoreh
  function nmlfi_hbar_star(w,lnRhoReh,Pstar,bfoldstar,xistar)
    implicit none
    real(kp) :: nmlfi_hbar_star
    real(kp), intent(in) :: lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar, xistar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF, xi, hbarend, hbarstar
    real(kp) :: epsOneEnd,potEnd

    real(kp) :: hbarmin,hbarmax

    type(transfert) :: nmlfiData
  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = 1._kp
!trick to use the std calfconst function not including the term in potend
    potEnd = 1._kp

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)
         
!calF(Vend=1). Missing part and HI correction added in find()
    nmlfiData%real1 = calF
    nmlfiData%real2 = Pstar
    nmlfiData%real3 = w
    
    hbarmin= hbarSmall
    hbarmax = hbarBig

    hbarstar = zbrent(find_nmlfi_hbar_star,hbarmin,hbarmax,tolzbrent,nmlfiData)    

    nmlfi_hbar_star = hbarstar        
    xi = nmlfi_xi_star(hbarstar,Pstar)
    
    if (present(bfoldstar)) then

       hbarend = nmlfi_parametric_hbar_endinf(xi)

       bfoldstar = -( nmlfi_parametric_efold_primitive(hbarstar,xi) &
            - nmlfi_parametric_efold_primitive(hbarend,xi) )
       
    end if

    if (present(xistar)) then
       xistar = xi
    endif
    
  end function nmlfi_hbar_star


  function find_nmlfi_hbar_star(hbar,nmlfiData)   
    implicit none
    real(kp) :: find_nmlfi_hbar_star
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: nmlfiData
    
    real(kp) :: xi, primStar,potStar,epsOneStar
    real(kp) :: hbarend, w
    real(kp) :: calF, effcalFplusNuEnd, primEnd, potEnd, Ps

    calF=nmlfiData%real1
    Ps = nmlfiData%real2
    w = nmlfiData%real3

    xi = nmlfi_xi_star(hbar,Ps)
    
    epsOneStar = nmlfi_parametric_epsilon_one(hbar,xi)

    potStar = nmlfi_norm_parametric_potential(hbar,xi)

    primStar = nmlfi_parametric_efold_primitive(hbar,xi)
    
    hbarend = nmlfi_parametric_hbar_endinf(xi)

    potEnd = nmlfi_norm_parametric_potential(hbarend,xi)
    
    primEnd = nmlfi_parametric_efold_primitive(hbarend,xi)

!the last term is coming from the convertion of JF rhoend -> EF rhoend
    effcalFplusNuEnd = calF - log(potEnd)/(3._kp+3._kp*w) &
         - (1._kp - 3._kp*w)/(6._kp+6._kp*w)*log(1._kp + hbarend*hbarend) &
         + primEnd

    
    find_nmlfi_hbar_star = find_reheat(primStar,effcalFplusNuEnd,w,epsOneStar,potStar)
  
  end function find_nmlfi_hbar_star




!returns k2 given potential parameters, scalar power, and lnRrad
  function nmlfi_hbar_rrad(lnRrad,Pstar,bfoldstar,xistar)    
    implicit none
    real(kp) :: nmlfi_hbar_rrad
    real(kp), intent(in) :: lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar,xistar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF, xi
    real(kp) :: epsOneEnd,potEnd
    real(kp) :: hbarmin, hbarmax, hbarstar, hbarend
    type(transfert) :: nmlfiData

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = 1._kp
!trick to use the std calfconst function
    potEnd = 1._kp
        
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)
         
!calF(Vend=1). Missing parts added in find()
    nmlfiData%real1 = calF
    nmlfiData%real2 = Pstar
    
    hbarmin= sqrt(epsilon(1._kp))
    hbarmax = hbarBig
     
    hbarstar = zbrent(find_nmlfi_hbar_rrad,hbarmin,hbarmax,tolzbrent,nmlfiData)

    nmlfi_hbar_rrad = hbarstar
    xi = nmlfi_xi_star(hbarstar,Pstar)
    
    if (present(bfoldstar)) then

       hbarend = nmlfi_parametric_hbar_endinf(xi)

       bfoldstar = -( nmlfi_parametric_efold_primitive(hbarstar,xi) &
            - nmlfi_parametric_efold_primitive(hbarend,xi) )
       
    end if

    if (present(xistar)) then
       xistar = xi
    endif

  end function nmlfi_hbar_rrad

  
  function find_nmlfi_hbar_rrad(hbar,nmlfiData)   
    implicit none
    real(kp) :: find_nmlfi_hbar_rrad
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: nmlfiData

    real(kp) :: xi, primStar,potStar,epsOneStar
    real(kp) :: hbarend
    real(kp) :: calF, effcalFplusNuEnd, primEnd, potEnd, Ps
    
    calF=nmlfiData%real1
    Ps = nmlfiData%real2

    xi = nmlfi_xi_star(hbar,Ps)
    
    epsOneStar = nmlfi_parametric_epsilon_one(hbar,xi)

    potStar = nmlfi_norm_parametric_potential(hbar,xi)

    primStar = nmlfi_parametric_efold_primitive(hbar,xi)
    
    hbarend = nmlfi_parametric_hbar_endinf(xi)

    potEnd = nmlfi_norm_parametric_potential(hbarend,xi)
    
    primEnd = nmlfi_parametric_efold_primitive(hbarend,xi)

    
!no last term as the JF rhoend is included in lnRrad
    effcalFplusNuEnd = calF - 0.25_kp*log(potEnd) + primEnd    
    
    find_nmlfi_hbar_rrad = find_reheat_rrad(primStar,effCalFplusNuEnd,epsOneStar,potStar)
  
  end function find_nmlfi_hbar_rrad
  


!returns hbar, and xistar, given calar power and lnRreh (non-standard,
!we need Pstar)
  function nmlfi_hbar_rreh(lnRreh,Pstar,bfoldstar,xistar)
    implicit none
    real(kp) :: nmlfi_hbar_rreh
    real(kp), intent(in) :: lnRreh,Pstar
    real(kp), intent(out), optional :: bfoldstar, xistar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF, xi
    real(kp) :: epsOneEnd,potEnd
    real(kp) :: hbarmin, hbarmax, hbarstar, hbarend
    type(transfert) :: nmlfiData

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = 1._kp
!trick to use the std calfconst function
    potEnd = 1._kp
    
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

!calF(Vend=1). Missing parts added in find()
    nmlfiData%real1 = calF
    nmlfiData%real2 = Pstar
    
    hbarmin= 0.1_kp
    hbarmax = hbarBig
         
    hbarstar = zbrent(find_nmlfi_hbar_rreh,hbarmin,hbarmax,tolzbrent,nmlfiData)

    nmlfi_hbar_rreh = hbarstar        
    xi = nmlfi_xi_star(hbarstar,Pstar)

    if (present(bfoldstar)) then

       hbarend = nmlfi_parametric_hbar_endinf(xi)

       bfoldstar = -( nmlfi_parametric_efold_primitive(hbarstar,xi) &
            - nmlfi_parametric_efold_primitive(hbarend,xi) )

    endif

    
    if (present(xistar)) then
       xistar = xi
    endif

  end function nmlfi_hbar_rreh

  
  function find_nmlfi_hbar_rreh(hbar,nmlfiData)   
    implicit none
    real(kp) :: find_nmlfi_hbar_rreh
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: nmlfiData

    real(kp) :: xi, primStar,potStar,epsOneStar
    real(kp) :: hbarend
    real(kp) :: calF,effcalFplusNuEnd, primEnd, potEnd, Ps
    
    calF=nmlfiData%real1
    Ps = nmlfiData%real2

    xi = nmlfi_xi_star(hbar,Ps)
    
    epsOneStar = nmlfi_parametric_epsilon_one(hbar,xi)

    potStar = nmlfi_norm_parametric_potential(hbar,xi)

    primStar = nmlfi_parametric_efold_primitive(hbar,xi)
    
    hbarend = nmlfi_parametric_hbar_endinf(xi)

    potEnd = nmlfi_norm_parametric_potential(hbarend,xi)
    
    primEnd = nmlfi_parametric_efold_primitive(hbarend,xi)
    
!The last term is present because we define lnRreh = lnRrad + 1/4
!ln(rhoendbarJF) and we have to call it inside find_nmlfi_hbar_rreh.
    effcalFplusNuEnd = calF - 0.5_kp*log(potEnd) + primEnd &
         -0.5_kp*log(1._kp + hbarend*hbarend)    
    
    find_nmlfi_hbar_rreh = find_reheat_rreh(primStar,effCalFplusNuEnd,potStar)

    
  end function find_nmlfi_hbar_rreh
  
!jordan frame
  function nmlfi_lnrhoreh_max(Pstar)
    implicit none
    real(kp) :: nmlfi_lnrhoreh_max
    real(kp), intent(in) :: Pstar

!trick to return rhoreh=rhoend    
    real(kp), parameter :: lnRrad=0._kp

    real(kp) :: xistar
    real(kp) :: hbarstar, hbarend
    real(kp) :: potStar, potEnd
    real(kp) :: epsOneStar, lnRhoEnd, epsOneEnd, lnOmega4End

    hbarstar = nmlfi_hbar_rrad(lnRrad,Pstar,xistar=xistar)
    
    epsOneStar = nmlfi_parametric_epsilon_one(hbarstar,xistar)
    potStar = nmlfi_norm_parametric_potential(hbarstar,xistar)

    hbarend = nmlfi_parametric_hbar_endinf(xistar)
    potEnd = nmlfi_norm_parametric_potential(hbarend,xistar)
!should be one
    epsOneEnd = nmlfi_parametric_epsilon_one(hbarend,xistar)

    lnOmega4End = 2._kp*log(1._kp+hbarend*hbarend)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'nmlfi_lnrhoreh_max: slow-roll violated!'
!Jordan frame
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar,lnOmega4End)
    
    nmlfi_lnrhoreh_max = lnRhoEnd

  end function nmlfi_lnrhoreh_max

  



end module nmlfireheat
