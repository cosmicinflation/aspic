!Reheating functions for Non-Minimal Large Field Inflation
!
!We here assume that the coupling constant (lambdabar in Encyclopedia
!Inflationaris) is free and we trade it for the potential
!normalization M. As such, the reheating equations are much simpler
!than for HI where lambdabar was fixed.
!
!Still, the reheating equation are non-standard due to the fact that
!lnRrad has to be evaluated in the Jordan Frame (see Encyclopedia
!Inflationaris). As such, expressed in terms of rhoend (Einstein
!Frame), there is an extra conformal factor in
!Omega4End=(1+hbarend^2)^2, see also SI. Therefore, energies are
!returned in the Jordan Frame and this matters for determining the
!maximal possible value of rhoreh.
!
!Mind the new definition:
!
! lnRreh = lnRrad + (1/4) ln(rhoendJF)
!
! to be consistent with lnRrad defined with JF energy densities and to
! ensure than all reheating-related energies are Jordan Frame-defined.
!
module nmlficomreh
  use infprec, only :  kp, tolkp, toldp, pi, transfert
  use inftools, only : zbrent
  use srreheat, only : display
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : pi, ln_rho_endinf, ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh

  use nmlficommon, only : nmlfi_x, nmlfi_hbar, hbarBig, hbarSmall
  use nmlficommon, only : nmlfi_parametric_epsilon_one, nmlfi_norm_parametric_potential
  use nmlficommon, only : nmlfi_parametric_efold_primitive
  
  implicit none

  private

!take as input the parametric field hbar  
  public nmlfi_lnrhoreh_max
  public nmlfi_hbar_star, nmlfi_hbar_rrad, nmlfi_hbar_rreh

!for convenience only
  public nmlfi_x_star, nmlfi_x_rrad, nmlfi_x_rreh

contains


!Convenience functions: inefficient wrapper taking as input the
!canonical field. For computations it is better to call the parametric
!functions directly, they work in terms of hbar rather than x to
!prevent unnecessary numerical inversions
  function nmlfi_x_star(xi,p,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nmlfi_x_star
    real(kp), intent(in) :: xi,p,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar, hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

    nmlfi_x_star = nmlfi_x(hbarstar,xi)
    
  end function nmlfi_x_star



  function nmlfi_x_rrad(xi,p,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nmlfi_x_rrad
    real(kp), intent(in) :: xi,p,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar,hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)

    nmlfi_x_rrad = nmlfi_x(hbarstar,xi)

  end function nmlfi_x_rrad

    
  function nmlfi_x_rreh(xi,p,xend,lnRreh,bfoldstar)
    implicit none
    real(kp) :: nmlfi_x_rreh
    real(kp), intent(in) :: xi,p,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar,hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi_hbar_rreh(xi,p,hbarend,lnRreh,bfoldstar)

    nmlfi_x_rreh = nmlfi_x(hbarstar,xi)
    
  end function nmlfi_x_rreh
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!the parametric reheating functions, the ones you should use !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function nmlfi_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nmlfi_hbar_star
    real(kp), intent(in) :: xi,p,hbarend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,primEnd,epsOneEnd,potEnd, lnOmega4End
    real(kp) :: mini,maxi,hbarstar
    type(transfert) :: nmlfiData
  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = nmlfi_parametric_epsilon_one(hbarend,xi,p)
    potEnd = nmlfi_norm_parametric_potential(hbarend,xi,p)

    primEnd = nmlfi_parametric_efold_primitive(hbarend,xi,p)
!ln(Omega^4) = ln[(1+ hbar^2)^2]
    lnOmega4End = 2._kp*log(1._kp + hbarend*hbarend)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd,lnOmega4End)
    
    nmlfiData%real1 = xi
    nmlfiData%real2 = p
    nmlfiData%real3 = w
    nmlfiData%real4 = calF + primEnd

    mini = hbarend
    maxi = hbarBig
    
    hbarstar = zbrent(find_nmlfi_hbar_star,mini,maxi,tolzbrent,nmlfiData)
    
    if (present(bfoldstar)) then
       bfoldstar = - (nmlfi_parametric_efold_primitive(hbarstar,xi,p) - primEnd)
    endif

  end function nmlfi_hbar_star


  function find_nmlfi_hbar_star(hbar,nmlfiData)
    implicit none
    real(kp) :: find_nmlfi_hbar_star
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: nmlfiData
    
    real(kp) :: xi, p, w, CalFplusprimEnd
    real(kp) :: primStar, potStar,epsOneStar

    xi =  nmlfiData%real1
    p =  nmlfiData%real2
    w = nmlfiData%real3
    CalFplusprimEnd = nmlfiData%real4

    primStar = nmlfi_parametric_efold_primitive(hbar,xi,p)
    epsOneStar = nmlfi_parametric_epsilon_one(hbar,xi,p)
    potStar = nmlfi_norm_parametric_potential(hbar,xi,p)
    
    find_nmlfi_hbar_star = find_reheat(primStar,CalFplusprimEnd,w,epsOneStar,potStar)    
    
  end function find_nmlfi_hbar_star
  


  function nmlfi_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nmlfi_hbar_rrad
    real(kp), intent(in) :: xi,p,hbarend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF
    real(kp) :: primEnd,epsOneEnd,potEnd
    real(kp) :: mini,maxi,hbarstar
    type(transfert) :: nmlfiData

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = nmlfi_parametric_epsilon_one(hbarend,xi,p)
    potEnd = nmlfi_norm_parametric_potential(hbarend,xi,p)

    primEnd = nmlfi_parametric_efold_primitive(hbarend,xi,p)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)
         
    nmlfiData%real1 = xi
    nmlfiData%real2 = p
    nmlfiData%real3 = calF + primEnd
    
    mini = hbarend
    maxi = hbarBig
     
    hbarstar = zbrent(find_nmlfi_hbar_rrad,mini,maxi,tolzbrent,nmlfiData)

    nmlfi_hbar_rrad = hbarstar
    
    if (present(bfoldstar)) then
       bfoldstar = -( nmlfi_parametric_efold_primitive(hbarstar,xi,p)- primEnd )
    end if

  end function nmlfi_hbar_rrad

  
  function find_nmlfi_hbar_rrad(hbar,nmlfiData)   
    implicit none
    real(kp) :: find_nmlfi_hbar_rrad
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: nmlfiData

    real(kp) :: xi, p
    real(kp) :: primStar,CalFplusprimEnd,potStar,epsOneStar
    

    xi = nmlfiData%real1
    p = nmlfiData%real2
    CalFplusprimEnd = nmlfiData%real3

    primStar = nmlfi_parametric_efold_primitive(hbar,xi,p)
    epsOneStar = nmlfi_parametric_epsilon_one(hbar,xi,p)
    potStar = nmlfi_norm_parametric_potential(hbar,xi,p)
    
    find_nmlfi_hbar_rrad = find_reheat_rrad(primStar,CalFplusprimEnd,epsOneStar,potStar)
  
  end function find_nmlfi_hbar_rrad
  


  function nmlfi_hbar_rreh(xi,p,hbarend,lnRreh,bfoldstar)
    implicit none
    real(kp) :: nmlfi_hbar_rreh
    real(kp), intent(in) :: xi,p,hbarend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF
    real(kp) :: primEnd,epsOneEnd,potEnd,lnOmega4End
    real(kp) :: mini,maxi,hbarstar
    type(transfert) :: nmlfiData

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = nmlfi_parametric_epsilon_one(hbarend,xi,p)
    potEnd = nmlfi_norm_parametric_potential(hbarend,xi,p)

    primEnd = nmlfi_parametric_efold_primitive(hbarend,xi,p)

!ln(Omega^4) = ln[(1+ hbar^2)^2]
    lnOmega4End = 2._kp*log(1._kp + hbarend*hbarend)
    
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd,lnOmega4End)

    nmlfiData%real1 = xi
    nmlfiData%real2 = p
    nmlfiData%real3 = calF + primEnd

    mini= hbarend
    maxi = hbarBig
         
    hbarstar = zbrent(find_nmlfi_hbar_rreh,mini,maxi,tolzbrent,nmlfiData)

    nmlfi_hbar_rreh = hbarstar        

    if (present(bfoldstar)) then
       bfoldstar = -( nmlfi_parametric_efold_primitive(hbarstar,xi,p) - primEnd )
    endif
    
  end function nmlfi_hbar_rreh

  
  function find_nmlfi_hbar_rreh(hbar,nmlfiData)   
    implicit none
    real(kp) :: find_nmlfi_hbar_rreh
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: nmlfiData

    real(kp) :: xi, p, primStar,potStar,CalFplusprimEnd
    

    xi = nmlfiData%real1
    p = nmlfiData%real2
    CalFplusprimEnd = nmlfiData%real3

    primStar = nmlfi_parametric_efold_primitive(hbar,xi,p)
    potStar = nmlfi_norm_parametric_potential(hbar,xi,p)        
    
    find_nmlfi_hbar_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
    
  end function find_nmlfi_hbar_rreh
  
!jordan frame
  function nmlfi_lnrhoreh_max(xi,p,hbarend,Pstar)
    implicit none
    real(kp) :: nmlfi_lnrhoreh_max
    real(kp), intent(in) :: xi,p,hbarend,Pstar

!trick to return rhoreh=rhoend    
    real(kp), parameter :: lnRrad=0._kp

    real(kp) :: xistar
    real(kp) :: hbarstar
    real(kp) :: potStar, potEnd
    real(kp) :: epsOneStar, lnRhoEnd, epsOneEnd, lnOmega4End

    potEnd = nmlfi_norm_parametric_potential(hbarend,xi,p)
!should be one
    epsOneEnd = nmlfi_parametric_epsilon_one(hbarend,xi,p)

    hbarstar = nmlfi_hbar_rrad(xi,p,hbarend,lnRrad,Pstar)
    
    epsOneStar = nmlfi_parametric_epsilon_one(hbarstar,xi,p)
    potStar = nmlfi_norm_parametric_potential(hbarstar,xi,p)
    lnOmega4End = 2._kp*log(1._kp+hbarend*hbarend)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'nmlfi_lnrhoreh_max: slow-roll violated!'

!Jordan frame
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar,lnOmega4End)
    
    nmlfi_lnrhoreh_max = lnRhoEnd

  end function nmlfi_lnrhoreh_max

  
end module nmlficomreh
