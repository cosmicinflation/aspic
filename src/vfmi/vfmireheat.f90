module vfmireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use vfmisr, only : vfmi_numacc_betamax, vfmi_x_potmax
  use vfmisr, only : vfmi_norm_potential, vfmi_epsilon_one
  use vfmisr, only : vfmi_x_trajectory, vfmi_efold_primitive, vfmi_x_endinf

  implicit none

!we require at minimum this amount of inflation before the potential
!becomes larger than huge(1.) >> super-planckian...

  real(kp), parameter :: efoldNumAccMin = 120._kp

  private

  public vfmi_x_star, vfmi_lnrhoreh_max
  public vfmi_x_rrad, vfmi_x_rreh
  
contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function vfmi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfold)
    implicit none
    real(kp) :: vfmi_x_star
    real(kp), intent(in) :: alpha,beta,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini, maxi, calF, x
    real(kp) :: primEnd, epsOneEnd, potEnd
    real(kp) :: betamax

    type(transfert) :: vfmiData
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    betamax = vfmi_numacc_betamax(efoldNumAccMin,alpha)
    if (beta.gt.betamax) then
       write(*,*)'vfmi_x_star: beta= betamax= ',beta,betamax
       stop 'beta too large!'
    endif


!should be one    
    epsOneEnd = vfmi_epsilon_one(xEnd,alpha,beta)
    potEnd = vfmi_norm_potential(xEnd,alpha,beta)
!should be zero
    primEnd = vfmi_efold_primitive(xEnd,alpha,beta)


    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    vfmiData%real1 = alpha
    vfmiData%real2 = beta
    vfmiData%real3 = w
    vfmiData%real4 = calF + primEnd

    if (alpha.le.2._kp) then

       mini = xEnd*(1._kp+epsilon(1._kp))
       maxi = vfmi_x_trajectory(-efoldNumAccMin,xend,alpha,beta)
       
    else
       mini = xEnd*(1._kp+epsilon(1._kp))
       maxi = vfmi_x_potmax(alpha,beta)*(1._kp-epsilon(1._kp))
    endif

    x = zbrent(find_vfmi_x_star,mini,maxi,tolFind,vfmiData)
    vfmi_x_star = x

    if (present(bfold)) then
       bfold = -(vfmi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function vfmi_x_star

  function find_vfmi_x_star(x,vfmiData)   
    implicit none
    real(kp) :: find_vfmi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: vfmiData

    real(kp) :: primStar,alpha,beta,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=vfmiData%real1
    beta=vfmiData%real2
    w = vfmiData%real3
    CalFplusPrimEnd = vfmiData%real4

    primStar = vfmi_efold_primitive(x,alpha,beta)
    epsOneStar = vfmi_epsilon_one(x,alpha,beta)
    potStar = vfmi_norm_potential(x,alpha,beta)
    find_vfmi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_vfmi_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function vfmi_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: vfmi_x_rrad
    real(kp), intent(in) :: alpha,beta,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    real(kp) :: betamax
    type(transfert) :: vfmiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    betamax = vfmi_numacc_betamax(efoldNumAccMin,alpha)
    if (beta.gt.betamax) then
       write(*,*)'vfmi_x_rrad: beta= betamax= ',beta,betamax
       stop 'beta too large!'
    endif
    
    epsOneEnd = vfmi_epsilon_one(xEnd,alpha,beta)
    potEnd = vfmi_norm_potential(xEnd,alpha,beta)
    primEnd = vfmi_efold_primitive(xEnd,alpha,beta)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    vfmiData%real1 = alpha
    vfmiData%real2 = beta
    vfmiData%real3 = calF + primEnd

    if (alpha.le.2._kp) then

       mini = xEnd*(1._kp+epsilon(1._kp))
       maxi = vfmi_x_trajectory(-efoldNumAccMin,xend,alpha,beta)
       
    else
       mini = xEnd*(1._kp+epsilon(1._kp))
       maxi = vfmi_x_potmax(alpha,beta)*(1._kp-epsilon(1._kp))
    endif

    x = zbrent(find_vfmi_x_rrad,mini,maxi,tolFind,vfmiData)
    vfmi_x_rrad = x

    if (present(bfold)) then
       bfold = -(vfmi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function vfmi_x_rrad

  function find_vfmi_x_rrad(x,vfmiData)   
    implicit none
    real(kp) :: find_vfmi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: vfmiData

    real(kp) :: primStar,alpha,beta,CalFplusPrimEnd,potStar,epsOneStar

    alpha=vfmiData%real1
    beta=vfmiData%real2
    CalFplusPrimEnd = vfmiData%real3

    primStar = vfmi_efold_primitive(x,alpha,beta)
    epsOneStar = vfmi_epsilon_one(x,alpha,beta)
    potStar = vfmi_norm_potential(x,alpha,beta)

    find_vfmi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_vfmi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function vfmi_x_rreh(alpha,beta,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: vfmi_x_rreh
    real(kp), intent(in) :: alpha,beta,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    real(kp) :: betamax
    type(transfert) :: vfmiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    betamax = vfmi_numacc_betamax(efoldNumAccMin,alpha)
    if (beta.gt.betamax) then
       write(*,*)'vfmi_x_rreh: beta= betamax= ',beta,betamax
       stop 'beta too large!'
    endif
    
    epsOneEnd = vfmi_epsilon_one(xEnd,alpha,beta)
    potEnd = vfmi_norm_potential(xEnd,alpha,beta)
    primEnd = vfmi_efold_primitive(xEnd,alpha,beta)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    vfmiData%real1 = alpha
    vfmiData%real2 = beta
    vfmiData%real3 = calF + primEnd

    if (alpha.le.2._kp) then

       mini = xEnd*(1._kp+epsilon(1._kp))
       maxi = vfmi_x_trajectory(-efoldNumAccMin,xend,alpha,beta)
       
    else
       mini = xEnd*(1._kp+epsilon(1._kp))
       maxi = vfmi_x_potmax(alpha,beta)*(1._kp-epsilon(1._kp))
    endif

    x = zbrent(find_vfmi_x_rreh,mini,maxi,tolFind,vfmiData)
    vfmi_x_rreh = x

    if (present(bfold)) then
       bfold = -(vfmi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function vfmi_x_rreh

  function find_vfmi_x_rreh(x,vfmiData)   
    implicit none
    real(kp) :: find_vfmi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: vfmiData

    real(kp) :: primStar,alpha,beta,CalFplusPrimEnd,potStar

    alpha=vfmiData%real1
    beta=vfmiData%real2
    CalFplusPrimEnd = vfmiData%real3

    primStar = vfmi_efold_primitive(x,alpha,beta)
    potStar = vfmi_norm_potential(x,alpha,beta)

    find_vfmi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_vfmi_x_rreh


  function vfmi_lnrhoreh_max(alpha,beta,xend,Pstar) 
    implicit none
    real(kp) :: vfmi_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    potEnd  = vfmi_norm_potential(xEnd,alpha,beta)
    epsOneEnd = vfmi_epsilon_one(xEnd,alpha,beta)
       
    x = vfmi_x_star(alpha,beta,xend,wrad,junk,Pstar)

    potStar = vfmi_norm_potential(x,alpha,beta)
    epsOneStar = vfmi_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'vfmi_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    vfmi_lnrhoreh_max = lnRhoEnd

  end function vfmi_lnrhoreh_max



end module vfmireheat
