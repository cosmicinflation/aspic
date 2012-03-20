!MSSM inflation reheating functions in the slow-roll approximations

module mssmireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use mssmisr, only : mssmi_epsilon_one, mssmi_epsilon_two, mssmi_epsilon_three
  use mssmisr, only : mssmi_norm_potential
  use mssmisr, only : mssmi_x_endinf, mssmi_efold_primitive
  use mssmisr, only : mssmi_x_epsilon1_min
  implicit none

  private

  public mssmi_x_star, mssmi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function mssmi_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: mssmi_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: mssmiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = mssmi_x_endinf(alpha,beta)
    epsOneEnd = mssmi_epsilon_one(xEnd,alpha,beta)
    potEnd = mssmi_norm_potential(xEnd,alpha,beta)
    primEnd = mssmi_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    mssmiData%real1 = alpha 
    mssmiData%real2 = beta 
    mssmiData%real3 = xEnd
    mssmiData%real4 = w
    mssmiData%real5 = calF + primEnd

    mini = xEnd

    if (alpha**2/beta>20._kp/9._kp) then
	maxi = mssmi_x_epsilon1_min(alpha,beta)!*(1._kp-epsilon(1._kp)) !local maximum of the potential
    else
	maxi=100._kp*mssmi_x_epsilon1_min(alpha,beta)
    endif



    x = zbrent(find_mssmi_x_star,mini,maxi,tolzbrent,mssmiData)
    mssmi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (mssmi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function mssmi_x_star

  function find_mssmi_x_star(x,mssmiData)   
    implicit none
    real(kp) :: find_mssmi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mssmiData

    real(kp) :: primStar,alpha,beta,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=mssmiData%real1
    beta=mssmiData%real2
    xEnd=mssmiData%real3
    w = mssmiData%real4
    CalFplusprimEnd = mssmiData%real5

    primStar = mssmi_efold_primitive(x,alpha,beta)
    epsOneStar = mssmi_epsilon_one(x,alpha,beta)
    potStar = mssmi_norm_potential(x,alpha,beta)

    find_mssmi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_mssmi_x_star



  function mssmi_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: mssmi_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = mssmi_x_endinf(alpha,beta)
    potEnd  = mssmi_norm_potential(xEnd,alpha,beta)
    epsOneEnd = mssmi_epsilon_one(xEnd,alpha,beta)

!   Trick to return x such that rho_reh=rho_end

    x = mssmi_x_star(alpha,beta,wrad,junk,Pstar)    
    potStar = mssmi_norm_potential(x,alpha,beta)
    epsOneStar = mssmi_epsilon_one(x,alpha,beta)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'mssmi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    mssmi_lnrhoend = lnRhoEnd

  end function mssmi_lnrhoend

  
end module mssmireheat
