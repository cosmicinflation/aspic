!Generalized Mixed inflation reheating functions in the slow-roll approximations

module gmireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use gmisr, only : gmi_epsilon_one, gmi_epsilon_two, gmi_epsilon_three
  use gmisr, only : gmi_norm_potential, gmi_x_endinf
  use gmisr, only :  gmi_efold_primitive
  implicit none

  private

  public gmi_x_star, gmi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function gmi_x_star(alpha,p,q,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: gmi_x_star
    real(kp), intent(in) :: alpha,p,q,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: gmiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = gmi_x_endinf(alpha,p,q)
    epsOneEnd = gmi_epsilon_one(xEnd,alpha,p,q)
    potEnd = gmi_norm_potential(xEnd,alpha,p,q)
    primEnd = gmi_efold_primitive(xEnd,alpha,p,q) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    gmiData%real1 = alpha
    gmiData%real2 = p
    gmiData%real3 = q
    gmiData%real4 = xEnd
    gmiData%real5 = w
    gmiData%real6 = calF + primEnd

    mini = xend
    maxi = xend*10._kp**(6._kp)

    x = zbrent(find_gmi_x_star,mini,maxi,tolzbrent,gmiData)
    gmi_x_star = x  



    if (present(bfoldstar)) then
       bfoldstar = - (gmi_efold_primitive(x,alpha,p,q) - primEnd)
    endif

  end function gmi_x_star

  function find_gmi_x_star(x,gmiData)   
    implicit none
    real(kp) :: find_gmi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gmiData

    real(kp) :: primStar,alpha,p,q,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=gmiData%real1
    p=gmiData%real2
    q=gmiData%real3
    xEnd=gmiData%real4
    w = gmiData%real5
    CalFplusprimEnd = gmiData%real6

    primStar = gmi_efold_primitive(x,alpha,p,q)
    epsOneStar = gmi_epsilon_one(x,alpha,p,q)
    potStar = gmi_norm_potential(x,alpha,p,q)

    find_gmi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_gmi_x_star



  function gmi_lnrhoend(alpha,p,q,Pstar) 
    implicit none
    real(kp) :: gmi_lnrhoend
    real(kp), intent(in) :: alpha,p,q,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = gmi_x_endinf(alpha,p,q)
    potEnd  = gmi_norm_potential(xEnd,alpha,p,q)
    epsOneEnd = gmi_epsilon_one(xEnd,alpha,p,q)

!   Trick to return x such that rho_reh=rho_end

    x = gmi_x_star(alpha,p,q,wrad,junk,Pstar)    
    potStar = gmi_norm_potential(x,alpha,p,q)
    epsOneStar = gmi_epsilon_one(x,alpha,p,q)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'gmi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    gmi_lnrhoend = lnRhoEnd

  end function gmi_lnrhoend

  
end module gmireheat
