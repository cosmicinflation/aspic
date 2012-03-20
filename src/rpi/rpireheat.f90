!R-R^p inflation reheating functions in the slow-roll approximations

module rpireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rpisr, only : rpi_epsilon_one, rpi_epsilon_two, rpi_epsilon_three
  use rpisr, only : rpi_norm_potential
  use rpisr, only : rpi_y_endinf, rpi_efold_primitive
  implicit none

  private

  public rpi_y_star, rpi_lnrhoend 

contains

!returns y =phi/Mp * sqrt(2/3) such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rpi_y_star(p,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rpi_y_star
    real(kp), intent(in) :: p,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp/1000000._kp
    real(kp) :: mini,maxi,calF,y
    real(kp) :: primEnd,epsOneEnd,yend,potEnd

    type(transfert) :: rpiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    yEnd = rpi_y_endinf(p)

    epsOneEnd = rpi_epsilon_one(yEnd,p)
    potEnd = rpi_norm_potential(yEnd,p)

    primEnd = rpi_efold_primitive(yEnd,p)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rpiData%real1 = p
    rpiData%real2 = w
    rpiData%real3 = calF + primEnd

    mini = yEnd

    if (p.eq.1._kp) then !Higgs Inflation Model (HI)

	maxi=100._kp !to avoid numerical explosion

    else

    maxi = log((2._kp*p-1._kp)/(p-1._kp))*(1._kp-1._kp*epsilon(1._kp)) !local maximum of the potential, numerical safety if inflation proceeds from the right to the left

    maxi = log((2._kp*p-1._kp)/(p-1._kp))*(1._kp+1._kp*epsilon(1._kp)) !local maximum of the potential, numerical safety if inflation proceeds from the left to the right

    endif


    y = zbrent(find_rpi_y_star,mini,maxi,tolzbrent,rpiData)
    rpi_y_star = y

    if (present(bfoldstar)) then
       bfoldstar = - (rpi_efold_primitive(y,p) - primEnd)
    endif

  end function rpi_y_star


  function find_rpi_y_star(y,rpiData)   
    implicit none
    real(kp) :: find_rpi_y_star
    real(kp), intent(in) :: y
    type(transfert), optional, intent(inout) :: rpiData

    real(kp) :: primStar,p,w,CalFplusprimEnd,potStar,epsOneStar

    p=rpiData%real1
    w = rpiData%real2
    CalFplusprimEnd = rpiData%real3

    primStar = rpi_efold_primitive(y,p)
    epsOneStar = rpi_epsilon_one(y,p)
    potStar = rpi_norm_potential(y,p)

    find_rpi_y_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rpi_y_star



  function rpi_lnrhoend(p,Pstar) 
    implicit none
    real(kp) :: rpi_lnrhoend
    real(kp), intent(in) :: p,Pstar

    real(kp) :: yEnd, potEnd, epsOneEnd
    real(kp) :: y, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    yEnd = rpi_y_endinf(p)


    potEnd  = rpi_norm_potential(yEnd,p)

    epsOneEnd = rpi_epsilon_one(yEnd,p)


!   Trick to return y such that rho_reh=rho_end

    y = rpi_y_star(p,wrad,junk,Pstar)  

 
    potStar = rpi_norm_potential(y,p)
    epsOneStar = rpi_epsilon_one(y,p)

    PRINT*,'rpi_lnrhoend   :ystar=',y,'  potStar=',potStar,'  epsOneStar=',epsOneStar

    
    if (.not.slowroll_validity(epsOneStar)) stop 'rpi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rpi_lnrhoend = lnRhoEnd

  end function rpi_lnrhoend

  
end module rpireheat
