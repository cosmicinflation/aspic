!R-R^p inflation reheating functions in the slow-roll approximations

module rpireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rpisr, only : rpi_epsilon_one, rpi_epsilon_two, rpi_epsilon_three
  use rpisr, only : rpi_norm_potential, rpi_x_potmax
  use rpisr, only : rpi_x_endinf, rpi_efold_primitive
  implicit none

  private

  public rpi_x_star, rpi_lnrhoend 

contains

!returns y =phi/Mp * sqrt(2/3) such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rpi_x_star(p,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rpi_x_star
    real(kp), intent(in) :: p,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=epsilon(1._kp)
    real(kp) :: mini,maxi,calF,y,yVmax
    real(kp) :: primEnd,epsOneEnd,yend,potEnd

    real(kp), parameter :: maxiMax = 100._kp

    type(transfert) :: rpiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    yEnd = rpi_x_endinf(p)

    epsOneEnd = rpi_epsilon_one(yEnd,p)
    potEnd = rpi_norm_potential(yEnd,p)

    primEnd = rpi_efold_primitive(yEnd,p)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rpiData%real1 = p
    rpiData%real2 = w
    rpiData%real3 = calF + primEnd
    

    if (p.eq.1._kp) then !Higgs Inflation Model (HI)

       mini = yEnd
       maxi=maxiMax !to avoid numerical explosion

    else

       yVmax = rpi_x_potmax(p)

       if (yend.lt.yVmax) then

          mini = yEnd
!local maximum of the potential, numerical safety if inflation proceeds from the right to the left
          maxi = yVmax - epsilon(1._kp)
!*(1._kp-1._kp*epsilon(1._kp))

       else
!local maximum of the potential, numerical safety if inflation proceeds from the left to the right
          mini = yVmax + epsilon(1._kp)
!*(1._kp+1._kp*epsilon(1._kp))
          maxi = maxiMax

          stop 'rpi_x_star: not implemented'

       endif

    endif


    y = zbrent(find_rpi_x_star,mini,maxi,tolzbrent,rpiData)

!consistency check
    if ((p.ne.1._kp).and.(yend.lt.yVmax).and.(y.le.yend)) stop 'rpi_x_star: out of numerical accuracy!'

    rpi_x_star = y

    if (present(bfoldstar)) then
       bfoldstar = - (rpi_efold_primitive(y,p) - primEnd)
    endif

  end function rpi_x_star


  function find_rpi_x_star(y,rpiData)   
    implicit none
    real(kp) :: find_rpi_x_star
    real(kp), intent(in) :: y
    type(transfert), optional, intent(inout) :: rpiData

    real(kp) :: primStar,p,w,CalFplusprimEnd,potStar,epsOneStar

    p=rpiData%real1
    w = rpiData%real2
    CalFplusprimEnd = rpiData%real3

    primStar = rpi_efold_primitive(y,p)
    epsOneStar = rpi_epsilon_one(y,p)
    potStar = rpi_norm_potential(y,p)

    find_rpi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rpi_x_star



  function rpi_lnrhoend(p,Pstar) 
    implicit none
    real(kp) :: rpi_lnrhoend
    real(kp), intent(in) :: p,Pstar

    real(kp) :: yEnd, potEnd, epsOneEnd
    real(kp) :: y, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    yEnd = rpi_x_endinf(p)


    potEnd  = rpi_norm_potential(yEnd,p)

    epsOneEnd = rpi_epsilon_one(yEnd,p)

!   Trick to return y such that rho_reh=rho_end

    y = rpi_x_star(p,wrad,junk,Pstar)  
     
    potStar = rpi_norm_potential(y,p)
    epsOneStar = rpi_epsilon_one(y,p)

    print *,'y yend',y,yend
    print *,'eps', rpi_epsilon_one(y,p),rpi_epsilon_one(yend,p)

    PRINT*,'rpi_lnrhoend   :ystar=',y,'  potStar=',potStar,'  epsOneStar=',epsOneStar

    
    if (.not.slowroll_validity(epsOneStar)) stop 'rpi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rpi_lnrhoend = lnRhoEnd

  end function rpi_lnrhoend

  
end module rpireheat
