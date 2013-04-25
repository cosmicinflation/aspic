!double well inflation reheating functions in the slow-roll approximations

module dwireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use dwisr, only : dwi_epsilon_one, dwi_epsilon_two, dwi_epsilon_three
  use dwisr, only : dwi_norm_potential
  use dwisr, only : dwi_x_endinf, dwi_efold_primitive
  implicit none

  private

  public dwi_x_star, dwi_lnrhoreh_max 

contains

!returns x =phi/phi0 such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function dwi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: dwi_x_star
    real(kp), intent(in) :: phi0,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: dwiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = dwi_x_endinf(phi0)

    epsOneEnd = dwi_epsilon_one(xEnd,phi0)
    potEnd = dwi_norm_potential(xEnd)

    primEnd = dwi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    dwiData%real1 = phi0
    dwiData%real2 = w
    dwiData%real3 = calF + primEnd

    mini = epsilon(1._kp)

    maxi = xEnd

   
    x = zbrent(find_dwi_x_star,mini,maxi,tolzbrent,dwiData)
    dwi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (dwi_efold_primitive(x,phi0) - primEnd)
    endif

  end function dwi_x_star


  function find_dwi_x_star(x,dwiData)   
    implicit none
    real(kp) :: find_dwi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: dwiData

    real(kp) :: primStar,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=dwiData%real1
    w = dwiData%real2
    CalFplusprimEnd = dwiData%real3

    primStar = dwi_efold_primitive(x,phi0)
    epsOneStar = dwi_epsilon_one(x,phi0)
    potStar = dwi_norm_potential(x)

    find_dwi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_dwi_x_star



  function dwi_lnrhoreh_max(phi0,Pstar) 
    implicit none
    real(kp) :: dwi_lnrhoreh_max
    real(kp), intent(in) :: phi0,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = dwi_x_endinf(phi0)


    potEnd  = dwi_norm_potential(xEnd)

    epsOneEnd = dwi_epsilon_one(xEnd,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = dwi_x_star(phi0,wrad,junk,Pstar)  

 
    potStar = dwi_norm_potential(x)
    epsOneStar = dwi_epsilon_one(x,phi0)

   ! PRINT*,'dwi_lnrhoreh_max   :xstar=',x,'  potStar=',potStar,'  epsOneStar=',epsOneStar

    
    if (.not.slowroll_validity(epsOneStar)) stop 'dwi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    dwi_lnrhoreh_max = lnRhoEnd

  end function dwi_lnrhoreh_max

  
end module dwireheat
