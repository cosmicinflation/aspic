!brane SUSY breaking model reheating functions in the slow-roll approximations

module bsusybireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use bsusybisr, only : bsusybi_epsilon_one, bsusybi_epsilon_two, bsusybi_epsilon_three
  use bsusybisr, only : bsusybi_norm_potential,bsusybi_efold_primitive
  implicit none

  private

  public bsusybi_x_star, bsusybi_lnrhoend,find_bsusybi_x_star

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function bsusybi_x_star(gammaBSUSYB,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: bsusybi_x_star
    real(kp), intent(in) :: gammaBSUSYB,xEnd,w,lnRhoReh,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: bsusybiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = bsusybi_epsilon_one(xEnd,gammaBSUSYB)
    potEnd = bsusybi_norm_potential(xEnd,gammaBSUSYB)
    primEnd = bsusybi_efold_primitive(xEnd,gammaBSUSYB)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    bsusybiData%real1 = gammaBSUSYB
    bsusybiData%real2 = w
    bsusybiData%real3 = calF + primEnd

    mini = xend
    maxi = 1._kp/(sqrt(6._kp)*(gammaBSUSYB-1._kp))* &
           log((sqrt(3._kp)-1._kp)/(1-gammaBSUSYB*sqrt(3._kp)))    !Value of x such that epsilon1=1

    x = zbrent(find_bsusybi_x_star,mini,maxi,tolFind,bsusybiData)
    bsusybi_x_star = x

    if (present(bfold)) then
       bfold = -(bsusybi_efold_primitive(x,gammaBSUSYB) - primEnd)
    endif

  end function bsusybi_x_star

  function find_bsusybi_x_star(x,bsusybiData)   
    implicit none
    real(kp) :: find_bsusybi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: bsusybiData

    real(kp) :: primStar,gammaBSUSYB,xend,w,CalFplusPrimEnd,potStar,epsOneStar

    gammaBSUSYB=bsusybiData%real1
    w = bsusybiData%real2
    CalFplusPrimEnd = bsusybiData%real3

    primStar = bsusybi_efold_primitive(x,gammaBSUSYB)
    epsOneStar = bsusybi_epsilon_one(x,gammaBSUSYB)
    potStar = bsusybi_norm_potential(x,gammaBSUSYB)


    find_bsusybi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)


  end function find_bsusybi_x_star


  function bsusybi_lnrhoend(gammaBSUSYB,xEnd,Pstar) 
    implicit none
    real(kp) :: bsusybi_lnrhoend
    real(kp), intent(in) :: gammaBSUSYB,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
         
    potEnd  = bsusybi_norm_potential(xEnd,gammaBSUSYB)
    epsOneEnd = bsusybi_epsilon_one(xEnd,gammaBSUSYB)
       
    x = bsusybi_x_star(gammaBSUSYB,xEnd,wrad,junk,Pstar)    
    potStar = bsusybi_norm_potential(x,gammaBSUSYB)
    epsOneStar = bsusybi_epsilon_one(x,gammaBSUSYB)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'bsusybi_lnrhoend: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    bsusybi_lnrhoend = lnRhoEnd

  end function bsusybi_lnrhoend

  
  
end module bsusybireheat
