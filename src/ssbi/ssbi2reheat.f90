!spontaneous symmetry breaking 2 reheating functions in the slow-roll approximations

module ssbi2reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : ln_rho_endinf, ln_rho_reheat
  use ssbicomreh, only : ssbi_x_star, ssbi_x_rrad, ssbi_x_rreh
  use ssbi2sr, only : ssbi2_epsilon_one, ssbi2_epsilon_two, ssbi2_epsilon_three
  use ssbi2sr, only : ssbi2_norm_potential
  use ssbi2sr, only : ssbi2_x_endinf, ssbi2_efold_primitive
  implicit none

  private

  public ssbi2_x_star, ssbi2_lnrhoreh_max
  public ssbi2_x_rrad, ssbi2_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi2_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: ssbi2_x_star
    real(kp), intent(in) :: alpha,beta,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi
   
    mini = epsilon(1._kp)
    maxi = ssbi2_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))

    ssbi2_x_star = ssbi_x_star(alpha,beta,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function ssbi2_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ssbi2_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssbi2_x_rrad
    real(kp), intent(in) :: alpha,beta,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = epsilon(1._kp)
    maxi = ssbi2_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))

    ssbi2_x_rrad = ssbi_x_rrad(alpha,beta,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function ssbi2_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function ssbi2_x_rreh(alpha,beta,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ssbi2_x_rreh
    real(kp), intent(in) :: alpha,beta,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = epsilon(1._kp)
    maxi = ssbi2_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))
        
    ssbi2_x_rreh = ssbi_x_rreh(alpha,beta,lnRreh,xend,mini,maxi,bfoldstar)

  end function ssbi2_x_rreh




  function ssbi2_lnrhoreh_max(alpha,beta,xend,Pstar) 
    implicit none
    real(kp) :: ssbi2_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = ssbi2_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssbi2_epsilon_one(xEnd,alpha,beta)

!    print*,'ssbi2_lnrhoreh_max:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssbi2_x_star(alpha,beta,xend,wrad,junk,Pstar)  


    potStar = ssbi2_norm_potential(x,alpha,beta)
    epsOneStar = ssbi2_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssbi2_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssbi2_lnrhoreh_max = lnRhoEnd

  end function ssbi2_lnrhoreh_max

  
end module ssbi2reheat
