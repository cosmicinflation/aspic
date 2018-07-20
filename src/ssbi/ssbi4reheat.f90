!spontaneous symmetry breaking 4 reheating functions in the slow-roll approximations

module ssbi4reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : ln_rho_endinf, ln_rho_reheat
  use ssbicomreh, only : ssbi_x_star, ssbi_x_rrad, ssbi_x_rreh
  use ssbi4sr, only : ssbi4_epsilon_one, ssbi4_epsilon_two, ssbi4_epsilon_three
  use ssbi4sr, only : ssbi4_norm_potential, ssbi4_x_potmax
  use ssbi4sr, only : ssbi4_x_endinf, ssbi4_efold_primitive
  implicit none

  private

  public ssbi4_x_star, ssbi4_lnrhoreh_max
  public ssbi4_x_rrad, ssbi4_x_rreh

contains



!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi4_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: ssbi4_x_star
    real(kp), intent(in) :: alpha,beta,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    mini = ssbi4_x_potmax(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssbi4_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))

    ssbi4_x_star = ssbi_x_star(alpha,beta,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function ssbi4_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ssbi4_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssbi4_x_rrad
    real(kp), intent(in) :: alpha,beta,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = ssbi4_x_potmax(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssbi4_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))
    
    ssbi4_x_rrad = ssbi_x_rrad(alpha,beta,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function ssbi4_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function ssbi4_x_rreh(alpha,beta,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ssbi4_x_rreh
    real(kp), intent(in) :: alpha,beta,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = ssbi4_x_potmax(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssbi4_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))
    
    ssbi4_x_rreh = ssbi_x_rreh(alpha,beta,lnRreh,xend,mini,maxi,bfoldstar)

  end function ssbi4_x_rreh


 

  function ssbi4_lnrhoreh_max(alpha,beta,xend,Pstar) 
    implicit none
    real(kp) :: ssbi4_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = ssbi4_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssbi4_epsilon_one(xEnd,alpha,beta)


!   Trick to return x such that rho_reh=rho_end

    x = ssbi4_x_star(alpha,beta,xend,wrad,junk,Pstar)  


    potStar = ssbi4_norm_potential(x,alpha,beta)
    epsOneStar = ssbi4_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssbi4_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssbi4_lnrhoreh_max = lnRhoEnd

  end function ssbi4_lnrhoreh_max

  
end module ssbi4reheat
