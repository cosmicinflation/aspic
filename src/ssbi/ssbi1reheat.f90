!spontaneous symmetry breaking 1 reheating functions in the slow-roll approximations

module ssbi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssbicomreh, only : ssbi_x_star, ssbi_x_rrad, ssbi_x_rreh
  use ssbi1sr, only : ssbi1_epsilon_one, ssbi1_epsilon_two, ssbi1_epsilon_three
  use ssbi1sr, only : ssbi1_norm_potential
  use ssbi1sr, only : ssbi1_x_endinf, ssbi1_efold_primitive

  implicit none

  private

  public ssbi1_x_star, ssbi1_lnrhoreh_max
  public ssbi1_x_rrad, ssbi1_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi1_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: ssbi1_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi
    real(kp) :: xend
   
    xEnd=ssbi1_x_endinf(alpha,beta)
    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)

    ssbi1_x_star = ssbi_x_star(alpha,beta,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function ssbi1_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ssbi1_x_rrad(alpha,beta,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssbi1_x_rrad
    real(kp), intent(in) :: alpha,beta,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi
    real(kp) :: xend

    xEnd=ssbi1_x_endinf(alpha,beta)
    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)
    
    ssbi1_x_rrad = ssbi_x_rrad(alpha,beta,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function ssbi1_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function ssbi1_x_rreh(alpha,beta,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ssbi1_x_rreh
    real(kp), intent(in) :: alpha,beta,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi
    real(kp) :: xend

    xEnd=ssbi1_x_endinf(alpha,beta)
    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)
    
    ssbi1_x_rreh = ssbi_x_rreh(alpha,beta,lnRreh,xend,mini,maxi,bfoldstar)

  end function ssbi1_x_rreh



  function ssbi1_lnrhoreh_max(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssbi1_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssbi1_x_endinf(alpha,beta)

    potEnd  = ssbi1_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssbi1_epsilon_one(xEnd,alpha,beta)

!    print*,'ssbi1_lnrhoreh_max:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssbi1_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssbi1_norm_potential(x,alpha,beta)
    epsOneStar = ssbi1_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssbi1_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssbi1_lnrhoreh_max = lnRhoEnd

  end function ssbi1_lnrhoreh_max

  
end module ssbi1reheat
