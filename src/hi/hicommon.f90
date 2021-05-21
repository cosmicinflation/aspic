! Common routines for Higgs inflation out of the large field approximation.
!
! With phi cannonically normalized field in the Einstein Frame:
!
! x = phi/Mg
!
! where Mg is the gravity scale in the Jordan Frame (Mg ~ Mpl).
!
! The Einstein frame parametric potential:
!
! W(x) = M^4 (hbar^2 - vbar^2)^2 / (1 + hbar^2)^2
!
! with hbar(x),  vbar^2 = xi*v^2/Mg^2 and  M^4 = lambda Mg^4/(4 xi^2)
!
! H is the Higgs field in the Jordan frame and hbar is defined as
!
! hbar = sqrt(xi) * H/Mg
!
! with the Higgs potential taken as
!
! V(H) = - lambda/4 ( H^2 - v^2 )^2
!
!


module hicommon
  use infprec, only : kp, pi, toldp, tolkp, transfert
  use inftools, only : zbrent
  use cosmopar, only : higgsVeV
  use specialinf, only : lambert
  implicit none

  real(kp), parameter :: sqr6 = sqrt(6._kp)

!for the record, this guy is 10^(-34), it starts being significant when
!multiplied by xi larger than 10^(34). Arf!
  
  real(kp), parameter :: vev2 = higgsVeV*higgsVeV
  
  private

  public vev2
  public hi_x, hi_hbar, hi_deriv_x, hi_deriv_second_x
  public hi_norm_parametric_potential, hi_norm_deriv_second_parametric_potential
  public hi_norm_deriv_parametric_potential, hi_parametric_hbar_endinf
  public hi_parametric_epsilon_one, hi_parametric_epsilon_two
  public hi_parametric_epsilon_three, hi_parametric_efold_primitive


contains


!returns the field value x from hbar and xi
  function hi_x(hbar,xi)
    implicit none
    real(kp) :: hi_x
    real(kp), intent(in) :: hbar, xi

    hi_x = log(hi_expx(hbar*hbar,xi))
   
  end function hi_x


!returns exp(x) as a function of hbar^2 and xi  
  function hi_expx(hbar2,xi)
    implicit none
    real(kp) :: hi_expx
    real(kp), intent(in) :: hbar2,xi
    
    hi_expx = ( sqrt(1._kp+hbar2) &
         / ( sqrt(1._kp + (1._kp+6._kp*xi)*hbar2) + sqrt(6._kp*xi*hbar2) ) )**sqr6 &
         * ( sqrt(1._kp + (1._kp + 6._kp*xi)*hbar2) &
         + sqrt((1._kp+6._kp*xi)*hbar2) )**(sqrt(6._kp+1._kp/xi))

  end function hi_expx
  
  

!returns hbar from the field value x by inverting hi_expx
  function hi_hbar(x,xi)
    implicit none
    real(kp) :: hi_hbar
    real(kp), intent(in) :: x,xi
    type(transfert) :: hiData
    real(kp), parameter :: tolFind = tolkp

    real(kp) :: hbar2
    
    real(kp) :: mini, maxi

    
    if (x.lt.0._kp) stop 'hi_hbar: x < 0!'

    mini = 0._kp
    maxi = huge(1._kp)

    hiData%real1 = exp(x)
    hiData%real2 = xi
    hbar2 = zbrent(find_hi_hbar2,mini,maxi,tolFind,hiData)

    hi_hbar = sqrt(hbar2)

  end function hi_hbar


  
  function find_hi_hbar2(hbar2,hiData)
    real(kp) :: find_hi_hbar2
    real(kp), intent(in) :: hbar2
    type(transfert), optional, intent(inout) :: hiData
    real(kp) :: expx

    expx = hiData%real1
    xi = hiData%real2
    find_hi_hbar2 = hi_expx(hbar2,xi) - expx

  end function find_hi_hbar2


!this is dx/dhbar
  function hi_deriv_x(hbar,xi)
    implicit none
    real(kp) :: hi_deriv_x
    real(kp), intent(in) :: hbar, xi

    hi_deriv_x = sqrt( (1._kp + (1._kp+6._kp*xi)*hbar*hbar ) &
         / (xi * (1._kp + hbar*hbar)**2._kp ) )
    
  end function hi_deriv_x


!this is d^x/dhbar^2  
  function hi_deriv_second_x(hbar,xi)
    implicit none
    real(kp) :: hi_deriv_second_x
    real(kp), intent(in) :: hbar, xi

    hi_deriv_second_x = - hbar*(1._kp + hbar*hbar +6._kp*xi*(hbar*hbar-1._kp)) &
         / (xi * (1._kp+hbar*hbar)**3._kp * sqrt(hi_deriv_x(hbar,xi)) )
    
  end function hi_deriv_second_x
  


!returns the normalized Einstein Frame potential W(hbar)/M^4 in terms of
!the parameter hbar
  function hi_norm_parametric_potential(hbar,xi)
    implicit none
    real(kp) :: hi_norm_parametric_potential
    real(kp), intent(in) :: hbar,xi

    hi_norm_parametric_potential = (hbar*hbar - xi*vev2)**2 / (1._kp &
         + hbar*hbar)**2


  end function hi_norm_parametric_potential



!returns the first derivative of the potential W with respect to hbar
  function hi_norm_deriv_parametric_potential(hbar,xi)
    implicit none
    real(kp) :: hi_norm_deriv_parametric_potential
    real(kp), intent(in) :: hbar,xi

    hi_norm_deriv_parametric_potential = 4._kp* (1._kp + xi*vev2) &
         * hbar*(hbar*hbar - xi*vev2) / (1._kp + hbar*hbar)**3._kp
    
  end function hi_norm_deriv_parametric_potential


!returns the 2nd derivative of the potential W with respect to hbar
  function hi_norm_deriv_second_parametric_potential(hbar,xi)
    implicit none
    real(kp) :: hi_norm_deriv_second_parametric_potential
    real(kp), intent(in) :: hbar,xi

    hi_norm_deriv_second_parametric_potential = -4._kp * (1._kp+xi*vev2) &
         * (xi*vev2 + 3._kp*hbar**4._kp - hbar*hbar*(3._kp + 5._kp*xi*vev2) ) &
         / (1._kp + hbar*hbar)**4._kp

  end function hi_norm_deriv_second_parametric_potential


  
!eps1 Einstein-Frame in terms of hbar
  function hi_parametric_epsilon_one(hbar,xi)
    real(kp) :: hi_parametric_epsilon_one
    real(kp), intent(in) :: hbar,xi
    real(kp) :: xivev2

    xivev2 = xi*vev2
    
    hi_parametric_epsilon_one = 8._kp*xi * (1._kp + xivev2)**2._kp &
         * hbar*hbar / (hbar*hbar - xivev2)**2._kp / (1._kp + (1._kp + 6._kp*xi) * hbar*hbar)
    
  end function hi_parametric_epsilon_one

  

!eps2 (EF) in terms of hbar
  function hi_parametric_epsilon_two(hbar,xi)
    implicit none
    real(kp) :: hi_parametric_epsilon_two
    real(kp), intent(in) :: hbar,xi
    real(kp) :: xivev2

    xivev2 = xi*vev2

    hi_parametric_epsilon_two =  (8._kp*(1._kp + hbar*hbar)*xi*(1._kp + xivev2) &
         *(hbar*hbar + 2._kp*hbar**4._kp*(1._kp + 6._kp*xi) + xivev2)) &
         / ((1._kp + hbar*hbar*(1._kp + 6._kp*xi))**2 * (hbar*hbar - xivev2)**2._kp)
        
  end function hi_parametric_epsilon_two


!eps3 (EF) in terms of hbar
  function hi_parametric_epsilon_three(hbar,xi)
    implicit none
    real(kp) :: hi_parametric_epsilon_three
    real(kp), intent(in) :: hbar,xi

    real(kp) :: xivev2

    xivev2 = xi*vev2
    
    hi_parametric_epsilon_three = (8._kp*hbar*hbar * xi * (1._kp + xivev2) &
         * (2._kp*hbar**8._kp*(1._kp + 6._kp*xi)**2._kp + 3._kp*xivev2 &
         - (1._kp + 12._kp*xi)*xivev2**2._kp + 2._kp*hbar**6._kp &
         * (1._kp + 6._kp*xi)**2._kp * (2._kp + xivev2) + 3._kp*hbar**4._kp &
         * (1._kp + 6._kp*xi)*(1._kp + 3._kp*xivev2) + hbar*hbar &
         * (1._kp + 2._kp*(5._kp + 21._kp*xi)*xivev2 - (1._kp + 6._kp*xi) &
         *xivev2**2._kp) ) ) &
         / ( (1._kp + hbar*hbar*(1._kp + 6._kp*xi))**2._kp &
         * (hbar*hbar - xivev2)**2._kp * ( hbar*hbar &
         + 2._kp*hbar**4._kp*(1._kp + 6._kp*xi) + xivev2 ) )

  end function hi_parametric_epsilon_three



!the solution of eps1=1 in terms of hbar (third order polynomial in hbar^2
!solution irreducible). This is the real root defined for all xi > 0
  function hi_parametric_hbar_endinf(xi)
    implicit none
    real(kp) :: hi_parametric_hbar_endinf
    real(kp), intent(in) :: xi

    real(kp), parameter :: onethird = 1._kp/3._kp

    real(kp) :: xivev2
    complex(kp) :: zend2


    xivev2 = xi*vev2
    
    zend2 = (4._kp - 8._kp*(1._kp + 6._kp*xi)*xivev2 - (cmplx(0._kp,2._kp) &
         *( cmplx(0._kp,-1._kp) + sqrt(3._kp)) * ((1._kp + 12._kp*xi)**2._kp &
         + 2._kp*(1._kp + 6._kp*xi)*(1._kp + 24._kp*xi)*xivev2 &
         + (1._kp + 6._kp*xi)*(1._kp + 30._kp*xi)*xivev2**2._kp)) &
         / (1._kp + 36._kp*xi + 216._kp*xi**2._kp + 3._kp*xivev2 + 18._kp*xi*xivev2 &
         - 432._kp*xi**2._kp*xivev2 - 2592._kp*xi**3._kp*xivev2 + 3._kp*xivev2**2._kp &
         - 72._kp*xi*xivev2**2._kp - 1404._kp*xi**2._kp*xivev2**2._kp &
         - 5184._kp*xi**3._kp*xivev2**2._kp + xivev2**3._kp - 54._kp*xi*xivev2**3._kp &
         - 756._kp*xi**2._kp*xivev2**3._kp - 2376._kp*xi**3._kp*xivev2**3._kp &
         + 6._kp*sqrt(6._kp)*(1._kp + 6._kp*xi)*(1._kp + xivev2) &
         * sqrt(xi*(-(xivev2*(1._kp + xivev2)**3._kp) - 24._kp*xi**3._kp*(4._kp &
         + 8._kp*xivev2 + xivev2**2._kp)**2._kp - 2._kp*xi*(1._kp + xivev2)**2._kp &
         * (1._kp + 20._kp*xivev2 + xivev2**2._kp) + 4._kp*xi**2._kp*(1._kp + xivev2) &
         * (-16._kp - 108._kp*xivev2 - 60._kp*xivev2**2._kp + 5._kp*xivev2**3._kp))))**onethird &
         + cmplx(0._kp,2._kp)*(cmplx(0._kp,1._kp) + sqrt(3._kp)) &
         * (1._kp + 36._kp*xi + 216._kp*xi**2._kp + 3._kp*xivev2 + 18._kp*xi*xivev2 &
         - 432._kp*xi**2._kp*xivev2 - 2592._kp*xi**3._kp*xivev2 + 3._kp*xivev2**2._kp &
         - 72._kp*xi*xivev2**2._kp - 1404._kp*xi**2._kp*xivev2**2._kp &
         - 5184._kp*xi**3._kp*xivev2**2._kp + xivev2**3._kp - 54._kp*xi*xivev2**3._kp &
         - 756._kp*xi**2._kp*xivev2**3._kp - 2376._kp*xi**3._kp*xivev2**3._kp &
         + 6._kp*sqrt(6._kp)*(1._kp + 6._kp*xi)*(1._kp + xivev2) &
         * sqrt(xi*(-(xivev2*(1._kp + xivev2)**3) - 24._kp*xi**3._kp &
         * (4._kp + 8._kp*xivev2 + xivev2**2._kp)**2._kp - 2._kp*xi*(1._kp + xivev2)**2._kp &
         * (1._kp + 20._kp*xivev2 + xivev2**2._kp) + 4._kp*xi**2._kp*(1._kp + xivev2) &
         * (-16._kp - 108._kp*xivev2 - 60._kp*xivev2**2._kp + 5._kp*xivev2**3._kp))))**onethird) &
         / (12._kp*(-1._kp - 6._kp*xi))

    hi_parametric_hbar_endinf = sqrt(real(zend2,kp))
    

  end function hi_parametric_hbar_endinf
    

  

!efold (Einstein Frame) primitive in terms of the parameter hbar
  function hi_parametric_efold_primitive(hbar,xi)
    implicit none
    real(kp) :: hi_parametric_efold_primitive
    real(kp), intent(in) :: hbar,xi
    type(transfert) :: hiData

    real(kp) :: xivev2

    xivev2 = xi*vev2
    
    hi_parametric_efold_primitive = 1._kp/(8._kp*xi*(1._kp+xivev2)) &
         * ( (1._kp + 6._kp*xi)*hbar*hbar - xivev2*log(hbar*hbar) &
         - 6._kp*xi*(1._kp + xivev2) * log(1._kp + hbar*hbar) )
     
  end function hi_parametric_efold_primitive


!returns hbar at bfold=-efolds before the end of inflation, ie N-Nend
  function hi_parametric_hbar_trajectory(bfold,hbarend,xi)
    implicit none
    real(kp), intent(in) :: bfold, hbarend,xi
    real(kp) :: hi_parametric_hbar_trajectory
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: rcipiData
    real(kp) :: mini, maxi, xivev2, hbar2
    
!there is an analytical solution (inverting
!hi_parametric_efold_primitive) but only for xivev2 = 0, let us do
!this for xivev2 << 1

    xivev2 = xi*vev2
    
    if (xivev2 .le. epsilon(1._kp)) then
       hi_parametric_hbar_trajectory = hi_approx_parametric_hbar_trajectory(bfold,hbarend,xi)
       return
    endif

    mini = hbarend2
    maxi = huge(1._kp)

    hiData%real1 = xi
    hiData%real2 = -bfold + hi_parametric_efold_primitive(hbarend,xi)

    hbar = zbrent(find_hi_parametric_hbar_trajectory,mini,maxi,tolFind,hiData)

  end function hi_parametric_hbar_trajectory
  
  function find_hi_parametric_hbar_trajectory(hbar,rcipiData)
    implicit none
    real(kp) :: find_hi_parametric_hbar_trajectory
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: hiData

    real(kp) :: xi,NplusNuend

    xi = rcipiData%real1
    NplusNuend = rcipiData%real2

    find_hi_parametric_hbar_trajectory = hi_parametric_efold_primitive(hbar,xi) - NplusNuend

  end function find_hi_parametric_hbar_trajectory


  
  
!analytical solution obtained for xivev2=0
  function hi_approx_parametric_hbar_trajectory(bfold,hbarend,xi)
    implicit none
    real(kp), intent(in) :: bfold, hbarend,xi
    real(kp) :: hi_approx_parametric_hbar_trajectory

    real(kp) :: xibfold, hbar2

    xibfold = xi*bfold
    hbarend2 = hbarend*hbarend
    
    hbar2 = -1._kp - 6._kp*xi/(1._kp + 6._kp*xi) &
         * lambert(- (1._kp + 1._kp/(6._kp*xi))*(1._kp+hbarend2) &
         * exp(-((1._kp + 6._kp*xi)*(1._kp*hbarend2) - 8._kp*xibfold )/(6._kp*xi)),-1)

    if (hbar2.lt.0._kp) then
       stop 'hi_approx_parametric_hbar_trajectory: hbar^2 < 0!'
    else    
       hi_approx_parametric_hbar_trajectory = sqrt(hbar2)
    endif
    
  end function hi_approx_parametric_hbar_trajectory



  
end module hicommon
