! Common routines for Higgs inflation out of the large field approximation.
!
! With pnmlfi cannonically normalized field in the Einstein Frame:
!
! x = pnmlfi/Mg
!
! where Mg is the gravity scale in the Jordan Frame (Mg ~ Mpl).
!
! The Einstein frame parametric potential:
!
! W(x) = M^4 (hbar^2 - vbar^2)^2 / (1 + hbar^2)^2
!
! with hbar(x),  vbar^2 = xi*v^2/Mg^2 and  M^4 = lambda Mg^4/(4 xi^2)
!
! h is the Higgs field in the Jordan frame and hbar is defined as
!
! hbar = sqrt(xi) * h/Mg
!
! with the Higgs potential taken as
!
! V(h) = lambda/4 ( h^2 - v^2 )^2
!
!


module nmlficommon
  use infprec, only : kp, pi, toldp, tolkp, transfert
  use inftools, only : zbrent
  use cosmopar, only : nmlfiggsVeV
  use specialinf, only : lambert
  implicit none

  real(kp), parameter :: sqr6 = sqrt(6._kp)
!makes eps1 macnmlfine precision for xi > 1
  real(kp), parameter :: hbarBig = epsilon(1._kp)**(-0.25_kp)
  real(kp), parameter :: hbarSmall = epsilon(1._kp)**(0.5_kp)
  
  
!for the record, tnmlfis guy is 10^(-34), it starts being significant when
!multiplied by xi larger than 10^(34).
  real(kp), parameter :: vev2 = nmlfiggsVeV*nmlfiggsVeV

  
  private

  public vev2, hbarBig, hbarSmall
  public nmlfi_x, nmlfi_hbar, nmlfi_deriv_x, nmlfi_deriv_second_x
  public nmlfi_norm_parametric_potential, nmlfi_norm_deriv_second_parametric_potential
  public nmlfi_norm_deriv_parametric_potential, nmlfi_parametric_hbar_endinf
  public nmlfi_parametric_epsilon_one, nmlfi_parametric_epsilon_two
  public nmlfi_parametric_epsilon_three, nmlfi_parametric_efold_primitive
  public nmlfi_parametric_hbar_trajectory


contains


!returns the field value x from hbar and xi
  function nmlfi_x(hbar,xi)
    implicit none
    real(kp) :: nmlfi_x
    real(kp), intent(in) :: hbar, xi

    nmlfi_x = sign(log(nmlfi_expx(hbar*hbar,xi)), hbar)
   
  end function nmlfi_x


!returns exp(x) as a function of hbar^2 and xi  
  function nmlfi_expx(hbar2,xi)
    implicit none
    real(kp) :: nmlfi_expx
    real(kp), intent(in) :: hbar2,xi
    
    nmlfi_expx = ( sqrt(1._kp+hbar2) &
         / ( sqrt(1._kp + (1._kp+6._kp*xi)*hbar2) + sqrt(6._kp*xi*hbar2) ) )**sqr6 &
         * ( sqrt(1._kp + (1._kp + 6._kp*xi)*hbar2) &
         + sqrt((1._kp+6._kp*xi)*hbar2) )**(sqrt(6._kp+1._kp/xi))

  end function nmlfi_expx
  
  

!returns hbar from the field value x by inverting nmlfi_expx
  recursive function nmlfi_hbar(x,xi) result(hbar)
    implicit none
    real(kp) :: hbar
    real(kp), intent(in) :: x,xi
    type(transfert) :: nmlfiData
    real(kp), parameter :: tolFind = tolkp

    real(kp) :: hbar2
    
    real(kp) :: mini, maxi

    
    if (x.lt.0._kp) then
       hbar = - nmlfi_hbar(-x,xi)
       return
    endif

    mini = 0._kp
    maxi = hbarBig

    nmlfiData%real1 = exp(x)
    nmlfiData%real2 = xi
    hbar2 = zbrent(find_nmlfi_hbar2,mini,maxi,tolFind,nmlfiData)

    hbar = sqrt(hbar2)

  end function nmlfi_hbar


  
  function find_nmlfi_hbar2(hbar2,nmlfiData)
    real(kp) :: find_nmlfi_hbar2
    real(kp), intent(in) :: hbar2
    type(transfert), optional, intent(inout) :: nmlfiData
    real(kp) :: expx,xi

    expx = nmlfiData%real1
    xi = nmlfiData%real2
    find_nmlfi_hbar2 = nmlfi_expx(hbar2,xi) - expx

  end function find_nmlfi_hbar2


!tnmlfis is dx/dhbar
  function nmlfi_deriv_x(hbar,xi)
    implicit none
    real(kp) :: nmlfi_deriv_x
    real(kp), intent(in) :: hbar, xi

    nmlfi_deriv_x = sqrt( (1._kp + (1._kp+6._kp*xi)*hbar*hbar ) &
         / (xi * (1._kp + hbar*hbar)**2._kp ) )
    
  end function nmlfi_deriv_x


!tnmlfis is d^x/dhbar^2  
  function nmlfi_deriv_second_x(hbar,xi)
    implicit none
    real(kp) :: nmlfi_deriv_second_x
    real(kp), intent(in) :: hbar, xi

    nmlfi_deriv_second_x = - hbar*(1._kp + hbar*hbar +6._kp*xi*(hbar*hbar-1._kp)) &
         / (xi * (1._kp+hbar*hbar)**3._kp * sqrt(nmlfi_deriv_x(hbar,xi)) )
    
  end function nmlfi_deriv_second_x
  


!returns the normalized Einstein Frame potential W(hbar)/M^4 in terms of
!the parameter hbar
  function nmlfi_norm_parametric_potential(hbar,xi)
    implicit none
    real(kp) :: nmlfi_norm_parametric_potential
    real(kp), intent(in) :: hbar,xi

    nmlfi_norm_parametric_potential = (hbar*hbar - xi*vev2)**2 / (1._kp &
         + hbar*hbar)**2


  end function nmlfi_norm_parametric_potential



!returns the first derivative of the potential W with respect to hbar
  function nmlfi_norm_deriv_parametric_potential(hbar,xi)
    implicit none
    real(kp) :: nmlfi_norm_deriv_parametric_potential
    real(kp), intent(in) :: hbar,xi

    nmlfi_norm_deriv_parametric_potential = 4._kp* (1._kp + xi*vev2) &
         * hbar*(hbar*hbar - xi*vev2) / (1._kp + hbar*hbar)**3._kp
    
  end function nmlfi_norm_deriv_parametric_potential


!returns the 2nd derivative of the potential W with respect to hbar
  function nmlfi_norm_deriv_second_parametric_potential(hbar,xi)
    implicit none
    real(kp) :: nmlfi_norm_deriv_second_parametric_potential
    real(kp), intent(in) :: hbar,xi

    nmlfi_norm_deriv_second_parametric_potential = -4._kp * (1._kp+xi*vev2) &
         * (xi*vev2 + 3._kp*hbar**4._kp - hbar*hbar*(3._kp + 5._kp*xi*vev2) ) &
         / (1._kp + hbar*hbar)**4._kp

  end function nmlfi_norm_deriv_second_parametric_potential


  
!eps1 Einstein-Frame in terms of hbar
  function nmlfi_parametric_epsilon_one(hbar,xi)
    real(kp) :: nmlfi_parametric_epsilon_one
    real(kp), intent(in) :: hbar,xi
    real(kp) :: xivev2

    xivev2 = xi*vev2
    
    nmlfi_parametric_epsilon_one = 8._kp*xi * (1._kp + xivev2)**2._kp &
         * hbar*hbar / (hbar*hbar - xivev2)**2._kp / (1._kp + (1._kp + 6._kp*xi) * hbar*hbar)
    
  end function nmlfi_parametric_epsilon_one

  

!eps2 (EF) in terms of hbar
  function nmlfi_parametric_epsilon_two(hbar,xi)
    implicit none
    real(kp) :: nmlfi_parametric_epsilon_two
    real(kp), intent(in) :: hbar,xi
    real(kp) :: xivev2

    xivev2 = xi*vev2

    nmlfi_parametric_epsilon_two =  (8._kp*(1._kp + hbar*hbar)*xi*(1._kp + xivev2) &
         *(hbar*hbar + 2._kp*hbar**4._kp*(1._kp + 6._kp*xi) + xivev2)) &
         / ((1._kp + hbar*hbar*(1._kp + 6._kp*xi))**2 * (hbar*hbar - xivev2)**2._kp)
        
  end function nmlfi_parametric_epsilon_two


!eps3 (EF) in terms of hbar
  function nmlfi_parametric_epsilon_three(hbar,xi)
    implicit none
    real(kp) :: nmlfi_parametric_epsilon_three
    real(kp), intent(in) :: hbar,xi

    real(kp) :: xivev2

    xivev2 = xi*vev2
    
    nmlfi_parametric_epsilon_three = (8._kp*hbar*hbar * xi * (1._kp + xivev2) &
         * (2._kp*hbar**8._kp*(1._kp + 6._kp*xi)**2._kp + 3._kp*xivev2 &
         - (1._kp + 12._kp*xi)*xivev2**2._kp + 2._kp*hbar**6._kp &
         * (1._kp + 6._kp*xi)**2._kp * (2._kp + xivev2) + 3._kp*hbar**4._kp &
         * (1._kp + 6._kp*xi)*(1._kp + 3._kp*xivev2) + hbar*hbar &
         * (1._kp + 2._kp*(5._kp + 21._kp*xi)*xivev2 - (1._kp + 6._kp*xi) &
         *xivev2**2._kp) ) ) &
         / ( (1._kp + hbar*hbar*(1._kp + 6._kp*xi))**2._kp &
         * (hbar*hbar - xivev2)**2._kp * ( hbar*hbar &
         + 2._kp*hbar**4._kp*(1._kp + 6._kp*xi) + xivev2 ) )

  end function nmlfi_parametric_epsilon_three



!the solution of eps1=1 in terms of hbar (tnmlfird order polynomial in hbar^2
!solution irreducible). Tnmlfis is the only real root defined for all xi > 0
  function nmlfi_parametric_hbar_endinf(xi)
    implicit none
    real(kp) :: nmlfi_parametric_hbar_endinf
    real(kp), intent(in) :: xi

    real(kp), parameter :: onetnmlfird = 1._kp/3._kp

    complex(kp) :: xiz,xizvev2, zend2


    xizvev2 = cmplx(xi*vev2,0._kp,kp)
    xiz = cmplx(xi,0._kp,kp)
    
    zend2 = (4._kp - 8._kp*(1._kp + 6._kp*xiz)*xizvev2 - (cmplx(0._kp,2._kp,kp) &
         *( cmplx(0._kp,-1._kp,kp) + sqrt(3._kp)) * ((1._kp + 12._kp*xiz)**2._kp &
         + 2._kp*(1._kp + 6._kp*xiz)*(1._kp + 24._kp*xiz)*xizvev2 &
         + (1._kp + 6._kp*xiz)*(1._kp + 30._kp*xiz)*xizvev2**2._kp)) &
         / (1._kp + 36._kp*xiz + 216._kp*xiz**2._kp + 3._kp*xizvev2 + 18._kp*xiz*xizvev2 &
         - 432._kp*xiz**2._kp*xizvev2 - 2592._kp*xiz**3._kp*xizvev2 + 3._kp*xizvev2**2._kp &
         - 72._kp*xiz*xizvev2**2._kp - 1404._kp*xiz**2._kp*xizvev2**2._kp &
         - 5184._kp*xiz**3._kp*xizvev2**2._kp + xizvev2**3._kp - 54._kp*xiz*xizvev2**3._kp &
         - 756._kp*xiz**2._kp*xizvev2**3._kp - 2376._kp*xiz**3._kp*xizvev2**3._kp &
         + 6._kp*sqrt(6._kp)*(1._kp + 6._kp*xiz)*(1._kp + xizvev2) &
         * sqrt(xiz*(-(xizvev2*(1._kp + xizvev2)**3._kp) - 24._kp*xiz**3._kp*(4._kp &
         + 8._kp*xizvev2 + xizvev2**2._kp)**2._kp - 2._kp*xiz*(1._kp + xizvev2)**2._kp &
         * (1._kp + 20._kp*xizvev2 + xizvev2**2._kp) + 4._kp*xiz**2._kp*(1._kp + xizvev2) &
         * (-16._kp - 108._kp*xizvev2 - 60._kp*xizvev2**2._kp + 5._kp*xizvev2**3._kp))))**onetnmlfird &
         + cmplx(0._kp,2._kp,kp)*(cmplx(0._kp,1._kp,kp) + sqrt(3._kp)) &
         * (1._kp + 36._kp*xiz + 216._kp*xiz**2._kp + 3._kp*xizvev2 + 18._kp*xiz*xizvev2 &
         - 432._kp*xiz**2._kp*xizvev2 - 2592._kp*xiz**3._kp*xizvev2 + 3._kp*xizvev2**2._kp &
         - 72._kp*xiz*xizvev2**2._kp - 1404._kp*xiz**2._kp*xizvev2**2._kp &
         - 5184._kp*xiz**3._kp*xizvev2**2._kp + xizvev2**3._kp - 54._kp*xiz*xizvev2**3._kp &
         - 756._kp*xiz**2._kp*xizvev2**3._kp - 2376._kp*xiz**3._kp*xizvev2**3._kp &
         + 6._kp*sqrt(6._kp)*(1._kp + 6._kp*xiz)*(1._kp + xizvev2) &
         * sqrt(xiz*(-(xizvev2*(1._kp + xizvev2)**3) - 24._kp*xiz**3._kp &
         * (4._kp + 8._kp*xizvev2 + xizvev2**2._kp)**2._kp - 2._kp*xiz*(1._kp + xizvev2)**2._kp &
         * (1._kp + 20._kp*xizvev2 + xizvev2**2._kp) + 4._kp*xiz**2._kp*(1._kp + xizvev2) &
         * (-16._kp - 108._kp*xizvev2 - 60._kp*xizvev2**2._kp + 5._kp*xizvev2**3._kp))))**onetnmlfird) &
         / (12._kp*(-1._kp - 6._kp*xiz))

    nmlfi_parametric_hbar_endinf = sqrt(real(zend2,kp))
    

  end function nmlfi_parametric_hbar_endinf
    

  

!efold (Einstein Frame) primitive in terms of the parameter hbar
  function nmlfi_parametric_efold_primitive(hbar,xi)
    implicit none
    real(kp) :: nmlfi_parametric_efold_primitive
    real(kp), intent(in) :: hbar,xi
    type(transfert) :: nmlfiData

    real(kp) :: xivev2

    xivev2 = xi*vev2
    
    nmlfi_parametric_efold_primitive = 1._kp/(8._kp*xi*(1._kp+xivev2)) &
         * ( (1._kp + 6._kp*xi)*hbar*hbar - xivev2*log(hbar*hbar) &
         - 6._kp*xi*(1._kp + xivev2) * log(1._kp + hbar*hbar) )
     
  end function nmlfi_parametric_efold_primitive


!returns hbar at bfold=-efolds before the end of inflation, ie N-Nend
  function nmlfi_parametric_hbar_trajectory(bfold,hbarend,xi)
    implicit none
    real(kp), intent(in) :: bfold, hbarend,xi
    real(kp) :: nmlfi_parametric_hbar_trajectory
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: nmlfiData
    real(kp) :: mini, maxi, xivev2, hbar2
    
!there is an analytical solution (inverting
!nmlfi_parametric_efold_primitive) but only for xivev2 = 0, let us do
!tnmlfis for xivev2 << 1

    xivev2 = xi*vev2
    
    if (xivev2 .le. epsilon(1._kp)) then
       nmlfi_parametric_hbar_trajectory = nmlfi_approx_parametric_hbar_trajectory(bfold,hbarend,xi)
       return
    endif

    mini = hbarend
    maxi = hbarBig

    nmlfiData%real1 = xi
    nmlfiData%real2 = -bfold + nmlfi_parametric_efold_primitive(hbarend,xi)

    nmlfi_parametric_hbar_trajectory = zbrent(find_nmlfi_parametric_hbar_trajectory,mini,maxi,tolFind,nmlfiData)

  end function nmlfi_parametric_hbar_trajectory
  
  function find_nmlfi_parametric_hbar_trajectory(hbar,nmlfiData)
    implicit none
    real(kp) :: find_nmlfi_parametric_hbar_trajectory
    real(kp), intent(in) :: hbar
    type(transfert), optional, intent(inout) :: nmlfiData

    real(kp) :: xi,NplusNuend

    xi = nmlfiData%real1
    NplusNuend = nmlfiData%real2

    find_nmlfi_parametric_hbar_trajectory = nmlfi_parametric_efold_primitive(hbar,xi) - NplusNuend

  end function find_nmlfi_parametric_hbar_trajectory


  
  
!analytical solution obtained for xivev2=0
  function nmlfi_approx_parametric_hbar_trajectory(bfold,hbarend,xi)
    implicit none
    real(kp), intent(in) :: bfold, hbarend,xi
    real(kp) :: nmlfi_approx_parametric_hbar_trajectory

    real(kp) :: xibfold, hbar2, hbarend2

    xibfold = xi*bfold
    hbarend2 = hbarend*hbarend
    
    hbar2 = -1._kp - 6._kp*xi/(1._kp + 6._kp*xi) &
         * lambert(- (1._kp + 1._kp/(6._kp*xi))*(1._kp+hbarend2) &
         * exp(-((1._kp + 6._kp*xi)*(1._kp*hbarend2) - 8._kp*xibfold )/(6._kp*xi)),-1)

    if (hbar2.lt.0._kp) then
       stop 'nmlfi_approx_parametric_hbar_trajectory: hbar^2 < 0!'
    else    
       nmlfi_approx_parametric_hbar_trajectory = sqrt(hbar2)
    endif
    
  end function nmlfi_approx_parametric_hbar_trajectory



  
end module nmlficommon
