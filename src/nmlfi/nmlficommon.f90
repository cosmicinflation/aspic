! Common routines for Non-Minimal Large Field Inflation
!
! V(hbar) = M^4[ hbar^p / (1+hbar^2)^2 ]
!
! x = phi/Mg ~ phi/Mpl, the Einstein frame field value in unit of the
!
! Jordan Frame gravity scale Mg
!
! (dx/dhbar)^2 = [ 1 + (1+6xi) hbar^2 ] / [ xi (1+hbar^2)^2]
!
! NB: hbar = sqrt(xi) phibar/Mg where phibar is the Jordan frame field
!
! The normalization of the potential is:
!
!  M^4 = lambdabar Mg^4 / xi^2
!


module nmlficommon
  use infprec, only : kp, pi, toldp, tolkp, transfert
  use inftools, only : zbrent
  use specialinf, only : lambert
  implicit none

!makes eps1 machine precision for xi > 1
  real(kp), parameter :: hbarBig = epsilon(1._kp)**(-0.25_kp)
  real(kp), parameter :: hbarSmall = epsilon(1._kp)**(0.5_kp)
      
  private

  public hbarBig, hbarSmall
  public nmlfi_x, nmlfi_hbar, nmlfi_deriv_x, nmlfi_deriv_second_x
  public nmlfi_norm_parametric_potential, nmlfi_norm_deriv_second_parametric_potential
  public nmlfi_norm_deriv_parametric_potential, nmlfi_parametric_hbar_endinf
  public nmlfi_parametric_epsilon_one, nmlfi_parametric_epsilon_two
  public nmlfi_parametric_epsilon_three, nmlfi_parametric_efold_primitive
  public nmlfi_parametric_hbar_trajectory, nmlfi_quartic_parametric_hbar_tracjectory


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
         / ( sqrt(1._kp + (1._kp+6._kp*xi)*hbar2) + sqrt(6._kp*xi*hbar2) ) )**sqrt(6._kp) &
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
  function nmlfi_norm_parametric_potential(hbar,xi,p)
    implicit none
    real(kp) :: nmlfi_norm_parametric_potential
    real(kp), intent(in) :: hbar,xi,p

    nmlfi_norm_parametric_potential = (hbar*hbar)**(p/2._kp) / (1._kp &
         + hbar*hbar)**2


  end function nmlfi_norm_parametric_potential



!returns the first derivative of the potential W with respect to hbar
  function nmlfi_norm_deriv_parametric_potential(hbar,xi,p)
    implicit none
    real(kp) :: nmlfi_norm_deriv_parametric_potential
    real(kp), intent(in) :: hbar,xi,p

    nmlfi_norm_deriv_parametric_potential = (p + (p-4._kp)*hbar*hbar)*hbar**(p-1._kp) &
         / (1._kp + hbar*hbar)**3._kp
    
  end function nmlfi_norm_deriv_parametric_potential


!returns the 2nd derivative of the potential W with respect to hbar
  function nmlfi_norm_deriv_second_parametric_potential(hbar,xi)
    implicit none
    real(kp) :: nmlfi_norm_deriv_second_parametric_potential
    real(kp), intent(in) :: hbar,xi

    nmlfi_norm_deriv_second_parametric_potential = hbar**(p-2._kp) * ( p*(p-1._kp) &
         + 2._kp*(p*(p-5._kp)-2._kp)*hbar*hbar +(p-4._kp)*(p-5._kp)*hbar**4._kp ) &
         / (1._kp + hbar*hbar)**4._kp

  end function nmlfi_norm_deriv_second_parametric_potential


  
!eps1 Einstein-Frame in terms of hbar
  function nmlfi_parametric_epsilon_one(hbar,xi,p)
    real(kp) :: nmlfi_parametric_epsilon_one
    real(kp), intent(in) :: hbar,xi,p

    
    nmlfi_parametric_epsilon_one = 0.5_kp*xi * ( p + (p-4._kp)*hbar*hbar)**2._kp &
         / (1._kp + (1._kp + 6._kp*xi)*hbar*hbar) / (hbar*hbar)
         
    
  end function nmlfi_parametric_epsilon_one

  

!eps2 (EF) in terms of hbar
  function nmlfi_parametric_epsilon_two(hbar,xi,p)
    implicit none
    real(kp) :: nmlfi_parametric_epsilon_two
    real(kp), intent(in) :: hbar,xi,p

    nmlfi_parametric_epsilon_two =  2_kp*xi * (1._kp + hbar*hbar)* &
         * (p + hbar*hbar*(4._kp + p + 12._kp*p*xi))  &
         /(1._kp + hbar*hbar*(1._kp + 6._kp*xi))**2 / (hbar*hbar)
        
  end function nmlfi_parametric_epsilon_two


!eps3 (EF) in terms of hbar
  function nmlfi_parametric_epsilon_three(hbar,xi,p)
    implicit none
    real(kp) :: nmlfi_parametric_epsilon_three
    real(kp), intent(in) :: hbar,xi,p

    
    nmlfi_parametric_epsilon_three = (2*(hbar*hbar*(-4._kp + p) + p)*xi &
         * (p + 3._kp*hbar*hbar*(p + 6._kp*p*xi) &
         + hbar**6._kp*(1._kp + 6._kp*xi)*(4._kp + p + 12._kp*p*xi) &
         + hbar**4._kp*(4._kp + 3._kp*p + 48._kp*xi + 36._kp*p*xi*(1._kp + 4._kp*xi)))) &
         /((hbar + hbar**3*(1 + 6._kp*xi))**2._kp*(p + hbar**2._kp*(4._kp + p + 12._kp*p*xi)))

  end function nmlfi_parametric_epsilon_three



!the solution of eps1=1 in terms of hbar. This is the only real root
!defined for all xi > 0 for p < p+ and for xi<ximax for p>p+
  function nmlfi_parametric_hbar_endinf(xi,p)
    implicit none
    real(kp) :: nmlfi_parametric_hbar_endinf
    real(kp), intent(in) :: xi,p

    real(kp) :: hbarend2
    
    hbarend2 = (sqrt( (1._kp + 2._kp*p*xi)*(1._kp+6._kp*p*xi) ) + p*xi*(p-4._kp) - 1._kp) &
         / (2._kp - xi*(p*(p-8._kp) + 4._kp))

    if (hbarend2.ge.0._kp) then
       nmlfi_parametric_hbar_endinf = sqrt(hbarend2)
    else
       write(*,*)'xi= p= ',xi,p
       write(*,*)'hbarend^2 = ',hbarend2
       stop 'nmlfi_parametric_hbar_endinf: hbarend^2 < 0'
    endif

  end function nmlfi_parametric_hbar_endinf
    

  

!efold (Einstein Frame) primitive in terms of the parameter hbar
  function nmlfi_parametric_efold_primitive(hbar,xi,p)
    implicit none
    real(kp) :: nmlfi_parametric_efold_primitive
    real(kp), intent(in) :: hbar,xi,p

    if (p.eq.4._kp) then

       nmlfi_parametric_efold_primitive = (1._kp + 6._kp*xi)*(1._kp+hbar*hbar)/(8._kp*xi) &
            - 0.75._kp*log(1._kp+hbar*hbar)
       
       return
    endif
    
    nmlfi_parametric_efold_primitive = (2._kp + 3._kp*xi*p)/(4._kp*xi*(p-.4._kp)) &
         * log(p+(p-4._kp)*hbar*hbar) - 0.75._kp*log(1._kp+hbar*hbar)
     
  end function nmlfi_parametric_efold_primitive


!returns hbar at bfold=-efolds before the end of inflation, ie N-Nend
  function nmlfi_parametric_hbar_trajectory(bfold,hbarend,xi,p)
    implicit none
    real(kp), intent(in) :: bfold, hbarend,xi
    real(kp) :: nmlfi_parametric_hbar_trajectory
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: nmlfiData
    
!there is an analytical solution only for p=4 (see next function)

    if (p.eq.4._kp) then
       nmlfi_parametric_hbar_trajectory = nmlfi_quartic_parametric_hbar_trajectory(bfold,hbarend,xi)
       return
    endif
    
    mini = hbarend
    maxi = hbarBig

    nmlfiData%real1 = xi
    nmlfiData%real2 = -bfold + nmlfi_parametric_efold_primitive(hbarend,xi,p)

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


  
  
!analytical solution for p=4
  function nmlfi_quartic_parametric_hbar_trajectory(bfold,hbarend,xi)
    implicit none
    real(kp), intent(in) :: bfold, hbarend,xi
    real(kp) :: nmlfi_quartic_parametric_hbar_trajectory

    real(kp) :: xibfold, hbar2, hbarend2

    xibfold = xi*bfold
    hbarend2 = hbarend*hbarend
    
    hbar2 = -1._kp - 6._kp*xi/(1._kp + 6._kp*xi) &
         * lambert(- (1._kp + 1._kp/(6._kp*xi))*(1._kp+hbarend2) &
         * exp(-((1._kp + 6._kp*xi)*(1._kp*hbarend2) - 8._kp*xibfold )/(6._kp*xi)),-1)

    if (hbar2.lt.0._kp) then
       stop 'nmlfi_quartic_parametric_hbar_trajectory: hbar^2 < 0!'
    else    
       nmlfi_quartic_parametric_hbar_trajectory = sqrt(hbar2)
    endif
    
  end function nmlfi_quartic_parametric_hbar_trajectory



  
end module nmlficommon
