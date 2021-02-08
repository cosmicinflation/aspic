!Common functions for the hybrid natural inflation 2 potential
!
!V(phi) = M**4 [ 1 + alpha cos(x) )
!
!x = phi/f
!
!For alpha < alphamax, inflation does not gracefully end and it stops
!at the extra parameter xend. If alpha > alphamax, then, xend has
!potential observable effects only for xend < xepsone

module hni2sr
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent

  use hnicommon, only : hni_check_params, hni_norm_potential
  use hnicommon, only : hni_norm_deriv_potential, hni_norm_deriv_second_potential
  use hnicommon, only : hni_epsilon_one, hni_epsilon_two, hni_epsilon_three
  use hnicommon, only : hni_x_epsoneunity, hni_alpha, hni_f, hni_x_potmin
  use hnicommon, only : find_hni_x_trajectory, hni_efold_primitive, hni_numacc_x_epsonenull
  use hnicommon, only : hni_numacc_xinimin
  
  implicit none

  private

  public hni2_norm_potential, hni2_norm_deriv_potential, hni2_norm_deriv_second_potential
  public hni2_epsilon_one, hni2_epsilon_two, hni2_epsilon_three
  public hni2_x_trajectory, hni2_check_params
  public hni2_xendmax, hni2_numacc_xendmin, hni2_numacc_xendmax
  public hni2_alphamax, hni2_fmin, hni2_numacc_efoldmax



contains
      

  
  function hni2_check_params(alpha,f)
    implicit none
    logical :: hni2_check_params
    real(kp), intent(in) :: alpha,f

    hni2_check_params = hni_check_params(alpha,f) .and. (alpha.lt.hni2_alphamax(f))

  end function hni2_check_params

  
  function hni2_alphamax(f)
    implicit none
    real(kp) :: hni2_alphamax
    real(kp), intent(in) :: f
    
    hni2_alphamax = hni_alpha(f)
    
  end function hni2_alphamax

  
  function hni2_fmin(alpha)
    implicit none
    real(kp) :: hni2_fmin
    real(kp), intent(in) :: alpha
    
    hni2_fmin = hni_f(alpha)
    
  end function hni2_fmin

  
!returns V/M**4
  function hni2_norm_potential(x,alpha,f)
    implicit none
    real(kp) :: hni2_norm_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: f

    hni2_norm_potential = hni_norm_potential(x,alpha,f)

  end function hni2_norm_potential


!returns the first derivative of the potential with respect to x=phi/f, divided by M**4
  function hni2_norm_deriv_potential(x,alpha,f)
    implicit none
    real(kp) :: hni2_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: f

    hni2_norm_deriv_potential = hni_norm_deriv_potential(x,alpha,f)

  end function hni2_norm_deriv_potential


!returns the second derivative of the potential with respect to x=phi/f, divided by M**4
  function hni2_norm_deriv_second_potential(x,alpha,f)
    implicit none
    real(kp) :: hni2_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: f

    hni2_norm_deriv_second_potential = hni_norm_deriv_second_potential(x,alpha,f)

  end function hni2_norm_deriv_second_potential


!epsilon1(x)
  function hni2_epsilon_one(x,alpha,f)    
    implicit none
    real(kp) :: hni2_epsilon_one
    real(kp), intent(in) :: x,alpha,f

    hni2_epsilon_one = hni_epsilon_one(x,alpha,f)    

  end function hni2_epsilon_one


!epsilon2(x)
  function hni2_epsilon_two(x,alpha,f)    
    implicit none
    real(kp) :: hni2_epsilon_two
    real(kp), intent(in) :: x,alpha,f

    hni2_epsilon_two = hni_epsilon_two(x,alpha,f) 

  end function hni2_epsilon_two

!epsilon3(x)
  function hni2_epsilon_three(x,alpha,f)    
    implicit none
    real(kp) :: hni2_epsilon_three
    real(kp), intent(in) :: x,alpha,f

    hni2_epsilon_three = hni_epsilon_three(x,alpha,f)

  end function hni2_epsilon_three


!return the maximal allowed valued for xend. If eps>1, this is the
!smallest root of eps1=1, oterwise this is the value at which the
!potential is minimal
  function hni2_xendmax(alpha,f)
    implicit none
    real(kp) :: hni2_xendmax
    real(kp), intent(in) :: alpha,f
    real(kp), dimension(2) :: xepsone

    if (alpha.lt.hni2_alphamax(f)) then

       hni2_xendmax = hni_x_potmin(alpha,f)

    else

       xepsone = hni_x_epsoneunity(alpha,f)
       
       hni2_xendmax = xepsone(1)

    endif
                  
  end function hni2_xendmax



  function hni2_numacc_xendmax(alpha,f)
    implicit none
    real(kp) :: hni2_numacc_xendmax
    real(kp), intent(in) :: alpha,f

    real(kp), dimension(2) :: xnumacc

    xnumacc = hni_numacc_x_epsonenull(alpha,f)

    hni2_numacc_xendmax = min(hni2_xendmax(alpha,f) - epsilon(1._kp) , xnumacc(2))

  end function hni2_numacc_xendmax

  

  function hni2_numacc_efoldmax(alpha,f)
    implicit none
    real(kp) :: hni2_numacc_efoldmax
    real(kp), intent(in) :: alpha,f

    real(kp) :: xinimin, xendmax    
    
    xendmax = hni2_numacc_xendmax(alpha,f)
    xinimin = hni_numacc_xinimin(alpha,f)

    hni2_numacc_efoldmax = -hni_efold_primitive(xendmax,alpha,f) &
         + hni_efold_primitive(xinimin,alpha,f)


  end function hni2_numacc_efoldmax
  

  
  function hni2_numacc_xendmin(efold,alpha,f)
    implicit none
    real(kp) :: hni2_numacc_xendmin
    real(kp), intent(in) :: efold,alpha,f

    real(kp) :: efoldMax
    real(kp) :: xinimin, xendmax
   
    
    xinimin = hni_numacc_xinimin(alpha,f)
    
    efoldMax = hni2_numacc_efoldmax(alpha,f)
    
    if (efold.gt.efoldMax) then
       write(*,*)'hni2_numacc_xendmin: not enough efolds!'
       write(*,*)'efold requested= efold maxi= ',efold,efoldMax       
       stop
    endif
    
    hni2_numacc_xendmin = hni2_x_trajectory(efold,xinimin,alpha,f)
    
  end function hni2_numacc_xendmin

  

!returns x at bfold=-efolds before the end of inflation
  function hni2_x_trajectory(bfold,xend,alpha,f)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,f
    real(kp) :: hni2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: hni2Data

    if (bfold.le.0._kp) then
       mini = epsilon(1._kp)
       maxi = xend
    else       
       mini = xend
       maxi = hni_x_potmin(alpha,f)
    endif   
    
    hni2Data%real1 = alpha
    hni2Data%real2 = f
    hni2Data%real3 = -bfold + hni_efold_primitive(xend,alpha,f)

    hni2_x_trajectory = zbrent(find_hni_x_trajectory,mini,maxi,tolFind,hni2Data)

  end function hni2_x_trajectory



end module hni2sr
