!slow-roll functions for the MSSM inflation potential
!
!V(phi) = M**4 [ x**2 - alpha x**6 + alpha,beta x**10 ]
!
!x = phi/Mp

module mssmisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : atan_ito_log
  implicit none

  private

  public  mssmi_norm_potential, mssmi_epsilon_one, mssmi_epsilon_two, mssmi_epsilon_three
  public  mssmi_x_endinf, mssmi_efold_primitive, mssmi_x_trajectory
  public  mssmi_norm_deriv_potential, mssmi_norm_deriv_second_potential
  public  mssmi_x_epsilon1_min,mssmi_alpha_min

 
contains
!returns V/M**4
  function mssmi_norm_potential(x,alpha,beta)
    implicit none
    real(kp) :: mssmi_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    mssmi_norm_potential = x**2-alpha*x**6+beta*x**10

  end function mssmi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function mssmi_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: mssmi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta

   mssmi_norm_deriv_potential = 2._kp*(x-3._kp*alpha*x**5+5._kp*beta*x**9)

  end function mssmi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function mssmi_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: mssmi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta

    mssmi_norm_deriv_second_potential = 2._kp*(1._kp-15._kp*alpha*x**4+45._kp*beta*x**8)

  end function mssmi_norm_deriv_second_potential



!epsilon_one(x)
  function mssmi_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: mssmi_epsilon_one
    real(kp), intent(in) :: x,alpha,beta
    
    mssmi_epsilon_one =(2._kp*(1._kp-3._kp*alpha*x**4+5._kp*beta*x**8)**2)/(x-alpha*x**5+beta*x**9)**2
    
  end function mssmi_epsilon_one


!epsilon_two(x)
  function mssmi_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: mssmi_epsilon_two
    real(kp), intent(in) :: x,alpha,beta
    
    mssmi_epsilon_two =(4._kp*(1._kp+4._kp*alpha*x**4+3._kp*alpha**2*x**8 + & 
                       beta*x**8*(-26._kp+5._kp*beta*x**8)))/(x-alpha*x**5+beta*x**9)**2
    
  end function mssmi_epsilon_two


!epsilon_three(x)
  function mssmi_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: mssmi_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
    
    mssmi_epsilon_three = (4._kp*(-1._kp+3._kp*alpha*x**4-5._kp*beta*x**8)* &
                          (-1._kp+3._kp*alpha**3*x**12+3._kp*alpha**2*x**8* &
                          (7._kp-5._kp*beta*x**8)+beta*x**8*(-87._kp-5._kp*beta* &
                          x**8*(-33._kp+beta*x**8))-3._kp*alpha*x**4*(-3._kp+ &
                          beta*x**8*(18._kp+5._kp*beta*x**8))))/((x-alpha*x**5+ &
                          beta*x**9)**2*(1._kp+4._kp*alpha*x**4+3._kp*alpha**2*x**8 + & 
                          beta*x**8*(-26._kp+5._kp*beta*x**8)))
    
  end function mssmi_epsilon_three

!Returns the position of the first local minimum of epsilon1
  function mssmi_x_epsilon1_min(alpha,beta)   
    implicit none
    real(kp) :: mssmi_x_epsilon1_min
    real(kp), intent(in) :: alpha,beta
    complex(kp) :: delta,BigDelta,sigma,BigSigma,x_eps2NUL
    
    if (alpha**2/beta<20._kp/9._kp) then

	delta=9._kp*alpha**4-156._kp*alpha**2*beta+736._kp*beta**2
	BigDelta=27._kp*alpha**8-11808._kp*alpha**4*beta**2+ &
                 153088._kp*alpha**2*beta**3-430336._kp*beta**4
	sigma=27._kp*alpha**6-702._kp*alpha**4*beta+6624._kp*alpha**2* &
              beta**2-12896._kp*beta**3+6._kp*sqrt(15._kp)*beta*sqrt(BigDelta)
	BigSigma=-6._kp*alpha**2+52._kp*beta+delta/(sigma**(1._kp/3._kp))+sigma**(1._kp/3._kp)

	x_eps2NUL=(1._kp/(2._kp*sqrt(15._kp)*beta)*(sqrt(BigSigma)-&
                  sqrt(156._kp*beta-18._kp*alpha**2-BigSigma &
                  -24._kp*sqrt(15._kp)*alpha*beta/(sqrt(BigSigma)))))**(0.25_kp)
    
	mssmi_x_epsilon1_min = real(x_eps2NUL,kp)

    else

	mssmi_x_epsilon1_min = (3._kp*alpha/(10._kp*beta)*(1._kp-sqrt(1._kp-20._kp*beta/(9._kp*alpha**2))))**(0.25_kp)

    endif   
    
  end function mssmi_x_epsilon1_min



!returns x at the end of inflation defined as epsilon1=1
  function mssmi_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: mssmi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mssmiData

    mini = 0._kp
    maxi = mssmi_x_epsilon1_min(alpha,beta)*(1._kp-epsilon(1._kp)) !Position of the first local minimum of epsilon1
  

    mssmiData%real1 = alpha
    mssmiData%real2 = beta
    
    mssmi_x_endinf = zbrent(find_mssmi_x_endinf,mini,maxi,tolFind,mssmiData)
   
   
  end function mssmi_x_endinf

  function find_mssmi_x_endinf(x,mssmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: mssmiData
    real(kp) :: find_mssmi_x_endinf
    real(kp) :: alpha,beta

    alpha = mssmiData%real1
    beta = mssmiData%real2

    find_mssmi_x_endinf = mssmi_epsilon_one(x,alpha,beta)-1._kp
   
  end function find_mssmi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function mssmi_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: mssmi_efold_primitive
    complex(kp) ::aplus,aminus,bplus,bminus

    if (beta.eq.0._kp.and.alpha.eq.0._kp) stop 'mssmi_efold_primitive: alpha=0 and beta=0!'

    if (alpha.lt.10._kp**(-4.).and.beta.lt.10._kp**(-4.)) then
    ! This is an approximated formula, valid when alpha<<1 , beta<<1 , and which leads to better numerical convergence
    mssmi_efold_primitive = real(x**2/4._kp)
    else

    aplus=-1.5_kp*alpha+0.5_kp*sqrt(complex(9._kp*alpha**2-20._kp*beta,kp))
    aminus=-1.5_kp*alpha-0.5_kp*sqrt(complex(9._kp*alpha**2-20._kp*beta,kp))
    bplus=(2._kp*aplus+alpha)/(aplus-aminus)
    bminus=(2._kp*aminus+alpha)/(aminus-aplus)

    mssmi_efold_primitive = real(x**2/20._kp+bplus/(10._kp*sqrt(aplus))*atan_ito_log(sqrt(aplus)*x**2) &
                            +bminus/(10._kp*sqrt(aminus))*atan_ito_log(sqrt(aminus)*x**2))


    endif

  end function mssmi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function mssmi_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, beta, xend
    real(kp) :: mssmi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mssmiData

  
    mini = xend

    if (alpha**2/beta>20._kp/9._kp) then
	maxi = mssmi_x_epsilon1_min(alpha,beta) !local maximum of the potential
    else
	maxi=100._kp*mssmi_x_epsilon1_min(alpha,beta)
    endif

    mssmiData%real1 = alpha
    mssmiData%real2 = beta
    mssmiData%real3 = -bfold + mssmi_efold_primitive(xend,alpha,beta)
    
    mssmi_x_trajectory = zbrent(find_mssmitraj,mini,maxi,tolFind,mssmiData)
       
  end function mssmi_x_trajectory

  function find_mssmitraj(x,mssmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: mssmiData
    real(kp) :: find_mssmitraj
    real(kp) :: alpha,beta,NplusNuend

    alpha= mssmiData%real1
    beta = mssmiData%real2
    NplusNuend = mssmiData%real3

    find_mssmitraj = mssmi_efold_primitive(x,alpha,beta) - NplusNuend
   
  end function find_mssmitraj


!Returns the prior alpha_min(beta) such that the first local minimum of epsilon1 is less than one
  function mssmi_alpha_min(beta)    
    implicit none
    real(kp), intent(in) :: beta   
    real(kp) :: mssmi_alpha_min
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mssmiData

    if (beta<0.000796037) then 
	mssmi_alpha_min = 0._kp
    else

    mini = 0._kp
    maxi = 10._kp
  

    mssmiData%real1 = beta

    mssmi_alpha_min = zbrent(find_mssmi_alpha_min,mini,maxi,tolFind,mssmiData)

    endif
   
  end function mssmi_alpha_min

  function find_mssmi_alpha_min(alpha,mssmiData)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: mssmiData
    real(kp) :: find_mssmi_alpha_min
    real(kp) :: beta

    beta= mssmiData%real1

    find_mssmi_alpha_min = mssmi_epsilon_one(mssmi_x_epsilon1_min(alpha,beta),alpha,beta)-1._kp
   
  end function find_mssmi_alpha_min


 


  
end module mssmisr
