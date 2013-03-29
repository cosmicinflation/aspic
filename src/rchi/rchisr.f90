!slow-roll functions for the radiatively corrected higgs inflation models
!
!V(phi) = M**4 [ 1 - 2 Exp( -2x/sqrt(6) ) + AI/(16 pi**2) x/sqrt(6) ]
!
!x = phi/Mp

module rchisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : polylog, lambert
  use inftools, only : zbrent
  implicit none

  private

  public  rchi_norm_potential, rchi_epsilon_one, rchi_epsilon_two, rchi_epsilon_three
  public  rchi_x_endinf, rchi_efold_primitive, rchi_x_trajectory
  public  rchi_norm_deriv_potential, rchi_norm_deriv_second_potential
 
contains
!returns V/M**4
  function rchi_norm_potential(x,AI)
    implicit none
    real(kp) :: rchi_norm_potential
    real(kp), intent(in) :: x,AI

    rchi_norm_potential = 1._kp-2._kp*exp(-sqrt((2._kp/3._kp))*x)+ &
                          (AI*x)/(16._kp*sqrt(6._kp)*acos(-1._kp)**2)

  end function rchi_norm_potential


!returns the first derivative of the potential with respect to x=phi/Mp, divided by M**4
  function rchi_norm_deriv_potential(x,AI)
    implicit none
    real(kp) :: rchi_norm_deriv_potential
    real(kp), intent(in) :: x,AI

   rchi_norm_deriv_potential = (64._kp*exp(-sqrt((2._kp/3._kp))*x)+ &
                               AI/acos(-1._kp)**2)/(16._kp*sqrt(6._kp))

  end function rchi_norm_deriv_potential



!returns the second derivative of the potential with respect to x=phi/Mp, divided by M**4
  function rchi_norm_deriv_second_potential(x,AI)
    implicit none
    real(kp) :: rchi_norm_deriv_second_potential
    real(kp), intent(in) :: x,AI

    rchi_norm_deriv_second_potential = -(4._kp/3._kp)*exp(-sqrt((2._kp/3._kp))*x)

  end function rchi_norm_deriv_second_potential




!epsilon_one(x)
  function rchi_epsilon_one(x,AI)    
    implicit none
    real(kp) :: rchi_epsilon_one
    real(kp), intent(in) :: x,AI
    
    !Approximated Formula commonly used in the litterature
    !rchi_epsilon_one = 4._kp/3._kp*(exp(-sqrt(2._kp/3._kp)*x)+ &
    !                   AI/(64._kp*acos(-1._kp)**2))**2

    rchi_epsilon_one = (3._kp*(AI*exp(sqrt(2._kp/3._kp)*x)+ &
                       64._kp*acos(-1._kp)**2)**2)/(-192._kp*acos(-1._kp)**2+ &
                       exp(sqrt(2._kp/3._kp)*x)*(96._kp*acos(-1._kp)**2+ &
                       sqrt(6._kp)*AI*x))**2




    
  end function rchi_epsilon_one


!epsilon_two(x)
  function rchi_epsilon_two(x,AI)    
    implicit none
    real(kp) :: rchi_epsilon_two
    real(kp), intent(in) :: x,AI

    !Approximated Formula commonly used in the litterature
    !rchi_epsilon_two = 4._kp*rchi_epsilon_one(x,AI)+8._kp/3._kp*exp(-sqrt(2._kp/3._kp)*x)

    rchi_epsilon_two = (4._kp*exp(sqrt(2._kp/3._kp)*x)*(3._kp*AI**2* &
                       exp(sqrt(2._kp/3._kp)*x)+6144._kp*acos(-1._kp)**4+ &
                       64._kp*AI*acos(-1._kp)**2*(6._kp+sqrt(6._kp)*x)))/ &
                       (-192._kp*acos(-1._kp)**2+exp(sqrt(2._kp/3._kp)*x)* &
                       (96._kp*acos(-1._kp)**2+sqrt(6._kp)*AI*x))**2




    
  end function rchi_epsilon_two


!epsilon_three(x)
  function rchi_epsilon_three(x,AI)    
    implicit none
    real(kp) :: rchi_epsilon_three
    real(kp), intent(in) :: x,AI
    
    rchi_epsilon_three = (12._kp*(AI*exp(sqrt(2._kp/3._kp)*x)+64._kp*acos(-1._kp)**2)* &
                         (3._kp*AI**3*exp(2._kp*sqrt(2._kp/3._kp)*x)+ &
                         2048._kp*acos(-1._kp)**4*(96._kp*acos(-1._kp)**2+ &
                         AI*(9._kp+sqrt(6._kp)*x))+32._kp*exp(sqrt(2._kp/3._kp)*x)* &
                         acos(-1._kp)**2*(3072._kp*acos(-1._kp)**4+32._kp* &
                         AI*acos(-1._kp)**2*(9._kp+2._kp*sqrt(6._kp)*x)+ &
                         AI**2*(18._kp+3._kp*sqrt(6._kp)*x+2._kp*x**2))))/ &
                         ((3._kp*AI**2*exp(sqrt(2._kp/3._kp)*x)+6144._kp*acos(-1._kp)**4+ &
                         64._kp*AI*acos(-1._kp)**2*(6._kp+sqrt(6._kp)*x))* &
                         (-192._kp*acos(-1._kp)**2+exp(sqrt(2._kp/3._kp)*x)* &
                         (96._kp*acos(-1._kp)**2+sqrt(6._kp)*AI*x))**2)
    
  end function rchi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function rchi_x_endinf(AI)
    implicit none
    real(kp), intent(in) :: AI
    real(kp) :: rchi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rchiData
    
    if (AI.gt.0) then !Uses the 0 branch of the Lambert Function
        rchi_x_endinf = 1._kp/sqrt(2._kp)-(16._kp*sqrt(6._kp)* &
                        acos(-1._kp)**2)/AI+sqrt(3._kp/2._kp)* &
                        lambert((64._kp*(3._kp+sqrt(3._kp))* &
                        exp(-(1._kp/sqrt(3._kp))+(32._kp*acos(-1._kp)**2)/AI)* &
                        acos(-1._kp)**2)/(3._kp*AI),0)
   
    elseif  (AI.lt.0)  then !Uses the -1 branch of the Lambert Function
    rchi_x_endinf = 1._kp/sqrt(2._kp)-(16._kp*sqrt(6._kp)* &
                        acos(-1._kp)**2)/AI+sqrt(3._kp/2._kp)* &
                        lambert((64._kp*(3._kp+sqrt(3._kp))* &
                        exp(-(1._kp/sqrt(3._kp))+(32._kp*acos(-1._kp)**2)/AI)* &
                        acos(-1._kp)**2)/(3._kp*AI),-1)

    end if

    if ( abs(AI) .lt. 1._kp ) then !In that case the expressions above become numerically too difficult to track down and a zbrent is used

    mini=sqrt(6._kp)*log(2._kp)/2._kp
    if (AI .lt. 0._kp ) then !maximum of the potential
      maxi=min(abs(sqrt(3._kp/2._kp)*log(abs(AI)/(64._kp*acos(-1._kp)**2))),100._kp)
    else
      maxi = max(abs(AI)*1000._kp/(48._kp*acos(-1._kp)**2)*sqrt(3._kp/2._kp)+10._kp,100._kp)
    endif
    rchiData%real1 = AI
    
    rchi_x_endinf = zbrent(find_rchi_x_endinf,mini,maxi,tolFind,rchiData)

    if  (AI.eq.0) then !In that case the expression is singluar

       rchi_x_endinf = sqrt(3._kp/2._kp)*log(2._kp+2._kp/sqrt(3._kp))

    endif

    end if

  end function rchi_x_endinf

  function find_rchi_x_endinf(x,rchiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: rchiData
    real(kp) :: find_rchi_x_endinf
    real(kp) :: AI

    AI = rchiData%real1

    find_rchi_x_endinf = rchi_epsilon_one(x,AI)-1._kp
   
  end function find_rchi_x_endinf



!this is integral(V(phi)/V'(phi) dphi)
  function rchi_efold_primitive(x,AI)
    implicit none
    real(kp), intent(in) :: x,AI
    real(kp) :: rchi_efold_primitive

     !Approximated Trajectory in  the vacuum dominated approximation
     !rchi_efold_primitive = 48._kp*acos(-1._kp)**2/AI*log(abs(1._kp+ & 
     !                      AI/(64._kp*acos(-1._kp)**2)*exp(sqrt(2._kp/3._kp)*x)))


    if (AI.eq.0._kp) then 
       rchi_efold_primitive = 0.25_kp*(3._kp*exp(sqrt(2._kp/3._kp)*x)-2._kp*sqrt(6._kp)*x)
    else

    rchi_efold_primitive = -sqrt(3._kp/2._kp)*x+48._kp*acos(-1._kp)**2/AI* &
                           (1._kp+AI/(32._kp*acos(-1._kp)**2)* &
                           (1._kp+sqrt(2._kp/3._kp)*x))*log(abs(1._kp+AI/(64._kp*acos(-1._kp)**2)* &
                           exp(sqrt(2._kp/3._kp)*x)))+1.5_kp*real(polylog(cmplx( &
                           -AI/(64._kp*acos(-1._kp)**2)*exp(sqrt(2._kp/3._kp)*x),0._kp), &
                           cmplx(2._kp,0._kp)),kp)

    endif



  end function rchi_efold_primitive




!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rchi_x_trajectory(bfold,xend,AI)
    implicit none
    real(kp), intent(in) :: bfold,AI,xend
    real(kp) :: rchi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rchiData

    mini=xend*(1._kp+epsilon(1._kp))
    if (AI .lt. 0._kp .and. abs(AI) .lt. 64._kp*acos(-1._kp)**2) then !maximum of the potential
      maxi=min(abs(sqrt(3._kp/2._kp)*log(abs(AI)/(64._kp*acos(-1._kp)**2))),xend+100._kp)
    else
      maxi = max(abs(AI)*1000._kp/(48._kp*acos(-1._kp)**2)*sqrt(3._kp/2._kp)+xend,xend+10._kp,100._kp)
    endif

    rchiData%real1 = AI
    rchiData%real2 = -bfold + rchi_efold_primitive(xend,AI)
    
    rchi_x_trajectory = zbrent(find_rchi_x_trajectory,mini,maxi,tolFind,rchiData)
       
  end function rchi_x_trajectory

  function find_rchi_x_trajectory(x,rchiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: rchiData
    real(kp) :: find_rchi_x_trajectory
    real(kp) :: AI,NplusNuend

    AI= rchiData%real1
    NplusNuend = rchiData%real2

    find_rchi_x_trajectory = rchi_efold_primitive(x,AI) - NplusNuend
   
  end function find_rchi_x_trajectory




  
end module rchisr
