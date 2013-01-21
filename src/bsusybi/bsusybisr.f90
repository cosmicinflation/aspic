!slow-roll functions for the brane SUSY breaking inflation potential
!
!V(phi) = M^4 [exp( sqrt(6) x ) + exp( sqrt(6) gamma x )]
!
!x = phi/Mp

module bsusybisr
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public bsusybi_norm_potential, bsusybi_norm_deriv_potential, bsusybi_norm_deriv_second_potential
  public bsusybi_epsilon_one, bsusybi_epsilon_two,bsusybi_epsilon_three
  public bsusybi_efold_primitive, bsusybi_x_trajectory
  public bsusybi_xendmax

 
contains
!returns V/M^4
  function bsusybi_norm_potential(x,gammaBSUSYB)
    implicit none
    real(kp) :: bsusybi_norm_potential
    real(kp), intent(in) :: x,gammaBSUSYB

    bsusybi_norm_potential = exp(sqrt(6._kp)*x)+exp(sqrt(6._kp)*gammaBSUSYB*x)

  end function bsusybi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function bsusybi_norm_deriv_potential(x,gammaBSUSYB)
    implicit none
    real(kp) :: bsusybi_norm_deriv_potential
    real(kp), intent(in) :: x,gammaBSUSYB

   bsusybi_norm_deriv_potential = sqrt(6._kp)*(exp(sqrt(6._kp)*x)+ &
        gammaBSUSYB*exp(sqrt(6._kp)*gammaBSUSYB*x))

  end function bsusybi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function bsusybi_norm_deriv_second_potential(x,gammaBSUSYB)
    implicit none
    real(kp) :: bsusybi_norm_deriv_second_potential
    real(kp), intent(in) :: x,gammaBSUSYB

    bsusybi_norm_deriv_second_potential = 6._kp*(exp(sqrt(6._kp)*x)+ &
         gammaBSUSYB**2*exp(sqrt(6._kp)*gammaBSUSYB*x))

  end function bsusybi_norm_deriv_second_potential

!epsilon1(x)
  function bsusybi_epsilon_one(x,gammaBSUSYB)    
    implicit none
    real(kp) :: bsusybi_epsilon_one
    real(kp), intent(in) :: x,gammaBSUSYB
    
    bsusybi_epsilon_one = (3._kp*(exp(sqrt(6._kp)*x) &
         + exp(sqrt(6._kp)*gammaBSUSYB*x)*gammaBSUSYB)**2) &
         /(exp(sqrt(6._kp)*x)+exp(sqrt(6._kp)*gammaBSUSYB*x))**2
    
  end function bsusybi_epsilon_one


!epsilon2(x)
  function bsusybi_epsilon_two(x,gammaBSUSYB)    
    implicit none
    real(kp) :: bsusybi_epsilon_two
    real(kp), intent(in) :: x,gammaBSUSYB
    
    bsusybi_epsilon_two = -((12._kp*exp(sqrt(6._kp)*(1._kp+gammaBSUSYB)*x)* &
         (-1._kp+gammaBSUSYB)**2)/(exp(sqrt(6._kp)*x)+ &
         exp(sqrt(6._kp)*gammaBSUSYB*x))**2)
    
  end function bsusybi_epsilon_two

!epsilon3(x)
  function bsusybi_epsilon_three(x,gammaBSUSYB)    
    implicit none
    real(kp) :: bsusybi_epsilon_three
    real(kp), intent(in) :: x,gammaBSUSYB
    
    bsusybi_epsilon_three = -((6._kp*(exp(sqrt(6._kp)*x)-exp(sqrt(6._kp)*gammaBSUSYB*x))* &
         (-1._kp+gammaBSUSYB)*(exp(sqrt(6._kp)*x)+exp(sqrt(6._kp)*gammaBSUSYB*x)* &
         gammaBSUSYB))/(exp(sqrt(6._kp)*x)+exp(sqrt(6._kp)*gammaBSUSYB*x))**2)
    
  end function bsusybi_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function bsusybi_efold_primitive(x,gammaBSUSYB)
    implicit none
    real(kp), intent(in) :: x,gammaBSUSYB
    real(kp) :: bsusybi_efold_primitive
    
     if (gammaBSUSYB.eq.1._kp) stop 'bsusybi_efold_primitive: gamma=1 is singular'

    bsusybi_efold_primitive = 1/sqrt(6._kp)*x-1._kp/(6._kp*gammaBSUSYB) *&
                              log(abs(1._kp+gammaBSUSYB*exp(sqrt(6._kp)*(gammaBSUSYB-1._kp)*x)))

  end function bsusybi_efold_primitive
 


!returns x at bfold=-efolds before the end of inflation
  function bsusybi_x_trajectory(bfold,xend,gammaBSUSYB)
    implicit none
    real(kp), intent(in) :: bfold,xend,gammaBSUSYB
    real(kp) :: bsusybi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: bsusybiData

    mini = xend
    maxi = 1._kp/(sqrt(6._kp)*(gammaBSUSYB-1._kp))* &
           log((sqrt(3._kp)-1._kp)/(1-gammaBSUSYB*sqrt(3._kp)))    !Value of x such that epsilon1=1

    bsusybiData%real1 = gammaBSUSYB
    bsusybiData%real2 = -bfold + bsusybi_efold_primitive(xend,gammaBSUSYB)
    
    bsusybi_x_trajectory = zbrent(find_bsusybi_x_trajectory,mini,maxi,tolFind,bsusybiData)
    
  end function bsusybi_x_trajectory

  function find_bsusybi_x_trajectory(x,bsusybiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: bsusybiData
    real(kp) :: find_bsusybi_x_trajectory
    real(kp) :: gammaBSUSYB,NplusPrimEnd

    gammaBSUSYB=bsusybiData%real1
    NplusPrimEnd = bsusybiData%real2

    find_bsusybi_x_trajectory = bsusybi_efold_primitive(x,gammaBSUSYB) - NplusPrimEnd
   
  end function find_bsusybi_x_trajectory

!Returns the maximum value of xend in order to realize the required -bdolstar e-folds.
  function bsusybi_xendmax(efold, gammaBSUSYB) 
    implicit none
    real(kp), intent(in) :: gammaBSUSYB, efold
    real(kp) :: bsusybi_xendmax
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini, maxi
    type(transfert) :: bsusybiData
   
    maxi = 1._kp/(sqrt(6._kp)*(gammaBSUSYB-1._kp))* &
         log((sqrt(3._kp)-1._kp)/(1-gammaBSUSYB*sqrt(3._kp)))    !Value of x such that epsilon1=1
    mini = - 10._kp*sqrt(6._kp) * efold + maxi

    bsusybiData%real1 = gammaBSUSYB
    bsusybiData%real2 = -efold
    
    bsusybi_xendmax = zbrent(find_bsusybi_xendmax,mini,maxi,tolFind,bsusybiData)

  end function bsusybi_xendmax

  function find_bsusybi_xendmax(x,bsusybiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: bsusybiData
    real(kp) :: find_bsusybi_xendmax,xmaxi
    real(kp) :: gammaBSUSYB,bfoldstar
!Value of x such that epsilon1=1

    xmaxi = 1._kp/(sqrt(6._kp)*(gammaBSUSYB-1._kp))* &
         log((sqrt(3._kp)-1._kp)/(1-gammaBSUSYB*sqrt(3._kp)))
    gammaBSUSYB=bsusybiData%real1
    bfoldstar = bsusybiData%real2

    find_bsusybi_xendmax = bfoldstar-(bsusybi_efold_primitive(x,gammaBSUSYB)- &
                           bsusybi_efold_primitive(xmaxi,gammaBSUSYB))
   
  end function find_bsusybi_xendmax



end module bsusybisr
