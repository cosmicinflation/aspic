!slow-roll functions for the radiative plateau quadratic inflation potential
!
!V(phi) = M**4 x**4 (1- 2(1-alpha) ln(x**2) + 2(beta+1) ln(x**2)**2 )
!
!x = phi/phi0

module rpqtisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent, QuarticRoots, Sort
  use specialinf, only :  cei

  implicit none

  private

  public rpqti_norm_potential, rpqti_norm_deriv_potential, rpqti_norm_deriv_second_potential
  public rpqti_epsilon_one, rpqti_epsilon_two, rpqti_epsilon_three
  public rpqti_x_endinf, rpqti_efold_primitive, rpqti_x_trajectory
  public rpqti_xstar_brackets

 
contains
!returns V/M**4
  function rpqti_norm_potential(x,phi0,alpha,beta)
    implicit none
    real(kp) :: rpqti_norm_potential
    real(kp), intent(in) :: x,phi0,alpha,beta

    rpqti_norm_potential = x**4*(1._kp-2._kp*(1._kp-alpha)*log(x**2)+ &
                            2._kp*(beta+1._kp)*log(x**2)**2)

  end function rpqti_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function rpqti_norm_deriv_potential(x,phi0,alpha,beta)
    implicit none
    real(kp) :: rpqti_norm_deriv_potential
    real(kp), intent(in) :: x,phi0,alpha,beta

   rpqti_norm_deriv_potential = 4._kp*x**3*(alpha+2._kp*log(x**2)*(alpha+beta+ &
                                (1._kp+beta)*log(x**2)))

  end function rpqti_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function rpqti_norm_deriv_second_potential(x,phi0,alpha,beta)
    implicit none
    real(kp) :: rpqti_norm_deriv_second_potential
    real(kp), intent(in) :: x,phi0,alpha,beta

    rpqti_norm_deriv_second_potential = 4._kp*x**2*(7._kp*alpha+4._kp*beta+ &
            2._kp*log(x**2)*(4._kp+3._kp*alpha+7._kp*beta+3._kp*(1._kp+beta)*log(x**2)))

  end function rpqti_norm_deriv_second_potential



!epsilon_one(x)
  function rpqti_epsilon_one(x,phi0,alpha,beta)    
    implicit none
    real(kp) :: rpqti_epsilon_one
    real(kp), intent(in) :: x,phi0,alpha,beta
    
    rpqti_epsilon_one = (8._kp*(alpha+2._kp*log(x**2)*(alpha+beta+(1._kp+beta)* &
                        log(x**2)))**2)/(phi0**2*(x+2._kp*x*log(x**2)*(-1._kp+alpha+ &
                        (1._kp+beta)*log(x**2)))**2)

 end function rpqti_epsilon_one


!epsilon_two(x)
  function rpqti_epsilon_two(x,phi0,alpha,beta)    
    implicit none
    real(kp) :: rpqti_epsilon_two
    real(kp), intent(in) :: x,phi0,alpha,beta
    
    rpqti_epsilon_two = (8._kp*(alpha*(-7._kp+4._kp*alpha)-4._kp*beta+2._kp*log(x**2)* &
                        (-4._kp-3._kp*beta+alpha*(4._kp+alpha+4._kp*beta)+log(x**2)* &
                        (5._kp+2._kp*alpha**2+alpha*(-1._kp+3._kp*beta)+beta*(7._kp+ &
                        4._kp*beta)+2._kp*(1._kp+beta)*log(x**2)*(-1._kp+2._kp*alpha+beta+ &
                        (1._kp+beta)*log(x**2))))))/(phi0**2*x**2*(1._kp+2._kp*log(x**2)* &
                        (-1._kp+alpha+(1._kp+beta)*log(x**2)))**2)

  end function rpqti_epsilon_two


!epsilon_three(x)
  function rpqti_epsilon_three(x,phi0,alpha,beta)    
    implicit none
    real(kp) :: rpqti_epsilon_three
    real(kp), intent(in) :: x,phi0,alpha,beta
    logical :: isnotalreadythere

    rpqti_epsilon_three = (8._kp*(alpha+2._kp*log(x**2)*(alpha+beta+(1._kp+beta)* &
                        log(x**2)))*(8._kp+alpha*(13._kp+2._kp*alpha*(-21._kp+8._kp*alpha)- &
                        24._kp*beta)+18._kp*beta+2._kp*log(x**2)*(-6._kp+alpha*(-31._kp+8._kp* &
                        alpha+6._kp*alpha**2)-23._kp*beta+24._kp*(-2._kp+alpha)*alpha* &
                        beta-24._kp*beta**2+log(x**2)*(-5._kp+2._kp*alpha*(-6._kp+alpha* &
                        (9._kp+alpha))-33._kp*beta+18._kp*alpha*(1._kp+alpha)*beta+24._kp* &
                        (-1._kp+alpha)*beta**2+2._kp*log(x**2)*(2._kp*(-3._kp+alpha*(8._kp+ &
                        (-1._kp+alpha)*alpha)+beta+alpha*(11._kp+2._kp*alpha)*beta+(7._kp+ &
                        6._kp*alpha)*beta**2+4._kp*beta**3)+(1._kp+beta)*log(x**2)*(10._kp- &
                        7._kp*alpha+6._kp*alpha**2+5._kp*(2._kp+alpha)*beta+6._kp*beta**2+2._kp* &
                        (1._kp+beta)*log(x**2)*(-2._kp+3._kp*alpha+beta+(1._kp+beta)*log(x**2)))))))) &
                        /(phi0**2*x**2*(1._kp+2._kp*log(x**2)*(-1._kp+alpha+(1._kp+beta)* &
                        log(x**2)))**2*(alpha*(-7._kp+4._kp*alpha)-4._kp*beta+2._kp*log(x**2)* &
                        (-4._kp-3._kp*beta+alpha*(4._kp+alpha+4._kp*beta)+log(x**2)*(5._kp+2._kp* &
                        alpha**2+alpha*(-1._kp+3._kp*beta)+beta*(7._kp+4._kp*beta)+2._kp*(1._kp+beta)* &
                        log(x**2)*(-1._kp+2._kp*alpha+beta+(1._kp+beta)*log(x**2))))))


  end function rpqti_epsilon_three

!returns the number of non-vanishing local extremas of eps1 and the ordered values of these extremas
  function rpqti_extremas_epsilon_one(alpha,beta)
    implicit none
    real(kp), dimension(1:5) :: rpqti_extremas_epsilon_one
    real(kp), intent(in) :: alpha,beta
    real(kp) :: a,b,c,d,e,phi0junk ! coefficients of the quartic polynomial
    real(kp), dimension(1:5) :: coefficients
    complex(kp), dimension(1:4) :: roots
    integer :: i,j,k ! number of distinct real roots that are not flat inflection points
    logical :: isnotalreadythere

    ! coefficients of a quartic polynomial a x^4 + b x^3 + c x^2 + d x + e =0
    a = 4._kp+8._kp*beta+4._kp*beta**2
    b = -4._kp+8._kp*alpha+8._kp*alpha*beta+4._kp*beta**2
    c = 10._kp-2._kp*alpha+4._kp*alpha**2+14._kp*beta+6._kp*alpha*beta+8._kp*beta**2
    d = -8._kp+8._kp*alpha+2._kp*alpha**2-6._kp*beta+8._kp*alpha*beta
    e = -7._kp*alpha+4._kp*alpha**2-4._kp*beta

    coefficients(1) = e
    coefficients(2) = d
    coefficients(3) = c
    coefficients(4) = b
    coefficients(5) = a

    CALL QuarticRoots(coefficients,roots)

    phi0junk = 1._kp

    k=0
    DO i =1, 4
        if (aimag(roots(i)) == 0._kp .and. rpqti_epsilon_three(exp(real(roots(i))/2._kp),phi0junk,alpha,beta) .ne. 0._kp) then! real root which is not a flat inflection point
            isnotalreadythere = .true.
            DO j=1,i-1 ! avoids multiple counting
                if ( real(roots(j)) == real(roots(i)) ) then
                    isnotalreadythere = .false.
                end if
            END DO
            if (isnotalreadythere) then
                k=k+1
                rpqti_extremas_epsilon_one(k+1)=exp(real(roots(i))/2._kp)
            endif
        end if
    END DO
    rpqti_extremas_epsilon_one(1) = real(k,kp)

    call Sort(rpqti_extremas_epsilon_one, 2, 1+k)

  end function rpqti_extremas_epsilon_one

!returns the minimal and maximal values of x to bracket xend
  function rpqti_xend_brackets(phi0,alpha,beta)
    implicit none
    real(kp), dimension(1:2) :: rpqti_xend_brackets
    real(kp), intent(in) :: phi0,alpha,beta
    real(kp) :: x1, x2, x3, x4, Delta, xmin_eps1, xmax_eps1, eps11, eps12, eps13, eps14, &
                mineps, maxeps, xa, xb
    integer :: Nepsroots, Nroots_lt_xa,i
    real(kp), dimension(1:5) :: extremas_epsilon_one

    extremas_epsilon_one = rpqti_extremas_epsilon_one(alpha,beta)
    Nepsroots = int(extremas_epsilon_one(1))
    x1 = extremas_epsilon_one(2)
    x2 = extremas_epsilon_one(3)
    x3 = extremas_epsilon_one(4)
    x4 = extremas_epsilon_one(5)

    Delta = beta**2+alpha*(alpha-2._kp)

    xmin_eps1 = epsilon(1._kp)
    xmax_eps1 = 1._kp/epsilon(1._kp)

    if (Delta <0._kp) then ! the potential is a monotonous function of the field
        if (Nepsroots == 0) then ! eps1 is a monotonous decreasing function of the field
            xmax_eps1 = 1._kp/epsilon(1._kp)
        else if (Nepsroots == 1) then ! this is a priori impossible
            write(*,*) 'one local extrema in epsilon_one, this is impossible'
        else if (Nepsroots ==2) then ! eps1 decreases, reaches a non-vanishing minimum, increases, reaches a local maximum, then goes to zero
            eps11 = rpqti_epsilon_one(x1,phi0,alpha,beta)
            if (eps11 .ge. 1._kp) then
                xmin_eps1 = x2
                xmax_eps1 = 1._kp/epsilon(1._kp)
            else if (eps11 <1._kp) then
                xmax_eps1 = x1
            end if
        else if (Nepsroots ==3) then ! this is a priori impossible
            write(*,*) '3 local extremas in epsilon_one, this is impossible'
        else if (Nepsroots ==4) then ! eps1 decreases, has a local min, a local max, a local min, a local max, then goes to 0
            eps11 = rpqti_epsilon_one(x1,phi0,alpha,beta)
            eps12 = rpqti_epsilon_one(x2,phi0,alpha,beta)
            eps13 = rpqti_epsilon_one(x3,phi0,alpha,beta)
            eps14 = rpqti_epsilon_one(x3,phi0,alpha,beta)
            if (eps11 <1._kp) then
                xmax_eps1 = x1
            else if (eps13 <1._kp) then
                xmax_eps1 = x3
                xmin_eps1 = x2
            else
                xmin_eps1 = x4
                xmax_eps1 = 1._kp/epsilon(1._kp)
            end if
        end if

    else if (Delta .ge. 0._kp) then ! the potential is flat at xa and xb, inflation proceeds at x<xa in this model
        xa = exp((-(alpha+beta)-sqrt(Delta))/(4._kp*(beta+1._kp)))
        xb = exp((-(alpha+beta)+sqrt(Delta))/(4._kp*(beta+1._kp)))
        Nroots_lt_xa =0
        do i =1,Nepsroots
            if (extremas_epsilon_one(i+1) < xa) then
                Nroots_lt_xa = Nroots_lt_xa +1
            end if
        end DO
        if (Nroots_lt_xa ==0) then ! before the flat point of the potential, eps1 monotonously decreases
            xmax_eps1 = xa
        else if (Nroots_lt_xa ==1) then ! this is a priori impossible
            write(*,*) 'one local extrema in epsilon_one before flat inflection point in the potential, this is impossible'
        else if (Nroots_lt_xa ==2) then ! eps1 decreases, reaches a local min, then a local max, then goes to 0 at xa
            eps11 = rpqti_epsilon_one(x1,phi0,alpha,beta)
            eps12 = rpqti_epsilon_one(x2,phi0,alpha,beta)
            if (eps11 <1._kp) then
                xmax_eps1 = x1
            else
                xmin_eps1 = x2
                xmax_eps1 = xa
            end if
        else if (Nroots_lt_xa ==3) then ! this is a priori impossible
            write(*,*) '3 local extrema in epsilon_one before flat inflection point in the potential, this is impossible'
        else if (Nroots_lt_xa ==4) then ! this is a priori impossible
            write(*,*) '3 local extrema in epsilon_one before flat inflection point in the potential, this is impossible'
        end if

    end if

    rpqti_xend_brackets(1) = xmin_eps1
    rpqti_xend_brackets(2) = xmax_eps1

  end function rpqti_xend_brackets

!returns the minimal and maximal values of x to bracket xstar
  function rpqti_xstar_brackets(phi0,alpha,beta)
    implicit none
    real(kp), dimension(1:2) :: rpqti_xstar_brackets
    real(kp), intent(in) :: phi0,alpha,beta
    real(kp) :: xend, Delta, xstar_min, xstar_max, xa

    xend = rpqti_x_endinf(phi0,alpha,beta)
    xstar_min = xend*(1._kp+epsilon(1._kp))
    Delta = beta**2+alpha*(alpha-2._kp)

    if (Delta <0._kp) then ! the potential is a monotonous function of the field
        xstar_max = 1._kp/epsilon(1._kp)
    else
        xa = exp((-(alpha+beta)-sqrt(Delta))/(4._kp*(beta+1._kp)))
        xstar_max = xa*(1._kp-epsilon(1._kp))
    end if

    rpqti_xstar_brackets(1) = xstar_min
    rpqti_xstar_brackets(2) = xstar_max

  end function rpqti_xstar_brackets


!returns the value for x=phi/phi0 defined as epsilon1=1, where inflation ends
  function rpqti_x_endinf(phi0,alpha,beta)
    implicit none
    real(kp), intent(in) :: phi0,alpha,beta
    real(kp) :: rpqti_x_endinf
    real(kp), dimension(1:2) :: xend_brackets
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rpqtiData

    xend_brackets = rpqti_xend_brackets(phi0,alpha,beta)
    mini = xend_brackets(1)
    maxi = xend_brackets(2)

    rpqtiData%real1 = phi0
    rpqtiData%real2 = alpha
    rpqtiData%real3 = beta
    
    rpqti_x_endinf = zbrent(find_rpqti_x_endinf,mini,maxi,tolFind,rpqtiData)

  end function rpqti_x_endinf

  function find_rpqti_x_endinf(x,rpqtiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: rpqtiData
    real(kp) :: find_rpqti_x_endinf
    real(kp) :: phi0,alpha,beta

    phi0 = rpqtiData%real1
    alpha = rpqtiData%real2
    beta = rpqtiData%real3

    find_rpqti_x_endinf = rpqti_epsilon_one(x,phi0,alpha,beta)-1._kp
   
  end function find_rpqti_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function rpqti_efold_primitive(x,phi0,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,phi0,alpha,beta
    real(kp) :: rpqti_efold_primitive
    complex(kp) :: Delta, A, xa, xb, expinta, expintb, result

    Delta = cmplx(beta**2+alpha*(alpha-2._kp),0._kp,kp)
    A = -(1._kp+beta)/sqrt(Delta)
    xa = exp((-(alpha+beta)-sqrt(Delta))/(4._kp*(beta+1._kp)))
    xb = exp((-(alpha+beta)+sqrt(Delta))/(4._kp*(beta+1._kp)))

    expinta =  cei(2._kp*log(x/xa))
    expintb =  cei(2._kp*log(x/xb))

    result =0.125_kp*phi0**2*(x**2-(1._kp-A)/2._kp*xa**2*(expinta)-(1._kp+A)/2._kp*xb**2*(expintb))

    rpqti_efold_primitive = real(result)

  end function rpqti_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rpqti_x_trajectory(bfold,xend,phi0,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, phi0,alpha,beta, xend
    real(kp) :: rpqti_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    real(kp), dimension(1:2) :: xstar_brackets
    type(transfert) :: rpqtiData

    xstar_brackets = rpqti_xstar_brackets(phi0,alpha,beta)
    mini = xstar_brackets(1)
    maxi = xstar_brackets(2)


    rpqtiData%real1 = phi0
    rpqtiData%real2 = alpha
    rpqtiData%real3 = beta
    rpqtiData%real4 = -bfold + rpqti_efold_primitive(xend,phi0,alpha,beta)
    
    rpqti_x_trajectory = zbrent(find_rpqti_x_trajectory,mini,maxi,tolFind,rpqtiData)
       
  end function rpqti_x_trajectory

  function find_rpqti_x_trajectory(x,rpqtiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: rpqtiData
    real(kp) :: find_rpqti_x_trajectory
    real(kp) :: phi0,alpha,beta,NplusNuend

    phi0= rpqtiData%real1
    alpha= rpqtiData%real2
    beta= rpqtiData%real3
    NplusNuend = rpqtiData%real4

    find_rpqti_x_trajectory = rpqti_efold_primitive(x,phi0,alpha,beta) - NplusNuend
   
  end function find_rpqti_x_trajectory

  
end module rpqtisr
