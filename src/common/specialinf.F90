module specialprec
  use infprec, only : kp
  implicit none
  integer, parameter :: dp = kp

end module specialprec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                        !
!                 MODULE FOR HYP_2F1                     !
!                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



MODULE HYP_2F1_MODULE
  use specialprec, only : dp,kp
  !--------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, PARAMETER :: PR=kp,IPR=KIND(1)
  REAL(PR)     :: EPS15=1.0D-15
  REAL(PR)     :: ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.50D0
  REAL(PR)     :: M_PI=3.14159265358979323846D0
  REAL(PR)     :: M_PI_2=1.57079632679489661923D0
  REAL(PR)     :: M_1_PI=0.31830988618379067154D0 
CONTAINS
  !
  FUNCTION INF_NORM(Z)
    COMPLEX(PR),INTENT(IN) :: Z
    REAL(PR)  :: INF_NORM
    INF_NORM=MAX(ABS(REAL(Z,PR)),ABS(AIMAG(Z)))
    RETURN
  END FUNCTION INF_NORM
  !
  FUNCTION TANZ(Z)
    COMPLEX(PR),INTENT(IN) :: Z
    COMPLEX(PR) :: TANZ
    TANZ=SIN(Z)/COS(Z)
    RETURN
  END FUNCTION TANZ
  !
  FUNCTION LOG1P(Z)
    COMPLEX(PR),INTENT(IN) :: Z
    REAL(PR) :: X,XP1,LOG1P_X
    REAL(PR) :: Y,YX,YX2,YX2P1,LOG1P_YX2
    REAL(PR) :: RE_LOG1P,IM_LOG1P
    COMPLEX(PR) :: LOG1P
    IF(INF_NORM(Z).LT.ONE) THEN
       X = REAL(Z,PR); XP1 = X+ONE
       IF(XP1.EQ.ONE) THEN
          LOG1P_X = X
       ELSE
          LOG1P_X = LOG(XP1)*X/(XP1-ONE)
       ENDIF
       Y = AIMAG(Z)
       YX = Y/XP1; YX2 = YX*YX; YX2P1 = YX2+ONE
       IF(YX2P1.EQ.ONE) THEN
          LOG1P_YX2 = YX2
       ELSE
          LOG1P_YX2 = LOG(YX2P1)*YX2/(YX2P1-ONE)
       ENDIF
       RE_LOG1P = LOG1P_X + HALF*LOG1P_YX2
       IM_LOG1P = ATAN2(Y,XP1)
       LOG1P = CMPLX(RE_LOG1P,IM_LOG1P,PR)
       RETURN
    ELSE
       LOG1P=LOG(ONE+Z)
       RETURN
    ENDIF
  END FUNCTION LOG1P
  !
  FUNCTION EXPM1(Z)
    COMPLEX(PR),INTENT(IN) :: Z
    REAL(PR) :: X,EXPM1_X,EXP_X,Y,SIN_HALF_Y
    REAL(PR) :: RE_EXPM1,IM_EXPM1
    COMPLEX(PR) :: EXPM1
    IF(INF_NORM(Z).LT.ONE) THEN
       X = REAL(Z,PR); EXP_X = EXP(X)
       Y = AIMAG(Z); SIN_HALF_Y=SIN(HALF*Y)
       IF(EXP_X.EQ.ONE) THEN
          EXPM1_X = X
       ELSE 
          EXPM1_X = (EXP_X-ONE)*X/LOG(EXP_X)
       ENDIF
       RE_EXPM1 = EXPM1_X-TWO*EXP_X*SIN_HALF_Y*SIN_HALF_Y 
       IM_EXPM1 = EXP_X*SIN(Y)
       EXPM1 = CMPLX(RE_EXPM1,IM_EXPM1,PR)
       RETURN
    ELSE
       EXPM1=EXP(Z)-ONE
       RETURN
    ENDIF
  END FUNCTION EXPM1
  !
END MODULE HYP_2F1_MODULE



module specialinf
  use specialprec, only : dp,kp
  implicit none

#ifdef NOF08
  interface atan
     module procedure atan_ito_log
  end interface atan

  interface atanh
     module procedure atanh_ito_log
  end interface atanh
#endif

interface expint
   module procedure expintkp
end interface expint

interface lambert
!   module procedure lambert_root
   module procedure lambert_iter
end interface lambert


contains


  SUBROUTINE xermsg(text1, text2, text3)
 	CHARACTER (LEN= *), INTENT(IN)  :: text1
 	CHARACTER (LEN= *), INTENT(IN)  :: text2
	CHARACTER (LEN= *), INTENT(IN)  :: text3
	WRITE(*, '(6a)') ' Error in call to ', text1, ' routine: ', text2, ' Mesg: ', text3
	return
  END SUBROUTINE xermsg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                        !
!                 SPECIAL FUNCTIONS                      !
!                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function log_otherbranchcut(z)
    implicit none
    complex(kp), intent(in) :: z
    complex(kp) :: log_otherbranchcut


    if (aimag(z)>0) then
       log_otherbranchcut  = log(z)
    else
       log_otherbranchcut  = log(z)+2*acos(-1._kp)*cmplx(0._kp,1._kp,kp)
    endif


  end function log_otherbranchcut

  function ellipticK(k2)
    use hyp_2f1_module, only : M_PI_2
    implicit none
    real(kp) :: ellipticK
    real(kp), intent(in) :: k2
    real(kp), save :: k2sav = -1._kp
    real(kp), save :: Ksav = -1._kp
!$omp threadprivate(k2sav, Ksav)
    if (k2.eq.k2sav) then
       ellipticK = Ksav
       return
    endif

    if ((k2.gt.1._kp).or.(k2.lt.0._kp)) then
       stop 'ellipticK: k2 > 1 or k2 < 0'
    endif
    
    ellipticK = M_PI_2 * hypergeom_2F1(0.5_kp,0.5_kp,1._kp,k2)
    k2sav = k2
    Ksav = ellipticK

  end function ellipticK



  function ellipticE(k2)
    use hyp_2f1_module, only : M_PI_2
    implicit none
    real(kp) :: ellipticE
    real(kp), intent(in) :: k2
    real(kp), save :: k2sav = -1._kp
    real(kp), save :: Esav = -1._kp
!$omp threadprivate(k2sav, Esav)

    if (k2.eq.k2sav) then
       ellipticE = Esav
       return
    endif

    if ((k2.gt.1._kp).or.(k2.lt.0._kp)) then
       stop 'ellipticE: k2 > 1 or k2 < 0'
    endif
    
    ellipticE = M_PI_2 * hypergeom_2F1(-0.5_kp,0.5_kp,1._kp,k2)
    k2sav = k2
    Esav = ellipticE

  end function ellipticE


  function hypergeom_2F1(a,b,c,z)
    implicit none
    real(kp), intent(in) :: a,b,c,z
    real(kp) :: hypergeom_2F1
    complex(dp) :: ac,bc,cc,zc
    complex(dp) :: res

    ac = cmplx(a,0,dp)
    bc = cmplx(b,0,dp)
    cc = cmplx(c,0,dp)
    zc = cmplx(z,0,dp)

    res = HYP_2F1(ac,bc,cc,zc)

    hypergeom_2F1  = real(res,kp)


  end function hypergeom_2F1

!Halley's iterative method
  function lambert_iter(x,n)
    implicit none
    real(kp) :: lambert_iter
    real(kp), intent(in) :: x
    integer, intent(in) :: n
    real(kp) :: w, expw, num, den
    real(kp) :: wzero, wnew
    real(kp), parameter :: xbig = 10._kp
    real(kp), parameter :: xsmall = 1._kp/xbig
    
    real(kp), parameter :: tolIter = 100._kp*epsilon(1._kp)

    integer :: count
    integer, parameter :: iterMax = 100000
    real(kp) :: error
    real(kp), parameter :: errorWarn = 1d-10


    if (x.eq.-exp(-1._kp)) then
       lambert_iter = -1._kp
       return
    endif

    select case (n)

    case (0)

       if (x.ge.xbig) then

          wzero = log(x) - log(log(x))

       elseif ( ((x.gt.-exp(-1._kp)).and.(x.lt.0._kp)) &
            .or. ((x.gt.0._kp).and.(x.lt.xbig)) ) then

          wzero = 0._kp

       elseif (x.eq.0._kp) then

          lambert_iter = 0._kp
          return

       else
          write(*,*)'x= ',x
          stop 'Lambert: x not defined in branch 0'
       endif

    case (-1)

       if ((x.gt.-exp(-1._kp)).and.(x.le.-xsmall)) then

          wzero = -2._kp

       elseif ((x.gt.-xsmall).and.(x.lt.0._kp)) then

          wzero = log(-x) - log(-log(-x))
                 
       else
          write(*,*)'x= ',x
          stop 'Lambert: x not defined in branch -1'
       endif

    case default
       stop 'Lambert function branch incorrect'

    end select
      


    w = wzero
    count = 0
    do
       expw = exp(w)
       num = w * expw - x
       den = expw*(w+1._kp) - (w+2._kp)*(w*expw-x)/(2._kp*w+2._kp)
              
       wnew = w - num/den

       error = abs(wnew-w)
       if (error.le.tolIter*abs(w)) exit
       count = count + 1

       if (count.gt.iterMax) then
          if (error.gt.errorWarn) then
             write(*,*)'x= w= deltaw= ',x,w,error
             write(*,*)'Lambert may be inaccurate!'
          endif
          exit
       endif

       w = wnew

    enddo

    lambert_iter = w

  end function lambert_iter


!lambert function of the nth-branch (n=0 or -1)
  REAL(kp) FUNCTION lambert_root(x,n) 
    implicit none
    REAL(kp), INTENT(IN)::x
    INTEGER, INTENT(IN)::n
    REAL(kp) ::y1,y2,y,L1,L2

    real(kp), parameter :: xSwitchTaylor = 0.0001_kp
    real(kp), parameter :: TargetPrec = 100._kp * epsilon(1._kp)

    IF(x<-1._kp/exp(1._kp)) THEN
       WRITE(*,*) 'Impossible to evaluate lambert Function: arg < -1/e'
    END IF
    IF(n==0) THEN
       y1=-1._kp
       y2=MAX(1._kp,x)
       IF((y1*exp(y1)-x)*(y2*exp(y2)-x)>0._kp) THEN
          WRITE(*,*) 'Wrong boundaries in lambert-0 Function calculation'
       END IF
!Uses a Taylor Series to avoid numerical errors
       IF (abs(x).lt.xSwitchTaylor) THEN 
          y=x-x**2+1.5_kp*x**3-8._kp/3._kp*x**4+125._kp/24._kp*x**5-54._kp/5._kp*x**6+ &
               16807._kp/720._kp*x**7-16384._kp/315._kp*x**8+531441._kp/4480._kp*x**9- &
               156250._kp/567._kp*x**(10)
       ELSE
!          DO WHILE(abs(y1-y2) > 1d-6)
          DO WHILE(abs(y1-y2)/abs(y1+y2) > TargetPrec)
             y=(y1+y2)/2._kp
             IF((y1*exp(y1)-x)*(y*exp(y)-x)<0._kp) THEN
                y2=y
             ELSE
                y1=y
             END IF
          END DO
       END IF
    END IF
    IF(n==-1) THEN
       IF(x>0._kp) THEN
          WRITE(*,*) 'Impossible to evaluate lambert-1 Function: x>0'
       END IF
       y1=-1.d6
       y2=-1._kp
       IF((y1*exp(y1)-x)*(y2*exp(y2)-x)>0._kp) THEN
          WRITE(*,*) 'Wrong boundaries in lambert-1 Function calculation'
       END IF
!       DO WHILE(abs(y1-y2)>1d-6)
       DO WHILE (abs(y1-y2)/abs(y1+y2) > TargetPrec)
          y=(y1+y2)/2._kp
          IF((y1*exp(y1)-x)*(y*exp(y)-x)<0._kp) THEN
             y2=y
          ELSE
             y1=y
          END IF
       END DO
!Uses an asymptotic solution to avoid numerical errors
       IF (-x .lt. epsilon(1._kp)) THEN 
          L1=log(-x)
          L2=log(-log(-x))
          y=L1-L2+L2/L1+L2*(-2._kp+L2)/(2._kp*L1**2)+ &
               L2*(6._kp-9._kp*L2+2._kp*L2**2)/(6._kp*L1**3)+ &
               L2*(-12._kp+36._kp*L2-22._kp*L2**2+3._kp*L2**3)/(12._kp*L1**4)
       ENDIF
    END IF
    lambert_root=y
  END FUNCTION lambert_root

!Exponential Integral function ExpInt(n,x)=Int_1^Infinity exp(-xt)/(t^n) dt
  FUNCTION expintkp(n,x) result(expint)
    use infprec, only : euler
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(kp), INTENT(IN) :: x
    REAL(kp) :: expint
    INTEGER, PARAMETER :: MAXIT=1000
    REAL(kp), PARAMETER :: EPS=epsilon(x),BIG=huge(x)*EPS

    !Evaluates the exponential integral En(x).  Parameters: MAXIT is
    !the maximum allowed number of iterations; EPS is the desired
    !relative error, not smaller than the machine precision; BIG is a
    !number near the largest representable floating-point number;
    !EULER (in nrtype) is Euler's constant γ.

    INTEGER :: i,nm1
    REAL(kp) :: a,b,c,d,del,fact,h

    if (.not.((n >= 0).and.(x>=-0._kp).and.((x > 0._kp .or. n > 1)))) then
       stop 'expint: incorrect parameter input!'
    endif
    if (n == 0) then !Special case.
       expint=exp(-x)/x
       RETURN
    end if

    nm1=n-1
    if (x == 0._kp) then !Another special case.
       expint=1.0_kp/nm1
    else if (x > 1._kp) then !Lentz's algorithm (§5.2).
       b=x+n
       c=BIG
       d=1.0_kp/b
       h=d
       do i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.0_kp
          d=1.0_kp/(a*d+b) !Denominators cannot be zero.
          c=b+a/c
          del=c*d
          h=h*del
          if (abs(del-1.0_kp) <= EPS) exit
!          print *,'test',i,abs(del-1.0_kp)
       end do
       if (i > MAXIT) stop 'expint: continued fraction failed'
       expint=h*exp(-x)
    else !Evaluate series.
       if (nm1 /= 0) then !Set first term.
          expint=1.0_kp/nm1
       else
          expint=-log(x)-EULER
       end if
       fact=1.0_kp
       do i=1,MAXIT
          fact=-fact*x/i
          if (i /= nm1) then
             del=-fact/(i-nm1)
          else !ψ(n) appears here.
             del=fact*(-log(x)-EULER+sum(1.0_kp/arthint(1,1,nm1)))
          end if
          expint=expint+del
          if (abs(del) < abs(expint)*EPS) exit
       end do
       if (i > MAXIT) stop 'expint: series failed'
    end if
  END FUNCTION expintkp


! Complex exponential Integral, by VV
  function cei(z) result(resultcei)
    complex(kp), intent(in) :: z
    complex(kp) :: resultcei
    complex(kp) :: Gamma0
    integer :: i,imax
    
    if (z.eq.cmplx(0._kp)) stop 'cei: z= 0 + i 0'

    imax = 1000

    ! This computes Gamma[0,-z] by means of a continued fraction
    Gamma0=1._kp/(-z+(2._kp*real(imax+1,kp)-1._kp))

    do i=imax,1,-1
       Gamma0 = -z+2._kp*real(i,kp)-1._kp-real(i,kp)**2/Gamma0
    end do
    Gamma0 = exp(z)/Gamma0

    resultcei = -Gamma0+0.5_kp*(log(z)-log(1._kp/z))-log(-z)

  end function cei



  FUNCTION arthint(first,increment,n)
    INTEGER, INTENT(IN) :: first,increment,n
    INTEGER, DIMENSION(n) :: arthint
    INTEGER :: k,k2,temp
    integer, parameter :: NPAR_ARTH = 16
    integer, parameter :: NPAR2_ARTH = 8
    if (n > 0) arthint(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arthint(k)=arthint(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arthint(k)=arthint(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arthint(k+1:min(k2,n))=temp+arthint(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arthint

 


  FUNCTION ei(x) !Exponential Integral Ei(x)=-Int_{-x}^Infinity exp(-t)/t dt

    use infprec, only : euler
    IMPLICIT NONE
    REAL(kp), INTENT(IN) :: x
    REAL(kp) :: ei
    INTEGER, PARAMETER :: MAXIT=10000
    REAL(kp), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
    INTEGER :: k
    REAL(kp) :: fact,prev,sm,term
    if (x<0._kp) then   	!If x < 0, uses the fact that Ei(-x)=-expint(1,-x)
       ei=-expint(1,-x)
    else if(x.eq.0._kp) then !If x=0, gives a numerical value for -Infinity
       ei=-1._kp/epsilon(1._kp)
    else if (isnan(x)) then
       stop 'ei(x): input x is NaN'
    else
       !Computes the exponential integral Ei(x) for x > 0.
       !Parameters: MAXIT is the maximum number of iterations allowed; EPS is the relative error,
       !or absolute error near the zero of Ei at x = 0.3725; FPMIN is a number near the smallest
       !representable floating-point number; EULER (in nrtype) is Euler's constant γ.

       !		call assert(x > 0.0, 'ei arg')
       if (.not.(x > 0._kp)) then
          write(*,*)'x= ',x
          stop 'ei: x<=0'
       endif
       if (x < FPMIN) then !Special case: avoid failure of convergence test
          ei=log(x)+EULER !because of underflow.
       else if (x <= -log(EPS)) then !Use power series.
          sm=0.0_kp
          fact=1.0_kp
          do k=1,MAXIT
             fact=fact*x/k
             term=fact/k
             sm=sm+term
             if (term < EPS*sm) exit
          end do
          if (k > MAXIT)  stop 'series failed in ei'

          ei=sm+log(x)+EULER
       else !Use asymptotic series.
          sm=0.0 !Start with second term.
          term=1.0
          do k=1,MAXIT
             prev=term
             term=term*k/x
             if (term < EPS) exit !Since final sum is greater than one, term itself
             if (term < prev) then !approximates the relative error.
                sm=sm+term !Still converging: add new term.
             else !Diverging: subtract previous term and exit.
                sm=sm-prev
                exit
             end if
          end do
          if (k > MAXIT) stop 'asymptotic failed in ei'
          ei=exp(x)*(1.0_kp+sm)/x
       end if
    end if
  END FUNCTION ei


! Gamma_inv denotes the entire inverse of the Gamma function.
! F(z) means 2F1(a,b,c,z) with the a, b, c and z given as inputs 
! in the routine.
!
! Elementary functions and standard constants 
! are defined in the module.
! See N.J.~Higham, ``Accuracy and Stability of Numerical Algorithms'',
! SIAM, Philadelphia, 1996 for expm1 implementation.
! log1p follows instantly.


RECURSIVE FUNCTION LOG_GAMMA_CPLX(Z) RESULT(RES)
!----------------------------------------------------------------------
! Logarithm of Gamma[z] and Gamma inverse function
! ------------------------------------------------
!
! For log[Gamma[z]],if z is not finite 
! or is a negative integer, the program 
! returns an error message and stops.
! The Lanczos method is used. Precision : ~ 1E-15
! The method works for Re[z]>0.5 .
! If Re[z]<=0.5, one uses the formula Gamma[z].Gamma[1-z]=Pi/sin(Pi.z)
! log[sin(Pi.z)] is calculated with the Kolbig method 
! (K.S. Kolbig, Comp. Phys. Comm., Vol. 4, p.221(1972)): 
! If z=x+iy and y>=0, log[sin(Pi.z)]=log[sin(Pi.eps)]-i.Pi.n, 
! with z=n+eps so 0<=Re[eps]< 1 and n integer.
! If y>110, log[sin(Pi.z)]=-i.Pi.z+log[0.5]+i.Pi/2 
! numerically so that no overflow can occur.
! If z=x+iy and y< 0, log[Gamma(z)]=[log[Gamma(z*)]]*, 
! so that one can use the previous formula with z*.
!
! For Gamma inverse, Lanczos method is also used 
! with Euler reflection formula.
! sin (Pi.z) is calculated as sin (Pi.(z-n)) 
! to avoid inaccuracy with z = n + eps 
! with n integer and |eps| as small as possible.
!
!
! Variables:
! ----------
! x,y: Re[z], Im[z]
! log_sqrt_2Pi,log_Pi : log[sqrt(2.Pi)], log(Pi).
! sum : Rational function in the Lanczos method
! log_Gamma_z : log[Gamma(z)] value.
! c : table containing the fifteen coefficients in the expansion 
! used in the Lanczos method.
! eps,n : z=n+eps so 0<=Re[eps]< 1 and n integer for Log[Gamma].
!         z=n+eps and n integer 
!         so |eps| is as small as possible for Gamma_inv.
! log_const : log[0.5]+i.Pi/2
! g : coefficient used in the Lanczos formula. It is here 607/128.
! z,z_m_0p5,z_p_g_m0p5,zm1 : argument of the Gamma function, 
! z-0.5, z-0.5+g, z-1 
! res: returned value
!----------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: Z
  INTEGER(IPR) :: N,I
  REAL(PR)     :: X,Y,LOG_SQRT_2PI,G,LOG_PI,M_LN2,C(0:14)
  COMPLEX(PR)  :: GAMMA_SUM,Z_M_0P5,Z_P_G_M0P5,ZM1
  COMPLEX(PR)  :: LOG_CONST,I_PI,EPS,LOG_SIN_PI_Z,RES
  !
  M_LN2=0.69314718055994530942D0; X=REAL(Z,PR); Y=AIMAG(Z)
  IF((Z.EQ.NINT(X)).AND.(X.LE.ZERO)) &
       STOP 'Z IS NEGATIVE INTEGER IN LOG_GAMMA_CPLX'
  IF(X.GE.HALF) THEN
     LOG_SQRT_2PI=0.91893853320467274177D0; G=4.7421875D0
     Z_M_0P5=Z-HALF; Z_P_G_M0P5=Z_M_0P5+G; ZM1=Z-ONE
     C=(/ 0.99999999999999709182D0,57.156235665862923517D0,       &
          -59.597960355475491248D0,  14.136097974741747174D0,     &
          -0.49191381609762019978D0, 0.33994649984811888699D-4,   &
          0.46523628927048575665D-4, -0.98374475304879564677D-4,  &
          0.15808870322491248884D-3, -0.21026444172410488319D-3,  &
          0.21743961811521264320D-3, -0.16431810653676389022D-3,  &
          0.84418223983852743293D-4, -0.26190838401581408670D-4,  &
          0.36899182659531622704D-5 /)

     GAMMA_SUM=C(0)
     DO I=1,14
        GAMMA_SUM=GAMMA_SUM+C(I)/(ZM1+I)
     ENDDO
     RES=LOG_SQRT_2PI+LOG(GAMMA_SUM)+Z_M_0P5*LOG(Z_P_G_M0P5) &
          -Z_P_G_M0P5
     RETURN
  ELSE IF(Y.GE.ZERO) THEN
     IF(X.LT.NINT(X)) THEN
        N=NINT(X)-1
     ELSE
        N=NINT(X)
     ENDIF
     LOG_PI=1.1447298858494002D0
     LOG_CONST=CMPLX(-M_LN2,M_PI_2,PR); I_PI=CMPLX(ZERO,M_PI,PR)
     EPS=Z-N
     IF(Y.GT.110.0D0) THEN
        LOG_SIN_PI_Z=-I_PI*Z+LOG_CONST
     ELSE
        LOG_SIN_PI_Z=LOG(SIN(M_PI*EPS))-I_PI*N
     ENDIF
     RES=LOG_PI-LOG_SIN_PI_Z-LOG_GAMMA_CPLX(ONE-Z);
     RETURN
  ELSE
     RES=CONJG(LOG_GAMMA_CPLX(CONJG(Z)))
     RETURN
  ENDIF
END FUNCTION LOG_GAMMA_CPLX
!
!----------------------------------------------------------------------
! Inverse of the Gamma function [1/Gamma](z)
! ------------------------------------------
! It is calculated with the Lanczos method for Re[z] >= 0.5 
! and is precise up to 10^{-15}.
! If Re[z] <= 0.5, one uses the formula 
! Gamma[z].Gamma[1-z] = Pi/sin (Pi.z).
! sin (Pi.z) is calculated as sin (Pi.(z-n)) to avoid inaccuracy,
! with z = n + eps with n integer and |eps| as small as possible.
! 
! Variables 
! ---------
! z : argument of the function
! x: Re[z]
! eps,n : z = n + eps with n integer and |eps| as small as possible.
! res: returned value
!----------------------------------------------------------------------
RECURSIVE FUNCTION GAMMA_INV(Z) RESULT(RES)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: Z
  INTEGER(IPR) :: N,I
  REAL(PR)     :: X,LOG_SQRT_2PI,G,C(0:14)
  COMPLEX(PR)  :: RES,GAMMA_SUM,Z_M_0P5,Z_P_G_M0P5,ZM1,EPS
  !
  X=REAL(Z,PR)
  IF(X.GE.HALF) THEN
     LOG_SQRT_2PI=0.91893853320467274177D0; G=4.7421875D0
     Z_M_0P5=Z-HALF; Z_P_G_M0P5=Z_M_0P5+G; ZM1=Z-ONE
     C=(/ 0.99999999999999709182D0,57.156235665862923517D0,       &
          -59.597960355475491248D0,  14.136097974741747174D0,     &
          -0.49191381609762019978D0, 0.33994649984811888699D-4,   &
          0.46523628927048575665D-4, -0.98374475304879564677D-4,  &
          0.15808870322491248884D-3, -0.21026444172410488319D-3,  &
          0.21743961811521264320D-3, -0.16431810653676389022D-3,  &
          0.84418223983852743293D-4, -0.26190838401581408670D-4,  &
          0.36899182659531622704D-5 /)

     GAMMA_SUM=C(0)
     DO I=1,14
        GAMMA_SUM=GAMMA_SUM+C(I)/(ZM1+I);
     ENDDO
     RES=EXP(Z_P_G_M0P5-Z_M_0P5*LOG(Z_P_G_M0P5)-LOG_SQRT_2PI) &
          /GAMMA_SUM
     RETURN
  ELSE
     X=REAL(Z,PR); N=NINT(X)
     EPS=Z-N
     IF(MOD(N,2).EQ.0) THEN
        RES=SIN(M_PI*EPS)*M_1_PI/GAMMA_INV (ONE-Z)
        RETURN
     ELSE
        RES=-SIN(M_PI*EPS)*M_1_PI/GAMMA_INV (ONE-Z)
        RETURN
     ENDIF
  ENDIF
END FUNCTION GAMMA_INV
!----------------------------------------------------------------------
!
! Calculation of H(z,eps) = [Gamma(z+eps)/Gamma(z) - 1]/eps, with e and
! ---------------------------------------------------------------------
! z complex so z,z+eps are not negative integers and 0 <= |eps|oo < 0.1
! ---------------------------------------------------------------------
! The function H(z,eps) = [Gamma(z+eps)/Gamma(z) - 1]/e is calculated 
! here with the Lanczos method.
! For the Lanczos method, the gamma parameter, denoted as g, 
! is 4.7421875 and one uses a sum of 15 numbers with the table c[15], 
! so that it is precise up to machine accuracy.
! The H(z,eps) function is used in formulas occuring in1-z and 1/z 
! transformations (see Comp. Phys. Comm. paper).
!
! One must have z and z+eps not negative integers as otherwise 
! it is clearly not defined.
! As this function is meant to be precise for small |eps|oo, 
! one has to have 0 <= |eps|oo < 0.1 .
! Indeed, a direct implementation of H(z,eps) with Gamma_inv or 
! log_Gamma for |eps|oo >= 0.1 is numerically stable.
! The returned function has full numerical accuracy 
! even if |eps|oo is very small.
!
! eps not equal to zero
! ---------------------
! If Re(z) >= 0.5 or Re(z+eps) >= 0.5, one clearly has Re(z) > 0.4 
! and Re(z+eps) > 0.4, 
! so that the Lanczos summation can be used for both Gamma(z) 
! and Gamma(z+eps).
! One then has:
! log[Gamma(z+eps)/Gamma(z)] = 
! (z-0.5) log1p[eps/(z+g-0.5)] + eps log(z+g-0.5+eps) - eps 
! + log1p[-eps \sum_{i=1}^{14} c[i]/((z-1+i)(z-1+i+eps)) 
! / (c[0] + \sum_{i=1}^{14} c[i]/(z-1+i))]
! H(z,eps) = expm1[log[Gamma(z+eps)/Gamma(z)]]/eps .
!
! If Re(z) < 0.5 and Re(z+eps) < 0.5, 
! Euler reflection formula is used for both Gamma(z) and Gamma(z+eps).
! One then has: 
! H(z+eps,-eps) = [cos(pi.eps) + sin(pi.eps)/tan(pi(z-n))].H(1-z,-eps) 
! + (2/eps).sin^2(eps.pi/2) - sin(pi.eps)/(eps.tan(pi.(z-n)))
! H(1-z,-eps) is calculated with the Lanczos summation 
! as Re(1-z) >= 0.5 and Re(1-z-eps) >= 0.5 .
! z-n is used in tan(pi.z) instead of z to avoid inaccuracies 
! due the finite number of digits of pi.
! H(z,eps) = H(z+eps,-eps)/(1 - eps.H(z+eps,-eps)) 
! provides the final result.
!
! eps equal to zero
! -----------------
! It is obtained with the previous case and eps -> 0 :
! If Re(z) >= 0.5, one has:
! H(z,eps) = (z-0.5)/(z+g-0.5) + log(z+g-0.5) - 1 -
! \sum_{i=1}^{14} c[i]/((z-1+i)^2)/(c[0]+\sum_{i=1}^{14} c[i]/(z-1+i))
!
! If Re(z) < 0.5, one has:
! H(z,0) = H(1-z,0) - pi/tan(pi.(z-n))
!
! Variables
! ---------
! z,eps: input variables of the function H(z,eps)
! g,c[15]: double and table of 15 doubles defining the Lanczos sum 
! so that it provides the Gamma function 
! precise up to machine accuracy.
! eps_pz,z_m_0p5,z_pg_m0p5,eps_pz_pg_m0p5,zm1,zm1_p_eps: 
! z+eps,z-0.5,z+g-0.5,z+eps+g-0.5,z-1,z-1+eps
! x,eps_px: real parts of z and z+eps.
! n,m: closest integer ot the real part of z, same for z+eps.
! sum_num,sum_den: \sum_{i=1}^{14} c[i]/((z-1+i)(z-1+i+eps)) 
! and (c[0] + \sum_{i=1}^{14} c[i]/(z-1+i)). 
! They appear respectively as numerator and denominator in formulas.
! Pi_eps,term,T1_eps_z: pi.eps, sin (pi.eps)/tan(pi.(z-n)), 
! [cos(pi.eps) + sin(pi.eps)/tan(pi(z-n))].H(1-z,-eps)
! sin_Pi_2_eps,T2_eps_z,T_eps_z: sin^2(eps.pi/2), 
! (2/eps).sin^2(eps.pi/2) - sin(pi.eps)/(eps.tan(pi.(z-n))), 
! H(z+eps,-eps)
! res: returned value
!----------------------------------------------------------------------
RECURSIVE FUNCTION GAMMA_RATIO_DIFF_SMALL_EPS(Z,EPS) RESULT(RES)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: Z,EPS
  INTEGER(IPR) :: N,M,I
  REAL(PR)     :: G,X,EPS_PX,C(0:14)
  COMPLEX(PR)  :: RES,SUM_NUM,SUM_DEN
  COMPLEX(PR)  :: EPS_PZ,Z_M_0P5,Z_PG_M0P5,EPS_PZ_PG_M0P5,ZM1
  COMPLEX(PR)  :: CI_ZM1_PI_INV,PI_EPS,TT,T1_EPS_Z,SIN_PI_2_EPS
  COMPLEX(PR)  :: ZM1_P_EPS,T2_EPS_Z,T_EPS_Z
  !
  G=4.74218750D0
  IF(INF_NORM(EPS).GT.0.1D0) &
       STOP 'ONE MUST HAVE |EPS|< 0.1 IN GAMMA_RATIO_DIFF_SMALL_EPS'
  EPS_PZ=Z+EPS; Z_M_0P5=Z-HALF; Z_PG_M0P5=Z_M_0P5+G
  EPS_PZ_PG_M0P5=Z_PG_M0P5+EPS; ZM1=Z-ONE; ZM1_P_EPS=ZM1+EPS
  X=REAL(Z,PR); EPS_PX=REAL(EPS_PZ,PR); N=NINT(X); M=NINT(EPS_PX)
  IF((Z.EQ.N).AND.(N.LE.0)) THEN
     STOP 'Z IS NEGATIVE INTEGER IN GAMMA_RATIO_DIFF_SMALL_EPS'
  ENDIF
  IF((EPS_PZ.EQ.M).AND.(M.LE.0)) THEN
     STOP 'Z+EPS IS NEGATIVE INTEGER IN GAMMA_RATIO_DIFF_SMALL_EPS'
  ENDIF
  C=(/ 0.99999999999999709182D0,57.156235665862923517D0,     &
       -59.597960355475491248D0,14.136097974741747174D0,     &
       -0.49191381609762019978D0,0.33994649984811888699D-4,  &
       0.46523628927048575665D-4,-0.98374475304879564677D-4, &
       0.15808870322491248884D-3,-0.21026444172410488319D-3, &
       0.21743961811521264320D-3,-0.16431810653676389022D-3, &
       0.84418223983852743293D-4,-0.26190838401581408670D-4, &
       0.36899182659531622704D-5 /)
  IF((X.GE.HALF).OR.(EPS_PX.GE.HALF)) THEN
     SUM_NUM=ZERO;SUM_DEN=C(0)
     DO I=1,14
        CI_ZM1_PI_INV=C(I)/(ZM1+I)
        SUM_NUM=SUM_NUM+CI_ZM1_PI_INV/(ZM1_P_EPS+I)
        SUM_DEN=SUM_DEN+CI_ZM1_PI_INV
     ENDDO
     IF(EPS.NE.ZERO) THEN
        RES=EXPM1(Z_M_0P5*LOG1P(EPS/Z_PG_M0P5) &
             +EPS*LOG(EPS_PZ_PG_M0P5)-EPS+LOG1P(-EPS*SUM_NUM/SUM_DEN))&
             /EPS
        RETURN
     ELSE
        RES=Z_M_0P5/Z_PG_M0P5 &
             +LOG(EPS_PZ_PG_M0P5)-ONE-SUM_NUM/SUM_DEN
        RETURN
     ENDIF
  ELSE
     IF(EPS.NE.ZERO) THEN
        PI_EPS=M_PI*EPS
        TT=SIN(PI_EPS)/TANZ(M_PI*(Z-N))
        T1_EPS_Z=(COS(PI_EPS)+TT)*& 
             GAMMA_RATIO_DIFF_SMALL_EPS(ONE-Z,-EPS)
        SIN_PI_2_EPS=SIN(M_PI_2*EPS)
        T2_EPS_Z=(TWO*SIN_PI_2_EPS*SIN_PI_2_EPS-TT)/EPS
        T_EPS_Z=T1_EPS_Z+T2_EPS_Z
        RES=(T_EPS_Z/(ONE-EPS*T_EPS_Z))
        RETURN
     ELSE
        RES=GAMMA_RATIO_DIFF_SMALL_EPS(ONE-Z,-EPS) &
             -M_PI/TANZ(M_PI*(Z-N))
        RETURN
     ENDIF
  ENDIF
END FUNCTION GAMMA_RATIO_DIFF_SMALL_EPS
!
!----------------------------------------------------------------------
! Calculation of G(z,eps) = [Gamma_inv(z) - Gamma_inv(z+eps)]/eps 
! ---------------------------------------------------------------
! with e and z complex
!---------------------
! The G(z,eps) function is used in formulas occuring in 1-z 
! and 1/z transformations (see Comp. Phys. Comm. paper).
! Several case have to be considered for its evaluation. 
! eps is considered equal to zero 
! if z+eps and z are equal numerically.
!
! |eps|oo > 0.1
! -------------
! A direct evaluation with the values Gamma_inv(z) 
! and Gamma_inv(z+eps) is stable and returned.
!
! |eps|oo <= 0.1 with z+eps and z numerically different
! -----------------------------------------------------
! If z is a negative integer, z+eps is not, 
! so that G(z,eps) = -Gamma_inv(z+eps)/eps, 
! for which a direct evaluation is precise and returned.
! If z+eps is a negative integer, z is not, 
! so that G(z,eps) = Gamma_inv(z)/eps, 
! for which a direct evaluation is precise and returned.
! If both of them are not negative integers, 
! one looks for the one of z and z+eps 
! which is the closest to a negative integer.
! If it is z, one returns H(z,eps).Gamma_inv(z+eps). 
! If it is z+eps, one returns H(z+eps,-eps).Gamma_inv(z).
! Both values are equal, so that one chooses the one 
! which makes the Gamma ratio Gamma(z+eps)/Gamma(z) 
! in H(z,eps) the smallest in modulus.
!
! z+eps and z numerically equal
! -----------------------------
! If z is negative integer, G(z,0) = (-1)^(n+1) n!, 
! where z = -n, n integer, which is returned.
! If z is not negative integer, one returns H(z,eps).Gamma_inv(z+eps)
!
! Variables
! ---------
! z,eps: input variables of the function G(z,eps)
! eps_pz,x,eps_px: z+eps,real parts of z and z+eps.
! n,m: closest integer ot the real part of z, same for z+eps.
! fact,k: (-1)^(n+1) n!, returned when z = -n, n integer 
! and z and z+eps identical numerically (eps ~ 0). 
! It is calculated with integer index k.
! is_z_negative_integer,is_eps_pz_negative_integer: 
! true if z is a negative integer, false if not, same for z+eps.
! z_neg_int_distance, eps_pz_neg_int_distance: 
! |z + |n||oo, |z + eps + |m||oo. 
! If |z + |n||oo < |z + eps + |m||oo, 
! z is closer to the set of negative integers than z+eps.
! Gamma_inv(z+eps) is then of moderate modulus 
! if Gamma_inv(z) is very small. 
! If z ~ n, H(z,eps) ~ -1/eps, 
! that so returning 
! G(z,eps) = H(z,eps).Gamma_inv(z+eps) here is preferred.
! Same for |z + |n||oo > |z + eps + |m||oo with z <-> z+eps.
!
!----------------------------------------------------------------------
FUNCTION GAMMA_INV_DIFF_EPS(Z,EPS)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: Z,EPS
  INTEGER(IPR) :: M,N,K
  REAL(PR)     :: X,EPS_PX,FACT
  REAL(PR)     :: Z_NEG_INT_DISTANCE
  REAL(PR)     :: EPS_PZ_NEG_INT_DISTANCE
  COMPLEX(PR)  :: GAMMA_INV_DIFF_EPS,EPS_PZ
!  COMPLEX(PR)  :: GAMMA_RATIO_DIFF_SMALL_EPS,GAMMA_INV
  LOGICAL      :: IS_Z_NEG_INT,IS_EPS_PZ_NEG_INT

  EPS_PZ=Z+EPS; X=REAL(Z,PR); EPS_PX=REAL(EPS_PZ,PR)
  N=NINT(X); M=NINT(EPS_PX)
  IS_Z_NEG_INT=(Z.EQ.N).AND.(N.LE.0)
  IS_EPS_PZ_NEG_INT=(EPS_PZ.EQ.M).AND.(M.LE.0)
  IF(INF_NORM(EPS).GT.0.10D0) THEN
     GAMMA_INV_DIFF_EPS = (GAMMA_INV (Z) - GAMMA_INV (EPS_PZ))/EPS
     RETURN
  ELSE IF(EPS_PZ.NE.Z) THEN 
     IF(IS_Z_NEG_INT) THEN
        GAMMA_INV_DIFF_EPS = (-GAMMA_INV (EPS_PZ)/EPS)
        RETURN
     ELSE IF(IS_EPS_PZ_NEG_INT) THEN
        GAMMA_INV_DIFF_EPS = (GAMMA_INV (Z)/EPS)
        RETURN
     ELSE
        Z_NEG_INT_DISTANCE = INF_NORM (Z + ABS (N))
        EPS_PZ_NEG_INT_DISTANCE = INF_NORM (EPS_PZ + ABS (M))
        IF(Z_NEG_INT_DISTANCE.LT.EPS_PZ_NEG_INT_DISTANCE) THEN
           GAMMA_INV_DIFF_EPS= &
                GAMMA_RATIO_DIFF_SMALL_EPS (Z,EPS)*GAMMA_INV (EPS_PZ)
           RETURN
        ELSE
           GAMMA_INV_DIFF_EPS= &
                GAMMA_RATIO_DIFF_SMALL_EPS (EPS_PZ,-EPS)*GAMMA_INV (Z)
           RETURN
        ENDIF
     ENDIF
  ELSE IF(IS_Z_NEG_INT.AND.IS_EPS_PZ_NEG_INT) THEN
     FACT = -ONE;K=-1
     DO WHILE (K.GE.N) 
        FACT=FACT*K
        K=K-1 
     ENDDO
     GAMMA_INV_DIFF_EPS = FACT
     RETURN
  ELSE
     GAMMA_INV_DIFF_EPS = &
          GAMMA_RATIO_DIFF_SMALL_EPS (Z,EPS)*GAMMA_INV (EPS_PZ)
     RETURN
  ENDIF
END FUNCTION GAMMA_INV_DIFF_EPS
!----------------------------------------------------------------------
!
! Calculation of Gamma_inv(1-m-eps)/eps of the A(z) polynomial in 1-z
! -------------------------------------------------------------------
! and 1/z transformations
! -----------------------
! This value occurs in A(z) in 1-z and 1/z transformations 
! (see Comp. Phys. Comm. paper) for m > 0.
! Both cases of 1-m-eps numerically negative integer 
! or not have to be considered
! 
! 1-eps-m and 1-m numerically different
! -------------------------------------
! One returns Gamma_inv(1-m-eps)/eps directly 
! as its value is accurate.
! To calculate Gamma_inv(1-m-eps), 
! one uses the value Gamma_inv(1-eps), 
! needed in considered transformations,
! and one uses the equality 
! Gamma_inv(1-m-eps) = Gamma_inv(1-eps) \prod_{i=1}^{m} (1-eps-i) 
! for m > 0.
! It is trivially demonstrated 
! from the equality Gamma(x+1) = x.Gamma(x). 
! One Gamma function evaluation is removed this way 
! from the calculation.
! 
! 1-eps-m and 1-m numerically equal
! ---------------------------------
! This implies that 1-m-eps is negative integer numerically.
! Here, eps~0, so that one returns the limit of Gamma_inv(1-m-eps)/eps
! for eps -> 0, which is (-1)^m (m-1)!
!
! Variables
! ---------
! m,eps: variable inputs of the function 
! (m,eps) -> Gamma_inv(1-m-eps)/eps
! Gamma_inv_one_meps: Gamma_inv(1-eps), 
! previously calculated and here recycled 
! to quickly calculate Gamma_inv(1-m-eps).
! one_meps: 1-eps
!----------------------------------------------------------------------
FUNCTION A_SUM_INIT(M,EPS,GAMMA_INV_ONE_MEPS)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  INTEGER(IPR),INTENT(IN) :: M
  COMPLEX(PR),INTENT(IN) :: EPS,GAMMA_INV_ONE_MEPS
  INTEGER(IPR) :: N,I
  REAL(PR)     :: FACT
  COMPLEX(PR)  :: A_SUM_INIT,ONE_MEPS
  COMPLEX(PR)  :: GAMMA_INV_ONE_MEPS_MM
  !
  ONE_MEPS = ONE - EPS
  IF(ONE_MEPS-M.NE.1-M) THEN
     GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS
     DO I=1,M
        GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS_MM*(ONE_MEPS-I)
     ENDDO
     A_SUM_INIT=GAMMA_INV_ONE_MEPS_MM/EPS
     RETURN
  ELSE
     FACT=ONE
     DO N=2,M-1
        FACT=FACT*N
     ENDDO
     IF(MOD(M,2).EQ.0) THEN
        A_SUM_INIT=FACT
     ELSE
        A_SUM_INIT=-FACT
     ENDIF
     RETURN
  ENDIF
END FUNCTION A_SUM_INIT
!
!----------------------------------------------------------------------
! Calculation of the log of Gamma_inv(1-m-eps)/eps
! ------------------------------------------------
! See previous function. 
! It is used in case Gamma_inv(1-m-eps)/eps might overflow.
!
! Variables
! ---------
! m,eps: variable inputs of the function 
! (m,eps) -> log[Gamma_inv(1-m-eps)/eps]
! one_meps_mm: 1-eps-m
! i_Pi: i.Pi
! log_fact: logarithm of (-1)^m (m-1)!, 
! here defined as log((m-1)!) + i.Pi if m is odd.
!----------------------------------------------------------------------
FUNCTION LOG_A_SUM_INIT(M,EPS)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  INTEGER(IPR),INTENT(IN) :: M
  COMPLEX(PR),INTENT(IN) :: EPS
  INTEGER(IPR) :: N
  REAL(PR)     :: LOG_FACT
  COMPLEX(PR)  :: ONE_MEPS_MM,LOG_A_SUM_INIT
!  COMPLEX(PR)  :: LOG_GAMMA_CPLX
  !
  ONE_MEPS_MM=ONE-EPS-M
  IF(ONE_MEPS_MM.NE.1-M) THEN
     LOG_A_SUM_INIT=(-LOG_GAMMA_CPLX(ONE_MEPS_MM) - LOG(EPS))
     RETURN
  ELSE
     LOG_FACT=ZERO
     DO N=2,M-1
        LOG_FACT=LOG_FACT + LOG(real(N,pr))
     ENDDO
     IF(MOD(M,2).EQ.0) THEN
        LOG_A_SUM_INIT=LOG_FACT
     ELSE
        LOG_A_SUM_INIT=CMPLX(LOG_FACT,M_PI,PR)
     ENDIF
     RETURN
  ENDIF
END FUNCTION LOG_A_SUM_INIT
!----------------------------------------------------------------------
! Calculation of the first term of the B(z) power series
! ------------------------------------------------------
! in the 1-z transformation, divided by (1-z)^m
! ----------------------------------------------
! In the 1-z transformation, 
! the power series B(z) = \sum_{n=0}^{+oo} \beta_n (1-z)^n occurs 
! (see Comp. Phys. Comm. paper).
! The first term \beta_0, divided by (1-z)^m, is calculated here. 
! m is the closest integer to Re(c-a-b) >= 0 and eps = c-a-b-m.
!
! One has to consider |eps|oo > 0.1 and |eps|oo <= 0.1, 
! where 1-m-eps and 1-m can be different or equal numerically, 
! leading to some changes in this last case.
!
! |eps|oo > 0.1
! -------------
! One has \beta_0/(1-z)^m = [(a)_m (b)_m Gamma_inv(1-eps) 
! Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) Gamma_inv(m+1)
! - (1-z)^eps Gamma_inv(a) Gamma_inv(b) Gamma_inv(1+m+eps)]
! [Gamma(c)/eps], stable in this regime for a direct evaluation.
!
! The values of Gamma(c), Gamma_inv(a+m+eps) 
! and Gamma_inv(b+m+eps) were already calculated and recycled here.
! Gamma_inv(m+1) is calculated as 1/(m!).
!
! Gamma_inv(1+m+eps) is calculated from Gamma_inv(1-eps), 
! using the equalities:
! Gamma_inv(1-m-eps) = Gamma_inv(1-eps) \prod_{i=1}^{m} (1-eps-i), 
! where the product is 1 by definition if m = 0,
! Gamma_inv(1+m+eps) = (-1)^m sin (pi.eps)
! /[pi.(eps+m).Gamma_inv(1-m-eps)] 
! from Euler reflection formula, Gamma(x+1) = x.Gamma(x) equality, 
! and m+eps no zero.
! This scheme is much faster than 
! to recalculate Gamma_inv(1+m+eps) directly.
! 
! |eps|oo <= 0.1
! --------------
! The \beta_0/(1-z)^m expression is rewritten 
! so that it contains no instabilities:
! \beta_0/(1-z)^m = Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) 
! [(G(1,-eps) Gamma_inv(m+1) + G(m+1,eps))
! - Gamma_inv(1+m+eps) (G(a+m,eps) Gamma_inv(b+m+eps) 
! + G(b+m,eps) Gamma_inv(a+m)) 
! - E(log(1-z),eps) Gamma_inv(a+m) Gamma_inv(b+m) Gamma_inv(1+m+eps)] 
! (a)_m (b)_m Gamma(c)
!
! E(log(1-z),eps) is [(1-z)^eps - 1]/eps 
! if 1-m-eps and 1-m are different numerically, 
! and log(1-z) otherwise (eps ~ 0).
! If 1-m-eps and 1-m are equal numerically, 
! Gamma_inv(1+m+eps) is numerically equal to Gamma_inv(1+m), 
! already calculated as 1/(m!).
! See |eps|oo > 0.1 case for data recycling of other values 
! or for 1-m-eps and 1-m different numerically.
!
!----------------------------------------------------------------------
! Variables
! ---------
! a,b,c,one_minus_z: a,b,c and 1-z parameters and arguments 
! of the 2F1(a,b,c,z) function.
! m,eps: closest integer to c-a-b, with Re(c-a-b) >= 0 
! and eps = c-a-b-m
! Gamma_c,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm, Gamma_inv_eps_pb_pm: 
! recycled values of Gamma(c), Gamma_inv(1-eps), 
! Gamma_inv(a+m+eps) and Gamma_inv(b+m+eps).
! inf_norm_eps,phase,a_pm,b_pm,one_meps,Pi_eps,Pi_eps_pm: 
! |eps|oo,(-1)^m,a+m,b+m,1-eps,pi.eps,pi.(eps+m)
! Gamma_inv_one_meps_mm,Gamma_inv_eps_pm_p1: 
! Gamma_inv(1-m-eps) and Gamma_inv(1+m+eps) 
! calculated with the recycling scheme.
! prod1: (a)_m (b)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) 
! x Gamma_inv(b+m+eps) Gamma_inv(m+1) in |eps|oo > 0.1 case.
! prod2: (1-z)^eps Gamma_inv(a) Gamma_inv(b) Gamma_inv(1+m+eps) 
! in |eps|oo > 0.1 case.
! Gamma_inv_mp1,prod_ab: Gamma_inv(m+1) calculated as 1/(m!) 
! and (a)_m (b)_m in |eps|oo <= 0.1 case.
! is_eps_non_zero: true if 1-m-eps and 1-m are different numerically,
! false if not.
! Gamma_inv_a_pm,Gamma_inv_b_pm,z_term: Gamma_inv(a+m),Gamma_inv(b+m),
! E(eps,log(1-z))
! prod1: Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) 
! x [(G(1,-eps) Gamma_inv(m+1) + G(m+1,eps)) in |eps|oo <= 0.1 case.
! prod2: Gamma_inv(1+m+eps) (G(a+m,eps) Gamma_inv(b+m+eps) 
! + G(b+m,eps) Gamma_inv(a+m))
! prod3: E(eps,log(1-z)) Gamma_inv(a+m) Gamma_inv(b+m) 
! Gamma_inv(1+m+eps) 
! res: returned \beta_0/(1-z)^m value in all cases.
!----------------------------------------------------------------------
FUNCTION B_SUM_INIT_PS_ONE(A,B,GAMMA_C,GAMMA_INV_ONE_MEPS, &
     GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM,MZP1,M,EPS)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  INTEGER(IPR),INTENT(IN) :: M
  COMPLEX(PR),INTENT(IN) :: A,B,GAMMA_C,GAMMA_INV_ONE_MEPS, &
       GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM,MZP1,EPS
  INTEGER(IPR) :: M_M1,N,I,PHASE
  REAL(PR)     :: INF_NORM_EPS,GAMMA_INV_MP1
  COMPLEX(PR)  :: A_PM,B_SUM_INIT_PS_ONE,PI_EPS,GAMMA_INV_ONE_MEPS_MM
  COMPLEX(PR)  :: B_PM,TMP1,TMP2
  COMPLEX(PR)  :: Z_TERM,PROD1,PROD2,PROD3,ONE_MEPS,PI_EPS_PM
  COMPLEX(PR)  :: GAMMA_INV_A_PM,PROD_AB,GAMMA_INV_B_PM
  COMPLEX(PR)  :: GAMMA_INV_EPS_PM_P1
!  COMPLEX(PR)  :: GAMMA_INV,GAMMA_INV_DIFF_EPS
  !       
  INF_NORM_EPS=INF_NORM(EPS); M_M1=M-1; A_PM=A+M; B_PM=B+M
  ONE_MEPS=ONE-EPS; PI_EPS=M_PI*EPS; PI_EPS_PM = M_PI*(EPS+M)
  IF(MOD(M,2).EQ.0) THEN
     PHASE = 1
  ELSE
     PHASE = -1
  ENDIF
  GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS
  DO I=1,M
     GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS_MM*(ONE_MEPS - I)
  ENDDO
  IF(INF_NORM_EPS.GT.0.10D0) THEN
     GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
          /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
     PROD1=GAMMA_INV_ONE_MEPS*GAMMA_INV_EPS_PA_PM*GAMMA_INV_EPS_PB_PM
     DO N=0,M_M1
        PROD1=PROD1*(A+N)*(B+N)/(N+ONE)
     ENDDO
     PROD2=GAMMA_INV(A)*GAMMA_INV(B)*GAMMA_INV_EPS_PM_P1*(MZP1**EPS)
     B_SUM_INIT_PS_ONE=GAMMA_C*(PROD1-PROD2)/EPS
     RETURN
  ELSE
     GAMMA_INV_MP1=ONE;PROD_AB=ONE
     DO N=0,M_M1
        GAMMA_INV_MP1 = GAMMA_INV_MP1/(N+ONE)
        PROD_AB = PROD_AB*(A+N)*(B+N)
     ENDDO
     IF(ONE_MEPS-M.NE.1-M) THEN
        Z_TERM=EXPM1(EPS*LOG(MZP1))/EPS
        GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
             /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
     ELSE
        Z_TERM=LOG(MZP1)
        GAMMA_INV_EPS_PM_P1 = GAMMA_INV_MP1
     ENDIF
     GAMMA_INV_A_PM=GAMMA_INV(A_PM);GAMMA_INV_B_PM=GAMMA_INV(B_PM)
     TMP1=ONE; TMP2=M+1;
     PROD1 = GAMMA_INV_EPS_PA_PM*GAMMA_INV_EPS_PB_PM    &
          *(GAMMA_INV_MP1*GAMMA_INV_DIFF_EPS(TMP1,-EPS) &
          +GAMMA_INV_DIFF_EPS(TMP2,EPS))
     PROD2 = GAMMA_INV_EPS_PM_P1 &
          *(GAMMA_INV_EPS_PB_PM*GAMMA_INV_DIFF_EPS(A_PM,EPS) &
          +GAMMA_INV_A_PM*GAMMA_INV_DIFF_EPS (B_PM,EPS))
     PROD3 = GAMMA_INV_A_PM*GAMMA_INV_B_PM*GAMMA_INV_EPS_PM_P1*Z_TERM
     B_SUM_INIT_PS_ONE=GAMMA_C*PROD_AB*(PROD1-PROD2-PROD3)
     RETURN
  ENDIF
END FUNCTION B_SUM_INIT_PS_ONE
!
!----------------------------------------------------------------------
! Calculation of the first term of the B(z) power series 
! ------------------------------------------------------
! in the 1/z transformation, divided by z^{-m}
!---------------------------------------------
! In the 1/z transformation, the power series 
! B(z) = \sum_{n=0}^{+oo} \beta_n z^{-n} occurs 
! (see Comp. Phys. Comm. paper).
! The first term \beta_0, divided by z^{-m}, is calculated here. 
! m is the closest integer to Re(b-a) >= 0 and eps = b-a-m.
!
! One has to consider |eps|oo > 0.1 and |eps|oo <= 0.1, 
! where 1-m-eps and 1-m can be different or equal numerically, 
! leading to some changes in this last case.
!
! |eps|oo > 0.1
! -------------
! One has \beta_0/z^{-m} = [(a)_m (1-c+a)_m Gamma_inv(1-eps) 
! Gamma_inv(a+m+eps) Gamma_inv(c-a) Gamma_inv(m+1)
! - (-z)^{-eps} (1-c+a+eps)_m Gamma_inv(a) Gamma_inv(c-a-eps) 
! Gamma_inv(1+m+eps)].[Gamma(c)/eps], 
! stable in this regime for a direct evaluation.
!
! The values of Gamma(c), Gamma_inv(c-a) and Gamma_inv(a+m+eps) 
! were already calculated and recycled here.
! Gamma_inv(m+1) is calculated as 1/(m!). 
! Gamma_inv(1+m+eps) is calculated from Gamma_inv(1-eps) 
! as in the 1-z transformation routine.
! 
! |eps|oo <= 0.1
! --------------
! The \beta_0/z^{-m} expression is rewritten 
! so that it contains no instabilities:
! \beta_0/z^{-m} = [((1-c+a+eps)_m G(1,-eps) - P(m,eps,1-c+a) 
! Gamma_inv(1-eps)) Gamma_inv(c-a) Gamma_inv(a+m+eps) Gamma_inv(m+1)
! + (1-c+a+eps)_m [G(m+1,eps) Gamma_inv(c-a) Gamma_inv(a+m+eps) 
! - G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)]
! - (G(c-a,-eps) - E(log(-z),-eps)) Gamma_inv(m+1+eps) 
! Gamma_inv(a+m)]] (a)_m Gamma(c)
!
! Definitions and method are the same 
! as in the 1-z transformation routine, except for P(m,eps,1-c+a).
! P(m,eps,s) = [(s+eps)_m - (s)_m]/eps 
! for eps non zero and has a limit for eps -> 0.
! Let n0 be the closest integer to -Re(s) for s complex. 
! A stable formula available for eps -> 0 for P(m,eps,s) is:
! P(m,eps,s) = (s)_m E(\sum_{n=0}^{m-1} L(1/(s+n),eps),eps) 
! if n0 is not in [0:m-1],
! P(m,eps,s) = \prod_{n=0, n not equal to n0}^{m-1} (s+eps+n) 
! + (s)_m E(\sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps),eps) 
! if n0 is in [0:m-1].
! L(s,eps) is log1p(s eps)/eps if eps is not zero, 
! and L(s,0) = s.
! This expression is used in the code.
!
! Variables
! ---------
! a,b,c,z: a,b,c and z parameters 
! and arguments of the 2F1(a,b,c,z) function.
! m,eps: closest integer to b-a, with Re(b-a) >= 0 and eps = b-a-m.
! Gamma_c,Gamma_inv_cma,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm: 
! recycled values of Gamma(c), Gamma_inv(c-a), Gamma_inv(1-eps) 
! and Gamma_inv(a+m+eps).
! inf_norm_eps,phase,cma,a_mc_p1,a_mc_p1_pm,cma_eps,eps_pa_mc_p1,a_pm: 
! |eps|oo,(-1)^m,c-a,1-c+a+m,c-a-eps,1-c+a+eps,a+m
! Gamma_inv_cma_meps,one_meps,Pi_eps,Pi_eps_pm: 
! Gamma_inv(c-a-eps),1-eps,pi.eps,pi.(eps+m)
! Gamma_inv_one_meps_mm,Gamma_inv_eps_pm_p1: Gamma_inv(1-m-eps) 
! and Gamma_inv(1+m+eps) calculated with the recycling scheme.
! prod1: (a)_m (1-c+a)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) 
! x Gamma_inv(c-a) Gamma_inv(m+1) in |eps|oo > 0.1 case.
! prod2: (-z)^{-eps} (1-c+a+eps)_m Gamma_inv(a) 
! x Gamma_inv(c-a-eps) Gamma_inv(1+m+eps) in |eps|oo > 0.1 case.
! n0: closest integer to -Re(1-c+a)
! is_n0_here: true is n0 belongs to [0:m-1], false if not.
! is_eps_non_zero: true if 1-m-eps and 1-m are different numerically, 
! false if not.
! Gamma_inv_mp1,prod_a,prod_a_mc_p1: 
! Gamma_inv(m+1) calculated as 1/(m!), 
! (a)_m and (1-c+a)_m in |eps|oo <= 0.1 case.
! prod_eps_pa_mc_p1_n0: 
! \prod_{n=0, n not equal to n0}^{m-1} (1-c+a+eps+n) 
! if n0 belongs to [0:m-1], 0.0 if not, in |eps|oo <= 0.1 case.
! prod_eps_pa_mc_p1: (1-c+a+eps)_m in |eps|oo <= 0.1 case.
! sum: \sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps) if 1-m-eps 
! and 1-m are different numerically, 
! \sum_{n=0, n not equal to n0}^{m-1} 1/(s+n) if not.
! a_pn,a_mc_p1_pn,eps_pa_mc_p1_pn: a+n,1-c+a+n,1-c+a+eps+n values 
! used in (a)_m, (1-c+a)_m and (1-c+a+eps)_m evaluations.
! sum_term,prod_diff_eps,z_term: 
! E(\sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps),eps), 
! P(m,eps,1-c+a), -E(-eps,log(-z))
! Gamma_inv_a_pm,Gamma_prod1: Gamma_inv(a+m), 
! Gamma_inv(c-a).Gamma_inv(a+m+eps)
! prod1: ((1-c+a+eps)_m G(1,-eps) 
! - P(m,eps,1-c+a) Gamma_inv(1-eps)) Gamma_inv(c-a) 
! x Gamma_inv(a+m+eps) Gamma_inv(m+1)
! prod_2a: Gamma_inv(c-a).Gamma_inv(a+m+eps).G(m+1,eps)
! prod_2b: G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)
! prod_2c: (G(c-a,-eps) 
! - E(log(-z),-eps)) Gamma_inv(m+1+eps) Gamma_inv(a+m)
! prod2: (1-c+a+eps)_m [G(m+1,eps) Gamma_inv(c-a) Gamma_inv(a+m+eps) 
! - G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)] 
! - (G(c-a,-eps) - E(log(-z),-eps)) 
! x Gamma_inv(m+1+eps) Gamma_inv(a+m)]]
! res: returned \beta_0/z^{-m} value in all cases.
!----------------------------------------------------------------------
FUNCTION B_SUM_INIT_PS_INFINITY(A,C,GAMMA_C,GAMMA_INV_CMA, &
     GAMMA_INV_ONE_MEPS,GAMMA_INV_EPS_PA_PM,Z,M,EPS)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  INTEGER(IPR),INTENT(IN) :: M
  COMPLEX(PR),INTENT(IN) :: A,C,GAMMA_C,GAMMA_INV_CMA,Z,EPS
  COMPLEX(PR),INTENT(IN) :: GAMMA_INV_ONE_MEPS,GAMMA_INV_EPS_PA_PM
  INTEGER(IPR) :: M_M1,I,N,N0,PHASE
  LOGICAL      :: IS_N0_HERE,IS_EPS_NON_ZERO
  REAL(PR)     :: INF_NORM_EPS,NP1,GAMMA_INV_MP1
  COMPLEX(PR)  :: B_SUM_INIT_PS_INFINITY,TMP1
  COMPLEX(PR)  :: CMA,A_MC_P1,A_MC_P1_PM,CMA_MEPS,EPS_PA_MC_P1,A_PM
  COMPLEX(PR)  :: GAMMA_INV_EPS_PM_P1,GAMMA_INV_CMA_MEPS,PI_EPS
  COMPLEX(PR)  :: PROD1,PROD2,A_PN,A_MC_P1_PN,ONE_MEPS
  COMPLEX(PR)  :: PROD_A,PROD_A_MC_P1,PROD_EPS_PA_MC_P1_N0,PI_EPS_PM
  COMPLEX(PR)  :: PROD_EPS_PA_MC_P1,SUM_N0,Z_TERM,SUM_TERM
  COMPLEX(PR)  :: PROD_DIFF_EPS,GAMMA_INV_A_PM,GAMMA_PROD1
  COMPLEX(PR)  :: PROD_2A,PROD_2B,PROD_2C
  COMPLEX(PR)  :: EPS_PA_MC_P1_PN,GAMMA_INV_ONE_MEPS_MM
!  COMPLEX(PR)  :: GAMMA_INV,GAMMA_INV_DIFF_EPS
  !
  INF_NORM_EPS=INF_NORM(EPS); CMA=C-A; A_MC_P1=A-C+ONE
  A_MC_P1_PM=A_MC_P1+M; CMA_MEPS=CMA-EPS; EPS_PA_MC_P1=EPS+A_MC_P1
  A_PM=A+M; M_M1=M-1; ONE_MEPS=ONE-EPS; PI_EPS=M_PI*EPS
  PI_EPS_PM=M_PI*(EPS+M); GAMMA_INV_CMA_MEPS=GAMMA_INV(CMA_MEPS)
  IF(MOD(M,2).EQ.0) THEN
     PHASE = 1
  ELSE
     PHASE = -1
  ENDIF
  GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS
  DO I=1,M
     GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS_MM*(ONE_MEPS - I)
  ENDDO
  IF(INF_NORM_EPS.GT.0.1D0) THEN
     GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
          /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
     PROD1 = GAMMA_INV_CMA*GAMMA_INV_EPS_PA_PM*GAMMA_INV_ONE_MEPS
     PROD2 = GAMMA_INV(A)*GAMMA_INV_CMA_MEPS*GAMMA_INV_EPS_PM_P1 &
          *((-Z)**(-EPS))
     DO N=0,M_M1
        A_PN=A+N; A_MC_P1_PN=A_MC_P1+N
        EPS_PA_MC_P1_PN=EPS+A_MC_P1_PN;NP1=N+ONE
        PROD1 = PROD1*A_PN*A_MC_P1_PN/NP1
        PROD2 = PROD2*EPS_PA_MC_P1_PN
     ENDDO
     B_SUM_INIT_PS_INFINITY = GAMMA_C*(PROD1-PROD2)/EPS
     RETURN
  ELSE
     N0=-NINT(REAL(A_MC_P1,PR))
     IS_EPS_NON_ZERO=ONE_MEPS-M.NE.1-M
     IS_N0_HERE=(N0.GE.0).AND.(N0.LT.M)     
     GAMMA_INV_MP1=ONE; PROD_A=ONE; PROD_A_MC_P1=ONE
     PROD_EPS_PA_MC_P1=ONE; SUM_N0=ZERO
     IF(IS_N0_HERE) THEN
        PROD_EPS_PA_MC_P1_N0 = ONE
     ELSE
        PROD_EPS_PA_MC_P1_N0 = ZERO
     ENDIF
     DO N=0,M_M1
        A_PN=A+N; A_MC_P1_PN=A_MC_P1+N
        EPS_PA_MC_P1_PN=EPS+A_MC_P1_PN; NP1=N+ONE
        PROD_A = PROD_A*A_PN
        PROD_A_MC_P1 = PROD_A_MC_P1*A_MC_P1_PN
        PROD_EPS_PA_MC_P1 = PROD_EPS_PA_MC_P1*EPS_PA_MC_P1_PN
        GAMMA_INV_MP1 = GAMMA_INV_MP1/NP1
        IF(N.NE.N0) THEN
           IF(IS_N0_HERE) THEN
              PROD_EPS_PA_MC_P1_N0=PROD_EPS_PA_MC_P1_N0 &
                   *EPS_PA_MC_P1_PN
           ENDIF
           IF(IS_EPS_NON_ZERO) THEN
              SUM_N0 = SUM_N0 + LOG1P(EPS/A_MC_P1_PN)
           ELSE
              SUM_N0 = SUM_N0 + ONE/A_MC_P1_PN
           ENDIF
        ENDIF
     ENDDO
     IF(IS_EPS_NON_ZERO) THEN
        GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
             /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
        SUM_TERM = EXPM1(SUM_N0)/EPS
        Z_TERM = EXPM1(-EPS*LOG(-Z))/EPS
     ELSE
        GAMMA_INV_EPS_PM_P1 = GAMMA_INV_MP1
        SUM_TERM = SUM_N0
        Z_TERM = -LOG(-Z)
     ENDIF
     PROD_DIFF_EPS = PROD_EPS_PA_MC_P1_N0 + PROD_A_MC_P1*SUM_TERM
     GAMMA_INV_A_PM = GAMMA_INV(A_PM)
     GAMMA_PROD1=GAMMA_INV_CMA*GAMMA_INV_EPS_PA_PM
     TMP1=ONE
     PROD1 = GAMMA_PROD1*GAMMA_INV_MP1*(GAMMA_INV_DIFF_EPS(TMP1,-EPS) &
          *PROD_EPS_PA_MC_P1 - GAMMA_INV_ONE_MEPS*PROD_DIFF_EPS)
     TMP1=M+1
     PROD_2A = GAMMA_PROD1*GAMMA_INV_DIFF_EPS(TMP1,EPS) 
     PROD_2B = GAMMA_INV_CMA*GAMMA_INV_EPS_PM_P1  &
          *GAMMA_INV_DIFF_EPS(A_PM,EPS)
     PROD_2C = GAMMA_INV_EPS_PM_P1*GAMMA_INV_A_PM &
          *(GAMMA_INV_DIFF_EPS(CMA,-EPS) + GAMMA_INV_CMA_MEPS*Z_TERM)
     PROD2 = PROD_EPS_PA_MC_P1*(PROD_2A - PROD_2B - PROD_2C)
     B_SUM_INIT_PS_INFINITY = GAMMA_C*PROD_A*(PROD1+PROD2)
     RETURN
  ENDIF
END FUNCTION B_SUM_INIT_PS_INFINITY
!
!----------------------------------------------------------------------
! Calculation of the derivative of the polynomial P(X) 
! ----------------------------------------------------
! testing power series convergence
! --------------------------------
! P(X) = |z(a+X)(b+X)|^2 - |(c+X)(X+1)|^2 
!      = \sum_{i=0}^{4} c[i] X^{i}, for |z| < 1.
! It is positive when the power series term modulus increases 
! and negative when it decreases, 
! so that its derivative provides information on its convergence 
! (see Comp. Phys. Comm. paper).
! Its derivative components cv_poly_der_tab[i] = (i+1) c[i+1] 
! for i in [0:3] 
! so that P'(X) = \sum_{i=0}^{3} cv_poly_der_tab[i] X^{i} 
! are calculated.
!
! Variables:
! ----------
! a,b,c,z: a,b,c and z parameters and arguments 
! of the 2F1(a,b,c,z) function.
! cv_poly_der_tab[3]: table of four doubles 
! containing the P'(X) components.
! mod_a2,mod_b2,mod_c2,mod_z2,R_a,Re_b,Re_c: |a|^2, |b|^2, |c|^2, 
! |z|^2, Re(a), Re(b), Re(c), with which P(X) can be expressed.
!----------------------------------------------------------------------
SUBROUTINE CV_POLY_DER_TAB_CALC(A,B,C,Z,CV_POLY_DER_TAB)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  REAL(PR),INTENT(OUT) :: CV_POLY_DER_TAB(0:3)
  REAL(PR)     :: MOD_A2,MOD_B2,MOD_C2,MOD_Z2
  REAL(PR)     :: RE_A,RE_B,RE_C,IM_A,IM_B,IM_C,RE_Z,IM_Z
  !
  RE_A=REAL(A,PR); IM_A=AIMAG(A); MOD_A2=RE_A*RE_A+IM_A*IM_A
  RE_B=REAL(B,PR); IM_B=AIMAG(B); MOD_B2=RE_B*RE_B+IM_B*IM_B
  RE_C=REAL(C,PR); IM_C=AIMAG(C); MOD_C2=RE_C*RE_C+IM_C*IM_C
  RE_Z=REAL(Z,PR); IM_Z=AIMAG(Z); MOD_Z2=RE_Z*RE_Z+IM_Z*IM_Z
  CV_POLY_DER_TAB(0)=TWO*((RE_A*MOD_B2+RE_B*MOD_A2)*MOD_Z2-RE_C-MOD_C2)
  CV_POLY_DER_TAB(1)=TWO*((MOD_A2+MOD_B2+4.0D0*RE_A*RE_B)*MOD_Z2 &
       -ONE-4.0D0*RE_C-MOD_C2)
  CV_POLY_DER_TAB(2)=6.0D0*((RE_A+RE_B)*MOD_Z2-RE_C-ONE)
  CV_POLY_DER_TAB(3)=4.0D0*(MOD_Z2-ONE)
END SUBROUTINE CV_POLY_DER_TAB_CALC
!
!----------------------------------------------------------------------
! Calculation of the derivative of the polynomial P(X) 
! ----------------------------------------------------
! testing power series convergence at one x value
! -----------------------------------------------
! P'(x) is calculated for a real x. 
! See P'(X) components calculation routine for definitions.
!----------------------------------------------------------------------
FUNCTION CV_POLY_DER_CALC(CV_POLY_DER_TAB,X)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  REAL(PR),INTENT(IN) :: X
  REAL(PR),INTENT(IN) :: CV_POLY_DER_TAB(0:3)
  REAL(PR) :: CV_POLY_DER_CALC
  !
  CV_POLY_DER_CALC=CV_POLY_DER_TAB(0)+X*(CV_POLY_DER_TAB(1) &
       +X*(CV_POLY_DER_TAB(2)+X*CV_POLY_DER_TAB(3)))
  RETURN
END FUNCTION CV_POLY_DER_CALC
!
!----------------------------------------------------------------------
! Calculation of an integer after which false convergence cannot occur
! --------------------------------------------------------------------
! See cv_poly_der_tab_calc routine for definitions.
! If P'(x) < 0 and P''(x) < 0 for x > xc, it will be so for all x > xc 
! as P(x) -> -oo for x -> +oo 
! and P(x) can have at most one maximum for x > xc. 
! It means that the 2F1 power series term modulus will increase 
! or decrease to 0 for n > nc, 
! with nc the smallest positive integer larger than xc.
!
! If P'(X) = C0 + C1.X + C2.X^2 + C3.X^3, 
! the discriminant of P''(X) is Delta = C2^2 - 3 C1 C3.
!
! If Delta > 0, P''(X) has two different real roots 
! and its largest root is -(C2 + sqrt(Delta))/(3 C3), 
! because C3 = 4(|z|^2 - 1) < 0.
! One can take xc = -(C2 + sqrt(Delta))/(3 C3) 
! and one returns its associated nc integer.
!
! If Delta <= 0, P''(X) has at most one real root, 
! so that P'(X) has only one root and then P(X) only one maximum.
! In this case, one can choose xc = nc = 0, which is returned.
!
! Variables
! ---------
! cv_poly_der_tab: table of four doubles 
! containing the P'(X) coefficients
! C1,C2,three_C3: cv_poly_der_tab[1], cv_poly_der_tab[2] 
! and 3.0*cv_poly_der_tab[3], so that P''(X) = C1 + 2.C2.x + three_C3.x^2
! Delta: discriminant of P''(X), equal to C2^2 - 3 C1 C3.
! largest_root: if Delta > 0, 
! P''(X) largest real root equal to -(C2 + sqrt(Delta))/(3 C3).
!----------------------------------------------------------------------
FUNCTION MIN_N_CALC(CV_POLY_DER_TAB)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  REAL(PR),INTENT(IN) :: CV_POLY_DER_TAB(0:3)
  INTEGER(IPR) :: MIN_N_CALC
  REAL(PR)     :: C1,C2,THREE_C3,DELTA,LARGEST_ROOT
  !
  C1=CV_POLY_DER_TAB(1); C2=CV_POLY_DER_TAB(2)
  THREE_C3=3.0D0*CV_POLY_DER_TAB(3); DELTA = C2*C2 - THREE_C3*C1
  IF(DELTA.LE.ZERO) THEN
     MIN_N_CALC = 0
     RETURN
  ELSE
     LARGEST_ROOT = -(C2 + SQRT (DELTA))/THREE_C3
     MIN_N_CALC = MAX(CEILING(LARGEST_ROOT),0)
     RETURN
  ENDIF
END FUNCTION MIN_N_CALC
!
!----------------------------------------------------------------------
! Calculation of the 2F1 power series converging for |z| < 1
! ----------------------------------------------------------
! One has 2F1(a,b,c,z) 
! = \sum_{n = 0}^{+oo} (a)_n (b)_n / ((c)_n n!) z^n,
! so that 2F1(a,b,c,z) = \sum_{n = 0}^{+oo} t[n] z^n, 
! with t[0] = 1 and t[n+1] = (a+n)(b+n)/((c+n)(n+1)) t[n] for n >= 0.
! If a or b are negative integers, 
! F(z) is a polynomial of degree -a or -b, evaluated directly.
! If not, one uses the test of convergence |t[n] z^n|oo < 1E-15 
! to truncate the series after it was checked 
! that false convergence cannot occur.
! Variables:
! ----------
! a,b,c,z: a,b,c and z parameters and arguments 
! of the 2F1(a,b,c,z) function. One must have here |z| < 1.
! term,sum: term of the 2F1 power series equal to t[n] z^n, 
! truncated sum at given n of the 2F1 power series.
! na,nb: absolute values of the closest integers to Re(a) and Re(b). 
! a = -na or b = -nb means one is in the polynomial case.
! cv_poly_der_tab: coefficients of the derivative 
! of the polynomial P(X) = |z(a+X)(b+X)|^2 - |(c+X)(X+1)|^2
! min_n: smallest integer after which false convergence cannot occur. 
! It is calculated in min_n_calc.
! possible_false_cv: always true if n < min_n. 
! If n >= min_n, it is true if P'(n) > 0. 
! If n >= min_n and P'(n) < 0, 
! it becomes false and remains as such for the rest of the calculation. 
! One can then check if |t[n] z^n|oo < 1E-15 to truncate the series.
!----------------------------------------------------------------------
FUNCTION HYP_PS_ZERO(A,B,C,Z)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  INTEGER(IPR) :: N,NA,NB,MIN_N
  COMPLEX(PR)  :: HYP_PS_ZERO,TERM
  LOGICAL :: POSSIBLE_FALSE_CV
  REAL(PR) :: CV_POLY_DER_TAB(0:3)
!  REAL(PR) :: CV_POLY_DER_CALC,MIN_N_CALC
  !
  NA = ABS(NINT(REAL(A,PR)))
  NB = ABS(NINT(REAL(B,PR)))
  TERM=ONE; HYP_PS_ZERO=ONE  
  IF(A.EQ.(-NA)) THEN
     DO N=0,NA-1
        TERM = TERM*Z*(A+N)*(B+N)/((N+ONE)*(C+N))
        HYP_PS_ZERO = HYP_PS_ZERO + TERM
     ENDDO
     RETURN
  ELSE IF(B.EQ.(-NB)) THEN
     DO N=0,NB-1
        TERM = TERM*Z*(A+N)*(B+N)/((N+ONE)*(C+N))
        HYP_PS_ZERO = HYP_PS_ZERO + TERM
     ENDDO
     RETURN
  ELSE
     CALL CV_POLY_DER_TAB_CALC(A,B,C,Z,CV_POLY_DER_TAB)
     POSSIBLE_FALSE_CV=.TRUE.
     MIN_N=MIN_N_CALC(CV_POLY_DER_TAB);N=0
     DO WHILE(POSSIBLE_FALSE_CV.OR.(INF_NORM(TERM).GT.EPS15))
        TERM = TERM*Z*(A+N)*(B+N)/((N+ONE)*(C+N))
        HYP_PS_ZERO = HYP_PS_ZERO + TERM
        IF(POSSIBLE_FALSE_CV.AND.(N.GT.MIN_N)) THEN
           POSSIBLE_FALSE_CV = &
                (CV_POLY_DER_CALC (CV_POLY_DER_TAB,real(N,pr)).GT.ZERO)
        ENDIF
        N=N+1 
     ENDDO
     RETURN
  ENDIF
END FUNCTION HYP_PS_ZERO
!
!----------------------------------------------------------------------
! Calculation of the 2F1 power series 
! -----------------------------------
! converging with the 1-z transformation
! --------------------------------------
! The formula for F(z) in the 1-z transformation holds:
! F(z) = (-1)^m (pi.eps)/sin (pi.eps) [A(z) + B(z)] 
! for eps not equal to zero, F(z) = (-1)^m [A(z) + B(z)] for eps = 0
! where m = |Re(c-a-b)], eps = c-a-b-m, 
! A(z) = \sum_{n=0}^{m-1} alpha[n] (1-z)^n, 
! B(z) = \sum_{n=0}^{+oo} beta[n] (1-z)^n, and:
!
! alpha[0] = [Gamma_inv(1-m-eps)/eps] Gamma_inv(a+m+eps) 
!          x Gamma_inv(b+m+eps) Gamma(c)
! [Gamma_inv(1-m-eps)/eps] is calculated in A_sum_init. 
! alpha[0] is calculated with log[Gamma] 
! if the previous expression might overflow, 
! and its imaginary part removed if a, b and c are real.
! alpha[n+1] = (a+n)(b+n)/[(n+1)(1-m-eps+n)] alpha[n], n in [0:m-2].
!
! beta[0] is defined in B_sum_init_PS_one function comments.
! gamma[0] = Gamma(c) (a)_m (b)_m (1-z)^m Gamma_inv(a+m+eps) 
!          x Gamma_inv(b+m+eps) Gamma_inv(m+1) Gamma_inv(1-eps)
!
! beta[n+1] = (a+m+n+eps)(b+m+n+eps)/[(m+n+1+eps)(n+1)] beta[n]
! + [(a+m+n)(b+m+n)/(m+n+1) - (a+m+n) - (b+m+n) - eps 
! + (a+m+n+eps)(b+m+n+eps)/(n+1)]
!             x gamma[n]/[(n+m+1+eps)(n+1+eps)], n >= 0.
! gamma[n+1] = (a+m+n)(b+m+n)/[(m+n+1)(n+1-eps)] gamma[n], n >= 0.
!
! B(z) converges <=> |1-z| < 1
! The test of convergence is |beta[n] (1-z)^n|oo < 1E-15 |beta[0]|oo
! for n large enough so that false convergence cannot occur.
!
! Variables
! ---------
! a,b,c,one_minus_z: a,b,c parameters 
! and 1-z from z argument of 2F1(a,b,c,z)
! m,phase,m_p1,eps,eps_pm,eps_pm_p1,
! a_pm,b_pm,one_meps,one_meps_pm: 
! |Re(c-a-b)], (-1)^m, m+1, c-a-b-m, 
! eps+m, eps+m+1, a+m, b+m, 1-eps, 1-eps-m
! eps_pa,eps_pb,eps_pa_pm,eps_pb_pm,Pi_eps,Gamma_c: 
! eps+a, eps+b, eps+a+m, eps+b+m, pi.eps, Gamma(c)
! Gamma_inv_eps_pa_pm,Gamma_inv_eps_pb_pm,Gamma_prod: 
! Gamma_inv(eps+a+m), Gamma_inv(eps+b+m), 
! Gamma(c).Gamma_inv(eps+a+m).Gamma_inv(eps+b+m)
! Gamma_inv_one_meps,A_first_term,A_sum,A_term: 
! Gamma_inv(1-eps), alpha[0], A(z), alpha[n] (1-z)^n
! pow_mzp1_m,B_first_term,prod_B,ratio: (1-z)^m, beta[0], 
! (a)_m (b)_m (1-z)^m, (a+n)(b+n)/(n+1) for n in [0:m-2].
! B_extra_term,B_term,B_sum,B_prec: 
! gamma[n], beta[n] (1-z)^n, B(z), 1E-15 |beta[0|oo
! cv_poly1_der_tab,cv_poly2_der_tab: P1'(X) and P2'(X) coefficients 
! of the potentials derivatives of P1(X) and P2(X) 
! defined in cv_poly_der_tab_calc with parameters 
! a1 = a, b1 = b, c1 = 1-m-eps, z1 = 1-z 
! and a2 = eps+b+m, b2 = eps+a+m,c2 = eps+m+1, z2 = 1-z.
! min_n: smallest integer after which false convergence cannot occur. 
! It is calculated in min_n_calc with both P1'(X) and P2'(X), 
! so one takes the largest integer coming from both calculations.
! possible_false_cv: always true if n < min_n. 
! If n >= min_n, it is true if P1'(n) > 0 or P2'(n) > 0. 
! If n >= min_n and P1'(n) < 0 and P2'(n) < 0, 
! it becomes false and remains as such for the rest of the calculation.
! One can then check if |beta[n] z^n|oo < 1E-15 to truncate the series.
! n,n_pm_p1,n_p1,a_pm_pn,b_pm_pn,eps_pm_p1_pn,n_p1_meps,eps_pa_pm_pn,
! eps_pb_pm_pn,eps_pm_pn: index of power series, n+m+1, n+1, 
! a+m+n, b+m+n, eps+m+n+1, n+1-eps, eps+a+m+n, eps+b+m+n, eps+m+n,
! prod1,prod2,prod3: (eps+a+m+n)(eps+b+m+n), 
! (eps+m+1+n)(n+1), (a+m+n)(b+m+n)
!----------------------------------------------------------------------
FUNCTION HYP_PS_ONE(A,B,C,MZP1)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,MZP1
  INTEGER(IPR) :: N,M,PHASE,M_M2,MIN_N,M_P1 
  REAL(PR)     :: B_PREC,N_P1,N_PM_P1
  COMPLEX(PR)  :: HYP_PS_ONE,EPS,EPS_PM,EPS_PM_P1,A_PM
  COMPLEX(PR)  :: B_PM,ONE_MEPS_MM,EPS_PA,EPS_PB,PI_EPS,GAMMA_PROD
  COMPLEX(PR)  :: EPS_PA_PM,EPS_PB_PM
  COMPLEX(PR)  :: A_SUM,A_TERM,ONE_MEPS
  COMPLEX(PR)  :: B_EXTRA_TERM,B_TERM,B_SUM,GAMMA_C,RATIO
  COMPLEX(PR)  :: A_PM_PN,B_PM_PN,EPS_PM_P1_PN,N_P1_MEPS
  COMPLEX(PR)  :: PROD1,PROD2,PROD3
  COMPLEX(PR)  :: EPS_PA_PM_PN,EPS_PB_PM_PN,EPS_PM_PN,PROD_B,POW_MZP1_M
  COMPLEX(PR)  :: GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM
  COMPLEX(PR)  :: GAMMA_INV_ONE_MEPS
!  COMPLEX(PR)  :: GAMMA_INV,A_SUM_INIT,LOG_GAMMA_CPLX,LOG_A_SUM_INIT,B_SUM_INIT_PS_ONE,MIN_N_CALC,CV_POLY_DER_CALC
  LOGICAL :: POSSIBLE_FALSE_CV
  REAL(PR) :: CV_POLY1_DER_TAB(0:3),CV_POLY2_DER_TAB(0:3)
  !
  M=NINT(REAL(C-A-B,PR)); M_M2=M-2; M_P1=M+1
  IF(MOD(M,2).EQ.0) THEN
     PHASE=1
  ELSE
     PHASE=-1
  ENDIF
  EPS=C-A-B-M; EPS_PM=EPS+M; EPS_PM_P1=EPS_PM+ONE; A_PM=A+M;B_PM=B+M
  ONE_MEPS=ONE-EPS; ONE_MEPS_MM=ONE_MEPS-M; EPS_PA=EPS+A; EPS_PB=EPS+B 
  PI_EPS=M_PI*EPS; EPS_PA_PM=EPS_PA+M; EPS_PB_PM=EPS_PB+M
  GAMMA_C=ONE/GAMMA_INV(C)
  GAMMA_INV_EPS_PA_PM=GAMMA_INV(EPS_PA_PM)
  GAMMA_INV_EPS_PB_PM=GAMMA_INV(EPS_PB_PM)
  GAMMA_PROD=GAMMA_C*GAMMA_INV_EPS_PA_PM*GAMMA_INV_EPS_PB_PM
  GAMMA_INV_ONE_MEPS=GAMMA_INV(ONE_MEPS)
  IF(M.EQ.0) THEN
     A_TERM=ZERO
  ELSE IF(INF_NORM(ONE_MEPS_MM &
       *(LOG(ONE + ABS(ONE_MEPS_MM))-ONE)).LT.300.0d0) THEN
     A_TERM=GAMMA_PROD*A_SUM_INIT(M,EPS,GAMMA_INV_ONE_MEPS)
  ELSE
     A_TERM=EXP(LOG_GAMMA_CPLX(C)-LOG_GAMMA_CPLX(EPS_PA_PM)&
          -LOG_GAMMA_CPLX(EPS_PB_PM)+LOG_A_SUM_INIT(M,EPS))
     IF((AIMAG(A).EQ.ZERO).AND.(AIMAG(B).EQ.ZERO)&
          .AND.(AIMAG(C).EQ.ZERO)) THEN
        A_TERM=REAL(A_TERM,PR)
     ENDIF
  ENDIF
  A_SUM=A_TERM
  POW_MZP1_M = MZP1**M
  B_TERM=B_SUM_INIT_PS_ONE(A,B,GAMMA_C,GAMMA_INV_ONE_MEPS, &
       GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM,MZP1,M,EPS)*POW_MZP1_M
  PROD_B=POW_MZP1_M
  DO N=0,M_M2
     RATIO=(A+N)*(B+N)/(N+ONE)
     A_TERM=A_TERM*MZP1*RATIO/(N+ONE_MEPS_MM)
     A_SUM=A_SUM+A_TERM
     PROD_B = PROD_B*RATIO
  ENDDO
  IF(M.GT.0) THEN
     PROD_B = PROD_B*(A+M-ONE)*(B+M-ONE)/real(M,pr)
  ENDIF
  B_EXTRA_TERM = PROD_B*GAMMA_PROD*GAMMA_INV_ONE_MEPS; B_SUM=B_TERM
  B_PREC=EPS15*INF_NORM(B_TERM)
  CALL CV_POLY_DER_TAB_CALC(A,B,ONE_MEPS_MM,MZP1,CV_POLY1_DER_TAB)
  CALL CV_POLY_DER_TAB_CALC(EPS_PB_PM,EPS_PA_PM,EPS_PM_P1,MZP1, &
       CV_POLY2_DER_TAB)
  MIN_N=MAX(MIN_N_CALC(CV_POLY1_DER_TAB),MIN_N_CALC(CV_POLY2_DER_TAB))
  POSSIBLE_FALSE_CV=.TRUE.; N=0
  DO WHILE(POSSIBLE_FALSE_CV.OR.(INF_NORM(B_TERM).GT.B_PREC))
     N_PM_P1=N+M_P1; N_P1=N+ONE; A_PM_PN=A_PM+N; B_PM_PN=B_PM+N
     EPS_PM_P1_PN=EPS_PM_P1+N; N_P1_MEPS=ONE_MEPS+N
     EPS_PM_PN=EPS_PM+N; EPS_PA_PM_PN=EPS_PA_PM+N 
     EPS_PB_PM_PN=EPS_PB_PM+N
     PROD1=EPS_PA_PM_PN*EPS_PB_PM_PN
     PROD2=EPS_PM_P1_PN*N_P1
     PROD3=A_PM_PN*B_PM_PN
     B_TERM = MZP1*(B_TERM*PROD1/PROD2+B_EXTRA_TERM*(PROD3/N_PM_P1 &
          -A_PM_PN-B_PM_PN-EPS+PROD1/N_P1)/(EPS_PM_P1_PN*N_P1_MEPS))
     B_SUM=B_SUM+B_TERM
     B_EXTRA_TERM=B_EXTRA_TERM*MZP1*PROD3/(N_PM_P1*N_P1_MEPS)
     IF(POSSIBLE_FALSE_CV.AND.(N.GT.MIN_N)) THEN
        POSSIBLE_FALSE_CV = &
             (CV_POLY_DER_CALC(CV_POLY1_DER_TAB,real(N,pr)).GT.ZERO).OR. &
             (CV_POLY_DER_CALC(CV_POLY2_DER_TAB,real(N,pr)).GT.ZERO)
     ENDIF
     N=N+1
  ENDDO
  IF(EPS.EQ.ZERO) THEN
     HYP_PS_ONE=PHASE*(A_SUM+B_SUM)
     RETURN
  ELSE
     HYP_PS_ONE=PHASE*(A_SUM+B_SUM)*PI_EPS/SIN(PI_EPS)
     RETURN
  ENDIF
END FUNCTION HYP_PS_ONE
!
!----------------------------------------------------------------------
! Calculation of the 2F1 power series 
! -----------------------------------
! converging with the 1/z transformation
! --------------------------------------
! The formula for F(z) in the 1/z transformation holds:
! F(z) = (-1)^m (pi.eps)/sin (pi.eps) [A(z) + B(z)] 
! for eps not equal to zero, 
! F(z) = (-1)^m [A(z) + B(z)] for eps = 0
! where m = |Re(b-a)], eps = b-a-m, 
! A(z) = \sum_{n=0}^{m-1} alpha[n] z^{-n}, 
! B(z) = \sum_{n=0}^{+oo} beta[n] z^{-n}, and:
!
! alpha[0] = [Gamma_inv(1-m-eps)/eps] Gamma_inv(c-a) 
!          x Gamma_inv(a+m+eps) Gamma(c)
! [Gamma_inv(1-m-eps)/eps] is calculated in A_sum_init. 
! alpha[0] is calculated with log[Gamma] 
! if the previous expression might overflow, 
! and its imaginary part removed if a, b and c are real.
! alpha[n+1] = (a+n)(1-c+a+n)/[(n+1)(1-m-eps+n)] alpha[n], 
! n in [0:m-2].
!
! beta[0] is defined in B_sum_init_PS_infinity function comments.
! gamma[0] = Gamma(c) (a)_m (1-c+a)_m z^{-m} Gamma_inv(a+m+eps) 
!          x Gamma_inv(c-a) Gamma_inv(m+1) Gamma_inv(1-eps)
!
! beta[n+1] = (a+m+n+eps)(1-c+a+m+n+eps)/[(m+n+1+eps)(n+1)] beta[n] 
! + [(a+m+n)(1-c+a+m+n)/(m+n+1) - (a+m+n) - (1-c+a+m+n) 
! - eps + (a+m+n+eps)(1-c+a+m+n+eps)/(n+1)]
! x gamma[n]/[(n+m+1+eps)(n+1+eps)], n >= 0.
! gamma[n+1] = (a+m+n)(b+m+n)/[(m+n+1)(n+1-eps)] gamma[n], n >= 0.
!
! B(z) converges <=> |z| > 1
! The test of convergence is |beta[n] z^{-n}|oo < 1E-15 |beta[0]|oo
! for n large enough so that false convergence cannot occur.
!
! Variables
! ---------
! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
! m,phase,m_p1,eps,a_mc_p1,one_meps,
! one_meps_pm,a_pm,a_mc_p1_pm,cma: |Re(b-a)], (-1)^m, m+1, b-a-m, 
! 1-c+a, 1-eps, 1-eps-m, a+m, 1-c+a+m, c-a
! eps_pa,eps_pm_p1,eps_pa_mc_p1_pm,Pi_eps,eps_pa_pm,eps_pm,Gamma_c: 
! eps+a, eps+m+1, eps+1-c+a+m, pi.eps, eps+a+m, eps+m, Gamma(c)
! Gamma_inv_eps_pa_pm,Gamma_inv_cma,z_inv,pow_mz_ma,
! Gamma_inv_one_meps,Gamma_prod: Gamma_inv(eps+a+m), Gamma_inv(c-a), 
! 1/z, (-z)^(-a), Gamma_inv(1-eps), 
! Gamma(c) Gamma_inv(c-a) Gamma_inv(eps+a+m)
! A_first_term,A_sum,A_term: alpha[0], A(z), alpha[n] z^{-n}
! pow_z_inv_m,B_first_term,prod_B,ratio: z^{-m}, beta[0], 
! (a)_m (1-c+a)_m z^{-m}, (a+n)(1-c+a+n)/(n+1) for n in [0:m-2].
! B_extra_term,B_term,B_sum,B_prec: 
! gamma[n], beta[n] z^{-n}, B(z), 1E-15 |beta[0|oo
! cv_poly1_der_tab,cv_poly2_der_tab: P1'(X) and P2'(X) coefficients 
! of the potentials derivatives of P1(X) and P2(X) 
! defined in cv_poly_der_tab_calc 
! with parameters a1 = a, b1 = 1-c+a, c1 = 1-m-eps, z1 = 1/z 
! and a2 = b, b2 = eps+1-c+a+m,c2 = eps+m+1, z2 = 1/z.
! min_n: smallest integer after which false convergence cannot occur. 
!        It is calculated in min_n_calc with both P1'(X) and P2'(X), 
! so one takes the largest integer coming from both calculations.
! possible_false_cv: always true if n < min_n. If n >= min_n, 
! it is true if P1'(n) > 0 or P2'(n) > 0. 
! If n >= min_n and P1'(n) < 0 and P2'(n) < 0, 
! it becomes false and remains as such for the rest of the calculation. 
! One can then check if |beta[n] z^n|oo < 1E-15 to truncate the series.
! n,n_pm_p1,n_p1,a_pm_pn,a_mc_p1_pm_pn,eps_pm_p1_pn,n_p1_meps,
! eps_pa_pm_pn,eps_pa_mc_p1_pm_pn,eps_pm_pn: 
! index of power series, n+m+1, n+1, a+m+n, 1-c+a+m+n, eps+m+n+1,
! n+1-eps, eps+a+m+n, eps+1-c+a+m+n, eps+m+n,
! prod1,prod2,prod3: (eps+a+m+n)(eps+1-c+a+m+n),
! (eps+m+1+n)(n+1), (a+m+n)(1-c+a+m+n)
!----------------------------------------------------------------------
FUNCTION HYP_PS_INFINITY(A,B,C,Z)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  INTEGER(IPR) :: N,M,PHASE,M_M2,MIN_N,M_P1
  REAL(PR)     :: B_PREC,N_P1,N_PM_P1
  COMPLEX(PR)  :: POW_Z_INV_M
  COMPLEX(PR)  :: HYP_PS_INFINITY,Z_INV,RATIO
  COMPLEX(PR)  :: EPS,A_MC_P1,ONE_MEPS,ONE_MEPS_MM,A_PM,A_MC_P1_PM
  COMPLEX(PR)  :: CMA,EPS_PA,EPS_PM_P1,EPS_PA_MC_P1_PM,PI_EPS
  COMPLEX(PR)  :: EPS_PA_PM,EPS_PM,GAMMA_C,GAMMA_INV_CMA,POW_MZ_MA
  COMPLEX(PR)  :: A_SUM,A_TERM
  COMPLEX(PR)  :: GAMMA_INV_EPS_PA_PM,GAMMA_INV_ONE_MEPS
  COMPLEX(PR)  :: PROD_B,B_EXTRA_TERM,B_TERM,B_SUM,PROD1
  COMPLEX(PR)  :: A_PM_PN,A_MC_P1_PM_PN,EPS_PM_P1_PN,N_P1_MEPS
  COMPLEX(PR)  :: PROD2,PROD3,GAMMA_PROD
  COMPLEX(PR)  :: EPS_PA_PM_PN,EPS_PA_MC_P1_PM_PN,EPS_PM_PN
!  COMPLEX(PR)  :: GAMMA_INV,A_SUM_INIT,LOG_GAMMA_CPLX,B_SUM_INIT_PS_INFINITY,MIN_N_CALC,CV_POLY_DER_CALC,LOG_A_SUM_INIT
  LOGICAL :: POSSIBLE_FALSE_CV
  REAL(PR) :: CV_POLY1_DER_TAB(0:3),CV_POLY2_DER_TAB(0:3)
  !
  M=NINT(REAL(B-A,PR)); M_M2=M-2;M_P1=M+1
  IF(MOD(M,2).EQ.0) THEN
     PHASE=1
  ELSE
     PHASE=-1
  ENDIF
  EPS=B-A-M; A_MC_P1=ONE-C+A; ONE_MEPS=ONE-EPS; ONE_MEPS_MM=ONE_MEPS-M
  A_PM=A+M; A_MC_P1_PM=A_MC_P1+M; CMA=C-A; EPS_PA=EPS+A
  EPS_PM=EPS+M; EPS_PM_P1=EPS_PM+ONE; EPS_PA_MC_P1_PM=EPS+A_MC_P1_PM
  PI_EPS=M_PI*EPS; EPS_PA_PM=EPS_PA+M
  GAMMA_C=ONE/GAMMA_INV(C); GAMMA_INV_EPS_PA_PM = GAMMA_INV(EPS_PA_PM)
  GAMMA_INV_ONE_MEPS = GAMMA_INV(ONE_MEPS)
  GAMMA_INV_CMA=GAMMA_INV(CMA); Z_INV=ONE/Z;POW_MZ_MA=(-Z)**(-A)
  GAMMA_PROD=GAMMA_C*GAMMA_INV_CMA*GAMMA_INV_EPS_PA_PM
  IF(M.EQ.0) THEN
     A_TERM=ZERO
  ELSE IF(INF_NORM(ONE_MEPS_MM &
       *(LOG(ONE + ABS(ONE_MEPS_MM))-ONE)).LT.300.0d0) THEN
     A_TERM=GAMMA_PROD*A_SUM_INIT(M,EPS,GAMMA_INV_ONE_MEPS)
  ELSE
     A_TERM=EXP(LOG_GAMMA_CPLX(C)-LOG_GAMMA_CPLX(CMA)-LOG_GAMMA_CPLX(B) &
          + LOG_A_SUM_INIT(M,EPS))
     IF((AIMAG(A).EQ.ZERO).AND.(AIMAG(B).EQ.ZERO).AND.     &
          (AIMAG(C).EQ.ZERO)) THEN
        A_TERM=REAL(A_TERM,PR)
     ENDIF
  ENDIF
  A_SUM=A_TERM
  POW_Z_INV_M=Z_INV**M
  B_TERM=B_SUM_INIT_PS_INFINITY(A,C,GAMMA_C,GAMMA_INV_CMA, &
       GAMMA_INV_ONE_MEPS,GAMMA_INV_EPS_PA_PM,Z,M,EPS)*POW_Z_INV_M
  PROD_B=POW_Z_INV_M
  DO N=0,M_M2
     RATIO=(A+N)*(A_MC_P1+N)/(N+ONE)
     A_TERM = A_TERM*Z_INV*RATIO/(N+ONE_MEPS_MM)
     A_SUM = A_SUM+A_TERM
     PROD_B = PROD_B*RATIO
  ENDDO
  IF (M.GT.0) THEN
     PROD_B=PROD_B*(A+M-ONE)*(A_MC_P1+M-ONE)/real(M,pr)
  ENDIF
  B_EXTRA_TERM = PROD_B*GAMMA_PROD*GAMMA_INV_ONE_MEPS
  B_SUM=B_TERM
  B_PREC=EPS15*INF_NORM(B_TERM)
  CALL CV_POLY_DER_TAB_CALC(A,A_MC_P1,ONE_MEPS_MM,Z_INV, &
       CV_POLY1_DER_TAB)
  CALL CV_POLY_DER_TAB_CALC(B,EPS_PA_MC_P1_PM,EPS_PM_P1, &
       Z_INV,CV_POLY2_DER_TAB)
  MIN_N=MAX(MIN_N_CALC(CV_POLY1_DER_TAB),MIN_N_CALC(CV_POLY2_DER_TAB))
  POSSIBLE_FALSE_CV=.TRUE.; N=0
  DO WHILE(POSSIBLE_FALSE_CV.OR.(INF_NORM(B_TERM).GT.B_PREC))
     N_PM_P1=N+M_P1; N_P1=N+ONE; A_PM_PN=A_PM+N
     A_MC_P1_PM_PN=A_MC_P1_PM+N; EPS_PM_P1_PN=EPS_PM_P1+N
     N_P1_MEPS=N_P1-EPS; EPS_PA_PM_PN=EPS_PA_PM+N
     EPS_PA_MC_P1_PM_PN=EPS_PA_MC_P1_PM+N; EPS_PM_PN=EPS_PM+N
     PROD1=EPS_PA_PM_PN*EPS_PA_MC_P1_PM_PN; PROD2=EPS_PM_P1_PN*N_P1
     PROD3=A_PM_PN*A_MC_P1_PM_PN
     B_TERM = Z_INV*(B_TERM*PROD1/PROD2+B_EXTRA_TERM*(PROD3/N_PM_P1 &
          -A_PM_PN-A_MC_P1_PM_PN-EPS+PROD1/N_P1)                    &
          /(EPS_PM_P1_PN*N_P1_MEPS))
     B_SUM=B_SUM+B_TERM
     B_EXTRA_TERM=B_EXTRA_TERM*Z_INV*PROD3/(N_PM_P1*N_P1_MEPS)
     IF(POSSIBLE_FALSE_CV.AND.(N.GT.MIN_N)) THEN
        POSSIBLE_FALSE_CV = (CV_POLY_DER_CALC( &
             CV_POLY1_DER_TAB,real(N,pr)).GT.ZERO).OR.(&
             CV_POLY_DER_CALC(CV_POLY2_DER_TAB,real(N,pr)).GT.ZERO)
     ENDIF
     N=N+1
  ENDDO
  IF(EPS.EQ.ZERO) THEN
     HYP_PS_INFINITY=PHASE*POW_MZ_MA*(A_SUM+B_SUM)
     RETURN
  ELSE
     HYP_PS_INFINITY=PHASE*POW_MZ_MA*(A_SUM+B_SUM)*PI_EPS &
          /SIN(PI_EPS)
     RETURN
  ENDIF
END FUNCTION HYP_PS_INFINITY
!
!----------------------------------------------------------------------
! Calculation of F(z) in transformation theory missing zones 
! ----------------------------------------------------------
! of the complex plane with a Taylor series
! -----------------------------------------
! If z is close to exp(+/- i.pi/3), no transformation in 1-z, z, 
! z/(z-1) or combination of them can transform z in a complex number 
! of modulus smaller than a given Rmax < 1 .
! Rmax is a radius for which one considers power series summation 
! for |z| > Rmax is too slow to be processed. One takes Rmax = 0.9 .
! Nevertheless, for Rmax = 0.9, 
! these zones are small enough to be handled 
! with a Taylor series expansion around a point z0 close to z 
! where transformation theory can be used to calculate F(z).
! One then chooses z0 to be 0.9 z/|z| if |z| < 1, and 1.1 z/|z| 
! if |z| > 1, 
! so that hyp_PS_zero or hyp_PS_infinity can be used 
! (see comments of these functions above).
! For this z0, F(z) = \sum_{n=0}^{+oo} q[n] (z-z0)^n, with:
! q[0] = F(z0), q[1] = F'(z0) = (a b/c) 2F1(a+1,b+1,c+1,z0)
! q[n+2] = [q[n+1] (n (2 z0 - 1) - c + (a+b+c+1) z0) 
! + q[n] (a+n)(b+n)/(n+1)]/(z0(1-z0)(n+2))
! As |z-z0| < 0.1, it converges with around 15 terms, 
! so that no instability can occur for moderate a, b and c.
! Convergence is tested 
! with |q[n] (z-z0)^n|oo + |q[n+1] (z-z0)^{n+1}|oo. 
! Series is truncated when this test is smaller 
! than 1E-15 (|q[0]|oo + |q[1] (z-z0)|oo).
! No false convergence can happen here 
! as q[n] behaves smoothly for n -> +oo.
!
! Variables
! ---------
! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
! abs_z,is_abs_z_small: |z|, true if |z| < 1 and false if not.
! z0,zc_z0_ratio,z0_term1,z0_term2: 0.9 z/|z| if |z| < 1, 
! and 1.1 z/|z| if |z| > 1, (z-z0)/(z0 (1-z0)), 
! 2 z0 - 1, c - (a+b+c+1) z0
! hyp_PS_z0,dhyp_PS_z0,prec: F(z0), F'(z0) calculated with 2F1 
! as F'(z0) = (a b/c) 2F1(a+1,b+1,c+1,z0), 
! precision demanded for series truncation 
! equal to 1E-15 (|q[0]|oo + |q[1] (z-z0)|oo).
! n,an,anp1,anp2,sum: index of the series, q[n] (z-z0)^n, 
! q[n+1] (z-z0)^{n+1}, q[n+2] (z-z0)^{n+2}, 
! truncated sum of the power series.
!----------------------------------------------------------------------
FUNCTION HYP_PS_COMPLEX_PLANE_REST(A,B,C,Z)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  INTEGER(IPR) :: N
  REAL(PR)     :: ABS_Z,PREC
  COMPLEX(PR)  :: HYP_PS_COMPLEX_PLANE_REST
  COMPLEX(PR)  :: Z0,ZC,ZC_Z0_RATIO,Z0_TERM1,Z0_TERM2
  COMPLEX(PR)  :: HYP_PS_Z0,DHYP_PS_Z0,AN,ANP1,ANP2
!  COMPLEX(PR)  :: HYP_PS_ZERO,HYP_PS_INFINITY
  !
  ABS_Z=ABS(Z)
  IF(ABS_Z.LT.ONE) THEN
     Z0=0.9D0*Z/ABS_Z; ZC=Z-Z0; ZC_Z0_RATIO=ZC/(Z0*(ONE-Z0))
     Z0_TERM1=TWO*Z0 - ONE; Z0_TERM2=C-(A+B+ONE)*Z0
     HYP_PS_Z0=HYP_PS_ZERO(A,B,C,Z0)
     DHYP_PS_Z0=HYP_PS_ZERO(A+ONE,B+ONE,C+ONE,Z0)*A*B/C 
  ELSE
     Z0=1.1D0*Z/ABS_Z; ZC=Z-Z0; ZC_Z0_RATIO=ZC/(Z0*(ONE-Z0))
     Z0_TERM1=TWO*Z0 - ONE; Z0_TERM2=C-(A+B+ONE)*Z0
     HYP_PS_Z0=HYP_PS_INFINITY(A,B,C,Z0)
     DHYP_PS_Z0=HYP_PS_INFINITY(A+ONE,B+ONE,C+ONE,Z0)*A*B/C 
  ENDIF
  AN=HYP_PS_Z0;ANP1=ZC*DHYP_PS_Z0;HYP_PS_COMPLEX_PLANE_REST=AN+ANP1
  PREC=EPS15*(INF_NORM(AN)+INF_NORM(ANP1)); N=0
  DO WHILE(INF_NORM(AN).GT.PREC)
     ANP2=ZC_Z0_RATIO*(ANP1*(N*Z0_TERM1-Z0_TERM2)+AN*ZC*(A+N)*(B+N) &
          /(N+ONE))/(N+TWO)
     HYP_PS_COMPLEX_PLANE_REST = HYP_PS_COMPLEX_PLANE_REST + ANP2
     N=N+1
     AN=ANP1
     ANP1=ANP2
  ENDDO
  RETURN
END FUNCTION HYP_PS_COMPLEX_PLANE_REST

!
!----------------------------------------------------------------------
! Calculation of F(z) for arbitrary z using previous routines
! -----------------------------------------------------------
! Firstly, it is checked if a,b and c are negative integers.
! If neither a nor b is negative integer but c is, 
! F(z) is undefined so that the program stops with an error message.
! If a and c are negative integers with c < a, 
! or b and c are negative integers with b < a, 
! or c is not negative integer integer but a or b is, 
! one is in the polynomial case.
! In this case, if |z| < |z/(z-1)| or z = 1, 
! hyp_PS_zero is used directly, as then |z| <= 2 
! and no instability arises with hyp_PS_zero 
! as long the degree of the polynomial is small (<= 10 typically).
! If not, one uses the transformation 
! F(z) = (1-z)^{-a} 2F1(a,c-b,c,z/(z-1)) if a is negative integer 
! or F(z) = (1-z)^{-b} 2F1(b,c-a,c,z/(z-1)) if b is negative integer 
! along with hyp_PS_zero.
! Indeed, 2F1(a,c-b,c,X) is a polynomial if a is negative integer, 
! and so is 2F1(b,c-a,c,X) if b is negative integer, 
! so that one has here |z/(z-1)| <= 2 
! and the stability of the method is the same 
! as for the |z| < |z/(z-1)| case.
! If one is in the non-polynomial case, one checks if z >= 1. 
! If it is, one is the cut of F(z) 
! so that z is replaced by z - 10^{-307}i.
! Then, using F(z) = 2F1(b,a,c,z) 
! and F(z) = (1-z)^{c-a-b} 2F1(c-a,c-b,c,z), 
! one replaces a,b,c parameters by combinations of them 
! so that Re(b-a) >= 0 and Re(c-a-b) >= 0.
! Exchanging a and b does not change convergence properties, 
! while having Re(c-a-b) >= 0 accelerates it 
! (In hyp_PS_zero, t[n] z^n ~ z^n/(n^{c-a-b}) for n -> +oo).
! If |1-z| < 1E-5, one uses hyp_PS_one 
! as the vicinity of the singular point z = 1 is treated properly.
! After that, one compares |z| and |z/(z-1)| 
! to R in {0.5,0.6,0.7,0.8,0.9}. 
! If one of them is smaller than R, 
! one uses hyp_PS_zero without transformation
! or with the transformation F(z) = (1-z)^{-a} 2F1(a,c-b,c,z/(z-1)).
! Then, if both of them are larger than 0.9, 
! one compares |1/z|, |(z-1)/z|, |1-z| and |1/(1-z)| 
! to R in {0.5,0.6,0.7,0.8,0.9}. 
! If one of them is found smaller than R, 
! with the condition that |c-b|oo < 5 for (z-1)/z transformation, 
! |a,b,c|oo < 5 for |1-z| transformation 
! and |a,c-b,c|oo < 5 for |1/(1-z)| transformation,
! the corresponding transformation is used. 
! If none of them was smaller than 0.9, 
! one is in the missing zones of transformation theory 
! so that the Taylor series of hyp_PS_complex_plane_rest is used.
!
! Variables
! ---------
! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
! Re_a,Re_b,Re_c,na,nb,nc,is_a_neg_int,is_b_neg_int,is_c_neg_int: 
! real parts of a,b,c, closest integers to a,b,c, 
! true if a,b,c is negative integers and false if not.
! zm1,z_over_zm1,z_shift: z-1, z/(z-1), z - 10^{-307}i in case z >= 1.
! ab_condition, cab_condition: true if Re(b-a) >= 0 and false if not, 
! true if Re(c-a-b) >= 0 and false if not.
! abs_zm1,abz_z,abs_z_inv,abs_z_over_zm1,abs_zm1_inv,abs_zm1_over_z: 
! |z-1|, |z|, |1/z|, |z/(z-1)|, |1/(z-1)|, |(z-1)/z|
! are_ac_small: true if |a|oo < 5 and |c|oo < 5, false if not.
! is_cmb_small: true if |c-b|oo < 5, false if not.
! are_abc_small: true if |a|oo < 5, |b|oo < 5 and |c|oo < 5, 
! false if not.
! are_a_cmb_c_small: true if |a|oo < 5, |c-b|oo < 5 and |c|oo < 5, 
! false if not.
! R_tab,R: table of radii {0.5,0.6,0.7,0.8,0.9}, one of these radii.
! res: returned result
!----------------------------------------------------------------------
RECURSIVE FUNCTION HYP_2F1(A,B,C,Z) RESULT(RES)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE
  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  INTEGER(IPR) :: NA,NB,NC,I
  REAL(PR)     :: RE_A,RE_B,RE_C,ABS_Z,ABS_ZM1,ABS_Z_OVER_ZM1
  REAL(PR)     :: ABS_ZM1_OVER_Z,ABS_ZM1_INV,R_TABLE(1:5),R,ABS_Z_INV
  COMPLEX(PR)  :: RES,Z_SHIFT
  COMPLEX(PR)  :: Z_OVER_ZM1,ZM1
!  COMPLEX(PR)  :: HYP_PS_ZERO,HYP_PS_ONE,HYP_PS_INFINITY,HYP_PS_COMPLEX_PLANE_REST
  LOGICAL      :: IS_A_NEG_INT,IS_B_NEG_INT,IS_C_NEG_INT
  LOGICAL      :: AB_CONDITION,CAB_CONDITION,ARE_A_CMB_C_SMALL
  LOGICAL      :: IS_CMB_SMALL,ARE_AC_SMALL,ARE_ABC_SMALL
  !
  RE_A=REAL(A,PR); RE_B=REAL(B,PR); RE_C=REAL(C,PR);
  NA=NINT(RE_A); NB=NINT(RE_B); NC=NINT(RE_C);
  IS_A_NEG_INT=A.EQ.NA.AND.NA.LE.0
  IS_B_NEG_INT=B.EQ.NB.AND.NB.LE.0
  IS_C_NEG_INT=C.EQ.NC.AND.NC.LE.0
  ZM1=Z-ONE
  IF(IS_C_NEG_INT) THEN
     ABS_Z=ABS(Z); Z_OVER_ZM1 = Z/ZM1
     ABS_Z_OVER_ZM1=ABS(Z_OVER_ZM1)
     IF(IS_A_NEG_INT.AND.(NC.LT.NA)) THEN
        IF((Z.EQ.ONE).OR.(ABS_Z.LT.ABS_Z_OVER_ZM1)) THEN
           RES=HYP_PS_ZERO(A,B,C,Z)
           RETURN
        ELSE
           RES=((-ZM1)**(-A))*HYP_PS_ZERO(A,C-B,C,Z_OVER_ZM1)
           RETURN
        ENDIF
     ELSE IF(IS_B_NEG_INT.AND.(NC.LT.NB)) THEN
        IF((Z.EQ.ONE).OR.(ABS_Z.LT.ABS_Z_OVER_ZM1)) THEN
           RES=HYP_PS_ZERO(A,B,C,Z)
           RETURN
        ELSE
           RES=((-ZM1)**(-B))*HYP_PS_ZERO(B,C-A,C,Z_OVER_ZM1)
           RETURN
        ENDIF
     ELSE
        STOP '2F1 UNDEFINED'
     ENDIF
  ENDIF
  IF(IS_A_NEG_INT) THEN
     ABS_Z=ABS(Z); Z_OVER_ZM1 = Z/ZM1
     ABS_Z_OVER_ZM1=ABS(Z_OVER_ZM1)
     IF((Z.EQ.ONE).OR.(ABS_Z.LT.ABS_Z_OVER_ZM1)) THEN
        RES=HYP_PS_ZERO(A,B,C,Z)
        RETURN
     ELSE
        RES=((-ZM1)**(-A))*HYP_PS_ZERO(A,C-B,C,Z_OVER_ZM1)
        RETURN
     ENDIF
  ELSE IF(IS_B_NEG_INT) THEN
     ABS_Z=ABS(Z); Z_OVER_ZM1 = Z/ZM1
     ABS_Z_OVER_ZM1=ABS(Z_OVER_ZM1)
     IF((Z.EQ.ONE).OR.(ABS_Z.LT.ABS_Z_OVER_ZM1)) THEN
        RES=HYP_PS_ZERO(A,B,C,Z)
        RETURN
     ELSE
        RES=((-ZM1)**(-B))*HYP_PS_ZERO(B,C-A,C,Z_OVER_ZM1)
        RETURN
     ENDIF
  ENDIF
  IF((REAL(Z,PR).GE.ONE).AND.(AIMAG(Z).EQ.ZERO)) THEN
     Z_SHIFT=CMPLX(REAL(Z,PR),-1.0D-307,PR)
     RES=HYP_2F1(A,B,C,Z_SHIFT)
     RETURN
  ENDIF
  AB_CONDITION = (RE_B.GE.RE_A)
  CAB_CONDITION = (RE_C.GE.RE_A + RE_B)
  IF ((.NOT.AB_CONDITION).OR.(.NOT.CAB_CONDITION)) THEN
     IF ((.NOT.AB_CONDITION).AND.(CAB_CONDITION)) THEN
        RES=HYP_2F1(B,A,C,Z)
        RETURN
     ELSE IF((.NOT.CAB_CONDITION).AND.(AB_CONDITION)) THEN
        RES=((-ZM1)**(C-A-B))*HYP_2F1(C-B,C-A,C,Z)
        RETURN
     ELSE
        RES=((-ZM1)**(C-A-B))*HYP_2F1(C-A,C-B,C,Z)
        RETURN 
     ENDIF
  ENDIF
  ABS_ZM1=ABS(ZM1)
  IF(ABS_ZM1.LT.1D-5) THEN 
     RES=HYP_PS_ONE (A,B,C,-ZM1)
     RETURN
  ENDIF
  ABS_Z=ABS(Z); ABS_Z_OVER_ZM1=ABS_Z/ABS_ZM1; ABS_Z_INV=ONE/ABS_Z
  ABS_ZM1_OVER_Z=ONE/ABS_Z_OVER_ZM1; ABS_ZM1_INV=ONE/ABS_ZM1
  IS_CMB_SMALL = INF_NORM(C-B).LT.5.0D0; 
  ARE_AC_SMALL = (INF_NORM(A).LT.5.0D0).AND.(INF_NORM(C).LT.5.0D0)
  ARE_ABC_SMALL = ARE_AC_SMALL.AND.(INF_NORM(B).LT.5.0D0)
  ARE_A_CMB_C_SMALL = ARE_AC_SMALL.AND.IS_CMB_SMALL
  R_TABLE=(/0.5D0,0.6D0,0.7D0,0.8D0,0.9D0/)
  DO I=1,5
     R=R_TABLE(I)
     IF(ABS_Z.LE.R) THEN 
        RES=HYP_PS_ZERO (A,B,C,Z)
        RETURN
     ENDIF
     IF(IS_CMB_SMALL.AND.(ABS_Z_OVER_ZM1.LE.R)) THEN
        RES=((-ZM1)**(-A))*HYP_PS_ZERO (A,C-B,C,Z/ZM1)
        RETURN
     ENDIF
  ENDDO
  DO I=1,5
     R=R_TABLE(I)
     IF(ABS_Z_INV.LE.R) THEN 
        RES=HYP_PS_INFINITY (A,B,C,Z)
        RETURN 
     ENDIF
     IF(IS_CMB_SMALL.AND.(ABS_ZM1_OVER_Z.LE.R)) THEN 
        RES=((-ZM1)**(-A))*HYP_PS_INFINITY (A,C-B,C,Z/ZM1)
        RETURN
     ENDIF
     IF(ARE_ABC_SMALL.AND.(ABS_ZM1.LE.R)) THEN 
        RES=HYP_PS_ONE (A,B,C,-ZM1)
        RETURN
     ENDIF
     IF(ARE_A_CMB_C_SMALL.AND.(ABS_ZM1_INV.LE.R)) THEN 
        RES=((-ZM1)**(-A))*HYP_PS_ONE (A,C-B,C,-ONE/ZM1)
        RETURN
     ENDIF
  ENDDO
  RES=HYP_PS_COMPLEX_PLANE_REST (A,B,C,Z)
  RETURN
END FUNCTION HYP_2F1



!
!----------------------------------------------------------------------
! Test of 2F1 numerical accuracy 
! ------------------------------
! using hypergeometric differential equation
! ------------------------------------------
! If z = 0, F(z) = 1 so that this value is trivially tested.
! To test otherwise if the value of F(z) is accurate, 
! one uses the fact that 
! z(z-1) F''(z) + (c - (a+b+1) z) F'(z) - a b F(z) = 0.
! If z is not equal to one, a relative precision test is provided 
! by |F''(z) + [(c - (a+b+1) z) F'(z) - a b F(z)]/[z(z-1)]|oo
! /(|F(z)|oo + F'(z)|oo + |F''(z)|oo).
! If z is equal to one, one uses |(c - (a+b+1)) F'(z) - a b F(z)|oo
! /(|F(z)|oo + F'(z)|oo + 1E-307).
! F'(z) and F''(z) are calculated using equalities 
! F'(z) = (a b/c) 2F1(a+1,b+1,c+1,z) 
! and F'(z) = ((a+1)(b+1)/(c+1)) (a b/c) 2F1(a+2,b+2,c+2,z).
!
! Variables
! ---------
! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
! F,dF,d2F: F(z), F'(z) and F''(z) calculated with hyp_2F1 
! using F'(z) = (a b/c) 2F1(a+1,b+1,c+1,z) 
! and F'(z) = ((a+1)(b+1)/(c+1)) (a b/c) 2F1(a+2,b+2,c+2,z).
!----------------------------------------------------------------------

FUNCTION TEST_2F1(A,B,C,Z,F)
  !--------------------------------------------------------------------
  USE HYP_2F1_MODULE
  IMPLICIT NONE

  COMPLEX(PR),INTENT(IN) :: A,B,C,Z
  REAL(PR)    :: TEST_2F1
  COMPLEX(PR) :: F,DF,D2F  
!  COMPLEX(PR) ::HYP_2F1
  IF(Z.EQ.ZERO) THEN
     TEST_2F1=INF_NORM(F-ONE)
     RETURN
  ELSE IF(Z.EQ.ONE) THEN
     DF = HYP_2F1(A+ONE,B+ONE,C+ONE,Z)*A*B/C
     TEST_2F1=INF_NORM((C-(A+B+ONE))*DF-A*B*F) &
          /(INF_NORM (F)+INF_NORM(DF)+1D-307)
     RETURN
  ELSE
     DF = HYP_2F1(A+ONE,B+ONE,C+ONE,Z)*A*B/C
     D2F = HYP_2F1(A+TWO,B+TWO,C+TWO,Z)*A*(A+ONE)*B*(B+ONE) &
          /(C*(C+ONE))
     TEST_2F1=INF_NORM(D2F+((C-(A+B+ONE)*Z)*DF-A*B*F)/(Z*(ONE-Z))) &
          /(INF_NORM(F)+INF_NORM(DF)+INF_NORM(D2F))
     RETURN
  ENDIF
END FUNCTION TEST_2F1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                        !
!                 INTEGRAL EVALUATORS                    !
!                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!DECK DLI
! Code converted using TO_F90 by Alan Miller
! Date: 2002-02-26  Time: 16:43:57
    FUNCTION dli(x) RESULT(fn_val)
!***BEGIN PROLOGUE  DLI
!***PURPOSE  Compute the logarithmic integral.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C5
!***TYPE      REAL (dp) (ALI-S, DLI-D)
!***KEYWORDS  FNLIB, LOGARITHMIC INTEGRAL, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
! DLI(X) calculates the REAL (dp) logarithmic integral
! for REAL (dp) argument X.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DEI, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DLI
      REAL (dp), INTENT(IN)  :: x
      REAL (dp)              :: fn_val
!***FIRST EXECUTABLE STATEMENT  DLI
      IF (x <= 0.d0) CALL xermsg('SLATEC', 'DLI', 'LOG INTEGRAL UNDEFINED FOR X <= 0')
      IF (x == 1.d0) CALL xermsg('SLATEC', 'DLI', 'LOG INTEGRAL UNDEFINED FOR X = 0')
      fn_val = dei(LOG(x))
      RETURN
    END FUNCTION dli


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                        !
!                  SERIES EVALUATORS                     !
!                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!DECK DCSEVL
    FUNCTION dcsevl(x, cs, n) RESULT(fn_val)
!***BEGIN PROLOGUE  DCSEVL
!***PURPOSE  Evaluate a Chebyshev series.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      REAL (dp) (CSEVL-S, DCSEVL-D)
!***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!  Evaluate the N-term Chebyshev series CS at X.  Adapted from
!  a method presented in the paper by Broucke referenced below.
!       Input Arguments --
!  X    value at which the series is to be evaluated.
!  CS   array of N terms of a Chebyshev series.  In evaluating
!       CS, only half the first coefficient is summed.
!  N    number of terms in array CS.
!***REFERENCES  R. Broucke, Ten subroutines for the manipulation of Chebyshev
!                 series, Algorithm 446, Communications of the A.C.M. 16,
!                 (1973) pp. 254-256.
!               L. Fox and I. B. Parker, Chebyshev Polynomials in Numerical
!                 Analysis, Oxford University Press, 1968, page 56.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900329  Prologued revised extensively and code rewritten to allow
!           X to be slightly outside interval (-1,+1).  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCSEVL
	REAL (dp), INTENT(IN)  :: x
	REAL (dp), INTENT(IN)  :: cs(:)
	INTEGER, INTENT(IN)    :: n
	REAL (dp)              :: fn_val
	REAL (dp)        :: b0, b1, b2, twox
	REAL (dp), SAVE  :: onepl
	LOGICAL, SAVE    :: first = .TRUE.
	INTEGER          :: i, ni
!***FIRST EXECUTABLE STATEMENT  DCSEVL
	IF (first) onepl = 1.0_dp + 2*EPSILON(0.0_dp)
	first = .false.
	IF (n < 1) CALL xermsg('SLATEC','DCSEVL', 'NUMBER OF TERMS <= 0')
	IF (n > 1000) CALL xermsg('SLATEC', 'DCSEVL', 'NUMBER OF TERMS > 1000')
	IF (ABS(x) > onepl) CALL xermsg('SLATEC','DCSEVL', 'X OUTSIDE THE INTERVAL (-1,+1)')
	b1 = 0.0_dp
	b0 = 0.0_dp
	twox = 2.0_dp * x
	DO  i = 1, n
  		b2 = b1
  		b1 = b0
 		 ni = n + 1 - i
 		 b0 = twox * b1 - b2 + cs(ni)
	END DO
	fn_val = 0.5_dp * (b0-b2)
	RETURN
   END FUNCTION dcsevl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                        !
!              POLYNOMIALS COMPUTATION                   !
!                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!DECK INITDS
FUNCTION initds(os, nos, eta) RESULT(ival)
!***BEGIN PROLOGUE  INITDS
!***PURPOSE  Determine the number of terms needed in an orthogonal
!            polynomial series so that it meets a specified accuracy.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      REAL (dp) (INITS-S, INITDS-D)
!***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
!             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!  Initialize the orthogonal series, represented by the array OS, so that
!  INITDS is the number of terms needed to insure the error is no larger than
!  ETA.  Ordinarily, ETA will be chosen to be one-tenth machine precision.
!             Input Arguments --
!   OS     REAL (dp) array of NOS coefficients in an orthogonal series.
!   NOS    number of coefficients in OS.
!   ETA    single precision scalar containing requested accuracy of series.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891115  Modified error message.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  INITDS
REAL (dp), INTENT(IN)  :: os(:)
INTEGER, INTENT(IN)    :: nos
REAL (dp), INTENT(IN)  :: eta
INTEGER                :: ival
INTEGER    :: i, ii
REAL (dp)  :: ERR
!***FIRST EXECUTABLE STATEMENT  INITDS
IF (nos < 1) CALL xermsg('SLATEC', 'INITDS', 'Number of coefficients < 1')
ERR = 0.0_dp
DO  ii = 1, nos
  i = nos + 1 - ii
  ERR = ERR + ABS(os(i))
  IF (ERR > eta) GO TO 20
END DO
20 IF (i == nos) CALL xermsg('SLATEC', 'INITDS',  &
    'Chebyshev series too short for specified accuracy')
ival = i
RETURN
END FUNCTION initds


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                        !
!           DOUBLONS, TO DELETE IF NOT USED              !
!                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!DECK DE1
FUNCTION de1(x) RESULT(fn_val)
!***BEGIN PROLOGUE  DE1
!***PURPOSE  Compute the exponential integral E1(X).
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C5
!***TYPE      REAL (dp) (E1-S, DE1-D)
!***KEYWORDS  E1 FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
! DE1 calculates the REAL (dp) exponential integral, E1(X), for positive
! REAL (dp) argument X and the Cauchy principal value for negative X.
! If principal values are used everywhere, then, for all X,
!    E1(X) = -Ei(-X)
! or
!    Ei(X) = -E1(-X).
! Series for AE10       on the interval -3.12500E-02 to  0.
!                                        with weighted error   4.62E-32
!                                         log weighted error  31.34
!                               significant figures required  29.70
!                                    decimal places required  32.18
! Series for AE11       on the interval -1.25000E-01 to -3.12500E-02
!                                        with weighted error   2.22E-32
!                                         log weighted error  31.65
!                               significant figures required  30.75
!                                    decimal places required  32.54
! Series for AE12       on the interval -2.50000E-01 to -1.25000E-01
!                                        with weighted error   5.19E-32
!                                         log weighted error  31.28
!                               significant figures required  30.82
!                                    decimal places required  32.09
! Series for E11        on the interval -4.00000E+00 to -1.00000E+00
!                                        with weighted error   8.49E-34
!                                         log weighted error  33.07
!                               significant figures required  34.13
!                                    decimal places required  33.80
! Series for E12        on the interval -1.00000E+00 to  1.00000E+00
!                                        with weighted error   8.08E-33
!                                         log weighted error  32.09
!                        approx significant figures required  30.4
!                                    decimal places required  32.79
! Series for AE13       on the interval  2.50000E-01 to  1.00000E+00
!                                        with weighted error   6.65E-32
!                                         log weighted error  31.18
!                               significant figures required  30.69
!                                    decimal places required  32.03
! Series for AE14       on the interval  0.          to  2.50000E-01
!                                        with weighted error   5.07E-32
!                                         log weighted error  31.30
!                               significant figures required  30.40
!                                    decimal places required  32.20
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891115  Modified prologue description.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  DE1
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val
REAL (dp)        :: eta, xmaxt
REAL (dp), SAVE  :: xmax
INTEGER, SAVE    :: ntae10, ntae11, ntae12, nte11, nte12, ntae13, ntae14
LOGICAL, SAVE    :: first = .TRUE.
REAL (dp), PARAMETER  :: ae10cs(50) = (/ +.3284394579616699087873844201881E-1_dp, &
 -.1669920452031362851476184343387E-1_dp, +.2845284724361346807424899853252E-3_dp,   &
 -.7563944358516206489487866938533E-5_dp, +.2798971289450859157504843180879E-6_dp,   &
 -.1357901828534531069525563926255E-7_dp, +.8343596202040469255856102904906E-9_dp,   &
 -.6370971727640248438275242988532E-10_dp, +.6007247608811861235760831561584E-11_dp, &
 -.7022876174679773590750626150088E-12_dp, +.1018302673703687693096652346883E-12_dp, &
 -.1761812903430880040406309966422E-13_dp, +.3250828614235360694244030353877E-14_dp, &
 -.5071770025505818678824872259044E-15_dp, +.1665177387043294298172486084156E-16_dp, &
 +.3166753890797514400677003536555E-16_dp, -.1588403763664141515133118343538E-16_dp, &
 +.4175513256138018833003034618484E-17_dp, -.2892347749707141906710714478852E-18_dp, &
 -.2800625903396608103506340589669E-18_dp, +.1322938639539270903707580023781E-18_dp, &
 -.1804447444177301627283887833557E-19_dp, -.7905384086522616076291644817604E-20_dp, &
 +.4435711366369570103946235838027E-20_dp, -.4264103994978120868865309206555E-21_dp, &
 -.3920101766937117541553713162048E-21_dp, +.1527378051343994266343752326971E-21_dp, &
 +.1024849527049372339310308783117E-22_dp, -.2134907874771433576262711405882E-22_dp, &
 +.3239139475160028267061694700366E-23_dp, +.2142183762299889954762643168296E-23_dp, &
 -.8234609419601018414700348082312E-24_dp, -.1524652829645809479613694401140E-24_dp, &
 +.1378208282460639134668480364325E-24_dp, +.2131311202833947879523224999253E-26_dp, &
 -.2012649651526484121817466763127E-25_dp, +.1995535662263358016106311782673E-26_dp, &
 +.2798995808984003464948686520319E-26_dp, -.5534511845389626637640819277823E-27_dp, &
 -.3884995396159968861682544026146E-27_dp, +.1121304434507359382850680354679E-27_dp, &
 +.5566568152423740948256563833514E-28_dp, -.2045482929810499700448533938176E-28_dp, &
 -.8453813992712336233411457493674E-29_dp, +.3565758433431291562816111116287E-29_dp, &
 +.1383653872125634705539949098871E-29_dp, -.6062167864451372436584533764778E-30_dp, &
 -.2447198043989313267437655119189E-30_dp, +.1006850640933998348011548180480E-30_dp, &
 +.4623685555014869015664341461674E-31_dp /)
REAL (dp), PARAMETER  :: ae11cs(60) = (/ +.20263150647078889499401236517381_dp,  &
 -.73655140991203130439536898728034E-1_dp, +.63909349118361915862753283840020E-2_dp,   &
 -.60797252705247911780653153363999E-3_dp, -.73706498620176629330681411493484E-4_dp,   &
 +.48732857449450183453464992488076E-4_dp, -.23837064840448290766588489460235E-5_dp,   &
 -.30518612628561521027027332246121E-5_dp, +.17050331572564559009688032992907E-6_dp,   &
 +.23834204527487747258601598136403E-6_dp, +.10781772556163166562596872364020E-7_dp,   &
 -.17955692847399102653642691446599E-7_dp, -.41284072341950457727912394640436E-8_dp,   &
 +.68622148588631968618346844526664E-9_dp, +.53130183120506356147602009675961E-9_dp,   &
 +.78796880261490694831305022893515E-10_dp, -.26261762329356522290341675271232E-10_dp, &
 -.15483687636308261963125756294100E-10_dp, -.25818962377261390492802405122591E-11_dp, &
 +.59542879191591072658903529959352E-12_dp, +.46451400387681525833784919321405E-12_dp, &
 +.11557855023255861496288006203731E-12_dp, -.10475236870835799012317547189670E-14_dp, &
 -.11896653502709004368104489260929E-13_dp, -.47749077490261778752643019349950E-14_dp, &
 -.81077649615772777976249734754135E-15_dp, +.13435569250031554199376987998178E-15_dp, &
 +.14134530022913106260248873881287E-15_dp, +.49451592573953173115520663232883E-16_dp, &
 +.79884048480080665648858587399367E-17_dp, -.14008632188089809829248711935393E-17_dp, &
 -.14814246958417372107722804001680E-17_dp, -.55826173646025601904010693937113E-18_dp, &
 -.11442074542191647264783072544598E-18_dp, +.25371823879566853500524018479923E-20_dp, &
 +.13205328154805359813278863389097E-19_dp, +.62930261081586809166287426789485E-20_dp, &
 +.17688270424882713734999261332548E-20_dp, +.23266187985146045209674296887432E-21_dp, &
 -.67803060811125233043773831844113E-22_dp, -.59440876959676373802874150531891E-22_dp, &
 -.23618214531184415968532592503466E-22_dp, -.60214499724601478214168478744576E-23_dp, &
 -.65517906474348299071370444144639E-24_dp, +.29388755297497724587042038699349E-24_dp, &
 +.22601606200642115173215728758510E-24_dp, +.89534369245958628745091206873087E-25_dp, &
 +.24015923471098457555772067457706E-25_dp, +.34118376888907172955666423043413E-26_dp, &
 -.71617071694630342052355013345279E-27_dp, -.75620390659281725157928651980799E-27_dp, &
 -.33774612157467324637952920780800E-27_dp, -.10479325703300941711526430332245E-27_dp, &
 -.21654550252170342240854880201386E-28_dp, -.75297125745288269994689298432000E-30_dp, &
 +.19103179392798935768638084000426E-29_dp, +.11492104966530338547790728833706E-29_dp, &
 +.43896970582661751514410359193600E-30_dp, +.12320883239205686471647157725866E-30_dp, &
 +.22220174457553175317538581162666E-31_dp /)
REAL (dp), PARAMETER  :: ae12cs(41) = (/ +.63629589796747038767129887806803_dp,  &
 -.13081168675067634385812671121135E+0_dp, -.84367410213053930014487662129752E-2_dp,   &
 +.26568491531006685413029428068906E-2_dp, +.32822721781658133778792170142517E-3_dp,   &
 -.23783447771430248269579807851050E-4_dp, -.11439804308100055514447076797047E-4_dp,   &
 -.14405943433238338455239717699323E-5_dp, +.52415956651148829963772818061664E-8_dp,   &
 +.38407306407844323480979203059716E-7_dp, +.85880244860267195879660515759344E-8_dp,   &
 +.10219226625855003286339969553911E-8_dp, +.21749132323289724542821339805992E-10_dp,  &
 -.22090238142623144809523503811741E-10_dp, -.63457533544928753294383622208801E-11_dp, &
 -.10837746566857661115340539732919E-11_dp, -.11909822872222586730262200440277E-12_dp, &
 -.28438682389265590299508766008661E-14_dp, +.25080327026686769668587195487546E-14_dp, &
 +.78729641528559842431597726421265E-15_dp, +.15475066347785217148484334637329E-15_dp, &
 +.22575322831665075055272608197290E-16_dp, +.22233352867266608760281380836693E-17_dp, &
 +.16967819563544153513464194662399E-19_dp, -.57608316255947682105310087304533E-19_dp, &
 -.17591235774646878055625369408853E-19_dp, -.36286056375103174394755328682666E-20_dp, &
 -.59235569797328991652558143488000E-21_dp, -.76030380926310191114429136895999E-22_dp, &
 -.62547843521711763842641428479999E-23_dp, +.25483360759307648606037606400000E-24_dp, &
 +.25598615731739857020168874666666E-24_dp, +.71376239357899318800207052800000E-25_dp, &
 +.14703759939567568181578956800000E-25_dp, +.25105524765386733555198634666666E-26_dp, &
 +.35886666387790890886583637333333E-27_dp, +.39886035156771301763317759999999E-28_dp, &
 +.21763676947356220478805333333333E-29_dp, -.46146998487618942367607466666666E-30_dp, &
 -.20713517877481987707153066666666E-30_dp, -.51890378563534371596970666666666E-31_dp /)
REAL (dp), PARAMETER  :: e11cs(29) = (/ -.16113461655571494025720663927566180E+2_dp,      &
 +.77940727787426802769272245891741497E+1_dp, -.19554058188631419507127283812814491E+1_dp,   &
 +.37337293866277945611517190865690209E+0_dp, -.56925031910929019385263892220051166E-1_dp,   &
 +.72110777696600918537847724812635813E-2_dp, -.78104901449841593997715184089064148E-3_dp,   &
 +.73880933562621681878974881366177858E-4_dp, -.62028618758082045134358133607909712E-5_dp,   &
 +.46816002303176735524405823868362657E-6_dp, -.32092888533298649524072553027228719E-7_dp,   &
 +.20151997487404533394826262213019548E-8_dp, -.11673686816697793105356271695015419E-9_dp,   &
 +.62762706672039943397788748379615573E-11_dp, -.31481541672275441045246781802393600E-12_dp, &
 +.14799041744493474210894472251733333E-13_dp, -.65457091583979673774263401588053333E-15_dp, &
 +.27336872223137291142508012748799999E-16_dp, -.10813524349754406876721727624533333E-17_dp, &
 +.40628328040434303295300348586666666E-19_dp, -.14535539358960455858914372266666666E-20_dp, &
 +.49632746181648636830198442666666666E-22_dp, -.16208612696636044604866560000000000E-23_dp, &
 +.50721448038607422226431999999999999E-25_dp, -.15235811133372207813973333333333333E-26_dp, &
 +.44001511256103618696533333333333333E-28_dp, -.12236141945416231594666666666666666E-29_dp, &
 +.32809216661066001066666666666666666E-31_dp, -.84933452268306432000000000000000000E-33_dp /)
REAL (dp), PARAMETER  :: e12cs(25) = (/ -.3739021479220279511668698204827E-1_dp,      &
 +.4272398606220957726049179176528E-1_dp, -.130318207984970054415392055219726_dp,     &
 +.144191240246988907341095893982137E-1_dp, -.134617078051068022116121527983553E-2_dp,   &
 +.107310292530637799976115850970073E-3_dp, -.742999951611943649610283062223163E-5_dp,   &
 +.453773256907537139386383211511827E-6_dp, -.247641721139060131846547423802912E-7_dp,   &
 +.122076581374590953700228167846102E-8_dp, -.548514148064092393821357398028261E-10_dp,  &
 +.226362142130078799293688162377002E-11_dp, -.863589727169800979404172916282240E-13_dp, &
 +.306291553669332997581032894881279E-14_dp, -.101485718855944147557128906734933E-15_dp, &
 +.315482174034069877546855328426666E-17_dp, -.923604240769240954484015923200000E-19_dp, &
 +.255504267970814002440435029333333E-20_dp, -.669912805684566847217882453333333E-22_dp, &
 +.166925405435387319431987199999999E-23_dp, -.396254925184379641856000000000000E-25_dp, &
 +.898135896598511332010666666666666E-27_dp, -.194763366993016433322666666666666E-28_dp, &
 +.404836019024630033066666666666666E-30_dp, -.807981567699845120000000000000000E-32_dp /)
REAL (dp), PARAMETER  :: ae13cs(50) = (/ -.60577324664060345999319382737747_dp,  &
 -.11253524348366090030649768852718E+0_dp, +.13432266247902779492487859329414E-1_dp,   &
 -.19268451873811457249246838991303E-2_dp, +.30911833772060318335586737475368E-3_dp,   &
 -.53564132129618418776393559795147E-4_dp, +.98278128802474923952491882717237E-5_dp,   &
 -.18853689849165182826902891938910E-5_dp, +.37494319356894735406964042190531E-6_dp,   &
 -.76823455870552639273733465680556E-7_dp, +.16143270567198777552956300060868E-7_dp,   &
 -.34668022114907354566309060226027E-8_dp, +.75875420919036277572889747054114E-9_dp,   &
 -.16886433329881412573514526636703E-9_dp, +.38145706749552265682804250927272E-10_dp,  &
 -.87330266324446292706851718272334E-11_dp, +.20236728645867960961794311064330E-11_dp, &
 -.47413283039555834655210340820160E-12_dp, +.11221172048389864324731799928920E-12_dp, &
 -.26804225434840309912826809093395E-13_dp, +.64578514417716530343580369067212E-14_dp, &
 -.15682760501666478830305702849194E-14_dp, +.38367865399315404861821516441408E-15_dp, &
 -.94517173027579130478871048932556E-16_dp, +.23434812288949573293896666439133E-16_dp, &
 -.58458661580214714576123194419882E-17_dp, +.14666229867947778605873617419195E-17_dp, &
 -.36993923476444472706592538274474E-18_dp, +.93790159936721242136014291817813E-19_dp, &
 -.23893673221937873136308224087381E-19_dp, +.61150624629497608051934223837866E-20_dp, &
 -.15718585327554025507719853288106E-20_dp, +.40572387285585397769519294491306E-21_dp, &
 -.10514026554738034990566367122773E-21_dp, +.27349664930638667785806003131733E-22_dp, &
 -.71401604080205796099355574271999E-23_dp, +.18705552432235079986756924211199E-23_dp, &
 -.49167468166870480520478020949333E-24_dp, +.12964988119684031730916087125333E-24_dp, &
 -.34292515688362864461623940437333E-25_dp, +.90972241643887034329104820906666E-26_dp, &
 -.24202112314316856489934847999999E-26_dp, +.64563612934639510757670475093333E-27_dp, &
 -.17269132735340541122315987626666E-27_dp, +.46308611659151500715194231466666E-28_dp, &
 -.12448703637214131241755170133333E-28_dp, +.33544574090520678532907007999999E-29_dp, &
 -.90598868521070774437543935999999E-30_dp, +.24524147051474238587273216000000E-30_dp, &
 -.66528178733552062817107967999999E-31_dp /)
REAL (dp), PARAMETER  :: ae14cs(64) = (/ -.1892918000753016825495679942820_dp, &
 -.8648117855259871489968817056824E-1_dp, +.7224101543746594747021514839184E-2_dp,   &
 -.8097559457557386197159655610181E-3_dp, +.1099913443266138867179251157002E-3_dp,   &
 -.1717332998937767371495358814487E-4_dp, +.2985627514479283322825342495003E-5_dp,   &
 -.5659649145771930056560167267155E-6_dp, +.1152680839714140019226583501663E-6_dp,   &
 -.2495030440269338228842128765065E-7_dp, +.5692324201833754367039370368140E-8_dp,   &
 -.1359957664805600338490030939176E-8_dp, +.3384662888760884590184512925859E-9_dp,   &
 -.8737853904474681952350849316580E-10_dp, +.2331588663222659718612613400470E-10_dp, &
 -.6411481049213785969753165196326E-11_dp, +.1812246980204816433384359484682E-11_dp, &
 -.5253831761558460688819403840466E-12_dp, +.1559218272591925698855028609825E-12_dp, &
 -.4729168297080398718476429369466E-13_dp, +.1463761864393243502076199493808E-13_dp, &
 -.4617388988712924102232173623604E-14_dp, +.1482710348289369323789239660371E-14_dp, &
 -.4841672496239229146973165734417E-15_dp, +.1606215575700290408116571966188E-15_dp, &
 -.5408917538957170947895023784252E-16_dp, +.1847470159346897881370231402310E-16_dp, &
 -.6395830792759094470500610425050E-17_dp, +.2242780721699759457250233276170E-17_dp, &
 -.7961369173983947552744555308646E-18_dp, +.2859308111540197459808619929272E-18_dp, &
 -.1038450244701137145900697137446E-18_dp, +.3812040607097975780866841008319E-19_dp, &
 -.1413795417717200768717562723696E-19_dp, +.5295367865182740958305442594815E-20_dp, &
 -.2002264245026825902137211131439E-20_dp, +.7640262751275196014736848610918E-21_dp, &
 -.2941119006868787883311263523362E-21_dp, +.1141823539078927193037691483586E-21_dp, &
 -.4469308475955298425247020718489E-22_dp, +.1763262410571750770630491408520E-22_dp, &
 -.7009968187925902356351518262340E-23_dp, +.2807573556558378922287757507515E-23_dp, &
 -.1132560944981086432141888891562E-23_dp, +.4600574684375017946156764233727E-24_dp, &
 -.1881448598976133459864609148108E-24_dp, +.7744916111507730845444328478037E-25_dp, &
 -.3208512760585368926702703826261E-25_dp, +.1337445542910839760619930421384E-25_dp, &
 -.5608671881802217048894771735210E-26_dp, +.2365839716528537483710069473279E-26_dp, &
 -.1003656195025305334065834526856E-26_dp, +.4281490878094161131286642556927E-27_dp, &
 -.1836345261815318199691326958250E-27_dp, +.7917798231349540000097468678144E-28_dp, &
 -.3431542358742220361025015775231E-28_dp, +.1494705493897103237475066008917E-28_dp, &
 -.6542620279865705439739042420053E-29_dp, +.2877581395199171114340487353685E-29_dp, &
 -.1271557211796024711027981200042E-29_dp, +.5644615555648722522388044622506E-30_dp, &
 -.2516994994284095106080616830293E-30_dp, +.1127259818927510206370368804181E-30_dp, &
 -.5069814875800460855562584719360E-31_dp /)
!***FIRST EXECUTABLE STATEMENT  DE1
IF (first) THEN
  eta = 0.1 * EPSILON(0.0_dp)
  ntae10 = initds(ae10cs, 50, eta)
  ntae11 = initds(ae11cs, 60, eta)
  ntae12 = initds(ae12cs, 41, eta)
  nte11 = initds(e11cs, 29, eta)
  nte12 = initds(e12cs, 25, eta)
  ntae13 = initds(ae13cs, 50, eta)
  ntae14 = initds(ae14cs, 64, eta)
  xmaxt = -LOG( TINY(0.0_dp) )
  xmax = xmaxt - LOG(xmaxt)
END IF
first = .false.
IF (x <= -1._dp) THEN
  IF (x <= -32._dp) THEN
    fn_val = EXP(-x) / x * (1._dp + dcsevl(64._dp/x+1._dp, ae10cs, ntae10))
    RETURN
  END IF
  IF (x <= -8._dp) THEN
    fn_val = EXP(-x) / x * (1._dp + dcsevl((64._dp/x+5._dp)/3._dp, ae11cs, ntae11))
    RETURN
  END IF
  IF (x <= -4._dp) THEN
    fn_val = EXP(-x) / x * (1._dp + dcsevl(16._dp/x+3._dp, ae12cs, ntae12))
    RETURN
  END IF
  fn_val = -LOG(-x) + dcsevl((2._dp*x+5._dp)/3._dp, e11cs, nte11)
  RETURN
END IF
IF (x <= 1.0_dp) THEN
  IF (x == 0._dp) CALL xermsg('SLATEC', 'DE1', 'X IS 0')
  fn_val = (-LOG(ABS(x)) - 0.6875_dp+x) + dcsevl(x, e12cs, nte12)
  RETURN
END IF
IF (x <= 4.0_dp) THEN
  fn_val = EXP(-x) / x * (1._dp + dcsevl((8._dp/x-5._dp)/3._dp, ae13cs, ntae13))
  RETURN
END IF
IF (x <= xmax) THEN
  fn_val = EXP(-x) / x * (1._dp + dcsevl(8._dp/x-1._dp, ae14cs, ntae14))
  RETURN
END IF
CALL xermsg('SLATEC', 'DE1', 'X SO BIG E1 UNDERFLOWS')
fn_val = 0._dp
RETURN
END FUNCTION de1



!DECK DEI
FUNCTION dei(x) RESULT(fn_val)
!***BEGIN PROLOGUE  DEI
!***PURPOSE  Compute the exponential integral Ei(X).
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C5
!***TYPE      REAL (dp) (EI-S, DEI-D)
!***KEYWORDS  EI FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
! DEI calculates the REAL (dp) exponential integral, Ei(X), for
! positive REAL (dp) argument X and the Cauchy principal value
! for negative X.  If principal values are used everywhere, then, for
! all X,
!    Ei(X) = -E1(-X)
! or
!    E1(X) = -Ei(-X).
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DE1
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   891115  Modified prologue description.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DEI
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val
!***FIRST EXECUTABLE STATEMENT  DEI
fn_val = -de1(-x)
RETURN
END FUNCTION dei


  function lambertw(x)
    implicit none
    real(kp) :: lambertw
    real(kp), intent(in) :: x
    real(kp) :: w,wTimesExpW,wPlusOneTimesExpW
    real(kp), parameter :: tol=1d-15

!rough estimation

    if ((x.le.500._kp).and.(x.ge.0._kp)) then
       w = 0.665_kp * (1._kp + 0.0195*log(x+1._kp)) &
            *log(x+1._kp) + 0.04_kp
    elseif (x.gt.500._kp) then
       w = log(x-4._kp) - (1._kp - 1._kp/log(x)) &
            *log(log(x))
    elseif (x.ge.-1._kp/exp(1._kp)) then
       w=-0.5
    else
       stop 'x<-1/e'
    endif

!recurrence
       
    do
       wTimesExpW = w*exp(w)
       wPlusOneTimesExpW = (w+1._kp)*exp(w)
       if (tol.gt.(abs((x-wTimesExpW)/wPlusOneTimesExpW))) then
          exit
       else
          w = w-(wTimesExpW-x) &
               /(wPlusOneTimesExpW-(w+2._kp)*(wTimesExpW-x) &
               /(2._kp*w+2._kp))
       endif
    enddo
    
    lambertw = w
    
  end function lambertw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                        !
!           CIRCULAR AND HYPERBOLIC FUNCTIONS            !
!                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this exists in recent fortran

!  function acosh(x)
!  USE nrutil, ONLY : nrerror
!    implicit none
!    real(kp) ::acosh
!    real(kp), intent(in) ::x
!	if (x<1.) call nrerror('specialinf/acosh: bad argument x<1. gives non real value for acosh(x)')
!	acosh=log(x+sqrt(x**2-1._kp))
!  end function acosh

!  function asinh(x)
!    implicit none
!    real(kp) ::asinh
!    real(kp), intent(in) ::x
!	asinh=log(x+sqrt(x**2+1._kp))
!  end function asinh

!  function atanh(x)
!  USE nrutil, ONLY : nrerror
!    implicit none
!    complex(kp) ::atanh
!    real(kp), intent(in) ::x
!	atanh=0.5_kp*log((1._kp+x)/(1._kp-x))
!  end function atanh

!ArcTan function expressed in terms of logarithmic functions, works
!for complex z. Normally supported in fortran08, already implemented
!in gfortran but not in ifort
  function atan_ito_log(x)
    implicit none
    complex(kp), intent(in) :: x
    complex(kp) :: atan_ito_log
    atan_ito_log=0.5_kp*(0._kp,1._kp)*log((1._kp-x*(0._kp,1._kp))/(1._kp+x*(0._kp,1._kp)))
  end function atan_ito_log


  function atanh_ito_log(x)
    implicit none
    complex(kp), intent(in) :: x
    complex(kp) :: atanh_ito_log
    atanh_ito_log=0.5_kp*log((1._kp+x)/(1._kp-x))
  end function atanh_ito_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                        !
!                      POLYLOGARITHM                     !
!                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function polylog(z,s) !polylogarithm functions
    implicit none
    complex(kp), intent(in) :: z,s
    complex(kp) :: polylog
    integer ::i,Nsum
    polylog=0._kp
    Nsum=1000

    !Uses inversion formulas for dilogarithms
    if (s.eq.2._kp.and.real(z,kp).gt.1._kp) then
       do i=1,Nsum
          polylog=polylog+(1._kp/z)**(real(i,kp))/(real(i,kp)**2)
       enddo
       polylog=acos(-1._kp)**2/3._kp-0.5_kp*log(z)**2-cmplx(0._kp,1._kp,kp)*acos(-1._kp)*log(z)-polylog
    elseif (s.eq.2._kp.and.real(z,kp).lt.-1._kp) then
       do i=1,Nsum
          polylog=polylog+(1._kp/z)**(real(i,kp))/(real(i,kp)**2)
       enddo
       polylog=-acos(-1._kp)**2/6._kp-0.5_kp*(log(-z))**2-polylog ! An error has been corrected from numerical recipe source (3->6)
    else
    do i=1,Nsum
      polylog=polylog+z**(real(i,kp))/(real(i,kp)**s)
    enddo
    endif
  end function polylog


end module specialinf
