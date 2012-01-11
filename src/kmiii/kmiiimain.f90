!test the reheating derivation from slow-roll
program kmiiimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use kmiiisr, only : kmiii_epsilon_one, kmiii_epsilon_two, kmiii_epsilon_three,kmiii_alphamin
  use kmiiireheat, only : kmiii_lnrhoend, kmiii_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  real(kp) :: alpha,beta,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::betamin,betamax,alphamin,alphamax

  Pstar = powerAmpScalar


!Calculates the prior space data

  betamin=0.1_kp
  betamax=10._kp
  call delete_file('kmiii_prior.dat')
  do i=0,1000
    beta=betamin+(betamax-betamin)*real(i,kp)/real(1000,kp)
    alpha=kmiii_alphamin(beta)
    call livewrite('kmiii_prior.dat',beta,alpha,beta*exp(1._kp)) !given beta, writes alphamin and alphamax between which alpha is allowed to vary.
  end do


!Calculates the Reaheating Constrains

  call delete_file('kmiii_predic.dat')
  call delete_file('kmiii_nsr.dat')

!  w = 1._kp/3._kp
  w=0._kp

   betamin=0.0005
   betamax=5._kp

   do i=0,20
     beta=betamin*(betamax/betamin)**(real(i,kp)/real(20,kp))
     alphamin=kmiii_alphamin(beta)*1.000001
     alphamax=beta*exp(1._kp)*0.999999
     do j=0,10
     alpha=alphamin+(alphamax-alphamin)*real(j,kp)/real(10,kp)



    lnRhoRehMin = lnRhoNuc
    lnRhoRehMax = kmiii_lnrhoend(alpha,beta,Pstar)
  

    print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax
  

  do k=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(k-1,kp)/real(npts-1,kp)

       xstar = kmiii_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)

       eps1 = kmiii_epsilon_one(xstar,alpha,beta)
       eps2 = kmiii_epsilon_two(xstar,alpha,beta)
       eps3 = kmiii_epsilon_three(xstar,alpha,beta)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

!       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'r=',r,'ns=',ns

       call livewrite('kmiii_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('kmiii_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  end do  

 end do



end program kmiiimain
