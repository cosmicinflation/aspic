!test the reheating derivation from slow-roll
program kmiiimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use kmiiisr, only : kmiii_epsilon_one, kmiii_epsilon_two, kmiii_epsilon_three
  use kmiiisr, only : kmiii_alphamin, kmiii_x_endinf, kmiii_x_endinf_appr
  use kmiiireheat, only : kmiii_lnrhoend, kmiii_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20, nalpha=20, nbeta=15
  integer :: nbetaprior, nbetaxend

  real(kp) :: alpha,beta,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::betamin,betamax,alphamin,alphamax,xendAppr,xend

  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                   !!!
!!!  Calculates the prior space data  !!!
!!!                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  betamin=10._kp**(10.)
  betamax=10._kp**(14.)
  nbetaprior=100
  call delete_file('kmiii_prior.dat')
  do i=0,nbetaprior
    beta=betamin*(betamax/betamin)**(real(i,kp)/real(nbetaprior,kp))! logarithmic step
    alpha=kmiii_alphamin(beta)
    call livewrite('kmiii_prior.dat',beta,alpha) !given beta, writes the minimum value of alpha for which inflation can sto by slow roll violation in the large field domain of the potential
  end do
  print*,'priors written'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                   !!!
!!!    Checking the  approached       !!!
!!!   analytical formula for xend     !!!
!!!                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   nbetaxend=100
   betamin=10._kp**(10.)
   betamax=10._kp**(14.)
   call delete_file('kmiii_xend.dat')
   do i=0,nbetaxend
     beta=betamin*(betamax/betamin)**(real(i,kp)/real(nbetaxend,kp))! logarithmic step
     alpha=kmiii_alphamin(beta)*10._kp
     !alpha=beta/1._kp
     xendAppr=kmiii_x_endinf_appr(alpha,beta)
     xend=kmiii_x_endinf(alpha,beta)
     call livewrite('kmiii_xend.dat',beta,xendAppr,xend)
   end do
  print*,'xend written'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                   !!!
!!!    Calculates the reheating       !!!
!!!  consistent slow roll predictions !!!
!!!                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('kmiii_predic.dat')
  call delete_file('kmiii_nsr.dat')

!  w = 1._kp/3._kp
  w=0._kp

   betamin=10._kp**(9.) 
   betamax=10._kp**(15.)

   do i=0,nbeta
     beta=betamin*(betamax/betamin)**(real(i,kp)/real(nbeta,kp))! logarithmic step
     alphamin=kmiii_alphamin(beta)*1.1_kp
     alphamax=beta/100._kp
     do j=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))! logarithmic step

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

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'r=',r,'ns=',ns

       call livewrite('kmiii_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('kmiii_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  end do  

 end do



end program kmiiimain
