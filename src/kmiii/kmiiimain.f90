!test the reheating derivation from slow-roll
program kmiiimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use kmiiisr, only : kmiii_epsilon_one, kmiii_epsilon_two, kmiii_epsilon_three
  use kmiiisr, only : kmiii_x_endinf_appr
  use kmiiireheat, only : kmiii_lnrhoreh_max, kmiii_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use kmiiisr, only : kmiii_norm_potential, kmiii_x_endinf
  use kmiiireheat, only : kmiii_x_rreh, kmiii_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 15, nalpha=20, nnu=30
  integer :: nbetaprior, nnuxend

  real(kp) :: alpha,beta,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::betamin,betamax,alphamin,alphamax,xendAppr,xend,nu,numin,numax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End

  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                   !!!
!!!    Checking the  approached       !!!
!!!   analytical formula for xend     !!!
!!!                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nnuxend=100
  numin=10._kp**(5.)
  numax=10._kp**(7.)
  call delete_file('kmiii_xend.dat')
  do i=0,nnuxend
     nu=numin*(numax/numin)**(real(i,kp)/real(nnuxend,kp))! logarithmic step
     alpha=nu**(5._kp/3._kp)
     beta=nu**(2._kp/3._kp)
     xendAppr=kmiii_x_endinf_appr(alpha,beta)
     xend=kmiii_x_endinf(alpha,beta)
     call livewrite('kmiii_xend.dat',nu,xendAppr,xend)
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

  call aspicwrite_header('kmiii',labeps12,labnsr,labbfoldreh,(/'calVs'/))
  
  !  w = 1._kp/3._kp
  w=0._kp

  numin=10._kp**(5.) 
  numax=10._kp**(7.)

  do i=0,nnu
     nu=numin*(numax/numin)**(real(i,kp)/real(nnu,kp))! logarithmic step
     alpha=nu**(5._kp/3._kp)
     beta=nu**(2._kp/3._kp)


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = kmiii_lnrhoreh_max(alpha,beta,Pstar)


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

           call livewrite('kmiii_predic.dat',nu,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('kmiii_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/nu/))
           
       end do

  end do

  call aspicwrite_end()
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                   !!!
!!!         Further Tests             !!!
!!!                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  nu = 10._kp**(6._kp)
  alpha=nu**(5._kp/3._kp)
  beta=nu**(2._kp/3._kp)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = kmiii_x_rrad(alpha,beta,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = kmiii_epsilon_one(xstar,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = kmiii_x_endinf(alpha,beta)
     eps1end =  kmiii_epsilon_one(xend,alpha,beta)
     VendOverVstar = kmiii_norm_potential(xend,alpha,beta)/kmiii_norm_potential(xstar,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = kmiii_x_rreh(alpha,beta,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = kmiii_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program kmiiimain
