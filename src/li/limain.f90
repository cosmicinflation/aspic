!test the reheating derivation from slow-roll
program limain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lisr, only : li_epsilon_one, li_epsilon_two, li_epsilon_three
  use lireheat, only : li_lnrhoreh_max, li_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use lisr, only : li_norm_potential, li_x_endinf, li_alphamin
  use lisr, only : li_x_epsoneunity, li_efold_primitive
  use lireheat, only : li_x_rreh, li_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : labeps12, labnsr, labbfoldreh
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::alphamin,alphamax

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  real(kp) :: DeltaNmax
  integer :: nalpha

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend
  real(kp), dimension(2) :: xepsones
  Pstar = powerAmpScalar

  call delete_file('li_predic.dat')
  call delete_file('li_nsr.dat')

  call aspicwrite_header('li',labeps12,labnsr,labbfoldreh,(/'alpha'/))
  
!!!!!!!!!!!!!!!!!!
!!!  Priors    !!!
!!!!!!!!!!!!!!!!!!

  nalpha=100
  alphamin=-1._kp
  alphamax=-0.2_kp

  call delete_file('li_alphamin.dat')
  do k=0,nalpha
     alpha = alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp))
     xepsones = li_x_epsoneunity(alpha)
     DeltaNmax = li_efold_primitive(xepsones(1),alpha)- &
                 li_efold_primitive(xepsones(2),alpha) 
     call livewrite('li_alphamin.dat',alpha,DeltaNmax)   
  end do



!!!!!!!!!!!!!!!!!!
!!!  alpha>0   !!!
!!!!!!!!!!!!!!!!!!

  npts = 20

  alphamin=0.003
  alphamax=2._kp
  nalpha=20


  !  w = 1._kp/3._kp
  w=0._kp

  do k=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp))

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = li_lnrhoreh_max(alpha,Pstar)

     print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = li_epsilon_one(xstar,alpha)
        eps2 = li_epsilon_two(xstar,alpha)
        eps3 = li_epsilon_three(xstar,alpha)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1*=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('li_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('li_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha/))

     end do

  end do


!!!!!!!!!!!!!!!!!!
!!!  alpha<0   !!!
!!!!!!!!!!!!!!!!!!

  npts = 5
  
  nalpha=100

  alphamin=-0.35
  alphamax=-0.1


  !  w = 1._kp/3._kp
  w=0._kp

  do k=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !logarithmic step
     alpha=-exp(log(-alphamin)/((log(-alphamin)/log(-alphamax))** &
          (real(k,kp)/real(nalpha,kp)))) !adapted step
     alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp))! arithmetic step

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = li_lnrhoreh_max(alpha,Pstar)

     print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


     do i=0,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i,kp)/real(npts-1,kp)

        xstar = li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = li_epsilon_one(xstar,alpha)
        eps2 = li_epsilon_two(xstar,alpha)
        eps3 = li_epsilon_three(xstar,alpha)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1*=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('li_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('li_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha/))
        
     end do

  end do

  call aspicwrite_end()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('li_predic_summarized.dat') 
  nalpha=1000
  alphamin=0.005
  alphamax=0.1
  w=0._kp
  do j=1,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
     lnRhoReh = lnRhoNuc
     xstarA = li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1A = li_epsilon_one(xstarA,alpha)
     eps2A = li_epsilon_two(xstarA,alpha)
     eps3A = li_epsilon_three(xstarA,alpha)
     nsA = 1._kp - 2._kp*eps1A - eps2A
     rA = 16._kp*eps1A
     lnRhoReh = li_lnrhoreh_max(alpha,Pstar)
     xstarB = li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1B = li_epsilon_one(xstarB,alpha)
     eps2B = li_epsilon_two(xstarB,alpha)
     eps3B = li_epsilon_three(xstarB,alpha)
     nsB = 1._kp - 2._kp*eps1B - eps2B
     rB =16._kp*eps1B
     call livewrite('li_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
  enddo



  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = -0.3
 

  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = li_x_rrad(alpha,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar


     eps1 = li_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = li_x_endinf(alpha)
     eps1end =  li_epsilon_one(xend,alpha)
     VendOverVstar = li_norm_potential(xend,alpha)/li_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = li_x_rreh(alpha,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar=', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar
     
  enddo





end program limain
