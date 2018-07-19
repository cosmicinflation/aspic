!test the reheating derivation from slow-roll
program rcqimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rcqisr, only : rcqi_epsilon_one, rcqi_epsilon_two, rcqi_epsilon_three
  use rcqireheat, only : rcqi_lnrhoreh_max, rcqi_x_star
  use infinout, only : delete_file, livewrite, has_not_shifted
  use srreheat, only : log_energy_reheat_ingev

  use rcqisr, only : rcqi_norm_potential, rcqi_x_endinf
  use rcqireheat, only : rcqi_x_rreh, rcqi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(1:6) ::alphavalues

  real(kp) :: alphamin,alphamax,eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nalpha

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  alphavalues(1)=(10._kp)**(-2.)
  alphavalues(2)=(10._kp)**(-0.7)
  alphavalues(3)=(10._kp)**(-0.6)
  alphavalues(4)=(10._kp)**(-0.55)
  alphavalues(5)=(10._kp)**(-0.5)
  alphavalues(6)=(10._kp)**(-0.485)

  alphamin=10.**(-2.0)
  alphamax=10.**(-0.4)



  Pstar = powerAmpScalar

  call delete_file('rcqi_predic.dat')
  call delete_file('rcqi_nsr.dat')

  call aspicwrite_header('rcqirad',labeps12,labnsr,labbfoldreh,(/'alpha'/))
  
  do j=1,1000
     w = 1._kp/3._kp
     alpha=alphamin+(alphamax-alphamin)*(real(j-1,kp)/real(1000,kp))



     lnRhoRehMin = lnRhoNuc
     xEnd = rcqi_x_endinf(alpha)
     lnRhoRehMax = rcqi_lnrhoreh_max(alpha,xend,Pstar)

     print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = rcqi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = rcqi_epsilon_one(xstar,alpha)
        eps2 = rcqi_epsilon_two(xstar,alpha)
        eps3 = rcqi_epsilon_three(xstar,alpha)

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha/))
        
        if (has_not_shifted(0.005_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
           cycle
        endif

        if ((eps1.lt.1e-5).or.(eps1.gt.0.1) &
             .and.(eps2.lt.0.1).and.(eps2.gt.0.15)) cycle


        call livewrite('rcqi_predic.dat',alpha,w,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('rcqi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        
        
     end do

  end do

  call aspicwrite_end()

  call aspicwrite_header('rcqi',labeps12,labnsr,labbfoldreh,(/'alpha'/))

  do j=1,1000
     w=0._kp

     alpha=alphamin+(alphamax-alphamin)*(real(j-1,kp)/real(1000,kp))

     lnRhoRehMin = lnRhoNuc
     xEnd = rcqi_x_endinf(alpha)
     lnRhoRehMax = rcqi_lnrhoreh_max(alpha,xend,Pstar)

     print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = rcqi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = rcqi_epsilon_one(xstar,alpha)
        eps2 = rcqi_epsilon_two(xstar,alpha)
        eps3 = rcqi_epsilon_three(xstar,alpha)

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha/))
        
        
         if (has_not_shifted(0.005_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
            cycle
         endif

           
         if ((eps1.lt.1e-5).or.(eps1.gt.0.1) &
              .and.(eps2.lt.0.1).and.(eps2.gt.0.15)) cycle


        call livewrite('rcqi_predic.dat',alpha,w,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('rcqi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

  call aspicwrite_end()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('rcqi_predic_summarized.dat') 
  nalpha=1000
  alphamin=10._kp**(-2.)
  alphamax=10._kp**(-0.3)
  w=0._kp
  do j=1,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
     xEnd = rcqi_x_endinf(alpha)
     lnRhoReh = lnRhoNuc
     xstarA = rcqi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
     eps1A = rcqi_epsilon_one(xstarA,alpha)
     eps2A = rcqi_epsilon_two(xstarA,alpha)
     eps3A = rcqi_epsilon_three(xstarA,alpha)
     nsA = 1._kp - 2._kp*eps1A - eps2A
     rA = 16._kp*eps1A
     lnRhoReh = rcqi_lnrhoreh_max(alpha,xend,Pstar)
     xstarB = rcqi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
     eps1B = rcqi_epsilon_one(xstarB,alpha)
     eps2B = rcqi_epsilon_two(xstarB,alpha)
     eps3B = rcqi_epsilon_three(xstarB,alpha)
     nsB = 1._kp - 2._kp*eps1B - eps2B
     rB =16._kp*eps1B
     call livewrite('rcqi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
  enddo

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 5e-2

  xEnd = rcqi_x_endinf(alpha)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rcqi_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rcqi_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rcqi_epsilon_one(xend,alpha)
     VendOverVstar = rcqi_norm_potential(xend,alpha)/rcqi_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rcqi_x_rreh(alpha,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rcqi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program rcqimain
