program saiii2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use saiii2sr, only : saiii2_epsilon_one, saiii2_epsilon_two, saiii2_epsilon_three
  use saiii2reheat, only : saiii2_lnrhoreh_max, saiii2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use saiii2sr, only : saiii2_norm_potential, saiii2_x_endinf
  use saiii2sr, only : saiii2_norm_deriv_potential, saiii2_norm_deriv_second_potential
  use saiii2sr, only : saiii2_numacc_mumin, saiii2_numacc_efoldmax
  use saiii2reheat, only : saiii2_x_rreh, saiii2_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use saiiicommon, only : saiii_x_potmax, beta1,beta2,beta3,beta0
  use saiiicommon, only : saiii_alpha_potneg, saiii_alpha_one, saiii_alpha_two, saiii_alpha_three
  use saiiicommon, only : saiii_x_epsoneunity, saiii_x_potzero
  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh



  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n,k,l
  integer :: npts = 15


  integer, parameter :: Na = 2
  real(kp) :: alphamin
  real(kp) :: alphamax

  integer, parameter :: Nb = 3
  real(kp) :: betamin
  real(kp) :: betamax

  integer, parameter :: Nmu=60
  real(kp) :: mumin = 0.1_kp
  real(kp) :: mumax = 100._kp

  real(kp) :: xmin = 0
  real(kp) :: xmax = 7

  real(kp) :: alpha,beta,mu,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax, efoldmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: V1,x,V,alpha1,alpha2,alpha3,alphan, alphav
  real(kp), dimension(3) :: xvzero
  real(kp), dimension(2) :: xeps
  real(kp) :: xvmax
  
  logical, parameter :: dumpExtra = .true.

  
  Pstar = powerAmpScalar


  if (dumpExtra) then
  
     call delete_file('saiii2_alpha1.dat')
     call delete_file('saiii2_alpha2.dat')
     call delete_file('saiii2_alpha3.dat')
     call delete_file('saiii2_alphan.dat')
     call delete_file('saiii2_alphav.dat')


     betamin = -2
     betamax = 4

     n=500

     do i=1,n
        beta = betamin + real(i-1,kp)*(betamax-betamin)/real(n-1,kp)        

        alpha2=0.
        alphan=0.
        alphav=0.

        if (beta.gt.beta0) then
           alpha2 = saiii_alpha_two(beta)

           call livewrite('saiii2_alpha2.dat',beta,alpha2)

        endif

        if ((beta.gt.beta3).and.(beta.lt.beta2)) then
           alphan = saiii_alpha_potneg(beta)

           call livewrite('saiii2_alphan.dat',beta,alphan)

        endif

        if ((beta.gt.-2._kp).and.(beta.lt.0._kp)) then


           alphav = -1._kp/(beta+2._kp)

           call livewrite('saiii2_alphav.dat',beta,alphav)

        endif

     enddo

     betamin = -4._kp
     betamax = -1._kp
     do i=1,n

        beta = betamin + real(i-1,kp)*(betamax-betamin)/real(n-1,kp) 

        alpha1 = saiii_alpha_one(beta)

        call livewrite('saiii2_alpha1.dat',beta,alpha1)

     enddo

     betamin = beta1
     betamax = -1._kp
     do i=1,n

        beta = betamin + real(i-1,kp)*(betamax-betamin)/real(n-1,kp) 

        alpha3 = saiii_alpha_three(beta)

        call livewrite('saiii2_alpha3.dat',beta,alpha3)

     enddo


     n=250
     
     call delete_file('saiii2_xV_1.dat')
     call delete_file('saiii2_xeps_1.dat')

     alphamin = -1.2_kp
     alphamax = -0.01_kp
     beta = -0.6

     do i=1,n

        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(n-1,kp) 

        xvmax = saiii_x_potmax(alpha,beta,mu)
        xvzero = saiii_x_potzero(alpha,beta,mu)


        call livewrite('saiii2_xV_1.dat',alpha,xvmax,xvzero(1),xvzero(2))

        xeps = saiii_x_epsoneunity(alpha,beta,mu=10._kp)

        call livewrite('saiii2_xeps_1.dat',alpha,xeps(1),xeps(2))

     enddo


     call delete_file('saiii2_xV_2.dat')
     call delete_file('saiii2_xeps_2.dat')

     alphamin = -1.2_kp
     alphamax = -0.01_kp
     beta = -0.1

     do i=1,n

        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(n-1,kp) 

        xvmax = saiii_x_potmax(alpha,beta,mu)
        xvzero = saiii_x_potzero(alpha,beta,mu)

        call livewrite('saiii2_xV_2.dat',alpha,xvmax,xvzero(1),xvzero(2))
           
        xeps = saiii_x_epsoneunity(alpha,beta,mu=10._kp)

        call livewrite('saiii2_xeps_2.dat',alpha,xeps(1),xeps(2))

     enddo


     call delete_file('saiii2_xV_3.dat')
     call delete_file('saiii2_xeps_3.dat')

     alphamin = 0.01_kp
     alphamax = 1.2_kp
     beta = 0.1

     do i=1,n

        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(n-1,kp) 

        xvmax = saiii_x_potmax(alpha,beta,mu)
        xvzero = saiii_x_potzero(alpha,beta,mu)


        call livewrite('saiii2_xV_3.dat',alpha,xvmax,xvzero(1),xvzero(2))

        xeps = saiii_x_epsoneunity(alpha,beta,mu=10._kp)

        call livewrite('saiii2_xeps_3.dat',alpha,xeps(1),xeps(2))

     enddo



     call delete_file('saiii2_xV_4.dat')
     call delete_file('saiii2_xeps_4.dat')

     alphamin = 0.01_kp
     alphamax = 1.2_kp
     beta = 0.6

     do i=1,n

        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(n-1,kp) 

        xvmax = saiii_x_potmax(alpha,beta,mu)
        xvzero = saiii_x_potzero(alpha,beta,mu)


        call livewrite('saiii2_xV_4.dat',alpha,xvmax,xvzero(1),xvzero(2),xvzero(3))

        xeps = saiii_x_epsoneunity(alpha,beta,mu=10._kp)

        call livewrite('saiii2_xeps_4.dat',alpha,xeps(1),xeps(2))

     enddo



  endif




  
  call delete_file('saiii2_predic.dat')
  call delete_file('saiii2_nsr.dat')


!$omp parallel sections &
!$omp default(shared) &
!$omp private(betamin,betamax,alphamin,alphamax) &  
!$omp private(alpha,beta,l,k,i,j,mu) &  
!$omp private(lnRhoRehMax,xend,lnRhoRehMin,lnRhoReh,bfoldstar) &
!$omp private(xstar,w,eps1,eps2,eps3,logErehGeV,Treh,ns,r)


!$omp section


  call aspicwrite_header('saiii2p',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','beta '/))

  w=0._kp


  betamin = 0.02_kp
  betamax = 0.9*beta2

  do l=1,Nb
     beta = betamin +  real(l-1,kp)*(betamax - betamin)/real(Nb-1,kp)

     alphamin = 1.01_kp*saiii_alpha_potneg(beta)
     alphamax = 5._kp*alphamin

     print *,'beta alphamin alphamax',beta, alphamin,alphamax

     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        mumin = saiii2_numacc_mumin(120._kp,alpha,beta)

        print *,'mumin= ',mumin

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = saiii2_numacc_efoldmax(alpha,beta,mu)
           print *,'alpha,beta,mu,efoldNumAccMax= ',alpha,beta,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = saiii2_x_endinf(alpha,beta,mu)

           lnRhoRehMax = saiii2_lnrhoreh_max(alpha,beta,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = saiii2_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = saiii2_epsilon_one(xstar,alpha,beta,mu)
              eps2 = saiii2_epsilon_two(xstar,alpha,beta,mu)
              eps3 = saiii2_epsilon_three(xstar,alpha,beta,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('saiii2_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('saiii2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,alpha,beta/))

           end do

        end do
     enddo
  enddo
  call aspicwrite_end()

!$omp section

  call aspicwrite_header('saiii2m',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','beta '/))

  w=0._kp


  betamin = 0.9*beta3
  betamax = -0.02_kp

  do l=1,Nb
     beta = betamin +  real(l-1,kp)*(betamax - betamin)/real(Nb-1,kp)

     alphamax = 1.01_kp*saiii_alpha_potneg(beta)
     alphamin = 5._kp*alphamax

     print *,'beta alphamin alphamax',beta, alphamin,alphamax

     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        mumin = saiii2_numacc_mumin(120._kp,alpha,beta)

        print *,'mumin= ',mumin

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = saiii2_numacc_efoldmax(alpha,beta,mu)
           print *,'alpha,beta,mu,efoldNumAccMax= ',alpha,beta,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = saiii2_x_endinf(alpha,beta,mu)

           lnRhoRehMax = saiii2_lnrhoreh_max(alpha,beta,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = saiii2_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = saiii2_epsilon_one(xstar,alpha,beta,mu)
              eps2 = saiii2_epsilon_two(xstar,alpha,beta,mu)
              eps3 = saiii2_epsilon_three(xstar,alpha,beta,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('saiii2_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('saiii2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,alpha,beta/))

           end do

        end do
     enddo
  enddo

  call aspicwrite_end()


!$omp end parallel sections



  write(*,*)
  write(*,*)'Testing Rrad/Rreh'



  lnRradmin=-42
  lnRradmax = 10

  beta = 0.01_kp
  alpha = 1._kp
  mu = 10._kp*saiii2_numacc_mumin(120._kp,alpha,beta)

  xend = saiii2_x_endinf(alpha,beta,mu)

  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = saiii2_x_rrad(alpha,beta,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = saiii2_epsilon_one(xstar,alpha,beta,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  saiii2_epsilon_one(xend,alpha,beta,mu)
     VendOverVstar = saiii2_norm_potential(xend,alpha,beta,mu) &
          /saiii2_norm_potential(xstar,alpha,beta,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = saiii2_x_rreh(alpha,beta,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = saiii2_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program saiii2main
