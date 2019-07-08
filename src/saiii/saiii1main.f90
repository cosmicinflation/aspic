program saiii1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use saiii1sr, only : saiii1_epsilon_one, saiii1_epsilon_two, saiii1_epsilon_three
  use saiii1reheat, only : saiii1_lnrhoreh_max, saiii1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use saiii1sr, only : saiii1_norm_potential, saiii1_x_endinf
  use saiii1sr, only : saiii1_norm_deriv_potential, saiii1_norm_deriv_second_potential
  use saiii1sr, only : saiii1_numacc_mumin, saiii1_numacc_efoldmax
  use saiii1reheat, only : saiii1_x_rreh, saiii1_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use saiiicommon, only : saiii_x_potmax, beta0,beta1,beta2,beta3
  use saiiicommon, only : saiii_alpha_two, saiii_alpha_one
  
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

  real(kp) :: V1,x,V


  
  Pstar = powerAmpScalar


  call delete_file('saiii1_potential.dat')
  call delete_file('saiii1_slowroll.dat')

  
  alpha = 0.3_kp
  beta = 0.9_kp
  mu=1._kp

  n=250

  
  do i=1,n
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)        

     V1 = saiii1_norm_potential(x,alpha,beta,mu)
        

     call livewrite('saiii1_potential.dat',x,V1)

     eps1 = saiii1_epsilon_one(x,alpha,beta,mu)
     eps2 = saiii1_epsilon_two(x,alpha,beta,mu)
     eps3 = saiii1_epsilon_three(x,alpha,beta,mu)

     call livewrite('saiii1_slowroll.dat',x,eps1,eps2,eps3)

  enddo

 
  call delete_file('saiii1_predic.dat')
  call delete_file('saiii1_nsr.dat')


!$omp parallel sections &
!$omp default(shared) &
!$omp private(betamin,betamax,alphamin,alphamax) &  
!$omp private(alpha,beta,l,k,i,j,mu) &  
!$omp private(lnRhoRehMax,xend,lnRhoRehMin,lnRhoReh,bfoldstar) &
!$omp private(xstar,w,eps1,eps2,eps3,logErehGeV,Treh,ns,r)


!$omp section
  
  call aspicwrite_header('saiii1p',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','beta '/))
  
  w=0._kp


  betamin = 0.05_kp
  betamax = 1.3_kp
  
  do l=1,Nb
     beta = betamin +  real(l-1,kp)*(betamax - betamin)/real(Nb-1,kp)

     if (beta.gt.beta0) then
        alphamax = 0.95*saiii_alpha_two(beta)
     else
        alphamax = 1.2_kp
     endif

     alphamin = alphamax/5._kp

     print *,'beta alphamin alphamax',beta, alphamin,alphamax
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        mumin = saiii1_numacc_mumin(120._kp,alpha,beta)

        print *,'mumin= ',mumin
        
        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = saiii1_numacc_efoldmax(alpha,beta,mu)
           print *,'alpha,beta,mu,efoldNumAccMax= ',alpha,beta,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = saiii1_x_endinf(alpha,beta,mu)

           lnRhoRehMax = saiii1_lnrhoreh_max(alpha,beta,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = saiii1_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = saiii1_epsilon_one(xstar,alpha,beta,mu)
              eps2 = saiii1_epsilon_two(xstar,alpha,beta,mu)
              eps3 = saiii1_epsilon_three(xstar,alpha,beta,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('saiii1_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('saiii1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,alpha,beta/))

           end do

        end do
     enddo
  enddo
  call aspicwrite_end()

!$omp section
  
  call aspicwrite_header('saiii1m',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','beta '/))
  
  w=0._kp


  betamin = -1.3_kp
  betamax = -0.05_kp
  
  do l=1,Nb
     beta = betamin +  real(l-1,kp)*(betamax - betamin)/real(Nb-1,kp)

     if (beta.le.-1._kp) then
        alphamin = 0.95*saiii_alpha_one(beta)
        alphamax = alphamin/5._kp
     else
        alphamin = -1.2_kp
        alphamax = -0.05_kp
     endif

     print *,'beta alphamin alphamax',beta, alphamin,alphamax
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        mumin = saiii1_numacc_mumin(120._kp,alpha,beta)

        print *,'mumin= ',mumin
        
        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = saiii1_numacc_efoldmax(alpha,beta,mu)
           print *,'alpha,beta,mu,efoldNumAccMax= ',alpha,beta,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = saiii1_x_endinf(alpha,beta,mu)

           lnRhoRehMax = saiii1_lnrhoreh_max(alpha,beta,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = saiii1_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = saiii1_epsilon_one(xstar,alpha,beta,mu)
              eps2 = saiii1_epsilon_two(xstar,alpha,beta,mu)
              eps3 = saiii1_epsilon_three(xstar,alpha,beta,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('saiii1_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('saiii1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

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

  beta = 0.9_kp
  alpha = 0.7
  mu = 10._kp*saiii1_numacc_mumin(120._kp,alpha,beta)

  xend = saiii1_x_endinf(alpha,beta,mu)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = saiii1_x_rrad(alpha,beta,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = saiii1_epsilon_one(xstar,alpha,beta,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  saiii1_epsilon_one(xend,alpha,beta,mu)
     VendOverVstar = saiii1_norm_potential(xend,alpha,beta,mu) &
          /saiii1_norm_potential(xstar,alpha,beta,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = saiii1_x_rreh(alpha,beta,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = saiii1_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program saiii1main
