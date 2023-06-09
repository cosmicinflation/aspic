program saiii3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use saiii3sr, only : saiii3_epsilon_one, saiii3_epsilon_two, saiii3_epsilon_three
  use saiii3reheat, only : saiii3_lnrhoreh_max, saiii3_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use saiii3sr, only : saiii3_norm_potential, saiii3_x_endinf
  use saiii3sr, only : saiii3_norm_deriv_potential, saiii3_norm_deriv_second_potential
  use saiii3reheat, only : saiii3_x_rreh, saiii3_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use saiiicommon, only : beta0, beta1, beta2, beta3
  use saiiicommon, only : saiii_alpha_one, saiii_alpha_two, saiii_alpha_three
  use saiiicommon, only : saiii_x_epsoneunity, saiii_x_potzero

  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n,k,l
  integer :: npts = 20


  integer, parameter :: Na = 2
  real(kp) :: alphamin
  real(kp) :: alphamax

  integer, parameter :: Nb = 3
  real(kp) :: betamin
  real(kp) :: betamax
  
  integer, parameter :: Nmu=59
  real(kp) :: mumin = 0.1_kp
  real(kp) :: mumax = 100._kp
  
  real(kp) :: xmin = 0
  real(kp) :: xmax = 10
  
  real(kp) :: alpha,beta,mu,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: V1,x,V
  real(kp), dimension(3) :: xvzero
  real(kp), dimension(2) :: xeps

  
  logical :: dumpExtra = .true.
  
  Pstar = powerAmpScalar


 

  if (dumpExtra) then


     call delete_file('saiii3_potential.dat')
     call delete_file('saiii3_slowroll.dat')

  
     n=500


     do i=1,n
        x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)        

        V = saiii3_norm_potential(x,alpha=1._kp,beta=2._kp,mu=1._kp)
        V1 = saiii3_norm_potential(x,alpha=-2._kp,beta=-1.2_kp,mu=1._kp)


        call livewrite('saiii3_potential.dat',x,V,V1)

        eps1 = saiii3_epsilon_one(x,alpha=1._kp,beta=2._kp,mu=40._kp)
        eps2 = saiii3_epsilon_two(x,alpha=1._kp,beta=2._kp,mu=40._kp)
        eps3 = saiii3_epsilon_three(x,alpha=1._kp,beta=2._kp,mu=40._kp)


        call livewrite('saiii3_slowroll.dat',x,eps1,eps2,eps3)

     enddo


     n=250

     call delete_file('saiii3_xV_1.dat')
     call delete_file('saiii3_xeps_1.dat')

     beta = -1.5
     alphamin = -3._kp
     alphamax = saiii_alpha_one(beta)

     
     do i=1,n

        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(n-1,kp) 

        xvzero = saiii_x_potzero(alpha,beta,mu)
        
        call livewrite('saiii3_xV_1.dat',alpha,xvzero(1),xvzero(2))

        xeps = saiii_x_epsoneunity(alpha,beta,mu=10._kp)

        call livewrite('saiii3_xeps_1.dat',alpha,xeps(1),xeps(2))

     enddo


     call delete_file('saiii3_xV_2.dat')
     call delete_file('saiii3_xeps_2.dat')


     beta = 1.5_kp
     alphamin = saiii_alpha_two(beta)
     alphamax = 3._kp

     do i=1,n

        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(n-1,kp) 

        xvzero = saiii_x_potzero(alpha,beta,mu)

        call livewrite('saiii3_xV_2.dat',alpha,xvzero(1),xvzero(2))

        xeps = saiii_x_epsoneunity(alpha,beta,mu=10._kp)
        print *,xeps,saiii3_epsilon_one(xeps(1),alpha,beta,mu=10._kp)
        call livewrite('saiii3_xeps_2.dat',alpha,xeps(1),xeps(2))

     enddo

  endif



!$omp parallel sections &
!$omp default(shared) &
!$omp private(betamin,betamax,alphamin,alphamax) &  
!$omp private(alpha,beta,l,k,i,j,mu) &  
!$omp private(lnRhoRehMax,xend,lnRhoRehMin,lnRhoReh,bfoldstar) &
!$omp private(xstar,w,eps1,eps2,eps3,logErehGeV,Treh,ns,r)


!$omp section
  
  call aspicwrite_header('saiii3p',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','beta '/))
  
  w=0._kp


  betamin = 1.5_kp*beta0
  betamax = 5._kp*betamin

  
  do l=1,Nb
     beta = betamin +  real(l-1,kp)*(betamax - betamin)/real(Nb-1,kp)

     alphamin = 1.5_kp*saiii_alpha_two(beta)
     alphamax = 5._kp*alphamin
     
     do k=1,Na
        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        print *,'alpha= beta= ',alpha,beta
        
        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))
           
           lnRhoRehMin = lnRhoNuc          
           
           xend = saiii3_x_endinf(alpha,beta,mu)

           
           lnRhoRehMax = saiii3_lnrhoreh_max(alpha,beta,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=1,npts

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = saiii3_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = saiii3_epsilon_one(xstar,alpha,beta,mu)
              eps2 = saiii3_epsilon_two(xstar,alpha,beta,mu)
              eps3 = saiii3_epsilon_three(xstar,alpha,beta,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1
         
              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,alpha,beta/))

           end do

        end do

     enddo
  enddo

  call aspicwrite_end()

!$omp section

  call aspicwrite_header('saiii3m',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','beta '/))

  betamax = -1.005_kp 
  betamin = -2._kp


  
  do l=1,Nb
     beta = betamin +  real(l-1,kp)*(betamax - betamin)/real(Nb-1,kp)

     alphamax = 1.02_kp*saiii_alpha_one(beta)

     if (beta.gt.beta1) then       
        alphamin = 0.98_kp*saiii_alpha_three(beta)
     else
        alphamin = 5._kp*alphamax
     endif
     
     do k=1,Na
        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        print *,'alpha= beta= ',alpha,beta
        
        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           lnRhoRehMin = lnRhoNuc

           xend = saiii3_x_endinf(alpha,beta,mu)

           lnRhoRehMax = saiii3_lnrhoreh_max(alpha,beta,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=1,npts

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = saiii3_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = saiii3_epsilon_one(xstar,alpha,beta,mu)
              eps2 = saiii3_epsilon_two(xstar,alpha,beta,mu)
              eps3 = saiii3_epsilon_three(xstar,alpha,beta,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1
           

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

  beta = 2.0_kp
  alpha = 3.0
  mu = 10._kp

  xend = saiii3_x_endinf(alpha,beta,mu)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = saiii3_x_rrad(alpha,beta,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = saiii3_epsilon_one(xstar,alpha,beta,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  saiii3_epsilon_one(xend,alpha,beta,mu)
     VendOverVstar = saiii3_norm_potential(xend,alpha,beta,mu) &
          /saiii3_norm_potential(xstar,alpha,beta,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     
     xstar = saiii3_x_rreh(alpha,beta,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = saiii3_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program saiii3main
