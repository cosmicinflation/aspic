program rclfi3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rclfi3sr, only : rclfi3_epsilon_one, rclfi3_epsilon_two, rclfi3_epsilon_three
  use rclfi3reheat, only : rclfi3_lnrhoreh_max, rclfi3_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use rclfi3sr, only : rclfi3_norm_potential, rclfi3_x_endinf, rclfi3_x_potzero
  use rclfi3sr, only : rclfi3_norm_deriv_potential, rclfi3_norm_deriv_second_potential
  use rclfi3sr, only : rclfi3_numacc_pmin, rclfi3_numacc_alphamin
  use rclfi3reheat, only : rclfi3_x_rreh, rclfi3_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use rclficommon, only : rclfi_alpha_zero
  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n,k,l
  integer :: npts = 15


  integer, parameter :: Na = 2
  real(kp) :: alphamin
  real(kp) :: alphamax

  integer, parameter :: Np = 2 
  real(kp) :: pmin
  real(kp) :: pmax
  
  integer, parameter :: Nmu=60
  real(kp) :: mumin = 0.01_kp
  real(kp) :: mumax = 2000._kp
  
  real(kp) :: xmin = 0
  real(kp) :: xmax = 8
  
  real(kp) :: p,alpha,mu,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax, efoldmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend



  
  Pstar = powerAmpScalar
 

!$omp parallel sections &
!$omp default(shared) &
!$omp private(alphamin,alphamax,pmin,pmax) &  
!$omp private(p,alpha,l,k,i,j,mu,xmin,xmax) &  
!$omp private(lnRhoRehMax,xend,lnRhoRehMin,lnRhoReh,bfoldstar) &
!$omp private(xstar,w,eps1,eps2,eps3,logErehGeV,Treh,ns,r)


!$omp section

  
  call aspicwrite_header('rclfi3pp',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))
  
  w=0._kp

  pmin = 4.1_kp
  pmax = 6._kp
  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)
     
     alphamin = 0.5_kp
     alphamax = 50._kp

     print *,'p alphamin alphamax',p, alphamin,alphamax
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)
        
        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           lnRhoRehMin = lnRhoNuc

           xend = rclfi3_x_endinf(p,alpha,mu)

           lnRhoRehMax = rclfi3_lnrhoreh_max(p,alpha,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi3_x_star(p,alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi3_epsilon_one(xstar,p,alpha,mu)
              eps2 = rclfi3_epsilon_two(xstar,p,alpha,mu)
              eps3 = rclfi3_epsilon_three(xstar,p,alpha,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1


              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,p,alpha/))

           end do

        end do
     enddo
  enddo
  call aspicwrite_end()

!$omp section

  call aspicwrite_header('rclfi3pm',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))

  w=0._kp



  pmin = 1.1_kp*rclfi3_numacc_pmin()
  pmax = 6._kp

  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)

     alphamax = 1.01_kp*rclfi_alpha_zero(p)
     alphamin = max(0.99_kp*rclfi3_numacc_alphamin(p),10._kp*alphamax)
     
     print *,'p alphamin alphamax',p, alphamin,alphamax
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           lnRhoRehMin = lnRhoNuc

           xend = rclfi3_x_endinf(p,alpha,mu)

           lnRhoRehMax = rclfi3_lnrhoreh_max(p,alpha,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi3_x_star(p,alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi3_epsilon_one(xstar,p,alpha,mu)
              eps2 = rclfi3_epsilon_two(xstar,p,alpha,mu)
              eps3 = rclfi3_epsilon_three(xstar,p,alpha,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1


              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,p,alpha/))

           end do

        end do
     enddo
  enddo
  call aspicwrite_end()

!$omp section

  call aspicwrite_header('rclfi3mp',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))

  w=0._kp

  
  pmin = 1.0_kp
  pmax = 3.5_kp

  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)

     alphamin = 1.1_kp*rclfi_alpha_zero(p)
     alphamax = 10._kp*alphamin
     
     print *,'p pmax alphamin alphamax',p, pmax, alphamin,alphamax
     alphamax = min(10._kp*alphamin,alphamax)
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        print *,'p= alpha= ',p,alpha
        

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           lnRhoRehMin = lnRhoNuc

           xend = rclfi3_x_endinf(p,alpha,mu)

           lnRhoRehMax = rclfi3_lnrhoreh_max(p,alpha,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi3_x_star(p,alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi3_epsilon_one(xstar,p,alpha,mu)
              eps2 = rclfi3_epsilon_two(xstar,p,alpha,mu)
              eps3 = rclfi3_epsilon_three(xstar,p,alpha,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('rclfi3_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('rclfi3_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,p,alpha/))

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

  p = 3.0_kp
  alpha = 4._kp
  mu = 10._kp

  xend = rclfi3_x_endinf(p,alpha,mu)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rclfi3_x_rrad(p,alpha,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rclfi3_epsilon_one(xstar,p,alpha,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rclfi3_epsilon_one(xend,p,alpha,mu)
     VendOverVstar = rclfi3_norm_potential(xend,p,alpha,mu) &
          /rclfi3_norm_potential(xstar,p,alpha,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rclfi3_x_rreh(p,alpha,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rclfi3_x_star(p,alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program rclfi3main
