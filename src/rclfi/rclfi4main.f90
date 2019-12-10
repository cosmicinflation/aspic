program rclfi4main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rclfi4sr, only : rclfi4_epsilon_one, rclfi4_epsilon_two, rclfi4_epsilon_three
  use rclfi4reheat, only : rclfi4_lnrhoreh_max, rclfi4_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use rclfi4sr, only : rclfi4_norm_potential, rclfi4_x_endinf
  use rclfi4sr, only : rclfi4_norm_deriv_potential, rclfi4_norm_deriv_second_potential
  use rclfi4reheat, only : rclfi4_x_rreh, rclfi4_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use rclficommon, only : rclfi_alpha_one
  
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
  real(kp) :: mumax = 10000._kp
  
  real(kp) :: xmin = 0
  real(kp) :: xmax = 8
  
  real(kp) :: alpha,p,mu,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax, efoldmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend



  
  Pstar = powerAmpScalar
 

!$omp parallel sections &
!$omp default(shared) &
!$omp private(alphamin,alphamax,pmin,pmax) &  
!$omp private(alpha,p,l,k,i,j,mu,xmin,xmax) &  
!$omp private(lnRhoRehMax,xend,lnRhoRehMin,lnRhoReh,bfoldstar) &
!$omp private(xstar,w,eps1,eps2,eps3,logErehGeV,Treh,ns,r)


!$omp section

  
  
  call aspicwrite_header('rclfi4p',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))
  
  w=0._kp

  pmin = 4.5_kp
  pmax = 6._kp
  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)
     
     alphamin = 0.99_kp*rclfi_alpha_one(p)
     alphamax = -0.8_kp

     print *,'p alphamin alphamax',p, alphamin,alphamax
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)
        
        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           lnRhoRehMin = lnRhoNuc

           xend = rclfi4_x_endinf(alpha,p,mu)

           print *,'xend',xend
           
           lnRhoRehMax = rclfi4_lnrhoreh_max(alpha,p,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi4_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi4_epsilon_one(xstar,alpha,p,mu)
              eps2 = rclfi4_epsilon_two(xstar,alpha,p,mu)
              eps3 = rclfi4_epsilon_three(xstar,alpha,p,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1


              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,alpha,p/))

           end do

        end do
     enddo
  enddo
  call aspicwrite_end()

!$omp section
  
  call aspicwrite_header('rclfi4m',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))

  w=0._kp



  pmin = 0.1_kp
  pmax = 3.0_kp

  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)

     alphamin = 0.3_kp
     alphamax = 0.99_kp*rclfi_alpha_one(p)

     
     print *,'p alphamin alphamax',p, alphamin,alphamax
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           lnRhoRehMin = lnRhoNuc

           xend = rclfi4_x_endinf(alpha,p,mu)

           lnRhoRehMax = rclfi4_lnrhoreh_max(alpha,p,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi4_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi4_epsilon_one(xstar,alpha,p,mu)
              eps2 = rclfi4_epsilon_two(xstar,alpha,p,mu)
              eps3 = rclfi4_epsilon_three(xstar,alpha,p,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,alpha,p/))

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
  alpha = 2._kp
  mu = 10._kp

  xend = rclfi4_x_endinf(alpha,p,mu)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rclfi4_x_rrad(alpha,p,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rclfi4_epsilon_one(xstar,alpha,p,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rclfi4_epsilon_one(xend,alpha,p,mu)
     VendOverVstar = rclfi4_norm_potential(xend,alpha,p,mu) &
          /rclfi4_norm_potential(xstar,alpha,p,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rclfi4_x_rreh(alpha,p,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rclfi4_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program rclfi4main
