program rclfi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rclfi2sr, only : rclfi2_epsilon_one, rclfi2_epsilon_two, rclfi2_epsilon_three
  use rclfi2reheat, only : rclfi2_lnrhoreh_max, rclfi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use rclfi2sr, only : rclfi2_norm_potential, rclfi2_x_endinf, rclfi2_x_potmax
  use rclfi2sr, only : rclfi2_norm_deriv_potential, rclfi2_norm_deriv_second_potential
  use rclfi2sr, only : rclfi2_numacc_mumin, rclfi2_numacc_efoldmax, rclfi2_numacc_pmax
  use rclfi2sr, only : rclfi2_numacc_alphamax
  use rclfi2reheat, only : rclfi2_x_rreh, rclfi2_x_rrad
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
  
  integer, parameter :: Nmu=50
  real(kp) :: mumin = 0.0001_kp
  real(kp) :: mumax = 1000._kp
  
  real(kp) :: xmin = 0
  real(kp) :: xmax = 8
  
  real(kp) :: alpha,p,mu,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax, efoldmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend



  
  Pstar = powerAmpScalar
 
  call aspicwrite_header('rclfi2pm',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))
  
  w=0._kp

  pmin = 4.1_kp
  pmax = 6._kp
  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)
     
     alphamax = 1.01_kp*rclfi_alpha_zero(p)


     alphamin = alphamax*10._kp

     print *,'p alphamin alphamax',p, alphamin,alphamax
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        mumin = rclfi2_numacc_mumin(120._kp,alpha,p)

        print *,'mumin= ',mumin
        if (mumin.gt.mumax) cycle
        
        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = rclfi2_numacc_efoldmax(alpha,p,mu)
           print *,'alpha,p,mu,efoldNumAccMax= ',alpha,p,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = rclfi2_x_endinf(alpha,p,mu)

           lnRhoRehMax = rclfi2_lnrhoreh_max(alpha,p,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi2_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi2_epsilon_one(xstar,alpha,p,mu)
              eps2 = rclfi2_epsilon_two(xstar,alpha,p,mu)
              eps3 = rclfi2_epsilon_three(xstar,alpha,p,mu)


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

  
  call aspicwrite_header('rclfi2mm',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))

  w=0._kp


  pmin = 0.1_kp
  pmax = 3.9_kp

  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)

     alphamax = -0.2_kp
     alphamin = alphamax*10._kp

     print *,'p alphamin alphamax',p, alphamin,alphamax

     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        mumin = rclfi2_numacc_mumin(120._kp,alpha,p)
        if (mumin.gt.mumax) cycle
        
        print *,'mumin= ',mumin

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = rclfi2_numacc_efoldmax(alpha,p,mu)
           print *,'alpha,p,mu,efoldNumAccMax= ',alpha,p,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = rclfi2_x_endinf(alpha,p,mu)

           lnRhoRehMax = rclfi2_lnrhoreh_max(alpha,p,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi2_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi2_epsilon_one(xstar,alpha,p,mu)
              eps2 = rclfi2_epsilon_two(xstar,alpha,p,mu)
              eps3 = rclfi2_epsilon_three(xstar,alpha,p,mu)


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

  
  call aspicwrite_header('rclfi2mp',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))

  w=0._kp

  
  pmin = 1.0_kp
  pmax = rclfi2_numacc_pmax()

  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)

     alphamin = 1.01_kp*rclfi_alpha_zero(p)
     alphamax = rclfi2_numacc_alphamax(p)
     
     print *,'p pmax alphamin alphamax',p, pmax, alphamin,alphamax
     alphamax = min(10._kp*alphamin,alphamax)
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        print *,'p= alpha= ',p,alpha
        
        mumin = rclfi2_numacc_mumin(120._kp,alpha,p)
        if (mumin.gt.mumax) cycle
        
        print *,'mumin= ',mumin

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = rclfi2_numacc_efoldmax(alpha,p,mu)
           print *,'alpha,p,mu,efoldNumAccMax= ',alpha,p,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = rclfi2_x_endinf(alpha,p,mu)

           lnRhoRehMax = rclfi2_lnrhoreh_max(alpha,p,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi2_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi2_epsilon_one(xstar,alpha,p,mu)
              eps2 = rclfi2_epsilon_two(xstar,alpha,p,mu)
              eps3 = rclfi2_epsilon_three(xstar,alpha,p,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('rclfi2_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('rclfi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,alpha,p/))

           end do

        end do
     enddo
  enddo
  call aspicwrite_end()
 

  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'


  
  lnRradmin=-42
  lnRradmax = 10

  p = 0.9_kp
  alpha = -0.7
  mu = 10._kp*rclfi2_numacc_mumin(120._kp,alpha,p)

  xend = rclfi2_x_endinf(alpha,p,mu)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rclfi2_x_rrad(alpha,p,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rclfi2_epsilon_one(xstar,alpha,p,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rclfi2_epsilon_one(xend,alpha,p,mu)
     VendOverVstar = rclfi2_norm_potential(xend,alpha,p,mu) &
          /rclfi2_norm_potential(xstar,alpha,p,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rclfi2_x_rreh(alpha,p,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rclfi2_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program rclfi2main
