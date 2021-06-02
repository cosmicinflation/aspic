program rclfi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rclfi1sr, only : rclfi1_epsilon_one, rclfi1_epsilon_two, rclfi1_epsilon_three
  use rclfi1reheat, only : rclfi1_lnrhoreh_max, rclfi1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use rclfi1sr, only : rclfi1_norm_potential, rclfi1_x_endinf, rclfi1_x_potmax, rclfi1_numacc_pmax
  use rclfi1sr, only : rclfi1_norm_deriv_potential, rclfi1_norm_deriv_second_potential
  use rclfi1sr, only : rclfi1_numacc_mumin, rclfi1_numacc_efoldmax, rclfi1_numacc_alphamax
  use rclfi1reheat, only : rclfi1_x_rreh, rclfi1_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use rclficommon, only : rclfi_alpha_one, rclfi_alpha_zero
  
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
  
  integer, parameter :: Nmu=100
  real(kp) :: mumin = 0.0001_kp
  real(kp) :: mumax = 1000._kp
  
  real(kp) :: xmin = 1e-6
  real(kp) :: xmax = 2
  
  real(kp) :: p,alpha,mu,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax, efoldmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: V1,x,V,V2


  
  Pstar = powerAmpScalar


  n=1000

  call delete_file('rclfi_alpha_zero.dat')
  call delete_file('rclfi_alpha_one.dat')

  pmin = -10
  pmax = 28
  
  do i=1,n
     p = pmin + + real(i-1,kp)*(pmax-pmin)/real(n-1,kp)        

     call livewrite('rclfi_alpha_zero.dat',p,rclfi_alpha_zero(p))
     call livewrite('rclfi_alpha_one.dat',p,rclfi_alpha_one(p))

  end do

  
  p=1.5_kp
  alpha = 8._kp
  call delete_file('rclfi1_potential.dat')
  call delete_file('rclfi1_slowroll.dat')
  
  do i=1,n
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)        

     V = rclfi1_norm_potential(x,alpha=p,alpha=p,mu=1._kp)
     
     V1 = rclfi1_norm_potential(x,alpha=-1.2_kp,p=-0.5_kp,mu=1._kp)


     call livewrite('rclfi1_potential.dat',x,V,V1)

     eps1 = rclfi1_epsilon_one(x,alpha=p,alpha=p,mu=100._kp)
     eps2 = rclfi1_epsilon_two(x,alpha=p,alpha=p,mu=100._kp)
     eps3 = rclfi1_epsilon_three(x,alpha=p,alpha=p,mu=100._kp)

     call livewrite('rclfi1_slowroll.dat',x,eps1,eps2,eps3)

  enddo

  
  call aspicwrite_header('rclfi1pm',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))
  
  w=0._kp

  pmin = 4.1_kp
  pmax = 6._kp
  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)
     
     alphamax = 1.01_kp*rclfi_alpha_one(p)


     alphamin = alphamax*10._kp

     print *,'p alphamin alphamax',p, alphamin,alphamax
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)

        mumin = rclfi1_numacc_mumin(120._kp,p,alpha)

        if (mumin.gt.mumax) cycle
        
        print *,'mumin= ',mumin
        
        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = rclfi1_numacc_efoldmax(p,alpha,mu)
           print *,'p,alpha,mu,efoldNumAccMax= ',p,alpha,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = rclfi1_x_endinf(p,alpha,mu)

           lnRhoRehMax = rclfi1_lnrhoreh_max(p,alpha,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi1_x_star(p,alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi1_epsilon_one(xstar,p,alpha,mu)
              eps2 = rclfi1_epsilon_two(xstar,p,alpha,mu)
              eps3 = rclfi1_epsilon_three(xstar,p,alpha,mu)


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

  
  call aspicwrite_header('rclfi1mm',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))

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

        mumin = rclfi1_numacc_mumin(120._kp,p,alpha)
        if (mumin.gt.mumax) cycle
        print *,'mumin= ',mumin

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = rclfi1_numacc_efoldmax(p,alpha,mu)
           print *,'p,alpha,mu,efoldNumAccMax= ',p,alpha,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = rclfi1_x_endinf(p,alpha,mu)

           lnRhoRehMax = rclfi1_lnrhoreh_max(p,alpha,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi1_x_star(p,alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi1_epsilon_one(xstar,p,alpha,mu)
              eps2 = rclfi1_epsilon_two(xstar,p,alpha,mu)
              eps3 = rclfi1_epsilon_three(xstar,p,alpha,mu)


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

  
  call aspicwrite_header('rclfi1mp',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))

  w=0._kp


  pmin = 0.5_kp
  pmax = 0.5_kp*rclfi1_numacc_pmax()
  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)

     alphamin = 1.01_kp*rclfi_alpha_one(p)
     alphamax = rclfi1_numacc_alphamax(p)
     
     print *,'p pmax alphamin alphamax',p, pmax, alphamin,alphamax
     alphamax = min(10._kp*alphamin,alphamax)

     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)
        
        mumin = rclfi1_numacc_mumin(120._kp,p,alpha)
        if (mumin.gt.mumax) cycle
        print *,'mumin= ',mumin,mumax

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = rclfi1_numacc_efoldmax(p,alpha,mu)
           print *,'p,alpha,mu,efoldNumAccMax= ',p,alpha,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = rclfi1_x_endinf(p,alpha,mu)

           lnRhoRehMax = rclfi1_lnrhoreh_max(p,alpha,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi1_x_star(p,alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi1_epsilon_one(xstar,p,alpha,mu)
              eps2 = rclfi1_epsilon_two(xstar,p,alpha,mu)
              eps3 = rclfi1_epsilon_three(xstar,p,alpha,mu)


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
 

  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'


  
  lnRradmin=-42
  lnRradmax = 10

  p = 0.9_kp
  alpha = -0.7
  mu = 10._kp*rclfi1_numacc_mumin(120._kp,p,alpha)

  xend = rclfi1_x_endinf(p,alpha,mu)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rclfi1_x_rrad(p,alpha,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rclfi1_epsilon_one(xstar,p,alpha,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rclfi1_epsilon_one(xend,p,alpha,mu)
     VendOverVstar = rclfi1_norm_potential(xend,p,alpha,mu) &
          /rclfi1_norm_potential(xstar,p,alpha,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rclfi1_x_rreh(p,alpha,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rclfi1_x_star(p,alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program rclfi1main
