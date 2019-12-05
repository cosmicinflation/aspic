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

  use rclficommon, only : rclfi_alpha_one
  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n,k,l
  integer :: npts = 15


  integer, parameter :: Na = 3
  real(kp) :: alphamin
  real(kp) :: alphamax

  integer, parameter :: Np = 2
  real(kp) :: pmin
  real(kp) :: pmax
  
  integer, parameter :: Nmu=10
  real(kp) :: mumin = 0.1_kp
  real(kp) :: mumax = 10000._kp
  
  real(kp) :: xmin = 0
  real(kp) :: xmax = 8
  
  real(kp) :: alpha,p,mu,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax, efoldmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: V1,x,V,V2


  
  Pstar = powerAmpScalar


  call delete_file('rclfi1_potential.dat')
  call delete_file('rclfi1_slowroll.dat')

  n=250

  goto 10
  
  do i=1,n
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)        

     V = rclfi1_norm_potential(x,alpha=0.8_kp,p=0.1_kp,mu=1._kp)
     
     V1 = rclfi1_norm_potential(x,alpha=-1.2_kp,p=-0.5_kp,mu=1._kp)


     call livewrite('rclfi1_potential.dat',x,V,V1)

     eps1 = rclfi1_epsilon_one(x,alpha=0.8_kp,p=0.1_kp,mu=20._kp)
     eps2 = rclfi1_epsilon_two(x,alpha=0.8_kp,p=0.1_kp,mu=20._kp)
     eps3 = rclfi1_epsilon_three(x,alpha=0.8_kp,p=0.1_kp,mu=20._kp)

     call livewrite('rclfi1_slowroll.dat',x,eps1,eps2,eps3)

  enddo

10 continue
 
  call delete_file('rclfi1_predic.dat')
  call delete_file('rclfi1_nsr.dat')

  
  call aspicwrite_header('rclfi1_pgt4_an',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))
  
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

        mumin = rclfi1_numacc_mumin(120._kp,alpha,p)

        print *,'mumin= ',mumin
        
        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = rclfi1_numacc_efoldmax(alpha,p,mu)
           print *,'alpha,p,mu,efoldNumAccMax= ',alpha,p,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = rclfi1_x_endinf(alpha,p,mu)

           lnRhoRehMax = rclfi1_lnrhoreh_max(alpha,p,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi1_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi1_epsilon_one(xstar,alpha,p,mu)
              eps2 = rclfi1_epsilon_two(xstar,alpha,p,mu)
              eps3 = rclfi1_epsilon_three(xstar,alpha,p,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('rclfi1_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('rclfi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,alpha,p/))

           end do

        end do
     enddo
  enddo
  call aspicwrite_end()


  call aspicwrite_header('rclfi1_plt4_an',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))

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

        mumin = rclfi1_numacc_mumin(120._kp,alpha,p)

        print *,'mumin= ',mumin

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = rclfi1_numacc_efoldmax(alpha,p,mu)
           print *,'alpha,p,mu,efoldNumAccMax= ',alpha,p,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = rclfi1_x_endinf(alpha,p,mu)

           lnRhoRehMax = rclfi1_lnrhoreh_max(alpha,p,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi1_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi1_epsilon_one(xstar,alpha,p,mu)
              eps2 = rclfi1_epsilon_two(xstar,alpha,p,mu)
              eps3 = rclfi1_epsilon_three(xstar,alpha,p,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('rclfi1_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('rclfi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/) &
                   ,(/mu,alpha,p/))

           end do

        end do
     enddo
  enddo
  call aspicwrite_end()

  
  call aspicwrite_header('rclfi1_plt4_ap',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha','p    '/))

  w=0._kp


  pmin = 3.0_kp
  pmax = 0.99_kp*rclfi1_numacc_pmax()

  
  
  do l=1,Np
     p = pmin +  real(l-1,kp)*(pmax - pmin)/real(Np-1,kp)

     alphamin = 1.01_kp*rclfi_alpha_one(p)
     alphamax = rclfi1_numacc_alphamax(p)
     
     print *,'p pmax alphamin alphamax',p, pmax, alphamin,alphamax
     alphamax = min(10._kp*alphamin,alphamax)
     
     do k=1,Na

        alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)
        
        mumin = rclfi1_numacc_mumin(120._kp,alpha,p)

        print *,'mumin= ',mumin

        do j=0,Nmu 

           mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

           efoldmax = rclfi1_numacc_efoldmax(alpha,p,mu)
           print *,'alpha,p,mu,efoldNumAccMax= ',alpha,p,mu,efoldmax

           lnRhoRehMin = lnRhoNuc

           xend = rclfi1_x_endinf(alpha,p,mu)

           lnRhoRehMax = rclfi1_lnrhoreh_max(alpha,p,mu,xend,Pstar)

           print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=npts,1,-1

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = rclfi1_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

              eps1 = rclfi1_epsilon_one(xstar,alpha,p,mu)
              eps2 = rclfi1_epsilon_two(xstar,alpha,p,mu)
              eps3 = rclfi1_epsilon_three(xstar,alpha,p,mu)


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)


              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('rclfi1_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('rclfi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

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
  mu = 10._kp*rclfi1_numacc_mumin(120._kp,alpha,p)

  xend = rclfi1_x_endinf(alpha,p,mu)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rclfi1_x_rrad(alpha,p,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rclfi1_epsilon_one(xstar,alpha,p,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rclfi1_epsilon_one(xend,alpha,p,mu)
     VendOverVstar = rclfi1_norm_potential(xend,alpha,p,mu) &
          /rclfi1_norm_potential(xstar,alpha,p,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rclfi1_x_rreh(alpha,p,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rclfi1_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program rclfi1main
