program rcpimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar

  use srreheat, only : log_energy_reheat_ingev
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat 
  
  use rcpisr, only : rcpi_epsilon_one, rcpi_epsilon_two, rcpi_epsilon_three
  use rcpisr, only : rcpi_norm_potential, rcpi_x_endinf, rcpi_efoldmax
  use rcpisr, only : rcpi_numacc_alphamin, rcpi_numacc_alphamax, rcpi_numacc_betamin
  
  use rcpireheat, only : rcpi_lnrhoreh_max, rcpi_x_star
  use rcpireheat, only : rcpi_x_rreh, rcpi_x_rrad

  use srflow, only : scalar_spectral_index, tensor_to_scalar_ratio

  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  
  implicit none
  
  real(kp) :: p,alpha,beta
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r
  real(kp) :: logErehGeV, Pstar, bfoldstar
  real(kp) :: lnRhoRehMin, lnRhoRehMax, w, efoldmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  integer :: i,j,k
  integer :: nalp,nbet,npts
  real(kp) :: betamin, betamax
  real(kp) :: alphamin, alphamax

  integer, parameter :: nvec = 3
  real(kp), dimension(nvec) :: betavec
  
  logical, parameter :: display = .true.
  real(kp), parameter :: efoldNum = 60._kp

  real(kp) :: a,b,mu
  
  w = 0._kp
  Pstar = powerAmpScalar 
  

  nalp = 500
  npts = 15

  
  
  call aspicwrite_header('rcpi',labeps12,labnsr,labbfoldreh,(/'alpha','beta ','p    '/))

 
  p=2._kp

  lnRhoRehMin = lnRhoNuc

  betamin = rcpi_numacc_betamin(p)

  betavec = (/betamin,betamin+0.05, betamin+0.1/)
  
  do j=1,nvec
     beta = betavec(j)

     alphamin = rcpi_numacc_alphamin(p,beta)
     alphamax = rcpi_numacc_alphamax(p,beta)


     print *,'alphamin, alphamax= ',beta, alphamin, alphamax
     
     do i=1,nalp
        alpha = alphamin + (i-1)*(alphamax-alphamin)/real(nalp-1,kp)

        efoldmax = rcpi_efoldmax(p,alpha,beta)
        

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif

        print *,' p alpha beta= ',p,alpha,beta

        xEnd = rcpi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcpi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcpi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcpi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcpi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcpi_epsilon_three(xstar,p,alpha,beta)                     

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1/))


           
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta,p/))
           
           if ((eps1.gt.0.01).or.(eps2.lt.-0.08).or.(eps2.gt.0.08)) cycle


        enddo

     enddo
  enddo


  p=4._kp

  betamin = rcpi_numacc_betamin(p)
  betavec = (/betamin,betamin+0.1, betamin+0.2/)

  lnRhoRehMin = lnRhoNuc
  
  do j=1,nvec
     beta = betavec(j)

     alphamin = rcpi_numacc_alphamin(p,beta)
     alphamax = rcpi_numacc_alphamax(p,beta)

     print *,'beta, alphamin, alphamax= ',beta, alphamin, alphamax

     do i=1,nalp
        alpha = alphamin + (i-1)*(alphamax-alphamin)/real(nalp-1,kp)
        
        efoldmax = rcpi_efoldmax(p,alpha,beta)

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,' p alpha beta= ',p,alpha,beta
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif

        print *,' p alpha beta= ',p,alpha,beta

        xEnd = rcpi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcpi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcpi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcpi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcpi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcpi_epsilon_three(xstar,p,alpha,beta)
          
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1,eps2/))

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta,p/))
           
           if ((eps1.gt.0.01).or.(eps2.lt.-0.08).or.(eps2.gt.0.08)) cycle


        enddo
     enddo
  enddo

  call aspicwrite_end()


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10

  p=2._kp
  alpha = -0.1_kp
  beta = 0.5_kp

  xEnd = rcpi_x_endinf(p,alpha,beta)

  if (rcpi_efoldmax(p,alpha,beta).lt.efoldNum) then
     stop 'efoldmax too small!'
  endif
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rcpi_x_rrad(p,alpha,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rcpi_epsilon_one(xstar,p,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rcpi_epsilon_one(xend,p,alpha,beta)
     VendOverVstar = rcpi_norm_potential(xend,p,alpha,beta)/rcpi_norm_potential(xstar,p,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rcpi_x_rreh(p,alpha,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rcpi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

  
end program rcpimain
