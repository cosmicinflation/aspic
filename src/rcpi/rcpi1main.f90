program rcpi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar

  use srreheat, only : log_energy_reheat_ingev
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : labeps12, labnsr, labbfoldreh
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  
  use rcpi1sr, only : rcpi1_epsilon_one, rcpi1_epsilon_two, rcpi1_epsilon_three
  use rcpi1sr, only : rcpi1_norm_potential, rcpi1_x_endinf, rcpi1_efoldmax
  use rcpi1reheat, only : rcpi1_lnrhoreh_max, rcpi1_x_star
  use rcpi1reheat, only : rcpi1_x_rreh, rcpi1_x_rrad

  use srflow, only : scalar_spectral_index, tensor_to_scalar_ratio
  
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

  logical, parameter :: display = .true.
  real(kp), parameter :: efoldNum = 60._kp
  
  w = 0._kp
  Pstar = powerAmpScalar
  
  nbet = 10
  

  nalp = 10

  npts = 20

  call aspicwrite_header('rcpi1',labeps12,labnsr,labbfoldreh,(/'alpha','beta ','p    '/))

  p=2._kp
  betamin = 0.1_kp
  betamax = 4._kp
  
  lnRhoRehMin = lnRhoNuc
  
  do j=1,nbet
     beta = betamin + (j-1)*(betamax-betamin)/real(nbet-1,kp)

     alphamin = -2*sqrt(beta)*(1._kp - 1._kp/real(nalp,kp))
     alphamax = 2*sqrt(beta)*(1._kp - 1._kp/real(nalp,kp))

     do i=1,nalp
        alpha = alphamin + (i-1)*(alphamax-alphamin)/real(nalp-1,kp)

        efoldmax = rcpi1_efoldmax(p,alpha,beta)
        
        print *,' p alpha beta= ',p,alpha,beta
        print *,'efoldmax =     ',efoldmax

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           write(*,*)
           cycle
        endif

        lnRhoRehMax = rcpi1_lnrhoreh_max(p,alpha,beta,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcpi1_x_star(p,alpha,beta,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcpi1_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcpi1_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcpi1_epsilon_three(xstar,p,alpha,beta)
          
           if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar)

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1,eps2/))

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta,p/))

        enddo

     enddo
  enddo


  p=4._kp
  betamin = 0.1_kp
  betamax = 4._kp
  
  lnRhoRehMin = lnRhoNuc
  
  do j=1,nbet
     beta = betamin + (j-1)*(betamax-betamin)/real(nbet-1,kp)

     alphamin = -2*sqrt(beta)*(1._kp - 1._kp/real(nalp,kp))
     alphamax = 2*sqrt(beta)*(1._kp - 1._kp/real(nalp,kp))

     do i=1,nalp
        alpha = alphamin + (i-1)*(alphamax-alphamin)/real(nalp-1,kp)

        efoldmax = rcpi1_efoldmax(p,alpha,beta)

        print *,' p alpha beta= ',p,alpha,beta
        print *,'efoldmax =     ',efoldmax

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           write(*,*)
           cycle
        endif
        
        lnRhoRehMax = rcpi1_lnrhoreh_max(p,alpha,beta,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcpi1_x_star(p,alpha,beta,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcpi1_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcpi1_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcpi1_epsilon_three(xstar,p,alpha,beta)
          
           if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar)

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1,eps2/))

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta,p/))

        enddo
     enddo
  enddo

  call aspicwrite_end()

  
end program rcpi1main
