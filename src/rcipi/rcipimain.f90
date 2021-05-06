program rcipimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar

  use srreheat, only : log_energy_reheat_ingev
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat 
  
  use rcipisr, only : rcipi_epsilon_one, rcipi_epsilon_two, rcipi_epsilon_three
  use rcipisr, only : rcipi_norm_potential, rcipi_x_endinf, rcipi_efoldmax
  use rcipisr, only : rcipi_numacc_alphamin, rcipi_numacc_alphamax, rcipi_numacc_betamin
  use rcipisr, only : rcipi_alpha_zero
  
  use rcipireheat, only : rcipi_lnrhoreh_max, rcipi_x_star
  use rcipireheat, only : rcipi_x_rreh, rcipi_x_rrad

  use srflow, only : scalar_spectral_index, tensor_to_scalar_ratio

  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  use infinout, only : livewrite, delete_file

  
  implicit none
  
  real(kp) :: p,alpha,beta
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r
  real(kp) :: logErehGeV, Pstar, bfoldstar
  real(kp) :: lnRhoRehMin, lnRhoRehMax, w, efoldmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  integer :: i,j,k,n
  integer :: nalp,nbet,npts
  real(kp) :: betamin, betamax
  real(kp) :: alphamin, alphamax, alphazero
  real(kp) :: alphanminmax

  integer, parameter :: nvec = 3
  real(kp), dimension(nvec) :: betavec
  
  logical, parameter :: display = .true.
  logical, parameter :: generateAllData = .true.
  
  real(kp), parameter :: efoldNum = 60._kp

  real(kp) :: x,xmin,xmax,V1,V2
  real(kp) :: eps1b, eps2b, eps3b
  
  w = 0._kp
  Pstar = powerAmpScalar 
   
  
  call delete_file('rcipi_alphazero.dat')
  call delete_file('rcipi_stability.dat')
  
  n=400
  xmin = 0._kp
  xmax = 6._kp

  p=2._kp
  
  do i=1,n

     beta = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)

     if (beta.le.p*p) then
        alpha = rcipi_alpha_zero(p,beta)     
        call livewrite('rcipi_alphazero.dat',beta,alpha,-alpha)
     endif

     alpha = 2._kp*sqrt(beta)
     
     call livewrite('rcipi_stability.dat',beta,alpha,-alpha)
        

  enddo

  
  call delete_file('rcipi1_potential.dat')
  call delete_file('rcipi1_slowroll.dat')

  
  n=1000

  xmin = 0.01_kp
  xmax = 3._kp
  
  do i=1,n
     
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)
     V1 = rcipi_norm_potential(x,p=2._kp,alpha=-2._kp,beta=2._kp)
     V2 = rcipi_norm_potential(x,p=2._kp,alpha=-2.5_kp,beta=2._kp)

     call livewrite('rcipi1_potential.dat',x,V1,V2)
     
     eps1 = rcipi_epsilon_one(x,p=2._kp,alpha=-2._kp,beta=2._kp)
     eps2 = rcipi_epsilon_two(x,p=2._kp,alpha=-2._kp,beta=2._kp)
     eps3 = rcipi_epsilon_three(x,p=2._kp,alpha=-2._kp,beta=2._kp)

     eps1b = rcipi_epsilon_one(x,p=2._kp,alpha=-2.5_kp,beta=2._kp)
     eps2b = rcipi_epsilon_two(x,p=2._kp,alpha=-2.5_kp,beta=2._kp)
     eps3b = rcipi_epsilon_three(x,p=2._kp,alpha=-2.5_kp,beta=2._kp)
     
     call livewrite('rcipi1_slowroll.dat',x,eps1,eps2,eps3,eps1b,eps2b,eps3b)
     
  enddo


  call delete_file('rcipi2_potential.dat')
  call delete_file('rcipi2_slowroll.dat')

  
  n=1000

  xmin = 0.01_kp
  xmax = 10._kp
  
  do i=1,n
     
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)
     V1 = rcipi_norm_potential(x,p=2._kp,alpha=0.5_kp,beta=2._kp)

     call livewrite('rcipi2_potential.dat',x,V1)
     
     eps1 = rcipi_epsilon_one(x,p=2._kp,alpha=0.5_kp,beta=2._kp)
     eps2 = rcipi_epsilon_two(x,p=2._kp,alpha=0.5_kp,beta=2._kp)
     eps3 = rcipi_epsilon_three(x,p=2._kp,alpha=0.5_kp,beta=2._kp)
     
     call livewrite('rcipi2_slowroll.dat',x,eps1,eps2,eps3)
     
  enddo



  if (generateAllData) then
  

!$omp parallel sections &
!$omp default(shared) &
!$omp private(p,alpha,beta,alphanminmax,alphazero,i,j,k,efoldmax) &
!$omp private(npts,nalp,nbet,betamin,betamax,alphamin,alphamax) &
!$omp private(xstar,w,eps1,eps2,eps3,logErehGeV,ns,r) &
!$omp private(lnRhoRehMax,xend,lnRhoRehMin,lnRhoReh,bfoldstar)


  

!$omp section

  call aspicwrite_header('rcipi1tune2',labeps12,labnsr,labbfoldreh,(/'alphapo','beta   ','p      '/))

 
  p=2._kp
  nbet = 3
  nalp = 500
  npts=30
  
  lnRhoRehMin = lnRhoNuc

  betamin = 1.0_kp
  betamax = 2.0_kp
  
  do j=1,nbet
     beta = betamin + (betamax-betamin)*(real(j-1,kp))/real(nbet-1,kp)

     alphazero = rcipi_alpha_zero(p,beta)
     
     alphamax = -alphazero - 0.01_kp
     alphamin = -alphazero
    
     
     do i=1,nalp
        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)
        
        efoldmax = rcipi_efoldmax(p,alpha,beta)
        

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif

        xEnd = rcipi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1/))


           
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           
           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha+alphazero,beta,p/))
           
        enddo

     enddo
  enddo


  call aspicwrite_end()

!$omp section

  call aspicwrite_header('rcipi1tune4',labeps12,labnsr,labbfoldreh,(/'alphapo','beta   ','p      '/))

 
  p=4._kp
  nbet = 3
  nalp = 500
  npts=30
  
  lnRhoRehMin = lnRhoNuc

  betamin = 7._kp
  betamax = 9._kp
  
  do j=1,nbet
     beta = betamin + (betamax-betamin)*(real(j-1,kp))/real(nbet-1,kp)

     alphazero = rcipi_alpha_zero(p,beta)
     
     alphamax = -alphazero - 0.0001_kp
     alphamin = -alphazero
    
     
     do i=1,nalp
        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)

        
        efoldmax = rcipi_efoldmax(p,alpha,beta)
        

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif


        xEnd = rcipi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1/))


           
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           
           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha+alphazero,beta,p/))
           
        enddo

     enddo
  enddo


  call aspicwrite_end()

!$omp section
  

  call aspicwrite_header('rcipi1',labeps12,labnsr,labbfoldreh,(/'alphanminmax','beta        ','p           '/))

 
  p=2._kp
  nbet = 3
  nalp = 500
  npts=15
  
  lnRhoRehMin = lnRhoNuc

  betamin = 0.06_kp
  betamax = 0.10_kp
  
  do j=1,nbet
     beta = betamin + (betamax-betamin)*(real(j-1,kp))/real(nbet-1,kp)

     alphazero = rcipi_alpha_zero(p,beta)
     

     alphamin = -2._kp*sqrt(beta) + epsilon(1._kp)
     alphamax = -alphazero - epsilon(1._kp)

!     alphamin = alphazero + epsilon(1._kp)
!     alphamax = 2._kp*sqrt(beta) - epsilon(1._kp)
    
     
     do i=1,nalp
        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)

        alphanminmax = (alpha-alphamin)/(alphamax-alphamin)
        
        efoldmax = rcipi_efoldmax(p,alpha,beta)
        

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif


        xEnd = rcipi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1/))

           
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           
           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alphanminmax,beta,p/))
           
        enddo

     enddo
  enddo


  p=3._kp
  nbet = 3
  nalp = 500
  npts=15
  
  lnRhoRehMin = lnRhoNuc

  betamin = 0.06_kp
  betamax = 0.10_kp
  
  
  do j=1,nbet
     beta = betamin + (betamax-betamin)*(real(j-1,kp))/real(nbet-1,kp)

     alphazero = rcipi_alpha_zero(p,beta)
     

     alphamin = -2._kp*sqrt(beta) + epsilon(1._kp)
     alphamax = -alphazero - epsilon(1._kp)

!     alphamin = alphazero + epsilon(1._kp)
!     alphamax = 2._kp*sqrt(beta) - epsilon(1._kp)
     
     
     do i=1,nalp
        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)

        alphanminmax = (alpha-alphamin)/(alphamax-alphamin)
        
        efoldmax = rcipi_efoldmax(p,alpha,beta)
        

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif


        xEnd = rcipi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1/))

           
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           
           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alphanminmax,beta,p/))
           
        enddo

     enddo
  enddo


  p=4._kp
  nbet = 3
  nalp = 500
  npts=15
  
  lnRhoRehMin = lnRhoNuc

  betamin = 0.06_kp
  betamax = 0.10_kp
 
  
  do j=1,nbet
     beta = betamin + (betamax-betamin)*(real(j-1,kp))/real(nbet-1,kp)

     alphazero = rcipi_alpha_zero(p,beta)
     

     alphamin = -2._kp*sqrt(beta) + epsilon(1._kp)
     alphamax = -alphazero - epsilon(1._kp)

!     alphamin = alphazero + epsilon(1._kp)
!     alphamax = 2._kp*sqrt(beta) - epsilon(1._kp)
     
     
     do i=1,nalp
        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)

        alphanminmax = (alpha-alphamin)/(alphamax-alphamin)
        
        efoldmax = rcipi_efoldmax(p,alpha,beta)
        

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif

        xEnd = rcipi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1/))

           
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           
           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alphanminmax,beta,p/))
           
        enddo

     enddo
  enddo
  

  call aspicwrite_end()


!$omp section

  call aspicwrite_header('rcipi2tune2',labeps12,labnsr,labbfoldreh,(/'alpha','beta ','p    '/))

 
  p=2._kp
  nalp = 20
  npts=15
  
  lnRhoRehMin = lnRhoNuc
  
  beta = 2._kp


  alphazero = rcipi_alpha_zero(p,beta)


  alphamin = -alphazero*0.99
  alphamax = alphazero*0.99

  !     alphamin = alphazero + epsilon(1._kp)
  !     alphamax = 2._kp*sqrt(beta) - epsilon(1._kp)


  do i=1,nalp
     alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)

     alphanminmax = (alpha-alphamin)/(alphamax-alphamin)

     efoldmax = rcipi_efoldmax(p,alpha,beta)


     if (efoldmax.lt.efoldNum) then
        write(*,*)'efoldmax too small!'
        print *,'efoldmax =     ',efoldmax
        write(*,*)
        cycle
     endif


     xEnd = rcipi_x_endinf(p,alpha,beta)
     lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

     do k=1,npts
        lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

        xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
        eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
        eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

        ns = scalar_spectral_index((/eps1,eps2/))
        r = tensor_to_scalar_ratio((/eps1/))


        if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r


        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta,p/))

     enddo

  enddo


  call aspicwrite_end()



!$omp section


  call aspicwrite_header('rcipi2tune4',labeps12,labnsr,labbfoldreh,(/'alpha','beta ','p    '/))
 
  p=4._kp
  nalp = 20
  npts=15
  
  lnRhoRehMin = lnRhoNuc
  
  beta = 8._kp


  alphazero = rcipi_alpha_zero(p,beta)


  alphamin = -alphazero*0.99
  alphamax = alphazero*0.99

  !     alphamin = alphazero + epsilon(1._kp)
  !     alphamax = 2._kp*sqrt(beta) - epsilon(1._kp)


  do i=1,nalp
     alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)

     alphanminmax = (alpha-alphamin)/(alphamax-alphamin)

     efoldmax = rcipi_efoldmax(p,alpha,beta)


     if (efoldmax.lt.efoldNum) then
        write(*,*)'efoldmax too small!'
        print *,'efoldmax =     ',efoldmax
        write(*,*)
        cycle
     endif


     xEnd = rcipi_x_endinf(p,alpha,beta)
     lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

     do k=1,npts
        lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

        xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
        eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
        eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

        ns = scalar_spectral_index((/eps1,eps2/))
        r = tensor_to_scalar_ratio((/eps1/))


        if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r


        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta,p/))

     enddo

  enddo


  call aspicwrite_end()

  
  

  
!$omp section
  
  call aspicwrite_header('rcipi22',labeps12,labnsr,labbfoldreh,(/'alphanminmax','beta        ','p           '/))

 
  p=2._kp
  nbet = 3
  nalp = 500
  npts=15
  
  lnRhoRehMin = lnRhoNuc

  betamin = 0.05_kp
  betamax = 0.15_kp
  
  do j=1,nbet
     beta = betamin + (betamax-betamin)*(real(j-1,kp))/real(nbet-1,kp)

     alphazero = rcipi_alpha_zero(p,beta)
     

     alphamin = -alphazero*0.9999
     alphamax = alphazero*0.9999

!     alphamin = alphazero + epsilon(1._kp)
!     alphamax = 2._kp*sqrt(beta) - epsilon(1._kp)

     
     do i=1,nalp
        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)

        alphanminmax = (alpha-alphamin)/(alphamax-alphamin)
        
        efoldmax = rcipi_efoldmax(p,alpha,beta)
        

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif


        xEnd = rcipi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1/))

           
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           
           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alphanminmax,beta,p/))
           
        enddo

     enddo
  enddo

  call aspicwrite_end()

!$omp section
  
  call aspicwrite_header('rcipi23',labeps12,labnsr,labbfoldreh,(/'alphanminmax','beta        ','p           '/))
  
  p=3._kp
  nbet = 3
  nalp = 1000
  npts=10
  
  lnRhoRehMin = lnRhoNuc

  betamin = 0.05_kp
  betamax = 0.15_kp
  
  do j=1,nbet
     beta = betamin + (betamax-betamin)*(real(j-1,kp))/real(nbet-1,kp)

     alphazero = rcipi_alpha_zero(p,beta)
     

     alphamin = -alphazero*0.9999
     alphamax = alphazero*0.9999

!     alphamin = alphazero + epsilon(1._kp)
!     alphamax = 2._kp*sqrt(beta) - epsilon(1._kp)

     
     do i=1,nalp
        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)

        alphanminmax = (alpha-alphamin)/(alphamax-alphamin)
        
        efoldmax = rcipi_efoldmax(p,alpha,beta)
        

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif


        xEnd = rcipi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1/))

           
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           
           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alphanminmax,beta,p/))
           
        enddo

     enddo
  enddo

  call aspicwrite_end()

!$omp section  
  
  call aspicwrite_header('rcipi24',labeps12,labnsr,labbfoldreh,(/'alphanminmax','beta        ','p           '/))

  p=4._kp
  nbet = 3
  nalp = 5000
  npts=10
  
  lnRhoRehMin = lnRhoNuc

  betamin = 0.05_kp
  betamax = 0.15_kp
  
  do j=1,nbet
     beta = betamin + (betamax-betamin)*(real(j-1,kp))/real(nbet-1,kp)

     alphazero = rcipi_alpha_zero(p,beta)
     

     alphamin = -alphazero*0.9999
     alphamax = alphazero*0.9999

!     alphamin = alphazero + epsilon(1._kp)
!     alphamax = 2._kp*sqrt(beta) - epsilon(1._kp)

     
     do i=1,nalp
        alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(nalp-1,kp)

        alphanminmax = (alpha-alphamin)/(alphamax-alphamin)
        
        efoldmax = rcipi_efoldmax(p,alpha,beta)
        

        if (efoldmax.lt.efoldNum) then
           write(*,*)'efoldmax too small!'
           print *,'efoldmax =     ',efoldmax
           write(*,*)
           cycle
        endif


        xEnd = rcipi_x_endinf(p,alpha,beta)
        lnRhoRehMax = rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar)

        do k=1,npts
           lnRhoReh = lnRhoRehMin + (k-1)*(lnRhoRehMax-lnRhoRehMin)/real(npts-1,kp)

           xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
           eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)
           eps2 = rcipi_epsilon_two(xstar,p,alpha,beta)
           eps3 = rcipi_epsilon_three(xstar,p,alpha,beta)                     

           ns = scalar_spectral_index((/eps1,eps2/))
           r = tensor_to_scalar_ratio((/eps1/))

           
           if (display) print *,'lnRhoReh= N*= ns= r= ',lnRhoReh,abs(bfoldstar),ns,r

           
           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alphanminmax,beta,p/))
           
        enddo

     enddo
  enddo
  call aspicwrite_end()
  
!$omp end parallel sections
   

  endif
  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10

  p=2._kp
  alpha = -0.1_kp
  beta = 0.5_kp

  npts = 20
  
  xEnd = rcipi_x_endinf(p,alpha,beta)

  if (rcipi_efoldmax(p,alpha,beta).lt.efoldNum) then
     stop 'efoldmax too small!'
  endif
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rcipi_x_rrad(p,alpha,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rcipi_epsilon_one(xstar,p,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rcipi_epsilon_one(xend,p,alpha,beta)
     VendOverVstar = rcipi_norm_potential(xend,p,alpha,beta)/rcipi_norm_potential(xstar,p,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rcipi_x_rreh(p,alpha,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

  
end program rcipimain
