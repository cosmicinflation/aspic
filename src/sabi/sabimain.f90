!test the reheating derivation from slow-roll
program sabimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sabisr, only : sabi_epsilon_one, sabi_epsilon_two,sabi_epsilon_three
  use sabireheat, only : sabi_lnrhoreh_max, sabi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use sabisr, only : sabi_norm_potential, sabi_x_endinf, sabi_x_trajectory
  use sabireheat, only : sabi_x_rreh, sabi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh


  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nn

  real(kp) :: alpha,n,w,bfoldstar,alphamin,alphamax,nmin,nmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(:), allocatable ::nvalues

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: xmin,xmax,V1,x

  Pstar = powerAmpScalar


  call delete_file('sabi_potential.dat')
  call delete_file('sabi_slowroll.dat')

  npts=100

  xmin = 1e-2_kp
  xmax = 5._kp
  
  do i=1,npts
     x = exp(log(xmin) + real(i-1,kp)*(log(xmax)-log(xmin))/real(npts-1,kp))

     V1 = sabi_norm_potential(x,alpha=0.1_kp,n=2._kp)


     call livewrite('sabi_potential.dat',x,V1)

     if (x.eq.0._kp) cycle
     
     eps1 = sabi_epsilon_one(x,alpha=0.1_kp,n=2._kp)
     eps2 = sabi_epsilon_two(x,alpha=0.1_kp,n=2._kp)
     eps3 = sabi_epsilon_three(x,alpha=0.1_kp,n=2._kp)
          
     call livewrite('sabi_slowroll.dat',x,eps1,eps2,eps3)
     
  enddo

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!
!!      n=1      !!
!!!!!!!!!!!!!!!!!!!

  call delete_file('sabi_predic_neq1.dat')
  call delete_file('sabi_nsr_neq1.dat')

  call aspicwrite_header('sabi',labeps12,labnsr,labbfoldreh,(/'alpha','n    '/))
  
  n=2._kp

  npts = 20
  nalpha = 500

  w=0._kp
  !  w = 1._kp/3._kp


  alphamin=0.05_kp
  alphamax=80._kp

  do k=0,nalpha
    alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp))

    print*,'alpha=',alpha,'n=',n

    lnRhoRehMin = lnRhoNuc
    xEnd = sabi_x_endinf(alpha,n)
    lnRhoRehMax = sabi_lnrhoreh_max(alpha,n,xend,Pstar)

    print *,lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


    do i=1,npts

      lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

      xstar = sabi_x_star(alpha,n,xend,w,lnRhoReh,Pstar,bfoldstar)

      eps1 = sabi_epsilon_one(xstar,alpha,n)
      eps2 = sabi_epsilon_two(xstar,alpha,n)
      eps3 = sabi_epsilon_three(xstar,alpha,n)

      print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1


      logErehGeV = log_energy_reheat_ingev(lnRhoReh)
      Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

      ns = 1._kp - 2._kp*eps1 - eps2
      r =16._kp*eps1

      print*, 'nS=',ns,'r=',r

      call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,n/))

      call livewrite('sabi_predic_neq1.dat',alpha,n,eps1,eps2,eps3,ns,r,Treh)

      call livewrite('sabi_nsr_neq1.dat',ns,r,abs(bfoldstar),lnRhoReh)


    end do


 end do


!!!!!!!!!!!!!!!!!!!
!!      n=5      !!
!!!!!!!!!!!!!!!!!!!

 call delete_file('sabi_predic_neq5.dat')
 call delete_file('sabi_nsr_neq5.dat')

  n=5._kp

  npts = 20
  nalpha = 500

  w=0._kp
  !  w = 1._kp/3._kp


  alphamin=0.05_kp
  alphamax=80._kp

  do k=0,nalpha
    alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp))

    print*,'alpha=',alpha,'n=',n


    lnRhoRehMin = lnRhoNuc
    xEnd = sabi_x_endinf(alpha,n)
    lnRhoRehMax = sabi_lnrhoreh_max(alpha,n,xend,Pstar)

    print *,lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


    do i=1,npts

      lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

      xstar = sabi_x_star(alpha,n,xend,w,lnRhoReh,Pstar,bfoldstar)

      eps1 = sabi_epsilon_one(xstar,alpha,n)
      eps2 = sabi_epsilon_two(xstar,alpha,n)
      eps3 = sabi_epsilon_three(xstar,alpha,n)

      print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1


      logErehGeV = log_energy_reheat_ingev(lnRhoReh)
      Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

      ns = 1._kp - 2._kp*eps1 - eps2
      r =16._kp*eps1

      print*, 'nS=',ns,'r=',r

      call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,n/))
      
      call livewrite('sabi_predic_neq5.dat',alpha,n,eps1,eps2,eps3,ns,r,Treh)

      call livewrite('sabi_nsr_neq5.dat',ns,r,abs(bfoldstar),lnRhoReh)


    end do


end do


  n=20._kp

  npts = 20
  nalpha = 500

  w=0._kp
  !  w = 1._kp/3._kp


  alphamin=0.05_kp
  alphamax=80._kp

  do k=0,nalpha
     alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp))

     print*,'alpha=',alpha,'n=',n


     lnRhoRehMin = lnRhoNuc
     xEnd = sabi_x_endinf(alpha,n)
     lnRhoRehMax = sabi_lnrhoreh_max(alpha,n,xend,Pstar)

     print *,lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = sabi_x_star(alpha,n,xend,w,lnRhoReh,Pstar,bfoldstar)

        eps1 = sabi_epsilon_one(xstar,alpha,n)
        eps2 = sabi_epsilon_two(xstar,alpha,n)
        eps3 = sabi_epsilon_three(xstar,alpha,n)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        print*, 'nS=',ns,'r=',r

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,n/))

        call livewrite('sabi_predic_neq5.dat',alpha,n,eps1,eps2,eps3,ns,r,Treh)

        call livewrite('sabi_nsr_neq5.dat',ns,r,abs(bfoldstar),lnRhoReh)


     end do


  end do


  call aspicwrite_end()
  


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 2._kp
  n = 2._kp
  xEnd = sabi_x_endinf(alpha,n)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sabi_x_rrad(alpha,n,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = sabi_epsilon_one(xstar,alpha,n)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  sabi_epsilon_one(xend,alpha,n)
     VendOverVstar = sabi_norm_potential(xend,alpha,n)/sabi_norm_potential(xstar,alpha,n)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = sabi_x_rreh(alpha,n,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = sabi_x_star(alpha,n,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo




end program sabimain
