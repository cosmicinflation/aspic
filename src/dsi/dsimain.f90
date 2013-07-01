!test the reheating derivation from slow-roll
program dsimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use dsisr, only : dsi_epsilon_one, dsi_epsilon_two, dsi_epsilon_three,dsi_xinimin,dsi_xendmin,dsi_xendmax,dsi_mumax
  use dsireheat, only : dsi_lnrhoreh_max, dsi_x_star
  use infinout, only : delete_file, livewrite, has_not_shifted
  use srreheat, only : log_energy_reheat_ingev

  use dsisr, only : dsi_norm_potential
  use dsireheat, only : dsi_x_rreh, dsi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,l
  integer :: npts,nmu,nxEnd,np

  real(kp) :: p,mu,xEnd,w,bfoldstar,q
  real(kp) :: mumin,mumax,xEndmin,xEndmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp), dimension(:), allocatable ::pvalues,qvalues,muminvalues,mumaxvalues
  integer, dimension(:), allocatable ::nmuvalues, nxEndvalues

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End

  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('dsi_predic.dat')
  call delete_file('dsi_nsr.dat')

  npts = 20

  w=0._kp
  !  w = 1._kp/3._kp

  np=3
  allocate(pvalues(1:np))
  pvalues(1)=2._kp
  pvalues(2)=3._kp
  pvalues(3)=4._kp

  allocate(qvalues(1:np))
  qvalues(1)=8._kp
  qvalues(2)=8._kp
  qvalues(3)=8._kp

  allocate(nmuvalues(1:np))
  nmuvalues(1)=15
  nmuvalues(2)=15
  nmuvalues(3)=15

  allocate(mumaxvalues(1:np))
  mumaxvalues(1)=dsi_mumax(60._kp,pvalues(1),qvalues(1))
  mumaxvalues(2)=dsi_mumax(60._kp,pvalues(2),qvalues(2))
  mumaxvalues(3)=dsi_mumax(60._kp,pvalues(3),qvalues(3))

  allocate(muminvalues(1:np))
  muminvalues(1)=mumaxvalues(1)*10.**(-5.)
  muminvalues(2)=mumaxvalues(2)*10.**(-6.)
  muminvalues(3)=mumaxvalues(3)*10.**(-7.)

  allocate(nxEndvalues(1:np))
  nxEndvalues(1)=50
  nxEndvalues(2)=50
  nxEndvalues(3)=50

  do l=1,np
     p=pvalues(l)
     q=qvalues(l)
     nmu=nmuvalues(l)
     nxEnd=nxEndvalues(l)
     mumin=muminvalues(l)
     mumax=mumaxvalues(l)


     do j=0,nmu
        mu=mumin*(mumax/mumin)**(real(j,kp)/real(nmu,kp))

        xEndmin=dsi_xendmin(70._kp,p,mu)
        xEndmin=dsi_xendmin(60._kp,p,mu)
        xEndmax=dsi_xendmax(p,mu,q)
        print*,'xEndmin=',xEndmin,'xEndmax=',xEndmax


        do k=1,nxEnd
           xEnd=xEndmin*(xEndmax/xEndmin)**(real(k,kp)/real(nxEnd,kp)) 


           lnRhoRehMin = lnRhoNuc
           lnRhoRehMax = dsi_lnrhoreh_max(p,mu,xEnd,Pstar)

           print *,'p=',p,'mu=',mu,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

           do i=1,npts

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = dsi_x_star(p,mu,xEnd,w,lnRhoReh,Pstar,bfoldstar)


              eps1 = dsi_epsilon_one(xstar,p,mu)
              eps2 = dsi_epsilon_two(xstar,p,mu)
              eps3 = dsi_epsilon_three(xstar,p,mu)


              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

              logErehGeV = log_energy_reheat_ingev(lnRhoReh)
              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
                 cycle
              endif

              call livewrite('dsi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('dsi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           end do

        end do

     end do

  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  p= 2
  mu= 0.1
  xend = 100
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = dsi_x_rrad(p,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = dsi_epsilon_one(xstar,p,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  dsi_epsilon_one(xend,p,mu)
     VendOverVstar = dsi_norm_potential(xend,p,mu)/dsi_norm_potential(xstar,p,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = dsi_x_rreh(p,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = dsi_x_star(p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program dsimain
