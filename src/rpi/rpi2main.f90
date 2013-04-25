!test the reheating derivation from slow-roll
program rpi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rpicommon, only : rpi_x_potmax
  use rpi1sr, only : rpi1_x_endinf
  use rpi2sr, only : rpi2_epsilon_one, rpi2_epsilon_two, rpi2_epsilon_three
  use rpi2reheat, only : rpi2_lnrhoreh_max, rpi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  integer, parameter :: Np=1
  real(kp), parameter :: pmin=1._kp
  real(kp), parameter :: pmax=1.5_kp

  integer, parameter :: Nyend = 1
  real(kp) :: yendMin, yendMax

  real(kp) :: p,w,bfoldstar,yend
  real(kp) :: lnRhoReh,ystar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('rpi2_predic.dat')
  call delete_file('rpi2_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Np 
    p=pmin+(pmax-pmin)*(real(j,kp)/Np)
    w=(1._kp-p)/(3._kp-1._kp)

    

    lnRhoRehMin = lnRhoNuc
    
    if (p.eq.1._kp) then
       yendMin = rpi1_x_endinf(p)
       yendMax = 2*yendMin
    else
       yendMin = rpi_x_potmax(p) + 10
       yendMax = 2*yendMin
    endif
    

    do k=0,Nyend
       yend = yendMin + (yendMax-yendMin)*(real(k,kp)/Nyend)

       print *,'p= yend= ',p,yend

       lnRhoRehMax = rpi2_lnrhoreh_max(p,yend,Pstar)

       print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

       do i=1,npts

          lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)



          ystar = rpi2_x_star(p,yend,w,lnRhoReh,Pstar,bfoldstar)


          print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'ystar=',ystar


          eps1 = rpi2_epsilon_one(ystar,p)
          eps2 = rpi2_epsilon_two(ystar,p)
          eps3 = rpi2_epsilon_three(ystar,p)


          logErehGeV = log_energy_reheat_ingev(lnRhoReh)


          Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


          ns = 1._kp - 2._kp*eps1 - eps2
          r =16._kp*eps1

          call livewrite('rpi2_predic.dat',p,yend,eps1,eps2,eps3,r,ns,Treh)

          call livewrite('rpi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

       end do

    end do

 enddo

end program rpi2main
