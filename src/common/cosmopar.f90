!Some cosmological parameters used for inflation. This module sets
!their fiducial values because observable predictions have only a weak
!dependence in this (may change with increased data accuracy).
module cosmopar
  use infprec, only : kp
  implicit none
  
  public

  real(kp), parameter :: pi=acos(-1._kp)

  real(kp), parameter :: HubbleSquareRootOf3OmegaRad = 7.4585d-63
  real(kp), parameter :: HubbleSquareRootOf2OmegaRad = sqrt(2._kp/3._kp)*HubbleSquareRootOf3OmegaRad

!q: number of relativistic entropic dof
!g: number of relativistic energy dof
!RelatDofRatio is: qo^(4/3)/go x greh/qreh^4/3
!If q=g this is  : (go/greh)^1/3
  real(kp), parameter :: RelatDofRatio = 1._kp


  real(kp), parameter :: lnMpcToKappa = 130.282_kp

  real(kp), parameter :: lnMpinGeV=42.334_kp

!1MeV
!!  real(kp), parameter :: lnRhoNuc = -196.97

!10MeV
!!  real(kp), parameter :: lnRhoNuc = -187.747

!Tnuc=(30/pi^2)^(1/4)*0.01 GeV
  real(kp), parameter :: lnRhoNuc = log((10._kp)**(-8)/(1.2209*10._kp**19 &
             /sqrt(8._kp*acos(-1._kp)))**4)

!100MeV
!!  real(kp), parameter :: lnRhoNuc = -178.55

!only used for reheating using slow-roll (libslowroll)
!COBE quadrupole moment
  real(kp), parameter :: QrmsOverT = 6e-6

!pivot scale in Mpc^-1
  real(kp), parameter :: kstar = 0.05_kp !Mpc^-1

!Best scalar amp for slow-roll (update with new constraints)
  real(kp), parameter :: powerAmpScalar = 2.165e-9

end module cosmopar
