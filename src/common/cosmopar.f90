!Some cosmological parameters used for inflation. This module sets
!their fiducial values because observable predictions have only an
!ultra weak dependence in this (may change with increased data
!accuracy like 21cm).
module cosmopar
  use infprec, only : kp,pi
  implicit none
  
  public

!does not depend on h but on Tcmb
  real(kp), parameter :: HubbleSquareRootOf3OmegaRad = 7.5437d-63
  real(kp), parameter :: HubbleSquareRootOf2OmegaRad = sqrt(2._kp/3._kp)*HubbleSquareRootOf3OmegaRad

!hard prior: all models are assumed to support, at least, this number (120)
!of efolds which allows for all possible reheating history
  real(kp), parameter :: efoldNum = 80._kp

!hard prior: all models inflating in a region with eps>epsMax are not counted
  real(kp), parameter :: epsilonMax = 0.2_kp

!q: number of relativistic entropic dof
!g: number of relativistic energy dof
!RelatDofRatio is: qo^(4/3)/go x greh/qreh^4/3
!If q=g this is  : (go/greh)^1/3
  real(kp), parameter :: RelatDofRatio = 1._kp


  real(kp), parameter :: lnMpcToKappa = 130.282_kp

  real(kp), parameter :: lnMpinGeV = 42.334_kp

!1MeV
!!  real(kp), parameter :: lnRhoNuc = -196.98

!10MeV
  real(kp), parameter :: lnRhoNuc = -187.77

!100MeV
!!  real(kp), parameter :: lnRhoNuc = -178.56


!pivot scale in Mpc^-1
  real(kp), parameter :: kstar = 0.05_kp !Mpc^-1

!scalar amp at kstar from planck 2013 + slow-roll second order
  real(kp), parameter :: powerAmpScalar = 2.2030e-09

!Effective COBE normalisation given the same amplitude as Pstar
  real(kp), parameter :: QrmsOverT = sqrt(powerAmpScalar/60._kp)

end module cosmopar
