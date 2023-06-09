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

!scalar amp at kstar from planck 2018 + BKP + slow-roll second order
  real(kp), parameter :: powerAmpScalar = 2.1127e-09

!Effective COBE normalisation given the same amplitude as Pstar
  real(kp), parameter :: QrmsOverT = sqrt(powerAmpScalar/60._kp)

!Higgs vacuum expectation value in reduced Planck Mass. The Higgs doublet Lagrangian is
!
!     L(S) = -(DS^+)(DS) - V(S)
!
!with
!
!     V(S)  = -mu^2 S^+S + lambda (S^+S)^2
!
!and
!     S = 1/sqrt(2) (H + v) [0,1]
!
!Higgs inflation scalar (canonically normalized) is h = H + v from which one gets
!
!     V(h) = lambda/4 (h^2 -v^2)^2 + Cte
!
!     v = mu/sqrt(lambda) = 246GeV
!     mH = sqrt(2 lambda) v = 125GeV
!
!
!v  
  real(kp), parameter :: HiggsVeV = 246._kp/exp(lnMpinGeV)
!mH
  real(kp), parameter :: HiggsMass = 125.1_kp/exp(lnMpinGeV)
!lambda = mH^2/(2 v^2)
  real(kp), parameter :: HiggsCoupling = 0.5_kp*(HiggsMass/HiggsVeV)**2
  
  
end module cosmopar
