!-----------------------------------------------------------!
!
!         Calculate HP1 Binding Energy
!
!            Started by Quinn 12/16/15
!
!
!  sign convention: EM and EU are more positive for favorable binding
!  Typical Values: EU=-1.52 and EM=0.01  

SUBROUTINE MC_bind(NT,BPM,IT1,IT2,AB,ABP,METH,EU,EM,DEBind,DEmu,mu,dx_mu)
use setPrecision
IMPLICIT NONE
INTEGER, intent(in) :: NT     ! Total number of beads in simulation
INTEGER, intent(in) :: BPM    ! Number of beads per monomer
INTEGER, intent(in) :: IT1    ! Start test bead
INTEGER, intent(in) :: IT2    ! Final test bead
INTEGER, intent(in) :: AB(NT)   ! Chemical identity (a.k.a. binding state)
INTEGER, intent(in) :: ABP(NT)  ! Test Chemical identity
INTEGER, intent(in) :: METH(NT) ! Methalation state (unerlyin chamical type)
DOUBLE PRECISION, intent(in) :: EU        ! Binding energy of Unemethalted state
DOUBLE PRECISION, intent(in) :: EM        ! Binding energy of methalated state
DOUBLE PRECISION, intent(in) :: mu   ! Chemical potential of HP1
DOUBLE PRECISION, intent(out) :: DEBind    ! Change in binding energy
DOUBLE PRECISION, intent(out) :: DEMu    ! Change in chemcial potential energy
double precision, intent(out) :: dx_mu ! -n_bound
INTEGER I      ! Index of bead being compared
DEBind=0.0_dp
DEMu=0.0_dp
Dx_mu=0.0_dp
DO I=IT1,IT2,BPM
    if(METH(I).EQ.1) then
        DEBind=DEBind+(-EM)*real(ABP(I)-AB(I))
        DEMu=DEMu+(-mu)*real(ABP(I)-AB(I))
        Dx_mu=Dx_mu-real(ABP(I)-AB(I))
    else
        DEBind=DEBind+(-EU)*real(ABP(I)-AB(I))
        DEMu=DEMu+(-mu)*real(ABP(I)-AB(I))
        Dx_mu=Dx_mu-real(ABP(I)-AB(I))
    endif
ENDDO

END
