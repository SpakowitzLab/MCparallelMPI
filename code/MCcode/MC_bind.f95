!-----------------------------------------------------------!
!
!         Calculate HP1 Binding Energy
!
!            Started by Quinn 12/16/15
!
!
!  sign convention: EM and EU are more positive for favorable binding
!  Typical Values: EU=-1.52 and EM=0.01  

SUBROUTINE MC_bind(NT,BPM,IT1,IT2,AB,ABP,METH,EU,EM,DEBind,mu)
use setPrecision
IMPLICIT NONE
INTEGER NT     ! Total number of beads in simulation
INTEGER BPM    ! Number of beads per monomer
INTEGER IT1    ! Start test bead
INTEGER IT2    ! Final test bead
INTEGER AB(NT)   ! Chemical identity (a.k.a. binding state)
INTEGER ABP(NT)  ! Test Chemical identity
INTEGER METH(NT) ! Methalation state (unerlyin chamical type)
DOUBLE PRECISION EU        ! Binding energy of Unemethalted state
DOUBLE PRECISION EM        ! Binding energy of methalated state
DOUBLE PRECISION  DEBind    ! Change in binding energy
INTEGER I      ! Index of bead being compared
DOUBLE PRECISION mu   ! Chemical potential of HP1

DEBind=0.0_dp

DO I=IT1,IT2,BPM
    if(METH(I).EQ.1) then
        DEBind=DEBind+(-mu-EM)*real(ABP(I)-AB(I))
    else
        DEBind=DEBind+(-mu-EU)*real(ABP(I)-AB(I))
        !print*, 'In MC_bind EU:',EU,' EM:',EM
    endif
ENDDO

if (abs(DEBind).gt.100000) then
    
    DO I=IT1,IT2,BPM
        print*, ABP(I), AB(I)
    ENDDO
    print*, "range:",IT1,IT2
    print*, "error in MC_bind"
    print*, "DEBind",DEBind
    print*, "mu", mu
    print*, "EM",EM,"   EU",EU
    stop 1
endif      


END
