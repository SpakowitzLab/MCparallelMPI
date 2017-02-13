!-----------------------------------------------------------!
!
!     Make energy Huge when polymer goes outside boundary
!
!            Started by Quinn 2/17/15
!
!
!  
! confineType  |  Discription
! _____________|_________________________________
!    0         |  No confinement
!    1         |  Betwene two plates in Z direction at 0 and LBox
!    2         |  Cube of size LBox**3,  range: 0-LBox
!    3         |  Circle of radius LBox, centered at LBox/2
!    4         |  Periodic, non-equal lengths
!    5         |  Cylinder with spherical caps

SUBROUTINE MC_confine(confineType, LBox, RP, NT, IT1, IT2, ECon,RCylinder,LCylinder)
use setPrecision


IMPLICIT NONE

DOUBLE PRECISION, intent(in) :: RCylinder
DOUBLE PRECISION, intent(in) :: LCylinder
INTEGER, intent(in) :: confineType  ! Specifier for type of confinement
DOUBLE PRECISION, intent(in) :: LBox(3) ! Side length of box
INTEGER, intent(in) :: NT     ! Total number of beads in simulation
DOUBLE PRECISION, intent(in) :: RP(NT,3)  ! Bead positions
INTEGER, intent(in) :: IT1    ! Start test bead
INTEGER, intent(in) :: IT2    ! Final test bead
INTEGER I      ! Index of bead being compared
DOUBLE PRECISION, intent(out) :: ECon
Double PRECISION r(3)
ECon=0.0_dp


if (confineType.EQ.0) then
    return
elseif(confineType.EQ.1) then
    ! Confinement only in the z-direction
    ! limits: 0 and LBox
    DO I=IT1,IT2
        if(RP(I,3)<0.0_dp) then
            ECon=9990000.0_dp
        elseif (RP(I,3)>LBox(3)) then
            ECon=9990000.0_dp
        endif
    ENDDO
elseif(confineType.EQ.2) then
    DO I=IT1,IT2
        if(RP(I,1)<0.0) then
            ECon=9990000.0_dp
        elseif(RP(I,1)>LBox(1)) then
            ECon=9990000.0_dp
        elseif(RP(I,2)<0.0) then
            ECon=9990000.0_dp
        elseif(RP(I,2)>LBox(2)) then
            ECon=9990000.0
        elseif(RP(I,3)<0.0) then
            ECon=9990000.0_dp
        elseif(RP(I,3)>LBox(3)) then
            ECon=9990000.0_dp
        endif
    ENDDO
elseif(confineType.EQ.3) then
    DO I=IT1,IT2
        if(((RP(I,1)-LBox(1)/2)**2 + (RP(I,2)-LBox(1)/2_dp)**2 + &
           (RP(I,3)-LBox(1)/2)**2).GT.dble(LBox(1)*LBox(1)*0.25_dp)) then
            ECon=9990000.0
            return
        endif
    Enddo    
elseif(confineType.EQ.4) then
    return
elseif(confineType.EQ.5) then ! Cylinder with spherical caps
    DO I=IT1,IT2
        if( (RP(I,2)-RCylinder)**2 + (RP(I,3)-RCylinder)**2 .GT. Rcylinder**2) then
            ECon=9990000.0
            return
        endif
        if( RP(I,1) .LT. RCylinder) then
            r(1) = RP(I,1) - RCylinder
            r(2) = RP(I,2) - RCylinder
            r(3) = RP(I,3) - RCylinder
            if (r(1)**2+r(2)**2+r(3)**2 .gt. RCylinder**2) then
                ECon=9990000.0
                return
            endif
        elseif ( RP(I,1) .GT. RCylinder+LCylinder) then
            r(1) = RP(I,1) - RCylinder + LCylinder
            r(2) = RP(I,2) - RCylinder
            r(3) = RP(I,3) - RCylinder
            if (r(1)**2+r(2)**2+r(3)**2 .gt. RCylinder**2) then
                ECon=9990000.0
                return
            endif
        endif
            
    enddo 
else 
   print*, "Undefined comfone Type"
   stop 1
endif




END
