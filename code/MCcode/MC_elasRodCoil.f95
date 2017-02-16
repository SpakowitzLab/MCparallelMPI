!---------------------------------------------------------------*

! Calculate the change in the bending energy for a displacement.
     
SUBROUTINE MC_eelasRodCoil(DEELAS,R,U,RP,UP,&
                    NT,NB,IB1,IB2,&
                    IT1,IT2,elasParam0,elasParam1,AB)
use elasMod
use setPrecision
IMPLICIT NONE
TYPE(elasParamType), intent(in) :: elasParam0
TYPE(elasParamType), intent(in) :: elasParam1
DOUBLE PRECISION, intent(in) :: R(NT,3)  ! Bead positions
INTEGER, intent(in) :: AB(NT)  ! Bead positions
DOUBLE PRECISION, intent(in) :: U(NT,3)  ! Tangent vectors
DOUBLE PRECISION, intent(in) :: RP(NT,3)  ! Bead positions
DOUBLE PRECISION, intent(in) :: UP(NT,3)  ! Tangent vectors
INTEGER, intent(in) :: NB                ! Number of beads in a polymer
INTEGER, intent(in) :: NT                ! Total number of beads
INTEGER, intent(in) :: IB1               ! Test bead position 1
INTEGER, intent(in) :: IT1               ! Index of test bead 1
INTEGER, intent(in) :: IB2               ! Test bead position 2
INTEGER, intent(in) :: IT2               ! Index of test bead 2

DOUBLE PRECISION, intent(out) :: DEELAS(3)   ! Change in ECOM      

!     Polymer properties

double precision EB
double precision EPAR
double precision EPERP
double precision GAM
double precision ETA

!     Variables for force and torque calculations

DOUBLE PRECISION DR(3),DRPAR,DRPERP(3)
DOUBLE PRECISION GI(3)

! Setup parameters

DEELAS(1)=0.0_dp
DEELAS(2)=0.0_dp
DEELAS(3)=0.0_dp

!     Calculate the change in the energy

if (AB(IT1-1).eq.0 .and. AB(IT1).eq.0) then
    EB=elasParam0%EB
    EPAR=elasParam0%EPAR
    EPERP=elasParam0%EPERP
    GAM=elasParam0%GAM
    ETA=elasParam0%ETA
else if (AB(IT1-1).eq.1 .and. AB(IT1).eq.1) then
    EB=   elasParam1%EB
    EPAR= elasParam1%EPAR
    EPERP=elasParam1%EPERP
    GAM=  elasParam1%GAM
    ETA=  elasParam1%ETA
else
    EB=   (elasParam1%EB   + elasParam0%EB   )/2 
    EPAR= (elasParam1%EPAR + elasParam0%EPAR )/2
    EPERP=(elasParam1%EPERP+ elasParam0%EPERP)/2
    GAM=  (elasParam1%GAM  + elasParam0%GAM  )/2
    ETA=  (elasParam1%ETA  + elasParam0%ETA  )/2
endif
if (IB1.NE.1) then
   
   DR(1)=R(IT1,1)-R(IT1-1,1)
   DR(2)=R(IT1,2)-R(IT1-1,2)
   DR(3)=R(IT1,3)-R(IT1-1,3)
   DRPAR=DR(1)*U(IT1-1,1)+DR(2)*U(IT1-1,2)+DR(3)*U(IT1-1,3)
   
   DRPERP(1)=DR(1)-DRPAR*U(IT1-1,1)
   DRPERP(2)=DR(2)-DRPAR*U(IT1-1,2)
   DRPERP(3)=DR(3)-DRPAR*U(IT1-1,3)
   !U1U2=U(IT1-1,1)*U(IT1,1)+U(IT1-1,2)*U(IT1,2)+U(IT1-1,3)*U(IT1,3)

   GI(1)=(U(IT1,1)-U(IT1-1,1)-ETA*DRPERP(1))
   GI(2)=(U(IT1,2)-U(IT1-1,2)-ETA*DRPERP(2))
   GI(3)=(U(IT1,3)-U(IT1-1,3)-ETA*DRPERP(3))

   DEELAS(1)=DEELAS(1)-0.5_dp*EB*(GI(1)**2+GI(2)**2+GI(3)**2) 
   DEELAS(2)=DEELAS(2)-0.5_dp*EPAR*(DRPAR-GAM)**2
   DEELAS(3)=DEELAS(3)-0.5_dp*EPERP*(DRPERP(1)**2+DRPERP(2)**2.+DRPERP(3)**2)

   DR(1)=RP(IT1,1)-R(IT1-1,1)
   DR(2)=RP(IT1,2)-R(IT1-1,2)
   DR(3)=RP(IT1,3)-R(IT1-1,3)
   DRPAR=DR(1)*U(IT1-1,1)+DR(2)*U(IT1-1,2)+DR(3)*U(IT1-1,3)
                  
   DRPERP(1)=DR(1)-DRPAR*U(IT1-1,1)
   DRPERP(2)=DR(2)-DRPAR*U(IT1-1,2)
   DRPERP(3)=DR(3)-DRPAR*U(IT1-1,3)
   !U1U2=U(IT1-1,1)*UP(IT1,1)+U(IT1-1,2)*UP(IT1,2)+U(IT1-1,3)*UP(IT1,3)

   GI(1)=(UP(IT1,1)-U(IT1-1,1)-ETA*DRPERP(1))
   GI(2)=(UP(IT1,2)-U(IT1-1,2)-ETA*DRPERP(2))
   GI(3)=(UP(IT1,3)-U(IT1-1,3)-ETA*DRPERP(3))

   DEELAS(1)=DEELAS(1)+0.5_dp*EB*(GI(1)**2+GI(2)**2+GI(3)**2)
   DEELAS(2)=DEELAS(2)+0.5_dp*EPAR*(DRPAR-GAM)**2
   DEELAS(3)=DEELaS(3)+0.5_dp*EPERP*(DRPERP(1)**2+DRPERP(2)**2+DRPERP(3)**2)
   
endif

if (AB(IT2).eq.0 .and. AB(IT2+1).eq.0) then
    EB=elasParam0%EB
    EPAR=elasParam0%EPAR
    EPERP=elasParam0%EPERP
    GAM=elasParam0%GAM
    ETA=elasParam0%ETA
else if (AB(IT2).eq.1 .and. AB(IT2+1).eq.1) then
    EB=   elasParam1%EB
    EPAR= elasParam1%EPAR
    EPERP=elasParam1%EPERP
    GAM=  elasParam1%GAM
    ETA=  elasParam1%ETA
else
    EB=   (elasParam1%EB   + elasParam0%EB   )/2 
    EPAR= (elasParam1%EPAR + elasParam0%EPAR )/2
    EPERP=(elasParam1%EPERP+ elasParam0%EPERP)/2
    GAM=  (elasParam1%GAM  + elasParam0%GAM  )/2
    ETA=  (elasParam1%ETA  + elasParam0%ETA  )/2
endif
if (IB2.NE.NB) then
   
   DR(1)=R(IT2+1,1)-R(IT2,1)
   DR(2)=R(IT2+1,2)-R(IT2,2)
   DR(3)=R(IT2+1,3)-R(IT2,3)
   DRPAR=DR(1)*U(IT2,1)+DR(2)*U(IT2,2)+DR(3)*U(IT2,3)
                  
   DRPERP(1)=DR(1)-DRPAR*U(IT2,1)
   DRPERP(2)=DR(2)-DRPAR*U(IT2,2)
   DRPERP(3)=DR(3)-DRPAR*U(IT2,3)
   !U1U2=U(IT2,1)*U(IT2+1,1)+U(IT2,2)*U(IT2+1,2)+U(IT2,3)*U(IT2+1,3)

   GI(1)=(U(IT2+1,1)-U(IT2,1)-ETA*DRPERP(1))
   GI(2)=(U(IT2+1,2)-U(IT2,2)-ETA*DRPERP(2))
   GI(3)=(U(IT2+1,3)-U(IT2,3)-ETA*DRPERP(3))
   
   DEELAS(1)=DEELAS(1)-0.5_dp*EB*(GI(1)**2.+GI(2)**2.+GI(3)**2.)
   DEELAS(2)=DEELAS(2)-0.5_dp*EPAR*(DRPAR-GAM)**2.
   DEELAS(3)=DEELAS(3)-0.5_dp*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)

   DR(1)=R(IT2+1,1)-RP(IT2,1)
   DR(2)=R(IT2+1,2)-RP(IT2,2)
   DR(3)=R(IT2+1,3)-RP(IT2,3)
   DRPAR=DR(1)*UP(IT2,1)+DR(2)*UP(IT2,2)+DR(3)*UP(IT2,3)
                  
   DRPERP(1)=DR(1)-DRPAR*UP(IT2,1)
   DRPERP(2)=DR(2)-DRPAR*UP(IT2,2)
   DRPERP(3)=DR(3)-DRPAR*UP(IT2,3)
   !U1U2=UP(IT2,1)*U(IT2+1,1)+UP(IT2,2)*U(IT2+1,2)+UP(IT2,3)*U(IT2+1,3)

   GI(1)=(U(IT2+1,1)-UP(IT2,1)-ETA*DRPERP(1))
   GI(2)=(U(IT2+1,2)-UP(IT2,2)-ETA*DRPERP(2))
   GI(3)=(U(IT2+1,3)-UP(IT2,3)-ETA*DRPERP(3))
   
   DEELAS(1)=DEELAS(1)+0.5_dp*EB*(GI(1)**2.+GI(2)**2+GI(3)**2)
   DEELAS(2)=DEELAS(2)+0.5_dp*EPAR*(DRPAR-GAM)**2.
   DEELAS(3)=DEELaS(3)+0.5_dp*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)

endif         

RETURN      
END

!---------------------------------------------------------------*
