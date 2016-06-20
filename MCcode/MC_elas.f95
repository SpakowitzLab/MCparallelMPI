!---------------------------------------------------------------*

! Calculate the change in the bending energy for a displacement.
     
SUBROUTINE MC_eelas(DEELAS,R,U,RP,UP,&
                    NT,NB,IP,IB1,IB2,&
                    IT1,IT2,EB,EPAR,EPERP,GAM,ETA)

use setPrecision
IMPLICIT NONE
DOUBLE PRECISION R(NT,3)  ! Bead positions
DOUBLE PRECISION U(NT,3)  ! Tangent vectors
DOUBLE PRECISION RP(NT,3)  ! Bead positions
DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
INTEGER NB                ! Number of beads in a polymer
INTEGER NT                ! Total number of beads

INTEGER IP  ! not in use  ! Test polymer 
INTEGER IB1               ! Test bead position 1
INTEGER IT1               ! Index of test bead 1
INTEGER IB2               ! Test bead position 2
INTEGER IT2               ! Index of test bead 2

DOUBLE PRECISION DEELAS(3)   ! Change in ECOM      

!     Polymer properties

DOUBLE PRECISION EB,EPAR,EPERP
DOUBLE PRECISION GAM,ETA

!     Variables for force and torque calculations

DOUBLE PRECISION DR(3),DRPAR,DRPERP(3),DRPERPM(3)
DOUBLE PRECISION FI(3),TI(3)
DOUBLE PRECISION U1U2,GI(3),DOTGU,HI(3)

! Setup parameters

DEELAS(1)=0.0_dp
DEELAS(2)=0.0_dp
DEELAS(3)=0.0_dp

!     Calculate the change in the energy

if (IB1.NE.1) then
   
   DR(1)=R(IT1,1)-R(IT1-1,1)
   DR(2)=R(IT1,2)-R(IT1-1,2)
   DR(3)=R(IT1,3)-R(IT1-1,3)
   DRPAR=DR(1)*U(IT1-1,1)+DR(2)*U(IT1-1,2)+DR(3)*U(IT1-1,3)
   
   DRPERP(1)=DR(1)-DRPAR*U(IT1-1,1)
   DRPERP(2)=DR(2)-DRPAR*U(IT1-1,2)
   DRPERP(3)=DR(3)-DRPAR*U(IT1-1,3)
   U1U2=U(IT1-1,1)*U(IT1,1)+U(IT1-1,2)*U(IT1,2)+U(IT1-1,3)*U(IT1,3)

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
   U1U2=U(IT1-1,1)*UP(IT1,1)+U(IT1-1,2)*UP(IT1,2)+U(IT1-1,3)*UP(IT1,3)

   GI(1)=(UP(IT1,1)-U(IT1-1,1)-ETA*DRPERP(1))
   GI(2)=(UP(IT1,2)-U(IT1-1,2)-ETA*DRPERP(2))
   GI(3)=(UP(IT1,3)-U(IT1-1,3)-ETA*DRPERP(3))

   DEELAS(1)=DEELAS(1)+0.5_dp*EB*(GI(1)**2+GI(2)**2+GI(3)**2)
   DEELAS(2)=DEELAS(2)+0.5_dp*EPAR*(DRPAR-GAM)**2
   DEELAS(3)=DEELaS(3)+0.5_dp*EPERP*(DRPERP(1)**2+DRPERP(2)**2+DRPERP(3)**2)
   
endif

if (IB2.NE.NB) then
   
   DR(1)=R(IT2+1,1)-R(IT2,1)
   DR(2)=R(IT2+1,2)-R(IT2,2)
   DR(3)=R(IT2+1,3)-R(IT2,3)
   DRPAR=DR(1)*U(IT2,1)+DR(2)*U(IT2,2)+DR(3)*U(IT2,3)
                  
   DRPERP(1)=DR(1)-DRPAR*U(IT2,1)
   DRPERP(2)=DR(2)-DRPAR*U(IT2,2)
   DRPERP(3)=DR(3)-DRPAR*U(IT2,3)
   U1U2=U(IT2,1)*U(IT2+1,1)+U(IT2,2)*U(IT2+1,2)+U(IT2,3)*U(IT2+1,3)

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
   U1U2=UP(IT2,1)*U(IT2+1,1)+UP(IT2,2)*U(IT2+1,2)+UP(IT2,3)*U(IT2+1,3)

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
