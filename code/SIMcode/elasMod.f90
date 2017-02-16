!   --------------------------------------------------------------
!
!    This module handles the SSWLC variables
!        By Quinn MacPherson ~ Winter 2017
!
!   --------------------------------------------------------------
Module elasMod
use setPrecision
IMPLICIT NONE

Type elasParamType   ! Structure for simulation variables of known size
    DOUBLE PRECISION PARA(10) ! Parameters for sswlc
        ! EB, EPAR, EPERP, GAM, ETA, ...

    DOUBLE PRECISION EB,EPAR,EPERP
    DOUBLE PRECISION GAM,ETA
    DOUBLE PRECISION L0       ! Equilibrium segment length
    DOUBLE PRECISION EPS       ! L0/2lp

end TYPE
contains
SUBROUTINE printElasParam(elasParam)
IMPLICIT NONE
TYPE(elasParamType), intent(in) :: elasParam

print*, "L0=",elasParam%L0
print*, "EPS=",elasParam%EPS
print*, "EB=",elasParam%EB
print*, "EPAR=",elasParam%EPAR
print*, "EPERP=",elasParam%EPERP
print*, "GAM=",elasParam%GAM

end subroutine
SUBROUTINE getpara(elasParam,EPS,L0)

IMPLICIT NONE
TYPE(elasParamType), intent(out) :: elasParam
DOUBLE PRECISION, intent(in) :: EPS   ! Elasticity l0/(2lp)
DOUBLE PRECISION, intent(in) :: L0       ! Equilibrium segment length
DOUBLE PRECISION PVEC(679,8)
INTEGER IND,CRS
DOUBLE PRECISION EB,EPAR,EPERP
DOUBLE PRECISION GAM,ETA
DOUBLE PRECISION M
DOUBLE PRECISION DEL
INTEGER I

!     Load in the parameters for the simulation

DEL=2.*EPS



!     Load the tabulated parameters

OPEN (UNIT=5,FILE='input/dssWLCparams',STATUS='OLD')
DO I=1,679
   READ(5,*) PVEC(I,1),PVEC(I,2),PVEC(I,3),PVEC(I,4),PVEC(I,5), &
             PVEC(I,6),PVEC(I,7),PVEC(I,8)
enddo 
CLOSE(5)

if (DEL.LT.PVEC(1,1)) then
   DEL=PVEC(1,1)
endif
if (DEL.GT.PVEC(679,1)) then
   DEL=PVEC(679,1)
endif

CRS=0
IND=1
do while (CRS.EQ.0)
   if (DEL.LE.PVEC(IND,1)) then
      CRS=1
   else
      IND=IND+1
   endif
enddo

I=2 
M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
EB=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

I=3 
M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
GAM=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

I=4
M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
EPAR=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

I=5
M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
EPERP=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

I=6
M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
ETA=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

!      I=7
!      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
!      XIU=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

!      I=8
!      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
!      DT=XIU*(M*(DEL-PVEC(IND,1))+PVEC(IND,I))

EB=EB/DEL
EPAR=EPAR*DEL/L0**2.
EPERP=EPERP*DEL/L0**2.
ETA=ETA*DEL/L0
GAM=L0*GAM
     
elasParam%EB=EB
elasParam%EPAR=EPAR
elasParam%EPERP=EPERP
elasParam%GAM=GAM
elasParam%ETA=ETA
elasParam%EPS=EPS
elasParam%L0=L0

elasParam%PARA(1)=EB
elasParam%PARA(2)=EPAR
elasParam%PARA(3)=EPERP
elasParam%PARA(4)=GAM
elasParam%PARA(5)=ETA
elasParam%PARA(6)=0
elasParam%PARA(7)=0
elasParam%PARA(8)=0
elasParam%PARA(9)=0
elasParam%PARA(10)=0

RETURN     
END
End MODULE      
!---------------------------------------------------------------*
