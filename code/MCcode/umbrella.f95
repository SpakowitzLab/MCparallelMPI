!---------------------------------------------------------------!
!
!     The subroutine calculates the reaction coordinate, Qout,
!     for the r values input as rIn. You can use this either to
!     calculate for current or proposed positions.
!
!     If you make I1<0 and I2<0 than it will do calculation based
!     only on mc%R.
!
!     Quinn started this subroutine on 1/11/17
!
! ---------------------------------------------------  
SUBROUTINE calcRxnQ(mc,md,Qout,I1,I2)
use simMod
use setPrecision
IMPLICIT NONE

TYPE(MCvar), intent(inout) :: mc   
TYPE(MCData), intent(inout) :: md
Double Precision, intent(out) :: Qout
integer, intent(in) :: I1
integer, intent(in) :: I2
double precision cen1(3)
double precision cen2(3)
integer skip  ! to save time
integer ii,npts


skip=15
cen1(1:3)=(/0.0, 0.0, 0.0/)
cen2(1:3)=(/0.0, 0.0, 0.0/)

npts=0
do ii = 1,mc%nBlockP2,skip
     if (ii.lt.I1 .or. ii.gt.I2) then
         cen1(1)=cen1(1)+md%R_P2(ii,1)
         cen1(2)=cen1(2)+md%R_P2(ii,2)
         cen1(3)=cen1(3)+md%R_P2(ii,3)
     else
         cen1(1)=cen1(1)+md%Rp_P2(ii,1)
         cen1(2)=cen1(2)+md%Rp_P2(ii,2)
         cen1(3)=cen1(3)+md%Rp_P2(ii,3)
     endif
     npts=npts+1
enddo
cen1=cen1/npts

npts=0
do ii = 2*mc%nBlockP2+1,mc%nBeadsP2,skip
     if (ii.lt.I1 .or. ii.gt.I2) then
         cen2(1)=cen2(1)+md%R_P2(ii,1)
         cen2(2)=cen2(2)+md%R_P2(ii,2)
         cen2(3)=cen2(3)+md%R_P2(ii,3)
     else
         cen2(1)=cen2(1)+md%Rp_P2(ii,1)
         cen2(2)=cen2(2)+md%Rp_P2(ii,2)
         cen2(3)=cen2(3)+md%Rp_P2(ii,3)
     endif
     npts=npts+1
enddo
cen2=cen2/npts

Qout=sqrt(sum((cen1-cen2)**2))
end subroutine

!---------------------------------------------------------------!
!
!     Calculate (lookup) umbrella energy
!     Quinn started this subroutine on 1/11/17
!
! ---------------------------------------------------  
SUBROUTINE calcUmbrellaE(mc,md,rxnQ,Energy,bin)
use simMod
use setPrecision
IMPLICIT NONE

TYPE(MCvar), intent(inout) :: mc   
TYPE(MCData), intent(inout) :: md
Double Precision, intent(in) :: rxnQ
Double Precision, intent(out) :: Energy
integer, intent(out) :: bin

if (rxnQ.lt.mc%minQ) then
    Energy=0.0
    bin=-1
    return
elseif (rxnQ.gt.mc%maxQ) then
    Energy=0.0
    bin=-1
    return
endif
bin=ceiling(mc%nUmbrellaBins*(rxnQ-mc%minQ)/(mc%maxQ-mc%minQ))
Energy=md%umbrellaV(bin)

end subroutine



