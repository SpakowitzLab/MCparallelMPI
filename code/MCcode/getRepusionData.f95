! *---------------------------------------------------------------*
Subroutine getRepusionData(mc)
use simMod
implicit none
TYPE(MCvar), intent(inout) :: mc
integer I

mc%nDiameters = 3
mc%nPhiValues = 3
mc%phiMax=0.3
mc%dmax=0.25
mc%phiMin=0.0
mc%dMin=0.01


ALLOCATE(mc%ERepusionData(mc%nPhiValues,mc%nDiameters))

OPEN (UNIT = 1, FILE = "../phiData", STATUS = 'OLD')      
DO I=1,mc%nPhiValues
    read(1,*)  mc%ERepusionData(I,:)
enddo 
CLOSE(1)

RETURN     
end subroutine      
!---------------------------------------------------------------*
