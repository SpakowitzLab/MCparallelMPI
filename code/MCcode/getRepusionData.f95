! *---------------------------------------------------------------*
Subroutine getRepusionData(mc)
use simMod
implicit none
TYPE(MCvar), intent(inout) :: mc
integer I

mc%nDiameters = 31
mc%nPhiValues = 200
mc%phiMax=1.0
mc%dmax=0.3
mc%phiMin=0.0
mc%dMin=0.01


ALLOCATE(mc%ERepusionData(mc%nPhiValues,mc%nDiameters))

OPEN (UNIT = 1, FILE = "input/mu", STATUS = 'OLD')      
DO I=1,mc%nPhiValues
    read(1,*)  mc%ERepusionData(I,:)
enddo 
CLOSE(1)

RETURN     
end subroutine      
!---------------------------------------------------------------*
