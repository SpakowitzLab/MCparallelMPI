! *---------------------------------------------------------------*
Subroutine getRepusionData(mc)
use simMod
implicit none
TYPE(MCvar), intent(inout) :: mc
integer I,J

mc%nDiameters = 31
mc%nPhiValues = 200
mc%phiMax=1.0
mc%dmax=0.3
mc%phiMin=0.0
mc%dMin=0.01


ALLOCATE(mc%exclusionMu(mc%nDiameters,mc%nPhiValues))

OPEN (UNIT = 1, FILE = "input/mu", STATUS = 'OLD')      
DO I=1,mc%nDiameters
    read(1,*)  mc%exclusionMu(I,:)
enddo 
CLOSE(1)

Do J=1,mc%nPhiValues
    Do I=1,mc%nDiameters 
        if (mc%exclusionMu(I,J)<0.0) then
            mc%exclusionMu(I,J)=0.0
        elseif (mc%exclusionMu(I,J)>100) then
            mc%exclusionMu(I,J)=100
        endif
        
    enddo
enddo

RETURN     
end subroutine      
!---------------------------------------------------------------*
