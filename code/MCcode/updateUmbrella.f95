!---------------------------------------------------------------!
!
!     Ajusts potential to insetivize moving to unexplored taritory.
!     Quinn started this subroutine on 1/11/17
!
! ---------------------------------------------------  
SUBROUTINE updateUmbrella(mc,md)
use simMod
use setPrecision
IMPLICIT NONE

TYPE(MCvar), intent(inout) :: mc   
TYPE(MCData), intent(inout) :: md

integer bin, total
double precision maxV
double precision minUmbrellaValue
character*16 fileName
double precision damping

md%umbrellaCounts=md%umbrellaCounts+1 ! incase of zero

total=sum(md%umbrellaCounts)
if (total.lt.1000) then
    if (total+mc%nOutside.lt.1000) then
        print*, "Something Probably when wrong with Umbrella Sampling"
        print*, "Total", total, "mc%nOutside", mc%nOutside
        stop 1
    else
        md%umbrellaCounts=md%umbrellaCounts-1 ! incase of zero
        print*, "Warning, mostly outside of umbrella range ..."
    endif
        
endif

damping=0.5
maxV=-100000.0
do bin=1,mc%nUmbrellaBins
     ! V = V + log(probability)
     md%umbrellaV(bin)=md%umbrellaV(bin)+&
         damping*log(real(md%umbrellaCounts(bin))/real(total)) 
     maxV=max(maxV,md%umbrellaV(bin))
enddo

md%umbrellaV=md%umbrellaV-maxV ! re-zero lowest
! Apply lower limit
minUmbrellaValue=-100.0  ! in kT
do bin=1,mc%nUmbrellaBins
    md%umbrellaV(bin)=max(minUmbrellaValue,md%umbrellaV(bin))
enddo

! Change current energy
call calcUmbrellaE(mc,md,mc%rxnQ,mc%EUmbrella,bin)

! save data for analisys
fileName='data/umbV'
call saveUmbrella(mc,md,fileName)

!rezero Counts
md%umbrellaCounts=0
mc%nOutside=0;
end subroutine
