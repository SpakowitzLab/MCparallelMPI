program main

use setPrecision
use simMod

IMPLICIT NONE
TYPE(MCvar) mc
TYPE(MCData) md
character*16 fileName ! file name to load from
character*16 paramsfile ! file name to load from
Integer I, IT1, IT2

fileName=trim('data/r1v1')
paramsfile=trim('input/params')
call MCvar_setParams(mc,paramsfile)
call MCvar_allocate(mc,md)
OPEN (UNIT = 5, FILE = fileName, STATUS = 'OLD')
Do I=1,mc%NT
   READ(5,*) md%R(I,1),md%R(I,2),md%R(I,3),md%AB(I)
enddo
CLOSE(5)
IT1=1
IT2=mc%NT
do I=1,mc%NBIN
     md%PHIA(I)=0.0_dp
     md%PHIB(I)=0.0_dp
enddo
do I=1,mc%NBIN
     md%Vol(I)=mc%del**3
enddo
call MC_int(mc,md,IT1,IT2,.True.)
mc%repSufix=''
fileName=trim('data/phi_from_r')
call MCVar_savePHI(mc,md,fileName)

end program
