Subroutine adaptCof(downSuccess,nPTReplicas,cof,N_average,&
                     lowerRepExe,upperRepExe,lowerCofRail,upperCofRail)
use setPrecision
implicit none

! inputs
integer nPTReplicas
integer downSuccess(nPTReplicas)
integer N_average
double precision lowerRepExe 
double precision upperRepExe 
double precision lowerCofRail
double precision upperCofRail

! input/output
double precision cof(nPTReplicas)

! internal variables
double precision, allocatable :: newCof(:)
integer rep
double precision successRate

allocate( newCof(1:nPTReplicas) )

newCof(1)=Cof(1)
do rep=2,nPTReplicas
    successRate=dble(downSuccess(rep))/dble(N_average) 
    if (Cof(rep).lt.Cof(rep-1)) then
        print*, "Error in adaptCof!"
        stop 1
    endif
    ! addapt spacing
    if (successRate.lt.lowerRepExe) then
        newCof(rep)=newCof(rep-1)+(Cof(rep)-Cof(rep-1))*0.95_dp
    elseif (successRate.gt.upperRepExe) then
        newCof(rep)=newCof(rep-1)+(Cof(rep)-Cof(rep-1))*1.05_dp
    else
        newCof(rep)=newCof(rep-1)+Cof(rep)-Cof(rep-1)
    endif
    ! inforce limits on addaption
    if ((newCof(rep)-newCof(rep-1)).lt.lowerCofRail) then
        newCof(rep)=newCof(rep-1)+lowerCofRail
    endif
    if ((newCof(rep)-newCof(rep-1)).gt.upperCofRail) then
        newCof(rep)=newCof(rep-1)+upperCofRail
    endif
    if (newCof(rep).lt.newCof(rep-1)) then
        print*, "Error in adaptCof!!!"
        stop 1
    endif
enddo
do rep=1,nPTReplicas
    if (newcof(rep).lt.0.00001) then
        print*, "Error in adaptCof! how did that happen?"
        print*, newCof
        print*, Cof
        stop 1
    endif
    cof(rep)=newcof(rep)
enddo
deallocate (newCof)
end subroutine
