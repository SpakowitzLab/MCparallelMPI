Subroutine save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                           cof,x,nodeNumber,N_average,nExchange,IND)
    use setPrecision
! Print Energy data
    IMPLICIT NONE
    LOGICAL isfile
    character*32 fullName
    integer nPTReplicas
    integer rep
    integer upSuccess(nPTReplicas)
    integer downSuccess(nPTReplicas)
    integer nodeNumber(nPTReplicas)
    double precision cof(nPTReplicas)
    double precision x(nPTReplicas)
    integer N_average
    integer nExchange
    integer IND
    double precision dxdcof
    fullName=  'data/repHistory'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif

    write(1,*) "~~~~~~~~~~~exchange: ",nExchange,", IND:",IND,"~~~~~~~~~~~~~~~~~~~~"
    write(1,*) " rep |  cof  |   x    |  up  | down |node| dxdcof|"
    do rep=1,nPTReplicas
        if (rep.ne.nPTReplicas) then
            dxdcof=(x(rep)-x(rep+1))*(cof(rep+1)-cof(rep))
        else
            dxdcof=0.0_dp
        endif
        write(1,"(I6,f8.3,f9.1,2f7.3,I5,f8.2)"), rep, cof(rep), x(rep), & 
                 real(upSuccess(rep))/real(N_average), &
                 real(downSuccess(rep))/real(N_average), nodeNumber(rep),&
                 dxdcof
    enddo  
    Close(1)
    fullName=  'data/cofData'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif
    write(1,*), IND, cof, nodeNumber
    Close(1)
end subroutine
