Subroutine replicaExchange(mc,md)
! This checks in with the mpi head node to 
! For parallel tempering of the form:  E=cof*x
! 1: Tell head node the x value
! 2: Recive replica assignment from head node
! 3: Recive assigned cof value
    use setPrecision
    use mpi
    use simMod
    IMPLICIT NONE
    integer, parameter :: nTerms=8  ! number of energy terms 
    integer (kind=4) id, ierror
    TYPE(MCvar), intent(inout) :: mc
    Type(MCData), intent(inout) :: md     ! system allocated data
    integer i  ! working intiger
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    integer (kind=4) nThreads
    character*16 iostr ! for handling sufix string
    integer status(MPI_STATUS_SIZE)  ! MPI status
    double precision cof(nTerms)
    double precision cofOld(nTerms)
    double precision x(nTerms)
    double precision test(5)
    logical isfile

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    if (nThreads.lt.3) return

    x(1)=mc%x_Chi
    x(2)=mc%x_mu
    x(3)=mc%x_Field
    x(4)=mc%x_couple
    x(5)=mc%x_kap
    x(6)=0.0_dp !x(6)=mc%EElas(1)/mc%Para(1)
    x(7)=0.0_dp !x(7)=mc%EElas(2)/mc%Para(2)
    x(8)=0.0_dp !x(8)=mc%EElas(3)/mc%Para(3)

    test(1)=mc%EChi/mc%Chi
    test(3)=mc%EField/mc%h_A
    test(4)=mc%ECouple/mc%HP1_Bind
    test(5)=mc%EKap/mc%Kap

    cofOld(1)=mc%chi      
    cofOld(2)=mc%mu     
    cofOld(3)=mc%h_A     
    cofOld(4)=mc%HP1_Bind
    cofOld(5)=mc%Kap    
    cofOld(6)=mc%Para(1)
    cofOld(7)=mc%Para(2)
    cofOld(8)=mc%Para(3) 

    do i=1,5
        if (i.eq.2) cycle ! doesn't work for mu
        if (abs(cofOld(I)*(test(I)-x(I))).lt.0.0001) cycle
        if (abs(cofOld(I)).lt.0.00000001) cycle
        inquire(file = "data/error", exist=isfile)
        if (isfile) then
            OPEN (UNIT = 1, FILE = "data/error", STATUS ='OLD', POSITION="append")
        else 
            OPEN (UNIT = 1, FILE = "data/error", STATUS = 'new')
        endif
        write(1,*), "Error in replicaExchange"
        write(1,*), "I",I," test",test(I)," x",x(I)," cof",cofOld(I)    
        print*, "Error in replicaExchange"
        print*, "I",I," test",test(I)," x",x(I)," cof",cofOld(I)
        close (1)
    enddo

    ! send number bound to head node
    dest=0
    call MPI_Send(x,nTerms,MPI_DOUBLE_PRECISION,dest,0,MPI_COMM_WORLD,mc%error)
    ! send umbrella info to head node
    if (mc%umbrellaOn) then
        call MPI_Send(md%umbrellaCounts,mc%nUmbrellaBins,MPI_INTEGER,&
                      dest,0,MPI_COMM_WORLD,mc%error)
        call MPI_Send(md%umbrellaV,mc%nUmbrellaBins,MPI_DOUBLE_PRECISION,&
                      dest,0,MPI_COMM_WORLD,mc%error)
        call MPI_Send(mc%IndUmbrella,1,MPI_INTEGER,&
                      dest,0,MPI_COMM_WORLD,mc%error)
        call MPI_Send(mc%nOutside,1,MPI_INTEGER,&
                      dest,0,MPI_COMM_WORLD,mc%error)
        call MPI_Send(mc%umbBin,1,MPI_INTEGER,&
                      dest,0,MPI_COMM_WORLD,mc%error)
    endif

    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    ! send IND to head node
    if (id.eq.1) then 
        call MPI_Send(mc%IND,1,MPI_INTEGER,dest,0,MPI_COMM_WORLD,mc%error)
    endif
    ! hear back on which replica and it's mu value
    source=0
    ! get new replica number
    call MPI_Recv(mc%rep,1,MPI_INTEGER,source,0, & 
                  MPI_COMM_WORLD,status,mc%error)
    ! get new mu value
    call MPI_Recv(cof,nTerms,MPI_DOUBLE_PRECISION,source,0,&
                  MPI_COMM_WORLD,status,mc%error)
    ! recv umbrella info to head node
    if (mc%umbrellaOn) then
        call MPI_Recv(md%umbrellaCounts,mc%nUmbrellaBins,MPI_INTEGER,&
                      source,0,MPI_COMM_WORLD,status,mc%error)
        call MPI_Recv(md%umbrellaV,mc%nUmbrellaBins,MPI_DOUBLE_PRECISION,&
                      source,0,MPI_COMM_WORLD,status,mc%error)
        i=mc%IndUmbrella
        call MPI_Recv(mc%IndUmbrella,1,MPI_INTEGER,&
                      source,0,MPI_COMM_WORLD,status,mc%error)
        if (i.ne.mc%IndUmbrella) then
            print*, "Problem with PT and Umbrella, IndUmbrella !="
            stop
        endif
        call MPI_Recv(mc%nOutside,1,MPI_INTEGER,&
                      source,0,MPI_COMM_WORLD,status,mc%error)
        call calcUmbrellaE(mc,md,mc%rxnQ,mc%EUmbrella,mc%umbBin)
    endif

    mc%chi      =cof(1)
    mc%mu       =cof(2)
    mc%h_A      =cof(3)
    mc%HP1_Bind =cof(4)
    mc%Kap      =cof(5)
    !mc%Para(1)  =cof(6)
    !mc%Para(2)  =cof(7)
    !mc%Para(3)  =cof(8)

    if (abs(mc%EChi-x(1)*CofOld(1)).gt.0.0000001_dp) then
        print*, "Error in replicaExchange"
        print*, "mc%EChi",mc%EChi,"x(1)*CofOld(1)",x(1)*CofOld(1)
        stop 1
    endif

    mc%EChi    =mc%EChi    +x(1)*(Cof(1)-CofOld(1)) 
    mc%EBind   =mc%EBind   +x(2)*(Cof(2)-CofOld(2))  
    mc%EField  =mc%EField  +x(3)*(Cof(3)-CofOld(3)) 
    mc%ECouple =mc%ECouple +x(4)*(Cof(4)-CofOld(4)) 
    mc%EKap    =mc%EKap    +x(5)*(Cof(5)-CofOld(5)) 
   ! mc%EElas(1)=mc%EElas(1)+x(6)*(Cof(6)-CofOld(6)) 
   ! mc%EElas(2)=mc%EElas(2)+x(7)*(Cof(7)-CofOld(7)) 
   ! mc%EElas(3)=mc%EElas(3)+x(8)*(Cof(8)-CofOld(8)) 

    if (abs(mc%EChi-x(1)*Cof(1)).gt.0.000001_dp) then
        print*, "Error in replicaExchange"
        print*, "mc%EChi",mc%EChi,"x(1)*Cof(1)",x(1)*Cof(1)
        stop 1
    endif

    ! change output file sufix
    write(iostr,"(I4)"), mc%rep
    iostr=adjustL(iostr)
    iostr=trim(iostr)
    iostr="v"//trim(iostr)
    iostr=trim(iostr)
    mc%repSufix=iostr


    ! keep track of which thread you are
    mc%id=int(id)
end Subroutine
