program main

!*****************************************************************************
!
!  This program is a MTMC with parallel tempering.
!  The MC code is in wlcsim 
!
!
!  Modified:   5/6/2016
!
!  Author:  Quinn MacPherson
!
!*****************************************************************************
  use mpi
  use setPrecision

  implicit none

  integer ( kind = 4 ) error
  integer ( kind = 4 ) id
  integer ( kind = 4 ) p


!
!  Initialize MPI.
!
  call MPI_Init ( error )
  if (error.ne.0) then
      print*, "MPI_Init", error
  endif
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
  if (error.ne.0) then
      print*, "MPI_Comm_size", error
  endif
!
!  Get the individual process ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
  if (error.ne.0) then
      print*, "MPI_Comm_rank", error
  endif
!
!  Print a message.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  MC with tempering using MPI'
    write ( *, '(a)' ) '  FORTRAN90/MPI version'
    write ( *, '(a)' ) '  This is a basic parallel tempering MC sim'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of threads being used is ', p
    write ( *, '(a,i8)' ) '  The number of temperatues is ', p-1
    
    
  end if
  if (p.gt.1) then
      call paraTemp ( p, id )
  else
      call singleCall()
  endif
!
!  Shut down MPI.
!
  call MPI_Finalize ( error )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
  end if

  stop
end
subroutine singleCall()
    use mersenne_twister
    implicit none
!   variable for random number generator seeding
    type(random_stat) rand_stat  ! state of random number chain
    integer Irand     ! Seed
    character*8 datedum  ! trash
    character*10 timedum ! trash
    character*5 zonedum  ! trash
    integer seedvalues(8) ! clock readings
    real urand(1)
    if (.false.) then ! set spedific seed
        Irand=7171
    else ! seed from clock
        call date_and_time(datedum,timedum,zonedum,seedvalues)
        Irand=-seedvalues(5)*1E7-seedvalues(6)*1E5-seedvalues(7)*1E3-seedvalues(8)
        Irand=mod(Irand,10000)
        print*, "Random Intiger seed:",Irand
    endif
    call random_setseed(Irand,rand_stat) ! random seed for head node
    print*, "calling single wlcsim"
    call wlcsim(rand_stat)
end subroutine
subroutine paraTemp ( p, id)

!*****************************************************************************
!
!
  use mpi
  use mersenne_twister
  use simMod
  use setPrecision

    implicit none

!   MPI variables
    integer ( kind = 4 ) dest   !destination id for messages
    integer ( kind = 4 ) source  !source id for messages
    integer ( kind = 4 ) error  ! error id for MIP functions
    integer ( kind = 4 ) id     ! which processor I am
    integer ( kind = 4 ) p ! number of threads
    integer ( kind = 4 ) status(MPI_STATUS_SIZE) ! MPI stuff
    integer nPTReplicas     ! number of replicas.

!   variable for random number generator seeding
    type(random_stat) rand_stat  ! state of random number chain
    integer Irand     ! Seed
    character*8 datedum  ! trash
    character*10 timedum ! trash
    character*5 zonedum  ! trash
    integer seedvalues(8) ! clock readings
    real urand(1)

!   worker node only variables
    double precision nan_dp !chemical potential

!   for head node use only variables
    integer rep ! physical replica number, for loops
    integer temp ! for castling
    logical keepGoing   ! set to false when NaN encountered
    integer, allocatable :: nodeNumber(:)  ! list of which nodes are which
    double precision, allocatable :: x(:)  ! sum of bound states 
    double precision, allocatable :: cof(:) ! mu or chi or whatever
    integer N_average      ! number of attempts since last average
    integer upSuccess(p-1)  ! number of successes since last average
    integer downSuccess(p-1) ! number of successes since last average
    integer nExchange ! total number of exchanges attemted
    type(MCvar) mc ! genaral symulation parameters
    character*16 fileName

    nExchange=0
    nPTReplicas=p-1 
    fileName='input/params'
!  -----------------------------------------------------------------
!
!      Head Node
!
!  -----------------------------------------------------------------
    if ( id == 0 ) then
        ! -----------------------------------------------
        !
        !   Generate thread safe random number seeds
        !
        !--------------------------------------------
        if (.false.) then ! set spedific seed
            Irand=7171
        else ! seed from clock
            call date_and_time(datedum,timedum,zonedum,seedvalues)
            Irand=-seedvalues(5)*1E7-seedvalues(6)*1E5-seedvalues(7)*1E3-seedvalues(8)
            Irand=mod(Irand,10000)
            print*, "Random Intiger seed:",Irand
        endif
        call random_setseed(Irand*(id+1),rand_stat) ! random seed for head node
        do dest=1,nPTReplicas ! send out the others
            call MPI_Send (Irand,1, MPI_INTEGER, dest,   0, &
                            MPI_COMM_WORLD,error )
        enddo  
        ! ------------------------- 
        !
        !   innitialize
        !
        ! --------------------------
        nPTReplicas = p-1;

        call MCvar_setParams(mc,fileName) ! so that the head thread knows the  parameters

        allocate( x(1:nPTReplicas))
        allocate( cof(1:nPTReplicas))
        allocate( nodeNumber(1:nPTReplicas))
        do rep=1,nPTReplicas
            upSuccess(rep)=0
            downSuccess(rep)=0
        enddo
        call PT_cofValues(cof,nPTReplicas)

        ! save cof values  
        OPEN (UNIT = 1, FILE = "data/cof", STATUS = 'NEW')
        Do rep=1,nPTReplicas
            WRITE(1,"(1f7.4)") cof(rep)
        enddo
        Close(1)

        N_average=0

        ! Initially replica numbers are same as nodes
        do rep=1,nPTReplicas
            nodeNumber(rep)=rep
        enddo
        !-----------------------------------------------------
        !
        !    Begin Main loop
        !
        !------------------------------------------------------
        keepGoing=.True.
        do while(keepGoing)
            ! give workers thier jobs
            do rep=1,nPTReplicas
                dest=nodeNumber(rep)
                call MPI_Send (rep,1, MPI_INTEGER, dest,   0, &
                                MPI_COMM_WORLD,error )
                call MPI_Send (cof(rep),1, MPI_DOUBLE_PRECISION, dest,   0, &
                                MPI_COMM_WORLD,error )
            enddo
            ! get results from workers
            
            do rep=1,nPTReplicas
                source=nodeNumber(rep)
                call MPI_Recv ( x(rep), 1, MPI_DOUBLE_PRECISION, source, 0, &
                               MPI_COMM_WORLD, status, error )
                if(x(rep).ne.x(rep)) then !test for NaN
                    keepGoing=.false.
                endif
            enddo
            
            ! do replica exchange here
            do rep=1,(nPTReplicas-1)
                call random_number(urand,rand_stat)
                if (exp((x(rep+1)-x(rep))*(cof(rep+1)-cof(rep))).gt.urand(1)) then 
                    temp=nodeNumber(rep)
                    nodeNumber(rep)=nodeNumber(rep+1)
                    nodeNumber(rep+1)=temp
                    upSuccess(rep)=upSuccess(rep)+1
                    downSuccess(rep+1)=downSuccess(rep+1)+1
                endif
            enddo

            source=1 
            call MPI_Recv (mc%Ind, 1, MPI_INTEGER, source, 0, &
                           MPI_COMM_WORLD, status, error )
            
            ! track/adapt acceptance rates
            N_average=N_average+1
            if (N_average.ge.mc%NRepAdapt) then
                call save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                                     cof,x,nodeNumber,N_average,nExchange,mc%IND)
                if (mc%IND.gt.mc%indStartRepAdapt) then ! insert input defined location here
                    call adaptCof(upSuccess,nPTReplicas,cof,N_average,&
                                   mc%lowerRepExe,mc%upperRepExe,mc%lowerCofRail,mc%upperCofRail)
                endif
                N_average=0
            endif
            do rep=1,nPTReplicas
                nExchange=nExchange+1
                upSuccess(rep)=0
                downSuccess(rep)=0
            enddo
        enddo

        deallocate (cof)
        deallocate (x)
    else
!  -----------------------------------------------------------------
!
!      Worker Node
!
!  -----------------------------------------------------------------

        source = 0
        dest = 0
        ! -----------------------------------------------
        !
        !   Generate thread safe random number chain: rand_stat
        !
        !--------------------------------------------
        
        call MPI_Recv ( Irand, 1, MPI_INTEGER, source, 0, &
                        MPI_COMM_WORLD, status, error )
        call random_setseed(Irand*(id+1),rand_stat) ! random seed for head node
        ! ------------------------------
        !
        ! call main simulation code
        !
        !  --------------------------------
        call wlcsim(rand_stat)
        nan_dp=0; nan_dp=nan_dp/nan_dp !NaN
        call MPI_Send(nan_dp,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,error)
    end if


    return
end
Subroutine PT_override(mc,md)
! Override initialization with parallel setup parameters
!  In particualar it changes: mc%AB, mc%rep, mc%mu, mc%repSufix
    use mpi
    use simMod
    Implicit none
    TYPE(MCvar) mc
    TYPE(MCData) md
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    integer (kind=4) id, nThreads,ierror
    integer (kind=4) error  ! error id for MIP functions
    character*16 iostrg    ! for file naming
    integer ( kind = 4 ) status(MPI_STATUS_SIZE) ! MPI stuff
    double precision cof

    
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    
    mc%rep=id    
    source=0
    dest=0
    !copy AB from replica 1 to others
    if (id.eq.1) then
        do dest=2,nThreads-1
            if(mc%simType.eq.1) then
                call MPI_Send (md%METH,mc%NT, MPI_INTEGER, dest,   0, &
                               MPI_COMM_WORLD,error )
            elseif(mc%simType.eq.0) then
                call MPI_Send (md%AB,mc%NT, MPI_INTEGER, dest,   0, &
                               MPI_COMM_WORLD,error )
            else
                print*, "Error in PT_override. simType doesn't exist."
                stop 1
            endif
        enddo
    else
        source=1
        if(mc%simType.eq.1) then
            call MPI_Recv (md%METH, mc%NT, MPI_INTEGER, source, 0, &
                           MPI_COMM_WORLD, status, error )
        elseif(mc%simType.eq.0) then
            call MPI_Recv (md%AB, mc%NT, MPI_INTEGER, source, 0, &
                           MPI_COMM_WORLD, status, error )
        else
            print*, "Error in PT_override. simType doesn't exist."
            stop 1
        endif
        source=0
    endif
    call MPI_Recv ( mc%rep, 1, MPI_INTEGER, source, 0, &
      MPI_COMM_WORLD, status, error )
    
    call MPI_Recv ( cof, 1, MPI_DOUBLE_PRECISION, source, 0, &
      MPI_COMM_WORLD, status, error ) 
   
    ! set cof
    if(mc%simType.eq.1) then
        mc%mu=cof
    elseif(mc%simType.eq.0) then
        mc%chi=cof
    else
        print*, "Error in PT_override. simType doesn't exist."
        stop 1
    endif

    write(iostrg,"(I4)"), mc%rep
    iostrg=adjustL(iostrg)
    iostrg=trim(iostrg)
    iostrg="v"//trim(iostrg)
    iostrg=trim(iostrg)
    mc%repSufix=iostrg
end Subroutine
Subroutine PT_cofValues(cof,nPTReplicas)
    use setPrecision
    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
    Implicit none
    Integer nPTReplicas 
    Double precision cof(nptReplicas)
    INteger rep
    ! IO variables
    character*16 fileName  ! file with parameters
    INTEGER :: PF   ! input file unit
    LOGICAL :: FILEEND=.FALSE. ! done reading file?
    CHARACTER*100 :: WORD ! keyword
    INTEGER :: NITEMS ! number of items on the line in the parameter file
    ! outputs
    Double precision gap
    Double precision minCof
    ! -----------------------
    !
    !  Read from file
    !
    !-------------------------
    fileName='input/RepSetting'
    PF=55
    OPEN(UNIT=PF,FILE=fileName,STATUS='OLD') 

    ! read in the keywords one line at a time
    DO 
       CALL READLINE(PF,FILEEND,NITEMS)
       IF (FILEEND.and.nitems.eq.0) EXIT

       ! skip empty lines
       IF (NITEMS.EQ.0) CYCLE

       ! Read in the keyword for this line
       CALL READA(WORD,CASESET=1)

       ! Skip any empty lines or any comment lines
       IF (WORD(1:1).EQ.'#') CYCLE

       SELECT CASE(WORD) ! pick which keyword
       CASE('MIN')
           Call READF(minCof)
       CASE('GAP')
           Call READF(gap)
       CASE DEFAULT
           print*, "Error in MCvar_setParams.  Unidentified keyword:", &
                   TRIM(WORD)
           stop 1
       ENDSELECT
    ENDDO

    print*, "MIN=",minCof," GAP=",GAP
    do rep=1,nPTReplicas
        !cof(rep)=2.0_dp-rep*0.08_dp  !over mu values
        cof(rep)=minCof+Gap*(rep-1)
    enddo
    print*, "cof values:"
    print*, cof
    Close(PF)
    return
end subroutine
Subroutine replicaExchange(mc)
! This checks in with the mpi head node to 
! For parallel tempering of the form:  E=cof*x
! 1: Tell head node the x value
! 2: Recive replica assignment from head node
! 3: Recive assigned cof value
    use setPrecision
    use mpi
    use simMod
    IMPLICIT NONE
    integer (kind=4) id, ierror
    TYPE(MCvar) mc
    integer i  ! working intiger
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    character*16 iostr ! for handling sufix string
    integer status(MPI_STATUS_SIZE)  ! MPI status
    double precision cof,x
    double precision mu_old
    double precision chi_old

    ! Calculate value conjagate to ajusted parameter
    if(mc%simType.eq.1) then
        mc%M=mc%EBind/(-1.0_dp*mc%mu)
        x=mc%M
        mu_old=mc%mu
    elseif(mc%simType.eq.0) then
        x=mc%EChi/(mc%Chi)  ! sum (Vol/V)*PHIA*PHIB
        chi_old=mc%chi
    else
        print*, "Error in PT_override. simType doesn't exist."
        stop 1
    endif


    if (x.ne.x) then
        print*, "Error in replicaExchange! NaN encountered"
        print*, "Chi:"
        print*, mc%Chi
        print*, "x:"
        print*, x
        stop 1
    endif

    ! send number bound to head node
    dest=0
    call MPI_Send(x,1,MPI_DOUBLE_PRECISION,dest,0,MPI_COMM_WORLD,mc%error)
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
    call MPI_Recv(cof,1,MPI_DOUBLE_PRECISION,source,0,&
                  MPI_COMM_WORLD,status,mc%error)

    ! update energy
    if(mc%simType.eq.1) then
        mc%mu=-cof
        mc%EBind=mc%EBind-mc%M*(mc%mu-mu_old)
    elseif(mc%simType.eq.0) then
        mc%chi=cof
        mc%EChi=mc%EChi+x*(mc%chi-chi_old)
    endif
    

    ! change output file sufix
    write(iostr,"(I4)"), mc%rep
    iostr=adjustL(iostr)
    iostr=trim(iostr)
    iostr="v"//trim(iostr)
    iostr=trim(iostr)
    mc%repSufix=iostr
end Subroutine
Subroutine save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                           cof,x,nodeNumber,N_average,nExchange,IND)
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


    fullName=  'data/repHistory'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif

    write(1,*) "~~~~~~~~~~~exchange: ",nExchange,", IND:",IND,"~~~~~~~~~~~~~~~~~~~~"
    write(1,*) " rep |  cof  |   x    |  up  | down |node"
    do rep=1,nPTReplicas
        write(1,"(I6,f8.3,f9.1,2f7.3,I4)"), rep, cof(rep), x(rep), & 
                 real(upSuccess(rep))/real(N_average), &
                 real(downSuccess(rep))/real(N_average), nodeNumber(rep)
    enddo  

    Close(1)
end subroutine
Subroutine PT_message(iostr) 
    use mpi
    implicit none
    integer (kind=4) id, nThreads,ierror
    character*16 iostr    ! for file naming

    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    print*, 'PT_message',id


end subroutine
