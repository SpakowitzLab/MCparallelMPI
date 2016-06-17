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

    integer ( kind = 4 ) dest   !destination id for messages
    integer ( kind = 4 ) source  !source id for messages
    integer ( kind = 4 ) error  ! error id for MIP functions
    integer ( kind = 4 ) id     ! which processor I am
    integer ( kind = 4 ) p ! number of threads
    integer ( kind = 4 ) status(MPI_STATUS_SIZE) ! MPI stuff
    integer nPTReplicas     ! number of replicas.

! now siulation variables
    type(MCvar) mc
    type(MCdata) md
    integer temp ! for castling
    logical keepGoing   ! set to false when NaN encountered
    integer rep ! physical replica number
    integer, allocatable :: nodeNumber(:)
    double precision, allocatable :: x(:)  ! sum of bound states
    double precision, allocatable :: cof(:)
    double precision nan_dp !chemical potential
      
!   variable for random number generator seeding
    type(random_stat) rand_stat  ! state of random number chain
    integer Irand     ! Seed
    character*8 datedum  ! trash
    character*10 timedum ! trash
    character*5 zonedum  ! trash
    integer seedvalues(8) ! clock readings
    real urand(1)

!   traking variables
    integer N_average
    integer upSuccess(p-1)
    integer downSuccess(p-1)
    integer NT
    integer nExchange
    nExchange=0

    nPTReplicas=p-1 
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
            !
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
            N_average=N_average+1
            if (N_average.gt.4999) then
                call save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                                     cof,x,nodeNumber,N_average,nExchange)
                N_average=0
            endif
            nExchange=nExchange+1
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
            call MPI_Send (md%METH,mc%NT, MPI_INTEGER, dest,   0, &
                           MPI_COMM_WORLD,error )
        enddo
    else
        source=1
        
        call MPI_Recv (md%METH, mc%NT, MPI_INTEGER, source, 0, &
                       MPI_COMM_WORLD, status, error )
        source=0
    endif
    call MPI_Recv ( mc%rep, 1, MPI_INTEGER, source, 0, &
      MPI_COMM_WORLD, status, error )
    
    call MPI_Recv ( cof, 1, MPI_DOUBLE_PRECISION, source, 0, &
      MPI_COMM_WORLD, status, error ) 
   
    ! set cof
    !mc%mu=cof
    mc%chi=cof

    write(iostrg,"(I4)"), mc%rep
    iostrg=adjustL(iostrg)
    iostrg=trim(iostrg)
    iostrg="v"//trim(iostrg)
    iostrg=trim(iostrg)
    mc%repSufix=iostrg
end Subroutine
Subroutine PT_cofValues(cof,nPTReplicas)
    use setPrecision
    Implicit none
    Integer nPTReplicas 
    Double precision cof(nptReplicas)
    INteger rep
    do rep=1,nPTReplicas
        !cof(rep)=2.0_dp-rep*0.08_dp  !over mu values
        cof(rep)=0.0001_dp+0.06_dp*rep
    enddo
end subroutine
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
    TYPE(MCvar) mc
    TYPE(MCData) md
    integer i  ! working intiger
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    character*16 iostr ! for handling sufix string
    integer status(MPI_STATUS_SIZE)  ! MPI status
    double precision cof,x
    double precision mu_old
    double precision chi_old

    ! Calculate value conjagate to ajusted parameter
    !mc%M=mc%EBind/(-1.0_dp*mc%mu)
    !x=mc%M
    !mu_old=mc%mu
    x=mc%EChi/(mc%Chi)  ! sum (Vol/V)*PHIA*PHIB
    chi_old=mc%chi


    if (x.ne.x) then
        print*, "Error in replicaExchange! NaN encountered"
        stop 1
    endif

    ! send number bound to head node
    dest=0
    call MPI_Send(x,1,MPI_DOUBLE_PRECISION,dest,0,MPI_COMM_WORLD,mc%error)
    ! hear back on which replica and it's mu value
    source=0
    ! get new replica number
    call MPI_Recv(mc%rep,1,MPI_INTEGER,source,0, & 
                  MPI_COMM_WORLD,status,mc%error)
    ! get new mu value
    call MPI_Recv(cof,1,MPI_DOUBLE_PRECISION,source,0,&
                  MPI_COMM_WORLD,status,mc%error)

    ! update energy
    !mc%mu=-cof
    !mc%EBind=mc%EBind-mc%M*(mc%mu-mu_old)
    mc%chi=cof
    mc%EChi=mc%EChi+x*(mc%chi-chi_old)
    

    ! change output file sufix
    write(iostr,"(I4)"), mc%rep
    iostr=adjustL(iostr)
    iostr=trim(iostr)
    iostr="v"//trim(iostr)
    iostr=trim(iostr)
    mc%repSufix=iostr
end Subroutine
Subroutine save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                           cof,x,nodeNumber,N_average,nExchange)
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


    fullName=  'data/repHistory'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif

    write(1,*) "~~~~~~~~~~~~~~~~~",nExchange,"~~~~~~~~~~~~~~~~~~~~~~~~"
    write(1,*) "  rep |  mu  |   M  |  up  | down |  node"
    do rep=1,nPTReplicas
        write(1,"(I7,f7.3,f7.1,2f7.3,I7)"), rep, -cof(rep), x(rep), & 
                 real(upSuccess(rep))/real(N_average), &
                 real(downSuccess(rep))/real(N_average), nodeNumber(rep)
        upSuccess(rep)=0
        downSuccess(rep)=0
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
