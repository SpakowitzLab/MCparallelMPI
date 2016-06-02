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
  call paraTemp ( p, id )
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
subroutine paraTemp ( p, id)

!*****************************************************************************
!
!
  use mpi
  use mersenne_twister
  use simMod

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
    double precision, allocatable :: M(:)  ! sum of bound states
    double precision, allocatable :: mu(:)
    double precision mu_value !chemical potential
      
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
        allocate ( mu(1:nPTReplicas))
        allocate ( M(1:nPTReplicas))
        allocate ( nodeNumber(1:nPTReplicas))
        !   mu(i) is the mu of the i'th physical replica,  The mu vectors doesn't change
        !   M(i) is the sum from the i'th physical replica
        !   nodeNumber(i) is the node running the i'th physical replica
        do rep=1,nPTReplicas
            mu(rep)=-2+rep*0.08
            upSuccess(rep)=0
            downSuccess(rep)=0
        enddo
        ! save mu values  
        OPEN (UNIT = 1, FILE = "data/mu", STATUS = 'NEW')
        Do rep=1,nPTReplicas
            WRITE(1,"(1f7.4)") mu(rep)
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
                call MPI_Send (mu(rep),1, MPI_DOUBLE_PRECISION, dest,   0, &
                                MPI_COMM_WORLD,error )
                
            enddo
            ! get results from workers
            
            do rep=1,nPTReplicas
                source=nodeNumber(rep)
                call MPI_Recv ( M(rep), 1, MPI_DOUBLE_PRECISION, source, 0, &
                               MPI_COMM_WORLD, status, error )
                if(M(rep).ne.M(rep)) then !test for NaN
                    keepGoing=.false.
                endif
            enddo
            !
            ! do replica exchange here
            do rep=1,(nPTReplicas-1)
                call random_number(urand,rand_stat)
                if (exp(-(M(rep+1)-M(rep))*(mu(rep+1)-mu(rep))).gt.urand(1)) then 
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
                                     mu,M,nodeNumber,N_average,nExchange)
                N_average=0
            endif
            nExchange=nExchange+1
        enddo

        deallocate (mu)
        deallocate (nodeNumber)
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
        print*, "Calling wlcsim for", mc%rep," node",id
        call wlcsim(rand_stat,id)
        print*, "Exiting from replica", mc%rep," node",id
        mu_value=0
        mu_value=1.0/mu_value !NaN
        call MPI_Send(mu_value,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,error)

    end if


    return
end
Subroutine PT_override(mc,md)
! Override initialization with parallel setup parameters
!  In particualar it changes: mc%AB, mc%rep, mc%mu, mc%repSufix
    use mpi
    use simMod
    TYPE(MCvar) mc
    TYPE(MCData) md
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    integer (kind=4) id, nThreads,ierror
    character*16 iostr    ! for file naming
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)

    mc%rep=id    
    source=0
    dest=0
    !copy AB from replica 1 to others
    if (id.eq.1) then
        do dest=2,nThreads-1
            !print*,"sending size ",mc%NT," to ",dest
            call MPI_Send (md%METH,mc%NT, MPI_INTEGER, dest,   0, &
                           MPI_COMM_WORLD,error )
        enddo
    else
        source=1
        
        print*,"receiving size ",mc%NT," at ",id
        call MPI_Recv (md%METH, mc%NT, MPI_INTEGER, source, 0, &
                       MPI_COMM_WORLD, status, error )
        source=0
    endif
    call MPI_Recv ( mc%rep, 1, MPI_INTEGER, source, 0, &
      MPI_COMM_WORLD, status, error )
    
    call MPI_Recv ( mc%mu, 1, MPI_DOUBLE_PRECISION, source, 0, &
      MPI_COMM_WORLD, status, error ) 
    
    write(iostr,"(I4)"), mc%rep
    iostr=adjustL(iostr)
    iostr=trim(iostr)
    iostr="v"//trim(iostr)
    iostr=trim(iostr)
    mc%repSufix=iostr
end Subroutine
Subroutine replicaExchange(mc,md)
! This checks in with the mpi head node to 
! 1: Tell head node the M value
! 2: Recive replica assignment from head node
! 3: Recive assigned mu value
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
    double precision mu_old

    ! Calculate M
    mc%M=0.0
    do i=1,mc%NT
       mc%M=mc%M+md%AB(i) 
    enddo
    mu_old=mc%mu
    ! send number bound to head node
    dest=0
    call MPI_Send(mc%M,1,MPI_DOUBLE_PRECISION,dest,0,MPI_COMM_WORLD,mc%error)
    ! hear back on which replica and it's mu value
    source=0
    ! get new replica number
    call MPI_Recv(mc%rep,1,MPI_DOUBLE_PRECISION,source,0, & 
                  MPI_COMM_WORLD,status,mc%error)
    ! get new mu value
    call MPI_Recv(mc%mu,1,MPI_DOUBLE_PRECISION,source,0,&
                  MPI_COMM_WORLD,status,mc%error)
    ! update energy
    mc%EBind=mc%EBind-mc%M*(mc%mu-mu_old) 
    ! change output file sufix
    write(iostr,"(I4)"), mc%rep
    iostr=adjustL(iostr)
    iostr=trim(iostr)
    iostr="v"//trim(iostr)
    iostr=trim(iostr)
    mc%repSufix=iostr
end Subroutine
Subroutine save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                           mu,M,nodeNumber,N_average,nExchange)
! Print Energy data
    IMPLICIT NONE
    LOGICAL isfile
    character*32 fullName
    integer nPTReplicas
    integer rep
    integer upSuccess(nPTReplicas)
    integer downSuccess(nPTReplicas)
    integer nodeNumber(nPTReplicas)
    double precision mu(nPTReplicas)
    double precision M(nPTReplicas)
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
        write(1,"(I7,f7.3,f7.1,2f7.3,I7)"), rep, mu(rep), M(rep), & 
                 real(upSuccess(rep))/real(N_average), &
                 real(downSuccess(rep))/real(N_average), nodeNumber(rep)
        upSuccess(rep)=0
        downSuccess(rep)=0
    enddo  

    Close(1)
end subroutine
Subroutine PT_message(iostr) 
    use mpi
    integer (kind=4) id, nThreads,ierror
    character*16 iostr    ! for file naming

    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    print*, 'PT_message',id


end subroutine
