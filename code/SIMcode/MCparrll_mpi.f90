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
    write ( *, '(a,i8)' ) '  The number of replicas is ', p-1
    
    
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
    !real urand(1) ! use this to generate a random number
    if (.false.) then ! set spedific seed
        Irand=7171
    else ! seed from clock
        call date_and_time(datedum,timedum,zonedum,seedvalues)
        Irand=int(-seedvalues(5)*1E7-seedvalues(6)*1E5&
                  -seedvalues(7)*1E3-seedvalues(8))
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
    double precision, allocatable :: xMtrx(:,:)  ! sum of bound states 
    double precision, allocatable :: cofMtrx(:,:) ! mu or chi or whatever
    double precision, allocatable :: s_vals(:) ! path parameter
    integer, parameter :: nTerms=8  ! number of energy terms 
    double precision x(nTerms) ! slice of xMtrx
    double precision cof(nTerms) ! slice of cofMtrx
    integer N_average      ! number of attempts since last average
    integer upSuccess(p-1)  ! number of successes since last average
    integer downSuccess(p-1) ! number of successes since last average
    integer nExchange ! total number of exchanges attemted
    type(MCvar) mc ! genaral symulation parameters
    character*16 fileName ! ouput filename
    double precision energy ! for deciding to accept exchange
    integer term ! for loopin over terms
    double precision h_path,chi_path ! functions

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
            Irand=int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
                      -seedvalues(7)*1E3-seedvalues(8))
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

        allocate( xMtrx(nPTReplicas,nTerms))
        allocate( cofMtrx(nPTReplicas,nTerms))
        allocate( nodeNumber(nPTReplicas))
        allocate( s_vals(nPTReplicas))

        do rep=1,nPTReplicas
            upSuccess(rep)=0
            downSuccess(rep)=0
            s_vals(rep)=0.01_dp*dble(rep)/dble(nPTReplicas)
        enddo

        do rep=1,nPTReplicas
            cofMtrx(rep,1)=chi_path(s_vals(rep))      
            cofMtrx(rep,2)=mc%mu     
            cofMtrx(rep,3)=h_path(s_vals(rep))     
            cofMtrx(rep,4)=mc%HP1_Bind
            cofMtrx(rep,5)=mc%EKap    
            cofMtrx(rep,6)=mc%Para(1)
            cofMtrx(rep,7)=mc%Para(2)
            cofMtrx(rep,8)=mc%Para(3) 
        enddo
        !call PT_cofValues(cof,nPTReplicas)

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
                cof=cofMtrx(rep,:)
                call MPI_Send (cof,nTerms, MPI_DOUBLE_PRECISION, dest,   0, &
                                MPI_COMM_WORLD,error )
            enddo
            ! get results from workers
            
            do rep=1,nPTReplicas
                source=nodeNumber(rep)
                call MPI_Recv ( x, nTerms, MPI_DOUBLE_PRECISION, source, 0, &
                               MPI_COMM_WORLD, status, error )
                xMtrx(rep,:)=x
                if(isnan(x(1))) then ! endo of program
                    keepGoing=.false.
                    return
                endif
            enddo
            
            ! do replica exchange
            do rep=1,(nPTReplicas-1)
                energy=0.0_dp
                do term=1,nTerms
                    energy=energy-(xMtrx(rep+1,term)-xMtrx(rep,term))*&
                                  (cofMtrx(rep+1,term)-cofMtrx(rep,term))
                enddo
                call random_number(urand,rand_stat)
                if (exp(-1.0_dp*energy).gt.urand(1)) then 
                    temp=nodeNumber(rep)
                    nodeNumber(rep)=nodeNumber(rep+1)
                    nodeNumber(rep+1)=temp
                    upSuccess(rep)=upSuccess(rep)+1
                    downSuccess(rep+1)=downSuccess(rep+1)+1
                    x=xMtrx(rep,:)
                    xMtrx(rep,:)=xMtrx(rep+1,:)
                    xMtrx(rep+1,:)=x
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

                if ((mc%IND.ge.mc%indStartRepAdapt).and. &
                    (mc%IND.lt.mc%indEndRepAdapt)) then ! insert input defined location here
                    call adaptCof(downSuccess,nPTReplicas,s_vals,N_average,&
                                   mc%lowerRepExe,mc%upperRepExe,& 
                                   mc%lowerCofRail,mc%upperCofRail,&
                                   mc%RepAnnealSpeed,mc%replicaBounds)
                    do rep=1,nPTReplicas
                        cofMtrx(rep,1)=chi_path(s_vals(rep))      
                        cofMtrx(rep,3)=h_path(s_vals(rep))     
                    enddo
                endif
                N_average=0
                do rep=1,nPTReplicas
                    upSuccess(rep)=0
                    downSuccess(rep)=0
                enddo
            endif
            nExchange=nExchange+1
        enddo

        deallocate(xMtrx)
        deallocate(cofMtrx)
        deallocate(nodeNumber)
        deallocate(s_vals)
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
        x(1)=nan_dp
        print*, "Node ",id," sending normal exit code."
        call MPI_Send(x,nTerms,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,error)
    end if


    print*, "Node ",id," exiting normally."
    return
end
function chi_path(s) result(chi)
    use setPrecision
    implicit none
    double precision, intent(in) :: s
    double precision chi
    double precision chi_max
    if (.false.) then
        chi_max=2.00_dp
        if (s.lt.0.5_dp) then
            chi=0.0_dp
        else
            chi=chi_max*2.0_dp*(s-0.5_dp)
        endif
    else
        chi=0.0_dp
    endif
end function chi_path
function h_path(s) result(h)
    use setPrecision
    implicit none
    double precision, intent(in) :: s
    double precision h
    double precision h_max
    if (.false.) then
        h_max=10.0_dp
        if (s.lt.0.5_dp) then
             h=h_max*s*2.0_dp
        else
             h=h_max*(1.0_dp-s)*2.0_dp
        endif
    else
        h=0.0_dp
    endif
end function h_path
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
    integer, parameter :: nTerms=8  ! number of energy terms 
    double precision cof(nTerms)

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nThreads,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
    
    !---------------------------------------------
    !
    !     Quenched Disorder must be same!
    !     Copy from replica 1 to others.
    !
    !----------------------------------------------
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
    endif

    !---------------------------------------------------
    !
    !   Receive instructions from head node
    !
    !----------------------------------------------------
    source=0
    dest=0
    call MPI_Recv ( mc%rep, 1, MPI_INTEGER, source, 0, &
      MPI_COMM_WORLD, status, error )
    
    call MPI_Recv ( cof, nTerms, MPI_DOUBLE_PRECISION, source, 0, &
      MPI_COMM_WORLD, status, error ) 
   
    mc%chi      =cof(1)
    mc%mu       =cof(2)
    mc%h_A      =cof(3)
    mc%HP1_Bind =cof(4)
    mc%EKap     =cof(5)
    !mc%Para(1)  =cof(6)
    !mc%Para(2)  =cof(7)
    !mc%Para(3)  =cof(8)

    write(iostrg,"(I4)"), mc%rep
    iostrg=adjustL(iostrg)
    iostrg=trim(iostrg)
    iostrg="v"//trim(iostrg)
    iostrg=trim(iostrg)
    mc%repSufix=iostrg

    ! keep track of which tread you are
    mc%id=int(id)
end Subroutine
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
    integer, parameter :: nTerms=8  ! number of energy terms 
    integer (kind=4) id, ierror
    TYPE(MCvar) mc
    integer i  ! working intiger
    integer (kind=4) dest ! message destination
    integer (kind=4) source ! message source
    character*16 iostr ! for handling sufix string
    integer status(MPI_STATUS_SIZE)  ! MPI status
    double precision cof(nTerms)
    double precision cofOld(nTerms)
    double precision x(nTerms)
    double precision test(5)

    x(1)=mc%x_Chi
    x(2)=mc%x_mu
    x(3)=mc%x_Field
    x(4)=mc%x_couple
    x(5)=mc%x_kap
    x(6)=0.0_dp !x(6)=mc%EElas(1)/mc%Para(1)
    x(7)=0.0_dp !x(7)=mc%EElas(2)/mc%Para(2)
    x(8)=0.0_dp !x(8)=mc%EElas(3)/mc%Para(3)

    test(1)=mc%EChi/mc%Chi
    test(2)=mc%EBind/mc%mu
    test(3)=mc%EField/mc%h_A
    test(4)=mc%ECouple/mc%HP1_Bind
    test(5)=mc%EKap/mc%EKap

    cofOld(1)=mc%chi      
    cofOld(2)=mc%mu     
    cofOld(3)=mc%h_A     
    cofOld(4)=mc%HP1_Bind
    cofOld(5)=mc%EKap    
    cofOld(6)=mc%Para(1)
    cofOld(7)=mc%Para(2)
    cofOld(8)=mc%Para(3) 

    do i=1,5
        if (isnan(test(I))) cycle
        if (abs(cofOld(I)).lt.0.0000001) cycle
        if (abs(test(I)-x(I)).lt.0.0000001) cycle
        print*, "Error in replicaExchange"
        print*, "I",I," test",test(I)," x",x(I)," cof",cofOld(I)
        stop 1
    enddo

    ! send number bound to head node
    dest=0
    call MPI_Send(x,nTerms,MPI_DOUBLE_PRECISION,dest,0,MPI_COMM_WORLD,mc%error)
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

    mc%chi      =cof(1)
    mc%mu       =cof(2)
    mc%h_A      =cof(3)
    mc%HP1_Bind =cof(4)
    mc%EKap     =cof(5)
    !mc%Para(1)  =cof(6)
    !mc%Para(2)  =cof(7)
    !mc%Para(3)  =cof(8)

    mc%EChi    =mc%EChi    +x(1)*(Cof(1)-CofOld(1)) 
    mc%EBind   =mc%EBind   +x(2)*(Cof(2)-CofOld(2))  
    mc%EField  =mc%EField  +x(3)*(Cof(3)-CofOld(3)) 
    mc%ECouple =mc%ECouple +x(4)*(Cof(4)-CofOld(4)) 
    mc%EKap    =mc%EKap    +x(5)*(Cof(5)-CofOld(5)) 
   ! mc%EElas(1)=mc%EElas(1)+x(6)*(Cof(6)-CofOld(6)) 
   ! mc%EElas(2)=mc%EElas(2)+x(7)*(Cof(7)-CofOld(7)) 
   ! mc%EElas(3)=mc%EElas(3)+x(8)*(Cof(8)-CofOld(8)) 

    ! change output file sufix
    write(iostr,"(I4)"), mc%rep
    iostr=adjustL(iostr)
    iostr=trim(iostr)
    iostr="v"//trim(iostr)
    iostr=trim(iostr)
    mc%repSufix=iostr


    ! keep track of which tread you are
    mc%id=int(id)
end Subroutine
