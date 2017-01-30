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
    write ( *, '(a)' ) '  WLC MC sim with tempering using MPI'
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
   !call simpleSim(rand_stat)
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
    integer ( kind = 4 ), intent(in) :: id     ! which processor I am
    integer ( kind = 4 ), intent(in) :: p ! number of threads
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
    double precision h_path,chi_path,mu_path,kap_path,HP1_Bind_path ! functions

!   Umbrella variables for head node
    double precision, allocatable :: umbrellaVMtrx(:,:)   
    integer, allocatable :: umbrellaCountsMtrx(:,:)
    integer, allocatable :: nOutsideVec(:)   
    integer, allocatable :: IndUmbrellaVec(:)
    integer, allocatable :: umbBinVec(:)
    INTEGER, ALLOCATABLE, DIMENSION(:):: umbrellaCounts
    Real(dp), ALLOCATABLE, DIMENSION(:):: umbrellaV
    

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
        !   Umbrella variables for head node
        allocate( umbrellaVMtrx(nPTReplicas,mc%nUmbrellaBins) )
        allocate( umbrellaCountsMtrx(nPTReplicas,mc%nUmbrellaBins) )
        allocate( nOutsideVec(nPTReplicas) )
        allocate( IndUmbrellaVec(nPTReplicas) )
        Allocate (umbBinVec(mc%nUmbrellaBins))
        Allocate (umbrellaCounts(mc%nUmbrellaBins))
        Allocate (umbrellaV(mc%nUmbrellaBins))

        do rep=1,nPTReplicas
            upSuccess(rep)=0
            downSuccess(rep)=0
            s_vals(rep)=mc%INITIAL_MAX_S*dble(rep)/dble(nPTReplicas)
        enddo

        do rep=1,nPTReplicas
            if (mc%PT_chi) then
                cofMtrx(rep,1)=chi_path(s_vals(rep))      
            else
                cofMtrx(rep,1)=mc%chi
            endif
            if (mc%PT_mu) then
                cofMtrx(rep,2)=mu_path(s_vals(rep))      
            else
                cofMtrx(rep,2)=mc%mu
            endif
            if (mc%PT_h) then
                cofMtrx(rep,3)=h_path(s_vals(rep))
            else
                cofMtrx(rep,3)=mc%h_A
            endif
            if (mc%PT_couple) then
                cofMtrx(rep,4)=HP1_Bind_path(s_vals(rep))
            else
                cofMtrx(rep,4)=mc%HP1_Bind
            endif
            if (mc%PT_Kap) then
                cofMtrx(rep,5)=kap_path(s_vals(rep))
            else
                cofMtrx(rep,5)=mc%KAP
            endif
            cofMtrx(rep,6)=mc%Para(1)
            cofMtrx(rep,7)=mc%Para(2)
            cofMtrx(rep,8)=mc%Para(3) 
        enddo

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
        mc%IND=0
        keepGoing=.True.
        do while(keepGoing)
            ! give workers thier jobs
            do rep=1,nPTReplicas
                dest=nodeNumber(rep)
                call MPI_Send (rep,1, MPI_INTEGER, dest,   0, &
                                MPI_COMM_WORLD,error )
                if (mc%restart.and.mc%IND.eq.0) then
                    source=dest
                    call MPI_Recv (cof, nTerms, MPI_DOUBLE_PRECISION, source, 0, &
                                   MPI_COMM_WORLD, status, error ) 
                    cofMtrx(rep,:)=cof
                else
                    cof=cofMtrx(rep,:)
                    call MPI_Send (cof,nTerms, MPI_DOUBLE_PRECISION, dest,   0, &
                                    MPI_COMM_WORLD,error )
                endif

                ! stuff for umbrella sampling
                if (mc%umbrellaOn .and. mc%IND.ne.0) then
                    umbrellaV = umbrellaVMtrx(rep,:) 
                    umbrellaCounts = umbrellaCountsMtrx(rep,:)
                    mc%nOutside = nOutsideVec(rep)
                    mc%IndUmbrella = IndUmbrellaVec(rep) 

                    call MPI_Send(umbrellaCounts,mc%nUmbrellaBins,MPI_INTEGER,&
                                  dest,0,MPI_COMM_WORLD,mc%error)
                    call MPI_Send(umbrellaV,mc%nUmbrellaBins,MPI_DOUBLE_PRECISION,&
                                  dest,0,MPI_COMM_WORLD,mc%error)
                    call MPI_Send(mc%IndUmbrella,1,MPI_INTEGER,&
                                  dest,0,MPI_COMM_WORLD,mc%error)
                    call MPI_Send(mc%nOutside,1,MPI_INTEGER,&
                                  dest,0,MPI_COMM_WORLD,mc%error)
                endif
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
                
                ! Receive Umbrella stuff
                if (mc%umbrellaOn) then
                    call MPI_Recv(umbrellaCounts,mc%nUmbrellaBins,MPI_INTEGER,&
                                  source,0,MPI_COMM_WORLD,status,mc%error)
                    call MPI_Recv(umbrellaV,mc%nUmbrellaBins,MPI_DOUBLE_PRECISION,&
                                  source,0,MPI_COMM_WORLD,status,mc%error)
                    call MPI_Recv(mc%IndUmbrella,1,MPI_INTEGER,&
                                  source,0,MPI_COMM_WORLD,status,mc%error)
                    call MPI_Recv(mc%nOutside,1,MPI_INTEGER,&
                                  source,0,MPI_COMM_WORLD,status,mc%error)
                    call MPI_Recv(mc%umbBin,1,MPI_INTEGER,&
                                  source,0,MPI_COMM_WORLD,status,mc%error)
                    umbrellaVMtrx(rep,:) = umbrellaV
                    umbrellaCountsMtrx(rep,:) = umbrellaCounts
                    nOutsideVec(rep) = mc%nOutside
                    IndUmbrellaVec(rep) = mc%IndUmbrella
                    umbBinVec(rep)=mc%umbBin
                endif

            enddo
            
            source=1 
            call MPI_Recv (mc%Ind, 1, MPI_INTEGER, source, 0, &
                           MPI_COMM_WORLD, status, error )



            ! do replica exchange
            do rep=1,(nPTReplicas-1)
                energy=0.0_dp
                do term=1,nTerms
                    energy=energy-(xMtrx(rep+1,term)-xMtrx(rep,term))*&
                                  (cofMtrx(rep+1,term)-cofMtrx(rep,term))
                enddo

                ! change in umbrella energy
                if (mc%umbrellaOn) then
                    if (umbBinVec(rep+1).ne.-1) then ! test to see if outside range
                        energy=energy &
                            +umbrellaVMtrx(rep,umbBinVec(rep+1))&
                            -umbrellaVMtrx(rep+1,umbBinVec(rep+1))
                    endif
                    if (umbBinVec(rep).ne.-1) then
                        energy=energy &
                            +umbrellaVMtrx(rep+1,umbBinVec(rep))&
                            -umbrellaVMtrx(rep,umbBinVec(rep))
                    endif
                endif

                call random_number(urand,rand_stat)
                if (exp(-1.0_dp*energy).gt.urand(1)) then 
                    if (mc%PTON) then
                        temp=nodeNumber(rep)
                        nodeNumber(rep)=nodeNumber(rep+1)
                        nodeNumber(rep+1)=temp
                    endif
                    upSuccess(rep)=upSuccess(rep)+1
                    downSuccess(rep+1)=downSuccess(rep+1)+1
                    x=xMtrx(rep,:)
                    xMtrx(rep,:)=xMtrx(rep+1,:)
                    xMtrx(rep+1,:)=x
                endif
            enddo

            
            ! track/adapt acceptance rates
            N_average=N_average+1
            if (N_average.ge.mc%NRepAdapt) then
                call save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                                     cofMtrx,xMtrx,nodeNumber,N_average,&
                                     nExchange,mc%IND,nTerms,s_vals)

                if ((mc%IND.ge.mc%indStartRepAdapt).and. &
                    (mc%IND.lt.mc%indEndRepAdapt)) then ! insert input defined location here
                    call adaptCof(downSuccess,nPTReplicas,s_vals,N_average,&
                                   mc%lowerRepExe,mc%upperRepExe,& 
                                   mc%lowerCofRail,mc%upperCofRail,&
                                   mc%RepAnnealSpeed,mc%replicaBounds)
                    do rep=1,nPTReplicas
                        if (mc%PT_chi) then
                            cofMtrx(rep,1)=chi_path(s_vals(rep))      
                        endif
                        if (mc%PT_mu) then
                            cofMtrx(rep,2)=mu_path(s_vals(rep))      
                        endif
                        if (mc%PT_h) then
                            cofMtrx(rep,3)=h_path(s_vals(rep))
                        endif
                        if (mc%PT_couple) then
                            cofMtrx(rep,4)=HP1_Bind_path(s_vals(rep))
                        endif
                        if (mc%PT_Kap) then
                            cofMtrx(rep,5)=kap_path(s_vals(rep))
                        endif
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
        !call simpleSim(rand_stat)
        nan_dp=0; nan_dp=nan_dp/nan_dp !NaN
        x(1)=nan_dp
        print*, "Node ",id," sending normal exit code."
        call MPI_Send(x,nTerms,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,error)
    end if


    print*, "Node ",id," exiting normally."
    return
end
