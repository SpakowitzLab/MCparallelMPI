Subroutine PT_override(mc,md)
! Override initialization with parallel setup parameters
!  In particualar it changes: mc%AB, mc%rep, mc%mu, mc%repSufix
    use mpi
    use simMod
    Implicit none
    TYPE(MCvar), intent(inout) :: mc
    TYPE(MCData), intent(inout) :: md
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
    if (nThreads.lt.3) then
        mc%repSufix="v1"
        mc%rep=1
        mc%id=int(id)
        print*, "No PT_override. Input values used."
        return
    endif
    
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
                call MPI_Send (md%AB,mc%NT+mc%nBeadsP2, MPI_INTEGER, dest,   0, &
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
            call MPI_Recv (md%AB, mc%NT+mc%nBeadsP2, MPI_INTEGER, source, 0, &
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
    mc%Kap      =cof(5)
    !mc%Para(1)  =cof(6)
    !mc%Para(2)  =cof(7)
    !mc%Para(3)  =cof(8)

    write(iostrg,"(I4)"), mc%rep
    iostrg=adjustL(iostrg)
    iostrg=trim(iostrg)
    iostrg="v"//trim(iostrg)
    iostrg=trim(iostrg)
    mc%repSufix=iostrg

    ! keep track of which thread you are
    mc%id=int(id)
end Subroutine
