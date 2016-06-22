program main

!*****************************************************************************80
!
!! MAIN is the main program for RING_MPI.
!
!  Discussion:
!
!  This program is a MTMC with parallel tempering
!  simulation for a single particle in a double well.
!  This problem can easily  be solved analytically.
! 
!
!
!  Modified:   5/6/2016
!
!  Author:  Quinn MacPherson
!
!
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
    write ( *, '(a)' ) ' basicPT_mpi:'
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
    write ( *, '(a)' ) 'basicPT_mpi'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
  end if

  stop
end
subroutine paraTemp ( p, id )

!*****************************************************************************80
!
!
  use mpi
  use mersenne_twister

  implicit none

  integer ( kind = 4 ) dest   !destination id for messages
  integer ( kind = 4 ) source  !source id for messages
  integer ( kind = 4 ) error  ! error id for MIP functions
  integer ( kind = 4 ) id     ! which processor I am
  integer ( kind = 4 ) j,i      !working inter for loops
  integer nPTReplicas     ! number of replicas.  Must equal p
  double precision, allocatable :: Energy(:)  ! energy of replics
  integer, allocatable ::  Temp(:)  ! temperatures of replicas
  integer, allocatable ::  TempID(:)  ! temperatures of replicas
  integer ( kind = 4 ) p ! number of threads
  integer ( kind = 4 ) status(MPI_STATUS_SIZE) ! MPI stuff
  type(random_stat), allocatable:: stat(:) !for random numer generator
  real urand(1) ! generic random number
  INTEGER IOStatus ! Status of file IO

  ! now siulation variables
  
  integer nr  ! number replica exchanges
  integer nbr  ! number steps between exchanges 
  double precision TRatio ! ratio between ajacent temperatures
  double precision a  ! variable in energy landscape
  double precision b  ! variable in energy landscape
  double precision x  ! position
  double precision xTrial  ! trial position
  double precision stepSize ! gaussian step size
  double precision E ! energy
  double precision T ! single Temperature
  double precision ETrial ! energy
  double precision holdT ! Temparary storage of Temp
  double precision holdTID ! Temparary storage of TempID
  integer step ! for counting steps
  integer numOne
  
!  !
!  ! Now make a structure
!  !
!  integer disp(3) ! for measureing DataStr
!  integer offsets(3)    ! for measureing DataStr 
!  integer blockCounts(3)    ! for measureing DataStr 
!  integer oldtypes(3)    ! for measureing DataStr 
!  integer MPI_Data    ! type of myData 
!  type DataStr
!    double precision :: T
!    double precision :: x
!    integer :: id
!  endtype
!  type (DataStr) myData
!    call mpi_get_address(myData%T, disp(0), error)
!    call mpi_get_address(myData%x, disp(2), error)
!    call mpi_get_address(myData%id, disp(3), error)
!    offsets(1)=0
!    offsets(2)=disp(2)-disp(1)
!    offsets(3)=disp(3)-disp(2)
!    blockcounts(1)=1;
!    blockcounts(2)=1;
!    blockcounts(3)=1;
!    oldtypes(1)=MPI_DOUBLE_PRECISION
!    oldtypes(2)=MPI_DOUBLE_PRECISION
!    oldtypes(3)=MPI_integer
!    call MPI_TYPE_CREATE_STRUCT(3,blockcounts,offsets,oldtypes,MPI_Data,error)
!    print*,'Struct_error=',error
!    call MPI_TYPE_COMMIT(MPI_Data,error)

  

    ! set variables
    nPTReplicas = p-1; ! make this an imput variable
    allocate ( Temp(1:nPTReplicas) )
    allocate ( Energy(1:nPTReplicas) )
    allocate ( stat(1:p) ) 
    allocate ( TempID(1:nPTReplicas))
    numOne =1;
    ! simulation variables
    nr=50000
    nbr=1000
    Tratio=2.0
    a=-7.0
    b=1.0
    x= 1.8 ! initial condition
    xTrial=x
    stepSize=0.4

    do j=1,p
        ! change the below line if you would like a random seed
        call random_setseed(7173*j,stat(j)) 
        ! this is a bit wastefull in that every thread generates
        ! every seed but only uses one of them.
    enddo
    !urand= random_number_f()

    if (nPTReplicas.ne.p-1) then
        print*, 'Error: p is not nPTReplicas!!!!!!!!!!'
    endif

    if ( id == 0 ) then
        OPEN(unit=1,file='data/xdata',IOSTAT=IOStatus,status='new')
        ! initialize temperature identities
        do j=1,nPTReplicas
            Temp(j)=Tratio**real(j-1)
            TempID(j)=j
        enddo
        do j=1, nr
            ! print*, j
            ! send out jobs
            do dest = 1,(p-1)
                T=Temp(dest)
                call MPI_Send ( T, 1, MPI_DOUBLE_PRECISION, dest,   0, &
                  MPI_COMM_WORLD,         error )
            enddo
            ! hear back from jobs
            do source = 1,(p-1)
                call MPI_Recv ( E, 1, MPI_DOUBLE_PRECISION, source, 0, &
                  MPI_COMM_WORLD, status, error )
                Energy(source)=E

                ! Here I bring the x data back to head node
                call MPI_Recv ( x, 1, MPI_DOUBLE_PRECISION, source, 0, &
                  MPI_COMM_WORLD, status, error )
                if (TempID(source).eq.1) then
                    WRITE(1,"(f10.5)") x
                endif
            enddo
            ! reorder jobs
            do i=1,(p-2)
                call random_number(urand,stat(id))
                if (exp((Energy(i+1)-Energy(i))*(1.0/Temp(i+1) - 1.0/Temp(i))).gt.urand(1)) then
                    !Swop Temperatures
                    holdT=Temp(i+1)
                    Temp(i+1)=Temp(i)
                    Temp(i)=holdT
                    holdTID=TempID(i+1)
                    TempID(i+1)=TempID(i)
                    TempID(i)=holdTID
                endif
            enddo


        enddo

        CLOSE(unit=1,IOSTAT=IOStatus)
!  Worker ID must receive first from ID-1, then send to ID+1.
!
    else

      source = 0
      dest = 0
 
      do j=1,nr
        call MPI_Recv ( T, 1, MPI_DOUBLE_PRECISION, source, 0, &
          MPI_COMM_WORLD, status, error )
        ! if calculating energy is expensive modify the code to only do it the first time  
        E=a*x*x+b*x**4
        do step=1,nbr

            ! calculate move
            call random_gauss(urand,stat(id))
            xTrial = x+stepSize*urand(1)
            
            ! calculate trial energy
            ETrial=a*xTrial*xTrial+b*xTrial**4

            call random_number(urand,stat(id))
            if (exp(-(ETrial-E)/T).gt.urand(1)) then
               ! print*, 'E',E,'ETrial',ETrial,'T',T
                E=ETrial
                x=xTrial
            endif
            
        enddo
        
       
        call MPI_Send ( E, 1, MPI_DOUBLE_PRECISION, dest,   0, &
          MPI_COMM_WORLD,         error )
        call MPI_Send ( x, 1, MPI_DOUBLE_PRECISION, dest,   0, &
          MPI_COMM_WORLD,         error )
      end do

    end if

    deallocate ( Temp )
    deallocate ( Energy )
   

  return
end

