!   --------------------------------------------------------------
!
!    This module if defines two strucutres.  MCvar is for simulation
!    parameters of a fixed size.  MCdata is for allocateable variables.
!    These strucutres contain most of the variables needed for a MC
!    simulation of a wlc polymer with beads interacting via a field.
!
!        By Quinn MacPherson ~ Spring 2016
!
!   --------------------------------------------------------------
Module simMod
    IMPLICIT NONE
    INTEGER, Parameter :: NmoveTypes = 7 ! ******* YOU MAY NEED TO CHAGE THIS ***
  Type MCvar   ! Structure for simulation variables of known size
!     Simulation parameters
    INTEGER NT                ! Total number of beads  NT=NP*N
    INTEGER NB                ! Number of beads in a polymer NB=N*G
    INTEGER N                 ! Number of monomers in a polymer
    INTEGER G                 ! Beads per monomer
    INTEGER NP                ! Number of polymers
    DOUBLE PRECISION LBOX     ! Box length (approximate)
    DOUBLE PRECISION DEL      ! Discretization size (approximate)
    DOUBLE PRECISION L0       ! Equilibrium segment length
    DOUBLE PRECISION V        ! Bead volume
    DOUBLE PRECISION FA       ! Fraction of A beads
    DOUBLE PRECISION LAM      ! Chemical correlation parameter
    DOUBLE PRECISION EPS      ! Elasticity l0/(2lp)
    DOUBLE PRECISION CHI      ! Chi parameter value (solvent-polymer)        
    DOUBLE PRECISION KAP      ! Incompressibility parameter
    DOUBLE PRECISION Fpoly    ! Fraction Polymer
    DOUBLE PRECISION EU       ! Energy of binding for unmethalated
    DOUBLE PRECISION EM       ! Energy of binding for methalated
    DOUBLE PRECISION HP1_Bind ! Energy of binding of HP1 to eachother
    DOUBLE PRECISION F_METH   ! Fraction methalated is using option 2
    DOUBLE PRECISION LAM_METH ! eigenvalue of methalation setup
    DOUBLE PRECISION mu       ! chemical potential of HP1
    INTEGER NBIN     ! Number of bins
    INTEGER NBINX    ! Number of bin on an edge
    DOUBLE PRECISION PARA(10) ! Parameters for sswlc
        ! EB, EPAR, EPERP, GAM, ETA, ...


!   Monte Carlo Variables (for addaptation)
    INTEGER moveTypes
    DOUBLE PRECISION MCAMP(NmoveTypes) ! Amplitude of random change
    DOUBLE PRECISION MinAMP(NmoveTypes) ! minium amplitude
    DOUBLE PRECISION MaxAMP(NmoveTypes) ! maximum amplitude
    INTEGER MOVEON(NmoveTypes)         ! Is the move active
    INTEGER WINDOW(NmoveTypes)         ! Size of window for bead selection
    INTEGER MAXWINDOW(NmoveTypes)         ! Size of window for bead selection
    DOUBLE PRECISION WA_ratio(NmoveTypes)         ! Size of window for bead selection
    INTEGER SUCCESS(NmoveTypes)        ! Number of successes
    DOUBLE PRECISION moveSlope(NmoveTypes) ! target for ratio of window to anmplitude
    DOUBLE PRECISION PDesire(NmoveTypes) ! desired hit rate     
    DOUBLE PRECISION PHit(NmoveTypes) ! hit rate 
    INTEGER NADAPT(NmoveTypes)        ! Nunber of steps before addapt


!   Energys
    DOUBLE PRECISION ENERGY   ! Total energy
    DOUBLE PRECISION Eint     ! running Eint
    DOUBLE PRECISION EELAS(3) ! Elastic force
    DOUBLE PRECISION ECHI     ! CHI energy
    DOUBLE PRECISION EKAP     ! KAP energy
    DOUBLE PRECISION EBind    ! binding energy

!   Switches
    INTEGER confineType       ! type of Boundary Conditions
    INTEGER setType           ! initial condition type
!    INTEGER INTON             ! interaction on

!   Move Variables 
    DOUBLE PRECISION DEELAS(3)   ! Change in bending energy
    DOUBLE PRECISION DEINT    ! Change in self energy
    DOUBLE PRECISION DEBind   ! Change in binding energy
    DOUBLE PRECISION ECon     ! Confinement Energy
    INTEGER NPHI  ! NUMBER o phi values that change

!   Parallel Tempering variables
    Character*16 repSufix    ! prefix for writing files
    integer rep  ! which replica am I
    integer (kind = 4) error  ! MPI error
    double precision M ! M=\sum_i \sigma_i   like ising magnitization
    integer NPT  ! number of steps before parallel tempering
    integer PTON

  end Type

  Type MCData  ! for Allocateable variables
!   Configuration Data
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R   ! Conformation of polymer chains
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U   ! Conformation of polymer chains 
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: RP !Test Bead positions - only valid from IT1 to IT2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: UP !Test target vectors - only valid from IT1 to IT2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PHIA ! Volume fraction of A
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PHIB ! Volume fraction of B
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: Vol ! Volume fraction of A
    INTEGER, ALLOCATABLE, DIMENSION(:):: AB            ! Chemical identity of beads
    INTEGER, ALLOCATABLE, DIMENSION(:):: ABP           ! Test Chemical identity of beads
    INTEGER, ALLOCATABLE, DIMENSION(:):: METH          ! Methalation state of beads
    DOUBLE PRECISION, Allocatable, Dimension(:):: DPHIA    ! Change in phi A
    DOUBLE PRECISION, Allocatable, Dimension(:):: DPHIB    ! Change in phi A
    INTEGER, Allocatable, Dimension(:) :: INDPHI           ! Indices of the phi
  end TYPE
contains
Subroutine MCvar_printDiscription(mc)
    IMPLICIT NONE
    TYPE(MCvar) mc
    print*, "Total number of beads, NT=", mc%NT
    print*, "Number of beads in a polymer, NB=", mc%NB
    print*, "Number of monomers in a polymer, N=", mc%N
    print*, "Number of polymers, NP=",mc%NP
    print*, "Number of bins", mc%NBIN
    print*, "Number of bins in x direction", mc%NBINX
    print*, "LBOX=", mc%LBOX
    print*, "spatial escritation DEL=",mc%DEL
    print*, "L0=", mc%L0
    print*, "bead volume V=", mc%V
    print*, "elasticity EPS =", mc%EPS
    print*, "persistance length =",(mc%L0/(2.0*mc%EPS))
    print*, "solvent-polymer CHI =",mc%CHI
    print*, "incompressibility KAP=",mc%KAP
    print*, "volume fraction polymer =", mc%Fpoly
    print*, " -energy of binding unmethalated ", mc%EU," more positive for favorable binding"
    print*, " -energy of binding methalated",mc%EM
    print*, "HP1_Binding energy parameter", mc%HP1_Bind
    print*, "fraction Methalated", mc%F_METH
    print*, "LAM_METH", mc%LAM_METH
    print*, "chemical potential of HP1", mc%mu
    
end Subroutine
Subroutine MCvar_allocate(mc,md)
    IMPLICIT NONE
    TYPE(MCvar) mc
    TYPE(MCData) md
    INTEGER NT  ! total number of beads
    INTEGER NBIN ! total number of bins
    NT=mc%NT
    NBIN=mc%NBIN
    
     
    ALLOCATE(md%R(NT,3))
    ALLOCATE(md%U(NT,3))
    Allocate(md%RP(NT,3))
    Allocate(md%UP(NT,3))
    ALLOCATE(md%AB(NT))   !Chemical identity aka binding state
    ALLOCATE(md%ABP(NT))   !Chemical identity aka binding state
    ALLOCATE(md%METH(NT)) !Underlying methalation profile 
    ALLOCATE(md%PHIA(NBIN))
    ALLOCATE(md%PHIB(NBIN))
    ALLOCATE(md%DPHIA(NBIN))
    ALLOCATE(md%DPHIB(NBIN))
    ALLOCATE(md%Vol(NBIN))
    Allocate(md%INDPHI(NBIN))

end subroutine
Subroutine MCvar_defaultAmp(mc,NSTEP)
    IMPLICIT NONE
    DOUBLE PRECISION, Parameter :: PI=3.141592654
    TYPE(MCvar) mc
    INTEGER MCTYPE ! Type of move
    INTEGER NSTEP ! number of steps per save point
!   ~~~~~~~~~~~~~~~~~~~
!    Edit the following variables for better performance
!   ~~~~~~~~~~~~~~~~~~~
    !  Monte-Carlo simulation parameters
    mc%MCAMP(1)=0.5*PI
    mc%MCAMP(2)=0.3*mc%L0
    mc%MCAMP(3)=0.5*PI
    mc%MCAMP(4)=0.5*PI
    mc%MCAMP(5)=0.5*PI
    mc%MCAMP(6)=5.0*mc%L0
    mc%MCAMP(7)=1.0
    !switches to turn on various types of moves
    mc%MOVEON(1)=1  ! crank-shaft move
    mc%MOVEON(2)=1  ! slide move
    mc%MOVEON(3)=1  ! pivot move
    mc%MOVEON(4)=1  ! rotate move
    mc%MOVEON(5)=0  ! full chain rotation
    mc%MOVEON(6)=0  ! full chain slide
    mc%MOVEON(7)=1  ! Change in Binding state
    
    !     Initial segment window for MC moves
    mc%WINDOW(1)=15 ! used to be N*G
    mc%WINDOW(2)=15 ! used to be N*G
    mc%WINDOW(3)=15 ! used to be N*G
    mc%WINDOW(4)=15 ! used to be N*G
    mc%WINDOW(5)=15 ! used to be N*G
    mc%WINDOW(6)=15 ! used to be N*G
    mc%WINDOW(7)=15 ! used to be N*G

    !    Maximum window size (large windows are expensive)
    mc%MAXWINDOW(1)=150 
    mc%MAXWINDOW(2)=150 
    mc%MAXWINDOW(3)=150 
    mc%MAXWINDOW(4)=150 
    mc%MAXWINDOW(5)=150 
    mc%MAXWINDOW(6)=150 
    mc%MAXWINDOW(7)=150 
    Do MCTYPE=1,mc%moveTypes
        mc%moveSlope(MCTYPE)=dble(mc%WINDOW(MCTYPE))/mc%MCAMP(MCTYPE)
    enddo

    mc%MINAMP(1)=0.1*PI
    mc%MINAMP(2)=0.2*mc%L0
    mc%MINAMP(3)=0.2*PI
    mc%MINAMP(4)=0.2*PI
    mc%MINAMP(5)=0.2*PI
    mc%MINAMP(6)=0.2*mc%L0
    mc%MINAMP(7)=1.0   !Amplitude of move 7 is irrelivant

    mc%MAXAMP(1)=1.0*PI
    mc%MAXAMP(2)=1.0*mc%L0
    mc%MAXAMP(3)=1.0*PI
    mc%MAXAMP(4)=1.0*PI
    mc%MAXAMP(5)=1.0*PI
    mc%MAXAMP(6)=0.1*mc%LBOX
    mc%MAXAMP(7)=1 !Amplitude of move 7 is irrelivant
     
    DO MCTYPE=1,mc%moveTypes
        mc%NADAPT(MCTYPE)=1000 ! addapt after at most 1000 steps
        if (NSTEP.LE.mc%NADAPT(MCTYPE)) then 
            mc%NADAPT(MCTYPE)=NSTEP ! addapt at least every save point
        endif
          
        mc%PDESIRE(MCTYPE)=0.5 ! Target
        mc%SUCCESS(MCTYPE)=0
    ENDDO
end subroutine
Subroutine MCvar_addapt(mc,MCTYPE,ISTEP,rand_stat)
! Run this after say 1000 move in order to improve performance
    !use mt19937, only : grnd
    use mersenne_twister 
    IMPLICIT NONE
    INTEGER ISTEP    ! Step number, need to initialize if 1
    TYPE(MCvar) mc
    INTEGER MCTYPE   ! Type of move
    Double Precision floatWindow  !like window but a floating point
    real urand(1)
    type(random_stat) rand_stat  !for random number generator
    mc%PHIT(MCTYPE)=real(mc%SUCCESS(MCTYPE))/real(mc%NADAPT(MCTYPE))
    if (mc%PHIT(MCTYPE).GT.mc%PDESIRE(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MCAMP(MCTYPE)*1.05
    else
       mc%MCAMP(MCTYPE)=mc%MCAMP(MCTYPE)*0.95
    endif
    
    if (ISTEP.eq.1) then
        ! initialize
        mc%WA_ratio(MCTYPE)=mc%WINDOW(MCTYPE)/mc%MCAMP(MCTYPE)
    endif
    !drift to chosen ratio
    mc%WA_ratio(MCTYPE)=mc%WA_ratio(MCTYPE)*0.96+0.04*mc%moveSlope(MCTYPE)

    ! probabalisticaly round window
    floatWindow=mc%WA_ratio(MCTYPE)*mc%MCAMP(MCTYPE)
    call random_number(urand,rand_stat) 
    if (urand(1).lt.(floatWindow-floor(floatWindow))) then
        mc%Window(MCTYPE)=CEILING(mc%WA_ratio(MCTYPE)*mc%MCAMP(MCTYPE))
    else
        mc%Window(MCTYPE)=floor(mc%WA_ratio(MCTYPE)*mc%MCAMP(MCTYPE))
    endif

    ! amplitude limits
    if (mc%MCAMP(MCTYPE).GT.mc%MAXAMP(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MAXAMP(MCTYPE)
    elseif (mc%MCAMP(MCTYPE).LT.mc%MINAMP(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MINAMP(MCTYPE)
    endif
    
    !window limits
    if (mc%WINDOW(MCTYPE).LT.1) then
       mc%WINDOW(MCTYPE)=1
    elseif (mc%WINDOW(MCTYPE).GT.mc%MAXWINDOW(MCTYPE)) then
       mc%WINDOW(MCTYPE)=mc%MAXWINDOW(MCTYPE)
    endif

    mc%WA_ratio(MCTYPE)=dble(mc%WINDOW(MCTYPE))/mc%MCAMP(MCTYPE)   
    mc%SUCCESS(MCTYPE)=0
end subroutine
Subroutine MCvar_recenter(mc,md)
!  Prevents drift in periodic BC
    IMPLICIT NONE
    TYPE(MCvar) mc
    TYPE(MCData) md
    INTEGER IB, I, J   ! Couners
    DOUBLE PRECISION R0(3)  ! Offset to move by
    IB=1
    DO I=1,mc%NP
       R0(1)=nint(md%R(IB,1)/mc%LBOX-0.5)*mc%LBOX
       R0(2)=nint(md%R(IB,2)/mc%LBOX-0.5)*mc%LBOX
       R0(3)=nint(md%R(IB,3)/mc%LBOX-0.5)*mc%LBOX
       if (abs(R0(1)*R0(2)*R0(3)) .gt. 0.0001) then
           DO J=1,mc%N
              md%R(IB,1)=md%R(IB,1)-R0(1)
              md%R(IB,2)=md%R(IB,2)-R0(2)
              md%R(IB,3)=md%R(IB,3)-R0(3)
              IB=IB+1
           Enddo
           print*, "Error in MCvar_recenter"
           print*, "Shouldn't have to recenter"
           print*, "Remove this error if repeating BC"
           stop 1
      endif
    enddo
end Subroutine
Subroutine MCvar_printEnergies(this)
! For realtime feedback on MC simulation
    IMPLICIT NONE
    TYPE(MCvar) this
    print*, "Eint:", this%Eint
    print*, "Bending energy", this%EELAS(1)
    print*, "Par compression energy", this%EELAS(2)
    print*, "Shear energy", this%EELAS(3)
    print*, "ECHI", this%ECHI
    print*, "EKAP", this%EKAP
    print*, "EBind", this%EBind
end subroutine
Subroutine MCvar_printWindowStats(mc)
! For realtime feedback on addaptation
    IMPLICIT NONE
    TYPE(MCvar) mc
    INTEGER I ! counter
    I=0
    print*, "Succes | MCAMP | ratio | WINDOW| Type "
    Do I=1,mc%moveTypes
        if (mc%MOVEON(i).eq.1) then
            write(*,"(3f8.2,2I8)"), mc%phit(i), mc%MCAMP(i), mc%WA_ratio(i) ,  mc%WINDOW(i), i
        endif
    enddo
end subroutine
Subroutine MCvar_loadAB(mc,md,fileName)
! Loads AB for file...has not been tested
    IMPLICIT NONE
    TYPE(MCvar) mc
    TYPE(MCData) md
    character*16 fileName ! file name to load from
    INTEGER IB, I, J ! counters
    OPEN (UNIT = 1, FILE = fileName, STATUS = 'OLD')      
    IB=1
    DO I=1,mc%NP
       DO J=1,mc%NB
          READ(1,"(I2)") md%AB(IB)
          IB=IB+1
          enddo
    enddo 
    CLOSE(1)
end subroutine
Subroutine MCvar_saveR(mc,md,fileName,repeatingBC)
! Writes R and AB to file for analysis
! Rx  Ry  Rz AB
    IMPLICIT NONE
    INTEGER repeatingBC  ! 1 for reapeating boundary conditions
    INTEGER I,J,IB  ! counters
    TYPE(MCvar) mc
    TYPE(MCData) md
    character*16 fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    fullName=trim(fullName)
    print*,"fullName: ",fullName
    OPEN (UNIT = 1, FILE = fullName, STATUS = 'NEW')
    IB=1
    if (repeatingBC.eq.1) then 
        Do I=1,mc%NP
            Do J=1,mc%NB
                    WRITE(1,"(3f7.2,I2)") , &
                          md%R(IB,1)-0.*nint(md%R(IB,1)/mc%LBOX-0.5)*mc%LBOX, &
                          md%R(IB,2)-0.*nint(md%R(IB,2)/mc%LBOX-0.5)*mc%LBOX, &
                          md%R(IB,3)-0.*nint(md%R(IB,3)/mc%LBOX-0.5)*mc%LBOX, & 
                          md%AB(IB)
                IB=IB+1
            enddo
        enddo
        print*, "Error in MCvar_saveR"
        print*, "Are you sure you want reapeating BC"
        stop 1
    else
        Do I=1,mc%NP
            Do J=1,mc%NB
                 WRITE(1,"(3f8.3,I2)") md%R(IB,1),md%R(IB,2),md%R(IB,3),md%AB(IB)
                IB=IB+1
            enddo
        enddo
    endif
    Close(1)
end subroutine
Subroutine MCVar_savePHI(mc,md,fileName)
! Saves PHIA and PHIB to file for analysis
    IMPLICIT NONE
    INTEGER I  ! counters
    TYPE(MCvar) mc
    TYPE(MCData) md
    character*16 fileName 
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    OPEN (UNIT = 1, FILE = fullName, STATUS = 'NEW')
    Do I=1,mc%NBIN
        WRITE(1,"(2f7.2)") md%PHIA(I),md%PHIB(I)
    enddo
    Close(1)
end subroutine
Subroutine MCvar_saveU(mc,md,fileName)
! Saves U to ASCII file for analisys
    IMPLICIT NONE
    INTEGER I,J,IB  ! counters
    TYPE(MCvar) mc
    TYPE(MCData) md
    character*16 fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    OPEN (UNIT = 1, FILE = fullName, STATUS = 'NEW')
    IB=1
    Do I=1,mc%NP
        Do J=1,mc%NB
            WRITE(1,"(3f8.3,2I2)") md%U(IB,1),md%U(IB,2),md%U(IB,3),md%AB(IB),md%METH(IB)
            IB=IB+1
        enddo
    enddo
    Close(1)
end subroutine
Subroutine MCvar_saveParameters(mc,fileName)
! Write a number of parameters ASCII variables to file for reccords
    IMPLICIT NONE
    TYPE(MCvar) mc
    character*16 fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    OPEN (UNIT =1, FILE = fullName, STATUS = 'NEW')
        WRITE(1,"(I8)") mc%NT ! 1 Number of beads in simulation
        WRITE(1,"(I8)") mc%N  ! 2 Number of monomers in a polymer
        WRITE(1,"(I8)") mc%NB ! 3 Number of beads in a polymer
        WRITE(1,"(I8)") mc%NP ! 4 Number of polymers in simulation
        WRITE(1,"(I8)") mc%NT ! 5 Number of beads in simulation
        WRITE(1,"(I8)") mc%G  ! 6 Number of beads per monomer

        WRITE(1,"(f10.5)") mc%L0    ! Equilibrium segment length 
        WRITE(1,"(f10.5)") mc%CHI  ! 8  initail CHI parameter value 
        WRITE(1,"(f10.5)") mc%Fpoly ! Fraction polymer
        WRITE(1,"(f10.5)") mc%LBOX  ! 10 Lenth of box
        WRITE(1,"(f10.5)") mc%EU    ! Energy unmethalated       
        WRITE(1,"(f10.5)") mc%EM    ! 12 Energy methalated
        WRITE(1,"(f10.5)") mc%HP1_Bind ! Energy of HP1 binding
        WRITE(1,"(f10.5)") (mc%L0/mc%EPS) ! 14 Khun lenth
        WRITE(1,"(A)") "-999"  ! for historic reasons
        WRITE(1,"(f10.5)") mc%F_METH  ! methalation fraction
        WRITE(1,"(f10.5)") mc%LAM_METH  ! methalation lambda
    CLOSE(1)
end subroutine
Subroutine MCvar_appendEnergyData(mc,fileName,lnNum)
! Print Energy data
    IMPLICIT NONE
    INTEGER lnNum
    TYPE(MCvar) mc
    LOGICAL isfile
    character*16 fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif
    WRITE(1,"(I5, 9f9.1)") lnNum, &
           mc%EELAS(1), mc%EELAS(2), mc%EELAS(3), mc%Eint, &
           mc%EKap, mc%ECHI, mc%EBind, mc%M, mc%HP1_Bind
    Close(1)
end subroutine
Subroutine MCvar_appendAdaptData(mc,fileName,lnNum)
! Appends MC move addaptation data to the file  
    IMPLICIT NONE
    INTEGER lnNum
    TYPE(MCvar) mc
    LOGICAL isfile
    character*16 fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif
    WRITE(1,"(I4,21f8.2)") lnNum,& 
          REAL(mc%WINDOW(1)),mc%MCAMP(1),mc%PHIT(1), &
          REAL(mc%WINDOW(2)),mc%MCAMP(2),mc%PHIT(2), &
          REAL(mc%WINDOW(3)),mc%MCAMP(3),mc%PHIT(3), &
          REAL(mc%MOVEON(4)),mc%MCAMP(4),mc%PHIT(4), &
          REAL(mc%MOVEON(5)),mc%MCAMP(5),mc%PHIT(5), &
          REAL(mc%MOVEON(6)),mc%MCAMP(6),mc%PHIT(6), &
          REAL(mc%MOVEON(7)),mc%MCAMP(7),mc%PHIT(7)
    Close(1)
end subroutine
Subroutine MCvar_writeBindary(mc,md,baceName)
!    This function writes the contence of the structures mc and md
!  to a binary file.  If you add more variables to md you need to 
!  a seperate write command for them as it is not possible to write
!  a structure with allocatables to a binar file.
!    The contence are stored in 
!     baceName//'R'
!     baceName//'U'
!     etc.
    IMPLICIT NONE
    INTEGER sizeOfType         ! for binary saving
    TYPE(MCvar) mc             ! to be save or filled
    TYPE(MCData) md             ! to be save or filled
    CHARACTER(LEN=16) baceName ! for example 'record/'
    CHARACTER(LEN=16) fileName ! fileName
    CHARACTER(LEN=16) rw       ! either 'read' or 'write'
    CHARACTER(LEN=16) sufix    ! end of file name
    LOGICAL exists    ! Does file already exist?

    !  ------parameters -----

    sizeOfType=SIZEOF(mc)
    sufix='parameters'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) mc
    close(1)    

    ! -------- R --------

    sizeOfType=SIZEOF(md%R)
    sufix='R'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) md%R
    close(1)

    ! -------- U --------

    sizeOfType=SIZEOF(md%U)
    sufix='U'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) md%U
    close(1)

    ! -------- AB --------

    sizeOfType=SIZEOF(md%AB)
    sufix='AB'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) md%AB
    close(1)

    ! -------- Vol --------

    sizeOfType=SIZEOF(md%Vol)
    sufix='Vol'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) md%Vol
    close(1)
end Subroutine

Subroutine MCvar_readBindary(mc,md,baceName)
! This function reads what MCvar_writeBinary writes and 
! stores it to mc and md.  Be sure to allocate md before 
! calling this command.
    IMPLICIT NONE
    INTEGER sizeOfType         ! for binary saving
    TYPE(MCvar) mc             ! to be save or filled
    TYPE(MCData) md             ! to be save or filled
    CHARACTER(LEN=16) baceName ! for example 'record/'
    CHARACTER(LEN=16) fileName ! fileName
    CHARACTER(LEN=16) sufix    ! end of file name
    LOGICAL exists    ! Does file already exist?

    !  ------parameters -----

    sizeOfType=SIZEOF(mc)
    sufix='parameters'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) mc
    close(1)    

    ! -------- R --------

    sizeOfType=SIZEOF(md%R)
    sufix='R'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%R
    close(1)

    ! -------- U --------

    sizeOfType=SIZEOF(md%U)
    sufix='U'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%U
    close(1)

    ! -------- AB --------

    sizeOfType=SIZEOF(md%AB)
    sufix='AB'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%AB
    close(1)

    ! -------- Vol --------

    sizeOfType=SIZEOF(md%Vol)
    sufix='Vol'
    fileName=trim(baceName) // trim(sufix)
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        print*, 'Error in MCvar_readBinary. File ',fileName,'does not exist'
        stop 1
    endif
    read(1,rec=1) md%Vol
    close(1)
end Subroutine
end module simMod
