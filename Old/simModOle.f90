Module simMod

  IMPLICIT NONE

  Type MCvar

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

!   Configuration Data
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R	 ! Conformation of polymer chains
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U	 ! Conformation of polymer chains
  
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PHIA ! Volume fraction of A
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PHIB ! Volume fraction of B
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: Vol ! Volume fraction of A
    INTEGER, ALLOCATABLE, DIMENSION(:):: AB            ! Chemical identity of beads
    INTEGER, ALLOCATABLE, DIMENSION(:):: ABP           ! Test Chemical identity of beads
    INTEGER, ALLOCATABLE, DIMENSION(:):: METH          ! Methalation state of beads

!   Monte Carlo Variables (for addaptation)
    Integer moveTypes ! Number of different MC move types
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: MCAMP ! Amplitude of random change
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: MinAMP ! minium amplitude
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: MaxAMP ! maximum amplitude
    INTEGER, ALLOCATABLE, DIMENSION(:):: MOVEON         ! Is the move active
    INTEGER, ALLOCATABLE, DIMENSION(:):: WINDOW         ! Size of window for bead selection
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: WA_ratio         ! Size of window for bead selection
    INTEGER, ALLOCATABLE, DIMENSION(:):: SUCCESS        ! Number of successes
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: moveSlope ! target for ratio of window to anmplitude
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PDesire ! desired hit rate 
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PHit ! hit rate 
    INTEGER, ALLOCATABLE, DIMENSION(:):: NADAPT        ! Nunber of steps before addapt


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
    INTEGER INTON             ! interaction on

!   Move Variables 
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: RP !Test Bead positions - only valid from IT1 to IT2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: UP !Test target vectors - only valid from IT1 to IT2
    DOUBLE PRECISION DEELAS(3)   ! Change in bending energy
    DOUBLE PRECISION DEINT    ! Change in self energy
    DOUBLE PRECISION DEBind   ! Change in binding energy
    DOUBLE PRECISION ECon     ! Confinement Energy
    DOUBLE PRECISION, Allocatable, Dimension(:):: DPHIA    ! Change in phi A
    DOUBLE PRECISION, Allocatable, Dimension(:):: DPHIB    ! Change in phi A
    INTEGER, Allocatable, Dimension(:) :: INDPHI           ! Indices of the phi
    INTEGER NPHI  ! NUMBER o phi values that change

  end Type
contains
Subroutine MCvar_allocate(this)
  IMPLICIT NONE
  TYPE(MCvar) this
  INTEGER NT  ! total number of beads
  INTEGER NBIN ! total number of bins
  INTEGER moveTypes ! number of MoveTypes
  NT=this%NT
  NBIN=this%NBIN
  moveTypes=this%moveTypes
  

  ALLOCATE(this%R(NT,3))
  ALLOCATE(this%U(NT,3))
  ALLOCATE(this%AB(NT))   !Chemical identity aka binding state
  ALLOCATE(this%ABP(NT))   !Chemical identity aka binding state
  ALLOCATE(this%METH(NT)) !Underlying methalation profile 
  ALLOCATE(this%PHIA(NBIN))
  ALLOCATE(this%PHIB(NBIN))
  ALLOCATE(this%DPHIA(NBIN))
  ALLOCATE(this%DPHIB(NBIN))
  ALLOCATE(this%Vol(NBIN))
  Allocate(this%INDPHI(NBIN))

  Allocate(this%MCAMP(moveTypes)) 
  Allocate(this%MinAMP(moveTypes)) 
  Allocate(this%MaxAMP(moveTypes)) 
  Allocate(this%MOVEON(moveTypes)) 
  Allocate(this%WINDOW(moveTypes)) 
  Allocate(this%WA_ratio(moveTypes)) 
  Allocate(this%SUCCESS(moveTypes)) 
  Allocate(this%PHIT(moveTypes)) 
  Allocate(this%moveSlope(moveTypes)) 
  Allocate(this%PDesire(moveTypes))
  Allocate(this%NADAPT(moveTypes))

  Allocate(this%RP(NT,3))
  Allocate(this%UP(NT,3))

  this%MCAMP(1)=0.5
!    this%NT = 100
!    this%NP = 1
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
      mc%MOVEON(3)=0  ! pivot move
      mc%MOVEON(4)=1  ! rotate move
      mc%MOVEON(5)=0  ! full chain rotation
      mc%MOVEON(6)=0  ! full chain slide
      mc%MOVEON(7)=1 ! Change in Binding state
    
      !     Initial segment window for MC moves
      mc%WINDOW(1)=15 ! used to be N*G
      mc%WINDOW(2)=15 ! used to be N*G
      mc%WINDOW(3)=15 ! used to be N*G
      mc%WINDOW(4)=15 ! used to be N*G
      mc%WINDOW(5)=15 ! used to be N*G
      mc%WINDOW(6)=15 ! used to be N*G
      mc%WINDOW(7)=15 ! used to be N*G
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
Subroutine MCvar_addapt(mc,MCTYPE,ISTEP)
    IMPLICIT NONE
    INTEGER ISTEP    ! Step number, need to initialize if 1
    TYPE(MCvar) mc   
    INTEGER MCTYPE   ! Type of move
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
    mc%WINDOW(MCTYPE)=CEILING(mc%WA_ratio(MCTYPE)*mc%MCAMP(MCTYPE)) 


    if (mc%MCAMP(MCTYPE).GT.mc%MAXAMP(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MAXAMP(MCTYPE)
       mc%WA_ratio(MCTYPE)=dble(mc%WINDOW(MCTYPE))/mc%MCAMP(MCTYPE)
       !WINDOW(MCTYPE)=WINDOW(MCTYPE)+1
    endif

    if (mc%MCAMP(MCTYPE).LT.mc%MINAMP(MCTYPE)) then
       mc%MCAMP(MCTYPE)=mc%MINAMP(MCTYPE)
       mc%WA_ratio(MCTYPE)=dble(mc%WINDOW(MCTYPE))/mc%MCAMP(MCTYPE)
       !WINDOW(MCTYPE)=WINDOW(MCTYPE)-1                
       !if ((MCTYPE.EQ.4).OR.(MCTYPE.EQ.5).OR.(MCTYPE.EQ.6)) then
       !   MOVEON(MCTYPE)=0
       !endif
    endif
    
    !window limits
    if (mc%WINDOW(MCTYPE).LT.1) then
       mc%WINDOW(MCTYPE)=1
       mc%WA_ratio(MCTYPE)=dble(mc%WINDOW(MCTYPE))/mc%MCAMP(MCTYPE)
    endif

    if (mc%WINDOW(MCTYPE).GT.mc%NB) then
       mc%WINDOW(MCTYPE)=mc%NB
       mc%WA_ratio(MCTYPE)=dble(mc%WINDOW(MCTYPE))/mc%MCAMP(MCTYPE)
    endif
   
    mc%SUCCESS(MCTYPE)=0
end subroutine
Subroutine MCvar_recenter(mc)
    IMPLICIT NONE
    TYPE(MCvar) mc
    INTEGER IB, I, J
    DOUBLE PRECISION R0(3)
    IB=1
    DO I=1,mc%NP
       R0(1)=nint(mc%R(IB,1)/mc%LBOX-0.5)*mc%LBOX
       R0(2)=nint(mc%R(IB,2)/mc%LBOX-0.5)*mc%LBOX
       R0(3)=nint(mc%R(IB,3)/mc%LBOX-0.5)*mc%LBOX
       if (abs(R0(1)*R0(2)*R0(3)) .gt. 0.0001) then
           DO J=1,mc%N
              mc%R(IB,1)=mc%R(IB,1)-R0(1)
              mc%R(IB,2)=mc%R(IB,2)-R0(2)
              mc%R(IB,3)=mc%R(IB,3)-R0(3)
              IB=IB+1
           Enddo
      endif
    enddo
end Subroutine
Subroutine MCvar_printEnergies(this)
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
    IMPLICIT NONE
    TYPE(MCvar) mc
    INTEGER I
    I=0
    print*, "Succes | MCAMP | ratio | WINDOW| Type "
    Do I=1,mc%moveTypes
        if (mc%MOVEON(i).eq.1) then
            write(*,"(3f8.2,2I8)"), mc%phit(i), mc%MCAMP(i), mc%WA_ratio(i) ,  mc%WINDOW(i), i
        endif
    enddo
end subroutine
Subroutine MCvar_loadAB(this,snapnm)
    IMPLICIT NONE
    TYPE(MCvar) this
    character*16 snapnm
    INTEGER IB, I, J
       OPEN (UNIT = 1, FILE = snapnm, STATUS = 'OLD')      
       IB=1
       DO 1 I=1,this%NP
          DO 2 J=1,this%NB
             READ(1,"(I2)") this%AB(IB)
             IB=IB+1
2         CONTINUE
1      CONTINUE 
       CLOSE(1)
end subroutine
Subroutine MCvar_saveR(this,snapnm,repeatingBC)
    IMPLICIT NONE
    INTEGER repeatingBC
    INTEGER I,J,IB  ! counters
    TYPE(MCvar) this
    character*16 snapnm
    OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
    IB=1
    if (repeatingBC.eq.1) then
        Do I=1,this%NP
            Do J=1,this%NB
                    WRITE(1,"(3f7.2,I2)") , &
                          this%R(IB,1)-0.*nint(this%R(IB,1)/this%LBOX-0.5)*this%LBOX, &
                          this%R(IB,2)-0.*nint(this%R(IB,2)/this%LBOX-0.5)*this%LBOX, &
                          this%R(IB,3)-0.*nint(this%R(IB,3)/this%LBOX-0.5)*this%LBOX, & 
                          this%AB(IB)
                IB=IB+1
            enddo
        enddo
    else
        Do I=1,this%NP
            Do J=1,this%NB
                 WRITE(1,"(3f8.3,I2)") this%R(IB,1),this%R(IB,2),this%R(IB,3),this%AB(IB)
                IB=IB+1
            enddo
        enddo
    endif
    Close(1)
end subroutine
Subroutine MCVar_savePHI(this,snapnm)
    IMPLICIT NONE
    INTEGER I  ! counters
    TYPE(MCvar) this
    character*16 snapnm 
    OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
    Do I=1,this%NBIN
        WRITE(1,"(2f7.2)") this%PHIA(I),this%PHIB(I)
    enddo
    Close(1)
end subroutine
Subroutine MCvar_saveU(this,snapnm)
    IMPLICIT NONE
    INTEGER I,J,IB  ! counters
    TYPE(MCvar) this
    character*16 snapnm
    OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
    IB=1
    Do I=1,this%NP
        Do J=1,this%NB
            WRITE(1,"(3f8.3,2I2)") this%U(IB,1),this%U(IB,2),this%U(IB,3),this%AB(IB),this%METH(IB)
            IB=IB+1
        enddo
    enddo
    Close(1)
end subroutine
Subroutine MCvar_saveParameters(this,snapnm)
    IMPLICIT NONE
    TYPE(MCvar) this
    character*16 snapnm
    OPEN (UNIT =1, FILE = snapnm, STATUS = 'NEW')
        WRITE(1,"(I8)") this%NT ! 1 Number of beads in simulation
        WRITE(1,"(I8)") this%N  ! 2 Number of monomers in a polymer
        WRITE(1,"(I8)") this%NB ! 3 Number of beads in a polymer
        WRITE(1,"(I8)") this%NP ! 4 Number of polymers in simulation
        WRITE(1,"(I8)") this%NT ! 5 Number of beads in simulation
        WRITE(1,"(I8)") this%G  ! 6 Number of beads per monomer

        WRITE(1,"(f10.5)") this%L0    ! Equilibrium segment length 
        WRITE(1,"(f10.5)") this%CHI  ! 8  initail CHI parameter value 
        WRITE(1,"(f10.5)") this%Fpoly ! Fraction polymer
        WRITE(1,"(f10.5)") this%LBOX  ! 10 Lenth of box
        WRITE(1,"(f10.5)") this%EU    ! Energy unmethalated       
        WRITE(1,"(f10.5)") this%EM    ! 12 Energy methalated
        WRITE(1,"(f10.5)") this%HP1_Bind ! Energy of HP1 binding
        WRITE(1,"(f10.5)") (this%L0/this%EPS) ! 14 Khun lenth
        WRITE(1,"(A)") "-999"  ! for historic reasons
        WRITE(1,"(f10.5)") this%F_METH  ! methalation fraction
        WRITE(1,"(f10.5)") this%LAM_METH  ! methalation lambda
    CLOSE(1)
end subroutine
Subroutine MCvar_appendAdaptData(this,snapnm,lnNum)
    IMPLICIT NONE
    INTEGER lnNum
    TYPE(MCvar) this
    LOGICAL isfile
    character*16 snapnm
    inquire(file = snapnm, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = snapnm, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = snapnm, STATUS = 'new')
    endif
    WRITE(1,"(I4,21f8.2)") lnNum,& 
          REAL(this%WINDOW(1)),this%MCAMP(1),this%PHIT(1), &
          REAL(this%WINDOW(2)),this%MCAMP(2),this%PHIT(2), &
          REAL(this%WINDOW(3)),this%MCAMP(3),this%PHIT(3), &
          REAL(this%MOVEON(4)),this%MCAMP(4),this%PHIT(4), &
          REAL(this%MOVEON(5)),this%MCAMP(5),this%PHIT(5), &
          REAL(this%MOVEON(6)),this%MCAMP(6),this%PHIT(6), &
          REAL(this%MOVEON(7)),this%MCAMP(7),this%PHIT(7)
    Close(1)
end subroutine
Subroutine MCvar_rwBindary(rw,mc,baceName)
IMPLICIT NONE
INTEGER, Parameter :: moveTypes = 7 ! ******* YOU MAY NEED TO CHAGE THIS ***
INTEGER sizeOfType         ! for binary saving
TYPE(MCvar) mc             ! to be save or filled
CHARACTER(LEN=16) baceName ! for example 'record/'
CHARACTER(LEN=16) fileName ! fileName
CHARACTER(LEN=16) rw       ! either 'read' or 'write'
CHARACTER(LEN=16) sufix    ! end of file name
LOGICAL exists    ! Does file already exist?


TYPE fixedMCvar
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
    Integer moveTypes ! Number of different MC move types
    DOUBLE PRECISION MCAMP(moveTypes) ! Amplitude of random change
    DOUBLE PRECISION MinAMP(moveTypes) ! minium amplitude
    DOUBLE PRECISION MaxAMP(moveTypes) ! maximum amplitude
    INTEGER MOVEON(moveTypes)         ! Is the move active
    INTEGER WINDOW(moveTypes)         ! Size of window for bead selection
    DOUBLE PRECISION WA_ratio(moveTypes)         ! Size of window for bead selection
    INTEGER SUCCESS(moveTypes)        ! Number of successes
    DOUBLE PRECISION moveSlope(moveTypes) ! target for ratio of window to anmplitude
    DOUBLE PRECISION PDesire(moveTypes) ! desired hit rate     
    DOUBLE PRECISION PHit(moveTypes) ! hit rate 
    INTEGER NADAPT(moveTypes)        ! Nunber of steps before addapt

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
    INTEGER INTON             ! interaction on

!   Move Variables 
    DOUBLE PRECISION DEELAS(3)   ! Change in bending energy
    DOUBLE PRECISION DEINT    ! Change in self energy
    DOUBLE PRECISION DEBind   ! Change in binding energy
    DOUBLE PRECISION ECon     ! Confinement Energy
    INTEGER NPHI  ! NUMBER o phi values that change

end type
type(fixedMCvar) f

if (moveTypes.ne.mc%moveTypes) then
    print*, "change moveTypes in MCvar_wrBinary"
    stop 1
endif


if (rw.eq."write") then

    f%NT           =mc%NT         
    f%NB           =mc%NB      
    f%N            =mc%N       
    f%G            =mc%G       
    f%NP           =mc%NP      
    f%LBOX         =mc%LBOX    
    f%DEL          =mc%DEL     
    f%L0           =mc%L0      
    f%V            =mc%V       
    f%FA           =mc%FA      
    f%LAM          =mc%LAM     
    f%EPS          =mc%EPS     
    f%CHI          =mc%CHI     
    f%KAP          =mc%KAP     
    f%Fpoly        =mc%Fpoly   
    f%EU           =mc%EU      
    f%EM           =mc%EM      
    f%HP1_Bind     =mc%HP1_Bind
    f%F_METH       =mc%F_METH  
    f%LAM_METH     =mc%LAM_METH
    f%mu           =mc%mu      
    f%NBIN         =mc%NBIN    
    f%NBINX        =mc%NBINX   
    f%PARA         =mc%PARA    
    f%moveTypes    =mc%moveTypes
    f%MCAMP        =mc%MCAMP    
    f%MinAMP       =mc%MinAMP   
    f%MaxAMP       =mc%MaxAMP   
    f%MOVEON       =mc%MOVEON   
    f%WINDOW       =mc%WINDOW   
    f%WA_ratio     =mc%WA_ratio 
    f%SUCCESS      =mc%SUCCESS  
    f%moveSlope    =mc%moveSlope
    f%PDesire      =mc%PDesire  
    f%PHit         =mc%PHit     
    f%NADAPT       =mc%NADAPT   
    f%ENERGY       =mc%ENERGY   
    f%Eint         =mc%Eint     
    f%EELAS        =mc%EELAS 
    f%ECHI         =mc%ECHI     
    f%EKAP         =mc%EKAP     
    f%EBind        =mc%EBind       
    f%confineType  =mc%confineType      
    f%setType      =mc%setType          
    f%INTON        =mc%INTON           
    f%DEELAS       =mc%DEELAS
    f%DEINT        =mc%DEINT    
    f%DEBind       =mc%DEBind   
    f%ECon         =mc%ECon     
    f%NPHI         =mc%NPHI        

    !  ------parameters -----

    sizeOfType=SIZEOF(f)
    sufix='parameters'
    fileName=trim(baceName) // trim(sufix)
    print*, "fileName: ", fileName
    inquire(file=fileName,exist=exists)
    if(exists) then
        open(unit=1,file=fileName, status='old', &
             form='unformatted',access='direct',recl=sizeOfType)
    else
        open(unit=1,file=fileName, status='new', &
             form='unformatted',access='direct',recl=sizeOfType)
    endif
    write(1,rec=1) f
    close(1)    

    ! -------- R --------

    sizeOfType=SIZEOF(R)
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
    write(1,rec=1) R
    close(1)

    ! -------- R --------

    sizeOfType=SIZEOF(R)
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
    write(1,rec=1) R

    ! -------- R --------

    sizeOfType=SIZEOF(R)
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
    write(1,rec=1) R


elseif(rw.eq.'read') then

    sizeOfType=SIZEOF(f)
    sufix='parameters'
    fileName=trim(baceName) // trim(sufix)
    open(unit=1,file=fileName, &
         form='unformatted',access='direct',recl=sizeOfType)
    read(1,rec=1) f
    close(1)

    mc%NT           =f%NT         
    mc%NB           =f%NB      
    mc%N            =f%N       
    mc%G            =f%G       
    mc%NP           =f%NP      
    mc%LBOX         =f%LBOX    
    mc%DEL          =f%DEL     
    mc%L0           =f%L0      
    mc%V            =f%V       
    mc%FA           =f%FA      
    mc%LAM          =f%LAM     
    mc%EPS          =f%EPS     
    mc%CHI          =f%CHI     
    mc%KAP          =f%KAP     
    mc%Fpoly        =f%Fpoly   
    mc%EU           =f%EU      
    mc%EM           =f%EM      
    mc%HP1_Bind     =f%HP1_Bind
    mc%F_METH       =f%F_METH  
    mc%LAM_METH     =f%LAM_METH
    mc%mu           =f%mu      
    mc%NBIN         =f%NBIN    
    mc%NBINX        =f%NBINX   
    mc%PARA         =f%PARA    
    mc%moveTypes    =f%moveTypes
    mc%MCAMP        =f%MCAMP    
    mc%MinAMP       =f%MinAMP   
    mc%MaxAMP       =f%MaxAMP   
    mc%MOVEON       =f%MOVEON   
    mc%WINDOW       =f%WINDOW   
    mc%WA_ratio     =f%WA_ratio 
    mc%SUCCESS      =f%SUCCESS  
    mc%moveSlope    =f%moveSlope
    mc%PDesire      =f%PDesire  
    mc%PHit         =f%PHit     
    mc%NADAPT       =f%NADAPT   
    mc%ENERGY       =f%ENERGY   
    mc%Eint         =f%Eint     
    mc%EELAS        =f%EELAS 
    mc%ECHI         =f%ECHI     
    mc%EKAP         =f%EKAP     
    mc%EBind        =f%EBind       
    mc%confineType  =f%confineType      
    mc%setType      =f%setType          
    mc%INTON        =f%INTON           
    mc%DEELAS       =f%DEELAS
    mc%DEINT        =f%DEINT    
    mc%DEBind       =f%DEBind   
    mc%ECon         =f%ECon     
    mc%NPHI         =f%NPHI        

else
    print*, "unkonw option. Please use 'read' or 'write'"
    stop 1
endif
end Subroutine

end module simMod
