!---------------------------------------------------------------*
      
Subroutine wlcsim
      
!     
!     This simulation tracks the dynamics of a single polymer
!     chain modeled as a discrete wormlike chain with bending
!     and stretching energy.
!     
!     Andrew Spakowitz
!     Written 9-2-13
!     
!     Edited by Quinn in 2016
!
!     Variables within the simulation

  use mt19937, only : grnd, sgrnd, rnorm, mt, mti
  use simMod

  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654 ! Value of pi
  INTEGER IOStatus              ! Status of file IO, may not be used
  INTEGER K                     ! Used in file IO
  INTEGER PTON                  ! Parallel Tempering on
  INTEGER NRABOVE                 ! Total number of replicas  
  INTEGER NRBELOW                 ! Total number of replicas  
  INTEGER NRNOW                 ! Total number of replicas  
  INTEGER NRMAX                 ! Total number of replicas  

  INTEGER NR                ! Replica Number
  INTEGER I,J,IB            ! Index
  INTEGER x                 ! Working interger
  INTEGER INDMAX            ! Maximum index in series
  INTEGER IND               ! Ind in series
  INTEGER TENS              ! Decimal of index
  character*4 fileind       ! Index of output
  character*4 repind        ! Replica index for output
  character*16 snapnm       ! File for output
  INTEGER INDEND            ! Restart index
  logical restart           ! Restart from previous?
  INTEGER WRTON

!     Simulation input variables
  
  INTEGER FRMFILE           ! Initial condition
  INTEGER BROWN             ! Include Brownian forces
  INTEGER INTON             ! Include polymer interactions
  INTEGER NSTEP             ! Number of MC steps between save
  INTEGER NINIT             ! Number of initialization MC steps
  INTEGER NNOINT            ! Number of initialization MC steps without interactions
  INTEGER NCHI              ! Number of savepoints between chi change
  INTEGER NPT               ! Number of savepoints between tempering
  INTEGER FRMCHEM           ! Initial chemical sequence
  INTEGER FRMMETH           ! Read methalation from file

      
!     Structure analysis
      
  DOUBLE PRECISION RCOM(3)  ! Center of mass
  DOUBLE PRECISION DELR(3)  ! Mag of gyration tensor
  DOUBLE PRECISION RCOM0(3) ! Init val RCOM
  DOUBLE PRECISION DELR0(3) ! Init val DELR
  DOUBLE PRECISION DRCOM    ! Change in RCOM
  DOUBLE PRECISION SIG(3,3)
  DOUBLE PRECISION COR
  
  INTEGER  SON               !calculate Structure Factors
  INTEGER, PARAMETER:: XNUM = 11
  INTEGER, PARAMETER:: KNUM = (XNUM)**3
  INTEGER NVEC              !number of times calculating SVEC
  DOUBLE PRECISION KVEC(KNUM)
  DOUBLE PRECISION SVEC(KNUM)
      
!     Variables for the random number generators

  INTEGER IDUM              ! Seed for the generator
  DOUBLE PRECISION MOM(6)

!     Simulation parameters
  
  DOUBLE PRECISION DCHI       ! Chi parameter value
  DOUBLE PRECISION CHI0         ! Initial Chi before ramping
  DOUBLE PRECISION DELCHI       ! Chi ramping rate
  INTEGER PTID              ! ID to pair up replicas for PT
  INTEGER ACCBELOW

  Type(MCvar) mc
  Type(MCData) md
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     Load in the parameters for the simulation
! 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  open (unit=5, file='input/input')
  read (unit=5, fmt='(4(/))')
  read (unit=5, fmt=*) PTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%N
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%G
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%LBOX
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%DEL
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%L0
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%V
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%FA
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%LAM
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FRMCHEM
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%EPS
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) CHI0
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) DELCHI
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%KAP
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) INDMAX
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FRMFILE
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) BROWN
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) INTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NCHI
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NPT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NNOINT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NINIT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NSTEP
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) WRTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) PTID
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRABOVE
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRBELOW
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRNOW
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRMAX
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%EU 
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%EM
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%HP1_Bind
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%F_METH
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%LAM_METH
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%Fpoly
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FRMMETH
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%mu

  close(5)
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!      choose simulation variables
! 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mc%EElas=0
  mc%CHI=CHI0

  mc%setType = 4 ! 4 for shereical
  mc%confineType = 3 ! 3 for sherical
  
  print*, 'EU', mc%EU, 'mu',mc%mu
  mc%NP=1; 
  mc%NT=mc%N*mc%NP*mc%G
  if (mc%confineType.eq.3) then
      mc%LBOX=(mc%V*mc%NT*6/(mc%Fpoly*PI))**(1./3.)
  else
      mc%LBOX=(mc%V*mc%NT/mc%Fpoly)**(1./3.)
  endif
  mc%NBINX=nint(mc%LBOX/mc%DEL)
  mc%NBIN=mc%NBINX**3.
  mc%LBOX = mc%NBINX*mc%DEL! used to be: DEL=LBOX/NBINX
  call getpara(mc%PARA,mc%EPS,mc%L0,mc%LBOX)
  mc%moveTypes=7
  print*,"Space descritation lenth, DEL=", mc%DEL
  print*,"Box lenth, LBOX=",mc%LBOX
  print*,"Monomer volume, V=",mc%V
  print*, "NT",mc%NT," NBIN",mc%NBIN," moveTypes",mc%moveTypes
  call MCvar_allocate(mc,md)

  print*, "Calculating Bin volumes" 
  mc%NBINX=nint(mc%LBOX/mc%DEL)
  if (mc%NBINX**3.ne.mc%NBIN) then
      print*, "error in MCsim. Wrong number of bins"
  endif
  call MC_caclVolume(mc%confineType,mc%NBINX,mc%DEL, mc%LBox, &
                     md%Vol)  ! calculate partial volumes

  print*, "setting up move types"
  call MCvar_defaultAmp(mc,NSTEP) ! set MC move aplitudes and windows

  INQUIRE (FILE = 'data/out1', exist = restart)
  if (.NOT.restart) then

    PRINT*, '-----new simulation-----'

!     Setup the initial condition
    print*, "setting initial position ..."
    mc%NB=mc%N*mc%G
    call initcond(md%R,md%U,md%AB,mc%NT,mc%NB,mc%NP,IDUM,FRMFILE,mc%PARA,mc%LBOX,mc%setType)

!     Load in AB sequence
    print *, "loading AB"
    IF (FRMCHEM.EQ.1) THEN
        snapnm='input/ab'
       call MCvar_loadAB(mc,md,snapnm)
    ELSE
       print*, "setting initial binding condition..."
       call initchem(md%AB,mc%NT,mc%N,mc%G,mc%NP,mc%FA,mc%LAM)
    ENDIF


    
!     Load methalation sequence
    IF (FRMMETH.EQ.1) THEN
        OPEN (UNIT = 2, FILE = 'input/meth', STATUS = 'OLD')
        ! more to come here ...
        CLOSE(2)
    ELSE
        print*, "setting initial chemical condition..."
        call initchem(md%METH,mc%NT,mc%N,mc%G,mc%NP,mc%F_METH,mc%LAM_METH)        
   ENDIF


    SON=0 ! calculate structure factor

    snapnm='data/r0'
    I=0;
    call MCvar_saveR(mc,md,snapnm,0)
   
    snapnm='data/params'
    call MCvar_saveParameters(mc,snapnm)

    snapnm='data/u0'
    call MCvar_saveU(mc,md,snapnm)

!     Open the output files
    INDEND = 0
    OPEN (UNIT = 1, FILE = 'data/out1', STATUS = 'NEW')
    OPEN (UNIT = 2, FILE = 'data/out2', STATUS = 'NEW')

 else

    PRINT*, '-----load simulation-----'
    mc%NB=mc%N*mc%G
    call load_from_old(md%R,md%U,md%AB,mc%CHI,DCHI,mc%NT,mc%NB,mc%NP,IDUM,FRMFILE,mc%PARA,mc%LBOX,INDEND)

 endif

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!              Begin simulation
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print*, 'Beginning simulation'
  IND=1

  !part 7 - PT
  OPEN(unit=1,file='data/ptnow',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  OPEN(unit=1,file='data/calcpnow',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  OPEN(unit=1,file='data/swapnow',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  OPEN(unit=1,file='data/swapend',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  mc%Eint=0.0
  DO WHILE ((IND+INDEND).LE.INDMAX)  ! INDEND is for sims. that are restarted
!     Parallel tempering + chi annealing
!     call chisched(IND,NR,CHI,DCHI,CHI0,DELCHI,NCHI,NPT,restart)

     if ((IND+INDEND).LE.NNOINT) then
         INTON=0
     else
         INTON=1
     endif
     
!     Perform a MC simulation
     SON=0 ! calculate structure factor
     call MCsim(mc,md,NSTEP,BROWN,INTON,KVEC,SVEC,SON, &
          PTON,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,NPT,PTID)
!     Save the conformation and the metrics
     TENS=nint(log10(1.*INDEND+IND)-0.4999)+1
     write (fileind,'(I4)'), INDEND+IND

     
     OPEN (UNIT = 1, FILE = 'data/out1', STATUS = 'OLD')
     DO K = 1,IND+INDEND-1
        READ(1,*)
     ENDDO
     WRITE(1,"(I4,5f20.3)") INDEND+IND,mc%EELAS(1),mc%EELAS(2),mc%EELAS(3),mc%Eint,mc%EKAP
     CLOSE(1)
      
     !part 2 - CHI
     OPEN (UNIT = 2, FILE = 'data/out2', STATUS = 'OLD')
     DO K = 1,IND+INDEND-1
        READ(2,*)
     ENDDO
     !  WRITE(2,"(I15,1f10.2)") INDPT,CHI
     WRITE(2,"(I4,1f20.8)") INDEND+IND,mc%CHI*mc%G
     CLOSE(2)

     !part 2.5 - adaptations
     snapnm='data/out3'
     I=INDEND+IND
     call MCvar_appendAdaptData(mc,snapnm,I)
     

     IF (WRTON.EQ.1) THEN
        !part 3 - R
        snapnm= 'data/r'//fileind((4-TENS+1):4)
        I=0
        call MCvar_saveR(mc,md,snapnm,I)
     ENDIF

     !part 6 - S(k)
     IF (SON.EQ.1) THEN
        snapnm= 'data/s'//fileind((4-TENS+1):4)
        OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
        DO K=1,KNUM
           WRITE(1,"(2f20.3)") KVEC(K),SVEC(K)
        ENDDO
        CLOSE(1)
     ENDIF

     
     PRINT*, '________________________________________'
     PRINT*, 'Time point ',IND+INDEND, ' out of', INDMAX
     call MCvar_printEnergies(mc)
     call MCvar_printWindowStats(mc)
         
     IND=IND+1
         
  ENDDO
  
END
      
!---------------------------------------------------------------*
