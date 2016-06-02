!---------------------------------------------------------------*
      
Subroutine wlcsim(rand_stat)
      
!     
!     This simulation tracks the dynamics of a single polymer
!     chain modeled as a discrete wormlike chain with bending
!     and stretching energy.
!     
!     Andrew Spakowitz
!     Written 9-2-13
!     
!     Edited by Shifan  prior to 2016
!     Edited heavily by Quinn in spring of 2016
!
!     Variables within the simulation

  use simMod
  use mersenne_twister  ! so that we know the size of rand_stat

  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654 ! Value of pi

  !Inputs
  type(random_stat) rand_stat ! state of random number generator
  integer ( kind = 4 ) id  ! thread id for PT

  ! miscellaneous
  INTEGER IOStatus          ! Status of file IO, may not be used
  INTEGER K                 ! Used in file IO

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
  INTEGER FRMCHEM           ! Initial chemical sequence
  INTEGER FRMMETH           ! Read methalation from file

! simulation data strucutres
  TYPE(MCvar) mc
  TYPE(MCData) md

!    historic useless things
 
  DOUBLE PRECISION trash

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     Load in the parameters for the simulation
! 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  open (unit=5, file='input/input')
  read (unit=5, fmt='(4(/))')
  read (unit=5, fmt=*) mc%PTON
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
  read (unit=5, fmt=*) mc%CHI
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) trash
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) mc%KAP
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) INDMAX
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FRMFILE
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) trash
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) INTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NCHI
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) trash
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NNOINT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NINIT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NSTEP
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) WRTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) trash
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) trash
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) trash
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) trash
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) trash
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
  mc%NPT=100
  mc%moveTypes=7
  mc%EElas(1)=0.0
  mc%EElas(2)=0.0
  mc%EElas(3)=0.0
  mc%Eint=0.0
  mc%EBind=0.0
  mc%EKap=0.0
  mc%ECHI=0.0
  mc%setType = 4 ! 4 for shereical
  mc%confineType = 3 ! 3 for sherical
  
  mc%NP=1; 
  mc%NT=mc%N*mc%NP*mc%G
  mc%NB=mc%N*mc%G
  if (mc%confineType.eq.3) then
      mc%LBOX=(mc%V*mc%NT*6/(mc%Fpoly*PI))**(1./3.)
  else
      mc%LBOX=(mc%V*mc%NT/mc%Fpoly)**(1./3.)
  endif
  mc%NBINX=nint(mc%LBOX/mc%DEL)
  mc%NBIN=mc%NBINX**3.
  mc%LBOX = mc%NBINX*mc%DEL! used to be: DEL=LBOX/NBINX
  call getpara(mc%PARA,mc%EPS,mc%L0,mc%LBOX)
  call MCvar_allocate(mc,md)

  print*, "Calculating Bin volumes" 
  mc%NBINX=nint(mc%LBOX/mc%DEL)
  if (mc%NBINX**3.ne.mc%NBIN) then
      print*, "error in MCsim. Wrong number of bins"
      stop 1
  endif
  call MC_caclVolume(mc%confineType,mc%NBINX,mc%DEL, mc%LBox, &
                     md%Vol,rand_stat)  ! calculate partial volumes

  print*, "setting up move types"
  call MCvar_defaultAmp(mc,NSTEP) ! set MC move aplitudes and windows

  INQUIRE (FILE = 'data/out1', exist = restart)
  if (.NOT.restart) then

    PRINT*, '-----new simulation-----'

!     Setup the initial condition
    print*, "setting initial position ..."
    call initcond(md%R,md%U,md%AB,mc%NT,mc%NB,mc%NP,FRMFILE,mc%PARA,mc%LBOX, &
                  mc%setType,rand_stat)

!     Load in AB sequence
    IF (FRMCHEM.EQ.1) THEN
        snapnm='input/ab'
        print *, "loading AB"
        call MCvar_loadAB(mc,md,snapnm)
    ELSE
        print*, "setting initial binding condition..."
        call initchem(md%AB,mc%NT,mc%N,mc%G,mc%NP,mc%FA,mc%LAM,rand_stat)
    ENDIF


    
!     Load methalation sequence
    IF (FRMMETH.EQ.1) THEN
        OPEN (UNIT = 2, FILE = 'input/meth', STATUS = 'OLD')
        ! more to come here ...
        CLOSE(2)
    ELSE
        print*, "setting initial chemical condition..."
        call initchem(md%METH,mc%NT,mc%N,mc%G,mc%NP,mc%F_METH,mc%LAM_METH,rand_stat)        
    ENDIF

    print*, snapnm
    if ( mc%PTON.eq.1) then
        call PT_override(mc,md)
    else
        mc%repSufix=''
    endif
    snapnm='data/r0'
    I=0;
    call MCvar_saveR(mc,md,snapnm,0)
   
    snapnm='data/params'
    call MCvar_saveParameters(mc,snapnm)

    snapnm='data/u0'
    call MCvar_saveU(mc,md,snapnm)

    INDEND = 0

 else

    PRINT*, '-----load simulation-----'
    snapnm='putBinaryFileNameHere'
    INDEND=0; INDEND=1/INDEND ! make this a variable of simmod before use
    call MCvar_readBindary(mc,md,snapnm)

 endif
    !call MCvar_printDiscription(mc)

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!              Begin simulation
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print*, 'Beginning simulation'
  IND=1

  DO WHILE ((IND+INDEND).LE.INDMAX)  ! INDEND is for sims. that are restarted

     call strength_schedule(IND+INDEND,mc%HP1_bind)
     if ((IND+INDEND).LE.NNOINT) then
         INTON=0
     else
         INTON=1
     endif
     
!    * Perform a MC simulation *
     call MCsim(mc,md,NSTEP,INTON,rand_stat)

!     Save the conformation and the metrics
     TENS=nint(log10(1.*INDEND+IND)-0.4999)+1
     write (fileind,'(I4)'), INDEND+IND

     !Save various energy contiributions to file 
     snapnm='data/out1'
     I=INDEND+IND
     call MCvar_appendEnergyData(mc,snapnm,I)
      

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
     
     PRINT*, '________________________________________'
     PRINT*, 'Time point ',IND+INDEND, ' out of', INDMAX
     call MCvar_printEnergies(mc)
     call MCvar_printWindowStats(mc)
         
     IND=IND+1
         
  ENDDO
  
END
Subroutine strength_schedule(ind,HP1_bind)
    Integer ind
    Double PRECISION HP1_bind
    Double precision maximum
    
    if(ind.lt.101) then
        HP1_bind=0
    elseif(ind.lt.111) then
        HP1_bind=maximum*0.1
    elseif(ind.lt.121) then
        HP1_bind=maximum*0.2
    elseif(ind.lt.131) then
        HP1_bind=maximum*0.3
    elseif(ind.lt.141) then
        HP1_bind=maximum*0.4
    elseif(ind.lt.151) then
        HP1_bind=maximum*0.5
    elseif(ind.lt.161) then
        HP1_bind=maximum*0.6
    elseif(ind.lt.171) then
        HP1_bind=maximum*0.7
    elseif(ind.lt.181) then
        HP1_bind=maximum*0.8
    elseif(ind.lt.191) then
        HP1_bind=maximum*0.9
    elseif(ind.lt.201) then
        HP1_bind=maximum*1.0
    elseif(ind.lt.211) then
        HP1_bind=maximum*1.3
    elseif(ind.lt.221) then
        HP1_bind=maximum*1.2
    elseif(ind.lt.231) then
        HP1_bind=maximum*1.1
    elseif(ind.lt.241) then
        HP1_bind=maximum*1.0
    elseif(ind.lt.251) then
        HP1_bind=maximum*0.9
    elseif(ind.lt.261) then
        HP1_bind=maximum*0.8
    elseif(ind.lt.271) then
        HP1_bind=maximum*0.7
    endif

end subroutine
!---------------------------------------------------------------*
