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

  use setPrecision
  use simMod
  use mersenne_twister  ! so that we know the size of rand_stat

  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654_dp ! Value of pi

  !Inputs
  type(random_stat) rand_stat ! state of random number generator

  ! miscellaneous
  INTEGER I            
  character*4 fileind       ! Index of output
  character*16 iostr       ! File for output
  logical restart           ! Restart from previous?

!     Simulation input variables
  
  INTEGER INTON             ! Include polymer interactions

! simulation data strucutres
  TYPE(MCvar) mc
  TYPE(MCData) md

  iostr='input/params'
  print*, "setting parameters from: ", iostr
  call MCvar_setParams(mc,iostr)
  call MCvar_allocate(mc,md)
 

  INQUIRE (FILE = 'data/out1', exist = restart)
  if (.NOT.restart) then

    PRINT*, '-----new simulation-----'
!    Calculate volume of bins
    if (mc%confineType.eq.3) then 
        print*, "Calculating Bin volumes"
        call MC_caclVolume(mc%confineType,mc%NBINX,mc%DEL, mc%LBox(1), &
                           md%Vol,rand_stat)  ! calculate partial volumes
        print*, "Done Calculating Bin volumes"
    else
        do I=1,mc%NBIN
             md%Vol(I)=mc%del**3
        enddo
    endif
    
!     Setup the initial condition
    print*, "setting initial position ..."
    call initcond(md%R,md%U,md%AB,mc%NT,mc%NB,mc%NP,mc%FRMFILE,mc%PARA,mc%LBOX, &
                  mc%setType,rand_stat)

!     Load in AB sequence
    IF (mc%FRMCHEM) THEN
        iostr='input/ab'
        print *, "loading AB"
        call MCvar_loadAB(mc,md,iostr)
    ELSE
        print*, "setting initial binding condition..."
        call initchem(md%AB,mc%NT,mc%N,mc%G,mc%NP,mc%FA,mc%LAM,rand_stat)
    ENDIF

    
!     Load methalation sequence
    IF (mc%FRMMETH) THEN
        OPEN (UNIT = 2, FILE = 'input/meth', STATUS = 'OLD')
        ! more to come here ...
        CLOSE(2)
    ELSE
        print*, "setting initial chemical condition..."
        call initchem(md%METH,mc%NT,mc%N,mc%G,mc%NP,mc%F_METH,mc%LAM_METH,rand_stat)        
    ENDIF

!      Get assignement from other threads
    if ( mc%PTON) then
        print*, "calling PT_overrid ..."
        call PT_override(mc,md)
    else
        mc%repSufix=''
    endif
    iostr='data/r0'
    I=0;
    print*, 'calling saveR...'
    call MCvar_saveR(mc,md,iostr,0)
   
    iostr='data/params'
    print*, "calling saveParameters..."
    call MCvar_saveParameters(mc,iostr)

    iostr='data/u0'
    print*, "calling saveU..."
    call MCvar_saveU(mc,md,iostr)


 else

    PRINT*, '-----load simulation-----'
    iostr='putBinaryFileNameHere'
    stop 1
    call MCvar_readBindary(mc,md,iostr)

 endif
 call MCvar_printDescription(mc)

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!              Begin simulation
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print*, 'Beginning simulation'
  mc%IND=1

  DO WHILE ((mc%IND).LE.mc%INDMAX) 

     if (mc%IND.LE.mc%NNOINT) then
         INTON=0
     else
         INTON=1
     endif
     ! for coupling schedule
     if (mc%UseSchedule) then
         call strength_schedule(mc)
     endif
     

!   * Perform a MC simulation *
    call MCsim(mc,md,mc%NSTEP,INTON,rand_stat)

!    Save the conformation and the metrics
    write (fileind,'(I4)'), mc%IND

    !Save various energy contiributions to file 
    iostr='data/out1'
    call MCvar_appendEnergyData(mc,iostr)
     
    !part 2.5 - adaptations
    iostr='data/out3'
    call MCvar_appendAdaptData(mc,iostr)

    !print Phi
    write(iostr,"(I6)"), mc%IND
    iostr='data/phi' // trim(adjustL(iostr))
    call MCVar_savePHI(mc,md,iostr)    

    !part 3 - R
    write(iostr,"(I6)"), mc%IND
    iostr='data/r' // trim(adjustL(iostr))
    call MCvar_saveR(mc,md,iostr,0)
    
    !part 4 - U
    write(iostr,"(I6)"), mc%IND
    iostr='data/u' // trim(adjustL(iostr))
    call MCvar_saveU(mc,md,iostr)


    PRINT*, '________________________________________'
    PRINT*, 'Time point ',mc%IND, ' out of', mc%INDMAX
    call MCvar_printEnergies(mc)
    call MCvar_printWindowStats(mc)
    !call MCvar_printPhi(mc,md)
    mc%IND=mc%IND+1    
  ENDDO
  
END
!---------------------------------------------------------------*
