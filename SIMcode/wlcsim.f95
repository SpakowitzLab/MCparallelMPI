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
  INTEGER IND, I            ! Ind in series
  character*4 fileind       ! Index of output
  character*16 iostr       ! File for output
  INTEGER INDEND            ! Restart index
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
  
  print*, "Calculating Bin volumes"
  !call MCvar_printDiscription(mc) 
  call MC_caclVolume(mc%confineType,mc%NBINX,mc%DEL, mc%LBox, &
                     md%Vol,rand_stat)  ! calculate partial volumes
  print*, "Done Calculating Bin volumes"

  INQUIRE (FILE = 'data/out1', exist = restart)
  if (.NOT.restart) then

    PRINT*, '-----new simulation-----'

    
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

    INDEND = 0

 else

    PRINT*, '-----load simulation-----'
    iostr='putBinaryFileNameHere'
    INDEND=0; INDEND=1/INDEND ! make this a variable of simmod before use
    call MCvar_readBindary(mc,md,iostr)

 endif

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!              Begin simulation
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print*, 'Beginning simulation'
  IND=1

  DO WHILE ((IND+INDEND).LE.mc%INDMAX)  ! INDEND is for sims. that are restarted

     if ((IND+INDEND).LE.mc%NNOINT) then
         INTON=0
     else
         INTON=1
     endif
     ! for coupling schedule
     if (.FALSE.) then
         call strength_schedule(IND+INDEND,mc%HP1_bind)
     endif
     

!   * Perform a MC simulation *
    call MCsim(mc,md,mc%NSTEP,INTON,rand_stat)

!    Save the conformation and the metrics
    write (fileind,'(I4)'), INDEND+IND

    !Save various energy contiributions to file 
    iostr='data/out1'
    I=INDEND+IND
    call MCvar_appendEnergyData(mc,iostr,I)
     
    !part 2.5 - adaptations
    iostr='data/out3'
    call MCvar_appendAdaptData(mc,iostr,I)
    

    !part 3 - R
    write(iostr,"(I6)"), I
    iostr='data/r' // trim(adjustL(iostr))
    call MCvar_saveR(mc,md,iostr,0)
    
    PRINT*, '________________________________________'
    PRINT*, 'Time point ',IND+INDEND, ' out of', mc%INDMAX
    call MCvar_printEnergies(mc)
    call MCvar_printWindowStats(mc)
    IND=IND+1    
  ENDDO
  
END
Subroutine strength_schedule(ind,HP1_bind)
    use setPrecision
    implicit none
    Integer ind
    Double PRECISION HP1_bind
    Double precision maximum
    maximum=-28.0_dp
    
    if(ind.lt.101) then
        HP1_bind=0.0_dp
    elseif(ind.lt.111) then
        HP1_bind=maximum*0.1_dp
    elseif(ind.lt.121) then
        HP1_bind=maximum*0.2_dp
    elseif(ind.lt.131) then
        HP1_bind=maximum*0.3_dp
    elseif(ind.lt.141) then
        HP1_bind=maximum*0.4_dp
    elseif(ind.lt.151) then
        HP1_bind=maximum*0.5_dp
    elseif(ind.lt.161) then
        HP1_bind=maximum*0.6_dp
    elseif(ind.lt.171) then
        HP1_bind=maximum*0.7_dp
    elseif(ind.lt.181) then
        HP1_bind=maximum*0.8_dp
    elseif(ind.lt.191) then
        HP1_bind=maximum*0.9_dp
    elseif(ind.lt.201) then
        HP1_bind=maximum*1.0_dp
    elseif(ind.lt.211) then
        HP1_bind=maximum*1.3_dp
    elseif(ind.lt.221) then
        HP1_bind=maximum*1.2_dp
    elseif(ind.lt.231) then
        HP1_bind=maximum*1.1_dp
    elseif(ind.lt.241) then
        HP1_bind=maximum*1.0_dp
    elseif(ind.lt.251) then
        HP1_bind=maximum*0.9_dp
    elseif(ind.lt.261) then
        HP1_bind=maximum*0.8_dp
    elseif(ind.lt.271) then
        HP1_bind=maximum*0.7_dp
    endif

end subroutine
!---------------------------------------------------------------*
