!---------------------------------------------------------------*
!
!
!      
!     This subroutine performs a Monte Carlo simulation on the 
!     polymer chain.
!
!    Quinn Made Changes to this file starting on 12/15/15    
!
      
SUBROUTINE MCsim(mc,md,NSTEP,INTON,rand_stat)

    !use mt19937, only : grnd, sgrnd, rnorm, mt, mti
    use mersenne_twister
    use simMod

    IMPLICIT NONE 

    DOUBLE PRECISION, PARAMETER :: PI=3.141592654 ! Value of pi

    INTEGER NSTEP             ! Number of MC steps
    INTEGER INTON             ! Include polymer interactions
    
!   Variables for the simulation
    
    INTEGER ISTEP             ! Current MC step index
    DOUBLE PRECISION PROB     ! Calculated test prob
    DOUBLE PRECISION TEST     ! Random test variable
    INTEGER IB                ! Test bead 
    INTEGER IP                ! Test polymer 
    INTEGER IB1               ! Test bead position 1
    INTEGER IT1               ! Index of test bead 1
    INTEGER IB2               ! Test bead position 2
    INTEGER IT2               ! Index of test bead 2

    INTEGER I,J,K
    DOUBLE PRECISION R0(3)
    
    
    INTEGER MCTYPE                    ! Type of MC move
    INTEGER initialize
    
    DOUBLE PRECISION EB,EPAR,EPERP
    DOUBLE PRECISION GAM,ETA
    DOUBLE PRECISION XIR,XIU
    DOUBLE PRECISION LHC      ! Length of HC int
    DOUBLE PRECISION VHC      ! HC strength
    DOUBLE PRECISION phiTot, phiTot2 ! for testing

! Things for random number generator
    real urnd(1) ! single random number
    type(random_stat) rand_stat
!   Load the input parameters
    Type(MCvar) mc 
    Type(MCData) md 
    EB=   mc%PARA(1)
    EPAR= mc%PARA(2)
    EPERP=mc%PARA(3)
    GAM=  mc%PARA(4)
    ETA=  mc%PARA(5)
    XIR=  mc%PARA(6)
    XIU=  mc%PARA(7)
    LHC=  mc%PARA(9)
    VHC=  mc%PARA(10)
! -------------------------------------
!
!   initialize densities and energies 
!
! -------------------------------------
    print*, "initializeing..." ! make sure this doesn't take too much time
    ! --- Binding Energy ---
    md%ABP=0 ! set entire array to zero
    !  Notide that ABP and AB are intensionally swapped below
    IT1=1; IT2=mc%NT
    call MC_bind(mc%NT,mc%G,IT1,IT2,md%ABP,md%AB,md%METH, &
                 mc%EU,mc%EM,mc%DEBind,mc%mu)
    if(abs(mc%EBind-mc%DEBind).gt.0.0001) then
        print*, "Warning. Integrated binding enrgy:", &
                mc%EBind," while absolute binding energy:", &
                mc%DEBind
    endif
    mc%EBind=mc%DEBind

    ! --- Elastic Energy ---
    call energy_elas(mc%DEELAS,md%R,md%U,mc%NT,mc%N,mc%NP,mc%Para)
    if(abs((mc%EElas(1)+  mc%EElas(2)+ mc%EElas(3))-& 
           (mc%DEElas(1)+mc%DEElas(2)+mc%DEElas(3))).gt.0.0001) then
        print*, "Warning. Integrated elastic enrgy:", &
                (mc%EElas(1)+mc%EElas(2)+mc%EElas(3)),&
                " while absolute elastic energy:", &
                (mc%DEElas(1)+mc%DEElas(2)+mc%DEElas(3))
    endif
    mc%EElas=mc%DEElas ! copy array

    ! --- Interaction Energy ---
    mc%DEINT=0.
    if (INTON.EQ.1) then
        ! initialize phi
        IT1=1
        IT2=mc%NT ! need to set up all beads
        initialize=1
        do I=1,mc%NBIN
             md%PHIA(I)=0
             md%PHIB(I)=0
        enddo
        call MC_int(mc%DEINT,md%R,md%AB,mc%NT,mc%NBIN, &
                    mc%V,mc%CHI,mc%KAP,mc%LBOX,mc%DEL,md%PHIA,md%PHIB,md%DPHIA,md%DPHIB, &
                    md%INDPHI,mc%NPHI,md%RP,IT1,IT2,mc%HP1_Bind,md%ABP,mc%confineType, &
                    md%Vol,initialize,mc%NBINX)
        do I=1,mc%NBIN
            phiTot=phiTot+(md%PHIA(I)+md%PHIB(I))*md%Vol(I)
        enddo
        print*, "phiTot", phiTot," NT:",mc%NT
        if(abs(mc%Eint-mc%DEINT).gt. 0.0001) then
             print*, "Warning. Intigrated interaction energy:", & 
                     mc%Eint,"  while absolute interaction energy:", &
                     mc%DEINT 
        endif
        mc%Eint=mc%DEINT
    !   call r_to_phi(R,AB,NT,N,NP,NBIN, &
    !        V,CHI,KAP,LBOX,DEL,PHIA,PHIB)
        do I=1,mc%NBIN
            if (md%Vol(I).eq.0.0) then
                Cycle
            endif
            if (md%PHIA(I) .ne. md%PHIA(I)) then
                write(*,"(A,I5,A)"), "PHIA(",I,")=NaN"
                write(*,"(A,I5,A,f8.4)"), "Vol(",I,")=",md%Vol(I)
                stop 1
            endif
            if (md%PHIB(I) .ne. md%PHIB(I)) then
                write(*,"(A,I5,A)"), "PHIB(",I,")=NaN"
                write(*,"(A,I5,A,f8.4)"), "Vol(",I,")=",md%Vol(I)
                stop 1
            endif
            if (md%Vol(I) .ne. md%Vol(I)) then
                write(*,"(A,I5,A)"), "Vol(",I,")=NaN"
                stop 1
            endif
        enddo

    else
        do I=1,mc%NBIN
             md%PHIA(I)=0.0
             md%PHIB(I)=0.0
        enddo
    endif
    print*, "done initializing"
  
! -------------------------------------
!
!   Begin Monte Carlo simulation
!
! -------------------------------------
    ISTEP=1
    DO WHILE (ISTEP.LE.NSTEP)
        
       DO MCTYPE=1,7  

          if (mc%MOVEON(MCTYPE).EQ.0) then
             cycle
          endif
          call MC_move(md%R,md%U,md%RP,md%UP,mc%NT,mc%N,mc%NP, &
                       IP,IB1,IB2,IT1,IT2,MCTYPE, & 
                       mc%MCAMP,mc%WINDOW,md%AB,md%ABP,mc%G,&
                       rand_stat)
          
!   Calculate the change in compression and bending energy
          if (MCTYPE.NE.7) then
              call MC_eelas(mc%DEELAS,md%R,md%U,md%RP,md%UP,mc%NT,mc%N, &
                            mc%NP,IP,IB1,IB2, & 
                            IT1,IT2,EB,EPAR,EPERP,GAM,ETA)
              !print*, mc%DEELAS(1), mc%DEELAS(2), mc%DEELAS(3)
          else
              mc%DEELAS(1)=0.0
              mc%DEELAS(2)=0.0
              mc%DEELAS(3)=0.0
          endif
!   Calculate the change in the binding energy
          if (MCTYPE.EQ.7) then

              !print*, 'MCsim says EM:',EM,'EU',EU
              call MC_bind(mc%NT,mc%G,IT1,IT2,md%AB,md%ABP,md%METH,mc%EU,mc%EM, &
                           mc%DEBind,mc%mu)
          else
              mc%DEBind=0.0
          endif
         
!   Calculate the change in the self-interaction energy (actually all
!   interation energy, not just self?)
          if (INTON.EQ.1) then
             initialize=0
             call MC_int(mc%DEINT,md%R,md%AB,mc%NT,mc%NBIN, &
                  mc%V,mc%CHI,mc%KAP,mc%LBOX,mc%DEL,md%PHIA,md%PHIB,md%DPHIA,md%DPHIB, &
                  md%INDPHI,mc%NPHI,md%RP,IT1,IT2,mc%HP1_Bind,md%ABP,& 
                  mc%confineType,md%Vol,initialize,mc%NBINX)
          endif

!   Calculate the change in confinement energy
          if (MCTYPE.NE.7) then
              call MC_confine(mc%confineType, mc%LBox, md%RP, mc%NT, & 
                              IT1,IT2,mc%ECon)
          else
              mc%ECon=0.0;
          endif


!   Change the position if appropriate
          mc%ENERGY=mc%DEELAS(1)+mc%DEELAS(2)+mc%DEELAS(3)+mc%DEINT+mc%DEBind+mc%ECon
          
          PROB=exp(-mc%ENERGY)
          call random_number(urnd,rand_stat)
          TEST=urnd(1)
          if (TEST.LE.PROB) then
             if(MCTYPE.EQ.7) then
                 DO I=IT1,IT2
                      md%AB(I)=md%ABP(I)
                 ENDDO
             else
                 DO I=IT1,IT2
                     md%R(I,1)=md%RP(I,1)
                     md%R(I,2)=md%RP(I,2)
                     md%R(I,3)=md%RP(I,3)
                     md%U(I,1)=md%UP(I,1)
                     md%U(I,2)=md%UP(I,2)
                     md%U(I,3)=md%UP(I,3)
                 enddo
             endif
             if (mc%ECon.ne.0.0) then
                 print*, "error in MCsim, out of bounds "
                 stop 1
             endif
             mc%EBind=mc%EBind+mc%DEBind
             mc%EELAS(1)=mc%EELAS(1)+mc%DEELAS(1)
             mc%EELAS(2)=mc%EELAS(2)+mc%DEELAS(2)
             mc%EELAS(3)=mc%EELAS(3)+mc%DEELAS(3)
             if (INTON.EQ.1) then
                DO I=1,mc%NPHI
                   J=md%INDPHI(I)
                   md%PHIA(J)=md%PHIA(J)+md%DPHIA(I)
                   md%PHIB(J)=md%PHIB(J)+md%DPHIB(I)  
                enddo
                mc%Eint=mc%Eint+mc%DEINT
             endif

             mc%SUCCESS(MCTYPE)=mc%SUCCESS(MCTYPE)+1
          endif
!   Adapt the amplitude of step every NADAPT steps

          !amplitude and window adaptations
          if (mod(ISTEP,mc%NADAPT(MCTYPE)).EQ.0) then  ! Addapt ever NADAPT moves
             call MCvar_addapt(mc,MCTYPE,ISTEP,rand_stat)
           
             ! move each chain back if drifted though repeated BC 
             call MCvar_recenter(mc,md)  ! You don't need to do this if there is confinement
          endif


       enddo ! End of movetype loop

       !  -----  Parallel tempering ----
       IF ((mc%PTON.EQ.1).and.(mod(ISTEP,mc%NPT)).eq.0) THEN
          call replicaExchange(mc,md)
       ENDIF
      
       ! seps in this subroutine
       ISTEP=ISTEP+1
    enddo ! end of ISTEP loop
    
    RETURN      
END
    
!-------------------------------------------------------------*
