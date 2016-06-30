!---------------------------------------------------------------!
      
!     
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!     
!     Corrections to force magnitude made 6-3-04.
!     
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!
!     Edited by Quinn in 2016
      
SUBROUTINE MC_int2(mc,md,I1,I2,I3,I4,initialize)
use simMod
use setPrecision
IMPLICIT NONE

!   iputs
TYPE(MCvar) mc   ! <---- Contains output
TYPE(MCData) md
LOGICAL initialize        ! if true calculate absolute energy
INTEGER I1                ! Test bead position 1
INTEGER I2                ! Test bead position 2
INTEGER I3                ! Test bead, first bead of second section
INTEGER I4                ! Test bead, second bead of second section

!   Internal variables
INTEGER I,J               ! For looping over bins
INTEGER IB                ! Bead index
INTEGER IB2               ! Index you are swapping with
INTEGER rrdr ! -1 if r, 1 if r+dr
DOUBLE PRECISION PHIPoly    ! Total fraction polymer
INTEGER IX(2),IY(2),IZ(2)      
DOUBLE PRECISION WX(2),WY(2),WZ(2)
DOUBLE PRECISION WTOT       ! total weight ascribed to bin
DOUBLE PRECISION RBIN(3)    ! bead position
INTEGER INDBIN              ! index of bin
INTEGER ISX,ISY,ISZ 
DOUBLE PRECISION VV         ! one of Vol
LOGICAL isA   ! The bead is of type A

! Copy so I don't have to type mc% everywhere
DOUBLE PRECISION LBOX(3)
DOUBLE PRECISION DEL
INTEGER NBINX(3)
DOUBLE PRECISION temp    !for speeding up code
LBOX=mc%LBOX
DEL=mc%DEL
NBINX=mC%NBINX

! -------------------------------------------------------------
!
!  Calculate change (or value if initialize) of phi for A and B
!
!--------------------------------------------------------------
if (initialize) then
    print*,"Error, Don't use MC_int2 for initialization"
    stop 1
endif

mc%NPHI=0
do IB=I1,I2
  IB2=IB+I3-I1
  !No need to do calculation if identities are the same
  if (md%AB(IB).eq.md%ABP(IB2)) cycle
  do rrdr=-1,1,2
   ! on initialize only add current position
   ! otherwise subract current and add new
   if (rrdr.eq.-1) then
       RBIN(1)=md%R(IB,1)
       RBIN(2)=md%R(IB,2)
       RBIN(3)=md%R(IB,3)
       isA=md%AB(IB).eq.1
   else     
       RBIN(1)=md%RP(IB,1)
       RBIN(2)=md%RP(IB,2)
       RBIN(3)=md%RP(IB,3)
       isA=md%ABP(IB).eq.1
   endif
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   SELECT CASE (mc%confineType)
   CASE (0) ! Box from 0-LBOX, Bins split by boundaries
       ! Periodic BC
       RBIN(1)=RBIN(1)-floor(RBIN(1)/LBOX(1))*LBOX(1)
       RBIN(2)=RBIN(2)-floor(RBIN(2)/LBOX(2))*LBOX(2)
       RBIN(3)=RBIN(3)-floor(RBIN(3)/LBOX(3))*LBOX(3)
  
       ! Binning  
       IX(1)=ceiling(RBIN(1)/DEL)
       IY(1)=ceiling(RBIN(2)/DEL)
       IZ(1)=ceiling(RBIN(3)/DEL)
       
       IX(2)=IX(1)-1
       IY(2)=IY(1)-1
       IZ(2)=IZ(1)-1
       
       ! Calculate the bin weighting
       WX(2)=(DEL*IX(1)-RBIN(1))/DEL   
       WX(1)=1.0_dp-WX(2)              
       WY(2)=(DEL*IY(1)-RBIN(2))/DEL   
       WY(1)=1.0_dp-WY(2)              
       WZ(2)=(DEL*IZ(1)-RBIN(3))/DEL   
       WZ(1)=1.0_dp-WZ(2)              
   
       ! Periodic BC on Bins:
       IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX(1))) * NBINX(1)
       IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX(1))) * NBINX(1)
       IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX(2))) * NBINX(2)
       IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX(2))) * NBINX(2)
       IZ(1)=IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBINX(3))) * NBINX(3)
       IZ(2)=IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBINX(3))) * NBINX(3)
   CASE (1)
       ! Periodic BC
       RBIN(1)=RBIN(1)-floor(RBIN(1)/LBOX(2))*LBOX(1)
       RBIN(2)=RBIN(2)-floor(RBIN(2)/LBOX(1))*LBOX(2)
  
       ! Binning  
       IX(1)=ceiling(RBIN(1)/DEL)
       IY(1)=ceiling(RBIN(2)/DEL)
       IZ(1)=nint(RBIN(3)/DEL)+1 ! Note 1.0 so that box centers are on half intigers 
       
       IX(2)=IX(1)-1
       IY(2)=IY(1)-1
       IZ(2)=IZ(1)-1
       
       ! Calculate the bin weighting
       WX(2)=(DEL*IX(1)-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
       WX(1)=1.0_dp-WX(2)              ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
       WY(2)=(DEL*IY(1)-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
       WY(1)=1.0_dp-WY(2)              ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
       WZ(2)=(DEL*IZ(1)-0.5_dp*DEL-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
       WZ(1)=1.0_dp-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)
   
       if ((WZ(1).lt.0).OR.(WZ(2).lt.0)) then
           print*, "negitive W"
           stop 1
       endif

       ! Periodic BC on Bins:
       IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX(1))) * NBINX(1)
       IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX(1))) * NBINX(1)
       IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX(2))) * NBINX(2)
       IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX(2))) * NBINX(2)
   CASE (2) ! Box confinement
       ! Binning  
       IX(1)=nint(RBIN(1)/DEL)+1 
       IY(1)=nint(RBIN(2)/DEL)+1
       IZ(1)=nint(RBIN(3)/DEL)+1 ! Note +1 because fortran starts a 1 
       
       IX(2)=IX(1)-1
       IY(2)=IY(1)-1
       IZ(2)=IZ(1)-1
        
       ! Calculate the bin weighting
       WX(2)=(DEL*IX(1)-0.5_dp*DEL-RBIN(1))/DEL
       WX(1)=1.0_dp-WX(2)             
       WY(2)=(DEL*IY(1)-0.5_dp*DEL-RBIN(2))/DEL 
       WY(1)=1.0_dp-WY(2)                           
       WZ(2)=(DEL*IZ(1)-0.5_dp*DEL-RBIN(3))/DEL    
       WZ(1)=1.0_dp-WZ(2)                  
   CASE (3) 
       ! Binning  
       IX(1)=nint(RBIN(1)/DEL)+1 
       IY(1)=nint(RBIN(2)/DEL)+1
       IZ(1)=nint(RBIN(3)/DEL)+1 ! Note +1 because fortran starts a 1 
           
       IX(2)=IX(1)-1
       IY(2)=IY(1)-1
       IZ(2)=IZ(1)-1
        
       ! Calculate the bin weighting
       WX(2)=(DEL*IX(1)-0.5_dp*DEL-RBIN(1))/DEL  
       WX(1)=1.0_dp-WX(2)                          
       WY(2)=(DEL*IY(1)-0.5_dp*DEL-RBIN(2))/DEL  
       WY(1)=1.0_dp-WY(2)                          
       WZ(2)=(DEL*IZ(1)-0.5_dp*DEL-RBIN(3))/DEL  
       WZ(1)=1.0_dp-WZ(2)                           
   CASE (4) ! Box from 0-LBOX, Bins split by boundaries
       ! Periodic BC
       RBIN(1)=RBIN(1)-floor(RBIN(1)/LBOX(1))*LBOX(1)
       RBIN(2)=RBIN(2)-floor(RBIN(2)/LBOX(2))*LBOX(2)
       RBIN(3)=RBIN(3)-floor(RBIN(3)/LBOX(3))*LBOX(3)
  
       ! Binning  
       IX(1)=ceiling(RBIN(1)/DEL)
       IY(1)=ceiling(RBIN(2)/DEL)
       IZ(1)=ceiling(RBIN(3)/DEL)
       
       IX(2)=IX(1)-1
       IY(2)=IY(1)-1
       IZ(2)=IZ(1)-1
       
       ! Calculate the bin weighting
       WX(2)=(DEL*IX(1)-RBIN(1))/DEL   
       WX(1)=1.0_dp-WX(2)              
       WY(2)=(DEL*IY(1)-RBIN(2))/DEL   
       WY(1)=1.0_dp-WY(2)              
       WZ(2)=(DEL*IZ(1)-RBIN(3))/DEL   
       WZ(1)=1.0_dp-WZ(2)              
   
       ! Periodic BC on Bins:
       IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX(1))) * NBINX(1)
       IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX(1))) * NBINX(1)
       IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX(2))) * NBINX(2)
       IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX(2))) * NBINX(2)
       IZ(1)=IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBINX(3))) * NBINX(3)
       IZ(2)=IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBINX(3))) * NBINX(3)
   END SELECT
   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (isA) then
       do ISX=1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX(1)+1))) CYCLE
          do ISY=1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX(2)+1))) CYCLE
             do ISZ=1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX(3)+1))) cycle
                WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)
                ! Generate list of which phi's change and by how much
                I=mc%NPHI
                do 
                   if (I.eq.0) then
                      mc%NPHI=mc%NPHI+1
                      md%INDPHI(mc%NPHI)=INDBIN
                      temp=rrdr*WTOT*mc%V/md%Vol(INDBIN)
                      md%DPHIA(mc%NPHI)=temp
                      md%DPHIB(mc%NPHI)=-temp
                      exit
                   elseif (INDBIN.EQ.md%INDPHI(I)) then
                      temp=rrdr*WTOT*mc%V/md%Vol(INDBIN)
                      md%DPHIA(I)=md%DPHIA(I)+temp
                      md%DPHIB(I)=md%DPHIB(I)-temp
                      exit
                   else
                      I=I-1
                   endif                     
                enddo
             enddo
          enddo
       enddo
   else
       do ISX=1,2
          if ((IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX(1)+1))) CYCLE
          do ISY=1,2
             if ((IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX(2)+1))) CYCLE
             do ISZ=1,2
                if ((IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX(3)+1))) cycle
                WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX(1)+(IZ(ISZ)-1)*NBINX(1)*NBINX(2)
                ! Generate list of which phi's change and by how much
                I=mc%NPHI
                do 
                   if (I.eq.0) then
                      mc%NPHI=mc%NPHI+1
                      md%INDPHI(mc%NPHI)=INDBIN
                      temp=rrdr*WTOT*mc%V/md%Vol(INDBIN)
                      md%DPHIA(mc%NPHI)=-temp
                      md%DPHIB(mc%NPHI)=temp
                      exit
                   elseif (INDBIN.EQ.md%INDPHI(I)) then
                      temp=rrdr*WTOT*mc%V/md%Vol(INDBIN)
                      md%DPHIA(I)=md%DPHIA(I)-temp
                      md%DPHIB(I)=md%DPHIB(I)+temp
                      exit
                   else
                      I=I-1
                   endif                     
                enddo
             enddo !ISZ
          enddo !ISY
       enddo !ISX 
   endif
 enddo ! loop over rrdr.  A.k.a new and old
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte change in (or abosolute) energy
!
!---------------------------------------------------------------------
mc%DEChi=0.0_dp
mc%DECouple=0.0_dp
mc%DEKap=0.0_dp
if (mc%simType.eq.0) then ! Melt Hamiltonian
    do I=1,mc%NPHI
        J=md%INDPHI(I)
        VV=md%Vol(J)
        if (VV.le.0.1_dp) CYCLE
        ! new
        mc%DEChi=mc%DEChi+VV*(mc%CHI/mc%V)*((md%PHIA(J)+md%DPHIA(I))*(md%PHIB(J)+md%DPHIB(I)))
        mc%DEKap=mc%DEKap+VV*(mc%KAP/mc%V)*((md%PHIA(J)+md%DPHIA(I)+md%PHIB(J)+md%DPHIB(I)-1.0_dp)**2)
        ! minus old
        mc%DEChi=mc%DEChi-VV*(mc%CHI/mc%V)*(md%PHIA(J)*md%PHIB(J))
        mc%DEKap=mc%DEKap-VV*(mc%KAP/mc%V)*((md%PHIA(J)+md%PHIB(J)-1.0_dp)**2)
        
    enddo
elseif(mc%simType.eq.1) then ! Chromatin Hamiltonian
    do I=1,mc%NPHI
        J=md%INDPHI(I)
        VV=md%Vol(J)
        if (VV.le.0.1_dp) CYCLE
        ! new ...
        PHIPoly=md%PHIA(J)+md%DPHIA(I)+md%PHIB(J)+md%DPHIB(I)
        mc%DEChi=mc%DEChi+VV*(mc%CHI/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
        mc%DECouple=mc%DECouple+VV*mc%HP1_Bind*(md%PHIA(J)+md%DPHIA(I))**2
        if(PHIPoly.GT.1.0_dp) then
           mc%DEKap=mc%DEKap+VV*(mc%KAP/mc%V)*(PHIPoly-1.0_dp)**2
        endif
        ! minus old
        PHIPoly=md%PHIA(J)+md%PHIB(J)
        mc%DEChi=mc%DEChi-VV*(mc%CHI/mc%V)*PHIPoly*(1.0_dp-PHIPoly)
        mc%DECouple=mc%DECouple-VV*mc%HP1_Bind*(md%PHIA(J))**2
        if(PHIPoly.GT.1.0_dp) then
           mc%DEKap=mc%DEKap-VV*(mc%KAP/mc%V)*(PHIPoly-1.0_dp)**2
        endif 
    enddo
endif

! Discount if interaction are only partial on
mc%DEChi=mc%DECHI*mc%CHI_ON
mc%DECouple=mc%DECouple*mc%Couple_ON
mc%DEKap=mc%DEKap*mc%KAP_ON

RETURN
END

!---------------------------------------------------------------!
