!---------------------------------------------------------------!
      
!     
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!     
!     Corrections to force magnitude made 6-3-04.
!     
!     Andrew Spakowitz
!     Written 6-29-04
      
      SUBROUTINE MC_int(DEINT,R,AB,NT,NBIN, &
                        V,CHI,KAP,LBOX,DEL,PHIA,PHIB,DPHIA,DPHIB, &
                        INDPHI,NPHI,RP,I1,I2,HP1_Bind,ABP,confineType, &
                        Vol,initialize,NBINX)
      IMPLICIT NONE

      DOUBLE PRECISION R(NT,3)  ! Conformation of polymer chains
      DOUBLE PRECISION RP(NT,3) ! Conformation of polymer chains
      INTEGER AB(NT)            ! Chemical identity of beads    
      INTEGER ABP(NT)           ! Proposed chemical identities
      INTEGER NT                ! Number of beads
      
!     Simulation input variables
      
      DOUBLE PRECISION V        ! Monomer volume
      DOUBLE PRECISION HP1_Bind ! HP1_Binding energy
      DOUBLE PRECISION CHI      ! Solvent polymer chi parameter value
      DOUBLE PRECISION KAP      ! Compressibility value
      DOUBLE PRECISION LBOX     ! Simulation box size (approximate)
      DOUBLE PRECISION DEL      ! Discretization size (approximate)
      INTEGER POLY              ! Polydisperse (step-growth stat)
      !DOUBLE PRECISION P        ! Degree of polymerization (polydisperse)
      INTEGER N                 ! Number of beads (monodisperse)
      
!     Variables for density calculation
      
      DOUBLE PRECISION PHIA(NBIN) ! Volume fraction of A
      DOUBLE PRECISION PHIB(NBIN) ! Volume fraction of B
      DOUBLE PRECISION DPHIA(NBIN) ! Delta Volume fraction of A
      DOUBLE PRECISION DPHIB(NBIN) ! Delta Volume fraction of B
      DOUBLE PRECISION PHIPoly    ! Total fraction polymer
      INTEGER INDPHI(NBIN)      ! Indices of the phi
      INTEGER NPHI              ! Number of phi values that change 
      DOUBLE PRECISION DEINT    ! Change in Self-interaction energy
      INTEGER I1                ! Test bead position 1
      INTEGER I2                ! Test bead position 2
      INTEGER I,J
      INTEGER IB
      INTEGER StopLoop
      INTEGER NBINX
      INTEGER rrdr ! -1 if r, 1 if r+dr

      INTEGER IX(2),IY(2),IZ(2)      
      DOUBLE PRECISION WX(2),WY(2),WZ(2)
      DOUBLE PRECISION WTOT
      DOUBLE PRECISION RBIN(3)
      INTEGER INDBIN
      INTEGER ISX,ISY,ISZ
      
      DOUBLE PRECISION Vol(NBIN)  ! Volume of bins 
      DOUBLE PRECISION VV         ! one of Vol
      INTEGER confineType
      INTEGER NBIN   ! Number of bins in simulation
      INTEGER initialize
      LOGICAL isA   ! The bead is of type A

      if ((initialize.ne.0).and.(initialize.ne.1)) then
          print*, "error in MC_int, initialize must be 1 or 0"
          stop 1
      endif

      do 10 I=1,NBIN
         DPHIA(I)=0.
         DPHIB(I)=0.
         INDPHI(I)=0
 10   continue
!      NBINX=nint(LBOX/DEL)
!      if (NBINX**3.ne.NBIN) then
!           print*, "error in MC_int. Wrong number of bins"
!      endif

      NPHI=0
      do 20 IB=I1,I2
        do 25 rrdr=(-1+2*initialize),1,2
         if ((rrdr.eq.-1).or.(initialize.eq.1)) then
             RBIN(1)=R(IB,1)
             RBIN(2)=R(IB,2)
             RBIN(3)=R(IB,3)
             isA=AB(IB).eq.1
         else     
             RBIN(1)=RP(IB,1)
             RBIN(2)=RP(IB,2)
             RBIN(3)=RP(IB,3)
             isA=ABP(IB).eq.1
         endif
         SELECT CASE (confineType)
         CASE (0) ! Box from 0-LBOX, Bins split by boundaries
             ! Periodic BC
             RBIN(1)=RBIN(1)-nint(RBIN(1)/LBOX-0.5)*LBOX
             RBIN(2)=RBIN(2)-nint(RBIN(2)/LBOX-0.5)*LBOX
             RBIN(3)=RBIN(3)-nint(RBIN(3)/LBOX-0.5)*LBOX
        
             ! Binning  
             IX(1)=nint(RBIN(1)/DEL+0.5)
             IY(1)=nint(RBIN(2)/DEL+0.5)
             IZ(1)=nint(RBIN(3)/DEL+0.5)
             
             IX(2)=IX(1)-1
             IY(2)=IY(1)-1
             IZ(2)=IZ(1)-1
             
             ! Calculate the bin weighting
             WX(2)=(DEL*IX(1)-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
             WX(1)=1-WX(2)                   ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
             WY(1)=(DEL*IY(1)-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
             WY(2)=1-WY(1)                   ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
             WZ(2)=(DEL*IZ(1)-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
             WZ(1)=1-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)
         
             ! Periodic BC on Bins:
             IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX)) * NBINX
             IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX)) * NBINX
             IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX)) * NBINX
             IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX)) * NBINX
             IZ(1)=IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBINX)) * NBINX
             IZ(2)=IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBINX)) * NBINX
         CASE (1)
             ! Periodic BC
             RBIN(1)=RBIN(1)-nint(RBIN(1)/LBOX-0.5)*LBOX
             RBIN(2)=RBIN(2)-nint(RBIN(2)/LBOX-0.5)*LBOX
        
             ! Binning  
             IX(1)=nint(RBIN(1)/DEL+0.5)
             IY(1)=nint(RBIN(2)/DEL+0.5)
             IZ(1)=nint(RBIN(3)/DEL+1.0) ! Note 1.0 so that box centers are on half intigers 
             
             IX(2)=IX(1)-1
             IY(2)=IY(1)-1
             IZ(2)=IZ(1)-1
             
             ! Calculate the bin weighting
             WX(2)=(DEL*IX(1)-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
             WX(1)=1-WX(2)                   ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
             WY(1)=(DEL*IY(1)-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
             WY(2)=1-WY(1)                   ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
             WZ(2)=(DEL*IZ(1)-0.5*DEL-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
             WZ(1)=1-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)
         
             if ((WZ(1).lt.0).OR.(WZ(2).lt.0)) then
                 print*, "negitive W"
                 stop 1
             endif

             ! Periodic BC on Bins:
             IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX)) * NBINX
             IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX)) * NBINX
             IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX)) * NBINX
             IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX)) * NBINX
         CASE (2) ! Box confinement
             ! Binning  
             IX(1)=nint(RBIN(1)/DEL+1.0)
             IY(1)=nint(RBIN(2)/DEL+1.0)
             IZ(1)=nint(RBIN(3)/DEL+1.0) ! Note 1.0 so that box centers are on half intigers 
             
             IX(2)=IX(1)-1
             IY(2)=IY(1)-1
             IZ(2)=IZ(1)-1
              
             ! Calculate the bin weighting
             WX(2)=(DEL*IX(1)-0.5*DEL-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
             WX(1)=1-WX(2)                   ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
             WY(1)=(DEL*IY(1)-0.5*DEL-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
             WY(2)=1-WY(1)                   ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
             WZ(2)=(DEL*IZ(1)-0.5*DEL-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
             WZ(1)=1-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)
         CASE (3)
             ! Binning  
             IX(1)=nint(RBIN(1)/DEL+1.0)
             IY(1)=nint(RBIN(2)/DEL+1.0)
             IZ(1)=nint(RBIN(3)/DEL+1.0) ! Note 1.0 so that box centers are on half intigers 
             !write(*,"(A,I4,A,I4,A,I4)"),"IX(1):",IX(1),"  IY(1):",IY(1),"  IZ(1)",IZ(1)
                 
             IX(2)=IX(1)-1
             IY(2)=IY(1)-1
             IZ(2)=IZ(1)-1
              
             ! Calculate the bin weighting
             WX(2)=(DEL*IX(1)-0.5*DEL-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
             WX(1)=1-WX(2)                   ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
             WY(1)=(DEL*IY(1)-0.5*DEL-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
             WY(2)=1-WY(1)                   ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
             WZ(2)=(DEL*IZ(1)-0.5*DEL-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
             WZ(1)=1-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL) 
             !write(*,"(A,f4.3,A,f4.3,A,f4.3)"),"WX(1):",WX(1),"  WY(1):",WY(1),"  WZ(1)",WZ(1)
         END SELECT
!   Add or Subtract volume fraction with weighting from each bin
         do 30 ISX=1,2
            do 40 ISY=1,2
               do 50 ISZ=1,2
                  if ( (IX(ISX).le.0).OR.(IX(ISX).ge.(NBINX+1)).OR. &
                       (IY(ISY).le.0).OR.(IY(ISY).ge.(NBINX+1)).OR. &
                       (IZ(ISZ).le.0).OR.(IZ(ISZ).ge.(NBINX+1)) ) then
                       ! interaction with boundary goes here
                       cycle
                  endif
                  WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                  INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX+(IZ(ISZ)-1)*NBINX**2
                  if (WTOT.gt.1.0) then
                      print*, "INDBIN", INDBIN
                      print*, "error in MC_int WTOT must be leq 1.0"
                      stop 1
                  endif

                  if (initialize.eq.1) then
                      ! Set all phi values on initialize
                      if (isA) then
                          PHIA(INDBIN)=PHIA(INDBIN)+WTOT*V !/Vol(INDBIN)
                      else
                          PHIB(INDBIN)=PHIB(INDBIN)+WTOT*V !/Vol(INDBIN)
                      endif
                  else
                      ! Generate list of which phi's change and by how much
                      I=1
                      StopLoop=0
                      do while (StopLoop.EQ.0)
                         if (INDBIN.EQ.INDPHI(I).AND.I.LE.NPHI) then
                            StopLoop=1
                            if (isA) then
                               DPHIA(I)=DPHIA(I)+rrdr*WTOT*V !/Vol(INDBIN)
                            else
                               DPHIB(I)=DPHIB(I)+rrdr*WTOT*V !/Vol(INDBIN)
                            endif
                         elseif (I.GT.NPHI) then
                            StopLoop=1
                            NPHI=NPHI+1
                            INDPHI(I)=INDBIN

                            if (isA) then
                               DPHIA(I)=DPHIA(I)+rrdr*WTOT*V !/Vol(INDBIN)
                            else
                               DPHIB(I)=DPHIB(I)+rrdr*WTOT*V !/Vol(INDBIN)
                            endif
                         else
                            I=I+1
                         endif                     
                      enddo
                  endif
 50            continue
 40         continue
 30      continue 
 25    continue ! loop over rrdr.  A.k.a new and old
 20   continue ! loop over IB  A.k.a. beads

      if (initialize.eq.1) then
          ! calculate absolute energy
          DEINT=0.0
          do I=1,NBIN
              if (Vol(I).eq.0.0) then
                  CYCLE
              endif
              VV=Vol(I)
!             DEINT=DEINT+(DEL**3.)*(CHI/V)*((PHIA(J))*(PHIB(J))-PHIA(J)*PHIB(J)) &
!                   +(DEL**3.)*(KAP/V)*((PHIA(J)+PHIB(J)-1.)**2.-(PHIA(J)+PHIB(J)-1.)**2.)
              ! new ...
              PHIPoly=PHIA(I)+PHIB(I)
              DEINT=DEINT&
                    +VV*(CHI/V)*PHIPoly*(1-PHIPoly)&   
                    +VV*HP1_Bind*(PHIA(I))**2
              if(PHIPoly.GT.1) then
                 DEINT=DEINT+VV*(KAP/V)*(PHIPoly-1.)**2
              endif
          enddo        
      

      else
          DEINT=0.
          do 100 I=1,NPHI
              J=INDPHI(I)
              VV=Vol(J)
!             DEINT=DEINT+(DEL**3.)*(CHI/V)*((PHIA(J)+DPHIA(I))*(PHIB(J)+DPHIB(I))-PHIA(J)*PHIB(J)) &
!                   +(DEL**3.)*(KAP/V)*((PHIA(J)+DPHIA(I)+PHIB(J)+DPHIB(I)-1.)**2.-(PHIA(J)+PHIB(J)-1.)**2.)
              ! new ...
              PHIPoly=PHIA(J)+DPHIA(I)+PHIB(J)+DPHIB(I)
              DEINT=DEINT&
                    +VV*(CHI/V)*PHIPoly*(1-PHIPoly)&   
                    +VV*HP1_Bind*(PHIA(J)+DPHIA(I))**2
              if(PHIPoly.GT.1) then
                 DEINT=DEINT+VV*(KAP/V)*(PHIPoly-1.)**2
              endif
              ! minus old
              PHIPoly=PHIA(J)+PHIB(J)
              DEINT=DEINT&
                    -VV*(CHI/V)*PHIPoly*(1-PHIPoly)&
                    -VV*HP1_Bind*PHIA(J)**2
              if(PHIPoly.GT.1) then
                 DEINT=DEINT-VV*(KAP/V)*(PHIPoly-1.)**2
              endif
              
 100      continue
      endif
      RETURN
      END
      
!---------------------------------------------------------------!
