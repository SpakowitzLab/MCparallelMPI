!  ----  testing new move ------

      PROGRAM test
      use mt19937, only : grnd, sgrnd, rnorm, mt, mti
      
      IMPLICIT NONE

      INTEGER I
      
      DO I=1,10,2
          print*, I

      ENDDO




!      INTEGER, PARAMETER :: NT = 20 !total number of beads
!      INTEGER IP
!      INTEGER IB1
!      INTEGER IB2
!      INTEGER IT1
!      INTEGER IT2
!      INTEGER TEMP
!      INTEGER win   !window, will be changed
!      INTEGER NP    !Number of polymers
!      INTEGER N     !Beads per polymer
!      INTEGER BPM   !Beads per nomomer
!      INTEGER AB(NT)
!      INTEGER I
!      INTEGER J
!
!      NP=1;
!      N=20;
!      BPM=4;
!      win=8;
!     
!
!      Do I=1,NT
!          AB(I)=0
!      ENDDO      
!
!      print "(20I2)", AB
!
!      DO I=1,500
!         IP=nint(0.5+grnd()*NP)
!         IB1=nint(0.5+grnd()*N)
!         IB2=IB1+nint((grnd()-0.5)*(2*win+1))
!
!         if (IB2.LT.1) then
!            IB2=1
!         endif
!         if (IB2.GT.N) then
!            IB2=N
!         endif
!
!         if (IB2.LT.IB1) then
!            TEMP=IB1
!            IB1=IB2
!            IB2=TEMP
!         endif        
!         IT1=N*(IP-1)+IB1
!         IT2=N*(IP-1)+IB2
!
!         !keep binding constant within monomers
!         IT1=IT1-MOD(IT1-1,BPM)
!         IT2=IT2-MOD(IT2-1,BPM)+BPM-1
!
!         DO J=IT1,IT2
!             AB(J)=1-AB(J)
!         ENDDO
!         print "(20I2)", AB
!      ENDDO


      END
