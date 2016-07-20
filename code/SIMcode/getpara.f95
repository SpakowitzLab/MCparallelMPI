! *---------------------------------------------------------------*
      
      SUBROUTINE getpara(PARA,EPS,L0,LBOX)
      

      DOUBLE PRECISION PARA(10)
      DOUBLE PRECISION DEL
      DOUBLE PRECISION PVEC(679,8)
      INTEGER IND,CRS
      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA
      DOUBLE PRECISION XIR,XIU
      DOUBLE PRECISION L0       ! Equilibrium segment length
      !DOUBLE PRECISION LHC      ! Length of HC int
      DOUBLE PRECISION VHC      ! HC strength
      DOUBLE PRECISION M
      !DOUBLE PRECISION DT
      INTEGER I
      DOUBLE PRECISION LBOX   ! Box length (approximate)
      !DOUBLE PRECISION L, LP
      DOUBLE PRECISION EPS   ! Elasticity l0/(2lp)
      
!     Load in the parameters for the simulation

      DEL=2.*EPS
  
      VHC=0.
      XIU=0.
      XIR=0.
  

!     Load the tabulated parameters
  
      OPEN (UNIT=5,FILE='input/dssWLCparams',STATUS='OLD')
      DO 10 I=1,679
         READ(5,*) PVEC(I,1),PVEC(I,2),PVEC(I,3),PVEC(I,4),PVEC(I,5), &
                   PVEC(I,6),PVEC(I,7),PVEC(I,8)
 10   CONTINUE 
      CLOSE(5)
      
      if (DEL.LT.PVEC(1,1)) then
         DEL=PVEC(1,1)
      endif
      if (DEL.GT.PVEC(679,1)) then
         DEL=PVEC(679,1)
      endif
      
      CRS=0
      IND=1
      do while (CRS.EQ.0)
         if (DEL.LE.PVEC(IND,1)) then
            CRS=1
         else
            IND=IND+1
         endif
      enddo
      
      I=2 
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      EB=M*(DEL-PVEC(IND,1))+PVEC(IND,I)
      
      I=3 
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      GAM=M*(DEL-PVEC(IND,1))+PVEC(IND,I)
      
      I=4
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      EPAR=M*(DEL-PVEC(IND,1))+PVEC(IND,I)
      
      I=5
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      EPERP=M*(DEL-PVEC(IND,1))+PVEC(IND,I)
      
      I=6
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      ETA=M*(DEL-PVEC(IND,1))+PVEC(IND,I)
      
      I=7
      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
      XIU=M*(DEL-PVEC(IND,1))+PVEC(IND,I)
      
!      I=8
!      M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
!      DT=XIU*(M*(DEL-PVEC(IND,1))+PVEC(IND,I))
      
      EB=EB/DEL
      EPAR=EPAR*DEL/L0**2.
      EPERP=EPERP*DEL/L0**2.
      ETA=ETA*DEL/L0
      GAM=L0*GAM
            
      PARA(1)=EB
      PARA(2)=EPAR
      PARA(3)=EPERP
      PARA(4)=GAM
      PARA(5)=ETA
      PARA(6)=0.
      PARA(7)=0.
      PARA(8)=LBOX
      PARA(9)=0.
      PARA(10)=0.

!      PRINT*, '>>Parameters == ', EPS, L0, DEL, PVEC(IND,1), PVEC(IND,2), PVEC(IND,3)      

      RETURN     
      END
      
!---------------------------------------------------------------*
