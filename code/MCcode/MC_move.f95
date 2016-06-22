!---------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!           
!    Quinn Made Changes to this file starting on 12/15/15
!       
!


! Find change in bead position for a crank-shaft type move
     
      SUBROUTINE MC_move(R,U,RP,UP,NT,NB,NP,IP,IB1,IB2,IT1,IT2,MCTYPE &
                        ,MCAMP,WINDOW,AB,ABP,BPM,rand_stat,winType)

      !use mt19937, only : grnd, sgrnd, rnorm, mt, mti
      use mersenne_twister      
      use setPrecision
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: PI=3.141592654 ! Value of pi 
      
      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION RP(NT,3)  ! Bead positions
      DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
      INTEGER NB                 ! Number of beads on a polymer
      INTEGER NP                ! Number of polymers
      INTEGER NT                ! Total beads in simulation 
      INTEGER BPM               ! Beads per monomer, aka G

      INTEGER IP                ! Test polymer 
      INTEGER IB1               ! Test bead position 1
      INTEGER IT1               ! Index of test bead 1
      INTEGER IB2               ! Test bead position 2
      INTEGER IT2               ! Index of test bead 2

      INTEGER I,J			! Test indices
! Things for random number generator
      type(random_stat) rand_stat  ! status of random number generator      
      real urand(3)  ! random vector
      real urnd(1) ! single random number
! Variables for the crank-shaft move

      DOUBLE PRECISION TA(3)    ! Axis of rotation
      DOUBLE PRECISION P1(3)    ! Point on rotation line
      DOUBLE PRECISION MAG      ! Magnitude of vector
      DOUBLE PRECISION ROT(4,4) ! Rotation matrix
            
      DOUBLE PRECISION ALPHA    ! Angle of move
      DOUBLE PRECISION BETA     ! Angle of move
      
!     MC adaptation variables
     
      INTEGER, PARAMETER :: moveTypes=7 ! Number of different move types 
      DOUBLE PRECISION MCAMP(moveTypes) ! Amplitude of random change      
      INTEGER MCTYPE            ! Type of MC move
      INTEGER winType
      DOUBLE PRECISION DR(3)    ! Displacement for slide move
      INTEGER TEMP
      Double precision WINDOW(moveTypes) ! Size of window for bead selection

! Variables for change of binding state move
      INTEGER AB(NT)            ! Chemical (binding) state
      INTEGER ABP(NT)          ! Underlying (methalation) state
      Double precision d1,d2  !for testing

!     Perform crank-shaft move (MCTYPE 1)
          
      if (MCTYPE.EQ.1) then
         call random_number(urand,rand_stat)
         IP=nint(0.5_dp+urand(1)*NP)
         IB1=nint(0.5_dp+urand(2)*NB)
         if (winType.eq.0) then
             IB2=IB1+nint((urand(3)-0.5_dp)*(2.0_dp*WINDOW(MCTYPE)+1.0))
         elseif (winType.eq.1) then 
             call random_number(urnd,rand_stat)
             IB2=IB1+(2.0_dp*nint(urand(3))-1.0)* &
                     nint(-1.0*log(urnd(1))*WINDOW(MCTYPE))
         endif

         if (IB2.LT.1) then
            IB2=1
         endif
         if (IB2.GT.NB) then
            IB2=NB
         endif

         if (IB2.LT.IB1) then
            TEMP=IB1
            IB1=IB2
            IB2=TEMP
         endif
         IT1=NB*(IP-1)+IB1
         IT2=NB*(IP-1)+IB2
         
         if (IB1.EQ.IB2.AND.IB1.EQ.1) then
            TA(1)=R(IT1+1,1)-R(IT1,1)
            TA(2)=R(IT1+1,2)-R(IT1,2)
            TA(3)=R(IT1+1,3)-R(IT1,3)
         elseif (IB1.EQ.IB2.AND.IB1.EQ.NB) then
            TA(1)=R(IT1,1)-R(IT1-1,1)
            TA(2)=R(IT1,2)-R(IT1-1,2)
            TA(3)=R(IT1,3)-R(IT1-1,3)
         elseif (IB1.EQ.IB2.AND.IB1.NE.1.AND.IB2.NE.NB) then
            TA(1)=R(IT1+1,1)-R(IT1-1,1)
            TA(2)=R(IT1+1,2)-R(IT1-1,2)
            TA(3)=R(IT1+1,3)-R(IT1-1,3)
         else
            TA(1)=R(IT2,1)-R(IT1,1)
            TA(2)=R(IT2,2)-R(IT1,2)
            TA(3)=R(IT2,3)-R(IT1,3)
         endif
         MAG=sqrt(TA(1)**2.+TA(2)**2.+TA(3)**2.)
         TA(1)=TA(1)/MAG
         TA(2)=TA(2)/MAG
         TA(3)=TA(3)/MAG
         P1(1)=R(IT1,1)
         P1(2)=R(IT1,2)
         P1(3)=R(IT1,3)	  
         
         call random_number(urnd,rand_stat)
         ALPHA=MCAMP(1)*(urnd(1)-0.5)
         
         ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
         ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
         ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
		 -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)
         
         ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
         ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
         ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
		 -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)
             
         ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
         ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
         ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
         ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
		 -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)
             
         DO 10 I=IT1,IT2
            RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
            RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
            RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
            UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
            UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
            UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)
            ABP(I)=AB(I)
 10      CONTINUE
        !  ------begining testing---------
        if(.false.) then
            ! This is a code block for testing
            if (abs(RP(IT1,1)-R(IT1,1)).gt.0.000001) then
                print*, "error in crank-shaft move"
                print*, RP(IT1,1), R(IT1,1)
                stop 1
            endif
            if (abs(RP(IT2,1)-R(IT2,1)).gt.0.000001) then
                print*, "error in crank-shaft move"
                print*, RP(IT1,1), R(IT1,1)
                stop 1
            endif
            if(IT1.ne.IT2) then
                d1=(R(IT1+1,1)-R(IT1,1))**2+&
                   (R(IT1+1,2)-R(IT1,2))**2+& 
                   (R(IT1+1,3)-R(IT1,3))**2
                d2=(RP(IT1+1,1)-RP(IT1,1))**2+&
                   (RP(IT1+1,2)-RP(IT1,2))**2+&
                   (RP(IT1+1,3)-RP(IT1,3))**2
                if (abs(d1-d2).gt.0.000001) then
                    print*, "error in crank-shaft move"
                    print*, "distance change in 1"
                    print*, "IT1",IT1," IT2",IT2
                    print*, d1,d2
                    stop 1
                endif
                d1=(R(IT2-1,1)-R(IT2,1))**2+&
                   (R(IT2-1,2)-R(IT2,2))**2+& 
                   (R(IT2-1,3)-R(IT2,3))**2
                d2=(RP(IT2-1,1)-RP(IT2,1))**2+&
                   (RP(IT2-1,2)-RP(IT2,2))**2+&
                   (RP(IT2-1,3)-RP(IT2,3))**2
                if (abs(d1-d2).gt.0.000001) then
                    print*, "error in crank-shaft move"
                    print*, "distance change in 2"
                    print*, d1,d2
                    stop 1
                endif
            endif
        endif
        ! --------end testing--------
        
!     Perform slide move (MCTYPE 2)
         
      elseif (MCTYPE.EQ.2) then
         call random_number(urand,rand_stat)
         IP=nint(0.5_dp+urand(1)*NP)
         IB1=nint(0.5_dp+urand(2)*NB)
         if (winType.eq.0) then
             IB2=IB1+nint((urand(3)-0.5_dp)*(2.0_dp*WINDOW(MCTYPE)+1.0))
         elseif (winType.eq.1) then 
             call random_number(urnd,rand_stat)
             IB2=IB1+(2.0_dp*nint(urand(3))-1.0)* &
                     nint(-1.0*log(urnd(1))*WINDOW(MCTYPE))
         endif
        
         if (IB2.LT.1) then
            IB2=1
         endif
         if (IB2.GT.NB) then
            IB2=NB
         endif

         if (IB2.LT.IB1) then
            TEMP=IB1
            IB1=IB2
            IB2=TEMP
         endif

         IT1=NB*(IP-1)+IB1
         IT2=NB*(IP-1)+IB2
                  
         call random_number(urand,rand_stat)
         DR(1)=MCAMP(2)*(urand(1)-0.5)
         DR(2)=MCAMP(2)*(urand(2)-0.5)
         DR(3)=MCAMP(2)*(urand(3)-0.5)
         
         DO 20 I=IT1,IT2
            RP(I,1)=R(I,1)+DR(1)
            RP(I,2)=R(I,2)+DR(2)
            RP(I,3)=R(I,3)+DR(3)
            UP(I,1)=U(I,1)
            UP(I,2)=U(I,2)
            UP(I,3)=U(I,3)
            ABP(I)=AB(I)
 20      CONTINUE    
             
!     Perform pivot move (MCTYPE 3)
             
      elseif (MCTYPE.EQ.3) then

         call random_number(urnd,rand_stat)
         IP=nint(0.5_dp+urnd(1)*NP)
         if (wintype.eq.0) then
             call random_number(urnd,rand_stat)
             IB1=nint(0.5+urnd(1)*(2.0_dp*WINDOW(MCTYPE)))
             if (IB1.LE.WINDOW(MCTYPE)) then
                IB2=IB1
                if (IB2.GT.NB) then
                   IB2=NB
                endif
                IB1=1
                IT1=NB*(IP-1)+IB1
                IT2=NB*(IP-1)+IB2
                P1(1)=R(IT2,1)
                P1(2)=R(IT2,2)
                P1(3)=R(IT2,3)	  
             else
                IB1=NB+WINDOW(MCTYPE)+1-IB1
                if (IB1.LT.1) then
                   IB1=1
                endif
                IB2=NB
                IT1=NB*(IP-1)+IB1
                IT2=NB*(IP-1)+IB2
                P1(1)=R(IT1,1)
                P1(2)=R(IT1,2)
                P1(3)=R(IT1,3)	  
             endif
         elseif(winType.eq.1) then
             call random_number(urnd,rand_stat)
             if (urnd(1).gt.0.5_dp) then
                call random_number(urnd,rand_stat)
                IB2=nint(-1.0_dp*log(urnd(1))*WINDOW(MCTYPE))+1
                if (IB2.GT.NB) then 
                    IB2=NB
                endif
                IB1=1
                IT1=NB*(IP-1)+IB1
                IT2=NB*(IP-1)+IB2
                P1(1)=R(IT2,1)
                P1(2)=R(IT2,2)
                P1(3)=R(IT2,3)
             else
                call random_number(urnd,rand_stat)
                IB1=NB-nint(-1.0_dp*log(urnd(1))*WINDOW(MCTYPE))
                if (IB1.LT.1) then
                   IB1=1
                endif
                IB2=NB
                IT1=NB*(IP-1)+IB1
                IT2=NB*(IP-1)+IB2
                P1(1)=R(IT1,1)
                P1(2)=R(IT1,2)
                P1(3)=R(IT1,3) 
             endif
         else
             print*, "Error, no other option"
             stop 1
         endif             
             
         call random_number(urand,rand_stat)
         ALPHA=2.*PI*urand(1)
         BETA=acos(2.*urand(2)-1.)
         TA(1)=sin(BETA)*cos(ALPHA)
         TA(2)=sin(BETA)*sin(ALPHA)
         TA(3)=cos(BETA)
             
         ALPHA=MCAMP(3)*(urand(3)-0.5)
             
         ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
         ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
         ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
		 -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)
             
         ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
         ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
         ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
		 -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)
             
         ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
         ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
         ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
         ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
		 -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)
             
         DO 30 I=IT1,IT2
            RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
            RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
            RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
            UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
            UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
            UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)		 
            ABP(I)=AB(I)
 30      CONTINUE

!     Perform rotate move (MCTYPE 4)
!     a.k.a. rotate a single 
             
      elseif (MCTYPE.EQ.4) then
         
         call random_number(urand,rand_stat)
         IP=nint(0.5+urand(1)*NP)
         IB1=nint(0.5+urand(2)*NB)
         IB2=IB1
         IT1=NB*(IP-1)+IB1
         IT2=NB*(IP-1)+IB2
         
         call random_number(urand,rand_stat)
         ALPHA=2.*PI*urand(1)
         BETA=acos(2.*urand(2)-1.)
         TA(1)=sin(BETA)*cos(ALPHA)
         TA(2)=sin(BETA)*sin(ALPHA)
         TA(3)=cos(BETA)
             
         ALPHA=MCAMP(4)*(urand(3)-0.5)
             
         ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
         ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
         ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
		 -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)
             
         ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
         ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
         ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
		 -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)
             
         ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
         ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
         ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
         ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
		 -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)
             
         I=IT1
         UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
         UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
         UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)
         RP(I,1)=R(I,1)
         RP(I,2)=R(I,2)
         RP(I,3)=R(I,3)
         ABP(I)=AB(I)

!     Perform a full chain rotation

      elseif (MCTYPE.EQ.5) then
         
         call random_number(urand,rand_stat)
         IP=nint(0.5_dp+urand(1)*NP)
         IB1=1
         IB2=NB
         IT1=NB*(IP-1)+IB1
         IT2=NB*(IP-1)+IB2
         
         ALPHA=2.0_dp*PI*urand(2)
         BETA=acos(2.0_dp*urand(3)-1.0_dp)
         TA(1)=sin(BETA)*cos(ALPHA)
         TA(2)=sin(BETA)*sin(ALPHA)
         TA(3)=cos(BETA)
         
         call random_number(urnd,rand_stat)
         ALPHA=MCAMP(5)*(urnd(1)-0.5)
             
         ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
         ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
         ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
		 -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)
             
         ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
         ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
         ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
		 -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)
             
         ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
         ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
         ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
         ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
		 -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)
             
         DO 40 I=IT1,IT2
            RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
            RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
            RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
            UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
            UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
            UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)		 
            ABP(I)=AB(I)
 40      CONTINUE

!     Perform full chain slide move (MCTYPE 6)
         
      elseif (MCTYPE.EQ.6) then
         
         call random_number(urnd,rand_stat)
         IP=nint(0.5_dp+urnd(1)*NP)
         IB1=1
         IB2=NB
         IT1=NB*(IP-1)+IB1
         IT2=NB*(IP-1)+IB2
         
         call random_number(urand,rand_stat)
         DR(1)=MCAMP(6)*(urand(1)-0.5_dp)
         DR(2)=MCAMP(6)*(urand(2)-0.5_dp)
         DR(3)=MCAMP(6)*(urand(3)-0.5_dp)
         
         DO 50 I=IT1,IT2
            RP(I,1)=R(I,1)+DR(1)
            RP(I,2)=R(I,2)+DR(2)
            RP(I,3)=R(I,3)+DR(3)
            UP(I,1)=U(I,1)
            UP(I,2)=U(I,2)
            UP(I,3)=U(I,3)
            ABP(I)=AB(I)
 50      CONTINUE    

      elseif (MCTYPE.EQ.7) then  
         ! Change AB (a.k.a HP1 binding type fore section of polymer)
         ! Move amplitude is ignored for this move type

         call random_number(urand,rand_stat)
         IP=nint(0.5_dp+urand(1)*NP)
         IB1=nint(0.5_dp+urand(2)*NB)
         if (winType.eq.0) then
             IB2=IB1+nint((urand(3)-0.5_dp)*(2.0_dp*WINDOW(MCTYPE)+1.0))
         elseif (winType.eq.1) then 
             call random_number(urnd,rand_stat)
             IB2=IB1+(2.0_dp*nint(urand(3))-1.0)* &
                     nint(-1.0*log(urnd(1))*WINDOW(MCTYPE))
         endif

         if (IB2.LT.1) then
            IB2=1
         endif
         if (IB2.GT.NB) then
            IB2=NB
         endif

         if (IB2.LT.IB1) then
            TEMP=IB1
            IB1=IB2
            IB2=TEMP
         endif
         IT1=NB*(IP-1)+IB1
         IT2=NB*(IP-1)+IB2

         !keep binding constant within monomers
         IT1=IT1-MOD(IT1-1,BPM)
         IT2=IT2-MOD(IT2-1,BPM)+BPM-1

         DO J=IT1,IT2
             ABP(J)=1-AB(J)
         ENDDO

         !This loop may not be necessary
         DO I=IT1,IT2
            RP(I,1)=R(I,1)
            RP(I,2)=R(I,2)
            RP(I,3)=R(I,3)
            UP(I,1)=U(I,1)
            UP(I,2)=U(I,2)
            UP(I,3)=U(I,3)
         ENDDO    


      endif
      RETURN      
      END
      
!---------------------------------------------------------------!      
