!---------------------------------------------------------------*
!
!           Makes Monti Carlo Moves
!           
!    Quinn Made Changes to this file starting on 12/15/15
!       
!


! Find change in bead position for a crank-shaft type move

SUBROUTINE MC_move(R,U,RP,UP,NT,NB,NP,IP,IB1,IB2,IT1,IT2,MCTYPE &
                  ,MCAMP,WINDOW,AB,ABP,BPM,rand_stat,winType &
                  ,IT3,IT4)

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
INTEGER IP2               ! Second Test polymer if applicable 
INTEGER IB1               ! Test bead position 1
INTEGER IT1               ! Index of test bead 1
INTEGER IB2               ! Test bead position 2
INTEGER IT2               ! Index of test bead 2
INTEGER IT3               ! Test bead position 3 if applicable
INTEGER IT4               ! Test bead position 4 if applicable

INTEGER I,J  ! Test indices
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

! variables for reptation move
!double precision dR(3) ! old displacement vector
double precision RperpOld(3)  ! portion of R in perpendicular direction old
double precision perpNew(3)  ! unit vector in new perp direction
double precision dRdotU
double precision RparaMag
double precision RperpMag
double precision upara  ! parallel component of other u 
double precision uperp  ! perpendicular component of other u
double precision utwist ! the rest of u
double precision Uvec(3) 
double precision twistDir(3) ! unit vector perp to upara and uperp
double precision dtemp  ! temparary storage of double
double precision vperp(3) ! unit vector allignled with rperpOld
double precision vperpdotU ! small round error correction

!     Perform crank-shaft move (MCTYPE 1)
    
if (MCTYPE.EQ.1) then
   call random_number(urand,rand_stat)
   IP=ceiling(urand(1)*NP)
   IB1=ceiling(urand(2)*NB)
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
       
   DO I=IT1,IT2
      RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
      RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
      RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
      UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
      UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
      UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)
      ABP(I)=AB(I)
  enddo
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
   IP=ceiling(urand(1)*NP)
   IB1=ceiling(urand(2)*NB)
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
   
   DO I=IT1,IT2
      RP(I,1)=R(I,1)+DR(1)
      RP(I,2)=R(I,2)+DR(2)
      RP(I,3)=R(I,3)+DR(3)
      UP(I,1)=U(I,1)
      UP(I,2)=U(I,2)
      UP(I,3)=U(I,3)
      ABP(I)=AB(I)
   enddo 
       
!     Perform pivot move (MCTYPE 3)
       
elseif (MCTYPE.EQ.3) then

   call random_number(urnd,rand_stat)
   IP=ceiling(urnd(1)*NP)
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
       
   DO I=IT1,IT2
      RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
      RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
      RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
      UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
      UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
      UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)
      ABP(I)=AB(I)
   enddo

!     Perform rotate move (MCTYPE 4)
!     a.k.a. rotate a single        
elseif (MCTYPE.EQ.4) then
   
   call random_number(urand,rand_stat)
   IP=ceiling(urand(1)*NP)
   IB1=ceiling(urand(2)*NB)
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
   IP=ceiling(urand(1)*NP)
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
       
   DO I=IT1,IT2
      RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
      RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
      RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
      UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
      UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
      UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3) 
      ABP(I)=AB(I)
   enddo

!     Perform full chain slide move (MCTYPE 6)
elseif (MCTYPE.EQ.6) then
   
   call random_number(urnd,rand_stat)
   IP=ceiling(urnd(1)*NP)
   IB1=1
   IB2=NB
   IT1=NB*(IP-1)+IB1
   IT2=NB*(IP-1)+IB2
   
   call random_number(urand,rand_stat)
   DR(1)=MCAMP(6)*(urand(1)-0.5_dp)
   DR(2)=MCAMP(6)*(urand(2)-0.5_dp)
   DR(3)=MCAMP(6)*(urand(3)-0.5_dp)
   
   DO I=IT1,IT2
      RP(I,1)=R(I,1)+DR(1)
      RP(I,2)=R(I,2)+DR(2)
      RP(I,3)=R(I,3)+DR(3)
      UP(I,1)=U(I,1)
      UP(I,2)=U(I,2)
      UP(I,3)=U(I,3)
      ABP(I)=AB(I)
   enddo

elseif (MCTYPE.EQ.7) then  
   ! Change AB (a.k.a HP1 binding type fore section of polymer)
   ! Move amplitude is ignored for this move type

   call random_number(urand,rand_stat)
   IP=ceiling(urand(1)*NP)
   IB1=ceiling(urand(2)*NB)
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

! chain flip move
elseif (MCTYPE.EQ.8) then
   call random_number(urand,rand_stat)
   IP=ceiling(urand(1)*NP)
   IB1=1
   IB2=NB
   IT1=NB*(IP-1)+IB1
   IT2=NB*(IP-1)+IB2
   DO I=0,NB-1
      RP(IT1+I,1)=R(IT2-I,1)
      RP(IT1+I,2)=R(IT2-I,2)
      RP(IT1+I,3)=R(IT2-I,3)
      UP(IT1+I,1)=-U(IT2-I,1)
      UP(IT1+I,2)=-U(IT2-I,2)
      UP(IT1+I,3)=-U(IT2-I,3)
      ABP(IT1+I)=AB(IT1+I)
   ENDDO    
! switch two chains
elseif(MCTYPE.EQ.9) then
   call random_number(urnd,rand_stat)
   IP=ceiling(urnd(1)*NP)
   call random_number(urnd,rand_stat)
   IP2=ceiling(urnd(1)*NP)
   ! Don't switch a chain with itself
   if (IP.eq.IP2) then
       IP2=IP-1
       if (IP2.eq.0) then
           IP2=2
       endif
   endif
   IT1=NB*(IP-1)+1
   IT2=NB*(IP-1)+NB
   IT3=NB*(IP2-1)+1
   IT4=NB*(IP2-1)+NB
   DO I=0,NB-1
      RP(IT1+I,1)=R(IT3+I,1)
      RP(IT1+I,2)=R(IT3+I,2)
      RP(IT1+I,3)=R(IT3+I,3)
      UP(IT1+I,1)=U(IT3+I,1)
      UP(IT1+I,2)=U(IT3+I,2)
      UP(IT1+I,3)=U(IT3+I,3)
      ABP(IT1+I)=AB(IT1+I)
      RP(IT3+I,1)=R(IT1+I,1)
      RP(IT3+I,2)=R(IT1+I,2)
      RP(IT3+I,3)=R(IT1+I,3)
      UP(IT3+I,1)=U(IT1+I,1)
      UP(IT3+I,2)=U(IT1+I,2)
      UP(IT3+I,3)=U(IT1+I,3)
      ABP(IT3+I)=AB(IT3+I)
   ENDDO
   IB1=0; IB1=IB1/IB1; IB2=IB1

! single bead reptation
elseif(MCTYPE.EQ.10) then
    call random_number(urnd,rand_stat)
    IP=ceiling(urnd(1)*NP)
    IT1=NB*(IP-1)+1
    IT2=NB*(IP-1)+NB
    ! move forward or backward
    call random_number(urnd,rand_stat)
    if (urnd(1).lt.0.5_dp) then
        dR(1)=R(IT1+1,1)-R(IT1,1)
        dR(2)=R(IT1+1,2)-R(IT1,2)
        dR(3)=R(IT1+1,3)-R(IT1,3)
        RparaMag=(dR(1)*u(IT1,1) + dR(2)*u(IT1,2) + dR(3)*u(IT1,3))
        dRdotU=dR(1)*U(IT1,1)+dR(2)*U(IT1,2)+dR(3)*U(IT1,3)
        RperpOld(1)=dR(1)-dRdotU*U(IT1,1) 
        RperpOld(2)=dR(2)-dRdotU*U(IT1,2) 
        RperpOld(3)=dR(3)-dRdotU*U(IT1,3) 
        RperpMag=sqrt(RperpOld(1)**2+RperpOld(2)**2+RperpOld(3)**2)
        upara=u(IT1,1)*u(IT1+1,1)+&
              u(IT1,2)*u(IT1+1,2)+&
              u(IT1,3)*u(IT1+1,3)
        ! Need to handle rownding error correctly
        ! generate a vector, vperp, that is perp to u
        if((RperpMag/RparaMag).lt.0.000000001_dp) then
            uperp=0.0_dp
        else
            vperp(1)=RperpOld(1)/RperpMag
            vperp(2)=RperpOld(2)/RperpMag
            vperp(3)=RperpOld(3)/RperpMag
            vperpdotU=vperp(1)*U(IT1,1)+vperp(2)*U(IT1,2)+vperp(3)*U(IT1,3)
            vperp(1)=vperp(1)-U(IT1,1)*vperpdotU
            vperp(2)=vperp(2)-U(IT1,2)*vperpdotU
            vperp(3)=vperp(3)-U(IT1,3)*vperpdotU
            mag=sqrt(vperp(1)**2+vperp(2)**2+vperp(3)**2)
            vperp(1)=vperp(1)/mag
            vperp(2)=vperp(2)/mag
            vperp(3)=vperp(3)/mag
            uperp=(vperp(1)*u(IT1+1,1)+&
                   vperp(2)*u(IT1+1,2)+&
                   vperp(3)*u(IT1+1,3))

        endif
        dtemp=1.0_dp-upara**2-uperp**2
        if (dtemp.lt.0.0_dp) then
            if (dtemp.gt.-0.0000001_dp) then
                utwist=0.0_dp
            else
                print*, "Error in MC_move. Imag Num. Forward"
                stop 1
            endif
        else
            utwist=sqrt(dtemp)
        endif
        
        !use Uvec for final bead orientation
        Uvec(1)=U(IT2,1)
        Uvec(2)=U(IT2,2)
        Uvec(3)=U(IT2,3)
        call random_perp(Uvec,perpNew,twistDir,rand_stat)
        ! perpNew is a unit vector
        UP(IT2,1)=U(IT2,1)*upara+perpNew(1)*uperp+twistDir(1)*utwist
        UP(IT2,2)=U(IT2,2)*upara+perpNew(2)*uperp+twistDir(2)*utwist
        UP(IT2,3)=U(IT2,3)*upara+perpNew(3)*uperp+twistDir(3)*utwist
        mag=sqrt(UP(IT2,1)**2+UP(IT2,2)**2+UP(IT2,3)**2)
        UP(IT2,1)=UP(IT2,1)/mag
        UP(IT2,2)=UP(IT2,2)/mag
        UP(IT2,3)=UP(IT2,3)/mag
        
        RP(IT2,1)=R(IT2,1)+U(IT2,1)*RparaMag+perpNew(1)*RperpMag
        RP(IT2,2)=R(IT2,2)+U(IT2,2)*RparaMag+perpNew(2)*RperpMag
        RP(IT2,3)=R(IT2,3)+U(IT2,3)*RparaMag+perpNew(3)*RperpMag
        
        DO I=IT1,IT2-1
           RP(I,1)=R(I+1,1)
           RP(I,2)=R(I+1,2)
           RP(I,3)=R(I+1,3)
           UP(I,1)=U(I+1,1)
           UP(I,2)=U(I+1,2)
           UP(I,3)=U(I+1,3)
        enddo

        !call test_equiv_forward(U,R,UP,RP,NT,IT1,IT2,RparaMag,RperpMag)

    else
        dR(1)=R(IT2,1)-R(IT2-1,1)
        dR(2)=R(IT2,2)-R(IT2-1,2)
        dR(3)=R(IT2,3)-R(IT2-1,3)
        RparaMag=(dR(1)*u(IT2,1) + dR(2)*u(IT2,2) + dR(3)*u(IT2,3))
        dRdotU=dR(1)*U(IT2,1)+dR(2)*U(IT2,2)+dR(3)*U(IT2,3)
        RperpOld(1)=dR(1)-dRdotU*U(IT2,1) 
        RperpOld(2)=dR(2)-dRdotU*U(IT2,2) 
        RperpOld(3)=dR(3)-dRdotU*U(IT2,3) 
        RperpMag=sqrt(RperpOld(1)**2+RperpOld(2)**2+RperpOld(3)**2)
        upara=u(IT2,1)*u(IT2-1,1)+&
              u(IT2,2)*u(IT2-1,2)+&
              u(IT2,3)*u(IT2-1,3)
        ! Need to handle rownding error correctly
        ! generate a vector, vperp, that is perp to u
        if((RperpMag/RparaMag).lt.0.000000001_dp) then
            uperp=0.0_dp
        else
            vperp(1)=RperpOld(1)/RperpMag
            vperp(2)=RperpOld(2)/RperpMag
            vperp(3)=RperpOld(3)/RperpMag
            vperpdotU=vperp(1)*U(IT2,1)+vperp(2)*U(IT2,2)+vperp(3)*U(IT2,3)
            vperp(1)=vperp(1)-U(IT2,1)*vperpdotU
            vperp(2)=vperp(2)-U(IT2,2)*vperpdotU
            vperp(3)=vperp(3)-U(IT2,3)*vperpdotU
            mag=sqrt(vperp(1)**2+vperp(2)**2+vperp(3)**2)
            vperp(1)=vperp(1)/mag
            vperp(2)=vperp(2)/mag
            vperp(3)=vperp(3)/mag
            uperp=(vperp(1)*u(IT2-1,1)+&
                   vperp(2)*u(IT2-1,2)+&
                   vperp(3)*u(IT2-1,3))

        endif
        dtemp=1.0_dp-upara**2-uperp**2
        if (dtemp.lt.0.0_dp) then
            if (dtemp.gt.-0.0000001_dp) then
                utwist=0.0_dp
            else
                print*, "Error in MC_move. Imag Num. Backward"
                stop 1
            endif
        else
            utwist=sqrt(dtemp)
        endif
        
        !use Uvec for final bead orientation
        Uvec(1)=U(IT1,1)
        Uvec(2)=U(IT1,2)
        Uvec(3)=U(IT1,3)

        call random_perp(Uvec,perpNew,twistDir,rand_stat)
        ! perpNew is a unit vector
        UP(IT1,1)=U(IT1,1)*upara+perpNew(1)*uperp+twistDir(1)*utwist
        UP(IT1,2)=U(IT1,2)*upara+perpNew(2)*uperp+twistDir(2)*utwist
        UP(IT1,3)=U(IT1,3)*upara+perpNew(3)*uperp+twistDir(3)*utwist
        mag=sqrt(UP(IT1,1)**2+UP(IT1,2)**2+UP(IT1,3)**2)
        UP(IT1,1)=UP(IT1,1)/mag
        UP(IT1,2)=UP(IT1,2)/mag
        UP(IT1,3)=UP(IT1,3)/mag

        RP(IT1,1)=R(IT1,1)-U(IT1,1)*RparaMag-perpNew(1)*RperpMag
        RP(IT1,2)=R(IT1,2)-U(IT1,2)*RparaMag-perpNew(2)*RperpMag
        RP(IT1,3)=R(IT1,3)-U(IT1,3)*RparaMag-perpNew(3)*RperpMag
        DO I=IT1+1,IT2
           RP(I,1)=R(I-1,1)
           RP(I,2)=R(I-1,2)
           RP(I,3)=R(I-1,3)
           UP(I,1)=U(I-1,1)
           UP(I,2)=U(I-1,2)
           UP(I,3)=U(I-1,3)
        enddo
    endif
    ! Testing
   ! if (abs(mag-1.0_dp).gt.0.000001_dp) then
   !     print*, "Error in MC_move. Bad U."
   !     stop 1
   ! endif
   ! if (ISNAN(UP(IT1,1))) then
   !     print*, "Error in MC_move, NAN encountered. UP(IT1,1)"
   !     print*, "upara",upara," uperp",uperp," utwist",utwist
   !     stop 1
   ! endif
   ! if (ISNAN(RP(IT1,1))) then
   !     print*, "Error in MC_move, NAN encountered. 1"
   !     print*, "dR", dR
   !     if (urnd(1).lt.0.5_dp) then
   !         print*, "U(IT1)", U(IT1,1),U(IT1,2),U(IT1,3)
   !         print*, "|U(IT1)|", abs(U(IT1,1)**2+U(IT1,2)**2+U(IT1,3)**2)
   !     else
   !         print*, "U(IT2)", U(IT2,1),U(IT2,2),U(IT2,3)
   !         print*, "|U(IT2)|", abs(U(IT2,1)**2+U(IT2,2)**2+U(IT2,3)**2)
   !     endif
   !     print*, "RparaMag",RparaMag,"RperpMag",RperpMag
   !     stop 1
   ! endif
   ! if (ISNAN(RP(IT2,1))) then
   !     print*, "Error in MC_move, NAN encountered. 2"
   !     stop 1
   ! endif
   ! if (ISNAN(R(IT1,1))) then
   !     print*, "Error in MC_move, NAN encountered. 3"
   !     stop 1
   ! endif
   ! if (ISNAN(R(IT2,1))) then
   !     print*, "Error in MC_move, NAN encountered. 4"
   !     stop 1
   ! endif
   ! if (abs(Uvec(1)**2+Uvec(2)**2+Uvec(3)**2-1.0_dp).gt.0.00001_dp) then
   !    print*, "Error in MC_move: Uvec not a unit vector"
   !    print*, Uvec(1)**2+Uvec(2)**2+Uvec(3)**2
   !    stop 1
   ! endif 
   ! if (abs(dR(1)**2+dR(2)**2+dR(3)**2 &
   !         -RparaMag**2-RperpMag**2).gt.0.0001_dp) then
   !     print*, "dR", dR
   !     print*, "RparaMag", RparaMag," RperpMag",RperpMag
   !     print*, dR(1)**2+dR(2)**2+dR(3)**2, " .ne. ",RparaMag**2+RperpMag**2
   !     print*, "urnd", urnd(1)
   !     print*, "RperpOld",RperpOld
   !     if (urnd(1).lt.0.5_dp) then
   !         print*, "U(IT1)", U(IT1,1),U(IT1,2),U(IT1,3)
   !         print*, "|U(IT1)|", abs(U(IT1,1)**2+U(IT1,2)**2+U(IT1,3)**2)
   !     else
   !         print*, "U(IT2)", U(IT2,1),U(IT2,2),U(IT2,3)
   !         print*, "|U(IT2)|", abs(U(IT2,1)**2+U(IT2,2)**2+U(IT2,3)**2)
   !     endif
   !     print*, "Error in MC_move, triangle mismatch"
   !     stop 1
   ! endif
   ! if (abs(upara**2+uperp**2 +utwist**2 - 1.0_dp).gt.0.00001_dp) then
   !     print*, "Error in MC_move, u vector incorrect"
   !     stop 1
   ! endif
    ! END Testing
    DO J=IT1,IT2
        ABP(J)=1-AB(J)
    ENDDO
endif


RETURN      
END
subroutine test_equiv_forward(U,R,UP,RP,NT,IT1,IT2,RparaMag,RperpMag)
use setPrecision
IMPLICIT NONE
! inputs
INTEGER NT,IT1,IT2
DOUBLE PRECISION R(NT,3)  ! Bead positions
DOUBLE PRECISION U(NT,3)  ! Tangent vectors
DOUBLE PRECISION RP(NT,3)  ! Bead positions
DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
DOUBLE PRECISION RparaMag, RperpMag

!defined
double precision drOld(3)
double precision drNew(3)
double precision drParOld, drParNew
double precision drPerpOld(3)
double precision drPerpNew(3)
double precision Eta
double precision GIOld(3)
double precision GINew(3)
Eta=1.89756278_dp

drOld(1)=R(IT1+1,1)-R(IT1,1)
drOld(2)=R(IT1+1,2)-R(IT1,2)
drOld(3)=R(IT1+1,3)-R(IT1,3)
DRPAROld=DROld(1)*U(IT1,1)+DROld(2)*U(IT1,2)+DROld(3)*U(IT1,3)
drNew(1)=RP(IT2,1)-RP(IT2-1,1)
drNew(2)=RP(IT2,2)-RP(IT2-1,2)
drNew(3)=RP(IT2,3)-RP(IT2-1,3)
DRPARNew=DRNew(1)*UP(IT2-1,1)+&
         DRNew(2)*UP(IT2-1,2)+&
         DRNew(3)*UP(IT2-1,3)
if (abs(drOld(1)**2+drOld(2)**2+drOld(3)**2&
      -(drNew(1)**2+drNew(2)**2+drNew(3)**2)).gt.0.000001) then
      print*, "drOld",drOld, " mag^2=",drOld(1)**2+drOld(2)**2+drOld(3)**2
      print*, "drNew",drNew, " mag^2=",drNew(1)**2+drNew(2)**2+drNew(3)**2
      print*, "Difference detected in test_equiv, 0"
      stop 1
endif

if (abs(drParOld-drParNew).gt.0.0000001_dp) then
    print*, "DRParOld",DRParOld,"DRParNew",DRParNew
    print*, "Difference detected in test_equiv, 1"
    stop 1
endif

drPerpOld(1)=drOld(1)-drParOld*U(IT1,1)
drPerpOld(2)=drOld(2)-drParOld*U(IT1,2)
drPerpOld(3)=drOld(3)-drParOld*U(IT1,3)
drPerpNew(1)=drNew(1)-drParNew*UP(IT2-1,1)
drPerpNew(2)=drNew(2)-drParNew*UP(IT2-1,2)
drPerpNew(3)=drNew(3)-drParNew*UP(IT2-1,3)

if (abs(drPerpOld(1)**2+drPerpOld(2)**2+drPerpOld(3)**2 &
      -(drPerpNew(1)**2+drPerpNew(2)**2+drPerpNew(3)**2)).gt.0.000001_dp) then
  print*, "drOld",sqrt(drOld(1)**2+drOld(2)**2+drOld(3)**2)
  print*, "drNew",sqrt(drNew(1)**2+drNew(2)**2+drNew(3)**2)
  print*, "dRparOld",dRparOld,"dRparNew",drParNew
  print*, "perp Old:", drPerpOld(1)**2+drPerpOld(2)**2+drPerpOld(3)**2
  print*, "perp New:", drPerpNew(1)**2+drPerpNew(2)**2+drPerpNew(3)**2
  print*, "RparaMag",RparaMag,"RperpMag",RperpMag
  print*, "Difference detected in test_equiv, 2"
  stop 1
endif

GIOld(1)=U(IT1+1,1)-U(IT1,1)-Eta*dRperpOld(1)
GIOld(2)=U(IT1+1,2)-U(IT1,2)-Eta*dRperpOld(2)
GIOld(3)=U(IT1+1,3)-U(IT1,3)-Eta*dRperpOld(3)
GINew(1)=UP(IT2,1)-UP(IT2-1,1)-Eta*dRperpNew(1)
GINew(2)=UP(IT2,2)-UP(IT2-1,2)-Eta*dRperpNew(2)
GINew(3)=UP(IT2,3)-UP(IT2-1,3)-Eta*dRperpNew(3)

if (abs(GIOld(1)**2+GIOld(2)**2+GIOld(3)**2&
      -(GINew(1)**2+GINew(2)**2+GINew(3)**2)).gt.0.000001_dp) then
  print*, "Difference detected in test_equiv, 3"
  stop 1
endif

return
end subroutine
subroutine random_perp(u,p,t,rand_stat)
use mersenne_twister      
use setPrecision
IMPLICIT NONE
DOUBLE PRECISION, PARAMETER :: PI=3.141592654 ! Value of pi 
type(random_stat) rand_stat  ! status of random number generator      
real urnd(1) ! single random number

double precision v(2) ! random 2-vec
double precision u(3) ! input
double precision p(3) ! output: random perpendicular to u
double precision t(3) ! orthogonal to p and u
double precision f

if (abs(u(1)**2+u(2)**2+u(3)**2-1.0_dp) .gt. 0.0000001_dp) then
    print*, u
    print*, "Error in random_perp, please give me a unit vector"
    stop 1
endif

call random_number(urnd,rand_stat)
v(1)=cos(2*PI*urnd(1))
v(2)=sin(2*PI*urnd(1))

if (u(3).gt.0.0) then
    f=1.0_dp/(1+u(3))
    p(1)=(u(3)+f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
    p(2)=(u(3)+f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
    p(3)=-1.0_dp*(u(2)*v(2)+u(1)*v(1))
else
    f=1.0_dp/(1-u(3))
    p(1)=(-u(3)+f*u(2)**2)*v(1) - u(2)*u(1)*v(2)*f
    p(2)=(-u(3)+f*u(1)**2)*v(2) - u(2)*u(1)*v(1)*f
    p(3)=(u(2)*v(2)+u(1)*v(1))

endif

t(1)=u(2)*p(3)-u(3)*p(2)
t(2)=u(3)*p(1)-u(1)*p(3)
t(3)=u(1)*p(2)-u(2)*p(1)

! random sign
call random_number(urnd,rand_stat)
if (urnd(1).lt.0.5_dp) then
    t(1)=-1.0_dp*t(1)
    t(2)=-1.0_dp*t(2)
    t(3)=-1.0_dp*t(3)
endif

! Testing
!if (abs(p(1)*u(1)+p(2)*u(2)+p(3)*u(3)).gt.0.000001_dp) then
!    print*, "Error in random_perp, 1"
!    stop 1
!endif
!if (abs(p(1)**2+p(2)**2+p(3)**2-1) .gt. 0.0000001_dp) then
!    print*, "Error in random_perp, 2"
!    stop 1
!endif
!if (abs(t(1)**2 + t(2)**2 + t(3)**2 -1).gt.0.000001_dp) then
!    print*, "Error in random_perp, 3"
!    stop 1
!endif
!if (abs(t(1)*p(1)+t(2)*p(2)+t(3)*p(3)).gt.0.0000001_dp) then
!    print*, "Error in random_perp, 4"
!    stop 1
!endif
!if (abs(t(1)*u(1)+t(2)*u(2)+t(3)*u(3)).gt.0.0000001_dp) then
!    print*, "Error in random_perp, 5"
!    stop 1
!endif
! END Testing

return
end subroutine
!---------------------------------------------------------------!      
