!! ---------------------------------------------------------------*
      
!     
!     This subroutine sets the initial condition for a rodCoil polymer.
!
!     
!     A linker segment is put between each monomer.
!     A=1, B=0, Link=1
!     Cuts off mid monomer at end of chain
!     Based on a code from Andrew Spakowitz written 4-16-04
!     By Quinn 2/16/17
!
      
SUBROUTINE initchemJie(AB,NT,NP,FA,LAM,rand_stat,numA,numB,numLink)

use mersenne_twister
use setPrecision
implicit none
INTEGER, intent(out) :: AB(NT)     ! Chemical identity of beads
INTEGER, intent(in) :: NP          ! Number of polymer chains
INTEGER, intent(in) :: NT          ! Total number of beads

INTEGER I,J,K,IB
real TEST(1)   ! changed to real by Quinn
type(random_stat), intent(inout) ::rand_stat    ! status of random number generator

DOUBLE PRECISION, intent(in) :: FA   ! Fraction of A beads
DOUBLE PRECISION, intent(in) :: LAM  ! Chemical correlation parameter
DOUBLE PRECISION PAA,PBB,PAB,PBA ! Chemical identity statistics


integer, intent(in) :: numA
integer, intent(in) :: numB
integer, intent(in) :: numLink
integer beadsPerPoly
integer ABtype

beadsPerPoly=NT/NP
if (NT.ne.beadsPerPoly*NP) then
    print*, "Error in initchemJi. NT",NT,"NP",NP,"bead/Poly",beadsPerPoly
endif

!		Translate LAM and FA to probabilities

PAA=FA*(1.0_dp-LAM)+LAM
PBB=FA*(LAM-1.0_dp)+1.0_dp
PBA=1.0_dp-PAA
PAB=1.0_dp-PBB

!		Determine the bead identities

IB=1
DO  I=1,NP
    call random_number(TEST,rand_stat)
    if (dble(TEST(1)).lt.FA) then
       ABtype=1
    else
       ABtype=0
    endif
    DO J=1,beadsPerPoly  ! maximum possible monomers is beads
        call random_number(TEST,rand_stat)

        if (ABtype.EQ.1) then
           if (TEST(1).LE.PAA) then
              ABtype=1
           else
              ABtype=0
           endif
        else
           if (TEST(1).LE.PAB) then
              ABtype=1
           else
              ABtype=0
           endif
        endif
         
        if (ABtype.EQ.1) then
            Do K=1,numA
                AB(IB)=1
                IB=IB+1
                if (IB>I*beadsPerPoly) then
                    GO To 10
                endif
            enddo
        else
            Do K=1,numB
                AB(IB)=0
                IB=IB+1
                if (IB>I*beadsPerPoly) then
                    GO To 10
                endif
            enddo
        endif
        DO K=1,numLink
           AB(IB)=1
           IB=IB+1
           if (IB>I*beadsPerPoly) then
               GO To 10
           endif
        enddo
    enddo
10  Continue ! next polymer
enddo 
RETURN     
END
      
!---------------------------------------------------------------*
