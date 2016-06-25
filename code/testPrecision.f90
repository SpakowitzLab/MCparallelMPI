program main
USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
use mersenne_twister

implicit none

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
REAL(dp) :: x
Double precision :: y
Double precision v
real z
INTEGER :: PF   ! input file unit
LOGICAL :: FILEEND=.FALSE. ! done reading file?
CHARACTER*100 :: WORD ! keyword
INTEGER :: NITEMS ! number of items on the line in the parameter file
character*16 fileName

!   variable for random number generator seeding
type(random_stat) rand_stat  ! state of random number chain
integer Irand     ! Seed
character*8 datedum  ! trash
character*10 timedum ! trash
character*5 zonedum  ! trash
integer seedvalues(8) ! clock readings
real urand(1)
integer ii

real A(1)
Double precision :: B,C

if (.false.) then ! set spedific seed
    Irand=7171
else ! seed from clock
    call date_and_time(datedum,timedum,zonedum,seedvalues)
    Irand=-seedvalues(5)*1E7-seedvalues(6)*1E5-seedvalues(7)*1E3-seedvalues(8)
    Irand=mod(Irand,10000)
    print*, "Random Intiger seed:",Irand
endif
call random_setseed(Irand,rand_stat) ! random seed for head node

fileName='../input/params'

x=0.1_dp
y=0.1_dp
z=0.1_dp
v=0.1234567890_dp
print*, 'x',x
print*, 'y',y
print*, 'z',z

print*, 'x*y',x*y
print*, 'x*z',x*z
 
PF=55
OPEN(UNIT=PF,FILE=fileName,STATUS='OLD') 

! read in the keywords one line at a time
DO 
   CALL READLINE(PF,FILEEND,NITEMS)
   IF (FILEEND.and.nitems.eq.0) EXIT

   ! skip empty lines
   IF (NITEMS.EQ.0) CYCLE

   ! Read in the keyword for this line
   CALL READA(WORD,CASESET=1)

   ! Skip any empty lines or any comment lines
   IF (WORD(1:1).EQ.'#') CYCLE

   SELECT CASE(WORD) ! pick which keyword
   CASE('V')
       Call READF(v)
   Endselect
enddo

print*, 'v',v

        
call random_number(urand,rand_stat)
print*, "random real", urand

A=0.51_dp
B=0.51_dp
C=0.51_dp

do ii=1,100000000
    call random_number(A,rand_stat)
    if (A(1).gt.B) B=A(1)
    if (A(1).lt.C) C=A(1)
enddo

print*, "minium found", C
print*, "maximum found", B
end program 
