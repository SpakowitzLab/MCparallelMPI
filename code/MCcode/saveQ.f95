!---------------------------------------------------------------!
!
!     Appends umbrellaV to file.
!     Appends counts to seperate file
!     First column of each file is savePoint.
!     Quinn started this subroutine on 1/11/17
!
! ---------------------------------------------------  
Subroutine saveQ(mc)
    use simMod
    use setPrecision
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    LOGICAL isfile
    character*16 fileName
    character*32 fullName
    fileName='data/rxnQ'
    fullName=  trim(fileName) // trim(mc%repSufix) 
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
        WRITE(1,*), "IND, Q, indUmbrella"
    endif

    WRITE(1,*) mc%IND, mc%rxnQ, mc%IndUmbrella
end subroutine
