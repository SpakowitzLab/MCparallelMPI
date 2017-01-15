!---------------------------------------------------------------!
!
!     Appends umbrellaV to file.
!     Appends counts to seperate file
!     First column of each file is savePoint.
!     Quinn started this subroutine on 1/11/17
!
! ---------------------------------------------------  
Subroutine saveUmbrella(mc,md,fileName)
    use simMod
    use setPrecision
    IMPLICIT NONE
    TYPE(MCvar), intent(in) :: mc
    TYPE(MCData), intent(in) :: md
    LOGICAL isfile
    character*16, intent(in) :: fileName
    character*32 fullName
    fullName=  trim(fileName) // trim(mc%repSufix) 
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
        WRITE(1,*), "IND, umbrellaV values ..."
    endif

    WRITE(1,*) mc%IND, md%umbrellaV
    Close(1)
    fullName=  'data/counts' // trim(mc%repSufix) 
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
        WRITE(1,*), "IND, umbrellaV values ..."
    endif

    WRITE(1,*) mc%IND, md%umbrellaCounts
    Close(1)
end subroutine
