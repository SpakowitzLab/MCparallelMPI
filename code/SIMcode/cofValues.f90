Subroutine PT_cofValues(cof,nPTReplicas)
    use setPrecision
    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
    Implicit none
    Integer nPTReplicas 
    Double precision cof(nptReplicas)
    INteger rep
    ! IO variables
    character*16 fileName  ! file with parameters
    INTEGER :: PF   ! input file unit
    LOGICAL :: FILEEND=.FALSE. ! done reading file?
    CHARACTER*100 :: WORD ! keyword
    INTEGER :: NITEMS ! number of items on the line in the parameter file
    ! outputs
    Double precision gap
    Double precision minCof
    ! -----------------------
    !
    !  Read from file
    !
    !-------------------------
    fileName='input/RepSetting'
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
       CASE('MIN')
           Call READF(minCof)
       CASE('GAP')
           Call READF(gap)
       CASE DEFAULT
           print*, "Error in MCvar_setParams.  Unidentified keyword:", &
                   TRIM(WORD)
           stop 1
       ENDSELECT
    ENDDO
    print*, "MIN=",minCof," GAP=",GAP
    do rep=1,nPTReplicas
        !cof(rep)=2.0_dp-rep*0.08_dp  !over mu values
        cof(rep)=minCof+Gap*(rep-1)
    enddo
    print*, "cof values:"
    print*, cof
    Close(PF)
    return
end subroutine
