!-------------------------------------------------------!
!
!       This is Quinn's Program to test things
!
!
      PROGRAM sandBox
      INTEGER NMT  !Total number of monomers
      INTEGER t  !index for loop
      NMT = 100000
 

      OPEN (UNIT=1, FILE = 'input/meth', STATUS = 'OLD')
      DO j=1,NMT/2
          write(1,"(I2)") 0
          write(1,"(I2)") 1
      ENDDO
      CLOSE(1)



      END
