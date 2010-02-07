SUBROUTINE Total(property, voll, total_property,NXm,NYm,NZm)  

USE COMDAT
INTEGER(4)   ,INTENT(IN)           :: NXm,NYm,NZm
REAL(8),                                  INTENT(OUT) :: total_property
REAL(8), DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1), INTENT(IN)  :: property, voll

total_property = 0
DO i=1,NXm;DO m=1,NYm;DO j=1,NZm
total_property = total_property + property(i,m,j)*voll(i,m,j)
ENDDO;ENDDO;ENDDO
ENDSUBROUTINE
!**************************************************************************





!**************************************************************************
SUBROUTINE ACC_1
USE COMDAT
  dt  = dt_o 
WRITE(*,*)
WRITE(*,*) '*************************************************'
WRITE(*,'(A,I5,A,f9.4,A,f9.4)') ' STEP =', step*subiteration, '  time_begin =', time_sec, '  dt =', dt

WRITE(18,*)
WRITE(18,*) '*************************************************'
WRITE(18,'(A,I5,A,f9.4,A,f9.4)') ' STEP =', step*subiteration, '  time_begin =', time_sec, '  dt =', dt

ENDSUBROUTINE

!**************************************************************************
SUBROUTINE ACC_2
USE COMDAT
  dt  = dt_o 
WRITE(*,*)
WRITE(*,*) '*************************************************'
WRITE(*,'(A,I5,A,f9.4,A,f9.4)') ' STEP =', step, '  time_begin =', time_sec, '  dt =', dt

WRITE(18,*)
WRITE(18,*) '*************************************************'
WRITE(18,'(A,I5,A,f9.4,A,f9.4)') ' STEP =', step, '  time_begin =', time_sec, '  dt =', dt

ENDSUBROUTINE
!**************************************************************************

