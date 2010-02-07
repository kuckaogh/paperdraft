SUBROUTINE Total(property, voll, total_property)

USE COMDAT
REAL(8),                                  INTENT(OUT) :: total_property
REAL(8), DIMENSION(0:NX+1,0:NY+1,0:NZ+1), INTENT(IN)  :: property, voll

total_property = 0
total_property = SUM( property(1:NX,1:NY,1:NZ)*voll(1:NX,1:NY,1:NZ) )

ENDSUBROUTINE
!**************************************************************************





!**************************************************************************
SUBROUTINE ACC_1
USE COMDAT

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

WRITE(*,*)
WRITE(*,*) '*************************************************'
WRITE(*,'(A,I5,A,f9.4,A,f9.4)') ' STEP =', step, '  time_begin =', time_sec, '  dt =', dt

WRITE(18,*)
WRITE(18,*) '*************************************************'
WRITE(18,'(A,I5,A,f9.4,A,f9.4)') ' STEP =', step, '  time_begin =', time_sec, '  dt =', dt

ENDSUBROUTINE
!**************************************************************************

SUBROUTINE ZERO_Vel(uu,vv,ww)
!Assign zero for uncomputed grid

USE COMDAT

REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1), INTENT(INOUT)    :: uu
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1), INTENT(INOUT)    :: vv
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ), INTENT(INOUT)    :: ww


 
DO i=0,NX; DO m=0,NY+1; DO j= NZ+1, NZ+1
 uu(i,m,j) = 0
ENDDO; ENDDO; ENDDO

DO i=0,NX; DO m=0,NY; DO j= NZ+1, NZ+1
 vv(i,m,j) = 0
ENDDO; ENDDO; ENDDO

DO i=0,NX; DO m=0,NY+1; DO j=NZ+1, NZ  ! 05/02/15
 ww(i,m,j) = 0
ENDDO; ENDDO; ENDDO

ENDSUBROUTINE

!**************************************************************************
