SUBROUTINE DENSITYSAL_R(rho,sal,NXm,NYm,NZm)


USE COMDAT
INTEGER(4)   ,INTENT(IN)           :: NXm,NYm,NZm
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)  :: sal
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT) :: rho


rho=0


DO i=0,NXm+1 !i=0,NX+1
DO m=0,NYm+1
DO j=0,NZm+1
!CALL STATECO(rho,sal,NX,NZ) ! this doesn't work
!rho(i,j) = den_w + 804.6*sal(i,j)
rho(i,m,j) = den_w + 0.758871487*sal(i,m,j)


ENDDO
ENDDO
ENDDo



ENDSUBROUTINE DENSITYSAL_R