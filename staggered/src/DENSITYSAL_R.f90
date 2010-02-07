SUBROUTINE DENSITYSAL_R(rho,sal)


USE COMDAT
REAL(8)   ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(IN)  :: sal
REAL(8)   ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(OUT) :: rho


rho=0


DO i=1-1,NX+1 !i=0,NX+1
DO m=0,NY+1
DO j=0,NZ+1
!CALL STATECO(rho,sal,NX,NZ) ! this doesn't work
!rho(i,j) = den_w + 804.6*sal(i,j)
rho(i,m,j) = den_w + 0.758871487*sal(i,m,j)


ENDDO
ENDDO
ENDDo



ENDSUBROUTINE DENSITYSAL_R