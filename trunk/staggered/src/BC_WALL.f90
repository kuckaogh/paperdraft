SUBROUTINE BC_wall(uu,vv,ww)          
!___________________________________________

USE COMDAT


!INTEGER(4) ,DIMENSION(0:NX+1,0:NY+1)        ,INTENT(IN)    :: topcell

REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(INOUT) :: uu
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(INOUT) :: vv
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(INOUT) :: ww

!_____________________________________________________________


!Assign zero for no-flux wall
!-------------------------------
DO m=0,NY+1
DO j=0,NZ+1 
 uu( 0,m,j) = 0
 uu(NX,m,j) = 0 
ENDDO
ENDDO


DO i=0,NX+1
DO j=0,NZ+1 
 vv(i, 0,j) = 0
 vv(i,NY,j) = 0
ENDDO
ENDDO


DO m=0,NY+1
DO i=0,NX+1  
 ww(i,m, 0) = 0
! ww(i,m,NZ) = 0
ENDDO
ENDDO


ENDSUBROUTINE BC_wall