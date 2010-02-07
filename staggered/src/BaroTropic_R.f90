SUBROUTINE BaroTropic_R(BaT,                       &   !OUT                                 
                      volsum_1)                  !IN                


USE COMDAT


REAL(8)    ,DIMENSION(0:NX+1,0:NY+1) ,INTENT(IN)    :: volsum_1
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)                :: volsum
!INTEGER(4) ,DIMENSION(0:NX+1,0:NY+1) ,INTENT(IN)    :: topcell

REAL(8)    ,DIMENSION(1:NX,1:NY,2)   ,INTENT(OUT)   :: BaT


BaT  = 0
volsum=(1/dz)*volsum_1

!******************************************
DO m=1,NY
DO i=1,NX-1 !i=1,NX-1

 BaT(i,m,1)= gravity*(volsum(i+1,m)-volsum(i,m))*dz/dx
 !BaT(i,m,1)= gravity*(volsum(i+1,m)-volsum(i,m))*dz*dy/dx

ENDDO
ENDDO
DO m=1,NY-1
DO i=1,NX !i=1,NX

 BaT(i,m,2)= gravity*(volsum(i,m+1)-volsum(i,m))*dz/dy
! BaT(i,m,2)= gravity*(volsum(i,m+1)-volsum(i,m))*dz*dx/dy

ENDDO
ENDDO
!******************************************

ENDSUBROUTINE BaroTropic_R

