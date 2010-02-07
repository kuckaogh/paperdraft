SUBROUTINE POISSON_Revised(pp2,                    & !OUT
                           Matrix,topcell,kp)        !IN

USE COMDAT

INTEGER(4) ,DIMENSION(0:NX+1,0:NY+1)        ,INTENT(IN)    :: topcell
REAL(8)    ,DIMENSION(1:NX  ,1:NY  ,1:NZ  ) ,INTENT(IN)    :: Matrix
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(INOUT) :: pp2
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: kp

dx2=dx*dx
dy2=dy*dy
dz2=dz*dz

DO j=0,NZ+1
DO m=0,NY+1
 pp2(0,m,j)=0
 pp2(NX+1,m,j)=0
ENDDO
ENDDO

DO i=0,NX+1
DO j=0,NZ+1
 pp2(i,0,j)=0
 pp2(i,NY+1,j)=0
ENDDO
ENDDO


DO i=0,NX+1 
DO m=0,NY+1
 pp2(i,m,0)=0
 pp2(i,m,NZ)=0   ! hydrostatic for the top layer
 pp2(i,m,NZ+1)=0
ENDDO
ENDDO

IF (NY==1) THEN

DO iSORPOS= 1, ii_SOR_Pos

 DO i=NX1B,NX1E !i=1,NX
 DO j=1,topcell(i,1) -1


   pp2(i,1,j)=1./kp(i,1,j)*  beta*( pp2(i+1,1,j)/dx2+pp2(i-1,1,j)/dx2   &
                                  + pp2(i,1,j+1)/dz2+pp2(i,1,j-1)/dz2   &
                                  - matrix(i,1,j) )  &           !- den_w/dt*matrix(i,m,j) )  &
                       + (1.-beta)*pp2(i,1,j)
 
 ENDDO
 ENDDO

ENDDO

ELSE
 
DO iSORPOS= 1, ii_SOR_Pos

 DO i=NX1B,NX1E !i=1,NX
 DO m=NY,1,-1
 DO j=1,topcell(i,m) -1


   pp2(i,m,j)=1./kp(i,m,j)*  beta*( pp2(i,m+1,j)/dy2+pp2(i,m-1,j)/dy2   &
                                  + pp2(i+1,m,j)/dx2+pp2(i-1,m,j)/dx2   &
                                  + pp2(i,m,j+1)/dz2+pp2(i,m,j-1)/dz2   &
                                  - matrix(i,m,j) )  &           !- den_w/dt*matrix(i,m,j) )  &
                       + (1.-beta)*pp2(i,m,j)
 
 ENDDO
 ENDDO
 ENDDO

ENDDO

ENDIF

ENDSUBROUTINE POISSON_Revised