SUBROUTINE MatrixUpdate(Matrix,                  & !OUT
                        pp2,kp,NXp,NYp,NZp)        !IN

USE COMDAT

!INTEGER(4) ,DIMENSION(0:NX+1,0:NY+1)        ,INTENT(IN)    :: topcell
INTEGER(4)                                     ,INTENT(IN)    :: NXp,NYp,NZp
REAL(8)    ,DIMENSION(1:NXp  ,1:NYp  ,1:NZp )  ,INTENT(OUT)   :: Matrix
REAL(8)    ,DIMENSION(0:NXp+1,0:NYp+1,0:NZp+1) ,INTENT(IN)    :: pp2
REAL(8)    ,DIMENSION(1:NXp  ,1:NYp  ,1:NZp)   ,INTENT(IN)    :: kp

Matrix=0

dx2=dx*dx
dy2=dy*dy
dz2=dz*dz

IF (NYp==1) THEN

DO iSORPOS= 1, ii_SOR_Pos

 DO i=1,NXp          !NX1B,NX1E !i=1,NX
 DO j=1,NZp          !-1


   Matrix(i,1,j) =    pp2(i+1,1,j)/dx2+pp2(i-1,1,j)/dx2   &
                    + pp2(i,1,j+1)/dz2+pp2(i,1,j-1)/dz2   &
                      -kp(i,1,j)*pp2(i,1,j)                     
 ENDDO
 ENDDO

ENDDO

ELSE
 
DO iSORPOS= 1, ii_SOR_Pos

 DO i=1,NXp !NX1B,NX1E !i=1,NX
 DO m=NYp,1,-1
 DO j=1,NZp !-1

   Matrix(i,m,j) =    pp2(i,m+1,j)/dy2+pp2(i,m-1,j)/dy2   &
                    + pp2(i+1,m,j)/dx2+pp2(i-1,m,j)/dx2   &
                    + pp2(i,m,j+1)/dz2+pp2(i,m,j-1)/dz2   &
                    - kp(i,m,j)*pp2(i,m,j) 
 
 ENDDO
 ENDDO
 ENDDO

ENDDO

ENDIF

ENDSUBROUTINE MatrixUpdate