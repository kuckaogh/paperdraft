SUBROUTINE POISSON_Revised_Multigrid(ppp,                                & !OUT
                                     Residual,kp,NXp,NYp,NZp)        !IN

USE COMDAT

INTEGER(4)                                ,INTENT(IN)   :: NXp,NYp,NZp
REAL(8),DIMENSION(1:NXp  ,1:NYp  ,1:NZp  ),INTENT(IN)   :: Residual
REAL(8),DIMENSION(0:NXp+1,0:NYp+1,0:NZp+1),INTENT(INOUT):: ppp
REAL(8),DIMENSION(1:NXp  ,1:NYp  ,1:NZp  ),INTENT(IN)   :: kp


dx2=dx*dx
dy2=dy*dy
dz2=dz*dz

IF (NYp==1) THEN

DO iSORPOS= 1, ii_SOR_Pos

 DO i=1,NXp; DO j=1,NZp 


   ppp(i,1,j)=1./kp(i,1,j)*  beta*( ppp(i+1,1,j)/dx2+ppp(i-1,1,j)/dx2   &
                                  + ppp(i,1,j+1)/dz2+ppp(i,1,j-1)/dz2   &
                                  - Residual(i,1,j) )  &           !- den_w/dt*Residual(i,m,j) )  &
                       + (1.-beta)*ppp(i,1,j)
 
 ENDDO;ENDDO

ENDDO

ELSE
 
DO iSORPOS= 1, ii_SOR_Pos

 DO i=1,NXp; DO m=NYp,1,-1; DO j=1,NZp


   ppp(i,m,j)=1./kp(i,m,j)*  beta*( ppp(i,m+1,j)/dy2+ppp(i,m-1,j)/dy2   &
                                  + ppp(i+1,m,j)/dx2+ppp(i-1,m,j)/dx2   &
                                  + ppp(i,m,j+1)/dz2+ppp(i,m,j-1)/dz2   &
                                  - Residual(i,m,j) )  &           !- den_w/dt*Residual(i,m,j) )  &
                       + (1.-beta)*ppp(i,m,j)
 
 ENDDO;ENDDO;ENDDO

ENDDO

ENDIF

ENDSUBROUTINE POISSON_Revised_Multigrid