SUBROUTINE BaroClinic_R(BaCX,BaCY,                 &   ! OUT                                 
                       den,NXm ,NYm ,NZm)                     ! IN                
USE COMDAT

INTEGER(4) ,INTENT(IN)           :: NXm,NYm,NZm
REAL(8)    ,DIMENSION(0:NXm  ,0:NYm  ,0:NZm+1)   ,INTENT(OUT)   :: BaCX, BaCY
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)   ,INTENT(IN)    :: den

BaCX  = 0
BaCY  = 0

!Compute g*sum[d rho/dx]: Baroclinic term
!******************************************

DO m=1,NYm !NY-1
DO i=1,NXm-1 !i=1,NX-1

  BaCX(i,m,NZm)  =  gravity*( den(i+1,m,NZm)-den(i,m,NZm) )*dz/dx



DO j=NZm-1,1,-1

  BaCX(i,m,j) = BaCX(i,m,j+1) + &
   & gravity*( den(i+1,m,j)-den(i,m,j) )*dz/dx

ENDDO
ENDDO
ENDDO


DO i=1,NXm !i=1,NX-1
DO m=1,NYm-1 !NY-1

  BaCY(i,m,NZm)  = &
   & gravity*( den(i,m+1,NZm)-den(i,m,NZm) )*dz/dy



DO j=NZm-1,1,-1

  BaCY(i,m,j) = BaCY(i,m,j+1) + &
   & gravity*( den(i,m+1,j)-den(i,m,j) )*dz/dy

ENDDO
ENDDO
ENDDO

IF (Bossi == 0) THEN

DO j=1,NZm
DO m=1,NYm !NY-1
DO i=1,NXm-1 !i=1,NX-1
  BaCX(i,m,j) = BaCX(i,m,j)/( den(i+1,m,j)+den(i,m,j) )
  
ENDDO;ENDDO;ENDDO

DO j=1,NZm  
DO i=1,NXm ! i=1,NX-1
DO m=1,NYm-1

  BaCY(i,m,j) = BaCY(i,m,j)/( den(i,m+1,j)+den(i,m,j) )  
   
ENDDO;ENDDO;ENDDO

ELSE

DO j=1,NZm
DO m=1,NYm !NY-1
DO i=1,NXm-1 !i=1,NX-1
  BaCX(i,m,j) = BaCX(i,m,j)/den_w
  
ENDDO;ENDDO;ENDDO

DO j=1,NZm  
DO i=1,NXm ! i=1,NX-1
DO m=1,NYm-1

  BaCY(i,m,j) = BaCY(i,m,j)/den_w  
   
ENDDO;ENDDO;ENDDO

ENDIF
!******************************************

ENDSUBROUTINE BaroClinic_R

