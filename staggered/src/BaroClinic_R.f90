SUBROUTINE BaroClinic_R(BaCX,BaCY,                 &   ! OUT                                 
                       den)                     ! IN                
USE COMDAT
REAL(8)    ,DIMENSION(0:NX  ,0:NY  ,0:NZ+1)   ,INTENT(OUT)   :: BaCX, BaCY
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)   ,INTENT(IN)    :: den






BaCX  = 0
BaCY  = 0

!Compute g*sum[d rho/dx]: Baroclinic term
!******************************************

DO m=1,NY !NY-1
DO i=1,NX-1 !i=1,NX-1

  BaCX(i,m,NZ)  =  gravity*( den(i+1,m,NZ)-den(i,m,NZ) )*dz/dx



DO j=NZ-1,1,-1

  BaCX(i,m,j) = BaCX(i,m,j+1) + &
   & gravity*( den(i+1,m,j)-den(i,m,j) )*dz/dx

ENDDO
ENDDO
ENDDO


DO i=1,NX !i=1,NX-1
DO m=1,NY-1 !NY-1

  BaCY(i,m,NZ)  = &
   & gravity*( den(i,m+1,NZ)-den(i,m,NZ) )*dz/dy



DO j=NZ-1,1,-1

  BaCY(i,m,j) = BaCY(i,m,j+1) + &
   & gravity*( den(i,m+1,j)-den(i,m,j) )*dz/dy

ENDDO
ENDDO
ENDDO

IF (Bossi == 0) THEN

DO j=1,NZ
DO m=1,NY !NY-1
DO i=1,NX-1 !i=1,NX-1
  BaCX(i,m,j) = BaCX(i,m,j)/( den(i+1,m,j)+den(i,m,j) )
  
ENDDO;ENDDO;ENDDO

DO j=1,NZ  
DO i=1,NX ! i=1,NX-1
DO m=1,NY-1

  BaCY(i,m,j) = BaCY(i,m,j)/( den(i,m+1,j)+den(i,m,j) )  
   
ENDDO;ENDDO;ENDDO

ELSE

DO j=1,NZ
DO m=1,NY !NY-1
DO i=1,NX-1 !i=1,NX-1
  BaCX(i,m,j) = BaCX(i,m,j)/den_w
  
ENDDO;ENDDO;ENDDO

DO j=1,NZ  
DO i=1,NX ! i=1,NX-1
DO m=1,NY-1

  BaCY(i,m,j) = BaCY(i,m,j)/den_w  
   
ENDDO;ENDDO;ENDDO

ENDIF
!******************************************

ENDSUBROUTINE BaroClinic_R

