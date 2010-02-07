SUBROUTINE HEIGHTVOL_colocated(volsumOUT,vollOUT,              &   !OUT
                     uu_IN,vv_IN,volsumIN,NXm,NYm,NZm)         !IN

USE COMDAT
INTEGER(4)   ,INTENT(IN)           :: NXm,NYm,NZm
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: uu_IN
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN)    :: vv_IN
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1)         ,INTENT(IN)    :: volsumIN
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1)         ,INTENT(OUT)   :: volsumOUT
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT)   :: vollOUT
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1,4)      :: in
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1)                :: SS
REAL(8)    ,DIMENSION(0:NXm  ,0:NYm+1,0:NZm+1)        :: uu
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm  ,0:NZm+1)        :: vv
 CALL ASSIGNVOL(vollOUT,volsumIN,NXm ,NYm ,NZm)

uu=0;vv=0
DO i=1,NXm-1;DO m=1,NYm;DO j=1,NZm
uu(i,m,j)= 0.5*(uu_IN(i,m,j)+uu_IN(i+1,m,j))
ENDDO;ENDDO;ENDDO

DO i=1,NXm;DO m=1,NYm-1;DO j=1,NZm
vv(i,m,j)= 0.5*(vv_IN(i,m,j)+vv_IN(i,m+1,j))
ENDDO;ENDDO;ENDDO



volsumOUT     = (1./dz) *volsumIN
vollOUT       = (1./dz) *vollOUT

in = 0
!**************************************************************************************

DO i=1,NXm ;DO m=1,NYm ;DO j=1,NZm

 IF (uu(i-1,m,j).GE.0) THEN
    in(i,m,j,1)=MIN(uu(i-1,m,j)*vollOUT(i-1,m,j)*dt/dx, 0.25*vollOUT(i-1,m,j))
 ELSE 
    in(i,m,j,1)=-MIN(-uu(i-1,m,j)*vollOUT(i,m,j)*dt/dx, 0.25*vollOUT(i,m,j))
 ENDIF
    
 IF (uu(i,m,j).GE.0) THEN
    in(i,m,j,2)=-MIN(uu(i,m,j)*vollOUT(i,m,j)*dt/dx, 0.25*vollOUT(i,m,j))
 ELSE 
    in(i,m,j,2)=MIN(-uu(i,m,j)*vollOUT(i+1,m,j)*dt/dx, 0.25*vollOUT(i+1,m,j))
 ENDIF

 !******************************** v ********************************************** problem
 
 IF (vv(i,m-1,j).GE.0) THEN
    in(i,m,j,3)=MIN(vv(i,m-1,j)*vollOUT(i,m-1,j)*dt/dy, 0.25*vollOUT(i,m-1,j))
 ELSE 
    in(i,m,j,3)=-MIN(-vv(i,m-1,j)*vollOUT(i,m,j)*dt/dy, 0.25*vollOUT(i,m,j))
 ENDIF
    
 IF (vv(i,m,j).GE.0) THEN
    in(i,m,j,4)=-MIN(vv(i,m,j)*vollOUT(i,m,j)*dt/dy, 0.25*vollOUT(i,m,j))
 ELSE 
    in(i,m,j,4)=MIN(-vv(i,m,j)*vollOUT(i,m+1,j)*dt/dy, 0.25*vollOUT(i,m+1,j))
 ENDIF

!********************************************************************************** 

ENDDO;ENDDO;ENDDO

SS = 0

DO i=1,NXm ;DO m=1,NYm ;DO j=1,NZm
  SS(i,m)=SS(i,m)+in(i,m,j,1)+in(i,m,j,2)+in(i,m,j,3)+in(i,m,j,4)
ENDDO;ENDDO;ENDDO

DO m=0,NYm+1 ;DO i=0,NXm+1
  volsumOUT(i,m) = volsumOUT(i,m) + SS(i,m)
ENDDO;ENDDO



!********************************************************************************************
! Assign Out

 volsumOUT  =  (dz) *volsumOUT 

 CALL ASSIGNVOL(vollOUT,volsumOUT,NXm ,NYm ,NZm)

ENDSUBROUTINE HEIGHTVOL_colocated

