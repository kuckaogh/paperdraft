SUBROUTINE ASSIGNVOL_R(vollOUT,volsumIN_1)

USE COMDAT


REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)        ,INTENT(IN)    :: volsumIN_1
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)                       :: volsumIN
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(OUT)   :: vollOUT



volsumIN= (1/dz)*volsumIN_1

vollOUT=0


DO m=0,NY+1
DO i=1-1,NX+1 !i = 0, NX+1

 IF(FLOOR(volsumIN(i,m))<NZ) THEN
  do j=1, FLOOR(volsumIN(i,m))
   vollOUT(i,m,j) =  1
  enddo
  do j=CEILING(volsumIN(i,m))+1, NZ
   vollOUT(i,m,j) = 0
  enddo
  IF (FLOOR(volsumIN(i,m)) .NE. CEILING(volsumIN(i,m))) THEN
  vollOUT(i,m,floor(volsumIN(i,m))+1) = volsumIN(i,m) - FLOOR(volsumIN(i,m))
  ENDIF
 ELSE
  do j=1, NZ-1
   vollOUT(i,m,j) =  1
  enddo
  do j = NZ, NZ
   vollOUT(i,m,j) = volsumIN(i,m) - NZ+1
  enddo
 ENDIF 
ENDDO
ENDDO




vollOUT = dz* vollOUT

ENDSUBROUTINE ASSIGNVOL_R