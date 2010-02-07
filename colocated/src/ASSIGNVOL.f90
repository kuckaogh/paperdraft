SUBROUTINE ASSIGNVOL(vollOUT,volsumIN_1,NXm ,NYm ,NZm)

USE COMDAT

INTEGER(4)   ,INTENT(IN)           :: NXm,NYm,NZm
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1)        ,INTENT(IN)    :: volsumIN_1
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1)                       :: volsumIN
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(OUT)   :: vollOUT



volsumIN= (1/dz)*volsumIN_1

vollOUT=0


DO m=0,NYm+1
DO i=0,NXm+1 !i = 0, NX+1

 IF(FLOOR(volsumIN(i,m))<NZm) THEN
  do j=1, FLOOR(volsumIN(i,m))
   vollOUT(i,m,j) =  1
  enddo
  do j=CEILING(volsumIN(i,m))+1, NZm
   vollOUT(i,m,j) = 0
  enddo
  IF (FLOOR(volsumIN(i,m)) .NE. CEILING(volsumIN(i,m))) THEN
  vollOUT(i,m,floor(volsumIN(i,m))+1) = volsumIN(i,m) - FLOOR(volsumIN(i,m))
  ENDIF
 ELSE
  do j=1, NZm-1
   vollOUT(i,m,j) =  1
  enddo
  do j = NZm, NZm
   vollOUT(i,m,j) = volsumIN(i,m) - NZm+1
  enddo
 ENDIF 
ENDDO
ENDDO




vollOUT = dz* vollOUT

ENDSUBROUTINE ASSIGNVOL