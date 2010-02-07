SUBROUTINE INIT_VOL(voll,volsum) 

USE COMDAT


REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)        ,INTENT(OUT)  :: volsum
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(OUT)  :: voll






!***** Assgin zero ************************

voll    = 0
volsum  = 0


!***** Assgin initial height (volsum,topcell) & volume (vol) ******
Do m=0,NY+1
DO i=0,NX+1
! volsum(i,m) = Am0*COS(pi*(dx*REAL(i-0.5))/10.0) + h0
! volsum(i,m) = Am0*COS(pi*(dx*REAL(i-0.5))/2.5) + h0

!  volsum(i,m) = h0
!  volsum(i,m) = Am0*COS(pi*(dx*REAL(i-0.5))/10.0) + h0
  volsum(i,m) = Am0*COS(pi*(dx*REAL(i-0.5))/8.1) + h0

!  volsum(i,m) = h0 +0.02*dz*(m-4)
! volsum(i,m) = Am0*COS(pi*(dy*REAL(m-0.5))/2.5) + h0
ENDDO
ENDDo


DO i=0,NX+1
DO m=0,NY+1
 volsum(i,m) = (1./dz)*volsum(i,m)
ENDDO
ENDDO

DO i=0,NX+1
DO m=0,NY+1
 IF(FLOOR(volsum(i,m))<NZ) THEN
  do j=1, FLOOR(volsum(i,m))
      voll(i,m,j) = 1
  enddo
  do j=CEILING(volsum(i,m))+1, NZ
      voll(i,m,j) = 0
  enddo
  IF (FLOOR(volsum(i,m)) .NE. CEILING(volsum(i,m))) THEN
   voll(i,m,FLOOR(volsum(i,m))+1) = volsum(i,m) - FLOOR(volsum(i,m))
  ENDIF 
 ELSE
  do j = 1, NZ-1
      voll(i,m,j) = 1
  enddo
  do j = NZ, NZ
      voll(i,m,j) = volsum(i,m) - NZ+1
  enddo   
 ENDIF
ENDDO
ENDDO

DO m=0,NY+1
DO i=0,NX+1
 NZ = NZ !CEILING(volsum(i,m))
ENDDO
ENDDO

volsum = dz* volsum
voll   = dz* voll

ENDSUBROUTINE INIT_VOL