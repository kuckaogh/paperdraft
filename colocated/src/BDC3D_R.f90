SUBROUTINE BDC3D_R(F,FX,FY,FZ,NXm,NYm,NZm)

USE COMDAT
INTEGER(4) ,INTENT(IN)           :: NXm,NYm,NZm
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)    :: f,fX,fY,fZ

DO i=0,NXm+1
!      DO i=0,NX+1
DO m=0,NZm+1
      F  (i,0,m)=F  (i,1,m)
      FX (i,0,m)=FX (i,1,m)
      FY (i,0,m)=FY (i,1,m)
      FZ (i,0,m)=FZ (i,1,m)
      
      F  (i,NYm+1,m)=F  (i,NYm,m)
      FX (i,NYm+1,m)=FX (i,NYm,m)
      FY (i,NYm+1,m)=FY (i,NYm,m)
      FZ (i,NYm+1,m)=FZ (i,NYm,m)
ENDDO
ENDDO

continue

DO J=0,NYm+1
DO m=0,NZm+1
      F  (0,J,m)=F  (1,J,m)
      FX (0,J,m)=FX (1,J,m)
      FY (0,J,m)=FY (1,J,m)
      FZ (0,J,m)=FZ (1,J,m)
      F  (NXm+1,J,m)=F  (NXm,J,m)
      FX (NXm+1,J,m)=FX (NXm,J,m)
      FY (NXm+1,J,m)=FY (NXm,J,m)
      FZ (NXm+1,J,m)=FZ (NXm,J,m)
ENDDO 
ENDDO
  
  DO J=0,NYm+1
  DO i=0,NXm+1
!      DO i=0,NX+1
      F  (i,J,0)=F  (i,J,1)
      FX (i,J,0)=FX (i,J,1)
      FY (i,J,0)=FY (i,J,1)
      FZ (i,J,0)=FZ (i,J,1)      
      F  (i,J,NZm+1)=F  (i,J,NZm)
      FX (i,J,NZm+1)=FX (i,J,NZm)
      FY (i,J,NZm+1)=FY (i,J,NZm)
      FZ (i,J,NZm+1)=FZ (i,J,NZm)
  ENDDO 
  ENDDO


      ENDSUBROUTINE BDC3D_R