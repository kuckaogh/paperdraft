SUBROUTINE BDC3D_R(F,FX,FY,FZ)

USE COMDAT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: f,fX,fY,fZ

DO i=1-1,NX+1
!      DO i=0,NX+1
DO m=0,NZ+1
      F  (i,0,m)=F  (i,1,m)
      FX (i,0,m)=FX (i,1,m)
      FY (i,0,m)=FY (i,1,m)
      FZ (i,0,m)=FZ (i,1,m)
      
      F  (i,NY+1,m)=F  (i,NY,m)
      FX (i,NY+1,m)=FX (i,NY,m)
      FY (i,NY+1,m)=FY (i,NY,m)
      FZ (i,NY+1,m)=FZ (i,NY,m)
ENDDO
ENDDO

DO J=0,NY+1
DO m=0,NZ+1
      F  (1-1,J,m)=F  (1,J,m)
      FX (1-1,J,m)=FX (1,J,m)
      FY (1-1,J,m)=FY (1,J,m)
      FZ (1-1,J,m)=FZ (1,J,m)
      F  (NX+1,J,m)=F  (NX,J,m)
      FX (NX+1,J,m)=FX (NX,J,m)
      FY (NX+1,J,m)=FY (NX,J,m)
      FZ (NX+1,J,m)=FZ (NX,J,m)
ENDDO 
ENDDO
  
  DO J=0,NY+1
  DO i=1-1,NX+1
!      DO i=0,NX+1
      F  (i,J,0)=F  (i,J,1)
      FX (i,J,0)=FX (i,J,1)
      FY (i,J,0)=FY (i,J,1)
      FZ (i,J,0)=FZ (i,J,1)      
      F  (i,J,NZ+1)=F  (i,J,NZ)
      FX (i,J,NZ+1)=FX (i,J,NZ)
      FY (i,J,NZ+1)=FY (i,J,NZ)
      FZ (i,J,NZ+1)=FZ (i,J,NZ)
  ENDDO 
  ENDDO


      ENDSUBROUTINE BDC3D_R