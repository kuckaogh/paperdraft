SUBROUTINE INIT3D(ff_initial,ff_strati,F,FX,FY,FZ,NX,NY,NZ,NXS_start,NXS_end,NYS,NZS,NZS2,dx,dy,dz)

IMPLICIT NONE
REAL(8)                                  ,INTENT(IN) :: dx,dy,dz,ff_initial,ff_strati
INTEGER(4)                               ,INTENT(IN) :: NX,NY,NZ,NXS_Start,NXS_end, NYS, NZS, NZS2
REAL(8) ,DIMENSION(0:NX+1,0:NY+1,0:nz+1) ,INTENT(OUT):: F,FX,FY,FZ
INTEGER i,j,k

F=0;FX=0;FY=0;FZ=0
    
    DO  I=1,NX
  DO  J=1,NY
    DO  K=1,NZS2 
      F(I,J,K)= REAL(NZS2-K)/REAL(NZS2-1)*ff_strati 
  ENDDO
  ENDDO
  ENDDO 
  
  DO  I=NXS_start,NXS_end 
      DO  J=1,NYS 
      do  k=1,NZS 

  F(I,J,k)=ff_initial 

      enddo
      ENDDO
      ENDDO

      DO  I=1,NX
      DO  J=1,NY
      do  k=1,NZ
      FX(I,J,k)=0.5*(F(I+1,J,k)-F(I-1,J,k))/DX
      FY(I,J,k)=0.5*(F(I,J+1,k)-F(I,J-1,k))/DY
      fz(I,J,k)=0.5*(F(I,J,k+1)-F(I,J,k-1))/DZ
      ENDDO
      ENDDO
      ENDDO
      
      RETURN
      ENDSUBROUTINE