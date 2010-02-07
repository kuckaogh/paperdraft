SUBROUTINE RCIP3D_R(F,FX,FY,FZ,U,V,W,DX,DY,DZ,DT,     &
                         NX_Start,NX_end,NX,NY,NZ)

!C     NX,NY,NZ in doloop and dimensions modified to extend to NX+1,NY+1,NZ+1 
!
!C     *******************************************************
!C     * A SUBROUTINE FOR CALCULATING 3-D ADVECTION EQUATION *
!C     * BASED ON AN INTERPOLATION OF RATIONAL FUNCTION.     *
!C     * FX: FIRST ORDER DERIVATIVE IN X-DIRECTION (3D ARRAY)*
!C     * FY: FIRST ORDER DERIVATIVE IN Y-DIRECTION (3D ARRAY)*
!C     * FZ: FIRST ORDER DERIVATIVE IN Z-DIRECTION (3D ARRAY)*
!C     * U:  ADVECTION SPEED IN X-DIRECTION (3D ARRAY)       *
!C     * V:  ADVECTION SPEED IN Y-DIRECTION (3D ARRAY)       *
!C     * W:  ADVECTION SPEED IN Z-DIRECTION (3D ARRAY)       *
!C     * DX: GRID SPACING IN X-DIRECTION                     *
!C     * DY: GRID SPACING IN Y-DIRECTION                     *
!C     * DZ: GRID SPACING IN Z-DIRECTION                     *
!C     * DT: TIME STEP                                       *
!C     *******************************************************
  
    IMPLICIT NONE

INTEGER(4) NX,NY,NZ,NX_Start,NX_End
INTEGER(4) i,j,k, IX, JY, KZ
REAL(8) dt,dx,dy,dz
REAL(8) F,FX,FY,FZ,U,V,W
REAL(8) FN,FXN,FYN,FZN
REAL(8) thetaCIP,au,av,aw,ddx,ddy,ddz
REAL(8) RXX,RYY,RZZ
REAL(8) FXS1,FYS1,FZS1
REAL(8) xi,yi,zi,sx,sy,sz
REAL(8) FXS0,FYS0,FZS0
REAL(8) SUMXY,SUMXZ,SUMYZ 
REAL(8) abs,c1,c2,c3,cxyz,fxdd,fydd,fzdd
REAL(8) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12
REAL(8) a13,a14,a15,a16,a17,a18,a19
DIMENSION  F(0:NX+1,0:NY+1,0:NZ+1),FX(0:NX+1,0:NY+1,0:NZ+1)
DIMENSION FY(0:NX+1,0:NY+1,0:NZ+1),FZ(0:NX+1,0:NY+1,0:NZ+1)
DIMENSION  U(0:NX+1,0:NY+1,0:NZ+1), V(0:NX+1,0:NY+1,0:NZ+1)
DIMENSION  W(0:NX+1,0:NY+1,0:NZ+1)
DIMENSION  FN(0:NX+1,0:NY+1,0:NZ+1),FXN(0:NX+1,0:NY+1,0:NZ+1)
DIMENSION FYN(0:NX+1,0:NY+1,0:NZ+1),FZN(0:NX+1,0:NY+1,0:NZ+1)
THETACIP=1.0 
RXX=1.0
RYY=1.0
RZZ=1.0

DO I=NX_Start-1,NX_end+1
DO J=0,NY+1
DO K=0,NZ+1
FN (I,J,K)=F (I,J,K)
FXN(I,J,K)=FX(I,J,K)
FYN(I,J,K)=FY(I,J,K)
FZN(I,J,K)=FZ(I,J,K)

ENDDO;ENDDO;ENDDO

DO 102 I=NX_Start,NX_end
DO 101 J=1,NY
DO 100 K=1,NZ

AU=U(I,J,K)
AV=V(I,J,K)
AW=W(I,J,K)

IF(AU.LE.0.0) THEN
DDX=DX
IX=1
ELSE 
DDX=-DX
IX=-1
ENDIF

IF(AV.LE.0.0) THEN
DDY=DY
JY=1
ELSE 
DDY=-DY
JY=-1
ENDIF

IF(AW.LE.0.0) THEN
DDZ=DZ
KZ=1
ELSE 
DDZ=-DZ
KZ=-1
ENDIF

XI=-AU*DT
YI=-AV*DT
ZI=-AW*DT

SX=(FN(I+IX,J,K)-FN(I,J,K))/DDX
SY=(FN(I,J+JY,K)-FN(I,J,K))/DDY
SZ=(FN(I,J,K+KZ)-FN(I,J,K))/DDZ

FXS1=FXN(I+IX,J,K)-SX
FYS1=FYN(I,J+JY,K)-SY
FZS1=FZN(I,J,K+KZ)-SZ

IF(ABS(FXS1).LE.1.0D-10) THEN
FXS1=1.0D-10
FXN(I,J,K)=SX
ENDIF

IF(ABS(FYS1).LE.1.0D-10) THEN
FYS1=1.0D-10
FYN(I,J,K)=SY
ENDIF

IF(ABS(FZS1).LE.1.0D-10) THEN
FZS1=1.0D-10
FZN(I,J,K)=SZ
ENDIF

FXDD=FXN(I+IX,J,K)*FXN(I,J,K)
FYDD=FYN(I,J+JY,K)*FYN(I,J,K)
FZDD=FZN(I,J,K+KZ)*FZN(I,J,K)

FXS0=SX-FXN(I,J,K)
C1=(ABS(FXS0/FXS1)-1.)/DDX
C1=C1*THETACIP*RXX
IF(FXDD.GT.0.0) C1=0.0

FYS0=SY-FYN(I,J,K)
C2=(ABS(FYS0/FYS1)-1.)/DDY
C2=C2*THETACIP*RYY
IF(FYDD.GT.0.0) C2=0.0

FZS0=SZ-FZN(I,J,K)
C3=(ABS(FZS0/FZS1)-1.)/DDZ
C3=C3*THETACIP*RZZ
IF(FZDD.GT.0.0) C3=0.0

A1=FXN(I,J,K)+C1*FN(I,J,K)
A2=FYN(I,J,K)+C2*FN(I,J,K)
A3=FZN(I,J,K)+C3*FN(I,J,K)

A17=((1.+C1*DDX)*(FXN(I+IX,J,K)-SX)+(FXN(I,J,K)-SX))/(DDX*DDX)
A18=((1.+C2*DDY)*(FYN(I,J+JY,K)-SY)+(FYN(I,J,K)-SY))/(DDY*DDY)
A19=((1.+C3*DDZ)*(FZN(I,J,K+KZ)-SZ)+(FZN(I,J,K)-SZ))/(DDZ*DDZ) 

A14=((1.+C1*DDX)*FN(I+IX,J,K)-FN(I,J,K)-A1*DDX)/(DDX*DDX)    &
            -A17*DDX
A15=((1.+C2*DDY)*FN(I,J+JY,K)-FN(I,J,K)-A2*DDY)/(DDY*DDY)    &
            -A18*DDY
A16=((1.+C3*DDZ)*FN(I,J,K+KZ)-FN(I,J,K)-A3*DDZ)/(DDZ*DDZ)    &
            -A19*DDZ
     
        SUMXY=FN(I,J,K)+A1*DDX+A2*DDY+A14*DDX*DDX+A15*DDY*DDY    &
             +A17*DDX*DDX*DDX+A18*DDY*DDY*DDY     
        SUMXZ=FN(I,J,K)+A1*DDX+A3*DDZ+A14*DDX*DDX+A16*DDZ*DDZ    &
             +A17*DDX*DDX*DDX+A19*DDZ*DDZ*DDZ     
        SUMYZ=FN(I,J,K)+A2*DDY+A3*DDZ+A15*DDY*DDY+A16*DDZ*DDZ    &
             +A18*DDY*DDY*DDY+A19*DDZ*DDZ*DDZ     
  
  A7=((1.+C1*DDX+C2*DDY)*FN(I+IX,J+JY,K)                       &
           -(1.+C2*DDY)*DDX*FXN(I,J+JY,K)-C1*DDX*FN(I,J+JY,K)    &
           +A1*DDX-SUMXY)/(DDX*DDX*DDY)                
  A8=((1.+C1*DDX+C2*DDY)*FN(I+IX,J+JY,K)                       &
           -(1.+C1*DDX)*DDY*FYN(I+IX,J,K)-C2*DDY*FN(I+IX,J,K)    &
           +A2*DDY-SUMXY)/(DDX*DDY*DDY)
  A4=((1.+C1*DDX+C2*DDY)*FN(I+IX,J+JY,K)                       &
           -A7*DDX*DDX*DDY-A8*DDX*DDY*DDY-SUMXY)/(DDX*DDY)
  
  A9=((1.+C1*DDX+C3*DDZ)*FN(I+IX,J,K+KZ)                       &
           -(1.+C3*DDZ)*DDX*FXN(I,J,K+KZ)-C1*DDX*FN(I,J,K+KZ)    &
           +A1*DDX-SUMXZ)/(DDX*DDX*DDZ)
       A10=((1.+C1*DDX+C3*DDZ)*FN(I+IX,J,K+KZ)                   &
           -(1.+C1*DDX)*DDZ*FZN(I+IX,J,K)-C3*DDZ*FN(I+IX,J,K)    &
           +A3*DDZ-SUMXZ)/(DDX*DDZ*DDZ)
  A5=((1.+C1*DDX+C3*DDZ)*FN(I+IX,J,K+KZ)                       & 
           -A9*DDX*DDX*DDZ-A10*DDX*DDZ*DDZ-SUMXZ)/(DDX*DDZ)
  
       A11=((1.+C2*DDY+C3*DDZ)*FN(I,J+JY,K+KZ)                   &
           -(1.+C3*DDZ)*DDY*FYN(I,J,K+KZ)-C2*DDY*FN(I,J,K+KZ)    &
           +A2*DDY-SUMYZ)/(DDY*DDY*DDZ)
       A12=((1.+C2*DDY+C3*DDZ)*FN(I,J+JY,K+KZ)                   &
           -(1.+C2*DDY)*DDZ*FZN(I,J+JY,K)-C3*DDZ*FN(I,J+JY,K)    &
           +A3*DDZ-SUMYZ)/(DDY*DDZ*DDZ)
  A6=((1.+C2*DDY+C3*DDZ)*FN(I,J+JY,K+KZ)                       &
           -A11*DDY*DDY*DDZ-A12*DDY*DDZ*DDZ-SUMYZ)/(DDY*DDZ)
  
  A13=(1.+C1*DDX+C2*DDY+C3*DDZ)*FN(I+IX,J+JY,K+KZ)-FN(I,J,K)     &
           -A1*DDX-A2*DDY-A3*DDZ-A4*DDX*DDY-A5*DDX*DDZ-A6*DDY*DDZ  & 
           -A7*DDX*DDX*DDY-A8*DDX*DDY*DDY-A9*DDX*DDX*DDZ           & 
           -A10*DDX*DDZ*DDZ-A11*DDY*DDY*DDZ-A12*DDZ*DDZ*DDY        &
           -A14*DDX*DDX-A15*DDY*DDY-A16*DDZ*DDZ                    &
           -A17*DDX*DDX*DDX-A18*DDY*DDY*DDY-A19*DDZ*DDZ*DDZ
     
  
  F(I,J,K)=FN(I,J,K)                                          &
           +A1*XI+A2*YI+A3*ZI+A4*XI*YI+A5*XI*ZI+A6*YI*ZI        &
           +A7*XI*XI*YI+A8*XI*YI*YI+A9*XI*XI*ZI                 &
           +A10*XI*ZI*ZI+A11*YI*YI*ZI+A12*ZI*ZI*YI              &
           +A13*XI*YI*ZI+A14*XI*XI+A15*YI*YI+A16*ZI*ZI          & 
           +A17*XI*XI*XI+A18*YI*YI*YI+A19*ZI*ZI*ZI
     
      CXYZ=1.+C1*XI+C2*YI+C3*ZI
      
  F(I,J,K)=F(I,J,K)/CXYZ   
   
  FX(I,J,K)=-C1*F(I,J,K)                      &
           +A1+A4*YI+A5*ZI                      &
           +2.*A7*XI*YI+A8*YI*YI+2.*A9*XI*ZI    &
           +A10*ZI*ZI                           &
           +A13*YI*ZI+2.*A14*XI                 &
           +3.*A17*XI*XI
  FX(I,J,K)=FX(I,J,K)/CXYZ   
  
  FY(I,J,K)=-C2*F(I,J,K)                      &
           +A2+A4*XI+A6*ZI                      &
           +A7*XI*XI+2.*A8*XI*YI+2.*A11*YI*ZI   &
           +A12*ZI*ZI                           &
           +A13*XI*ZI+2.*A15*YI                 &
           +3.*A18*YI*YI
  FY(I,J,K)=FY(I,J,K)/CXYZ   
  
  FZ(I,J,K)=-C3*F(I,J,K)                      & 
           +A3+A5*XI+A6*YI                      & 
           +A9*XI*XI+2.*A10*XI*ZI+A11*YI*YI     &
           +2.*A12*ZI*YI                        &
           +A13*XI*YI+2.*A16*ZI                 &  
           +3.*A19*ZI*ZI
  FZ(I,J,K)=FZ(I,J,K)/CXYZ   

2000  CONTINUE
100 ENDDO
101 ENDDO
102 ENDDO

      DO I=NX_Start,NX_end
      DO J=1,NY
      DO K=1,NZ
        FXN(I,J,K) = FX (I,J,K)
        FYN(I,J,K) = FY (I,J,K)
        FZN(I,J,K) = FZ (I,J,K)
      ENDDO;ENDDO;ENDDO
       
      DO  I=NX_Start,NX_end
      DO  J=1,NY
      DO  K=1,NZ
        FX(I,J,K) = FXN (I,J,K)                                    &
                  - FXN (I,J,K)*(U(I+1,J,K)-U(I-1,J,K))*0.5*DT/DX  &
                  - FYN (I,J,K)*(V(I+1,J,K)-V(I-1,J,K))*0.5*DT/DX  &
                  - FZN (I,J,K)*(W(I+1,J,K)-W(I-1,J,K))*0.5*DT/DX
     
        FY(I,J,K) = FYN (I,J,K)                                    &
                  - FXN (I,J,K)*(U(I,J+1,K)-U(I,J-1,K))*0.5*DT/DY  &
                  - FYN (I,J,K)*(V(I,J+1,K)-V(I,J-1,K))*0.5*DT/DY  &
                  - FZN (I,J,K)*(W(I,J+1,K)-W(I,J-1,K))*0.5*DT/DY

        FZ(I,J,K) = FZN (I,J,K)                                    &
                  - FXN (I,J,K)*(U(I,J,K+1)-U(I,J,K-1))*0.5*DT/DZ  &
                  - FYN (I,J,K)*(V(I,J,K+1)-V(I,J,K-1))*0.5*DT/DZ  &
                  - FZN (I,J,K)*(W(I,J,K+1)-W(I,J,K-1))*0.5*DT/DZ
     
      ENDDO;ENDDO;ENDDO 

  END



