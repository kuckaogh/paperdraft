SUBROUTINE RCIP2D(F,FX,FY,U,V,NX,NY,DX,DY,DT)

!C     *******************************************************
!C     * A SUBROUTINE FOR CALCULATING 2-D ADVECTION EQUATION *
!C     * BASED ON AN INTERPOLATION OF RATIONAL FUNCTION.     *
!C     * F:  DEPENDENT VARIABLE (2D ARRAY)                   *
!C     * FX: FIRST ORDER DERIVATIVE IN X-DIRECTION (2D ARRAY)*
!C     * FY: FIRST ORDER DERIVATIVE IN Y-DIRECTION (2D ARRAY)*
!C     * U:  ADVECTION SPEED IN X-DIRECTION (2D ARRAY)       *
!C     * V:  ADVECTION SPEED IN Y-DIRECTION (2D ARRAY)       *
!C     * DX: GRID SPACING IN X-DIRECTION                     *
!C     * DY: GRID SPACING IN Y-DIRECTION                     *
!C     * DT: TIME STEP                                       *
!C     *******************************************************

IMPLICIT NONE
INTEGER NX,NY
INTEGER i,j,ISX,JSY
REAL(8) dt,dx,dy
REAL(8) F,FX,FY,U,V
REAL(8) FN,FXN,FYN
REAL(8) theta,au,av,rxd,ryd,ddx,ddy,xi,yi,sx,sy,fxj
REAL(8) fyi,abs,cx1,c1,c2,cy2,b1,b2,b3,b4,fxdd,fydd
REAL(8) a1,a2,a3,a4,a5,a6,a7,a8,a9


DIMENSION F(0:NX+1,0:NY+1),FX(0:NX+1,0:NY+1),FY(0:NX+1,0:NY+1)
DIMENSION U(0:NX+1,0:NY+1),V(0:NX+1,0:NY+1)
DIMENSION FN(0:NX+1,0:NY+1),FXN(0:NX+1,0:NY+1),FYN(0:NX+1,0:NY+1)


THETA=0
THETA=1.0


FN=F
FXN=FX
FYN=FY


DO 101 I=1,NX
DO 100 J=1,NY
AU=U(I,J)
AV=V(I,J)
RXD=1.
RYD=1.
IF(AU.LE.0.0) THEN
DDX=DX
ISX=1
ELSE 
DDX=-DX
ISX=-1
ENDIF
IF(AV.LE.0.0) THEN
DDY=DY
JSY=1
ELSE 
DDY=-DY
JSY=-1
ENDIF
XI=-AU*DT
YI=-AV*DT

SX=(FN(I+ISX,J)-FN(I,J))/DDX
FXJ=FXN(I+ISX,J)-SX
IF(ABS(FXJ).LE.1.0D-10) THEN
FXJ=1.0D-10
FXN(I,J)=SX
ENDIF

SY=(FN(I,J+JSY)-FN(I,J))/DDY
FYI=FYN(I,J+JSY)-SY
IF(ABS(FYI).LE.1.0D-10) THEN
FYI=1.0D-10
FYN(I,J)=SY
ENDIF

FXDD=FXN(I+ISX,J)*FXN(I,J)
FYDD=FYN(I,J+JSY)*FYN(I,J)

B1=SX-FXN(I,J)
B2=FXJ
C1=(ABS(B1/B2)-1.)/DDX
IF(FXDD.LE.0.0) THEN
C1=C1*THETA*RXD
ELSE
C1=0.0
ENDIF
CX1=1.+C1*DDX

B3=SY-FYN(I,J)
B4=FYI
C2=(ABS(B3/B4)-1.)/DDY
IF(FYDD.LE.0.0) THEN
C2=C2*THETA*RYD
ELSE
C2=0.0
ENDIF
CY2=1.+C2*DDY

A1=FXN(I,J)+C1*FN(I,J)
A2=FYN(I,J)+C2*FN(I,J)
A7=(CY2*FYN(I,J+JSY)-SY*CY2+FYN(I,J)-SY)/DDY/DDY
A6=(CX1*FXN(I+ISX,J)-SX*CX1+FXN(I,J)-SX)/DDX/DDX
A9=(CY2*FN(I,J+JSY)-FN(I,J)-A2*DDY)/DDY/DDY-A7*DDY
A8=(CX1*FN(I+ISX,J)-FN(I,J)-A1*DDX)/DDX/DDX-A6*DDX
A3=CY2*FXN(I,J+JSY)/DDY+CX1*FYN(I+ISX,J)/DDX                &
          +C1*FN(I,J+JSY)/DDY+C2*FN(I+ISX,J)/DDX                  & 
          +(FN(I,J)-(1.+C1*DDX+C2*DDY)*FN(I+ISX,J+JSY))/DDX/DDY   &
          +A6*DDX*DDX/DDY+A7*DDY*DDY/DDX+A8*DDX/DDY+A9*DDY/DDX    
A5=CX1*FYN(I+ISX,J)/DDX/DDX+C2*FN(I+ISX,J)/DDX/DDX          &
          -A2/DDX/DDX-A3/DDX                                      
A4=CY2*FXN(I,J+JSY)/DDY/DDY+C1*FN(I,J+JSY)/DDY/DDY          &
          -A1/DDY/DDY-A3/DDY        

F(I,J)=FN(I,J)+A1*XI+A2*YI+A3*XI*YI+A4*XI*YI*YI             &
              +A5*XI*XI*YI+A6*XI*XI*XI+A7*YI*YI*YI+A9*YI*YI       &
              +A8*XI*XI
F(I,J)=F(I,J)/(1.+C1*XI+C2*YI)    
FX(I,J)=A1+A3*YI+A4*YI*YI+2.*A5*XI*YI+3.*A6*XI*XI           &
               +2.*A8*XI-C1*F(I,J)
FX(I,J)=FX(I,J)/(1.+C1*XI+C2*YI)    
FY(I,J)=A2+A3*XI+A5*XI*XI+2.*A4*XI*YI+3.*A7*YI*YI           &
               +2.*A9*YI-C2*F(I,J)
FY(I,J)=FY(I,J)/(1.+C1*XI+C2*YI)    

100 ENDDO
101 ENDDO

      DO I=1,NX
      DO J=1,NY
        FXN(I,J) = FX(I,J)
        FYN(I,J) = FY(I,J)
        
      ENDDO;ENDDO

      DO I=1,NX
      DO J=1,NY
        FX(I,J) = FXN (I,J)                                         &
                          -FXN (I,J)*(U(I+1,J)-U(I-1,J))*0.5*DT/DX  &
                          -FYN (I,J)*(V(I+1,J)-V(I-1,J))*0.5*DT/DX
        FY(I,J) = FYN (I,J)                                         &
                          -FXN (I,J)*(U(I,J+1)-U(I,J-1))*0.5*DT/DY  &
                          -FYN (I,J)*(V(I,J+1)-V(I,J-1))*0.5*DT/DY

      ENDDO;ENDDO

ENDSUBROUTINE RCIP2D


