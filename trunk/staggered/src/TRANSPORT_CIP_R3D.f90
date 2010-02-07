SUBROUTINE TRANSPORT_CIP_R3D(ffOUT,ffXOUT,ffYOUT,ffZOUT,        & ! OUT
                             ffIN, ffXIN, ffYIN, ffZIN,         & ! IN
                              uu,vv,ww) !,              & ! IN
                         !KdifX,KdifZ,                & ! IN
                         !NX,NZ,dx,dz,dt)               ! IN common
!__________________________________________________________________
!surface BC is difficult to set, use ff=0 at topcell 

USE COMDAT



REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(OUT)   :: ffOUT,ffXOUT,ffYOUT,ffZOUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: ffIN, ffXIN ,ffYIN ,ffZIN

REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(IN)    :: uu
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(IN)    :: vv
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(IN)    :: ww

REAL(8)    ,DIMENSION(0:NX+1,0:NZ+1)    :: Ua2d,Wa2d

!REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: voll
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: fftemp,ffXtemp,ffYtemp,ffZtemp
REAL(8)    ,DIMENSION(0:NX+1,0:NZ+1)    :: fftemp2d,ffXtemp2d,ffZtemp2d
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: difX,difY,difZ,Ua,Va,Wa
REAL(8) :: Max_ff,Max_ffX, Max_ffY, Max_ffZ, Max_all
!__________________________________________________________________________
fftemp=0;ffXtemp=0;ffYtemp=0;ffZtemp=0
difX=0;difY=0;difZ=0;Ua=0;Va=0;Wa=0



fftemp = ffIN ; ffXtemp= ffXIN ; ffYtemp= ffYIN ; ffZtemp= ffZIN





DO i=1,NX
DO m=1,NY
DO j=1,NZ   !NZ

IF (uu(i-1,m,j)*uu(i,m,j)>0) THEN
Ua(i,m,j)=transp_x*SIGN(1.0 , uu(i,m,j))*MAX(ABS(uu(i-1,m,j)),ABS(uu(i,m,j)))+(1-transp_x)*0.5*(uu(i-1,m,j)+uu(i,m,j))
ELSE 
Ua(i,m,j)=0.5*(uu(i-1,m,j)+uu(i,m,j))
ENDIF

IF (vv(i,m-1,j)*vv(i,m,j)>0) THEN
Va(i,m,j)=transp_y*SIGN(1.0 , vv(i,m,j))*MAX(ABS(vv(i,m-1,j)),ABS(vv(i,m,j)))+(1-transp_y)*0.5*(vv(i,m-1,j)+vv(i,m,j))
ELSE 
Va(i,m,j)=0.5*(vv(i,m-1,j)+vv(i,m,j))
ENDIF

ENDDO
ENDDO
ENDDO

DO i=1,NX
DO m=1,NY
DO j=1,NZ-1   !NZ

IF (ww(i,m,j-1)*ww(i,m,j)>0) THEN
Wa(i,m,j)=transp_z*SIGN(1.0 , ww(i,m,j))*MAX(ABS(ww(i,m,j-1)),ABS(ww(i,m,j)))+(1-transp_z)*0.5*(ww(i,m,j-1)+ww(i,m,j))
ELSE 
Wa(i,m,j)=0.5*(ww(i,m,j-1)+ww(i,m,j))
ENDIF

ENDDO
ENDDO
ENDDO

!---------------------------  	


DO  i=1,NX
DO  m=1,NY
DO  j=1,NZ
 IF (fftemp(i,m,j)<0) fftemp(i,m,j)=0
 IF (fftemp(i,m,j)>ff_initial) fftemp(i,m,j)=ff_initial !05/02/08
ENDDO
ENDDO
ENDDO


      
 CALL BDC3D_R(fftemp,ffXtemp,ffYtemp,ffZtemp)
!---------------------------

!---------------------------------------------------------
 DO i=1,NX
 DO m=1,NY
 
!   DifX(i,m,NZ) = 1./dx*(   KdifX*(fftemp(i+1,m,j)-fftemp(i,m,j))/dx    &
!                      - KdifX*(fftemp(i,m,j)-fftemp(i-1,m,j))/dx     )
!
!   DifY(i,m,NZ) = 1./dy*(   KdifY*(fftemp(i,m+1,j)-fftemp(i,m,j))/dy    &
!                      - KdifY*(fftemp(i,m,j)-fftemp(i,m-1,j))/dy     )
!
!   DifZ(i,m,NZ) = 1./dz*(                                    &
!                      - KdifZ*(fftemp(i,m,j)-fftemp(i,m,j-1))/dz     )
 
 Do j=1,NZ !NZ-1
   DifX(i,m,j) = 1./dx*(   KdifX*(fftemp(i+1,m,j)-fftemp(i,m,j))/dx    &
                        - KdifX*(fftemp(i,m,j)-fftemp(i-1,m,j))/dx     )
                        
   DifY(i,m,j) = 1./dy*(   KdifY*(fftemp(i,m+1,j)-fftemp(i,m,j))/dy    &
                        - KdifY*(fftemp(i,m,j)-fftemp(i,m-1,j))/dy     )                     

   DifZ(i,m,j) = 1./dz*(   KdifZ*(fftemp(i,m,j+1)-fftemp(i,m,j))/dz    &
                        - KdifZ*(fftemp(i,m,j)-fftemp(i,m,j-1))/dz     )
 ENDDO
 ENDDO
 ENDDO                     
!---------------------------------------------------------

 DO i=1,NX
 DO m=1,NY
 DO j=1,NZ
  fftemp(i,m,j)=fftemp(i,m,j)+dt*( DifX(i,m,j)+DifY(i,m,j)+DifZ(i,m,j) )
 ENDDO
 ENDDO
 ENDDO
 
 !--------------------------- 	
DO  i=1,NX
DO  m=1,NY
DO  j=1,NZ
 IF (fftemp(i,m,j)<0) fftemp(i,m,j)=0
 IF (fftemp(i,m,j)>ff_initial) fftemp(i,m,j)=ff_initial !05/02/08
ENDDO
ENDDO
ENDDO
!---------------------------
 
!    CALL   BC_SURFACE_F_R(fftemp,                   & !INOUT
!                        topcell,                  & !IN
!                        voll)                      !IN common    
 
    CALL BDC3D_R(fftemp,ffXtemp,ffYtemp,ffZtemp) 
      Max_ff =MAXVAL(ABS(fftemp))
      Max_ffX=MAXVAL(ABS(ffXtemp));Max_ffY=MAXVAL(ABS(ffYtemp));Max_ffZ=MAXVAL(ABS(ffZtemp))
      Max_all=MAX(Max_ffX,Max_ffY,Max_ffZ)    
    WRITE(18,'(f8.1,4X,E9.2,4X,E9.2,4X,E9.2)') Max_ff,Max_ffX,Max_ffY,Max_ffZ
    WRITE(*,'(f8.1,4X,E9.2,4X,E9.2,4X,E9.2)') Max_ff,Max_ffX,Max_ffY,Max_ffZ
   
    IF (CIP2D==0.OR.NY>1) THEN
     CALL RCIP3D_R(fftemp,ffXtemp,ffYtemp,ffZtemp,Ua,Va,Wa,dx,dy,dz,dt,1,NX,NX,NY,NZ)
    
     
    ELSEIF(NY==1.and.CIP2D>0) THEN
     fftemp2d  = fftemp(:,1,:)
     ffXtemp2d = ffXtemp(:,1,:)
     ffZtemp2d = ffZtemp(:,1,:)
     Ua2d      = Ua(:,1,:)
     Wa2d      = Wa(:,1,:)
     
     CALL RCIP2D(fftemp2d,ffXtemp2d,ffZtemp2d,Ua2d,Wa2d,NX,NZ,dx,dz,dt)
     
     fftemp(:,1,:) = fftemp2d
     ffXtemp(:,1,:)= ffXtemp2d
     ffZtemp(:,1,:) = ffZtemp2d
     ffYtemp = 0
        
     ELSE
     WRITE(*,*) 'check CIP2D & NY';STOP
     
    ENDIF

        

DO  i=1,NX
DO  m=1,NY
DO  j=1,NZ
 IF (fftemp(i,m,j)<0) fftemp(i,m,j)=0
 IF (fftemp(i,m,j)>ff_initial) fftemp(i,m,j)=ff_initial !05/02/08
ENDDO
ENDDO
ENDDO

ffOUT  = fftemp
ffXOUT = ffXtemp
ffYOUT = ffYtemp
ffZOUT = ffZtemp
  
CONTINUE

ENDSUBROUTINE TRANSPORT_CIP_R3D