SUBROUTINE HEIGHTONLY_R(volsum,                                   &   !OUT
                      vollIN_1,                                   &   !IN
                      uu,vv,volsum_old_1)                             !IN

USE COMDAT
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(IN)    :: uu
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(IN)    :: vv
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)        ,INTENT(IN)    :: volsum_old_1
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(IN)    ::     vollIN_1
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)                       :: volsum_old
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)                ::     vollIN



REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)        ,INTENT(OUT)   :: volsum

REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)               ::   volsumOUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)        ::     vollTEMP


REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1,4)      :: in
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1)               :: SS

volsum_old  = (1./dz) *volsum_old_1 
vollIN      = (1./dz) *vollIN_1

volsumOUT   = volsum_old
 vollTEMP   = vollIN

 volsum       = 0


in = 0

!**************************************************************************************

DO i=1,NX !i = 1, NX
DO m=1,NY
DO j=1,NZ

 IF (uu(i-1,m,j).GE.0) THEN

    in(i,m,j,1)=MIN(uu(i-1,m,j)*vollTEMP(i-1,m,j)*dt/dx, 0.25*vollTEMP(i-1,m,j))

 ELSE 
    
    in(i,m,j,1)=-MIN(-uu(i-1,m,j)*vollTEMP(i,m,j)*dt/dx, 0.25*vollTEMP(i,m,j))

 ENDIF
    
 IF (uu(i,m,j).GE.0) THEN

    in(i,m,j,2)=-MIN(uu(i,m,j)*vollTEMP(i,m,j)*dt/dx, 0.25*vollTEMP(i,m,j))

   ELSE 

    in(i,m,j,2)=MIN(-uu(i,m,j)*vollTEMP(i+1,m,j)*dt/dx, 0.25*vollTEMP(i+1,m,j))
    
 ENDIF
 !goto 245 !skip problem area
 !******************************** v ********************************************** problem
 
 IF (vv(i,m-1,j).GE.0) THEN

    in(i,m,j,3)=MIN(vv(i,m-1,j)*vollTEMP(i,m-1,j)*dt/dy, 0.25*vollTEMP(i,m-1,j))

 ELSE 
    
    in(i,m,j,3)=-MIN(-vv(i,m-1,j)*vollTEMP(i,m,j)*dt/dy, 0.25*vollTEMP(i,m,j))

 ENDIF
    
 IF (vv(i,m,j).GE.0) THEN

    in(i,m,j,4)=-MIN(vv(i,m,j)*vollTEMP(i,m,j)*dt/dy, 0.25*vollTEMP(i,m,j))

 ELSE 

    in(i,m,j,4)=MIN(-vv(i,m,j)*vollTEMP(i,m+1,j)*dt/dy, 0.25*vollTEMP(i,m+1,j))
    
 ENDIF

!********************************************************************************** 
245 continue
ENDDO
ENDDO
ENDDO

SS = 0
DO i=1,NX !i=1,NX
DO m=1,NY
DO j=1,NZ
SS(i,m)=SS(i,m)+in(i,m,j,1)+in(i,m,j,2)+in(i,m,j,3)+in(i,m,j,4)
ENDDO
ENDDO
ENDDO

DO m=0,NY+1
DO i=1-1,NX+1 !i=0,NX+1
   volsumOUT(i,m) = volsum_old(i,m) + SS(i,m)
ENDDO 
ENDDO



!********************************************************************************************

! Assign Out

 volsum=volsumOUT
 
 CONTINUE


volsum      =  (dz) *volsum 

ENDSUBROUTINE HEIGHTONLY_R

