SUBROUTINE SIMPLE_R_S(ppOUT,uuOUT,vvOUT,wwOUT,          & !OUT
                      ppIN, uuIN, vvIN, wwIN,           & !IN
                      topcell,domain,kp)                  !IN
! skip right boundary cell, no residuls are computed at NX1E 


USE COMDAT

INTEGER(4) ,DIMENSION(0:NX+1,0:NY+1)        ,INTENT(IN)    :: topcell
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: ppIN
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(OUT)   :: ppOUT

REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(IN)    :: uuIN
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(IN)    :: wwIN
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(IN)    :: vvIN
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(OUT)   :: uuOUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(OUT)   :: wwOUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(OUT)   :: vvOUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: kp

REAL(8)    ,DIMENSION(0:nx+1,0:NY+1,0:NZ+1)    :: PPtemp
REAL(8)    ,DIMENSION(1:NX  ,1:NY  ,1:NZ  )    :: Matrix
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: pp2,PPc,PPo
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1)    :: uutemp
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  )    :: wwtemp
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1)    :: vvtemp
INTEGER(2), INTENT(IN) :: domain
INTEGER(4)             :: beginNX, endNX

ppOUT=0
uuOUT=0
wwOUT=0
vvOUT=0
pp2=0
PPc=0

uutemp=uuIN
wwtemp=wwIN
vvtemp=vvIN
PPtemp=ppIN

PPo=PPtemp

!!!!!! moving backward one cell produce perfect results
Select CASE (domain)
CASE(0)
  beginNX=NX1B;endNX=NX1E
CASE(1) !left domain, skip right
  beginNX=NX1B;endNX=NX1E-1
CASE(2) !right domain, skip left
  beginNX=NX1B+1;endNX=NX1E
ENDSELECT

IF((ANH==1).OR.(ANH==2)) THEN ! NonHydrostatic without Hydrostatic Correction and Adjustable NonHydrostatic.

    LOOP_1: DO iSOR=1,ii_SOR 
    !************************************************************************************

    matrix = 0

        DO m=1,NY
        DO i=beginNX,endNX !NX1B,NX1E !i=NX1B,NX1E !i=1,nx
        DO j=1,topcell(i,m) -1
        matrix(i,m,j) =           (                                            &
                                      (1./dx)*(uutemp(i,m,j)-uutemp(i-1,m,j))  &
                                     +(1./dy)*(vvtemp(i,m,j)-vvtemp(i,m-1,j))  &
                                     +(1./dz)*(wwtemp(i,m,j)-wwtemp(i,m,j-1))  & 
                                  )
        ENDDO;ENDDO;ENDDO

        err=MAXVAL(ABS(matrix))


        IF (iSOR==1) THEN
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
        ENDIF 
 


        Matrix=Matrix*den_w/dt

        CALL      POISSON_Revised(pp2,                            & !OUT
                                   Matrix,topcell,kp)               !IN


        ! Enforcing p=0 >= surface cell  
        !------------------
        DO m=0,NY+1
        DO i=beginNX-1,endNX+1
        DO j=NZ,NZ+1
               PPc(i,m,j)=0
            PPtemp(i,m,j)=0
        ENDDO;ENDDO ;ENDDO

        PPtemp=PPtemp+pp2
        PPc=PPc+pp2

        !--------------------------------------------------------------------------------
        !(4) Correct uu,vv,ww
        DO m=1,NY
        DO i=beginNX,endNX-1  !DO i=NX1B,NX1E-1 !i=1,nx-1
        DO j=1,topcell(i,m)
	         !uutemp(i,m,j)=uuIN(i,m,j) -dt/dx*(PPc(i+1,m,j)-PPc(i,m,j))/(0.5*(den(i+1,m,j)+den(i,m,j)))
	         uutemp(i,m,j)=uuIN(i,m,j) -dt/dx*(PPc(i+1,m,j)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO

        DO m=1,NY-1
        DO i=beginNX,endNX !DO i=NX1B,NX1E !i=1,nx
        DO j=1,topcell(i,m)
	         vvtemp(i,m,j)=vvIN(i,m,j) -dt/dy*(PPc(i,m+1,j)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO
        
        DO m=1,NY
        DO i=beginNX,endNX        
        DO j=1,topcell(i,m)-1 
	         wwtemp(i,m,j) = wwIN(i,m,j)-dt/dz*(PPc(i,m,j+1)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO
 
        DO m=1,NY
        DO i=beginNX,endNX       
        DO j=topcell(i,m),topcell(i,m)
!             wwtemp(i,m,j) =  wwtemp(i,m,j-1)
             wwtemp(i,m,j)   = wwtemp(i,m,j-1) -  (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                               -  (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO
    
!        CALL       BC_SURFACE_R(uutemp,vvtemp,wwtemp,              & !INOUT
!                                  w_kineticIN,topcell)               !IN common

!************************************************************************************
 IF (ISOR_ErrMax==211) THEN

        IF ((err<=errmax).and.(iSOR>=ii_SOR_Min)) then
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
            EXIT LOOP_1
        ENDIF
 
 ELSEIF (ISOR_ErrMax==210) THEN

        IF (iSOR>=ii_SOR_Min) then
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
            EXIT LOOP_1
        ENDIF 
 
 ELSEIF (ISOR_ErrMax==201) THEN

        IF ((err<=errmax)) then
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
            EXIT LOOP_1
        ENDIF 
 
 ELSEIF (ISOR_ErrMax==200) THEN

        IF ((err<=errmax).or.(iSOR>=ii_SOR_Min)) then
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
            EXIT LOOP_1
        ENDIF
 ELSE
    WRITE(*,*) 'check ISOR_ErrMax';STOP
 ENDIF
!************************************************************************************

    ENDDO LOOP_1
    
    
    IF( ANH == 2) THEN !Nonhydrostatic with Hydrostatic Correction -> Adjustable NonHydrostatic

        DO j=1,topcell(i,m) !-1
        DO m=1,NY
        DO i=beginNX,endNX	 
             wwtemp(i,m,j)   = wwtemp(i,m,j-1) -  (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                               -  (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO
        
    ENDIF

ELSEIF (ANH == 0) THEN ! Fully Hydrostatic

        DO j=1,topcell(i,m) !-1
        DO m=1,NY
        DO i=beginNX,endNX	 
             wwtemp(i,m,j)   = wwtemp(i,m,j-1) -  (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                               -  (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO

ENDIF    





! Enforcing p=0 >= surface cell  
!------------------
DO m=0,NY+1
DO i=NX1B,NX1E !i=1,NX
       PPtemp(i,m,NZ+1)=0
ENDDO
ENDDO

ppOUT=PPtemp
uuOUT=uutemp
vvOUT=vvtemp
wwOUT=wwtemp

CONTINUE


ENDSUBROUTINE SIMPLE_R_S