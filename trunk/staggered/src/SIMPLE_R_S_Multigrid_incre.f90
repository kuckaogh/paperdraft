SUBROUTINE SIMPLE_R_S_Multigrid_incre(ppOUT,uuOUT,vvOUT,wwOUT,      & !OUT
                                ppIN, uuIN, vvIN, wwIN)              !IN
! skip right boundary cell, no residuls are computed at NX 
USE COMDAT


REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: ppIN
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(OUT)   :: ppOUT

REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(IN)    :: uuIN
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(IN)    :: wwIN
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(IN)    :: vvIN
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(OUT)   :: uuOUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(OUT)   :: wwOUT
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(OUT)   :: vvOUT
REAL(8)    ,DIMENSION(1:NXm ,1:NYm ,1:NZm )     :: kp
REAL(8)    ,DIMENSION(1:NXm2,1:NYm2,1:NZm2)     :: kp2
REAL(8)    ,DIMENSION(1:NXm3,1:NYm3,1:NZm3)     :: kp3
REAL(8)    ,DIMENSION(1:NXm4,1:NYm4,1:NZm4)     :: kp4
REAL(8)    ,DIMENSION(1:NXm5,1:NYm5,1:NZm5)     :: kp5
REAL(8)    ,DIMENSION(0:nx+1,0:NY+1,0:NZ+1)    :: PPtemp
REAL(8)    ,DIMENSION(1:NX  ,1:NY  ,1:NZ  )    :: Residual, Resi, Matrix
REAL(8)    ,DIMENSION(1:NXm2,1:NYm2,1:NZm2)    :: Resi2
REAL(8)    ,DIMENSION(1:NXm3,1:NYm3,1:NZm3)    :: Resi3
REAL(8)    ,DIMENSION(1:NXm4,1:NYm4,1:NZm4)    :: Resi4
REAL(8)    ,DIMENSION(1:NXm5,1:NYm5,1:NZm5)    :: Resi5
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1)    :: ppp,PPc
REAL(8)    ,DIMENSION(0:NXm2+1,0:NYm2+1 ,0:NZm2+1)  :: ppp2
REAL(8)    ,DIMENSION(0:NXm3+1,0:NYm3+1 ,0:NZm3+1)  :: ppp3
REAL(8)    ,DIMENSION(0:NXm4+1,0:NYm4+1 ,0:NZm4+1)  :: ppp4
REAL(8)    ,DIMENSION(0:NXm5+1,0:NYm5+1 ,0:NZm5+1) :: ppp5
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1)    :: uutemp
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  )    :: wwtemp
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1)    :: vvtemp

INTEGER(4)             :: idum


uuOUT=0;wwOUT=0;vvOUT=0
ppOUT=0;PPc=0
ppp=0;ppp2=0;ppp3=0;ppp4=0;ppp5=0
Resi=0;Resi2=0;Resi3=0;Resi4=0;Resi5=0

uutemp=uuIN;wwtemp=wwIN;vvtemp=vvIN
PPtemp=ppIN

!!!!!! moving backward one cell produce perfect results



IF((ANH==1).OR.(ANH==2)) THEN ! NonHydrostatic without Hydrostatic Correction and Adjustable NonHydrostatic.
        ! find kp
        !---------------------------------------------------
        !NXm=NX; NYm=NY; NZm=NZ
        CALL INIT_KP_All(kp, NXm, NYm, NZm)        
        !NXm2=FLOOR(REAL(NXm)/2); NYm2=MAX(FLOOR(REAL(NYm)/2),1); NZm2=FLOOR(REAL(NZm)/2)
        IF (Multilevel>1) CALL INIT_KP_All(kp2,NXm2,NYm2,NZm2)
        !NXm3=FLOOR(REAL(NXm2)/2);NYm3=MAX(FLOOR(REAL(NYm2)/2),1);NZm3=FLOOR(REAL(NZm2)/2)
        IF (Multilevel>2) CALL INIT_KP_All(kp3,NXm3,NYm3,NZm3)         
        IF (Multilevel>3) CALL INIT_KP_All(kp4,NXm4,NYm4,NZm4)
        IF (Multilevel>4) CALL INIT_KP_All(kp5,NXm5,NYm5,NZm5)
        !---------------------------------------------------

   
    LOOP_1: DO iSOR=1,ii_SOR 
!************************************************************************************
    ! Find Residual
   !(============================================================================
  If (SurfaceCorrection==0) THEN
        DO m=1,NY;DO i=1,NX;DO j=NZ,NZ
           wwtemp(i,m,j)=wwtemp(i,m,j-1) - (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                         - (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO 
  ENDIF

    Residual = 0
        DO m=1,NY;DO j=1,NZ;DO i=1,NX !i=1,nx
  
        Residual(i,m,j) = (   (1./dx)*(uutemp(i,m,j)-uutemp(i-1,m,j))   &
                             +(1./dy)*(vvtemp(i,m,j)-vvtemp(i,m-1,j))   &
                             +(1./dz)*(wwtemp(i,m,j)-wwtemp(i,m,j-1))  )
        ENDDO;ENDDO;ENDDO
   !(============================================================================  

        err=MAXVAL(ABS(Residual)) 

        IF (iSOR==1) THEN
            WRITE(18,*) err, iSOR
            WRITE(*,*) err, iSOR
        ENDIF 
 
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
 
        Resi=Residual*den_w/dt
!ppp=0;ppp2=0;ppp3=0;ppp4=0;ppp5=0


        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)
        !-----------------------------------------------
        CALL MatrixUpdate(Matrix,              & !OUT
                          ppp,kp,NXm,NYm,NZm)    !IN 
        Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------

IF(Multilevel==1 .or. err<errLevel1) GOTO 915
        
        !Coarse2 Grid (
        CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        
        CALL AveMatrix2SmallP(ppp2,NXm2,NYm2,NZm2,ppp,NXm,NYm,NZm)
        DO idum = 1, i_multi2
         CALL POISSON_Revised_Multigrid(ppp2,                          & !INOUT
                                        Resi2,kp2,NXm2,NYm2,NZm2)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp,NXm,NYm,NZm,ppp2,NXm2,NYm2,NZm2)
        !)

        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)
        !-----------------------------------------------
        CALL MatrixUpdate(Matrix,              & !OUT
                          ppp,kp,NXm,NYm,NZm)    !IN 
        Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------

IF(Multilevel==2 .or. err<errLevel2) GOTO 915

        !Coarse3 Grid (
        CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        
        CALL AveMatrix2SmallP(ppp3,NXm3,NYm3,NZm3,ppp2,NXm2,NYm2,NZm2)
        DO idum = 1, i_multi3
         CALL POISSON_Revised_Multigrid(ppp3,                          & !INOUT
                                        Resi3,kp3,NXm3,NYm3,NZm3)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp2,NXm2,NYm2,NZm2,ppp3,NXm3,NYm3,NZm3)
        !)

        !Coarse2 Grid (

        !CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        !CALL AveMatrix2SmallP(ppp2,NXm2,NYm2,NZm2,ppp,NXm,NYm,NZm)
        DO idum = 1, i_multi2
         CALL POISSON_Revised_Multigrid(ppp2,                          & !INOUT
                                       Resi2,kp2,NXm2,NYm2,NZm2)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp,NXm,NYm,NZm,ppp2,NXm2,NYm2,NZm2)
        !)

        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)
        !-----------------------------------------------
        CALL MatrixUpdate(Matrix,              & !OUT
                          ppp,kp,NXm,NYm,NZm)    !IN 
        Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------

IF(Multilevel==3 .or. err<errLevel3) GOTO 915

        !Coarse4 Grid (
        CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        CALL AveMatrix2Small(Resi4,NXm4,NYm4,NZm4,Resi3,NXm3,NYm3,NZm3)
        
        CALL AveMatrix2SmallP(ppp4,NXm4,NYm4,NZm4,ppp3,NXm3,NYm3,NZm3)
        DO idum = 1, i_multi4
         CALL POISSON_Revised_Multigrid(ppp4,                          & !INOUT
                                       Resi4,kp4,NXm4,NYm4,NZm4)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp3,NXm3,NYm3,NZm3,ppp4,NXm4,NYm4,NZm4)
        !)
        
        !Coarse3 Grid (
        !CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        !CALL AveMatrix2SmallP(ppp3,NXm3,NYm3,NZm3,ppp2,NXm2,NYm2,NZm2)
        DO idum = 1, i_multi3
         CALL POISSON_Revised_Multigrid(ppp3,                          & !INOUT
                                        Resi3,kp3,NXm3,NYm3,NZm3)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp2,NXm2,NYm2,NZm2,ppp3,NXm3,NYm3,NZm3)
        !)

        !Coarse2 Grid (
        !CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        !CALL AveMatrix2SmallP(ppp2,NXm2,NYm2,NZm2,ppp,NXm,NYm,NZm)
        DO idum = 1, i_multi2
        CALL POISSON_Revised_Multigrid(ppp2,                          & !INOUT
                                       Resi2,kp2,NXm2,NYm2,NZm2)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp,NXm,NYm,NZm,ppp2,NXm2,NYm2,NZm2)
        !)

        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)
        !-----------------------------------------------
        CALL MatrixUpdate(Matrix,              & !OUT
                          ppp,kp,NXm,NYm,NZm)    !IN 
        Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------        

IF(Multilevel==4 .or. err<errLevel4) GOTO 915

        !Coarse5 Grid (
        CALL AveMatrix2Small(Resi2,NXm2,NYm2,NZm2,Resi,NXm,NYm,NZm)
        CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        CALL AveMatrix2Small(Resi4,NXm4,NYm4,NZm4,Resi3,NXm3,NYm3,NZm3)       
        CALL AveMatrix2Small(Resi5,NXm5,NYm5,NZm5,Resi4,NXm4,NYm4,NZm4)
        
        CALL AveMatrix2SmallP(ppp5,NXm5,NYm5,NZm5,ppp4,NXm4,NYm4,NZm4)
        DO idum = 1, i_multi5
         CALL POISSON_Revised_Multigrid(ppp5,                          & !INOUT
                                        Resi5,kp5,NXm5,NYm5,NZm5)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp4,NXm4,NYm4,NZm4,ppp5,NXm5,NYm5,NZm5)
        !)

        !Coarse4 Grid (
        !CALL AveMatrix2Small(Resi4,NXm4,NYm4,NZm4,Resi3,NXm3,NYm3,NZm3)
        !CALL AveMatrix2SmallP(ppp4,NXm4,NYm4,NZm4,ppp3,NXm3,NYm3,NZm3)
        DO idum = 1, i_multi4
         CALL POISSON_Revised_Multigrid(ppp4,                          & !INOUT
                                        Resi4,kp4,NXm4,NYm4,NZm4)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp3,NXm3,NYm3,NZm3,ppp4,NXm4,NYm4,NZm4)
        !)
        
        !Coarse3 Grid (
        !CALL AveMatrix2Small(Resi3,NXm3,NYm3,NZm3,Resi2,NXm2,NYm2,NZm2)
        !CALL AveMatrix2SmallP(ppp3,NXm3,NYm3,NZm3,ppp2,NXm2,NYm2,NZm2)
        DO idum = 1, i_multi3
         CALL POISSON_Revised_Multigrid(ppp3,                          & !INOUT
                                        Resi3,kp3,NXm3,NYm3,NZm3)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp2,NXm2,NYm2,NZm2,ppp3,NXm3,NYm3,NZm3)
        !)        

        !Coarse2 Grid (
        DO idum = 1, i_multi2
         CALL POISSON_Revised_Multigrid(ppp2,                          & !INOUT
                                        Resi2,kp2,NXm2,NYm2,NZm2)        !IN
        ENDDO
        CALL CloneMatrix2Large(ppp,NXm,NYm,NZm,ppp2,NXm2,NYm2,NZm2)
        !)
              
        !Fine Grid (
        DO idum = 1, i_multi1
         CALL POISSON_Revised_Multigrid(ppp,                      & !INOUT
                                        Resi,kp,NXm,NYm,NZm)        !IN
        ENDDO
        !)

        !-----------------------------------------------
        !CALL MatrixUpdate(Matrix,              & !OUT
        !                  ppp,kp,NXm,NYm,NZm)    !IN 
        !Resi = Resi - Matrix                                               
        PPtemp=PPtemp+ppp; PPc=PPc+ppp
        !-----------------------------------------------
        !Residual = Resi*dt/den_w

915     CONTINUE
   !)============================================================================ 
        !(4) Correct uu,vv,ww
        DO m=1,NY;DO i=1,NX-1;DO j=1,NZ
           uutemp(i,m,j)=uuIN(i,m,j) -dt/dx*(PPc(i+1,m,j)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO
        DO m=1,NY-1;DO i=1,NX;DO j=1,NZ
           vvtemp(i,m,j)=vvIN(i,m,j) -dt/dy*(PPc(i,m+1,j)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO        
        DO m=1,NY;DO i=1,NX;DO j=1,NZ-1 
           wwtemp(i,m,j)=wwIN(i,m,j)-dt/dz*(PPc(i,m,j+1)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO
  If (SurfaceCorrection==0) THEN
        DO m=1,NY;DO i=1,NX;DO j=NZ,NZ
           wwtemp(i,m,j)=wwtemp(i,m,j-1) - (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                         - (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO 
  ELSE 
        DO m=1,NY;DO i=1,NX;DO j=NZ,NZ 
           wwtemp(i,m,j)=wwIN(i,m,j)-dt/dz*(PPc(i,m,j+1)-PPc(i,m,j))/den_w
        ENDDO;ENDDO;ENDDO 
  ENDIF      
   !)============================================================================ 
    ENDDO LOOP_1

    
    
    IF( ANH == 2) THEN !Nonhydrostatic with Hydrostatic Correction -> Adjustable NonHydrostatic

        DO j=1,NZ !-1
        DO m=1,NY
        DO i=1,NX   
             wwtemp(i,m,j)   = wwtemp(i,m,j-1) -  (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                               -  (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO
        
    ENDIF

ELSEIF (ANH == 0) THEN ! Fully Hydrostatic

        DO j=1,NZ !-1
        DO m=1,NY
        DO i=1,NX   
             wwtemp(i,m,j)   = wwtemp(i,m,j-1) -  (uutemp(i,m,j)-uutemp(i-1,m,j))*dz/dx   &
                                               -  (vvtemp(i,m,j)-vvtemp(i,m-1,j))*dz/dy
        ENDDO;ENDDO;ENDDO

ENDIF    

ppOUT=PPtemp
uuOUT=uutemp
vvOUT=vvtemp
wwOUT=wwtemp

CONTINUE

ENDSUBROUTINE SIMPLE_R_S_Multigrid_incre