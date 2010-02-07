SUBROUTINE UPREDICT_CENTRAL   (uustar,vvstar,wwstar,                    & !OUT
                              uu_old,vv_old,ww_old,uu_m,vv_m,ww_m,      & !IN
                              pp,BaCX,BaCY,BaT,den)               !IN  
!_______________________________________________________________________________________

USE COMDAT


REAL(8)    ,DIMENSION(0:nx+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: pp
REAL(8)    ,DIMENSION(0:NX  ,0:NY  ,0:NZ+1) ,INTENT(IN)    :: BaCX,BaCY
REAL(8)    ,DIMENSION(1:NX  ,1:NY  ,2)      ,INTENT(IN)    :: BaT
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ,INTENT(IN)    :: den
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(IN)    :: uu_m,uu_old
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(IN)    :: vv_m,vv_old
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(IN)    :: ww_m,ww_old

REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(OUT)   :: uustar
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(OUT)   :: vvstar
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(OUT)   :: wwstar

REAL(8)                                        :: viscX,viscZ
REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1)    :: uuA,uuD,uutemp
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1)    :: vvA,vvD,vvtemp
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  )    :: wwA,wwD,wwtemp
!________________________________________________________________
viscX=visc
viscZ=visc
uutemp=uu_m
vvtemp=vv_m
wwtemp=ww_m

 CALL       BC_wall(uutemp,vvtemp,wwtemp)
 CALL       VelocityBC(uutemp,vvtemp,wwtemp)



!-------------------------------Advection----------------------------------------------------------
!------------------------uuA---------------------
DO m=1,NY
DO i=1,NX-1 !i=1,NX-1
DO j=1,NZ !-1 ! 05/02/15

 uuA(i,m,j)=(1./dx)*(((uutemp(i+1,m,j)+uutemp(i,m,j))/2)**2-((uutemp(i,m,j)+uutemp(i-1,m,j))/2)**2)  &
  &      +(1./dy)*( ((uutemp(i,m+1,j)+uutemp(i,m,j))/2)*((vvtemp(i+1,m,j)+vvtemp(i,m,j))/2)      &
  &               -((uutemp(i,m,j)+uutemp(i,m-1,j))/2)*((vvtemp(i+1,m-1,j)+vvtemp(i,m-1,j))/2))  &
  &    +(1./dz)*( ((uutemp(i,m,j+1)+uutemp(i,m,j))/2)*((wwtemp(i+1,m,j)+wwtemp(i,m,j))/2)      &
  &               -((uutemp(i,m,j)+uutemp(i,m,j-1))/2)*((wwtemp(i+1,m,j-1)+wwtemp(i,m,j-1))/2))

ENDDO
ENDDO
ENDDO
!------------------------vvA---------------------
DO m=1,NY-1
DO i=1,NX !i=1,NX
DO j=1,NZ !-1 ! 05/02/15

 vvA(i,m,j)=(1./dy)*(((vvtemp(i,m+1,j)+vvtemp(i,m,j))/2)**2-((vvtemp(i,m,j)+vvtemp(i,m-1,j))/2)**2)  &
  &      +(1./dx)*(  ((vvtemp(i+1,m,j)+vvtemp(i,m,j))/2)*((uutemp(i,m+1,j)+uutemp(i,m,j))/2)       &
  &               -((vvtemp(i,m,j)+vvtemp(i-1,m,j))/2)*((uutemp(i-1,m+1,j)+uutemp(i-1,m,j))/2))  &
  &    +(1./dz)*( ((vvtemp(i,m,j+1)+vvtemp(i,m,j))/2)*((wwtemp(i,m+1,j)+wwtemp(i,m,j))/2)      &
  &               -((vvtemp(i,m,j)+vvtemp(i,m,j-1))/2)*((wwtemp(i,m+1,j-1)+wwtemp(i,m,j-1))/2))

ENDDO
ENDDO
ENDDO
!------------------------wwA---------------------
DO m=1,NY
DO i=1,NX !i=1,NX
DO j=1,NZ -1
  
 wwA(i,m,j)=(1./dz)*(((wwtemp(i,m,j+1)+wwtemp(i,m,j))/2)**2-((wwtemp(i,m,j)+wwtemp(i,m,j-1))/2)**2)    &
  &      +(1./dy)*( ((vvtemp(i,m,j+1)+vvtemp(i,m,j))/2)*((wwtemp(i,m,j)+wwtemp(i,m+1,j))/2)            &
  &               -((vvtemp(i,m-1,j+1)+vvtemp(i,m-1,j))/2)*((wwtemp(i,m,j)+wwtemp(i,m-1,j))/2))      &
  &      +(1./dx)*( ((uutemp(i,m,j+1)+uutemp(i,m,j))/2)*((wwtemp(i,m,j)+wwtemp(i+1,m,j))/2)            &
  &               -((uutemp(i-1,m,j+1)+uutemp(i-1,m,j))/2)*((wwtemp(i,m,j)+wwtemp(i-1,m,j))/2))

ENDDO
ENDDO
ENDDO

DO m=1,NY
DO i=1,NX !i=1,NX
DO j=NZ,NZ
  
 wwA(i,m,j)=(1./dz)*(((wwtemp(i,m,j)+wwtemp(i,m,j))/2)**2-((wwtemp(i,m,j)+wwtemp(i,m,j-1))/2)**2)    &
  &      +(1./dy)*( ((vvtemp(i,m,j+1)+vvtemp(i,m,j))/2)*((wwtemp(i,m,j)+wwtemp(i,m+1,j))/2)            &
  &               -((vvtemp(i,m-1,j+1)+vvtemp(i,m-1,j))/2)*((wwtemp(i,m,j)+wwtemp(i,m-1,j))/2))      &
  &      +(1./dx)*( ((uutemp(i,m,j+1)+uutemp(i,m,j))/2)*((wwtemp(i,m,j)+wwtemp(i+1,m,j))/2)            &
  &               -((uutemp(i-1,m,j+1)+uutemp(i-1,m,j))/2)*((wwtemp(i,m,j)+wwtemp(i-1,m,j))/2))

ENDDO
ENDDO
ENDDO

!Ghost Points for d2u/dy2
!----------------------------------
!DO i=1, NX-1
! uutemp(i,NZ+1)  = 0            + 1./3.*uutemp(i,NZ-1)-2*uutemp(i,NZ) 
! uutemp(i,0)     = 8./3.*uuwall + 1./3.*uutemp(i,2)   -2*uutemp(i,1)
!ENDDO
!DO j=1, NZ-1
! wwtemp(NX+1,j) = 0             + 1./3.*wwtemp(NX-1,j)-2*wwtemp(NX,j)
! wwtemp(0,j)    = 0             + 1./3.*wwtemp(2,j)   -2*wwtemp(1,j)
!ENDDO

!CALL       BC_SURFACE_R(uutemp,vvtemp,wwtemp,              & !INOUT
!                      w_kineticIN,topcell)                        !IN common 


DO m=1,NY
DO i=1,NX-1 !i=1,NX-1
DO j=1,NZ !-1 ! 05/02/15
 uuD(i,m,j)=   1./dx*(   viscX*(uutemp(i+1,m,j)-uutemp(i,m,j))/dx       &
                       + viscX*(uutemp(i-1,m,j)-uutemp(i,m,j))/dx   )   &
            +  1./dy*(   viscX*(uutemp(i,m+1,j)-uutemp(i,m,j))/dy       &
                       + viscX*(uutemp(i,m-1,j)-uutemp(i,m,j))/dy   )   &
            +  1./dz*(   viscZ*(uutemp(i,m,j+1)-uutemp(i,m,j))/dz       &
                       + viscZ*(uutemp(i,m,j-1)-uutemp(i,m,j))/dz   ) 
ENDDO
ENDDO
ENDDO


DO m=1,NY-1
DO i=1,NX !i=1,NX
DO j=1,NZ !-1 ! 05/02/15
 vvD(i,m,j)=   1./dx*(   viscX*(vvtemp(i+1,m,j)-vvtemp(i,m,j))/dx       &
                       + viscX*(vvtemp(i-1,m,j)-vvtemp(i,m,j))/dx   )   &
            +  1./dy*(   viscX*(vvtemp(i,m+1,j)-vvtemp(i,m,j))/dy       &
                       + viscX*(vvtemp(i,m-1,j)-vvtemp(i,m,j))/dy   )   &
            +  1./dz*(   viscZ*(vvtemp(i,m,j+1)-vvtemp(i,m,j))/dz       &
                       + viscZ*(vvtemp(i,m,j-1)-vvtemp(i,m,j))/dz   ) 
ENDDO
ENDDO
ENDDO

DO m=1,NY
DO i=1,NX !i=1,NX
DO j=1,NZ -1

 wwD(i,m,j)=   1./dx*(   viscX*(wwtemp(i+1,m,j)-wwtemp(i,m,j))/dx       &
                     +   viscX*(wwtemp(i-1,m,j)-wwtemp(i,m,j))/dx   )   &
            +  1./dy*(   viscX*(wwtemp(i,m+1,j)-wwtemp(i,m,j))/dy       &
                     +   viscX*(wwtemp(i,m-1,j)-wwtemp(i,m,j))/dy   )   &
            +  1./dz*(   viscZ*(wwtemp(i,m,j+1)-wwtemp(i,m,j))/dz       &
                     +   viscZ*(wwtemp(i,m,j-1)-wwtemp(i,m,j))/dz   ) 

ENDDO
ENDDO
ENDDO

DO m=1,NY
DO i=1,NX !i=1,NX
DO j=NZ,NZ

 wwD(i,m,j)=   1./dx*(   viscX*(wwtemp(i+1,m,j)-wwtemp(i,m,j))/dx       &
                     +   viscX*(wwtemp(i-1,m,j)-wwtemp(i,m,j))/dx   )   &
            +  1./dy*(   viscX*(wwtemp(i,m+1,j)-wwtemp(i,m,j))/dy       &
                     +   viscX*(wwtemp(i,m-1,j)-wwtemp(i,m,j))/dy   )   &
            +  1./dz*(   viscZ*(wwtemp(i,m,j  )-wwtemp(i,m,j))/dz       &
                     +   viscZ*(wwtemp(i,m,j-1)-wwtemp(i,m,j))/dz   ) 
ENDDO
ENDDO
ENDDO

DO m=1,NY
DO i=1,NX-1 !i=1,NX-1
DO j=1,NZ !-1  ! 05/02/15

 uustar(i,m,j)=uu_old(i,m,j)+dt*(-uuA(i,m,j)+uuD(i,m,j)-BaT(i,m,1)-WBaC*BaCX(i,m,j))    &
  &         -dt/dx*(pp(i+1,m,j)-pp(i,m,j))/(0.5*(den(i+1,m,j)+den(i,m,j)))  

                                         
ENDDO
ENDDO
ENDDO

DO m=1,NY-1
DO i=1,NX !i=1,NX
DO j=1,NZ !-1 ! 05/02/15

 vvstar(i,m,j)=vv_old(i,m,j)+dt*(-vvA(i,m,j)+vvD(i,m,j)-BaT(i,m,2)-WBaC*BaCY(i,m,j))    &
  &         -dt/dy*(pp(i,m+1,j)-pp(i,m,j))/(0.5*(den(i,m+1,j)+den(i,m,j)))  

                                            
ENDDO
ENDDO
ENDDO
  
DO m=1,NY
DO i=1,NX !i=1,NX
DO j=1,NZ -1

 wwstar(i,m,j)=ww_old(i,m,j)+dt*(-wwA(i,m,j)+wwD(i,m,j))   -dt*(1.-WBaC)*( (0.5*(den(i,m,j+1)+den(i,m,j))-den_w)/den_w*gravity )    & 
  &        -dt/dz*(pp(i,m,j+1)-pp(i,m,j))/(0.5*(den(i,m,j+1)+den(i,m,j)))

ENDDO
ENDDO
ENDDO

DO m=1,NY
DO i=1,NX !i=1,NX
DO j=NZ,NZ

 wwstar(i,m,j)=ww_old(i,m,j)+dt*(-wwA(i,m,j)+wwD(i,m,j))   -dt*(1.-WBaC)*( (0.5*(den(i,m,j+1)+den(i,m,j+1))-den_w)/den_w*gravity )    & 
  &        -dt/dz*(pp(i,m,j+1)-pp(i,m,j))/(0.5*(den(i,m,j+1)+den(i,m,j+1)))

ENDDO
ENDDO
ENDDO

CONTINUE

ENDSUBROUTINE UPREDICT_CENTRAL