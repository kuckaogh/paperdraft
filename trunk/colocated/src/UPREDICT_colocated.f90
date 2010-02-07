SUBROUTINE UPREDICT_colocated(uu_OUT,vv_OUT,ww_OUT,                     & !OUT
                              uu_old,vv_old,ww_old,uu_m,vv_m,ww_m,      & !IN
                              pp,BaT,den,NXm,NYm,NZm,wall)                !IN  
!_______________________________________________________________________________________

USE COMDAT

INTEGER(4),INTENT(IN)           :: NXm,NYm,NZm,wall
!REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1),INTENT(IN)  :: pp
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)             :: pp
!REAL(8)   ,DIMENSION(0:NXm  ,0:NYm  ,0:NZm+1),INTENT(IN)  :: BaCX,BaCY
REAL(8)   ,DIMENSION(1:NXm  ,1:NYm  ,2)      ,INTENT(IN)  :: BaT
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1),INTENT(IN)  :: den
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1),INTENT(IN)  :: uu_m,uu_old
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1),INTENT(IN)  :: vv_m,vv_old
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1),INTENT(IN)  :: ww_m,ww_old
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1),INTENT(OUT) :: uu_OUT
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1),INTENT(OUT) :: vv_OUT
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1),INTENT(OUT) :: ww_OUT
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: uutemp
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: vvtemp
REAL(8)   ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)  :: wwtemp
REAL(8)                                        :: viscX,viscZ
REAL(8)   ,DIMENSION(1:NXm+1,1:NYm+1,1:NZm+1)  :: uuA,uuD
REAL(8)   ,DIMENSION(1:NXm+1,1:NYm+1,1:NZm+1)  :: vvA,vvD
REAL(8)   ,DIMENSION(1:NXm+1,1:NYm+1,1:NZm+1)  :: wwA,wwD
!________________________________________________________________
viscX=visc
viscZ=visc
uutemp=uu_m
vvtemp=vv_m
wwtemp=ww_m

 CALL       BC_wall_colocated   (uutemp,vvtemp,wwtemp,NXm,NYm,NZm,wall) 
 CALL       VelocityBC_colocated(uutemp,vvtemp,wwtemp,NXm,NYm,NZm) 

!-------------------------------Advection----------------------------------------------------------
!------------------------uuA---------------------
DO m=1,NYm;DO i=1,NXm;DO j=1,NZm 

uuA(i,m,j) =(                                                                    &
               dy*dz*(    (0.5*( uutemp(i+1,m,j)+uutemp(i,m,j) ))                &
                         *(0.5*( uutemp(i+1,m,j)+uutemp(i,m,j) ))                &   
                       -  (0.5*( uutemp(i-1,m,j)+uutemp(i,m,j) ))                &
                         *(0.5*( uutemp(i-1,m,j)+uutemp(i,m,j) ))                &
                     )                                                           &
             + dx*dz*(    (0.5*( uutemp(i,m+1,j)+uutemp(i,m,j) ))                &
                         *(0.5*( vvtemp(i,m+1,j)+vvtemp(i,m,j) ))                &   
                       -  (0.5*( uutemp(i,m-1,j)+uutemp(i,m,j) ))                &
                         *(0.5*( vvtemp(i,m-1,j)+vvtemp(i,m,j) ))                &
                     )                                                           &                       
             + dx*dy*(    (0.5*( uutemp(i,m,j+1)+uutemp(i,m,j) ))                &
                         *(0.5*( wwtemp(i,m,j+1)+wwtemp(i,m,j) ))                &   
                       -  (0.5*( uutemp(i,m,j-1)+uutemp(i,m,j) ))                &
                         *(0.5*( wwtemp(i,m,j-1)+wwtemp(i,m,j) ))                &
                     )                                                           &
             )*(1./(dx*dy*dz))                                                          
ENDDO;ENDDO;ENDDO
!------------------------vvA---------------------
DO m=1,NYm;DO i=1,NXm;DO j=1,NZm 

vvA(i,m,j) =(                                                                    &
               dy*dz*(    (0.5*( vvtemp(i+1,m,j)+vvtemp(i,m,j) ))                &
                         *(0.5*( uutemp(i+1,m,j)+uutemp(i,m,j) ))                &   
                       -  (0.5*( vvtemp(i-1,m,j)+vvtemp(i,m,j) ))                &
                         *(0.5*( uutemp(i-1,m,j)+uutemp(i,m,j) ))                &
                     )                                                           &
             + dx*dz*(    (0.5*( vvtemp(i,m+1,j)+vvtemp(i,m,j) ))                &
                         *(0.5*( vvtemp(i,m+1,j)+vvtemp(i,m,j) ))                &   
                       -  (0.5*( vvtemp(i,m-1,j)+vvtemp(i,m,j) ))                &
                         *(0.5*( vvtemp(i,m-1,j)+vvtemp(i,m,j) ))                &
                     )                                                           &                       
             + dx*dy*(    (0.5*( vvtemp(i,m,j+1)+vvtemp(i,m,j) ))                &
                         *(0.5*( wwtemp(i,m,j+1)+wwtemp(i,m,j) ))                &   
                       -  (0.5*( vvtemp(i,m,j-1)+vvtemp(i,m,j) ))                &
                         *(0.5*( wwtemp(i,m,j-1)+wwtemp(i,m,j) ))                &
                     )                                                           &
             )*(1./(dx*dy*dz))                                                          
ENDDO;ENDDO;ENDDO
!------------------------wwA---------------------
DO m=1,NYm;DO i=1,NXm;DO j=1,NZm 

wwA(i,m,j) =(                                                                    &
               dy*dz*(    (0.5*( wwtemp(i+1,m,j)+wwtemp(i,m,j) ))                &
                         *(0.5*( uutemp(i+1,m,j)+uutemp(i,m,j) ))                &   
                       -  (0.5*( wwtemp(i-1,m,j)+wwtemp(i,m,j) ))                &
                         *(0.5*( uutemp(i-1,m,j)+uutemp(i,m,j) ))                &
                     )                                                           &
             + dx*dz*(    (0.5*( wwtemp(i,m+1,j)+wwtemp(i,m,j) ))                &
                         *(0.5*( vvtemp(i,m+1,j)+vvtemp(i,m,j) ))                &   
                       -  (0.5*( wwtemp(i,m-1,j)+wwtemp(i,m,j) ))                &
                         *(0.5*( vvtemp(i,m-1,j)+vvtemp(i,m,j) ))                &
                     )                                                           &                       
             + dx*dy*(    (0.5*( wwtemp(i,m,j+1)+wwtemp(i,m,j) ))                &
                         *(0.5*( wwtemp(i,m,j+1)+wwtemp(i,m,j) ))                &   
                       -  (0.5*( wwtemp(i,m,j-1)+wwtemp(i,m,j) ))                &
                         *(0.5*( wwtemp(i,m,j-1)+wwtemp(i,m,j) ))                &
                     )                                                           &
             )*(1./(dx*dy*dz))                                                          
ENDDO;ENDDO;ENDDO

DO m=1,NYm;DO i=1,NXm;DO j=1,NZm 

 uuD(i,m,j)=   1./dx*(   viscX*(uutemp(i+1,m,j)-uutemp(i,m,j))/dx       &
                       + viscX*(uutemp(i-1,m,j)-uutemp(i,m,j))/dx   )   &
            +  1./dy*(   viscX*(uutemp(i,m+1,j)-uutemp(i,m,j))/dy       &
                       + viscX*(uutemp(i,m-1,j)-uutemp(i,m,j))/dy   )   &
            +  1./dz*(   viscZ*(uutemp(i,m,j+1)-uutemp(i,m,j))/dz       &
                       + viscZ*(uutemp(i,m,j-1)-uutemp(i,m,j))/dz   ) 
ENDDO;ENDDO;ENDDO

DO m=1,NYm;DO i=1,NXm;DO j=1,NZm 

 vvD(i,m,j)=   1./dx*(   viscX*(vvtemp(i+1,m,j)-vvtemp(i,m,j))/dx       &
                       + viscX*(vvtemp(i-1,m,j)-vvtemp(i,m,j))/dx   )   &
            +  1./dy*(   viscX*(vvtemp(i,m+1,j)-vvtemp(i,m,j))/dy       &
                       + viscX*(vvtemp(i,m-1,j)-vvtemp(i,m,j))/dy   )   &
            +  1./dz*(   viscZ*(vvtemp(i,m,j+1)-vvtemp(i,m,j))/dz       &
                       + viscZ*(vvtemp(i,m,j-1)-vvtemp(i,m,j))/dz   ) 
ENDDO;ENDDO;ENDDO

DO m=1,NYm;DO i=1,NXm;DO j=1,NZm 

 wwD(i,m,j)=   1./dx*(   viscX*(wwtemp(i+1,m,j)-wwtemp(i,m,j))/dx       &
                     +   viscX*(wwtemp(i-1,m,j)-wwtemp(i,m,j))/dx   )   &
            +  1./dy*(   viscX*(wwtemp(i,m+1,j)-wwtemp(i,m,j))/dy       &
                     +   viscX*(wwtemp(i,m-1,j)-wwtemp(i,m,j))/dy   )   &
            +  1./dz*(   viscZ*(wwtemp(i,m,j+1)-wwtemp(i,m,j))/dz       &
                     +   viscZ*(wwtemp(i,m,j-1)-wwtemp(i,m,j))/dz   ) 
ENDDO;ENDDO;ENDDO

PP(0,:,:)     = PP(1,:,:)
PP(NXm+1,:,:) = PP(NXm,:,:)
PP(:,0,:)     = PP(:,1,:)
PP(:,NYm+1,:) = PP(:,NYm,:)
PP(:,:,0)     = PP(:,:,1)
PP(:,:,NZm+1) = PP(:,:,NZm)

DO m=1,NYm;DO i=1,NXm;DO j=1,NZm  

 uu_OUT(i,m,j)=uu_old(i,m,j)+dt*(-uuA(i,m,j)+uuD(i,m,j)-BaT(i,m,1) )    &
  &         -dt/dx*0.5*(pp(i+1,m,j)-pp(i-1,m,j))/den(i,m,j)  
                                         
ENDDO;ENDDO;ENDDO

DO m=1,NYm;DO i=1,NXm;DO j=1,NZm 

 vv_OUT(i,m,j)=vv_old(i,m,j)+dt*(-vvA(i,m,j)+vvD(i,m,j)-BaT(i,m,2) )    &
  &         -dt/dy*0.5*(pp(i,m+1,j)-pp(i,m-1,j))/den(i,m,j)   
                                            
ENDDO;ENDDO;ENDDO
  
DO m=1,NYm;DO i=1,NXm;DO j=1,NZm 

 ww_OUT(i,m,j)=ww_old(i,m,j)+dt*(-wwA(i,m,j)+wwD(i,m,j))-dt*( ( den(i,m,j)-den_w )/den_w*gravity  )   & 
  &         -dt/dz*0.5*(pp(i,m,j+1)-pp(i,m,j-1))/den(i,m,j) 

ENDDO;ENDDO;ENDDO


CONTINUE

 CALL       BC_wall_colocated   (uu_OUT,vv_OUT,ww_OUT,NXm,NYm,NZm,wall) 
 CALL       VelocityBC_colocated(uu_OUT,vv_OUT,ww_OUT,NXm,NYm,NZm) 
 
ENDSUBROUTINE UPREDICT_colocated