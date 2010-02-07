SUBROUTINE BaroTropic_colocated(BaT,                       &   !OUT                                 
                      volsum_1,NXm ,NYm)                  !IN                


USE COMDAT

INTEGER(4)   ,INTENT(IN)           :: NXm,NYm
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1) ,INTENT(IN)    :: volsum_1
REAL(8)    ,DIMENSION(0:NXm+1,0:NYm+1)                :: volsum
REAL(8)    ,DIMENSION(1:NXm,1:NYm,2)   ,INTENT(OUT)   :: BaT


BaT  = 0
volsum=(1/dz)*volsum_1

volsum(0,:)=volsum(1,:)
volsum(NXm+1,:)=volsum(NXm,:)
volsum(:,0)=volsum(:,1)
volsum(:,NYm+1)=volsum(:,NYm)

!******************************************
DO m=1,NYm
DO i=1,NXm !i=1,NX-1

 BaT(i,m,1)= gravity*0.5*(volsum(i+1,m)-volsum(i-1,m))*dz/dx
 !BaT(i,m,1)= gravity*(volsum(i+1,m)-volsum(i,m))*dz*dy/dx

ENDDO
ENDDO
DO m=1,NYm
DO i=1,NXm !i=1,NX

 BaT(i,m,2)= gravity*0.5*(volsum(i,m+1)-volsum(i,m-1))*dz/dy
! BaT(i,m,2)= gravity*(volsum(i,m+1)-volsum(i,m))*dz*dx/dy

ENDDO
ENDDO
!******************************************
CONTINUE
ENDSUBROUTINE BaroTropic_colocated

