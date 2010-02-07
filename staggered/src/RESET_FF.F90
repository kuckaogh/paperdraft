SUBROUTINE RESET_FF(total_ref, total,ff,ff_i)
USE COMDAT
REAL(8)                                    ,INTENT(IN)    :: total_ref,total
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1),INTENT(INOUT) :: ff
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1),INTENT(IN)    :: ff_i
REAL(8)                                                   :: ffd,tr

tr=(total/total_ref)**conservation_power

 DO i=1,NX
 DO m=1,NY
 DO j=1,NZ !-2 !1,NZ
   !ff(i,m,j) = ff(i,m,j) + (ff(i,m,j)-ff_0(i,m,j))*(1.0d0-total/total_ref)
   ffd = ff(i,m,j)-ff_i(i,m,j)
   ff(i,m,j) = ff_i(i,m,j) + ffd *( 1.0  - sign(1.0,ffd)*( tr - 1.0 )  )
 ENDDO
 ENDDO
 ENDDO
 

!ff(1:NX,1:NY,1:NZ) = ff(1:NX,1:NY,1:NZ) + (ff(1:NX,1:NY,1:NZ)-ff_0(1:NX,1:NY,1:NZ))*(1.0d0-totalsalin/totalsalin_ori)

ENDSUBROUTINE RESET_FF
