SUBROUTINE OUTPUT_colocated(voll,volsum,uu,vv,ww,pp,ff,den,NXm,NYm,NZm)

USE COMDAT
INTEGER(4)   ,INTENT(IN)           :: NXm,NYm,NZm
REAL(8)      ,DIMENSION(0:NXm+1,0:NYm+1)        ,INTENT(IN) :: volsum
REAL(8)      ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN) :: den
REAL(8)      ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)             :: pp,ff
REAL(8)      ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1)             :: x, y, z
REAL(8)      ,DIMENSION(0:NXm+1,0:NYm+1,0:NZm+1) ,INTENT(IN) :: voll
REAL(8)      ,DIMENSION(0:NXm+1,0:NYm+1)                    :: volsum_new
REAL(8)       uu(0:NXm,0:NYm+1,0:NZm+1), ww(0:NXm+1,0:NYm+1,0:NZm), vv(0:NXm+1,0:NYm,0:NZm+1)

x=0
y=0
z=0
!NZ=20
!Plot Average
!-----------------------------------------------------------------------

!volsum_new=volsum*dz-h0
volsum_new=volsum-h0

WRITE(15,100)

WRITE(15,101) step*dt, NXm, NYm, NZm

  
WRITE(15,102) (((den(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( pp(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( ff(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( voll(i,m,j), i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( (uu(i-1,m,j)+uu(i,m,j))/2, i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( (vv(i,m-1,j)+vv(i,m,j))/2, i=1,nxm),m=1,NYm),j=1,NZm)
WRITE(15,102) ((( (ww(i,m,j-1)+ww(i,m,j))/2, i=1,nxm),m=1,NYm),j=1,NZm)

DO j=1,NZm ;DO m=1,NYm ;DO i=1,NXm
WRITE(15,102) 0.5*( (uu(i,m,j+1)-uu(i,m,j-1))/dz - (ww(i+1,m,j)-ww(i-1,m,j))/dx ) 
ENDDO;ENDDO;ENDDO

100   format('VARIABLES = x, y, z, den, p, f, lvol, u, v, w, curl')
101   format(' zone t="',f9.4,'", i=',i5, ', j=',i5,', k=',i5,', f=block, VARSHARELIST=([1,2,3]=1)')
102   format(E15.6)
!NZ=25
ENDSUBROUTINE OUTPUT_colocated