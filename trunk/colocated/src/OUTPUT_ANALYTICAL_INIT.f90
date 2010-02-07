SUBROUTINE OUTPUT_ANALYTICAL_INIT(NX,NY,NZ,dx,dz,dt,step)
!SUBROUTINE OUTPUT_ANALYTICAL_INIT(h0,ff,vol,volsum,den,NX,NZ,dx,dz,uu,ww,dt,step)

IMPLICIT NONE

INTEGER NX,NY,NZ,step
INTEGER i,m,j

REAL(8)                         ,INTENT(IN) ::dt,dx,dz
REAL(8) ,DIMENSION(0:NX+1,0:NY+1,0:NZ+1) ::x, z,x_new, z_new

REAL(8) h0,xs(0:NX)

REAL(8), PARAMETER   :: pi=3.141592653589793238462643383279502884197169399375 
REAL(8)              :: Kapa,Omega,Length,Am,Time,gravity,Density
REAL(8) ,DIMENSION(1:NX,1:NY,1:NZ) :: u,w,p,Eta
x=0
z=0
x_new=0
z_new=0
xs = 0
u=0;w=0;p=0;Eta=0

!Plot Average
!-----------------------------------------------------------------------
DO i=1,NX
DO m=1,NY
DO j=1,NZ
          x(i,m,j) = dx*(i-1)+dx/2                    ! Assigning x, y for plotting
          z(i,m,j) = dz*(j-1)+dz/2
ENDDO
ENDDO
ENDDO

Density=1000.0
gravity=9.8
Time = step*dt
h0 = 10.0
Length = 10.0
Am=0.1
Kapa =pi/Length
Omega = (gravity*Kapa*TANH(Kapa*h0))**0.5
z=z-h0


DO i=1,NX
DO m=1,NY
DO j=1,NZ

u(i,m,j)=Omega*Am*SIN(Kapa*x(i,m,j))*SIN(Omega*Time)*COSH(Kapa*(z(i,m,j)+h0))/SINH(Kapa*h0)
w(i,m,j)=-Omega*Am*COS(Kapa*x(i,m,j))*SIN(Omega*Time)*SINH(Kapa*(z(i,m,j)+h0))/SINH(Kapa*h0)
Eta(i,m,j)=Am*COS(Kapa*x(i,m,j))*COS(Omega*Time)

ENDDO
ENDDO
ENDDO

DO i=1,NX
DO m=1,NY
DO j=1,NZ
p(i,m,j)=-density*gravity*Eta(i,m,j) + density*gravity*Am* COSH(Kapa*(z(i,m,j)+h0))/COSH(Kapa*h0)*COS(Kapa*x(i,m,j))*COS(Omega*Time)
ENDDO
ENDDO
ENDDO



!NZ=20
WRITE(153,100)
WRITE(153,101) step*dt, NX, NY, NZ
WRITE(153,102) ((( x(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(153,102) ((( z(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(153,102) ((( p(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(153,102) ((( u(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
WRITE(153,102) ((( w(i,m,j), i=1,nx),m=1,NY),j=1,NZ) 
WRITE(153,102) ((( Eta(i,m,j), i=1,nx),m=1,NY),j=1,NZ)
!NZ=25
100   format('VARIABLES = x, z, p, u, w, height')
101   format(' zone t="',f9.4,'", i=',i5, ', j=',i5,', k=',i5,', f=block')
102   format(E15.6)


!---------------------------------------------------------------



! Surface Plot
!---------------------------------------
DO i=1,NX
          xs(i) = dx*(i-1)+dx/2                    
ENDDO

WRITE(156,109)
WRITE(156,110) step*dt, NX
WRITE(156,*) 'DATAPACKING=POINT'
DO i=1,NX
 WRITE(156,111) xs(i), xs(i)/h0, Am*COS(Kapa*xs(i))*COS(Omega*Time)
ENDDO

109   format('VARIABLES = "x","x/h", "surface"')
110   format(' zone t="',f9.4,'", i=',i5, ', j=1, k=1, zonetype=ordered,  DATAPACKING=POINT')
111   format(E15.6,E15.6,E15.6,E15.6)

ENDSUBROUTINE OUTPUT_ANALYTICAL_INIT