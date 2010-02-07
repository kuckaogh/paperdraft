SUBROUTINE CloneMatrix2Large(M_L_OUT,NX_L,NY_L,NZ_L,M_s_IN,NX_s,NY_s,NZ_s)

INTEGER(4)                                   ,INTENT(IN) :: NX_s,NY_s,NZ_s,NX_L,NY_L,NZ_L
REAL(8),DIMENSION(0:NX_s+1,0:NY_s+1,0:NZ_s+1),INTENT(IN) :: M_s_IN
REAL(8),DIMENSION(0:NX_L+1,0:NY_L+1,0:NZ_L+1),INTENT(OUT):: M_L_OUT
INTEGER(4)                                               :: i,m,j,i2,m2,j2
REAL(8),DIMENSION(0:NX_s+1,0:NY_s+1,0:NZ_s+1)            :: M_s
REAL(8),DIMENSION(0:NX_L+1,0:NY_L+1,0:NZ_L+1)            :: M_L

M_s=M_s_IN
M_L=0



IF (NX_s*2==NX_L) THEN

    IF (NY_s==1) THEN
     DO i=1,NX_L; DO j=1,NZ_L
      i2 = CEILING(REAL(i)/2); j2 = CEILING(REAL(j)/2)
      M_L(i,1,j) = M_s(i2,1,j2) 
     ENDDO; ENDDO
    ELSE
     DO i=1,NX_L; DO m=1,NY_L; DO j=1,NZ_L
      i2 = CEILING(REAL(i)/2);  m2 = CEILING(REAL(m)/2);  j2 = CEILING(REAL(j)/2)
      M_L(i,m,j) = M_s(i2,m2,j2) 
     ENDDO; ENDDO; ENDDO
    ENDIF

ELSE

    IF (NY_s==1) THEN
    
    M_s(0     ,1,1:NZ_s)=M_s(1   ,1,1:NZ_s)
    M_s(NX_s+1,1,1:NZ_s)=M_s(NX_s,1,1:NZ_s)
    
    M_s(1:NX_s,1,0     )=M_s(1:NX_s,1,1   )
    M_s(1:NX_s,1,NZ_s+1)=M_s(1:NX_s,1,NZ_s)
    
    M_s(0     ,1,0     )=    M_S(1,   1,1   )
    M_s(NX_s+1,1,0     )=    M_s(NX_s,1,1   )
    M_s(0     ,1,NZ_s+1)=    M_s(1,   1,NZ_s)
    M_s(NX_s+1,1,NZ_s+1)=    M_s(NX_s,1,NZ_s)

    
     DO i=0,NX_L+1,2; DO j=0,NZ_L+1,2
      i2 = i/2; j2 = j/2
      M_L(i,1,j) = M_s(i2,1,j2) 
     ENDDO; ENDDO
     DO i=0,NX_L+1,2; DO j=1,NZ_L,2
      M_L(i,1,j) = 0.5*(M_L(i,1,j-1)+M_L(i,1,j+1))
     ENDDO; ENDDO     
     DO i=1,NX_L,2; DO j=1,NZ_L
      M_L(i,1,j) = 0.5*(M_L(i-1,1,j)+M_L(i+1,1,j))
     ENDDO; ENDDO
    ELSE
     WRITE(*,*) 'CloneMatrix2Large error'; STOP       
    ENDIF
ENDIF     

M_L_OUT = 0
M_L_OUT(1:NX_L,1:NY_L,1:NZ_L) = M_L(1:NX_L,1:NY_L,1:NZ_L)

ENDSUBROUTINE CloneMatrix2Large