SUBROUTINE AveMatrix2Small(M_s,NX_s,NY_s,NZ_s,M_I,NX_I,NY_I,NZ_I) 

INTEGER(4)                                ,INTENT(IN) :: NX_I,NY_I,NZ_I
INTEGER(4)                                ,INTENT(IN) :: NX_s,NY_s,NZ_s
REAL(8)   ,DIMENSION(1:NX_I,1:NY_I,1:NZ_I),INTENT(IN) :: M_I
REAL(8)   ,DIMENSION(1:NX_s,1:NY_s,1:NZ_s),INTENT(OUT):: M_s
INTEGER(4)                                            :: i,m,j

M_s=0

IF (NX_s*2==NX_I) THEN

    IF (NY_I==1) THEN
      DO i=2,NX_I,2; DO j=2,NZ_I,2
       M_s(i/2,1,j/2) = 0.25*( M_I(i-1,1,j-1)+M_I(i-1,1,j)+M_I(i,1,j-1)+M_I(i,1,j) ) 
      ENDDO; ENDDO
    ELSE
      DO i=2,NX_I,2; DO j=2,NZ_I,2; DO m=2,NY_I,2
       M_s(i/2,m/2,j/2) = &
                 0.125*( M_I(i-1,m  ,j-1)+M_I(i-1,m  ,j)+M_I(i,m  ,j-1)+M_I(i,m  ,j)  &
                        +M_I(i-1,m-1,j-1)+M_I(i-1,m-1,j)+M_I(i,m-1,j-1)+M_I(i,m-1,j) )                                                                   
      ENDDO; ENDDO; ENDDO
    ENDIF

ELSE

    IF (NY_I==1) THEN
      DO i=1,NX_s ; DO j=1,NZ_s
       M_s(i,1,j) = 0.25*( 0.25*M_I(2*i-1,1,2*j-1)+0.5*M_I(2*i  ,1,2*j-1)+0.25*M_I(2*i+1,1,2*j-1) &
                         + 0.5 *M_I(2*i-1,1,2*j  )+    M_I(2*i  ,1,2*j  )+0.5* M_I(2*i+1,1,2*j  ) &
                         + 0.25*M_I(2*i-1,1,2*j+1)+0.5*M_I(2*i  ,1,2*j+1)+0.25*M_I(2*i+1,1,2*j+1) ) 
       
        
      ENDDO; ENDDO
    ELSE
      DO i=1,NX_s ; DO m=1,NY_s ; DO j=1,NZ_s
       M_s(i,1,j) = 0.125*(                                                                    &
       0.125*M_I(2*i-1,2*m-1,2*j-1)+0.25*M_I(2*i  ,2*m-1,2*j-1)+0.125*M_I(2*i+1,2*m-1,2*j-1)   &
     + 0.25* M_I(2*i-1,2*m-1,2*j  )+0.5* M_I(2*i  ,2*m-1,2*j  )+0.25* M_I(2*i+1,2*m-1,2*j  )   &
     + 0.125*M_I(2*i-1,2*m-1,2*j+1)+0.25*M_I(2*i  ,2*m-1,2*j+1)+0.125*M_I(2*i+1,2*m-1,2*j+1)   & 
     + 0.25* M_I(2*i-1,2*m  ,2*j-1)+0.5* M_I(2*i  ,2*m  ,2*j-1)+0.25* M_I(2*i+1,2*m  ,2*j-1)   &
     + 0.5*  M_I(2*i-1,2*m  ,2*j  )+     M_I(2*i  ,2*m  ,2*j  )+0.5*  M_I(2*i+1,2*m  ,2*j  )   &
     + 0.25* M_I(2*i-1,2*m  ,2*j+1)+0.5* M_I(2*i  ,2*m  ,2*j+1)+0.25* M_I(2*i+1,2*m  ,2*j+1)   &
     + 0.125*M_I(2*i-1,2*m+1,2*j-1)+0.25*M_I(2*i  ,2*m+1,2*j-1)+0.125*M_I(2*i+1,2*m+1,2*j-1)   &
     + 0.25* M_I(2*i-1,2*m+1,2*j  )+0.5* M_I(2*i  ,2*m+1,2*j  )+0.25* M_I(2*i+1,2*m+1,2*j  )   &
     + 0.125*M_I(2*i-1,2*m+1,2*j+1)+0.25*M_I(2*i  ,2*m+1,2*j+1)+0.125*M_I(2*i+1,2*m+1,2*j+1)   )                              
      ENDDO; ENDDO; ENDDO
    ENDIF

ENDIF

ENDSUBROUTINE AveMatrix2Small