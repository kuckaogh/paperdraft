SUBROUTINE VelocityBC(uuB,vvB,wwB)   !INOUT  
!_______________________________________________________________________________________

USE COMDAT

REAL(8)    ,DIMENSION(0:NX  ,0:NY+1,0:NZ+1) ,INTENT(INOUT)   :: uuB
REAL(8)    ,DIMENSION(0:NX+1,0:NY  ,0:NZ+1) ,INTENT(INOUT)   :: vvB
REAL(8)    ,DIMENSION(0:NX+1,0:NY+1,0:NZ  ) ,INTENT(INOUT)   :: wwB

!________________________________________________________________

IF     (SideSlip==0) THEN !Nonslip

        !Ghost Points for du/dy ! slip or non-slip
        !---------------------------------- uuB --------
        DO i=1,NX-1     ! i=1, NX-1
        DO m=0,NY+1
        uuB(i,m,NZ+1)  = 0 - uuB(i,m,NZ) 
        uuB(i,m,0)     = 0 - uuB(i,m,1)
        ENDDO
        ENDDO
        DO i=1,NX-1 ! i=1, NX-1
        DO j=0,NZ+1
        uuB(i,NY+1,j)  = 0 - uuB(i,NY,j) 
        uuB(i,0,j)     = 0 - uuB(i,1,j)
        ENDDO
        ENDDO
        !---------------------------------- vvB --------
        DO m=1,NY-1
        DO i=1-1,NX+1 ! i=1, NX
        vvB(i,m,NZ+1)  = 0 - vvB(i,m,NZ) 
        vvB(i,m,0)     = 0 - vvB(i,m,1)
        ENDDO
        ENDDO
        DO m=1,NY-1
        DO j=0,NZ+1
        vvB(NX+1,m,j)  = 0 - vvB(NX,m,j) 
        vvB(0,m,j)     = 0 - vvB(1,m,j)
        ENDDO
        ENDDO

        !---------------------------------- wwB --------
        DO j=1,NZ-1
        DO m=0,NY+1
        wwB(NX+1,m,j) =0 - wwB(NX,m,j)
        wwB(0,m,j)    =0 - wwB(1,m,j)
        ENDDO
        ENDDO
        DO j=1,NZ-1
        DO i=0,NX+1
        wwB(i,NY+1,j) =0 - wwB(i,NY,j)
        wwB(i,0,j)    =0 - wwB(i,1,j)
        ENDDO
        ENDDO

ELSEIF (SideSlip==1) THEN !TopSlip                              

        !Ghost Points for du/dy ! slip or non-slip
        !---------------------------------- uuB --------
        DO i=1,NX-1     ! i=1, NX-1
        DO m=0,NY+1
        uuB(i,m,NZ+1)  =     uuB(i,m,NZ) 
        uuB(i,m,0)     = 0 - uuB(i,m,1)
        ENDDO
        ENDDO
        DO i=1,NX-1 ! i=1, NX-1
        DO j=0,NZ+1
        uuB(i,NY+1,j)  = 0 - uuB(i,NY,j) 
        uuB(i,0,j)     = 0 - uuB(i,1,j)
        ENDDO
        ENDDO
        !---------------------------------- vvB --------
        DO m=1,NY-1
        DO i=1-1,NX+1 ! i=1, NX
        vvB(i,m,NZ+1)  =     vvB(i,m,NZ) 
        vvB(i,m,0)     = 0 - vvB(i,m,1)
        ENDDO
        ENDDO
        DO m=1,NY-1
        DO j=0,NZ+1
        vvB(NX+1,m,j)  = 0 - vvB(NX,m,j) 
        vvB(0,m,j)     = 0 - vvB(1,m,j)
        ENDDO
        ENDDO
        !---------------------------------- wwB --------
        DO j=1,NZ-1
        DO m=0,NY+1
        wwB(NX+1,m,j) = 0 - wwB(NX,m,j)
        wwB(0,m,j)    = 0 - wwB(1,m,j)
        ENDDO
        ENDDO
        DO j=1,NZ-1
        DO i=0,NX+1
        wwB(i,NY+1,j) = 0 - wwB(i,NY,j)
        wwB(i,0,j)    = 0 - wwB(i,1,j)
        ENDDO
        ENDDO

ELSEIF (SideSlip==2) THEN !Top and Side Slip
                                       
        !Ghost Points for du/dy ! slip or non-slip
        !---------------------------------- uuB --------
        DO i=1,NX-1     ! i=1, NX-1
        DO m=0,NY+1
        uuB(i,m,NZ+1)  =     uuB(i,m,NZ) 
        uuB(i,m,0)     = 0 - uuB(i,m,1)
        ENDDO
        ENDDO
        DO i=1,NX-1 ! i=1, NX-1
        DO j=0,NZ+1
        uuB(i,NY+1,j)  =     uuB(i,NY,j) 
        uuB(i,0,j)     =     uuB(i,1,j)
        ENDDO
        ENDDO
        !---------------------------------- vvB --------
        DO m=1,NY-1
        DO i=1-1,NX+1 ! i=1, NX
        vvB(i,m,NZ+1)  =     vvB(i,m,NZ) 
        vvB(i,m,0)     = 0 - vvB(i,m,1)
        ENDDO
        ENDDO
        DO m=1,NY-1
        DO j=0,NZ+1
        vvB(NX+1,m,j)  = 0 - vvB(NX,m,j) 
        vvB(0,m,j)     = 0 - vvB(1,m,j)
        ENDDO
        ENDDO
        !---------------------------------- wwB --------
        DO j=1,NZ-1
        DO m=0,NY+1
        wwB(NX+1,m,j) = 0 - wwB(NX,m,j)
        wwB(0,m,j)    = 0 - wwB(1,m,j)
        ENDDO
        ENDDO
        DO j=1,NZ-1
        DO i=0,NX+1
        wwB(i,NY+1,j) =     wwB(i,NY,j)
        wwB(i,0,j)    =     wwB(i,1,j)
        ENDDO
        ENDDO

ELSEIF (SideSlip==3) THEN !All Slip

        !Ghost Points for du/dy ! slip or non-slip
        !---------------------------------- uuB --------
        DO i=1,NX-1     ! i=1, NX-1
        DO m=0,NY+1
        uuB(i,m,NZ+1)  =  uuB(i,m,NZ) 
        uuB(i,m,0)     =  uuB(i,m,1)
        ENDDO
        ENDDO
        DO i=1,NX-1 ! i=1, NX-1
        DO j=0,NZ+1
        uuB(i,NY+1,j) =  uuB(i,NY,j) 
        uuB(i,0,j)    =  uuB(i,1,j)
        ENDDO
        ENDDO
        !---------------------------------- vvB --------
        DO m=1,NY-1
        DO i=0,NX+1 
        vvB(i,m,NZ+1)  =  vvB(i,m,NZ) 
        vvB(i,m,0)     =  vvB(i,m,1)
        ENDDO
        ENDDO
        DO m=1,NY-1
        DO j=0,NZ+1
        vvB(NX+1,m,j)  =  vvB(NX,m,j) 
        vvB(0,m,j)     =  vvB(1,m,j)
        ENDDO
        ENDDO
        !---------------------------------- wwB --------
        DO j=1,NZ-1
        DO m=0,NY+1
        wwB(NX+1,m,j) =  wwB(NX,m,j)
        wwB(0,m,j)    =  wwB(1,m,j)
        ENDDO
        ENDDO
        DO j=1,NZ-1
        DO i=0,NX+1
        wwB(i,NY+1,j) =  wwB(i,NY,j)
        wwB(i,0,j)    =  wwB(i,1,j)
        ENDDO
        ENDDO

ELSE

        WRITE(*,*) 'Slip Model Error.....';STOP 

ENDIF



ENDSUBROUTINE VelocityBC 