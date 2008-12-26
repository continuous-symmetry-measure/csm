      PROGRAM MIR_SYMMETRY

      PARAMETER (NN      = 500)
      PARAMETER (N_POINT = 1000000)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*100 CHARBUFF
      INTEGER*8 I0
      REAL*4 FZ_XY(NN,NN,NN)
      DIMENSION FX(NN), FY_X(NN,NN)
      DIMENSION XX(NN), YY(NN), ZZ(NN)
      DIMENSION PX(N_POINT),PY(N_POINT),PZ(N_POINT)
      DIMENSION NUMBERS(N_POINT)
      DIMENSION NP(21)
      DIMENSION DIR_COS(3)

      COMMON /RRAND/ I0

      DATA I0/762399453125/
      DATA NP/10,18,32,56,100,178,316,562,1000,
     *1778,3162,5622,10000,17782,31622,56234,100000,
     *177820,316220,562340,1000000/


c----------------------------------------------------------


      READ *, DIR_COS
      
      ABSV = DIR_COS(1)**2+DIR_COS(2)**2+DIR_COS(3)**2
      IF(ABS(ABSV-1.0D0).GT.1.D-10) THEN
         PRINT *,'WRONG DIR. COS. VECTOR'
         STOP
      END IF
      

      CALL GETARG(3, CHARBUFF)
      OPEN(UNIT=10,FILE=CHARBUFF,
     *     STATUS='UNKNOWN',FORM='FORMATTED')


      CALL PDF(NN,XX,YY,ZZ,FX,FY_X,FZ_XY)

      DO I=1,21
         N = NP(I)

         I0 = 762399453125

         CALL GEN(N,NN,XX,YY,ZZ,FX,FY_X,FZ_XY,PX,PY,PZ)

         CALL CENTER(N,PX,PY,PZ)

         CALL MIR_SYMM(N, PX, PY, PZ, NUMBERS, SYM, DIR_COS)

         WRITE(10,*) N, SYM
         ENDFILE(10)

      END DO

      CLOSE(10)

      STOP
      END


C*********************************************************************
C*********************************************************************



      SUBROUTINE TDRAND(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*8 I, IA, IB, IT
      COMMON /RRAND/ I
      DATA IA, IB, IT /513, 29741096258473, 140737488355328/
      DATA T /1.40737488355328D14/
      I=IA*I+IB
      I=I-(I/IT)*IT
      X=DBLE(I)/T

      RETURN
      END


C*********************************************************************
C*********************************************************************


      SUBROUTINE PDF(NN,X,Y,Z,FX,FY_X,FZ_XY)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*100 CHARBUFF
      INTEGER*4 nr1, nr2, nx, ny, nz
      REAL*4 xmax, xmin, ymax, ymin, zmax, zmin
      REAL*4 FZ_XY(NN,NN,NN)

      DIMENSION X(NN),Y(NN),Z(NN)
      DIMENSION W_XY(NN,NN),W_X(NN)
      DIMENSION FX(NN), FY_X(NN,NN)


c*******   INPUT    ***************


        CALL GETARG(1,CHARBUFF)
      OPEN(UNIT=3,FILE=CHARBUFF,
     *STATUS='OLD',ACCESS='DIRECT',RECL=44, FORM='UNFORMATTED')
      nrec = 1
      read (3,rec=nrec) nr1, nr2, nx, ny, nz,
     *      zmin, zmax, ymin, ymax, xmin, xmax
     
     
      dx = (xmax-xmin)/dble(nx-1)
      dy = (ymax-ymin)/dble(ny-1)
      dz = (zmax-zmin)/dble(nz-1)
      CLOSE(3)     
      

      CALL GETARG(2, CHARBUFF)
      OPEN(UNIT=3,FILE=CHARBUFF,
     *STATUS='OLD',ACCESS='DIRECT',RECL=nx * 4, FORM='UNFORMATTED')
      NREC = 1
      DO K=1,nz
         DO J=1,ny
            read (3,rec=nrec)(FZ_XY(I,J,K),i=1,nx)
            NREC = NREC+1
         END DO
      END DO
      CLOSE(3)

      DO I=1,NX
         X(I) = XMIN+DX*(I-1)
      END DO
      DO J=1,NY
         Y(J) = YMIN+DY*(J-1)
      END DO
      DO K=1,NZ
         Z(K) = ZMIN+DZ*(K-1)
      END DO


c************   normal 3d test   *********


c      dx = 0.01d0
c      dy = 0.01d0
c      dz = 0.01d0
c
c      DO I=1,NN
c         x(i) = -2.5d0+dx*(i-1)
c         y(i) = -2.5d0+dy*(i-1)
c         z(i) = -2.5d0+dz*(i-1)
c      end do
c
c      DO I=1,NN
c         DO J=1,NN
c            DO K=1,NN
c               FZ_XY(I,J,K)=exp(-(x(i)**2+y(j)**2+z(k)**2)-z(k)*Y(J))
c            END DO
c         END DO
c      END DO


c**********************
C     2D PDF W(X,Y) 
c**********************  


      DX = X(2)-X(1)
      DY = Y(2)-y(1)
      DZ = Z(2)-Z(1)

      DO I=1,NN
         DO J=1,NN
            W_XY(I,J) = 0.0D0
         END DO
      END DO

      DO I=1,NN
         DO J=1,NN
            DO K=1,NN
               W_XY(I,J) = W_XY(I,J) + FZ_XY(I,J,K)*DZ
            END DO
         END DO
      END DO


c**********************
C     1D PDF W(X) 
c********************** 


      DO I=1,NN
         W_X(I) = 0.0D0
      END DO

      DO I=1,NN
         DO J=1,NN
            W_X(I) = W_X(I) + W_XY(I,J)*DY
         END DO
      END DO


c*********************************
C     Conditional PDF W(Z|X,Y) 
c*********************************


      DO K=1,NN
         DO J=1,NN
            DO I=1,NN
               FZ_XY(I,J,K) = FZ_XY(I,J,K)/W_XY(I,J)
            END DO
         END DO
      END DO


c*********************************
C     Conditional PDF W(Y|X) 
c********************************* 


      DO J=1,NN
         DO I=1,NN
            W_XY(I,J) = W_XY(I,J)/W_X(I)
         END DO
      END DO


c****************************************
C     Conditional integral PDF  F(X)
c**************************************** 


      FX(1) = 0.0D0
      DO I=2,NN
         FX(I) = FX(I-1)+0.5D0*(W_X(I-1)+W_X(I))*DX
      END DO

      DO I=1,NN
         FX(I) = FX(I)/FX(NN)
      END DO


c****************************************
C     Conditional integral PDF F(Y|X)
c**************************************** 


      DO I=1,NN
         FY_X(I,1) = 0.0D0
         DO J=2,NN
            FY_X(I,J) = FY_X(I,J-1)+
     *      0.5D0*(W_XY(I,J-1)+W_XY(I,J))*DY
         END DO
      END DO

      DO I=1,NN
         DO J=1,NN
            FY_X(I,J) = FY_X(I,J)/FY_X(I,NN)
         END DO
      END DO



c*********************************************
C     Conditional integral PDF F(Z|X,Y)
c********************************************* 


      DO I=1,NN
         DO J=1,NN
            W1 = FZ_XY(I,J,1)
            FZ_XY(I,J,1) = 0.0D0
            DO K=2,NN
               W2 = FZ_XY(I,J,K)
               FZ_XY(I,J,K) = FZ_XY(I,J,K-1)+0.5D0*(W1+W2)*DZ
               W1 = W2
            END DO
         END DO
      END DO

      DO I=1,NN
         DO J=1,NN
            DO K=1,NN
               FZ_XY(I,J,K) = FZ_XY(I,J,K)/FZ_XY(I,J,NN)
            END DO
         END DO
      END DO


      RETURN
      END


C*********************************************************************
C*********************************************************************


      SUBROUTINE GEN(N_POINT,NN,XX,YY,ZZ,FX,FY_X,FZ_XY,PPX,PPY,PPZ)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      REAL*4 FZ_XY(NN,NN,NN)

      DIMENSION XX(NN), YY(NN), ZZ(NN), FX(NN), FY(NN), FZ(NN)
      DIMENSION FY_X(NN,NN)
      DIMENSION PPX(N_POINT),PPY(N_POINT),PPZ(N_POINT) 


c*******   First  random number    ***************


      DO KK=1,100
         CALL TDRAND(AX)
      END DO

      DO KK=1,N_POINT
         CALL TDRAND(AX)
         DO I=1,NN
            IF(AX.LT.FX(I)) THEN
               M=I
               GO TO 1
            END IF
         END DO
    1    CONTINUE
         PX = XX(M-1)+(XX(M)-XX(M-1))/(FX(M)-FX(M-1))*(AX-FX(M-1))

         I1 = (PX-XX(1))/(XX(2)-XX(1))
         I1 = I1+1
         I2 = I1+1

         DO J=1,NN
            FY(J) = (FY_X(I1,J)*(XX(I2)-PX)+FY_X(I2,J)*(PX-XX(I1)))/
     *              (XX(I2)-XX(I1))
         END DO


c*******   Second  random number    ***************


         CALL TDRAND(AY)
         CALL TDRAND(AY)
         CALL TDRAND(AY)
         CALL TDRAND(AY)

         DO J=1,NN
            IF(AY.LT.FY(J)) THEN
               M=J
               GO TO 2
            END IF
         END DO
    2    CONTINUE
         PY = YY(M-1)+(YY(M)-YY(M-1))/(FY(M)-FY(M-1))*(AY-FY(M-1))

         J1 = (PY-YY(1))/(YY(2)-YY(1))
         J1 = J1+1
         J2 = J1+1

         DO K=1,NN
            FZ1 = (FZ_XY(I1,J1,K)*(XX(I2)-PX)+FZ_XY(I2,J1,K)*
     *            (PX-XX(I1)))/(XX(I2)-XX(I1))
            FZ2 = (FZ_XY(I1,J2,K)*(XX(I2)-PX)+FZ_XY(I2,J2,K)*
     *            (PX-XX(I1)))/(XX(I2)-XX(I1))
            FZ(K) = FZ1+(FZ2-FZ1)*(PY-YY(J1))/(YY(J2)-YY(J1))
         END DO


c*******  Third  random number    ***************


         CALL TDRAND(AZ)
         CALL TDRAND(AZ)
         CALL TDRAND(AZ)
         CALL TDRAND(AZ)

         DO K=1,NN
            IF(AZ.LT.FZ(K)) THEN
               M=K
               GO TO 3
            END IF
         END DO
    3    CONTINUE
         PZ = ZZ(M-1)+(ZZ(M)-ZZ(M-1))/(FZ(M)-FZ(M-1))*(AZ-FZ(M-1))

         PPX(KK) = PX
         PPY(KK) = PY
         PPZ(KK) = PZ

      END DO

      RETURN
      END


C*********************************************************************
C*********************************************************************


      SUBROUTINE CENTER(N_POINT,PX,PY,PZ)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION PX(N_POINT),PY(N_POINT),PZ(N_POINT)

      XC = 0.0D0
      YC = 0.0D0
      ZC = 0.0D0
      DO I=1,N_POINT
         XC = XC+PX(I)
         YC = YC+PY(I)
         ZC = ZC+PZ(I)
      END DO
      XC = XC/DBLE(N_POINT)
      YC = YC/DBLE(N_POINT)
      ZC = ZC/DBLE(N_POINT)

      SIG = 0.0D0
      DO I=1,N_POINT
         SIG = SIG+(PX(I)-XC)**2+(PY(I)-YC)**2+(PZ(I)-ZC)**2
      END DO
      SIG  = SQRT(SIG/DBLE(N_POINT))

      DO I=1,N_POINT
         PX(I) = (PX(I)-XC)/SIG
         PY(I) = (PY(I)-YC)/SIG
         PZ(I) = (PZ(I)-ZC)/SIG
      END DO

      RETURN
      END 


C*********************************************************************
C*********************************************************************


      SUBROUTINE MIR_SYMM(N_POINT, PX, PY, PZ, NUMBERS, SYM, DIR_COS)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION PX(N_POINT),PY(N_POINT),PZ(N_POINT)
      DIMENSION NUMBERS(N_POINT)
      DIMENSION CX(2), CY(2), CZ(2)
      DIMENSION DIR_COS(3)


      DO I=1,N_POINT
         NUMBERS(I) = I
      END DO

      NNP = N_POINT
      SYM = 0.0D0

      sum = 0.0d0
      
      

      DO I=1,N_POINT,2

         CALL TDRAND(SSS)

         NNN   = SSS*NNP+1.0D0
         NN    = NUMBERS(NNN)
         CX(1) = PX(NN)
         CY(1) = PY(NN)
         CZ(1) = PZ(NN)

         NNP=NNP-1
         DO J=NNN,NNP
            NUMBERS(J)=NUMBERS(J+1)
         END DO

         R2MIN = 1.0D10
         DO II=1,NNP
            JJ = NUMBERS(II)
	  
	  
	    SK_MUL = PX(JJ)*DIR_COS(1)+PY(JJ)*DIR_COS(2)+
     *                PZ(JJ)*DIR_COS(3)
            X = PX(JJ)-2.0D0*SK_MUL*DIR_COS(1)
            Y = PY(JJ)-2.0D0*SK_MUL*DIR_COS(2)	  
            Z = PZ(JJ)-2.0D0*SK_MUL*DIR_COS(3)	  	    
            R2 = (X-CX(1))**2+(Y-CY(1))**2+(Z-CZ(1))**2
            IF(R2.LT.R2MIN) THEN
               R2MIN = R2
               KK    = JJ
               NNN   = II
            END IF
         END DO
         CX(2) = PX(KK)
         CY(2) = PY(KK)
         CZ(2) = PZ(KK)
         NNP = NNP-1
         DO II=NNN,NNP
            NUMBERS(II)=NUMBERS(II+1)
         END DO
         CALL CALC_MIR_SYMM(CX, CY, CZ, SIG2, DIR_COS)
         SYM = SYM+2.0d0*SIG2

      END DO
      

      SYM = 1.0D2*SYM/DBLE(N_POINT)
      

      RETURN
      END


C*********************************************************************
C*********************************************************************


      SUBROUTINE CALC_MIR_SYMM(CX, CY, CZ, SIG2, DIR_COS)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION CX(2), CY(2), CZ(2)
      DIMENSION DIR_COS(3)
      
      
      SK_MUL = CX(2)*DIR_COS(1)+CY(2)*DIR_COS(2)+
     *         CZ(2)*DIR_COS(3)
      X = CX(2)-2.0D0*SK_MUL*DIR_COS(1)
      Y = CY(2)-2.0D0*SK_MUL*DIR_COS(2)	  
      Z = CZ(2)-2.0D0*SK_MUL*DIR_COS(3)	  	          

      XC = 0.5D0*(CX(1)+X)
      YC = 0.5D0*(CY(1)+Y)
      ZC = 0.5D0*(CZ(1)+Z)

      SIG2 = (CX(1)-XC)**2+(CY(1)-YC)**2+(CZ(1)-ZC)**2

      RETURN
      END
