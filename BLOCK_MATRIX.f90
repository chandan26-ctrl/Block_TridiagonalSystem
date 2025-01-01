!****************************************************************************                                                                         
! *   FILE         = BLOCK_MATRIX.F90                                        *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *                                                                        
!**************************************************************************** 
      PROGRAM BLOCK_MATRIX
      IMPLICIT NONE
      INTEGER, PARAMETER :: M=8
      DOUBLE PRECISION :: A(M,M), B(M), X(M)
      INTEGER :: I, N
      
       

      A(1,:) = (/1,10,0,0,0,0,0,0/)
      A(2,:) = (/5,2,8,0,0,0,0,0/)
      A(3,:) = (/0,13,4,6,0,0,0,0/)
      A(4,:) = (/0,0,1,1,7,0,0,0/)
      A(5,:) = (/0,0,0,4,3,5,0,0/)
      A(6,:) = (/0,0,0,0,1,6,3,0/)
      A(7,:) = (/0,0,0,0,0,4,9,1/)
      A(8,:) =(/0,0,0,0,0,0,2,5/)

      B(:) = (/3,2,5,7,8,9,4,1/)
      ! SIZE OF BLOCK

      N=2  

      CALL SOLVE(A,B,M,N,X)
     
      DO I =1, M
         PRINT*, X(I)
      END DO


      END PROGRAM

      SUBROUTINE SOLVE(A,B, M, N, X)
      IMPLICIT NONE
      INTEGER :: M, N
      DOUBLE PRECISION :: A(M,M), B(M), X(M)
      
      INTEGER :: ROWS, COLS, REM, NBLK
      INTEGER :: I,J,P, K
      INTEGER :: KM1, KM2, KI, KJ
      DOUBLE PRECISION, ALLOCATABLE :: B1(:,:), X1(:,:), D(:,:,:), C1(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: C(:,:,:), BL(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Q(:,:,:), G(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: MI(:,:), IM(:,:) , M2(:,:)
      
      ROWS=M      
      REM= MOD(ROWS, N)
 
      IF  (REM .NE. 0) THEN
          PRINT*, 'MATRIX CANNOT BE DIVIDED IN BLOCKS'
      END IF

      !NO. OF BLOCKS

      NBLK=ROWS/N

      ALLOCATE( B1(N, NBLK), C1 (N, NBLK))
      ALLOCATE( X1(N, NBLK))
      ALLOCATE( D(N, N, NBLK), BL(N,N,NBLK), C(N,N,NBLK))
      ALLOCATE( Q(N,N,NBLK), G(N,N,NBLK))
      ALLOCATE( MI(N,N), IM(N,N), M2(N,N))
      
      
      
      P=1
      DO J = 1 ,NBLK
        DO I =1, N
           B1(I,J)=B(P)
            P=P+1
        END DO
      END DO

      DO I = 1, N
         DO J = 1, NBLK
          X1(I,J) = 0.0D0
          C1(I,J) =  0.0D0
         END DO
      END DO

      DO I = 1,N
        DO J = 1,N
          DO K= 1, NBLK
            D(I,J,K) =  0.0D0
            Q(I,J,K) =  0.0D0
            G(I,J,K) =  0.0D0
          END DO
        END DO
      END DO

      DO I = 1,N
        DO J = 1,N
          DO K= 1, NBLK-1
            C(I,J,K) = 0.0D0
            BL(I,J,K) = 0.0D0
          END DO
        END DO
      END DO

 
      DO K =1, NBLK-1
         KM1 =(1+(K-1)*N)
         KM2 =(N+(K-1)*N)
         KI=0 
         DO I =KM1,KM2
           KI =KI+1
           KJ=0
           DO J=KM1, KM2
             KJ=KJ+1
             D(KI,KJ,K) =A(I,J)
           END DO
         END DO
         KI=0
       

         DO I = KM1+N, KM2+N  
            KI = KI +1
            KJ=0
            DO J =KM1, KM2
              KJ = KJ +1
              BL(KI,KJ,K) = A(I,J)
            END DO
         END DO

         KI=0
       
         DO I = KM1, KM2 
            KI = KI +1
            KJ=0
            DO J =KM1+N, KM2+N
              KJ = KJ +1
              C(KI,KJ,K) = A(I,J)
            END DO
         END DO

         KM1 = KM1+N
         KM2 = KM2+N
      END DO


      KM1 = 1+(NBLK-1)*N
      KM2 = N+(NBLK-1)*N

      KI=0
       
      DO I= KM1,KM2
         KI=KI+1
         KJ=0
         DO J = KM1,KM2
            KJ=KJ+1
            D(KI,KJ,NBLK)= A(I,J)
         END DO
      END DO

      DO I = 1,N
         DO J =1,N
           Q(I,J,1)= D(I,J,1)
         END DO
      END DO

      DO I = 1,N
         DO J =1,N
            MI(I,J) =Q(I,J,1)
         END DO
      END DO


       
      CALL INVERSE(MI, IM, N)
        
      DO I = 1,N
        DO J =1,N
          MI(I,J) =C(I,J,1)
        END DO
      END DO

      CALL MULTIAB (IM , MI , N, N, M2)

      DO I = 1,N
        DO J =1,N
          G(I,J,1)= M2(I,J)
        END DO
      END DO

      DO K = 2, NBLK-1
        
        DO I =1,N
          DO J = 1,N
             IM(I,J) = BL(I,J,K-1)
             MI(I,J) = G(I,J,K-1)
          END DO
        END DO

        CALL MULTIAB(IM, MI, N, N, M2)
         
        DO I =1,N
          DO J = 1,N
            Q(I,J,K) = D(I,J,K)-M2(I,J)
            MI(I,J)=Q(I,J,K)
          END DO
        END DO
        
        CALL INVERSE(MI, IM, N)

        DO I = 1,N
          DO J =1,N
            MI(I,J)=C(I,J,K)
          END DO
        END DO
          
        CALL MULTIAB (IM , MI , N, N, M2)
        
        
        DO I = 1,N
          DO J = 1,N
            G(I,J,K) = M2(I,J)
          END DO
        END DO

      END DO


      DO I =1,N
        DO J = 1,N
           IM(I,J) = BL(I,J,NBLK-1)
           MI(I,J) = G(I,J,NBLK-1)
        END DO
      END DO
         
      CALL MULTIAB(IM, MI, N, N, M2)

      DO I = 1,N
        DO J = 1,N
           Q(I,J,NBLK) = D(I,J,NBLK)-M2(I,J)
        END DO
      END DO
       
      DO I = 1,N
        DO J = 1,N
           MI(I,J)= Q(I,J,1)
        END DO
      END DO

      CALL INVERSE(MI, IM, N)
         
      DO I =1,N
         MI(I,1)=B1(I,1)
      END DO
         
      CALL MULTIAB(IM, MI, N, 1, M2)
          
      DO I =1, N
         C1(I,1)= M2(I,1)
      END DO


      DO K = 2, NBLK

        DO I = 1,N
          DO J = 1,N
             IM (I,J)= BL(I,J,K-1)
          END DO
          MI(I,1) = C1(I,K-1)
        END DO

        CALL MULTIAB(IM, MI, N, 1, M2)
            
        DO I = 1, N
          DO J=1,N
             IM(I,J) = Q(I,J,K)
          END DO
        END DO
              
        CALL INVERSE(IM, MI, N)
             
        DO I =1,N
           IM(I,1) = B1(I,K)-M2(I,1)
        END DO

        CALL MULTIAB(MI, IM, N, 1, M2)

        DO I = 1,N
           C1(I,K)=M2(I,1)
        END DO

      END DO

      DO I =1,N
         X1(I,NBLK)= C1(I,NBLK)
      END DO

      DO K=(NBLK-1), 1, -1

         DO I =1, N
           DO J =1,N
              MI(I,J) = G(I,J,K)
           END DO
           IM(I,1) = X1(I,K+1)
        END DO

            
        CALL MULTIAB(MI, IM, N, 1, M2)

        DO I =1, N
           X1(I,K) =C1(I,K)-M2(I,1)
        END DO

      END DO

      J=0
      DO K=1, NBLK
        DO I =1,N
           J = J+1
           X(J)= X1(I,K)
        END DO
      END DO



     END SUBROUTINE

     SUBROUTINE INVERSE(A,C,N)
     IMPLICIT NONE 
     INTEGER N
     DOUBLE PRECISION A(N,N), C(N,N), Y(N), Z(N)
     DOUBLE PRECISION L(N,N), U(N,N), B(N), D(N), X(N)
     DOUBLE PRECISION COEFF, SUMM
     INTEGER I, J, K

     L=0.0
     U=0.0
     B=0.0

! STEP 1: FORWARD ELIMINATION
     DO K=1, N-1
      DO I=K+1,N
        COEFF=A(I,K)/A(K,K)
        L(I,K) = COEFF
        DO J=K+1,N
          A(I,J) = A(I,J)-COEFF*A(K,J)
        END DO
      END DO
     END DO

! STEP 2: PREPARE L AND U MATRICES 

     DO I=1,N
        L(I,I) = 1.0
     END DO
! U MATRIX IS THE UPPER TRIANGULAR PART OF A
     DO J=1,N
       DO I=1,J
         U(I,J) = A(I,J)
      END DO
     END DO

     D(1)=0.0
     D(2)=0.0
     D(3)=1.0
     DO I=1,N
        Y(I)=0.0
     END DO
! STEP 3A: SOLVE LD=B USING THE FORWARD SUBSTITUTION
     DO I=1,N
        SUMM =0.0
      DO J=1,N
        SUMM=SUMM-L(I,J)*Y(J)
      END DO
       Y(I)=(D(I)+SUMM)/L(I,I)
     END DO
! STEP 3B: SOLVE UX=D USING THE BACK SUBSTITUTION
     DO I=1,N
        D(I)=0.0
     END DO
     DO I=N,1,-1
       SUMM=0.0
     DO J=N,1,-1
       SUMM= SUMM-U(I,J)*D(J)
     END DO
       D(I)=(Y(I)+SUMM)/U(I,I)
     END DO
    


! STEP 3: COMPUTE COLUMNS OF THE INVERSE MATRIX C
     DO K=1,N
        B(K)=1.0
        D(1) = B(1)
! STEP 3A: SOLVE LD=B USING THE FORWARD SUBSTITUTION
     DO I=2,N
         D(I)=B(I)
     DO J=1,I-1
         D(I) = D(I) - L(I,J)*D(J)
     END DO
     END DO
! STEP 3B: SOLVE UX=D USING THE BACK SUBSTITUTION
     X(N)=D(N)/U(N,N)
       DO I = N-1,1,-1
          X(I) = D(I)
          DO J=N,I+1,-1
            X(I)=X(I)-U(I,J)*X(J)
          END DO
           X(I) = X(I)/U(I,I)
     END DO
! STEP 3C: FILL THE SOLUTIONS X(N) INTO COLUMN K OF C
     DO I=1,N
        C(I,K) = X(I)
     END DO
       B(K)=0.0
     END DO
     END SUBROUTINE INVERSE


     SUBROUTINE MULTIAB(A,B,N,M,X)
     IMPLICIT NONE
     INTEGER :: I,J,K,N,M
     DOUBLE PRECISION :: A(N,N), B(N,M), X(N,M)
      
     DO I =1, N
      DO J =1,M
         X(I,J) = 0.0D0
        DO K =1,N
         X(I,J) =X(I,J) +A(I,K)*B(K,J)
        END DO
      END DO
     END DO

     END SUBROUTINE
      

