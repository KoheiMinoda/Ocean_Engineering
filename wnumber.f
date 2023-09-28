      SUBROUTINE WNUMBER(X,Wk0)
      IMPLICIT NONE
      DOUBLE PRECISION X , a , p
      INTEGER i
      DOUBLE PRECISION Wk0
      DOUBLE PRECISION c(9) , b(6) , sum
      DATA c/1.D0 , -0.33333372D0 , -0.01109668D0 , 0.01726435D0 ,
     &     0.01325580D0 , -0.00116594D0 , 0.00829006D0 , -0.01252603D0 ,
     &     0.00404923D0/
      DATA b/0.000000122D0 , 0.073250017D0 , -0.009899981D0 ,
     &     0.002640863D0 , -0.000829239D0 , -0.000176411D0/
      SAVE b , c
      IF ( X.LT.0.0D0 ) STOP 'error at wnumber'
      IF ( X.LE.2.0D0 ) THEN
         a = 0.5D0*X
         p = 1.0D0
         sum = c(1)
         DO i = 2 , 9
            p = p*a
            sum = sum + c(i)*p
         ENDDO
         Wk0 = SQRT(X)/sum
      ELSE
         a = 0.5D0*X*EXP(4.0D0-2.0D0*X)
         p = 1.0D0
         sum = b(1)
         DO i = 2 , 6
            p = p*a
            sum = sum + b(i)*p
         ENDDO
         Wk0 = X + sum
      ENDIF
      END
