

*********************************                                       
*   Subroutine for integration                                          
********************************                                        
      SUBROUTINE GADAP(A0,B0,F,EPS,SUM)                                 
      COMMON/GADAP1/ NUM,IFU                                            
      EXTERNAL F                                                        
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)     
    1 FORMAT(16H GADAP:I TOO BIG)                                       
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)          
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.3   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)    
      F2(1)=F(0.5*(A0+B0))  
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)    
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(A(I)+B(I)-W1)   
      F2(I+1)=F3(I) 
      F3(I+1)=F(B(I)-A(I+2)+W1) 
      F1(I+2)=F(U2) 
      F2(I+2)=F2(I) 
      F3(I+2)=F(B(I+2)+A(I+2)-U2)   
      F1(I+3)=F(A(I)+A(I+2)-W1) 
      F2(I+3)=F1(I) 
      F3(I+3)=F(W1) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   


      SUBROUTINE GADAPu(A0,B0,F,EPS,SUM)                                 
      COMMON/GADAPu1/ NUM,IFU                                            
      EXTERNAL F                                            
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)     
*1     FORMAT(16H GADAPu:I TOO BIG)              
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)          
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.3   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)    
      F2(1)=F(0.5*(A0+B0))  
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)    
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(A(I)+B(I)-W1)   
      F2(I+1)=F3(I) 
      F3(I+1)=F(B(I)-A(I+2)+W1) 
      F1(I+2)=F(U2) 
      F2(I+2)=F2(I) 
      F3(I+2)=F(B(I+2)+A(I+2)-U2)   
      F1(I+3)=F(A(I)+A(I+2)-W1) 
      F2(I+3)=F1(I) 
      F3(I+3)=F(W1) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   

      SUBROUTINE GADAP2(A0,B0,FL,FU,F,EPS,SUM)  
      COMMON/GADAP_2/ NUM,IFU    
      EXTERNAL F,FL,FU  
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      X=0.5*(1+C)*A0+0.5*(1-C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F1(1)=FGADAP(X,AY,BY,F,EPS)   
      X=0.5*(A0+B0) 
      AY=FL(X)  
      BY=FU(X)  
      F2(1)=FGADAP(X,AY,BY,F,EPS)   
      X=0.5*(1-C)*A0+0.5*(1+C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F3(1)=FGADAP(X,AY,BY,F,EPS)   
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      X=A(I)+B(I)-W1    
      AY=FL(X)  
      BY=FU(X)  
      F1(I+1)=FGADAP(X,AY,BY,F,EPS) 
      F2(I+1)=F3(I) 
      X=B(I)-A(I+2)+W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+1)=FGADAP(X,AY,BY,F,EPS) 
      X=U2  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+2)=FGADAP(X,AY,BY,F,EPS) 
      F2(I+2)=F2(I) 
      X=B(I+2)+A(I+2)-U2    
      AY=FL(X)  
      BY=FU(X)  
      F3(I+2)=FGADAP(X,AY,BY,F,EPS) 
      X=A(I)+A(I+2)-W1  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+3)=FGADAP(X,AY,BY,F,EPS) 
      F2(I+3)=F1(I) 
      X=W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+3)=FGADAP(X,AY,BY,F,EPS) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   
      FUNCTION FGADAP(X,A0,B0,F,EPS)    
      COMMON/GADAP_2/ NUM,IFU    
      EXTERNAL F    
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(X,0.5*(1+C)*A0+0.5*(1-C)*B0)  
      F2(1)=F(X,0.5*(A0+B0))    
      F3(1)=F(X,0.5*(1-C)*A0+0.5*(1+C)*B0)  
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(X,A(I)+B(I)-W1) 
      F2(I+1)=F3(I) 
      F3(I+1)=F(X,B(I)-A(I+2)+W1)   
      F1(I+2)=F(X,U2)   
      F2(I+2)=F2(I) 
      F3(I+2)=F(X,B(I+2)+A(I+2)-U2) 
      F1(I+3)=F(X,A(I)+A(I+2)-W1)   
      F2(I+3)=F1(I) 
      F3(I+3)=F(X,W1)   
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  FGADAP=SUM    
      EPS=EPS/RED   
      RETURN    
      END   




      SUBROUTINE GADAPS2(A0,B0,FL,FU,F,EPS,SUM)  
      COMMON/GADAPS_2/ NUM,IFU    
      EXTERNAL F,FL,FU  
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      X=0.5*(1+C)*A0+0.5*(1-C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F1(1)=FGADAPS(X,AY,BY,F,EPS)   
      X=0.5*(A0+B0) 
      AY=FL(X)  
      BY=FU(X)  
      F2(1)=FGADAPS(X,AY,BY,F,EPS)   
      X=0.5*(1-C)*A0+0.5*(1+C)*B0   
      AY=FL(X)  
      BY=FU(X)  
      F3(1)=FGADAPS(X,AY,BY,F,EPS)   
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      X=A(I)+B(I)-W1    
      AY=FL(X)  
      BY=FU(X)  
      F1(I+1)=FGADAPS(X,AY,BY,F,EPS) 
      F2(I+1)=F3(I) 
      X=B(I)-A(I+2)+W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+1)=FGADAPS(X,AY,BY,F,EPS) 
      X=U2  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+2)=FGADAPS(X,AY,BY,F,EPS) 
      F2(I+2)=F2(I) 
      X=B(I+2)+A(I+2)-U2    
      AY=FL(X)  
      BY=FU(X)  
      F3(I+2)=FGADAPS(X,AY,BY,F,EPS) 
      X=A(I)+A(I+2)-W1  
      AY=FL(X)  
      BY=FU(X)  
      F1(I+3)=FGADAPS(X,AY,BY,F,EPS) 
      F2(I+3)=F1(I) 
      X=W1  
      AY=FL(X)  
      BY=FU(X)  
      F3(I+3)=FGADAPS(X,AY,BY,F,EPS) 
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  RETURN    
      END   
      FUNCTION FGADAPS(X,A0,B0,F,EPS)    
      COMMON/GADAPS_2/ NUM,IFU    
      EXTERNAL F    
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300) 
    1 FORMAT(16H GADAP:I TOO BIG)   
      DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)  
      IF(EPS.LT.1.0E-8) EPS=1.0E-8  
      RED=1.4   
      L=1   
      I=1   
      SUM=0.    
      C=SQRT(15.)/5.    
      A(1)=A0   
      B(1)=B0   
      F1(1)=F(X,0.5*(1+C)*A0+0.5*(1-C)*B0)  
      F2(1)=F(X,0.5*(A0+B0))    
      F3(1)=F(X,0.5*(1-C)*A0+0.5*(1+C)*B0)  
      IFU=3 
      S(1)=  DSUM(F1(1),F2(1),F3(1),A0,B0)  
  100 CONTINUE  
      L=L+1 
      N(L)=3    
      EPS=EPS*RED   
      A(I+1)=A(I)+C*(B(I)-A(I)) 
      B(I+1)=B(I)   
      A(I+2)=A(I)+B(I)-A(I+1)   
      B(I+2)=A(I+1) 
      A(I+3)=A(I)   
      B(I+3)=A(I+2) 
      W1=A(I)+(B(I)-A(I))/5.    
      U2=2.*W1-(A(I)+A(I+2))/2. 
      F1(I+1)=F(X,A(I)+B(I)-W1) 
      F2(I+1)=F3(I) 
      F3(I+1)=F(X,B(I)-A(I+2)+W1)   
      F1(I+2)=F(X,U2)   
      F2(I+2)=F2(I) 
      F3(I+2)=F(X,B(I+2)+A(I+2)-U2) 
      F1(I+3)=F(X,A(I)+A(I+2)-W1)   
      F2(I+3)=F1(I) 
      F3(I+3)=F(X,W1)   
      IFU=IFU+6 
      IF(IFU.GT.5000) GOTO 130  
      S(I+1)=  DSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))  
      S(I+2)=  DSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))  
      S(I+3)=  DSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))  
      SS=S(I+1)+S(I+2)+S(I+3)   
      I=I+3 
      IF(I.GT.300)GOTO 120  
      SOLD=S(I-3)   
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100  
      SUM=SUM+SS    
      I=I-4 
      N(L)=0    
      L=L-1 
  110 CONTINUE  
      IF(L.EQ.1) GOTO 130   
      N(L)=N(L)-1   
      EPS=EPS/RED   
      IF(N(L).NE.0) GOTO 100    
      I=I-1 
      L=L-1 
      GOTO 110  
  120 CONTINUE
C      WRITE(6,1)    
 130  FGADAPS=SUM    
      EPS=EPS/RED   
      RETURN    
      END   





      SUBROUTINE simpson(A0,B0,F,n,SUM) 
C ********************************************************************** 
C PURPOSE - INTEGRATE A FUNCTION F(X) 
C METHOD - STUPID
C USAGE - CALL simpson(a0,b0,f,n,sum)
C PARAMETERS A0 - LOWER LIMIT (INPUT,REAL) 
C B0 - UPPER LIMIT (INPUT,REAL) 
C F - FUNCTION F(X) TO BE INTEGRATED. MUST BE 
C SUPPLIED BY THE USER. (INPUT,REAL FUNCTION) 
C n - NUMBER OF DIVISIONS 
C SUM - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL) 
C PRECISION - SINGLE 
C REQ'D PROG'S - F 
C
* FIU 
C....................................................................... 
*      implicit real*8(a-h,o-x)
*      integer*8 n
      external f
*      tn = n
      width = (b0-a0)/n
      sumi = 0.0
      do i = 1,n
      x = a0 + float(i)*width 
      sumi = sumi + f(x)* width
      enddo
     
      sum = sumi
      return
      end

***************************************************************************
*
* $Id
*
* $Log
*
* #include "gen/pilot.h"
*#if defined(CERNLIB_DOUBLE)
*      SUBROUTINE RADMUL
*     1 (F,N,A,B,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,NFNEVL,IFAIL)
*      CHARACTER NAME*(*)
*      PARAMETER (NAME = 'RADMUL')
*      CALL MTLPRT(NAME,'D120',
*     +'not available on this machine - see documentation')
*      RETURN
*      END

*      SUBROUTINE DADMUL
*     1 (F,N,A,B,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,NFNEVL,IFAIL)
*#include "gen/imp64.inc"

*#else
*      SUBROUTINE DADMUL
*     1 (F,N,A,B,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,NFNEVL,IFAIL)
* #include "gen/imp128.inc"
*     CHARACTER NAME*(*)
*      PARAMETER (NAME = 'DADMUL')
*      CALL MTLPRT(NAME,'D120',
*     +'not available on this machine - see documentation')
*      RETURN
*      END

      SUBROUTINE RADMUL
     1 (F,N,A,B,MINPTS,MAXPTS,EPS,WK,IWK,RESULT,RELERR,NFNEVL,IFAIL)
*#endif
 
      LOGICAL LDV
 
      DIMENSION A(*),B(*),WK(*)
      DIMENSION CTR(15),WTH(15),WTHL(15),Z(15)
      DIMENSION W(2:15,5),WP(2:15,3)
 
      PARAMETER (R1 = 1, HF = R1/2)
 
      PARAMETER (XL2 =  0.35856 85828 00318 073D0)
      PARAMETER (XL4 =  0.94868 32980 50513 796D0)
      PARAMETER (XL5 =  0.68824 72016 11685 289D0)
 
      PARAMETER (W2 =  980*R1/6561, W4 = 200*R1/19683)
      PARAMETER (WP2 =  245*R1/486, WP4 = 25*R1/729)
 
      DATA (W(N,1),W(N,3),N=2,15)
     1/-0.193872885230909911D+00,  0.518213686937966768D-01,
     2 -0.555606360818980835D+00,  0.314992633236803330D-01,
     3 -0.876695625666819078D+00,  0.111771579535639891D-01,
     4 -0.115714067977442459D+01, -0.914494741655235473D-02,
     5 -0.139694152314179743D+01, -0.294670527866686986D-01,
     6 -0.159609815576893754D+01, -0.497891581567850424D-01,
     7 -0.175461057765584494D+01, -0.701112635269013768D-01,
     8 -0.187247878880251983D+01, -0.904333688970177241D-01,
     9 -0.194970278920896201D+01, -0.110755474267134071D+00,
     A -0.198628257887517146D+01, -0.131077579637250419D+00,
     B -0.198221815780114818D+01, -0.151399685007366752D+00,
     C -0.193750952598689219D+01, -0.171721790377483099D+00,
     D -0.185215668343240347D+01, -0.192043895747599447D+00,
     E -0.172615963013768225D+01, -0.212366001117715794D+00/
 
      DATA (W(N,5),W(N+1,5),N=2,14,2)
     1/ 0.871183254585174982D-01,  0.435591627292587508D-01,
     2  0.217795813646293754D-01,  0.108897906823146873D-01,
     3  0.544489534115734364D-02,  0.272244767057867193D-02,
     4  0.136122383528933596D-02,  0.680611917644667955D-03,
     5  0.340305958822333977D-03,  0.170152979411166995D-03,
     6  0.850764897055834977D-04,  0.425382448527917472D-04,
     7  0.212691224263958736D-04,  0.106345612131979372D-04/
 
      DATA (WP(N,1),WP(N,3),N=2,15)
     1/-0.133196159122085045D+01,  0.445816186556927292D-01,
     2 -0.229218106995884763D+01, -0.240054869684499309D-01,
     3 -0.311522633744855959D+01, -0.925925925925925875D-01,
     4 -0.380109739368998611D+01, -0.161179698216735251D+00,
     5 -0.434979423868312742D+01, -0.229766803840877915D+00,
     6 -0.476131687242798352D+01, -0.298353909465020564D+00,
     7 -0.503566529492455417D+01, -0.366941015089163228D+00,
     8 -0.517283950617283939D+01, -0.435528120713305891D+00,
     9 -0.517283950617283939D+01, -0.504115226337448555D+00,
     A -0.503566529492455417D+01, -0.572702331961591218D+00,
     B -0.476131687242798352D+01, -0.641289437585733882D+00,
     C -0.434979423868312742D+01, -0.709876543209876532D+00,
     D -0.380109739368998611D+01, -0.778463648834019195D+00,
     E -0.311522633744855959D+01, -0.847050754458161859D+00/
 
      RESULT=0
      ABSERR=0
      IFAIL=3
      IF(N .LT. 2 .OR. N .GT. 15) RETURN
      IF(MINPTS .GT. MAXPTS) RETURN
 
      IFNCLS=0
      LDV=.FALSE.
      TWONDM=2**N
      IRGNST=2*N+3
      IRLCLS=2**N+2*N*(N+1)+1
      ISBRGN=IRGNST
      ISBRGS=IRGNST
      IF(MAXPTS .LT. IRLCLS) RETURN
      DO 10 J = 1,N
      CTR(J)=(B(J)+A(J))*HF
   10 WTH(J)=(B(J)-A(J))*HF
 
   20 RGNVOL=TWONDM
      DO 30 J = 1,N
      RGNVOL=RGNVOL*WTH(J)
   30 Z(J)=CTR(J)
      SUM1=F(N,Z)
 
      DIFMAX=0
      SUM2=0
      SUM3=0
      DO 40 J = 1,N
      Z(J)=CTR(J)-XL2*WTH(J)
      F2=F(N,Z)
      Z(J)=CTR(J)+XL2*WTH(J)
      F2=F2+F(N,Z)
      WTHL(J)=XL4*WTH(J)
      Z(J)=CTR(J)-WTHL(J)
      F3=F(N,Z)
      Z(J)=CTR(J)+WTHL(J)
      F3=F3+F(N,Z)
      SUM2=SUM2+F2
      SUM3=SUM3+F3
      DIF=ABS(7*F2-F3-12*SUM1)
      DIFMAX=MAX(DIF,DIFMAX)
      IF(DIFMAX .EQ. DIF) IDVAXN=J
   40 Z(J)=CTR(J)
 
      SUM4=0
      DO 70 J = 2,N
      J1=J-1
      DO 60 K = J,N
      DO 50 L = 1,2
      WTHL(J1)=-WTHL(J1)
      Z(J1)=CTR(J1)+WTHL(J1)
      DO 50 M = 1,2
      WTHL(K)=-WTHL(K)
      Z(K)=CTR(K)+WTHL(K)
   50 SUM4=SUM4+F(N,Z)
   60 Z(K)=CTR(K)
   70 Z(J1)=CTR(J1)
 
      SUM5=0
      DO 80 J = 1,N
      WTHL(J)=-XL5*WTH(J)
   80 Z(J)=CTR(J)+WTHL(J)
   90 SUM5=SUM5+F(N,Z)
      DO 100 J = 1,N
      WTHL(J)=-WTHL(J)
      Z(J)=CTR(J)+WTHL(J)
      IF(WTHL(J) .GT. 0) GO TO 90
  100 CONTINUE
 
      RGNCMP=RGNVOL*(WP(N,1)*SUM1+WP2*SUM2+WP(N,3)*SUM3+WP4*SUM4)
      RGNVAL=W(N,1)*SUM1+W2*SUM2+W(N,3)*SUM3+W4*SUM4+W(N,5)*SUM5
      RGNVAL=RGNVOL*RGNVAL
      RGNERR=ABS(RGNVAL-RGNCMP)
      RESULT=RESULT+RGNVAL
      ABSERR=ABSERR+RGNERR
      IFNCLS=IFNCLS+IRLCLS
 
      IF(LDV) THEN
  110  ISBTMP=2*ISBRGN
       IF(ISBTMP .GT. ISBRGS) GO TO 160
       IF(ISBTMP .LT. ISBRGS) THEN
        ISBTPP=ISBTMP+IRGNST
        IF(WK(ISBTMP) .LT. WK(ISBTPP)) ISBTMP=ISBTPP
       ENDIF
       IF(RGNERR .GE. WK(ISBTMP)) GO TO 160
       DO 130 K = 0,IRGNST-1
  130  WK(ISBRGN-K)=WK(ISBTMP-K)
       ISBRGN=ISBTMP
       GO TO 110
      ENDIF
  140 ISBTMP=(ISBRGN/(2*IRGNST))*IRGNST
      IF(ISBTMP .GE. IRGNST .AND. RGNERR .GT. WK(ISBTMP)) THEN
       DO 150 K = 0,IRGNST-1
  150  WK(ISBRGN-K)=WK(ISBTMP-K)
       ISBRGN=ISBTMP
       GO TO 140
      ENDIF
 
  160 WK(ISBRGN)=RGNERR
      WK(ISBRGN-1)=RGNVAL
      WK(ISBRGN-2)=IDVAXN
      DO 170 J = 1,N
      ISBTMP=ISBRGN-2*J-2
      WK(ISBTMP+1)=CTR(J)
  170 WK(ISBTMP)=WTH(J)
      IF(LDV) THEN
       LDV=.FALSE.
       CTR(IDVAX0)=CTR(IDVAX0)+2*WTH(IDVAX0)
       ISBRGS=ISBRGS+IRGNST
       ISBRGN=ISBRGS
       GO TO 20
      ENDIF
      RELERR=ABSERR/ABS(RESULT)
      IF(ISBRGS+IRGNST .GT. IWK) IFAIL=2
      IF(IFNCLS+2*IRLCLS .GT. MAXPTS) IFAIL=1
      IF(RELERR .LT. EPS .AND. IFNCLS .GE. MINPTS) IFAIL=0
      IF(IFAIL .EQ. 3) THEN
       LDV=.TRUE.
       ISBRGN=IRGNST
       ABSERR=ABSERR-WK(ISBRGN)
       RESULT=RESULT-WK(ISBRGN-1)
       IDVAX0=WK(ISBRGN-2)
       DO 190 J = 1,N
       ISBTMP=ISBRGN-2*J-2
       CTR(J)=WK(ISBTMP+1)
  190  WTH(J)=WK(ISBTMP)
       WTH(IDVAX0)=HF*WTH(IDVAX0)
       CTR(IDVAX0)=CTR(IDVAX0)-WTH(IDVAX0)
       GO TO 20
      ENDIF
      NFNEVL=IFNCLS
      RETURN
      END


      subroutine qgaus(func,a,b,ss)
      real a,b,ss,func
      external func
      integer j
      real dx,xm,xr,w(5),x(5)
      save w,x
      data w/0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 
     &       0.0666713443/
      data x/0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666,
     &       0.9739065285/

      xm = 0.5*(b+a)
      xr = 0.5*(b-a)
      ss = 0.0
      do j =1,5
         dx = xr*x(j)
         ss = ss + w(j)*(func(xm+dx)+func(xm-dx))
      enddo
      ss = xr*ss
      return
      end
