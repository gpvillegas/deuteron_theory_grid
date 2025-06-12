*      program test
*      common/par/pi,pm,dm
*      dimension q(3),p2(3),p1(3)
*      real je_re, je_im
*      dimension je_re(2,2),je_im(2,2)
*      complex je
*      dimension je(2,2)
*      pi = acos(-1.0)
*      pm = 0.938279
*      dm = 1.875
*      kn = 1
*      q2 = 1.0
*      q(1) = 0.0; q(2) = 0.0
*      q(3)   = 1.5 
*      p1(1) = 0.0; p1(2) = 0.0; p1(3) = 0.0
*      p2(1) = 0.0; p2(2) = 0.0
*      p2(3) = q(3)
*      mu = 0
*      call  j_e(p2,p1,je,mu)
*      je_re = real(je)
*      je_im = imag(je)
      
*      print *,je_re
*      print *,je_im
*      print *
*      print *,je
*      end



      subroutine j_e(p2,p1,je,mu)
*************************************************************************************
* Calculates electromagnetic current of electron
*  \bar u(s2,p2)|\gamma^\mu |u(p1,s1)  = j_re/im(s2,s1), where
*  1,1 = u2,u1, 1,2=u2,d1,  2,1=d2,u1, 2,2=d2,d1
*
*
* p1 - p1(3) three momentum of the initial nucleon
* p2 - p2(3) three momentum of the final nucleon 
* mu - 0,1,2,3, - components
* je - 2x2 complex matrix out - returns re and im part of the current
*
*
*
* Misak Sargsian
* FIU, Miami, FL
* July, 2005
* 11-July-2005 - last updated
*
****************************************************************************************

 
      implicit none
      complex, dimension(2,2) :: je
      real, dimension(2,2) :: je_re, je_im
      real, dimension(3) :: p2,p1,datark
      real, dimension(3) :: p2xp1
      real, dimension(2,2) :: amre,amim
      real :: p2_2,p1_2,p2p1,p2xp1m
      real :: em,e2,e1,e2m,e1m,e21m,fc
      integer :: x=1,y=2,z=3,mu
      em  = 0.000511
      datark= 0.0
      call prode(p2,p2,datark,p2_2,0)
      e2 = sqrt(em**2 + p2_2)
      call prode(p1,p1,datark,p1_2,0)
      e1 = sqrt(em**2 + p1_2)

*************************************************
*      Scalar products
*************************************************
      call prode(p2,p1,datark,p2p1,0)

*************************************************
*      Vector Products
*************************************************
      call prode(p2,p1,p2xp1,p2xp1m,1)


      e2m  = e2+em
      e1m  = e1+em
      e21m = e2m*e1m
      fc   = sqrt(e21m)


*********************************************
*  In the matrix (s2,s1)
!************ mu = 0 ************************
      if(mu.eq.0)then
       amre(1,1) = 1 + p2p1/e21m; amre(1,2) = p2xp1(y)/e21m
       amre(2,1) =    -amre(1,2); amre(2,2) =     amre(1,1)
       
       je_re = fc * amre

       amim(1,1)= p2xp1(z)/e21m;  amim(1,2) = p2xp1(x)/e21m
       amim(2,1)=     amim(1,2);  amim(2,2) =     -amim(1,1)

       je_im = fc * amim
       elseif(mu.eq.x)then 
       amre(1,1) = p2(x)/e2m + p1(x)/e1m; amre(1,2)=p2(z)/e2m-p1(z)/e1m
       amre(2,1) =            -amre(1,2); amre(2,2)=          amre(1,1)

       je_re = fc * amre

       amim(1,1) = -p2(y)/e2m + p1(y)/e1m; amim(1,2) = 0.0
       amim(2,1) = 0.0;                    amim(2,2) = - amim(1,1)

       je_im = fc * amim
      elseif(mu.eq.y)then          
       amre(1,1)  = p2(y)/e2m + p1(y)/e1m;  amre(1,2) = 0.0
       amre(2,1)  = 0.0;                    amre(2,2) = amim(1,1)

       je_re = fc * amre
       amim(1,1) = p2(x)/e2m - p1(x)/e1m; amim(1,2)=-p2(z)/e2m+p1(z)/e1m
       amim(2,1) =             amim(1,2); amim(2,2)=          -amim(1,1)

       je_im = fc * amim
      elseif(mu.eq.z)then          
       amre(1,1) = p2(z)/e2m + p1(z)/e1m; amre(1,2)=-p2(x)/e2m+p1(x)/e1m
       amre(2,1) =            -amre(1,2); amre(2,2)=           amre(1,1)
       je_re = fc * amre

       amim(1,1)  = 0.0;       amim(1,2)  = p2(y)/e2m - p1(y)/e1m;  
       amim(2,1)  = amim(1,2); amim(2,2)  = 0.0

       je_im = fc * amim
       endif

       je = cmplx(je_re,je_im)
      return
      end


      subroutine prode(a,b,c,cm,i)
******************************************************************
*      calculating the different products of  vectors a and b
*
*      i = 0  calculates the scalar product gives out only cm 
*      i = 1  calculates the vector product
*      i = 2  calculates the circular product
*
*      c(1),c(2),c(3) - x,y,z component of the product vector
*      cm   - magnitude of the product vectors
*
*   Misak Sargsian
*   February, 2005 
*   FIU, Miami
*******************************************************************
      implicit none
      integer i
      real a,b,c,cm
      dimension a(3),b(3),c(3)
      c(1) = 0.0; c(2)=0.0; c(3) = 0.0
      cm   = 0.0
      if(i.eq.0)then
      cm = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
      elseif(i.eq.1)then
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1) 
      cm   = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
      elseif(i.eq.2)then
      c(1) = a(2)*b(3) + a(3)*b(2)
      c(2) = a(3)*b(1) + a(1)*b(3)
      c(3) = a(1)*b(2) + a(2)*b(1) 
      cm   = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
      endif
      return
      end

 
