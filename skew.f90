      subroutine skew
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension t(2,2)
      do i=1,intnum
         nd=numint(i)
	   nd1=id(2*(nd-1)+1)
	   nd2=id(2*(nd-1)+2)
         ns=ninsk(i)
         xa=xang(ns)*pi/180.0d0
         t(1,1)=dcos(xa)
         t(2,1)=-dsin(xa)
         t(1,2)=dsin(xa)
         t(2,2)=dcos(xa)
         ci=t(1,1)*drf(nd1)+t(1,2)*drf(nd2)
         cj=t(2,1)*drf(nd1)+t(2,2)*drf(nd2)
         drf(nd1)=ci
         drf(nd2)=cj
         ai=t(1,1)**2*sk(nd1)+t(1,2)**2*sk(nd2)
         aj=t(2,1)**2*sk(nd1)+t(2,2)**2*sk(nd2)
         sk(nd1)=ai
         sk(nd2)=aj
      enddo
      return
      end