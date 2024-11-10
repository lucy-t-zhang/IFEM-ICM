      subroutine smaterj(wto,toc,xmi,xmj,dxmj,ddxmj)
      implicit real*8 (a-h,o-z) 
      include 'common'
      dimension xmj(3),xmi(3),toc(3,3),&
     	dli(3,6),ddli(3,6,6),dxmj(3,6),ddxmj(3,6,6)
      cc=0.0d0
      do i=1,2
         do j=1,2
            cc=cc+toc(i,j)*toc(i,j)
         enddo
	enddo
      cc=cc+1.0d0
      xmi(1)=toc(1,1)+toc(2,2)+1.0d0
      xmi(2)=0.5d0*(xmi(1)**2-cc)
      xmi(3)=toc(1,1)*toc(2,2)-toc(1,2)*toc(2,1)
      xmi1=dexp(-x13*dlog(xmi(3)))
      xmi2=dexp(-x23*dlog(xmi(3)))
      xmi4=dexp(-x43*dlog(xmi(3)))
      xmi5=dexp(-x53*dlog(xmi(3)))
      xmi7=dexp(-x73*dlog(xmi(3)))
      xmi8=dexp(-x83*dlog(xmi(3)))
      xmj(1)=xmi(1)*xmi1
      xmj(2)=xmi(2)*xmi2
      xmj(3)=dsqrt(xmi(3))
      wto=rc1*(xmj(1)-3.0d0)+rc2*(xmj(2)-3.0d0)+&
          0.5d0*rk*(xmj(3)-1.0d0)**2
      do i=1,3
         do j=1,4
            dli(i,j)=0.0d0
            do m=1,4
               ddli(i,m,j)=0.0d0
            enddo
		enddo
	enddo
      dli(3,1)=2.0d0*toc(2,2)
      dli(3,2)=2.0d0*toc(1,1)
      dli(3,3)=-2.0d0*toc(2,1)
      dli(3,4)=2.0d0*(toc(1,1)*toc(2,2)-toc(2,1)*toc(1,2))
      dli(1,1)=2.0d0
      dli(1,2)=2.0d0
      dli(1,4)=2.0d0
      dli(2,1)=2.0d0+2.0d0*toc(2,2)
      dli(2,2)=2.0d0+2.0d0*toc(1,1)
      dli(2,3)=-2.0d0*toc(1,2)
      dli(2,4)=2.0d0*(toc(2,2)+toc(1,1))
      ddli(3,1,2)=4.0d0
      ddli(3,2,1)=4.0d0
      ddli(3,3,3)=-2.0d0
      ddli(2,1,2)=4.0d0
      ddli(2,2,1)=4.0d0
      ddli(2,3,3)=-2.0d0
      do i=1,3
         dxmj(1,i)=dli(1,i)*xmi1-xmi(1)*xmi4*dli(3,i)*x13
         dxmj(2,i)=dli(2,i)*xmi2-xmi(2)*xmi5*dli(3,i)*x23
         dxmj(3,i)=0.5d0/dsqrt(xmi(3))*dli(3,i)
         do j=1,3
            ddxmj(1,i,j)=-x13*xmi4*(dli(1,i)*dli(3,j)+&
                dli(1,j)*dli(3,i)+ddli(3,i,j)*xmi(1))+&
                x49*dli(3,j)*dli(3,i)*xmi7*xmi(1)
            ddxmj(2,i,j)=ddli(2,i,j)*xmi2-xmi5*x23*&
               (dli(2,i)*dli(3,j)+dli(2,j)*dli(3,i)+&
               ddli(3,i,j)*xmi(2))+dli(3,j)*dli(3,i)*xmi8*xmi(2)*x109         
            ddxmj(3,i,j)=-0.25d0*dli(3,i)*dli(3,j)/xmi(3)/dsqrt(xmi(3))+&
                0.5d0*ddli(3,i,j)/dsqrt(xmi(3))
		enddo
	enddo
      dxmj(1,4)=dli(1,4)*xmi1-xmi(1)*xmi4*dli(3,4)*x13
      dxmj(2,4)=dli(2,4)*xmi2-xmi(2)*xmi5*dli(3,4)*x23
      dxmj(3,4)=0.5d0/dsqrt(xmi(3))*dli(3,4)
      return
      end