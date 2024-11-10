      subroutine spiola(ocpp,xmj,dxmj)
      implicit real*8 (a-h,o-z) 
      include 'common'
      dimension xmj(3),dxmj(3,6)
      do i=1,4
         btos(i)=rc1*dxmj(1,i)+rc2*dxmj(2,i)+
     $        rk*(xmj(3)-1.0d0)*dxmj(3,i)
         tos(i)=btos(i)+ocpp*(bpre-cpre)*dbpre(i)
      enddo
      return
      end