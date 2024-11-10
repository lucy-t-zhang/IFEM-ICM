      subroutine sboc(obc,ocpp,ocuu,ocup,xmj,dxmj,ddxmj)
      implicit real*8 (a-h,o-z) 
      include 'common'
      dimension obc(6,6),xmj(3),dxmj(3,6),ddxmj(3,6,6),&
          ocuu(6,6),ocup(6)
      ocpp=-1.0d0/rk
      do i=1,3
         ocup(i)=-ocpp*dbpre(i)
         do j=1,3
            obc(i,j)=rc1*ddxmj(1,i,j)+rc2*ddxmj(2,i,j)+&
                rk*dxmj(3,i)*dxmj(3,j)+rk*ddxmj(3,i,j)*(xmj(3)-1.0d0)
            ocuu(i,j)=obc(i,j)+ocpp*dbpre(i)*dbpre(j)+&
               ocpp*(bpre-cpre)*ddbpre(i,j)
		enddo
	enddo
      return
      end