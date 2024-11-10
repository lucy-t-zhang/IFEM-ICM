      subroutine antiskew
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension t(2,2)
      do i=1,intnum
		nd=numint(i)
		ns=ninsk(i)
		xa=xang(ns)*pi/180.0d0
		t(1,1)=dcos(xa)
		t(2,1)=dsin(xa)
		t(1,2)=-dsin(xa)
		t(2,2)=dcos(xa)
		ai=t(1,1)*dui(1,nd)+t(1,2)*dui(2,nd)
		aj=t(2,1)*dui(1,nd)+t(2,2)*dui(2,nd)
		dui(1,nd)=ai
		dui(2,nd)=aj
      enddo
      return
      end