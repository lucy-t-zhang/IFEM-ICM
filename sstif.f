      subroutine sstif(ocpp,ocuu,ocup,ne,w,toxj)
      implicit real*8 (a-h,o-z)
      include 'common'
      dimension ocuu(6,6),ocup(6),tram(2),toxj(2,2)
      do i=1,2
		xac(i)=0.0d0
          do k=1,nis(2)
			ntem=nea(ne,k)
              xac(i)=xac(i)+h(k)*acm(i,ntem,1)
		enddo
	enddo
	do i=1,2
		xac(i)=xac(i)-fbacco(i)
	enddo
      do ni=1,nis(2)
		nu1=2*(nea(ne,ni)-1)+1
          nv1=2*(nea(ne,ni)-1)+2
          drf(id(nu1))=drf(id(nu1))-w*sdensi*h(ni)*xac(1)
          drf(id(nv1))=drf(id(nv1))-w*sdensi*h(ni)*xac(2)
          sk(id(nu1))=sk(id(nu1))+w*sdensi*h(ni)*h(ni)/beta/dt
          sk(id(nv1))=sk(id(nv1))+w*sdensi*h(ni)*h(ni)/beta/dt
      enddo
      do ni=1,nis(ig)
		do i=1,2
			nu1=2*(nea(ne,ni)-1)+i
			call scalfu(fu,i,ni)
              drf(id(nu1))=drf(id(nu1))-fu*w
              call scalkuu(xkuu,ocuu,i,ni)
              sk(id(nu1))=sk(id(nu1))+w*xkuu*alpha/beta*dt
		enddo
	enddo	
      do i=1,nump
		np1=nndtem+idp(neap(ne,i))
          call scalfp(fp,ocpp,i)
          drf(np1)=drf(np1)-fp*w
          call scalkpp(xkpps,ocpp,i)
          sk(np1)=sk(np1)+xkpps*w
      enddo
      return
      end        