      subroutine surfree(nes)
      implicit real*8 (a-h,o-z)
      include 'common'
      do lx=1,nint
		r=xg(lx,nint)
		call belement(r,det,nes)
		w=wgt(lx,nint)*det
		ydis=0.0d0
		yddis=0.0d0
		do k=1,niss
			ntem=nfsnd(nes,k)
              ydis=ydis+hfs(k)*dis(2,ntem,1)
              yddis=yddis+hfsp(k)*dis(2,ntem,1)/det
		enddo
		do i=1,niss
			ni=2*nnd+nfsnd(nes,i)
			sk(ni)=sk(ni)+alpha*dt/beta*xmat*w*
     $			hfs(i)*hfs(i)
              drf(ni)=drf(ni)-w*hfs(i)*ydis*xmat
		enddo
      enddo
      return
      end