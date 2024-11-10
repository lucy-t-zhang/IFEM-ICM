      subroutine surele(w,nes,det)
      implicit real*8 (a-h,o-z)
      include 'common'
      call suracc(nes,det)
      do i=1,niss
         nu1=2*nnd+neas(nes,i)
         drf(nu1)=drf(nu1)-
     $        w*(xvms-xv2)*hfs(i)
         if (ncons .ne. 1) then
            drf(nu1)=drf(nu1)-w*(xdsur*xv1)*hfs(i)
         endif
         if (nsmoos .eq. 4) then
            drf(nu1)=drf(nu1)-w*hfsp(i)/det*
     $           xdsur*xuk(1,1)
         endif
         if (nsmoos .eq. 1) then
            drf(nu1)=drf(nu1)-w*(xvms-xv2+xdsur*xv1)*
     $           hfsp(i)/det*upws
         endif
         sk(nu1)=sk(nu1)+w*hfs(i)*hfs(i)
         if (ncons .ne. 1) then
               sk(nu1)=sk(nu1)+xv1*
     $              w*hfsp(i)/det*hfs(i)*alpha*dt/beta
         endif
         if (nsmoos .eq. 1) then
               sk(nu1)=sk(nu1)+
     $              w*hfsp(i)/det*upws*hfs(i)
               sk(nu1)=sk(nu1)+xv1*
     $              w*hfsp(i)/det*hfsp(i)/det*upws*alpha*dt/beta
         endif
         if (nsmoos .eq. 4) then
               sk(nu1)=sk(nu1)+xuk(1,1)*w*
     $              hfsp(i)/det*hfsp(i)/det*alpha*dt/beta
               sk(nu1)=sk(nu1)+xdsur*
     $              hfsp(i)/det*alpha*dt/beta*xuk(1,1)*
     $              w*xdsur*hfsp(i)/det
         endif
	enddo
      return
      end