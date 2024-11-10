      subroutine surfluid(w,ne)
      implicit real*8 (a-h,o-z)
      dimension bdup(2,9)
      include 'common'
      tem1=w*fdensi
      tem3=w/fbulk
      tem5=fdensi*w/falph
      call acc(ne)
      if (nsmoo .eq. 1) then
         do i=1,nis(ig)
            bdup(1,i)=(bd(1,i)*defv(1)+bd(2,i)*defv(2))*upw
         enddo
      endif
      bpres=0.0d0
      bpresv=0.0d0       
      do i=1,nump
		ni=neap(ne,i)
          bpres=bpres+epc(ni,1)*hp(i)
          bpresv=bpresv+epcv(ni,1)*hp(i)
      enddo
      xelep=xdve(1,1)+xdve(2,2)+bpresv/fbulk
      do i=1,nump
            xkpp(i,ne)=xkpp(i,ne)-tem3*hp(i)*hp(i)/beta/dt
      enddo
      if ((nsmoo .ne. 5) .and. (ncon .ne. 1)) then
         do ni=1,nis(ig)
               xfsk(ni)=xfsk(ni)+
     $              tem1*h(ni)*(defv(1)*bd(1,ni)+
     $              defv(2)*bd(2,ni))
               xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+
     $              tem1*h(ni)*(defv(1)*bd(1,ni)+defv(2)*bd(2,ni))  
	   enddo
      endif
      do ni=1,nis(ig)
            xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+h(ni)*h(ni)*tem1/beta/dt
            xfsk(ni)=xfsk(ni)+h(ni)*h(ni)*tem1/beta/dt
      enddo
      do ni=1,nis(ig)
		xfsk(ni)=xfsk(ni)+tem5*(bd(1,ni)*bd(1,ni)+
     $           bd(2,ni)*bd(2,ni))
          xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+
     $           tem5*(bd(1,ni)*bd(1,ni)+bd(2,ni)*bd(2,ni))
      enddo
      if (nsmoo .eq. 4) then
c     isotropic
         if (niso .eq. 1) then
            do ni=1,nis(ig)
                  xfsk(ni)=xfsk(ni)+xuk(1,1)*
     $                 tem1*(bd(1,ni)*bd(1,ni)+
     $                 bd(2,ni)*bd(2,ni))
                  xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+
     $                 tem1*xuk(1,1)*(bd(1,ni)*bd(1,ni)+
     $                 bd(2,ni)*bd(2,ni))
            enddo
         else
c     anisotropic
            do ni=1,nis(ig)
				xfsk(ni)=xfsk(ni)+tem1*(bd(1,ni)*(xuk(1,1)*
     $                 bd(1,ni)+xuk(1,2)*bd(2,ni))+bd(2,ni)*
     $                 (xuk(2,1)*bd(1,ni)+xuk(2,2)*bd(2,ni)))
                  xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+
     $                 tem1*(bd(1,ni)*(xuk(1,1)*bd(1,ni)+xuk(1,2)*
     $                 bd(2,ni))+bd(2,ni)*(xuk(2,1)*bd(1,ni)+
     $                 xuk(2,2)*bd(2,ni)))   
            enddo
         endif
      endif
      if (nsmoo .eq. 1) then
         do ni=1,nis(ig)
               xfsk(ni)=xfsk(ni)+
     $              tem1*bdup(1,ni)*(h(nj)*xdve(1,1)+
     $              defv(1)*bd(1,ni)+defv(2)*bd(2,ni))
               xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+
     $              tem1*bdup(1,ni)*(h(nj)*xdve(2,2)+
     $              defv(1)*bd(1,ni)+defv(2)*bd(2,ni))
         enddo         
         do ni=1,nis(ig)
               xfsk(ni+nis(ig))=xfsk(ni+nis(ig))+
     $              bdup(1,ni)*h(ni)*tem1/beta/dt
         enddo
         do ni=1,nis(ig)
               xfsk(ni)=xfsk(ni)+bdup(1,ni)*h(ni)*tem1/beta/dt  
         enddo
       endif
      do nj=1,nis(ig)
         ntem=nea(ne,nj)
         do j=1,n2link(ntem)
            nu2=nod2link(j,ntem)
            nu2tem=2*(nu2-1)+2
            tem4=w*coe2link(j,ntem)
            nu2t=nnni(nod2link(j,ntem))
            if (nu2t .gt. 0) then
               nu2tem=nu2t+2*nnd
            endif
            sk(nu2tem)=sk(nu2tem)-fdensi*
     $              tem4*h(ni)*h(ni)*xdve(1,2)
            sk(nu2tem)=sk(nu2tem)+tem4*alpha*dt/beta*
     $              bd(2,ni)*(h(ni)*xac(1)*fdensi-
     $              bpres*bd(1,ni)+
     $              tem5*(bd(1,ni)*xdve(1,1)+
     $              bd(2,ni)*xdve(1,2)))
            sk(nu2tem)=sk(nu2tem)+tem4*alpha*dt/beta*
     $              bd(1,ni)*bd(2,ni)*bpres
            if ((nsmoo .ne. 5) .and. (ncon .ne. 1)) then
                  sk(nu2tem)=sk(nu2tem)+
     $                 fdensi*tem4*alpha*dt/beta*h(ni)*
     $                 (-defv(1)*xdve(1,2)*bd(1,ni)-
     $                 defv(2)*xdve(1,2)*bd(2,ni))
            endif
            if (nsmoo .eq. 1) then
               do ni=1,nis(ig)
                  nu1=nea(ne,ni)
				if (nu1 .eq. nu2tem) then
                  sk(nu1)=sk(nu1)-fdensi*
     $                 tem4*bdup(1,ni)*h(nj)*xdve(1,2)
                  sk(nu1)=sk(nu1)+tem4*alpha*dt/beta*
     $                 bd(2,nj)*bdup(1,ni)*(xac(1)*fdensi)
                  sk(nu1)=sk(nu1)+
     $                 tem4*alpha*dt/beta*bdup(1,ni)*fdensi*
     $                 (-defv(1)*xdve(1,2)*bd(1,nj)-
     $                 defv(2)*xdve(1,2)*bd(2,nj))
                  sk(nu1)=sk(nu1)-
     $                 bd(2,ni)*tem4*alpha*dt/beta*bdup(1,nj)*
     $                 (fdensi*xac(1))
				endif
               enddo
            endif               
ccccccccccccccccccccccccccccccccccc
c     Artificial viscosity
ccccccccccccccccccccccccccccccccccc
            if (nsmoo .eq. 4) then
               if (niso .eq. 1) then
                  do ni=1,nis(ig)
                     nu1=nea(ne,ni)
                     xy1=xuk(1,1)*tem4*alpha*dt/beta*fdensi
				if (nu1 .eq. nu2tem) then
                     sk(nu1)=sk(nu1)-xy1*bd(2,ni)*
     $                    (bd(1,nj)*xdve(1,1)+bd(2,nj)*xdve(1,2))     
                     sk(nu1)=sk(nu1)-xy1*
     $                    xdve(1,2)*(bd(1,ni)*bd(1,nj)+
     $                    bd(2,ni)*bd(2,nj))
                     sk(nu1)=sk(nu1)+xy1*
     $                    (bd(1,ni)*xdve(1,1)+bd(2,ni)*
     $                    xdve(1,2))*bd(2,nj)
				endif
                enddo
               else
                  do ni=1,nis(ig)
                     nu1=nea(ne,ni)
                     xy2=tem4*alpha*dt/beta*fdensi
				if (nu1 .eq. nu2tem) then
                     sk(nu1)=sk(nu1)+xy2*
     $                    bd(2,nj)*(bd(1,ni)*(xuk(1,1)*xdve(1,1)+
     $                    xuk(1,2)*xdve(1,2))+bd(2,ni)*(xuk(2,1)*
     $                    xdve(1,1)+xuk(2,2)*xdve(1,2)))     
                     sk(nu1)=sk(nu1)-xy2*
     $                       bd(2,ni)*(bd(1,nj)*(xdve(1,1)*xuk(1,1)+
     $                    xdve(1,2)*xuk(1,2))+bd(2,nj)*(xdve(1,1)*
     $                    xuk(2,1)+xdve(1,2)*xuk(2,2)))     
                     sk(nu1)=sk(nu1)-xy2*
     $                    xdve(1,2)*(bd(1,ni)*(xuk(1,1)*bd(1,nj)+
     $                    xuk(1,2)*bd(2,nj))+bd(2,ni)*(xuk(2,1)*
     $                    bd(1,nj)+xuk(2,2)*bd(2,nj)))     
				endif
                 enddo
               endif
            endif
        enddo
      enddo
      do i=1,nis(ig)
         xfu(i)=xfu(i)-(tem1*h(i)*xac(1)-w*bpres*bd(1,i)+
     $        tem5*(bd(1,i)*xdve(1,1)+bd(2,i)*xdve(1,2)))
         xfu(i+nis(ig))=xfu(i+nis(ig))-(tem1*h(i)*xac(2)-
     $        w*bpres*bd(2,i)+tem5*(bd(1,i)*xdve(2,1)+
     $           bd(2,i)*xdve(2,2)))
      enddo
      do i=1,nump
         xfp(i,ne)=xfp(i,ne)-w*hp(i)*xelep
      enddo
      if (nsmoo .eq. 1) then
         do i=1,nis(ig)
            xfu(i)=xfu(i)-tem1*bdup(1,i)*xac(1)
            xfu(i+nis(ig))=xfu(i+nis(ig))-tem1*bdup(1,i)*xac(2)
         enddo
      endif
      if (nsmoo .eq. 4) then
         do i=1,nis(ig)
            if (niso .eq. 1) then
               xfu(i)=xfu(i)-tem1*xuk(1,1)*(bd(1,ni)*
     $              xdve(1,1)+bd(2,ni)*xdve(1,2))
               xfu(i+nis(ig))=xfu(i+nis(ig))-tem1*xuk(1,1)*
     $              (bd(1,ni)*xdve(2,1)+bd(2,ni)*xdve(2,2))
            else
               xfu(i)=xfu(i)-tem1*(bd(1,ni)*(xuk(1,1)*
     $              xdve(1,1)+xuk(1,2)*xdve(1,2))+bd(2,ni)*
     $              (xuk(2,1)*xdve(1,1)+xuk(2,2)*xdve(1,2)))
               xfu(i+nis(ig))=xfu(i+nis(ig))-tem1*(bd(1,ni)*
     $              (xuk(1,1)*xdve(2,1)+xuk(1,2)*xdve(2,2))+
     $              bd(2,ni)*(xuk(2,1)*xdve(2,1)+
     $              xuk(2,2)*xdve(2,2)))
            endif
         enddo
      endif
      return
      end