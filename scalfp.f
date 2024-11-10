      subroutine scalfp(fp,ocpp,i)
      implicit real*8 (a-h,o-z) 
      include 'common'
      fp=-ocpp*(bpre-cpre)*hp(i)
      return
      end