      subroutine cal_boost_force
      include 'moldy.h'
c
      integer*4 i, j, ii, ihe, ijk
      real*8 xijFe, yijFe, zijFe, xijHe, yijHe, zijHe
      real*8 dispmaxFe, dispmaxHe
      integer*4 idya1, iikey

c sum all Fe displacement
      sumdisp = 0.0
      do i = 1, idya
         j = idya_n(i)
         idj = id(j)
         if(idj .eq. 1) sumdisp = sumdisp + displace(i)
      end do
      if(nacmd .eq.3000)then
        if(sumdisp.lt.rcutq-1.0 .or. sumdisp.ge.rcutq) then
          rcutq = sumdisp + 0.1
        endif
      endif
c calculate the boost energy and force
      do i = 1, nn1
         fbstx(i) = 0.0
         fbsty(i) = 0.0
         fbstz(i) = 0.0
      end do
c check the sumdisp
      diffdisp = sumdisp - rcutq
      vboost = 0.0
      if(sumdisp.le.rcutq)then
        vboost = sigma*(1.0-(sumdisp/rcutq)**2)
        forcefactor = 2.0*sigma*(sumdisp)/rcutq**2
        if(vboost .gt. vebcric)then
          vcfact = vboost/vebcric
          forcefactor = forcefactor/vcfact
          vboost = vebcric
        endif
        do i = 1,idya
           j = idya_n(i)
           xij = x0(j)-x0f(j)
           yij = y0(j)-y0f(j)
           zij = z0(j)-z0f(j)
           fbstx(j) = forcefactor * xij/displace(i)
           fbsty(j) = forcefactor * yij/displace(i)
           fbstz(j) = forcefactor * zij/displace(i)
        end do
      endif
      najumpsum = najumpsum + 1
      if(nacmd .gt.0.and.mod(najumpsum,2000).eq.0)then
         najump = 0
      endif
      write(*,*) nacmd, najump, vboost, rcutq, sigma
      if(najump.eq.0 .and. vboost.eq.0 .and. 
     &   sumdisp.ge.rcutq+rinc)then
          rcutq = rcutq + rinc
          sigma = sigma + sigmainc
          najump = 1
          najumpsum = 0
      endif
     
c output the system boost energy and the force 
      write(*,*)  'sumdisp = ', sumdisp, vboost, rcutq, sigma
c      do i = 1, nn1
c         if(fbstx(i).ne.0.0 .or. fbsty(i).ne.0.0 .or. 
c     &      fbstz(i).ne.0.0)then
c            write(*,*) 'force- ', i, fbstx(i), fbsty(i), fbstz(i)
c         endif
c      end do
c calculate the time
      timeboost = timeboost + exp(vboost/bk/(temprq/2.0))
      write(*,*) 'timeboost = ',timeboost
      return
      end subroutine 
