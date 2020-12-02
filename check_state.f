************************************************************************
      subroutine check_state
      include 'moldy.h'
c
c     Local variables

      do i = 1, nn1
         if(pbcx)then
           if(x0(i) .ge. 0.5) x0(i) = x0(i) - 1.0
           if(x0(i) .lt.-0.5) x0(i) = x0(i) + 1.0
         endif
         if(pbcy)then
           if(y0(i) .ge. 0.5) y0(i) = y0(i) - 1.0
           if(y0(i) .lt.-0.5) y0(i) = y0(i) + 1.0
         endif
         if(pbcz)then
           if(z0(i) .ge. 0.5) z0(i) = z0(i) - 1.0
           if(z0(i) .lt.-0.5) z0(i) = z0(i) + 1.0
         endif
      end do

c calculate the distance to the fixed position
      distmax = 0.0
      do j = 1, idya
         i = idya_n(j)
         xij = (x0(i) - x0f(i))*b0(1,1)
         yij = (y0(i) - y0f(i))*b0(2,2)
         zij = (z0(i) - z0f(i))*b0(3,3)

         if(xij > b0(1,1)/2.0d0) xij = xij - b0(1,1)
         if(xij <-b0(1,1)/2.0d0) xij = xij + b0(1,1)
         if(yij > b0(2,2)/2.0d0) yij = yij - b0(2,2)
         if(yij <-b0(2,2)/2.0d0) yij = yij + b0(2,2)
         if(zij > b0(3,3)/2.0d0) zij = zij - b0(3,3)
         if(zij <-b0(3,3)/2.0d0) zij = zij + b0(3,3)

         dist = sqrt(xij*xij + yij*yij + zij*zij)
         !dist = sqrt(yij*yij + zij*zij)
         if(dist .gt. 0.4) write(*,*) 'dist = ',i, dist
         if(distmax .le. dist) then
           distmax = dist
           idjump = i
         endif
      end do
c for vacancy case, max distance > 1.3
      if(distmax .ge. rstop)then
        write(*,*) 'jump with distance ', distmax, timeboost, idjump
c        vacposx = x0f(idjump)
c        vacposy = y0f(idjump)
c        vacposz = z0f(idjump)
        timeboost = 0.0
        nacmd = 0
        do_boost = .false.
        jumptime = jumptime + 1
        call rite
      endif
c output the boost atomic position
c Ning Gao: write the OVITO readin format
      if(mod(nst,100).eq.0)then
        open(unit=21,file='mdyn_region.psi',position='append')
        write(21,'(a)') 'ITEM: TIMESTEP'
        write(21,'(i9)') nst
        write(21,'(a)') 'ITEM: NUMBER OF ATOMS'
        write(21,'(i6)') idya
        write(21,'(a)') 'ITEM: BOX BOUNDS'
        write(21,'(2(f18.9))') -b0(1,1)/2.0d0, b0(1,1)/2.0d0
        write(21,'(2(f18.9))') -b0(2,2)/2.0d0, b0(2,2)/2.0d0
        write(21,'(2(f18.9))') -b0(3,3)/2.0d0, b0(3,3)/2.0d0
        write(21,'(a)') 'ITEM: ATOMS id x y z itype pe ke ic'
        do 12211 j = 1, idya
           iij = idya_n(j)
           xiij = x0(iij)*b0(1,1)
           yiij = y0(iij)*b0(2,2)
           ziij = z0(iij)*b0(3,3)
           write(21,'(i6,2x,3(f12.5),2x,i2,2(f12.5),1x,i3)')
     &     iij,xiij,yiij,ziij,id(iij),atpe(iij),ekin(iij),
     &     ibelongclusterid(j)
12211   continue
        close(21)
      endif

      return
      end
