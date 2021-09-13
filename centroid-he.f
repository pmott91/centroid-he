      implicit real*8 (a-h, o-z)

c --- parameters are compile-time constants.

c --- this parameter is the number of "beads" in a "necklace".

c --- this parameter is the number of iterations in a block.  only one
c     snapshot per block is ever recorded; this is to allow the beads
c     plenty of time to randomize between snapshots.

c --- this stores the coordinates of the beads.  the first index is the
c     atom number (there are two "real" atoms in this system).  the second
c     index is the bead number for a given atom.  the third index is 1, 2,
c     or 3 for the x, y, or z coordinates, respectively.

      dimension path(2, 10000, 3)

c --- a common block sets up global variables that can be shared between
c     subroutines.  every subroutine that wants access to these variables
c     must have the same common statement near the beginning.

      common /rancom/ rscale, rstate(6)

      dimension rseed(6)

      common /params/ force, beta, zmass, delta

c --- read three simulation parameters from the input stream:

c     niter = # of iterations (each iteration consists of NBLOCK steps)
c     temp = temperature in kelvin
c     dist = distance between centroids in bohr

      write(6,*)'Enter niter, temp, dist, NBEADS, NBLOCK, delta, iskip:'
      read (5, *) niter, temp, dist, NBEADS, NBLOCK, delta, iskip

c --- set up random number generator.

      rscale=1.0d0/4294967088.0d0

      call rcheck

      do i=1, 6
         rseed(i)=12345.0d0
      end do

c      open (3, file='skip')

c      read (3, *) iskip

      do i=1, iskip
         call rskip(rseed)
      end do

      close (3)

      do j=1, 6
         rstate(j)=rseed(j)
      end do

c --- delta is the random motion step size.  a bead can move up to
c     plus or minus (delta/2.0) in each direction (x, y, and z).

c ---      delta=0.35d0

c --- set up the initial bead locations.

      do i=1, NBEADS

         do j=1, 3
            path(1, i, j)=0.0d0
            path(2, i, j)=0.0d0
         end do

c ------ the centers of mass of the two necklaces are on the x axis,
c        separated by the chosen distance.

         path (1, i, 1)=-0.5d0*dist
         path (2, i, 1)=0.5d0*dist

      end do

c --- calculate some physical quantities for the simulation.

      beta=1.0d0/(temp*3.1668329d-6)

      zmass=4.00260325d0*1822.8885d0

      force=zmass*dble(NBEADS)/beta**2

      write (6, *) 'beta = ', beta
      write (6, *) 'force = ', force
      write (6, *) ''

c --- initialize the potential energy curve.

      call vinit(r2min, bin)

c --- initialize the interatomic force curve

      call finit(r2min, bin)

c --- warm-up the simulation.

      trej=0.0

      do n=1, 1000

         if (mod(n, 100).eq.0) write (6, *) 'WARMUP: ', n

c ------ beginning of the loop over a single block.

         nrej=0

         do nb=1, NBLOCK*NBEADS

            call move(1, path, nrej,NBEADS)

            call move(2, path, nrej,NBEADS)

c ------ end of the loop over a single block.

         end do

         trej=trej+dble(nrej)

c --- end of the warm-up phase.

      end do

      trej=trej/(2.0d0*1000.0d0*dble(NBLOCK)*dble(NBEADS))

      write (6, *) ''
      write (6, *) 'warmup rejection fraction = ', trej

c --- this counts the number of rejected moves.

      trej=0.0

c --- start the simulation.

      do n=1, niter

         nrej=0

c ------ beginning of the loop over a single block.

         do nb=1, NBLOCK*NBEADS

            call move(1, path, nrej,NBEADS)

            call move(2, path, nrej,NBEADS)

c ------ end of the loop over a single block.

         end do

         trej=trej+dble(nrej)

c ------ compute a few different things at the end of every block.

         rmax=-1.0d0
         rmin=10000.0d0

c ------ compute the mean potential energy and mean interatomic force of the pair of atoms.

         v=0.0d0
         f=0.0d0

         do i=1, NBEADS

            xx=path(2, i, 1)-path(1, i, 1)
            yy=path(2, i, 2)-path(1, i, 2)
            zz=path(2, i, 3)-path(1, i, 3)

            r2=xx*xx+yy*yy+zz*zz

c --------- also remember the smallest and largest distance.

            if (r2.lt.rmin) rmin=r2
            if (r2.gt.rmax) rmax=r2

            v=v+vlu(r2)
            f=f+flu(r2)
            
c --------- and record a snapshot of both sets of beads.

c --------- a non-blank in column 6 indicates a continuation of the
c           preceding line in Fortran...

c --------- replace the comment characters with spaces to activate these
c           write statements.

           write (11, 9100) path(1, i, 1), path(1, i, 2),
     +                      path(1, i, 3)
           write (12, 9100) path(2, i, 1), path(2, i, 2),
     +                      path(2, i, 3)

9100        format (f9.5, 2x, f9.5, 2x, f9.5)

            if (i.eq.1) then
            rout1=sqrt(r2)
            v1=vlu(r2)
            f1=flu(r2)
            fx1=f1*(xx/rout1)
            fy1=f1*(yy/rout1)
            fz1=f1*(zz/rout1)
            end if
            if (i.eq.1+NBEADS/4)then
            rout2=sqrt(r2)
            v2=vlu(r2)
            f2=flu(r2)
            fx2=f2*(xx/rout2)
            fy2=f2*(yy/rout2)
            fz2=f2*(zz/rout2)
            end if
            if (i.eq.1+NBEADS/2) then
            rout3=sqrt(r2)
            v3=vlu(r2)
            f3=flu(r2)
            fx3=f3*(xx/rout3)
            fy3=f3*(yy/rout3)
            fz3=f3*(zz/rout3)
            end if
            if (i.eq.1+3*NBEADS/4)then
            rout4=sqrt(r2)
            v4=vlu(r2)
            f4=flu(r2)
            fx4=f4*(xx/rout4)
            fy4=f4*(yy/rout4)
            fz4=f4*(zz/rout4)
            end if


         end do

         write (18, 1800) rout1, rout2, rout3, rout4
1800     format (4(f15.10, 1x))
         write (19,1900) n, v1*27.2112d0*8065.54d0, 
     +   v2*27.2112d0*8065.54d0,v3*27.2112d0*8065.54d0, v4*27.2112d0*8065.54d0
         write (23, 1900) n, f1, f2, f3, f4
1900     format(i9, 4(2x, 1pe15.8))
         write (20, 1800) fx1, fx2, fx3, fx4
         write (21, 1800) fy1, fy2, fy3, fy4
         write (22, 1800) fz1, fz2, fz3, fz4
         
         call flush(18)

c ------ this pads each snapshot with a blank line at the end.

        write (11, *) ''
        write (12, *) ''

         v=v/dble(NBEADS)
         f=f/dble(NBEADS)

c ------ write the information out to unit 7.

c ------ the numerical factors multiplying v convert it to convenient
c        units, called wavenumbers.

         write (7, 7001) n, v*27.2112d0*8065.54d0, sqrt(rmin), 
     +                   sqrt(rmax)
7001     format (i9, 3(2x, 1pe15.8))

c ------ compute the average x position for the left atom.  (this
c        should remain constant during the simulation.)

         x1=0.0d0
         x2=0.0d0
         x3=0.0d0
         x4=0.0d0

         do i=1, NBEADS
            x1=x1+path(1, i, 1)
         end do

         x1=x1/dble(NBEADS)
        
         do i=1, NBEADS
         write (31, 3100) path(1, i, 1)
         write (32, 3100) path(1, i, 2)
         write (33, 3100) path(1, i, 3)
         write (41, 3100) path(2, i, 1)
         write (42, 3100) path(2, i, 2)
         write (43, 3100) path(2, i, 3)
 3100            format (f10.5)
         end do 
c ------ compute the mean squared, cubic, and quartic deviation
c        of the left atoms beads along the x axis.

c ------ the squared and quartic deviation help us decide whether the
c        distribution can be modeled by a Gaussian.

c ------ the cubic deviation would be sensitive to any asymmetry in the
c        distribution.

         do i=1, NBEADS
            x2=x2+(path(1, i, 1)-x1)**2
            x3=x3+(path(1, i, 1)-x1)**3
            x4=x4+(path(1, i, 1)-x1)**4
         end do

         x2=x2/dble(NBEADS)
         x3=x3/dble(NBEADS)
         x4=x4/dble(NBEADS)

c ------ write the information out to unit 8.

         write (8, 8001) n, x2, x3, x4, 'X'
8001     format (i9, 3(2x, 1pe15.8), 2x, a1)

c ------ repeat for the y direction.

         x1=0.0d0
         x2=0.0d0
         x3=0.0d0
         x4=0.0d0

         do i=1, NBEADS
            x1=x1+path(1, i, 2)
         end do

         x1=x1/dble(NBEADS)

         do i=1, NBEADS
            x2=x2+(path(1, i, 2)-x1)**2
            x3=x3+(path(1, i, 2)-x1)**3
            x4=x4+(path(1, i, 2)-x1)**4
         end do

         x2=x2/dble(NBEADS)
         x3=x3/dble(NBEADS)
         x4=x4/dble(NBEADS)

         write (8, 8001) n, x2, x3, x4, 'Y'

c ------ repeat for the z direction.

         x1=0.0d0
         x2=0.0d0
         x3=0.0d0
         x4=0.0d0

         do i=1, NBEADS
            x1=x1+path(1, i, 3)
         end do

         x1=x1/dble(NBEADS)

         do i=1, NBEADS
            x2=x2+(path(1, i, 3)-x1)**2
            x3=x3+(path(1, i, 3)-x1)**3
            x4=x4+(path(1, i, 3)-x1)**4
         end do

         x2=x2/dble(NBEADS)
         x3=x3/dble(NBEADS)
         x4=x4/dble(NBEADS)

         write (8, 8001) n, x2, x3, x4, 'Z'

      end do

c --- compute the overall fraction of moves that are rejected.

      trej=trej/(dble(NBLOCK)*dble(NBEADS)*dble(2*niter))
 
      write (6, *) 'rejection fraction = ', trej
      write (50,*) trej
c --- end of the entire simulation.

      stop
      end

c --- this subroutine computes the quantity V(Q) that appears in the
c     probability function.

      subroutine prob(force, nn, ds1f, ds1b, ds2f, ds2b, dr1, dr2, p)

      implicit real*8 (a-h, o-z)

      p=0.5d0*force*(ds1f+ds1b+ds2f+ds2b)

      p=p+(vlu(dr1)+vlu(dr2))/dble(nn)

      return
      end

c --- this is the He-He potential, evaluated from the lookup table.

      double precision function vlu(r2)

      implicit real*8 (a-h, o-z)

      parameter (NVBINS=40000)

      common /potcom/ v(2, NVBINS)

      if (r2.gt.9.0d0) then
 
c ------ if the distance is greater than 3 bohr, use the lookup table.

         ibin=idint(20.0d0*r2-180.0d0)+1

         vlu=v(1, ibin)+v(2, ibin)*r2

         return

      else

c ------ otherwise use the actual function.

         vlu=hfdbhe(r2)

         return

      end if

      end

c --- this is the He-He force, evaluateed from the lookup table.

      double precision function flu(r2)
      
      implicit real*8 (a-h, o-z)

      parameter (NVBINS=40000)

      common /forcom/ f(2, NVBINS)

      if (r2.gt.9.0d0) then

c --- if the distance is greater tan 3 bohr, use the lookup table.

         ibin=idint(20.0d0*r2-180.0d0)+1
      
         flu=f(1, ibin)+f(2, ibin)*r2

         return
  
      else

c --- otherwise use the actual function.

          flu=hfdbheforce(r2)

          return

       end if

       end

c --- this subroutine carries out a Monte Carlo move.

      subroutine move(id, path, nrej,NBEADS)

      implicit real*8 (a-h, o-z)

c --- id = 1 to move left atom; id = 2 to move right atom.

      dimension path(2, 10000, 3)

      common /rancom/ rscale, rstate(6)

      common /params/ force, beta, zmass, delta

c --- first choose two beads to jigger.

120   call rstep(rstate, z, rscale)

      k1=idint(z*NBEADS)+1

      call rstep(rstate, z, rscale)
 
      k2=idint(z*NBEADS)+1

c --- make sure that we chose two different beads.

      if (k1.eq.k2) goto 120

c --- put them in order so that k1 < k2.

      if (k1.gt.k2) then
         kk1=k1
         k1=k2
         k2=kk1
      end if

c --- calculate the indices of adjacent beads: f = forward, b = backward

      k1b=k1-1
      k1f=k1+1
      k2b=k2-1
      k2f=k2+1

c --- account for "wrap around".

      if (k1.eq.1) k1b=NBEADS
      if (k1.eq.NBEADS) k1f=1

      if (k2.eq.1) k2b=NBEADS
      if (k2.eq.NBEADS) k2f=1

c --- calculate the harmonic link lengths (squared).

      ds1b=(path(id, k1, 1)-path(id, k1b, 1))**2+
     +     (path(id, k1, 2)-path(id, k1b, 2))**2+
     +     (path(id, k1, 3)-path(id, k1b, 3))**2

      ds1f=(path(id, k1, 1)-path(id, k1f, 1))**2+
     +     (path(id, k1, 2)-path(id, k1f, 2))**2+
     +     (path(id, k1, 3)-path(id, k1f, 3))**2

      ds2b=(path(id, k2, 1)-path(id, k2b, 1))**2+
     +     (path(id, k2, 2)-path(id, k2b, 2))**2+
     +     (path(id, k2, 3)-path(id, k2b, 3))**2

      ds2f=(path(id, k2, 1)-path(id, k2f, 1))**2+
     +     (path(id, k2, 2)-path(id, k2f, 2))**2+
     +     (path(id, k2, 3)-path(id, k2f, 3))**2

c --- check for the case in which k1 and k2 are only different by one;
c     in this case, there are only three "links" to deal with instead
c     of four.

      if (k2.eq.k1+1) then

         ds2b=0.0d0

      end if

      if (k1.eq.1.and.k2.eq.NBEADS) then

         ds1b=0.0d0

      end if

c --- calculate the distance between the left and right atoms.

      dr1=(path(id, k1, 1)-path(3-id, k1, 1))**2+
     +    (path(id, k1, 2)-path(3-id, k1, 2))**2+
     +    (path(id, k1, 3)-path(3-id, k1, 3))**2

      dr2=(path(id, k2, 1)-path(3-id, k2, 1))**2+
     +    (path(id, k2, 2)-path(3-id, k2, 2))**2+
     +    (path(id, k2, 3)-path(3-id, k2, 3))**2

c --- calculate the probability for this configuration.

      call prob(force, NBEADS, ds1f, ds1b, ds2f, ds2b, dr1, dr2, prob1)

c --- choose a random displacement for jiggering the beads.

      call rstep(rstate, dx, rscale)
      call rstep(rstate, dy, rscale)
      call rstep(rstate, dz, rscale)

      dx=(dx-0.5d0)*2.0d0*delta
      dy=(dy-0.5d0)*2.0d0*delta
      dz=(dz-0.5d0)*2.0d0*delta

c --- save the old positions, in case we reject the move.

      x1o=path(id, k1, 1)
      y1o=path(id, k1, 2)
      z1o=path(id, k1, 3)

      x2o=path(id, k2, 1)
      y2o=path(id, k2, 2)
      z2o=path(id, k2, 3)

c --- jigger the beads, keeping the centroid fixed.

      path(id, k1, 1)=path(id, k1, 1)+dx
      path(id, k1, 2)=path(id, k1, 2)+dy
      path(id, k1, 3)=path(id, k1, 3)+dz

      path(id, k2, 1)=path(id, k2, 1)-dx
      path(id, k2, 2)=path(id, k2, 2)-dy
      path(id, k2, 3)=path(id, k2, 3)-dz

c --- calculate the new harmonic link lengths (squared).

      ds1b=(path(id, k1, 1)-path(id, k1b, 1))**2+
     +     (path(id, k1, 2)-path(id, k1b, 2))**2+
     +     (path(id, k1, 3)-path(id, k1b, 3))**2

      ds1f=(path(id, k1, 1)-path(id, k1f, 1))**2+
     +     (path(id, k1, 2)-path(id, k1f, 2))**2+
     +     (path(id, k1, 3)-path(id, k1f, 3))**2

      ds2b=(path(id, k2, 1)-path(id, k2b, 1))**2+
     +     (path(id, k2, 2)-path(id, k2b, 2))**2+
     +     (path(id, k2, 3)-path(id, k2b, 3))**2

      ds2f=(path(id, k2, 1)-path(id, k2f, 1))**2+
     +     (path(id, k2, 2)-path(id, k2f, 2))**2+
     +     (path(id, k2, 3)-path(id, k2f, 3))**2

c --- check for the case in which k1 and k2 are only different by one;
c     in this case, there are only three "links" to deal with instead
c     of four.

      if (k2.eq.k1+1) then

         ds2b=0.0d0

      end if

      if (k1.eq.1.and.k2.eq.NBEADS) then

         ds1b=0.0d0
         
      end if

c --- calculate the new distance between the left and right atoms.

      dr1=(path(id, k1, 1)-path(3-id, k1, 1))**2+
     +    (path(id, k1, 2)-path(3-id, k1, 2))**2+
     +    (path(id, k1, 3)-path(3-id, k1, 3))**2

      dr2=(path(id, k2, 1)-path(3-id, k2, 1))**2+
     +    (path(id, k2, 2)-path(3-id, k2, 2))**2+
     +    (path(id, k2, 3)-path(3-id, k2, 3))**2

c --- calculate the probability for this new configuration.

      call prob(force, NBEADS, ds1f, ds1b, ds2f, ds2b, dr1, dr2, prob2)

c --- check the probabilities to decide whether to accept the move.

c --- if prob2 <= prob1 then we will accept the move.  this is because
c     the real probability is given by exp(-beta*prob), so lower values
c     of prob are really higher values of the actual probability.

c --- so we only need to do additional work if prob2 > prob1.

      if (prob2.gt.prob1) then

         test=exp(-beta*(prob2-prob1))

         call rstep(rstate, z, rscale)

         if (z.gt.test) then

c --------- here we reject the move and reset the beads.

            nrej=nrej+1

            path(id, k1, 1)=x1o
            path(id, k1, 2)=y1o
            path(id, k1, 3)=z1o
   
            path(id, k2, 1)=x2o
            path(id, k2, 2)=y2o
            path(id, k2, 3)=z2o

         end if

      end if

      return
      end
