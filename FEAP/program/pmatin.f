c$Id:$
      subroutine pmatin(tx,d,ul,xl,tl,s,p,idl,ie,iedof,lie,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1.  Modify the write 2006 record to avoid too much  19/06/2007
c           information in the record
c       2.  Increase number of allowable dofs to 30         05/07/2007
c       3.  Add 'plm' for element multiplier partition      20/07/2007
c       4.  Correct format 2006 -- remove 'a1/2x'           20/08/2007
c       5.  Add global equation values for element          27/03/2008
c       6.  Automatically set 'geqnum' and 'nadd'           23/04/2008
c           Add include for 'pglob1.h', allocate for np(258)
c       7.  Check that geqnum > 0 to allocate array         02/05/2008
c       8.  Correct setting of idl array for ndf > 13       09/04/2009
c       9.  Modify type inputs to accept 2 words for types  14/01/2010
c           of standard feap elements (e.g., couple thmech)
c      10.  Cast 'nh1' & 'nh3' to 'int' to assign to 'ie'   18/06/2011
c      11.  Dimension uelnum(1)                             01/05/2011
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Data input routine for material parameters

c      Inputs:
c         tx          - Option identifier
c         prt         - Flag, print input data if true
c         prth        - Flag, print title/header if true

c      Scratch:
c         ul(*)       - Local nodal solution vector
c         xl(*)       - Local nodal coordinate vector
c         tl(*)       - Local nodal temperature vector
c         s(*)        - Element array
c         p(*)        - Element vector
c         idl(*)      - Local degree of freedom integer data
c         lie(ndf,*)  - Local dof assignment array

c      Outputs:
c         ie(nie,*)   - Element descriptor data
c         iedof(*,*,*) - Element nodal assembly data
c         d(ndd,*)    - Material parameters generated by elements
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'chdata.h'
      include  'comfil.h'
      include  'eldata.h'
      include  'erotas.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'iodata.h'
      include  'iosave.h'
      include  'pdata3.h'
      include  'pglob1.h'
      include  'pointer.h'
      include  'rigid1.h'
      include  'rigid2.h'
      include  'strnum.h'
      include  'sdata.h'

      logical   pcomp,pinput,tinput,vinput,errck,prt,prth,doflg
      logical   setval,palloc
      character tx(2)*15,mtype*69,etype*20
      integer   i,j, ii,il,is, geqold
      integer   lie(ndf,*), ie(nie,*),iedof(ndf,nen,*),idl(*)
      real*8    d(ndd,*),ul(*),xl(*),tl(*),s(*),p(*), td(50)
      real*8    uelnum(1)

      save

c     Data input for material set ma

      errck = vinput(yyy(16:30),15,td,1)
      ma     = nint(td(1))
      geqold = geqnum
      if(nummat.eq.1) ma = max(1,ma)
      if(ma.le.0 .or. ma.gt.nummat) then
        if(ior.gt.0) then
          write(iow,3000) ma
          write(iow,3001)
          write(ilg,3000) ma
          write(ilg,3001)
          call plstop()
        else
          write(*,3000) ma,' Reinput mate,ma'
        endif
        return
      endif

      if(prt) then
        call prtitl(prth)
        write(iow,2000)
      endif

c     Set material identifier

      do j = 3,80
        if(xxx(j:j).eq.' ' .or. xxx(j:j).eq.',') then
          do i = j+1,80
            if(xxx(i:i).eq.' ' .or. xxx(i:i).eq.',') go to 300
          end do ! i
        endif
      end do ! j
      i = 69
300   mtype = xxx(i+1:70)

c     Input record for element type selection

301   if(ior.lt.0) then
        write(*,2004)
        call pprint('   >')
      endif
      errck = tinput(tx,2,td,min(ndf+2,14))
      if(errck) go to 301

c     Set flag to save material properties

      inquire(file=fmtl, opened = errck)
      if(errck) then
        lmate = .true.
      else
        lmate = .false.
      endif

c     Set material type for standard and user elements

      if(.not.pcomp(tx(1),'user',4)) then
        call pelnum(tx,iel,errck)
      else
        errck = vinput(tx(2),15,uelnum,1)
        iel   = nint(uelnum(1))
      endif

      if(ie(nie-1,ma).ne.0 .and. iel.eq.0) then
        iel = ie(nie-1,ma)
      else
        ie(nie-2,ma) = nint(td(1))
        if(ie(nie-2,ma).le.0) ie(nie-2,ma) = ma
      endif

c     Set print head

      if(ma.eq.ie(nie-2,ma)) then
        etype = ' '
      else
        write(etype,'(a16,i4)') 'Element Material',ie(nie-2,ma)
      endif

c     Set idl for first group of dof's

      do j = 1,min(ndf,12)
        idl(j) = nint(td(j+1))
      end do ! j

c     For large number of dof's input additional records and set idl

      if(ndf.gt.12) then
        il = 12
        do ii = 1,(ndf+2)/16
          is = il+1
          il = min(is+15,ndf)
302       errck = pinput(td,16)
          if(errck) go to 302
          do j = is,il
            idl(j) = nint(td(j-is+1))
          end do ! j
        end do ! ii
      endif

c     Check to see if degree of freedoms to be reassigned

      do i = 1,ndf
        if(idl(i).ne.0) go to 303
      end do ! i

c     Reset all zero inputs

      do i = 1,ndf
        idl(i) = i
      end do ! i

303   ie(nie-1,ma) = iel

c     Set flags for number of history and stress terms

      mct  = 0
      nh1  = 0
      nh2  = 0
      nh3  = 0
      nlm  = 0
      plm  = 0
      nge  = 0
      istv = 0
      istp = 0
      istc = 1

c     Output information

      if(prt) then
        if(iel.gt.0) then
          write(iow,2001) ma,iel,etype,(j,idl(j),j=1,ndf)
        else
          write(iow,2002) ma,tx(1)(1:5),etype,(j,idl(j),j=1,ndf)
        endif
        write(iow,2003) mtype
      else
        if(iel.gt.0) then
          write(iow,2001) ma,iel,etype
        else
          write(iow,2002) ma,tx(1)(1:5),etype
        endif
        write(iow,2003) mtype
      endif

c     Save material sets to file for some later uses

      if(lmate) then
        write(iwd,2005) ma
        if(iel.gt.0) then
          write(iwd,2006) ' user ',iel,(nint(td(j)),j=1,ndf+1)
        else
          if(pcomp(tx(2),'     ',5)) then
            tx(2) = 'elmt '
          endif
          write(iwd,2007) tx(1)(1:5),tx(2)(1:5),(nint(td(j)),j=1,ndf+1)
        endif
      endif

c     Obtain inputs from element routine

      do j = 1,nen+1
        do i = 1,ndf
          lie(i,j) = i
        end do ! i
      end do ! j
      rotyp = 0

c     Set default element plot type to mesh dimension

      pstyp = ndm

c     Get material input data element library (isw = 1)

      call elmlib(d(1,ma),ul,xl,lie,tl,s,p,ndf,ndm,nst,iel,1)

c     Check if dof values set by each node

      doflg = .false.    ! Check if mixed conditions exist
      do i = 1,ndf
        do j = 1,nen
          if(lie(i,j+1).eq.0) then
            doflg = .true.
          endif
        end do ! j
      end do ! i

c     Set assembly information

      if(doflg) then
        do i = 1,ndf
          do j = 1,nen
            if(lie(i,j+1).gt.0) then
              iedof(i,j,ma) = idl(i)
            else
              iedof(i,j,ma) = 0
            endif
          end do ! j
        end do ! i
      else
        do i = 1,ndf
          if(lie(i,1).gt.0) then
            do j = 1,nen
              iedof(i,j,ma) = idl(i)
            end do ! j
          else
            do j = 1,nen
              iedof(i,j,ma) = 0
            end do ! j
          endif
        end do ! i
      endif

c     Set element plot shape type (can be reset by element library)

      ie(1,ma) = pstyp

c     Set element rigid body option

      if(nrmatl.gt.0) then
        ie(2,ma)       = nrmatl
        nrbody         = max(nrbody,nrmatl)
        rbtype(nrmatl) = 0
        rbcen(nrmatl)  = 0
      endif

c     Set number of history terms

      if(nh1.ne.0) then
        ie(nie,ma) = int(nh1)
      else
        ie(nie,ma) = mct
      endif
      ie(nie-5,ma) = int(nh3)
      if(nlm.gt.0) then
        ie(nie-8,ma) = nlm
        nadd         = max(nadd,nlm)
      endif
      if(plm.gt.0) then
        ie(nie-9,ma) = plm
      endif
      if(nge.gt.0) then
        ie(nie-10,ma) = nge
        nadd          = max(nadd,nge)
        geqnum        = max(geqnum,nge)
        if(np(258).eq.0 .and. geqnum.gt.0) then
          setval = palloc( 258,'GUVAL',geqnum, 2)
        endif
      endif

c     Set maximum number of element plot variables

      npstr = max(npstr,istv)

c     Set rotational update type

      if(rotyp.ne.0) then
        ie(nie-6,ma) = rotyp
      endif
      lmate = .false.

c     Allocate storage for global equations

      if(np(258).eq.0 .and. (geqnum.gt.0 .or. geqnum.gt.geqold)) then
        setval = palloc( 258,'GUVAL', geqnum, 2)
      endif

c     Formats

2000  format('   M a t e r i a l    P r o p e r t i e s')

2001  format(/5x,'Material  Number',i4,' for User Element Type',i3,5x,/
     &        5x,a,/:/
     &        5x,'Degree of Freedom Assignments    Local    Global'/
     &       37x,'Number',4x,'Number'/(31x,2i10))

2002  format(/5x,'Material  Number',i4,' for Element Type: ',a,5x,/
     &        5x,a,/:/
     &        5x,'Degree of Freedom Assignments    Local    Global'/
     &       37x,'Number',4x,'Number'/(31x,2i10))

2003  format(5x,a)

2004  format(' Input: Elmt type, Id No., dof set')

2005  format('MATErial',i8)
2006  format(a5,15i4:/(16i4))
2007  format(a5,1x,a5,1x,14i4:/(16i4))

3000  format(' *ERROR* PMATIN: Illegal material number: ma=',i5:,a)
3001  format('                 Check value of NUMMAT on control record')

      end
