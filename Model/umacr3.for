      subroutine umacr3(lct,ctl)


c__________________________________________________TCG. 12.05.2003___71
c
c     Advance crackfront 
c
c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters
c         prt       - Flag, output if true
c
c        common blocks 
c        /tets/      Latent tets block
c        /tip/       Crack tip data block
c
c
c    Actions   
c            (i)   update History vector with crack data
c            (ii)  Advance crack front
c
c
c
c  Declare variable types
c     
c     IXpoint     Poiter of IX-Array 
c     IXlen       Length of IX-Array 
c     IXpre       Precision of IX-Array 
c     IXflag      Flag of IX-Array 
c     Xpoint      Poiter of X-Array (Ref. Nodal coordinates)
c     Xlen        Length of X-Array 
c     Xpre        Precision of X-Array 
c     Xflag       Flag of X-Array 
c     Hpoint      Poiter of H-Array 
c     Hlen        Length of H-Array 
c     Hpre        Precision of H-Array 
c     Hflag       Flag of H-Array 
c     ns          No. of shared faces
c     ng          global node no. of tet
c     elmtnr      Element no.
c     ns          No. of shared faces
c     sflag       Shareflag (lineNo. in cTip, which shares with tet)
c                 advancing location of the crack tip
c     flag        .false. no shared faces ; .true. shared faces
c     Ae          Area of discontinuity in ref config
c     Ve          Volume in ref. config
c     rnode_flag  Real value representing the node_flag
c     n1=n2     Length of n1 and n2 History
c     n3         Length of n3 History
c     iel        element no. stored in common block eldata.h
c
c     /tets/      Latent tets block
c        nTets    No. of lines in cTets
c        cTets    Data as following typical line
c 
c           Elmt    node1  node2  node3  node4  id1  id2  id3  id4
c
c
c     /tip/      Crack tip data block
c        ncTip    No. of lines in cTip
c        cTip    Data as following typical line
c 
c           node1  node2  node3  id1  id2  id3 
c
c
c     /node_id/   Node id data block
c        nid    node id 
c 
c           nid = +1 .... node is in Omega plus
c           nid =  0  ... node is in Omega minus or does not belong to crossed element
c____________________________________________________________________71


      implicit  none

      include  'umac1.h'
	include  'comblk.h'
	include  'pointer.h'
	include  'cdata.h'   
	include  'sdata.h' 
      include  'eldata.h'
      include  'hdata.h'
      include  'hdatam.h'
      include  'iofile.h'

      include  'cTets.h'
      include  'cTip.h'
      include  'cGeom.h'
      include  'nodes_id.h'
      include  'nLength.h'
      include  'Visco.h'

	integer pN0,pAe,pVe,pLcn,pVepVe,pXbar,pcini,pnodeflag,pcn
	integer pconp, pconm, pXl2d, psig, pfi, pNelem, pT0, PNold
	integer  i,j,k, elmtnr, ns, idtet(4)

	logical abbr

      logical   pcomp, prt
      character lct*15
      logical Xflag
	integer Xpoint, Xlen, Xpre
      logical IXflag
	integer IXpoint, IXlen, IXpre
      logical dflag
	integer dpoint, dlen, dpre
      logical hflag
	integer hpoint, hlen, hpre, nn
      real*8  ctl(3)
	integer sflag(4)
	logical flag
	integer ntets_old, fn_unit, restart_unit
	real*8 N0(3), fp(3)
	integer n1, n2, n3      
      integer HistpElmt, HLocpElmt

      fn_unit=98
      restart_unit = 99

      if(pcomp(uct,'mac3',4)) then      ! Usual    form
	  uct='GROWTH'
	  iniflag=.false.
      elseif(urest.eq.1) then           ! Read  restart data
	  write(*,*)
	  write(*,*)'  Read from file: user_restart'

        OPEN(UNIT=restart_unit,FILE='user_restart')
	  write(*,*)'         Crack tip data'
c iniflag
	  read(restart_unit,3000) iniflag

c inimodel
	  read(restart_unit,3001) imo
       
c ctip
	  read(restart_unit,3001) nctip
	   do i= 1,nctip
	     read(restart_unit,3002) (cTip_point(i,j),j=1,3)
	   enddo
	   do i= 1,nctip
	     read(restart_unit,3003) (cTip(i,j),j=1,6)
	   enddo
	   do i= 1,nctip
	     read(restart_unit,3001) mTi(i)
	   enddo

c cmat
	  read(restart_unit,3001) ncmat
	  do i= 1, ncmat
	     read(restart_unit,3001) cmat(i)
	  enddo
c nid
	  write(*,*)'         Node id list'
        do i=1,numnp
	    read(restart_unit,3004) nid(i)	       
	  enddo
c cGeom 
	  write(*,*)'         Crack geometry'
        read(restart_unit,3001) nsurf
        do i=1,nsurf
	    read(restart_unit,3006) ncorn(i),matnc(i),
	1                            (cnode(i,j),j=1,12)        
	  enddo
	  write(*,*)
	  Close(restart_unit)

      elseif(urest.eq.2) then           ! Write restart data

	  write(*,*)
	  write(*,*)'  Save to file: user_restart'

        OPEN(UNIT=restart_unit,FILE='user_restart')
	  write(*,*)'         Crack tip data'
c iniflag
	  write(restart_unit,3000) iniflag
c inimodel
	  write(restart_unit,3001) imo
c ctip
	  write(restart_unit,3001) nctip
	   do i= 1,nctip
	     write(restart_unit,3002) (cTip_point(i,j),j=1,3)
	   enddo
	   do i= 1,nctip
	     write(restart_unit,3003) (cTip(i,j),j=1,6)
	   enddo
	   do i= 1,nctip
	     write(restart_unit,3001) mTi(i)
	   enddo
c cmat
	  write(restart_unit,3001) ncmat
	  do i= 1, ncmat
	     write(restart_unit,3001) cmat(i)
	  enddo

c nid
	  write(*,*)'         Node id list'
         do i=1,numnp
	     write(restart_unit,3004) nid(i)	       
	   enddo
c cGeom 
	  write(*,*)'         Crack geometry'
         write(restart_unit,3001) nsurf
         do i=1,nsurf
	     write(restart_unit,3006) ncorn(i), matnc(i),
     1                             (cnode(i,j),j=1,12)        
	   enddo
	  Close(restart_unit)
	  write(*,*)

      else                              ! Perform user operation
	 


c set elementhistory pointer for elem11
	 pAe       = 0
	 pVe       = 1
	 pVepVe    = 2   
       pN0       = 3
	 pcini     = 6  
	 pcn       = 7   
	 pconp     = 8   
	 pconm     = 9   
	 pLcn      = 10   
	 pXbar     = 26  
	 pXL2d     = 30  
	 pSig      = 38  
	 pFi       = 44  
	 pNelem    = 53  
	 pT0       = 355  
	 pNold     = 356
       !pT0n     = 360  !CJM
       !pT0t     = 361  !CJM
	 pnodeflag = 0   ! not used for elmt11
      
 
 

c read history pointer from hdatam.h
       n1 = nhmax
	 n2 = nhmax
	 n3 = nh3max

c      write(*,*) 'n1, n2, n3', n1,n2,n3
c      write(iow,*) 'n1, n2, n3', n1,n2,n3
c pointer of n3 array
       nn = 2*n1 + n3




c  get pointer and length of the History     
	 call pgetd ('H  ',hpoint,hlen,hpre,hflag)
c  get pointer and length of the X     
	 call pgetd ('X ',xpoint,xlen,xpre,xflag)
c  get pointer and length of the IX     
	 call pgetd ('IX ',ixpoint,ixlen,ixpre,ixflag)
c  get pointer and length of the D     
	 call pgetd ('D ',dpoint,dlen,dpre,dflag)

c initialize cTip
	 if (not(iniflag)) then
c initialize cTip
        call cTip_ini (ndm,numel,nen,nen1,hr(xpoint),xlen,
	1                  mr(ixpoint),ixlen)
c initialize node id
        do i = 1,numnp
	    nid(i) = 0
        enddo
c initialize fn file
        OPEN(UNIT=fn_unit,FILE='fn')
	  write (fn_unit,*) 'BOUNdary'
	  write(fn_unit,*)
	  write(fn_unit,*) 'END'
	  Close(fn_unit)

c compute neighbouring elements
        if (chl.ne.0.0d0) then
	    write(*,*) 'Compute neighbour elements......'
          call neighbour_11(hr(dpoint),dlen,n1,n3,nen1,pXbar,pNelem,
	1                  pVe,hr(hpoint),hlen,mr(ixpoint),ixlen)
        endif  
	  iniflag=.true.
c initialize cGeom
        nsurf=0
	  call plot_ctip11

        goto 1000
	 endif
      
c  check crack-initialization criterion
c  compute non-local stress and deformation
       call asseble_ctets_11(hr(dpoint),dlen,n1,n3,nen1,pCini,pN0,
	1                     pVe,pXbar,pSig,pFi,pNelem,pT0,pNold,
     2                     hr(hpoint),hlen,mr(ixpoint),ixlen)


       if (ntets.ne.0) then

c remove duplicate lines in ctet if necessary
 	  do j= 1,ntets
         do i =j+1,ntets
	     if (cTets(i,1).eq.cTets(j,1)) then
	       cTets(j,1) = 0
	     endif
	   enddo
	  enddo

c rebuild cTtets (remove zero line and count active tets)
       ntets_old = nTets  
	 nTets = 0
	 do i= 1,ntets_old
	  if(cTets(i,1).ne.0) then
	   nTets = nTets + 1
	   do k= 1,9
	    cTets(nTets,k) = cTets(i,k)
	   enddo
	  endif
	 enddo       

	call plot_ctets11

c check if tets share a crackTip facet and compute node id
1      abbr=.true.  
       do i = 1, ntets
        
c retrieve global element number
	  elmtnr = cTets(i,1)    !necessary for cini and nid of quadratic elem

c get number of element nodes
        call get_nodes_elem(elmtnr,nen,nen1,mr(ixpoint),ixlen, nel)

        ns=0
	  sflag=0
        do j= 1,ncTip
c compare current tet and triangle if there have shared faces
c and predefine tet id 
	    call compare10(i,j,flag)
	    if (flag) then
	      abbr=.false.
	      ns = ns+1
	      sflag(ns) = j
	    endif
	  enddo
	  if (ns.gt.4) then
	    write(*,*) 'CRITICAL WARNING ** ns too large in umacr3'
	    write(iow,*) 'CRITICAL WARNING ** ns too large in umacr3'
	  endif
c compared the current tet with all triangles

c introduce crack in tets and advance crack front
       if (sflag(1).NE.0) then ! tet is at the front
	  if(xflag.and.ixflag.and.hflag) then
	      call advance_11(i,ns,sflag, 
	1        ndm,nen,numel,nen1,hr(xpoint),xlen,mr(ixpoint),ixlen,
     2	    hr(hpoint),hlen,n1,n3,idtet,N0,fp,
     3        pN0,pAe,pVe,pLcn,pVepVe,pXbar,pcini,pnodeflag,pcn,
     4        pconp,pconm,pXl2d,PNold,pNelem)
            
c re-allocate history varibales to the subelemnet's Gausspoints
            HistpElmt= n1+n1+n3        ! Length of history vector per element
            HLocpElmt= (elmtnr-1)*HistpElmt   !begin of the element history in the history vector
	      call reallo_hist_variables(hr(hpoint+ HLocpElmt),HistpElmt,n1,n3)
        else
	   write (*,*) 'cannot read common block!!'
	  endif
c include enhancement information (stored  in nid)

c	  call plot_ctip11

c complete enhancement information
c for corner nodes
        do j = 2, 5
	    nid(cTets(i,j)) = cTets(i,j+4)
	  enddo
c for mideside nodes

        if (nel.eq.10) then
          call id_mideside_11(elmtnr, ndm,nen,nen1,hr(xpoint),xlen,
	1                mr(ixpoint),ixlen, N0,fp)
	  endif


c remove line form cTets
        cTets(i,1) = 0
       endif
	 enddo

c rebuild cTtets (remove zero line and count active tets)
       ntets_old = nTets  
	 nTets = 0
	 do i= 1,ntets_old
	  if(cTets(i,1).ne.0) then
	   nTets = nTets + 1
	   do k= 1,9
	    cTets(nTets,k) = cTets(i,k)
	   enddo
	  endif
	 enddo
	
	call plot_ctets11
	call plot_ctip11




c do until no further latent tet is on the tip
       if (not(abbr)) goto 1

c write some data for information only
       write(*,5000)  nsurf,ntets,ncTip
       write(iow,5000)  nsurf,ntets,ncTip

c initialize latent tet data
      cTets = 0
	ntets = 0
      
c modify nid in view of closed crack at the tip 
      call close_front_11(mr(ixpoint),ixlen,hr(hpoint),hlen,
	1                    n1,n3,pcini,nen1)

c write fn file
       OPEN(UNIT=fn_unit,FILE='fn')
	 write (fn_unit,*) 'BOUNdary'
	 do i = 1,numnp
	  if (nid(i).ne.0)  write(fn_unit,1001) i 
     !!   if (nid(i).ne.0)  write(fn_unit,1002) i, 
     !!1  boun_og(i,1), boun_og(i,2), boun_og(i,3) 
       enddo
	 write(fn_unit,*)
	 write(fn_unit,*) 'END'
	 Close(fn_unit)

      endif

1001  format(i6, ',,0,0,0,0,0,0')
1002  format(i6, ',,',i1,',',i1,',',i1,',0,0,0')      

2001  format('** New boundary conditions of node ', i4, ' **' )

3000  format(5x, L1/)
3001  format(5x, i6/)
3002  format(5x, e12.5, e12.5, e12.5/)
3003  format(5x, i6, i6, i6, i6, i6, i6/)
3004  format(5x, i6/)
3005  format(5x, i6, i6, i6, i6, i6, i6, i6, i6, i6/)
3006  format(5x, i6, i6, e12.5, e12.5, e12.5, e12.5, e12.5, 
     1               e12.5, e12.5, e12.5, e12.5,e12.5, e12.5, e12.5/)
4000  format('  *ERROR* Element of unknown order')

5000  format(
     1 5x, '------------------------------------------------------'/
     2  10x,'Number of discontinuities:             ',i8/
     3  10x,'Number of tets meeting failure crit.:  ',i8/
     4  10x,'Number of facets at crack tip:         ',i8/
     5  10x,'Numbers must not exceed dimension of associated'/
     6  10x,'arrays defined in cGeom.h, cTets.h and cTip.h  '/
     7 5x, '------------------------------------------------------')
     
	endif

	
      

1000  end

	subroutine compare10(tetln,fln,flag)
c__________________________________________________TCG. 24.11.2002___71
c
c     compare if tetraeder shares tria 
c
c
c  Declare variable types
c     
c     tetln   line in ctets
c     fln     line in ctip
c     elemtnr element no
c     tet     Global node No. of tetraeder
c     tria    Global node No. of triangle
c     idtria  Id of triangle nodes (+1 or -1)
c     idtet   Id of tetraeder nodes (+1 , -1 or 0)
c     flag    .true. tria is part of tet
c             .false. tria is not part of tet
c
c____________________________________________________________________71
c
      implicit none

	logical flag

	integer tet(4), tria(3), cont, i, j, idtria(3)
	integer tetln, fln

      include  'cTets.h'
      include  'cTip.h'



c retrieve global nodes of current element 
      do j= 1,4
	  tet(j) =  cTets(tetln,1+j)
	enddo
c retrieve  global nodes and id of current triangle
      do j= 1,3
	  tria(j) =  cTip(fln,j)
	  idtria(j) =  cTip(fln,j+3)
	enddo

      cont = 0

      do i= 1,4
	  do j= 1,3 
	   if (tet(i).eq.tria(j)) then
	     if (abs(idtria(j)).ne.0) then
	       cTets(tetln,5+i) = idtria(j)/abs(idtria(j))
	     else
	       cTets(tetln,5+i) = idtria(j)
	     endif
	     cont = cont +1
	   endif
	  enddo 
	enddo
  
      if (cont.eq.3) then
	  flag = .true.
	else
	  flag =.false.
	endif
   
      end




      subroutine  advance_11(tetln,ns,sflag,ndm,nen,numel,nen1,
	1       x,xlen,ix,ixlen, hist,hlen,n1,n3,idtet,N0,fp,
     2       pN0,pAe,pVe,pLcn,pVepVe,pXbar,pcini,pnodeflag,pcn,
     3	   pconp,pconm,pXl2d,pNold,pNelem)
c__________________________________________________TCG. 19.11.2002___71
c
c  Declare variable types
c
c     tetln       current line in ctets
c     ndm         Number of dimension
c     numel       Number of elements
c     numnp       Number of nodal points
c     sflag       shareflag
c     ntet        nodes of tetraeder
c     idtet       ID of nodes of tetraeder
c     x           Ref. coordinate of nodes
c     xlen        Length of x-vector
c     ix          Element-connectivity
c     ixlen       Length of ix-vector
c     hist        history vector
c     hlen        Length of history vector
c     h1len       Length of n1 (and n2) history vector per elment
c     h3len       Length of n3 history vector per elment
c     ns          No. of shared faces
c     xc          center of elment
c     xc_curr     center of current element
c     nl          noin local measure (e.g. 3 elements)
c     chl         characzeristic length, i.e. nl*V_e^1/3
c     cn          No. of corners of disc.
c____________________________________________________________________71
      implicit none

	integer pN0,pAe,pVe,pLcn,pVepVe,pXbar,pcini,pnodeflag,pcn
	integer pconp, pconm, pXl2d, pNold, pNelem
	
	integer  xlen, hlen, nn, ns, tetln, cn, tn1
	integer  i, k, j, ndm, ixlen, n1,n3, nen, numel
      integer  ix(ixlen), nodenr, elmtnr,  nen1
	integer  idtet(4), idtemp(4), ntet(4), sflag(4)
	integer  ncTip_old
	real*8   x(xlen), hist(hlen), fp(3)
	real*8   N0(3), Ae, Ve, xl(3,4)
	real*8   VepVe, Lcn(4,4), rconp, rconm, xl2d(2,4)
	real*8   xn3d(3,4),dN


      include  'cGeom.h'
      include  'cTets.h'
      include  'cTip.h'
      include  'eldata.h'
      include  'iofile.h'



c pointer of n3 array
      tn1=2*n1
      nn = tn1 + n3

c retrieve global element number
	 elmtnr = cTets(tetln,1)    
c retrieve tet nodes and id
      do i= 1,4
	 ntet(i) = cTets(tetln,i+1)
	 idtet(i) = cTets(tetln,i+5)
	enddo
c compute average frontpoint
      fp =0.0d0
      do i= 1,ns
	  do j= 1,3
	    fp(j) = fp(j) + cTip_point(sflag(i),j)
	  enddo
	enddo
	fp = fp/Dfloat(ns)

c read direction vectors
	do k= 1,3
	  N0(k) = Hist(1+(elmtnr-1)*nn+tn1+pN0-1+k)
	enddo

cc modify normal if it is too far away from previous
cc non-uniqueness of principle stress criterion!
c      call mod_princ_dir_11(ndm, elmtnr, n1,n3,nen1, 
c     1         pN0, pcini, pNelem, pNold,
c     2         pAe, hlen, hist, ixlen, ix, N0)

c retrieve reference coordinates   
      do j = 1,4
	  nodenr = ix((elmtnr-1)*nen1 + j)
	  do k= 1,ndm
	    xl(k,j) = x((nodenr-1)*ndm + k) 
	  enddo
	enddo	   
	
c modify Normal if bad geometry criterion is satisfied
      call mod_N0_11(ndm, nen1, xlen, x, hlen, hist, ixlen, ix, 
	1            elmtnr, ns, sflag, pN0, pcini, nn, tn1, N0)

c predictor step	
	call cut_tetra_11(elmtnr,xl,N0,fp,Ae,Ve,VepVe,Lcn,xl2d,
	1				xn3d,idtemp, cn,rconp,rconm)
      
      if ((cn.lt.3).or.(cn.gt.4)) then
	   write(*,*) 
     1		'CRITICAL WARNING ** illegal discontinuity'
	   write(iow,*) 
     1		'CRITICAL WARNING ** illegal discontinuity'

      endif

c corrector step (smooth cracksurface)
      call smooth_11 (cn, xn3d, fp, Ve,N0,dN)
	 
c      write (*,*) 'Change',dn
c      write (iow,*) 'Change',dn

c modify Normal if bad geometry criterion is satisfied
      call mod_N0_11(ndm, nen1, xlen, x, hlen, hist, ixlen, ix, 
	1           elmtnr, ns, sflag, pN0, pcini, nn, tn1, N0)
c recompute cut
	call cut_tetra_11(elmtnr,xl,N0,fp,Ae,Ve,VepVe,Lcn,xl2d,
	1				xn3d,idtemp, cn,rconp,rconm)
c update idtet
      do i= 1,4
       if (idtet(i).eq.0) idtet(i) = idtemp(i)
	enddo

c write individuell crack data into history 
      Hist(1+(elmtnr-1)*nn+tn1+pAe) = Ae             !discontinuity area
      Hist(1+(elmtnr-1)*nn+tn1+pVe) = Ve             !elementvolume
      Hist(1+(elmtnr-1)*nn+tn1+pVepVe) = VepVe       !volumefraction
      Hist(1+(elmtnr-1)*nn+tn1+pcn) = cn             !No. of corners of disc.
      Hist(1+(elmtnr-1)*nn+tn1+pconp) = rconp        !connectivity in omega+
      Hist(1+(elmtnr-1)*nn+tn1+pconm) = rconm        !connectivity in omega-

	k= pLcn
	do i= 1,4  
	  do j= 1,4                                     !Volume coordinates of
	    Hist(1+(elmtnr-1)*nn+tn1+k)= Lcn(j,i)      !disc. corner nodes
	    k= k+1
	  enddo
	enddo 
	k= pXl2d
	do i= 1,4  
	  do j= 1,2                                     !2d coordinates of
	    Hist(1+(elmtnr-1)*nn+tn1+k)= xl2d(j,i)     !disc. corner nodes
	    k= k+1
	  enddo
      enddo 

c write modified direction vectors
	do k= 1,3
	  Hist(1+(elmtnr-1)*nn+tn1+pN0-1+k) = N0(k)
	enddo


c store xn3d in crack_geom
      nsurf = nsurf+1
	if (nsurf.gt.50000) then
	  write(*,*) 
     1		'CRITICAL WARNING ** nsurf too large in advance_11'
	  write(iow,*) 
     1		'CRITICAL WARNING ** nsurf too large in advance_11'
	endif
c get material no. of current tet
      matnc(nsurf) = int(ix((elmtnr-1)*nen1 + nen1))

	k = 1
	ncorn(nsurf) = cn
      do i= 1,4
	  do j= 1,3
	    cnode(nsurf,k) = xn3d(j,i)
	    k= k+1
	  enddo
      enddo    

c store idetet in cTets
      do i= 1,4
	  cTets(tetln,i+5)=idtet(i)
	enddo     

c advance crack front

      call add_triag10 (ns,sflag,fp,ntet,xl,N0,idtet)

      
c sort cTip (remove zero line and count active facets)
      ncTip_old = ncTip  
	ncTip = 0
	do i= 1,ncTip_old
	  if(cTip(i,1).ne.0) then
	   ncTip = ncTip + 1
	   do k= 1,6
	    cTip(ncTip,k) = cTip(i,k)
	   enddo
	   do k= 1,3
          cTip_point(ncTip,k) = cTip_point(i,k)
	   enddo
	  endif
	enddo

c set cini falg 
      Hist(1+(elmtnr-1)*nn+tn1+pcini) = 1.0d0
	 
 	end




	subroutine add_triag10 (ns,sflag,fp,ntet,xl,N0,idtet)
c__________________________________________________TCG. 19.11.2002___71
c
c  add tet to crack front which shares faces with front
c  effect  ->  (1) remove all shared sides from the front
c          ->  (2) add remaining side to front if id changes
c
c
c  sflag      number of line of cTip with share a tet side
c  ntet       node number of tet
c  idtet      ID of tet nodes (+1,-1)
c  ns         No. of shared faces
c  fp         frontpoint
c____________________________________________________________________71
      implicit none


	integer sflag(4), ntet(4), idtet(4)
	integer trian(4,3)
	integer i,j, ns, nnew
	integer res(4,6)
	real*8 fp(3),xl(3,4)
	real*8 N0(3), fp_new(4,3)

      include  'cTets.h'
      include  'cTip.h'


c retrieve trianglenodes and id
      do i = 1,ns
	  do j= 1,3
	    trian(i,j) = cTip(sflag(i),j)
	  enddo
	enddo

c remove shared faces from cTip
      do j=1,ns
c	  cTip(sflag(j),1) = 0
	enddo

     
	if (ns.ne.4) then 
      
        call compare_triangle (ns,trian,fp,ntet,idtet,xl,N0,
	1                        nnew, res, fp_new)  


cc add new frontfacets
        do i= 1,nnew
	     ncTip = ncTip + 1
	     do j= 1,6
	      cTip(ncTip,j)       = res(i,j)
	     enddo
	     do j= 1,3
	      cTip_point(ncTip,j) = fp_new(i,j)
	     enddo
	  enddo
	endif

      end




      subroutine compare_triangle (n,tria,fp,ntet,idtet,xl,N0,
     1                        ns, res, fp_new)
c__________________________________________________TCG. 19.11.2002___71
c
c
c  compare shared tetsides
c
c  n     No. of triangles to be checked
c  res   sides in sorted order
c  ftet  tetraeder face
c  new_flag old flag --> .false.  tetside is at the old front
c                        .true.   tetside is at the new front
c  fp         frontpoint
c
c____________________________________________________________________71
      implicit none

	include 'tdata.h'
      include 'iofile.h'

      logical new_flag, n1, n2, n3, flag
     
      integer n, ns,ii,jj
	integer tria(4,3), ntet(4), idtet(4), res(4,6)
	integer ftet(4,6), np 
	integer i,j,k,l
	real*8 fp(3), xl(3,4), xtria(4,3,3), m(3), mu, temp
	real*8 N0(3), eps, fp_new(4,3)
	real*8 x(3), y(3), z(3)

	data eps/1.0d-8/

c generate tet faces
      
      ftet(1,1) = ntet(1)
      ftet(1,2) = ntet(2)
      ftet(1,3) = ntet(3)
      ftet(1,4) = idtet(1)
      ftet(1,5) = idtet(2)
      ftet(1,6) = idtet(3)
	do i= 1,3
	  xtria(1,i,1) = xl(i,1)
	  xtria(1,i,2) = xl(i,2)
	  xtria(1,i,3) = xl(i,3)
	enddo


      ftet(2,1) = ntet(1)
      ftet(2,2) = ntet(2)
      ftet(2,3) = ntet(4)
      ftet(2,4) = idtet(1)
      ftet(2,5) = idtet(2)
      ftet(2,6) = idtet(4)
	do i= 1,3
	  xtria(2,i,1) = xl(i,1)
	  xtria(2,i,2) = xl(i,2)
	  xtria(2,i,3) = xl(i,4)
	enddo

      ftet(3,1) = ntet(2)
      ftet(3,2) = ntet(3)
      ftet(3,3) = ntet(4)
      ftet(3,4) = idtet(2)
      ftet(3,5) = idtet(3)
      ftet(3,6) = idtet(4)
	do i= 1,3
	  xtria(3,i,1) = xl(i,2)
	  xtria(3,i,2) = xl(i,3)
	  xtria(3,i,3) = xl(i,4)
	enddo

      ftet(4,1) = ntet(1)
      ftet(4,2) = ntet(3)
      ftet(4,3) = ntet(4)
      ftet(4,4) = idtet(1)
      ftet(4,5) = idtet(3)
      ftet(4,6) = idtet(4)
	do i= 1,3
	  xtria(4,i,1) = xl(i,1)
	  xtria(4,i,2) = xl(i,3)
	  xtria(4,i,3) = xl(i,4)
	enddo

      res=0
      ns=0
	fp_new = 0.0d0

      do j= 1,4
c initialize
       new_flag=.true.
c check differnet id
	 If(Abs(ftet(j,4)+ftet(j,5)+ftet(j,6)).ne.3) then ! ftet is crossed by disc
	   do i= 1,n
	     n1 = (tria(i,1).eq.ftet(j,1)).or.(tria(i,1).eq.ftet(j,2))
	1			 .or.(tria(i,1).eq.ftet(j,3))
	     n2 = (tria(i,2).eq.ftet(j,1)).or.(tria(i,2).eq.ftet(j,2))
     2        	 .or.(tria(i,2).eq.ftet(j,3))
	     n3 = (tria(i,3).eq.ftet(j,1)).or.(tria(i,3).eq.ftet(j,2))
     3    		 .or.(tria(i,3).eq.ftet(j,3))
c check coincidence with old triangles -->true
	     IF(n1.and.n2.and.n3) then
		  new_flag=.false.
	     endif
	   enddo
	   if (new_flag) then
	     ns = ns+1
	     do k= 1,3
	       res(ns,k) = ftet(j,k)
	     enddo
c nodes are on the frontline -> multiply by 2
	     do k= 4,6
c	       res(ns,k) = 2*ttim*ftet(j,k)
	       res(ns,k) = ftet(j,k)
	     enddo

c compute new frontpoint
           np=0
           do k= 1,3  !lines
c compute direction vectors of lines 
            do l= 1,3   
	       if(k.eq.1) then  !line 1
              m(l) = xtria(j,l,2) - xtria(j,l,1) 
	        flag =  not(ftet(j,5).eq.ftet(j,4)) 
	       elseif(k.eq.2) then  !line 2
              m(l) = xtria(j,l,3) - xtria(j,l,2) 	       
	        flag =  not(ftet(j,6).eq.ftet(j,5)) 
	       elseif(k.eq.3) then  !line 3
              m(l) = xtria(j,l,1) - xtria(j,l,3) 
	        flag =  not(ftet(j,4).eq.ftet(j,6)) 
	       endif
	      enddo

            if (flag) then   ! consider only lines with different id
c compute crosspoint
              temp = m(1)*N0(1) + m(2)*N0(2) + m(3)*N0(3) 
	        If (Abs(temp).gt.eps) then ! plane is not parallel to m
	          mu = (N0(1)*(fp(1)-xtria(j,1,k)) + 
	1                N0(2)*(fp(2)-xtria(j,2,k)) + 
	1                N0(3)*(fp(3)-xtria(j,3,k)))/temp  
                 mu = dMax1(0.0001d0,mu)
                 mu = dMin1(0.9999d0,mu)
	           if ((mu.eq.0.0001d0).or.(mu.eq.0.9999d0)) then
	        write(*,*)   'WARNING --> Closest point approximation'
	        write(iow,*)   'WARNING --> Closest point approximation'
	           endif
	           np = np + 1
	           do l= 1,3
	             fp_new(ns,l) = fp_new(ns,l) + 
     1                          0.5d0*(xtria(j,l,k) + mu*m(l))
	           enddo
	        endif
	      endif
		 enddo
	     if(np.ne.2) then ! no intersection -> take center  
	       write(*,*)   'WARNING --> Center point approximation'
	       write(iow,*)   'WARNING --> Center point approximation'
	       fp_new=0.0d0
	       do k= 1,3
	         do l= 1,3
	           fp_new(ns,l) = fp_new(ns,l) + 
     1             0.3333333333d0*xtria(j,l,k)
	         enddo
	       enddo
           endif
c  check result
           do ii= 1,3
	       x(ii) = xtria(j,ii,3) - xtria(j,ii,1)
	       y(ii) = xtria(j,ii,2) - xtria(j,ii,1)
	       z(ii) = fp_new(ns,ii) - xtria(j,ii,1)
		 enddo 
		 temp = (-x(3)*y(2) + x(2)*y(3))*z(1) + 
     1		     (x(3)*y(1) - x(1)*y(3))*z(2) + 
     2            (-x(2)*y(1) + x(1)*y(2))*z(3) 
	     if(abs(temp).gt.eps) then
	       write(*,*)     'ERROR --> fp not on facet'
	       write(iow,*)   'ERROR --> fp not on facet'
             write(iow,1000) ((xl(ii,jj),ii=1,3),jj=1,4)
             write(iow,1001) (fp_new(ns,ii),ii=1,3),(fp(ii),ii=1,3), 
     1						(N0(ii),ii=1,3)
	     endif
	   endif
	 endif

	enddo

c define fp_new
	do i= 1, ns
         if ((fp_new(i,1).eq.0.0d0).and.(fp_new(i,2).eq.0.0d0).and.
     1       (fp_new(i,3).eq.0.0d0)) then
	    write(*,*)   'WARNING --> took old frontpoint'
	    write(iow,*)   'WARNING --> took old frontpoint'
          do l= 1,3
           fp_new(i,l) = fp(l)
          enddo
         endif
      enddo
1000  format(
     1  10x,'Tet data   '/
     2  10x,'node 1      ',e12.5,e12.5,e12.5/
     2  10x,'node 2      ',e12.5,e12.5,e12.5/
     2  10x,'node 3      ',e12.5,e12.5,e12.5/
     2  10x,'node 4      ',e12.5,e12.5,e12.5/)
1001  format(
     2  10x,'Plane data              '/
     2  10x,'node_old    ',e12.5,e12.5,e12.5/
     2  10x,'node_new    ',e12.5,e12.5,e12.5/
     2  10x,'normal      ',e12.5,e12.5,e12.5/)
      end



      subroutine cTip_ini (ndm,numel,nen,nen1,x,xlen,
	1                  ix,ixlen)
c__________________________________________________TCG. 18.04.2003___71
c
c
c  CrackTip initialization defined by ctip.dat
c
c    1				input model (-2= multiple automatic, -1= single automatic, 0= manuell, 1= no. of lines, 2,3= circle, )
c
c model=-2:      
c    blank line
c
c model=-1:
c    blank line
c
c model=0:
c    1               no of facets
c
c    1390 1391 1417  0 0 0  4.3654 3.3111 18.5233
c
c
c model=1:
c    1				no. of lines
c
c    222.4,-10,80		a(1)
c    222.4, 110,80	b(1)
c    1,0,0            n(1)
c
c    1.0d-6           tol  
c
c model=3:
c
c    0.0, 0.0,   154	     m
c    0.0, 0.0,   1.0	     Nc
c    200				r
c    1.0                 alpha (N = alpha Nc +(1-alpha)Nr)
c
c    2.0d0              tol
c
c
c  materials where the crack is allowed to propagate
c
c    ncmat
c
c    cmat1
c    cmat2
c    ....
c
c____________________________________________________________________71
      implicit none

	logical cflag

      include  'cTip.h'
      include  'iofile.h'

	integer xlen
	integer i, k, j, l, ndm, ixlen, nen, numel
      integer ix(ixlen), nodenr(4),  nen1, ncirc
	real*8 x(xlen), xl(3,4), facecor(3,3,4)
	integer io, faceno(3,4)
      integer itip, fid(3), ityp, nline
	real*8 a(3,100), b(3,100), n(3,100)
	real*8 max_tol, fp(3)
	real*8 cc(3,100),nc(3,100),r(100), alpha(100)
	integer mat




	itip=38

c read crackTip data from cTip.dat
      OPEN(UNIT=itip,FILE='CTip.dat',status='OLD',IOSTAT=io)
	IF(io.eq.0) THEN
c read inputtype 
     	  read (itip,*) ityp
	  read (itip,*) 
c automatic -2
        if (ityp.eq.-2) then 
	    write (iow,*) 'Multiple automatic crack tip initialisation'
	    imo = -2

c automatic -1
        elseif (ityp.eq.-1) then 
	    write (iow,*) 'Single automatic crack tip initialisation'
	    imo = -1

c manuell
        elseif (ityp.eq.0) then 
	    imo = 0
	    read (itip,*) ncTip
	    read (itip,*) 
	    do j= 1, ncTip
            read(itip,*) (cTip(j,k),k=1,6), (cTip_point(j,k),k=1,3) 
	    enddo
c number of lines
        elseif (ityp.eq.1) then 
	    imo = 1
	    read (itip,*) nline
	    write (iow,*) 'no. of lines', nline              ! no. of circles
	    read (itip,*) 
	    do j= 1, nline
	      read (itip,*) (a(i,j),i=1,3)
	      read (itip,*) (b(i,j),i=1,3)
	      read (itip,*) (n(i,j),i=1,3)
	      read (itip,*) 

	      write (iow,*) 'Point a',(a(i,j),i=1,3)
	      write (iow,*) 'Point b',(b(i,j),i=1,3)
	      write (iow,*) 'Normal',(n(i,j),i=1,3)
          enddo
c circle a
	  elseif (ityp.eq.2) then
	    imo = 2
	    read (itip,*) ncirc              ! no. of circles
	    write (iow,*) 'no. of circles', ncirc              ! no. of circles
	    read (itip,*) 
	    do j= 1, ncirc
	      read (itip,*) (cc(i,j),i=1,3)    ! center of circle
	      read (itip,*) (nc(i,j),i=1,3)    ! normal of circle
	      read (itip,*) r(j)               ! radius of circle
	      read (itip,*) alpha(j)           ! crack tip normal of circle (not necessary)

	      write (iow,*) 'center', (cc(i,j),i=1,3)    ! center of circle
	      write (iow,*) 'normal', (nc(i,j),i=1,3)    ! normal of circle
	      write (iow,*) 'radius', r(j)               ! radius of circle
	      write (iow,*) 'normal', alpha(j)           ! crack tip normal of circle (not necessary)
          enddo
c circle b
	  elseif (ityp.eq.3) then
	    imo = 3
	    read (itip,*) ncirc              ! no. of circles
	    write (iow,*) 'no. of circles', ncirc              ! no. of circles
	    read (itip,*) 
	    do j= 1, ncirc
	      read (itip,*) (cc(i,j),i=1,3)    ! center of circle
	      read (itip,*) (nc(i,j),i=1,3)    ! normal of circle
	      read (itip,*) r(j)               ! radius of circle
	      read (itip,*) alpha(j)           ! crack tip normal of circle (not necessary)

	      write (iow,*) 'center', (cc(i,j),i=1,3)    ! center of circle
	      write (iow,*) 'normal', (nc(i,j),i=1,3)    ! normal of circle
	      write (iow,*) 'radius', r(j)               ! radius of circle
	      write (iow,*) 'normal', alpha(j)           ! crack tip normal of circle (not necessary)
          enddo
	  else
c type not known
          write (1001,*)
	  endif
c read accuracy
	  if (ityp.gt.0) then
	    read (itip,*) max_tol   ! tolerance 
          write (iow,*) 'Tolerance',	 max_tol 
	  endif            
c read materials where the crack is allowed to propagate
	  read (itip,*) 
     	  read (itip,*) ncmat
	  write (iow,*) 'no. of crack mat.', ncmat             
	  read (itip,*) 
	    do j= 1, ncmat
            read(itip,*) cmat(j) 
c store automatic initialisation in cmat
c		  if(ityp.eq.-1) cmat(j) = -cmat(j)           
c		  if(ityp.eq.-2) cmat(j) = -cmat(j)           
	      write (iow,*) 'crack mat.', cmat(j) 
	    enddo
        CLOSE(itip)
      ELSE
	    goto 1000
      ENDIF
      if (ityp.eq.0) then 
c all done
      elseif (ityp.eq.1) then 
c number of lines
c loop over lines
        do l= 1,nline 
c get mesh data
         do i= 1,numel
c retrieve node no. and lagrangian coordinates
c          do j = 1,nen
          do j = 1,4
	     nodenr(j) = ix((i-1)*nen1 + j)
	     do k= 1,ndm
	       xl(k,j) = x((nodenr(j)-1)*ndm + k) 
	     enddo
          enddo
c filter out crack tip facets
c set up tet faces
          faceno(1,1) = nodenr(1)
          faceno(2,1) = nodenr(2)
          faceno(3,1) = nodenr(3)
	    do j=1,3
            facecor(j,1,1) = xl(j,1)
            facecor(j,2,1) = xl(j,2)
            facecor(j,3,1) = xl(j,3)
	    enddo

          faceno(1,2) = nodenr(1)
          faceno(2,2) = nodenr(2)
          faceno(3,2) = nodenr(4)
	    do j=1,3
            facecor(j,1,2) = xl(j,1)
            facecor(j,2,2) = xl(j,2)
            facecor(j,3,2) = xl(j,4)
	    enddo

          faceno(1,3) = nodenr(2)
          faceno(2,3) = nodenr(3)
          faceno(3,3) = nodenr(4)
	    do j=1,3
            facecor(j,1,3) = xl(j,2)
            facecor(j,2,3) = xl(j,3)
            facecor(j,3,3) = xl(j,4)
	    enddo

          faceno(1,4) = nodenr(1)
          faceno(2,4) = nodenr(3)
          faceno(3,4) = nodenr(4)
	    do j=1,3
            facecor(j,1,4) = xl(j,1)
            facecor(j,2,4) = xl(j,3)
            facecor(j,3,4) = xl(j,4)
	    enddo

c compare faces with line
          do j= 1,4
c lines no. k 
	      call face_cross_line (n(1,l),a(1,l),b(1,l),max_tol,
	1                     facecor(1,1,j), cflag, fid, fp)
            if (cflag) then
c write nodes in ctip
	  	    ncTip = ncTip + 1
	        do k=1,3
                cTip(ncTip,k) = faceno(k,j) 
                cTip(ncTip,k+3) = fid(k) 
	          cTip_point(ncTip,k) = fp(k)
	        enddo
	      endif
	    enddo

	   enddo
	  enddo
      elseif (ityp.eq.2) then 
c circle a
c loop over circles
        do l= 1,ncirc 
c get mesh data
         do i= 1,numel
c retrieve node no. and lagrangian coordinates
c          do j = 1,nen
          do j = 1,4
	     nodenr(j) = ix((i-1)*nen1 + j)
	     do k= 1,ndm
	       xl(k,j) = x((nodenr(j)-1)*ndm + k) 
	     enddo
          enddo
c filter out crack tip facets
c set up tet faces
          faceno(1,1) = nodenr(1)
          faceno(2,1) = nodenr(2)
          faceno(3,1) = nodenr(3)
	    do j=1,3
            facecor(j,1,1) = xl(j,1)
            facecor(j,2,1) = xl(j,2)
            facecor(j,3,1) = xl(j,3)
	    enddo

          faceno(1,2) = nodenr(1)
          faceno(2,2) = nodenr(2)
          faceno(3,2) = nodenr(4)
	    do j=1,3
            facecor(j,1,2) = xl(j,1)
            facecor(j,2,2) = xl(j,2)
            facecor(j,3,2) = xl(j,4)
	    enddo

          faceno(1,3) = nodenr(2)
          faceno(2,3) = nodenr(3)
          faceno(3,3) = nodenr(4)
	    do j=1,3
            facecor(j,1,3) = xl(j,2)
            facecor(j,2,3) = xl(j,3)
            facecor(j,3,3) = xl(j,4)
	    enddo

          faceno(1,4) = nodenr(1)
          faceno(2,4) = nodenr(3)
          faceno(3,4) = nodenr(4)
	    do j=1,3
            facecor(j,1,4) = xl(j,1)
            facecor(j,2,4) = xl(j,3)
            facecor(j,3,4) = xl(j,4)
	    enddo

c compare faces with circle
          do j= 1,4
	      call face_cross_circlea (cc(1,l),nc(1,l),r(l),alpha(l),
	1                    max_tol,facecor(1,1,j), cflag, fid, fp)
            if (cflag) then
c write nodes in ctip
	  	    ncTip = ncTip + 1
	        do k=1,3
                cTip(ncTip,k) = faceno(k,j) 
                cTip(ncTip,k+3) = fid(k) 
	          cTip_point(ncTip,k) = fp(k)
	        enddo
	      endif
	    enddo

	   enddo
	  enddo
      elseif (ityp.eq.3) then 
c circle b
c loop over circles
        do l= 1,ncirc 
c get mesh data
         do i= 1,numel
c retrieve node no. and lagrangian coordinates
c          do j = 1,nen
          do j = 1,4
	     nodenr(j) = ix((i-1)*nen1 + j)
	     do k= 1,ndm
	       xl(k,j) = x((nodenr(j)-1)*ndm + k) 
	     enddo
          enddo
	    mat = ix((i-1)*nen1 + 11)
	    if (mat.eq.1) then
c filter out crack tip facets
c set up tet faces
            faceno(1,1) = nodenr(1)
            faceno(2,1) = nodenr(2)
            faceno(3,1) = nodenr(3)
	      do j=1,3
              facecor(j,1,1) = xl(j,1)
              facecor(j,2,1) = xl(j,2)
              facecor(j,3,1) = xl(j,3)
	      enddo

            faceno(1,2) = nodenr(1)
            faceno(2,2) = nodenr(2)
            faceno(3,2) = nodenr(4)
	      do j=1,3
              facecor(j,1,2) = xl(j,1)
              facecor(j,2,2) = xl(j,2)
              facecor(j,3,2) = xl(j,4)
	      enddo

            faceno(1,3) = nodenr(2)
            faceno(2,3) = nodenr(3)
            faceno(3,3) = nodenr(4)
	      do j=1,3
              facecor(j,1,3) = xl(j,2)
              facecor(j,2,3) = xl(j,3)
              facecor(j,3,3) = xl(j,4)
	      enddo

            faceno(1,4) = nodenr(1)
            faceno(2,4) = nodenr(3)
            faceno(3,4) = nodenr(4)
	      do j=1,3
              facecor(j,1,4) = xl(j,1)
              facecor(j,2,4) = xl(j,3)
              facecor(j,3,4) = xl(j,4)
	      enddo

c compare faces with circle
            do j= 1,4
	        call face_cross_circleb (cc(1,l),nc(1,l),r(l),alpha(l),
	1                     max_tol,facecor(1,1,j), cflag, fid, fp)
              if (cflag) then
c write nodes in ctip
	  	      ncTip = ncTip + 1
	          do k=1,3
                  cTip(ncTip,k) = faceno(k,j) 
                  cTip(ncTip,k+3) = fid(k) 
	            cTip_point(ncTip,k) = fp(k)
	          enddo
	        endif
	      enddo
          endif
	   enddo
	 enddo
	endif

1001  format('** Cracktip model not implemented!')
1000  end





      subroutine face_cross_line (n0,Xb,XXb,max_tol,xt,cflag,
	1                              fid, fp)
c__________________________________________________TCG. 25.11.2002___71
c
c
c  compute if two sides triangle are crossed by a line 
c
c
c   n0        normal onto initial crack surface
c   Xb        first point of line
c   XXb       secound point of line
c   max_tol   toleranz value
c   xt        nodes of actuell triangle
c   cflag     cut flag   .true. -> line crosses actuell triangle
c                        .false. -> line does not cross actuell triangle
c   fid       node id (positive and negative side)
c   fp        frontpoint
c
c
c____________________________________________________________________71

      implicit none

	integer i,j,k, cont, fid(3)

	real*8 n0(3), Xb(3), xxb(3), xt(3,3),   det
	real*8 r,s, max_tol, temp, fp(3)
	real*8 ab, rhs(2), a(3), b(3), Yb(3)
	real*8 YbmXb(3), eps, vect(3), smax, rmax
	real*8 xtma(3)

	logical cflag


      cont = 0
c compute direction vector of line a-b
      a = XXb - Xb
c normalize 
      rmax= sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      a = a/rmax

	fp=0.0d0
      
	do i= 1,3    ! lines
       do k= 1,3   
	    if(i.eq.1) then  !line 1
            b(k) = xt(k,2) - xt(k,1) 
	      Yb(k) = xt(k,1)
	    elseif(i.eq.2) then  !line 2
            b(k) = xt(k,3) - xt(k,2) 	       
	      Yb(k) = xt(k,2)
	    elseif(i.eq.3) then  !line 3
            b(k) = xt(k,1) - xt(k,3) 
	      Yb(k) = xt(k,3)
	    endif
	 enddo
c normalize 
      smax= sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
      b = b/smax
c compute n.m
          ab = a(1)*b(1) + a(2)*b(2) + a(3)*b(3) 
          det = 1.0d0 - ab*ab
c compute rhs
          YbmXb = Yb - Xb
          rhs(1) =  YbmXb(1)*b(1) + YbmXb(2)*b(2) + YbmXb(3)*b(3)
          rhs(2) =  YbmXb(1)*a(1) + YbmXb(2)*a(2) + YbmXb(3)*a(3)
	    If (abs(det).gt.1.0d-8) then !line not parallel to n
c compute r and s
            r = (-ab*rhs(1) + rhs(2))/det
            s = (-rhs(1) + ab*rhs(2))/det
	      if ((r.gt.0.0d0).and.(r.lt.rmax).and.
     1          (s.gt.0.0d0).and.(s.lt.smax)) then
c minimum distance
              vect = Yb-Xb + s*b - r*a
	        eps= sqrt(vect(1)*vect(1) + vect(2)*vect(2) + 
	1                  vect(3)*vect(3))
	        If (eps.lt.max_tol) then  ! faceside crosses line 
                cont = cont +1
	          do k= 1,3
	            fp(k) = fp(k) + 0.5d0*(Yb(k) + s*b(k))
	          enddo
	        endif
	      endif
	    endif
      enddo   
	
	if (cont.eq.2) then
	  cflag=.true.  
	else
	  cflag=.false.
	endif
c compute face id if triangle is crossed
      if (cflag) then
        do i= 1,3
	   do j= 1,3
	     xtma(j) = xt(j,i) - fp(j)
	   enddo
	   temp = xtma(1)*N0(1) + xtma(2)*N0(2) + xtma(3)*N0(3) 
	   if (temp.gt.0) then
c	     fid(i) = -2
	     fid(i) = 0
	   else
c	     fid(i) = +2
	     fid(i) = 0
	   endif
	  enddo
      endif

      end





      subroutine face_cross_circlea (c,n,r,alpha,max_tol,xt,cflag,
	1                              fid, fp)
c__________________________________________________TCG. 18.06.2003___71
c
c
c  compute if two sides triangle are crossed by a circle
c
c
c   c         Center of circle
c   n         Normal of circle
c   r         Radius of circle
c   alpha     normal factor at cracktip of circle
c   max_tol   tolerance value
c   xt        nodes of actuell triangle
c   cflag     cut flag   .true. -> line crosses actuell triangle
c                        .false. -> line does not cross actuell triangle
c   fid       node id (positive and negative side)
c   fp        frontpoint
c
c
c____________________________________________________________________71

      implicit none

	integer i,j,k, cont, fid(3)

	real*8 c(3), n(3), r, xt(3,3)
	real*8 max_tol, temp, fp(3)
	real*8 a(3), eps, alpha
      real*8 ac, cc, aa, lambda(2), x, xtma(3)
	real*8 a1(3), a2(3), a1n, a2n
	real*8 xb(3), xbxb, axb, cxb, Nr(3), n0(3)

	logical cflag, plane

      data eps/1.0d-8/
      cont = 0

c normalize 

	fp=0.0d0
      
	do i= 1,3    ! lines
       do k= 1,3   
	    if(i.eq.1) then      !line 1
	      xb(k)= xt(k,1)
            a(k) = xt(k,2) - xt(k,1) 
	      a1(k)= xt(k,1) - c(k)
	      a2(k)= xt(k,2) - c(k)
	    elseif(i.eq.2) then  !line 2
	      xb(k)= xt(k,2)
            a(k) = xt(k,3) - xt(k,2) 	       
	      a1(k)= xt(k,2) - c(k)
	      a2(k)= xt(k,3) - c(k)
	    elseif(i.eq.3) then  !line 3
	      xb(k)= xt(k,3)
            a(k) = xt(k,1) - xt(k,3) 
	      a1(k)= xt(k,3) - c(k)
	      a2(k)= xt(k,1) - c(k)
	    endif
	 enddo
c compute inproducts
       a1n = a1(1)*n(1)  + a1(2)*n(2)  + a1(3)*n(3)  
       a2n = a2(1)*n(1)  + a2(2)*n(2)  + a2(3)*n(3)  
       ac  = a(1)*c(1)   + a(2)*c(2)   + a(3)*c(3) 
       cc  = c(1)*c(1)   + c(2)*c(2)   + c(3)*c(3) 
       aa  = a(1)*a(1)   + a(2)*a(2)   + a(3)*a(3) 
       axb = a(1)*xb(1)  + a(2)*xb(2)  + a(3)*xb(3) 
       cxb = c(1)*xb(1)  + c(2)*xb(2)  + c(3)*xb(3) 
       xbxb= xb(1)*xb(1) + xb(2)*xb(2) + xb(3)*xb(3) 

c plane
       plane = (((a1n*a2n).lt.0.0d0).or.
	1          (abs(a1n)+abs(a2n)).lt.max_tol)
        if(plane) then   !linesegment is in plane of interesst
	 
c solve quadratic eq.
	   if (aa.gt.eps) then                 ! length of element side .gt. zero 
c compute root
          temp = (axb-ac)*(axb-ac) - aa*(cc - 2.0d0*cxb + xbxb - r*r)
	     if (abs(max_tol)+temp.gt.0.0d0) then    ! solutions exist
	      if (temp.lt.0.0d0) temp = 0.0d0
	      temp = sqrt(temp)
            lambda(1) = (ac-axb + temp)/aa
            lambda(2) = (ac-axb - temp)/aa
c find solution of interesst
	      x = -1.0d0       ! no solution
            if((lambda(1).gt.0.0d0).and.(lambda(1).lt.1.0d0)) then 
		   if((lambda(2).gt.0.0d0).and.(lambda(2).lt.1.0d0)) then  
	         x = (lambda(1) + lambda(2))/2.0d0 
	       else
	          x = lambda(1)
	       endif
            elseif((lambda(2).gt.0.0d0).and.(lambda(2).lt.1.0d0)) then 
		   if((lambda(1).gt.0.0d0).and.(lambda(1).lt.1.0d0)) then  
	        x = (lambda(1) + lambda(2))/2.0d0 
	       else
	        x = lambda(2)
	       endif
	      endif
c compute data 
            if (x.gt.0.0d0) then
	        cont = cont + 1
	        do k= 1,3
	         fp(k) = fp(k) + xb(k) + x*a(k)
	        enddo
	      endif
	     endif
          endif
         endif
      enddo   
	
	if ((cont.eq.2).or.(cont.eq.3)) then
	  cflag=.true.  
	  fp = fp/dfloat(cont)
c compute face id if triangle is crossed
c compute normal at cracktip
        Nr = fp - c
	  temp = sqrt(Nr(1)*Nr(1) + Nr(2)*Nr(2) +Nr(3)*Nr(3)) 
	  Nr = Nr/temp
        N0 = alpha*N + (1.0d0-alpha)*Nr
        if (cflag) then
          do i= 1,3
	     do j= 1,3
	       xtma(j) = xt(j,i) - fp(j)
	     enddo
	     temp = xtma(1)*N0(1) + xtma(2)*N0(2) + xtma(3)*N0(3) 
	     if (temp.gt.0) then
c	       fid(i) = -2
	       fid(i) = 0
	     else
c	       fid(i) = +2
	       fid(i) = 0
	     endif
	    enddo
        endif
	else
	  cflag=.false.
	endif

      end






      subroutine face_cross_circleb (c,n,r,alpha,max_tol,xt,cflag,
	1                              fid, fp)
c__________________________________________________TCG. 23.06.2003___71
c
c
c  compute if two sides triangle are crossed by a circle
c
c
c   c         Center of circle
c   n         Normal of circle
c   r         Radius of circle
c   alpha     normal factor at cracktip of circle
c   max_tol   tolerance value
c   xt        nodes of actuell triangle
c   cflag     cut flag   .true. -> line crosses actuell triangle
c                        .false. -> line does not cross actuell triangle
c   fid       node id (positive and negative side)
c   fp        frontpoint
c
c
c____________________________________________________________________71

      implicit none

	integer i,j,k, cont, fid(3)

	real*8 c(3), n(3), r, xt(3,3)
	real*8 max_tol, temp, fp(3)
	real*8 a(3), eps, alpha
      real*8 lambda, xtma(3)
	real*8 an, bn, tt1(3), tt2(3), ntt1, ntt2
	real*8 xb(3), Nr(3), n0(3)

	logical cflag

      data eps/1.0d-8/

      cont = 0
	fp=0.0d0

c normalize 

      
	do i= 1,3    ! lines
       do k= 1,3   
	    if(i.eq.1) then      !line 1
	      xb(k)= xt(k,1)
            a(k) = xt(k,2) - xt(k,1) 
	    elseif(i.eq.2) then  !line 2
	      xb(k)= xt(k,2)
            a(k) = xt(k,3) - xt(k,2) 	       
	    elseif(i.eq.3) then  !line 3
	      xb(k)= xt(k,3)
            a(k) = xt(k,1) - xt(k,3) 
	    endif
	 enddo
c compute inproducts
       an = a(1)*n(1)  + a(2)*n(2)  + a(3)*n(3)  
       bn = (xb(1)-c(1))*n(1)+(xb(2)-c(2))*n(2)+(xb(3)-c(3))*n(3)  

c compute lambda
       if (abs(an).gt.eps) then
	  lambda = -bn/an 
	  if ((lambda.gt.0.0d0).and.(lambda.lt.1.0d0)) then
	    tt1 = xb - c
	    ntt1 = sqrt(tt1(1)*tt1(1) + tt1(2)*tt1(2) + tt1(3)*tt1(3)) 
	    tt2 = xb + a - c
	    ntt2 = sqrt(tt2(1)*tt2(1) + tt2(2)*tt2(2) + tt2(3)*tt2(3)) 

	    temp= (1.0d0-lambda)*ntt1 + lambda*ntt2 
	   if ((abs(temp).lt.r+max_tol).and.
	1       (abs(temp).gt.r-max_tol))  then
	     cont=cont+1
	     fp = fp + xb + lambda*a
         endif
	  endif
	 endif
      enddo   
	
	if (cont.eq.2) then
	  cflag=.true.  
	  fp = fp/2.0d0
c compute face id if triangle is crossed
c compute normal at cracktip
        Nr = fp - c
	  temp = sqrt(Nr(1)*Nr(1) + Nr(2)*Nr(2) +Nr(3)*Nr(3)) 
	  Nr = Nr/temp
        N0 = alpha*N + (1.0d0-alpha)*Nr
        if (cflag) then
          do i= 1,3
	     do j= 1,3
	       xtma(j) = xt(j,i) - fp(j)
	     enddo
	     temp = xtma(1)*N0(1) + xtma(2)*N0(2) + xtma(3)*N0(3) 
	     if (temp.gt.0) then
c	       fid(i) = -2
	       fid(i) = 0
	     else
c	       fid(i) = +2
	       fid(i) = 0
	     endif
	    enddo
        endif
	else
	  cflag=.false.
	endif

      end

      subroutine smooth_11 (cn, xn3d, fp, Ve, N0, dN)
c__________________________________________________TCG. 09.08.2004___71
c
c
c  smooth cracksurface
c
c
c      Ae    elementvolume
c      cn    no. of corners at discontinuity
c      xn3d  coordinates of corners
c      fp    point on discontinuity
c      N0    Normal of discontinuity
c      R2    square nonlocal radius
c      an    no. of nodes to be included in the approximation
c      anode nodes to be included in the approximation
c      cc    centerof node cloud
c      nfit  order of surface to be fitted to data
c
c
c____________________________________________________________________71
	implicit none

      include  'nLength.h'
      include  'cGeom.h'
      include  'iofile.h'

	integer i,j,k, cn, an, itemp, info

	real*8 xn3d(3,4), R2, anode(3,1010), temp, delta(3)
      real*8 fp(3), xi(3), cc(3), Ixx, Iyy, Izz, Ixy, Ixz, Iyz
	real*8 mm(3,3), II(3), N0(3), Ve, QQ(3,3), tempv(3)
	real*8 para(6), dzdx, dzdy, dN, Nold(3), fpb(3)

c define nonlocal radius
      if (ichl.eq.1) then
	 R2 = chr**2
	else
       R2 = chr**2*Ve**0.6666666667d0
	endif


c step I  Define nodes to be included
      an = 0
      do i= 1, nsurf
	  do j = 1, ncorn(nsurf)
	    do k= 1,3
	      xi(k) = cnode(i,(j-1)*3 + k)
	    enddo
	    delta = xi - fp
	    temp = delta(1)**2 + delta(2)**2 + delta(3)**2
	    if(temp.lt.R2) then 
	      an = an + 1
	      if (an.le.1000) then
	        do k= 1,3
	          anode(k,an) = xi(k)
	        enddo
	      else 	        
	        write (*,*)
	        write (*,1000) 
	        write (*,*)
	        write (iow,*)
	        write (iow,1000)
	        write (iow,*)
	      endif
	    endif
	  enddo
	enddo



	if (an.gt.8) then ! perform analysis only if enough nodes available
c add nodes of actual discontinuity
      do i= 1, cn
	  an = an + 1
	  do k= 1,3
	    anode(k,an) = xn3d(k,i)
	  enddo
	enddo

c add frontpoint
	an = an + 1
	 do k= 1,3
	   anode(k,an) = fp(k)
	enddo

c step II define local coordinate system
c (a) compute center
      cc = 0.0d0

	do i= 1, an
	  do j= 1,3
	    cc(j) = cc(j) + anode(j,i)
	  enddo
	enddo
	cc = cc/dfloat(an)


c shift nodes
	do i= 1, an
	  do j= 1,3
	    anode(j,i) = anode(j,i) - cc(j)
	  enddo
	enddo


c (b) compute eigenvectors and values
      Ixx =0.0d0
      Iyy =0.0d0
      Izz =0.0d0
      Ixy =0.0d0
      Ixz =0.0d0
      Iyz =0.0d0

	do i= 1,an
	  Ixx = Ixx + anode(1,i)*anode(1,i)
	  Iyy = Iyy + anode(2,i)*anode(2,i)
	  Izz = Izz + anode(3,i)*anode(3,i)
	  Ixy = Ixy + anode(1,i)*anode(2,i)
	  Ixz = Ixz + anode(1,i)*anode(3,i)
	  Iyz = Iyz + anode(2,i)*anode(3,i)
	enddo

	mm(1,1) = Ixx
	mm(1,2) = Ixy
	mm(1,3) = Ixz
	mm(2,1) = Ixy
	mm(2,2) = Iyy
	mm(2,3) = Iyz
	mm(3,1) = Ixz
	mm(3,2) = Iyz
	mm(3,3) = Izz

      call eigen_11(mm,II,itemp)
      call sort_11(mm,II)

c (c) compute coordinates in local frame
  
      do i= 1,3
	  do j= 1,3
	    QQ(i,j) = mm(j,i)
	  enddo
	enddo



	do i= 1,an
	  tempv = 0.0d0
	  do j= 1,3
	    do k= 1,3
	      tempv(j) = tempv(j) + QQ(j,k)*anode(k,i)
	    enddo
	  enddo
	  do j= 1,3
	    anode(j,i) = tempv(j)
	  enddo
	enddo


c step III fit data to surface
c get frontpoint in local frame
      do i= 1,3
	  fpb(i) = anode(i,an)
	enddo

      if (nfit.eq.1)    then  ! linear surface
        call fit_lin_11(an, anode, para, info) 
c step IV redefine normal
        Nold = N0
        if (info.eq.0) then
	    dzdx = para(2) 
     	    dzdy = para(3) 
	    temp = sqrt(dzdx**2 + dzdy**2 + 1.0d0)

          N0(1) = -dzdx/temp
	    N0(2) = -dzdy/temp
	    N0(3) = 1.0d0/temp
        endif
	elseif(nfit.eq.2) then  ! quadratic surface
        call fit_quad_11(an, anode, para, info)       
c step IV redefine normal
        Nold = N0
        if (info.eq.0) then
	    dzdx = para(2) + 2.0d0*para(4)*fpb(1) + para(6)*fpb(2)
     	    dzdy = para(3) + 2.0d0*para(5)*fpb(2) + para(6)*fpb(1)
	    temp = sqrt(dzdx**2 + dzdy**2 + 1.0d0)

          N0(1) = -dzdx/temp
	    N0(2) = -dzdy/temp
	    N0(3) = 1.0d0/temp
        endif
	else
	  write (*,*)
	  write (*,1002) 
	  write (iow,*)
	  write (iow,1002)
      endif


c step IV redefine normal
      Nold = N0
      if (info.eq.0) then
c rotate N0 into global frame
	  tempv = 0.0d0
	  do i= 1,3
	    do j= 1,3
	      tempv(i) = tempv(i) + QQ(j,i)*N0(j)
	    enddo
	  enddo
	  N0 = tempv
      endif
      dn = sqrt((Nold(1)-N0(1))**2 + (Nold(2)-N0(2))**2 + 
	1          (Nold(3)-N0(3))**2)
	else
	endif

1000  format('==============================================='/
     1       'WARNING --> More than 1000 neighbouring nodes'/
     3       '===============================================')
1001  format(e13.5,e13.5,e13.5,e13.5,e13.5,e13.5)
1002  format('==============================================='/
     1       'ERROR --> Approximation surface not available'/
     3       '===============================================')
	end

      subroutine sort_11(n, lamb)
c__________________________________________________TCG. 09.08.2004___71
c
c  sort eigenvalues
c
c__________________________________________________TCG. 09.08.2004___71

      implicit none

	logical abbr

	integer i, j

	real*8 n(3,3), lamb(3), temp

       
c sort values 
      
      abbr=.false.
      do while(not(abbr))
	  abbr =.true.
	  do i=1,2
	    if (lamb(i).lt.lamb(i+1)) then
	      temp      = lamb(i)
	      lamb(i)   = lamb(i+1)
	      lamb(i+1) = temp
	      do j= 1,3
	        temp     = n(j,i)
	        n(j,i)   = n(j,i+1)
	        n(j,i+1) = temp
	      enddo
	      abbr      =.false.
	    endif
	  enddo
	enddo

	end

      subroutine fit_quad_11(an, anode, r, info) 
c__________________________________________________TCG. 09.08.2004___71
c
c  fit quadratic surface to data
c  z = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y
c
c least square problem leading to a linear symmetric system for the
c coefficients a0, ..., a5
c
c   an     no. of nodes
c   anode  coorsinates of nodes
c   r      coefficients
c   info   info about the system of eqns
c          info=0 is ok
c
c__________________________________________________TCG. 09.08.2004___71
      implicit none

      integer an, i, j, kpvt(6),info

	real*8 anode(3,an), k(6,6), r(6)
	real*8 x4(an), x3(an), x3y(an), x2(an), x2y(an), x2y2(an)
      real*8 x(an), xy(an), xy2(an), xy3(an) 
	real*8 y(an), y2(an), y3(an), y4(an), z(an), xz(an), x2z(an)
	real*8 yz(an), y2z(an), xyz(an)

	real*8 sx4, sx3, sx3y, sx2, sx2y, sx2y2, sx, sxy, sxy2, sxy3 
	real*8 sy, sy2, sy3, sy4, sz, sxz, sx2z, syz, sy2z, sxyz

c compute scalars
      do i= 1,an
	  x(i)   = anode(1,i)
	  y(i)   = anode(2,i)
	  z(i)   = anode(3,i)
	  x2(i)  = x(i)*x(i)
	  x3(i)  = x2(i)*x(i)
	  x4(i)  = x3(i)*x(i)
	  y2(i)  = y(i)*y(i)
	  y3(i)  = y2(i)*y(i)
	  y4(i)  = y3(i)*y(i)
	  xy(i)  = x(i)*y(i)
	  x2y(i) = x2(i)*y(i)
	  x3y(i) = x3(i)*y(i)
	  xy2(i) = xy(i)*y(i)
	  xy3(i) = xy2(i)*y(i)
	  x2y2(i)= x2(i)*y2(i)
	  xz(i)  = x(i)*z(i)
	  yz(i)  = y(i)*z(i)
	  x2z(i) = x2(i)*z(i)
	  y2z(i) = y2(i)*z(i)
	  xyz(i) = xy(i)*z(i)

	enddo

c compute summs

	sx   = 0.0d0
	sy   = 0.0d0
	sz   = 0.0d0
	sx2  = 0.0d0
	sx3  = 0.0d0
	sx4  = 0.0d0
	sy2  = 0.0d0
	sy3  = 0.0d0
	sy4  = 0.0d0
	sxy  = 0.0d0
	sx2y = 0.0d0
	sx3y = 0.0d0
	sxy2 = 0.0d0
	sxy3 = 0.0d0
	sx2y2= 0.0d0
	sxz  = 0.0d0
	syz  = 0.0d0
	sx2z = 0.0d0
	sy2z = 0.0d0
	sxyz = 0.0d0

	do i= 1,an
	sx   = sx + x(i)
	sy   = sy + y(i)
	sz   = sz + z(i)
	sx2  = sx2 + x2(i)
	sx3  = sx3 + x3(i)
	sx4  = sx4 + x4(i)
	sy2  = sy2 + y2(i)
	sy3  = sy3 + y3(i)
	sy4  = sy4 + y4(i)
	sxy  = sxy + xy(i)
	sx2y = sx2y+ x2y(i)
	sx3y = sx3y+ x3y(i)
	sxy2 = sxy2+ xy2(i)
	sxy3 = sxy3+ xy3(i)
	sx2y2= sx2y2+x2y2(i)
	sxz  = sxz + xz(i)
	syz  = syz + yz(i)
	sx2z = sx2z+ x2z(i)
	sy2z = sy2z+ y2z(i)
	sxyz = sxyz+ xyz(i)
	enddo

c assemble matrix
      k(1,1) = dfloat(an)
	k(1,2) = sx
	k(1,3) = sy
	k(1,4) = sx2
	k(1,5) = sy2
	k(1,6) = sxy

	k(2,2) = sx2
	k(2,3) = sxy
	k(2,4) = sx3
	k(2,5) = sxy2
	k(2,6) = sx2y

	k(3,3) = sy2
	k(3,4) = sx2y
	k(3,5) = sy3
	k(3,6) = sxy2

	k(4,4) = sx4
	k(4,5) = sx2y2
	k(4,6) = sx3y

	k(5,5) = sy4
	k(5,6) = sxy3

	k(6,6) = sx2y2

c sym
      do i= 1,6
	  do j= i+1,6
	    k(j,i) = k(i,j)
	  enddo
	enddo

c assemble vector
      r(1) = sz
      r(2) = sxz
      r(3) = syz
      r(4) = sx2z
      r(5) = sy2z
      r(6) = sxyz


c write
c Set counter i
c      OPEN(UNIT=47,FILE='kk')
c	do i= 1,6
c	  write(47,1000) (k(i,j),j=1,6)
c	enddo
c	CLOSE(47)
c
c      OPEN(UNIT=47,FILE='rr')
c	do i= 1,6
c	  write(47,1000) r(i)
c	enddo
c	CLOSE(47)


c solve linear system using linpack
c factorize
      call dsifa(k,6,6,kpvt,info)
c solve
      if (info.eq.0) then
        call dsisl(k,6,6,kpvt,r)
      endif

c      OPEN(UNIT=47,FILE='para')
c	do i= 1,6
c	  write(47,1000) r(i)
c	enddo
c	CLOSE(47)


1000  format(e13.5,e13.5,e13.5,e13.5,e13.5,e13.5)
      end



      subroutine fit_lin_11(an, anode, r, info) 
c__________________________________________________TCG. 09.08.2004___71
c
c  fit linear surface to data
c  z = a0 + a1*x + a2*y 
c
c least square problem leading to a linear symmetric system for the
c coefficients a0, ..., a2
c
c   an     no. of nodes
c   anode  coorsinates of nodes
c   r      coefficients
c   info   info about the system of eqns
c          info=0 is ok
c
c__________________________________________________TCG. 09.08.2004___71
      implicit none

      integer an, i, j, kpvt(3),info

	real*8 anode(3,an), k(3,3), r(3)
	real*8 x2(an),x(an), xy(an),y(an), y2(an), z(an), xz(an),yz(an)
	real*8 sx2, sx, sxy, sy, sy2, sz, sxz, syz

c compute scalars
      do i= 1,an
	  x(i)   = anode(1,i)
	  y(i)   = anode(2,i)
	  z(i)   = anode(3,i)
	  x2(i)  = x(i)*x(i)
	  y2(i)  = y(i)*y(i)
	  xy(i)  = x(i)*y(i)
	  xz(i)  = x(i)*z(i)
	  yz(i)  = y(i)*z(i)
	enddo

c compute summs

	sx   = 0.0d0
	sy   = 0.0d0
	sz   = 0.0d0
	sx2  = 0.0d0
	sy2  = 0.0d0
	sxy  = 0.0d0
	sxz  = 0.0d0
	syz  = 0.0d0

	do i= 1,an
	sx   = sx + x(i)
	sy   = sy + y(i)
	sz   = sz + z(i)
	sx2  = sx2 + x2(i)
	sy2  = sy2 + y2(i)
	sxy  = sxy + xy(i)
	sxz  = sxz + xz(i)
	syz  = syz + yz(i)
	enddo

c assemble matrix
      k(1,1) = dfloat(an)
	k(1,2) = sx
	k(1,3) = sy

	k(2,2) = sx2
	k(2,3) = sxy

	k(3,3) = sy2

c sym
      do i= 1,3
	  do j= i+1,3
	    k(j,i) = k(i,j)
	  enddo
	enddo

c assemble vector
      r(1) = sz
      r(2) = sxz
      r(3) = syz



c solve linear system using linpack
c factorize
      call dsifa(k,3,3,kpvt,info)
c solve
      if (info.eq.0) then
        call dsisl(k,3,3,kpvt,r)
      endif



1000  format(e13.5,e13.5,e13.5,e13.5,e13.5,e13.5)
      end


      subroutine compare_triangle_sav (n,tria,fp,ntet,idtet,xl,N0,
     1                        ns, res, fp_new)
c__________________________________________________TCG. 19.11.2002___71
c
c
c  compare shared tetsides
c
c  n     No. of triangles to be checked
c  res   sides in sorted order
c  ftet  tetraeder face
c  new_flag old flag --> .false.  tetside is at the old front
c                        .true.   tetside is at the new front
c  fp         frontpoint
c
c____________________________________________________________________71
      implicit none

	include 'tdata.h'
      include 'iofile.h'

      logical new_flag, n1, n2, n3
     
      integer n, ns,ii,jj
	integer tria(4,3), ntet(4), idtet(4), res(4,6)
	integer ftet(4,6), np 
	integer i,j,k,l
	real*8 fp(3), xl(3,4), xtria(4,3,3), m(3), mu, temp
	real*8 N0(3), eps, fp_new(4,3), di(3), delta
	real*8 min
	real*8 x(3), y(3), z(3)

	data eps/1.0d-8/
	data delta/1.0d-4/

c generate tet faces
      
      ftet(1,1) = ntet(1)
      ftet(1,2) = ntet(2)
      ftet(1,3) = ntet(3)
      ftet(1,4) = idtet(1)
      ftet(1,5) = idtet(2)
      ftet(1,6) = idtet(3)
	do i= 1,3
	  xtria(1,i,1) = xl(i,1)
	  xtria(1,i,2) = xl(i,2)
	  xtria(1,i,3) = xl(i,3)
	enddo


      ftet(2,1) = ntet(1)
      ftet(2,2) = ntet(2)
      ftet(2,3) = ntet(4)
      ftet(2,4) = idtet(1)
      ftet(2,5) = idtet(2)
      ftet(2,6) = idtet(4)
	do i= 1,3
	  xtria(2,i,1) = xl(i,1)
	  xtria(2,i,2) = xl(i,2)
	  xtria(2,i,3) = xl(i,4)
	enddo

      ftet(3,1) = ntet(2)
      ftet(3,2) = ntet(3)
      ftet(3,3) = ntet(4)
      ftet(3,4) = idtet(2)
      ftet(3,5) = idtet(3)
      ftet(3,6) = idtet(4)
	do i= 1,3
	  xtria(3,i,1) = xl(i,2)
	  xtria(3,i,2) = xl(i,3)
	  xtria(3,i,3) = xl(i,4)
	enddo

      ftet(4,1) = ntet(1)
      ftet(4,2) = ntet(3)
      ftet(4,3) = ntet(4)
      ftet(4,4) = idtet(1)
      ftet(4,5) = idtet(3)
      ftet(4,6) = idtet(4)
	do i= 1,3
	  xtria(4,i,1) = xl(i,1)
	  xtria(4,i,2) = xl(i,3)
	  xtria(4,i,3) = xl(i,4)
	enddo

      res=0
      ns=0
	fp_new = 0.0d0

      do j= 1,4
c initialize
       new_flag=.true.
c check differnet id
	 If(Abs(ftet(j,4)+ftet(j,5)+ftet(j,6)).ne.3) then ! ftet is crossed by disc
	   do i= 1,n
	     n1 = (tria(i,1).eq.ftet(j,1)).or.(tria(i,1).eq.ftet(j,2))
	1			 .or.(tria(i,1).eq.ftet(j,3))
	     n2 = (tria(i,2).eq.ftet(j,1)).or.(tria(i,2).eq.ftet(j,2))
     2        	 .or.(tria(i,2).eq.ftet(j,3))
	     n3 = (tria(i,3).eq.ftet(j,1)).or.(tria(i,3).eq.ftet(j,2))
     3    		 .or.(tria(i,3).eq.ftet(j,3))
c check coincidence with old triangles -->true
	     IF(n1.and.n2.and.n3) then
		  new_flag=.false.
	     endif
	   enddo
	   if (new_flag) then
	     ns = ns+1
	     do k= 1,3
	       res(ns,k) = ftet(j,k)
	     enddo
c nodes are on the frontline -> multiply by 2
	     do k= 4,6
c	       res(ns,k) = 2*ttim*ftet(j,k)
	       res(ns,k) = ftet(j,k)
	     enddo

c compute new frontpoint
           np=0
           do k= 1,3  !lines
c compute direction vectors of lines
            do l= 1,3   
	       if(k.eq.1) then  !line 1
              m(l) = xtria(j,l,2) - xtria(j,l,1) 
	       elseif(k.eq.2) then  !line 2
              m(l) = xtria(j,l,3) - xtria(j,l,2) 	       
	       elseif(k.eq.3) then  !line 3
              m(l) = xtria(j,l,1) - xtria(j,l,3) 
	       endif
	      enddo

c compute crosspoint
            temp = m(1)*N0(1) + m(2)*N0(2) + m(3)*N0(3) 
	      If (Abs(temp).gt.eps) then ! plane is not parallel to m
	        mu = (N0(1)*(fp(1)-xtria(j,1,k)) + 
	1              N0(2)*(fp(2)-xtria(j,2,k)) + 
	1              N0(3)*(fp(3)-xtria(j,3,k)))/temp  
	        if ((mu.gt.0.0d0).and.(mu.lt.1.0d0)) then
	         np = np + 1
	         do l= 1,3
	           fp_new(ns,l) = fp_new(ns,l) + 
     1                          0.5d0*(xtria(j,l,k) + mu*m(l))
	         enddo
	        endif
	      endif
	     enddo
	     if(np.ne.2) then ! no intersection -> find closest point 
	       min = 1.0d16
	       do k= 1,3 ! nodes
	         temp = abs((xtria(j,1,k)-fp(1))*N0(1)+
     1					(xtria(j,2,k)-fp(2))*N0(2)+
     2					(xtria(j,3,k)-fp(3))*N0(3))
	         if(temp.lt.min) then
	           min = temp
c  and move point a little into interior
                 if (k.eq.1) then
	             do l= 1,3
	               di(l) = 0.5d0*(xtria(j,l,2) + xtria(j,l,3)) 
	1								- xtria(j,l,1)
	             enddo	             
	           elseif (k.eq.2) then
	             do l= 1,3
	               di(l) = 0.5d0*(xtria(j,l,1) + xtria(j,l,3)) 
	1								- xtria(j,l,2)
	             enddo	             
	           elseif (k.eq.3) then
	             do l= 1,3
	               di(l) = 0.5d0*(xtria(j,l,1) + xtria(j,l,2)) 
	1								- xtria(j,l,3)
	             enddo	             
	           endif
	           do l= 1,3
	             fp_new(ns,l) = xtria(j,l,k) + delta*di(l)
	           enddo
	         endif
	       enddo
           endif
c  check result
           if (np.ne.2) then
	       write(*,*)   'WARNING --> Closest point approximation'
	       write(iow,*)   'WARNING --> Closest point approximation'
	     endif
           do ii= 1,3
	       x(ii) = xtria(j,ii,3) - xtria(j,ii,1)
	       y(ii) = xtria(j,ii,2) - xtria(j,ii,1)
	       z(ii) = fp_new(ns,ii) - xtria(j,ii,1)
		 enddo 
		 temp = (-x(3)*y(2) + x(2)*y(3))*z(1) + 
     1		     (x(3)*y(1) - x(1)*y(3))*z(2) + 
     2            (-x(2)*y(1) + x(1)*y(2))*z(3) 
	     if(abs(temp).gt.eps) then
	       write(*,*)     'ERROR --> fp not on facet'
	       write(iow,*)   'ERROR --> fp not on facet'
             write(iow,1000) ((xl(ii,jj),ii=1,3),jj=1,4)
             write(iow,1001) (fp_new(ns,ii),ii=1,3),(fp(ii),ii=1,3), 
     1						(N0(ii),ii=1,3)
	     endif
	   endif
	 endif

	enddo

c define fp_new
	do i= 1, ns
         if ((fp_new(i,1).eq.0.0d0).and.(fp_new(i,2).eq.0.0d0).and.
     1       (fp_new(i,3).eq.0.0d0)) then
          do l= 1,3
           fp_new(i,l) = fp(l)
          enddo
         endif
      enddo
1000  format(
     1  10x,'Tet data   '/
     2  10x,'node 1      ',e12.5,e12.5,e12.5/
     2  10x,'node 2      ',e12.5,e12.5,e12.5/
     2  10x,'node 3      ',e12.5,e12.5,e12.5/
     2  10x,'node 4      ',e12.5,e12.5,e12.5/)
1001  format(
     2  10x,'Plane data              '/
     2  10x,'node_old    ',e12.5,e12.5,e12.5/
     2  10x,'node_new    ',e12.5,e12.5,e12.5/
     2  10x,'normal      ',e12.5,e12.5,e12.5/)
      end


	subroutine auto_cini_11 (elem,N0,curmat)
c__________________________________________________TCG. 19.09.2005___71
c
c
c  automatic crackt tip initialisation
c
c   elem         elemt no.
c   N0           Normal to failure surface
c   sol          node number of tip facet
c   xb           tip point coordinates
c   curmat       current material
c____________________________________________________________________71
      implicit none

      include  'iofile.h'
	include  'comblk.h'
	include  'sdata.h'
      include  'cTip.h'
      
	integer elem,i,j,k,sol(3),ntria(3),nodenr(4), curmat

	real*8 N0(3), xl(ndm,4), a(3), b(3), in, xb(3), xtemp(3)
	real*8 nf(3), min

      logical Xflag
	integer Xpoint, Xlen, Xpre
      logical IXflag
	integer IXpoint, IXlen, IXpre


c write message
         write(*,*)   'Automatic crack initialisation in Elmt.',elem
         write(iow,*) 'Automatic crack initialisation in Elmt.',elem


c  get pointer and length of the X     
	 call pgetd ('X ',xpoint,xlen,xpre,xflag)
c  get pointer and length of the IX     
	 call pgetd ('IX ',ixpoint,ixlen,ixpre,ixflag)

c retrieve node no. and lagrangian coordinates
      do j = 1,4
	  nodenr(j) = mr(ixpoint + (elem-1)*nen1 + j -1)
	  do k= 1,ndm
	    xl(k,j) = hr(xpoint+(nodenr(j)-1)*ndm + k - 1) 
	  enddo
      enddo

      min = 1.0d20
	do i= 1,4
c compte normal to facet
        if (i.eq.1) then
          do j= 1,3
	      a(j) = xl(j,2) - xl(j,1)
	      b(j) = xl(j,3) - xl(j,1)
	      xtemp(j) = 0.33333d0*(xl(j,1)+xl(j,2)+xl(j,3))
	    enddo
	    ntria(1) = nodenr(1)
	    ntria(2) = nodenr(2)
	    ntria(3) = nodenr(3)
	  elseif(i.eq.2) then
          do j= 1,3
	      a(j) = xl(j,2) - xl(j,1)
	      b(j) = xl(j,4) - xl(j,1)
	      xtemp(j) = 0.33333d0*(xl(j,1)+xl(j,2)+xl(j,4))
	    enddo
	    ntria(1) = nodenr(1)
	    ntria(2) = nodenr(2)
	    ntria(3) = nodenr(4)
	  elseif(i.eq.3) then
          do j= 1,3
	      a(j) = xl(j,3) - xl(j,1)
	      b(j) = xl(j,4) - xl(j,1)
	      xtemp(j) = 0.33333d0*(xl(j,1)+xl(j,3)+xl(j,4))
	    enddo
	    ntria(1) = nodenr(1)
	    ntria(2) = nodenr(3)
	    ntria(3) = nodenr(4)
	  elseif(i.eq.4) then
          do j= 1,3
	      a(j) = xl(j,3) - xl(j,2)
	      b(j) = xl(j,4) - xl(j,2)
	      xtemp(j) = 0.33333d0*(xl(j,2)+xl(j,3)+xl(j,4))
	    enddo
	    ntria(1) = nodenr(2)
	    ntria(2) = nodenr(3)
	    ntria(3) = nodenr(4)
	  endif
        Nf(1) = -a(3)*b(2) + a(2)*b(3)
	  Nf(2) =  a(3)*b(1) - a(1)*b(3)
	  Nf(3) = -a(2)*b(1) + a(1)*b(2)
c compute inner product
        !in = Nf(1)*N0(1)+Nf(2)*N0(2)+Nf(3)*N0(3)
        in = abs(Nf(1)*N0(1)+Nf(2)*N0(2)+Nf(3)*N0(3))
	  if (in.lt.min) then
	    sol = ntria
	    xb  = xtemp
	    min = in
	  endif

      enddo
c write solution into ctip data

	ncTip = ncTip + 1
	mTi(ncTip)=curmat
	do k=1,3
        cTip(ncTip,k) = sol(k) 
        cTip(ncTip,k+3) = 0 
	  cTip_point(ncTip,k) = xb(k)
	enddo
	end




      subroutine debug_history_11(hist, hlen)
      implicit none
      
      integer hlen, i
      real*8 hist(hlen)
      
      do i= 1,hlen
      if(abs(hist(i)-3.530460559031873D-004).lt.1.0d-8) then
            write(*,*) 'hit'
      endif
      enddo
     
      end