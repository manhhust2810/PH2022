	program CNA
c	This common neighbor analysis (CNA) program to dectect crystalline structure
c	The program was writen by:
c	                 Le Van Vinh, PhD
c	        Department of Computational Physics
c	      Hanoi University of Science & Technology
	integer i,j,na,n
	integer cp(9999)
	real x2(9999),y2(9999),z2(9999)
	real r,rx,ry,rz,lx,ly,lz
	
	call readfile(na,lx,ly,lz,x2,y2,z2,cp)
	write(*,*)'na=',na
	write(*,*)'lx=',lx,'  ly=',ly,'  lz=',lz
	write(*,*)x2(na),'	',y2(na),'	',z2(na),'	',cp(na)
	lx=lx*0.529177
	ly=ly*0.529177
	lz=lz*0.529177
	do i=1,na
	  x2(i)=x2(i)*0.529177
	  y2(i)=y2(i)*0.529177
	  z2(i)=z2(i)*0.529177
	enddo
	call dcluster(na,lx,ly,lz,x2,y2,z2,cp)
	write(*,*)'finished!'
	end
c	***********************************
c	read an input file of coordinations
	subroutine readfile(na,lx,ly,lz,x2,y2,z2,cp)
	integer i,na,n
	integer cp(9999)
	real lx,ly,lz
	real x2(9999),y2(9999),z2(9999)
c	the number of atoms in the model: na
	open(1,file='al2o3.dat')
	read(1,*)lx,ly,lz
	read(1,*)na
	do i=1,na
	  read(1,*)n,x2(i),y2(i),z2(i),cp(i)
	enddo
	endfile(1)
	end	
c	***********************************
c	the periodic boundary condition!
	real function pbc(l,xa,xb)
	real l,l2,xa,xb,r
	l2=l/2.0
	r=xb-xa
	if (r.ge.l2) then
	  r=-l+r
	elseif (r.le.-l2) then
	  r=l+r
	endif
	pbc=r
	return
	end
c	***********************************
	subroutine dcluster(na,lx,ly,lz,x2,y2,z2,cp)
	integer i,j,j1,j2,k,k1,k2,na,n,n1,n2,nal,ncr,non,ntol,nf
	integer ic12,ic11,ic10,ic9,ic8,ic7,ic6,ic5,ic4,ic3,ic2,ic1
	integer n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16
	integer o4,o5,o6,o45,o46,o56,o456,cpf(9999),cph(9999)
	integer cp(9999),neib(9999),id(10),h421(9999),h422(9999),cpol2(9999)
	integer h442(9999),h662(9999),cp1(9999),pol(9999),cp2(9999)
	integer fnei(9999,20),h(9999,20),l(9999,20),m(9999,20),cl(9999)
	integer pol1(9999),pol2(9999),cp3(9999),id1(9999,2),h552(9999)
	integer h332(9999),cr(9999),h555(9999),h544(9999),h433(9999)
	integer neibO(9999),neibM(9999),fneiM(9999,20),fneiO(9999,20)
	integer nm1(400),nm2(400),nfcc(9999),cpO(9999),nmi(9999)
	real lx,ly,lz,rc,rx,ry,rz,rij,a0,lv,x,y,z,xi,yi,zi,xa,ya,za,r1,r2
	real x1,y1,z1,rcM
	real x2(9999),y2(9999),z2(9999)

	rc=3.63
	rcM=2.35
	do i=1,na
	  neib(i)=0
	  neibO(i)=0
	  neibM(i)=0
	  h421(i)=0
	  h422(i)=0
	  h442(i)=0
	  h662(i)=0
	  h552(i)=0
	  h332(i)=0
	  h555(i)=0
	  h544(i)=0
	  h433(i)=0
	  cp1(i)=0
	  cp2(i)=0
	  cp3(i)=0
	  pol(i)=0
	  pol1(i)=0
	  pol2(i)=0
	  cpol2(i)=0
	  cr(i)=0
	  nfcc(i)=0
	  cpO(i)=0
	  neibM(i)=0
	  neibO(i)=0
	  nmi(i)=0
	  cpf(i)=0
	  cph(i)=0
	enddo
	nal=0
	do i=1,na
	  if(cp(i).eq.1)then
	    nal=nal+1
	  endif
	enddo
	write(*,*)'nAl=',nal
c	--------------
c	calculate coodinate number of Metal and Oxygen:
c	for Metal:
	do i=1,nal
	  n=0
	  do j=nal+1,na
	    rx=pbc(lx,x2(j),x2(i))
	    ry=pbc(ly,y2(j),y2(i))
	    rz=pbc(lz,z2(j),z2(i))
	    rij=sqrt(rx*rx+ry*ry+rz*rz)
	    if (rij.lt.rcM) then
	          neibM(i)=neibM(i)+1
	          n=n+1
	          fneiM(i,n)=j
		  endif
	  enddo
	enddo
c	for Oxyen:
	do i=nal+1,na
	  n=0
	  do j=1,nal
	    rx=pbc(lx,x2(j),x2(i))
	    ry=pbc(ly,y2(j),y2(i))
	    rz=pbc(lz,z2(j),z2(i))
	    rij=sqrt(rx*rx+ry*ry+rz*rz)
	    if (rij.lt.rcM) then
	          neibO(i)=neibO(i)+1
	          n=n+1
	          fneiO(i,n)=j
		  endif
	  enddo
	enddo

c	--------------
c	Calculate CNA:
	do i=nal+1,na
	  n=0
	  neib(i)=0
	  do j=nal+1,na
	    if (i.ne.j) then
	      rx=pbc(lx,x2(j),x2(i))
	      ry=pbc(ly,y2(j),y2(i))
	      rz=pbc(lz,z2(j),z2(i))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
	      if (rij.lt.rc) then
	          neib(i)=neib(i)+1
	          n=n+1
	          fnei(i,n)=j
		  endif
	    endif
	  enddo
	enddo
c	-----------	
	do i=nal+1,na
	  do k=1,neib(i)
	    j=fnei(i,k)
	    n=0
	    do k1=1,neib(j)
	      j1=fnei(j,k1)
		  do k2=1,neib(i)
	        j2=fnei(i,k2)
	        if(j1.eq.j2)then
	          n=n+1
	          id(n)=j1
			endif
		  enddo
		enddo
c	    ----------
	    h(i,k)=n
c	    ----------
	    n1=0
		do k1=1,n-1
	      do k2=k1+1,n
	        j1=id(k1)
	        j2=id(k2)
	        rx=pbc(lx,x2(j2),x2(j1))
	        ry=pbc(ly,y2(j2),y2(j1))
	        rz=pbc(lz,z2(j2),z2(j1))
	        rij=sqrt(rx*rx+ry*ry+rz*rz)
		    if(rij.lt.rc)then
	          n1=n1+1
			endif
		  enddo
		enddo
	    l(i,k)=n1
c	    -----------
		n2=0
		do k1=1,n
	      n1=0
		  do k2=1,n
	        j1=id(k1)
	        j2=id(k2)
	        if(j1.ne.j2)then
	          rx=pbc(lx,x2(j2),x2(j1))
	          ry=pbc(ly,y2(j2),y2(j1))
	          rz=pbc(lz,z2(j2),z2(j1))
	          rij=sqrt(rx*rx+ry*ry+rz*rz)
		      if(rij.lt.rc)then
	            n1=n1+1
			  endif
		    endif
		  enddo
		  if(n2.lt.n1)then
	        n2=n1
		  endif
		enddo	  
		m(i,k)=n2  
c	    ----------
	  enddo
	enddo
c	------------
	do i=nal+1,na
	  do j=1,neib(i)
	    if(h(i,j).eq.4)then
	      if(l(i,j).eq.2)then
	        if(m(i,j).eq.1)then
	          h421(i)=h421(i)+1
			endif
	        if(m(i,j).eq.2)then
	          h422(i)=h422(i)+1
			endif
		  endif
		  if(l(i,j).eq.4)then
	        if(m(i,j).eq.2)then
			  h442(i)=h442(i)+1
			endif	      
		  endif
		elseif(h(i,j).eq.6)then
	      if(l(i,j).eq.6)then
	        if(m(i,j).eq.2)then
			  h662(i)=h662(i)+1
			endif
		  endif
		elseif(h(i,j).eq.5)then
	      if(l(i,j).eq.5)then
	        if(m(i,j).eq.2)then
	          h552(i)=h552(i)+1
			endif
		  endif
		elseif(h(i,j).eq.3)then
	      if(l(i,j).eq.3)then
	        if(m(i,j).eq.2)then
	          h332(i)=h332(i)+1
			endif
		  endif
		elseif(h(i,j).eq.5)then
	      if(l(i,j).eq.5)then
	        if(m(i,j).eq.5)then
	          h555(i)=h555(i)+1
			endif
		  endif
		elseif(h(i,j).eq.5)then
	      if(l(i,j).eq.4)then
	        if(m(i,j).eq.4)then
	          h544(i)=h544(i)+1
			endif
		  endif
		elseif(h(i,j).eq.4)then
	      if(l(i,j).eq.3)then
	        if(m(i,j).eq.3)then
	          h433(i)=h433(i)+1
			endif
		  endif
		endif
	  enddo
	enddo
c	------------
	open(2,file='Icohedral.txt')
	  write(2,*)'i:	h421:	h552:	h422:	h332:	h433:'
	  do i=nal+1,na
	    write(2,*)i,'	',h421(i),'	',h552(i),'	',h422(i),'	',h332(i),
     +'	',h433(i)
	  enddo
	endfile(2)
c	------------
	ic1=0
	ic2=0
	ic3=0
	ic4=0
	ic5=0
	ic6=0
	ic7=0
	ic8=0
	ic9=0
	ic10=0
	ic11=0
	ic12=0
	do i=nal+1,na
	select case(h552(i))
	  case(12)
	    ic12=ic12+1
	  case(11)
	    ic11=ic11+1
	  case(10)
	    ic10=ic10+1
	  case(9)
	    ic9=ic9+1
	  case(8)
	    ic8=ic8+1
	  case(7)
	    ic7=ic7+1
	  case(6)
	    ic6=ic6+1
	  case(5)
	    ic5=ic5+1
	  case(4)
	    ic4=ic4+1
	  case(3)
	    ic3=ic3+1
	  case(2)
	    ic2=ic2+1
	  case(1)
	    ic1=ic1+1
	endselect
	enddo

	n=na-nal

	open(2,file='Icodist.txt')
	write(2,*)'Distribution of h552:'
	write(2,*)'ico12:','	',ic12,'	',ic12*1./n
	write(2,*)'ico11:','	',ic11,'	',ic11*1./n
	write(2,*)'ico10:','	',ic10,'	',ic10*1./n
	write(2,*)'ico9:','	',ic9,'	',ic9*1./n
	write(2,*)'ico8:','	',ic8,'	',ic8*1./n
	write(2,*)'ico7:','	',ic7,'	',ic7*1./n
	write(2,*)'ico6:','	',ic6,'	',ic6*1./n
	write(2,*)'ico5:','	',ic5,'	',ic5*1./n
	write(2,*)'ico4:','	',ic4,'	',ic4*1./n
	write(2,*)'ico3:','	',ic3,'	',ic3*1./n
	write(2,*)'ico2:','	',ic2,'	',ic2*1./n
	write(2,*)'ico1:','	',ic1,'	',ic1*1./n

	ic1=0
	ic2=0
	ic3=0
	ic4=0
	ic5=0
	ic6=0
	ic7=0
	ic8=0
	ic9=0
	ic10=0
	ic11=0
	ic12=0
	do i=nal+1,na
	select case(h332(i))
	  case(12)
	    ic12=ic12+1
	  case(11)
	    ic11=ic11+1
	  case(10)
	    ic10=ic10+1
	  case(9)
	    ic9=ic9+1
	  case(8)
	    ic8=ic8+1
	  case(7)
	    ic7=ic7+1
	  case(6)
	    ic6=ic6+1
	  case(5)
	    ic5=ic5+1
	  case(4)
	    ic4=ic4+1
	  case(3)
	    ic3=ic3+1
	  case(2)
	    ic2=ic2+1
	  case(1)
	    ic1=ic1+1
	endselect
	enddo
	write(2,*)
	write(2,*)'Distribution of h332:'
	write(2,*)'ico12:','	',ic12,'	',ic12*1./n
	write(2,*)'ico11:','	',ic11,'	',ic11*1./n
	write(2,*)'ico10:','	',ic10,'	',ic10*1./n
	write(2,*)'ico9:','	',ic9,'	',ic9*1./n
	write(2,*)'ico8:','	',ic8,'	',ic8*1./n
	write(2,*)'ico7:','	',ic7,'	',ic7*1./n
	write(2,*)'ico6:','	',ic6,'	',ic6*1./n
	write(2,*)'ico5:','	',ic5,'	',ic5*1./n
	write(2,*)'ico4:','	',ic4,'	',ic4*1./n
	write(2,*)'ico3:','	',ic3,'	',ic3*1./n
	write(2,*)'ico2:','	',ic2,'	',ic2*1./n
	write(2,*)'ico1:','	',ic1,'	',ic1*1./n
	write(2,*)
	write(2,*)'i:	h552:	h332:	neib:'
	do i=nal+1,na
	  write(2,*)i,'	',h552(i),'	',h332(i),'	',neib(i)
	enddo
	endfile(2)
c	------------
c	calculate number of crystalline atoms: 
c	cp1=1 => core-fcc
c	cp1=2 => shell- fcc
c	cp1=3 => core-hcp
c	cp1=4 => shell-hcp
c	cp1=5 => core-bcc
c	cp1=6 => shell-bcc
c	cp1=7 => core-ico
c	cp1=8 => shell-ico
	do i=nal+1,na
	  if(h421(i).eq.12)then
	    cp1(i)=1
	    cpf(i)=1
	    do k=1,neib(i)
	      j=fnei(i,k)
		  cpf(j)=2
		  if((cp1(j).ne.1).and.(cp1(j).ne.3))then
		    cp1(j)=2
		  endif
		enddo
	  endif
	  if(h421(i).eq.6)then
	    if(h422(i).eq.6)then
	      cp1(i)=3
	      cph(i)=3
		  do k=1,neib(i)
	        j=fnei(i,k)
	        cph(j)=4
			if((cp1(j).ne.1).and.(cp1(j).ne.3).and.(cp1(j).ne.2))then
		      cp1(j)=4
		    endif
		  enddo
		endif
	  endif
	  if(h442(i).eq.6)then
	    if(h662(i).eq.8)then
	      cp1(i)=5
	      do k=1,neib(i)
	        j=fnei(i,k)
	        if((cp1(j).ne.1).and.(cp1(j).ne.3).and.(cp1(j).ne.5))then
			  cp1(j)=6
	        endif
		  enddo
		endif
	  endif
	  if(h552(i).eq.12)then
	    cp1(i)=7
	    do k=1,neib(i)
	      j=fnei(i,k)
	      if((cp1(j).ne.1).and.(cp1(j).ne.3).and.(cp1(j).ne.6))then
	        if((cp1(j).ne.2).and.(cp1(j).ne.4))then
		      cp1(j)=8
	        endif
		  endif
		enddo
	  endif
	enddo
c	-----------
	k1=0
	k2=0
	j1=0
	j2=0
	k=0
	j=0
	n=0
	n1=0
	do i=nal+1,na
	  if(cp1(i).eq.1)then
	    k1=k1+1
	  elseif(cp1(i).eq.2)then
	    k2=k2+1
	  elseif(cp1(i).eq.3)then
	    j1=j1+1
	  elseif(cp1(i).eq.4)then
	    j2=j2+1
	  elseif(cp1(i).eq.5)then
	    k=k+1
	  elseif(cp1(i).eq.6)then
	    j=j+1
	  elseif(cp1(i).eq.7)then
	    n=n+1
	  elseif(cp1(i).eq.8)then
	    n1=n1+1
	  endif
	enddo
c	  --------
	open(2,file='numcrystatoms.txt')
	  write(2,*)'number of core fcc atoms:',k1
	  write(2,*)' '
	  write(2,*)'number of shell fcc atoms:',k2
	  write(2,*)' '
	  write(2,*)'number of core hcp atoms:',j1
	  write(2,*)' '
	  write(2,*)'number of shell hcp atoms:',j2
	  write(2,*)' '
	  write(2,*)'number of core bcc atoms:',k
	  write(2,*)' '
	  write(2,*)'number of shell bcc atoms:',j
	  write(2,*)' '
	  write(2,*)'number of core ICO atoms:',n
	  write(2,*)' '
	  write(2,*)'number of shell ICO atoms:',n1
	  write(2,*)' '
	  write(2,*)'number of crystalline atoms:',k1+k2+j1+j2+k+j
	  write(2,*)' '
	  write(2,*)'number of amorphous atoms:',na-(k1+k2+j1+j2+k+j+n+n1)
	  write(2,*)'i:	cp1:	neib:'
	  do i=1,na
	    write(2,*)i,'	',cp1(i),'	',neib(i)
	  enddo
	endfile(2)
c	------------
	open(2,file='CrysOxy.txt')
	write(2,*)'i:	cp1:	neibO: n4:	n5:	n6:'
	do i=nal+1,na
	  if(cp1(i).eq.1)then
	    n4=0
	    n5=0
	    n6=0
		do k=1,neibO(i)
	      j=fneiO(i,k)
	      select case(neibM(j))
	      case(4)
	        n4=n4+1
	      case(5)
	        n5=n5+1
	      case(6)
	        n6=n6+1
		  endselect
		enddo
	  write(2,*)i,'	',cp1(i),'	',neibO(i),'	',n4,'	',n5,'	',n6 
	  endif
	enddo
	do i=nal+1,na
	  if(cp1(i).eq.2)then
	    n4=0
	    n5=0
	    n6=0
		do k=1,neibO(i)
	      j=fneiO(i,k)
	      select case(neibM(j))
	      case(4)
	        n4=n4+1
	      case(5)
	        n5=n5+1
	      case(6)
	        n6=n6+1
		  endselect
		enddo
	  write(2,*)i,'	',cp1(i),'	',neibO(i),'	',n4,'	',n5,'	',n6 
	  endif
	enddo
	do i=nal+1,na
	  if(cp1(i).eq.3)then
	    n4=0
	    n5=0
	    n6=0
		do k=1,neibO(i)
	      j=fneiO(i,k)
	      select case(neibM(j))
	      case(4)
	        n4=n4+1
	      case(5)
	        n5=n5+1
	      case(6)
	        n6=n6+1
		  endselect
		enddo
	  write(2,*)i,'	',cp1(i),'	',neibO(i),'	',n4,'	',n5,'	',n6 
	  endif
	enddo
	do i=nal+1,na
	  if(cp1(i).eq.4)then
	    n4=0
	    n5=0
	    n6=0
		do k=1,neibO(i)
	      j=fneiO(i,k)
	      select case(neibM(j))
	      case(4)
	        n4=n4+1
	      case(5)
	        n5=n5+1
	      case(6)
	        n6=n6+1
		  endselect
		enddo
	  write(2,*)i,'	',cp1(i),'	',neibO(i),'	',n4,'	',n5,'	',n6 
	  endif
	enddo
	do i=nal+1,na
	  if(cp1(i).eq.7)then
	    n4=0
	    n5=0
	    n6=0
		do k=1,neibO(i)
	      j=fneiO(i,k)
	      select case(neibM(j))
	      case(4)
	        n4=n4+1
	      case(5)
	        n5=n5+1
	      case(6)
	        n6=n6+1
		  endselect
		enddo
	  write(2,*)i,'	',cp1(i),'	',neibO(i),'	',n4,'	',n5,'	',n6 
	  endif
	enddo
	endfile(2)
c	------------
c	write a file of poly-crystalline atoms:
	do i=1,na
	  if(cpf(i).eq.1)then
	    cp2(i)=1
	    do k=1,neibO(i)
	      j=fneiO(i,k)
		  cp2(j)=5
		enddo
	  elseif(cpf(i).eq.2)then
	    cp2(i)=2
	    do k=1,neibO(i)
	      j=fneiO(i,k)
		  cp2(j)=5
		enddo
	  endif
	  if(cph(i).eq.3)then
	    cp2(i)=3
	    do k=1,neibO(i)
	      j=fneiO(i,k)
		  cp2(j)=6
		enddo
	  elseif(cph(i).eq.4)then
	    cp2(i)=4
	    do k=1,neibO(i)
	      j=fneiO(i,k)
		  cp2(j)=6
		enddo
	  endif
	enddo
c	---
	open(2,file='xyzCrystalMO.txt')
	  write(2,*)'lx=',lx
	  k=0
	  do i=1,na
	    if(cp2(i).ne.0)then
		  k=k+1
		  write(2,*)k,'	',x2(i),'	',y2(i),'	',z2(i),'	',cp2(i)
	    endif
	  enddo
	endfile(2)
c	----
	open(2,file='xyzCrystalO.txt')
	  write(2,*)'lx=',lx
	  k=0
	  do i=1,na
	    if(cp1(i).ne.0)then
		  k=k+1
		  write(2,*)k,'	',x2(i),'	',y2(i),'	',z2(i),'	',cp1(i)
	    endif
	  enddo
	endfile(2)
c	-------------	
	open(2,file='IntersticeM.txt')
	n4=0
	n5=0
	n6=0
	do i=1,nal
	  if(cp2(i).eq.5)then
	    selectcase(neibM(i))
	    case(4)
	      n4=n4+1
	    case(5)
	      n5=n5+1
	    case(6)
	      n6=n6+1
	    case(7)
	      n6=n6+1
		endselect
	  endif
	enddo	
	  write(2,*)'fcc:'
	  write(2,*)'ZAl-O=4:	',n4
	  write(2,*)'ZAl-O=5:	',n5
	  write(2,*)'ZAl-O=6:	',n6
c	hcp:
	n4=0
	n5=0
	n6=0
	do i=1,nal
	  if(cp2(i).eq.6)then
	    selectcase(neibM(i))
	    case(4)
	      n4=n4+1
	    case(5)
	      n5=n5+1
	    case(6)
	      n6=n6+1
	    case(7)
	      n6=n6+1
		endselect
	  endif
	enddo	
	  write(2,*)'hcp:'
	  write(2,*)'ZAl-O=4:	',n4
	  write(2,*)'ZAl-O=5:	',n5
	  write(2,*)'ZAl-O=6:	',n6
	endfile(2)
c	------------
c	------------
c	calculate fcc-hcp crystalline clusters:
	do i=1,na
	  pol2(i)=0
	enddo
	do i=1,na
	  if((cp1(i).eq.1).or.(cp1(i).eq.3))then
	    pol2(i)=i
	    do k=1,neib(i)
	      j=fnei(i,k)
	      pol2(j)=i
		enddo
	  endif
	enddo
	do j1=1,12
	 do i=1,na
	  n=na+1-i
	  if((cp1(n).eq.1).or.(cp1(n).eq.3))then
	    k1=n
		do k=1,neib(n)
	      j=fnei(n,k)
	      k2=pol2(j)
	      if(k2.le.k1)then
		    k1=k2
		  endif
	    enddo
	    pol2(n)=k1
	    do k=1,neib(n)
	      j=fnei(n,k)
	      pol2(j)=k1
		enddo
	  endif
	 enddo
	enddo
c	------------
c	write a file of polyAlO4:
	open(2,file='polyfcc_hcp.txt')
c	-----------
	do i=1,na
	  id1(i,1)=pol2(i)
	  id1(i,2)=i
	  cl(i)=0
	enddo
	do i=1,na-1
	  do j=i+1,na
	    if(id1(i,1).gt.id1(j,1))then
	      k1=id1(i,1)
	      id1(i,1)=id1(j,1)
	      id1(j,1)=k1
	      
		  k2=id1(i,2)
	      id1(i,2)=id1(j,2)
	      id1(j,2)=k2
		endif
	  enddo
	enddo
	k1=0
	do i=1,na-1
	  if(id1(i,1).eq.0)then
	    j1=0
	  endif
	  if(id1(i,1).ne.id1(i+1,1))then
	    j1=j1+1

		k1=k1+1
	    cl(k1)=i
	    write(*,*)'k1=',k1,'	',cl(k1)
	  endif
	  cpol2(id1(i+1,2))=j1
	enddo
	k1=k1+1
	cl(k1)=na
	
	  write(2,*)'pol2[',cl(1),']=',id1(1,1),'	cl=',cl(1)
	  do i=2,k1
	    k2=cl(i)-cl(i-1)
	    write(2,*)'pol2[',cl(i),']=',id1(cl(i),1),'	cl=',k2
	  enddo
	  write(2,*)'i=		pol4:		cp:	cpol2:'
c	--------------------
	endfile(2)
c	----------------
c	------------
c	write a file of poly-crystalline atoms:
	open(2,file='polyfcc_hcp.txt')
	  write(2,*)'lx=',lx
	  k=0
	  do i=1,na
	    if(pol2(i).ne.0)then
		  k=k+1
	      write(2,*)k,'	',pol2(i)
	    endif
	  enddo
	endfile(2)
c	------------
c	+++++++++++++++++
	n=0
	do i=nal+1,na
	  n4=0
	  n5=0
	  n6=0
	  do k=1,neibO(i)
	    j=fneiO(i,k)
	    select case(neibM(j))
	    case(4)
	      n4=1
	    case(5)
	      n5=1
	    case(6)
	      n6=1
		endselect
	  enddo
	  if(n4.eq.1)then
	    if(n5.eq.0)then
	      if(n6.eq.0)then
	        cpO(i)=4
		    n=n+1
		  elseif(n6.eq.1)then
	        cpO(i)=46
		    n=n+1
		  endif
		elseif(n5.eq.1)then
	      if(n6.eq.0)then
	        cpO(i)=45
		    n=n+1
		  elseif(n6.eq.1)then
	        cpO(i)=456
		    n=n+1
		  endif	     
		endif
	  elseif(n4.eq.0)then
	    if(n5.eq.0)then
	      if(n6.eq.1)then
	        cpO(i)=6
		    n=n+1
		  endif
		elseif(n5.eq.1)then
	      if(n6.eq.0)then
	        cpO(i)=5
		    n=n+1
		  elseif(n6.eq.1)then
	        cpO(i)=56
		    n=n+1
		  endif	     
		endif
	  endif
	enddo
	write(*,*)'nO=',n
c	----------
c	calculate correlation Al-O and O-O:
c	----------
	open(2,file='correlationNeiborO.txt')
c	---cpO=4
	n7=0
	n8=0
	n9=0
	n10=0
	n11=0
	n12=0
	n13=0
	n14=0
	n15=0
	n16=0
	do i=nal+1,na
	  if(cpO(i).eq.4)then
	    if(cp1(i).ne.0)then
	      n7=n7+1
	    endif
	    select case(neib(i))
	    case(8)
	      n8=n8+1
	    case(9)
	      n9=n9+1
	    case(10)
	      n10=n10+1
	    case(11)
	      n11=n11+1
	    case(12)
	      n12=n12+1
	    case(13)
	      n13=n13+1
	    case(14)
	      n14=n14+1
	    case(15)
	      n15=n15+1
	    case(16)
	      n16=n16+1
	    end select
	  endif
	enddo
	n=n8+n9+n10+n11+n12+n13+n14+n15+n16
	n1=n8*8+n9*9+n10*10+n11*11+n12*12+n13*13+n14*14+n15*15+n16*16
c	  write to file:
	  write(2,*)'ZAl-O=4:'
	  write(2,*)'ncrys:','	',n7,'	',n7*1./n
	  write(2,*)8,'	','n8=',n8,'	',n8*1./n
	  write(2,*)9,'	','n9=',n9,'	',n9*1./n
	  write(2,*)10,'	','n10=',n10,'	',n10*1./n
	  write(2,*)11,'	','n11=',n11,'	',n11*1./n
	  write(2,*)12,'	','n12=',n12,'	',n12*1./n
	  write(2,*)13,'	','n13=',n13,'	',n13*1./n
	  write(2,*)14,'	','n14=',n14,'	',n14*1./n
	  write(2,*)15,'	','n15=',n15,'	',n15*1./n
	  write(2,*)16,'	','n16=',n16,'	',n16*1./n
	  write(2,*)'nO-O=','	',n
	  write(2,*)'Zo-o=','	',n1*1./n
c	------cpO=5
	n7=0
	n8=0
	n9=0
	n10=0
	n11=0
	n12=0
	n13=0
	n14=0
	n15=0
	n16=0
	do i=nal+1,na
	  if(cpO(i).eq.5)then
	    if(cp1(i).ne.0)then
	      n7=n7+1
	    endif
	    select case(neib(i))
	    case(8)
	      n8=n8+1
	    case(9)
	      n9=n9+1
	    case(10)
	      n10=n10+1
	    case(11)
	      n11=n11+1
	    case(12)
	      n12=n12+1
	    case(13)
	      n13=n13+1
	    case(14)
	      n14=n14+1
	    case(15)
	      n15=n15+1
	    case(16)
	      n16=n16+1
	    end select
	  endif
	enddo
	n=n8+n9+n10+n11+n12+n13+n14+n15+n16
	n1=n8*8+n9*9+n10*10+n11*11+n12*12+n13*13+n14*14+n15*15+n16*16
c	  write to file:
	  write(2,*)'ZAl-O=5:'
	  write(2,*)'ncrys:','	',n7,'	',n7*1./n
	  write(2,*)8,'	','n8=',n8,'	',n8*1./n
	  write(2,*)9,'	','n9=',n9,'	',n9*1./n
	  write(2,*)10,'	','n10=',n10,'	',n10*1./n
	  write(2,*)11,'	','n11=',n11,'	',n11*1./n
	  write(2,*)12,'	','n12=',n12,'	',n12*1./n
	  write(2,*)13,'	','n13=',n13,'	',n13*1./n
	  write(2,*)14,'	','n14=',n14,'	',n14*1./n
	  write(2,*)15,'	','n15=',n15,'	',n15*1./n
	  write(2,*)16,'	','n16=',n16,'	',n16*1./n
	  write(2,*)'nO-O=','	',n
	  write(2,*)'Zo-o=','	',n1*1./n

c	----cpO=6
	n7=0
	n8=0
	n9=0
	n10=0
	n11=0
	n12=0
	n13=0
	n14=0
	n15=0
	n16=0
	do i=nal+1,na
	  if(cpO(i).eq.6)then
	    if(cp1(i).ne.0)then
	      n7=n7+1
	    endif
	    select case(neib(i))
	    case(8)
	      n8=n8+1
	    case(9)
	      n9=n9+1
	    case(10)
	      n10=n10+1
	    case(11)
	      n11=n11+1
	    case(12)
	      n12=n12+1
	    case(13)
	      n13=n13+1
	    case(14)
	      n14=n14+1
	    case(15)
	      n15=n15+1
	    case(16)
	      n16=n16+1
	    end select
	  endif
	enddo
	n=n8+n9+n10+n11+n12+n13+n14+n15+n16
	n1=n8*8+n9*9+n10*10+n11*11+n12*12+n13*13+n14*14+n15*15+n16*16
c	  write to file:
	  write(2,*)'ZAl-O=6:'
	  write(2,*)'ncrys:','	',n7,'	',n7*1./n
	  write(2,*)8,'	','n8=',n8,'	',n8*1./n
	  write(2,*)9,'	','n9=',n9,'	',n9*1./n
	  write(2,*)10,'	','n10=',n10,'	',n10*1./n
	  write(2,*)11,'	','n11=',n11,'	',n11*1./n
	  write(2,*)12,'	','n12=',n12,'	',n12*1./n
	  write(2,*)13,'	','n13=',n13,'	',n13*1./n
	  write(2,*)14,'	','n14=',n14,'	',n14*1./n
	  write(2,*)15,'	','n15=',n15,'	',n15*1./n
	  write(2,*)16,'	','n16=',n16,'	',n16*1./n
	  write(2,*)'nO-O=','	',n
	  write(2,*)'Zo-o=','	',n1*1./n
c	----cpO=45
	n7=0
	n8=0
	n9=0
	n10=0
	n11=0
	n12=0
	n13=0
	n14=0
	n15=0
	n16=0
	do i=nal+1,na
	  if(cpO(i).eq.45)then
	    if(cp1(i).ne.0)then
	      n7=n7+1
	    endif
	    select case(neib(i))
	    case(8)
	      n8=n8+1
	    case(9)
	      n9=n9+1
	    case(10)
	      n10=n10+1
	    case(11)
	      n11=n11+1
	    case(12)
	      n12=n12+1
	    case(13)
	      n13=n13+1
	    case(14)
	      n14=n14+1
	    case(15)
	      n15=n15+1
	    case(16)
	      n16=n16+1
	    end select
	  endif
	enddo
	n=n8+n9+n10+n11+n12+n13+n14+n15+n16
	n1=n8*8+n9*9+n10*10+n11*11+n12*12+n13*13+n14*14+n15*15+n16*16
c	  write to file:
	  write(2,*)'ZAl-O=45:'
	  write(2,*)'ncrys:','	',n7,'	',n7*1./n
	  write(2,*)8,'	','n8=',n8,'	',n8*1./n
	  write(2,*)9,'	','n9=',n9,'	',n9*1./n
	  write(2,*)10,'	','n10=',n10,'	',n10*1./n
	  write(2,*)11,'	','n11=',n11,'	',n11*1./n
	  write(2,*)12,'	','n12=',n12,'	',n12*1./n
	  write(2,*)13,'	','n13=',n13,'	',n13*1./n
	  write(2,*)14,'	','n14=',n14,'	',n14*1./n
	  write(2,*)15,'	','n15=',n15,'	',n15*1./n
	  write(2,*)16,'	','n16=',n16,'	',n16*1./n
	  write(2,*)'nO-O=','	',n
	  write(2,*)'Zo-o=','	',n1*1./n
c	----cpO=46
	n7=0
	n8=0
	n9=0
	n10=0
	n11=0
	n12=0
	n13=0
	n14=0
	n15=0
	n16=0
	do i=nal+1,na
	  if(cpO(i).eq.46)then
	    if(cp1(i).ne.0)then
	      n7=n7+1
	    endif
	    select case(neib(i))
	    case(8)
	      n8=n8+1
	    case(9)
	      n9=n9+1
	    case(10)
	      n10=n10+1
	    case(11)
	      n11=n11+1
	    case(12)
	      n12=n12+1
	    case(13)
	      n13=n13+1
	    case(14)
	      n14=n14+1
	    case(15)
	      n15=n15+1
	    case(16)
	      n16=n16+1
	    end select
	  endif
	enddo
	n=n8+n9+n10+n11+n12+n13+n14+n15+n16
	n1=n8*8+n9*9+n10*10+n11*11+n12*12+n13*13+n14*14+n15*15+n16*16
c	  write to file:
	  write(2,*)'ZAl-O=46:'
	  write(2,*)'ncrys:','	',n7,'	',n7*1./n
	  write(2,*)8,'	','n8=',n8,'	',n8*1./n
	  write(2,*)9,'	','n9=',n9,'	',n9*1./n
	  write(2,*)10,'	','n10=',n10,'	',n10*1./n
	  write(2,*)11,'	','n11=',n11,'	',n11*1./n
	  write(2,*)12,'	','n12=',n12,'	',n12*1./n
	  write(2,*)13,'	','n13=',n13,'	',n13*1./n
	  write(2,*)14,'	','n14=',n14,'	',n14*1./n
	  write(2,*)15,'	','n15=',n15,'	',n15*1./n
	  write(2,*)16,'	','n16=',n16,'	',n16*1./n
	  write(2,*)'nO-O=','	',n
	  write(2,*)'Zo-o=','	',n1*1./n
c	----cpO=56
	n7=0
	n8=0
	n9=0
	n10=0
	n11=0
	n12=0
	n13=0
	n14=0
	n15=0
	n16=0
	do i=nal+1,na
	  if(cpO(i).eq.56)then
	    if(cp1(i).ne.0)then
	      n7=n7+1
	    endif
	    select case(neib(i))
	    case(8)
	      n8=n8+1
	    case(9)
	      n9=n9+1
	    case(10)
	      n10=n10+1
	    case(11)
	      n11=n11+1
	    case(12)
	      n12=n12+1
	    case(13)
	      n13=n13+1
	    case(14)
	      n14=n14+1
	    case(15)
	      n15=n15+1
	    case(16)
	      n16=n16+1
	    end select
	  endif
	enddo
	n=n8+n9+n10+n11+n12+n13+n14+n15+n16
	n1=n8*8+n9*9+n10*10+n11*11+n12*12+n13*13+n14*14+n15*15+n16*16
c	  write to file:
	  write(2,*)'ZAl-O=56:'
	  write(2,*)'ncrys:','	',n7,'	',n7*1./n
	  write(2,*)8,'	','n8=',n8,'	',n8*1./n
	  write(2,*)9,'	','n9=',n9,'	',n9*1./n
	  write(2,*)10,'	','n10=',n10,'	',n10*1./n
	  write(2,*)11,'	','n11=',n11,'	',n11*1./n
	  write(2,*)12,'	','n12=',n12,'	',n12*1./n
	  write(2,*)13,'	','n13=',n13,'	',n13*1./n
	  write(2,*)14,'	','n14=',n14,'	',n14*1./n
	  write(2,*)15,'	','n15=',n15,'	',n15*1./n
	  write(2,*)16,'	','n16=',n16,'	',n16*1./n
	  write(2,*)'nO-O=','	',n
	  write(2,*)'Zo-o=','	',n1*1./n
c	----cpO=456
	n7=0
	n8=0
	n9=0
	n10=0
	n11=0
	n12=0
	n13=0
	n14=0
	n15=0
	n16=0
	do i=nal+1,na
	  if(cpO(i).eq.456)then
	    if(cp1(i).ne.0)then
	      n7=n7+1
	    endif
	    select case(neib(i))
	    case(8)
	      n8=n8+1
	    case(9)
	      n9=n9+1
	    case(10)
	      n10=n10+1
	    case(11)
	      n11=n11+1
	    case(12)
	      n12=n12+1
	    case(13)
	      n13=n13+1
	    case(14)
	      n14=n14+1
	    case(15)
	      n15=n15+1
	    case(16)
	      n16=n16+1
	    end select
	  endif
	enddo
	n=n8+n9+n10+n11+n12+n13+n14+n15+n16
	n1=n8*8+n9*9+n10*10+n11*11+n12*12+n13*13+n14*14+n15*15+n16*16
c	  write to file:
	  write(2,*)'ZAl-O=456:'
	  write(2,*)'ncrys:','	',n7,'	',n7*1./n
	  write(2,*)8,'	','n8=',n8,'	',n8*1./n
	  write(2,*)9,'	','n9=',n9,'	',n9*1./n
	  write(2,*)10,'	','n10=',n10,'	',n10*1./n
	  write(2,*)11,'	','n11=',n11,'	',n11*1./n
	  write(2,*)12,'	','n12=',n12,'	',n12*1./n
	  write(2,*)13,'	','n13=',n13,'	',n13*1./n
	  write(2,*)14,'	','n14=',n14,'	',n14*1./n
	  write(2,*)15,'	','n15=',n15,'	',n15*1./n
	  write(2,*)16,'	','n16=',n16,'	',n16*1./n
	  write(2,*)'nO-O=','	',n
	  write(2,*)'Zo-o=','	',n1*1./n
c	  -------------
	endfile(2)
c	+++++++++++++++++
	ncr=0
	non=0
	ntol=0
	do i=nal+1,na
	  do j=nal+1,na
		if(i.ne.j)then
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
	      if((rij.gt.3.50).and.(rij.lt.3.90))then
			nmi(i)=nmi(i)+1
	        ntol=ntol+1
	        if((cp1(i).ne.0).or.(cp1(j).ne.0))then
                ncr=ncr+1
	        endif
			if((cp1(i).eq.0).and.(cp1(j).eq.0))then
			  non=non+1
			endif
	      endif
	   endif
	  enddo
	enddo
	o4=0
	o5=0
	o6=0
	o45=0
	o46=0
	o56=0
	o456=0
	n4=0
	n5=0
	n6=0
	n7=0
	n8=0
	n9=0
	n10=0
	ic1=0
	ic2=0
	ic3=0
	ic4=0
	ic5=0
	ic6=0
	ic7=0
	do i=nal+1,na
	  selectcase(cpO(i))
	  case(4)
	    o4=o4+nmi(i)
	    n4=n4+1
	    if(cp1(i).ne.0)then
	      ic1=ic1+1
		endif
	  case(5)
	    o5=o5+nmi(i)
	    n5=n5+1
	    if(cp1(i).ne.0)then
	      ic2=ic2+1
		endif
	  case(6)
	    o6=o6+nmi(i)
	    n6=n6+1
	    if(cp1(i).ne.0)then
	      ic3=ic3+1
		endif
	  case(45)
	    o45=o45+nmi(i)
	    n7=n7+1
	    if(cp1(i).ne.0)then
	      ic4=ic4+1
		endif
	  case(46)
	    o46=o46+nmi(i)
	    n8=n8+1
	    if(cp1(i).ne.0)then
	      ic5=ic5+1
		endif
	  case(56)
	    o56=o56+nmi(i)
	    n9=n9+1
	    if(cp1(i).ne.0)then
	      ic6=ic6+1
		endif
	  case(456)
	    o456=o456+nmi(i)
	    n10=n10+1
	    if(cp1(i).ne.0)then
	      ic7=ic7+1
		endif
	  endselect
	enddo
	n=o4+o5+o6+o45+o46+o56+o456	
c	
	open(3,file='distfcc.txt')
	  write(3,*)'ncr=','	',ncr
	  write(3,*)'non=','	',non
	  write(3,*)'ntol=','	',ntol
	  write(3,*)'Distribution of intermediate O-O:'
	  write(3,*)'Type:	O-O pair:	number:		%:'
	  write(3,*)'O4:','	',o4,'	',n4,'	',o4*1./n4,'	',o4*1./n
	  write(3,*)'O5:','	',o5,'	',n5,'	',o5*1./n5,'	',o5*1./n
	  write(3,*)'O6:','	',o6,'	',n6,'	',o6*1./n6,'	',o6*1./n
	  write(3,*)'O45:','	',o45,'	',n7,'	',o45*1./n7,'	',o45*1./n
	  write(3,*)'O46:','	',o46,'	',n8,'	',o46*1./n8,'	',o46*1./n
	  write(3,*)'O56:','	',o56,'	',n9,'	',o56*1./n9,'	',o56*1./n
	  write(3,*)'O456:','	',o456,'	',n10,'	',o456*1./n10,'	',o456*1./n
	  write(3,*)'Ototal:','	',o4+o5+o6+o45+o46+o56+o456
	  write(3,*)'	'

	  write(3,*)'Number of crys in ZAL-O:'
	  write(3,*)'Type:	Ncrys:	NZAl-O:		%:'
	  write(3,*)'O4:','	',ic1,'	',n4,'	',ic1*1./n4
	  write(3,*)'O5:','	',ic2,'	',n5,'	',ic2*1./n5
	  write(3,*)'O6:','	',ic3,'	',n6,'	',ic3*1./n6
	  write(3,*)'O45:','	',ic4,'	',n7,'	',ic4*1./n7
	  write(3,*)'O46:','	',ic5,'	',n8,'	',ic5*1./n8
	  write(3,*)'O56:','	',ic6,'	',n9,'	',ic6*1./n8
	  write(3,*)'O456:','	',ic7,'	',n10,'	',ic7*1./n9
	  n1=ic1+ic2+ic3+ic4+ic5+ic6+ic7
	  n2=na-nal-n1
	  write(3,*)'crys:		non-crys:'
	  write(3,*)ncr*1./n1,'	',non*1./n2
	endfile(3)
c	+++++++++++++++++

c	========
	end
c	*********************************
