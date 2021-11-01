	program PRDF
c	This MD program uses the Sutton-Chen embedded-atom potentials
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
	call neighbor(na,lx,ly,lz,x2,y2,z2,cp)
c	call PRDFs(na,lx,ly,lz,x2,y2,z2,cp)
	call distangle(na,lx,ly,lz,x2,y2,z2,cp)
	call angleMOM(na,lx,ly,lz,x2,y2,z2,cp)
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
	if (r.gt.l2) then
	  r=-l+r
	elseif (r.lt.-l2) then
	  r=l+r
	endif
	pbc=r
	return
	end

c	***********************************
	subroutine neighbor(na,lx,ly,lz,x2,y2,z2,cp)
	integer i,j,na,nal,n2,n3,n4,n5,n6,n7,n8,n
	integer no2,no3,no4,no5,no6,no7,no
	integer g4(100),g5(100),g6(100)
	integer cp(9999),neib(9999),neia(9999),fnei(9999,20)
	real lx,ly,lz,rc,rx,ry,rz,rij,r1,rb
	real x2(9999),y2(9999),z2(9999),r4(9999,4),r5(9999,5),r6(9999,7)

	write(*,*)cp(na),'	',lx,'	',x2(na)
	
	nal=0
	rc=2.61
	do i=1,na
	  if(cp(i).eq.1)then
	    nal=nal+1
	  endif
	enddo
	write(*,*)'nAl=',nal

	do i=1,na
	  neib(i)=0
	  neia(i)=0
	enddo
	do i=1,nal
	  n=0
	  do j=i+1,na
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
	      if (rij.lt.rc) then
	        if (cp(j).eq.2) then
	          neib(i)=neib(i)+1
	          n=n+1
	          fnei(i,n)=j
			endif
		  endif
	  enddo
	enddo

	do i=nal+1,na
	  n=0
	  do j=1,nal
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
	      if (rij.lt.rc) then
	          neia(i)=neia(i)+1
c	          n=n+1
c	          fnei(i,n)=j
		  endif
	  enddo
	enddo
c	----------------------
c	calculate distribution	
	no2=0
	no3=0
	no4=0
	no5=0
	no6=0
	no7=0
	do i=nal+1,na
	  select case(neia(i))
	  case(2)
	    no2=no2+1
	  case(3)
	    no3=no3+1
	  case(4)
	    no4=no4+1
	  case(5)
	    no5=no5+1
	  case(6)
	    no6=no6+1
	  case(7)
	    no7=no7+1
	  end select
	enddo
	no=no2+no3+no4+no5+no6+no7
c	----------------------
c	calculate distribution	
	n3=0
	n4=0
	n5=0
	n6=0
	n7=0
	n8=0
	do i=1,nal
	  select case(neib(i))
	  case(3)
	    n3=n3+1
	  case(4)
	    n4=n4+1
	  case(5)
	    n5=n5+1
	  case(6)
	    n6=n6+1
	  case(7)
	    n7=n7+1
	  case(8)
	    n8=n8+1
	  end select
	enddo
	n=n3+n4+n5+n6+n7+n8
c	-------------------
	open(2,file='neigbor.txt')
	  write(2,*)'For Al:'
	  write(2,*)'n3=',n3,'	',n3*1./n
	  write(2,*)'n4=',n4,'	',n4*1./n
	  write(2,*)'n5=',n5,'	',n5*1./n
	  write(2,*)'n6=',n6,'	',n6*1./n
	  write(2,*)'n7=',n7,'	',n7*1./n
	  write(2,*)'n8=',n8,'	',n8*1./n
	  write(2,*)'nAl=',n
	  write(2,*)'For O:'
	  write(2,*)'no2=',no2,'	',no2*1./no
	  write(2,*)'no3=',no3,'	',no3*1./no
	  write(2,*)'no4=',no4,'	',no4*1./no
	  write(2,*)'no5=',no5,'	',no5*1./no
	  write(2,*)'no6=',no6,'	',no6*1./no
	  write(2,*)'no7=',no7,'	',no7*1./no
	  write(2,*)'nO=',no
	  do i=1,na
	    write(2,*)i,'	',neib(i),'	',neia(i)
	  enddo

	endfile(2)
c	+++++++++++++++++
c	calculate the distribution of the bond-length of units:
	n4=0
	n5=0
	n6=0
	do i=1,nal
	  select case(neib(i))
	  case(4)
	  do k=1,neib(i)
	    j=fnei(i,k)
	    rx=pbc(lx,x2(i),x2(j))
	    ry=pbc(ly,y2(i),y2(j))
	    rz=pbc(lz,z2(i),z2(j))
	    rij=sqrt(rx*rx+ry*ry+rz*rz)
	    r4(i,k)=rij
	    n4=n4+1
	  enddo
	  case(5)
	  do k=1,neib(i)
	    j=fnei(i,k)
	    rx=pbc(lx,x2(i),x2(j))
	    ry=pbc(ly,y2(i),y2(j))
	    rz=pbc(lz,z2(i),z2(j))
	    rij=sqrt(rx*rx+ry*ry+rz*rz)
	    r5(i,k)=rij
	    n5=n5+1
	  enddo
	  case(6)
	  do k=1,neib(i)
	    j=fnei(i,k)
	    rx=pbc(lx,x2(i),x2(j))
	    ry=pbc(ly,y2(i),y2(j))
	    rz=pbc(lz,z2(i),z2(j))
	    rij=sqrt(rx*rx+ry*ry+rz*rz)
	    r6(i,k)=rij
	    n6=n6+1
	  enddo
	  case(7)
	  do k=1,neib(i)
	    j=fnei(i,k)
	    rx=pbc(lx,x2(i),x2(j))
	    ry=pbc(ly,y2(i),y2(j))
	    rz=pbc(lz,z2(i),z2(j))
	    rij=sqrt(rx*rx+ry*ry+rz*rz)
	    r6(i,k)=rij
	    n6=n6+1
	  enddo
	  end select
	enddo
	do i=1,100
	  g4(i)=0
	  g5(i)=0
	  g6(i)=0
	enddo

	n4=0
	n5=0
	n6=0	
	do i=1,60
	  r1=1.5+0.02*(i-1)
	  rb=r1+0.02
	  do j=1,nal
	    select case(neib(j))
	    case(4)
	      do k=1,neib(j)
	        if((r4(j,k).gt.r1).and.(r4(j,k).le.rb))then
	          g4(i)=g4(i)+1
	          n4=n4+1
			endif
		  enddo
	    case(5)
	      do k=1,neib(j)
	        if((r5(j,k).gt.r1).and.(r5(j,k).le.rb))then
	          g5(i)=g5(i)+1
			  n5=n5+1
			endif
		  enddo
	    case(6)
	      do k=1,neib(j)
	        if((r6(j,k).gt.r1).and.(r6(j,k).le.rb))then
	          g6(i)=g6(i)+1
			  n6=n6+1
			endif
		  enddo
	    case(7)
	      do k=1,neib(j)
	        if((r6(j,k).gt.r1).and.(r6(j,k).le.rb))then
	          g6(i)=g6(i)+1
			  n6=n6+1
			endif
		  enddo
	    end select
	  enddo
	enddo
	open(3,file='bldistriAl_O.txt')
	do i=1,60
	  r1=1.5+0.02*(i-1)+0.01
	  write(3,*)r1,'	',g4(i)*1./n4,'	',g5(i)*1./n5,'	',g6(i)*1./n6
	enddo
	endfile(3)
c	+++++++++++++++++
c	-------------------
	end
c	***********************************
	subroutine distangle(na,lx,ly,lz,x2,y2,z2,cp)
	parameter(pi=3.14159265)
	integer i,j,na,nal,no,n,k,l,m,n4,n5,n6,nt,pb
	integer cp(9999),neib(9999),pba4(120),pba5(120),pba6(120)
	integer fnei(9999,10)
	real lx,ly,lz,rc,rx,ry,rz,rij,rx1,ry1,rz1,rim,dag,r1,r2,x
	real x2(9999),y2(9999),z2(9999)
	real angle(9999,16)

	rc=2.23

	nal=0
	do i=1,na
	  if(cp(i).eq.1)then
	    nal=nal+1
	  endif
	enddo
	write(*,*)'nAl=',nal	

	do i=1,nal
	  neib(i)=0
	  do j=1,16
	    angle(i,j)=0.0
	  enddo
	enddo
	do i=1,nal
	  n=0
	  do j=i+1,na
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
	      if (rij.lt.rc) then
	        if (cp(j).eq.2) then
	          neib(i)=neib(i)+1
	          n=n+1
	          fnei(i,n)=j
			endif
		  endif
	  enddo
	enddo
c	calculate angle(i,no)
	do i=1,nal
	  select case(neib(i))
	  case(4)
	    no=0
		do k=1,3
	      j=fnei(i,k)
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
		  do l=k+1,4
	        m=fnei(i,l)
	        rx1=pbc(lx,x2(i),x2(m))
	        ry1=pbc(ly,y2(i),y2(m))
	        rz1=pbc(lz,z2(i),z2(m))
	        rim=sqrt(rx1*rx1+ry1*ry1+rz1*rz1)
	        no=no+1
	        x=(rx*rx1+ry*ry1+rz*rz1)/rij/rim
			angle(i,no)=Acos(x)*180.0/pi
		  enddo
	    enddo
	  case(5)
	    no=0
		do k=1,4
	      j=fnei(i,k)
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
		  do l=k+1,5
	        m=fnei(i,l)
	        rx1=pbc(lx,x2(i),x2(m))
	        ry1=pbc(ly,y2(i),y2(m))
	        rz1=pbc(lz,z2(i),z2(m))
	        rim=sqrt(rx1*rx1+ry1*ry1+rz1*rz1)
	        no=no+1
	        x=(rx*rx1+ry*ry1+rz*rz1)/rij/rim
			angle(i,no)=Acos(x)*180.0/pi
		  enddo
	    enddo
	  case(6)
	    no=0
		do k=1,5
	      j=fnei(i,k)
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
		  do l=k+1,6
	        m=fnei(i,l)
	        rx1=pbc(lx,x2(i),x2(m))
	        ry1=pbc(ly,y2(i),y2(m))
	        rz1=pbc(lz,z2(i),z2(m))
	        rim=sqrt(rx1*rx1+ry1*ry1+rz1*rz1)
	        no=no+1
	        x=(rx*rx1+ry*ry1+rz*rz1)/rij/rim
			angle(i,no)=Acos(x)*180.0/pi
		  enddo
	    enddo

	  end select
	enddo

c	------------
c	calculate the angle distribution
	dag=3.0
	n4=0
	n5=0
	n6=0
	do i=1,120
	  pba4(i)=0
	  pba5(i)=0
	  pba6(i)=0
	enddo
	do i=1,nal
	  select case(neib(i))
	  case(4)
	    do j=1,120
	      r1=(j-1)*dag
	      r2=r1+dag
	      do k=1,16
	        if ((angle(i,k).gt.r1).and.(angle(i,k).le.r2)) then
	          pba4(j)=pba4(j)+1
	          n4=n4+1
	        endif
	      enddo
		enddo
	  case(5)
	    do j=1,120
	      r1=(j-1)*dag
	      r2=r1+dag
	      do k=1,16
	        if ((angle(i,k).gt.r1).and.(angle(i,k).le.r2)) then
	          pba5(j)=pba5(j)+1
	          n5=n5+1
	        endif
	      enddo
		enddo
	  case(6)
	    do j=1,120
	      r1=(j-1)*dag
	      r2=r1+dag
	      do k=1,16
	        if ((angle(i,k).gt.r1).and.(angle(i,k).le.r2)) then
	          pba6(j)=pba6(j)+1
	          n6=n6+1
	        endif
	      enddo
		enddo
	  end select
	enddo
	nt=n4+n5+n6
c	-----------
	open(6,file='distangle0.txt')
c
	  write(6,*)'angle:	pba:'
	  do j=1,120
	    pb=pba4(j)+pba5(j)+pba6(j)
		r1=(j-1)*dag+1.5
		write(6,*)r1,'	',pb,'	',pb*1.0/nt
	  enddo
c
	  write(6,*)'angle:	pba4:'
	  do j=1,120
		r1=(j-1)*dag+1.5
	    write(6,*)r1,'	',pba4(j),'	',1.*pba4(j)/n4
	  enddo
c
	  write(6,*)'angle:	pba5:'
	  do j=1,120
		r1=(j-1)*dag+1.5
	    write(6,*)r1,'	',pba5(j),'	',1.*pba5(j)/n5
	  enddo
c
	  write(6,*)'angle:	pba6:'
	  do j=1,120
		r1=(j-1)*dag+1.5
	    write(6,*)r1,'	',pba6(j),'	',1.*pba6(j)/n6
	  enddo
	endfile(6)
	
c	finished	
	end
c	***********************************
c	***********************************
	subroutine PRDFs(na,lx,ly,lz,x2,y2,z2,cp)
	parameter (pi=3.14159265)
	integer i,j,k,na,kk,n,naa,nbb,ncc,kkk,fac
	integer cp(9999),nn(800)
	integer n1(800),n2(800),n3(800)
	
	real x2(9999),y2(9999),z2(9999)
	real g1(800),g2(800),g3(800),g4(800),g5(800),g6(800)
	real gg(800),S(800),gN(800),gX(800),Sn(800),Sx(800)
	real c1,c2,b1,b2,f1,f2,sumN,sumX
	real r0,r1,r2,rij,dr,rx,ry,rz,lx,ly,lz,vol,r,sum,h,Dn,Df
	
	n=300
	r0=1.0
	dr=0.02

	vol=lx*ly*lz
	naa=0
	nbb=0
	ncc=0
	write(*,*)'dr=',dr

	do i=1,na
	  if (cp(i).eq.1)then
	    naa=naa+1
	  elseif(cp(i).eq.2)then
	    nbb=nbb+1
	  endif
	enddo
	write(*,*)'na=',naa,'nb=',nbb

	c1=1.0*naa/na
	c2=1.0*nbb/na
	write(*,*)'c1=',c1,' c2=',c2

	b1=0.00003449
	b2=0.00005805

	f1=8.144
	f2=14.41

	Dn=c1*b1+c2*b2
	Dn=Dn*Dn

	Df=c1*f1+c2*f2
	Df=Df*Df

	do i=1,n
	  n1(i)=0
	  n2(i)=0
	  n3(i)=0
	  nn(i)=0
	enddo
	do i=1,na-1
	  do j=1+i,na
	    rx=pbc(lx,x2(i),x2(j))
	    ry=pbc(ly,y2(i),y2(j))
	    rz=pbc(lz,z2(i),z2(j))
	    rij=sqrt(rx*rx+ry*ry+rz*rz)
c	    calcualte for G(r):
          do k=1,n
	      r1=r0+(k-1)*dr
	      r2=r1+dr
		  if ((rij .gt. r1).and.(rij.le.r2)) then
		    nn(k)=nn(k)+1
		  endif 
	    enddo
c	  options: kkk=1-AlO, kkk=2-AlAl, kkk=3-OO
	    kkk=0
	    if((cp(i).eq.1).and.(cp(j).eq.2)) then
	      kkk=1
	    elseif((cp(i).eq.2).and.(cp(j).eq.1)) then
	      kkk=1
	    elseif((cp(i).eq.1).and.(cp(j).eq.1)) then
	      kkk=2
	    elseif((cp(i).eq.2).and.(cp(j).eq.2)) then
	      kkk=3
	    endif

	    if (kkk.eq.1) then
	        do k=1,n
	          r1=r0+(k-1)*dr
	          r2=r1+dr
			  if ((rij .gt. r1).and.(rij.le.r2)) then
			     n1(k)=n1(k)+1
			  endif 
			enddo
	    elseif(kkk.eq.2)then
	        do k=1,n
	          r1=r0+(k-1)*dr
	          r2=r1+dr
			  if ((rij .gt. r1).and.(rij.le.r2)) then
			     n2(k)=n2(k)+1
			  endif 
			enddo
	    elseif(kkk.eq.3)then
	        do k=1,n
	          r1=r0+(k-1)*dr
	          r2=r1+dr
			  if ((rij .gt. r1).and.(rij.le.r2)) then
			     n3(k)=n3(k)+1
			  endif 
			enddo
		endif
	  enddo
	enddo

	open(7,file='PRDFs0.txt')
	do i=1,n
	  r=r0+(i-1)*dr+0.01
	  g1(i)=vol*n1(i)/4.0/pi/dr/naa/nbb/(r**2)
	  g2(i)=vol*n2(i)/2.0/pi/dr/naa/naa/(r**2)
	  g3(i)=vol*n3(i)/2.0/pi/dr/nbb/nbb/(r**2)
	  write(7,*)r,'	',g1(i),'	',g2(i),'	',g3(i)
	enddo
	endfile(7)

c	---------------------
c	calculate the Gn,x(r):
	do i=1,n
	  gN(i)=c1*c1*b1*b1*g2(i)+c2*c2*b2*b2*g3(i)+2*c1*c2*b1*b2*g1(i)
	  gN(i)=gN(i)/Dn

	  gX(i)=c1*c1*f1*f1*g2(i)+c2*c2*f2*f2*g3(i)+2*c1*c2*f1*f2*g1(i)
	  gX(i)=gX(i)/Df
	enddo

	open(9,file='Grtotal.txt')
	do i=1,n
	  r=r0+(i-1)*dr+0.01
	  gg(i)=vol*nn(i)/2.0/pi/dr/na/na/(r**2)
	  write(9,*)r,'	',gg(i),'	',gN(i),'	',gX(i)
	enddo
	endfile(9)
c	---------------------------
c	calculate structure factor:
	write(*,*)'calculate the structure factor...'
	do k=10,200
	  h=k*0.1
	  sum=(gg(1)-1)*1.01*1.01*sin(h*1.01)/(h*1.01)
	  sumN=(gN(1)-1)*1.01*1.01*sin(h*1.01)/(h*1.01)
	  sumX=(gX(1)-1)*1.01*1.01*sin(h*1.01)/(h*1.01)

	  fac=2
	  do i=2,n-1
	    if (fac.eq.2) then
	      fac=4
	    else
	      fac=2
		endif
	    r=r0+(i-1)*dr+0.01
		sum=sum+fac*(gg(i)-1)*r*r*sin(h*r)/(h*r)	  
		sumN=sumN+fac*(gN(i)-1)*r*r*sin(h*r)/(h*r)	  
		sumX=sumX+fac*(gX(i)-1)*r*r*sin(h*r)/(h*r)	  
	  enddo
	  r=r0+(n-1)*dr+0.01
	  sum=sum+(gg(n)-1)*r*r*sin(h*r)/(h*r)
	  sum=dr*sum/3.0
	  S(k)=1+4*pi*3000*sum/vol

	  sumN=sumN+(gN(n)-1)*r*r*sin(h*r)/(h*r)
	  sumN=dr*sumN/3.0
	  Sn(k)=1+4*pi*3000*sumN/vol

	  sumX=sumX+(gX(n)-1)*r*r*sin(h*r)/(h*r)
	  sumX=dr*sumX/3.0
	  Sx(k)=1+4*pi*3000*sumX/vol
	enddo
	open(9,file='Sfactor.txt')
	do k=10,200
	  write(9,*)k*0.1,'	',S(k),'	',Sn(k),'	',Sx(k)
	enddo
	endfile(9)
c	--------------------------
	end
c	***********************************
	subroutine angleMOM(na,lx,ly,lz,x2,y2,z2,cp)
	parameter(pi=3.14159265)
	integer i,j,na,nal,no,n,k,l,m,n2,n3,n4,n5,n6,nt,pb
	integer cp(9999),neib(9999),pba2(120),pba3(120),pba4(120)
	integer fnei(9999,20),pba5(120),pba6(120)
	real lx,ly,lz,rc,rx,ry,rz,rij,rx1,ry1,rz1,rim,dag,r1,r2,x
	real x2(9999),y2(9999),z2(9999)
	real angle(9999,20)

	rc=2.39

	nal=0
	do i=1,na
	  if(cp(i).eq.1)then
	    nal=nal+1
	  endif
	enddo
	write(*,*)'nAl=',nal	

	do i=1,nal
	  neib(i)=0
	  do j=1,20
	    angle(i,j)=0.0
	  enddo
	enddo
		

	do i=nal+1,na
	  n=0
	  do j=1,nal
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
	      if (rij.lt.rc) then
	          neib(i)=neib(i)+1
	          n=n+1
	          fnei(i,n)=j
		  endif
	  enddo
	enddo
c	calculate angle(i,no)
	do i=nal+1,na
	  select case(neib(i))
	  case(2)
	    no=0
		do k=1,1
	      j=fnei(i,k)
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
		  do l=k+1,2
	        m=fnei(i,l)
	        rx1=pbc(lx,x2(i),x2(m))
	        ry1=pbc(ly,y2(i),y2(m))
	        rz1=pbc(lz,z2(i),z2(m))
	        rim=sqrt(rx1*rx1+ry1*ry1+rz1*rz1)
	        no=no+1
	        x=(rx*rx1+ry*ry1+rz*rz1)/rij/rim
			angle(i,no)=Acos(x)*180.0/pi
		  enddo
	    enddo
	  case(3)
	    no=0
		do k=1,2
	      j=fnei(i,k)
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
		  do l=k+1,3
	        m=fnei(i,l)
	        rx1=pbc(lx,x2(i),x2(m))
	        ry1=pbc(ly,y2(i),y2(m))
	        rz1=pbc(lz,z2(i),z2(m))
	        rim=sqrt(rx1*rx1+ry1*ry1+rz1*rz1)
	        no=no+1
	        x=(rx*rx1+ry*ry1+rz*rz1)/rij/rim
			angle(i,no)=Acos(x)*180.0/pi
		  enddo
	    enddo
	  case(4)
	    no=0
		do k=1,3
	      j=fnei(i,k)
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
		  do l=k+1,4
	        m=fnei(i,l)
	        rx1=pbc(lx,x2(i),x2(m))
	        ry1=pbc(ly,y2(i),y2(m))
	        rz1=pbc(lz,z2(i),z2(m))
	        rim=sqrt(rx1*rx1+ry1*ry1+rz1*rz1)
	        no=no+1
	        x=(rx*rx1+ry*ry1+rz*rz1)/rij/rim
			angle(i,no)=Acos(x)*180.0/pi
		  enddo
	    enddo
	  case(5)
	    no=0
		do k=1,4
	      j=fnei(i,k)
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
		  do l=k+1,5
	        m=fnei(i,l)
	        rx1=pbc(lx,x2(i),x2(m))
	        ry1=pbc(ly,y2(i),y2(m))
	        rz1=pbc(lz,z2(i),z2(m))
	        rim=sqrt(rx1*rx1+ry1*ry1+rz1*rz1)
	        no=no+1
	        x=(rx*rx1+ry*ry1+rz*rz1)/rij/rim
			angle(i,no)=Acos(x)*180.0/pi
		  enddo
	    enddo
	  case(6)
	    no=0
		do k=1,5
	      j=fnei(i,k)
	      rx=pbc(lx,x2(i),x2(j))
	      ry=pbc(ly,y2(i),y2(j))
	      rz=pbc(lz,z2(i),z2(j))
	      rij=sqrt(rx*rx+ry*ry+rz*rz)
		  do l=k+1,6
	        m=fnei(i,l)
	        rx1=pbc(lx,x2(i),x2(m))
	        ry1=pbc(ly,y2(i),y2(m))
	        rz1=pbc(lz,z2(i),z2(m))
	        rim=sqrt(rx1*rx1+ry1*ry1+rz1*rz1)
	        no=no+1
	        x=(rx*rx1+ry*ry1+rz*rz1)/rij/rim
			angle(i,no)=Acos(x)*180.0/pi
		  enddo
	    enddo
	  end select
	enddo

c	------------
c	calculate the angle distribution
	dag=3.0
	n2=0
	n3=0
	n4=0
	n5=0
	n6=0
	do i=1,120
	  pba2(i)=0
	  pba3(i)=0
	  pba4(i)=0
	  pba5(i)=0
	  pba6(i)=0
	enddo
c	note: i=nal+1!!!
	do i=nal+1,na
	  select case(neib(i))
	  case(2)
	    do j=1,120
	      r1=(j-1)*dag
	      r2=r1+dag
	      do k=1,16
	        if ((angle(i,k).gt.r1).and.(angle(i,k).le.r2)) then
	          pba2(j)=pba2(j)+1
	          n2=n2+1
	        endif
	      enddo
		enddo
	  case(3)
	    do j=1,120
	      r1=(j-1)*dag
	      r2=r1+dag
	      do k=1,16
	        if ((angle(i,k).gt.r1).and.(angle(i,k).le.r2)) then
	          pba3(j)=pba3(j)+1
	          n3=n3+1
	        endif
	      enddo
		enddo
	  case(4)
	    do j=1,120
	      r1=(j-1)*dag
	      r2=r1+dag
	      do k=1,16
	        if ((angle(i,k).gt.r1).and.(angle(i,k).le.r2)) then
	          pba4(j)=pba4(j)+1
	          n4=n4+1
	        endif
	      enddo
		enddo
	  case(5)
	    do j=1,120
	      r1=(j-1)*dag
	      r2=r1+dag
	      do k=1,20
	        if ((angle(i,k).gt.r1).and.(angle(i,k).le.r2)) then
	          pba5(j)=pba5(j)+1
	          n5=n5+1
	        endif
	      enddo
		enddo
	  case(6)
	    do j=1,120
	      r1=(j-1)*dag
	      r2=r1+dag
	      do k=1,20
	        if ((angle(i,k).gt.r1).and.(angle(i,k).le.r2)) then
	          pba6(j)=pba6(j)+1
	          n6=n6+1
	        endif
	      enddo
		enddo
	  end select
	enddo
	nt=n2+n3+n4+n5+n6
c	-----------
	open(6,file='AngleMOM.txt')
c
	  write(6,*)'angle:	pba:'
	  do j=1,120
	    pb=pba2(j)+pba3(j)+pba4(j)+pba5(j)+pba6(j)
		r1=(j-1)*dag+1.5
		write(6,*)r1,'	',pb,'	',pb*1.0/nt
	  enddo
c
	  write(6,*)'angle:	pba2:'
	  do j=1,120
		r1=(j-1)*dag+1.5
	    write(6,*)r1,'	',pba2(j),'	',1.*pba2(j)/n2
	  enddo
c
	  write(6,*)'angle:	pba3:'
	  do j=1,120
		r1=(j-1)*dag+1.5
	    write(6,*)r1,'	',pba3(j),'	',1.*pba3(j)/n3
	  enddo
c
	  write(6,*)'angle:	pba4:'
	  do j=1,120
		r1=(j-1)*dag+1.5
	    write(6,*)r1,'	',pba4(j),'	',1.*pba4(j)/n4
	  enddo
c
	  write(6,*)'angle:	pba5:'
	  do j=1,120
		r1=(j-1)*dag+1.5
	    write(6,*)r1,'	',pba5(j),'	',1.*pba5(j)/n5
	  enddo
c
	  write(6,*)'angle:	pba6:'
	  do j=1,120
		r1=(j-1)*dag+1.5
	    write(6,*)r1,'	',pba6(j),'	',1.*pba6(j)/n6
	  enddo
	endfile(6)

c	finished	
	end
c	***********************************
