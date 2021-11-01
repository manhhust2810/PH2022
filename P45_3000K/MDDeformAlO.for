	program MD_DeformAl2O3
c	This MD program
c	The program was writen by:
c	                 Le Van Vinh, PhD
c	        Department of Computational Physics
c	      Hanoi University of Technology
	integer i,j,na,n,kk
	integer cp(9999)
	real x2(9999),y2(9999),z2(9999),vx(9999),vy(9999),vz(9999)
	real ep(9999),fx(9999),fy(9999),fz(9999),mass(9999)
	real fp(9999),fbx(9999),fby(9999),fbz(9999),vol(9999)
	real lx,ly,lz,h,h2,temp,tmpkin,pr,pr0,ekin,epot,etot
	real dx,dy,dz,delP,kb,fpp,ox,oy,oz,epsi,lz0,Vo

c	a0=4.38/0.529177
	kb=8.617343*1e-5

	h=1.0/0.024189
	h2=h*h

	temp=3000.0
	delP=0.0
	pr0=45.0


	call readfile(na,lx,ly,lz,x2,y2,z2,cp)
	write(*,*)'lx=',lx,'	ly=',ly,'	lz=',lz
	write(*,*)na,' ',x2(na),' ',y2(na),' ',z2(na),' ',cp(na)
	lz0=lz
	Vo=lz0*lz0*lz0/na
	write(*,*)'Volume atom=',Vo
c	---------------
c	initial velocity
	do i=1,na
	  if (cp(i).eq.1) then
	    mass(i)=26.9815*1836.15
	    vol(i) =Vo
	  elseif (cp(i).eq.2) then
	    mass(i)=15.999*1836.15
	    vol(i)=Vo
	  endif
	  vx(i)=0.0
	  vy(i)=0.0
	  vz(i)=0.0

	  fx(i)=0.0
	  fy(i)=0.0
	  fz(i)=0.0

	enddo
c	--------------
c	Make a file:
	open(6,file='energyAlO3000K.txt')
	write(6,*)'temp:	  etot:	  ekin:	  epot:'
c	------------------------------------
c	The loops
	kk=0
	do while (kk.lt.100000)
c	  -------------
c	  displacement:
c	  -------------
c	  temp=3500.0-kk*0.01
c	  -------------
	  ekin=0.0
	  do i=1,na
	    x2(i)=x2(i)+h*vx(i)+h2*fx(i)/2.0/mass(i)
	    y2(i)=y2(i)+h*vy(i)+h2*fy(i)/2.0/mass(i)
	    z2(i)=z2(i)+h*vz(i)+h2*fz(i)/2.0/mass(i)
	    
		if(x2(i).gt.lx)then
	       x2(i)=x2(i)-lx
	    endif
	    if(x2(i).lt.0)then
	      x2(i)=x2(i)+lx
	    endif
	    if(y2(i).gt.ly)then
	      y2(i)=y2(i)-ly
		endif
	    if(y2(i).lt.0)then
	      y2(i)=y2(i)+ly
		endif
	    if(z2(i).gt.lz)then
	      z2(i)=z2(i)-lz
		endif
	    if(z2(i).lt.0)then
	      z2(i)=z2(i)+lz
	    endif

	    fbx(i)=fx(i)
	    fby(i)=fy(i)
	    fbz(i)=fz(i)
	  enddo
c	  -----------------
c	  calculate the force acting on the atoms:	  
	  call calcforce(na,lx,ly,lz,x2,y2,z2,cp,fx,fy,fz,ep,fp,
     +	ox,oy,oz)
	  do i=1,na
	    vx(i)=vx(i)+h*(fx(i)+fbx(i))/2.0/mass(i)
	    vy(i)=vy(i)+h*(fy(i)+fby(i))/2.0/mass(i)
	    vz(i)=vz(i)+h*(fz(i)+fbz(i))/2.0/mass(i)
	    ekin=ekin+mass(i)*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
	  enddo
	
	  ekin=ekin/2.0

c	  rescale the velocity:	  
	  call temp_scal(na,temp,tmpkin,vx,vy,vz,ekin)
c
c	  if(mod(kk,5).eq.0)then
	    call pr_scal(na,lx,ly,lz,x2,y2,z2,ekin,pr,pr0,fp)
c	  endif
c	  ----------------
	  if (mod(kk,500).eq.0) then
c	    calculate the pressure:
c	    fpp=0.0
c	    do i=1,na
c	      fpp=fpp+fp(i)
c	    enddo
c	    fpp=2.0*fpp
c		pr=(2.0*ekin+fpp)/3.0
c	    pr=29392.92*pr/lx/ly/lz
c	    ------------
	    epot=0.0
	    do i=1,na
	      epot=epot+ep(i)
	    enddo
	    epot=2.*epot
	    etot=ekin+epot
c	  -------------
	    write(*,*)kk,' ep=',epot,'	lx=',lx,' ',tmpkin
	    write(6,*)tmpkin,'	',etot,'	',ekin,'	',epot,'	',lx
	  endif
c	  ------------

	  selectcase(kk)
	  case (20000)
    	    open(4,file='xyz1.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (40000)
    	    open(4,file='xyz2.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (60000)
    	    open(4,file='xyz3.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (80000)
    	    open(4,file='xyz4.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (100000)
    	    open(4,file='xyz5.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (120000)
    	    open(4,file='xyz6.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (140000)
    	    open(4,file='xyz7.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (160000)
    	    open(4,file='xyz8.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (180000)
    	    open(4,file='xyz9.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (200000)
    	    open(4,file='xyz10.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (220000)
    	    open(4,file='xyz11.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (240000)
    	    open(4,file='xyz12.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (260000)
    	    open(4,file='xyz13.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (280000)
    	    open(4,file='xyz14.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  case (300000)
    	    open(4,file='xyz15.txt')
	    write(4,*)lx,'	',ly,'	',lz
	    write(4,*)na
	    do i=1,na
	       write(4,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	    enddo
	    write(4,*)'temp=',temp
	    endfile(4)
	  endselect
	  kk=kk+1
c	---------
	enddo
c	close the file:
	endfile(6)	
c	------------------
	open(3,file='xyzAlO3000K.txt')
	  write(3,*)lx,'	',ly,'	',lz
	  write(3,*)na
	  do i=1,na
	     write(3,*)i,' ',x2(i),' ',y2(i),' ',z2(i),'	',cp(i)
	  enddo
	endfile(3)

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
c	***********************************
c	Calculated the force1=ZiZj/rij^2 between two atoms with a distance R
	Real function forc1(rij,cp1,cp2)
c	parameter(pi=3.14159265)
	integer cp1,cp2
	real zp(2)
	real rij,r,aph,ff
	  zp(1)=1.4175
	  zp(2)=-0.9450
	  aph=0.13
	  r=rij*aph
	  ff=zp(cp1)*zp(cp2)/rij/rij
	  forc1=ff*(erfc(r)+2.0*aph*rij*exp(-r*r)/1.772454)
      return
	End

c	***************************************
c	Calculated the F(r)=-6*Ci*Cj/rij^-7 +D*Exp([Ai+Aj-rij]/[Bi+Bj]) between two atoms with a distance R
	Real function forc2(rij,cp1,cp2)
	integer cp1,cp2
	real A(2),B(2),C(2)
	real rij,r,D,ff,Aij,Bij
c	  parameters of Buckingham potential:1-Al, 2-Si, 3-O;
	  A(1)=1.483814
	  A(2)=3.442138
	  B(1)=0.064251
	  B(2)=0.260782
	  C(1)=4.849199
	  C(2)=11.933349
	  D=0.0008433
c	  options: k=1-AlAl, k=2-AlN, k=3-NN
	  r=rij*rij
	  r=r*r*r
	  r=r*rij
	  Bij=B(cp1)+B(cp2)
	  Aij=A(cp1)+A(cp2)-rij
	  ff=-6.0*C(cp1)*C(cp2)/r
	  ff=ff+D*exp(Aij/Bij)
	  forc2=ff
	return
	End
c	***************************************
c	calculate Phi(rij)
	real function phi1(rij,cp1,cp2)
c	parameter(pi=3.14159265)
	integer cp1,cp2
	real zp(3)
	real rij,r,aph,rc
	  zp(1)=1.4175
	  zp(2)=-0.9450
	  aph=0.13
	  r=rij*aph
	  phi1=0.5*zp(cp1)*zp(cp2)*erfc(r)/rij
	return
	end
c	*************************************
	real function phi2(rij,cp1,cp2)
	integer cp1,cp2
	real A(2),B(2),C(2)
	real rij,r,D,ff,Bij,Aij
c	  parameters of Buckingham potential:1-Al, 2-Si, 3-O;
	  A(1)=1.483814
	  A(2)=3.442138
	  B(1)=0.064251
	  B(2)=0.260782
	  C(1)=4.849199
	  C(2)=11.933349
	  D=0.0008433
c	  options: k=1-AlAl, k=2-AlN, k=3-NN
	  r=rij*rij
	  r=r*r*r
	  Bij=B(cp1)+B(cp2)
	  Aij=A(cp1)+A(cp2)-rij
	  ff=-C(cp1)*C(cp2)/r
	  ff=ff+D*Bij*exp(Aij/Bij)
	  phi2=ff
	return
	end
c	*************************************
c	calculate the force acting on the atoms	
	subroutine calcforce(na,lx,ly,lz,x2,y2,z2,cp,fx,fy,fz,ep,fp,
     +	ox,oy,oz)
	parameter (rc=12.2)
	integer i,j,na
	integer cp(9999)
	real x2(9999),y2(9999),z2(9999),fx(9999),fy(9999),fz(9999)
	real ep(9999),fp(9999)
	real rij,rx,ry,rz,lx,ly,lz,f,fc,uc,oz,oy,ox
	ox=0.0
	oy=0.0
	oz=0.0
	do i=1,na
	  fx(i)=0.0
	  fy(i)=0.0
	  fz(i)=0.0
	  ep(i)=0.0
	  fp(i)=0.0
	enddo
	fc=0.0
	uc=0.0
	do i=1,na-1
	  do j=i+1,na
	    rx=pbc(lx,x2(i),x2(j))
	    ry=pbc(ly,y2(i),y2(j))
	    rz=pbc(lz,z2(i),z2(j))
	    rij=sqrt(rx*rx+ry*ry+rz*rz)
	    if (rij.lt.rc) then
	      f=forc1(rij,cp(i),cp(j))/rij
	      f=f+forc2(rij,cp(i),cp(j))/rij
	      fc=forc1(rc,cp(i),cp(j))
	      fc=fc+forc2(rc,cp(i),cp(j))
	      f=f-fc/rc
	      fx(i)=fx(i)-f*rx
	      fx(j)=fx(j)+f*rx
	      fy(i)=fy(i)-f*ry
	      fy(j)=fy(j)+f*ry
	      fz(i)=fz(i)-f*rz
	      fz(j)=fz(j)+f*rz
c	      ----------------
c	      calculate the potential:
	      ep(i)=ep(i)+phi1(rij,cp(i),cp(j))
	      ep(i)=ep(i)+phi2(rij,cp(i),cp(j))
	      uc=phi1(rc,cp(i),cp(j))
	      uc=uc+phi2(rc,cp(i),cp(j))
	      ep(i)=ep(i)-uc+fc*(rij-rc)
c	      -----------------
c	      calculate the part of pressure:
	      fp(i)=fp(i)+f*rij*rij
c	      fp(i)=fp(i)+phi1(rij,cp(i),cp(j))
c	      fp(i)=fp(i)-phi1(rc,cp(i),cp(j))+fc*(rij-rc)
c	      calculate the stress:
c	      oz=oz+f*rz*rz
c	      ox=ox+f*rx*rx
c	      oy=oy+f*ry*ry
	    endif
	  enddo
	enddo
	end
c	***********************************
	real function erfc(x)
	parameter(p=0.3275911,a1=0.254829592,a2=-0.284496736,
     +          a3=1.421413741,a4=-1.453152027,a5=1.061405429)
	real T,x,h,T2
	
	T=1./(1.+p*x)
	T2=T*T
	h=T*(a1+a3*T2+a5*T2*T2)
	h=h+a2*T2+a4*T2*T2
	erfc=h*exp(-x*x)
	return
	end
c	***********************************
c	Adjust to Temperature
	subroutine temp_scal(na,temp,tmpkin,vx,vy,vz,ekin)
	integer i,na
	real temp,tmpkin,ekin,kt,fact
	real vx(9999),vy(9999),vz(9999)
	kt=8.617343
	tmpkin=2.0*ekin*27.212
	tmpkin=tmpkin*100000.0/3.0/kt/na
	if(tmpkin .gt. temp) then
	   tmpkin=tmpkin-0.01
	else
	   tmpkin=tmpkin+0.01
	endif
	fact=sqrt(temp/tmpkin)
	do i=1,na
	  vx(i)=fact*vx(i)
	  vy(i)=fact*vy(i)
	  vz(i)=fact*vz(i)
	enddo
	end
c	************************************
c	Adjust to Pressure!
	subroutine pr_scal(na,lx,ly,lz,x2,y2,z2,ekin,pr,pr0,fp)
	integer i,na
	real x2(9999),y2(9999),z2(9999)
	real fp(9999)
	real lx,ly,lz,ekin,pr,pr0,dp,fpp,fact
	fpp=0.0
	dp=0.0001
	do i=1,na
	  fpp=fpp+fp(i)
	enddo
	fpp=2.0*fpp
	pr=2.0*ekin/3.0
	pr=pr+fpp/3.0
	pr=pr*29392.92
	pr=pr/lx/ly/lz

	if (pr.le.pr0) then
	  fact=1.0-dp
	else
	  fact=1+dp
	endif
	do i=1,na
	  x2(i)=x2(i)*fact
	  y2(i)=y2(i)*fact
	  z2(i)=z2(i)*fact
	enddo
	lx=lx*fact
	ly=ly*fact
	lz=lz*fact
	end
