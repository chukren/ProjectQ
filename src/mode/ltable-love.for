C  Output from Public domain Ratfor, version 1.0
c  open(90,file='layer.love') for wkern.for, read-engi.for disperion curves ...   
c  recorded the layers of eigenfunctions and the phase velocity... 
c  => w  vphase gcom number-layer stratl endl           
cc
cc to calc low oder l<10
cc
c  => 10.0d  change to 2.0d
c  
c  => wdiff > 0.001
c
      subroutine ltable(iout,ioeig,mtype)
      parameter(nlmax=700)
      implicit real*8(a-h,o-z)
      integer*4 ls0,m1,kind,null
      real*4 tlen,dcda,dcdb,dcdrh,y1,dy1dr,dum,vphase,gcom,qmod,fl4,sqll
     *,ekin, wdiff,w4,zero,y2,dy2dr,rn4,vp4,vs4,rh4,q4,df4,cmax4,dcdxi
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl, fl1,fl
     *2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      common/eifx/a(5,nlmax),dcda(nlmax),dcdb(nlmax),dcdrh(nlmax),dcdxi(nlmax), y1
     *(nlmax),dy1dr(nlmax),y2(nlmax),dy2dr(nlmax),dum(3345)
      common/stdep/ls,maxlayr
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
c Ying      
c      dimension xl1(0:600),xl2(0:600),d1(0:600),d2(0:600),k1(0:600),k2(0
c     *:600)
      dimension xl1(0:4800),xl2(0:4800),d1(0:4800),d2(0:4800),k1(0:4800),k2(0
     *:4800)
c Ying      
c      dimension xlpred(0:600),dlpred(0:600)
      dimension xlpred(0:4800),dlpred(0:4800)
      character*2 ichar(2)
c      data nmx/600/,inss/5/,ichar/' R',' L'/,zero/0./,null/0/
      data nmx/4800/,inss/5/,ichar/' R',' L'/,zero/0./,null/0/
      chz=2.0d0*pi
      stepf=1.d0
      print *,'enter eps'
      read(*,*) eps
      eps1=eps
      eps2=eps
      modtot=0
      call steps(eps)
      print *,'frequency spacing is computed as 1.0/timelength'
      print *,'if used for synthetics, time length must be dt * 2**n'
      print *,'enter power of two (e.g. 8):'
      read *,npow
      print *,'enter dt (sec):'
      read *,dt
      ttlen=dt*(2**npow)
      tlen=ttlen
      print *,'enter upper frequency limit (0 for Nyquist) in Hz:'
      read *,fnyquist
      if(.not.(fnyquist.le.0.))goto 23000
      fnyquist=1.0d0/(2.0d0*dt)
23000 continue
      print *,'enter highest mode number (fundamental=0):'
      read *,nbran
      print *,'enter min, max Vphase (km/s):'
      read(*,*) cmin,cmax
      print *,'the modes file contains eigenvectors at the top of the mo
     *del'
      print *,'too many layers will use up lots of disk space'
      print *,'how many layers (from top) to eigenfunction file?'
      read *,maxlayr
      df=1.0d0/ttlen
      nfreq=fnyquist/df
      write(iout,80) cmin,cmax,df,nfreq
80    format(/,'Phase velocity bounds:',2g12.4,' km/s',//, 'Frequency sp
     *acing',g14.6,' Hz,  number of fft bins:',i6)
      write(iout,100) eps,eps1
100   format(/,'Integration precision =',g12.4,'  root precision =', g12
     *.4,///,3x,'mode',5x,'order', 7x,'f(hz)',10x,'t(secs)',4x,'phase vel(km/s)',6x,'grp vel(km/s)', 8x,'q',13x,'raylquo',/) 
      print *,'Frequency step=',df,' Hz'
      print *,'Highest frequency=',nfreq*df,' Hz'
      print *,'When used with newtons, this modesfile will be suitable f
     *or'
      print *,'the inversion of time series of length',tlen,tlen/2.,tlen
     */4.
      print *,'or',tlen/8.,' seconds. Change time length if necessary'
      wmin=df*chz
      wmax=nfreq*wmin
      if(.not.(nbran.lt.0.or.nbran.gt.nmx))goto 23002
      nbran=nmx
23002 continue
      open(90,file='layer.love')
      write(90,*)'skipped  w  vphase gcom number-layer stratl endl'
      open(10,file='dummyeig',form='unformatted')
      read(10) rn4,in4,n4
      df4=df
      cmax4=cmax
      rewind ioeig
      write(ioeig) rn4,in4,n4,nfreq,nbran,df4,cmax4
      do23004 i=1,n4 
      read(10) rn4,vp4,vs4,rh4,q4
      write(ioeig) rn4,vp4,vs4,rh4,q4
23004 continue
23005 continue
      close(10)
      write(6,*) 'Computation progress'
      write(6,*) 'If mmax exceeds nbran=',nbran,' it is truncated'
      write(6,*)
c chukren
c      write(6,10)
10    format(5x,'i',10x,'Hz  mmin  mmax')
      do23006 i=1,nfreq
      xlpred(i)=0.0
23006 continue
23007 continue
      do23008 i=1,nfreq
      dlpred(i)=0.0
23008 continue
23009 continue
c Ying       
c      do23010 ifreq=1,nfreq 
      do23010 ifreq=1,nfreq/2 
      wrad=ifreq*wmin
      fhz=wrad/chz
      flmin=6371.0d0*wrad/cmax - 0.5d0
c     if(.not.(flmin.lt.10.0d0))goto 23012
c     flmin=10.0d0
c Ying
      if(.not.(flmin.lt.2.0d0))goto 23012
      flmin=2.0d0
23012 continue
      flmax=6371.0d0*wrad/cmin - 0.5d0
c      if(.not.(flmax.lt.10.0d0))goto 23014
c      flmax=10.0d0
c  Ying
      if(.not.(flmax.lt.2.0d0))goto 23014
      flmax=2.0d0
23014 continue
      if(.not.(flmax.le.flmin))goto 23016
      goto 23010
23016 continue
      knsw=1
      maxo=inss
      call ldetqn(flmin,wrad,kmax,dmax,0,nerr)
      if(.not.(nerr.gt.0))goto 23018
      print *,'ERROR in ldetqn, frequency=',flmin
      kmax=mmax
23018 continue
      call ldetqn(flmax,wrad,kmin,dmin,0,nerr)
      if(.not.(nerr.gt.0))goto 23020
      print *,'ERROR in ldetqn, frequency=',flmax
      kmin=mmin
23020 continue
      mmax=kmax
      mmin=kmin+1
c chukren
c      write(6,20) ifreq,fhz,mmin,mmax
20    format(i6,f12.6,2i6)
      if(.not.(mmax.lt.mmin))goto 23022
      goto 23010
23022 continue
      if(.not.(mmax.gt.nbran))goto 23024
      mmax=nbran
23024 continue
      do23026 m=mmin,mmax 
      xl1(m)=flmin
      k1(m)=kmax
      d1(m)=dmax
      xl2(m)=flmax
      k2(m)=kmin
      d2(m)=dmin
23026 continue
23027 continue
      do23028 m=mmin,mmax 
23030 if(.not.(k1(m).ne.k2(m)+1))goto 23031
      fltry=0.5d0*(xl1(m)+xl2(m))
      call ldetqn(fltry,wrad,ktry,dtry,0,nerr)
      if(.not.(nerr.gt.0))goto 23032
      print *,'WARNING: nerr=1 on ldetqn, ktry may be corrupt'
23032 continue
      do23034 mm=ktry+1,mmax 
      if(.not.(xl2(mm).gt.fltry))goto 23036
      xl2(mm)=fltry
      k2(mm)=ktry
      d2(mm)=dtry 
23036 continue
23034 continue
23035 continue
      do23038 mm=mmin,ktry 
      if(.not.(xl1(mm).lt.fltry))goto 23040
      xl1(mm)=fltry
      k1(mm)=ktry
      d1(mm)=dtry 
23040 continue
23038 continue
23039 continue
      goto 23030
23031 continue
23028 continue
23029 continue
      do23042 m=mmin,mmax 
      knsw=0
      maxo=8
c Ying      
c      if(.not.(xl2(m).lt.10.0d0))goto 23044
      if(.not.(xl2(m).lt.2.0d0))goto 23044
      goto 23042
23044 continue
      call rootf(wrad,flroot,detroot,xl1(m),xl2(m),d1(m),d2(m),eps)
c Ying      
c      if(.not.(flroot.lt.10.0d0))goto 23046
      if(.not.(flroot.lt.2.0d0))goto 23046
      goto 23042
23046 continue
      call detqn(wrad,kroot,droot,1,nerr)
      wdiff=(wrad-wray*wn)/wrad
      gcom=vn*cg/1000.d0
      qmod=0.0
      if(.not.(qinv .gt. 0.0))goto 23048
      qmod=1.0/qinv
23048 continue
      vphase=6371.d0*wrad/(fl+0.5d0)
      w4=wrad
      ekin=1.0
      fl4=fl
      sqll=sfl3
      ls00=max(ls,nsl-maxlayr)
      ls0=ls00
      m1=m+1
      kind=3-mtype
      if(.not.(nerr.gt.0))goto 23050
      write(iout,190) m,ichar(mtype),flroot,fhz,1.0/fhz,vphase,gcom,qmod
     *,wdiff,'warn: detqn error'
23050 continue
      if(.not.(wdiff.gt.0.001))goto 23052
      write(iout,190) m,ichar(mtype),flroot,fhz,1.0/fhz,vphase,gcom,qmod
     *,wdiff,'skip: inaccurate'
190   format(i5,a2,f10.2,6g16.7,1x,a)
      goto 23042
23052 continue
      isig=+1
      if(.not.(y1(ls).lt.0.))goto 23054
      isig=-1
      do23056 i=ls00,n 
      y1(i)=-y1(i)
      dy1dr(i)=-dy1dr(i)
      if(.not.(kind.eq.1))goto 23058
      goto 23056
23058 continue
      y2(i)=-y2(i)
      dy2dr(i)=-dy2dr(i)
23056 continue
23057 continue
23054 continue
      write(iout,200) m,ichar(mtype),flroot,fhz,1.0/fhz,vphase,gcom,qmod
     *,wdiff,isig
200   format(i5,a2,f10.2,6g16.7,i5)
      gcom=gcom*1000.0
      vphase=vphase*1000.0
      if(.not.(qmod.gt.0.))goto 23060
      qmod=1.0/qmod
23060 continue
      write(ioeig) w4,vphase,gcom,qmod,ekin,wdiff,fl4,sqll,ls0,m1,kind
      if(.not.(kind.eq.1))goto 23062
c      write(*,*)
c      write(*,*)'love engifunction : from to:ls00= n= ',ls00,n
      number = n-ls00 +1
c Ying      
       if(.not.(wdiff.gt.0.001)) write(90,91)w4,vphase,gcom,number,ls00,n
91    format(1x,e12.4,1x,f12.5,1x,f12.5,1x,i5,i5,i5)      
      do23064 i=ls00,n
      write(ioeig) y1(i),dy1dr(i),dcda(i),dcdb(i),dcdrh(i),dcdxi(i)
23064 continue
23065 continue
      goto 23063
23062 continue
c      write(*,*)
c      write(*,*)'rayleigh engifunction : from to:ls00= n= ',ls00,n
      number = n-ls00 +1
c Ying      
       if(.not.(wdiff.gt.0.001)) write(90,91)w4,vphase,gcom,number,ls00,n
      do23066 i=ls00,n
      write(ioeig) y1(i),dy1dr(i),y2(i),dy2dr(i),dcda(i),dcdb(i),dcdrh(i
     *),dcdxi(i)
23066 continue
23067 continue
23063 continue
23042 continue
23043 continue
23010 continue
23011 continue
      write(ioeig) zero,zero,zero,zero,zero,zero,zero,zero,null,null,nul
     *l
      close(90)
      return
      end
      subroutine ldetqn(fldum,wdum,kdum,ddum,ifeif,nerr)
      implicit real*8(a-h,o-z)
      real*4 f4,d4
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl, fl1,fl
     *2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
c      data nmx/600/,inss/5/
      data nmx/4800/,inss/5/
      l=fl
      fl=fldum
      fl1=fl+1.d0
      fl2=fl+fl1
      fl3=fl*fl1
      sfl3=dsqrt(fl3)
      call detqn(wdum,kdum,ddum,ifeif,nerr)
      f4=fl
      d4=ddum
      return
      end
      subroutine rootf(wrad,flroot,detroot,x1,x2,d1,d2,tol)
      implicit real*8 (a-h,o-z)
      parameter (itmax=100,eps=3.e-8)
      a=x1
      b=x2
      fa=d1
      fb=d2
      if(.not.(fb*fa.gt.0.))goto 23068
      stop 'root must be bracketed for rootf.'
23068 continue
      fc=fb
      do23070 iter=1,itmax 
      if(.not.(fb*fc.gt.0.))goto 23072
      c=a
      fc=fa
      d=b-a
      e=d
23072 continue
      if(.not.(abs(fc).lt.abs(fb)))goto 23074
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
23074 continue
      tol1=2.*eps*abs(b)+0.5*tol
      xm=.5*(c-b)
      if(.not.(abs(xm).le.tol1 .or. fb.eq.0.))goto 23076
      flroot=b
      return
23076 continue
      if(.not.(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)))goto 23078
      s=fb/fa
      if(.not.(a.eq.c))goto 23080
      p=2.*xm*s
      q=1.-s
      goto 23081
23080 continue
      q=fa/fc
      r=fb/fc
      p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
      q=(q-1.)*(r-1.)*(s-1.)
23081 continue
      if(.not.(p.gt.0.))goto 23082
      q=-q
23082 continue
      p=abs(p)
      if(.not.(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))))goto 23084
      e=d
      d=p/q
      goto 23085
23084 continue
      d=xm
      e=d
23085 continue
      goto 23079
23078 continue
      d=xm
      e=d
23079 continue
      a=b
      fa=fb
      if(.not.(abs(d) .gt. tol1))goto 23086
      b=b+d
      goto 23087
23086 continue
      b=b+sign(tol1,xm)
23087 continue
      call ldetqn(b,wrad,kroot,detroot,0,nerr)
      fb=detroot
23070 continue
23071 continue
      print *, 'rootf exceeding maximum iterations.'
      flroot=b
      return
      end
      subroutine partials(y,dw,i)
      parameter(nlmax=700)
      implicit real*8 (a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl, fl1,fl
     *2,fl3,sfl3,nord,ll,kount,knsw,ifanis,iback
      common r(nlmax),fmu(nlmax),flam(nlmax),qshear(nlmax),qkappa(nlmax), xa2(nlmax)
     *,xlam(nlmax),rho(nlmax),qro(3,nlmax),g(nlmax),qg(3,nlmax), fcon(nlmax),fspl(3,
     *nlmax),lcon(nlmax),lspl(3,nlmax),ncon(nlmax), nspl(3,nlmax),ccon(nlmax),cspl(3
     *,nlmax),acon(nlmax),aspl(3,nlmax)
      dimension y(4),dw(4)
      a=acon(i)
      f=fcon(i)
      c=ccon(i)
      xl=lcon(i)
      xn=ncon(i)
      d=rho(i)
      z=r(i)
      z2=z*z
      y1=y(1)
      y3=y(3)
      y1r=y1/z
      y3r=y3/z
      t1=2.0d0*y1-fl3*y3
      y2=c*y(2)+f*t1/z
      y4=xl*(y(4)+y1r-y3r)
      vsv=xl/d
      if(.not.(xl.gt.0.0d0))goto 23088
      vsv=dsqrt(vsv)
23088 continue
      vph=sqrt(a/d)
      t2=t1*t1
      t3=a-xl-xl
      dw(3)=-wsq*d*z2*(y1*y1+fl3*y3*y3)+z2*y2*y2/c+ (a-xn-f*f/c)*t2+(fl-
     *1.0d0)*fl3*(fl+2.0d0)*xn*y3*y3
      if(.not.(xl.gt.0.0d0))goto 23090
      dw(3)= fl3*z2*y4*y4/xl+dw(3)
23090 continue
      dw(3)=dw(3)-2.0d0*d*z*y1*g(i)*t1
      dw(3)=0.5d0*w*dw(3)/d
      dw(1)=((z*y2+2.0d0*f*xl*t1/t3)**2)/c+a*t2*(1.0d0-f*f*a/(c*t3*t3))
      dw(1)=w*dw(1)/vph
      if(.not.(xl.gt.0.0d0))goto 23092
      dw(2)=fl3*z2*y4*y4/xl-4.0d0*f*xl*z*y(2)*t1/t3+2.0d0*xn*(fl3*y3*(y1
     *- y3)-y1*t1)
      dw(2)=w*dw(2)/vsv
      goto 23093
23092 continue
      dw(2)=0.0d0
23093 continue
      dw(4)=w*xn*(fl3*y3*(y1-y3)-y1*t1)
      return
      end
