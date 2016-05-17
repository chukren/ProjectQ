program love

!     this is an adaptation of Guy Masters program tor1.f

!     differences are:
! (1) modes are computed at fixed frequency for non-integer angular
!     order
! (2) output is conform that given by surfc.f - at this stage anisotropy
!     is ignored (although it is still used in the computations, the
!     model output and partial derivatives are for vph and vsv)
!     the eigenvector scales to produce unit kinetic energy as
!     defined by Dahlen, GJ, 59:19-42,1979, and unlike Aki&Richards,
!     who omit a factor l(l+1) under the integral. Units are MKS on
!     output (8/92 G.N.).

! compile:  f77 love.f ltable.r
! source text of several subroutines is ratfor, in ltable.r

      character*256 filnam
      print *,'enter input model file: '
      read(*,'(a256)') filnam
      open(7,file=filnam,status='old',form='formatted')
      print *,'enter output file: '
      read(*,'(a256)') filnam
      open(8,file=filnam,form='formatted')
      print *,'enter eigfn output file: '
      read(*,'(a256)') filnam
      open(9,file='dummyeig',form='unformatted')
      call model(7,8,9)
      close(7)            
      close(9)
      open(9,file=filnam,form='unformatted')
      call ltable(8,9,2)
      close(9)
      close(8)
      stop
end

subroutine detqn(wdim,knt,det,ifeif,nerr)
      parameter(nlmax=700)
!**** supervises the integration of the equations,it returns the value
!**** of the secular determinant as det and the count of zero crossings.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      real*4 dcda,dcdb,dcdxi,dcdrh,y1,dy1dr,dum,y2,dy2dr
      common/eifx/a(5,nlmax),dcda(nlmax),dcdb(nlmax),dcdrh(nlmax),dcdxi(nlmax),
     +      y1(nlmax),dy1dr(nlmax),y2(nlmax),dy2dr(nlmax),dum(3345)
      common r(nlmax),fmu(nlmax),qshear(nlmax),qkappa(nlmax),
     + rho(nlmax),qro(3,nlmax),lcon(nlmax),lspl(3,nlmax),ncon(nlmax),nspl(3,nlmax)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/stdep/ls,maxlayr
      dimension ass(6)
!  nerr deals with negative group velocities (which violate assumptions)
      nerr=0
      w=wdim/wn
      wsq=w*w
      kount=0
      fct=0.d0
      if(tref.gt.0.d0) fct=2.d0*dlog(tref*wdim)/pi 
      ass(1)=1.d0
      do 5 ind=2,6
    5   ass(ind)=0.d0
      q=0.d0
      ls=nocp1
      call startl(nocp1,nsl,fmu,ls,q)
      if(ls.ne.nocp1) call tps(ls,ass,ifeif)
      if(ifeif.eq.1) goto 10
      call tprop(ls,nsl,ass)
      det=ass(2)/dsqrt(ass(1)*ass(1)+ass(2)*ass(2))
      if(knsw.ne.1) return
      knt=kount-2             
      if(l.eq.1) return
      knt=knt+1
      irem=mod(knt,2)
      if(irem.eq.0.and.det.lt.0.d0) return
      if(irem.ne.0.and.det.gt.0.d0) return
      knt=knt+1
      return
   10 call etprop(ls,nsl,ass)
      if(ass(3).lt.0.) then
          ass(3)=-ass(3)
          nerr=1
      endif
      if(ass(6)/ass(3).lt.0.) then
          ass(6)=-ass(6)
          nerr=nerr+2
      endif
      det=ass(2)/dsqrt(ass(1)*ass(1)+ass(2)*ass(2))
      rnrm=1.d0/(w*dsqrt(ass(3)))
      cg=(fl+0.5d0)*ass(4)/(w*ass(3))
      if(nerr.eq.1.or.nerr.eq.3) cg=-cg
      qinv=ass(5)/(wsq*ass(3))
      wray=dsqrt(ass(6)/ass(3))
      if(nerr.eq.1.or.nerr.eq.2) wray=-wray
!  compute partials (isotropy assumed, follows Takeuchi and Saito 1972)
      do 11 i=1,n
      y1(i)=0.
      dy1dr(i)=0.
      dcdrh(i)=0.
      dcdxi(i)=0.
      dcda(i)=0.
11    dcdb(i)=0.
      ekin=wsq*ass(3)
      vf=w/(fl+0.5d0)
      x1=vf*vf/(ekin*cg)
      do 12 i=ls,nsl
      z1=a(1,i)**2
      z2=a(2,i)**2
      rr=r(i)**2
      dcdrh(i)=0.5d0*x1*(-wsq*rho(i)*rr*z1+rr*z2/lcon(i)+
     +  (fl-1.0d0)*(fl+2.0d0)*ncon(i)*z1)/rho(i)
      dcda(i)=0.
      dcdb(i)=x1*(rr*z2/lcon(i)+(fl-1.0d0)*(fl+2.0d0)*ncon(i)*z1)*
     +  sqrt(rho(i)/lcon(i))
      dcdxi(i)=0.5*x1*(fl-1.0)*(fl+2.0)*ncon(i)*z1
12    continue
!   transform a2 to da1/dr
      do 15 i=ls,nsl
   15   a(2,i)=a(1,i)/r(i)+a(2,i)/(lcon(i)*(1.d0+qshear(i)*fct))
!   normalize to unite kinetic energy (but in funny units)
      do 20 i=ls,nsl
        do 20 j=1,2
   20     a(j,i)=a(j,i)*rnrm
!   normalize to MKS, such that omega**2 x I1 = 1.0 Nm
!   note I1 includes factor fl3 (as in Dahlen, GJ, 59:19-43, 1979 and
!   unlike Aki & Richards). Remove division by sfl3 to conform to A&R.
!   note dcdb is in units 1/rn since it is an integrand
!   the factors below have been obtained as follows:
!   7.78d-10=1/sqrt(wn*wn*rhobar*rn)
!   1.22d-16=7.78d-10/rn
!   1.94d-7=vn/(rhobar*rn)
!   1.56d-7=1/rn
      do  i=ls,nsl
          y1(i)   =a(1,i)*7.78808387d-10/sfl3
          dy1dr(i)=a(2,i)*1.22242723d-16/sfl3
          dcdrh(i)=dcdrh(i)*1.949574d-7
          dcdb(i) =dcdb(i)*1.569612d-7
          dcdxi(i)=dcdxi(i)*1.0752d-3
      enddo 
      
      return
end

      subroutine tps(i,a,ifeif)
      parameter(nlmax=700)
!*** toroidal mode start soln using sph bessel fns.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(nlmax),fmu(nlmax),qshear(nlmax),qkappa(nlmax),
     + rho(nlmax),qro(3,nlmax),lcon(nlmax),lspl(3,nlmax),ncon(nlmax),nspl(3,nlmax)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension a(2)
      fu=fmu(i)*(1.d0+qshear(i)*fct)
      sqk=wsq*rho(i)/fu
      zsq=r(i)*r(i)*sqk
      call bfs(fl,zsq,eps,fp)
      a(1)=r(i)
      a(2)=fu*(fp-1.d0)       
!      if(ifeif.ne.1) return
!      qrmu=fu*qshear(i)
!      z=dsqrt(zsq)
!      fnor=r(i)*r(i)/dsqrt(sqk)
!      fi1=bfint(l,z)
!      fi2=z*(z*z+fp*(fp+1)-fl3)/2.d0
!      fi3=z*fp+fi2-fl3*fi1
!      fi4=fi3-z+fi1
!      a(3)=fi2*rho(i)*fnor/sqk
!      a(4)=fi1*ncon(i)*fnor
!      a(5)=qrmu*(fi4+(fl3-1)*fi1)*fnor
!      a(6)=(lcon(i)*(fi4+fi1)+(fl3-2)*ncon(i)*fi1)*fnor
      return
      end

      subroutine etprop(jf,jl,f)
      parameter(nlmax=700)
!*** propagates f from jf to jl - toroidal modes ***
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(nlmax),fmu(nlmax),qshear(nlmax),qkappa(nlmax),
     + rho(nlmax),qro(3,nlmax),lcon(nlmax),lspl(3,nlmax),ncon(nlmax),nspl(3,nlmax)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/a(5,nlmax),kdum(5129)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension h(6,10),s(6),f(6)
      fl3m2=fl3-2.d0
      maxo1=maxo-1
      y=r(jf)
      qy=1.d0/y+dsqrt(dabs(rho(jf)*wsq/fmu(jf)-fl3/(y*y)))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2) 
      if(i.eq.jl) return
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qll=1.d0+qshear(iq)*fct
      qx=qy
      qy=1.d0/y+dsqrt(dabs(rho(i)*wsq/fmu(i)-fl3/(y*y)))
      q=dmax1(qx,qy)
      del=step(maxo)/q
      dxs=0.d0   
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q,maxo1)
        dxs=dx
        do 20 nf=1,6
   20     s(nf)=f(nf)
        do 40 ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
          nn=ll
          if(ifanis.ne.0) nn=(ncon(iq)+
     +        t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          h(3,ni)=ro*f(1)*f(1)*z*z
          z=1.d0/z
          h(1,ni)=z*f(1)+f(2)/ll
          h(2,ni)=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
          t1=(h(1,ni)/z-f(1))**2
          t2=fl3m2*f(1)*f(1)
          h(4,ni)=nn*f(1)*f(1)
          h(5,ni)=(t1+t2)*ll*qshear(iq)
          h(6,ni)=ll*t1+t2*nn
   40     call rkdot(f,s,h,6,ni)
        x=y
        if(y.ne.r(i)) goto 15
      goto 10
      end          

      subroutine tprop(jf,jl,f)
      parameter(nlmax=700)
!*** propagates f from jf to jl - toroidal modes ***
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(nlmax),fmu(nlmax),qshear(nlmax),qkappa(nlmax),
     + rho(nlmax),qro(3,nlmax),lcon(nlmax),lspl(3,nlmax),ncon(nlmax),nspl(3,nlmax)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/a(5,nlmax),kdum(5129)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension h(2,10),s(2),f(2)
      fl3m2=fl3-2.d0
      maxo1=maxo-1
      y=r(jf)
      qy=1.d0/y+dsqrt(dabs(rho(jf)*wsq/fmu(jf)-fl3/(y*y)))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.eq.jl) return                                     
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qll=1.d0+qshear(iq)*fct
      qx=qy
      qy=1.d0/y+dsqrt(dabs(rho(i)*wsq/fmu(i)-fl3/(y*y)))
      q=dmax1(qx,qy)
      del=step(maxo)/q
      dxs=0.d0   
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q,maxo1)
        dxs=dx
        s(1)=f(1)
        s(2)=f(2)
        do 40 ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
          nn=ll
          if(ifanis.ne.0) nn=(ncon(iq)+
     +      t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          z=1.d0/z
          h(1,ni)=z*f(1)+f(2)/ll
          h(2,ni)=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
   40     call rkdot(f,s,h,2,ni)
        if(knsw.eq.1) then
          if(s(2)*f(2).le.0.d0) then
            if(f(2).eq.0.d0) then
              tes=-s(2)*s(1)
            else
              tes=f(2)*s(1)-s(2)*f(1)
            endif
            if(tes.lt.0.d0) kount=kount+1
            if(tes.gt.0.d0) kount=kount-1
          end if
        end if
        x=y
        if(y.ne.r(i)) goto 15
      goto 10
      end          

      subroutine model(iin,iout,ioeig)
      parameter(nlmax=700)
      implicit real*8(a-h,o-z)
      integer*4 ititle(20),n4,in4
      real*4 rn4,rh4,vp4,vs4,q4
      real*8 lcon,ncon,lspl,nspl
      common r(nlmax),fmu(nlmax),qshear(nlmax),qkappa(nlmax),
     + rho(nlmax),qro(3,nlmax),lcon(nlmax),lspl(3,nlmax),ncon(nlmax),nspl(3,nlmax)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/vpv(nlmax),vph(nlmax),vsv(nlmax),vsh(nlmax),eta(nlmax),
     + 	wrk(2564),kdum
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      data bigg,tau/6.6723d-11,1.d3/,rhobar/5515.d0/
      pi=3.14159265358979d0
      read(iin,100) (ititle(i),i=1,20)
  100 format(20a4)
      read(iin,*) ifanis,tref,ifdeck
      if(ifdeck.eq.0) go to 1000
!*** card deck model ***
      read(iin,*) n,nic,noc
! added if-block to allow for isotropic models too. SL Aug 93
      if (ifanis.eq.0) then
	 do 101 i=1,n
101      read(iin,*) r(i),rho(i),vpv(i),vsv(i),qkappa(i),qshear(i)
      else
	 do 102 i=1,n
102      read(iin,*) r(i),rho(i),vpv(i),vsv(i),
     +     qkappa(i),qshear(i),vph(i),vsh(i),eta(i)
      endif
      go to 2000
!*** polynomial model ***
 1000 read(iin,*) nreg,nic,noc,rx
      rx=rx*tau
      n=0
      knt=0
      jj=5
      if(ifanis.ne.0) jj=8
      do 10 nn=1,nreg
      read(iin,*) nlay,r1,r2
      r1=r1*tau
      r2=r2*tau
      dr=(r2-r1)/float(nlay-1)
      do 15 i=1,nlay
      n=n+1
   15 r(n)=r1+dr*float(i-1)
      do 20 j=1,jj
      read(iin,*) (wrk(i),i=1,5)
      do 20 i=1,nlay
      ind=knt+i
      rt=r(ind)/rx
      val=wrk(1)+rt*(wrk(2)+rt*(wrk(3)+rt*(wrk(4)+rt*wrk(5))))
      if(j.eq.1) rho(ind)=val*tau
      if(j.eq.2) vpv(ind)=val*tau
      if(j.eq.3) vsv(ind)=val*tau
      if(j.eq.4) qkappa(ind)=val
      if(j.eq.5) qshear(ind)=val
      if(ifanis.eq.0) goto 20
      if(j.eq.6) vph(ind)=val*tau
      if(j.eq.7) vsh(ind)=val*tau
      if(j.eq.8) eta(ind)=val
   20 continue
   10 knt=knt+nlay
 2000 if(ifanis.ne.0) go to 3000
      do 25 i=1,n
      vph(i)=vpv(i)
      vsh(i)=vsv(i)
   25 eta(i)=1.d0
 3000 if(iout.lt.0) goto 30
!*** write out model ***
      write(iout,900) (ititle(k),k=1,20),tref
  900 format(1x,20a4,' ref per =',f6.1,' secs',///,2x,'level',
     1 4x,'radius(m)',5x,'rho(kg/m3)',2x,'vpv(m/s)',4x,'vph(m/s)',4x,
     2 'vsv(m/s)',4x,'vsh(m/s)',4x,'eta',9x,'qmu ',8x,'qkap',/)
      write(iout,905) (i,r(i),rho(i),vpv(i),vph(i),vsv(i),vsh(i),
     1 eta(i),qshear(i),qkappa(i),i=1,n)
  905 format(3x,i3,f12.1,5f12.2,f12.5,2f12.2)
  30  continue
      if(r(n).lt.6.e6.or.rho(n).lt.500.0.or.vpv(n).lt.100.0) then
	write(6,*) 'WARNING!!!! Is your input in MKS units???'
	write(iout,*) 'WARNING!!!! Is your input in MKS units???'
! if not MKS units --> make MKS units:  (SL Aug 93)
	do 32 i=1,n
	r(i) = r(i)*tau
	rho(i) = rho(i)*tau
	vpv(i) = vpv(i)*tau
	vsv(i) = vsv(i)*tau
	vph(i) = vph(i)*tau
	vsh(i) = vsh(i)*tau
32      continue
	write(6,*) 'corrected for WARNING as follows:'
	write(6,*) 'r,rho,vpv,vsv,vph,vsh multiplied by 1000 to'
	write(6,*) 'convert to MKS units.'
	write(iout,*) 'corrected for WARNING as follows:'
	write(iout,*) 'r,rho,vpv,vsv,vph,vsh multiplied by 1000 to'
	write(iout,*) 'convert to MKS units.'
      endif
   60 nsl=n
      if(vsv(nsl).gt.0.d0) go to 70
   65 nsl=nsl-1
      if(vsv(nsl).le.0.d0) go to 65
   70 nicp1=nic+1
      nocp1=noc+1
      nslp1=nsl+1
      n4=n
      in4=nsl
      rn4=r(n)
      write(ioeig) rn4,in4,n4
      do 75 i=1,n
      rn4=r(i)
      vp4=vph(i)
      vs4=vsv(i)
      rh4=rho(i)
      q4=qshear(i)
75    write(ioeig) rn4,vp4,vs4,rh4,q4
!*** normalise and spline ***
      rn=r(n)
      gn=pi*bigg*rhobar*rn
      vn2=gn*rn
      vn=dsqrt(vn2)
      wn=vn/rn
      do 45 i=1,n
      r(i)=r(i)/rn
      if(i.gt.1.and.dabs(r(i)-r(i-1)).lt.1.d-7) r(i)=r(i-1)
      if(qshear(i).gt.0.d0) qshear(i)=1.d0/qshear(i)
      if(qkappa(i).gt.0.d0) qkappa(i)=1.d0/qkappa(i)
      rho(i)=rho(i)/rhobar
      lcon(i)=rho(i)*vsv(i)*vsv(i)/vn2
      ncon(i)=rho(i)*vsh(i)*vsh(i)/vn2
   45 fmu(i)=lcon(i)
      call drspln(1,n,r,rho,qro,wrk)
      call drspln(1,n,r,lcon,lspl,wrk)
      if(ifanis.eq.0) goto 80
      call drspln(1,n,r,ncon,nspl,wrk)
80      tref=0.5d0*tref/pi
      return
      end        


      subroutine startl(jf,jl,v,ls,q)
      parameter(nlmax=700)
!*** finds start level between jf and jl using velocityv and ang. ord. l.
!*** upon entry q is the value of the exponent at r(jf) or at the turning
!*** point(q=0) depending on previous calls to startl. upon exit q is the
!*** value of the exponent at the starting level ls.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(nlmax),fmu(nlmax),qshear(nlmax),qkappa(nlmax),
     + rho(nlmax),qro(3,nlmax),lcon(nlmax),lspl(3,nlmax),ncon(nlmax),nspl(3,nlmax)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension rrlog(nlmax),p(nlmax),v(nlmax)
      save rrlog,p,ifirst
      data ifirst/1/
      if(ifirst.ne.1) goto 5
      ifirst=0
      vertno=-dlog(eps)
      do 1 i=3,n
    1 rrlog(i)=.5d0*dlog(r(i)/r(i-1))
    5 do 10 j=jf,jl
      pp=fl3-wsq*r(j)*r(j)*rho(j)/v(j)
      if(pp.le.0.d0) goto 15
   10 p(j)=dsqrt(pp)
   15 p(j)=0.d0
   20 k=j
      j=j-1
      if(j.le.jf) go to 25
      q=q+rrlog(k)*(p(j)+p(k))
      if(q.lt.vertno) go to 20
      ls=j
      return
   25 ls=jf
      return
      end

      subroutine baylis(q,maxo1)
!    baylis returns the coefficients for rks integration.
!    see e. baylis shanks(1966 a. m. s.) and references therein for the
!    coefficients. the eight runge-kutta-shanks formulae are (1-1) (2-2)
!    (3-3) (4-4) (5-5) (6-6) (7-7) (8-10). for orders greater than 4 the
!    formulae are approximate rather than exact so incurring less roundoff.
      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,i
      ds=q*dabs(dx)
      do 10 j=1,maxo1
      if(ds.gt.step(j)) go to 10
      i=j
      go to 15
   10 continue
      i=maxo
   15 c(1)=0.d0
      go to (1,2,3,4,5,6,7,8),i
    1 b(1)=dx
      return
    2 c(2)=dx
      b(1)=dx
      b(2)=.5d0*dx
      b(3)=1.d0
      return
    3 c(2)=.5d0*dx
      c(3)=dx
      b(1)=c(2)
      b(2)=-dx
      b(3)=-2.d0
      b(4)=.16666666666667d0*dx
      b(5)=4.d0
      b(6)=1.d0
      return
    4 c(2)=.01d0*dx
      c(3)=.6d0*dx
      c(4)=dx
      b(1)=c(2)
      b( 2)=-.17461224489790d+02*dx
      b( 3)=-.10343618513324d+01
      b( 4)= .59691275167780d+02*dx
      b( 5)=-.10140620414448d+01
      b( 6)= .30814908546230d-01
      b( 7)=-.25555555555556d+01*dx
      b( 8)=-.11165449632656d+01
      b( 9)=-.22568165070006d+00
      b(10)=-.49077733860351d-01
      return
    5 c( 2)= 1.1111111111111d-04*dx
      c( 3)= 3.0d-01*dx
      c( 4)= 7.5d-01*dx
      c( 5)= dx
      b( 1)=c(2)
      b( 2)=-.40470000000000d+03*dx
      b( 3)=-.10007412898443d+01
      b( 4)= .25301250000000d+04*dx
      b( 5)=-.10004446420631d+01
      b( 6)= .74107010523195d-03
      b( 7)=-.11494333333333d+05*dx
      b( 8)=-.10004929965491d+01
      b( 9)= .52629261224803d-03
      b(10)=-.12029545422812d-03
      b(11)= .92592592592593d-01*dx
      b(12)= .00000000000000d+00
      b(13)= .47619047619048d+01
      b(14)= .42666666666667d+01
      b(15)= .77142857142857d+00
      return
    6 c(2)=3.3333333333333d-03*dx
      c(3)=.2d0*dx
      c(4)=.6d0*dx
      c(5)=9.3333333333333d-01*dx
      c(6)=dx
      b( 1)=c(2)
      b( 2)=-.58000000000000d+01*dx
      b( 3)=-.10344827586207d+01
      b( 4)= .64600000000000d+02*dx
      b( 5)=-.10216718266254d+01
      b( 6)= .30959752321982d-01
      b( 7)=-.62975802469136d+03*dx
      b( 8)=-.10226149961576d+01
      b( 9)= .24906685695466d-01
      b(10)=-.37737402568887d-02
      b(11)=-.54275714285714d+04*dx
      b(12)=-.10225567867765d+01
      b(13)= .25375487829097d-01
      b(14)=-.31321559234596d-02
      b(15)= .12921040478749d-03
      b(16)= .53571428571429d-01*dx
      b(17)= .00000000000000d+00
      b(18)= .61868686868687d+01
      b(19)= .77777777777778d+01
      b(20)= .40909090909091d+01
      b(21)=-.38888888888889d+00
      return
    7 c(2)=5.2083333333333d-03*dx
      c(3)=1.6666666666667d-01*dx
      c(4)=.5d0*dx
      c(5)=dx
      c(6)=8.3333333333333d-01*dx
      c(7)=dx
      b( 1)=c(2)
      b( 2)=-.25000000000000d+01*dx
      b( 3)=-.10666666666667d+01
      b( 4)= .26166666666667d+02*dx
      b( 5)=-.10421204027121d+01
      b( 6)= .61228682966918d-01
      b( 7)=-.64500000000000d+03*dx
      b( 8)=-.10450612653163d+01
      b( 9)= .51262815703925d-01
      b(10)=-.77519379844961d-02
      b(11)=-.93549382716049d+02*dx
      b(12)=-.10450293206756d+01
      b(13)= .48394546673620d-01
      b(14)=-.11877268228307d-01
      b(15)=-.39590894094358d-03
      b(16)= .35111904761905d+03*dx
      b(17)=-.10446476812124d+01
      b(18)= .52479782656724d-01
      b(19)=-.71200922221468d-02
      b(20)=-.61029361904114d-03
      b(21)= .27463212856852d-02
      b(22)= .46666666666667d-01*dx
      b(23)= .57857142857143d+01
      b(24)= .78571428571429d+01
      b(25)= .00000000000000d+00
      b(26)= b(23)
      b(27)= .10000000000000d+01
      return
    8 c(2)=.14814814814815d0*dx
      c(3)=.22222222222222d0*dx
      c(4)=.33333333333333d0*dx
      c(5)= .5d0*dx
      c(6)=.66666666666667d0*dx
      c(7)=.16666666666667d0*dx
      c(8)=dx
      c(9)=.83333333333333d0*dx
      c(10)=dx
      b( 1)=c(2)
      b( 2)= .55555555555556d-01*dx
      b( 3)= .30000000000000d+01
      b( 4)= .83333333333333d-01*dx
      b( 5)= .00000000000000d+00
      b( 6)= .30000000000000d+01
      b( 7)= .12500000000000d+00*dx
      b( 8)= .00000000000000d+00
      b( 9)= .00000000000000d+00
      b(10)= .30000000000000d+01
      b(11)= .24074074074074d+00*dx
      b(12)= .00000000000000d+00
      b(13)=-.20769230769231d+01
      b(14)= .32307692307692d+01
      b(15)= .61538461538461d+00
      b(16)= .90046296296295d-01*dx
      b(17)= .00000000000000d+00
      b(18)=-.13881748071980d+00
      b(19)= .24832904884319d+01
      b(20)=-.21182519280206d+01
      b(21)= .62467866323908d+00
      b(22)=-.11550000000000d+02*dx
      b(23)=-.35064935064935d+00
      b(24)= .50389610389610d+01
      b(25)=-.28398268398268d+01
      b(26)= .52813852813853d+00
      b(27)=-.34632034632035d+01
      b(28)=-.44097222222222d+00*dx
      b(29)=-.14173228346457d+00
      b(30)= .53385826771654d+01
      b(31)=-.35905511811023d+01
      b(32)= .70866141732284d-01
      b(33)=-.45354330708661d+01
      b(34)=-.31496062992126d-01
      b(35)= .18060975609756d+01*dx
      b(36)=-.54692775151925d-01
      b(37)= .47967589466576d+01
      b(38)=-.22795408507765d+01
      b(39)= .48615800135044d-01
      b(40)=-.34031060094530d+01
      b(41)=-.40513166779204d-01
      b(42)= .48615800135044d+00
      b(43)= .48809523809524d-01*dx
      b(44)= .65853658536585d+00
      b(45)= .66341463414634d+01
      b(46)= .52682926829268d+01
      i=10
      return
      end

      subroutine steps(eps)
!*** computes 8 dimensionless step sizes for rks integration
      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      ps=dlog(eps)
      fac=1.d0
      do 2 n=1,8
      fn=n+1
      fac=fac*fn
      x=(dlog(fac)+ps)/fn
      x=dexp (x)
      s=x
      do 1 i=1,n
    1 s=x*dexp(-s/fn)
    2 step(n)=s
      return
      end

      subroutine drspln(i1,i2,x,y,q,f)
      implicit real*8(a-h,o-z)
!   rspln computes cubic spline interpolation coefficients
!   for y(x) between grid points i1 and i2 saving them in q.  the
!   interpolation is continuous with continuous first and second
!   derivitives.  it agrees exactly with y at grid points and with the
!   three point first derivitives at both end points (i1 and i2).
!   x must be monotonic but if two successive values of x are equal
!   a discontinuity is assumed and seperate interpolation is done on
!   each strictly monotonic segment.  the arrays must be dimensioned at
!   least - x(i2), y(i2), q(3,i2), and f(3,i2).  f is working storage
!   for rspln.
!                                                     -rpb
      dimension x(i2),y(i2),q(3,i2),f(3,i2),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0.d0/
      j1=i1+1
      y0=0.d0
!   bail out if there are less than two points total.
      if(i2-i1)13,17,8
 8    a0=x(j1-1)
!   search for discontinuities.
      do 3 i=j1,i2
      b0=a0
      a0=x(i)
      if(a0-b0)3,4,3
 3    continue
 17   j1=j1-1
      j2=i2-2
      go to 5
 4    j1=j1-1
      j2=i-3
!   see if there are enough points to interpolate (at least three).
 5    if(j2+1-j1)9,10,11
!   only two points.  use linear interpolation.
 10   j2=j2+2
      y0=(y(j2)-y(j1))/(x(j2)-x(j1))
      do 15 j=1,3
      q(j,j1)=yy(j)
 15   q(j,j2)=yy(j)
      go to 12
!   more than two points.  do spline interpolation.
 11   a0=0.d0
      h=x(j1+1)-x(j1)
      h2=x(j1+2)-x(j1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
!   calculate derivitive at near end.
      b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
      b1=b0
!   explicitly reduce banded matrix to an upper banded matrix.
      do 1 i=j1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=q(3,i)
 1    b0=f(3,i)
!   take care of last two rows.
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
!   calculate derivitive at far end.
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
!   solve upper banded matrix by reverse iteration.
      do 2 j=j1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
 2    i=k
      q(1,i)=b1
!   fill in the last point with a linear extrapolation.
 9    j2=j2+2
      do 14 j=1,3
 14   q(j,j2)=yy(j)
!   see if this discontinuity is the last.
 12   if(j2-i2)6,13,13
!   no.  go back for more.
 6    j1=j2+2
      if(j1-i2)8,8,7
!   there is only one point left after the latest discontinuity.
 7    do 16 j=1,3
 16   q(j,i2)=yy(j)
!   fini.
 13   return
      end

      subroutine dsplin(n,x,y,q,f)
      implicit real*8(a-h,o-z)
      dimension x(3),y(3),q(3,3),f(3,3),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0.d0/
      a0=0.d0
      j2=n-2
      h=x(2)-x(1)
      h2=x(3)-x(1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
      b0=(y(1)*(h-h2)+y(2)*h2-y(3)*h)/y0
      b1=b0
      do 5 i=1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=q(3,i)
    5 b0=f(3,i)
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
      do 10 j=1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
   10 i=k
      q(1,i)=b1
      do 15 j=1,3
   15 q(j,n)=yy(j)
      return
      end

      subroutine rkdot(f,s,h,nvec,ni)
!*** performs dot product with rks coefficients ***
      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension s(nvec),f(nvec),h(nvec,10)
      goto (1,2,3,4,5,6,7,8,9,10),ni
    1 do 21 j=1,nvec
   21 f(j)=s(j)+b(1)*h(j,1)
      return
    2 do 22 j=1,nvec
   22 f(j)=s(j)+b(2)*(h(j,1)+b(3)*h(j,2))
      return
    3 do 23 j=1,nvec
   23 f(j)=s(j)+b(4)*(h(j,1)+b(5)*h(j,2)+b(6)*h(j,3))
      return
    4 do 24 j=1,nvec
   24 f(j)=s(j)+b(7)*(h(j,1)+b(8)*h(j,2)+b(9)*h(j,3)+b(10)*h(j,4))
      return
    5 do 25 j=1,nvec
   25 f(j)=s(j)+b(11)*(h(j,1)+b(12)*h(j,2)+b(13)*h(j,3)+b(14)*h(j,4)+
     +b(15)*h(j,5))
      return
    6 do 26 j=1,nvec
   26 f(j)=s(j)+b(16)*(h(j,1)+b(17)*h(j,2)+b(18)*h(j,3)+b(19)*h(j,4)+
     +b(20)*h(j,5)+b(21)*h(j,6))
      return
    7 do 27 j=1,nvec
   27 f(j)=s(j)+b(22)*(h(j,1)+b(23)*h(j,3)+b(24)*h(j,4)+b(25)*h(j,5)+
     +b(26)*h(j,6)+b(27)*h(j,7))
      return
    8 do 28 j=1,nvec
   28 f(j)=s(j)+b(28)*(h(j,1)+b(29)*h(j,3)+b(30)*h(j,4)+b(31)*h(j,5)+
     +b(32)*h(j,6)+b(33)*h(j,7)+b(34)*h(j,8))
      return
    9 do 29 j=1,nvec
   29 f(j)=s(j)+b(35)*(h(j,1)+b(36)*h(j,3)+b(37)*h(j,4)+b(38)*h(j,5)+
     +b(39)*h(j,6)+b(40)*h(j,7)+b(41)*h(j,8)+b(42)*h(j,9))
      return
   10 do 30 j=1,nvec
   30 f(j)=s(j)+b(43)*(h(j,1)+h(j,10)+b(44)*(h(j,4)+h(j,6))+
     +b(45)*h(j,5)+b(46)*(h(j,7)+h(j,9)))
      return
      end

      subroutine entry(w,imax,kei)
      implicit real*8(a-h,o-z)
!      common/mtab/we(1200),de(1200),ke(1200),wtry(600),bm(600),um(600)
      common/mtab/we(9600),de(9600),ke(4800),wtry(4800),bm(4800),um(4800)
      call detqn(w,kei,dei,0,nerr)
      indx=min0(max0(2*(kei-ke(1)),1),imax)
      if(indx.eq.1.and.we(1).lt.w) goto 10
      if(indx.eq.imax.and.we(imax).gt.w) goto 10
      if(kei.ne.ke(indx)) goto 5
      if(we(indx).gt.w) goto 10
      indx=indx+1
      if(we(indx).lt.w) goto 10
      return
    5 we(indx)=w
      ke(indx)=kei
      de(indx)=dei
      indx=indx+1
   10 we(indx)=w
      ke(indx)=kei
      de(indx)=dei
      return
      end

subroutine bfs(l,xsq,eps,fp)
!  this routine calculates spherical bessel function of the ist kind.
!  fp is equivalent to (r*dj/dr)/j
!  where r is radius and j is the sbf of order l and argument x=k*r
!  the technique employs the continued fraction approach
!  described in w. lentz's article in applied qptics, vol.15, #3, 1976
      implicit real*8(a-h,o-z)
!  28-8-92: l and lp1 added to real*8 statement (GN)
      real*8 numer,nu,l,lp1
      lp1=l+1
      if(xsq.le.0.d0) goto 10
      x=dsqrt(xsq)
      rx=2.0d0/x
      nu=lp1-0.5d0
      rj=nu*rx
      rx=-rx 
      nu=nu+1
      denom=nu*rx
      numer=denom+1.0d0/rj
      rj=rj*numer/denom
    5   nu=nu+1
        rx=-rx
        a3=nu*rx
        denom=a3+1.d0/denom
        numer=a3+1.d0/numer
        ratio=numer/denom
        rj=rj*ratio
        if(dabs(dabs(ratio)-1.d0).gt.eps) goto 5
      fp=rj*x-lp1
      return
!  series solution
   10 f=1.d0
      fp=l
      a=1.d0
      b=l+lp1
      c=2.d0
      d=l+2.d0
   15   a=-a*xsq/(c*(b+c))
        f=f+a
        fp=fp+a*d
        if(dabs(a*d).lt.eps) goto 20
        c=c+2.d0
        d=d+2.d0
        goto 15
   20 fp=fp/f
      return
      end
