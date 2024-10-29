      program main
      implicit none
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,   hsf0,   tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
      real*8 ebeam,sumsib,sumsigi,sumsigt,sumsibs,sumsigsi,sumsigst,sumsibc,sumsigci,sumsigct,sib,sigi,sigt,del(5,3),urand8
      real*8 zmin,zmax,pt2min,pt2max,pt2up,amhh,amp,amp2,m2th  
      real*8 pi
      integer*8 iy
      integer nev,i,j,ipoint,ii,iz,ipt,ipol,iphi,inev
      character*19 file_name
        data zmin/0.22d0/,zmax/0.77d0/,pt2min/0.015d0/,amhh/0.1395675d0/,amp/0.938272d0/,m2th/1.1645277295317797d0/
      pi=atan(1d0)*4d0
        amp2=amp**2
      iy=12345_ 8
        nev=-2
       open(8,file='bins.txt')
        do ipoint=1,9
        read(8,*)ii,ebeam,q2,x
         if(ipoint.eq.1.or.ipoint.eq.5.or.ipoint.eq.9)then   
        print*,ii,ebeam,q2,x
        do iz=1,3   
        write (file_name,"(i2.2,'.z',i0,'.dat')") ipoint, iz
        open(16,file=trim(file_name))

        z=zmin+5d-1*dble(iz-1)*(zmax-zmin)
        pt2up=(q2*z/2d0/amp/x)**2-amhh**2
        if(amhh**2 + amp2 - (q2**2*z)/(2.*amp2*x**2) - (q2*(-1 + x + z))/x.lt.m2th)
     -   pt2up=(-(amp2**3*x**3) + q2**2*(amhh**2*x*(-1 + z) - m2th*x*z + 
     -       q2*(-1 + x)*(-1 + z)*z) + 
     -    2*amp2**2*x**2*((-amhh**2 + m2th)*x + q2*(-1 + x + z)) - 
     -    amp2*x*((amhh**2 - m2th)**2*x**2 + 
     -       q2**2*(1 + x**2 + 2*x*(-1 + z) + (-3 + z)*z) + 
     -       2*q2*x*(amhh**2*(1 + x - z) + m2th*(-1 + x + z))))/(q2*x*(q2 + 4*amp2*x**2))
c        print*,'pt',pt2up,pt2max,pt2
        pt2max=min(0.999d0*pt2up,1.33d0)
        do ipt=1,10
        pt=sqrt(pt2min+dble(ipt-1)/9d0*(pt2max-pt2min))
      un=1
      le=0d0
      etal=0d0
      etat=1d0
      nev=1e7
      sumsib=0d0
      sumsigi=0d0
      sumsigt=0d0
      sumsibs=0d0
      sumsigsi=0d0
      sumsigst=0d0
      sumsibc=0d0
      sumsigci=0d0
      sumsigct=0d0
      do j=1,nev
      phih=2d0*pi*urand8(iy)
      pheta=2d0*pi*urand8(iy)
      call haprad3(ebeam,iy,1,sib,sigi,sigt) 
c      print*,sigi      
c      print*,sigt      
c
c            stop
      sumsib=sumsib+4d0*pi**2*sib/dble(nev)
      sumsigi=sumsigi+4d0*pi**2*sigi/dble(nev)
      sumsigt=sumsigt+4d0*pi**2*sigt/dble(nev)
      sumsibs=sumsibs+4d0*pi**2*sin(phih-pheta)*sib/dble(nev)
      sumsigsi=sumsigsi+4d0*pi**2*sin(phih-pheta)*sigi/dble(nev)
      sumsigst=sumsigst+4d0*pi**2*sin(phih-pheta)*sigt/dble(nev)
      sumsibc=sumsibc+4d0*pi**2*sin(phih+pheta)*sib/dble(nev)
      sumsigci=sumsigci+4d0*pi**2*sin(phih+pheta)*sigi/dble(nev)
      sumsigct=sumsigct+4d0*pi**2*sin(phih+pheta)*sigt/dble(nev)
      enddo
c           print*,1,ipoint,iz,ipt,z,sumsibs/sumsib,sumsigs/sumsig,sumsigs*sumsib/sumsig/sumsibs
      write(16,'(18g14.5)')q2,x,z,pt**2,sumsibs/sumsib,sumsigsi/sumsigi,sumsigst/sumsigt,sumsibc/sumsib,sumsigci/sumsigi,sumsigct/sumsigt
      write(*,'(i5,18g14.5)')ipoint,q2,x,z,pt**2,sumsibs/sumsib,sumsigsi/sumsigi,sumsigst/sumsigt,sumsibc/sumsib,sumsigci/sumsigi,sumsigct/sumsigt
        enddo
        close(16)
        enddo
        endif
        enddo
      end
      
      subroutine haprad3(ebeam,iy,nev,sib,sigi,sigt)  
      implicit none
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
          real*8 ebeam,born,borntest,sib,deltavr,vacpol,sig,sigi,sigt,sphi,fspen
          real*8 sh,slsh,xh,slxh,slm,llm,px2
          real*8 URAND8,sum,mu,rh,ta,phik,sigrf,sigrfsphi,sigrfsphiold,sigrfd0,sigex,sigexd0,sigexta,sex,am(3),bm(3),wrk(500),re,rere,rein,reex,ot,otr, phiar(4),tar(6)
          real*8 sum1in,sum2in,sum1ex,sum2ex,ta1,ta2,tam,rhophik,slx,cta1,cta2,xi,rhota
          integer nev,i,j,j1,mir,ma,id,iph,ico
          integer*8 iy
          external sigrfd0,sigrfsphi,sigrfsphiold,sigexd0,sigexta
      
      call haprad_init(ebeam)
      
      sib=borntest(x,z,q2,pt,ebeam,un,etal,le,etat,phih,pheta)
c      print'(4i3,3g12.5)',int(un),int(le),int(etal),int(etat),phih/pi*180,pheta/pi*180,sib  
c     print*,sib
      sib=born()
c       print*,'sib',sib
         if(abs(sib).lt.1d-23)return
c        stop
      px2=mpi2-q2-2d0*vm+mp2+(1d0-z)*sx 
      sh=s-q2-v1
      slsh=sqrt(sh**2-4d0*ml2*px2)
      xh=xx+q2-v2
      slxh=sqrt(xh**2-4d0*ml2*px2)
      slm=sqrt(q2*(q2+4d0*ml2))
      llm=log((q2+slm)**2/4d0/ml2/q2)/slm
      rhmax=px2-(mp+mpi)**2  

      deltavr=2d0*((q2+2d0*ml2)*llm-1d0)*log(rhmax/ml/sqrt(px2))+sh*log((sh+slsh)**2/4d0/ml2/px2)/slsh/2d0
     .        +xh*log((xh+slxh)**2/4d0/ml2/px2)/slxh/2d0+sphi(sh,xh,q2+2d0*ml2,ml2,ml2,px2)-2d0+(1.5d0*q2+4d0*ml2)*llm
     .        -(q2+2d0*ml2)/slm*(slm**2*llm**2/2d0+2d0*fspen(2d0/(q2/slm+1d0))-pi**2/2d0)
      ta1=-q2/s
      ta2=q2/xx
        rein=0d0
        reex=0d0
        if(nev.ge.1)then
      tam=(ta1+ta2)/2d0
       cta1=(s+sls)*(q2*s+2d0*ml2*sx+sls*(sls*tam+sqrt((q2+s*tam)**2+4d0*ml2*(q2+tam*(sx-mp2*tam)))))
     .       /2d0/ml2/(2d0*mp2*q2+sqrt(lq)*sls+s*sx)
       slx=sqrt(xx**2-4d0*ml2*mp2) 
       cta2=(slx+xx)*(sx*xx+sqrt(lq)*slx-2*mp2*q2)/2d0/mp2/(2d0*ml2*sx-q2*xx+slx*(slx*tam
     .  +sqrt(4d0*ml2*(q2+tam*(sx-mp2*tam))+(q2-tam*xx)**2)))
      sum1in=0d0
      sum2in=0d0
      sum1ex=0d0
      sum2ex=0d0
      do i=1,nev
       xi=ml2*(2d0*mp2*q2+sqrt(lq)*sls+s*sx)/(s+sls)*cta1**(urand8(iy))
       ta=(xi-ml2*l1/xi-s*q2-2d0*sx*ml2)/ls 
       rhota=sls/log(cta1)/sqrt((q2+s*ta)**2+4d0*ml2*(q2+ta*(sx-mp2*ta))) 
      phik=2d0*atan(sqrt(lq*((q2+s*ta)**2+4*ml2*(q2+ta*(sx-mp2*ta))))*tan(pi*urand8(iy))/
     -  (s*sx*ta+q2*(sp+2d0*mp2*ta)+2d0*sqrt(l1*(q2+ta*(sx-mp2*ta)))))
      if(phik.lt.0d0)phik=2d0*pi+phik
      rhophik=sqrt(lq*((q2+s*ta)**2+4d0*ml2*(q2+ta*(sx-mp2*ta))))/
     -  (2d0*pi*(s*sx*ta + q2*(sp + 2*mp2*ta)-2d0*sqrt(l1*(q2+ta*(sx-mp2*ta)))*cos(phik)))
      rh=urand8(iy)*rhmax
      sum1in=sum1in+1d0/dble(nev)*sigrf(rh,ta,phik)/rhophik/rhota*rhmax
      sum1ex=sum1ex+1d0/dble(nev)*sigex(ta,phik)/rhophik/rhota
       xi=(slx+xx)*(sx*xx+sqrt(lq)*slx-2*mp2*q2)/4d0/mp2*cta2**(urand8(iy)-1d0)!urand8(iy)
       ta=(xi-ml2*l1/xi+xx*q2-2d0*sx*ml2)/slx**2 
       rhota=slx/log(cta2)/sqrt((q2-xx*ta)**2+4d0*ml2*(q2+ta*(sx-mp2*ta))) 
      phik=2d0*atan(sqrt(lq*(4d0*ml2*(q2+ta*(sx-mp2*ta))+(q2-ta*xx)**2))*tan(pi*urand8(iy))/
     -  (2d0*q2*(sp-mp2*ta)+2d0*sqrt(l1*(q2+ta*(sx-mp2*ta)))+sx*ta*xx))
      if(phik.lt.0d0)phik=2d0*pi+phik
      rhophik=sqrt(lq*(4d0*ml2*(q2+ta*(sx-mp2*ta))+(q2-ta*xx)**2))/
     -  (2d0*pi*(q2*(sp-2*mp2*ta)+sx*ta*xx-2d0*sqrt(l1*(q2+ta*(sx-mp2*ta)))*cos(phik)))
      rh=urand8(iy)*rhmax

      sum2in=sum2in+1d0/dble(nev)*sigrf(rh,ta,phik)/rhophik/rhota*rhmax
      sum2ex=sum2ex+1d0/dble(nev)*sigex(ta,phik)/rhophik/rhota
       enddo
            rein=sum1in+sum2in
            reex=sum1ex+sum2ex
        elseif(nev.eq.0)then
        am(1)=0d0    
        am(2)=tamin    
        am(3)=0d0    
        bm(1)=rhmax 
        bm(2)=tamax    
        bm(3)=2d0*pi
      id=1
      mir=100000
      ma=10*mir
            ot=1d-3 
      phiar(1)=0.d0
      phiar(2)=0.01d0*pi
      phiar(3)=2.d0*pi-0.01d0*pi
      phiar(4)=2.d0*pi
      tar(1)=tamin
      tar(2)=ta1-0.15d0*(ta1-tamin)
      tar(3)=ta1+0.15d0*(ta2-ta1)
      tar(4)=ta2-0.15d0*(ta2-ta1)
      tar(5)=ta2+0.15d0*(tamax-ta2)
      tar(6)=tamax

      do iph=1,3
      do ico=1,5
      am(2)=tar(ico)
      bm(2)=tar(ico+1)
      am(3)=phiar(iph)
      bm(3)=phiar(iph+1)
       if(am(2).gt.bm(2))write(*,*)' am(2)<bm(2)'
       if(am(1).gt.bm(1))write(*,*)' am(1)<bm(2)'
 
          call d01fce(3,am,bm,mir,ma,sigrfd0,ot,otr,500,wrk,rere,id)
        re=re+rere
        print*,iph,ico,re,rere
      enddo  
      enddo
        elseif(nev.eq.-1)then
        print*,'simps'
      call simpux(0d0,2d0*pi,150,1d-4,sigrfsphi,re)
        print*,'simps'
c        re=-re
        elseif(nev.eq.-2)then
c        print*,'simpsold'
c      call simpux(0d0,2d0*pi,150,1d-4,sigrfsphiold,re)
c        print*,'simpsold'
c        re=-re
c        print*,'exclu'
      id=1
      mir=6000
      ma=10*mir
            ot=1d-3 
      phiar(1)=0.d0
      phiar(2)=0.01d0*pi
      phiar(3)=2.d0*pi-0.01d0*pi
      phiar(4)=2.d0*pi
      tar(1)=tamin
      tar(2)=ta1-0.15d0*(ta1-tamin)
      tar(3)=ta1+0.15d0*(ta2-ta1)
      tar(4)=ta2-0.15d0*(ta2-ta1)
      tar(5)=ta2+0.15d0*(tamax-ta2)
      tar(6)=tamax

c      do iph=1,3
c      do ico=1,5
c      am(1)=tar(ico)
c      bm(1)=tar(ico+1)
c      am(2)=phiar(iph)
c      bm(2)=phiar(iph+1)
c       if(am(2).gt.bm(2))write(*,*)' am(2)<bm(2)'
c       if(am(1).gt.bm(1))write(*,*)' am(1)<bm(2)'
c      print*,'exlim',am(1),bm(1),am(2),bm(2)      
cx        pause    
c          call d01fce(2,am,bm,mir,ma,sigexd0,ot,otr,500,wrk,rere,id)
c        re=re+rere
c        print*,iph,ico,re,rere
c      enddo  
c      enddo
c 
        rere=0d0
c      print*,alfa**3/(256.d0*pi**4*sqly*amp)*ys*sx*barn/1.5435330503633793d0
      do ico=1,5
      call simpsx(tar(ico),tar(ico+1),150,1d-4,sigexta,re)
      rere=rere+re
      print*,'rrr',rere,re
      enddo
         sex=-alpha**3/(256.d0*pi**4*sqrt(lq)*mp)*sx**2/s*barn/pl*rere
      print*,sex,sib,sex/sib
       return     
       endif
         rein=2d0*pi*barn*alpha**3*s*sx**2/64d0/pi**2/mp/pl/ls/sqrt(lq)*rein
      sigi=(alpha/pi*(deltavr+vacpol(q2))+1d0)*sib+rein
         sigt=sigi-alpha**3/(256.d0*pi**4*sqrt(lq)*mp)*sx**2/s*barn/pl*reex
cx        print*,'rc',sig/sib-1d0
c       print*,j,sig-sib,sig/sib
        return
      end
      
      subroutine haprad_init(ebeam)
      implicit none
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
          real*8 ebeam
          integer ikey
      data ikey/0/
      pi=atan(1d0)*4d0
      data alpha/.729735d-2/,barn/0.389379d+6/,mp/.938272d0/,mn/.939565d0/,mpi/0.1395675d0/,ml/0.5110d-03/!,ml/0.105257d0/
      mp2=mp**2  
      mn2=mn**2  
      ml2=ml**2  
      mpi2=mpi**2  
        if(ikey.eq.0)call gridsidis_sf

        if(ikey.eq.0)call gridexc_sf
         ikey=1
        s=2d0*ebeam*mp
        ls=s**2-4d0*mp2*ml2
        sls=sqrt(ls)
        sx=q2/x
        xx=s-sx
        sp=s+xx
        lq=sx**2+4d0*mp2*q2
        l1=q2*(s*xx-mp2*q2)-ml2*lq   
        ph0=z*sx/2d0/mp
        pl=sqrt(ph0**2-pt**2-mpi2)
        v1=ph0*s/mp-pl*(s*sx+ 2d0*mp2*q2)/mp/sqrt(lq)-2d0*pt*sqrt(l1/lq)*cos(phih)
        v2=ph0*xx/mp-pl*(xx*sx-2d0*mp2*q2)/mp/sqrt(lq)-2d0*pt*sqrt(l1/lq)*cos(phih)
        vp=(v1+v2)/2d0
        vm=(v1-v2)/2d0
          vm=(ph0*sx-pl*sqrt(lq))/2d0/mp
        eh=-pt*sqrt(l1)*sin(phih)/2d0
        tamin=(sx-sqrt(lq))/2d0/mp2
        tamax=(sx+sqrt(lq))/2d0/mp2
        tamin=-q2/mp2/tamax
        tmb(1)=q2-2d0*ml2
        tmb(2)=(s*xx-mp2*q2)/2d0
        tmb(3)=(v1*v2-mpi2*q2)/2d0
        tmb(4)=(s*v2+xx*v1-z*q2*sx)/2d0
        tmb(5)=2d0*le*s*eh/sls
        tmb(6)=-sp*eh
        tmb(7)=le*s/4d0/sls*(lq*vp-sp*sx*(z*q2+vm)) 
        tmb(8)=-2d0*vp*eh
        tmb(9)=le/2d0/sls*(s*(q2*(z*sx*vp-mpi2*sp)+vm*(s*v2-xx*v1))+2d0*ml2*(4d0*mp2*vm**2+lq*mpi2-z*sx**2*(z*Q2+2d0*vm) ))
       call sfsidis(x,z,q2,pt,un,etal,etat,cos(pheta-phih),sin(pheta-phih),hsf0)
      end

      subroutine gridsidis_sf
      implicit none
      integer nx,nz,nq2,iq2,iz,ix,narg(3)
      parameter(nx=20)
      parameter(nz=20)
      parameter(nq2=20)
      common/sfgridsidis/x_net,z_net,q2_net,fuur,fllr,futf1tperpr,futh1r
     .,futh1tpr,fuucos2phir,fltr,fulsin2phir,fltcosphir,fulsinphir
     .,fllcosphir,fltcos2phir,fuucosphir,futsinphisr,futsin2phi1r
     .,futsin2phi2r,rarg
      real*8 x_net(nx),z_net(nz),q2_net(nq2),fuur(nq2,nz,nx),fllr(nq2,nz,nx),futf1tperpr(nq2,nz,nx),futh1r(nq2,nz,nx)
     .,futh1tpr(nq2,nz,nx),fuucos2phir(nq2,nz,nx),fltr(nq2,nz,nx),fulsin2phir(nq2,nz,nx),fltcosphir(nq2,nz,nx),fulsinphir(nq2,nz,nx)
     .,fllcosphir(nq2,nz,nx),fltcos2phir(nq2,nz,nx),fuucosphir(nq2,nz,nx),futsinphisr(nq2,nz,nx),futsin2phi1r(nq2,nz,nx)
     .,futsin2phi2r(nq2,nz,nx),rarg(nx+nz+nq2)
       print*,'---1'
        open(44,file='sidis.grid') 
       do ix=1,nx
        do iz=1,nz
          do iq2=1,nq2
        
       read(44,*)x_net(ix),z_net(iz),q2_net(iq2),fuur(iq2,iz,ix),fllr(iq2,iz,ix),futf1tperpr(iq2,iz,ix),futh1r(iq2,iz,ix)
     .,futh1tpr(iq2,iz,ix),fuucos2phir(iq2,iz,ix),fltr(iq2,iz,ix),fulsin2phir(iq2,iz,ix),fltcosphir(iq2,iz,ix),fulsinphir(iq2,iz,ix)
     .,fllcosphir(iq2,iz,ix),fltcos2phir(iq2,iz,ix),fuucosphir(iq2,iz,ix),futsinphisr(iq2,iz,ix),futsin2phi1r(iq2,iz,ix)
     .,futsin2phi2r(iq2,iz,ix)
        enddo
        enddo
        enddo
        close(44)
       do iq2=1,nq2
        rarg(iq2)=q2_net(iq2)
       enddo
       do iz=1,nz
        rarg(nq2+iz)=z_net(iz)
       enddo
       do ix=1,nx
        rarg(ix+nz+nq2)=x_net(ix)
       enddo
      end

      subroutine gridexc_sf
      implicit none
      integer nt,nw,nq2,it,iw,iq2,narg(3)
      parameter(nq2=18)
      parameter(nw=47)
      parameter(nt=61)
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2

      common/sfgridex/q2net,wnet,tnet,f11n,f22n,f33n,f44n,f55n,f66n,
     .f12rn,f12in,f13rn,f13in,f14rn,
     .f14in,f15rn,f15in,f16rn,f16in,
     .f23rn,f23in,f24rn,f24in,f25rn,
     .f25in,f26rn,f26in,f34rn,f34in,
     .f35rn,f35in,f36rn,f36in,f45rn,
     .f45in,f46rn,f46in,f56rn,f56in,rargex

      real*8 q2net(nq2),wnet(nw),tnet(nt),f11n(nt,nw,nq2),f22n(nt,nw,nq2),f33n(nt,nw,nq2),f44n(nt,nw,nq2),f55n(nt,nw,nq2),f66n(nt,nw,nq2),
     .f12rn(nt,nw,nq2),f12in(nt,nw,nq2),f13rn(nt,nw,nq2),f13in(nt,nw,nq2),f14rn(nt,nw,nq2),
     .f14in(nt,nw,nq2),f15rn(nt,nw,nq2),f15in(nt,nw,nq2),f16rn(nt,nw,nq2),f16in(nt,nw,nq2),
     .f23rn(nt,nw,nq2),f23in(nt,nw,nq2),f24rn(nt,nw,nq2),f24in(nt,nw,nq2),f25rn(nt,nw,nq2),
     .f25in(nt,nw,nq2),f26rn(nt,nw,nq2),f26in(nt,nw,nq2),f34rn(nt,nw,nq2),f34in(nt,nw,nq2),
     .f35rn(nt,nw,nq2),f35in(nt,nw,nq2),f36rn(nt,nw,nq2),f36in(nt,nw,nq2),f45rn(nt,nw,nq2),
     .f45in(nt,nw,nq2),f46rn(nt,nw,nq2),f46in(nt,nw,nq2),f56rn(nt,nw,nq2),f56in(nt,nw,nq2),rargex(nq2+nw+nt)
      real*8 a1r,a1i,a2r,a2i,a3r,a3i,a4r,a4i,a5r,a5i,a6r,a6i
      real*8 f1r,f1i,f2r,f2i,f3r,f3i,f4r,f4i,f5r,f5i,f6r,f6i
      real*8 w,t,q2       
        open(44,file='exclu.grid')
         read(44,*)   
       do iq2=1,nq2
        do iw=1,nw
          do it=1,nt
         read(44,*)q2net(iq2),wnet(iw),tnet(it),a1r,a1i,a2r,a2i,a3r,a3i,a4r,a4i,a5r,a5i,a6r,a6i
            w=wnet(iw)*1d-3
           q2=q2net(iq2)
            t=mpi2-q2-2*(((mpi2+w**2-mn2)*(W**2-q2-mp2))/(4.*W**2) - 
     -     Sqrt(-mpi2+(W**2-mn2+mpi2)**2/(4.*W**2))*Cos(pi*tnet(it)/180d0)*
     -      Sqrt((W**2-mp2-q2)**2/(4.*W**2) + q2))
        f1r= a1r+(2d0*a6r*q2-(a3r-a4r)*(q2+t-mpi2))/2d0/(w-mp)+a4r*(w-mp)
        f1i= a1i+(2d0*a6i*q2-(a3i-a4i)*(q2+t-mpi2))/2d0/(w-mp)+a4i*(w-mp)
        f2r=-a1r+(2d0*a6r*q2-(a3r-a4r)*(q2+t-mpi2))/2d0/(w+mp)+a4r*(w+mp)
        f2i=-a1i+(2d0*a6i*q2-(a3i-a4i)*(q2+t-mpi2))/2d0/(w+mp)+a4i*(w+mp)
        f3r=a3r-a4r+a2r*(w-mp)+(a2r-2d0*a5r)*q2/2d0/(w+mp)    
        f3i=a3i-a4i+a2i*(w-mp)+(a2i-2d0*a5i)*q2/2d0/(w+mp)
        f4r=a3r-a4r-a2r*(w+mp)-(a2r-2d0*a5r)*q2/2d0/(w-mp)    
        f4i=a3i-a4i-a2i*(w+mp)-(a2i-2d0*a5i)*q2/2d0/(w-mp)
        f5r=(4*(a1r+(a4r-a6r)*(w-mp))*(q2+(w+mp)**2) - 
     -    2*a2r*(4*q2*w**2+(mp2+q2-w**2)**2) + 
     -    2*(w**2+ mpi2-mn2)*(2*(a3r-a4r)*(w+mp)-2*a5r*q2+ 
     -    a2r*(q2+2*w**2-2*mp2))+(q2+t-mpi2)*
     -     (4*(a3r-a4r)*w-2*a5r*(mp2+q2-w**2) + 
     -       a2r*(3*(mp2+q2)+5*w**2)))/(8.*w)    
        f5i=(4*(a1i+(a4i-a6i)*(w-mp))*(q2+(w+mp)**2) - 
     -    2*a2i*(4*q2*w**2+(mp2+q2-w**2)**2) + 
     -    2*(w**2+ mpi2-mn2)*(2*(a3i-a4i)*(w+mp)-2*a5i*q2+ 
     -    a2i*(q2+2*w**2-2*mp2))+(q2+t-mpi2)*
     -     (4*(a3i-a4i)*w-2*a5i*(mp2+q2-w**2)+ 
     -       a2i*(3*(mp2+q2)+5*w**2)))/(8.*w)
        f6r=(-4*a4r*(t-mn2+mp2)*w+4*a4r*mp*(mp2+q2-w**2) - 
     -    4*(q2+(mp-w)**2)*(a1r+a6r*(mp+w))+2*(w**2-mn2+mpi2)*
     -     (2*a4r*mp+2*a5r*q2+2*a3r*(w-mp)+a2r*(2*mp2-q2-2*w**2))+ 
     -    2*a2r*((mp2+q2)**2+w**2*(w**2-2*mp2+2*q2))-(mpi2-q2-t)*
     -     (4*a3r*w+2*a5r*(mp2+q2-w**2)-a2r*(3*mp2+3*q2+5*w**2)))/8./w    
        f6i=(-4*a4i*(t-mn2+mp2)*w+4*a4i*mp*(mp2+q2-w**2) - 
     -    4*(q2+(mp-w)**2)*(a1i+a6i*(mp+w))+2*(w**2-mn2+mpi2)*
     -     (2*a4i*mp+2*a5i*q2+2*a3i*(w-mp)+a2i*(2*mp2-q2-2*w**2))+ 
     -    2*a2i*((mp2+q2)**2+w**2*(w**2-2*mp2+2*q2))-(mpi2-q2-t)*
     -     (4*a3i*w+2*a5i*(mp2+q2-w**2)-a2i*(3*mp2+3*q2+5*w**2)))/8./w
c            if(q2.gt.1d-2.and.tnet(it).gt.1d0)then
c           print*,it,iw,iq2,q2,w,t,tnet(it)
c           print*,'const',mp,mn,mpi
c           print'(13g12.5)',a1r,a1i,a2r,a2i,a3r,a3i,a4r,a4i,a5r,a5i,a6r,a6i
c           print'(13g12.5)',f1r,f1i,f2r,f2i,f3r,f3i,f4r,f4i,f5r,f5i,f6r,f6i
c            stop
c            endif
        f11n(it,iw,iq2)=f1r*f1r+f1i*f1i
        f12rn(it,iw,iq2)=f1r*f2r+f1i*f2i
        f12in(it,iw,iq2)=f1i*f2r-f1r*f2i
        f13rn(it,iw,iq2)=f1r*f3r+f1i*f3i
        f13in(it,iw,iq2)=f1i*f3r-f1r*f3i
        f14rn(it,iw,iq2)=f1r*f4r+f1i*f4i
        f14in(it,iw,iq2)=f1i*f4r-f1r*f4i
        f15rn(it,iw,iq2)=f1r*f5r+f1i*f5i
        f15in(it,iw,iq2)=f1i*f5r-f1r*f5i
        f16rn(it,iw,iq2)=f1r*f6r+f1i*f6i
        f16in(it,iw,iq2)=f1i*f6r-f1r*f6i
        f22n(it,iw,iq2)=f2r*f2r+f2i*f2i    
        f23rn(it,iw,iq2)=f2r*f3r+f2i*f3i
        f23in(it,iw,iq2)=f2i*f3r-f2r*f3i
        f24rn(it,iw,iq2)=f2r*f4r+f2i*f4i
        f24in(it,iw,iq2)=f2i*f4r-f2r*f4i
        f25rn(it,iw,iq2)=f2r*f5r+f2i*f5i
        f25in(it,iw,iq2)=f2i*f5r-f2r*f5i
        f26rn(it,iw,iq2)=f2r*f6r+f2i*f6i
        f26in(it,iw,iq2)=f2i*f6r-f2r*f6i
        f33n(it,iw,iq2)=f3r*f3r+f3i*f3i    
        f34rn(it,iw,iq2)=f3r*f4r+f3i*f4i
        f34in(it,iw,iq2)=f3i*f4r-f3r*f4i
        f35rn(it,iw,iq2)=f3r*f5r+f3i*f5i
        f35in(it,iw,iq2)=f3i*f5r-f3r*f5i
        f36rn(it,iw,iq2)=f3r*f6r+f3i*f6i
        f36in(it,iw,iq2)=f3i*f6r-f3r*f6i
        f44n(it,iw,iq2)=f4r*f4r+f4i*f4i    
        f45rn(it,iw,iq2)=f4r*f5r+f4i*f5i
        f45in(it,iw,iq2)=f4i*f5r-f4r*f5i
        f46rn(it,iw,iq2)=f4r*f6r+f4i*f6i
        f46in(it,iw,iq2)=f4i*f6r-f4r*f6i
        f55n(it,iw,iq2)=f5r*f5r+f5i*f5i    
        f56rn(it,iw,iq2)=f5r*f6r+f5i*f6i
        f56in(it,iw,iq2)=f5i*f6r-f5r*f6i
        f66n(it,iw,iq2)=f6r*f6r+f6i*f6i    
        enddo
        enddo
        enddo
          do it=1,nt
           rargex(it)=tnet(it) 
        enddo
        do iw=1,nw
           rargex(nt+iw)=wnet(iw) 
        enddo
       do iq2=1,nq2
           rargex(nt+nw+iq2)=q2net(iq2) 
        enddo
      end
      

      real*8 function borntest(x,z,q2,pt,ebeam,un,etal,le,etat,phih,pheta)
      implicit none
      common/const/alpha,pi,barn,mp,mp2,mpi,mpi2,me,me2
      real*8 alpha,pi,barn,mp,mp2,mpi,mpi2,me,me2
      real*8 pt,ebeam,un,etal,le,etat,phih,pheta,s,y,ga,eps
      real*8 x,z,q2
      common/fstrf/fuu,fll,futf1tperp,futh1,futh1tp,fuucos2phi,flt,fulsin2phi,
     . fltcosphi, fulsinphi, fllcosphi, fltcos2phi, fuucosphi, futsinphis,futsin2phi
      real*8 fuu,fll,futf1tperp,futh1,futh1tp,fuucos2phi,flt,fulsin2phi,
     . fltcosphi, fulsinphi, fllcosphi, fltcos2phi, fuucosphi, futsinphis,futsin2phi
        s=2d0*ebeam*mp
        y=q2/x/s
        ga=2d0*mp*x/sqrt(q2)
       eps=(1d0-y-ga**2*y**2/4d0)/(1d0-y+y**2/2d0+ga**2*y**2/4d0)
       borntest=2d0*pi*barn*alpha**2*y**2/(2d0*x*y*q2*(1d0-eps))*(1d0+ga**2/2d0/x)*(un*(fuu+sqrt(2d0*eps*(1d0+eps))*cos(phih)*fuucosphi
     .+eps*cos(2d0*phih)*fuucos2phi)
     .+etal*sqrt(2d0*eps*(1d0+eps))*sin(phih)*fulsinphi+etal*eps*sin(2d0*phih)*fulsin2phi+etal*le*(sqrt(1-eps**2)*fll
     .+sqrt(2d0*eps*(1d0-eps))*cos(phih)*fllcosphi)+etat*(sin(phih-pheta)*futf1tperp+eps*sin(phih+pheta)*futh1
     .+eps*sin(3d0*phih-pheta)*futh1tp+sqrt(2d0*eps*(1d0+eps))*(sin(pheta)*futsinphis+sin(2d0*phih-pheta)*futsin2phi)
     .+le*(sqrt(1-eps**2)*cos(phih-pheta)*flt+sqrt(2d0*eps*(1d0-eps))*(cos(pheta)*fltcosphi+cos(2d0*phih-pheta)*fltcos2phi))))
cv        print*,'b1=',borntest 
       end

      real*8 function born()
      implicit none
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
      real*8 sum
      integer i
         sum=0d0 
         do i=1,9
         sum=sum+tmb(i)*hsf0(i)
         enddo
        born=2d0*pi*barn*alpha**2*s*sx**2/8d0/mp/pl/q2**2/ls*sum
c        print*,'b1=',born
       end

      double precision function vacpol(t)
      implicit none
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      double precision t,am2(3),suml,a2,sqlmi,allmi,aaa,bbb,ccc,sumh
      integer i
c
c    am2 : squared masses of charge leptons
c
      data am2/0.26110d-6,0.111637d-1,3.18301d0/

      suml=0.d0
      do i=1,3
	 a2=2.d0*am2(i)
	 sqlmi=sqrt(t*t+2.d0*a2*t)
	 allmi=log((sqlmi+t)/(sqlmi-t))/sqlmi
       suml=suml+2.d0*(t+a2)*allmi/3.d0-10.d0/9.d0+4.d0*a2*(1.d0-a2*allmi)/3.d0/t
       enddo 
      if(t.lt.1.d0)then
	aaa = -1.345d-9
	bbb = -2.302d-3
	ccc =  4.091d0
      elseif(t.lt.64.d0)then
	aaa = -1.512d-3
	bbb = -2.822d-3
	ccc =  1.218d0
      else
	aaa = -1.1344d-3
	bbb = -3.0680d-3
	ccc =  9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.d0+ccc*t))*2.d0*pi/alpha

      vacpol=suml+sumh

      return
      end

      real*8 function sigrfd0(ndim,arg)
      implicit none
      integer ndim
      real*8 arg(15),sigrf
      sigrfd0=sigrf(arg(1),arg(2),arg(3))
        print*,arg(1),arg(2),arg(3),sigrfd0
      return
      end

c      subroutine simpsx(a,b,np,ep,func,res)
c      subroutine simptx(a,b,np,ep,func,res)
c      subroutine simpux(a,b,np,ep,func,res)

      real*8 function sigrfsphi(phikk)
      implicit none
      common/phik/phik  
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      real*8 sigrfsta,phik,phikk
      external sigrfsta   
      phik=phikk  
      call simptx(tamin,tamax,150,1d-5,sigrfsta,sigrfsphi)
      return
      end

      real*8 function sigrfsta(taa)
      implicit none
      common/ta/ta
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
      real*8 sigrfsrh,ta,taa,phik
      external sigrfsrh   
      ta=taa  
      call simpsx(rhmax,1d-10,150,1d-6,sigrfsrh,sigrfsta)
      sigrfsta=-sigrfsta  
      return
      end

      real*8 function sigrfsphiold(phikk)
      implicit none
      common/phik/phik  
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      real*8 sigrfstaold,phik,phikk,ta1,ta2,tar(6),ep,re
      integer i
      external sigrfstaold   
      phik=phikk
            ep=1d-12
      ta1=-q2/s  
      ta2=q2/xx  
      tar(1)=tamin
      tar(2)=ta1-0.15d0*(ta1-tamin)
      tar(3)=ta1+0.15d0*(ta2-ta1)
      tar(4)=ta2-0.15d0*(ta2-ta1)
      tar(5)=ta2+0.15d0*(tamax-ta2)
      tar(6)=tamax
      sigrfsphiold=0d0
       do i=1,5
      call simptx(log(x+tar(i))+ep,log(x+tar(i+1))-ep,100,1d-3,sigrfstaold,re)
      sigrfsphiold=sigrfsphiold+re
c        print*,i,sigrfsphiold,re    
       enddo 
c        stop
c      call simptx(tamin,tamax,150,1d-3,sigrfstaold,sigrfsphiold)
      return
      end

      real*8 function sigrfstaold(taa)
      implicit none
      common/ta/ta
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
      real*8 sigrfsr,ta,taa,phik
      common/mu/mu  
      real*8 mu
      external sigrfsr   
      ta=exp(taa)-x  
        mu=ph0/mp+pl*(2d0*ta*mp2-sx)/mp/sqrt(lq)-2d0*mp*pt*cos(phih-phik)*sqrt((ta-tamin)*(tamax-ta)/lq)
c        mu=(ph0/mp+pl*(2d0*ta*mp2-sx)/mp/sqrt(lq)-pt*cos(phih-phik)*sqrt(1d0-(2d0*ta*mp2-sx)**2/lq))/mp
      call simpsx(rhmax/(1d0+ta-mu),1d-10,150,1d-4,sigrfsr,sigrfstaold)
      sigrfstaold=-sigrfstaold*(x+ta)  
      return
      end

      real*8 function sigrfsr(r)
      implicit none
      common/phik/phik  
      common/ta/ta
      real*8 sigrf,r,ta,phik  
      common/mu/mu  
      real*8 mu
      sigrfsr=(1d0+ta-mu)*sigrf(r*(1d0+ta-mu),ta,phik)
      return
      end

      real*8 function sigrfsrh(rh)
      implicit none
      common/phik/phik  
      common/ta/ta
      real*8 sigrf,rh,ta,phik  
      sigrfsrh=sigrf(rh,ta,phik)
      return
      end


       
      real*8 function sigrf(rh,ta,phik)
      implicit none
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
      real*8 rh,r,ta,phik 
      real*8 mu,z1,tmr(9,4),xr,zr,q2r,ptr,lqr,hsf(9),sum
      real*8 etaq,etak,etak1,etaph,etalr,etatr,cdph,sdph  
      integer i,j
        
        mu=ph0/mp+pl*(2d0*ta*mp2-sx)/mp/sqrt(lq)-2d0*mp*pt*cos(phih-phik)*sqrt((ta-tamin)*(tamax-ta)/lq)
c        mu=(ph0+pl*(2d0*ta*mp2-sx)/sqrt(lq)-pt*cos(phih-phik)*sqrt(1d0-(2d0*ta*mp2-sx)**2/lq))/mp
        z1=(q2*sp+ta*(s*sx+2d0*mp2*q2)-2d0*mp*sqrt((ta-tamin)*(tamax-ta)*l1)*cos(phik))/lq
        call tails(ta,phik,tmr,mu,z1)
          r=rh/(1d0+ta-mu) 
         xr=(q2+r*ta)/(sx-r)
         zr=z*sx/(sx-r)
         q2r=q2+r*ta
c        vm=(ph0*sx-pl*sqrt(lq))/2d0/mp
            lqr=(sx-r)**2+4d0*mp2*q2r   
c         ptr=sqrt(ph0**2-(z*sx*(sx-r)-2d0*mp2*(2d0*vm-mu*r))**2/4d0/mp2/((sx-r)**2+4d0*mp2*q2r)-mpi2)
c            print*,'1',ptr
c         ptr=sqrt(ph0**2-(2d0*mp2*mu*r+2d0*mp*pl*sqrt(lq)-z*sx*r)**2/lqr/4d0/mp2-mpi2)
c            print*,'2',ptr
         ptr=sqrt((z**2*sx**2-(2d0*mp2*mu*r+2d0*mp*pl*sqrt(lq)-z*sx*r)**2/lqr)/4d0/mp2-mpi2)
            if(isnan(ptr))then
                sigrf=0d0
                return
              endif  
c            print*,'2',ptr
c                stop
            if(isnan(ptr))print*,ph0**2-(z*sx*(sx-r)-2d0*mp2*(2d0*vm-mu*r))**2/4d0/mp2/((sx-r)**2+4d0*mp2*q2r)-mpi2
            if(etal**2+etat**2.gt.5d-1)then
********** scalar product vector decomposition
 
           etaq=-etal*sqrt(lq)/2d0/mp
           etak=-r*(etal*(sx-2d0*mp2*ta)+2d0*etat*mp2*sqrt((tamax-ta)*(ta-tamin))*cos(phik-pheta))/2d0/mp/sqrt(lq)
           etak1=-(etal*(s*sx+2d0*mp2*q2)+2d0*mp*etat*sqrt(l1)*cos(pheta))/2d0/mp/sqrt(lq)
           etaph=-etal*pl-etat*pt*cos(phih-pheta)
***********************************8
 
         etalr=(etal*(lq+r*(2d0*ta*mp2-sx))-2d0*etat*mp2*r*cos(pheta-phik)*sqrt((ta-tamin)*(tamax-ta)))/sqrt(lq*lqr)
c         etalr=2d0*mp*(etak-etaq)/sqrt(lqr)
          etatr=sqrt(etat**2+4d0*mp2*r*sqrt((tamax-ta)*(ta-tamin))*(mp2*r*sqrt((tamax-ta)*(ta-tamin))*(etal**2-etat**2*cos(phik-pheta)**2)
     -          +etal*etat*cos(phik-pheta)*(lq+r*(2d0*ta*mp2-sx)))/lq/lqr)  
c            print*,'4',etalr,etatr

        cdph=(-etaph+((etak-etaq)*(r*sx*z-2d0*sqrt(lq)*mp*pl-2*mp2*mu*r))/((r-sx)**2+4*mp2*(q2+r*ta)))/etatr/ptr
               if(isnan(cdph))print*,cdph,rh,ta,phik
               if(isnan(cdph))stop
 
c                print*,cdph
cx        cdph=(-etaph+(etalr*(r*sx*z-2d0*sqrt(lq)*mp*pl-2d0*mp2*mu*r))/2d0/mp/sqrt((r-sx)**2+4*mp2*(q2+r*ta)))/etatr/ptr
        sdph=-(-(((sx*(2*s*(Vm+q2*z)-sx*(V1+q2*z)))+r**2*(mu*s-V1+sx*z*z1)+ 
     -       r*(2*sx*V1 + s*(2*sx*ta*z-mu*sx-2*Vm)+sx*z*(q2-sx*z1))-2*mp2*(q2*(mu*r+2*Vp)+r*(2*ta*V1+mu*r*z1-2*Vm*z1)))*
     -         (-etat*(2*ml2*(sx*(sx-r)+2*mp2*(2*q2+r*ta))+2*mp2*q2*(q2+r*z1)-s*(r*s*ta+q2*(r+2*xx)-r*sx*z1))*Sin(pheta)
     -        +(2*mp*r*(etaq*(4*ml2*mp2-s**2)+etak1*(2*mp2*q2+s*sx))*Sqrt((tamax-ta)*(ta-tamin))*Sin(phik))/Sqrt(lq))) + 
     -      2*(etak1*((sx-r)**2+4*mp2*(q2+r*ta))+(etak-etaq)*(s*(sx-r)+2*mp2*(q2+r*z1)))*
     -       (pt*(2*ml2*(sx*(sx-r)+2*mp2*(2*q2+r*ta))+2*mp2*q2*(q2+r*z1)-s*(r*s*ta+q2*(r+2*xx)-r*sx*z1))*Sin(phih) + 
     -       (mp*r*Sqrt((tamax-ta)*(ta-tamin))*(s*(V1*xx+q2*sx*z-s*V2)-2*mp2*q2*V1+2*ml2*(sx**2*z-4*mp2*Vm))*Sin(phik))/Sqrt(lq)))/
     -      (4.*etatr*ptr*((sx-r)**2+4*mp2*(q2+r*ta))**1.5*sqrt(l1)*((q2*s*(r+xx)-mp2*(q2+r*z1)**2+r*s*(s*ta+(r-sx)*z1))
     -      /((r-sx)**2+4*mp2*(q2+r*ta))-ml2))
 
c               print*,'c1',cdph,sdph,cdph**2+sdph**2
c                stop
c        cdph=max(-1d0,min(1d0,cdph))
c        sdph=sign(sqrt(1d0-cdph**2),sdph)
        sdph=max(-1d0,min(1d0,sdph))
        cdph=sign(sqrt(1d0-sdph**2),cdph)
c                print*,'c2',cdph,sdph,cdph**2+sdph**2
        call sfsidis(xr,zr,q2r,ptr,un,etalr,etatr,cdph,sdph,hsf)
        else
        call sfsidis(xr,zr,q2r,ptr,un,etal,etat,cos(pheta-phih),sin(pheta-phih),hsf)
        endif
        sum=0d0
        
       do i=1,9
         sum=sum+(hsf(i)/q2r**2-hsf0(i)/q2**2)*tmr(i,1)/r
       do j=2,4
         sum=sum+hsf(i)/q2r**2*tmr(i,j)*r**(j-2)
       enddo
c        print*,i,sum
c        stop
       enddo
c        print'(23g14.5)',hsf0
c        print'(23g14.5)',hsf
c        print*,rh,ta,phik
c        print*
c        pause
        
c        stop
        sigrf=-sum/(1d0+ta-mu)
c        print*,sigrf,rh,ta,phik
        if(isnan(sigrf))pause
       end
       
      real*8 function sigexd0(ndim,arg)
      implicit none
      real*8 ta,phik,arg(15),sigex
      integer ndim
            ta=arg(1)
            phik=arg(2)
            sigexd0=sigex(ta,phik)
            return
c            print*,ta,phik,sigexd0
c            pause
            end

      real*8 function sigexta(ta1)
      implicit none
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      common/ta1/ta
      real*8 ta,ta1,sigexphik,phiar(4),re
      integer i
      external sigexphik
      ta=ta1
      phiar(1)=0.d0
      phiar(2)=0.01d0*pi
      phiar(3)=2.d0*pi-0.01d0*pi
      phiar(4)=2.d0*pi
        sigexta=0d0    
        do i=1,3
      call simpux(phiar(i),phiar(i+1),150,1d-4,sigexphik,re)
        sigexta=sigexta+re    
c       print*,i,sigexta,re     
        enddo
c       print*     
       end
     
      real*8 function sigexphik(phik)
      implicit none
      common/ta1/ta
      real*8 ta,phik,sigex  
            sigexphik=sigex(ta,phik)
            end


      real*8 function sigex(ta,phik)
      implicit none
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
      real*8 ta,phik,mu,z1,tmr(9,4),r,q2r,w2r,hsfex(9),cdph,sdph
      real*8 lqr,etaq,etak,etak1,etaph,ptr,etalr,etatr
      integer i,j
        mu=ph0/mp+pl*(2d0*ta*mp2-sx)/mp/sqrt(lq)-2d0*mp*pt*cos(phih-phik)*sqrt((ta-tamin)*(tamax-ta)/lq)
        z1=(q2*sp+ta*(s*sx+2d0*mp2*q2)-2d0*mp*sqrt((ta-tamin)*(tamax-ta)*l1)*cos(phik))/lq
        call tails(ta,phik,tmr,mu,z1)
            r=(mp2+mpi2-mn2-q2+(1d0-z)*sx-2d0*vm)/(1d0+ta-mu)
            q2r=q2+ta*r
            w2r=sx-q2+mp2-(ta+1d0)*r
            if(etal**2+etat**2.gt.5d-1)then

         lqr=(sx-r)**2+4d0*mp2*q2r   
         etalr=(etal*(lq+r*(2d0*ta*mp2-sx))-2d0*etat*mp2*r*cos(pheta-phik)*sqrt((ta-tamin)*(tamax-ta)))/sqrt(lq*lqr)
          etatr=sqrt(etat**2+4d0*mp2*r*sqrt((tamax-ta)*(ta-tamin))*(mp2*r*sqrt((tamax-ta)*(ta-tamin))*(etal**2-etat**2*cos(phik-pheta)**2)
     -          +etal*etat*cos(phik-pheta)*(lq+r*(2d0*ta*mp2-sx)))/lq/lqr)

            if(etatr.gt.1d-12)then
           etaq=-etal*sqrt(lq)/2d0/mp
           etak=-r*(etal*(sx-2d0*mp2*ta)+2d0*etat*mp2*sqrt((tamax-ta)*(ta-tamin))*cos(phik-pheta))/2d0/mp/sqrt(lq)
           etak1=-(etal*(s*sx+2d0*mp2*q2)+2d0*mp*etat*sqrt(l1)*cos(pheta))/2d0/mp/sqrt(lq)
           etaph=-etal*pl-etat*pt*cos(phih-pheta)
         ptr=sqrt((z**2*sx**2-(2d0*mp2*mu*r+2d0*mp*pl*sqrt(lq)-z*sx*r)**2/lqr)/4d0/mp2-mpi2)
        cdph=(-etaph+((etak-etaq)*(r*sx*z-2d0*sqrt(lq)*mp*pl-2*mp2*mu*r))/((r-sx)**2+4*mp2*(q2+r*ta)))/etatr/ptr
        sdph=-(-(((sx*(2*s*(Vm+q2*z)-sx*(V1+q2*z)))+r**2*(mu*s-V1+sx*z*z1)+ 
     -       r*(2*sx*V1 + s*(2*sx*ta*z-mu*sx-2*Vm)+sx*z*(q2-sx*z1))-2*mp2*(q2*(mu*r+2*Vp)+r*(2*ta*V1+mu*r*z1-2*Vm*z1)))*
     -         (-etat*(2*ml2*(sx*(sx-r)+2*mp2*(2*q2+r*ta))+2*mp2*q2*(q2+r*z1)-s*(r*s*ta+q2*(r+2*xx)-r*sx*z1))*Sin(pheta)
     -        +(2*mp*r*(etaq*(4*ml2*mp2-s**2)+etak1*(2*mp2*q2+s*sx))*Sqrt((tamax-ta)*(ta-tamin))*Sin(phik))/Sqrt(lq))) + 
     -      2*(etak1*((sx-r)**2+4*mp2*(q2+r*ta))+(etak-etaq)*(s*(sx-r)+2*mp2*(q2+r*z1)))*
     -       (pt*(2*ml2*(sx*(sx-r)+2*mp2*(2*q2+r*ta))+2*mp2*q2*(q2+r*z1)-s*(r*s*ta+q2*(r+2*xx)-r*sx*z1))*Sin(phih) + 
     -       (mp*r*Sqrt((tamax-ta)*(ta-tamin))*(s*(V1*xx+q2*sx*z-s*V2)-2*mp2*q2*V1+2*ml2*(sx**2*z-4*mp2*Vm))*Sin(phik))/Sqrt(lq)))/
     -      (4.*etatr*ptr*((sx-r)**2+4*mp2*(q2+r*ta))**1.5*sqrt(l1)*((q2*s*(r+xx)-mp2*(q2+r*z1)**2+r*s*(s*ta+(r-sx)*z1))
     -      /((r-sx)**2+4*mp2*(q2+r*ta))-ml2))
 
c               print*,'c1',cdph,sdph,cdph**2+sdph**2
        sdph=max(-1d0,min(1d0,sdph))
        cdph=sign(sqrt(1d0-sdph**2),cdph)


        call sfexclus(q2r,w2r,mpi2-q2-2d0*vm+(mu-ta)*r,un,etalr,etatr,cdph,sdph,hsfex)
            else
        call sfexclus(q2r,w2r,mpi2-q2-2d0*vm+(mu-ta)*r,un,etal,etat,cos(pheta-phih),sin(pheta-phih),hsfex)
            endif
            else
        call sfexclus(q2r,w2r,mpi2-q2-2d0*vm+(mu-ta)*r,un,etal,etat,cos(pheta-phih),sin(pheta-phih),hsfex)
            endif
            sigex=0d0
            do i=1,9
            do j=1,4
            sigex=sigex+hsfex(i)*tmr(i,j)/q2r**2*r**(j-2)
            enddo     
            enddo     

            sigex=sigex/(1d0+ta-mu)
            end

        subroutine sfexclus(q2,w2,t,un,etal,etat,cdph,sdph,hsfex)
            implicit none
      integer nt,nw,nq2,it,iw,iq2,nargex(3),i
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      parameter(nq2=18)
      parameter(nw=47)
      parameter(nt=61)
      data nargex/nt,nw,nq2/
      common/sfgridex/q2net,wnet,tnet,f11n,f22n,f33n,f44n,f55n,f66n,
     .f12rn,f12in,f13rn,f13in,f14rn,
     .f14in,f15rn,f15in,f16rn,f16in,
     .f23rn,f23in,f24rn,f24in,f25rn,
     .f25in,f26rn,f26in,f34rn,f34in,
     .f35rn,f35in,f36rn,f36in,f45rn,
     .f45in,f46rn,f46in,f56rn,f56in,rargex

      real*8 q2net(nq2),wnet(nw),tnet(nt),f11n(nt,nw,nq2),f22n(nt,nw,nq2),f33n(nt,nw,nq2),f44n(nt,nw,nq2),f55n(nt,nw,nq2),f66n(nt,nw,nq2),
     .f12rn(nt,nw,nq2),f12in(nt,nw,nq2),f13rn(nt,nw,nq2),f13in(nt,nw,nq2),f14rn(nt,nw,nq2),
     .f14in(nt,nw,nq2),f15rn(nt,nw,nq2),f15in(nt,nw,nq2),f16rn(nt,nw,nq2),f16in(nt,nw,nq2),
     .f23rn(nt,nw,nq2),f23in(nt,nw,nq2),f24rn(nt,nw,nq2),f24in(nt,nw,nq2),f25rn(nt,nw,nq2),
     .f25in(nt,nw,nq2),f26rn(nt,nw,nq2),f26in(nt,nw,nq2),f34rn(nt,nw,nq2),f34in(nt,nw,nq2),
     .f35rn(nt,nw,nq2),f35in(nt,nw,nq2),f36rn(nt,nw,nq2),f36in(nt,nw,nq2),f45rn(nt,nw,nq2),
     .f45in(nt,nw,nq2),f46rn(nt,nw,nq2),f46in(nt,nw,nq2),f56rn(nt,nw,nq2),f56in(nt,nw,nq2),rargex(nq2+nw+nt),a2,a30,a31,a3,w2arg,wcor,q2cor,q2arg

      real*8 f11,f22,f33,f44,f55,f66,
     .f12r,f12i,f13r,f13i,f14r,
     .f14i,f15r,f15i,f16r,f16i,
     .f23r,f23i,f24r,f24i,f25r,
     .f25i,f26r,f26i,f34r,f34i,
     .f35r,f35i,f36r,f36i,f45r,
     .f45i,f46r,f46i,f56r,f56i
      real*8 q2,w2,t,sx,zz,lq,thcm,pt,un,etal,etat,cdph,sdph,hsfex(9)
      real*8 arg(3),dfint
      real*8 w00p,w11p,w22,w22r,w01r,w01i,w02r,w02i,w12r,w12i
      real*8 r0,r1,r2,r3,r4,r5,w
      data a2/1.15d0/,a30/-1.23d0/,a31/0.16d0/
      thcm=acos(((w2+mpi2-mn2)*(w2-q2-mp2)+2d0*w2*(t+q2-mpi2))/sqrt(((w2+mpi2-mn2)**2-4d0*mpi2*w2)*((w2-q2-mp2)**2+4d0*q2*w2)))
            if(q2.gt.5.0) then
c      print*,'Warning: Q2>5 GeV^2 in exclusive model!'
c      print*,'Using extrapolation from MAID2003'
      q2cor=(5.d0**a2)/(q2**a2)
      q2arg=5.d0
      else
      q2cor=1.d0
      q2arg=q2
      endif

            if(sqrt(w2).lt.1.07784) return
      if(sqrt(w2).gt.2.0) then
c      print*,'Warning: W>2 GeV in exclusive model!'
c      print*,'Using extrapolation from MAID2003 (A.Browman PRL35,Cornell)'
      a3=a30+a31*thcm
      if(thcm.lt.50.0) a3=a30+a31*50.d0
      if(thcm.gt.100.0) a3=a30+a31*100.d0
      wcor=(2.d0**a3)/(w2**(a3/2d0))
      w2arg=4.d0
      else
      w2arg=w2
      wcor=1.d0
      endif

        arg(1)=180d0/pi*thcm
        arg(2)=sqrt(w2arg)*1000d0
        arg(3)=q2arg
            sx=w2+q2-mp2
            lq=sx**2+4d0*mp2*q2
            zz=1d0+(t+mp2-mn2)/sx
            pt=sqrt(((mp2+sx+ t-mn2)*(sx*(mpi2-q2-t)+q2*(mp2+sx+t-mn2))-mpi2*sx**2-mp2*((mpi2+q2+ t)**2-4*mpi2*t))/lq)  
        f11=wcor*q2cor*dfint(3,arg,nargex,rargex,f11n)   
        f22=wcor*q2cor*dfint(3,arg,nargex,rargex,f22n)   
        f33=wcor*q2cor*dfint(3,arg,nargex,rargex,f33n)   
        f44=wcor*q2cor*dfint(3,arg,nargex,rargex,f44n)   
        f55=wcor*q2cor*dfint(3,arg,nargex,rargex,f55n)   
        f66=wcor*q2cor*dfint(3,arg,nargex,rargex,f66n)   
        f12r=wcor*q2cor*dfint(3,arg,nargex,rargex,f12rn)   
        f12i=wcor*q2cor*dfint(3,arg,nargex,rargex,f12in)
        f13r=wcor*q2cor*dfint(3,arg,nargex,rargex,f13rn)   
        f13i=wcor*q2cor*dfint(3,arg,nargex,rargex,f13in)   
        f14r=wcor*q2cor*dfint(3,arg,nargex,rargex,f14rn)   
        f14i=wcor*q2cor*dfint(3,arg,nargex,rargex,f14in)   
        f15r=wcor*q2cor*dfint(3,arg,nargex,rargex,f15rn)   
        f15i=wcor*q2cor*dfint(3,arg,nargex,rargex,f15in)   
        f16r=wcor*q2cor*dfint(3,arg,nargex,rargex,f16rn)   
        f16i=wcor*q2cor*dfint(3,arg,nargex,rargex,f16in)   
        f23r=wcor*q2cor*dfint(3,arg,nargex,rargex,f23rn)
        f23i=wcor*q2cor*dfint(3,arg,nargex,rargex,f23in)   
        f24r=wcor*q2cor*dfint(3,arg,nargex,rargex,f24rn)   
        f24i=wcor*q2cor*dfint(3,arg,nargex,rargex,f24in)   
        f25r=wcor*q2cor*dfint(3,arg,nargex,rargex,f25rn)   
        f25i=wcor*q2cor*dfint(3,arg,nargex,rargex,f25in)   
        f26r=wcor*q2cor*dfint(3,arg,nargex,rargex,f26rn)   
        f26i=wcor*q2cor*dfint(3,arg,nargex,rargex,f26in)   
        f34r=wcor*q2cor*dfint(3,arg,nargex,rargex,f34rn)   
        f34i=wcor*q2cor*dfint(3,arg,nargex,rargex,f34in)
        f35r=wcor*q2cor*dfint(3,arg,nargex,rargex,f35rn)   
        f35i=wcor*q2cor*dfint(3,arg,nargex,rargex,f35in)   
        f36r=wcor*q2cor*dfint(3,arg,nargex,rargex,f36rn)   
        f36i=wcor*q2cor*dfint(3,arg,nargex,rargex,f36in)   
        f45r=wcor*q2cor*dfint(3,arg,nargex,rargex,f45rn)   
        f45i=wcor*q2cor*dfint(3,arg,nargex,rargex,f45in)   
        f46r=wcor*q2cor*dfint(3,arg,nargex,rargex,f46rn)   
        f46i=wcor*q2cor*dfint(3,arg,nargex,rargex,f46in)   
        f56r=wcor*q2cor*dfint(3,arg,nargex,rargex,f56rn)   
        f56i=wcor*q2cor*dfint(3,arg,nargex,rargex,f56in)   
      
        w=sqrt(w2)
        r0=mp2*(q2+t-mpi2)+q2*(q2 + t-2d0*mn2+mpi2)+(mpi2+q2-t)*w2 
        r1=(w+mp)**2+q2    
        r2=(w-mp)**2+q2
        r3=(w+mn)**2-mpi2
        r4=(w-mn)**2-mpi2
        r5=(mpi2-mn2)*(mp2+q2)+(mn2+mp2+mpi2-q2-2d0*t)*w2-w2**2
        w00p=un*((f11*r1*r3*(mp-w)**2+f22*r2*r4*(mp+w)**2)/2.+2*q2*(r3/r1*f55+r4/r2*f66)*w2+ 
     -       (r5*(lq*(w2-mp2)*f12r-4*q2*w2*f56r))/r1/r2)/w2+
     -  2*etat*pt*sdph*(lq*(w2-mp2)*f12i+4*q2*w2*f56i)/sqrt(lq)/w 
        w22=(2*etat*f12i*sqrt(lq)*pt*sdph*(w2-mp2))/w+un*(f11*r1*r3*(w-mp)**2+f22*r2*r4*(mp+W)**2+2*f12r*r5*(w2-mp2))/2d0/w2
        w11p= (un*(r1*(w-mp)**2*(r4*f44+4*w*f14r)+r2*(w+mp)**2*(r3*f33+4*w*f23r)+2*r5*(mp2-w2)*f34r)
     -       +2*etat*sdph*(r5*(f14i*r1*(mp-W)**2-f23i*r2*(mp+w)**2)+lq*(r4*f24i-r3*f13i-4*w*f12i+2*pt**2*w*f34i)*(w2-mp2))/sqrt(lq)/pt)/2./w2    
        w01r=pt*sqrt(q2)*un*(r1*(mp-W)*(r4*f46r+2*w*f16r)-r2*(mp+w)*(r3*f35r+2*w*f25r)+r5*(f45r*(w-mp)+ (mp + W)*f36r))/sqrt(lq)/w + 
     -  (etat*sqrt(q2)*sdph*(r5*(f16i*r1*(mp-w)+f25i*r2*(mp+w))-lq*(f15i*r3*(mp-w)+f26i*r4*(mp+w)+2*pt**2*W*(f45i*(mp-W)+f36i*(mp+w)))))/r1/r2/w
        w01i=pt*Sqrt(Q2)*un*(r1*(w-mp)*(r4*f46i+2*w*f16i)+r2*(mp+w)*(f35i*r3+2*f25i*w)+r5*(f45i*(mp-w)-f36i*(mp+w)))/sqrt(lq)/w 
     -      +etat*sqrt(q2)*sdph*(r5*(f16r*r1*(mp-w)+f25r*r2*(mp+w))-lq*(f15r*r3*(mp-w)+f26r*r4*(mp+w)+2*pt**2*w*(f45r*(mp-w)+f36r*(mp+w))))/r1/r2/W
        w02r=2*etal*pt*Sqrt(Q2)*(f25i*r2*(mp+w)-f16i*r1*(mp-w))/Sqrt(lq) + 
     -  cdph*etat*Sqrt(Q2)*((f16i*r1*r5-f15i*lq*r3)*(w-mp)+r2*(f26i*r1*r4-f25i*r5)*(mp+w))/r1/r2/w    
        w02i=2*etal*pt*Sqrt(Q2)*(f25r*r2*(mp+w)-f16r*r1*(mp-w))/Sqrt(lq) + 
     -  (cdph*etat*Sqrt(Q2)*((f16r*r1*r5-f15r*lq*r3)*(w-mp)+r2*(f26r*r1*r4-f25r*r5)*(mp+w)))/r1/r2/w
        w12r=-etal*pt**2*(f14i*r1*(mp-w)**2+f23i*r2*(mp+w)**2)/w + 
     -  cdph*etat*pt*(r2*(mp+w)*(f23i*r5*(mp+w)+r1*(mp-w)*(f24i*r4-4*f12i*w))-f14i*r1*r5*(mp-w)**2-f13i*lq*r3*(mp2-w2))/2d0/sqrt(lq)/w2  
        w12i=etal*(-(r1*r2*(f11*r1*r3*(mp-w)**2+f22*r2*r4*(mp+w)**2))+ 
     -       2*(-(lq*pt**2*W*(f14r*r1*(mp-w)**2+f23r*r2*(mp+w)**2))+ 
     -          f12r*lq*r5*(mp2-w2)))/2./r1/r2/w2 + 
     -  (cdph*etat*pt*(-f14r*r1*r5*(mp-W)**2+r2*(f23r*r5*(mp+W)**2+f24r*r1*r4*(mp2-w2))+f13r*lq*r3*(w2-mp2)))/
     -   (2.*Sqrt(lq)*W**2)
        hsfex(1)=w22/8d0/pi/alpha
        hsfex(2)=(4d0*(q2*w00p+sqrt(q2/lq)*r0/pt*w01r)+r0**2/lq*w11p)/lq/8d0/pi/alpha
        hsfex(3)=w11p/8d0/pi/alpha 
        hsfex(4)=-(r0/lq*w11p+2d0*sqrt(q2/lq)*w01r/pt)/8d0/pi/alpha
        hsfex(5)=Sqrt(q2/lq)*w01i/pt/4d0/pi/alpha 
        hsfex(6)=-(2d0*sqrt(q2)*w02r+r0*w12r/pt/sqrt(lq))/pt/lq/4d0/pi/alpha
        hsfex(7)=(2d0*sqrt(q2)*w02i+r0*w12i/pt/sqrt(lq))/pt/lq/4d0/pi/alpha
        hsfex(8)=w12r/sqrt(lq)/pt**2/4d0/pi/alpha
        hsfex(9)=-w12i/sqrt(lq)/pt**2/4d0/pi/alpha

            
 
            end

       
       subroutine tails(ta,phik,tmr,mu,z1)
       implicit none
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
            real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      common/var/x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb,hsf0,tamin,tamax,rhmax
          real*8 x,z,q2,pt,un,le,etal,etat,phih,pheta,s,ls,sls,sx,xx,sp,lq,l1,ph0,pl,v1,v2,vp,vm,eh,tmb(9),hsf0(9),tamin,tamax,rhmax
        real*8 ta,phik,ekr,mu,z1,z2,bid,bip,bi21,bis,bir,fir,tmr(9,4)
        ekr=-mp*sqrt((ta-tamin)*(tamax-ta)*l1/lq)*sin(phik)/2d0
        z2=z1-ta
       bid=1d0/z1/z2
       bip=1d0/z2+1d0/z1
       bi21=1d0/z1**2
       bis=1d0/z2**2+1d0/z1**2
       bir=1d0/z2**2-1d0/z1**2
       fir=ml2*bis-(q2+2d0*ml2)*bid
c      print*,'fir',fir,ml2*bis,q2+2d0*ml2,bid  
        tmr(1,1)=4d0*fir*(q2-2d0*ml2)
        tmr(1,2)=4d0*ta*fir
        tmr(1,3)=-4d0-2d0*ta**2*bid
        tmr(1,4)=0d0
        tmr(2,1)=2d0*fir*(s*xx-mp2*q2)  
        tmr(2,2)=(sp*sx*bip+2d0*ml2*sp*bir+2d0*(sx-2d0*mp2*ta)*fir-bid*sp**2*ta)/2d0
        tmr(2,3)=((4d0*ml2+ta*(2d0*ta*mp2-sx))*bid-sp*bip+4d0*mp2)/2d0
        tmr(2,4)=0d0
        tmr(3,1)=4d0*fir*(v1*v2-mpi2*q2)/2d0
        tmr(3,2)=2d0*((mu*vm-ta*mpi2)*fir+vp*(mu*ml2*bir+bip*vm-ta*vp*bid))
        tmr(3,3)=(2d0*mu**2*ml2 + ta*(mpi2*ta - mu*vm))*bid-mu*vp*bip+2d0*mpi2
        tmr(3,4)=0d0
        tmr(4,1)=4d0*fir*(s*v2+xx*v1-z*q2*sx)/2d0
        tmr(4,2)=(s*v1-v2*xx)*bip+ml2*(mu*sp+2d0*vp)*bir-2d0*ta*sp*vp*bid+((mu-2d0*ta*z)*sx+2d0*vm)*fir
        tmr(4,3)=((8d0*mu*ml2+ta*((2d0*ta*z-mu)*sx-2d0*vm))*bid-(mu*sp+2d0*vp)*bip + 4d0*sx*z)/2d0
        tmr(4,4)=0d0
c        print*,'tm1',tmr(1,1),tmr(1,2),tmr(1,3)
c        print*,'tm2',tmr(2,1),tmr(2,2),tmr(2,3)
c        print*,'tm3',tmr(3,1),tmr(3,2),tmr(3,3)
c        print*,'tm4',tmr(4,1),tmr(4,2),tmr(4,3)
c        stop
        tmr(5,1)=8d0*fir*le*s*eh/sls
        tmr(5,2)=le*s/l1/sls*(eh*(2d0*(sx*(q2+4d0*ml2)+2d0*ta*(s*xx-2d0*mp2*(2d0*ml2+q2)))*fir
     .             +q2*(sp*(sx*bip+2d0*ml2*bir)-ta*(4*s*xx+sx**2)*bid)) 
     .             +2d0*ekr*(ml2*(sx*(s*v2-v1*xx-q2*sp*z)+4d0*mp2*q2*vp)*bir
     .             +((q2+4d0*ml2)*(4d0*mp2*vm-sx**2*z)+sp*(s*v2-v1*xx))*fir))
cccccccccccccccc (B2) contribution 
     .             +2d0*le*ml2/l1/sls*(eh*(2d0*(2d0*ml2*lq+(q2+s*ta)*(2d0*mp2*q2+s*sx))*bi21-sx*lq*bip
     .             +(2d0*q2*xx*sx+ta*sx*(2d0*s**2-sp**2)+4d0*mp2*q2*(ta*s-q2)-4d0*ml2*lq)*bid) 
     .             +2d0*ekr*(sx*(z*sp*q2+v1*xx-s*v2)-4d0*mp2*q2*vp)*(bid*xx-bi21*s))  
c        tmr(5,3)=le*s/l1/sls*(eh*(8d0*ml2*(ta*(mp2*ta-sx)-q2)*bi21+(q2*(sp+4d0*mp2*ta)+2d0*ta*s*sx)*bip 
c     .           +ta*(4d0*ml2*(2d0*ta*mp2-sx)+q2*(sx-4d0*s)-2d0*ta*s**2)*bid)
c     .           +2d0*ekr*(2d0*ml2*(sx*(2d0*z*q2+2d0*vm+(ta*z-mu)*sx)-4d0*mp2*(mu*q2+ta*vm))*bi21
c     .           +ta*(2d0*ml2*(z*sx**2-4d0*mp2*vm)-2d0*mp2*q2*v1+s*(z*sx*q2-s*v2+v1*xx))*bid))
ccccccccccccccccccccccccc hat theta_53
c     .           +2d0*bi21*le*s/l1/sls*(ekr*(2d0*(mu*q2+ta*v1)*(s*xx-mp2*q2)+(q2+s*ta)*(z*q2*sx-s*v2-v1*xx))-eh*(q2+s*ta)**2)
        tmr(5,3)=(le*s*(eh*(2*lq - bip*q2*sp + 
     -         bid*ta*((4*ml2 + q2)*(-sx + 2*mp2*ta) - 2*ta*(-(mp2*q2) + s*xx))) + 
     -      ekr*(bip*(-(lq*vp) + sp*sx*(vm + q2*z)) + 
     -         bid*ta*(sp*(-(s*v2) + v1*xx) - (4*ml2 + q2)*(4*mp2*vm - sx**2*z)))))/
     -  (l1*sls) 
cccccccccccccccc (B2) contribution 
     .           +2d0*le*ml2/l1/sls*(eh*(2d0*((q2+2d0*ml2)*(2d0*ta*mp2+xx)-(ta*xx+2d0*ml2)*s)*bi21-lq*bip 
     .           +(4d0*ml2*(sx-2d0*ta*mp2)+2d0*q2*s+ta*(s**2+xx**2))*bid) 
     .           +2d0*ekr*((2d0*ml2*(z*sx**2-4d0*mp2*vm)+2d0*mp2*q2*v2+xx*(v1*xx-s*v2-q2*sx*z))*bi21
     .           +(2d0*ml2*(4d0*mp2*vm-z*sx**2)+2d0*mp2*q2*v1+s*(s*v2-v1*xx-q2*sx*z))*bid))
        tmr(5,4)=0d0
        tmr(6,1)=-4d0*fir*sp*eh
        tmr(6,2)=(eh*((4d0*mp2*q2*(q2+4d0*ml2)-sx**2*(q2-4d0*ml2)-8d0*q2*s*xx)*(sx*bip+2d0*ml2*bir-ta*sp*bid) 
     .            +2d0*sp*(2d0*ta*(2d0*mp2*(q2+2d0*ml2)-s*xx)-sx*(q2+4d0*ml2))*fir) 
     .            +2d0*ekr*sp*(ml2*(sx*(z*sp*q2-s*v2+v1*xx)-4d0*mp2*q2*vp)*bir 
     .            +((q2+4d0*ml2)*(z*sx**2-4d0*mp2*vm)+sp*(v1*xx-s*v2))*fir))/2d0/l1
        tmr(6,3)=(2d0*eh*((2d0*q2*(s*xx-2d0*mp2*q2)-ta*sx*(sx**2+3d0*s*xx-4d0*ml2*mp2)-(q2+2d0*ml2)*sx**2)*bip 
     .            +ml2*(2d0*ta*(2d0*mp2*(q2+2d0*ml2)-s*xx)-sx*(q2+4d0*ml2))*bir-q2*sp*fir
     .            +sp*(ta**2*(sx**2+2d0*s*xx-2d0*mp2*(q2+4d0*ml2))+2d0*ta*sx*(q2+2*ml2)-4d0*ml2*q2)*bid+sp*sx**2)  
     .            +ekr*(((q2+4d0*ml2)*(z*sx**2-4d0*mp2*vm)+sp*(v1*xx-s*v2))*(2*ml2*bir+sx*bip) 
     .            +2d0*ml2*(sx*(z*q2*sp-s*v2+v1*xx)-4d0*mp2*q2*vp)*bis 
     .            +(4d0*ta*(mp2*q2*(4d0*s*vm+sx*(v2-vm))+2d0*s*xx*(s*v2-v1*xx)+2d0*ml2*sp*(4d0*mp2*vm - z*sx**2))
     .            +sx*(3d0*ta*sx+2d0*(q2-2d0*ml2))*(s*v2-v1*xx-z*q2*sp)+8d0*mp2*q2*(q2-2d0*ml2)*vp)*bid))/2d0/l1
        tmr(6,4)=(eh*(((q2+4d0*ml2)*sx+2d0*ta*(s*xx-2d0*mp2*(q2+2d0*ml2)))*bip+sp*(ta*q2*bid-2*sx)) 
     .           +ekr*(((q2+4d0*ml2)*(4d0*mp2*vm-sx**2*z)+sp*(s*v2-v1*xx))*bip
     .           +bid*ta*(4*mp2*q2*vp+sx*(s*v2-v1*xx-z*q2*sp))))/2d0/l1
        tmr(7,1)=fir*le*s/sls*(lq*vp-sp*sx*(z*q2+vm)) 
        tmr(7,2)=le*s/2d0/sls*(q2*(4d0*mp2*vm-sx**2*z)*bip+ml2*(mu*lq-2d0*sx*(z*q2+vm))*bir 
     .           +(2d0*(4d0*ta*mp2-sx)*vp+(mu-2d0*ta*z)*sp*sx-2d0*s*v2+2d0*v1*xx)*fir
     .           +ta*(q2*(z*sp*sx-4d0*mp2*vp)+sx*(v1*xx-s*v2))*bid)
cccccccccccccccc (B2) contribution 
     .           +le*ml2/sls*((4*mp2*(ta*s*vm-q2*vp)-sx**2*(ta*z*s+z*q2+v1)+mu*lq*s)*bi21
     .           +bid*(4*mp2*(q2*vp-ta*vm*xx)+sx**2*(ta*z*xx+v2-z*q2)-mu*lq*xx))
        tmr(7,3)=le*s/4d0/sls*((sx*(4d0*z*q2+2d0*vm-mu*sx)-8d0*mp2*mu*q2)*bip+2d0*ml2*(4d0*mu*ta*mp2+2d0*vm-sx*(mu+2d0*ta*z))*bir 
     .           +2d0*(2d0*vp-mu*sp)*fir+ta*(4d0*(sx-2d0*ta*mp2)*vp+sp*(sx*(2d0*ta*z-mu)-2d0*vm))*bid)
cccccccccccccccc (B2) contribution 
     .           +le*ml2/sls*(bi21*(2d0*mp2*(mu*(q2+s*ta)-2d0*ta*vp)+sx*(s*(ta*z-2d0*mu)+2d0*vp-z*q2)+(mu-ta*z)*sx**2)
     .           +bid*(2d0*mp2*(mu*(q2-ta*xx)+2d0*ta*vp)+sx*((2d0*mu-ta*z)*s-2d0*vp-q2*z)-mu*sx**2))
        tmr(7,4)=le*s/4d0/sls*(bip*((mu+2d0*ta*z)*sx-2d0*vm-4d0*mu*ta*mp2)+ta*(mu*sp-2d0*vp)*bid)
cccccccccccccccc (B2) contribution 
     .           +le*ml2/sls*(bi21*(2*mu*ta*mp2+mu*xx-ta*z*sx-v2)+bid*(2*mu*ta*mp2+v1-mu*s-sx*ta*z))
        tmr(8,1)=-8d0*fir*vp*eh
        tmr(8,2)=(eh*((q2*sx*(sx*vp-2d0*s*v2)-2d0*vm*(2d0*l1+q2*s*sx))*bip-2d0*ml2*(2d0*mu*l1+q2*sp*vp)*bir 
     .            +vp*(2d0*ml2*(2d0*ta*(2d0*mp2*(q2+2d0*ml2)-s*xx)-(q2+4d0*ml2)*sx)*bis 
     .            +(4d0*ml2*((3d0*q2+4d0*ml2)*sx+ta*(2d0*s*xx-4d0*mp2*(3d0*q2+2d0*ml2)-sx**2))
     .            +q2*(ta*(12d0*s*xx+sx**2)+2d0*q2*(sx-6d0*ta*mp2)))*bid))
     .            +2d0*ekr*vp*(((q2+4d0*ml2)*(z*sx**2-4d0*mp2*vm)+sp*(v1*xx-s*v2))*fir
     .            +ml2*(sx*(v1*xx-s*v2+z*sp*q2)-4d0*mp2*q2*vp)*bir))/l1
        tmr(8,3)=(eh*((2d0*mu*(q2-2d0*ml2)*q2*sp+ta*(2d0*(q2+8d0*ml2)*sx*vp+q2*(mu*sp*sx+2d0*s*v1-2d0*v2*xx)) 
     .           -2d0*ta**2*(4d0*(q2+4d0*ml2)*vp*mp2-sp*(s*v1+v2*xx)))*bid 
     .           +2d0*ml2*mu*(2d0*ta*(2d0*(q2+2d0*ml2)*mp2-s*xx)-(q2+4d0*ml2)*sx)*bir+4d0*sp*sx*vm-2d0 *mu*ml2*q2*sp*bis 
     .           +bip*(2d0*ta*(v2*xx**2-s**2*v1)+8*ml2*(2*mp2*ta-sx)*vm+q2*(4d0*mu*(s*xx-2d0*mp2*q2)-sx*(mu*sx+2d0*vm))))  
     .           +2d0*ekr*(((q2+4*ml2)*(z*sx**2-4*mp2*vm)+sp*(v1*xx-s*v2))*(vm*bip+mu*ml2*bir) 
     .           +mu*ml2*(sx*(z*sp*q2+v1*xx-s*v2)-4*mp2*q2*vp)*bis+bid*(mu*(q2-2d0*ml2)*(4*mp2*q2*vp+sx*(s*v2-v1*xx-z*sp*q2))
     .           +ta*(sp*sx*v2*(v1+vp)+2d0*vm*vp*((sx-4d0*s)*xx+2d0*mp2*(3d0*q2+8d0*ml2))
     .           +z*sx*(q2*(v2*xx-s*v1)-(8d0*ml2+q2)*sx*vp)))))/2d0/l1
        tmr(8,4)=mu*(eh*(bip*((q2+4d0*ml2)*sx+2d0*ta*(s*xx-2*(q2+2*ml2)*mp2))+sp*(ta*q2*bid-2*sx)) 
     .           +ekr*(((q2+4*ml2)*(4*mp2*vm-sx**2*z)+sp*(s*v2-v1*xx))*bip+ta*(4d0*q2*vp*mp2+sx*(s*v2-v1*xx-q2*sp*z))*bid))/2d0/l1
        tmr(9,1)=2d0*fir*le/sls*(s*(q2*(z*sx*vp-mpi2*sp)+vm*(s*v2-xx*v1))+2d0*ml2*(4d0*mp2*vm**2+lq*mpi2-z*sx**2*(z*Q2+2d0*vm)))
        tmr(9,2)=le*s/sls*(q2*sx*(z*vm-mpi2)*bip+ml2*(q2*(mu*z*sx-2d0*mpi2)+(mu*sx-2d0*vm)*vm)*bir 
     .           +ta*(q2*(mpi2*sp-z*sx*vp)+vm*(v1*xx-s*v2))*bid+(2d0*vm*(2d0*mu*s-vp)+2d0*ta*(z*sx*vp-mpi2*sp)-mu*(v1+vm)*sx)*fir)
cccccccccccccccc (B2) contribution 
     .           +2d0*le*ml2/sls*(2d0*ml2*(2d0*mpi2*(2d0*ta*mp2-sx)+2d0*(sx*z-2d0*mp2*mu)*vm+z*(mu-ta*z)*sx**2)*bis 
     .           +sx*(z*q2*(mu*s-vp)+vm*((mu+ta*z)*s-v1)-mpi2*(q2+s*ta))*bi21 
     .           +(sx*((mpi2-vm*z)*(ta*xx+3d0*q2+8d0*ml2)+(vm+q2*z)*(v2-mu*xx))
     .           +2d0*(4d0*mp2*(mu*vm-mpi2*ta)+(ta*z-mu)*z*sx**2)*(q2+2d0*ml2))*bid)  
        tmr(9,3)=le*s/2d0/sls*((2d0*(2d0*mpi2*q2+vm**2)-mu*sx*(2d0*z*q2+vm))*bip+ml2*(mu*((2d0*ta*z-mu)*sx+2d0*vm)-4d0*mpi2*ta)*bir
     .            +mu*(2d0*vp-mu*sp)*fir+ta*(2d0*ta*(mpi2*sp-sx*vp*z)+mu*sx*(v1+vm)+2d0*vm*(vp-2d0*mu*s))*bid)
cccccccccccccccc (B2) contribution 
     .            +le*ml2/sls*(4d0*ml2*(mpi2+mp2*mu**2-mu*z*sx)*bis
     .            +(2d0*mpi2*(ta*xx-q2)+sx*(mu*(ta*z-mu)*s-2d0*ta*z*vp+mu*z*q2+mu*v1)+2d0*vm*(v2-mu*xx))*bi21
     .            +bid*(mu*((q2+2*ml2)*(5d0*z*sx-4d0*mp2*mu)+(mu-ta*z)*sx*xx)-2d0*mpi2*(ta*s+3d0*q2+4d0*ml2)
     .            +sx*(2d0*z*ta*vp-mu*v2-2d0*ml2*mu*z)+2d0*vm*(mu*s-v1)))
        tmr(9,4)=le*s*((2*ta*(2*mpi2-mu*sx*z)+mu*(mu*sx-2*vm))*bip+mu*ta*(mu*sp-2*vp)*bid)/4d0/sls
cccccccccccccccc (B2) contribution 
     .            +le*ml2/sls*(bi21*(mu*(ta*z*sx+mu*xx-v2)-2*mpi2*ta)+bid*(mu*(ta*z*sx+v1-mu*s)-2*ta*mpi2))        
        return
        end
       
      subroutine sfsidis(x,z,q2,pt,un,etal,etat,cdph,sdph,hsf)
      implicit none
      common/const/alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      real*8 alpha,pi,barn,mp,mp2,mn,mn2,mpi,mpi2,ml,ml2
      integer nx,nz,nq2,narg(3)
      parameter(nx=20)
      parameter(nz=20)
      parameter(nq2=20)
      data narg/nx,nz,nq2/
      common/sfgridsidis/x_net,z_net,q2_net,fuur,fllr,futf1tperpr,futh1r
     .,futh1tpr,fuucos2phir,fltr,fulsin2phir,fltcosphir,fulsinphir
     .,fllcosphir,fltcos2phir,fuucosphir,futsinphisr,futsin2phi1r
     .,futsin2phi2r,rarg
      real*8 x_net(nx),z_net(nz),q2_net(nq2),fuur(nq2,nz,nx),fllr(nq2,nz,nx),futf1tperpr(nq2,nz,nx),futh1r(nq2,nz,nx)
     .,futh1tpr(nq2,nz,nx),fuucos2phir(nq2,nz,nx),fltr(nq2,nz,nx),fulsin2phir(nq2,nz,nx),fltcosphir(nq2,nz,nx),fulsinphir(nq2,nz,nx)
     .,fllcosphir(nq2,nz,nx),fltcos2phir(nq2,nz,nx),fuucosphir(nq2,nz,nx),futsinphisr(nq2,nz,nx),futsin2phi1r(nq2,nz,nx)
     .,futsin2phi2r(nq2,nz,nx),rarg(nx+nz+nq2)
      real*8 arg(3),dfint
      common/fstrf/fuu,fll,futf1tperp,futh1,futh1tp,fuucos2phi,flt,fulsin2phi,
     . fltcosphi, fulsinphi, fllcosphi, fltcos2phi, fuucosphi, futsinphis,futsin2phi
      real*8 fuu, fll, futf1tperp, futh1, futh1tp, fuucos2phi, flt, fulsin2phi,
     . fltcosphi, fulsinphi, fllcosphi, fltcos2phi, fuucosphi, futsinphis,futsin2phi
      real*8 hsf(9)
      real*8 q2,pt,un,etal,etat,cdph,sdph,x,z
      real*8 avp,avk,avkg
      real*8 avpt,avllpt,avutptf1tperp,avutpth1,avutpth1tp,avuucos2phipt
      real*8 ph0,pl,c1,sx,lq,l2,l3,vm,eta1,eta2
        data avp,avk,avkg/2d-1,2.5d-1,1.9d-1/
        arg(1)=q2
        arg(2)=z
        arg(3)=x
c        arg(1)=q2_net(10)
c        arg(2)=z_net(1)
c        arg(3)=x_net(1)
c        print*,'s1',fuur(10,1,1)
c        print*,'s2',dfint(3,arg,narg,rarg,fuur)
c            stop
cccccccccccccccc pochemy v sechenii fltcos2phi -> futsin2phi       
       avpt=avp+avk*z**2
       fuu=exp(-pt**2/avpt)/pi/avpt*x*dfint(3,arg,narg,rarg,fuur)

       avllpt=avp+avkg*z**2
       fll=exp(-pt**2/avllpt)/pi/avllpt*x*dfint(3,arg,narg,rarg,fllr)
       
       
       avutptf1tperp=avp+avk*0.19d0/(avk+0.19d0)*z**2
       futf1tperp=2d0*pt*exp(-pt**2/avutptf1tperp)/(pi*avutptf1tperp)**1.5*x*dfint(3,arg,narg,rarg,futf1tperpr)
       avutpth1=avp*1.5d0/(avp+1.5d0)+avk*z**2
       futh1=2d0*pt*exp(-pt**2/avutpth1)/(pi*avutpth1)**1.5*x*dfint(3,arg,narg,rarg,futh1r)

       avutpth1tp=avp*1.5d0/(avp+1.5d0)+avk*z**2*1.8d-1/(avk+1.8d-1)
       futh1tp=4d0*pt**3*exp(-pt**2/avutpth1tp)/3d0/pi**1.5/avutpth1tp**2.5*x*dfint(3,arg,narg,rarg,futh1tpr)

       avuucos2phipt=avp*1.5d0/(avp+1.5d0)+avk*z**2*1.9d-1/(avk+1.9d-1)
       fuucos2phi=pt**2*exp(-pt**2/avuucos2phipt)/pi/avuucos2phipt**2*x*dfint(3,arg,narg,rarg,fuucos2phir)
       
       flt=2d0*pt*exp(-pt**2/avllpt)/(pi*avllpt)**1.5*x*dfint(3,arg,narg,rarg,fltr)
       
       fulsin2phi=pt**2*exp(-pt**2/avutpth1)/pi/avutpth1**2*x*dfint(3,arg,narg,rarg,fulsin2phir)

       fltcosphi=exp(-pt**2/avllpt)/pi/avllpt*x*dfint(3,arg,narg,rarg,fltcosphir)

       fulsinphi=2d0*pt*exp(-pt**2/avutpth1)/(pi*avutpth1)**1.5*x*dfint(3,arg,narg,rarg,fulsinphir)

       fllcosphi=2d0*pt*exp(-pt**2/avllpt)/(pi*avllpt)**1.5*x*dfint(3,arg,narg,rarg,fllcosphir)

       fltcos2phi=pt**2*exp(-pt**2/avllpt)/pi/avllpt**2*x*dfint(3,arg,narg,rarg,fltcos2phir)

       fuucosphi=2d0*pt*exp(-pt**2/avpt)/(pi*avpt)**1.5*x*dfint(3,arg,narg,rarg,fuucosphir)
       
       futsinphis=(1d0-pt**2/avutpth1)*exp(-pt**2/avutpth1)*x*dfint(3,arg,narg,rarg,futsinphisr)  

       futsin2phi=pt**2*exp(-pt**2/avutptf1tperp)/pi/avutptf1tperp**2*x*dfint(3,arg,narg,rarg,futsin2phi1r)
     .+pt**2*avutptf1tperp*exp(-pt**2/avutpth1tp)/pi/avutpth1tp**3*x*dfint(3,arg,narg,rarg,futsin2phi2r)

cccccccccccccccccccccccc y prokudina ccccccccccccccccccc     
c       fltcos2phi -> futsin2phi
c         i vse sf *x
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       print*,arg,pt
c       print*,fuu,fll,futf1tperp,futh1,futh1tp
c       print*,fuucos2phi,flt,fulsin2phi,fltcosphi,fulsinphi
c       print*,fllcosphi,fltcos2phi,fuucosphi,futsinphis,futsin2phi
c       stop
        sx=q2/x
        ph0=z*sx/2d0/mp
        pl=sqrt(ph0**2-pt**2-mpi**2)
        c1=4d0*mp*pl*(q2+2d0*x*mp**2)/q2**2
        lq=sx**2+4d0*mp**2*q2
        vm=(sx*ph0-sqrt(lq)*pl)/2d0/mp
        l2=vm**2+mpi**2*q2
        l3=vm+z*q2
        eta1=etat*cdph
        eta2=etat*sdph
c         futsinphis=0d0 
        hsf(1)=c1*(un*(fuu-fuucos2phi)+eta2*(futh1tp-futf1tperp-futh1)) 
        hsf(2)=4d0*c1/lq**2/pt**2*((l3**2*sx**2-l2*lq)*(un*fuu-eta2*futf1tperp)
     .+(l3**2*sx**2+l2*lq)*(un*fuucos2phi+eta2*(futh1-futh1tp))+2d0*pt*sqrt(q2*lq)*l3*sx*(un*fuucosphi+eta2*(futsinphis-futsin2phi)))
        hsf(3)=2d0*c1/pt**2*(un*fuucos2phi+eta2*(futh1-futh1tp)) 
        hsf(4)=-2d0*c1/lq/pt*(2d0*l3*sx/pt*(un*fuucos2phi+eta2*(futh1-futh1tp))+sqrt(q2*lq)*(un*fuucosphi+eta2*(futsinphis-futsin2phi))) 
        hsf(5)=2d0*c1*eta2/pt*sqrt(q2/lq)*(fltcosphi-fltcos2phi) 
        hsf(6)=4d0*c1/lq/pt*((l3*sx/sqrt(lq)/pt*(eta1*(futh1+futh1tp)+etal*fulsin2phi))
     .+sqrt(q2)*(etal*fulsinphi+eta1*(futsinphis+futsin2phi)))
        hsf(7)=-4d0*c1/lq/pt*(l3*sx/sqrt(lq)/pt*(etal*fll+eta1*flt)+sqrt(q2)*(eta1*(fltcosphi+fltcos2phi)+etal*fllcosphi))
        hsf(8)=-2d0*c1/sqrt(lq)/pt**2*(eta1*(futh1+futh1tp)+etal*fulsin2phi) 
        hsf(9)=2d0*c1/sqrt(lq)/pt**2*(etal*fll+eta1*flt)
       end

       
       real*8 function sphi(s1,s2,s3,m12,m22,m32)
      implicit none
      real*8 s1,s2,s3,sl1,sl2,sl3,m12,m22,m32,fspen
      sl1=sqrt(s1**2-4d0*m12*m32)
      sl2=sqrt(s2**2-4d0*m22*m32)
      sl3=sqrt(s3**2-4d0*m12*m22)
      sphi=-s3*(-fspen(1d0+((sl2-s2)*(s2*(s3+sl3)-2d0*m22*s1))/4d0/m22/m32/sl3) 
     -         -fspen(1d0-((sl2+s2)*(s2*(s3+sl3)-2d0*m22*s1))/4d0/m22/m32/sl3) 
     -         +fspen(1d0+((sl1-s1)*(s2*(s3+sl3)-2d0*m22*s1))/2d0/m32/sl3/(s3+sl3))
     -         +fspen(1d0-((sl1+s1)*(s2*(s3+sl3)-2d0*m22*s1))/2d0/m32/sl3/(s3+sl3)) 
     -         +log((s1-sl1)/(s1+sl1))**2/4d0-log((s2-sl2)/(s2+sl2))**2/4d0)/sl3
      return
      end

      
****************** fspens *************************************

      double precision function fspens(x)
c
c    spence function
c
      implicit real*8(a-h,o-z)
      f=0.d0
      a=1.d0
      an=0.d0
      tch=1.d-16
  1   an=an+1.d0
      a=a*x
      b=a/an**2
      f=f+b
      if(b-tch)2,2,1
  2   fspens=f
      return
      end

      double precision function fspen(x)
      implicit real*8(a-h,o-z)
      data f1/1.644934d0/
      if(x)8,1,1
  1   if(x-.5d0)2,2,3
    2 fspen=fspens(x)
      return
    3 if(x-1d0)4,4,5
    4 fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
      return
    5 if(x-2d0)6,6,7
    6 fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
    7 fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
      return
    8 if(x+1d0)10,9,9
   9  fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
      return
  10  fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
      return
      end
      
       
      function dfint(narg,arg,nent,ent,table)
      implicit double precision (a-h,o-z)
      dimension arg(narg),nent(narg),ent(60),table(8000)
      dimension d(narg),ncomb(narg),ient(narg)
      kd=1
      m=1
      ja=1
	 do 5 i=1,narg
      ncomb(i)=1
      jb=ja-1+nent(i)
	 do 2 j=ja,jb
      if (arg(i).le.ent(j)) go to 3
    2 continue
      j=jb
    3 if (j.ne.ja) go to 4
      j=j+1
    4 jr=j-1
      d(i)=(ent(j)-arg(i))/(ent(j)-ent(jr))
      ient(i)=j-ja
      kd=kd+ient(i)*m
      m=m*nent(i)
    5 ja=jb+1
      dfint=0.d0
   10 fac=1.d0
      iadr=kd
      ifadr=1
	 do 15 i=1,narg
      if (ncomb(i).eq.0) go to 12
      fac=fac*(1.d0-d(i))
      go to 15
   12 fac=fac*d(i)
      iadr=iadr-ifadr
   15 ifadr=ifadr*nent(i)
      dfint=dfint+fac*table(iadr)
      il=narg
   40 if (ncomb(il).eq.0) go to 80
      ncomb(il)=0
      if (il.eq.narg) go to 10
      il=il+1
	 do 50	k=il,narg
   50 ncomb(k)=1
      go to 10
   80 il=il-1
      if(il.ne.0) go to 40
      return
      end
      FUNCTION	 URAND8(IY)
      IMPLICIT	 NONE
      REAL(8)	 URAND8
      INTEGER(8) IY
C     INTEGER(8),PARAMETER :: M=2**63
      INTEGER(8),PARAMETER :: A =3622009729038561421_8,
     *			      M2=2_8**62,
     *			      C =1949127854270302053_8
      REAL(8),PARAMETER    :: S =1.0842021724855044E-19
C     URAND8.FOR FORTRAN PROGRAM
C
C     Генератор Псевдослучайных Чисел
C     равномерно распределенных на
C     интервале (0.,1.)
C
C     INTERFACE TO REAL*8 FUNCTION URAND8(IY)
C     INTEGER(8) IY
C     END
C
C     Использован алгоритм из книги
C     Дж.Форсайт, M.Mалькольм, K.Mоулер
C     MAШИHHЫE METOДЫ MATEMATИЧECKИX BЫЧИCЛEHИЙ
C
C     Y    = A*Y  +  C (mod M),
C      i+1	i
C
C     здесь Y	  - положительное целое из интервала (0 , M)
C	    A,C,M - большие положительные целые константы,
C		    значения которых выбраны следующим образом:
C		    M = 2**63
C		    A = 8*int[(M/2)*(п/32)] + 5 =
C		      = 3622009729038561421 = (десятичное значение)
C		      = 03243f6a8885a308dh    (шестнадцатиричное значение)
C		    C = 2*int{(M/4)*[1-1/sqrt(3)]} + 1 =
C		      = 1949127854270302053 =
C		      = 01b0cb174df99c765h
C	    S	  - нормировочная константа 1/M
C
C     На выходе значение REAL*8 функции URAND8 равно  S*Y


C     PRINT *,S
      IY=IY*A+C
      IF(IY.LT.0)IY=(IY+M2)+M2
      URAND8=DFLOAT(IY)*S
      END

**************************simps
      subroutine simpu(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit double precision(a-h,o-z)
      dimension f(7),p(5)
      h=sign(h1,b1-a1)
      s=sign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=abs(reps1)
      aeps=abs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.d0
      
c      print*,'c ',x,f(1)
      
    4 x0=x
      if((x0+4.d0*h-b)*s) 5,5,6
    6 h=(b-x0)/4.d0
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=abs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.d0
      
c      print*,'d'
c      write(*,'(2f19.15)') x,f(k)
      
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.d0*f(3)+f(5))*2.d0*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=abs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=abs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.d0) 17,14,14
   17 h=2.d0*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.d0
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.d0
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

CDECK  ID>, SIMPS.
      subroutine simpsz(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
c simps
c a1,b1 -the limits of integration
c h1 -an initial step of integration
c reps1,aeps1 - relative and absolute precision of integration
c funct -a name of function subprogram for calculation of integrand +
c x - an argument of the integrand
c ai - the value of integral
c aih- the value of integral with the step of integration
c aiabs- the value of integral for module of the integrand
c this subrogram calculates the definite integral with the relative or
c absolute precision by simpson+s method with the automatical choice
c of the step of integration
c if aeps1    is very small(like 1.e-17),then calculation of integral
c with reps1,and if reps1 is very small (like 1.e-10),then calculation
c of integral with aeps1
c when aeps1=reps1=0. then calculation with the constant step h1
c
      implicit double precision(a-h,o-z)
      dimension f(7),p(5)
      h=sign(h1,b1-a1)
      s=sign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=abs(reps1)
      aeps=abs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.d0
    4 x0=x
      if((x0+4.d0*h-b)*s) 5,5,6
    6 h=(b-x0)/4.d0
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=abs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.d0
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.d0*f(3)+f(5))*2.d0*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=abs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=abs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.d0) 17,14,14
   17 h=2.d0*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.d0
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.d0
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpsx(a,b,np,ep,func,res)
      implicit double precision (a-h,o-z)
      external func
      step=(b-a)/dble(np)
      call simpsz(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      return
      end

      subroutine simptx(a,b,np,ep,func,res)
      implicit double precision (a-h,o-z)
      external func
      step=(b-a)/dble(np)
      call simpt(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      return
      end

      subroutine simpux(a,b,np,ep,func,res)
      implicit double precision (a-h,o-z)
      external func
      step=(b-a)/dble(np)
      call simpu(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      return
      end

      subroutine simpt(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit double precision(a-h,o-z)
      dimension f(7),p(5)
      h=sign(h1,b1-a1)
      s=sign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=abs(reps1)
      aeps=abs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.d0
      
c      print*,'a ',x,f(1)
      
    4 x0=x
      if((x0+4.d0*h-b)*s) 5,5,6
    6 h=(b-x0)/4.d0
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=abs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.d0
      
c      print*,'b ',x,f(k)
      
      
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.d0*f(3)+f(5))*2.d0*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=abs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=abs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.d0) 17,14,14
   17 h=2.d0*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.d0
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.d0
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end


*******************************d01fce
      subroutine d01fce(ndim, a, b, minpts, maxpts, functn, eps,
     * acc, lenwrk, wrkstr, finval, ifail)
      implicit double precision(a-h,o-z)
c     mark 8 release. nag copyright 1979.
c
c     adaptive multidimensional integration subroutine
c
c     *********  parameters for d01fce ****************************
c
c      input parameters
c
c     ndim    integer number of variables, must exceed 1 but
c	  not exceed 15.
c
c     a       real array of lower limits, with dimension ndim
c
c     b       real array of upper limits, with dimension ndim
c
c     minpts  integer minimum number of integrand values to be
c	  allowed, which must not exceed maxpts.
c
c     maxpts  integer maximum number of integrand values to be
c	  allowed, which must be at least
c	  2**ndim+2*ndim**2+2*ndim+1.
c
c     functn  externally declared user defined real function
c	  integrand. it must have parameters (ndim,z),
c	  where z is a real array of dimension ndim.
c
c     eps     real required relative accuracy, must be greater
c	  than zero
c
c     lenwrk  integer length of array wrkstr, must be at least
c	  2*ndim+4.
c
c     ifail   integer nag failure parameter
c	  ifail=0 for hard fail
c	  ifail=1 for soft fail
c
c      output parameters
c
c     minpts  integer number of integrand values used by the
c	  routine
c
c     wrkstr  real array of working storage of dimension (lenwrk).
c
c     acc     real estimated relative accuracy of finval
c
c     finval  real estimated value of integral
c
c     ifail   ifail=0 for normal exit, when estimated relative
c	    less integaccuracy rand values used.
c
c      ifail=1 if ndim.lt.2, ndim.gt.15, minpts.gt.maxpts,
c	    maxpts.lt.2**ndim+2*ndim*(ndim+1)+1, eps.le.0
c	    or lenwrk.lt.2*ndim+4.
c
c      ifail=2 if maxpts was too small for d01fce to obtain the
c	    required relative accuracy eps.  in this
c	    case d01fce returns a value of finval
c	    with estimated relative accuracy acc.
c
c      ifail=3 if lenwrk too small for maxpts integrand
c	    values.  in this case d01fce returns a
c	    value of finval with estimated accuracy
c	    acc using the working storage
c	    available, but acc will be greater
c	    than eps.
c
c     **************************************************************
c
c     .. scalar arguments ..
      double precision eps, finval, acc
      integer ifail, lenwrk, maxpts, minpts, ndim
c     .. array arguments ..
      double precision a, b, wrkstr
      dimension a(ndim), b(ndim), wrkstr(lenwrk)
c     .. function arguments ..
      double precision functn
c     ..
c     .. local scalars ..
      character*8 srname
      double precision 
     * abserr, df1, df2, difmax, f1, f2, f3, f4, half, lamda2,
     * lamda4, lamda5, one, ratio, rgncmp, rgnerr, rgnert, rgnval,
     * rgnvlt, rgnvol, rlndim, sum1, sum2, sum3, sum4, sum5, two,
     * twondm, weit1, weit2, weit3, weit4, weit5, weitp1, weitp2,
     * weitp3, weitp4, zero
      integer dvaxes, dvaxis, dvflag, funcls, ierror, j, k, maxaxs,
     * mxrgns, pointr, rgncls, rulcls, sbrgns, subrgn, subtmp,
     * tpontp, tpontr
c     .. local arrays ..
      dimension center(15), dif(15), oldcnt(15), width(15), z(15)
      integer dvcntl(15), dvcntr(15)
c     .. function references ..
      double precision x02aae
      integer p01aae, x02bbe
c     ..
      data srname /'  d01fce'/
      data zero, one, two, half /0.d0, 1.d0, 2.d0, 0.5d0/
c
c   subroutine initialisation and parameter checking
c
      if (ndim.lt.2 .or. ndim.gt.15) go to 560
      if (minpts.gt.maxpts) go to 560
      if (eps.le.zero) go to 560
      if (lenwrk.lt.2*ndim+4) go to 560
      funcls = 0
      finval = zero
      abserr = zero
      twondm = two**ndim
      rgnvol = twondm
      dvflag = 1
      fffff1 = dble(x02bbe(one))
      fffff2 = 1.d0/x02aae(0.0d0)
      maxaxs = int(min(fffff1,fffff2))
c     maxaxs = int(amin1(float(x02bbe(one)),1.0/x02aae(0.0d0)))
      maxaxs = (maxaxs-ndim)/(ndim+1)
      mxrgns = lenwrk/(2*ndim+4)
      sbrgns = 0
      rgnvlt = zero
      rgnert = zero
      do 20 j=1,ndim
       center(j) = (a(j)+b(j))*half
       dif(j) = zero
       width(j) = (b(j)-a(j))*half
       dvcntl(j) = 1
       dvcntr(j) = 1
       oldcnt(j) = center(j)
       rgnvol = rgnvol*width(j)
   20 continue
c
c   end subroutine initialisation
c   basic rule initialisation
c
      rulcls = 2**ndim + 2*ndim*ndim + 2*ndim + 1
      funcls = rulcls
      if (maxpts.lt.rulcls) go to 560
      rlndim = ndim
      lamda2 = sqrt(9.d0/70.d0)
      lamda4 = sqrt(9.d0/10.d0)
      lamda5 = sqrt(9.d0/19.d0)
      weit1 = (12824.d0-9120.d0*rlndim+400.d0*rlndim*rlndim)/19683.d0
      weit2 = 980.d0/6561.d0
      weit3 = (1820.d0-400.d0*rlndim)/19683.d0
      weit4 = 200.d0/19683.d0
      weit5 = 6859.d0/19683.d0/twondm
      weitp1 = (729.d0-950.d0*rlndim+50.d0*rlndim**2)/729.d0
      weitp2 = 245.d0/486.d0
      weitp3 = (265.d0-100.d0*rlndim)/1458.d0
      weitp4 = 25.d0/729.d0
      ratio = (lamda2/lamda4)**2
c
c   end basic rule initialisation
      go to 100
c   divide subregion with largest error and prepare to use
c   basic rule on each portion
c
   40 subrgn = 1
      pointr = wrkstr(1)
      rgncls = rulcls
      rgnvol = twondm
      tpontr = pointr + 2
      do 60 j=1,ndim
       tpontr = tpontr + 2
       center(j) = wrkstr(tpontr-1)
       width(j) = wrkstr(tpontr)
       dvcntr(j) = 1
       dvcntl(j) = 1
       oldcnt(j) = center(j)
       rgnvol = rgnvol*width(j)
   60 continue
      dvaxes = wrkstr(pointr+2)
      if (dvaxes.lt.0) go to 600
   80 dvaxis = dvaxes
      dvaxes = dvaxis/(ndim+1)
      dvaxis = dvaxis - (ndim+1)*dvaxes
      dvcntl(dvaxis) = 2*dvcntl(dvaxis)
      rgncls = rgncls*2
      if (dvaxes.gt.0) go to 80
      if (funcls+rgncls.gt.maxpts) go to 580
      if (rgncls/rulcls+sbrgns-1.gt.mxrgns) dvflag = 2
      funcls = funcls + rgncls
c      print *,funcls
      abserr = abserr - wrkstr(pointr)
      finval = finval - wrkstr(pointr+1)
c
c   begin basic rule
  100 do 120 j=1,ndim
       z(j) = center(j)
  120 continue
      sum1 = functn(ndim,z)
      sum2 = zero
      sum3 = zero
      do 140 j=1,ndim
       z(j) = center(j) - lamda2*width(j)
       f1 = functn(ndim,z)
       z(j) = center(j) + lamda2*width(j)
       f2 = functn(ndim,z)
       z(j) = center(j) - lamda4*width(j)
       f3 = functn(ndim,z)
       z(j) = center(j) + lamda4*width(j)
       f4 = functn(ndim,z)
       sum2 = sum2 + f1 + f2
       sum3 = sum3 + f3 + f4
       df1 = f1 + f2 - two*sum1
       df2 = f3 + f4 - two*sum1
       dif(j) = dif(j) + abs(df1-ratio*df2)
       z(j) = center(j)
  140 continue
      sum4 = zero
      do 200 j=2,ndim
       z(j-1) = center(j-1) - lamda4*width(j-1)
       do 160 k=j,ndim
	  z(k) = center(k) - lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k) + lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k)
  160  continue
       z(j-1) = center(j-1) + lamda4*width(j-1)
       do 180 k=j,ndim
	  z(k) = center(k) - lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k) + lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k)
  180  continue
       z(j-1) = center(j-1)
  200 continue
      sum5 = zero
      do 220 j=1,ndim
       z(j) = center(j) - lamda5*width(j)
  220 continue
  240 do 260 j=2,ndim
       if (z(j-1).lt.center(j-1)+width(j-1)) go to 280
       z(j-1) = center(j-1) - lamda5*width(j-1)
       z(j) = z(j) + two*lamda5*width(j)
  260 continue
      if (z(ndim).gt.center(ndim)+width(ndim)) go to 300
  280 sum5 = sum5 + functn(ndim,z)
      z(1) = z(1) + two*lamda5*width(1)
      go to 240
  300 rgnval = rgnvol*(weit1*sum1+weit2*sum2+weit3*sum3+weit4*
     * sum4+weit5*sum5)
      rgncmp = rgnvol*(weitp1*sum1+weitp2*sum2+weitp3*sum3+weitp4*
     * sum4)
      rgnerr = abs(rgnval-rgncmp)
c
c   end basic rule
c   store results of basic rule application
c
      rgnvlt = rgnvlt + rgnval
      rgnert = rgnert + rgnerr
      finval = finval + rgnval
      abserr = abserr + rgnerr
      if (dvflag.eq.0) go to 340
      if (dvflag.eq.2) go to 500
      pointr = mxrgns + sbrgns*(2*ndim+3) + 1
      sbrgns = sbrgns + 1
      wrkstr(sbrgns) = pointr
      subrgn = sbrgns
      tpontr = pointr + 2
      do 320 j=1,ndim
       tpontr = tpontr + 2
       wrkstr(tpontr-1) = center(j)
       wrkstr(tpontr) = width(j)
  320 continue
  340 wrkstr(pointr) = rgnert
      wrkstr(pointr+1) = rgnvlt
c   determine axis along which fourth difference is largest
      difmax = zero
      do 380 j=1,ndim
       if (difmax.gt.dif(j)) go to 360
       difmax = dif(j)
       dvaxis = j
  360	    dif(j) = zero
  380 continue
      tpontr = pointr + 2*(dvaxis+1)
      wrkstr(tpontr) = width(dvaxis)*half
      wrkstr(tpontr-1) = center(dvaxis) - wrkstr(tpontr)
      if (dvflag.ne.2) go to 400
      dvaxes = wrkstr(pointr+2)
      if (dvaxes.gt.maxaxs) dvaxes = -1
      dvaxis = dvaxis + (ndim+1)*dvaxes
  400 wrkstr(pointr+2) = dvaxis
      if (dvflag.eq.1) go to 460
c   determine the position in the parially ordered list of
c   the subregion which replaces most recently divided subregion
  420 subtmp = 2*subrgn
      if (subtmp.gt.sbrgns) go to 480
      tpontr = wrkstr(subtmp)
      if (subtmp.eq.sbrgns) go to 440
      tpontp = wrkstr(subtmp+1)
      if (wrkstr(tpontr).ge.wrkstr(tpontp)) go to 440
      subtmp = subtmp + 1
      tpontr = tpontp
  440 if (rgnert.ge.wrkstr(tpontr)) go to 480
      wrkstr(subtmp) = pointr
      wrkstr(subrgn) = tpontr
      subrgn = subtmp
      go to 420
c   when working storage is not used up, determine the
c   position in the partially ordered list for the description
c   of other portion(s) of most recently divided subregion
  460 subtmp = subrgn/2
      if (subtmp.lt.1) go to 480
      tpontr = wrkstr(subtmp)
      if (rgnert.le.wrkstr(tpontr)) go to 480
      wrkstr(subtmp) = pointr
      wrkstr(subrgn) = tpontr
      subrgn = subtmp
      go to 460
  480 rgnvlt = zero
      rgnert = zero
      if (dvflag.eq.2) go to 540
      dvflag = 1 - dvflag
c   count to determine the next part of the recently divided
c   subregion for application of the basic rule
  500 center(1) = center(1) + two*width(1)
      dvcntr(1) = dvcntr(1) + 1
      do 520 j=2,ndim
       if (dvcntr(j-1).le.dvcntl(j-1)) go to 100
       dvcntr(j-1) = 1
       center(j-1) = oldcnt(j-1)
       dvcntr(j) = dvcntr(j) + 1
       center(j) = center(j) + two*width(j)
  520 continue
      if (dvcntr(ndim).le.dvcntl(ndim)) go to 100
      center(ndim) = oldcnt(ndim)
      if (dvflag.eq.2) go to 340
c
c   end ordering of basic rule results
c   make checks for possible termination of routine
c
  540 acc = abserr/abs(finval)
      if (acc.gt.eps .or. funcls.lt.minpts) go to 40
c
c   loop back to apply basic rule
c
c   termination point, set ifail and return
c
      ierror = 0
      go to 620
  560 ierror = 1
      go to 620
  580 ierror = 2
      go to 620
  600 ierror = 3
  620 minpts = funcls
      ifail = p01aae(ifail,ierror,srname)
      return
      end

      double precision function x02aae(x)
      implicit double precision(a-h,o-z)
c     nag copyright 1975
c     mark 4.5 release
c+self,if=ibm.
cc     for ibm/360/370/3090
c      data z/z3380000000000000/
c      x02aae = z
c     for sun
      data z/1.1d-16/
      x02aae = z
c     * eps *
c     returns the value eps where eps is the smallest
c     positive
c     number such that 1.0 + eps > 1.0
c     the x parameter is not used
c     for icl 1900
c     x02aae = 2.0**(-37.0)
c+self,if=pc.
c     for pdp11
c      x02aae=2.d0**(-23.d0)
c+self.

      return
      end
c
      integer  function x02bbe(x)
      implicit double precision(a-h,o-z)
c     nag copyright 1975
c     mark 4.5 release
*     real x
c     * maxint *
c     returns the largest integer representable on the computer
c     the x parameter is not used
c     for icl 1900
c      x02bbe = 8388607
c     for ibm,sun,vax,ibm pc/386/486
       x02bbe = 2147483647
c   for pdp11
c     x02bbe=32767
      return
      end

      integer function p01aae(ifail, error, srname)
c     mark 1 release.  nag copyright 1971
c     mark 3 revised
c     mark 4a revised, ier-45
c     mark 4.5 revised
c     mark 7 revised (dec 1978)
c     returns the value of error or terminates the program.
      integer error, ifail, nout
      character*8 srname
c     test if no error detected
      if (error.eq.0) go to 20
c     determine output unit for message
      call x04aae (0,nout)
c     test for soft failure
      if (mod(ifail,10).eq.1) go to 10
c     hard failure
      write (nout,99999) srname, error
c     stopping mechanism may also differ
      stop
c     soft fail
c     test if error messages suppressed
   10 if (mod(ifail/10,10).eq.0) go to 20
      write (nout,99999) srname, error
   20 p01aae = error
      return
99999 format (1h0, 38herror detected by nag library routine , a8,
     * 11h - ifail = , i5//)
      end
      subroutine x04aae(i,nerr)
c     mark 7 release. nag copyright 1978
c     mark 7c revised ier-190 (may 1979)
c     if i = 0, sets nerr to current error message unit number
c     (stored in nerr1).
c     if i = 1, changes current error message unit number to
c     value specified by nerr.
c
c     *** note ***
c     this routine assumes that the value of nerr1 is saved
c     between calls.  in some implementations it may be
c     necessary to store nerr1 in a labelled common
c     block /ax04aa/ to achieve this.
c
c     .. scalar arguments ..
      integer i, nerr
c     ..
c     .. local scalars ..
      integer nerr1
c     ..
      data nerr1 /5/
      if (i.eq.0) nerr = nerr1
      if (i.eq.1) nerr1 = nerr
      return
      end

      subroutine exclusive_model(q2m,wm,csthcm,st,sl,stt,stl,stlp)
      implicit none
      double precision st,sl,stt,stl,stlp
      integer nq,nw,nt,iq,iw,it,nc,narg(3)
      parameter(nq=18)
      parameter(nw=47)
      parameter(nt=61)
      double precision q2,w,csthcm,dfint,ee,th_cm,degrad,
     &q2_pn(nq),w_pn(nw),th_cm_pn(nt),
     &ft_cs(nq,nw,nt),fl_cs(nq,nw,nt),ftt_cs(nq,nw,nt),
     &ftl_cs(nq,nw,nt),ftlp_cs(nq,nw,nt),
     &arg(3),rarg(1000),a2,a30,a31,a3,wcor,q2cor
      double precision q2m,wm
      data q2_pn/0.0d0,0.3d0,0.6d0,0.9d0,1.2d0,1.5d0,1.8d0,2.1d0,
     &2.4d0,2.7d0,3.0d0,3.3d0,3.6d0,3.9d0,4.2d0,4.5d0,4.8d0,5.0d0/,
     &w_pn/1.08d0,1.10d0,1.12d0,1.14d0,1.16d0,1.18d0,1.20d0,1.22d0,
     &1.24d0,1.26d0,1.28d0,1.30d0,1.32d0,1.34d0,1.36d0,1.38d0,1.40d0,
     &1.42d0,1.44d0,1.46d0,1.48d0,1.50d0,1.52d0,1.54d0,1.56d0,1.58d0,
     &1.60d0,1.62d0,1.64d0,1.66d0,1.68d0,1.70d0,1.72d0,1.74d0,1.76d0,
     &1.78d0,1.80d0,1.82d0,1.84d0,1.86d0,1.88d0,1.90d0,1.92d0,1.94d0,
     &1.96d0,1.98d0,2.00d0/,
     &th_cm_pn/0.d0,3.d0,6.d0,9.d0,12.d0,15.d0,18.d0,21.d0,24.d0,27.d0,30.d0,
     &33.d0,36.d0,39.d0,42.d0,45.d0,48.d0,51.d0,54.d0,57.d0,60.d0,63.d0,66.d0,
     &69.d0,72.d0,75.d0,78.d0,81.d0,84.d0,87.d0,90.d0,93.d0,96.d0,99.d0,102.d0,
     &105.d0,108.d0,111.d0,114.d0,117.d0,120.d0,123.d0,126.d0,129.d0,132.d0,
     &135.d0,138.d0,141.d0,144.d0,147.d0,150.d0,153.d0,156.d0,159.d0,162.d0,
     &165.d0,168.d0,171.d0,174.d0,177.d0,180.d0/,
     &nc/0/,narg/nq,nw,nt/,degrad/57.29577952d0/,
     &a2/1.15d0/,a30/-1.23d0/,a31/0.16d0/
      common/exlusive/rarg,ft_cs,fl_cs,ftt_cs,ftl_cs,ftlp_cs
c Init
      st=0.d0
      sl=0.d0
      stt=0.d0
      stl=0.d0
      stlp=0.d0
c new variables
      q2=q2m
      w=wm
      th_cm=acos(csthcm)*degrad
c Check kinematics
      if(q2.lt.0.0) then
      print*,'Warning: Q2<0 in exclusive model!'
      print*,'Using Q2=0'
      q2=0.d0
      endif
      if(q2.gt.5.0) then
c      print*,'Warning: Q2>5 GeV^2 in exclusive model!'
c      print*,'Using extrapolation from MAID2003'
      q2cor=(5.d0**a2)/(q2**a2)
      q2=5.d0
      else
      q2cor=1.d0
      endif
      if(w.lt.1.07784) return
      if(w.gt.2.0) then
c      print*,'Warning: W>2 GeV in exclusive model!'
c      print*,'Using extrapolation from MAID2003 (A.Browman PRL35,Cornell)'
      a3=a30+a31*th_cm
      if(th_cm.lt.50.0) a3=a30+a31*50.d0
      if(th_cm.gt.100.0) a3=a30+a31*100.d0
      wcor=(2.d0**a3)/(w**a3)
      w=2.d0
      else
      wcor=1.d0
      endif
      if(abs(csthcm).gt.1.0) return
c Read data from file
      if(nc.eq.0) then
      open(44,file='pi_n_maid.dat',status='old')
       do iq=1,nq
        do iw=1,nw
         do it=1,nt
      read(44,*) ft_cs(iq,iw,it),fl_cs(iq,iw,it),
     &ftt_cs(iq,iw,it),ftl_cs(iq,iw,it),ftlp_cs(iq,iw,it)
         enddo
        enddo
       enddo
      close(44)
       do iq=1,nq
	rarg(iq)=q2_pn(iq)
       enddo
       do iw=1,nw
	rarg(iw+nq)=w_pn(iw)
       enddo
       do it=1,nt
	rarg(it+nw+nq)=th_cm_pn(it)
       enddo
	nc=nc+1
      endif
c Extract interpolated cross section
c      if(q2.ge.q2_pn(1).and.w.ge.w_pn(1).and.th_cm.ge.th_cm_pn(1)
c     &.and.q2.le.q2_pn(nq).and.w.le.w_pn(nw)
c     &.and.th_cm.le.th_cm_pn(nt)) then

	arg(1)=q2
	arg(2)=w
	arg(3)=th_cm

      st=dfint(3,arg,narg,rarg,ft_cs)*wcor*q2cor
      sl=dfint(3,arg,narg,rarg,fl_cs)*wcor*q2cor
      stt=dfint(3,arg,narg,rarg,ftt_cs)*wcor*q2cor
      stl=dfint(3,arg,narg,rarg,ftl_cs)*wcor*q2cor
      stlp=0.d0
c      stlp=dfint(3,arg,narg,rarg,ftlp_cs)*wcor*q2cor

c	else
c      st=0d0
c      sl=0d0
c      stt=0d0
c      stl=0d0
c      stlp=0d0
c	endif
	return
	end
