c***********************************************************************
c        G R O W T H : 3-D gravity inversion           10/09/2018      *
c                                                                      *
c        Compilation with Intel Visual Fortran,                        *
c            (Quickwin Application Project)                            *
c                                                                      *
c  Aim:      3-D gravity inversion searching for anomalous structures  *
c            with prescribed positive and negative density contrasts   *
c            in an exploratory process                                 *
c                                                                      *
c  Author:    Antonio G. Camacho    (antonio_camacho@mat.ucm.es)       *
c             Instituto de Astronomia y Geodesia (CSIC-UCM)            *
c             Facultad de CC. Matematicas                              *
c             Ciudad Universitaria,     28040 Madrid (Spain)           *
c                                                                      *
c  Reference: Camacho et al. "The 3-D gravity inversion package        *
c             GROWTH2.0 and its application to Tenerife Island, Spain" *
c               Computers & Geosciences 37 (2011) 621–633              *
c               See also the file GROWTH.TXT                           *
c ---------------------------------------------------------------------*
c                                                                      *
c  INPUT                                                               *
c  -----                                                               *
c  file "GRA.txt": coordinates, altitudes, observed gravity values     *
c        and relative errors for the gravity stations :                *
c         xp,yp,zp : planar coordinates (metres) and altitude (metres) *
c         dg: observed gravity values (in microgals)                   *
c         eg: stand.dev. values for the observed gravity               *
c         terco: applied terrain correction for unit density           *
c                                                                      *
c  OUTPUT                                                              *
c  ------                                                              *
c  file "MOD.txt": Resulting 3-D model of the subsurface anomalous     *
c        structures as given by means of the 3-D grid of cells         *
c        Values for several characteristic parameters of the model     *
c                                                                      *
c  file "FIL.txt": regional, local, calculated and residual gravity    *
c        values resulting from inversion                               *
c                                                                      *
c ---------------------------------------------------------------------*
c  DIMENSIONS                                                          *
c  ----------                                                          *
c   ms:  max. number of gravity stations                               *
c   mc:  max. number of prismatic cells for the 3-D grid               *
c   mi:  max. number of horizontal levels                              *
c   mp:  number of steps to estimate a growth trend                    *       
c***********************************************************************

      use IFQWIN
      use IFPORT

c  Dimensions

      parameter (ms=5000,mc=99000,mi=50,mh=20,mp=300,mr=20)
      type (xycoord) xy
      dimension 
     & xc(mc),yc(mc),zc(mc),tx(mc),ty(mc),tz(mc),df(mc),far(mc),
     & at(mc),ireg(mc),vol(mc),dec(mc),
     & rr(ms),g(ms),gc(ms),x(ms),y(ms),z(ms),re(ms),t(ms),w(ms),w0(ms),
     & xn(ms),yn(ms),zn(ms), 
     & jca(mi),nus(mh),  de(2),zdm(2),vdm(2),szdm(2),
     & jes(mh),jdi(mh),jd(mh),jz(mh),apla(mh),is(mh),cov(mr),swc(mr)
     &   , wme(mh,2),nwe(mh,2),gpe(mp)
      integer m(mc),nada
      integer*2 u1,u2,v1,v2
      real*8 gxx,gyy,gxg,gyg,gxz,gyz,gzz,gzg,ryy,ryg,ryz,rzz,rzg,
     & px,py,pz,f1g,f11,srg,su,dg,f  
      character*50 obs,mod,fil
      character*120 texto
      character*9 hoy
      character*1 ts1,ts2,ts3,ts4,ts5,ts6,bl,as
      data bl/' '/,as/'*'/

      call seed(2016)
      nada=setbkcolorrgb(#f0f0f0)
      call clearscreen($gclearscreen)
      call date(hoy)
      nada=initializefonts()
      nada=settextcolorrgb(#000000)

c     Reading file names and default parameters

      ka=8000                     ! fator de forma exponencial fijo 8000
      fam=4.
      sla=3.
      fla=40.
      uw=0
      na0=20
      seli=40
      blun=3.5
      esp=2.
      call dial1(obs,mod,fil,kp,na0,fam,sla,fla,uw)

      open(6,file='out.dat')
      
c     Reading position of bechmarks 

      open(1,file=obs)
      ax=9.d9
      bx=-ax
      ay=ax
      by=-ay
      i=0
   71 read(1,*,err=72,end=72) xp,yp,zp
      if(xp.eq.0.and.yp.eq.0.and.zp.eq.0) go to 72
      i=i+1
      if(xp.lt.ax) ax=xp
      if(xp.gt.bx) bx=xp
      if(yp.lt.ay) ay=yp
      if(yp.gt.by) by=yp
      go to 71
   72 nms=i   
      close(1)
      db=(by-ay)*(bx-ax)/10000.

c     Reading anomaly data      
      
   70 open(1,file=obs)
      xd=0.
      yd=0.
      zd=0.
      gm=0.
      gmi=9.d9
      ns=0
   73 read(1,*,err=75,end=75) xp,yp,zp,dg,eg
      if(xp.eq.0.and.yp.eq.0.and.zp.eq.0) go to 75
      if(ns.gt.0) then
        do 74 j=1,ns
        px=xp-x(j)
        py=yp-y(j)
        d2=px*px+py*py
        if(d2.le.db) go to 73
   74   continue       
      endif
      ns=ns+1
      x(ns)=xp
      y(ns)=yp
      z(ns)=zp  
      g(ns)=dg
      t(ns)=1
      w0(ns)=eg
      xd=xd+xp
      yd=yd+yp
      zd=zd+zp
      gm=gm+dg
      if(dg.lt.gmi) gmi=dg
      if(ns.lt.ms) go to 73
             close(1)
             db=db*1.2
             go to 70
   75 if(ns.gt.1) go to 76
        write(texto,'(a)') ' >> Error in anomaly data file'
        call moveto(150,200,xy)
        call outgtext(texto)
        stop
   76 close(1)
      
c    Normalized coordinates
      
      gm=gm/ns
      xd=xd/ns
      yd=yd/ns
      zd=zd/ns
      ixm=nint(xd/10.)*10.
      iym=nint(yd/10.)*10.
      izm=nint(zd/10.)*10.
      fe=(bx-ax+by-ay)/2./30000.
      techo=-9.d9     
      do 77 i=1,ns
      x(i)=x(i)-ixm
      y(i)=y(i)-iym
      z(i)=z(i)-izm
      zn(i)=z(i)/fe
      if(zn(i).gt.techo) techo=zn(i)    
      xn(i)=x(i)/fe
   77 yn(i)=y(i)/fe
      ax=(ax-ixm)/fe
      bx=(bx-ixm)/fe
      ay=(ay-iym)/fe
      by=(by-iym)/fe

c    Read geometrical parameters for the 3D grid

      pp=500
      base=-12000
   10 iex=3
      call dial2(ns,izm,fe,techo,base,pp)

c   Creating a 3D grid of cells 

      avm=pp*pp*pp*pp*1.0d-10   
      nc=0
      kr=0
      nada=setfont('t''Courier''h17w7')
   80 ncr=nc
      kr=kr+1
      zs=techo
      ax=ax-4000
      ay=ay-4000
      bx=bx+4000
      by=by+4000
   82 dz=pp*0.8
      call especapa(zs,avm,ms,ns,xn,yn,zn,dz) 
      if(dz.gt.1) go to 81
         zs=zs-pp*0.5
         if(zs.gt.base) go to 82
   81 if(kr.eq.2) dz=dz*1.7   ! 1.7
      if(kr.eq.3) dz=dz*2.0   ! 2.0
      call celdascapa(ax,bx,ay,by,zs,dz,avm,ms,ns,xn,yn,zn,
     & mc,nc,xc,yc,zc,tx,ty,tz,df,seli,ncr,ireg)
      if(nc.lt.0) then 
        write(*,'(////80x,a,i6)') '***** Warning: Num. of cells >',mc
        write(*,'(/84x,a,i6)') '>> Try a smaller length for cells'
        go to 10
      endif
      zs=zs-dz                      
      if(zs.gt.base) go to 82  
      if(kr.lt.iex) go to 80   
      call dial3(ns,izm,fe,techo,base,pp,nc,ire)
      if(ire.eq.1) go to 10
      do 5 i=1,nc
      xc(i)=xc(i)*fe
      yc(i)=yc(i)*fe
      zc(i)=zc(i)*fe
      tx(i)=tx(i)*fe
      ty(i)=ty(i)*fe
    5 tz(i)=tz(i)*fe
      siz=nint(pp*fe)

c  Some initial preparation for Cells and Gravity data
    
      ztop=-9.d9
      srr=0.
      d=0.
      disg=0
      do 6 i=1,ns
      d=d+x(i)*x(i)+y(i)*y(i)+z(i)*z(i)     
      if(z(i).gt.ztop) ztop=z(i)
      if(kp.eq.0) w0(ns)=1.
      if(kp.eq.1) w0(ns)=1./sqrt(w0(ns)+0.01)
      srr=srr+w0(ns)
      disg=disg+(g(i)-gm)**2
    6 if(t(i).gt.1) g(i)=g(i)/t(i)
      disg=sqrt(disg/ns)
      gmi=(gmi+gm)/2.  
      srr=srr/ns
      fe=sqrt(d/ns)/8000.
      tic=2000.*fe     
      call round(tic,r)   
      tic=r
      do 9 j=1,nc
    9 vol(j)=tx(j)*ty(j)*tz(j)
      i=0
      avm=siz/50. 
      fsiz=siz/fe
      fsl=1.e13/fsiz/fsiz/fsiz
      zme=-2500*fe
      ex0=0.
      zdm(1)=0
      zdm(2)=0 
      na1=200 
      
c  Reading options and parameters for inversion

      call step(ms,ns,x,y,t,dr)
      dr=20000*fe/sqrt(ns*1.)
      kwl=4*dr
      gm0=gmi
      call dial4(obs,mod,fil,kp,na0,kwl,nc,ns,fla,de,
     & ig0s,gm0,igxys,gx0,gy0,ict,jla,sla,jden,fam,uw)
c      uw=fam/3.
      if(jla.eq.1) sla=0.
      if(jla.eq.1) fam=0.
      if(jden.eq.1) de(1)=-1
      if(jden.eq.1) de(2)= 1      
      fls=fla*fla*fe*fe*1.e2      !*1.e3 /fsiz*900.  * avm**esp  
      fas=fam
      
c  Calculating design matrix A

      r=1.e6*fe*fe
      do 11 j=1,nc
      sv=0.
      c=0.
      sp=0
      do 12 i=1,ns
      xd=x(i)-xc(j)
      yd=y(i)-yc(j)
      zd=z(i)-zc(j)                             ! spherical aprox.
      d=xd*xd+yd*yd+zd*zd+10.*fe 
      av2=zd*zd/d/d/d
      av=sqrt(av2)*6.672e-3*vol(j)
      av=av/avm
      c=c+av*av
      sv=sv+av2
   12 continue
      dec(j)=c/ns                                     ! atraccion normalizada
      c=3.d-7/sqrt(sv/ns)/fe/fe                       ! sensit. of cells
      if(c.gt.seli.and.ireg(j).eq.0) ireg(j)=1
   11 continue   
             
c  Terrain correction option

      if(ict.eq.1) then
       open(2, file=obs)
       do 13 i=1,ns
       read(2,*,err=92) xp,yp,zp,dg,eg,terco
   13  rr(i)=41.9214*zp-terco
       close(2)
       p2=kwl*kwl
       
       do 15 i=1,ns
       c=0
       ne=0
       do 14 j=1,ns
       xd=x(i)-x(j)
       yd=y(i)-y(j)
       d=xd*xd+yd*yd
       if(d.gt.p2) go to 14
       c=c+rr(j)
       ne=ne+1
   14  continue
   15  t(i)=rr(i)-c/ne
      endif
      
c     Initial trend     
      
      do 16 i=1,ns
      g(i)=g(i)-gm0-gx0*x(i)/1000.-gy0*y(i)/1000.
   16 w0(i)=w0(i)/srr
      krk=0

c     Optional Stratification

      nes=0
      bot=izm-8000*fe 
      if(jla.eq.0) go to 20
      nes=8
      q1n=520
      q2n=-12
      q1p=0
      q2p=40
      fas=4.
      ra=2.0
      call dial5(nes,fas,jind,bot,ra,q1n,q2n,q1p,q2p)
      bot=bot-izm
      if(jbot.eq.1) jind=0
      if(nes.gt.mh) nes=mh
      do 17 i=1,mh
      jdi(i)=de(2)
   17 apla(i)=fas
      kes=1  
      jca(1)=zc(1)+tz(1)/2.
      do 19 j=2,nc
      c=zc(j)-tz(j)/2. 
      do 18 k=1,kes
      r=(jca(k)-c)
      if(r.lt.2) go to 19 
   18 continue  
      kes=kes+1
      jca(kes)=c
   19 continue
      fapl=fsiz**0.33/fe*1.e-5
      
c  ------Graphics and headings ------------------------------------------
   20 nada=setcolorrgb(#f0f0f0)                                         !
      nada=rectangle($gfillinterior,10,60,600,305)                      ! cuadro total
      nada=setcolorrgb(#f8f8f8)                                         !
      nada=rectangle($gfillinterior,1265,330,1295,345)                  !
      nada=setcolorrgb(#000000)                                         !
      nada=initializefonts()                                            !
      write(texto,'(a)') 'GROWTH - 3D Gravity inversion'                !
      nada=setfont('t''Times New Roman''h23w11')                        !
      call moveto(100,30,xy)                                            !
      call outgtext(texto)                                              !
      nada=setfont('t''Courier''h17w7')                                 !
      write(texto,'(a)') hoy                                            !
      call moveto(1210,40,xy)                                           !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'v. 3.0'                                       !
      call moveto(510,35,xy)                                            !
      call outgtext(texto)                                              !
      ix=1183                                                           !
      write(texto,'(a,i4,a)') 'Mean reliev:',izm,' m'                   !
      call moveto(ix,330,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Sea level'                                    !
      call moveto(ix,352,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Neg. density contr.'                          !
      call moveto(ix,376,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Pos. density contr.'                          !
      call moveto(ix,399,xy)                                            !
      call outgtext(texto)                                              !
      call setgtextrotation(900)                                        !
      write(texto,'(a)') 'A.G. Camacho 2020'                            !
      call moveto(1030,300,xy)                                          !
      call outgtext(texto)                                              !
      write(texto,'(a,i5,a)') 'Tic=',nint(tic),' m'                     !
      call moveto(600,320,xy)                                           !
      call outgtext(texto)                                              !
      call moveto(600,600,xy)                                           !
      call outgtext(texto)                                              !
      call setgtextrotation(000)                                        !
      write(texto,'(a,7x,a)') ' Number of',                             !
     -'Dens.   Resid. Offset   Regional     Depth cdm'                  !
      call moveto(60,150,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a,a)') 'Steps Cells Out  ',                         !
     -'(kg/m3) (uGal)  (uGal)  (uGal/km)      - (m) +'                  !
      call moveto(55,165,xy)                                            !
      call outgtext(texto)                                              !
      call ventana(45,187,550,205)                                      !results
      call ventana(620,470,1030,640)                                    !WE
      call ventana(1060,470,1470,640)                                   !SN      
      call ventana(620,20,1030,420)                                     !horiz
      call ventana(30,469,570,680)                                      !profun vs. mass
      fdib=80                                                           !
c      if(ns.gt.1000) fdib=70                                            !
c      if(ns.gt.2000) fdib=60                                            !

      j=500+(izm/fe-techo)/fdib
      j=500
      nada=setcolorrgb(#dfffdf)                                         !verde
      nada=rectangle($gfillinterior,620,20,1030,420)                    !hor
      nada=rectangle($gfillinterior,620,j,1030,640)                     !WE
      nada=rectangle($gfillinterior,1060,j,1470,640)                    !SN

      ix=1140                                                           !
      nada=rectangle($gfillinterior,ix,335,ix+30,425)                   !muestra
      nada=setcolorrgb(#ffffff)                                         !
      nada=rectangle($gfillinterior,ix,315,ix+30,335)                   !muestra      
      nada=setcolorrgb(#f08000)                                         !
      nada=rectangle($gfillinterior,ix+10,380,ix+20,385)                !blue
      nada=setcolorrgb(#0040f0)                                         !
      nada=rectangle($gfillinterior,ix+10,405,ix+20,410)                !red
      nada=setcolorrgb(#e5e5e5)                                         !
      nada=rectangle($gfillinterior,30,j,570,680)                       !mass

      u1=500+izm/fe/fdib                                                !sea levels
      nada=setcolorrgb(#f0f000)                                         !
      call moveto(620,u1,xy)                                            !
      nada=lineto(1030,u1)                                              !
      call moveto(1050,u1,xy)                                           !
      nada=lineto(1460,u1)                                              !
      call moveto(1140,360,xy)                                          !
      nada=lineto(1170,360)                                             !
      call moveto(30,u1,xy)                                             !
      nada=lineto(570,u1)                                               !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i5,a)') nint(tic),' m'                              !
      call moveto(1347,350,xy)                                          !
      call outgtext(texto)                                              !
      call moveto(1350,375,xy)                                          ! escala
      nada=lineto(1400,375)                                             !
      call moveto(1360,370,xy)                                          !
      nada=lineto(1360,375)                                             !
      u1=1360+tic/fe/fdib                                               !
      call moveto(u1,370,xy)                                            !
      nada=lineto(u1,375)                                               !
      j=izm/tic+2                                                       !
      do 21 i=1,20                                                      !
      v1=20+(i-1)*tic/fe/fdib                                           !
      call moveto(620,v1,xy)                                            !
      if(v1.lt.420) nada=lineto(624,v1)                                 !left horiz
      call moveto(1030,v1,xy)                                           !
      if(v1.lt.420) nada=lineto(1026,v1)                                !right horiz
      v1=500+((i-j)*tic+izm)/fe/fdib                                    !
      call moveto(620,v1,xy)                                            !
      if(v1.lt.640) nada=lineto(624,v1)                                 !
      call moveto(1030,v1,xy)                                           !
      if(v1.lt.640) nada=lineto(1025,v1)                                !
      call moveto(1060,v1,xy)                                           !
      if(v1.lt.640) nada=lineto(1064,v1)                                !
      call moveto(1470,v1,xy)                                           !
      if(v1.lt.640) nada=lineto(1465,v1)                                !
      call moveto(570,v1,xy)                                            !
      if(v1.lt.680) nada=lineto(564,v1)                                 !right mass
      call moveto(30,v1,xy)                                             !
      if(v1.lt.680) nada=lineto(35,v1)                                  !left mass
      u1=620+(i-1)*tic/fe/fdib                                          !  
      call moveto(u1,20,xy)                                             !top horiz
      if(u1.lt.1030) nada=lineto(u1,25)                                 !
      call moveto(u1,420,xy)                                            !bot horiz
      if(u1.lt.1030) nada=lineto(u1,415)                                !
      call moveto(u1,470,xy)                                            !top WE
      if(u1.lt.1030) nada=lineto(u1,475)                                !
      u1=1060+(i-1)*tic/fe/fdib                                         !
      call moveto(u1,470,xy)                                            !top NS
      if(u1.lt.1470) nada=lineto(u1,475)                                !
   21 continue                                                          !
      lx=29                                                             !
      ly=240                                                            !
      mx=522                                                            !
      my=38                                                             !
      write(texto,'(a)') 'Filling grid cells'                           !
      call moveto(lx+162,ly+my+50,xy)                                   !
      call outgtext(texto)                                              !
      call setgtextrotation(900)                                        !
      write(texto,'(a)') 'Dens.Scale'                                   !
      call moveto(lx-22,ly+my+35,xy)                                    !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Fit'                                          !
      call moveto(lx-22,ly+my+110,xy)                                   !
      call outgtext(texto)                                              !
      call setgtextrotation(000)                                        !
      v1=ly+160                                                         !
      call moveto(lx,v1,xy)                                             ! eje hor offset
      nada=lineto(570,v1)                                               !
      call moveto(lx,v1+5,xy)                                           ! eje vert offset
      nada=lineto(lx,v1-60)                                             !
      nada=rectangle($gfillinterior,16,v1-60,17,v1-40)                  ! muestra fit
      write(texto,'(a)') '0%'                                           !
      call moveto(lx-2,ly+my+43,xy)                                     !
      call outgtext(texto)                                              !
      write(texto,'(a)') '20%'                                          !
      call moveto(lx+mx*2/7-10,ly+my+43,xy)                             !
      call outgtext(texto)                                              !
      write(texto,'(a)') '40%'                                          !
      call moveto(lx+mx*4/7-10,ly+my+43,xy)                             !
      call outgtext(texto)                                              !
      write(texto,'(a)') '60%'                                          !
      call moveto(lx+mx*6/7-10,ly+my+43,xy)                             !
      call outgtext(texto)                                              !
      nada=setcolorrgb(#f0f0f0)                                         !
      u1=lx-1                                                           !
      v1=ly-my+38                                                       !
      u2=lx+mx+20                                                       !
      v2=ly+my+42                                                       !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        ! dens
      nada=setcolorrgb(#000000)                                         !
      nada=rectangle($gborder,u1,v1,u2,v1)                              ! dens recuadro
      nada=rectangle($gborder,u1,v2,u2,v2)                              !
      nada=rectangle($gborder,u1,v1,u1,v2)                              !
      u1=lx+mx*2./7.                                                    !
      nada=rectangle($gborder,u1,v1,u1,v2)                              !
      u1=lx+mx*4./7.                                                    !
      nada=rectangle($gborder,u1,v1,u1,v2)                              !
      u1=lx+mx*6./7.                                                    !
      nada=rectangle($gborder,u1,v1,u1,v2)                              !
      call ventana(545, 82,572, 94)                                     !
      call ventana(104,112,131,124)                                     !
      call ventana(237,112,268,124)                                     !
      call ventana(391,112,423,124)                                     !
      call ventana(546,112,576,124)                                     !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a,i6,i5,5x,a,i7,4x,a,i5,5x,a,i4)')                  !
     & 'Num.data:',nms,ns,'Num.cells:',nc,'Cell side(m):',nint(siz),    !
     & 'Rand:',na0                                                      !
      call moveto(25,80,xy)                                             !
      call outgtext(texto)                                              !
      write(texto,'(a,i4,4x,a,f5.1,4x,a,f5.1,4x,a,f5.1)')               !
     & 'Smoothness:',nint(fla),'Flattening:',fas,                       !
     & 'Lateral.Coef:',sla,'Up-Weighting:',uw                           !
      if(jla.ne.0) write(texto,'(a,i4,4x,a,f5.1,4x,a,i5,4x,a,f5.1)')    !
     & 'Smoothness:',nint(fla),'Flattening:',fas,                       !
     & 'Numb. layers:',nes,'Up-Weighting:',uw                           !
      call moveto(25,110,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Plan view'                                    !
      call moveto(800,430,xy)                                           !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Elevation view W - E'                         !
      call moveto(750,650,xy)                                           !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Elevation view S - N'                         !
      call moveto(1200,650,xy)                                          !
      call outgtext(texto)                                              !
      j=410                                                             !
      write(texto,'(a)') 'Offset     WE     SN'                         !
      call moveto(405,j-8,xy)                                           !
      call outgtext(texto)                                              !
      nada=setcolorrgb(#d000d0)                                         !
      nada=rectangle($gfillinterior,450,j,470,j+1)                      !
      nada=setcolorrgb(#00d0d0)                                         !
      nada=rectangle($gfillinterior,500,j,520,j+1)                      !
      nada=setcolorrgb(#00c000)                                         !
      nada=rectangle($gfillinterior,550,j,570,j+1)                      !
      j=660                                                             !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'c.d.m. zm-      zm+'                          !
      call moveto(55,j-8,xy)                                            !
      call outgtext(texto)                                              !
      nada=setcolorrgb(#0040f0)                                         !
      nada=rectangle($gfillinterior,130,j,150,j+1)                      !
      nada=setcolorrgb(#f08000)                                         !
      nada=rectangle($gfillinterior,190,j,210,j+1)                      !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'Initial steps . . . . . .'                    !
      call moveto(60,189,xy)                                            !
      call outgtext(texto)                                              !
c-----------------------------------------------------------------------!
     
c     Adjust of bottom for stratified model
      
      if(jla.eq.0)  go to 25    
      do 22 i=1,nes
      c=(i-1.)/(nes-1.)
   22 jz(i)=ztop-(exp(c*ra)-1.)*(ztop-bot)/(exp(ra)-1.)+izm
      if(jind.eq.1) then
        call dial6(de(2),mh,nes,jz,jdi,apla)
        if(jdi(9).eq.0) nes=8
        if(jdi(8).eq.0) nes=7
        if(jdi(7).eq.0) nes=6
        if(jdi(6).eq.0) nes=5
        if(jdi(5).eq.0) nes=4
        if(jdi(4).eq.0) nes=3
        if(jdi(3).eq.0) nes=2
        if(jdi(2).eq.0) nes=1
      endif
      do 24 i=1,nes
      jd(i)=jdi(i)
      jz(i)=jz(i)-izm
      r=9.d9  
      do 23 k=1,kes
      c=abs(jz(i)-jca(k))
      if(c.ge.r) go to 23   
      kk=k
      r=c   
   23 continue
      jz(i)=jca(kk) 
      jes(i)=kk  
      is(i)=kk  
   24 continue

c    More initial values

   25 krk=krk+1
      ig0=ig0s
      igxy=igxys
      do 26 i=1,ns
      gc(i)=0.
   26 w(i)=w0(i)
      nbf=0
      nb0=0
      if(nes.gt.0) then 
       do 27 i=1,nes 
       wme(i,1)=0
       wme(i,2)=0
       nwe(i,1)=0
       nwe(i,2)=0
   27  jz(i)=jca(jes(i))
      endif   
      peme=0
      npe=0
      do 29 j=1,nc
      at(j)=dec(j)  *exp((zme-zc(j))/fe/8000*uw) 
      df(j)=0.
      m(j)=1
      if(nes.le.1) go to 29  
      do 28 k=1,nes-1
   28 if(nes.ge.k.and.zc(j).lt.jz(k)) m(j)=k+1
   29 continue
      ryz=0.
      rzz=0.
      rzg=0.
      pxs=0
      pys=0
      pzs=0
      nu=0
      ncf=0
      mm=0
      sm=0.
      nk=0
      gxmi=9.d9
      gymi=9.d9
      fit0=999.d9
      fs1=9.d9
      fsa=9.d9
      fl=fls
      szdm(1)=0
      szdm(2)=0
      vdm(1)=0.1
      vdm(2)=0.1
      do 47 i=1,mp
   47 gpe(i)=0.
      pen=0.
      nbum=0
      denmi=9.d9
         kre=0
         tkre=1
         
c-----------------------------------------------------------------------!
c       Step by step growth of the anom anomalous structure             !
c-----------------------------------------------------------------------!

   30 jend=1 
      ncc=nc/na1
      c=nu+1
      if(nu.lt.na1) ncc=nc/c
      st0=1.
      if(sla.gt.0.01) st0=1.-(c*sla/fsl)**1    ! lateral smoothing  **2
      if(st0.le.0.001) go to 50
      if(fam.gt.0.01) ex0=fam*exp(zme/fe/ka)
      
   31 em=0.
      fit=fit0*1.0005 
      if(nu.lt.50) fit=fit0*1.1
      s11=0.
      s1x=0.
      s1y=0.
      s1g=0.
      sxx=0.
      sxy=0.
      syy=0.
      sxg=0.
      syg=0.
      s1z=0.
      sxz=0.
      syz=0.
      szz=0.
      szg=0.
      do 32 i=1,ns
      p2=w(i)
      xd=x(i)/1000.
      yd=y(i)/1000.
      dg=g(i)
      s11=s11+p2
      s1x=s1x+xd*p2
      s1y=s1y+yd*p2
      s1g=s1g+dg*p2
      sxx=sxx+xd*xd*p2
      sxy=sxy+xd*yd*p2
      sxg=sxg+xd*dg*p2
      syy=syy+yd*yd*p2
      syg=syg+yd*dg*p2
      if(ict.eq.0) go to 32
      tt=t(i)
      s1z=s1z+tt*p2
      sxz=sxz+xd*tt*p2
      syz=syz+yd*tt*p2
      szg=szg+tt*dg*p2
      szz=szz+tt*tt*p2
   32 continue
      js=0
      do 39 j1=1,ncc
      j=j1
      if(ncc.lt.nc) call RANDOM_NUMBER(ranval)
      if(ncc.lt.nc) j=(nc-1)*ranval+1
      if(m(j).eq.0) go to 39
      cmi=0.
      do 38 k=1,2                                            ! 1: neg 2 pos
      if(fsa.lt.abs(de(k))) go to 38
      if(fam.gt.0.01) cmi=fam*exp(zc(j)/fe/ka)-ex0           ! downward continuous increase
        st=st0
        if(ireg(j).eq.1) st=1. 
        db=2*k-3 + cmi                                       ! downward increase  
        db=db*st                                             ! lateral smooth
        if(k.eq.1.and.db.ge.0) go to 38
        if(k.eq.2.and.db.le.0) go to 38
      p2=1.  
       if(nes.ge.1) then                                      ! layered model
          jc=m(j)+k-2                                         ! superficie que pasaria
          if(jc.le.0.or.jc.gt.nes) go to 38 
          kk=2*k-3  
          db=kk*jd(jc)
          c=(zc(j)-jz(jc))*kk
          if(c.lt.0) go to 38                 
            if(k.eq.1) r=q1n*fe-jz(jc)*q2n/100. 
            if(k.eq.2) r=q1p*fe-jz(jc)*q2p/100. 
            p2=exp(fapl*(c-r)*apla(jc))   
        endif
      if(k.eq.1.and.db.gt.0) db=0.
      if(k.eq.2.and.db.lt.0) db=0.
      if(db.eq.0.) go to 38
      px=pxs
      py=pys
      pz=pzs
      dg=dgs
      srr=0.
      srg=0.
      srx=0.
      sry=0.
      srz=0
      s1r=0.
      do 33 i=1,ns
      xd=x(i)-xc(j)
      yd=y(i)-yc(j)
      zd=z(i)-zc(j)                             ! spherical aprox.
      d=xd*xd+yd*yd+zd*zd+10.*fe 
      av2=zd*zd/d/d/d
      av=sqrt(av2)*6.672e-3*vol(j)

      re(i)=gc(i)+av*db    !+A(i,j)*db     
      zb=re(i)*w(i)
      s1r=s1r+zb
      srg=srg+zb*g(i)
      srx=srx+zb*x(i)/1000.
      sry=sry+zb*y(i)/1000.
      if(ict.eq.1) srz=srz+zb*t(i)
   33 srr=srr+zb*re(i)
      sv=sm+at(j)*db*db *p2    ! uGal*uGal
      su=srr+fl*sv
      f11=s11*su-s1r*s1r
      f1x=s1x*su-s1r*srx
      f1y=s1y*su-s1r*sry
      f1g=s1g*su-s1r*srg
      fxx=sxx*su-srx*srx
      fxy=sxy*su-srx*sry
      fxg=sxg*su-srx*srg
      fyy=syy*su-sry*sry
      fyg=syg*su-sry*srg
      if(ict.ne.1) go to 34
        f1z=s1z*su-s1r*srz
        fxz=sxz*su-srx*srz
        fyz=syz*su-sry*srz
        fzz=szz*su-srz*srz
        fzg=szg*su-srz*srg
   34 gxx=f11*fxx-f1x*f1x
      gxy=f11*fxy-f1x*f1y
      gxg=f11*fxg-f1x*f1g
      gyy=f11*fyy-f1y*f1y
      gyg=f11*fyg-f1y*f1g
      if(ict.ne.1) go to 35
       gxz=f11*fxz-f1x*f1z
       gyz=f11*fyz-f1y*f1z
       gzz=f11*fzz-f1z*f1z
       gzg=f11*fzg-f1z*f1g
   35 ryy=gyy*gxx-gxy*gxy
      ryg=gyg*gxx-gxy*gxg
      if(ict.ne.1) go to 36
       ryz=gyz*gxx-gxy*gxz
       rzz=gzz*gxx-gxz*gxz
       rzg=gzg*gxx-gxz*gxg
   36  f=rzz*ryy-ryz*ryz
      if(ict.eq.1.and.f.ne.0.) pz=(rzg*ryy-ryz*ryg)/f
      if(igxy.eq.1) py=(ryg-pz*ryz)/ryy
      if(igxy.eq.1) px=(gxg-py*gxy-pz*gxz)/gxx
      if(ig0. eq.1) dg=(f1g-px*f1x-py*f1y-pz*f1z)/f11
      f=(srg-dg*s1r-px*srx-py*sry-pz*srz)/su
      if(f.lt.0.9999) go to 38
      sv=f*f*sv*fl
      em=0
      do 37 i=1,ns
      c=g(i)-dg-px*x(i)/1000.-py*y(i)/1000.-pz*t(i)-re(i)*f        ! uGal 
   37 em=em+c*c*w(i)                                   ! uGal*uGal (residuo medio al cuadrado)
      c=sv+em 
      if(c.ge.fit) go to 38
      fit=c
      fs=f
      dgs=dg
      pxs=px
      pys=py
      pzs=pz
      js=j
      dbs=db
      ks=k
      pes=p2
      svs=sv
      sts=st
   38 continue
   39 continue
      
c     Growth does stop
      
      if(js.eq.0) go to 50
      if(fs.le.1.0) go to 50                         
      if(nu.eq.0.and.fs.lt.1.4) go to 90
      
c     Selected values    

      nu=nu+1
      if(fit.lt.fit0) fit0=fit
      if(df(js).eq.0) ncf=ncf+1
      do 41 i=1,ns
         xd=x(i)-xc(js)
         yd=y(i)-yc(js)
         zd=z(i)-zc(js)                       
         d=xd*xd+yd*yd+zd*zd+10.*fe 
         av=zd/sqrt(d*d*d)*6.672e-3*vol(js)
   41 gc(i)=gc(i)+av*dbs   !A(i,js)*dbs
      jc=m(js)+ks-2
      m(js)=0
      sm=sm+at(js)*dbs*dbs  *pes  
      peme=peme+pes
      npe=npe+1    
      
      if(ireg(js).eq.0) then                                 ! local
         c=vol(js)*abs(dbs)                                  ! mass
         vdm(ks)=vdm(ks)+c 
         szdm(ks)=szdm(ks)+zc(js)*c
         zdm(ks)=szdm(ks)/vdm(ks)
      endif
      
      if(nes.gt.0.and.ireg(js).eq.0) then                    ! local
         m(js)=jc+ks-1 
         fl=fls/peme*npe    
         wme(jc,ks)=wme(jc,ks)+pes
         nwe(jc,ks)=nwe(jc,ks)+1
      endif
      df(js)=df(js)+dbs
      far(js)=fs
      if(df(js).eq.0) far(js)=0
      if(far(js).gt.99999.) far(js)=99999.
      
c     Modelled and residual values 
      
      ems=0.
      c=0
      do 42 i=1,ns
      tt=g(i)-x(i)/1000.*pxs-y(i)/1000.*pys-t(i)*pzs-dgs
      re(i)= tt-gc(i)
      rr(i)=abs(tt-fs*gc(i))
      c=c+w(i)
   42 ems=ems+w(i)*rr(i)*rr(i) 
      ems=sqrt(ems/c)
      tt=fs
      call dmedian(ms,ns,rr,em)
      em=em/0.6745*blun*tt*tt *tt
      mm=0
      dg=0.
      do 43 i=1,ns
   43 if(w(i).gt.dg) dg=w(i)
      do 44 i=1,ns
      c=rr(i)/em
      if(c.le.1.) go to 44
       c=c*c-1.
       mm=mm+1
       c=dg/(1.+c*c)
       if(c.lt.w(i)) w(i)=c
   44 continue
      
      if(nu.eq.1) fs1=fs
      fsa=fs
      do 45 i=1,mp-1
   45 gpe(i)=gpe(i+1)
      gpe(mp)=dgs
      c=0.
      do 46 i=1,mp-1
   46 c=c+(gpe(i+1)-gpe(i))       
      c=c/mp
      if(nu.le.mp) c=0.

          kre=kre+1
          igxy=0
          ig0=0
          if(kre.ge.tkre) then
            if(igxys.eq.1) igxy=1
            if(ig0s .ne.0) ig0=1
            kre=0
            tkre=tkre
          endif     
       if(c*pen.lt.0.and.ig0.eq.1) tkre=tkre+nc/1000.     ! se detiene la regional     
       if(c*pen.lt.0.and.ig0.eq.1) tkre=tkre+nc/1000.     ! se detiene la regional

c      if(c*pen.lt.0.and.ig0.eq.1) igxy=0     ! se detiene la regional     
c      if(c*pen.lt.0.and.ig0.eq.1) ig0=0      ! se detiene la regional

      pen=c
      
      if(abs(fs*sts).lt.denmi) denmi=abs(fs*sts) 
           c=(zdm(1)+zdm(2))/2.+fam*20*fe
           if(c.gt.zme) zme=c
       if(krk.lt.2) then 
          if(nu.ge.nc/100.) then
            call str(fe,mc,nc,xc,yc,zc,tx,ty,tz,df,ireg,fdib,pen,pep,0)
            na1=na0
            go to 25
          endif 
          go to 30
      endif    
c---------Graphics-------------------------------------------------------
      write(texto,'(i7,i6,i4,i8,i7,i9,2i6,i8,i6)')                      !
     & nu,ncf,nm,nint(fs*(2*ks-3)*sts),nint(ems),                       !
     & nint(dgs+gm0),nint(gx0+pxs),nint(gy0+pys)                        !
     & ,nint(zdm(1)+izm),nint(zdm(2)+izm)                               !
      if(nu.le.200) gxmi=pxs                                            !
      if(nu.le.200) gymi=pys                                            !
      if(nu.le.200) gmi =dgs                                            !
      if(nu.le.200) zmi =(zdm(1)+zdm(2))/2./fe                          !
      nada=setcolorrgb(#ffffff)                                         !
      nada=rectangle($gfillinterior,45,189,520,205)                     !  blanco en ventana
      nada=setcolorrgb(#000000)                                         !
      call moveto(40,189,xy)                                            !
      call outgtext(texto)                                              !
      ix=xc(js)/fe/fdib+830                                             !
      if(ix.lt.622.or.ix.gt.1028) go to 30                              !
      iy=220-yc(js)/fe/fdib                                             !
      if(iy.lt. 22.or.iy.gt.418) go to 30                               !
      if(ks.eq.1) nada=setcolorrgb(#f08000)                             ! blue
      if(ks.eq.2) nada=setcolorrgb(#0040f0)                             ! red
      if(ireg(js).eq.1) then                                            ! regional
         if(ks.eq.1) nada=setcolorrgb(#f0d0b0)                          !
         if(ks.eq.2) nada=setcolorrgb(#b0d0f0)                          !
      endif                                                             !
      if(jend.gt.1) nada=setcolorrgb(#000500)                           !
      i=tx(js)/2.2/fe/fdib                                              !
      u1=ix-i                                                           !pixel cell
      if(u1.lt.622) u1=622                                              !
      u2=ix+i                                                           !
      if(u2.gt.1028) u2=1028                                            !
      j=ty(js)/2.2/fe/fdib                                              !
      v1=iy-j                                                           !
      if(v1.lt.22) v1=22                                                !
      v2=iy+j                                                           !
      if(v2.gt.418) v2=418                                              !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !horiz
      iz=-zc(js)/fe/fdib+500                                            !
      k=tz(js)/2.2/fe/fdib                                              !
      v1=iz-k                                                           !
      v2=iz+k                                                           !
      if(iz.lt.640) nada=rectangle($gfillinterior,u1,v1,u2,v2)          !WE
      iy=1265+yc(js)/fe/fdib                                            !
      u1=iy-i                                                           !
      u2=iy+i                                                           !
      if(iz.lt.640) nada=rectangle($gfillinterior,u1,v1,u2,v2)          !SN
      
      u1=lx+ncf/0.7/nc*mx                                               ! dens
      if(u1.gt.570) go to 30                                            !
      v1=ly+40                                                          !
      if(ks.eq.1) call moveto(u1,v1+3,xy)                               !
      if(ks.eq.2) call moveto(u1,v1-3,xy)                               !
      v2=v1+((2*ks-3)*sts)*my                                           !
      nada=lineto(u1,v2)                                                !
         u2=u1+1                                                        !
         nada=setcolorrgb(#f08000)                                      !cdm neg azul
         v1=ly+260+(-zdm(1))/fe/fdib                                 !
         if(zdm(1).ne.0.) nada=rectangle($gfillinterior,u1,v1,u2,v1+1)  !
         nada=setcolorrgb(#0040f0)                                      !cdm npos rojo
         v1=ly+260+(-zdm(2))/fe/fdib                                 !
         if(zdm(2).ne.0.) nada=rectangle($gfillinterior,u1,v1,u2,v1+1)  !
         
C                     j=500+(izm/fe-techo)/fdib
         
         
         if(fitg.eq.0) fitg=fit                                         !
         nada=setcolorrgb(#000000)                                      ! fit negro
         v1=ly+150-fit/fitg*200.                                        !
         if(u1.lt.600) nada=rectangle($gfillinterior,u1,v1,u2,v1+1)     !
         
         nada=setcolorrgb(#d000d0)                                      ! offset line rosa
         v1=ly+160-(dgs-gmi)/fe/250.                                    !
         if(u1.lt.600) nada=rectangle($gfillinterior,u1,v1,u2,v1+1)     !
         nada=setcolorrgb(#00d0d0)                                      ! offset line amar
         v1=ly+160-(pxs-gxmi)/10.                                       !
         if(u1.lt.600) nada=rectangle($gfillinterior,u1,v1,u2,v1+1)     !
         nada=setcolorrgb(#00c000)                                      ! offset line verde
         v1=ly+160-(pys-gymi)/10.                                       !
         if(u1.lt.600) nada=rectangle($gfillinterior,u1,v1,u2,v1+1)     !

!-----------------------------------------------------------------------
      
      go to 30

c----------------------------------------------------------------------!
c              End of model growth                                     !
c----------------------------------------------------------------------!

   50 fs=fsa
      fit0=sqrt(fit0) 
      if(nu.le.300) go to 90
     
    
c----------------------------------------------------------------------!
c           End of process                                             !
c----------------------------------------------------------------------!

      if(npe.gt.0) peme=peme/npe
      if(nes.gt.0) then
        do 51 i=1,nes
        jd(i)=jd(i)*fs 
   51   jz(i)=jca(is(i))
      endif 
      c=0.
      em=0
      do 52 i=1,ns
      gc(i)=gc(i)*fs
      re(i)=g(i)-dgs-(x(i)/1000.*pxs+y(i)/1000.*pys+t(i)*pzs)-gc(i)
      em=em+re(i)
   52 w0(i)=abs(re(i))
      call dmedian(ms,ns,w0,syg)
      syg=syg/0.6745
      call cov2(ms,ns,x,y,re,w,dr,mr,cov,swc)  ! dr/1000.
      covc=(cov(1)+cov(2))/2.
      zb=9.d9
      do 53 j=1,nc
      if(df(j).le.0.or.ireg(j).eq.1) go to 53 
      if(zc(j).lt.zb) zb=zc(j)-tz(j)/1.5
   53 continue 
      ddm=-fam*exp(zb/fe/ka)+ex0
      
c   Writing resulting model 

      open(1,file=mod)
      kk=0
      write(1,'(4i7)') nc,nint(siz),nint(dr)
      js=0
      ks=0
      z1=9.d9
      z2=-z1
      srx=0.
      sry=0.
      fxx=0.
      fyy=0.
      fxz=0.
      f1g=0
      fyz=0.
      fxg=9.d9
      fyg=-fxg
      gxz=0.
      gyz=0.
      px=0.
      py=0.
      s1r=0.
      srr=0.
      srg=0.
      smr=0
      ncc=0
      do 59 j=1,nc
      df(j)=df(j)*fs 
      dd=0
      if(fam.gt.0.01) dd=-fam*exp(zc(j)/fe/ka)+ex0 
      k=-sqrt(dec(j))*avm*1.e4                               ! m.c. atraction ugal for density=10 Tm/m3
      if(ireg(j).eq.0.and.dd.gt.ddm) dd=ddm
      dd=dd*fs
      if(ireg(j).eq.1) dd=0
      if(nes.eq.0) go to 56
      dd=0.
      do 55 i=1,nes
   55 if(zc(j).lt.jz(i)) dd=dd+jd(i)
   56 dec(j)=df(j)+dd
      
      r=(tx(j)+ty(j)+tz(j))/2.9
      if(ireg(j).eq.1.and.df(j).ne.0.) smr=smr+vol(j)        ! tamaño regional
      if(ireg(j).eq.1) go to 57                              !  regional
      k=-k
      if(df(j).ne.0) then
        ncc=ncc+1
        if(zt.gt.z2) z2=zt
        if(zb.lt.z1) z1=zb
        f1g=f1g+r                                           ! mean side
        if(r.lt.fxg) fxg=r                                  ! min. side
        if(r0.gt.fyg) fyg=r                                 ! max. side
        s1r=s1r+vol(j)                                      ! volume
        dg=vol(j)*df(j)                                     ! mass
        f1z=abs(dg)                                         ! absol mass
        srg=srg+f1z                                         ! addition absol. mass
        if(zc(j).gt.0.) srr=srr+f1z                         ! addition terrain mass
      endif
      if(df(j).lt.0) then
        js=js+1
        srx=srx+f1z                                         ! addition neg. mass
        fxx=fxx+vol(j)                                      ! addition neg. volume
        px=px+f1z*zc(j)                                     ! mean depth
        gxz=gxz+zc(j)*f1z*vol(j)                            ! weighted depth
        fxz=fxz+vol(j)*f1z                                  ! addition neg. weights
      endif
      if(df(j).gt.0) then
        ks=ks+1
        sry=sry+f1z                                         ! addition pos. mass
        fyy=fyy+vol(j)                                      ! addition posit. volumw
        py=py+f1z*zc(j)                                     ! mean depth
        gyz=gyz+f1z*zc(j)*vol(j)                            ! weighted depth
        fyz=fyz+vol(j)*f1z                                  ! addition pos. weights
      endif

   57 write(1,'(i7,2i8,3i7,3i7,f8.1)') 
     -nint(xc(j)+ixm),nint(yc(j)+iym),nint(zc(j)+izm),
     -nint(tx(j)),nint(ty(j)),nint(tz(j)),
     -nint(df(j)),nint(dec(j)),k,far(j)
   59 continue
      
      ncf=ks+js                                            ! pos and neg cells
      f1g=f1g/ncf                                          ! mean side
      if(js.gt.0) mnd=srx/fxx                              ! mean neg. dens.
      if(ks.gt.0) mpd=sry/fyy                              ! mean pos. dens.
      mde=izm+(px+py)/srg                                  ! mean depth
      rzz=izm+(gxz+gyz)/(fxz+fyz)                          ! weighted mean depth
      if(srx.gt.0.) px=px/srx+izm                          ! weighted neg. depth
      if(sry.gt.0.) py=py/sry+izm                          ! weighted pos. depth
      if(fxz.gt.0.) gxz=gxz/fxz+izm
      if(fyz.gt.0.) gyz=gyz/fyz+izm

      j=99999
      k=-j

c  Writing modelled and residual values

      open(2,file=obs)
      open(3,file=fil)
      mm=0
      nms=0
      dtg=0.
      j1=9999
      j2=-j1
      gm0=gm0+dgs
      gx0=gx0+pxs
      gy0=gy0+pys
      disr=0
      disl=0
   61 if(ict.eq.0) read(2,*,end=63) xp,yp,zp,dg,eg
      if(ict.eq.1) read(2,*,end=63) xp,yp,zp,dg,eg,terco
      if(xp.eq.0.and.yp.eq.0.and.zp.eq.0) go to 63
      xp=xp-ixm
      yp=yp-iym
      greg=gm0+xp/1000.*gx0+yp/1000.*gy0+pzs*(41.9214*zp-terco)
      zp=zp-izm
      gloc=dg-greg
      dtg=dtg+gloc*gloc
      nms=nms+1
      gcal=0
      do 62 j=1,nc
      xd=xp-xc(j)
      yd=yp-yc(j)
      zd=zp-zc(j)                        
      d=xd*xd+yd*yd+zd*zd+10.*fe 
      av2=zd*zd/d/d/d
      av=sqrt(av2)*6.672e-3*vol(j)
   62 gcal=gcal+av*df(j)
      gres=gloc-gcal
      disl=disl+gloc*gloc
      disr=disr+(greg-gm0)**2
      gxg=abs(gres)/syg
      if(gxg.gt.blun) mm=mm+1
      ts1=bl
      ts2=bl
      ts3=bl
      ts4=bl
      ts5=bl 
      ts6=bl
      if(gxg.gt.1.0*blun) ts1=as
      if(gxg.gt.1.5*blun) ts2=as
      if(gxg.gt.2.0*blun) ts3=as
      if(gxg.gt.2.5*blun) ts4=as      
      if(gxg.gt.3.0*blun) ts5=as      
      if(gxg.gt.4.0*blun) ts6=as            
      write(3,203) xp+ixm,yp+iym,greg,gloc,gcal,gres,gxg,
     -  ts1,ts2,ts3,ts4,ts5,ts6
  203 format(f8.0,5f9.0,f9.2,2x,6a1)
      if(gres.gt.j2) j2=gres
      if(gres.lt.j1) j1=gres
      go to 61
   63 close(3)
      close(2)
      dtg=sqrt(dtg/nms)
      disl=sqrt(disl/nms)
      disr=sqrt(disr/nms)
          
c  Writing some characteristic results in MOD.txt

      dgs=0.
      if(pys.ne.0) dgs= atan2(pxs,pys)*180./3.141592
      if(dgs.lt.0.) dgs=dgs+360.
      if(dgs.gt.180.) dgs=dgs-180.
      sne=0
      spo=0
      if(nes.eq.0) go to 65
      snne=0
      snpo=0
      do 64 i=1,nes
      snne=snne+nwe(i,1)
      snpo=snpo+nwe(i,2)
      sne=sne+wme(i,1)
   64 spo=spo+wme(i,2)
      sne=sne/snne
      spo=spo/snpo 
   65 down=fam*exp(0./fe/ka)-fam*exp(-1000./fe/ka)
      mde=srg/s1r/2.
      write(1,204) -srx/fxx,sry/fyy,srg/s1r,de(1),de(2),
     - px,py,zme+izm,z2+izm,z1+izm,
     - gxz,gyz,rzz,
     - js,ks,ncf,
     - srx,sry,srg,srr,
     - nint(fxg),nint(fyg),nint(f1g),siz,
     - nms,fe,blun,mm,mm*100./ns,
     - disg,disr,disl,ems,
     - fla,uw,fs1,fs,kp,na0,
     - gm0,ig0s,gx0,igxy,gy0,igxy,
     - sqrt(pxs*pxs+pys*pys),dgs,ixm,iym,izm,
     - pzs*1000,ict,kwl,
     - nint(em/ns),j1,j2,ems,syg,fit0,av*100,
     - covc,nint(dr), dtg,ems/dtg,
     - sla,fam,nint(down*fs),
     - jla,jden
  204 format(25('-'),'(-)-------(+)-----(anom)-----(min)-----(max)-----'
     -/'Mean Dens.(kg/m3) :',f10.0,f9.0,5x,f6.0,3x,f6.0,4x,f6.0  
     -/'Depth of cells(m) :',5(f10.0)/'Weighted depths   :',3(f10.0)
     -/'Number of cells   :',3(4x,i6),
     -/'Anomal. mass (kg) :',1x,3e10.3,' (rel.',e10.3,')'
     -/'Side of cells (m) :'
     -/  4x,'Min=',i5,6x,'Max=',i5,6x,'Average=',i5,6x,'Ref='f6.0   
     -/'Num.of stations=',i4,4x,'Scale fact. =',f9.3
     -/  4x,'Outliers(>',f3.1,')=',i3,' (',f4.1,'%)'
     -/'Gravity anomaly dispersion (uGal)'
     -/  4x,'Obs.=',f7.0,3x,'Reg.=',f7.0,3x,'Loc.=',f7.0,3x,'Res.=',f7.0
     -/'Modelling parameters'
     -/  4x,'Lambda fac=',f7.1,6x,'Up weighting=',f4.1
     -/  4x,'Scale fact: in=',f8.0,' fin=',f6.2
     -/  4x,'A priori data weighting Y(1)/N(0):',i1
     -/  4x,'Random exploration coef.:',i3,'   (1 sistem.)'
     -/'Trend =',f9.0,'(',i1,') uGal + (',
     -      2(f7.0,'(',i1,')'),' ) uGal/km'  
     -/4x,'Slope=',f6.0,' uGal/km    Azim.= N',f5.0,' E   (1:adj,0:fix)'
     -/  4x,'Incremental coord. from (xm,ym,zm)=',2i9,i6,' m'
     -/'Additional terrain density=',f6.0,' (',i1,') kg/m3'
     -/  4x,'Terrain filtering radius=',i6,' m' 
     -/'Residues (uGal):'/4x,'Mean=',i4,4x,'Min.=',i6,3x,'Max.=',i5
     -/  4x,'Weight. StDv=',f7.1,4x,'Consist. Stdv=',f7.1,
     -/  4x,'Fit.=',f8.0,' uGal'
     -/  4x,'Mean Sensitivity for cells=',f6.1,' uGal for 100 kg/m3'
     -/  4x,'Autocorrelation =',f6.2,5x,'Correlat. step=',i6,' m' 
     -/  4x,'StDv Data =',f7.0,' uGal     SD Resid/SD Data =',f5.2
     -/'Dens. pattern:'
     -/  4x,'Lateral smoothing (1-10)=',f5.2
     -/  4x,'Flattening=',f4.1,6x,'Downward Incr.:',i4,' kg/m3 each km'
     -/  4x,'Sub-horiz.discontin. (Yes(1)/No(0)):',i2
     -/  4x,'Density contrast fit (Yes(1)/No(0)):',i2)
      if(jden.eq.0) write(1,205) de(1),de(2)
  205 format(4x,'A priori max.densit.contrast=',2f6.0,' kg/m3')
      if(jla.eq.1) write(1,206) ra,bot+izm,
     - nint(q1n),nint(q2n),nint(q1p),nint(q2p),
     - nes,(((jz(i)+izm)/10)*10,i=1,9),
     - (jd(i),i=1,9),(apla(i),i=1,9),
     - (wme(i,1)/nwe(i,1),i=1,9),(wme(i,2)/nwe(i,2),i=1,9),
     - (nwe(i,1)+nwe(i,2),i=1,9), 
     - sne,spo,peme
  206 format(    
     - / 9x,'Layers model:   Shape=',f4.2,'  Bottom=',f8.0
     - / 9x,'Weight.model:  q1n=',i3,' q2n=',i3,'  q1p=',i3,' q2p=',i3
     - / 9x,'Num. of layers=',i2,
     - / 9x,'Depths(m)',9i7/9x,'Density  ',9i7/9x,'Flattening ',9f7.1
     - / 9x,'Neg.weig:',9f7.2/  9x,'Pos.weig:',9f7.2
     - / 9x,'Num.cells',9i7
     - / 9x,'Aver.weights:   Neg=',f5.2,'  Pos=',f5.2,'  Tot=',f5.2 )
      write(1,'(a,a9,a)') '-------------- Date:',hoy,'-----------------'
      close(1)
      
c     Calculating mass instability
      
      tinp=0
      tinn=0
      sum=9000.*fe
      sum=sum*sum
      do 78 i=1,nc
      if(ireg(i).eq.1) go to 78                         ! regional
      if(df(i).lt.0) go to 78
      c=tx(i)/2.1
      do 79 j=1,nc
      if(ireg(j).eq.1) go to 79                         ! regional
      if(df(j).gt.0) go to 79
      if(df(i).eq.0.and.df(j).eq.0) go to 79
      xd=abs(xc(i)-xc(j))
      if(xd.gt.c) go to 79
      yd=abs(yc(i)-yc(j))
      if(yd.gt.c) go to 79
      zd=zc(i)-zc(j)                                   ! i arriba   j abajo
      if(zd.le.1) go to 79
      if(zd.gt.1.4*tz(i)) go to 79
      if((dec(i)-dec(j)).lt.mde/2.) go to 79
      if(df(i).gt.0) tinp=tinp+c*c
      if(df(j).lt.0) tinn=tinn+c*c
c-----------------------------------------------------------------------
c      ix=xc(j)/fe/fdib+830                                              !
c      iy=220-yc(j)/fe/fdib                                              !
c      if(df(i).gt.0) nada=setcolorrgb(#000090)                          !
c      if(df(j).lt.0) nada=setcolorrgb(#b00000)                          !
c      u1=ix-1                                                           !pixel cell
c      u2=ix+2                                                           !
c      v1=iy-1                                                           !
c      v2=iy+2                                                           !
c      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !horiz
c      iz=-zc(j)/fe/fdib+500                                             !
c      v1=iz-1                                                           !
c      v2=iz+2                                                           !
c      if(iz.lt.640) nada=ellipse($gfillinterior,u1,v1,u2,v2)            !WE
c      iy=1275+yc(j)/fe/fdib                                             !
c      u1=iy-1                                                           !
c      u2=iy+2                                                           !
c      if(iz.lt.640) nada=ellipse($gfillinterior,u1,v1,u2,v2)            !SN
c-----------------------------------------------------------------------
   79 continue
   78 continue       
      tinp=100.*tinp/sum  
      tinn=100.*tinn/sum  
      call str(fe,mc,nc,xc,yc,zc,tx,ty,tz,df,ireg,fdib,pen,pep,1)
 
c----------------------- Final drawing ----------------------------------
      ix=1110                                                           !
      iy=810                                                            !
      call ventana(ix-50,iy-105,ix+360,iy+50)                           !
      nada=setcolorrgb(#000000)                                         !
      call moveto(ix,iy,xy)                                             !
      nada=lineto(ix+350,iy)                                            !
      call moveto(ix,iy+40,xy)                                          !
      nada=lineto(ix,iy-80)                                             !
      write(texto,'(a,f6.0)') 'Autocorrel. of Residues    Step(m)=',dr  !
      call moveto(ix+65,iy-100,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,f5.2)') 'Autocorrel =',covc                       !
      call moveto(ix+226,iy-80,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,f9.0)') 'Horiz.Dist.'                             !
      call moveto(ix+280,iy+15,xy)                                      !
      call outgtext(texto)                                              !
      do 84 i=1,7                                                       !
      v1=iy-80*(i-3)*0.25                                               !
      call moveto(ix-5,v1,xy)                                           !
      nada=lineto(ix,v1)                                                !
      write(texto,'(f5.2)') -0.75+i*0.25                                !
      call moveto(ix-42,v1-10,xy)                                       !
      call outgtext(texto)                                              !
      u1=ix+i*tic/fe/fdib                                               !
      call moveto(u1,iy,xy)                                             !
      nada=lineto(u1,iy+6)                                              !
   84 continue                                                          !
      u1=ix+dr*0.5/fe/fdib                                              !
      v1=iy-cov(1)*50                                                   !
      nada=setcolorrgb(#009000)                                         !
      call moveto(u1,v1,xy)                                             !
      do 83 i=1,mr                                                      !
      u1=ix+(i-0.5)*dr/fe/fdib                                          !
      v1=iy-cov(i)*80                                                   !
   83 nada=ellipse($gfillinterior,u1-3,v1-3,u1+3,v1+3)                  !

      call moveto(1350,375,xy)                                          !
      nada=lineto(1400,375)                                             !
      call moveto(1360,370,xy)                                          !
      nada=lineto(1360,375)                                             !
      u1=1360+tic/fe/fdib                                               !
      call moveto(u1,370,xy)                                            !
      nada=lineto(u1,375)                                               !

      nada=setcolorrgb(#f8f8f8)                                         !
      nada=rectangle($gfillinterior,130,708,290,750)                    !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(2f14.2,a)') -pen,pep,'  ......  Mass/depth slope'   !
      call moveto(82,710,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(2f14.1,a,6x,a)') tinn,tinp,                         !
     & '  ......  Mass instability %','- - - - - - ->'                  !
      call moveto(82,730,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)') '>>> Completed process <<<'                    !
      call moveto(290,800,xy)                                           !
      call outgtext(texto)                                              !
      write(texto,'(a)') ' ===> Press (Enter)'                          !
      call moveto(300,840,xy)                                           !
      call outgtext(texto)                                              !
      
      ix=720                                                            !
      write(texto,'(a)') ' Some suggestions:'                           !
      call moveto(680,ix+20,xy)                                         !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Flattening:'                                  !
      call moveto(840,ix-10,xy)                                         !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Up-Weight :'                                  !
      call moveto(840,ix+10,xy)                                         !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Lateral Co:'                                  !
      call moveto(840,ix+30,xy)                                         !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Random Coe:'                                  !
      call moveto(840,ix+50,xy)                                         !
      call outgtext(texto)                                              !
      write(texto,'(a,13x,a)') 'Smoothness:' ,'<- - -'                  !
      call moveto(840,ix+70,xy)                                         !
      call outgtext(texto)                                              !
      if(pep.lt.2.or.pen.lt.2.or.tinn.gt.10.or.tinp.gt.10) then         !
        write(texto,'(a)') 'Increase'                                   !
        call moveto(933,ix-10,xy)                                       !
        call outgtext(texto)                                            !
      endif                                                             !
      if(zdm(1).lt.1.5*zdm(2)) then                                     !
        write(texto,'(a)') 'Increase'                                   !
        call moveto(933,ix+10,xy)                                       !
        call outgtext(texto)                                            !
      endif                                                             !
        if(denmi.gt.1.5) write(texto,'(a)') 'Increase'                  !
        if(denmi.le.1.5) write(texto,'(a)') 'Decrease'                  !
        call moveto(933,ix+30,xy)                                       !
        call outgtext(texto)                                            !
      if(na0.gt.6) then                                                 !
        write(texto,'(a,i4,a)')  'Decrease'                             !
        call moveto(933,ix+50,xy)                                       !
        call outgtext(texto)                                            !
      endif                                                             !
      if(covc.gt.0.05) then                                             !
        write(texto,'(a,i4,a)')  'Decrease'                             !
        call moveto(933,ix+70,xy)                                       !
        call outgtext(texto)                                            !
      endif                                                             !
      if(covc.lt.0.01) then                                             !
        write(texto,'(a,i4,a)')  'Increase'                             !
        call moveto(933,ix+70,xy)                                       !
        call outgtext(texto)                                            !
      endif                                                             !

      ix=1085                                                           !
      iy=810                                                            !
      nada=setcolorrgb(#f8f8f8)                                         !
      nada=rectangle($gfillinterior,ix+120, 80,ix+210, 97)              !
      nada=rectangle($gfillinterior,ix+120,117,ix+390,159)              !
      nada=rectangle($gfillinterior,ix+120,180,ix+390,290)              !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a,f8.0,a)') 'Misfit value.....',fit0,' uGal'        !
      call moveto(ix,80,xy)                                             !
      call outgtext(texto)                                              !
      write(texto,'(a)')  'Tot.        Neg.       Posit.'               !
      call moveto(ix+160,100,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(a,3e12.2)') 'Anom. mass(kg)...',srg,srx,sry         !
      call moveto(ix,120,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a,3f12.0)') 'Aver. depth(m)...',zme+izm,px,py       !
      call moveto(ix,140,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Tot.   Reg.   Loc.  Resid.  Res/Loc'          !
      call moveto(ix+145,160,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(a,2x,4f7.0,f7.2)') 'Dispersion(uGal).',             !
     & disg,disr,disl,ems,  ems/dtg                                     !
      call moveto(ix,180,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a,1x,f8.2)') 'Average weight...',peme               !
      call moveto(ix,200,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a,3x,e9.2)') 'External cells...',smr*fs             !
      call moveto(ix,220,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a,i9,a,i5,a,i5,a)') 'Regional trend...',nint(gm0),  !
     & ' uGal + (',nint(gx0),',',nint(gy0),') uGal/km '                 !
      call moveto(ix,240,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a,f7.0,a)') 'Addit.terr.dens..',pzs*1000,' kg/m3'   !
      call moveto(ix,260,xy)                                            !
      call outgtext(texto)                                              !

c-----------------------------------------------------------------------!
      go to 99

C  Some messages for error or warning

   90 nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)')  '*** Error: Unexpected execution end !!!'
      call moveto(650,720,xy)
      call outgtext(texto)
      if(nu.le.1) write(texto,'(a)')
     &'Verify data suitability (units, magnitudes) and model hypothesis'
      call moveto(660,740,xy)
      call outgtext(texto)

      go to 99
   92 write(texto,'(3a)') '*** Error reading file ',obs,' !!'
      call moveto(650,720,xy)
      call outgtext(texto)
      go to 99
   99 read(*,'(i1)') i
      stop
      end

      
c***********************************************************************
c  attrac(+uGal) of a right prism, parallel to axis, upon coord origin *
c       xc,yc: coordinates (meters) of the prism center                *
c       dx,dy: sides (meters) of the prism                             *
c       zb,zt: vertical coordinates of prism bottom and top            *
c       d: density (kg/m3)        av: vertical attraction (uGal) (+z)  *
c***********************************************************************
      subroutine VAP(xc,yc,dx,dy,zb,zt,d,av)
      dimension x(4),y(4),s(4),z(2)
         zc=(zt+zb)/2.                                  ! spherical aprox.
         a=xc*xc+yc*yc+zc*zc
              r=(dx+dy+zb-zt)/2.9
              r=r*r 
              if(a.lt.r) a=r                             ! ojo!!!!     
         av=dx*dy*(zb-zt)*zc/sqrt(a*a*a) *d*6.672e-3
         if(zc.ne.0.) return
      a=dx/2.
      b=dy/2.
      x(1)=xc+a
      x(2)=xc-a
      x(3)=x(1)
      x(4)=x(2)
      y(1)=yc-b
      y(2)=y(1)
      y(3)=yc+b
      y(4)=y(3)
      z(1)=zb
      z(2)=zt
      do 1 i=1,4
      s(i)=0.
      if(x(i).eq.0.or.y(i).eq.0.) go to 1
      call PIC(x(i),y(i),z,s(i))
    1 continue
      av=(s(3)+s(2)-s(1)-s(4))*d*6.672e-3
      return

      end
      subroutine PIC(xa,ya,z,av)
      dimension z(2)
      x= abs(xa)
      y= abs(ya)
      xx=x*x
      yy=y*y
      av=0.
      i=1
      do 1 k=1,2
      i=-i
      zz=z(k)*z(k)
      r= sqrt(xx+yy+zz)
      rx= sqrt(xx+zz)
      ry= sqrt(yy+zz)
      if(z(k).eq.0.) go to 1
      av=av+(z(k)* atan(z(k)*r/x/y)-abs(z(k))*1.570796327)*i
    1 av=av+(x* log((y+r)/rx)+y* log((x+r)/ry))*i
      if(xa.lt.0.) av=-av
      if(ya.lt.0.) av=-av
      return
      end

c***********************************************************************
c     Autocorrelation of n<m data v(x,y) for a distance dr             c
c***********************************************************************
      subroutine cov2(m,n,x,y,v,w,dr,mr,cov,sw)
      dimension x(m),y(m),v(m),w(m),cov(mr),sw(mr)
        
      vm=0.
      suw=0.
      do 1 i=1,n
      suw=suw+w(i)
    1 vm=vm+v(i)*w(i)
      vm=vm/suw
      sd2=0
      do 2 i=1,n
      vi=v(i)-vm
    2 sd2=sd2+vi*vi*w(i)
      sd2=sd2/suw
      
      do 7 k=1,mr
      cov(k)=0.
    7 sw(k)=0. 
      dr2=dr*dr
     
      do 3 i=1,n-1
      do 4 j=i+1,n
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      d=xx*xx+yy*yy
      do 5 k=1,mr
      if(d.ge.dr2*(k-1).and.d.le.dr2*k) then
        cov(k)=cov(k)+(v(j)-vm)*(v(i)-vm)*w(j)*w(i)
        sw(k)=sw(k)   +w(j)*w(i)
      endif
    5 continue  
    4 continue
    3 continue
      do 6 k=1,mr
    6 if(sw(k).ne.0) cov(k)=cov(k)/sw(k)/sd2
      
      return
      end
      
c***********************************************************************
c   median of absolute values for n<m different numbers x(i)           c
c***********************************************************************
      subroutine dmedian(m,n,x,xa)
      dimension x(m)
      xi=x(1)
      xs=xi
      do 1 i=1,n
      xx=x(i)
      if(xx.lt.xi)  xi=xx
      if(xx.gt.xs)  xs=xx
        do 6 j=i+1,n
    6   if(xx.eq.x(j)) x(j)=x(j)*1.00001
    1 continue
      j=0
    4 xa=xs
      do 2 i=1,n
      xx=x(i)
      if(xx.lt.xa.and.xx.gt.xi) xa=xx
    2 continue
      j=j+1
      if(j.eq.n) go to 9
      xi=xa
      xb=xi
      do 3 i=1,n
      xx=x(i)
      if(xx.gt.xb.and.xx.lt.xs) xb=xx
    3 continue
      j=j+1
      xa=(xa+xb)/2.
      xs=xb
      if(j.lt.n) go to 4
    9 return
      end
      
c***********************************************************************
c     Execution stops when number of data is higher than dimensions    *
c***********************************************************************
      subroutine dim(m)
      write(*,200) m
  200 format(/6x,'*** Error: data size > max. dimension',i6,' !!!'/
     -10x,'Reduce the data side or change dimension in the code'/)
      stop
      end
      
c***********************************************************************
c     Sensitivity calculus
c***********************************************************************
      subroutine at2(ms,ns,x,y,z,xx,yy,zz,dx,dy,dz,sen,av,zp)
      dimension x(ms),y(ms),z(ms)
      av=0.
      t=(dx+dy+dz)/3
      t2=t*t
      zp=0
      sp=0 
      do 1 i=1,ns
      zr=zz-z(i)
      xr=xx-x(i)
      yr=yy-y(i)
      d=xr*xr+yr*yr+zr*zr+500.
      if(d.lt.t2) d=t2
      av=av+zr*zr/d/d/d
      d=1/d/d
      sp=sp+d
    1 zp=zp+z(i)*d
      av=av/ns
      sen=3.d-7/sqrt(av)
      t2=dx*dy*dz
      av=av*t2*t2  
      zp=zp/sp
      return
      end

c***********************************************************************
c     New cells for a layer
c           zs: top of the layer    tz: width of the layer 
c***********************************************************************
      subroutine celdascapa(ax,bx,ay,by,zs,tz,avm,ms,ns,x,y,z,
     -mc,nc,xc,yc,zc,dx,dy,dz,at,seli,ncr,ireg)
      dimension xc(mc),yc(mc),zc(mc),dx(mc),dy(mc),dz(mc),
     -x(ms),y(ms),z(ms),ireg(mc)
      character*50 texto
      real*4 at(mc)
      zi=zs-tz                           
      zz=zs-tz/2.
      tx=tz                               ! easting size
      ty1=tz/3                            ! minimo atraccion  /4
      ty2=tz*3                            ! maximo atraccion  *4
      tt=tz*tz 
      kk=nc
      
      do 7 i=1,999                       
      xx=ax+tx*(i-1)
      if(xx.gt.bx) go to 9
      ty=tz
      yr=ay
      do 6 j=1,999                        
      if(yr.ge.by) go to 7
      c=ty
      l=0
    3 yy=yr+c/2
      call at2(ms,ns,x,y,z,xx,yy,zz,tx,c,tz,sen,am,zp)
      if(sen.gt.seli*10.) go to 6             
      if((zp-100).lt.zs.or.am.le.0) go to 6   
      av=am/avm
      if(abs(av-1.).le.0.1) go to 4         ! 0.2 
      c=c/av**0.41
      if(c.gt.ty2) go to 6            !si tamaño en Y supera limites no vale
      if(c.lt.ty1) go to 6            !si es muy en Y pequeño no vale
      l=l+1
      if(l.le.10) go to 3
      go to 6
    4 if(c.gt.9999) go to 6
    5 ty=c                               ! northing size

c   Eliminacion celdas superpuestas

      if(ncr.eq.0) go to 12
      do 11 k=1,ncr
      if(abs(xx-xc(k)).gt.(dx(k)+tx)/2.6) go to 11     ! 2.6
      if(abs(yy-yc(k)).gt.(dy(k)+ty)/2.6) go to 11     ! 2.6
      if(abs(zz-zc(k)).gt.(dz(k)+tz)/2.6) go to 11     ! 2.6
      go to 6
   11 continue
   12 kk=kk+1
      xc(kk)=xx
      yc(kk)=yy
      zc(kk)=zz
      dx(kk)=tx
      dy(kk)=ty
      dz(kk)=tz
      at(kk)=av
      ireg(kk)=0
      if(ncr.gt.1) ireg(kk)=1
      if(kk.lt.mc) go to 6 
         kk=-9999
         go to 9
    6 yr=yr+ty
    7 continue
      
    9 nc=kk
      return
      end
      
c***********************************************************************
c   Width of a layer
c***********************************************************************
      subroutine especapa(zs,avm,ms,ns,x,y,z,tz)
      dimension x(ms),y(ms),z(ms)
      l=0
      c=tz                           
    1 zz=zs-c*0.5                    
      k=0
      ama=-9.d9
      ame=0
      do 2 i=1,ns
      call at2(ms,ns,x,y,z,x(i)*1.5,y(i)*1.5,zz,c,c,c,se,aa,zp)
      if((zp-zz).le.100) go to 2
      if(aa.gt.ama) ama=aa
      ame=ame+aa
      k=k+1
    2 continue
      if(k.le.1) go to 9
      am=ame/k/avm
      c=c/am**0.13
      l=l+1
      if(abs(am-1.).gt.0.1.and.l.lt.12) go to 1
      
      if(c.gt.tz) tz=c
    9 return
      end
      
c***********************************************************************
c     Step for autocorrelation                                         c
c***********************************************************************
      subroutine step(m,n,x,y,d,dr)
      dimension x(m),y(m),d(m)
      do 1 i=1,n
      d4=9.e9
      d3=9.e9
      d2=9.e9
      de=9.e9
      do 2 j=1,n
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      dd= xx*xx+yy*yy
      if(dd.gt.d4.or.dd.le.1.) go to 2
      if(dd.lt.de) de=dd
      if(dd.lt.d2.and.dd.gt.de) d2=dd
      if(dd.lt.d3.and.dd.gt.d2) d3=dd
      if(dd.lt.d4.and.dd.gt.d3) d4=dd
    2 continue
      de=sqrt(de)
      d2=sqrt(d2)
      d3=sqrt(d3)
      d4=sqrt(d4)
    1 d(i)=(de+d2+d3+d4)/4.
      call dmedian(m,n,d,dr)
      return
      end
      
c****************************************************************
c     Square window                                             *
c****************************************************************
      subroutine ventana(i1,j1,i2,j2)
      use IFQWIN
      integer*2 x1,x2,y1,y2
      nada=setcolorrgb(#909090)
      x1=i1-2
      y1=j1-2
      x2=i1-1
      y2=j2+1
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      x2=i2+1
      y2=j1-1
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      nada=setcolorrgb(#fefefe)
      x1=i1-1
      y1=j2+2 
      x2=i2+2
      y2=j2+2
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      y1=j1
      nada=rectangle($gfillinterior,x2,y1,x2,y2)
      nada=setcolorrgb(#ffffff)
      x1=i1
      x2=i2
      y2=j2
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      return
      end
      
c**********************************************************************
c        Rounding                                                     *
c**********************************************************************    
      subroutine round(c,r)
      j=1
      if(c.lt.0) j=-1
      r=abs(c)
      i=log10(r)                                               
      r=j*nint(r/10**i)*10**i  
      return 
      end
      
c***********************************************************************
c     Dialog interface 1 (data files)                                  *
c***********************************************************************
      subroutine dial1(obs,mod,fil,kp,na,fam,sla,fla,uw)
      use dflogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto,obs,mod,fil
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Growth, dlg)

      retlog = DlgSet(dlg,IDC_EDIT_e1 ,"Gra.dat" )
      retlog = DlgSet(dlg,IDC_EDIT_e3 ,"Mod.dat" )
      retlog = DlgSet(dlg,IDC_EDIT_e4 ,"Fil.dat" )

      retlog = DlgSet(dlg,IDC_EDIT_e12,"-100" )
      retlog = DlgSet(dlg,IDC_EDIT_e13,"100")
      retlog = DlgSet(dlg,IDC_EDIT_e8,"0.")
      retlog = DlgSet(dlg,IDC_EDIT_e9,"0")
      retlog = DlgSet(dlg,IDC_EDIT_e10,"0")      
      write(texto,'(i3)') na            
      retlog = DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i3)') nint(fla)      
      retlog = DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(f3.1)') sla
      retlog = DlgSet(dlg,IDC_EDIT_e11,texto)      
      write(texto,'(f3.1)') fam      
      retlog = DlgSet(dlg,IDC_EDIT_e14,texto) 
      write(texto,'(f3.1)') uw     
      retlog = DlgSet(dlg,IDC_EDIT_e15,texto) 
      retlog = DlgSet(dlg,IDC_CHECK_c2,.false.)
      retlog = DlgSet(dlg,IDC_CHECK_c4,.true.)
      retlog = DlgSet(dlg,IDC_CHECK_c5,.false.)
      retlog = DlgSet(dlg,IDC_CHECK_c3,.true.)      
      retlog = DlgSet(dlg,IDC_CHECK_c1,.false.)
      retlog = DlgSet(dlg,IDC_CHECK_c6,.true.)

      retlog=DlgSet(dlg,IDC_TEXT_t5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t7 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t8 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t9 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t10 ,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_TEXT_t11 ,.false.,dlg_enable) 
      retlog=DlgSet(dlg,IDC_TEXT_t12,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t13,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t14,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_EDIT_e5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e7 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e8 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e9 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e10,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e11,.false.,dlg_enable)            
      retlog=DlgSet(dlg,IDC_EDIT_e12,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_EDIT_e14,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e15,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_EDIT_e32,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e33,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e34,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e12,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e13,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_BOX_b4  ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_BOX_b5  ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_BOX_b6  ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_BOX_b7  ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_CHECK_c2,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_CHECK_c3,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_CHECK_c4,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_CHECK_c5,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_CHECK_c1,.true.,dlg_enable)
      retlog=DlgSet(dlg,IDC_CHECK_c6,.false.,dlg_enable)
      retlog = DlgSetSub(dlg,IDCANCEL,fin)

      retint = DlgModal( dlg )

      nada=DlgGet(dlg,IDC_EDIT_e1, obs )
      nada=DlgGet(dlg,IDC_EDIT_e3, mod )
      nada=DlgGet(dlg,IDC_EDIT_e4, fil )
      kp=0
      retlog=DlgGet( dlg, IDC_CHECK_c1, nada)
      if(nada) kp=1

      return
      end      
           
      
c***********************************************************************
c     Dialog interface 2 (geometrical parameters)                      *
c***********************************************************************
      subroutine dial2(ns,izm,fe,techo,base,pp)
      use dialogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto
      external fin
      type(dialog) dlg
      retlog=DlgInit(IDD_Grid3D, dlg)
      retlog = DlgSet(dlg,IDC_TEXT_t4,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e4,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_CHECK_c1,.false.)
      write(texto,'(i7)') nint(techo*fe+izm)
      nada=DlgSet(dlg,IDC_EDIT_e1,texto)
      write(texto,'(i7)') nint(base*fe+izm)
      nada=DlgSet(dlg,IDC_EDIT_e2,texto)
      write(texto,'(i6)') nint(pp*fe)
      nada=DlgSet(dlg,IDC_EDIT_e3,texto)
      write(texto,'(i6)') ns
      nada=DlgSet(dlg,IDC_EDIT_e5,texto)
      retlog=DlgSetSub(dlg,IDCANCEL,fin)
      retint = DlgModal( dlg )
      nada=DlgGet(dlg,IDC_EDIT_e1,texto)
      read(texto,*) k
      techo=(k-izm)/fe
      nada=DlgGet(dlg,IDC_EDIT_e2,texto)
      read(texto,*) k
      base=(k-izm)/fe
      nada=DlgGet(dlg,IDC_EDIT_e3,texto)
      read(texto,*) k
      pp=k/fe

      return
      end
c***********************************************************************
c     Dialog interface 3 ( continue / re-run )                         *
c***********************************************************************
      subroutine dial3(ns,izm,fe,techo,base,pp,nc,ire)
      use dialogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto
      external fin
      type(dialog) dlg
              ! supress compiler warnings for unreferenced arguments
      if(.not.DlgInit(IDD_Grid3D,dlg)) then
      write(*,*) "error: resource not found,  vale?"
      else
          
      retlog=DlgSet(dlg,IDC_CHECK_c1,.false.)
      retlog=DlgSetSub(dlg,IDCANCEL,fin)
      write(texto,'(i7)') nint(techo*fe+izm)
      nada=DlgSet(dlg,IDC_EDIT_e1,texto)
      write(texto,'(i7)') nint(base*fe+izm)
      nada=DlgSet(dlg,IDC_EDIT_e2,texto)
      write(texto,'(i6)') nint(pp*fe)
      nada=DlgSet(dlg,IDC_EDIT_e3,texto)   
      write(texto,'(i6)') nc
      nada=DlgSet(dlg,IDC_EDIT_e4,texto) 
      write(texto,'(i6)') ns
      nada=DlgSet(dlg,IDC_EDIT_e5,texto)

      ire=0
      retint = DlgModal( dlg )
      retlog=DlgGet(dlg,IDC_CHECK_c1,nada )
      if(nada) ire=1
      end if
      return
      end
      
c***********************************************************************
c     Dialog interface 4  (inversion parameters)                       *
c***********************************************************************
      subroutine dial4(obs,mod,fil,kp,na,kwl,nc,ns,fla,
     & de,ig0,gm0,igxy,gx0,gy0,ict,jla,sla,jden,fam,uw)
      use dflogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto,obs,mod,fil
      dimension de(2)
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Growth, dlg)
                                                    
      retlog = DlgSet(dlg,IDC_EDIT_e1,obs )
      retlog = DlgSet(dlg,IDC_EDIT_e3,mod )
      retlog = DlgSet(dlg,IDC_EDIT_e4,fil )

      retlog = DlgSet(dlg,IDC_TEXT_t1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t4,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_BOX_b1 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_BOX_b2 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e4,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c3,.true.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c4,.true.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c5,.true.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c6,.true.,dlg_enable)      
      retlog = DlgSet(dlg,IDC_CHECK_c7,.true.,dlg_enable)
      write(texto,'(i6)') ns
      retlog = DlgSet(dlg,IDC_EDIT_e30,texto)
      write(texto,'(i6)') kwl
      retlog = DlgSet(dlg,IDC_EDIT_e7,texto)
      write(texto,'(i6)') nc
      retlog = DlgSet(dlg,IDC_EDIT_e31,texto)
      write(texto,'(i3)') na            
      retlog = DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i3)') nint(fla)      
      retlog = DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(f3.1)') sla
      retlog = DlgSet(dlg,IDC_EDIT_e11,texto)      
      write(texto,'(f3.1)') fam      
      retlog = DlgSet(dlg,IDC_EDIT_e14,texto) 
      write(texto,'(f3.1)') uw     
      retlog = DlgSet(dlg,IDC_EDIT_e15,texto) 
      retlog = DlgSet(dlg,IDC_EDIT_e12,"-100" )
      retlog = DlgSet(dlg,IDC_EDIT_e13,"100")
      write(texto,'(i8)') nint(gm0)
      retlog = DlgSet(dlg,IDC_EDIT_e8,texto)
      retlog = DlgSet(dlg,IDC_EDIT_e9,"0.")
      retlog = DlgSet(dlg,IDC_EDIT_e10,"0.")
      retlog = DlgSet(dlg,IDC_CHECK_c1,.false.)
      if(kp.eq.1) retlog = DlgSet(dlg,IDC_CHECK_c1,.true.)
      retlog = DlgSet(dlg,IDC_CHECK_c2,.false.)      
      retlog = DlgSet(dlg,IDC_CHECK_c3,.true.)
      retlog = DlgSet(dlg,IDC_CHECK_c4,.true.)
      retlog = DlgSet(dlg,IDC_CHECK_c5,.false.)
      retlog = DlgSet(dlg,IDC_CHECK_c6,.true.)
      retlog = DlgSetSub(dlg,IDCANCEL,fin)
      
      retint = DlgModal( dlg )
                                                  
      retlog = DlgGet( dlg, IDC_EDIT_e5,texto)
      read (texto,*) na
      if(na.lt.1) na=1
      retlog = DlgGet( dlg, IDC_EDIT_e6, texto)
      read(texto,*) fla
      retlog = DlgGet( dlg, IDC_EDIT_e7,texto)
      read (texto,*) kwl
      retlog = DlgGet( dlg, IDC_EDIT_e9,texto)
      read(texto,*) gx0
      retlog = DlgGet( dlg, IDC_EDIT_e10,texto)
      read(texto,*) gy0
      retlog = DlgGet( dlg, IDC_EDIT_e12, texto)
      read(texto,*) de(1)
      retlog = DlgGet( dlg, IDC_EDIT_e13, texto)
      read(texto,*) de(2)
      retlog = DlgGet( dlg, IDC_EDIT_e8, texto)
      read(texto,*) gm0
      retlog = DlgGet( dlg, IDC_EDIT_e11, texto)
      read(texto,*) sla
      retlog = DlgGet( dlg, IDC_EDIT_e14, texto)
      read(texto,*) fam
      retlog = DlgGet( dlg, IDC_EDIT_e15, texto)
      read(texto,*) uw 
      ig0=0
      retlog = DlgGet( dlg, IDC_CHECK_c3, nada)
      if(nada) ig0=1
      igxy=0
      retlog = DlgGet( dlg, IDC_CHECK_c4, nada)
      if(nada) igxy=1
      ict=0
      retlog = DlgGet( dlg, IDC_CHECK_c5, nada)
      if(nada) ict=1
      jla=0
      retlog = DlgGet( dlg, IDC_CHECK_c2, nada)
      if(nada) jla=1
      jden=0
      retlog = DlgGet( dlg, IDC_CHECK_c6, nada)
      if(nada) jden=1

      return
      end
c***********************************************************************
c     Dialog interface 5                                               *
c***********************************************************************
      subroutine dial5(nes,ap,jind,bot,ra,q1n,q2n,q1p,q2p)
      use dflogm
      include 'growth.fd'
      logical retlog,nada
      external fin
      character*50 texto
      type(dialog) dlg
      retlog=DlgInit( IDD_Model,dlg)
      write(texto,'(f3.1)') ra
      retlog = DlgSet(dlg,IDC_EDIT1,texto)
      write(texto,'(f8.0)') bot
      retlog = DlgSet(dlg,IDC_EDIT3,texto)
      write(texto,'(i2)') nes
      retlog = DlgSet(dlg,IDC_EDIT4,texto)
      write(texto,'(f4.1)') ap
      retlog = DlgSet(dlg,IDC_EDIT5,texto)
      write(texto,'(f4.0)') q1n
      retlog = DlgSet(dlg,IDC_EDIT6,texto)
      write(texto,'(f4.0)') q2n
      retlog = DlgSet(dlg,IDC_EDIT7,texto)
      write(texto,'(f4.0)') q1p
      retlog = DlgSet(dlg,IDC_EDIT8,texto)
      write(texto,'(f4.0)') q2p
      retlog = DlgSet(dlg,IDC_EDIT9,texto)
      retlog = DlgSet(dlg,IDC_CHECK1,.false.)
      retlog = DlgSetSub(dlg,IDCANCEL,fin)

      retint = DlgModal( dlg )

      retlog = DlgGet( dlg, IDC_EDIT1,texto)
      read (texto,*) ra
      retlog = DlgGet( dlg, IDC_EDIT3,texto)
      read (texto,*) bot
      retlog = DlgGet( dlg, IDC_EDIT4,texto)
      read (texto,*) nes
      retlog = DlgGet( dlg, IDC_EDIT5,texto)
      read (texto,*) ap
      retlog = DlgGet( dlg, IDC_EDIT6,texto)
      read (texto,*) q1n
      retlog = DlgGet( dlg, IDC_EDIT7,texto)
      read (texto,*) q2n
      retlog = DlgGet( dlg, IDC_EDIT8,texto)
      read (texto,*) q1p
      retlog = DlgGet( dlg, IDC_EDIT9,texto)
      read (texto,*) q2p
      jind=0 
      retlog = DlgGet( dlg, IDC_CHECK1, nada)
      if(nada) jind=1

      call DlgUninit(dlg )
      return
      end
c***********************************************************************
c     Dialog interface 6  (a priori depth for layers)                  *
c***********************************************************************
      subroutine dial6(den,mh,nes,jz,jd,apla)
      use dflogm
      include 'growth.fd'
      logical retlog,nada
      external fin
      character*50 texto
      dimension jd(mh),jz(mh),apla(mh)
      type(dialog) dlg
      retlog=DlgInit( IDD_Layers,dlg)
      write(texto,'(i6)') jz(1)
      retlog = DlgSet(dlg,IDC_EDIT1,texto)
      write(texto,'(i6)') jz(2)
      retlog = DlgSet(dlg,IDC_EDIT2,texto)
      write(texto,'(i6)') jz(3)
      retlog = DlgSet(dlg,IDC_EDIT3,texto)
      write(texto,'(i6)') jz(4)
      retlog = DlgSet(dlg,IDC_EDIT4,texto)
      write(texto,'(i6)') jz(5)
      retlog = DlgSet(dlg,IDC_EDIT5,texto)
      write(texto,'(i6)') jz(6)
      retlog = DlgSet(dlg,IDC_EDIT6,texto)
      write(texto,'(i6)') jz(7)
      retlog = DlgSet(dlg,IDC_EDIT20,texto)
      write(texto,'(i6)') jz(8)
      retlog = DlgSet(dlg,IDC_EDIT23,texto)
      write(texto,'(i6)') jz(9)
      retlog = DlgSet(dlg,IDC_EDIT26,texto)

      write(texto,'(i4)') nint(den)
      retlog = DlgSet(dlg,IDC_EDIT7,"0")
      if(nes.ge.1) retlog = DlgSet(dlg,IDC_EDIT7,texto)
      retlog = DlgSet(dlg,IDC_EDIT8,"0")
      if(nes.ge.2) retlog = DlgSet(dlg,IDC_EDIT8,texto)
      retlog = DlgSet(dlg,IDC_EDIT9,"0")
      if(nes.ge.3) retlog = DlgSet(dlg,IDC_EDIT9,texto)
      retlog = DlgSet(dlg,IDC_EDIT10,"0")
      if(nes.ge.4) retlog = DlgSet(dlg,IDC_EDIT10,texto)
      retlog = DlgSet(dlg,IDC_EDIT11,"0")
      if(nes.ge.5) retlog = DlgSet(dlg,IDC_EDIT11,texto)
      retlog = DlgSet(dlg,IDC_EDIT12,"0")
      if(nes.ge.6) retlog = DlgSet(dlg,IDC_EDIT12,texto)
      retlog = DlgSet(dlg,IDC_EDIT21,"0")
      if(nes.ge.7) retlog = DlgSet(dlg,IDC_EDIT21,texto)
      retlog = DlgSet(dlg,IDC_EDIT24,"0")
      if(nes.ge.8) retlog = DlgSet(dlg,IDC_EDIT24,texto)
      retlog = DlgSet(dlg,IDC_EDIT27,"0")
      if(nes.ge.9) retlog = DlgSet(dlg,IDC_EDIT27,texto)

      write(texto,'(f5.1)') apla(1)
      retlog = DlgSet(dlg,IDC_EDIT13,texto)
      write(texto,'(f5.1)') apla(2)
      retlog = DlgSet(dlg,IDC_EDIT14,texto)
      write(texto,'(f5.1)') apla(3)
      retlog = DlgSet(dlg,IDC_EDIT15,texto)
      write(texto,'(f5.1)') apla(4)
      retlog = DlgSet(dlg,IDC_EDIT16,texto)
      write(texto,'(f5.1)') apla(5)
      retlog = DlgSet(dlg,IDC_EDIT17,texto)
      write(texto,'(f5.1)') apla(6)
      retlog = DlgSet(dlg,IDC_EDIT18,texto)
      write(texto,'(f5.1)') apla(7)
      retlog = DlgSet(dlg,IDC_EDIT22,texto)
      write(texto,'(f5.1)') apla(8)
      retlog = DlgSet(dlg,IDC_EDIT25,texto)
      write(texto,'(f5.1)') apla(9)
      retlog = DlgSet(dlg,IDC_EDIT28,texto)

      retlog = DlgSetSub(dlg,IDCANCEL,fin)

      retint = DlgModal( dlg )

      retlog = DlgGet( dlg, IDC_EDIT1,texto)
      read (texto,*) jz(1)
      retlog = DlgGet( dlg, IDC_EDIT2,texto)
      read (texto,*) jz(2)
      retlog = DlgGet( dlg, IDC_EDIT3,texto)
      read (texto,*) jz(3)
      retlog = DlgGet( dlg, IDC_EDIT4,texto)
      read (texto,*) jz(4)
      retlog = DlgGet( dlg, IDC_EDIT5,texto)
      read (texto,*) jz(5)
      retlog = DlgGet( dlg, IDC_EDIT6,texto)
      read (texto,*) jz(6)
      retlog = DlgGet( dlg, IDC_EDIT20,texto)
      read (texto,*) jz(7)
      retlog = DlgGet( dlg, IDC_EDIT23,texto)
      read (texto,*) jz(8)
      retlog = DlgGet( dlg, IDC_EDIT26,texto)
      read (texto,*) jz(9)

      retlog = DlgGet( dlg, IDC_EDIT7,texto)
      read (texto,*) jd(1)
      retlog = DlgGet( dlg, IDC_EDIT8,texto)
      read (texto,*) jd(2)
      retlog = DlgGet( dlg, IDC_EDIT9,texto)
      read (texto,*) jd(3)
      retlog = DlgGet( dlg, IDC_EDIT10,texto)
      read (texto,*) jd(4)
      retlog = DlgGet( dlg, IDC_EDIT11,texto)
      read (texto,*) jd(5)
      retlog = DlgGet( dlg, IDC_EDIT12,texto)
      read (texto,*) jd(6)
      retlog = DlgGet( dlg, IDC_EDIT21,texto)
      read (texto,*) jd(7)
      retlog = DlgGet( dlg, IDC_EDIT24,texto)
      read (texto,*) jd(8)
      retlog = DlgGet( dlg, IDC_EDIT27,texto)
      read (texto,*) jd(9)

      retlog = DlgGet( dlg, IDC_EDIT13,texto)
      read (texto,*) apla(1)
      retlog = DlgGet( dlg, IDC_EDIT14,texto)
      read (texto,*) apla(2)
      retlog = DlgGet( dlg, IDC_EDIT15,texto)
      read (texto,*) apla(3)
      retlog = DlgGet( dlg, IDC_EDIT16,texto)
      read (texto,*) apla(4)
      retlog = DlgGet( dlg, IDC_EDIT17,texto)
      read (texto,*) apla(5)
      retlog = DlgGet( dlg, IDC_EDIT18,texto)
      read (texto,*) apla(6)
      retlog = DlgGet( dlg, IDC_EDIT22,texto)
      read (texto,*) apla(7)
      retlog = DlgGet( dlg, IDC_EDIT25,texto)
      read (texto,*) apla(8)
      retlog = DlgGet( dlg, IDC_EDIT28,texto)
      read (texto,*) apla(9)

      call DlgUninit(dlg )
      return
      end
     
c***********************************************************************
      subroutine fin
      use dialogm
      include 'Growth.fd'
      stop
      end
      
c************************************************************************
c        Distribution:   anomalous mass  vs  depth                      *
c************************************************************************
      subroutine str(fe,mc,nc,xc,yc,zc,tx,ty,tz,df,ireg,fdib,pen,pep,k)
      use IFQWIN
      use IFPORT
      parameter(mi=200)
      type (xycoord) xy
      character*120 texto
      dimension xc(mc),yc(mc),zc(mc),tx(mc),ty(mc),tz(mc),
     & df(mc),ireg(mc), zca(mi),nca(mi),npo(mi),nne(mi),sz(mi)
     
c   Num de niveles, cotas y numero de celdas no nulas
      
      niv=0
      do 2 j=1,nc       
      if(ireg(j).ne.0) go to 2 
      if(niv.eq.0) then
        niv=niv+1
        zca(niv)=zc(j)
        sz(niv)=tz(j)
        nca(niv)=nca(niv)+1        
        if(df(j).gt.0) npo(niv)=npo(niv)+1
        if(df(j).lt.0) nne(niv)=nne(niv)+1
        go to 2
      endif
      do 1 i=1,niv
      d=abs(zc(j)-zca(i))
      if(d.gt.2) go to 1
        nca(i)=nca(i)+1
        if(df(j).gt.0) npo(i)=npo(i)+1
        if(df(j).lt.0) nne(i)=nne(i)+1
        go to 2
    1 continue 
       niv=niv+1    
       zca(niv)=zc(j)
       sz(niv)=tz(j)
       nca(niv)=nca(niv)+1
       if(df(j).gt.0) npo(niv)=npo(niv)+1
       if(df(j).lt.0) nne(niv)=nne(niv)+1
    2 continue       
      
c     Ajuste de pendientes          
      
      iy=nc/50./niv
      sx=0 
      sy=0 
      sx2=0
      sy2=0
      sxy=0
      sp=0 
      do 3 i=1,niv
      if(nca(i).lt.iy) go to 3
      zd=nne(i)*300./nca(i)   
      if(zd.eq.0) go to 3
      r=1.*nca(i)
      sy=sy+zd*r 
      c=-zca(i)/fe/fdib
      sx=sx+ c*r 
      sy2=sy2+zd*zd*r
      sx2=sx2+ c* c*r
      sxy=sxy+ c*zd*r
      sp=sp+r 
   3  continue 
        pen=(sp*sxy-sx*sy)/(sp*sx2-sx*sx)
        cn=sy/sp-pen*sx/sp 
      sx=0 
      sy=0 
      sx2=0
      sy2=0
      sxy=0
      sp=0 
      do 4 i=1,niv    
      if(nca(i).lt.iy) go to 4
      zd=npo(i)*300./nca(i) 
      if(zd.eq.0.) go to 4
      r=1.*nca(i)
      sy=sy+zd*r   
      c=-zca(i)/fe/fdib
      sx=sx+ c*r          
      sy2=sy2+zd*zd*r     
      sx2=sx2+ c* c*r     
      sxy=sxy+ c*zd*r     
      sp=sp+r             
    4 continue            
        pep=(sp*sxy-sx*sy)/(sp*sx2-sx*sx)
        cp=sy/sp-pep*sx/sp   

      if(k.eq.0) return  
c---------------- Graphics
      write(texto,'(a,8x,a,10x,a)')                                     !
     & 'Negat. Mas.','Posit. Mas.','Combin. Mas.'                       !
      call moveto(80,470,xy)                                            !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      sx=0                                                              !
      sy=0                                                              !
      sx2=0                                                             !
      sy2=0                                                             !
      sxy=0                                                             !
      sp=0                                                              !
      do 5 i=1,niv                                                      !
      if(nca(i).lt.iy) go to 5                                          !
      r=sz(i)/fe/fdib/2.-1.                                             !
      c=-zca(i)/fe/fdib                                                 !
      v2=500+c+r                                                        !
      v1=500+c-r                                                        !
      if(v2.gt.850) go to 5                                             !
      zd=nne(i)*250./nca(i)                                             !
        nada=setcolorrgb(#ff0000)                                       ! azul
        u2=160-zd                                                       !
        if(zd.ne.0) nada=rectangle($gfillinterior,160,v1,u2,v2)         !
      zd=npo(i)*250./nca(i)                                             !
        nada=setcolorrgb(#0000ff)                                       ! rojo
        u2=200+zd                                                       !
        if(zd.ne.0) nada=rectangle($gfillinterior,200,v1,u2,v2)         !
      zd=(npo(i)-nne(i))*250./nca(i)                                    !
        nada=setcolorrgb(#00a000)                                       ! verde
        u2=400+zd                                                       !
        if(zd.ne.0) nada=rectangle($gfillinterior,400,v1,u2,v2)         !
    5 continue                                                          !
         nada=setcolorrgb(#000000)                                      !
         iy=-5                                                          !
         ix2=cn+iy*pen                                                  !
         call moveto(160-ix2,500+iy,xy)                                 !
         iy=80                                                          !
         u1=cn+iy*pen                                                   !
         nada=lineto(160-u1,500+iy)                                     !
         
         iy=-5                                                          !
         ix2=cp+iy*pep                                                  !
         call moveto(200+ix2,500+iy,xy)                                 !
         iy=80                                                          !
         u1=cp+iy*pep                                                   !
         nada=lineto(200+u1,500+iy)                                     !
      return
      end
c***********************************************************************