c***********************************************************************
c                         V I E W                          12/01/2020  *
c                                                                      *
c       Graphic presentation for GROWTH gravity inversion.             *
c                                                                      *
c       Compiled with Intel Visual Fortran,                            *
c            (Quickwin Application Project)                            *
c                                                                      *
c  Aim:      Visualization of 3-D models coming from the program       *
c            GROWTH.FOR of gravity inversion.                          *
c            Display of vertical profiles and horizontal sections of   *
c            the model.                                                *
c            Display of maps of anomaly, altitudes, residuals, regional*
c            trend, local anomaly, model anomaly.                      *
c  Author:   Antonio G. Camacho  (antonio_camacho@mat.ucm.es)          *
c            Instituto de Astronomia y Geodesia (CSIC-UCM)             *
c            Facultad de CC. Matematicas                               *
c            Ciudad Universitaria,     28040 Madrid (Spain)            *
c  Reference: Camacho et al. "The 3-D gravity inversion package        *
c             GROWTH2.0 and its application to Tenerife Island, Spain" *
c               Computers & Geosciences 37 (2011) 621–633              *
c               See also the file GROWTH.TXT                           *
c ---------------------------------------------------------------------*
c                                                                      *
c  INPUT                                                               *
c  -----                                                               *
c  file "GRA.txt": coordinates, altitudes, observed gravity values     *
c        and relative errors for the gravity stations                  *
c                                                                      *
c  file "MOD.txt": 3-D model of the anomalous density contrast         *
c        distribution as obtained by means of GROWTH.FOR               *
c                                                                      *
c  file "MAP.BLN": base map for overlapping upon the figures           *
c                                                                      *
c  file "FIL.txt": regional,local,adjusted and residual gravity values *
c                                                                      *
c                                                                      *
c  OUTPUT                                                              *
c  ------                                                              *
c  several graphic views on the screen (anomaly, altitudes, residues..)*
c                                                                      *
c  files "OUT.txt" for residues, calculated gravity anomaly and        *
c       mass and density distributions,                                *
c  file "MAT.txt" for a 3D matrix for anomalous density values         *
c                                                                      *
c  grid files "*.GRD" (and base map files "*.bln") for further drawing *
c  of sections and profiles of the model into Surfer (Golden Software) *
c                                                                      *
c ---------------------------------------------------------------------*
c  Dimensions:  ms: maximum number of gravity stations                 *
c               mc: maximum number of cells for the model              *
c               mb: maximum number of points for the base map file     *
c***********************************************************************
      use IFQWIN   !msflib
      TYPE (wxycoord) wxy
      TYPE (xycoord) xy
      TYPE (rccoord) curpos

c  Dimensions

      parameter (ms=500000,mc=230000,mb=6000,mt=5000,npc=20,mg=200)
      character*20 fobs,fmod,fres,fmap
      character*120 texto,fname
      character*9 hoy
      dimension g(ms),gc(ms),gr(ms), r(npc,mt) ,sb(mc)
      integer zs(mc),zi(mc), ij(mt),jb(mc), col(30), xp(mb),yp(mb),
     -  ii(mc),lr(mc),      
     - x(ms),y(ms),z(ms),ieg(ms), jp(mb),nada ,dbm(mc),
     - xb(mc),yb(mc),dx(mc),dy(mc),db(mc),eb(mc),cb(mc),sen(mc),fq(mc)
      data nc/23/,
     -fobs/'Gra.dat'/,fmod/'mod.dat'/,fres/'fil.dat'/,fmap/'map.bln'/

      col( 1)=#a01000     ! b g r
      col( 2)=#c03000
      col( 3)=#d05000
      col( 4)=#e06000
      col( 5)=#f07000  ! 7
      col( 6)=#f09000
      col( 7)=#f0b500  !5
      col (8)=#f0d000
      col( 9)=#f0e020
      col(10)=#f0f060
      col(11)=#f0f0b0
      col(12)=#f0f0f0  ! blanco
      col(13)=#80f0f0  
      col(14)=#00e7f0
      col(15)=#00d6f0  ! 3
      col(16)=#00c0f0  ! 4
      col(17)=#00a6f0  ! 5
      col(18)=#0085f0  ! 6    
      col(19)=#0060f0  ! 7
      col(20)=#0030f0
      col(21)=#0000e0
      col(22)=#0000c0
      col(23)=#0000a0
      nada=settextcolorrgb(#000000)
      nada=setbkcolorrgb(#ffffff)
      nada=initializefonts()

c  Reading file names

      call clearscreen($clearscreen)
      call dial1(fobs,fmod,fres,fmap)

c   Reading gravity data  GRA.txt

      open(1,file=fobs)
      ax=9.d9
      bx=-ax
      ay=ax
      by=-ay
      zt=0.
      j=0
      cxm=0
      cym=0
      czm=0
    2 read(1,*,err=95,end=6) xx,yy,zz,gg,ee
      if(xx.eq.0.and.yy.eq.0.and.zz.eq.0) go to 6
      j=j+1
      x(j)=xx
      y(j)=yy
      z(j)=zz
      cxm=cxm+xx   
      cym=cym+yy
      czm=czm+zz      
      if(xx.lt.ax) ax=x(j)
      if(xx.gt.bx) bx=x(j)
      if(yy.lt.ay) ay=y(j)
      if(yy.gt.by) by=y(j)
      g(j)=gg
      gr(j)=0.
      ieg(j)=ee
      zt=zt+z(j)
      go to 2
    6 close(1)
      n=j
      if(n.eq.0) stop
      zt=zt/n
      cxm=cxm/n
      cym=cym/n
      czm=czm/n
      d=0
      do 49 j=1,n
      xx=x(j)-cxm
      yy=y(j)-cym
      zz=z(j)-czm 
      d=d+xx*xx+yy*yy+zz*zz 
   49 continue
      fe=sqrt(d/n)/8000.

c  Reading anomalous model MOD.txt

      open(2,file=fmod)
      read(2,*,err=96) nb,siz,dr 
      if(nb.gt.mc) call out(mc)
      tb=0.
      do 1 i=1,nb
      read(2,*) xb(i),yb(i),zz,dx(i),dy(i),a3,dd,dm,ee, ff
      zs(i)=zz+a3/2.
      zi(i)=zz-a3/2.
      lr(i)=0
      fq(i)=0.
c     if(dd.ne.0) fq(i)=ff*dd/abs(dd)
        fq(i)=dd
      db(i)=dd
      dbm(i)=dm
      xx=xb(i)
      yy=yb(i)
      zz=(zs(i)+zi(i))/2.
      eb(i)=ee
      call SST(ms,n,x,y,z,xx,yy,zz,avm)
      sen(i)=avm*1.e9*fe*fe
      tb=tb+(a1+a2+a3)/3.
    1 continue
      if(nb.gt.0) tb=tb/nb
      close(2)

c  Limits for pictures

      tic=(bx-ax+by-ay)/2./10.
      i=log10(tic)
      tic=nint(tic/10**i)*10**i
      if(tic.lt.1.) tic=1.
      tam=(bx-ax)*1.1
      tt=(by-ay)*1.1
      if(tt.gt.tam) tam=tt
      xx=(bx+ax)/2.
      yy=(by+ay)/2.
      ax=xx-tam/2.
      bx=xx+tam/2.
      ay=yy-tam/2.
      by=yy+tam/2.
      if(n.gt.ms) then
      write(*,*) n,'>', ms
      stop
      endif

c  Reading base map file  MAP.BLN

      open(3,file=fmap)
      k=1
    4 read(3,*,end=5) j
      do 3 i=k,(j+k-1)
      read(3,*) xp(i),yp(i)
      jp(i)=0
      if(i.ge.mb) go to 5
    3 continue
      jp(k)=1
      k=k+j
      go to 4
    5 np=k-1
      close(3)

c  Number of data items

      call dial2(fobs,fmod,fres,fmap,n,nb,np)
      call date(hoy)

c   Options

   10 call clearscreen($clearscreen)
      nada=setcolorrgb(#000050)
      nada=setwindow(.true.,0.,0.,1000.,720.)
      call setviewport(0,0,1000,720)
      call titulo(hoy)
      nada=initializefonts()
      write(*,'(26x,a)') '   1  Limits of figures, scale   '
      write(*,'(26x,a)') '          DATA FILE              '
      write(*,'(26x,a)') '   2  Benchmark elevations       '
      write(*,'(26x,a)') '   3  Gravity anomaly            '
      write(*,'(26x,a)') '   4  Anomaly Stdev.             '
      write(*,'(26x,a)') '          FILE OF RESIDUALS      '
      write(*,'(26x,a)') '   5  Regional trend             '
      write(*,'(26x,a)') '   6  "Observed" local anomaly   '
      write(*,'(26x,a)') '   7  Modelled local anomaly     '
      write(*,'(26x,a)') '   8  Residual anomaly           '
      write(*,'(26x,a)') '          MODEL FILE             '
      write(*,'(26x,a)') '   9  Profiles: density contrast '
      write(*,'(26x,a)') '  10  Profiles: domain sensitivity '
      write(*,'(26x,a)') '  11  Mass, vol. and density vs. depth'
      write(*,'(26x,a)') '  12  Mass of isolated bodies    '
      write(*,'(26x,a)') '  13  3-D matrix for matlab  (MAT.txt)'
      write(*,'(26x,a)') '  14  Discontinuity surface        '
      write(*,'(26x,a)') '  15  Profiles: model sensitivity  '
      write(*,'(26x,a)') '  16  Perspect of the model        '
      write(*,'(26x,a)') '  17  Remove isolated cells        '
      write(*,'(26x,a)') '  18  Local/Regional cells         '
      write(*,'(26x,a)') '  19  Pressure Ton/m2, loading maps'
      write(*,'(26x,a)') '  20  Write (out.txt) filled cells '

      write(*,'(26x,a)') '                                 '
      write(*,'(26x,a)') '  99  Exit                       '
      write(*,'(//40x,a\)') ' >> Option ? '
      read(*,'(i2)') iop
      if(iop.eq.99) go to 99
      if(iop.le.0.or.iop.gt.20) go to 10

c  Option 1: limits of figures, axis

      if(iop.eq.1) then
      call clearscreen($clearscreen)
      call titulo(hoy)
      write(*,'(20x,a)') ' Input 0 for default '
      write(*,'(//20x,a\)') ' >> Factor for size change (e.g. 1.2) ? '
      read(*,'(f4.1)') dd
      if(dd.ne.0.) then
       xx=(ax+bx)/2.
       ax=xx+(ax-xx)*dd
       bx=xx+(bx-xx)*dd
       yy=(ay+by)/2.
       ay=yy+(ay-yy)*dd
       by=yy+(by-yy)*dd
      endif
      write(*,'(/20x,a,f9.0,a\)') ' >> Limit AX (def.', ax,') ? '
      read(*,'(f8.0)') xx
      if(xx.eq.0.) go to 8
      write(*,'(20x,a,f9.0,a\)') ' >> Limit BX (def.',bx,') ? '
      ax=xx
      read(*,'(f8.0)') xx
      bx=xx
    8 write(*,'(/20x,a,f9.0,a\)') ' >> Limits AY (def.',ay,') ? '
      read(*,'(f9.0)') xx
      if(xx.eq.0.) go to 7
      write(*,'( 20x,a,f9.0,a\)') ' >> Limits BY (def.',by,') ? '
      ay=xx
      read(*,'(f9.0)') xx
      by=xx
    7 tam=bx-ax
      if((by-ay).gt.tam) tam=by-ay
      xx=(bx+ax)/2.
      yy=(by+ay)/2.
      ax=xx-tam/2.
      bx=xx+tam/2.
      ay=yy-tam/2.
      by=yy+tam/2.
       tic=(bx-ax+by-ay)/2./10.
       i=log10(tic)
       tic=nint(tic/10**i)*10**i
       if(tic.lt.1.) tic=1.
      write(*,'(/20x,a,f6.0,a\)')
     -' >> Axis unit (def. ',tic,' m ) ? '
      read(*,'(f8.0)') xx
      if(xx.ne.0.) tic=xx
      go to 10
      endif

c     Option 2:  map of data altitudes

      if(iop.eq.2) then
      do 9 i=1,n
    9 gc(i)=z(i)
      tic=-tic
      fname='Benchmark elevations (m)'
      call cmap(ms,ax,bx,ay,by,nc,n,x,y,gc,ieg,tic,fe,col,dr,
     -mb,np,xp,yp,jp,fname,hoy,iop,fmod)
      read(*,'(i1)') i
      go to 10
      endif

c     Option 3:  map of data anomaly

      if(iop.eq.3) then
      fname='"Obs." Gravity Anomaly'
      call cmap(ms,ax,bx,ay,by,nc,n,x,y,g,ieg,tic,fe,col,dr,
     -mb,np,xp,yp,jp,fname,hoy,iop,fmod)
      read(*,'(i1)') i
      go to 10
      endif

c     Option 4:  map of comparative data errors

      if(iop.eq.4) then
      do 11 i=1,n
   11 gc(i)=ieg(i)
      fname='Data St.Dv.'
      call cmap(ms,ax,bx,ay,by,nc,n,x,y,gc,ieg,tic,fe,col,dr,
     -mb,np,xp,yp,jp,fname,hoy,iop,fmod)
      read(*,'(i1)') i     
      go to 10
      endif

c     Option 5:  map of regional anomaly

      if(iop.eq.5) then
      open(1,file=fres)
      do 12 i=1,n
      read(1,*,err=97) xx,yy,gc(i)
   12 gc(i)=gc(i)+gr(i)
      close(1)
      fname='Regional Grav. Anomaly'
      call cmap(ms,ax,bx,ay,by,nc,n,x,y,gc,ieg,tic,fe,col,dr,
     -mb,np,xp,yp,jp,fname,hoy,iop,fmod)
      read(*,'(i1)') i
      go to 10
      endif

c     Option 6: map of local anomaly

      if(iop.eq.6) then
      open(1,file=fres)
      do 13 i=1,n
      read(1,*,err=97) xx,yy,gg,gc(i)
   13 gc(i)=gc(i)-gr(i) 
      close(1)
      fname='Obs. Local Grav. Anomaly'
      call cmap(ms,ax,bx,ay,by,nc,n,x,y,gc,ieg,tic,fe,col,dr,
     -mb,np,xp,yp,jp,fname,hoy,iop,fmod)
      read(*,'(i1)') i
      go to 10
      endif

c     Option 7: map of model anomaly

      if(iop.eq.7) then
      open(1,file=fres)
      do 14 i=1,n
      read(1,*,err=97) xx,yy,gg,tt,gc(i)
   14 gc(i)=gc(i)-gr(i)
      close(1)
      fname='Modelled Local Anomaly'
      call cmap(ms,ax,bx,ay,by,nc,n,x,y,gc,ieg,tic,fe,col,dr,
     -mb,np,xp,yp,jp,fname,hoy,iop,fmod)
      read(*,'(i1)') i
      go to 10
      endif

c     Option 8: map of residual anomaly

      if(iop.eq.8) then
      open(1,file=fres)
      do 15 i=1,n
      read(1,*,err=97) xx,yy,zz,gg,tt,gc(i),c
   15 ieg(i)=1.             ! 5./sqrt(c+0.001)
      close(1)
      fname='Inversion Residuals'
      call cmap(ms,ax,bx,ay,by,nc,n,x,y,gc,ieg,tic,fe,col,dr,
     -mb,np,xp,yp,jp,fname,hoy,iop,fmod)
      read(*,'(i1)') i
   60 nada=setcolorrgb(#a0ffa0)
      nada=rectangle($gfillinterior,800,672,823,687)
      c=0.
      do 58 i=1,n
   58 c=c+gc(i)
      c=c/n
      tt=0
      do 59 i=1,n
      rl=gc(i)-c 
   59 tt=tt+rl*rl
      tt=sqrt(tt/n) 
      call settextposition(45,70,curpos)  
      write(*,'(8x,a4,i4,4x,a5,f7.0,4x,a4,f6.0)')
     -'Num=',n,'Mean=',c,'SD=',tt
      endif

c     Option 9: vertical and horizontal sections of the density model

      if(iop.eq.9) then
      write(*,'(/30x,a\)') '   >> Interpolation NO(0)/yes(1) ? '
      read(*,'(i1)') is
      write(*,'(/30x,a\)') ' >> Stratified background: YES(0)/no(1) ? '
      read(*,'(i1)') j
      do 17 i=1,nb
      cb(i)=dbm(i)          
      if(j.eq.1) cb(i)=db(i)
   17 continue
      call profiles(iop,is,fe,tb,mc,nb,xb,yb,zi,zs,dx,dy,eb,
     -cb,lr,ax,ay,tam,nc,col,zt,czm,tic,mb,np,xp,yp,jp,ms,n,x,y,z)
      endif

c     Option 10:  sections of the domain sensitibity

      if(iop.eq.10) then
      write(*,'(//30x,a\)')  ' >> Mas (0), depth (1), horiz (2) ? '
      read(*,'(i3)') k
      iop=k+30
      is=1
      do 16 i=1,nb
   16 cb(i)=abs(eb(i))
      call clearscreen($clearscreen)
      call profiles(iop,is,fe,tb,mc,nb,xb,yb,zi,zs,dx,dy,eb,
     -cb,lr,ax,ay,tam,nc,col,zt,czm,tic,mb,np,xp,yp,jp,ms,n,x,y,z)
      endif

c     Option 11:  mass and density distributions vs. depth

      if(iop.eq.11) then
      rl=50
      write(*,'(//30x,a\)')  ' >> Domain limit (0-100) ? '
      read(*,'(i3)') k
      call clearscreen($clearscreen)
      nada=setcolorrgb(#000050)
      nada=setfont('t''Arial''h12w6i')
      write(texto,'(a)') 'GROWTH 2020'
      call moveto(40,7,xy)
      call outgtext(texto)
      write(texto,'(a9)') hoy
      call moveto(750,8,xy)
      nada=setfont('t''Arial''h14w7i')
      call outgtext(texto)
      if(k.ne.0) rl=k
      ninv=0
      do 23 i=1,nb
      sb(i)=9999.
      if(db(i).eq.0) go to 23
      
      do 19 j=1,nb
      if(zi(j).ne.zs(i)) go to 18
      if((db(j)-10).le.db(i)) go to 18
      if((xb(i)+dx(i)/2).le.xb(j)) go to 18
      if((xb(i)-dx(i)/2).ge.xb(j)) go to 18
      if((yb(i)+dy(i)/2).le.yb(j)) go to 18
      if((yb(i)-dy(i)/2).ge.yb(j)) go to 18
      ninv=ninv+1   
   18 continue
      if(zi(i).ne.zs(j)) go to 19
      if((db(i)-10).le.db(j)) go to 19
      if((xb(j)+dx(j)/2).le.xb(i)) go to 19
      if((xb(j)-dx(j)/2).ge.xb(i)) go to 19
      if((yb(j)+dy(j)/2).le.yb(i)) go to 19
      if((yb(j)-dy(j)/2).ge.yb(i)) go to 19
      ninv=ninv+1   
   19 continue

      xx=xb(i)
      yy=yb(i)
      zz=(zs(i)+zi(i))/2.
      call SST(ms,n,x,y,z,xx,yy,zz,avm)
      sb(i)=avm*fe*fe*1.e9
   23 continue
      ee=abs(eb(1))
      ij(1)=zi(1)
      tt=zs(1)
      bb=zi(1)
      dd=dx(1)*dy(1)
      ps=zs(1)-zi(1)
      tot=ps*ps*ps
      nsd=0
      do 21 i=2,nb
      if(zs(i).gt.tt) tt=zs(i)
      if(zi(i).lt.bb) bb=zi(i)
      ee=ee+abs(eb(i))
      zz=dx(i)*dy(i)
      if(zz.lt.dd) dd=zz
      c=zs(i)-zi(i)
      tot=tot+c*zz
      if(c.lt.ps) ps=c
      do 22 j=1,nsd
      if(zi(i).eq.ij(j)) go to 21
   22 continue
      nsd=nsd+1
      if(nsd.gt.mt) then
      write(*,'(/6x,a,i3)')  ' **** Error: number of depth levels >',mt
      write(*,'(/20x,a)')  ' ===> Press (Enter)'
      read(*,'(i1)') kk
      go to 10
      endif
      ij(nsd)=zi(i)
   21 continue
      tot=tot/nsd
      ee=ee/nb
      dd=sqrt(dd)
      eg=ee*ee
      gc(1)=0.
      nx=tam/dd
      a1=9.d9
      a2=-a1
      tm=0
      tv=0
      do 26 k=1,nsd                  !values for each sheet
      do 24 l=1,npc
      gc(l)=0
   24 r(l,k)=0.
      zz=ij(k)
      sc=0
      sn=0
      sp=0
      do 27 l=1,nb
      ee=eb(l)*eb(l)/eg
      if(zi(l).ne.ij(k)) go to 25
       if(db(l).eq.0.or.sb(l).gt.rl) go to 27
       c=dx(l)*dy(l)*(zs(l)-zi(l))           ! cell volume
       dd=c*db(l)                            ! cell mass
       r(9,k)=r(9,k)+c                       ! total volume
       r(1,k)=r(1,k)+abs(dd)                 ! total absolute mass
       if(db(l).gt.0) r(2,k)=r(2,k)+dd       ! posit mass
       if(db(l).gt.0) r(10,k)=r(10,k)+c      ! posit volume
       if(db(l).lt.0) r(3,k)=r(3,k)+dd       ! negative mass
       if(db(l).lt.0) r(11,k)=r(11,k)+c      ! negative volume
       r(4,k)=r(4,k)+dd                      ! total mass
   25 zr=0.
      c=0.
      yy=zz+ps
      if(zs(l).lt.yy) yy=zs(l)
      xx=zz
      if(zi(l).gt.xx) xx=zi(l)
      if(yy.gt.xx) c=c+(yy-xx)
        yy=zz+ps*2.
        if(zs(l).lt.yy) yy=zs(l)
        xx=zz+ps
        if(zi(l).gt.xx) xx=zi(l)
        if(yy.gt.xx) c=c+(yy-xx)*0.5
      yy=zz
      if(zs(l).lt.yy) yy=zs(l)
      xx=zz-ps
      if(zi(l).gt.xx) xx=zi(l)
      if(yy.gt.xx) c=c-(yy-xx)
        yy=zz-ps
        if(zs(l).lt.yy) yy=zs(l)
        xx=zz-ps*2
        if(zi(l).gt.xx) xx=zi(l)
        if(yy.gt.xx) c=c-(yy-xx)*0.5
      if(c.eq.0.) go to 27
      zr=c*db(l)/ee
      do 28 i=1,nx
      xx=ax+dd*i
      c=dx(l)/2.
      if(xx.gt.(xb(l)+c).or.xx.lt.(xb(l)-c)) go to 28
      do 29 j=1,nx
      yy=ay+dd*j
      c=dy(l)/2.
      if(yy.gt.(yb(l)+c).or.yy.lt.(yb(l)-c)) go to 29
      r(8,k)=r(8,k)+zr
   29 continue
   28 continue
   27 continue
      if(r(1,k).gt.tm) tm=r(1,k)
      if(r(9,k).gt.tv) tv=r(9,k)
      if(r(8,k).gt.a2) a2=r(8,k)
      if(r(8,k).lt.a1) a1=r(8,k)
      if(r(9,k).ne.0.) r(5,k)=r(1,k)/r(9,k)/800.
      if(r(10,k).ne.0.) r(6,k)=r(2,k)/r(10,k)/800.
      if(r(11,k).ne.0.) r(7,k)=r(3,k)/r(11,k)/800.
   26 continue

      if((a2-a1).gt.0.) r(8,1)=r(8,1)/(a2-a1)
      do 30 k=1,nsd
      if(k.gt.1.and.(a2-a1).gt.0.) r(8,k)=r(8,k)/(a2-a1)
      r(1,k)=r(1,k)/tm
      r(2,k)=r(2,k)/tm
      r(3,k)=r(3,k)/tm
      r(4,k)=r(4,k)/tm
      r(9,k)=r(9,k)/tv
      r(10,k)=r(10,k)/tv
   30 r(11,k)=r(11,k)/tv

      nada=setcolorrgb(#000000)
      nada=setfont('t''Arial''h16w8e')
      write(texto,'(a)')
     -'Interfaces, Mass, Volume & Density contrasts   vs.  depth (m)'
      call moveto(250,20,xy)
      call outgtext(texto)
      nada=setfont('t''Courier''h16w8')
      write(texto,'(a,i6)') 
     -'Total:   Positive:   Negative:   Differ:          Invers:',ninv
      call moveto(300,50,xy)
      call outgtext(texto)
      nada=setcolorrgb(#a0a0a0)
      nada=rectangle($gborder,352,52,353,64)
      nada=setcolorrgb(#6060f0)
      nada=rectangle($gborder,448,52,449,64)
      nada=setcolorrgb(#f08080)
      nada=rectangle($gborder,544,52,545,64)
      nada=setcolorrgb(#40f040)
      nada=rectangle($gborder,626,52,627,64)
      nada=setcolorrgb(#000000)
      write(texto,'(a,8x,a,e8.2,a,10x,a,e8.2,a)')
     -'Mean density contrast (kg/m3)',
     -'Anomalous volume (',tv/100,' m3) ',
     -'Anomalous mass (',tm/100,' kg)'
      call moveto(120,85,xy)
      call outgtext(texto)
      write(texto,'(a)') '-800    -400       0       400      800'
      call moveto(50,105,xy)
      call outgtext(texto)
      write(texto,'(a)') '-100    -50        0        50     100'
      call moveto(365,105,xy)
      call outgtext(texto)
      call moveto(675,105,xy)
      call outgtext(texto)

      c=(tt-bb)/40.
      bb=bb-c
      tt=tt+c
      c=(tt-bb)/12.*fe
      i=log10(c)
      a1=nint((tt*fe+czm)/10**i)*10**i
      c=nint(c/10**i)*10**i
      mx=300
      my=480

      call setviewport(60,120,60+mx,120+my)
      nada=setwindow(.true.,-1.,bb,1.,tt)
      nada=setcolorrgb(#e0e0e0)
      nada=rectangle($gfillinterior,0,0,mx,my)
      nada=setcolorrgb(#005000)
      zz=zs(1)
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(1,zz)
      zz=ij(1)
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(1,zz)
      yr=(zz+zs(1))/2.
      do 31 k=2,nsd
      xr=(ij(k)+ij(k-1))/2.
      zz=ij(k)
      nada=setcolorrgb(#005000)
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(1,zz)                             !divisoria capa
      nada=setcolorrgb(#a0a0a0)
      call moveto_w(dble(r(5,k-1)),dble(yr),wxy)
      nada=lineto_w(r(5,k),xr)                          ! dens total
      nada=setcolorrgb(#6060f0)
      call moveto_w(dble(r(6,k-1)),dble(yr),wxy)
      nada=lineto_w(r(6,k),xr)                          ! dens posit
      nada=setcolorrgb(#f08080)
      call moveto_w(dble(r(7,k-1)),dble(yr),wxy)
      nada=lineto_w(r(7,k),xr)                          ! dens negat
   31 yr=xr
      call ejes(mx,my,tt,bb,c,a1,czm,fe)

      call setviewport(70+mx,120,70+mx+mx,120+my)
      nada=setwindow(.true.,-1.,bb,1.,tt)
      nada=setcolorrgb(#e0e0e0)
      nada=rectangle($gfillinterior,0,0,mx,my)
      nada=setcolorrgb(#005000)
      zz=zs(1)
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(1,zz)
      zz=ij(1)
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(1,zz)
      yr=(zz+zs(1))/2.
      do 32 k=2,nsd
      xr=(ij(k)+ij(k-1))/2.
      zz=ij(k)
      nada=setcolorrgb(#005000)
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(1,zz)
      nada=setcolorrgb(#a0a0a0)
      call moveto_w(dble(r(9,k-1)),dble(yr),wxy)
      nada=lineto_w(r(9,k),xr)                          ! volum total
      nada=setcolorrgb(#6060f0)
      call moveto_w(dble(r(10,k-1)),dble(yr),wxy)
      nada=lineto_w(r(10,k),xr)                          ! volum total
      nada=setcolorrgb(#f08080)
      call moveto_w(dble(-r(11,k-1)),dble(yr),wxy)
      nada=lineto_w(-r(11,k),xr)                          ! volum total
   32 yr=xr
      call ejes(mx,my,tt,bb,c,a1,czm,fe)

      call setviewport(80+mx+mx,120,80+mx+mx+mx,120+my)
      nada=setwindow(.true.,-1.,bb,1.,tt)
      nada=setcolorrgb(#e0e0e0)
      nada=rectangle($gfillinterior,0,0,mx,my)
      nada=setcolorrgb(#005000)
      zz=zs(1)
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(1,zz)
      zz=ij(1)
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(1,zz)
      yr=(zz+zs(1))/2.
        kk=0
      do 33 k=2,nsd
      xr=(ij(k)+ij(k-1))/2.
      zz=ij(k)
      nada=setcolorrgb(#005000)
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(1,zz)
      nada=setcolorrgb(#f0f080)
      call moveto_w(dble(r(9,k-1)),dble(yr),wxy)
      nada=setcolorrgb(#a0a0a0)
      call moveto_w(dble(r(1,k-1)),dble(yr),wxy)
      nada=lineto_w(r(1,k),xr)                          ! masa absol
      nada=setcolorrgb(#40f040)
      call moveto_w(dble(r(4,k-1)),dble(yr),wxy)
      nada=lineto_w(r(4,k),xr)                          ! masa difer
         if(r(4,k).lt.r(4,k-1).and.kk.eq.0) then
         kk=1
         call moveto_w(dble(r(2,k-1)),dble(yr),wxy)
         nada=setcolorrgb(#00ffff)
         nada=rectangle_w($gfillinterior,0,yr,r(4,k-1),yr)
         endif
         if(r(4,k).ge.r(4,k-1)) kk=0
      nada=setcolorrgb(#6060f0)
      call moveto_w(dble(r(2,k-1)),dble(yr),wxy)
      nada=lineto_w(r(2,k),xr)                          ! masa posit
      nada=setcolorrgb(#f08080)
      call moveto_w(dble(r(3,k-1)),dble(yr),wxy)
      nada=lineto_w(r(3,k),xr)                          ! masa negat
   33 yr=xr
      call ejes(mx,my,tt,bb,c,a1,czm,fe)

      call setviewport(0,0,1000,720)
      nada=setwindow(.true.,0.,0.,1000.,720.)
      ps=my/(tt-bb)/fe
      xx=c*ps
      ps=112+(tt*fe+czm-a1)*ps
      do 34 i=1,30
      yy=a1-c*(i-1)
      zz=(yy-czm)/fe
      if(zz.lt.bb) go to 34
      write(texto,'(i6)') nint(yy)
      call moveto(6,nint(ps+(i-1)*xx),xy)
      call outgtext(texto)
   34 continue
      call settextposition(42,60,curpos)
      write(*,'(a\)') ' >> Output to file OUT.txt: y(1)/N(0) ? '
      read(*,'(i1)') i
      if(i.ne.1) go to 10
      open(1,file='out.txt')
      write(1,'(a,e8.2,a/a)')
     -' Depth     Anomalous mass (',tm/100,'kg)      Density contrast ',
     -' m asl  Total  Positive Negative  Differ  Mean  Positive  Negat'
      do 35 k=1,nsd
   35 write(1,'(f7.0,4f8.2,4f8.1)') ij(k)*fe+czm,(r(l,k)*100,l=1,4),
     -(r(l,k)*1000,l=5,8)
      close(1)
      go to 10
      endif

c     Option 12:  calculating the mass of isolated anomalous bodies

      if(iop.eq.12) then
      call clearscreen($clearscreen)
      call titulo(hoy)
      write(*,'(/8x,a,5f9.0,f6.0)') ' Limits of the model =',ax
     -,bx,ay,by,zt-tam/2,zt
   20 write(*,'(/4x,a\)')
     -' >> Dens. contrast :  neg(1), pos(2), tot(3), exit (0) ? '
      read(*,'(i1)') j
      if(j.eq.0) go to 10
      if(j.ne.1.and.j.ne.2.and.j.ne.3) go to 20
      write(*,'(/4x,a)')
     -' >> Limits (1), total (0) ? '
      read(*,'(i1)') is
      if(is.eq.1) then
        write(*,'(/4x,a)')
     -  ' >> Limits of the body:  ax, bx, ay, by, zbot, ztop ? '
        read(*,*) xx,a1,yy,a2,bb,tt
      endif
      xr=0
      yr=0
      zz=0
      c=0.
      ps=0
      k=0  
      do 36 i=1,nb
      if(j.eq.3.and.db(i).eq.0) go to 36
      if(j.eq.1.and.db(i).ge.0) go to 36
      if(j.eq.2.and.db(i).le.0) go to 36
      if(is.eq.1) then
      if((xb(i)-dx(i)/2).lt.xx.or.(xb(i)+dx(i)/2).gt.a1) go to 36
      if((yb(i)-dy(i)/2).lt.yy.or.(yb(i)+dy(i)/2).gt.a2) go to 36
      if(zi(i).lt.bb.or.zs(i).gt.tt) go to 36
      endif
      p=1.*abs(db(i))*(zs(i)-zi(i))*dx(i)*dy(i)
      ps=ps+abs(db(i))
      k=k+1
      c=c+p
      xr=xr+xb(i)*p
      yr=yr+yb(i)*p
      zz=zz+(zi(i)+zs(i))/2.*p
   36 continue
      ps=ps/k
      if(c.eq.0.) go to 20
      write(*,'(8x,a,e9.3,a)') ' -->  Anomalous mass =',c,' kg'
      a2=c/ps/1.d9
      write(*,'(8x,a,f6.0,a,e9.3,a)') ' -->  Volume (',ps,')=',a2,' km3'
      xr=xr/c
      yr=yr/c
      zz=zz/c
      write(*,'(8x,a,2f9.0,f6.0)') ' -->  Center of mass =',xr,yr,zz
      go to 20
      endif

c     Option 13: calculating a 3-D matrix

      if(iop.eq.13) then
      open(1,file='MAT.txt',status='unknown')
      open(2,file='Out.dat')
      call clearscreen($clearscreen)
      call titulo(hoy)
c      rl=50
c      write(*,'(/10x,a\)') ' >> Domain limit: 0,..,100 (def. 50) ? '
c      read(*,*) yy
c      if(yy.ne.0.) rl=yy
      write(*,'(//a)') '    Limits for the matrix'
      write(*,'(/10x,a,f9.0,a\)') ' >> Xmin (def.', ax,') ? '
      read(*,'(f8.0)') xx
      if(xx.ne.0.) ax=xx
      write(*,'(10x,a,f9.0,a\)') ' >> Xmax (def.',bx,') ? '
      read(*,'(f8.0)') xx
      if(xx.ne.0.) bx=xx
      write(*,'(/10x,a,f9.0,a\)') ' >> Ymin (def.',ay,') ? '
      read(*,'(f9.0)') xx
      if(xx.ne.0.) ay=xx
      write(*,'( 10x,a,f9.0,a\)') ' >> Ymax (def.',by,') ? '
      read(*,'(f9.0)') xx
      if(xx.ne.0.) by=xx
      bb=zt-tam/2
      write(*,'(/10x,a,f9.0,a\)') ' >> Zmin (def.',bb,') ? '
      read(*,'(f9.0)') xx
      if(xx.ne.0.) bb=xx
      tt=zt
      write(*,'( 10x,a,f9.0,a\)') ' >> Zmax (e.g.',tt,') ? '
      read(*,'(f9.0)') xx
      tt=xx
   37 ps=(by-ay+bx-ax)/2./60.
      write(*,'(/10x,a,f6.0,a\)') ' >> Grid step (def.',ps,') ? '
      read(*,*) xx
      if(xx.ne.0.) ps=xx
      write(*,'(/30x,a\)') ' >> Stratified background: YES(0)/no(1) ? '
      read(*,'(i1)') kst
      nx=(bx-ax)/ps+1
      ny=(by-ay)/ps+1
      nz=(tt-bb)/ps+1
      if(nx.gt.mt) then
      write(*,200) mt
  200 format(/6x,'*** Error: max. dimension ',i7,' is surpased !!')
      go to 37
      endif
      write(*,'(40x,3(3x,a,i3))') 'nx=',nx,'ny=',ny,'nz=',nz
      write(1,'(f9.0)') ax
      bx=ax+(nx-1)*ps
      write(1,'(f9.0)') bx
      write(1,'(f9.0)') ay
      by=ay+(ny-1)*ps
      write(1,'(f9.0)') by
      write(1,'(f9.0)') bb
      tt=bb+(nz-1)*ps
      write(1,'(f9.0)') tt
      write(1,'(i4)') nx,ny,nz
      do 42 k=1,nz
      zz=bb+ps*(k-1)
      ns=0
      tt=0 
      do 38 i=1,nb
      if(fq(i).eq.0) go to 38    
      c=zs(i)-zi(i)
      if(zz.gt.(zs(i)+c).or.zz.lt.(zi(i)-c)) go to 38
      ns=ns+1
      tt=tt+(dx(i)+zs(i)-zi(i))/2.
      jb(ns)=i
   38 continue
      if(ns.gt.0) tt=tt/ns 
      do 41 j=1,ny
      yy=ay+ps*(j-1)
      do 39 i=1,nx
      xx=ax+ps*(i-1)
      c=0
      if(ns.gt.0.and.kst.eq.1) 
     - call am(mc,xb,yb,dx,dy,zs,zi,fq,jb,ns,xx,yy,zz,tt,c)
      if(ns.gt.0.and.kst.eq.0) 
     - call am(mc,xb,yb,dx,dy,zs,zi,db,jb,ns,xx,yy,zz,tt,c)

c      write(2,'(4i9,f6.0,i6)') nint(xx),nint(yy),nint(zz),nint(c),tt,ns
c                             c=c*(1.+ps/abs(zz))
c                  if(c.lt.0) c=0
   39 ij(i)=c
   41 write(1,'(140i5)') (ij(i),i=1,nx)
      write(*,'(62x,i5,f9.0,i6)') k,zz,ns
   42 write(1,'(a)') ' '
      close(1)
      close(2)
      write(*,'(20x,a)') ' >> Press (Enter)'
      read(*,'(i1)') i
      go to 10
      endif

c     Option 14: map of discontinuity surface

      if(iop.eq.14) then
      call clearscreen($clearscreen)
      open(1,file='out.txt')
      nada=setfont('t''Times New Roman''h18w9e')
      write(texto,'(a40,i1)') 'Depth (m) of discontinuity interface '
      call moveto(120,8,xy)
      call outgtext(texto)
      nada=setcolorrgb(#000050)
      nada=setfont('t''Arial''h12w6i')
      write(texto,'(a)') 'GROWTH 2020'
      call moveto(18,8,xy)
      call outgtext(texto)
      write(texto,'(a9)') hoy
      nada=setfont('t''Arial''h13w7i')
      call moveto(506,8,xy)
      call outgtext(texto)
      nada=setcolorrgb(#000000)
      nada=setfont('t''Courier''h16w8')
      call moveto(24,80,xy)
      write(texto,'(a30)') '.'
      call outgtext(texto)
      call settextposition(25,80,curpos)
      write(*,'(a\)') ' >> Density for discont. ? '
      read(*,*) dl
c      call settextposition(27,80,curpos)
c      write(*,'(a\)') ' >> Include topogr. Y(1)/N(0) ? '
c      read(*,'(i1)') kt
      call settextposition(20,30,curpos)
      write(*,'(a)') ' Please, wait ...'
      open(1,file='out.txt')
      ng=0
      dn=0  
      do 44 i=1,nb
      if(sen(i).lt.40) go to 44
      d1=dbm(i)
      zz=(zs(i)+zi(i))/2.
      if(d1.lt.dl) go to 44
      xx=dx(i)/1.8
      yy=dy(i)/1.8
        k=0
      do 43 j=1,nb
      if(sen(j).lt.40) go to 43
      if(abs(zi(j)-zs(i)).gt.2) go to 43
      if(abs(xb(j)-xb(i)).gt.xx) go to 43
      if(abs(yb(j)-yb(i)).gt.yy) go to 43
      d2=dbm(j)
      zz=(zs(j)+zi(j))/2.
        k=1
      if(d2.ge.dl) go to 43
       ng=ng+1
       ii(ng)=j
       dn=dn+zi(j)
       write(1,'(4i8)') xb(j),yb(j),zi(j),1
   43 continue
c       if(k.eq.0) then
c         ng=ng+1
c         ii(ng)=i
c         dn=dn+zs(i)
c         write(1,'(4i8)') xb(i),yb(i),zs(i),2
c       endif  
   44 continue
      close(1)
      if(ng.gt.0) dn=dn/ng
      call settextposition(30,84,curpos)
      write(*,'(a,f7.0,a)') 'Mean depth =',dn*fe+czm,' m (see out.txt)'
      call grid(mc,ng,ax,ay,tam,tic,xb,yb,zi,dx,dy,ii,nc,col,fe,
     -mb,np,xp,yp,jp)
      close(1)
      endif

c     Option 15:  sections of the model sensitibity

      if(iop.eq.15) then
      is=0
      do 45 i=1,nb
   45 cb(i)=abs(eb(i))
      call profiles(iop,is,fe,tb,mc,nb,xb,yb,zi,zs,dx,dy,eb,
     -cb,lr,ax,ay,tam,nc,col,zt,czm,tic,mb,np,xp,yp,jp,ms,n,x,y,z)
      endif

c     Option 16:  perpective view

      if(iop.eq.16) then
      do 47 i=1,nb
      c=db(i)
      zz=(zs(i)+zi(i))/2.
      cb(i)=c
   47 continue
      call perspec(ax,ay,tam,fe,mc,nb,xb,yb,dx,dy,zi,zs,ii,cb,eb,
     -ms,n,x,y,z, mb,np,xp,yp,jp,tic)
      endif

c     Option 17:  delete isolated cells

      if(iop.eq.17) then
      call clean(mc,nb,xb,yb,dx,dy,zi,zs,db,cb)
      close(1)
      endif

c     Option 18: Local/Regional cells

      if(iop.eq.18) then
      call clearscreen($clearscreen)
      call settextposition(20,21,curpos)
      write(*,'(a\)') ' >> Limiting value (def. 60.) ?  '
      read(*,'(f4.0)') pn
      if(pn.eq.0) pn=60.
      p=pn/fe*1.e6 
      write(*,'(/20x,a,f4.0,a,e7.2,a)') 
     -' >> Wait (',pn,') (sensit. ',p,' Tm/uGal ) .....'
      nl=0
      nr=0
      do 50 j=1,n
   50 gr(j)=0
      do 48 i=1,nb
      lr(i)=0     
      tt=(zs(i)+zi(i))/2. 
      dd=0 
      do 53 j=1,n
      xx=x(j)-xb(i)
      yy=y(j)-yb(i)
      zz=z(j)- tt 
      d2=xx*xx+yy*yy+zz*zz
   53 dd=dd+zz*zz/d2/d2/d2              !atrac normal
      dd=sqrt(dd/n) *1000.*6.672e-3     ! atrac media uGal por tonelada normal
      dd=1.e-6/dd                       ! sensit Ton por uGal
      if(dd.gt.pn) lr(i)=1              ! regional
      if(lr(i).eq.0) nl=nl+1 
      if(lr(i).eq.0) go to 48
      nr=nr+1 
      c=fe*dx(i)*dy(i)*(zs(i)-zi(i))*db(i)*6.672e-3
      do 46 j=1,n 
      xx=x(j)-xb(i)
      yy=y(j)-yb(i)
      zz=z(j)- tt 
      d2=xx*xx+yy*yy+zz*zz
   46 gr(j)=gr(j)+c*zz/sqrt(d2*d2*d2)
   48 continue

      write(*,'(/23x,a,3i7)') 
     -' Number of Total, Local and Regional cells:',nb,nl,nr
      write(*,'(//20x,a)') ' >> Press (Enter)'
      read(*,'(i1)') i
      close(1)
      endif

c     Option 19: map of pressure and loading

      if(iop.eq.19) then
      call clearscreen($clearscreen)
      open(1,file='out.txt')
      call settextposition(20,20,curpos)
      write(*,'(a\)') 'Profundidad maxima masa pesante (ej.-5000 m) ? '
      read(*,*) pma
      call settextposition(30,20,curpos)
      write(*,'(a)') 'Pressure (Ton/m2)....Please, wait ...'
      pi=3.141592653589793D0
      dmu=10.
      sig=0.25
      u1=9.8066/4./pi/dmu/1.e7 
      np=100
      pp=tam/np
      do 52 i=1,np
      xx=ax+(i-1)*pp
      do 52 j=1,np
      yy=ay+(j-1)*pp
      x1=xx-pp/2.
      x2=xx+pp/2.
      y1=yy-pp/2.
      y2=yy+pp/2.
      pr=0
      cdz=0
      do 51 k=1,nb
      if(db(k).eq.0) go to 51
      if(zi(k).lt.pma) go to 51
      zz=zs(k)-zi(k)
      vo=zz*dx(k)*dy(k)
      xr=xb(k)-xx
      yr=yb(k)-yy
      zr=(zs(k)+zi(k))/2.-0.
      zz=zr*zr
      d2=xr*xr+yr*yr+zz
      d=sqrt(d2)
      cdz=cdz-u1/d*(2.*(1.-sig)+zz/d2)*vo*db(k)
      xi=xb(k)-dx(k)/2. 
      if(x1.gt.xi) xi=x1
      xf=xb(k)+dx(k)/2. 
      if(x2.lt.xf) xf=x2
      if(xi.ge.xf) go to 51
      yi=yb(k)-dy(k)/2. 
      if(y1.gt.yi) yi=y1
      yf=yb(k)+dy(k)/2. 
      if(y2.lt.yf) yf=y2
      if(yi.ge.yf) go to 51
      pr=pr+db(k)*zz*(xf-xi)*(yf-yi)
   51 continue
   52 write(1,'(2f9.0,f11.0,f9.0)') xx,yy,pr/pp/pp/1000.,cdz
      close(1)
      endif

c    Option 20: Writing cells
      
      if(iop.eq.20) then
      call clearscreen($clearscreen)
      write(*,'(///4x,a\)') '>> Local(1), Regional (2), Both (0) ? '
      read(*,'(i1)') jr
      write(*,'(//4x,a\)') '>> Filled(1), Total (0) ? '
      read(*,'(i1)') k
      np=0
      do 57 i=1,nb
      if(db(i).eq.0.and.k.eq.1) go to 57        
      if(lr(i).eq.0.and.jr.eq.2) go to 57 
      if(lr(i).eq.1.and.jr.eq.1) go to 57 
      np=np+1 
   57 continue
      open(1,file='out.txt')
      write(1,*) np
      do 54 i=1,nb
      if(db(i).eq.0.and.k.eq.1) go to 54        
      if(lr(i).eq.0.and.jr.eq.2) go to 54 
      if(lr(i).eq.1.and.jr.eq.1) go to 54 
      zz=(zi(i)+zs(i))/2.       
      write(1,'(i7,i8,i7,3i6,3i8)') xb(i),yb(i),nint(zz),
     -dx(i),dy(i),zs(i)-zi(i),db(i),dbm(i),sen(i)
   54 continue 
      close(1)
      endif

c
      go to 10

C  End and some messages for error or warning

      stop
   97 write(*,'(//20x,a//)') '  --> Error reading "fil.dat" !! '
      stop
   96 write(*,'(//20x,a//)') '  --> Error reading "mod.dat" !! '
      stop
   95 write(*,'(//20x,a//)') '  --> Error reading "gra.dat" !! '
   99 STOP
      END
c********************************************************
c  Determining a colour k corresponding to a value dn
c********************************************************
      subroutine acol(nc,ad,eq,dn,k)
      c=eq/2
      do 1 i=1,nc
      a=abs(ad+(i-1)*eq-dn)
      if(a.gt.c) go to 1
      k=i
      go to 2
    1 continue
      if(dn.le.(ad-0.5*eq)) k=0
      if(dn.gt.(ad+(nc-0.5)*eq)) k=0
    2 return
      end
c*********************************************************
c   Coloured map of a planar data distribution x,y,g
c*********************************************************
      subroutine cmap(ms,ax,bx,ay,by,nc,n,x,y,g,ig,tic,fe,col,
     -dr,mb,np,xp,yp,jp,fname,hoy,iop,fmod)
      use IFQWIN !msflib
      TYPE (xycoord) xy
      TYPE (wxycoord) wxy
      TYPE (rccoord) curpos
      parameter (mr=60)
      character*120 texto,fname,tex
      character*20 fmod
      character*4 un,ung,unm
      character*9 hoy
      integer nada,x(ms),y(ms), jp(mb), ig(ms),
     - col(30), ij(30), xp(mb),yp(mb)
      dimension g(ms),cv(mr),w(mr)
      data ung/'uGal'/,unm/'mts '/
      tam=bx-ax
      if((by-ay).gt.tam) tam=by-ay
      ag=9.d9
      bg=-ag
      gm=0.
      sp=0
         gx=0
         gy=0
         lx=0
         ly=0
      do 1 i=1,n
      if(g(i).lt.ag) ag=g(i)
      if(g(i).gt.bg) bg=g(i)
      c=1./ig(i)/ig(i)
      sp=sp+c
         do 15 j=1,n     
         t=g(i)-g(j)
         tx=x(i)-x(j)
         if(abs(tx).gt.12000*fe) then
         lx=lx+1   
         gx=t/tx+gx  
         endif 
         ty=y(i)-y(j)
         if(abs(ty).gt.12000*fe) then
         ly=ly+1   
         gy=t/ty+gy  
         endif 
   15    continue
    1 gm=gm+g(i)*c
      gm=gm/sp
         gx=gx/lx*1000.
         gy=gy/ly*1000.
      sg=0.
      do 13 i=1,n
      c=g(i)-gm
   13 sg=sg+c*c/ig(i)/ig(i)
      c=0
      if(sp.gt.0) c=sg/sp
      sg=0
      if(c.gt.0.) sg=sqrt(c)
      call clearscreen($clearscreen)
      nada=setcolorrgb(#000050)
      nada=setfont('t''Arial''h12w6i')
      write(texto,'(a)') 'GROWTH 2020'
      call moveto(28,8,xy)
      call outgtext(texto)
      write(texto,'(a9)') hoy
      nada=setfont('t''Arial''h13w7i')
      call moveto(480,8,xy)
      call outgtext(texto)
      nada=setcolorrgb(#000000)
      nada=setfont('t''Times New Roman''h18w9e')
      write(texto,'(a30)')  fname
      call moveto(170,8,xy)
      call outgtext(texto)

c   Title and colour scale

      if(tic.gt.0.) un=ung
      if(tic.lt.0.) un=unm
      nada=setfont('t''Courier''h16w8')
      nada=setcolorrgb(#000000)
      write(texto,'(a,i6)') ' Number=',n
      call moveto(770,330,xy)
      call outgtext(texto)
      write(texto,'(a,f8.0,1x,a4)') ' Mean=',gm,un
      call moveto(770,345,xy)
      call outgtext(texto)
      write(texto,'(a,f8.0,1x,a4)') ' S.D.=',sg,un
      call moveto(770,360,xy)
      call outgtext(texto)
      write(texto,'(a)') un
      ly=30
      call moveto(650,ly,xy)
      call outgtext(texto)
      eq=(bg-ag)/(nc-1)
      do 2 k=1,nc
      c=ag+(k-1)*eq
      if(tic.gt.0.) write(texto,'(f8.0)') c
      if(tic.lt.0.) write(texto,'(f8.2)') c
      call moveto(583,ly+k*16,xy)
      nada=setcolorrgb(#000000)
      call outgtext(texto)
      nada=setcolorrgb(col(k))
      c=16*k+ly
    2 nada=rectangle($gfillinterior,650,c,680,c+16)
      nada=setcolorrgb(#000000)
      nada=rectangle($gborder,650,ly+16,680,ly+16+nc*16)
      lx=20
      ly=40
      mx=530

c   Setting limits. Drawing axes

      call settextposition(22,74,curpos)
      if(tic.lt.0.) tic=-tic
      nada=setfont('t''Courier''h16w8')
      write(texto,'(a,f6.0,a)')  'Axis tick interval: ',tic,' m'
      call moveto(200,580,xy)
      call outgtext(texto)
      write(texto,'(a,2f9.0)') 'X limits (m):',ax,bx
      call moveto(200,600,xy)
      call outgtext(texto)
      write(texto,'(a,2f9.0)') 'Y limits (m):',ay,by
      call moveto(200,620,xy)
      call outgtext(texto)
      write(texto,'(2f8.0,a9)')  gx,gy,'ugal/km'
      call moveto(200,650,xy)
      call outgtext(texto)
      call setviewport(lx,ly,lx+mx,ly+mx)
      nada=setcolorrgb(#d0d0d0)
      nada=rectangle($gfillinterior,0,0,mx,mx)
      nada=setcolorrgb(#000000)
      nada=rectangle($gborder,0,0,mx,mx)
      nada=setwindow(.true.,ax,ay,ax+tam,ay+tam)
      c=tic/5
      do 3 i=1,60
      call moveto_w(dble(ax),dble(ay+i*tic),wxy)
      nada=lineto_w(ax+c,ay+i*tic)
      call moveto_w(dble(ax+tam),dble(ay+i*tic),wxy)
      nada=lineto_w(ax+tam-c,ay+i*tic)
      call moveto_w(dble(ax+i*tic),dble(ay),wxy)
      nada=lineto_w(ax+i*tic,ay+c)
      call moveto_w(dble(ax+i*tic),dble(ay+tam),wxy)
    3 nada=lineto_w(ax+i*tic,ay+tam-c)

c   Drawing the coloured distribution for data

      d=tam/1400.
      mx=tam/100.
      if(n.gt.10) mx=tam/60
      if(n.gt.50) mx=tam/100
      if(n.gt.300) mx=tam/150
      if(n.gt.2000) mx=tam/300
      if(n.gt.6000) mx=tam/500

      if(mx.eq.0) mx=1
      do 5 i=1,n
      do 4 j=1,nc
      c=ag+(j-1)*eq
      if(g(i).gt.(c-0.5*eq).and.g(i).le.(c+0.5*eq)) k=j
    4 continue
      nada=setcolorrgb(col(k))
      nada=ellipse_w($gfillinterior,x(i)-mx,y(i)-mx,x(i)+mx,y(i)+mx)
      nada=setcolorrgb(#000000)
    5 nada=ellipse_w($gfillinterior,x(i)-d,y(i)-d,x(i)+d,y(i)+d)

c   Drawing the base map

      call cont(mb,np,xp,yp,jp)

      if(iop.eq.7) then
      open(1,file=fmod)
      read(1,*) nb
      do 21 i=1,nb
   21 read(1,*)
      call setviewport(0,0,1000,720)
      nada=setfont('t''Arial''h16w7')
      read(1,'(/////////a80)') tex
      write(texto,'(a)') tex
      call moveto(580,500,xy)
      call outgtext(texto)
      read(1,'(a45)') tex
      write(texto,'(a)') tex
      call moveto(580,520,xy)
      call outgtext(texto)
      close(1)
      endif

c  Histogram

      c=bg-ag
      if(c.eq.0.) go to 99
      eq=c/4
      call setviewport(0,0,1000,720)
      lx=750
      write(texto,'(a)') 'Num.benchmarks'
      nada=setfont('t''Courier''h14w7')
      call moveto(lx-50,75,xy)
      call outgtext(texto)
      write(texto,'(5i8)')
     -nint(ag),nint(ag+eq),nint(ag+2*eq),nint(ag+3*eq),nint(bg)
      call moveto(lx-50,300,xy)
      call outgtext(texto)
      write(texto,'(a)') un
      call moveto(lx+200,310,xy)
      call outgtext(texto)
      mx=15
      eq=c/mx
      k=0
      do 7 j=1,mx
      ij(j)=0
      c=ag+eq*j
      d=c-eq
      do 8 i=1,n
    8 if(g(i).ge.d.and.g(i).lt.c) ij(j)=ij(j)+1
    7 if(ij(j).gt.k) k=ij(j)
      i=k/12.
      t=(i+1)*12
      k=t/12
      do 6 i=1,13
      j=(13-i)*k
      write(texto,'(i5)') j
      call moveto(lx-38,68+17*i,xy)
    6 call outgtext(texto)
      call setviewport(lx,92,lx+220,298)
      nada=setwindow(.true.,ag,0,bg,t)
      nada=setcolorrgb(#e0e0e0)
      nada=rectangle_w($gfillinterior,ag,0,bg,t)
      nada=setcolorrgb(#2060f0)
      do 9 j=1,mx
      if(ij(j).eq.0) go to 9
      c=ag+eq*j
      d=c-eq
      nada=rectangle_w($gfillinterior,d+eq/10,0,c-eq/10,ij(j))
    9 continue
      nada=setcolorrgb(#f00000)
      call moveto_w(dble(gm),0.d0,wxy)
      nada=lineto_w(gm,t)
      nada=setcolorrgb(#00f000)
      call moveto_w(0.d0,0.d0,wxy)
      nada=lineto_w(0.,t)
      nada=setcolorrgb(#000000)
      nada=rectangle_w($gborder,ag,0,bg,t)
      k=t/12
      do 11 i=1,13
      j=(13-i)*k
      call moveto_w(dble(ag),dble(j),wxy)
   11 nada=lineto_w(ag+eq/10,j)
      eq=(bg-ag)/4
      do 12 i=1,5
      call moveto_w(dble(ag+i*eq),dble(0),wxy)
   12 nada=lineto_w(ag+i*eq,k/5)
      if(iop.eq.7.or.iop.eq.3) go to 99

c   Autocorrelation analysis

      call cov2(ms,n,x,y,g,ig,mr,dr,cv,w,a,b,c,sv2,em)
      lx=590
      ly=430
      tx=410
      ty=230
      call setviewport(0,0,1000,720)
      nada=setcolorrgb(#e0e0e0)
      nada=rectangle($gfillinterior,lx,ly,lx+tx,ly+ty)
      nada=setcolorrgb(#000000)
      nada=rectangle($gborder,lx,ly,lx+tx,ly+ty)
      nada=setfont('t''Courier''h14w7')
      call setviewport(lx-30,ly,lx,ly+ty)       ! eje vertical
      nada=setwindow(.true.,0,-0.5,30,1.)
      px=10*tic
      do 41 i=1,10
      yy=-0.6+(i-1)*0.2
      call moveto_w(dble(25),dble(yy),wxy)
      nada=setcolorrgb(#000000)
      nada=lineto_w(30,yy)
      call moveto_w(0.d0,dble(yy+0.02),wxy)
      write(texto,'(f3.1)') yy
      nada=setcolorrgb(#f00000)
   41 call outgtext(texto)
      call setviewport(lx,ly+ty,lx+tx,ly+ty+30)  ! eje horizontal
      nada=setwindow(.true.,0.,0.,px,30)
      nada=setcolorrgb(#0000f0)
      do 42 i=1,12
      j=(i-1)*tic
      nada=setcolorrgb(#000000)
      call moveto_w(dble(j),25.d0,wxy)
      nada=lineto_w(j,30.)
      call moveto_w(dble(j-0.48*tic),25.d0,wxy)
      write(texto,'(i5)') j
      nada=setcolorrgb(#f00000)
   42 call outgtext(texto)
      write(texto,'(a)') 'Mutual distance (m)'
      call moveto_w(5.d0*tic,12.d0,wxy)
      call outgtext(texto)
      call setviewport(lx,ly,lx+tx,ly+ty)
      nada=setwindow(.true.,0,-0.5,px,1.)
      call moveto_w(0.d0,0.d0,wxy)
      nada=lineto_w(px,0.)
      nada=setcolorrgb(#009000)
      px=tam/200.
      py=0.02
      do 43 i=1,mr
      xx=i-0.5
      xx=xx*dr! tic
      fk=fbx(xx,a,b,c)
      nada=ellipse_w($gfillinterior,xx-px,cv(i)-py,xx+px,cv(i)+py)
   43 continue
      nada=setcolorrgb(#0000f0)
      fk=fbx(0.,a,b,c)
      call moveto_w(0.d0,dble(fk),wxy)
      do 44 i=1,2*mr
      xx=(i-0.5)/2
      fk=fbx(xx,a,b,c)
   44 nada=lineto_w(xx*dr,fk)
      nada=setcolorrgb(#009000)
      nada=ellipse($gfillinterior,96,12,100,15)
      nada=setcolorrgb(#0000f0)
      call moveto(96,22,xy)
      nada=lineto(99.,29.)
      nada=setcolorrgb(#000000)
      write(texto,'(a)') 'Empirical autocorrel.'
      call moveto(103,5,xy)
      call outgtext(texto)
      write(texto,'(a)') 'Analytic.autocovar.func.'
      call moveto(103,18,xy)
      call outgtext(texto)
      write(texto,'(a)') 'y=a Jo(cx) e(-bx)'
      call moveto(115,30,xy)
      call outgtext(texto)
      write(texto,'(2(a,f4.2),2(a,f5.2))')
     -'a=',a,' (C1=',cv(1),')  b=',b,'  c=',c
      call moveto(115,42,xy)
      call outgtext(texto)
      write(texto,'(a,f6.0,a)')
     -'correlat.step =',dr,' m'
      call moveto(115,54,xy)
      call outgtext(texto)
      write(texto,'(a,f5.2)') 'rms.resid=',em
      call moveto(115,66,xy)
      call outgtext(texto)

   99 call setviewport(0,0,600,750)
      write(texto,'(a)') '>> Press (Enter)'
      call moveto(400,700,xy)
      nada=setcolorrgb(#000000)
      nada=setfont('t''Courier''h16w8')
      call outgtext(texto)
      return
      end
c****************************************************
c     Drawing a base map
c****************************************************
      subroutine cont(mb,np,xp,yp,jp)
      use msflib
      TYPE (rccoord) curpos
      TYPE (wxycoord) wxy
      integer nada,jp(mb),xp(mb),yp(mb)
      nada=setcolorrgb(#000000)
      do 1 i=1,np
      if(jp(i).eq.1) call moveto_w(dble(xp(i)),dble(yp(i)),wxy)
      if(jp(i).eq.0) nada=lineto_w(xp(i),yp(i))
    1 continue
      call settextposition(26,11,curpos)
      return
      end
c****************************************************
c    Combining basic colours
c****************************************************
      integer*4 function irva(r,v,a)
      integer r,v,a
      irva=ishl( ishl( a, 8 ) .or. v, 8 ) .or. r
      return
      end
c*************************************************
c    Determining a topographical profile
c*************************************************
      subroutine pertop(ms,x,y,z,n,xx,yy,zz,pp,sp,zt)
      integer x(ms),y(ms),z(ms)
      fr=500.   ! 5.
      zz=0.
      sp=0.
      np=0
      do 1 i=1,n
      xr=abs(x(i)-xx)
      if(xr.gt.(fr*pp)) go to 1
      yr=abs(y(i)-yy)
      if(yr.gt.(fr*pp)) go to 1
      dd= xr*xr+yr*yr+1.
      pe=1./dd
      zz=zz+pe*z(i)
      sp=sp+pe
      np=np+1
    1 continue
      if(sp.ne.0.) zz=zz/sp
      if(sp.eq.0.) zz=zt
      sp=np
      return
      end
c****************************************************
c   Interpolating and smoothing from a 3-D grid
c****************************************************
      subroutine am(mc,xb,yb,dx,dy,zs,zi,db,jb,ns,xx,yy,zz,tm,dn)
      integer zs(mc),zi(mc),jb(mc),xb(mc),yb(mc),db(mc),dx(mc),dy(mc)
      ar=8*tm*tm*tm
      dn=0.
      sp=0.
      de=0.
      x2=xx+tm
      x1=xx-tm
      y2=yy+tm
      y1=yy-tm
      z2=zz+tm
      z1=zz-tm
      do 1 i=1,ns
      j=jb(i)
    2 c=dx(j)/2.
      bx=xb(j)+c
      if(bx.le.x1) go to 1
      if(x2.lt.bx) bx=x2
      ax=xb(j)-c
      if(ax.ge.x2) go to 1
      if(x1.gt.ax) ax=x1
      c=dy(j)/2.
      by=yb(j)+c
      if(by.le.y1) go to 1
      if(y2.lt.by) by=y2
      ay=yb(j)-c
      if(ay.ge.y2) go to 1
      if(y1.gt.ay) ay=y1
      bz=zs(j)
      if(bz.le.z1) go to 1
      if(z2.lt.bz) bz=z2
      az=zi(j)
      if(az.ge.z2) go to 1
      if(z1.gt.az) az=z1
      c=(bx-ax)*(by-ay)*(bz-az)
      de=de+c*db(j)
      sp=sp+c
    1 continue
      dn=de/ar 
      if(dn.eq.0.and.sp.gt.0.) dn=0.01
   99 return
      end
c***********************************************************************
c     Autocorrelation analysis
c***********************************************************************
      subroutine cov2(m,n,x,y,v,iev,mr,dr, cv,w,a,b,c,sv2,em)
      real*4 cv(mr),w(mr),v(m)
      integer iev(m),x(m),y(m)
      dr2=dr*dr
      do 1 k=1,mr
      w(k)=0.
    1 cv(k)=0.
      vm=0.
      sp=0
      do 5 i=1,n
      p=1./iev(i)/iev(i)
      sp=sp+p
    5 vm=vm+v(i)*p
      vm=vm/sp
      sv2=0.
      do 6 i=1,n
      vi=v(i)-vm
      p=1./iev(i)/iev(i)
    6 sv2=sv2+vi*vi*p
      sv2=sv2/sp
      do 2 i=1,n-1
      vi=v(i)-vm
      do 3 j=i+1,n
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      c=xx*xx+yy*yy
      k=nint(c/dr2+0.5)     
      if(k.gt.mr) go to 3
      p=1./iev(i)/iev(j)
      cv(k)=cv(k)+(v(j)-vm)*vi*p
      w(k)=w(k)+p
    3 continue
    2 continue
      do 4 k=1,mr
      if(w(k).eq.0.) cv(k)=0
      if(w(k).gt.0.) cv(k)=cv(k)/w(k)/sv2
    4 continue  
      call busca(mr,cv,w,1,a,b,c,em)
   99 return
      end
c**************************************************************
c     Bessel fuction
c**************************************************************
      FUNCTION BESSJ0(X)
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     *   -.2073370639D-5,.2093887211D-6/,
     *Q1,Q2,Q3,Q4,Q5/-.1562499995D-1,.1430488765D-3,-.6911147651D-5,
     *    .7621095161D-6,-.934945152D-7/,
     *R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,651619640.7D0,
     *    -11214424.18D0,77392.33017D0,-184.9052456D0/,
     *S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,
     *    9494680.718D0,59272.64853D0,267.8532712D0,1.D0/
      IF(ABS(X).LT.8.)THEN
        y=x*x
        BESSJ0=(R1+y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        ax=abs(x)
        z=8./ax
        y=z*z
        xx=ax-.785398164
        BESSJ0=SQRT(.636619772/ax)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*Q5)))))
      ENDIF
      RETURN
      END
c**************************************************************
c     Function  expon-Bessel
c**************************************************************
      function fbx(x,a,b,c)
      fbx=0.
      t=x/10.
      if(b*t.le.170.) fbx=a*exp(-b*t)*BESSJ0(c*t)
      return
      end
c***********************************************************************
c     Fit of an analytical covariance function
c***********************************************************************
      subroutine busca(mr,cv,we,kf,a,b,c,er)
      dimension cv(mr),we(mr)
      er=9.d9
      a=0.5
      b=2.0
      c=6.0
      da=0.1
      db=0.4
      dc=1.
    5 do 1 i=1,11
      ae=a+(i-6)*da
      if(ae.le.0.01.or.ae.gt.0.99) go to 1
      do 2 j=1,11
      be=b+(j-6)*db
      if(be.lt.0.) go to 2
      do 3 k=1,11
      ce=c+(k-6)*dc
      if(ce.le.0.) go to 3
      dc0=2.405/ce*50.
      ere=0.
      sp=0.
      do 4 l=1,mr
      fk=0.
      r=l-0.5
      if(r.gt.dc0) go to 6
      if(kf.eq.1) fk=fbx(r,ae,be,ce)
      d=cv(l)-fk
      ere=ere+we(l)*d*d
      sp=sp+we(l)
    4 continue
    6 ere=ere/sp
      if(ere.ge.er) go to 3
        er=ere
        a1=ae
        b1=be
        c1=ce
    3 continue
    2 continue
    1 continue
      dc=dc/3
      db=db/3
      da=da/3
      a=a1
      b=b1
      c=c1
      if(er.gt.0.) er=sqrt(er)
      if(da.gt.0.001) go to 5
   99 return
      end
c***********************************************************************
c      Draw a title
c***********************************************************************
      subroutine titulo(hoy)
      use IFQWIN !msflib
      TYPE (xycoord) xy
      TYPE (rccoord) curpos
      character*120 texto
      character*9 hoy
      nada=setcolorrgb(#e0e0e0)
      nada=rectangle($gfillinterior,100,30,600,180)
      nada=setcolorrgb(#000000)
      nada=rectangle($gborder,100,30,600,180)
      nada=initializefonts()
      nada=setfont('t''Arial''h20w10e')
      write(texto,'(a)') 'GROWTH  Gravity Inversion'
      call moveto(200,70,xy)
      call outgtext(texto)
      write(texto,'(a)') 'VIEW:  user interface'
      nada=setfont('t''Times New Roman''h18w9e')
      call moveto(230,110,xy)
      call outgtext(texto)
      write(texto,'(a)') 'Camacho 2020'
      nada=setfont('t''Times New Roman''h12w6')
      call moveto(525,168,xy)
      nada=setcolorrgb(#909090)
      call outgtext(texto)
      write(texto,'(a)') hoy
      nada=setcolorrgb(#000000)
      nada=setfont('t''Courier''h16w8')
      call moveto(300,150,xy)
      call outgtext(texto)
      call settextposition(15,0,curpos)
      return
      end
c********************************************************
c     Stop when number of data surpases dimension
c********************************************************
      subroutine out(m)
      write(*,200) m
  200 format(/6x,'*** Error: max. dimension ',i7,' is surpased !!')
      stop
      end
c*********************************************************
c     Axis
c*********************************************************
      subroutine ejes(mx,my,tt,bb,c,a1,czm,fe)
      use IFQWIN !msflib
      TYPE (wxycoord) wxy
      nada=setcolorrgb(#000000)                      ! bordes y ticks
      nada=rectangle($gborder,0,0,mx,my)
      yy=(tt-bb)/150.
      do 1 i=1,30
      xx=-1.25+i*0.25
      call moveto_w(dble(xx),dble(tt),wxy)
      nada=lineto_w(xx,tt-yy)
      call moveto_w(dble(xx),dble(bb),wxy)
      nada=lineto_w(xx,bb+yy)
      zz=(a1-c*(i-1)-czm)/fe
      if(zz.lt.bb) go to 1
      call moveto_w(dble(-1.),dble(zz),wxy)
      nada=lineto_w(-0.98,zz)
    1 continue
      call moveto_w(dble(0.),dble(tt),wxy)
      nada=lineto_w(0,bb)
      return
      end
c***********************************************************************
c     Dialog interface 1                                               *
c***********************************************************************
      subroutine dial1(obs,mod,fil,map)
      use dialogm
      include 'View.fd'
      logical retlog,nada
      character*20 obs,mod,fil,map
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_View, dlg)
      retlog = DlgSet(dlg,IDC_EDIT_e1 ,obs )
      retlog = DlgSet(dlg,IDC_EDIT_e2 ,mod )
      retlog = DlgSet(dlg,IDC_EDIT_e3 ,fil )
      retlog = DlgSet(dlg,IDC_EDIT_e4 ,map )
      retlog = DlgSetSub(dlg,IDCANCEL,fin)
      retlog = DlgSet(dlg,IDC_EDIT_e6 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t6 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e7 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t7 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e8 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t8 ,.false.,dlg_enable)

      retint = DlgModal( dlg )
                                                   ! input parameters
      nada = DlgGet( dlg, IDC_EDIT_e1, obs )
      nada = DlgGet( dlg, IDC_EDIT_e2, mod )
      nada = DlgGet( dlg, IDC_EDIT_e3, fil )
      nada = DlgGet( dlg, IDC_EDIT_e4, map )

      call DlgUninit(dlg)
      return
      end
c***********************************************************************
c     Dialog interface 2                                               *
c***********************************************************************
      subroutine dial2(obs,mod,fil,map,n,nb,np)
      use dialogm
      include 'View.fd'
      logical retlog
      character*20 obs,mod,fil,map,texto
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_View, dlg)
      retlog = DlgSet(dlg,IDC_EDIT_e1 ,obs )
      retlog = DlgSet(dlg,IDC_EDIT_e2 ,mod )
      retlog = DlgSet(dlg,IDC_EDIT_e3 ,fil )
      retlog = DlgSet(dlg,IDC_EDIT_e4 ,map )
      retlog = DlgSetSub(dlg,IDCANCEL,fin)
      write(texto,'(i5)') n
      retlog = DlgSet(dlg,IDC_EDIT_e6 ,texto )
      write(texto,'(i6)') nb
      retlog = DlgSet(dlg,IDC_EDIT_e7 ,texto )
      write(texto,'(i5)') np
      retlog = DlgSet(dlg,IDC_EDIT_e8 ,texto )
      retlog = DlgSet(dlg,IDC_EDIT_e1 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t1 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e2 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t2 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e3 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t3 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e4 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t4 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_BOX_b1  ,.false.,dlg_enable)

      retint = DlgModal( dlg )
                                                   ! input parameters
      call DlgUninit(dlg)
      return
      end
c***********************************************************************
c     Dialog interface 3                                               *
c***********************************************************************
      subroutine dial3(nn,sr,igr,ibl,rl,rs)
      use dialogm
      include 'View.fd'
      logical retlog,nada
      character*20 texto
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Interp, dlg)

      nn=200
      write(texto,'(i3)') nn
      retlog = DlgSet(dlg,IDC_EDIT3,texto)
      write(texto,'(f6.0)') sr
      retlog = DlgSet(dlg,IDC_EDIT4,texto)
      retlog = DlgSet(dlg,IDC_EDIT5,'60.')
      retlog = DlgSetSub(dlg,IDCANCEL,fin)

      retint = DlgModal( dlg )
                                                   ! input parameters
      nada = DlgGet( dlg, IDC_EDIT3, texto)
      read(texto,*) nn
      nada = DlgGet( dlg, IDC_EDIT4, texto)
      read(texto,*) sr
      nada = DlgGet( dlg, IDC_EDIT5, texto)
      read(texto,*) rl
      nada = DlgGet( dlg, IDC_EDIT6, texto)
      read(texto,*) rs
      igr=0
c      retlog=DlgGet( dlg, IDC_CHECK2, nada)
c      if(nada) igr=1
      ibl=0
c      retlog=DlgGet( dlg, IDC_CHECK3, nada)
c      if(nada) ibl=1

      call DlgUninit(dlg)
      return
      end
c***********************************************************************
c     Dialog interface 4                                               *
c***********************************************************************
      subroutine dial4(ad,eq,id0,fta)
      use dialogm
      include 'View.fd'
      logical retlog,nada
      character*20 texto
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Mapa, dlg)
      write(texto,'(f7.0)') ad
      retlog = DlgSet(dlg,IDC_EDIT1,texto)
      write(texto,'(f7.2)') eq
      retlog = DlgSet(dlg,IDC_EDIT2,texto)
      write(texto,'(f7.0)') ad
      retlog = DlgSet(dlg,IDC_EDIT3,texto)
      write(texto,'(f3.1)') 0.3
      retlog = DlgSet(dlg,IDC_EDIT4,texto)

      retint = DlgModal( dlg )
                                                   ! input parameters
      nada = DlgGet( dlg, IDC_EDIT1, texto)
      read(texto,*) ad
      nada = DlgGet( dlg, IDC_EDIT2, texto)
      read(texto,*) eq
      nada = DlgGet( dlg, IDC_EDIT3, texto)
      read(texto,*) id0
      nada = DlgGet( dlg, IDC_EDIT4, texto)
      read(texto,*) fta

      call DlgUninit(dlg)
      return
      end


c***********************************************************************
      subroutine fin()
      use dialogm
      include 'View.fd'
      stop
      end
c**********************************************************************
      subroutine grid(m,n,ax,ay,tam,tic,x,y,v,dx,dy,ii,nc,col,
     -fe,mb,np,xp,yp,jp)
      use msflib
      TYPE (xycoord) xy
      TYPE (wxycoord) wxy
      TYPE (rccoord) curpos
      real*8 xdp,ydp
      character*120 texto
      character*1 key
      integer nada,jp(mb),col(30),xp(mb),yp(mb),
     -x(m),y(m),v(m),dx(m),dy(m),ii(m)
      p=tam/n
      ag=9.d9
      bg=-ag
      do 1 i1=1,n
      i=ii(i1) 
      d=v(i)
      if(d.lt.ag) ag=d
      if(d.gt.bg) bg=d
    1 continue
      eq=(bg-ag)/(nc-1)
      ly=40
      do 2 k=1,nc
      c=ag+(k-1)*eq
      write(texto,'(f8.0)') c
      call moveto(583,ly+k*16,xy)
      nada=setcolorrgb(#000000)
      call outgtext(texto)
      nada=setcolorrgb(col(k))
      c=16*k+ly
    2 nada=rectangle($gfillinterior,650,c,680,c+16)
      continue
      nada=setcolorrgb(#000000)
      nada=rectangle($gborder,650,ly+16,680,ly+16+nc*16)
      lx=20
      mx=550
      call setviewport(lx,ly,lx+mx,ly+mx)
      nada=setcolorrgb(#d0d0d0)
      nada=rectangle($gfillinterior,0,0,mx,mx)
      nada=setcolorrgb(#000000)
      nada=rectangle($gborder,0,0,mx,mx)
      nada=setwindow(.true.,ax,ay,ax+tam,ay+tam)
      do 4 i1=1,n
      i=ii(i1)
      xc=x(i)
      yc=y(i)
      rx=dx(i)/3.
      ry=dy(i)/3. 
      if(ry.lt.rx) rx=ry
      d=v(i)
      if(d.ge.bg) k1=nc
      if(d.le.ag) k1=1
      do 3 k=1,nc
      c=ag+(k-1.5)*eq
      if(d.ge.c.and.d.le.c+eq) k1=k
    3 continue
      nada=setcolorrgb(col(k1))
      nada=rectangle_w($gfillinterior,xc-rx,yc-rx,
     -xc+rx,yc+rx)
    4 continue
      c=tic/5
      nada=setcolorrgb(#000000)
      do 5 i=1,60
      call moveto_w(dble(ax),dble(ay+i*tic),wxy)
      nada=lineto_w(ax+c,ay+i*tic)
      call moveto_w(dble(ax+tam),dble(ay+i*tic),wxy)
      nada=lineto_w(ax+tam-c,ay+i*tic)
      call moveto_w(dble(ax+i*tic),dble(ay),wxy)
      nada=lineto_w(ax+i*tic,ay+c)
      call moveto_w(dble(ax+i*tic),dble(ay+tam),wxy)
    5 nada=lineto_w(ax+i*tic,ay+tam-c)
      nada=rectangle($gborder,0,0,mx,mx)
      call cont(mb,np,xp,yp,jp)
      call settextposition(39,20,curpos)
      write(*,'(a,i6,a)') ' Axis tick interval:',nint(tic),' m'
      xdp=ax+tam/2.
      ydp=ay+tam/2.
      call settextposition(34,80,curpos)
      write(*,'(a)') ' >> Cursor: N, S, E, W, q(quit) '
      nada=setcolorrgb(#00f000)
   10 call settextposition(36,86,curpos)
      write(*,'(a,2i8,a)') ' (',nint(xdp),nint(ydp),')'
c     call getimage_w(xdp-p,ydp,xdp+p,ydp,c1)
c     call getimage_w(xdp,ydp+p,xdp,ydp-p,c2)
      call moveto_w(xdp-p,ydp,wxy)
      nada=lineto_w(xdp+p,ydp)
      call moveto_w(xdp,ydp+p,wxy)
      nada=lineto_w(xdp,ydp-p)
      key=GETCHARQQ()
c     call putimage_w(xdp-p,ydp,c1,$gpset)
c     call putimage_w(xdp,ydp+p,c2,$gpset)
      if(key.eq.'q') go to 99
      if(key.eq.'s') ydp=ydp-p*2
      if(key.eq.'n') ydp=ydp+p*2
      if(key.eq.'w') xdp=xdp-p*2
      if(key.eq.'e') xdp=xdp+p*2
      go to 10

   99 return
      end
c******************************************************
c     Drawing profiles of a 3D distribution
c******************************************************
      subroutine profiles(iop,is,fe,tb,mc,nb,xb,yb,zi,zs,dx,dy,
     -eb,cb,lr,ax,ay,tam,nc,col,zt,czm,tic,mb,np,xp,yp,jp,ms,n,x,y,z)
      use IFQWIN !msflib
      TYPE (wxycoord) wxy
      TYPE (xycoord) xy
      TYPE (rccoord) curpos

c  Dimensions

      parameter (mt=5000)
      character*120 texto,fname
      character*9 hoy
      integer zs(mc),zi(mc),ij(mt),jb(mc), col(30), xp(mb),yp(mb),
     -    x(ms),y(ms),z(ms), jp(mb),nada , 
     -    xb(mc),yb(mc),dx(mc),dy(mc),cb(mc),lr(mc),eb(mc)

   30 call clearscreen($clearscreen)
      nada=setcolorrgb(#000050)
      nada=setfont('t''Arial''h12w6i')
      write(texto,'(a)') 'GROWTH 2020'
      call moveto(6,1,xy)
      call outgtext(texto)
      call titulo(hoy)
      sr=nint(tb/10.)*5.*fe
      pp=tam/150.
      rl=60
      if(is.eq.1) call dial3(nn,sr,igr,ibl,rl,rs)
      sr=sr/fe
      ad=0.
      if(iop.eq.15) ad=9.d9
      bd=0.
      do 1 i=1,nb
      if(cb(i).eq.0) go to 1    
      dd=cb(i)
      zz=(zi(i)+zs(i))/2.
      if(eb(i).lt.1) go to 1  
      if(dd.lt.ad) ad=dd
      if(dd.gt.bd) bd=dd
    1 continue
          ad=ad-(bd-ad)/60.
          bd=bd+(bd-ad)/60
  	    if(iop.ge.30) ad=0
          if(iop.ge.30) bd=200
      eq=(bd-ad)/(nc-1.)
      call dial4(ad,eq,id0,rs)
      if(is.eq.1) then
       nz=nn/2.+nn/10.
       ps=tam/nn
       p=ps/2.
       nx=tam/ps
       ny=tam/ps
       is=is+igr
      endif
      call clearscreen($clearscreen)
      nada=setcolorrgb(#000050)
      nada=setfont('t''Arial''h16w8i')
      write(texto,'(a)') 'GROWTH 2020'
      call moveto(40,7,xy)
      call outgtext(texto)
      write(texto,'(a9)') hoy
      call moveto(680,8,xy)
      nada=setfont('t''Arial''h18w9i')
      call outgtext(texto)
      nada=setcolorrgb(#000000)
      if(iop.eq.9) write(texto,'(a)')
     -' Sections across the 3D contrast density model'
      if(iop.eq.15) write(texto,'(a)')
     -' Sections across the sensitivity for the 3D model'
      if(iop.ge.30) write(texto,'(a)')
     -' Sections across the sensitivity for the domain'
      call moveto(280,6,xy)
      nada=setfont('t''Times New Roman''h18w9e')
      call outgtext(texto)
      nada=setfont('t''Courier''h18w8')
      call settextposition(4,1,curpos)
      if(iop.eq.9) write(texto,'(a)') ' Dens.contrast'
      if(iop.eq.15) write(texto,'(a)') ' Cell Sensitivity'
      if(iop.ge.30) write(texto,'(a)') 'Domain Sensitivity'
      call moveto(30,40,xy)
      call outgtext(texto)
      if(iop.eq.9) write(texto,'(4x,a)') 'kg/m3'
      if(iop.eq.15) write(texto,'(5x,a)') 'uGal'
      if(iop.eq.30) write(texto,'(a)') ' 10**5 Tm/uGal'
      if(iop.eq.31) write(texto,'(3x,a)') ' km/uGal'
      if(iop.eq.32) write(texto,'(3x,a)') ' km/uGal'
      call moveto(30,60,xy)
      call outgtext(texto)
      do 2 k=1,nc
      c=(k-1)*eq+ad
      write(texto,'(i7)') nint(c+id0-ad)
      call moveto(20,70+18*k,xy)
      nada=setcolorrgb(#000000)
      call outgtext(texto)
      nada=setcolorrgb(col(k))
    2 nada=rectangle($gfillinterior,90,70+18*k,110,70+18*(k+1))
      nada=setcolorrgb(#000000)
      nada=rectangle($gborder,89,70+(18-1),111,70+18*(1+nc))
      xc=ax+tam/2.
      yc=ay+tam/2.
      write(texto,'(a,i6,a)') 'Axis tic interval:',nint(tic),' m'
      call moveto(110,570,xy)
      call outgtext(texto)
      call settextposition(39,7,curpos)
      write(*,213) xc,yc,czm
  213 format(4x,'Center:  ',2f9.0,f7.0)
      write(*,214) 'Limits:',ax,ay,zt+tam/30
      write(*,214) '',ax+tam,ay+tam,zt-tam/3
  214 format(10x,a7,2x,2f9.0,f8.0)
      nli=20

      mx =480   ! 380
      l1x=180    ! 230
      l2x=mx/2+190    ! 440
      l1y=55
      l2y=mx+80      ! 460 
      call setviewport(l2x,l1y,l2x+mx,l1y+mx)
      nada=setcolorrgb(#e0e0e0)
      nada=rectangle($gfillinterior,0,0,mx,mx)
      nada=setcolorrgb(#000000)
      nada=rectangle($gborder,0,0,mx,mx)
      nada=setwindow(.true.,ax,ay,ax+tam,ay+tam)
      c=tic/5.
      do 3 i=1,40
      call moveto_w(dble(ax),dble(ay+i*tic),wxy)
      nada=lineto_w(ax+c,ay+i*tic)
      call moveto_w(dble(ax+tam),dble(ay+i*tic),wxy)
      nada=lineto_w(ax+tam-c,ay+i*tic)
      call moveto_w(dble(ax+i*tic),dble(ay),wxy)
      nada=lineto_w(ax+i*tic,ay+c)
      call moveto_w(dble(ax+i*tic),dble(ay+tam),wxy)
    3 nada=lineto_w(ax+i*tic,ay+tam-c)
      call cont(mb,np,xp,yp,jp)

      nada=setcolorrgb(#000000)
      call setviewport(l2x,l2y,l2x+mx,l2y+mx*(1/3.1+1/10.))
      nada=setwindow(.true.,ax,zt+tam/8,ax+tam,zt-tam/3.1)
      nada=setcolorrgb(#e0e0e0)
      nada=rectangle_w($gfillinterior,ax,zt,ax+tam,zt-tam/3.1)
      nada=setcolorrgb(#000000)
      call moveto_w(dble(ax),dble(zt+tam/8),wxy)
      nada=lineto_w(ax,zt-tam/3.1)
      nada=lineto_w(ax+tam,zt-tam/3.1)
      nada=lineto_w(ax+tam,zt+tam/8)
      i=zt/tic
      zz=(i+1)*tic
      do 5 i=1,20
      call moveto_w(dble(ax+tam),dble(zz-i*tic),wxy)
      nada=lineto_w(ax+tam-c,zz-i*tic)
      call moveto_w(dble(ax),dble(zz-i*tic),wxy)
    5 nada=lineto_w(ax+c,zz-i*tic)

      call setviewport(l1x,l1y,l1x+mx*(1/8.+1/3.1),l1y+mx)
      nada=setwindow(.true.,zt-tam/3.1,ay,zt+tam/8.,ay+tam)
      nada=setcolorrgb(#e0e0e0)
      nada=rectangle_w($gfillinterior,zt-tam/3.1,ay,zt,ay+tam)
      nada=setcolorrgb(#000000)
      call moveto_w(dble(zt+tam/8.),dble(ay+tam),wxy)
      nada=lineto_w(zt-tam/3.1,ay+tam)
      nada=lineto_w(zt-tam/3.1,ay)
      nada=lineto_w(zt+tam/8.,ay)
      do 6 i=1,20
      call moveto_w(dble(zz-i*tic),dble(ay),wxy)
      nada=lineto_w(zz-i*tic,ay+c)
      call moveto_w(dble(zz-i*tic),dble(ay+tam),wxy)
    6 nada=lineto_w(zz-i*tic,ay+tam-c)
 
   40 call setviewport(0,0,1000,900)
      nada=setcolorrgb(#f2fff2)  
      nada=rectangle($gfillinterior,20,mx+200,390,mx+300)
      call settextposition(44,5,curpos)
      write(*,'(a\)') ' >> exit 0, hor 1, WE 2, SN 3, obli 4 ? '
      read(*,'(i1)') ir
      if(ir.gt.4.or.ir.lt.0) go to 40
      if(ir.eq.0) go to 99
      if(ir.eq.1) write(*,'(4x,a\)') ' >> z value (m)  horiz section ? '
      if(ir.eq.2) write(*,'(4x,a\)') ' >> Y value (UTM)   WE profile ? '
      if(ir.eq.3) write(*,'(4x,a\)') ' >> X value (UTM)   NS profile ? '
      if(ir.eq.4) write(*,'(4x,a\)') ' >> azimut N (degs) ob profile ? '
      read(*,*) gg
      if(is.eq.2) then
        write(*,'(4x,a\)') ' >> Grid  file (a20) ? '
        read(*,'(a)') fname
        open(1,file=fname)
        if(ibl.eq.1) then
        write(*,'(4x,a\)') ' >> Base map file (a20) ? '
        read(*,'(a)') fname
        open(2,file=fname)
        endif
      endif
      if(ir.eq.1) then
c  horizontal section
      nada=setcolorrgb(#ffffff)
      nada=rectangle($gfillinterior,600,25,800,45)    
      nada=setcolorrgb(#000000)
      write(texto,'(a,i7,a)') ' Horizontal',nint(gg),' m'
      call moveto(600,30,xy)
      call outgtext(texto)
      call setviewport(l2x,l1y,l2x+mx,l1y+mx)
      nada=setcolorrgb(#d0d0d0)
      nada=rectangle($gfillinterior,0,0,mx,mx)
      nada=setwindow(.true.,ax,ay,ax+tam,ay+tam)
      zz=gg

      if(is.eq.0) then
       do 7 i=1,nb
       if(zz.lt.zi(i).or.zz.gt.zs(i)) go to 7
       dn=cb(i)
       tt=zs(i)
       bb=zi(i)
       call acol(nc,ad,eq,dn,k)
       nada=setcolorrgb(col(k))
       if(k.gt.nc) nada=setcolorrgb(col(2*nc-k))
       xx=xb(i)
       yy=yb(i)
       call SST(ms,n,x,y,z,xx,yy,zz,avm)
       if(eb(i).lt.0.and.rs.eq.0) go to 7
       xr=dx(i)/2.
       yr=dy(i)/2.
       if(eb(i).lt.0) xr=xr*rs
       if(eb(i).lt.0) yr=yr*rs
       nada=rectangle_w($gfillinterior,xx-xr,yy-yr,xx+xr,yy+yr)
    7  continue
      endif

      if(is.ge.1) then
       ns=0
       do 8 i=1,nb
       if(zi(i).gt.(zz+2*sr).or.zs(i).lt.(zz-2*sr)) go to 8
       ns=ns+1
       jb(ns)=i
    8  continue
       xx=ax
       yy=ay
       if(is.eq.2) write(1,227) nx+1,ny+1,xx,xx+ps*nx*fe,yy,yy+ps*ny*fe
  227  format('DSAA'/2i5/2f9.0/2f9.0/' -5000 5000')
       do 9 j=1,ny+1
       yy=ay+ps*(j-1)
       l=0
      do 11 i=1,nx+1
      xx=ax+ps*(i-1)
      call SST(ms,n,x,y,z,xx,yy,zz,avm)
             c=1.e-5/avm
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+1000,av2)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-1000,av3)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+500,av4)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-500,av5)
      av4=av4*2
      av5=av5*2
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy+1000,zz,av2)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy-1000,zz,av3)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy+1000,zz,av4)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy-1000,zz,av5)
      if(iop.ge.31)
     -c=4.e8*(abs(avm-av2)+abs(avm-av3)+abs(avm-av4)+abs(avm-av5))*fe*fe 

      if(iop.ge.30) go to 12
      if(c.lt.rl) c=9999.
      if(c.ge.rl) call am(ms,xb,yb,dx,dy,zs,zi,cb,jb,ns,xx,yy,zz,sr,c)
   12 call acol(nc,ad,eq,c,k)
      nada=setcolorrgb(col(k))
      if(k.gt.nc) nada=setcolorrgb(col(2*nc-k))
      if(is.ne.2) go to 11
      l=l+1
      if(l.gt.mt) call out(mt)
      ij(l)=c+id0-ad
      if(l.lt.10) go to 11
      write(1,'(10i7)') (ij(l),l=1,10)
      l=0
   11 if(k.ge.1) nada=rectangle_w($gfillinterior,xx-p,yy-p,xx+p,yy+p)
      if(l.gt.0) write(1,'(10i7)') (ij(i),i=1,l)
    9 if(is.eq.2) write(1,'(a)') ' '
      endif
      nada=setcolorrgb(#000000)
      nada=rectangle_w($gborder,ax,ay,ax+tam,ay+tam)
      c=tic/5
      do 13 i=1,40
      call moveto_w(dble(ax),dble(ay+i*tic),wxy)
      nada=lineto_w(ax+c,ay+i*tic)
      call moveto_w(dble(ax+tam),dble(ay+i*tic),wxy)
      nada=lineto_w(ax+tam-c,ay+i*tic)
      call moveto_w(dble(ax+i*tic),dble(ay),wxy)
      nada=lineto_w(ax+i*tic,ay+c)
      call moveto_w(dble(ax+i*tic),dble(ay+tam),wxy)
   13 nada=lineto_w(ax+i*tic,ay+tam-c)
      call cont(mb,np,xp,yp,jp)
      close(1)
      if(ibl.eq.1) close(2)
      go to 40
      endif
c  WE vertical profile
      if(ir.eq.2) then
      call setviewport(l2x,l1y,l2x+mx,l1y+mx)
      nada=setwindow(.true.,ax,ay,ax+tam,ay+tam)
      nada=setcolorrgb(#00f000)
      yy=gg
      nada=rectangle_w($gfillinterior,ax,yy-tam/500,ax+tam,yy+tam/500)
      c=mx*(1/3.1+1/10.)
      call setviewport(l2x,l2y,l2x+mx,l2y+c)
      nada=setwindow(.true.,ax,zt-tam/3.1,ax+tam,zt+tam/8.)
      nada=setcolorrgb(#d0d0d0)
      nada=rectangle_w($gfillinterior,ax,zt+tam/8.,ax+tam,zt-tam/3.1)

      if(is.eq.0) then
      yr=yy
      do 14 i=1,nb
      rr=dy(i)/2.
      if((yr+rr).lt.yb(i).or.(yr-rr).gt.yb(i)) go to 14
      xx=xb(i)
      xr=dx(i)/2.
      tt=zs(i)
      bb=zi(i)
      dn=cb(i)
      call acol(nc,ad,eq,dn,k)
      nada=setcolorrgb(col(k))
      if(k.gt.nc) nada=setcolorrgb(col(2*nc-k))
       call SST(ms,n,x,y,z,xx,yy,(bb+tt)/2.,avm)
       if(eb(i).lt.0.and.rs.eq.0) go to 14
       if(eb(i).lt.0) xr=xr*rs
       if(eb(i).lt.0) tt=tt-(tt-bb)/2.*(1-rs)
       if(eb(i).lt.0) bb=bb+(tt-bb)/2.*(1-rs)
      nada=rectangle_w($gfillinterior,xx-xr,bb,xx+xr,tt)
   14 continue
      ns=20
      endif
      
      if(is.ge.1) then
      ns=0
      do 15 i=1,nb
      zz=tb
      if(dy(i).gt.zz) zz=dy(i)
      yr=abs(yb(i)-yy)-zz/2.
      if(yr.le.sr) then
        ns=ns+1
        jb(ns)=i
      endif
   15 continue
      zz=(zt-tam/3.1)*fe+czm
      xx=ax
      if(is.eq.2) write(1,227) nx+1,nz+1,xx,xx+ps*nx*fe,zz,zz+ps*nz*fe
      do 18 j=1,nz+1
      zz=zt-tam/3.1+ps*(j-1)
      l=0
      do 17 i=1,nx+1
      xx=ax+ps*(i-1)
      call SST(ms,n,x,y,z,xx,yy,zz,avm)
                    c=1.e-5/avm
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+1000,av2)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-1000,av3)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+500,av4)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-500,av5)
      av4=av4*2
      av5=av5*2
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy+1000,zz,av2)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy-1000,zz,av3)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy+1000,zz,av4)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy-1000,zz,av5)
      if(iop.ge.31)
     -c=4.e8*(abs(avm-av2)+abs(avm-av3)+abs(avm-av4)+abs(avm-av5))*fe*fe 

      if(iop.ge.30) go to 16
      if(c.lt.rl) c=9999.
      if(c.gt.rl) call am(ms,xb,yb,dx,dy,zs,zi,cb,jb,ns,xx,yy,zz,sr,c)
   16 call acol(nc,ad,eq,c,k)
      nada=setcolorrgb(col(k))
      if(k.gt.nc) nada=setcolorrgb(col(nc*2-k))
      if(is.ne.2) go to 17
      l=l+1
      if(l.gt.mt) call out(mt)
      ij(l)=c+id0-ad
      if(l.lt.10) go to 17
      write(1,'(10i7)') (ij(l),l=1,10)
      l=0
   17 if(k.ge.1) nada=rectangle_w($gfillinterior,xx-p,zz-p,xx+p,zz+p)
      if(l.gt.0) write(1,'(10i7)') (ij(i),i=1,l)
   18 if(is.eq.2) write(1,'(a)') ' '
      ns=30
      endif
      
      xr=ax
      zr=zt
      c=tam/ns
      if(is.eq.2.and.ibl.eq.1) write(2,*) ns+1
      do 19 i=1,ns+1
      xx=ax+(i-1)*c
      call pertop(ms,x,y,z,n,xx,yy,zz,pp,rr,zr)
      nada=setcolorrgb(#ffffff)
      nada=rectangle_w($gfillinterior,xr,(zz+zr)/2.,xx,zt+tam/8)
      if(is.eq.2.and.ibl.eq.1) write(2,201) xx*fe+cxm,zz*fe+czm
  201 format(f9.0,f8.1,2x,f9.0)
      nada=setcolorrgb(#00f000)
      if(rr.ne.0.) call moveto_w(dble(xr),dble(zr),wxy)
      if(rr.ne.0.) nada=lineto_w(xx,zz)
      xr=xx
   19 zr=zz
      nada=setcolorrgb(#f0f000)
      call moveto_w(dble(ax),dble(0),wxy)
      nada=lineto_w(ax+tam,0)
      c=tic/5
      i=zt/tic
      zz=(i+10)*tic
      nada=setcolorrgb(#000000)
      do 20 i=1,30
      zr=zz-i*tic
      call moveto_w(dble(ax+tam),dble(zr),wxy)
      nada=lineto_w(ax+tam-c,zr)
      call moveto_w(dble(ax),dble(zr),wxy)
   20 nada=lineto_w(ax+c,zr)
      if(is.eq.2.and.ibl.eq.1) write(2,*) 2
      if(is.eq.2.and.ibl.eq.1) write(2,201) ax+tam,0.
      if(is.eq.2.and.ibl.eq.1) write(2,201) ax,0.
      call moveto_w(dble(ax),dble(zt+tam/8),wxy)
      nada=lineto_w(ax,zt-tam/3.1)
      nada=lineto_w(ax+tam,zt-tam/3.1)
      nada=lineto_w(ax+tam,zt+tam/8)
      if(ibl.eq.1) close(2)
      close(1)
      write(texto,'(a,i7,a)') '   WE  Y=',nint(gg),' UTM'
      nada=setcolorrgb(#000000)
      call setviewport(0,0,1000,900)
      call moveto(590,570,xy)
      call outgtext(texto)
      go to 40
      endif
      if(ir.eq.3) then
c  SN vertical profile
      nada=setcolorrgb(#ffffff)
      nada=rectangle($gfillinterior,200,26,350,46)    
      nada=setcolorrgb(#000000)
      write(texto,'(a,i6,a)') '  NS  X=',nint(gg),' UTM'
      call moveto(200,30,xy)
      call outgtext(texto)
      call setviewport(l2x,l1y,l2x+mx,l1y+mx)
      nada=setwindow(.true.,ax,ay,ax+tam,ay+tam)
      nada=setcolorrgb(#00f000)
      xx=gg
      nada=rectangle_w($gfillinterior,xx-tam/500,ay,xx+tam/500,ay+tam)
      c=mx*(1/8.+1./3.1)
      call setviewport(l1x,l1y,l1x+c,l1y+mx)
      nada=setwindow(.true.,zt-tam/3.1,ay,zt+tam/8.,ay+tam)
      nada=setcolorrgb(#d0d0d0)
      nada=rectangle_w($gfillinterior,zt-tam/3.1,ay,zt+tam/8,ay+tam)

      if(is.eq.0) then
      do 21 i=1,nb
      rr=dx(i)/2.
      if((xx+rr).lt.xb(i)) go to 21
      if((xx-rr).gt.xb(i)) go to 21
      yy=yb(i)
      yr=dy(i)/2.
      tt=zs(i)
      bb=zi(i)
      dn=cb(i)
      call acol(nc,ad,eq,dn,k)
      nada=setcolorrgb(col(k))
      if(k.gt.nc) nada=setcolorrgb(col(2*nc-k))
      call SST(ms,n,x,y,z,xx,yy,(bb+tt)/2.,avm)
       if(eb(i).lt.0.and.rs.eq.0) go to 21
       if(eb(i).lt.0) yr=yr*rs
       if(eb(i).lt.0) tt=tt-(tt-bb)/2.*(1-rs)
       if(eb(i).lt.0) bb=bb+(tt-bb)/2.*(1-rs)
      nada=rectangle_w($gfillinterior,bb,yy-yr,tt,yy+yr)
   21 continue
      ns=20
      endif
      
      if(is.ge.1) then
      ns=0
      do 22 i=1,nb
      zz=tb
      if(dx(i).gt.zz) zz=dx(i)
      xr=abs(xb(i)-xx)-zz/2.
      if(xr.le.sr) then
      ns=ns+1
      jb(ns)=i
      endif
   22 continue
      zz=(zt-tam/3.1)*fe+czm
      yy=ay*fe+cym
      if(is.eq.2) write(1,227) ny+1,nz+1,yy,yy+ny*ps,zz,zz+ps*nz
      do 25 j=1,nz+1
      zz=zt-tam/3.1+ps*(j-1)
      l=0
      do 24 i=1,ny+1
      yy=ay+ps*(i-1)
      call SST(ms,n,x,y,z,xx,yy,zz,avm)
                    c=1.e-5/avm
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+1000,av2)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-1000,av3)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+500,av4)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-500,av5)
      av4=av4*2
      av5=av5*2
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy+1000,zz,av2)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy-1000,zz,av3)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy+1000,zz,av4)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy-1000,zz,av5)
      if(iop.ge.31)
     -c=4.e8*(abs(avm-av2)+abs(avm-av3)+abs(avm-av4)+abs(avm-av5))*fe*fe 

      if(iop.ge.30) go to 23
      if(c.lt.rl) c=9999.
      if(c.ge.rl) call am(ms,xb,yb,dx,dy,zs,zi,cb,jb,ns,xx,yy,zz,sr,c)
   23 call acol(nc,ad,eq,c,k)
      nada=setcolorrgb(col(k))
      if(k.gt.nc) nada=setcolorrgb(col(2*nc-k))
      if(is.ne.2) go to 24
      l=l+1
      if(l.gt.mt) call out(mt)
      ij(l)=c+id0-ad
      if(l.lt.10) go to 24
      write(1,'(10i7)') (ij(l),l=1,10)
      l=0
   24 if(k.ge.1) nada=rectangle_w($gfillinterior,zz-p,yy-p,zz+p,yy+p)
      if(l.gt.0) write(1,'(10i7)') (ij(i),i=1,l)
   25 if(is.eq.2) write(1,'(a)') ' '
      ns=30
      endif
      yr=ay
      zr=zt
      c=tam/ns
      if(is.eq.2.and.ibl.eq.1) write(2,*) ns+1
      do 26 i=1,ns+1
      yy=ay+(i-1)*c
      call pertop(ms,x,y,z,n,xx,yy,zz,pp,rr,zr)
      nada=setcolorrgb(#ffffff)
      nada=rectangle_w($gfillinterior,(zr+zz)/2.,yr,zt+tam/8,yy)
      nada=setcolorrgb(#00f000)
      if(is.eq.2.and.ibl.eq.1) write(2,201) yy,zz
      if(rr.ne.0.) call moveto_w(dble(zr),dble(yr),wxy)
      if(rr.ne.0.) nada=lineto_w(zz,yy)
      yr=yy
   26 zr=zz
      nada=setcolorrgb(#f0f000)
      call moveto_w(dble(0),dble(ay),wxy)
      nada=lineto_w(0,ay+tam)
      nada=setcolorrgb(#000000)
      call moveto_w(dble(zt+tam/8.),dble(ay+tam),wxy)
      nada=lineto_w(zt-tam/3.1,ay+tam)
      nada=lineto_w(zt-tam/3.1,ay)
      nada=lineto_w(zt+tam/8.,ay)
      i=zt/tic
      zz=(i+10)*tic
      c=tic/5
      do 27 i=1,30
      zr=zz-i*tic
      call moveto_w(dble(zr),dble(ay),wxy)
      nada=lineto_w(zr,ay+c)
      call moveto_w(dble(zr),dble(ay+tam),wxy)
   27 nada=lineto_w(zr,ay+tam-c)
      if(is.eq.2.and.ibl.eq.1) write(2,*) 2
      if(is.eq.2.and.ibl.eq.1) write(2,201) ay+tam,0.
      if(is.eq.2.and.ibl.eq.1) write(2,201) ay,0.
      close(1)
      if(ibl.eq.1) close(2)
      go to 40
      endif
      if(ir.eq.4) then
c  oblique vertical profiles
      xc=ax+tam/2.
      yc=ay+tam/2.
      write(*,'(4x,a\)') ' >> X,Y center ? '
      read(*,*) xc1,yc1
      if(xc1.ne.0.and.yc1.ne.0) then
      xc=xc1
      yc=yc1
      endif
      rr=gg
      if(gg.lt.0) rr=gg+360
      if(gg.gt.180) rr=gg-180
      a1=sin(gg*0.01745329)
      a2=cos(gg*0.01745329)
c  WE oblique vertical profiles
      if(45.le.rr.and.rr.le.135.) then
      gg=a2/a1
      call setviewport(l2x,l1y,l2x+mx,l1y+mx)
      nada=setwindow(.true.,ax,ay,ax+tam,ay+tam)
      xx=ax
      yy=yc+(xx-xc)*gg
      nada=setcolorrgb(#00f000)
      call moveto_w(dble(xx),dble(yy),wxy)
      xx=ax+tam
      yy=yc+(xx-xc)*gg
      nada=lineto_w(xx,yy)
      c=mx*(1/3.1+1/10.)
      call setviewport(l2x,l2y,l2x+mx,l2y+c)
      nada=setwindow(.true.,ax,zt-tam/3.1,ax+tam,zt+tam/8.)
      nada=setcolorrgb(#d0d0d0)
      nada=rectangle_w($gfillinterior,ax,zt-tam/3.1,ax+tam,zt+tam/8.)

      if(is.eq.0) then
      do 28 i=1,nb
      xx=xb(i)-xc
      yy=yb(i)-yc
      if(abs(xx).eq.0.) xx=0.1
      tt=yy/xx
      c=(dx(i)+dy(i))/3.7
      if(c.eq.0.) go to 28
      rr=sqrt(xx*xx+yy*yy)
      if(rr.ge.c) then
      c=rr/c
      a2=(c*tt+1.)/(c-tt)
      a1=(c*tt-1.)/(c+tt)
      if(a1.gt.gg.or.gg.ge.a2) go to 28
      endif
      xx=xb(i)
      yy=yb(i)
      xr=dx(i)/2.
      tt=zs(i)
      bb=zi(i)
      dn=cb(i)
      call acol(nc,ad,eq,dn,k)
      nada=setcolorrgb(col(k))
      if(k.gt.nc) nada=setcolorrgb(col(2*nc-k))
       call SST(ms,n,x,y,z,xx,yy,(bb+tt)/2.,avm)
       c=1.e9*avm*fe*fe
       if(eb(i).lt.0.and.rs.eq.0) go to 28
       if(eb(i).lt.0) xr=xr*rs
       if(eb(i).lt.0) tt=tt-(tt-bb)/2.*(1-rs)
       if(eb(i).lt.0) bb=bb+(tt-bb)/2.*(1-rs)
      nada=rectangle_w($gfillinterior,xx-xr,bb,xx+xr,tt)
   28 continue
      ns=20
      endif
      if(is.ge.1) then
      ns=0
      dn=2.*sr
      do 29 i=1,nb
      xr=xb(i)-xc
      yr=yb(i)-yc
      dd=abs(yr*a1-xr*a2)
      if(dd.le.dn) then
      ns=ns+1
      jb(ns)=i
      endif
   29 continue
      xx=ax
      yy=tam
      zz=zt
      if(is.eq.2) write(1,227) nn+1,nz+1,xx,xx+yy,zz-yy/2,zz+yy/10
      do 33 j=1,nz+1
      zz=zt-tam/3.1+ps*(j-1)
      l=0
      do 32 i=1,nn+1
      xx=ax+ps*(i-1)
      yy=yc+(xx-xc)*gg
      call SST(ms,n,x,y,z,xx,yy,zz,avm)
                    c=1.e-5/avm
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+1000,av2)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-1000,av3)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+500,av4)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-500,av5)
      av4=av4*2
      av5=av5*2
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy+1000,zz,av2)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy-1000,zz,av3)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy+1000,zz,av4)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy-1000,zz,av5)
      if(iop.ge.31)
     -c=4.e8*(abs(avm-av2)+abs(avm-av3)+abs(avm-av4)+abs(avm-av5))*fe*fe 

      if(iop.ge.30) go to 31
      if(c.lt.rl) c=9999.
      if(c.ge.rl) call am(ms,xb,yb,dx,dy,zs,zi,cb,jb,ns,xx,yy,zz,sr,c)
   31 call acol(nc,ad,eq,c,k)
      nada=setcolorrgb(col(k))
      if(k.gt.nc) nada=setcolorrgb(col(2*nc-k))
      if(is.ne.2) go to 32
      l=l+1
      ij(l)=c+id0-ad
      if(l.lt.10) go to 32
      write(1,'(10i7)') (ij(l),l=1,10)
      l=0
   32 if(k.ge.1) nada=rectangle_w($gfillinterior,xx-p,zz-p,xx+p,zz+p)
      if(l.gt.0) write(1,'(10i7)') (ij(i),i=1,l)
   33 if(is.eq.2) write(1,'(a)') ' '
      ns=30
      endif
      xr=ax
      yr=yc+(xr-xc)*gg
      zr=zt
      c=tam/ns
      if(is.eq.2.and.ibl.eq.1) write(2,*) ns+1
      do 34 i=1,ns+1
      xx=ax+(i-1)*c
      yy=yc+(xx-xc)*gg
      call pertop(ms,x,y,z,n,xx,yy,zz,pp,rr,zr)
      nada=setcolorrgb(#ffffff)
      nada=rectangle_w($gfillinterior,xr,(zr+zz)/2.,xx,zt+tam/8)
      nada=setcolorrgb(#00f000)
      if(is.eq.2.and.ibl.eq.1) write(2,201) xx,zz,yy
      if(rr.ne.0.) call moveto_w(dble(xr),dble(zr),wxy)
      if(rr.ne.0.) nada=lineto_w(xx,zz)
      xr=xx
   34 zr=zz
      nada=setcolorrgb(#f0f000)
      call moveto_w(dble(ax),0.d0,wxy)
      nada=lineto_w(ax+tam,0)
      c=tic/5
      i=zt/tic
      zz=(i+10)*tic
      nada=setcolorrgb(#000000)
      do 35 i=1,30
      zr=zz-i*tic
      call moveto_w(dble(ax+tam),dble(zr),wxy)
      nada=lineto_w(ax+tam-c,zr)
      call moveto_w(dble(ax),dble(zr),wxy)
   35 nada=lineto_w(ax+c,zr)
      if(is.eq.2.and.ibl.eq.1) write(2,*) 2
      if(is.eq.2.and.ibl.eq.1) write(2,201) (ax+tam)*fe+cxm,0.
      if(is.eq.2.and.ibl.eq.1) write(2,201) ax*fe+cxm,0.
      nada=setcolorrgb(#000000)
      call moveto_w(dble(ax),dble(zt+tam/8),wxy)
      nada=lineto_w(ax,zt-tam/3.1)
      nada=lineto_w(ax+tam,zt-tam/3.1)
      nada=lineto_w(ax+tam,zt+tam/8)
      if(ibl.eq.1) close(2)
      close(1)
      call settextposition(37,75,curpos)
      write(*,'(a,f6.0,a)') 'Oblique of Azim',rr,'     '
      go to 40
      endif
c  SN oblique vertical profiles
      if(45.gt.rr.or.rr.gt.135.) then
      call settextposition(3,25,curpos)
      write(*,218) rr
  218 format(' Oblique Az',f6.0,'    ')
      gg=a1/a2
      call setviewport(l2x,l1y,l2x+mx,l1y+mx)
      nada=setwindow(.true.,ax,ay,ax+tam,ay+tam)
      yy=ay
      xx=xc+(yy-yc)*gg
      nada=setcolorrgb(#00f000)
      call moveto_w(dble(xx),dble(yy),wxy)
      yy=ay+tam
      xx=xc+(yy-yc)*gg
      nada=lineto_w(xx,yy)
      c=mx*(1/8.+1./3.1)
      call setviewport(l1x,l1y,l1x+c,l1y+mx)
      nada=setwindow(.true.,zt-tam/3.1,ay,zt+tam/8,ay+tam)
      nada=setcolorrgb(#d0d0d0)
      nada=rectangle_w($gfillinterior,zt-tam/3.1,ay,zt+tam/8,ay+tam)

      if(is.eq.0) then
      do 36 i=1,nb
      xx=xb(i)-xc
      yy=yb(i)-yc
      if(abs(yy).eq.0.) yy=0.1
      tt=xx/yy
      c=(dx(i)+dy(i))/3.7
      if(c.eq.0.) go to 36
      rr=sqrt(xx*xx+yy*yy)
      if(rr.ge.c) then
      c=c/rr
      a1=(tt-c)/(1.+c*tt)
      a2=(tt+c)/(1.-c*tt)
      if(a1.gt.gg.or.gg.ge.a2) go to 36
      endif
      xx=xb(i)
      yy=yb(i)
      yr=dy(i)/2.
      tt=zs(i)
      bb=zi(i)
      dn=cb(i)
      call acol(nc,ad,eq,dn,k)
      nada=setcolorrgb(col(k))
      if(k.gt.nc) nada=setcolorrgb(col(2*nc-k))
       if(eb(i).lt.0.and.rs.eq.0) go to 36
       if(eb(i).lt.0) yr=yr*rs
       if(eb(i).lt.0) tt=tt-(tt-bb)/2.*(1-rs)
       if(eb(i).lt.0) bb=bb+(tt-bb)/2.*(1-rs)
      nada=rectangle_w($gfillinterior,bb,yy-yr,tt,yy+yr)
   36 continue
      ns=20
      endif
      if(is.ge.1) then
      ns=0
      dn=2.*sr
      do 37 i=1,nb
      xr=xb(i)-xc
      yr=yb(i)-yc
      dd=abs(yr*a1-xr*a2)
      if(dd.le.dn) then
      ns=ns+1
      jb(ns)=i
      endif
   37 continue
      xx=ay*fe+cym
      yy=tam*fe
      zz=zt*fe+czm
      if(is.eq.2) write(1,227) nn+1,nz+1,xx,xx+yy,zz-yy/2,zz+yy/10
      do 41 j=1,nz+1
      zz=zt-tam/3.1+ps*(j-1)
      l=0
      do 39 i=1,nn+1
      yy=ay+ps*(i-1)
      xx=xc+(yy-yc)*gg
      call SST(ms,n,x,y,z,xx,yy,zz,avm)
                    c=1.e-5/avm
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+1000,av2)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-1000,av3)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz+500,av4)
      if(iop.eq.31) call SST(ms,n,x,y,z,xx,yy ,zz-500,av5)
      av4=av4*2
      av5=av5*2
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy+1000,zz,av2)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy-1000,zz,av3)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx-1000,yy+1000,zz,av4)
      if(iop.eq.32) call SST(ms,n,x,y,z,xx+1000,yy-1000,zz,av5)
      if(iop.ge.31)
     -c=4.e8*(abs(avm-av2)+abs(avm-av3)+abs(avm-av4)+abs(avm-av5))*fe*fe 

      if(iop.ge.30) go to 38
      if(c.lt.rl) c=9999.
      if(c.ge.rl) call am(ms,xb,yb,dx,dy,zs,zi,cb,jb,ns,xx,yy,zz,sr,c)
   38 call acol(nc,ad,eq,c,k)
      nada=setcolorrgb(col(k))
      if(k.gt.nc) nada=setcolorrgb(col(2*nc-k))
      if(is.ne.2) go to 39
      l=l+1
      ij(l)=c+id0-ad
      if(l.lt.10) go to 39
      write(1,'(10i7)') (ij(l),l=1,10)
      l=0
   39 if(k.ge.1) nada=rectangle_w($gfillinterior,zz-p,yy-p,zz+p,yy+p)
      if(l.gt.0) write(1,'(10i7)') (ij(i),i=1,l)
   41 if(is.eq.2) write(1,'(a)') ' '
      ns=30
      endif
      yr=ay
      xr=xc+(yr-yc)*gg
      zr=zt
      c=tam/ns
      if(is.eq.2.and.ibl.eq.1) write(2,*) ns+1
      do 42 i=1,ns+1
      yy=ay+(i-1)*c
      xx=xc+(yy-yc)*gg
      call pertop(ms,x,y,z,n,xx,yy,zz,pp,rr,zr)
      nada=setcolorrgb(#ffffff)
      nada=rectangle_w($gfillinterior,(zr+zz)/2.,yr,zt+tam/8,yy)
      nada=setcolorrgb(#00f000)
      if(is.eq.2.and.ibl.eq.1) write(2,201)
     - yy*fe+cym,zz*fe+czm,xx*fe+cxm
      if(rr.ne.0.) call moveto_w(dble(zr),dble(yr),wxy)
      if(rr.ne.0.) nada=lineto_w(zz,yy)
      yr=yy
   42 zr=zz
      nada=setcolorrgb(#f0f000)
      call moveto_w(0.d0,dble(ay),wxy)
      nada=lineto_w(-czm/fe,ay+tam)
      nada=setcolorrgb(#000000)
      call moveto_w(dble(zt+tam/8.),dble(ay+tam),wxy)
      nada=lineto_w(zt-tam/3.1,ay+tam)
      nada=lineto_w(zt-tam/3.1,ay)
      nada=lineto_w(zt+tam/8.,ay)
      i=zt/tic
      zz=(i+10)*tic
      c=tic/5
      do 43 i=1,30
      zr=zz-i*tic-czm/fe
      call moveto_w(dble(zr),dble(ay),wxy)
      nada=lineto_w(zr,ay+c)
      call moveto_w(dble(zr),dble(ay+tam),wxy)
   43 nada=lineto_w(zr,ay+tam-c)
      if(is.eq.2.and.ibl.eq.1) write(2,*) 2
      if(is.eq.2.and.ibl.eq.1) write(2,201) (ay+tam)*fe+cym,0.
      if(is.eq.2.and.ibl.eq.1) write(2,201) ay*fe+cym,0.
      close(1)
      if(ibl.eq.1) close(2)
      go to 40
      endif
      go to 40
      endif
   99 return
      end
c**********************************************************************
c     Mean sensitivity of the network for a mass located en XX,YY,ZZ
c**********************************************************************
      subroutine SST(ms,ns,x,y,z,xx,yy,zz,s)
      integer x(ms),y(ms),z(ms)
      r=0
      do 2 i=1,ns
    2 r=r+x(i)-x(1)+y(i)-y(1)
      r=r/ns/2.5    !/2.
      rr=r*r
      sp=0
      am=0.
      do 1 i=1,ns
      xr=xx-x(i)
      yr=yy-y(i)
      zr=zz-z(i)
      d2=xr*xr+yr*yr+zr*zr+1.
      p=1
      if(d2.lt.rr) p=d2/rr
      p=p*p
      aa=zr*zr/d2/d2/d2
      am=am+aa*p 
      sp=sp+p
    1 continue
      am=sqrt(am/sp)
      s=am*1000.*6.672e-3    ! uGal media por tonelada
      return
      end

c*********************************************************************
c      Perspectiva 3D del modelo
C*********************************************************************
      subroutine perspec(ax,ay,tam,fe,
     -mc,nc,xc,yc,sx,sy,zb,zt,ij,pr,sen,
     -ms,ns,xs,ys,zs,
     -mb,nb,xb,yb,jb,tic)

      use msflib
      use msfwin
      use portlib
      TYPE (rccoord) curpos
      TYPE (wxycoord) wxy
      TYPE (xycoord) xy
      character*1 key /'A'/
      integer*4 col(70),xb(mb),yb(mb),jb(mb),ij(mc),sx(mb),sy(mb)
     -,xs(ms),ys(ms),zs(ms),xc(mc),yc(mc),zb(mc),zt(mc),pr(mc),sen(mc)
      integer*2 nada
      character*120 texto,tex2

      call setviewport(0,0,1024,768)

      pi=3.141592653589793D0
      dera=pi/180.

      ncol=32
      col( 1)=#400000     ! b g r
      col( 2)=#600000
      col( 3)=#800000
      col( 4)=#900000
      col( 5)=#a00000
      col( 6)=#b00000
      col( 7)=#c00000
      col( 8)=#d00000
      col( 9)=#e00000
      col(10)=#f00000
      col(11)=#f02000
      col(12)=#f04000
      col(13)=#f06000
      col(14)=#f08000
      col(15)=#f09000
      col(16)=#f0a000
      col(17)=#f0b000
      col(18)=#f0c000
      col(19)=#f0c020
      col(20)=#f0d020
      col(21)=#f0e020
      col(22)=#f0f000
      col(23)=#f0f020
      col(24)=#f0f040
      col(25)=#f0f060
      col(26)=#f0f080
      col(27)=#f0f0a0
      col(28)=#f0f0b0
      col(29)=#f0f0c0
      col(30)=#f0f0d0
      col(31)=#f0f0e0
      col(32)=#f0f0f0

      col(33)=#000040     ! b g r
      col(34)=#000060
      col(35)=#000080
      col(36)=#000090
      col(37)=#0000a0
      col(38)=#0000b0
      col(39)=#0000c0
      col(40)=#0000d0
      col(41)=#0000e0
      col(42)=#0000f0
      col(43)=#0020f0
      col(44)=#0040f0
      col(45)=#0060f0
      col(46)=#0080f0
      col(47)=#0090f0
      col(48)=#00a0f0
      col(49)=#00b0f0
      col(50)=#00c0f0
      col(51)=#20c0f0
      col(52)=#20d0f0
      col(53)=#20e0f0
      col(54)=#00f0f0
      col(55)=#20f0f0
      col(56)=#40f0f0
      col(57)=#60f0f0
      col(58)=#80f0f0
      col(59)=#a0f0f0
      col(60)=#b0f0f0
      col(61)=#c0f0f0
      col(62)=#d0f0f0
      col(63)=#e0f0f0
      col(64)=#f0f0f0

c -------------  Colores y limites -------------

      y1=9.d9                          ! Limites modelo
      y2=-y1
      do 2 i=1,nc
      ij(i)=0
      if(pr(i).eq.0.) go to 2
      c=yc(i)-sy(i)/2.
      if(c.le.y1) y1=c
      c=yc(i)+sy(i)/2.
      if(c.ge.y2) y2=c
    2 continue
      c=(y2-y1)/ncol/3.
      if(c.lt.2) c=2
      y1=y1-c
      y2=y2+c
           y1=ay
           y2=ay+tam
      py=(y2-y1)/ncol                    ! pasos de color
      z0=0                           ! Techo del modelo en estaciones
      do 3 i=1,ns
    3 z0=z0+zs(i)
      z0=z0/ns

c ------------- Perpectiva por temas -------

   29 call clearscreen($clearscreen)
      write(texto,'(a,f7.0,a))') 'Perspective view      Tick=',tic*fe,
     -'m       Rotation: S,D,W,E     Exit:9' 
      nada=setfont('t''Courier''h14w7')
      call moveto(100,10,xy)
      call outgtext(texto)

      call dial5(al,be,del,nl,ip,in,svh,bot)
      del=del/fe
      ly=250
      lx=450
      c=500/tam
      fx=c
      fy=c*svh
      xg=ax+tam/2
      yg=ay+tam/2
      zg=z0
                                                 ! Ejes y cubo limite
  10  continue
      sa=sin(al*dera)
      ca=cos(al*dera)
      sb=sin(be*dera)
      cb=cos(be*dera)
      nada=setcolorrgb(#ffffff)
      nada=rectangle($gfillinterior,10,30,1024,768)
      nada=setcolorrgb(#000000)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay+tam,z0,
     -ixp,iyp)
      call moveto(ixp,iyp,xy)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay,z0,
     -ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay,
     -z0-bot*tam,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay+tam,
     -z0-bot*tam,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay+tam,z0,
     -ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam,ay+tam,
     -z0,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam,ay+tam,
     -z0-bot*tam,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay+tam,
     -z0-bot*tam,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay,
     -z0-bot*tam,ixp,iyp)
      call moveto(ixp,iyp,xy)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam,ay,
     -z0-bot*tam,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam,ay+tam,
     -z0-bot*tam,ixp,iyp)
      nada=lineto(ixp,iyp)

      d=tam/(nl-1)                       ! paso grid 3D
      d2x=d/2*fx
      d2y=d/2*fy
      do 4 j=1,nl                       ! para cada plano de Y
      yy=ay+tam-(j-1)*d
      nce=0
      do 5 l=1,nc
      if(pr(l).eq.0) go to 5
      if(sen(l).gt.40) go to 5
      if(zt(l).gt.(z0-del)) go to 5
      if(zt(l).lt.(z0-bot*tam)) go to 5
      if(in.eq.0.and.pr(l).le.0) go to 5
      if(ip.eq.0.and.pr(l).ge.0) go to 5
      if(abs(yy-yc(l)).gt.sy(l)/2.) go to 5
      nce=nce+1
      ij(nce)=l
    5 continue
      if(nce.eq.0) go to 4
      do 9 k=1,nl                       ! para cada nivel de Z
      zz=z0-(k-1)*d
      if(zz.lt.(z0-bot*tam)) go to 9
      do 8 i=1,nl                       ! para cada nivel de X
      xx=ax+(i-1)*d
      ll=0
      do 6 le=1,nce
      l=ij(le)
      if(zz.gt.zt(l).or.zz.lt.zb(l)) go to 6
      dx=xc(l)-xx
      if(abs(dx).gt.sx(l)/2.) go to 6
      ll=l
    6 continue
      if(ll.eq.0) go to 8
      k1=(yy-y1)/py                       ! asigna color
      if(pr(ll).gt.0) k1=ncol+k1
      cx=xx+(yy-ay)*fpx
      cz=zz+(yy-ay)*fpy
      nada=setcolorrgb(col(k1))
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,xx,yy,zz,
     -ixp,iyp)
      nada=rectangle($gfillinterior,ixp-d2x,iyp-d2y,ixp+d2x,iyp+d2y)
    8 continue
    9 continue
    4 continue

      do 1 i=1,nl        ! ejes con color
      cy=d*(i-1)
      k1=(ay+cy-y1)/py
      if(k1.le.0) k1=1
      if(k1.gt.ncol) k1=ncol
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,
     -ax,ay+cy,z0,ixp,iyp)
      nada=setcolorrgb(col(k1+ncol))
      if(ip.eq.1) nada=rectangle($gfillinterior,ixp-5,iyp,ixp,iyp+5)
      cx=ax+tam+fpx*cy
      cz=z0+fpy*cy
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,
     -ax+tam,ay+cy,z0,ixp,iyp)
      if(k1.lt.0) k1=1
      if(k1.gt.ncol) k1=ncol
      nada=setcolorrgb(col(k1))
      if(in.eq.1) nada=rectangle($gfillinterior,ixp,iyp,ixp+5,iyp+5)
    1 continue

      nada=setcolorrgb(#000000)
      r=tam/50
      do 13 i=1,20
      c=tic*i
      if(c.gt.bot*tam) go to 13
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay,
     -z0-c,ixp,iyp)
      call moveto(ixp,iyp,xy)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax-r,ay,
     -z0-c,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam,ay,
     -z0-c,ixp,iyp)
      call moveto(ixp,iyp,xy)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam+r,ay,
     -z0-c,ixp,iyp)
      nada=lineto(ixp,iyp)
   13 continue
      do 33 i=1,20
      c=tic*i
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+c,ay,
     -z0,ixp,iyp)
      call moveto(ixp,iyp,xy)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+c,ay,
     -z0-r,ixp,iyp)
      if(c.le.tam) nada=lineto(ixp,iyp)
   33 continue

      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam,ay,
     -z0-tam*bot,ixp,iyp)
      call moveto(ixp,iyp,xy)  
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam,ay,
     -z0,ixp,iyp)
      nada=lineto(ixp,iyp)  ! vertical delantero derecho
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam,ay+tam,
     -z0,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax+tam,ay,
     -z0,ixp,iyp)
      call moveto(ixp,iyp,xy)  
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay,
     -z0,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay,
     -z0-bot*tam,ixp,iyp)
      nada=lineto(ixp,iyp)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay,
     -z0,ixp,iyp)
      call moveto(ixp,iyp,xy)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,ax,ay+tam,
     -z0,ixp,iyp)
      nada=lineto(ixp,iyp)


      nada=setcolorrgb(#00a000)        ! contorno map.bln
      k=0
      do 11 i=1,nb
      if(xb(i).lt.ax.or.xb(i).gt.(ax+tam)) k=0
      if(xb(i).lt.ax.or.xb(i).gt.(ax+tam)) go to 11
      if(yb(i).lt.ay.or.yb(i).gt.(ay+tam)) go to 11
      xx=xb(i)
      yy=yb(i) 
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,
     -xx,yy,z0,ixp,iyp)
      if(jb(i).eq.1.or.k.eq.0) call moveto(ixp,iyp,xy)
      if(jb(i).eq.0) nada=lineto(ixp,iyp)
      k=1 
   11 continue

      nada=setcolorrgb(#005000)
      r=tam/400*fx                         ! estaciones dato
      if(ns.gt.400) r=r/2.
      do 12 i=1,ns
      xx=xs(i)
      yy=ys(i)
      zz=zs(i)
      call pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,
     -xx,yy,zz,ixp,iyp)
   12 nada=ellipse($gborder,ixp-r,iyp-r/3,ixp+r,iyp+r/3)
       
      key=GETCHARQQ()
      if(key.eq.'9') go to 99
      if(key.eq.'s') be=be-5
      if(key.eq.'d') be=be+5
      if(key.eq.'w') al=al+5
      if(key.eq.'e') al=al-5

      if(al.lt.-90.) al=-90 
      if(al.gt.90.) al=90 
      go to 10      
   
   99 return
      end
c***********************************************************************
c     Dialog interface 5                                               *
c***********************************************************************
      subroutine dial5(al,be,del,nl,ip,in,svh,bot)
      use dialogm
      include 'View.fd'
      logical retlog,nada
      character*20 texto
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Perspec, dlg)

      retlog = DlgSet(dlg,IDC_EDIT1 ,'15')
      retlog = DlgSet(dlg,IDC_EDIT2 ,'15')
      retlog = DlgSet(dlg,IDC_EDIT3 ,' 0')
      retlog = DlgSet(dlg,IDC_EDIT4 ,'100')
      retlog = DlgSet(dlg,IDC_EDIT5 ,'1.1')
      retlog = DlgSet(dlg,IDC_EDIT6 ,'0.6')
      retlog = DlgSet(dlg,IDC_CHECK1,.true.)
      retlog = DlgSet(dlg,IDC_CHECK2,.true.)

      retint = DlgModal( dlg )

      nada = DlgGet( dlg, IDC_EDIT1, texto)
      read(texto,*) al
      nada = DlgGet( dlg, IDC_EDIT2, texto)
      read(texto,*) be
      nada = DlgGet( dlg, IDC_EDIT3, texto)
      read(texto,*) del
      nada = DlgGet( dlg, IDC_EDIT4, texto)
      read(texto,*) nl
      nada = DlgGet( dlg, IDC_EDIT5, texto)
      read(texto,*) svh
      nada = DlgGet( dlg, IDC_EDIT6, texto)
      read(texto,*) bot
      ip=0
      retlog = DlgGet( dlg, IDC_CHECK1, nada)
      if(nada) ip=1
      in=0
      retlog = DlgGet( dlg, IDC_CHECK2, nada)
      if(nada) in=1

      call DlgUninit(dlg)
      return
      end

c***********************************************************************
c      Perspectiva 3D de un punto
C*********************************************************************
      subroutine pers(sa,ca,sb,cb,lx,ly,fx,fy,xg,yg,zg,xx,yy,zz,ixp,iyp)
      ixp=lx+( (xx-xg)*ca +(yy-yg)*sa)*fx
      iyp=ly-(-(xx-xg)*sa*sb+(yy-yg)*ca*sb+(zz-zg)*cb)*fy

      return
      end
c********************************************************************
      subroutine clean(mc,nb,xb,yb,dx,dy,zi,zs,db,cb)
      integer zs(mc),zi(mc),xb(mc),yb(mc),dx(mc),dy(mc),db(mc),cb(mc)
      n=0
      do 1 i=1,nb
      cb(i)=0
      if(db(i).eq.0) go to 1 
      n=n+1
      nn=0
      do 2 j=1,nb
      if(db(j).eq.0.or.(db(j)*db(i)).lt.0) go to 2
      if(abs(xb(i)-xb(j)).gt.(dx(i)+dx(j))*0.7) go to 2
      if(abs(yb(i)-yb(j)).gt.(dy(i)+dy(j))*0.7) go to 2
      if(abs(zi(j)-zi(i)).gt.(zs(j)-zi(j)+zs(i)-zi(i))*0.7) go to 2  
      nn=nn+1
    2 continue
      if(nn.gt.4) cb(i)=db(i)
    1 continue      
      nn=0
      do 3 i=1,nb
      if(cb(i).ne.0) nn=nn+1
    3 db(i)=cb(i)   
      write(*,'(/30x,a/40x,3(a,i6)/30x,a)') 
     -'Number of filled cells:','  Initial=',n,'   Removed=',n-nn,
     -'   Resulting=',nn,'See out.txt'
      read(*,'(i1)') n
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