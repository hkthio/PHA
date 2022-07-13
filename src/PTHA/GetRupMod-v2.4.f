C    
C     -- Par file = Event-scaling-segment-top-bottom-asperity-epsilon.par
C     -- SI units, L and W in km, input slip rates in mm/yr
      implicit none
      integer MAX_L,MAX_W,MAX_GRIDS,MAX_SEG
      parameter(MAX_L=210,MAX_W=45,MAX_GRIDS=MAX_L*MAX_W,MAX_SEG=15)

      character*80 version
C     --index parameters --
      integer i,i2,j,k,l,m,n,kl,jj,kk,idx,iii
      integer nevents

C     -- grid parameters --
      integer nsfRup
      integer n_flt,ipt(MAX_L,MAX_W,2),nz(MAX_L,2),igr
      integer ix(MAX_GRIDS,2),iz(MAX_GRIDS,2),ngrid(2),n_len(2),n_wid(MAX_L,2)
      real lon(MAX_GRIDS,2),lat(MAX_GRIDS,2),dep(MAX_GRIDS,2),str(MAX_GRIDS,2)
      real dip(MAX_GRIDS,2),rak(MAX_GRIDS,2),slp,sfL(MAX_GRIDS,2)
      real sfW(MAX_GRIDS,2)
      character*80 infile(2),name(2),subfile(2)

C     -- subfault parameters --
      real sx(4,MAX_GRIDS),sy(4,MAX_GRIDS)

C     -- length distribution parameters --
      integer iLstart(MAX_SEG,2),iLend(MAX_SEG,2),iLs,iLe,iLtaper,nL,nLen
      real Rate(MAX_SEG),Coupl(Max_SEG)

C     -- logic tree parameters --
      integer nTop,nBot,lBot
      real zSplay(3),zTrench(3)
      real zMed(3),zBot(MAX_L,3),zBot1(3),zTop(MAX_L,3),zTop1(3)
      real WtTop(3), WtBot(3),wtLen(3,1)

C     -- bottom parameters
      real dist,lat1,lat2,lon1,lon2,dep1,dep2,dist1,dist2,dbot,mindist,mindz
      real lon_w(1000,5),lat_w(1000,5),top_w(1000,5),bot_w(1000,5),vel_w(1000,5)
      real d_bot
      integer i_bot(MAX_L,5),icur,n_dep,n_w(5)
      character*80 widthfile

C     -- source parameters
      character*80 sfnm,mxfile,cnum*4,out1
      data sfnm /'Cascadia'/
      real scale_slip(MAX_GRIDS,2),ratio,DA,DA_target,amom,amag
      real slip(MAX_GRIDS,2),area,dD,aMu

C     -- scaling variability
      integer ii,n_Scale
      real aM(3),bM(3),sigM(3),WtScale(3)
C     -- aleatory parameters --
      integer ieps
      real pMag(5)

C     -- slip variability --
      integer ls(6),nAsper
      real DmaxDav,DminDav
      real pAsper

C     -- Cumulative rates --
      real Cumul_Slip(MAX_GRIDS,2)

      real Wt,Prob,rate_event
      character*80 outfile
      real PI, conv

C     -- site coordinates for joint hazard --
      real stlo, stla,atop,abot

C     ----------------------- Data block and definitions ------------------------
C     -- set the parameters for slip models --
      data zSplay   / 0.0, 0.0, 0.0/
      data zTrench  / 1.0, .5, 0.0/
      data zMed     /  5.0, 10.0, 15. /
C     data zBot     / 10.0, 15.0, 20.0/
C     data iLstart  /  1,  1,  1,  1,  1,  1, 1,  1,  1,  1,  1,  1, 1, 1, 1
C    1                 1,  1,  1,  1,  1,  1, 1,  1,  1,  1,  1,  1, 1, 1, 1
C     data iLend    / 35, 11,  9,  6, 11,  8,35, 11,  9,  6, 11,  8,35,11, 9
C    1                37, 33, 26, 19, 33, 25 /
      data name     /'slab','splay'/
      data WtTop    / 0.33, 0.34, 0.33/
      data WtBot    / 0.33, 0.34, 0.33/
C     data Rate     / 0.00390, 0.0002, 0.0004, 0.0008, 0.0002, 0.0004/
      data Rate     / .01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01/
      data pMag     /.0,  .0, 1.0,  .0, .0/
C     data pMag     /.0, .28, .44, .28, .0/
C     data pMag     /.07, .24, .38, .24, .07/
C     -- Scaling relations Goda, Skarlatoudis and Murotani --
c     data aM       /2.553, 3.685, 3.806/
c     data bM       /1.283, 1.000, 1.000/
c     data sigM     /0.241, 0.176, 0.187/
c     data WtScale  /0.333, 0.333, 0.334/
C     -- Scaling relations Papazachos, Leonard, Hanks_Bakun
      data aM       /3.402, 3.990, 3.980/
      data bM       /1.219, 1.000, 1.000/
      data sigM     /0.190, 0.176, 0.037/
      data WtScale  /0.333, 0.333, 0.334/

      n_Scale=3
      PI=acos(-1.0)
      conv=PI/180.

1     format(a80)
C     -- number of faults
      n_flt=1
C     -- read general input file --
      read(*,1) version
      read(*,1) name(1)
      read(*,1) sfnm
      print*,name(1),sfnm
C     -- we assume two source files with i_invall form --
      infile(1)='i_invall-'//name(1)
      subfile(1)='o_subfaults-'//name(1)
      read(*,*) nTop
      if(nTop .gt. 0) then
         read(*,*) (zTop1(i),i=1,nTop)
         read(*,*) (wtTop(i),i=1,nTop)
      else
         read(*,*) (wtTop(i),i=1,-nTop)
      endif

      read(*,*) nBot
      if(nBot .gt. 0) then
C        -- common bottom depth --
         lBot=1
         read(*,*) (zBot1(i),i=1,nBot)
         read(*,*) (wtBot(i),i=1,nBot)
      else if(nBot .eq. 0) then 
         lBot=-1
         widthfile='i_width-'//name(1)
C        -- read top and bottom contours --
         open(9,file=widthfile)
         read(9,*) nBot
         do j=1,nBot
            read(9,*) n_w(j),wtTop(j)
            do i=1,n_w(j)
C              read(9,*,end=9) lon_w(i,j),lat_w(i,j),top_w(i,j),bot_w(i,j),vel_w(i,j)
               read(9,*) lat_w(i,j),lon_w(i,j),bot_w(i,j)
            enddo
         enddo
9        continue
         close(9)
         n_dep=i-1
      else
         lBot=1
         read(*,*) (wtBot(i),i=1,-nBot)
      endif

      read(*,*) nLen
      read(*,*) (iLstart(i,1),i=1,nLen)
      read(*,*) (iLend(i,1) ,i=1,nLen)
      read(*,*) (wtLen(i,1) ,i=1,nLen)
      read(*,*) (Rate(i) ,i=1,nLen)
      read(*,*) (Coupl(i) ,i=1,nLen)

C     -- set segment top and bottom --
C     -- if variable zTop --
      if(nTop .lt. 0) then
         nTop=-nTop
         do i=1,nTop
            read(*,*) (zTop(j,i) ,j=1,nLen)
         enddo
      else
         do i=1,nTop
            do j=1,nLen
               zTop(j,i)=zTop1(i)
            enddo
         enddo
      endif
    
      if(nBot .lt. 0) then
         nBot=-nBot
         do i=1,nBot
            read(*,*) (zBot(j,i) ,j=1,nLen)
         enddo
      else
         do i=1,nBot
            do j=1,nLen
               zBot(j,i)=zBot1(i)
            enddo
         enddo
      endif
C     -- 

      do i=1,nLen
         iLstart(i,1)= iLstart(i,1)+1
         iLend(i,1)  = iLend(i,1)+1
         Rate(i)=Coupl(i)*Rate(i)/1000.
      enddo

      read(*,*) aMu
      read(*,*) nAsper, DmaxDav
C     read(*,*) stlo,stla
      print*, nBot,nTop,nLen
C     ----------------------------------

C     -- taper between length segments --
      iLtaper =  4
C     -- min and max displacement relative to Dav --
C     -- nAsper is 1/(fraction of apserity vs rupture) --
C     nAsper = 3
C     DmaxDav = 2.0
      DminDav = (nAsper-DmaxDav)/(nAsper-1)
      print*,'DmaxDav DminDav',DmaxDav,DminDav

C     ----------------------- Begin program -------------------------------------
C     -- read input files --
      do k=1,n_flt
         open(1,file=infile(k))
         do i=1,MAX_GRIDS
C           print*,k,i
            read(1,*,end=99) lon(i,k),lat(i,k),dep(i,k),str(i,k),dip(i,k),
     1                       rak(i,k),slp,sfL(i,k),sfW(i,k),ix(i,k),iz(i,k)
            ix(i,k)=ix(i,k)+1
            iz(i,k)=iz(i,k)+1
C           -- find rows --
            if(i .eq. 1) then
               j=1
               n_wid(j,k)=1
               icur=ix(i,k)
            else if (ix(i,k) .ne. ix(i-1,k)) then
               j=j+1
               n_wid(j,k)=1
            else
               n_wid(j,k)=n_wid(j,k)+1
            endif
            n_len(k)=j

            ipt(ix(i,k),iz(i,k),k)=i
C           print*,ix(i,k),iz(i,k)
C           -- The following assumes that the last element is the deepest --
            nz(ix(i,k),k)=iz(i,k)
            Cumul_Slip(i,k)=0.0
         enddo
99       continue
         ngrid(k)=i-1
         print*,'Fault ',k,'- read ',ngrid(k),'subfaults, along strike',n_len(k)
         close(1)
      enddo



C     -- compute termination depth for different epistemic models
C     -- only for first grid right now --
C     -- find nearest 2 points on thetermination contour and interpolate the depth
      do kk=1,nLen
C        -- loop over shallow termination/splay --
         iLs=iLstart(kk,1)
         iLe=iLend(kk,1)
         nL=iLe-iLs+1
         do j=iLs,iLe
            if(lBot .eq. 1) then
               do m=1,nBot
                  i_bot(j,m) = 1
                  do i=1,n_wid(j,1)
                     ii=ipt(j,i,1)
                     if(dep(ii,1) .le. zBot(kk,m) ) i_bot(j,m) = i
                  enddo
                  print*,"Bottom: ",j,m,i_bot(j,m),dep(ipt(j,i_bot(j,m),1),1),n_wid(j,1),zBot(j,m)
               enddo
            else if (lBot .eq. -1) then
               do m=1,nBot
                  mindist=99999.
                  do i=1,n_wid(j,1)
                     l=ipt(j,i,1)
C                    print*, l,lon(l,1),lat(l,1)
                     do k=1,n_w(m)
                        dist=sqrt(((lon(l,1)-lon_w(k,m))*sin((90.-lat(l,1))*conv))**2+
     1                       (lat(l,1)-lat_w(k,m))**2)
C                       print*,k,lon_w(k,m),lat_w(k,m),bot_w(k,m),dist
                        if(dist .lt. mindist) then
                           ii=k
                           jj=l
c                          dist2=mindist
c                          lat2=lat1
c                          lon2=lon1
c                          dep2=dep1
                           dist1=dist
                           lat1=lat_w(k,m)
                           lon1=lon_w(k,m)
                           dep1=bot_w(k,m)
                           mindist=dist
                        endif
                     enddo
                  enddo
C                 Need to find the scond nearest point here
                  if(ii .gt. 1 .and. ii .lt. n_w(m)) then
                     dist=sqrt(((lon(l,1)-lon_w(ii-1,m))*sin((90.-lat(l,1))*conv))**2+
     1                       (lat(l,1)-lat_w(ii-1,m))**2)
                     dist2=sqrt(((lon(l,1)-lon_w(ii+1,m))*sin((90.-lat(l,1))*conv))**2+
     1                       (lat(l,1)-lat_w(ii+1,m))**2)
                     if(dist .lt. dist2) then
                        dist2=dist
                        dep2=bot_w(ii-1,m)
                     else
                        dep2=bot_w(ii+1,m)
                     endif
                  else
                     i2=n_w(m)-1
                     if(ii .eq. 1) i2=2
                     dist2=sqrt(((lon(l,1)-lon_w(i2,m))*sin((90.-lat(l,1))*conv))**2+
     1                       (lat(l,1)-lat_w(i2,m))**2)
                     dep2=bot_w(i2,m)
                  endif

                  d_bot=(dist1*dep2+dist2*dep1)/(dist1+dist2)
                  mindz=9999.
                  do i=1,n_wid(j,1)
                     l=ipt(j,i,1)
                     if(abs(dep(l,1)-d_bot) .lt. mindz) then
                        mindz=abs(dep(l,1)-d_bot)
                        i_bot(j,m)=i
                     endif
                  enddo
               print*, "Bottom:",i_bot(j,m),d_bot,dist1,dist2,dep1,dep2
               enddo
            else
               print*,'lBot=',lBot,'is not defined'
               stop
            endif
         enddo
      enddo
 
C     -- compute termination points for every column --
      

C     -- compute maximum width for rupture --

C     -- termination depths completed

      open(2,file='Events.table')
      open(7,file='Events_extended.table')
      idx=index(sfnm,' ')-1
      out1='i_multimux-'//sfnm(1:idx)
      open(3,file=out1)
      write(3,'(a80)') sfnm


      nevents=0
C     -- loop over scaling relations --
      do ii=1,n_Scale
C        -- loop over length --
         do j=1,nLen
C           -- loop over shallow termination/splay --
            iLs=iLstart(j,1)
            iLe=iLend(j,1)
            nL=iLe-iLs+1
            print*,'Length=',j,iLs,iLe
            do k=1,nTop
C              -- loop over depth termination --
               print*,'Top=',k
               do l=1,nBot
C                 -- set counters for area, cumulative slip and slip x area
                  print*,'Bottom=',l
                  area=0.0
C                 -- set slip to zero everywhere
                  do i = 1,n_flt
                     do jj = 1,ngrid(i)
                        slip(jj,i)=0.0
                     enddo
                  enddo
C                 -- determine the base average slip for every element --
                  print*,'Start average slip'
                  do jj=iLs,iLe
C                    print*, j,k,l,jj,nz(jj,1)
C                    -- compute average slip for along-strike position --
                     igr=ipt(jj,i_bot(jj,l),1)
                     zMed(l)=.50*dep(igr,1)
                     zBot1(l)=dep(igr,1)
                     print*,"Row: ",jj,igr,i_bot(jj,l),zBot1(l),zMed(l),n_wid(jj,1),dep(igr,1)
                     do kl=1,i_bot(jj,l)
                        igr=ipt(jj,kl,1)
C                       -- temp fix to set depth taper
C                       -- set slip according to depth --
C                       print*,jj,kl
                        if(dep(igr,1) .lt. zTop(j,k)) then
C                          -- taper towards the trench --
C                          slip(igr,1)=1.0-(zTop(j,k)-dep(igr,1))/(zTop(j,k)-zTrench(k))
                           slip(igr,1)=0.0
C                          if(zTrench(k) .ne. 0.0) then
C                             area=area+sfL(igr,1)*sfW(igr,1)
C                          endif
C                          print*,0,igr
                        elseif(dep(igr,1) .lt. zMed(l)) then
C                          -- between trench and mid rupture --
                           area=area+sfL(igr,1)*sfW(igr,1)
                           slip(igr,1)=1.0
C                          print*,1,igr
                        elseif(dep(igr,1) .lt. zBot1(l)) then
C                          -- taper towards the bottom --
                           slip(igr,1)=1.0-(dep(igr,1)-zMed(l))/(zBot1(l)-zMed(l))
                           area=area+sfL(igr,1)*sfW(igr,1)
C                          print*,2,igr
                        else
C                          print*,3,igr,dep(igr,1),zBot1(l)
                           slip(igr,1)=0.0
                        endif
C                       -- add the splay fault --
                        if(zSplay(k) .ne. 0.0 .and. kl .le. nz(jj,2)) then
                           igr=ipt(jj,kl,2)
                           if(jj .ge.   iLstart(j,2) .and. jj .le. iLend(j,2)) then
                              area=area+sfL(igr,2)*sfW(igr,2)
                              slip(igr,2)=zSplay(k)
                           endif
                        endif
                     enddo
                  enddo
                  print*,'Done Average slip'
C                 -- end of rupture loop --
C                 -- we now have area and scaled momentum DA to work with --
C                 -- iterate over aleatory uncertainty --
                  do n=1,nAsper
                     ls(n)=iLs+(n-1)*nL/nAsper
                  enddo
                  ls(nAsper+1)=iLe
                  dD=(DmaxDav-DminDav)/iLtaper
C                 -- now include slip variability --
                  print*,'Start slip variability'
                  
                  pAsper=1./nAsper
                  do n=1,nAsper
                     DA=0.0
C                    print*,'Asperity: ',n,ls(n),ls(n+1)
                     do kk=1,n_flt
                        do m=iLStart(j,kk),iLend(j,kk)
C                          print*,m,nz(m,kk)
                           do kl=1,i_bot(m,l)
                              igr=ipt(m,kl,kk)
                              if(m .le. (ls(n)-iLtaper)) then
                                 scale_slip(igr,kk)=slip(igr,kk)*DminDav
                              elseif(m .le. ls(n)) then
                                 scale_slip(igr,kk)=slip(igr,kk)*(DmaxDav - (ls(n)-m)*dD)
                              elseif(m .gt. ls(n) .and. m .le. ls(n+1)) then
                                 scale_slip(igr,kk)=slip(igr,kk)*DmaxDav
C                                print*,'Max slip: ',n,m
                              elseif(m .lt. ls(n+1)+iLtaper) then 
                                 scale_slip(igr,kk)=slip(igr,kk)*(DmaxDav - (m-ls(n+1))*dD)
                              else
                                 scale_slip(igr,kk)=slip(igr,kk)*DminDav
                              endif
                              if(m .lt. 5) then 
                                 scale_slip(igr,kk)=m*scale_slip(igr,kk)/5
                              else if (m .gt. (n_len(1)-5)) then
                                 scale_slip(igr,kk)=(n_len(1)-m)*scale_slip(igr,kk)/5
                              endif
                              DA=DA+scale_slip(igr,kk)*sfL(igr,kk)*sfW(igr,kk)
                           enddo
                        enddo
                     enddo
C                    -- iterate over aleatory uncertainty --
                     print*,'Start aleatory uncertainty'
                     do ieps=0,0,1
C                    do ieps=-1,1,1
                        nevents=nevents+1

                        amag=aM(ii)+bM(ii)*log10(area)+sigM(ii)*ieps
                        amom=10.**(1.5*amag+9.1)
                        DA_target=amom/(aMu*1e3*1e3)
                        ratio=DA_target/DA

                        Wt=wtLen(j,1)*WtScale(ii)*WtTop(k)*WtBot(l)
C                       Prob=pMag(ieps+3)*pAsper
C                       -- no aleatory magnitude --
                        Prob=pAsper*pMag(ieps+3)
                        rate_event=Wt*Prob*Rate(j)/(DA_target/area)
                        print*,"Rate",wtLen(j,1),WtScale(ii),WtTop(k),WtBot(l),pAsper,pMag(ieps+3),Rate(j),DA_target,area
                        print*,"Total",wtLen(j,1)*WtScale(ii)*WtTop(k)*WtBot(l)*pAsper*pMag(ieps+3)*Rate(j)/(DA_target/area)
                        write(outfile,'("Event-",i1.1,"-",i2.2,"-",3(i1,
     1                                  "-"),i2.2,".par")') ii,j,k,l,n,ieps+2 
C                       print*, outfile
C                       write(2,10) outfile,amag,aL,aW,aD,eps,nL,nstepL
                        open(1,file=outfile)
                        

                        nsfRup=0
                        minDist=9999.
                        do kk=1,n_flt
                           do m=iLStart(j,kk),iLend(j,kk)
                              do kl=1,i_bot(m,1)
                                 igr=ipt(m,kl,kk)
                                 if(scale_slip(igr,kk) .ne. 0.0)  then
                                    nsfRup=nsfRup+1
                                    dist=111.*sqrt(((stlo-lon(igr,1))*sin((90.-stla)*conv))**2+(stla-lat(igr,1))**2)
                                    dist=sqrt(dep(igr,1)**2+dist**2)
                                    if(dist < minDist) then
                                       minDist=dist
                                       abot=dep(ipt(m,i_bot(m,1),kk),1)
                                       atop=dep(ipt(m,1,kk),1)
                                    endif
                                 endif
                              enddo
                           enddo
                        enddo

                        write(3,'(i5,x,a60,x,4(x,f8.1),x,e12.6)')
     1                            nsfRup,outfile,amag,DA_target/area,area,area,
     2                            rate_event
                        write(2,'("Event-",i1.1,"-",i2.2,"-",3(i1,"-"),i2.2,x,e12.6,x,f5.3,x,f15.1,2(x,i1),2(x,f5.1))') ii,j,k,l,n,
     1                                         ieps+2,rate_event ,amag,mindist,n,j,atop,abot
                        write(7,'("Event-",i1.1,"-",i2.2,"-",3(i1,"-"),i2.2,x,e12.6,x,f5.3,x,e12.6,x,f7.2,x,f10.2)') 
     1                            ii,j,k,l,n,ieps+2,rate_event,amag,amom,DA_target/area,area
C                       print*,DA,DA_target,area,ratio,amag
                        do kk=1,n_flt
                           do m=iLStart(j,kk),iLend(j,kk)
                              do kl=1,i_bot(m,1)
                                 igr=ipt(m,kl,kk)
                                 write(cnum,'(i4.4)') igr-1
C                                print*,m,kl,igr
                                 mxfile=sfnm(1:idx)//"-"//cnum//".grd"
                                 if(scale_slip(igr,kk) .ne. 0.0) then
                                    write(3,'(a60,1x,3(e12.6,1x),i4,2(x,f6.3))')
     1                                        mxfile,ratio*scale_slip(igr,kk),
     2                                        rate_event,amag,igr
                                    write(1,'(f7.3,x,i4.4,x,a80,2(f9.4,x),f4.1)') ratio*scale_slip(igr,kk),igr-1,name(kk),
     1                                                                      lon(igr,1),lat(igr,1),dep(igr,1)
                                    Cumul_Slip(igr,kk)=Cumul_Slip(igr,kk)+ratio*scale_slip(igr,kk)*rate_event
                                 endif
                              enddo
                           enddo
                        enddo
                        close(1)
                     enddo
C                   -- close Aleatory loop --
                  enddo
C                 -- close SlipVar loop --
               enddo
C              -- close Bottom loop --
            enddo
C           -- close Top loop --
         enddo
C        -- close Length loop --
      enddo
C     -- close Scaling loop --
      print*,nevents
      close(2)
      close(3)
      close(7)
      open(3,file='Rates.xy')
      do k=1,2
         do j=i,ngrid(k)
            write(3,'(f12.7,x,f11.7,x,e12.6)') lon(j,k),lat(j,k),Cumul_Slip(j,k)
         enddo
      enddo
      close(3)
10    format(a20,x,f4.2,x,f7.2,x,f6.2,x,f4.1,x,f4.1,2(x,i3))
      end
