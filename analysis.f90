      Program analyze
      implicit none
!====  Analysis of simulations for phase population studies ====
!By: Evan Thilo, Dr Eric Shea-Brown, University of Washington, 7/2008
      
      integer(kind=8) a
      integer(kind=8) histogram
      integer(kind=8) histogram2
      integer(kind=8) maxtime
      integer(kind=8) numosc
      integer(kind=8) numsims
      integer(kind=8) winsizerange
      
      parameter(a=2)
      parameter(histogram=6001)
      parameter(histogram2=141)
      parameter(maxtime=800000)
      parameter(numosc=2)
      parameter(numsims=10)
      parameter(winsizerange=10)
      
      character(len=200) :: params
      character(len=100) input_filename
      character(len=220) hstname
      character(len=220) filename
      
      double precision bursttime(numsims),combursttime(numsims)
      double precision bstmatrix(maxtime,7),c
      double precision coeffv1(numsims),coeffv2(numsims),coeffv1exp,coeffv2exp,corrwidth
      double precision covar(numsims),covarexp,covarvector,cv21(numsims),cv22(numsims)
      double precision cv21exp,cv22exp,dt,fano1(numsims),fano2(numsims)
      double precision fano1exp,fano2exp,fr1(numsims),fr2(numsims),fr1exp,fr2exp
      double precision ISI1,ISI2,mu,oscdata(maxtime,a)
      double precision osc1(maxtime),osc2(maxtime),osc1exp,osc2exp
      double precision osc1count,osc2count,osc1spiketotal,osc2spiketotal
      double precision osc1vector(maxtime),osc2vector(maxtime)
      double precision restcovar(numsims),restrho(numsims)
      double precision restvar1(numsims),restvar2(numsims)
      double precision restrhoexp(winsizerange),restvar1exp(winsizerange)
      double precision restvar2exp(winsizerange),restcovarexp(winsizerange)
      double precision bursttimeexp(winsizerange),combursttimeexp(winsizerange)
      double precision rho(numsims),rhoexp,sigma
      double precision stdvfr1(numsims),stdvfr2(numsims)
      double precision stdvcoeffv1(numsims),stdvcoeffv2(numsims)
      double precision stdvfano1(numsims),stdvfano2(numsims)
      double precision stdvrho(numsims),stdvvar1(numsims),stdvvar2(numsims)
      double precision stdvbursttime(numsims),stdvcombursttime(numsims)
      double precision stdvcovar(numsims),stdvcv21(numsims),stdvcv22(numsims)
      double precision stdvrestrho(numsims),stdvrestcovar(numsims)
      double precision stdvrestvar1(numsims),stdvrestvar2(numsims)
      double precision stdvfr1exp,stdvfr2exp,stdvcoeffv1exp,stdvcoeffv2exp
      double precision stdvfano1exp,stdvfano2exp,stdvrhoexp
      double precision stdvvar1exp,stdvvar2exp,stdvcovarexp
      double precision stdvcv21exp,stdvcv22exp
      double precision stdvrestrhoexp(winsizerange),stdvrestcovarexp(winsizerange)
      double precision stdvrestvar1exp(winsizerange),stdvrestvar2exp(winsizerange)
      double precision stdvbursttimeexp(winsizerange)
      double precision stdvcombursttimeexp(winsizerange)
      double precision tmax,tpad
      double precision var1vector(maxtime),var2vector(maxtime)
      double precision var1(numsims),var2(numsims)
      double precision var1exp,var2exp
      double precision winend,winsize,winstart,winsteplength
      
      integer(kind=8) auto1(histogram),auto2(histogram),correlationloop,distribution
      integer(kind=8) model_choice,modifier_choice
      integer(kind=8) runtime,runtimeend(8),runtimestart(8),runtimetotal
      integer(kind=8) spkcount1,spkcount2,xcor1(histogram),xcor2(histogram)
      integer(kind=8) burstxcor(histogram2),winsizevector(winsizerange)
      integer(kind=8) pad
      
      integer(kind=8) Aprobe,bstranks(maxtime,4),cvprobe1,cvprobe2,d,g,h
      integer(kind=8) i,ilast,ilastold,ISIcount,j,jlast,jlastold,k,x,z
      integer(kind=8) lastspikerank1,lastspikerank2
      integer(kind=8) p,probe1,probe2,r1,r2,ref,r
      integer(kind=8) tpadstartrank1,tpadstartrank2
      integer(kind=8) u,w,winrank,wintotal,winstep,Xprobe
      
      logical nospikes1,nospikes2,killall
      
      double precision simstart,simend

!      write (*,990) winsize,fr1exp,stdvfr1exp,fr2exp,stdvfr2exp,&
!                   &coeffv1exp,stdvcoeffv1exp,coeffv2exp,stdvcoeffv2exp,&
!                   &fano1exp,stdvfano1exp,fano2exp,stdvfano2exp,&
!                   &rhoexp,stdvrhoexp,&
!                   &var1exp,stdvvar1exp,var2exp,stdvvar2exp,&
!                   &covarexp,stdvcovarexp,&
!                   &cv21exp,stdvcv21exp,cv22exp,stdvcv22exp,&
!                   &restrhoexp(k),stdvrestrhoexp(k),restcovarexp(k),stdvrestcovarexp(k),&
!                   &restvar1exp(k),stdvrestvar1exp(k),restvar2exp(k),stdvrestvar2exp(k),&
!                   &bursttimeexp(k),stdvbursttimeexp(k),&
!                   &combursttimeexp(k),stdvcombursttime(k)

  989 format(7F11.2)
  990 format(F5.0,   F12.4,   F9.4,    F12.4,   F9.4,&
            &F9.4,   F12.7,   F9.4,    F12.7,&
            &F11.4,  F12.7,   F11.4,   F12.7,&
            &F10.4,  F12.7,&
            &F10.4,  F12.7,   F10.4,   F12.7,&
            &F10.4,  F12.7,&
            &F10.4,  F12.7,   F10.4,   F12.7,&
            &F10.4,  F12.7,   F10.4,   F12.7,&
            &F10.4,  F12.7,   F10.4,   F12.7,&
            &F10.4,  F12.7,&
            &F10.4,  F12.7)
  991 format(' mu=',F5.2,' sigma=',F5.2,' c=',F5.2,' tmax=',F8.0,' simsize=',F6.0)
  992 format('runtime= ',I3,'   wincount= ',I2,'   mu= ',F7.4,'   sig= ',F7.4,'   c= ',F7.4,'   dt= ',F7.4,'   tmax= ',F8.0)
  993 format(F20.4)
  994 format(F1.0,F11.2)
  995 format(F20.4)
  996 format(4I10)
  997 format(I10)
  998 format(F12.3,F12.3,I6,I6,F7.1,I6)
  999 format(2F8.3)         
     

!****************************************************************************
!*********************Call Arguments, Assign Parameters**********************
!****************************************************************************
      
      call getarg(1,params)
      read(params,*),mu
      call getarg(2,params)
      read(params,*),sigma
      call getarg(3,params)
      read(params,*),c
      call getarg(4,params)
      read(params,*),model_choice
      call getarg(5,params)
      read(params,*),modifier_choice
      call getarg(6,params)
      read(params,*),distribution
      call getarg(7,params)
      read(params,*),correlationloop
      call getarg(8,params)
      read(params,*),runtime
      call getarg(9,params)
      read(params,*),dt
      call getarg(10,params)
      read(params,*),tmax
!      call getarg(10,params)
!      read(params,*),corrwidth
      corrwidth=(histogram-1)/2
      call getarg(11,params)
      read(params,*),input_filename
      
!winsize range is the number of different window sizes to try out
      winsizevector=[1,2,4,8,16,32,64,128,256,1000]
!      winsizevector=[1,2,4,8,16,32]
!pad the window so it doesn't start at the beginning
      tpad = 200.d0
      
!****************************************************************************
!**********************Generate Spike-time Vectors***************************
!****************************************************************************
!Pull the data out of the file and into matrix oscdata

!      open(unit=10,file="oscdata.dat",status='replace')
!      open(unit=11,file="oscvectors.dat",status='replace')
!      open(unit=12,file="bstmatrix.dat",status='replace')
!      open(unit=13,file="bstranks.dat",status='replace')
      
      if (distribution .eq. 3) then
         filename=input_filename
         open(unit=9,file=filename)
         if (correlationloop .gt. 0) then
            hstname='hst_' // input_filename
            open(unit=13,file=hstname,status='replace')
         end if
         if (runtime .eq. 1) then
            call date_and_time(values=runtimestart)
         end if
      else
         filename='code/corr/data/' // input_filename 
         open(unit=9,file=filename)
         if (correlationloop .gt. 0) then
            hstname='code/corr/hst_' // input_filename
            open(unit=13,file=hstname,status='replace')
         end if
         if (runtime .eq. 1) then
            call date_and_time(values=runtimestart)
         end if
      end if
      
      do 100 i = 1,maxtime
      read(9,994,end=40) (oscdata(i,j),j=1,2)
  100 end do
  
!these keep track of array rank in osc1 and osc2
  40  r1=1
      r2=1
      Aprobe=0
      Xprobe=0
      tpadstartrank1=0
      tpadstartrank2=0
      lastspikerank1=0
      lastspikerank2=0

      print *,'tmax=',tmax,'     tpad=',tpad
!separate the spikes into 2 separate lists osc1, osc2
      do 101 i=1,size(oscdata(:,1))
         if (oscdata(i,1) .le. 1.5) then
            osc1(r1)= oscdata(i,2)
            if ((osc1(r1) .ge. tpad) .and. (tpadstartrank1 .eq. 0)) then
               tpadstartrank1=r1
            end if
            if (((osc1(r1) .ge. tmax-tpad) .and. (lastspikerank1 .eq. 0)).or.(osc1(r1).eq.0)) then
            
               lastspikerank1=r1-1
            end if
            r1=r1+1
         else
            osc2(r2)= oscdata(i,2)
            if ((osc2(r2) .ge. tpad) .and. (tpadstartrank2 .eq. 0)) then
               tpadstartrank2=r2
            end if
            if (((osc2(r2) .ge. tmax-tpad) .and. (lastspikerank2 .eq. 0)).or.(osc2(r2).eq.0)) then
               lastspikerank2=r2-1
            end if
            r2=r2+1
         end if
         if (oscdata(i,1) .eq. 0) exit
  101 end do
!      print *,'lastspikerank1=',lastspikerank1,'    lastspikerank2=',lastspikerank2
!      write(11,995) osc1
!****************************************************************************
!**************************Auto/Cross Correlations***************************
!****************************************************************************
      if (correlationloop .eq. 1) then 
         do 126 i = 1,size(auto1)
            auto1(i)=0
            xcor1(i)=0
            auto2(i)=0
            xcor2(i)=0
  126    end do
         do 103 ref = tpadstartrank1,lastspikerank1
!Autocor of osc1
!Set the starting corr position
            do while (abs(osc1(ref)-osc1(Aprobe)) .gt. corrwidth)
               Aprobe=Aprobe+1
            end do
!walk down spiketimes, record the ISI
            do 104 k = Aprobe,size(osc1)
               if (abs(osc1(ref)-osc1(k)) .lt. corrwidth) then
                  d=dnint(osc1(ref)-osc1(k))+corrwidth
                  auto1(d)=auto1(d)+1
               else
                  exit
               end if
  104       end do
!Xcor, osc1 as ref
            do while (abs(osc1(ref)-osc2(Xprobe)) .gt. corrwidth)
               Xprobe=Xprobe+1
            end do
            
            do 105 k = Xprobe, size(osc2)
               if (abs(osc1(ref)-osc2(k)) .lt. corrwidth) then
                  d=dnint(osc1(ref)-osc2(k))+corrwidth
                  xcor1(d)=xcor1(d)+1
               else
                  exit
               end if
  105       end do

  103    end do
  
         Aprobe=0
         Xprobe=0
         do 106 ref = tpadstartrank2,lastspikerank2
!Autocor of osc2
!Set the starting corr position
            do while (abs(osc2(ref)-osc2(Aprobe)) .gt. corrwidth)
               Aprobe=Aprobe+1
            end do
            
!walk down spiketimes, record the ISI
            do 107 k = Aprobe,size(osc2)
               if (abs(osc2(ref)-osc2(k)) .lt. corrwidth) then
                  d=dnint(osc2(ref)-osc2(k))+corrwidth
                  auto2(d)=&
                     &auto2(d)+1
               else
                  exit
               end if
  107       end do
!Xcor, osc2 as ref
!Set the starting corr position
            do while (abs(osc2(ref)-osc1(Xprobe)) .gt. corrwidth)
               Xprobe=Xprobe+1
            end do
            
!walk down spiketimes, record the ISI
            do 108 k = Xprobe, size(osc1)
               if (abs(osc2(ref)-osc1(k)) .lt. corrwidth) then
                  d=dnint(osc2(ref)-osc1(k))+corrwidth
                  xcor2(d)=&
                     &xcor2(d)+1
               else
                  exit
               end if
  108       end do
  106    end do
         write(13,*) 'auto1, xcor1, auto2, xcor2'
         do 109 i=1,size(auto1)
            write(13,996) auto1(i),xcor1(i),auto2(i),xcor2(i)
 
  109    end do
      end if


!****************************************************************************
!*********************Burst Correlations, restricted rho*********************
!****************************************************************************
      if ((correlationloop .eq. 2).or.(modifier_choice .eq. 5)) then
         do 182 k=1,maxtime
            bstmatrix(k,1)=osc1(k)
            bstmatrix(k,2)=osc2(k)
            do 185 i=3,7
               bstmatrix(k,i)=0.0d0
  185       end do
!                     write(12,989) bstmatrix(k,:)
  182    end do
  
         
         do 181 k=1,size(winsizevector)
            winsize=winsizevector(k)
            winsteplength=winsize
            pad=(histogram2-1)/2
            tpadstartrank1=0
            tpadstartrank2=0
            lastspikerank1=0
            lastspikerank2=0
            do 186 p=1,numsims
               winrank=0
               bursttime(p)=0
               
               do 192 h=1,maxtime
                  do 193 i=1,4
                     bstranks(h,i)=0
  193             end do
  192          end do
               
               combursttime(p)=0
               simstart=(p-1)*(tmax/numsims)
               simend=(p)*(tmax/numsims)
               
               !clear the window vectors
               if (modifier_choice.eq.5) then
                  do 191 i=1,maxtime
                     osc1vector(i)=0.d0
                     osc2vector(i)=0.d0
  191             end do
               end if
               !DEFINE WHERE THE 1ST AND LAST RANKS ARE FOR BOTH OSCs
               do while ((osc1(tpadstartrank1).lt.simstart+tpad).and.&
                         &(osc1(tpadstartrank1+1).gt.osc1(tpadstartrank1)))
                  tpadstartrank1=tpadstartrank1+1
               end do
               do while ((osc2(tpadstartrank2).lt.simstart+tpad).and.&
                         &(osc2(tpadstartrank2+1).gt.osc2(tpadstartrank2)))
                  tpadstartrank2=tpadstartrank2+1
               end do
               
               do while (osc1(lastspikerank1+1).lt.simend)
                  if (osc1(lastspikerank1+1).lt.osc1(lastspikerank1)) exit 
                  lastspikerank1=lastspikerank1+1
               end do
               do while (osc2(lastspikerank2+1).lt.simend)
                  if (osc2(lastspikerank2+1).lt.osc2(lastspikerank2)) exit
                  lastspikerank2=lastspikerank2+1
               end do
               
               
               if (correlationloop.eq.2) then
                  do 102 h = 1,size(burstxcor)
                     burstxcor(h)=0
  102             end do
               end if
               
               
               
!               write(*,911),k,p
!  911          format('winsize=',I3,'   simulation number=',I3)
!               write(*,910),simstart,simend
!  910          format('simstart=',F11.2,'   simend=',F11.2)
!               write(*,909),tpadstartrank1,osc1(tpadstartrank1),tpadstartrank2,osc2(tpadstartrank2)
!  909          format('1stosc1=',I8,'  osc1(1st)=',F11.2,'  1stosc2=',I8,'  osc2(1st)=',F11.2)
!               write(*,908),lastspikerank1,osc1(lastspikerank1),lastspikerank2,osc2(lastspikerank2)
!  908          format('lastosc1=',I8,'  osc1(last)=',F11.2,'  lastosc2=',I8,'  osc2(last)=',F11.2)
 

!make a list of the start/finish ranks of all bursts for osc1
               i=tpadstartrank1
               ilast=-1
               g=0
               killall=.false.
               if (osc1(1).eq.0) then
                  killall=.true.
               end if
               do while ((ilast.lt. lastspikerank1-1).and.(killall.eqv..false.))
                  if (ilast.ne.-1)then
                     i=ilast
                  end if
                  ilastold=ilast
                  call findaburst(i,ilast,osc1,lastspikerank1,maxtime,killall)
                  if (ilast.eq.ilastold) exit
                  g=g+1
                  bstranks(g,1)=i
                  bstranks(g,2)=ilast
                  bursttime(p)=bursttime(p)+((osc1(ilast)-osc1(i)))
!                  write(*,923) osc1(i),i,osc1(ilast),ilast,bursttime(p)
!   923            format(' osc1(i)= ',F11.3,' i= ',I8,' osc1(ilast)= '&
!                  &,F11.3,' ilast= ',I8,' bursttime(p)= ',F11.3)
               end do
               bursttime(p)=bursttime(p)/(simend-simstart-tpad)

!make a list of the start/finish ranks of all bursts for osc2
               j=tpadstartrank2
               jlast=-1
               h=0
               killall=.false.
               if (osc2(1).eq.0) then
                  killall=.true.
               end if
               do while ((jlast.lt. lastspikerank2-1).and.(killall.eqv..false.))
                  if (jlast.ne.-1)then
                     j=jlast
                  end if
                  jlastold=jlast
                  call findaburst(j,jlast,osc2,lastspikerank2,maxtime,killall)
                  if (jlast.eq.jlastold) exit
                  h=h+1
                  bstranks(h,3)=j
                  bstranks(h,4)=jlast
               end do


!               do 149 r=1,11000
!                  write(13,833),bstranks(r,:),bstmatrix(r,:)
!   149         end do
!   833            format(4I5,7F11.3)


!walk down the lists and analyze common bursting regions               
               z=1
               x=1
               do 150 z=1,g
                  do 151 x=1,h
                     i=bstranks(z,1)
                     ilast=bstranks(z,2)
                     j=bstranks(x,3)
                     jlast=bstranks(x,4)
                     if (((osc1(i).le.osc2(j)).and.(osc1(ilast).gt.osc2(j))).or.&
                        &((osc2(j).lt.osc1(i)).and.(osc2(jlast).ge.osc1(i)))) then
                        call winanalysis(osc1vector,osc2vector,winrank,&
                              &i,ilast,j,jlast,bstmatrix,winsize,maxtime,&
                              &winsteplength,dt,simend)
                        combursttime(p)=combursttime(p)+(min(osc1(ilast),osc2(jlast))-max(osc1(i),osc2(j)))
                     end if
   151            end do
   150         end do
               combursttime(p)=combursttime(p)/(simend-simstart-tpad)
         
               if (modifier_choice.eq.5) then
                  do 183 g=1,30000
!                     write(12,989) bstmatrix(g,:)
 183              end do
            
!this sums all the spikes that we count using windows (because windows overlap, spikes are counted multiple times)
                  osc1spiketotal = 0
                  osc2spiketotal = 0
                  osc1exp        = 0.0d0
                  osc2exp        = 0.0d0
                  restvar1(p)    = 0.0d0
                  restvar2(p)    = 0.0d0
                  restcovar(p)   = 0.0d0
                  restrho(p)     = 0.0d0
              
                  do 177 i = 1,winrank
                     if (osc1vector(i).gt.0) then
                        !write(*,999),osc1vector(i),osc2vector(i)
                     end if
                     osc1spiketotal=osc1spiketotal+osc1vector(i)
                     osc2spiketotal=osc2spiketotal+osc2vector(i)
  177             end do

!                  print *,'osc1spiketotal=',osc1spiketotal,'   osc2spiketotal=',osc2spiketotal
                  if ((osc1spiketotal .gt. 2).and.(osc2spiketotal.gt.2)) then
                     osc1exp = osc1spiketotal/winrank
                     osc2exp = osc2spiketotal/winrank
!                  print *,'osc1spiketotal=',osc1spiketotal,'   osc2spiketotal=',osc2spiketotal
!                  print *,'winrank=',winrank
!                  print *,'osc1exp=',osc1exp,'   osc2exp=',osc2exp
               
                     do 178 i = 1,winrank
                        restvar1(p)  = restvar1(p) + (osc1vector(i)-osc1exp)**2
                        restvar2(p)  = restvar2(p) + (osc2vector(i)-osc2exp)**2
  178                end do
  
                     restvar1(p)     = restvar1(p)/winrank
                     restvar2(p)     = restvar2(p)/winrank
            
                     do 179 i = 1,winrank
                        restcovar(p) = restcovar(p) + (osc1vector(i)-osc1exp)*(osc2vector(i)-osc2exp)
  179                end do
  
                     restcovar(p)    = restcovar(p)/winrank
                     
                     if ((restvar1(p).ne.0).and.(restvar2(p).ne.0)) then
                        restrho(p)      = restcovar(p)/sqrt(restvar1(p)*restvar2(p))
                     end if
!                     print *, 'winrank=',winrank
!                     print *, 'restrho=',restrho(p)

  
                  end if
               end if
  186       end do
         bursttimeexp(k)   = sum(bursttime)   / numsims
         combursttimeexp(k)= sum(combursttime)/ numsims
         restrhoexp(k)     = sum(restrho)     / numsims
         restvar1exp(k)    = sum(restvar1)    / numsims
         restvar2exp(k)    = sum(restvar2)    / numsims
         restcovarexp(k)   = sum(restcovar)   / numsims
!         print *,'windowsize=',winsizevector(k)
!         print*,'restvar1exp=',restvar1exp(k),'   restvar2exp=',restvar2exp(k)
         do 187 j=1,numsims
            stdvrestrho(j)     = ((restrho(j)-restrhoexp(k))**2)
            stdvrestvar1(j)    = ((restvar1(j)-restvar1exp(k))**2)
            stdvrestvar2(j)    = ((restvar2(j)-restvar2exp(k))**2)
            stdvrestcovar(j)   = ((restcovar(j)-restcovarexp(k))**2)
            stdvbursttime(j)   = ((bursttime(j)-bursttimeexp(k))**2)
            stdvcombursttime(j)= ((combursttime(j)-combursttimeexp(k))**2)
  187    end do
  
         stdvrestrhoexp(k)     = sqrt(sum(stdvrestrho)     /numsims) / sqrt(dble(numsims))
         stdvrestvar1exp(k)    = sqrt(sum(stdvrestvar1)    /numsims) / sqrt(dble(numsims))
         stdvrestvar2exp(k)    = sqrt(sum(stdvrestvar2)    /numsims) / sqrt(dble(numsims))
         stdvrestcovarexp(k)   = sqrt(sum(stdvrestcovar)   /numsims) / sqrt(dble(numsims))
         stdvbursttimeexp(k)   = sqrt(sum(stdvbursttime)   /numsims) / sqrt(dble(numsims))
         stdvcombursttimeexp(k)= sqrt(sum(stdvcombursttime)/numsims) / sqrt(dble(numsims))
         
         bursttimeexp(k)=bursttimeexp(k)
         stdvbursttimeexp(k)=stdvbursttimeexp(k)
         
         combursttimeexp(k)=combursttimeexp(k)
         stdvcombursttimeexp(k)=stdvcombursttimeexp(k)
  
  181    end do
      end if
      
!***********************************************************************************
!************************Individual Simulation Statistics***************************
!***********************************************************************************

      write(*,991) mu,sigma,c,tmax,tmax/numsims
      print *, '+/- designates the standard error of the mean'
      print *, 'win_______&
               &fr1_______+/-_________fr2______+/-_____&
               &coeffv1_____+/-_____coeffv2______+/-________&
               &fano1_______+/-_________fano2______+/-_________&
               &rho_______+/-________&
               &var1_______+/-________var2________+/-_______&
               &covar______+/-________&
               &cv21_______+/-________cv22________+/-______&
               &restrho_____+/-_____restcovar_____+/-______&
               &restvar1_____+/-______restvar2_____+/-____&
               &bursttime_____+/-___combursttime_____+/-______'
      
!Loop over all the possible window sizes
      winsteplength=.25
      do 110 k=1,size(winsizevector)
         g=1
         h=1
         winsize = winsizevector(k)
         do 111 j=1,numsims
            simstart=(j-1)*(tmax/numsims)
            simend=(j)*(tmax/numsims)
              
            !            print *,'simstart=',simstart,'simend=',simend
            
            wintotal = dnint(((simend-simstart)-tpad)/(.25*winsize))-3
            do while (((3+wintotal)*(.25*winsize)) .gt. ((simend-simstart)-tpad))
               wintotal=wintotal-1
            end do
            
!            print *, 'winsize=',winsize,'wintotal=',wintotal
            do 112 i=1,maxtime
               osc1vector(i)=0
               osc2vector(i)=0
  112       end do

            winstart = (.25*(1-1))*winsize + tpad + simstart
!            print *, 'winstart=',winstart
            
!count the number of spikes within each window range  
            probe1=1
            probe2=1
             
            do 113 winstep=1,wintotal
               winstart    = (.25*(winstep-1))*winsize + tpad + simstart
               winend      = winstart + winsize
               osc1count   = 0.0d0
               osc2count   = 0.0d0
               
               do while (osc1(probe1) .gt. winstart)
                  if (probe1 .eq. 1) exit
                  probe1=probe1-1
               end do
               if (osc1(probe1) .lt. osc1(probe1-1)) then
                  probe1=probe1-1
               end if
               do while (osc1(probe1) .le. winstart)
                  if (osc1(probe1) .le. osc1(probe1-1)) exit
                  probe1=probe1+1
               end do
               do while (osc1(probe1) .le. winend)
                  if (osc1(probe1) .le. osc1(probe1-1)) exit
                  osc1count=osc1count+1
                  probe1=probe1+1
               end do
               

               do while (osc2(probe2) .gt. winstart)
                  if (probe2 .eq. 1) exit
                  probe2=probe2-1
               end do
               if (osc2(probe2) .lt. osc2(probe2-1)) then
                  probe2=probe2-1
               end if
               do while (osc2(probe2) .le. winstart)
                  if (osc2(probe2) .le. osc2(probe2-1)) exit
                  probe2=probe2+1
               end do
               do while (osc2(probe2) .le. winend)
                  if (osc2(probe2) .le. osc2(probe2-1)) exit
                  osc2count=osc2count+1
                  probe2=probe2+1
               end do
            
!write the total spike count for that window size into vector entries
               osc1vector(winstep)=osc1count
               osc2vector(winstep)=osc2count
!               write(*,998) winstart,winend,winstep,wintotal,winsize,k
  113       end do
!            print *,'hello!'
!            write(12,999) osc1vector,osc2vector

!calculate firing rate (using the number of ACTUAL spikes inside one simulation)
            osc1spiketotal = 0
            osc2spiketotal = 0
            do while (osc1(g) .lt. simstart+tpad)
               g=g+1
               if (osc1(g) .le. osc1(g-1)) exit
            end do
            do while ((osc1(g) .gt. simstart+tpad) .and. (osc1(g) .le. winend))
               osc1spiketotal=osc1spiketotal+1
               g=g+1
            end do
            
            do while (osc2(h) .lt. simstart+tpad)
               h=h+1
               if (osc2(h) .le. osc2(h-1)) exit
            end do
            do while ((osc2(h) .gt. simstart+tpad) .and. (osc2(h) .le. winend))
               osc2spiketotal=osc2spiketotal+1
               h=h+1
            end do
            nospikes1=.false.
            if (osc1spiketotal .le. 2) then
               nospikes1=.true.
            end if
            nospikes2=.false.
            if (osc2spiketotal .le. 2) then
               nospikes2=.true.
            end if
!            print *,nospikes1,nospikes2
!            print *, 'osc1total= ',osc1spiketotal,'osc2total= ',osc2spiketotal,'SIMSIZE',winend-simstart
            fr1(j) = 1000*osc1spiketotal/((winend-simstart)-tpad)
            fr2(j) = 1000*osc2spiketotal/((winend-simstart)-tpad)
!this sums all the spikes that we count using windows (because windows overlap, spikes are counted multiple times)
            osc1spiketotal = 0
            osc2spiketotal = 0
            osc1exp = 0.0d0
            osc2exp = 0.0d0
            if (nospikes1 .eqv. .false.) then
               do 114 i = 1,maxtime
                  osc1spiketotal=osc1spiketotal+osc1vector(i)
  114          end do
               osc1exp = osc1spiketotal/wintotal
            end if
            if (nospikes2 .eqv. .false.) then
               do 115 i = 1,maxtime 
                  osc2spiketotal=osc2spiketotal+osc2vector(i)
  115          end do
               osc2exp = osc2spiketotal/wintotal
            end if
            
!            print *,osc1exp,osc2exp
!            print *,osc1spiketotal,osc2spiketotal

!***********************************************************************************
!**********************Compiled Cross-Simulation Statistics*************************
!***********************************************************************************
            var1(j)    = 0.0d0
            var2(j)    = 0.0d0
            covar(j)   = 0.0d0
            coeffv1(j) = 0.0d0
            coeffv2(j) = 0.0d0
            rho(j)     = 0.0d0
            fano1(j)   = 0.0d0
            fano2(j)   = 0.0d0
            cv21(j)    = 0.0d0
            cv22(j)    = 0.0d0
            
!Make the independent statistics
            if (nospikes1 .eqv. .false.) then
               do 116 i = 1,wintotal
                  var1(j)  = var1(j) + (osc1vector(i)-osc1exp)**2
  116          end do
               var1(j)    = var1(j)/wintotal
               coeffv1(j) = sqrt(var1(j))/osc1exp
               fano1(j)   = var1(j)/osc1exp
! cv2 calculations
               cvprobe1=1
               ISIcount=0
               do while (osc1(cvprobe1-1) .lt. simstart+tpad)
                  cvprobe1=cvprobe1+1
               end do
               do while ((osc1(cvprobe1-1) .gt. simstart+tpad) .and. (osc1(cvprobe1+1) .le. simend))
                  ISI1=osc1(cvprobe1+1)-osc1(cvprobe1)
                  ISI2=osc1(cvprobe1)-osc1(cvprobe1-1)
                  cv21(j)=cv21(j) + (2*abs(ISI1-ISI2))/(ISI1+ISI2)
                  cvprobe1=cvprobe1+1
                  ISIcount=ISIcount+1
               end do
               cv21(j)=cv21(j)/ISIcount
            end if

            if (nospikes2 .eqv. .false.) then
               do 117 i = 1,wintotal
                  var2(j)  = var2(j) + (osc2vector(i)-osc2exp)**2
  117          end do
               var2(j)    = var2(j)/wintotal
               coeffv2(j) = sqrt(var2(j))/osc2exp
               fano2(j)   = var2(j)/osc2exp
! cv2 calculations
               cvprobe2=1
               ISIcount=0
               do while (osc2(cvprobe2-1) .lt. simstart+tpad)
                  cvprobe2=cvprobe2+1
               end do
               do while ((osc2(cvprobe2-1) .gt. simstart+tpad) .and. (osc2(cvprobe2+1) .le. simend))
                  ISI1=osc2(cvprobe2+1)-osc2(cvprobe2)
                  ISI2=osc2(cvprobe2)-osc2(cvprobe2-1)
                  cv22(j)=cv22(j) + (2*abs(ISI1-ISI2))/(ISI1+ISI2)
                  cvprobe2=cvprobe2+1
                  ISIcount = ISIcount+1
               end do
               cv22(j) = cv22(j)/ISIcount
            end if

!Interdependent statistics
            if ((nospikes1 .eqv. .false.) .and. (nospikes2 .eqv. .false.)) then
               do 118 i = 1,wintotal
                  covar(j) = covar(j) + (osc1vector(i)-osc1exp)*(osc2vector(i)-osc2exp)
  118          end do
               covar(j)   = covar(j)/wintotal
               rho(j)     = covar(j)/sqrt(var1(j)*var2(j))
            end if
            
  111    end do
         
!         print *,'hello!'
  
!after all sims have been completed for 1 window size (loop 333), you'll have stat vectors each of length(numsims) for the winsize

         fr1exp     = sum(fr1)     / numsims
         fr2exp     = sum(fr2)     / numsims
         coeffv1exp = sum(coeffv1) / numsims
         coeffv2exp = sum(coeffv2) / numsims
         fano1exp   = sum(fano1)   / numsims
         fano2exp   = sum(fano2)   / numsims
         rhoexp     = sum(rho)     / numsims
         var1exp    = sum(var1)    / numsims
         var2exp    = sum(var2)    / numsims
         covarexp   = sum(covar)   / numsims
         cv21exp    = sum(cv21)    / numsims
         cv22exp    = sum(cv22)    / numsims

         do 119 j=1,numsims
            stdvfr1(j)     = ((fr1(j)-fr1exp)**2)
            stdvfr2(j)     = ((fr2(j)-fr2exp)**2)
            stdvcoeffv1(j) = ((coeffv1(j)-coeffv1exp)**2)
            stdvcoeffv2(j) = ((coeffv2(j)-coeffv2exp)**2)
            stdvfano1(j)   = ((fano1(j)-fano1exp)**2)
            stdvfano2(j)   = ((fano2(j)-fano2exp)**2)
            stdvrho(j)     = ((rho(j)-rhoexp)**2)
            stdvvar1(j)    = ((var1(j)-var1exp)**2)
            stdvvar2(j)    = ((var2(j)-var2exp)**2)
            stdvcovar(j)   = ((covar(j)-covarexp)**2)
            stdvcv21(j)    = ((cv21(j)-cv21exp)**2)
            stdvcv22(j)    = ((cv22(j)-cv22exp)**2)
  119    end do
         
         stdvfr1exp     = sqrt(sum(stdvfr1)     /numsims) / sqrt(dble(numsims))
         stdvfr2exp     = sqrt(sum(stdvfr2)     /numsims) / sqrt(dble(numsims))
         stdvcoeffv1exp = sqrt(sum(stdvcoeffv1) /numsims) / sqrt(dble(numsims))
         stdvcoeffv2exp = sqrt(sum(stdvcoeffv2) /numsims) / sqrt(dble(numsims))
         stdvfano1exp   = sqrt(sum(stdvfano1)   /numsims) / sqrt(dble(numsims))
         stdvfano2exp   = sqrt(sum(stdvfano2)   /numsims) / sqrt(dble(numsims))
         stdvrhoexp     = sqrt(sum(stdvrho)     /numsims) / sqrt(dble(numsims))
         stdvvar1exp    = sqrt(sum(stdvvar1)    /numsims) / sqrt(dble(numsims))
         stdvvar2exp    = sqrt(sum(stdvvar2)    /numsims) / sqrt(dble(numsims))
         stdvcovarexp   = sqrt(sum(stdvcovar)   /numsims) / sqrt(dble(numsims))
         stdvcv21exp    = sqrt(sum(stdvcv21)    /numsims) / sqrt(dble(numsims))
         stdvcv22exp    = sqrt(sum(stdvcv22)    /numsims) / sqrt(dble(numsims))
         
         write (*,990) winsize,fr1exp,stdvfr1exp,fr2exp,stdvfr2exp,&
                      &coeffv1exp,stdvcoeffv1exp,coeffv2exp,stdvcoeffv2exp,&
                      &fano1exp,stdvfano1exp,fano2exp,stdvfano2exp,&
                      &rhoexp,stdvrhoexp,&
                      &var1exp,stdvvar1exp,var2exp,stdvvar2exp,&
                      &covarexp,stdvcovarexp,&
                      &cv21exp,stdvcv21exp,cv22exp,stdvcv22exp,&
                      &restrhoexp(k),stdvrestrhoexp(k),restcovarexp(k),stdvrestcovarexp(k),&
                      &restvar1exp(k),stdvrestvar1exp(k),restvar2exp(k),stdvrestvar2exp(k),&
                      &bursttimeexp(k),stdvbursttimeexp(k),&
                      &combursttimeexp(k),stdvcombursttime(k)
                         
  110 end do
      if (distribution .eq. 3) then
         if (runtime .eq. 1) then
            filename='data/analysisruntime.dat'
            open(unit=10,file=filename,access='append')
         end if
      else
         if (runtime .eq. 1) then
            filename='code/corr/analysisruntime.dat'
            open(unit=10,file=filename,access='append')
         end if
      end if
  
      if (runtime .eq. 1) then
         call date_and_time(values=runtimeend)
         runtimetotal=(runtimeend(7)+runtimeend(6)*60+runtimeend(5)*&
            &60**2+runtimeend(3)*24*60**2) - (runtimestart(7)+runtimestart(6)*&
            &60+runtimestart(5)*60**2+runtimestart(3)*24*60**2)
         write(10,992) runtimetotal,size(winsizevector),mu,sigma,c,dt,tmax
      end if
      stop
      end program analyze










!************************************************************************************
!***********************************Subroutines**************************************
!************************************************************************************

!******probe for a burst, return the first rank and last rank of the burst***********
      subroutine findaburst(k,klast,osc,lastspikerank,maxtime,killall)  
      integer(kind=8) k,n,m,lastspikerank,klast
      integer(kind=8) maxtime
      double precision osc(maxtime)
      logical killall
!      print *,'burst probe'
      
      do while (k .le. lastspikerank-1)
         m=k+1
         n=1
!         print *, 'm=',m-1,'   time=',osc(m-1)
         do while ((osc(m)-osc(m-1) .le. 20) .and. (m .le. lastspikerank-1))
!            print *, 'm=',m,'   time=',osc(m)
            m=m+1
            n=n+1
         end do
         if (n .ge. 10) exit
         k=k+1
         if (k .eq. lastspikerank-1) then
            killall=.true.
         end if
!         print *, 'no burst'
      end do
      klast=k+n-1
!      print *,'klast=',klast,'   k=',k,'   n=',n,'   m=',m,'   killall=',killall,'   lastspikerank=',lastspikerank
      return
      end subroutine findaburst

!*****************************Xcorr of common burst regions****************************
      subroutine bstXcor(burstxcor,i,ilast,j,jlast,bstmatrix,pad,maxtime,histogram2)
      
      integer(kind=8) histogram2
      integer(kind=8) pad,burstxcor(histogram2)
      integer(kind=8) ref,Xprobe
      integer(kind=8) d,k
      integer(kind=8) i,ilast,j,jlast
      integer(kind=8) maxtime
      double precision bstmatrix(maxtime,7)
!      print *,'Xcor'
      ref=i
      Xprobe=j
!set the starting position of the reference to be +70ms ahead of the beginnings of both bursts
      do while (((bstmatrix(ref,1).le.bstmatrix(i,1)+dble(pad))&
               &.or.(bstmatrix(ref,1).le.bstmatrix(j,2)+dble(pad)))&
               &.and.(bstmatrix(ref,1).le.bstmatrix(ilast,1)-dble(pad)))
         ref=ref+1
      end do
!      print *, 'START OF COMMON REGION+PAD=',bstmatrix(ref),' lastspike=',bstmatrix(ilast)
!if the bursts overlap, you just found the start of the overlap + the 70ms pad.
      do while ((bstmatrix(ref,1).le.bstmatrix(ilast,1)-dble(pad))&
               &.and.(bstmatrix(ref,1).le.bstmatrix(jlast,2)-dble(pad)))
      
!set the starting position of the Xprobe to be within the pad of the reference spike
         do while (abs(bstmatrix(ref,1)-bstmatrix(Xprobe,2)) .gt. pad)
            Xprobe=Xprobe+1
         end do

         do k=Xprobe,size(bstmatrix(:,1))
            if (abs(bstmatrix(ref,1)-bstmatrix(k,2)).le.pad) then
               d=dnint(bstmatrix(ref,1)-bstmatrix(k,2))+pad
               burstxcor(d)=burstxcor(d)+1
!               print *, 'IM RECORDING CORRELATIONS!!'
            else
               exit
            end if
         end do

         ref=ref+1
      end do
      if (bstmatrix(ilast,1).ge.bstmatrix(jlast,2)) then
         j=jlast
!         print *,'NEW j VALUE=',j,'osc2(j)=',bstmatrix(j,2)
      else
         i=ilast
!         print *,'NEW i VALUE=',i,'osc1(i)=',bstmatrix(i,1)
      end if
!      print *, 'END OF COMMON REGION-PAD=',bstmatrix(ref,1)
!      print *, 'osc1last=',bstmatrix(ilast,1),' osc2last=',bstmatrix(jlast,2)
      return
      end subroutine bstXcor

!*****************************Window Analysis****************************
      subroutine winanalysis(osc1vector,osc2vector,winrank,&
                             &i,ilast,j,jlast,bstmatrix,winsize,&
                             &maxtime,winsteplength,dt,simend)
      integer(kind=8) maxtime,probe1,probe2
      integer(kind=8) i,ilast,j,jlast,k
      integer(kind=8) winrank,wintotal,winstep
      double precision bstmatrix(maxtime,7)
      double precision osc1vector(maxtime),osc2vector(maxtime) 
      double precision Tstart,Tend
      double precision winend,winsize,winstart,winsteplength
      double precision dt,simend
      
!      print *,'window analysis'
!define the timewindow of analysis      
      Tstart=max(bstmatrix(i,1),bstmatrix(j,2))
      Tend=min(bstmatrix(ilast,1),bstmatrix(jlast,2),simend)

      wintotal = (Tend-Tstart)/(winsteplength*winsize)
!      print *,'Tstart=',Tstart,'Tend=',Tend
!      print *,'i=',i,'ilast=',ilast,'j=',j,'jlast=',jlast
!      print *,'winrank=',winrank,'winsteplength=',winsteplength
!      print *, 'winsize=',winsize,'wintotal=',wintotal


!count the number of spikes within each window range   
      probe1=i
      probe2=j
      do 113 winstep=1,wintotal   
         winstart    = (winsteplength*(winstep-1))*winsize + Tstart
         winend      = winstart + winsize
         
!        keep track of window times:
         bstmatrix(winrank+winstep,5)=winstart
         
!         print *,'winstart=',winstart,'winend=',winend,'probe1=',probe1,'probe2=',probe2
!         print *,'bstmatrix(probe1,1)=',bstmatrix(probe1,1),'bstmatrix(probe2,2)=',bstmatrix(probe2,2)
!set the probe each time to right inside the window
         do while (bstmatrix(probe1,1) .gt. winstart)
            if (probe1 .eq. 1) exit
            probe1=probe1-1
         end do
         !this protects from arriving at the very end of the spiketime list
         if (bstmatrix(probe1,1) .lt. bstmatrix(probe1-1,1)) then
            probe1=probe1-1
         end if
         !this puts the probe right into the window
         do while (bstmatrix(probe1,1) .lt. winstart)
            if (bstmatrix(probe1,1) .le. bstmatrix(probe1-1,1)) exit
            probe1=probe1+1
         end do
         !this will write the spikes into the window vector
         do while (bstmatrix(probe1,1) .le. winend)
            if (bstmatrix(probe1,1) .le. bstmatrix(probe1-1,1)) exit
            osc1vector(winstep+winrank)=osc1vector(winstep+winrank)+1
            bstmatrix(probe1,3)=bstmatrix(probe1,1)
            probe1=probe1+1
         end do
               
         do while (bstmatrix(probe2,2) .gt. winstart)
            if (probe2 .eq. 1) exit
            probe2=probe2-1
         end do
         if (bstmatrix(probe2,2) .lt. bstmatrix(probe2-1,2)) then
            probe2=probe2-1
         end if
         do while (bstmatrix(probe2,2) .lt. winstart)
            if (bstmatrix(probe2,2) .le. bstmatrix(probe2-1,2)) exit
            probe2=probe2+1
         end do
         do while (bstmatrix(probe2,2) .le. winend)
           if (bstmatrix(probe2,2) .le. bstmatrix(probe2-1,2)) exit
           osc2vector(winstep+winrank)=osc2vector(winstep+winrank)+1
           bstmatrix(probe2,4)=bstmatrix(probe2,2)
           probe2=probe2+1
         end do
         bstmatrix(winstep+winrank,6)=osc1vector(winstep+winrank)
         bstmatrix(winstep+winrank,7)=osc2vector(winstep+winrank)
  113 end do
!set the winrank to the new value
      winrank=winrank+wintotal
!      print *,winrank

      return
      end subroutine winanalysis
