cc seperate
cc          love.out : disp.love.0  disp.love.1
cc          layer.love: layer.love.0  layer.love.2
cc          eigf.love:  eigf.love.0 eigf.love.1


cc     new  layer.love.?? disp.love.??  eigf.love.??

cc     open love.out, layer.love, eigf.love
cc               read love.out --- one line 
cc               read layer.love --- one line
cc               1) which mode  ??  n   how many lines ?? ist-iend
cc               2) write this line to layer.love.??  disp.love.?? 
cc               3) read eigf.love (ist-iend) lines , write to  eigf.love.??  

cc           
cc note : write total number at layer.love.??


         parameter(modemax=300,layermax=650)

         character*70 filelay(modemax), filedisp(modemax), fileeigf(modemax)
         integer id(modemax), il(modemax), ie(modemax)
         integer*4 m 
         real*4 dcda,dcdb,dcdrh
         character*2 ichar,charL,charR
         real*4  WW,DW,UU,DU,VV,DV,dcdxi 
         real*4 w4,vphase,gcom,qmod,ekin,wdiff,fl4,sqll 
         integer*4 ls0,m1,kind
          data  charL/' L'/,charR/' R'/

c         immax = 31 
c         iheadout = 653
c         iomax = 1918

cc         immax = upper-mod + 1  n=0  means the id is id(1) il(1) ie(1)          

        write(*,*)'love wave (1)  rayleigh (2)'
        read(*,*)is_love

        if(is_love.eq.1) then
         write(*,*)'highest mode in love.in :'
         read(*,*)mod_upp
         immax = mod_upp+1
          write(*,*)'the headers lines before love.out:'
          read(*,*)iheadout
          write(*,*)'the lines of modes love.out:'
          read(*,*)iomax
        else
         write(*,*)'highest mode in rayleigh.in :'
         read(*,*)mod_upp
         immax = mod_upp+1
          write(*,*)'the headers lines before rayleigh.out:'
          read(*,*)iheadout
          write(*,*)'the lines of modes rayleigh.out:'
          read(*,*)iomax
        endif

cc open new files
         do im = 1,immax 
          if = im - 1 
          if(is_love.eq.1) then
                 if(if.lt.10)then
                 write(filedisp(im), '(''disp.love.'',i1)' )if
                 write(filelay(im),  '(''layer.love.'',i1)' )if 
                 write(fileeigf(im), '(''eigf.love.'',i1)' )if 
                 endif
                 if(if.gt.9)then
                 write(filedisp(im), '(''disp.love.'',i2)' )if
                 write(filelay(im),  '(''layer.love.'',i2)' )if 
                 write(fileeigf(im), '(''eigf.love.'',i2)' )if 
                 endif
                 if(if.gt.99)then
                 write(filedisp(im), '(''disp.love.'',i3)' )if
                 write(filelay(im),  '(''layer.love.'',i3)' )if 
                 write(fileeigf(im), '(''eigf.love.'',i3)' )if 
                 endif
          else
                 if(if.lt.10)then
                 write(filedisp(im), '(''disp.rayleigh.'',i1)' )if
                 write(filelay(im),  '(''layer.rayleigh.'',i1)' )if 
                 write(fileeigf(im), '(''eigf.rayleigh.'',i1)' )if 
                 endif
                 if(if.gt.9)then
                 write(filedisp(im), '(''disp.rayleigh.'',i2)' )if
                 write(filelay(im),  '(''layer.rayleigh.'',i2)' )if 
                 write(fileeigf(im), '(''eigf.rayleigh.'',i2)' )if 
                 endif
                 if(if.gt.99)then
                 write(filedisp(im), '(''disp.rayleigh.'',i3)' )if
                 write(filelay(im),  '(''layer.rayleigh.'',i3)' )if 
                 write(fileeigf(im), '(''eigf.rayleigh.'',i3)' )if 
                 endif
          endif       


           id(im) = 100 + im
           il(im) = 500 + im
           ie(im) = 1000 + im


          
          open(id(im),file=filedisp(im))
          open(il(im),file=filelay(im))
          open(ie(im),file=fileeigf(im),form='unformatted')
          
          write(il(im),1)'skipped  w  vphase gcom number-layer stratl endl'  
1         format(a)          
         enddo  


cc open input files
cc  love.out 

       if(is_love.eq.1) then
         open(10,file='love.out')
         open(20,file='layer.love')
         open(30,file='eigf.love',form='unformatted')
       else  
         open(10,file='rayleigh.out')
         open(20,file='layer.rayleigh')
         open(30,file='eigf.rayleigh',form='unformatted')
       endif  


cc read headers 
         do i = 1,iheadout  
         read(10,*)
         enddo
         read(20,*)
         read(30)rn4,in4,n4,nfreq,nbran,df4,cmax4
         do j=1,immax
         write(ie(j))rn4,in4,n4,nfreq,nbran,df4,cmax4
         enddo
         do i = 1,n4
         read(30)rn4,vp4,vs4,rh4,q4 
           do j=1,immax
           write(ie(j))rn4,vp4,vs4,rh4,q4
           enddo
         enddo

cc read line by line 

         do 2000 io = 1, iomax 
         read(10,9)m,ichar,flroot,fhz,period,vphase,gcom,qmod,wdiff
9         format(i5,a2,f10.2,6g16.7)         

         write(*,*)'io= order= ichar',io,m,ichar,flroot

         read(20,8)w4,vphase,gcom,number,ls00,n
8         format(1x,e12.4,1x,f12.5,1x,f12.5,1x,i5,i5,i5)         

         write(id(m+1),*)period,vphase
         write(il(m+1),8)w4,vphase,gcom,number,ls00,n

c         write(*,*)'number = ',number

cc  the mode eigf id j

         j = m+1

         read(30) w4,vphase,gcom,qmod,ekin,wdiff,fl4,sqll,ls0,m1,kind 
         write(ie(j)) w4,vphase,gcom,qmod,ekin,wdiff,fl4,sqll,ls0,m1,kind

c         write(*,*)w4,vphase,gcom,qmod,ekin,wdiff,fl4,sqll,ls0,m1,kind


         do 1000 ilay = 1,number
           if(ichar.eq.charL) then 
                read(30) WW,DW,dcda,dcdb,dcdrh,dcdxi
                write(ie(j)) WW,DW,dcda,dcdb,dcdrh,dcdxi
           endif           
           if(ichar.eq.charR) then 
                read(30) UU,DU,VV,DV,dcda,dcdb,dcdrh,dcdxi
                write(ie(j)) UU,DU,VV,DV,dcda,dcdb,dcdrh,dcdxi
           endif           

1000      continue

2000      continue


cc   close all files


         do im = 1,immax 
          write(*,*)'closed :',id(im),il(im),ie(im)
          close(id(im))
          close(il(im))
          close(ie(im))
          enddo



         close(10)
         close(20)
         close(30)




         stop
         end


