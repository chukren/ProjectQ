c
c   read eigf.out(love)  r_eigf.out(rayleigh) to plot mode depth profiles of
c       W , dcdb dcdrho 
c       U , V, dcda , dcdb dcdrho
c 
c   Tony  s = U*Plm + U*Blm + W*Clm
c   Guust s = U*Plm + l*(l+1)*U*Blm  +  l*(l+1)*W*Clm
c 
c units:   
c         displacement in meters
c         R in meters
c         rho in kg/m^3
c
c
       parameter(pi2 = 3.1415926*2,nlmax=1000,modmax=1000)
       integer*4 ls0,m1,kind
       real*4 w4,vphase,gcom,qmod,ekin,wdiff,fl4,sqll
       real*4 y1,dy1dr,y2,dy2dr,dcda,dcdb,dcdrh,dcdxi
       dimension nlay(modmax),nstart(modmax),nend(modmax),radius(nlmax),beta(nlmax)
       character*70 filelay,fileeigf,filedisp,filetest,file_model


       write(*,*)'model file name:'
       read(*,'(a)')file_model
       write(*,*)'1:love wave  2:rayleigh wave '
       read(*,*)is_love
       write(*,*)'order n  from 0 :'
       read(*,*)no
       write(*,*)'freqency in Hz (0.01 for 100s):'
       read(*,*)freq

       ! write(*,*)'reading model file:', file_model ! no output chukren

       if(is_love.eq.1) then
        if(no.lt.10) write(filelay,  '(''layer.love.'',i1)' )no  
        if(no.gt.9) write(filelay,  '(''layer.love.'',i2)' )no  
        if(no.lt.10) write(fileeigf,  '(''eigf.love.'',i1)' )no  
        if(no.gt.9) write(fileeigf,  '(''eigf.love.'',i2)' )no  
        if(no.lt.10) write(filedisp,  '(''disp.love.'',i1)' )no  
        if(no.gt.9) write(filedisp,  '(''disp.love.'',i2)' )no  
       else
        if(no.lt.10) write(filelay,  '(''layer.rayleigh.'',i1)' )no  
        if(no.gt.9) write(filelay,  '(''layer.rayleigh.'',i2)' )no  
        if(no.lt.10) write(fileeigf,  '(''eigf.rayleigh.'',i1)' )no  
        if(no.gt.9) write(fileeigf,  '(''eigf.rayleigh.'',i2)' )no  
        if(no.lt.10) write(filedisp,  '(''disp.rayleigh.'',i1)' )no  
        if(no.gt.9) write(filedisp,  '(''disp.rayleigh.'',i2)' )no  
       endif 

       open(10,file=filedisp)
       read(10,*)itotal
       !write(*,*)itotal ! no output chukren
       close(10)

c read recorded layer number for diff  freq 

         w_need = freq*pi2

         w4_old = 0

       open(90,file=filelay)
         read(90,*)
         do i =1,itotal
           read(90,91)w4,vphase,gcom,nlay(i),nstart(i),nend(i)
           if(w_need.ge.w4_old.and.w_need.lt.w4) mod = i
           w4_old = w4
         enddo 
         close(90) 

91      format(1x,e12.4,1x,f12.5,1x,f12.5,1x,i5,i5,i5)     

       
       open(40,file=file_model)
       read(40,*)
       read(40,*)
       read(40,*)number
       do i = 1, number
       read(40,*)radius(i),on,on,beta(i)
       enddo
       close(40)
       !write(*,*)'layer=',number ! no output chukren

       !write(*,*)'1: reading eigenfunction',fileeigf ! no output
       !chukren

       open(10,file=fileeigf,form='unformatted')

       !write(*,*)'2: reading eigenfunction',fileeigf ! no output
       !chukren

       if(is_love.eq.1) then        

       if(no.lt.10) write(filetest,  '(''function.W.'',i1)' )no  
       if(no.gt.9) write(filetest,  '(''function.W.'',i2)' )no  
       open(20,file=filetest)
       !write(*,*)' love wave eigen function', filetest ! no output

       if(no.lt.10) write(filetest,  '(''function.dW.'',i1)' )no  
       if(no.gt.9) write(filetest,  '(''function.dW.'',i2)' )no  
       open(21,file=filetest)
       !write(*,*)' love wave eigen function', filetest ! no output

       else 

       if(no.lt.10) write(filetest,  '(''function.U.'',i1)' )no  
       if(no.gt.9) write(filetest,  '(''function.U.'',i2)' )no  
       open(20,file=filetest)
       !write(*,*)' rayleigh wave eigen function : ', filetest ! no
       !output chukren

       if(no.lt.10) write(filetest,  '(''function.dU.'',i1)' )no  
       if(no.gt.9) write(filetest,  '(''function.dU.'',i2)' )no  
       open(21,file=filetest)
       !write(*,*)' rayleigh wave eigen function : ', filetest ! no
       !output chukren

       if(no.lt.10) write(filetest,  '(''function.V.'',i1)' )no  
       if(no.gt.9) write(filetest,  '(''function.V.'',i2)' )no  
       open(30,file=filetest)
       !write(*,*)' rayleigh wave eigen function : ', filetest ! no
       !output chukren

       if(no.lt.10) write(filetest,  '(''function.dV.'',i1)' )no  
       if(no.gt.9) write(filetest,  '(''function.dV.'',i2)' )no  
       open(31,file=filetest)

       !write(*,*)' rayleigh wave eigen function : ', filetest ! no
       !output chukren
       endif


       read(10)rn4,in4,n4,nfreq,nbran,df4,cmax4
       ! write(*,*)rn4,in4,n4,nfreq,nbran,df4,cmax4 ! no output chukren

       ! write(*,*)'n4=',n4

       do i = 1,n4
        read(10)rn4,vp4,vs4,rh4,q4
c        write(*,*)rn4,vp4,vs4,rh4,q4
       enddo 


cc love wave        
       if(is_love.eq.1)then

        do k = 1,mod

         read(10) w4,vphase,gcom,qmod,ekin,wdiff,fl4,sqll,ls0,m1,kind
         !write(*,*)'T=',(pi2/w4),'vphase=',vphase,'gcom=',gcom,'fl4=',fl4 ! no output chukren


          fact = sqrt((fl4+1)*fl4/(vphase*gcom))
          factor = w4 * 6371000 *fact

           if(k.eq.mod) write(20,*)'function    depth'
           if(k.eq.mod) write(21,*)'function    depth'

c           write(*,*)'k= nlay(k)=',k,nlay(k)

           add =0
            do i = 1,nlay(k)
              read(10) y1,dy1dr,dcda,dcdb,dcdrh,dcdxi
              if(k.eq.mod) then
                ir = nstart(k)+i-1 

cc r^2 included in dcdb
              add = add + dcdb*(radius(ir)-radius(ir-1))*beta(ir)/vphase

               y1 = y1 * factor*1e9 
               dy1dr = dy1dr * factor*1e9 

               write(20,5)y1,(6371000-radius(ir))/1000
               write(21,5)dy1dr,(6371000-radius(ir))/1000

              endif
           enddo

        enddo   

cc rayleigh        
       else 

        do k = 1,mod
        read(10) w4,vphase,gcom,qmod,ekin,wdiff,fl4,sqll,ls0,m1,kind
        write(*,*) 'T=',(pi2/w4),'vphase=',vphase,'gcom=',gcom,'fl4=',fl4
cc change norm        
           fact = sqrt((fl4+1)*fl4/(vphase*gcom))
           fact2 = sqrt(1.0/(vphase*gcom))
           factor = w4 * 6371000 *fact
           factor2 = w4 * 6371000 *fact2
           if(k.eq.mod) write(20,*)'function    depth'
           if(k.eq.mod) write(30,*)'function    depth'
           if(k.eq.mod) write(21,*)'function    depth'
           if(k.eq.mod) write(31,*)'function    depth'

           add =0
           do i =1, nlay(k)
           read(10) y1,dy1dr,y2,dy2dr,dcda,dcdb,dcdrh,dcdxi
           if(k.eq.mod) then
                ir = nstart(k)+i-1 

cc r^2 included in dcdb
              add = add + dcdb*(radius(ir)-radius(ir-1))*beta(ir)/vphase
c              write(*,*)ir,radius(ir),beta(ir),vphase

cc change norm and direction for V              
cc no direction change

              y1 = y1 * factor2*1e9 
              y2 = y2 * factor  *1e9 
              dy1dr = dy1dr * factor2*1e9 
              dy2dr = dy2dr * factor  *1e9 

           write(20,5)y1,(6371000-radius(ir))/1000
           write(30,5)y2,(6371000-radius(ir))/1000
           write(21,5)dy1dr,(6371000-radius(ir))/1000
           write(31,5)dy2dr,(6371000-radius(ir))/1000
          endif       
          enddo
        enddo

       endif 

       !    write(*,*)'dc/db integreal with depth ',add   ! no ouput
       !    chukren 

       close(10)
       close(20)
       close(21)
       if(is_love .eq. 2) then
       close(30)
       close(31)
       endif
3      format(6(e12.4,1x))
4      format(8(e12.4,1x))
5        format(2(e12.4,1x))
       stop
       end
