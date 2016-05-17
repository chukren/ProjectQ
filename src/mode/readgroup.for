	  parameter(modemax=50,layermax=650)

	  integer*4 m 
	  character*2 ichar,charL,charR
	  real*4 vphase,gcom,qmod,wdiff,flroot 
	   data  charL/' L'/,charR/' R'/

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

cc  love.out 

	if(is_love.eq.1) then
	  open(10,file='love.out')
	  open(20,file='love.grp')
	  write(*,*)'output love.grp'
	else  
	  open(10,file='rayleigh.out')
	  open(20,file='rayleigh.grp')
	  write(*,*)'output rayleigh.grp'
	endif  


cc read headers 
	  do i = 1,iheadout  
	  read(10,*)
	  enddo

	  do 2000 io = 1, iomax 
	  read(10,9)m,ichar,flroot,fhz,period,vphase,gcom,qmod,wdiff
	  if(m.eq.0) write(20,*)period,(gcom)
2000	  continue
	  
9         format(i5,a2,f10.2,6g16.7)	  

	  close(10)

          close(20)


	  stop
	  end


