module infinout

use infprec, only : kp

private

interface livewrite
   module procedure sp_livewrite, kp_livewrite
end interface

interface allwrite
   module procedure sp_allwrite, kp_allwrite
end interface

interface binallwrite
   module procedure sp_binallwrite, kp_binallwrite
end interface

integer, parameter :: reclUnit = 4

public delete_file
public livewrite, allwrite, binallwrite
 
contains


  subroutine delete_file(name)
    implicit none
    character(len=*) :: name
    logical :: isthere

    inquire(file=name,exist=isthere)

    if (isthere) then
       open(unit=10,file=name)
       close(unit=10,status='delete')
    endif

  end subroutine delete_file


  subroutine sp_livewrite(name,x,a,b,c,d,e,f,g,h,i,j)
    implicit none
    character(len=*) :: name    
    real :: x,a
    real, optional :: b,c,d,e,f,g,h,i,j
      
    open(10,file=name,position='append',status='unknown')
    
    if (.not.present(b)) then
       write(10,100) x,a         
    elseif (.not.present(c)) then           
       write(10,100) x,a,b         
    elseif (.not.present(d)) then             
       write(10,100) x,a,b,c                    
    elseif (.not.present(e)) then         
       write(10,100) x,a,b,c,d                    
    elseif (.not.present(f)) then
       write(10,100) x,a,b,c,d,e            
    elseif (.not.present(g)) then
       write(10,100) x,a,b,c,d,e,f            
    elseif (.not.present(h)) then
       write(10,100) x,a,b,c,d,e,f,g  
    elseif (.not.present(i)) then
       write(10,100) x,a,b,c,d,e,f,g,h       
    elseif (.not.present(j)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i
    else
       write(10,100) x,a,b,c,d,e,f,g,h,i,j 
    endif
    
    close(10)

100 format(8(ES25.16E3))

  end subroutine sp_livewrite


  subroutine kp_livewrite(name,x,a,b,c,d,e,f,g,h,i,j)
    implicit none
    character(len=*) :: name    
    real(kp) :: x,a
    real(kp), optional :: b,c,d,e,f,g,h,i,j
    
    open(10,file=name,position='append',status='unknown')

    if (.not.present(b)) then
       write(10,100) x,a         
    elseif (.not.present(c)) then           
       write(10,100) x,a,b         
    elseif (.not.present(d)) then             
       write(10,100) x,a,b,c                    
    elseif (.not.present(e)) then         
       write(10,100) x,a,b,c,d                    
    elseif (.not.present(f)) then
       write(10,100) x,a,b,c,d,e            
    elseif (.not.present(g)) then
       write(10,100) x,a,b,c,d,e,f            
    elseif (.not.present(h)) then
       write(10,100) x,a,b,c,d,e,f,g  
    elseif (.not.present(i)) then
       write(10,100) x,a,b,c,d,e,f,g,h       
    elseif (.not.present(j)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i
    else
       write(10,100) x,a,b,c,d,e,f,g,h,i,j 
    endif
    
    close(10)
    
100 format(8(ES25.16E3))      
    
  end subroutine kp_livewrite


  subroutine sp_allwrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: j,npts
    real :: x(:),a(:)
    real, optional :: b(:),c(:),d(:),e(:),f(:),g(:)
    
    npts=ubound(x,1)
      
    if (ubound(a,1).ne.npts) then
       write(*,*)'WARNING: vectors length differ'
    endif
    
    write(*,*)'__write: save in ',name
    open(10,file=name,status='unknown')
    
    if (.not.present(b)) then
       do j=1,npts      
          write(10,100) x(j),a(j)
       enddo
    elseif (.not.present(c)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j)
       enddo
    elseif (.not.present(d)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j)            
       enddo
    elseif (.not.present(e)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j),d(j)            
       enddo
    elseif (.not.present(f)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j),d(j),e(j)            
       enddo
    elseif (.not.present(g)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
       enddo
    else
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
       enddo
    endif
    
    close(10)
    
100 format(8(ES25.16E3))      

  end subroutine sp_allwrite
  


  subroutine kp_allwrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: j,npts
    real(kp) :: x(:),a(:)
    real(kp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)

    npts=ubound(x,1)
      
    if (ubound(a,1).ne.npts) then
       write(*,*)'WARNING: vectors length differ'
    endif
    
    write(*,*)'__write: save in ',name
    open(10,file=name,status='unknown')
    
    if (.not.present(b)) then
       do j=1,npts      
          write(10,100) x(j),a(j)
       enddo
    elseif (.not.present(c)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j)
       enddo
    elseif (.not.present(d)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j)            
       enddo
    elseif (.not.present(e)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j),d(j)            
       enddo
    elseif (.not.present(f)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j),d(j),e(j)            
       enddo
    elseif (.not.present(g)) then
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
       enddo
    else
       do j=1,npts      
          write(10,100) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
       enddo
    endif
    
    close(10)
    
100 format(8(ES25.16E3))      
    
  end subroutine kp_allwrite



  subroutine sp_binallwrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: j,npts
    real :: x(:),a(:)
    real, optional :: b(:),c(:),d(:),e(:),f(:),g(:)
    
    integer :: datarecl
    integer :: recnum


    recnum = 0
    npts=ubound(x,1)
      
    if (ubound(a,1).ne.npts) then
       write(*,*)'WARNING: vectors length differ'
    endif
    
    write(*,*)'__write: save in ',name

    
    if (.not.present(b)) then
       datarecl=2*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j)
       enddo
    elseif (.not.present(c)) then
       datarecl=3*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j)
       enddo
    elseif (.not.present(d)) then
       datarecl=4*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j)            
       enddo
    elseif (.not.present(e)) then
       datarecl=5*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j)            
       enddo
    elseif (.not.present(f)) then
       datarecl=6*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j)            
       enddo
    elseif (.not.present(g)) then
       datarecl=7*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
       enddo
    else
       datarecl=8*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=datarecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
       enddo
    endif
    
    close(10)
    
  end subroutine sp_binallwrite
  


  subroutine kp_binallwrite(name,x,a,b,c,d,e,f,g)
    implicit none
    character(*) :: name
    integer :: j,npts
    real(kp) :: x(:),a(:)
    real(kp), optional :: b(:),c(:),d(:),e(:),f(:),g(:)

    integer :: recnum
    integer :: datarecl

    npts=ubound(x,1)
    recnum=0

    if (ubound(a,1).ne.npts) then
       write(*,*)'WARNING: vectors length differ'
    endif
    
    write(*,*)'__write: save in ',name
    
    if (.not.present(b)) then
       datarecl = 4*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum=recnum+1
          write(10,rec=recnum) x(j),a(j)
       enddo
    elseif (.not.present(c)) then
       datarecl = 6*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum =recnum+1
          write(10,rec=recnum) x(j),a(j),b(j)
       enddo
    elseif (.not.present(d)) then
       datarecl = 8*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum =recnum+1
          write(10,rec=recnum) x(j),a(j),b(j),c(j)            
       enddo
    elseif (.not.present(e)) then
       datarecl = 10*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum = recnum + 1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j)            
       enddo
    elseif (.not.present(f)) then
       datarecl = 12*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum = recnum + 1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j)            
       enddo
    elseif (.not.present(g)) then
       datarecl = 14*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum = recnum + 1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j),f(j)       
       enddo
    else
       datarecl = 16*reclUnit
       open(10,file=name,status='unknown',form='unformatted',access='direct',recl=dataRecl)
       do j=1,npts      
          recnum = recnum + 1
          write(10,rec=recnum) x(j),a(j),b(j),c(j),d(j),e(j),f(j),g(j)            
       enddo
    endif
    
    close(10)
        
  end subroutine kp_binallwrite


end module infinout


