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


integer :: npar, neps, npl, nreh

integer, parameter :: reclUnit = 4

integer, parameter :: lenmax = 80
character(len=lenmax), dimension(:), allocatable :: header

logical, save :: reset = .true.

public delete_file
public livewrite, allwrite, binallwrite

public cleanwrite_header, cleanwrite_data
public has_shifted, has_not_shifted 


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


  

  subroutine sp_livewrite(name,x,a,b,c,d,e,f,g,h,i,j,k,l,m,n)
    implicit none
    character(len=*) :: name    
    real :: x,a
    real, optional :: b,c,d,e,f,g,h,i,j,k,l,m,n
      
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
    elseif (.not.present(k)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i,j
    elseif (.not.present(l)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i,j,k
    elseif (.not.present(m)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i,j,k,l
    elseif (.not.present(n)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i,j,k,l,m
    else
       write(10,100) x,a,b,c,d,e,f,g,h,i,j,k,l,m,n
    endif
    
    close(10)

 100 format(15(ES25.16E3))

  end subroutine sp_livewrite


  subroutine kp_livewrite(name,x,a,b,c,d,e,f,g,h,i,j,k,l,m,n)
    implicit none
    character(len=*) :: name    
    real(kp) :: x,a
    real(kp), optional :: b,c,d,e,f,g,h,i,j,k,l,m,n
    
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
    elseif (.not.present(k)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i,j
    elseif (.not.present(l)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i,j,k
    elseif (.not.present(m)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i,j,k,l
    elseif (.not.present(n)) then
       write(10,100) x,a,b,c,d,e,f,g,h,i,j,k,l,m
    else
       write(10,100) x,a,b,c,d,e,f,g,h,i,j,k,l,m,n
    endif
    
    close(10)
    
100 format(15(ES25.16E3))
    
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




  function has_not_shifted(maxshift,x,y,z)
    implicit none
    logical :: has_not_shifted
    real(kp), intent(in) :: maxshift,x
    real(kp), intent(in), optional :: y,z
    
    has_not_shifted = .not.has_shifted(maxshift,x,y,z)

  end function has_not_shifted


  
  function has_shifted(maxshift,x,y,z)
    implicit none
    logical :: has_shifted
    real(kp), intent(in) :: maxshift,x
    real(kp), intent(in), optional :: y,z

    real(kp), save :: xold = tiny(1._kp)
    real(kp), save :: yold = tiny(1._kp)
    real(kp), save :: zold = tiny(1._kp)
       
    real(kp) :: d2

    has_shifted = .true.

    d2 = (x-xold)**2

    if (present(y)) then
       d2 = d2 + (y-yold)**2
    endif

    if (present(z)) then
       d2 = d2 + (z-zold)**2
    end if
    

    if (d2 == 0) then
       has_shifted = .false.
       return
    endif

    if (sqrt(d2).le.maxshift) then
       has_shifted = .false.
       return
    endif

    xold = x
    if (present(y)) yold = y
    if (present(z)) zold = z

  end function has_shifted


  
  subroutine int2char_adjtrim(icount,clen,ccount)
    implicit none
    integer, intent(in) :: icount
    integer, intent(in) :: clen
    character(len=clen), intent(out) :: ccount

    character(len=lenmax) :: numToStrg
    character(len=clen) :: strg

    if (clen.gt.lenmax) then
       stop 'int2char: clen > lenIoMax!'
    endif

    write(numToStrg,*) icount
    strg = trim(adjustl(numToStrg))    
    ccount(1:clen) = strg(1:clen)
  end subroutine int2char_adjtrim


  
  subroutine cleanwrite_header(asname,texpar,texeps,texpl,texreh)
    implicit none
    character(len=*), intent(in) :: asname
    character(len=lenmax), dimension(:), intent(in) :: texpar, texeps
    character(len=lenmax), dimension(:), intent(in) :: texpl, texreh

    integer :: i, ncols
    character :: c
    
    if (allocated(header)) stop 'cleanwrite_header: already created!'
    if (allocated(srfilenames).or.allocated(plfilenames)) then
       stop 'cleanwrite_header: filenames already set!'
    endif

    npar = size(texpar,1)
    neps = size(texeps,1)
    npl = size(texpl,1)
    nreh = size(texreh,1)

    ncols = npar + neps + npl + nreh
    allocate(header(ncols))

    do i=1,npar
       header(i) = '#'//trim(adjustl(texpar(i)))
    enddo
    do i=1,neps
       header(npar+i) = '#'//trim(adjustl(texeps(i)))
    enddo
    do i=1,npl
       header(npar+neps+i) = '#'//trim(adjustl(texpl(i)))
    enddo
    do i=1,nreh
       header(npar+neps+npl+i) = '#'//trim(adjustl(texreh(i)))
    enddo
           
       
  end subroutine cleanwrite_header



  subroutine cleanwrite_data(param, epsv, pl, reh)
    implicit none
!varying model parameters, from slowest to fastest
    real(kp), dimension(:), intent(in) :: param
!slow-roll parameters: epsV1, epsV2 etc...
    real(kp), dimension(:), intent(in) :: epsv
!powerlaw parameters: ns, r, alphaS etc...
    real(kp), dimension(:), intent(in) :: pl
!reheating parameters, Treh, Ereh, Delta Nstar...
    real(kp), dimension(:), intent(in) :: reh

    character(len=lenmax) :: srfilename, plfilename
    
    real(kp), save :: stored
    integer, save :: count

    integer, parameter :: clen = 2
    character(len=clen), save :: ccount

    integer, parameter :: slu = 300
    integer, parameter :: plu = 301
    
    logical :: sane

    
    sane = (size(param,1).eq.npar).and.(size(epsv,1).eq.neps) &
         .and.(size(pl,1).eq.npl).and.(size(reh,1).eq.nreh)

    if (.not.sane) stop 'cleanwrite_data: sizes input do no match header!'

    if (reset) then
       count = 0
       stored = 0._kp
       reset = .false.
    endif
    
    srfilename = trim(asname)//'_slowroll_'
    plfilename = trim(asname)//'_powerlaw_'
        
    if ((param(npar-1).ne.stored)) then

       count = count + 1
       call int2char_adjtrim(count,clen,ccount)

       open(sru,file=srfilename//ccount//'.dat',status='unknown')
       write(sru,*) (header(i),i=1,npar-1)
       write(sru,*) (param(i),i=1,npar-1)
       write(sru,*) header(npar),(header(npar+i),i=1,neps) &
            ,(header(npar+neps+npl+j),j=1,nreh)
       write(sru,*) param(npar),(epsv(i),i=1,neps),(reh(j),j=1,nreh)

       open(plu,file=plfilename//ccount//'.dat',status='unknown')
       write(plu,*) (header(i),i=1,npar-1)
       write(plu,*) (param(i),i=1,npar-1)
       write(plu,*) header(npar),(header(npar+i),i=1,neps) &
            ,(header(npar+neps+npl+j),j=1,nreh)
       write(plu,*) param(npar),(pl(i),i=1,neps),(reh(j),j=1,nreh)

       close(sru,plu)

    else

       open(sru,file=srfilename//ccount//'.dat',status='append')
       write(sru,*) param(npar),(epsv(i),i=1,neps),(reh(j),j=1,nreh)

       open(plu,file=plfilename//ccount//'.dat',status='append')
       write(plu,*) param(npar),(pl(i),i=1,neps),(reh(j),j=1,nreh)

       close(sru,plu)

    endif

    stored = param(npar-1)


  end subroutine cleanwrite_data


  subroutine cleanwrite_end()
    implicit none

    reset = .true.
    
  end subroutine cleanwrite_end



end module infinout


