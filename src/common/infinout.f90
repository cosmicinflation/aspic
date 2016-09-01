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



integer, parameter :: lenhead = 30
character(len=lenhead), dimension(:), allocatable :: header

integer, parameter :: srun = 300
integer, parameter :: plun = 301
integer, parameter :: lenmax = 80
character(len=lenmax) :: srfileprefix, plfileprefix

character(len=*), dimension(2), parameter :: labeps12 = &
     (/'$\epsilon_1$','$\epsilon_2$'/)
character(len=*), dimension(3) , parameter :: labeps123 = &
     (/'$\epsilon_1$','$\epsilon_2$','$\epsilon_3$'/)
character(len=*), dimension(2), parameter :: labnsr = &
     (/'$n_\mathrm{S}$','$r$           '/)
character(len=*), dimension(2), parameter :: labbfoldreh = &
     (/'$\Delta N_*$            ','$\ln(\rho_\mathrm{reh})$'/)


logical, save :: reset = .true.


public delete_file
public livewrite, allwrite, binallwrite
public has_shifted, has_not_shifted

public labeps12, labeps123, labnsr, labbfoldreh
public aspicwrite_header, aspicwrite_data, aspicwrite_end



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


  
  subroutine aspicwrite_header(asname,labeps,labpl,labreh,labpar)
    implicit none
    character(len=*), intent(in) :: asname
    character(len=*), dimension(:), intent(in) :: labeps, labpl, labreh
    character(len=*), dimension(:), intent(in), optional :: labpar
    
    integer :: i, j, ncols
    
    if (allocated(header)) stop 'cleanwrite_header: already created!'

    srfileprefix = trim(asname)//'_slowroll_'
    plfileprefix = trim(asname)//'_powerlaw_'
       
    neps = size(labeps,1)
    npl = size(labpl,1)
    nreh = size(labreh,1)
    
    if (present(labpar)) then
       npar = size(labpar,1)
    else
       npar = 0
    endif
    
    ncols = neps + npl + nreh + npar
    allocate(header(ncols))
    
    do i=1,neps
       header(i) = trim(adjustl(labeps(i)))
    enddo
    do i=1,npl
       header(neps+i) = trim(adjustl(labpl(i)))
    enddo
    do i=1,nreh
       header(neps+npl+i) = trim(adjustl(labreh(i)))
    enddo

    if (present(labpar)) then
       header(neps+npl+nreh+1) = trim(adjustl(labpar(1)))
       do i=2,npar
          header(neps+npl+nreh+i) = '# '//trim(adjustl(labpar(i)))
       enddo
    endif

    header(1)='# '//trim(header(1))
    header(neps+1)='# '//trim(header(neps+1))
    
    
    reset = .true.
       
  end subroutine aspicwrite_header



  subroutine aspicwrite_data(epsv, pl, reh, param)
    implicit none
!slow-roll parameters: epsV1, epsV2 etc...
    real(kp), dimension(:), intent(in) :: epsv
!powerlaw parameters: ns, r, alphaS etc...
    real(kp), dimension(:), intent(in) :: pl
!reheating parameters, Treh, Ereh, Delta Nstar...
    real(kp), dimension(:), intent(in) :: reh
!varying model parameters, from fastest to slowest: only param(1) is
!varying in the output file. There are as many outout files as
!different values of params(2) x params(3) x ...
    real(kp), dimension(:), intent(in), optional :: param
    
    real(kp), save :: stored
    integer, save :: count

    integer, parameter :: clen = 2
    character(len=clen), save :: ccount

    integer :: i,j,ipar
    logical :: sane, isopened
    character(len=lenmax) :: filename
    real(kp) :: current
    
    sane = (size(epsv,1).eq.neps) &
         .and.(size(pl,1).eq.npl).and.(size(reh,1).eq.nreh)

    if (present(param)) then
       sane = sane.and.(size(param,1).eq.npar)
    else
       sane = sane.and.(npar.eq.0)
    endif
    
    if (.not.sane) stop 'cleanwrite_data: sizes input do no match header!'

    if (reset) then
       count = 0
       stored = 0._kp
       current = -1._kp
       reset = .false.
    endif
    
    if (npar.gt.1) current = param(2)

    ipar = neps+npl+nreh
    
    if (current.ne.stored) then

       count = count + 1
       call int2char_adjtrim(count,clen,ccount)

       inquire(unit=srun,opened=isopened)
       if (isopened) close(srun)
       filename = trim(srfileprefix)//trim(ccount)//'.dat'
       open(srun,file=trim(filename),status='replace',position='append')

       inquire(unit=plun,opened=isopened)
       if (isopened) close(plun)
       filename = trim(plfileprefix)//trim(ccount)//'.dat'
       open(plun,file=trim(filename),status='replace',position='append')

       if (npar.gt.1) then
          do i=2,npar           
             write(srun,*) trim(header(i+ipar)),param(i)
             write(plun,*) trim(header(i+ipar)),param(i)
          enddo          
       endif
       
       if (npar.ge.1) then
          write(srun,*) (header(i),i=1,neps),(header(neps+npl+j),j=1,nreh),header(ipar+1)
          write(srun,*) (epsv(i),i=1,neps),(reh(j),j=1,nreh),param(1)

          write(plun,*) (header(neps+i),i=1,npl),(header(neps+npl+j),j=1,nreh),header(ipar+1)
          write(plun,*) (pl(i),i=1,neps),(reh(j),j=1,nreh),param(1)

       else
          write(srun,*) (header(i),i=1,neps),(header(neps+npl+j),j=1,nreh)
          write(srun,*) (epsv(i),i=1,neps),(reh(j),j=1,nreh)
          
          write(plun,*) (header(neps+i),i=1,npl),(header(neps+npl+j),j=1,nreh)
          write(plun,*) (pl(i),i=1,neps),(reh(j),j=1,nreh)
       endif
                
    else

       if (npar.ge.1) then
          write(srun,*) (epsv(i),i=1,neps),(reh(j),j=1,nreh),param(1)
          write(plun,*) (pl(i),i=1,neps),(reh(j),j=1,nreh),param(1)
       else
          write(srun,*) (epsv(i),i=1,neps),(reh(j),j=1,nreh)
          write(plun,*) (pl(i),i=1,neps),(reh(j),j=1,nreh)
       endif

    endif

    stored = current


  end subroutine aspicwrite_data


  subroutine aspicwrite_end()
    implicit none

    reset = .true.
    close(srun)
    close(plun)
    
  end subroutine aspicwrite_end



end module infinout


