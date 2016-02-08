module displine
  use infprec, only : kp, pi
  use bspline, only : dbsnak, dbsint, dbsval
  
  implicit none


  private

  interface di_spline_x
     module procedure di_x
  end interface di_spline_x

  interface di_spline_k2
     module procedure di_k2
  end interface di_spline_k2

!precomputed modulus field correspondance
  character(len=*), parameter :: tablefile = 'fieldvalues.sav'

  integer, parameter :: order = 3
  integer, save :: bcoefNum = 0
  integer, save :: ndata = 0

  real(kp), save :: xmin, xmax, k2min, k2max
  real(kp), dimension(:), allocatable :: xknots, xbcoef
  real(kp), dimension(:), allocatable :: k2knots, k2bcoef

!max number of records in precomputed data file
  integer, parameter :: nrecmax = 1000
  integer, parameter :: ncols = 2

!set it to false to use binary preprocessed data
  logical, parameter :: usePP = .true.

!set it to true to generate binary preprocessed data in fieldvalues.pp
!from data file
  logical, parameter :: createPP = .false.


  public k2min, k2max, xmin, xmax
  public di_spline_x, di_spline_k2
  public di_check_splines, di_free_splines,di_set_splines


contains
  
  
  subroutine di_preprocessed_data(xdata,k2data)
    implicit none
    real(kp), dimension(nrecmax,ncols) :: buffer
    real(kp), dimension(:), allocatable, intent(inout) :: xdata, k2data
    integer :: i,nrec

    buffer = -1._kp

#define PPDATA(fooidx,fook2,foox) \
    buffer(fooidx,1) = fook2 ; \
    buffer(fooidx,2) = foox
#include "fieldvalues.pp"
#undef PPDATA

    do i=1,nrecmax
       if (buffer(i,1).eq.-1._kp) exit
    enddo

    write(*,*)
    write(*,*)'di_preprocessed_data:'
    write(*,*)'number of records: ',i-1

    nrec = i-1

    allocate(xdata(nrec),k2data(nrec))
    k2data = buffer(:,1)
    xdata = buffer(:,2)

    ndata = nrec
    k2min = k2data(1)
    k2max = k2data(ndata)
    xmin = xdata(ndata)
    xmax = xdata(1)
    
  end subroutine di_preprocessed_data



  subroutine di_readump_data(filename,xdata,k2data)
    implicit none

    character(len=*), intent(in) :: filename
    real(kp), dimension(:), allocatable, intent(inout) :: xdata, k2data

    integer, parameter :: nunit = 110
    logical, parameter :: dumpPPdata = .true.

    integer :: ioerr, nrec,i

    real(kp), dimension(ncols), save :: statbuffer
    real(kp), dimension(nrecmax,ncols) :: buffer
    
    write(*,*)'reading field values table...'


    open(unit=nunit,file=filename,status='old')   

    do i=1,nrecmax
       read(nunit,*,iostat=ioerr) statbuffer         

       buffer(i,1) = statbuffer(1)
       buffer(i,2) = statbuffer(2)

       if (ioerr.ne.0) exit

    enddo

    close(nunit)

    write(*,*)
    write(*,*)'di_read_table:'
    write(*,*)'number of records: ',i-1

    nrec = i-1

    allocate(xdata(nrec),k2data(nrec))
    k2data = buffer(:,1)
    xdata = buffer(:,2)

    ndata = nrec
    k2min = k2data(1)
    k2max = k2data(ndata)
    xmin = xdata(ndata)
    xmax = xdata(1)

    if (createPP) then
       open(unit=nunit,file='fieldvalues.pp',status='new')
       do i=1,nrec
          write(nunit,*)'PPDATA(',i,',',buffer(i,1),',',buffer(i,2),')'
       enddo
       close(nunit)
    endif

  end subroutine di_readump_data



  function di_check_splines()
    implicit none
    logical :: di_check_splines

    di_check_splines = allocated(k2knots).or.allocated(xknots)
  
!sanity checks
    if (di_check_splines) then      
       if (.not.allocated(k2bcoef)) stop 'di_check_splines: k2bcoef not found!'
       if (.not.allocated(xbcoef)) stop 'di_check_splines: xbcoef not found!'
    endif

  end function  di_check_splines
  


  subroutine di_free_splines()
    implicit none

    if (di_check_splines()) then
       deallocate(k2knots,xknots)
       deallocate(k2bcoef,xbcoef)
    endif

  end subroutine di_free_splines



  subroutine di_set_splines
    implicit none

    real(kp), dimension(:), allocatable :: xdata, k2data

    if (di_check_splines()) stop 'di_set_splines: splines already set!'

    if (usePP) then
       call di_preprocessed_data(xdata, k2data)
    else
       call di_readump_data(tablefile, xdata, k2data)
    endif

    allocate(k2knots(ndata+order),xknots(ndata+order))
    allocate(k2bcoef(ndata),xbcoef(ndata))

    call dbsnak(ndata,k2data,order,k2knots)
    call dbsint(ndata,k2data,xdata,order,k2knots,k2bcoef)

    call dbsnak(ndata,xdata(ndata:1:-1),order,xknots)
    call dbsint(ndata,xdata(ndata:1:-1),k2data(ndata:1:-1),order,xknots,xbcoef)
    
    deallocate(xdata,k2data)

  end subroutine di_set_splines


  function di_x(k2)
    implicit none
    real(kp) :: di_x
    real(kp), intent(in) :: k2        

    if (k2.eq.1._kp) then
       di_x = 0._kp
       return
    endif

    if ((k2.le.k2min) .or. (k2.gt.k2max)) then      
       write(*,*)'k2= k2min= k2max= ',k2, k2min, k2max
       stop 'di_k2: x out of range!'
    endif

    di_x = dbsval(k2,order,k2knots,ndata,k2bcoef)

  end function di_x



  function di_k2(x)
    implicit none
    real(kp) :: di_k2
    real(kp), intent(in) :: x

    if (x.eq.0._kp) then
       di_k2 = 1._kp
       return
    endif

    if ((x.le.xmin).or.(x.gt.xmax)) then
       write(*,*)'x= xmin= xmax= ',x, xmin, xmax
       stop 'di_k2: x out of range!'
    endif
    
    di_k2 = dbsval(x,order,xknots,ndata,xbcoef)

  end function di_k2
  

end module displine
