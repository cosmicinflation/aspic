module inftools
  use infprec, only : kp,transfert  
  implicit none

  private

  public easydverk, tunedverk
  public zbrent, newton
  public quarticroots, sort, selectsort



contains
 
 

  subroutine easydverk(n,fcn,x,y,xend,tol,extradata)
    implicit none
     
    integer :: n,ind
    real(kp) :: x,y(n),xend,tol,c(24),w(n,9)    
    type(transfert), optional,intent(inout) :: extradata    

!    external :: fcn

    include 'inftools.h'

   
!    ind=1    
    ind=2
    c = 0._kp
    c(3) = epsilon(1._kp)    
 
    call dverk(n,fcn,x,y,xend,tol,ind,c,n,w,extradata)  
    
    if (ind.ne.3) then
       write(*,*) 'easydverk: stop ind = ',ind
       write(*,*) 'try to reduce tolerance'
       stop
    endif    
  end subroutine easydverk
 
 

  subroutine tunedverk(n,fcn,x,y,xend,tol,extradata)
    implicit none
     
    integer :: n,ind
    real(kp) :: x,y(n),xend,tol,c(24),w(n,9)    
    type(transfert), optional,intent(inout) :: extradata    

!    external :: fcn

    include 'inftools.h'

   

    ind=2
    c = 0._kp
    c(3) = epsilon(1._kp)
    c(4) = epsilon(1._kp)
 
    call dverk(n,fcn,x,y,xend,tol,ind,c,n,w,extradata)  
    
    if (ind.ne.3) then
       write(*,*) 'unsanedverk: stop ind = ',ind
       write(*,*) 'desesperate accuracy unreachable...'
       write(*,*) 'try tuning c(4) and c(6) in inftools:unsanedverk!'
       stop
    endif    
  end subroutine tunedverk
 
  




  subroutine diydverk(n,fcn,x,y,xend,tol,ind,c,extradata)
    implicit none
     
    integer :: n,ind
    real(kp) :: x,y(n),xend,tol,c(24),w(n,9)    
    type(transfert), optional, intent(inout) :: extradata    

!    external :: fcn

    include 'inftools.h'

    ind=2
    c = 0._kp

    c(3) = epsilon(1._kp)    
    c(4) = 0.001_kp

    call dverk(n,fcn,x,y,xend,tol,ind,c,n,w,extradata)  
    
    if ((ind.ne.3).and.(ind.ne.7)) then
       write(*,*) 'diydverk: stop ind = ',ind
       stop
    endif    
  end subroutine diydverk




!  function easyfixpnf(xstart,f,fjac,extradata)    
!    use hompack, only : rhodum, rhojacdum,fixpnf    
!    implicit none
   
!    real(kp) :: xstart, easyfixpnf
!    type(transfert), optional :: extradata


!    real(kp), parameter :: tolZero = 1e-5

!    integer, parameter :: n = 1
!    integer :: iflag, trace, nfe
!    real(kp) :: arcre,arcae
!    real(kp) :: ansre,ansae
!    real(kp) :: arclen

!    real(kp), dimension(n) :: a
!    real(kp), dimension(n+1) :: y,yp,yold,ypold

!    real(kp), dimension(1:N,1:N+2) :: qr
!    real(kp), dimension(n) :: alpha
!    real(kp), dimension(n+1) :: tz,w,wp,z0,z1

!    real(kp), dimension(8) :: sspar
   
!    integer, dimension(n+1) :: pivot

!    include 'hominfpack.h'

!    iflag = -1
!    trace = -1
!    sspar = -1

!    arcre = -1.
!    arcae = -1.
!    ansre = tolZero
!    ansae = tolZero

    
!    y(n+1) = xstart

!    call fixpnf(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,NFE, &
!         ARCLEN,YP,YOLD,YPOLD,QR,ALPHA,TZ,PIVOT,W,WP,Z0,Z1,SSPAR, &
!         f,fjac,rhodum,rhojacdum,extradata)

!    if (iflag.ne.1) then
!       write(*,*)'easyzero: iflag = ',iflag
!       stop
!    endif

!    easyfixpnf = y(n+1)

!  end function easyfixpnf







  subroutine dverk(n,fcn,x,y,xend,tol,ind,c,nw,w,extradata)
    implicit none
    integer :: n, ind, nw, k
    real(kp) :: x, y(n), xend, tol, c(*), w(nw,9), temp
    type(transfert), optional, intent(inout) :: extradata
    

!
!***********************************************************************
!                                                                      *
! note added 11/14/85.                                                 *
!                                                                      *
! if you discover any errors in this subroutine, please contact        *
!                                                                      *
!        kenneth r. jackson                                            *
!        department of computer science                                *
!        university of toronto                                         *
!        toronto, ontario,                                             *
!        canada   m5s 1a4                                              *
!                                                                      *
!        phone: 416-978-7075                                           *
!                                                                      *
!        electroni! mail:                                              *
!        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
!        csnet:  krj@toronto                                           *
!        arpa:   krj.toronto@csnet-relay                               *
!        bitnet: krj%toronto@csnet-relay.arpa                          *
!                                                                      *
! dverk is written in fortran 66.                                      *
!                                                                      *
! the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
! set for a  vax  in  double  precision.  they  should  be  reset,  as *
! described below, if this program is run on another machine.          *
!                                                                      *
! the c array is declared in this subroutine to have one element only, *
! although  more  elements  are  referenced  in this subroutine.  this *
! causes some compilers to issue warning messages.  there is,  though, *
! no  error  provided  c is declared sufficiently large in the calling *
! program, as described below.                                         *
!                                                                      *
! the following external statement  for  fcn  was  added  to  avoid  a *
! warning  message  from  the  unix  f77 compiler.  the original dverk *
! comments and code follow it.                                         *
!                                                                      *
!***********************************************************************
!

!EXTRADATA

!might be dangerous (xlf90) with optional argument, rather use explicit
!interface
!      external fcn
    include 'inftools.h'

    if (present(extradata)) extradata%update = .false.
!
!***********************************************************************
!                                                                      *
!     purpose - this is a runge-kutta  subroutine  based  on  verner's *
! fifth and sixth order pair of formulas for finding approximations to *
! the solution of  a  system  of  first  order  ordinary  differential *
! equations  with  initial  conditions. it attempts to keep the global *
! error proportional to  a  tolerance  specified  by  the  user.  (the *
! proportionality  depends  on the kind of error control that is used, *
! as well as the differential equation and the range of integration.)  *
!                                                                      *
!     various options are available to the user,  including  different *
! kinds  of  error control, restrictions on step sizes, and interrupts *
! which permit the user to examine the state of the  calculation  (and *
! perhaps make modifications) during intermediate stages.              *
!                                                                      *
!     the program is efficient for non-stiff systems.  however, a good *
! variable-order-adams  method  will probably be more efficient if the *
! function evaluations are very costly.  such a method would  also  be *
! more suitable if one wanted to obtain a large number of intermediate *
! solution values by interpolation, as might be the case  for  example *
! with graphical output.                                               *
!                                                                      *
!                                    hull-enright-jackson   1/10/76    *
!                                                                      *
!***********************************************************************
!                                                                      *
!     use - the user must specify each of the following                *
!                                                                      *
!     n  number of equations                                           *
!                                                                      *
!   fcn  name of subroutine for evaluating functions - the  subroutine *
!           itself must also be provided by the user - it should be of *
!           the following form                                         *
!              subroutine fcn(n, x, y, yprime)                         *
!              integer n                                               *
!              double precision x, y(n), yprime(n)                     *
!                      *** etc ***                                     *
!           and it should evaluate yprime, given n, x and y            *
!                                                                      *
!     x  independent variable - initial value supplied by user         *
!                                                                      *
!     y  dependent variable - initial values of components y(1), y(2), *
!           ..., y(n) supplied by user                                 *
!                                                                      *
!  xend  value of x to which integration is to be carried out - it may *
!           be less than the initial value of x                        *
!                                                                      *
!   tol  tolerance - the subroutine attempts to control a norm of  the *
!           local  error  in  such  a  way  that  the  global error is *
!           proportional to tol. in some problems there will be enough *
!           damping  of  errors, as well as some cancellation, so that *
!           the global error will be less than tol. alternatively, the *
!           control   can   be  viewed  as  attempting  to  provide  a *
!           calculated value of y at xend which is the exact  solution *
!           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
!           is proportional to tol.  (the norm  is  a  max  norm  with *
!           weights  that  depend on the error control strategy chosen *
!           by the user.  the default weight for the k-th component is *
!           1/max(1,abs(y(k))),  which therefore provides a mixture of *
!           absolute and relative error control.)                      *
!                                                                      *
!   ind  indicator - on initial entry ind must be set equal to  either *
!           1  or  2. if the user does not wish to use any options, he *
!           should set ind to 1 - all that remains for the user to  do *
!           then  is  to  declare c and w, and to specify nw. the user *
!           may also  select  various  options  on  initial  entry  by *
!           setting ind = 2 and initializing the first 9 components of *
!           c as described in the next section.  he may also  re-enter *
!           the  subroutine  with ind = 3 as mentioned again below. in *
!           any event, the subroutine returns with ind equal to        *
!              3 after a normal return                                 *
!              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
!              -1, -2, or -3 after an error condition (see below)      *
!                                                                      *
!     c  communications vector - the dimension must be greater than or *
!           equal to 24, unless option c(1) = 4 or 5 is used, in which *
!           case the dimension must be greater than or equal to n+30   *
!                                                                      *
!    nw  first dimension of workspace w -  must  be  greater  than  or *
!           equal to n                                                 *
!                                                                      *
!     w  workspace matrix - first dimension must be nw and second must *
!           be greater than or equal to 9                              *
!                                                                      *
!     the subroutine  will  normally  return  with  ind  =  3,  having *
! replaced the initial values of x and y with, respectively, the value *
! of xend and an approximation to y at xend.  the  subroutine  can  be *
! called  repeatedly  with new values of xend without having to change *
! any other argument.  however, changes in tol, or any of the  options *
! described below, may also be made on such a re-entry if desired.     *
!                                                                      *
!     three error returns are also possible, in which  case  x  and  y *
! will be the most recently accepted values -                          *
!     with ind = -3 the subroutine was unable  to  satisfy  the  error *
!        requirement  with a particular step-size that is less than or *
!        equal to hmin, which may mean that tol is too small           *
!     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
!        probably  means  that the requested tol (which is used in the *
!        calculation of hmin) is too small                             *
!     with ind = -1 the allowed maximum number of fcn evaluations  has *
!        been  exceeded,  but  this  can only occur if option c(7), as *
!        described in the next section, has been used                  *
!                                                                      *
!     there are several circumstances that will cause the calculations *
! to  be  terminated,  along with output of information that will help *
! the user determine the cause of  the  trouble.  these  circumstances *
! involve  entry with illegal or inconsistent values of the arguments, *
! such as attempting a normal  re-entry  without  first  changing  the *
! value of xend, or attempting to re-enter with ind less than zero.    *
!                                                                      *
!***********************************************************************
!                                                                      *
!     options - if the subroutine is entered with ind = 1, the first 9 *
! components of the communications vector are initialized to zero, and *
! the subroutine uses only default values  for  each  option.  if  the *
! subroutine  is  entered  with ind = 2, the user must specify each of *
! these 9 components - normally he would first set them all  to  zero, *
! and  then  make  non-zero  those  that  correspond to the particular *
! options he wishes to select. in any event, options may be changed on *
! re-entry  to  the  subroutine  -  but if the user changes any of the *
! options, or tol, in the course of a calculation he should be careful *
! about  how  such changes affect the subroutine - it may be better to *
! restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
! program  -  the information is available to the user, but should not *
! normally be changed by him.)                                         *
!                                                                      *
!  c(1)  error control indicator - the norm of the local error is  the *
!           max  norm  of  the  weighted  error  estimate  vector, the *
!           weights being determined according to the value of c(1) -  *
!              if c(1)=1 the weights are 1 (absolute error control)    *
!              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
!                 control)                                             *
!              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
!                 (relative  error  control,  unless abs(y(k)) is less *
!                 than the floor value, abs(c(2)) )                    *
!              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
!                 (here individual floor values are used)              *
!              if c(1)=5 the weights are 1/abs(c(k+30))                *
!              for all other values of c(1), including  c(1) = 0,  the *
!                 default  values  of  the  weights  are  taken  to be *
!                 1/max(1,abs(y(k))), as mentioned earlier             *
!           (in the two cases c(1) = 4 or 5 the user must declare  the *
!           dimension of c to be at least n+30 and must initialize the *
!           components c(31), c(32), ..., c(n+30).)                    *
!                                                                      *
!  c(2)  floor value - used when the indicator c(1) has the value 3    *
!                                                                      *
!  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
!           to be abs(c(3)) - otherwise it uses the default value      *
!              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
!           where dwarf is a very small positive  machine  number  and *
!           rreb is the relative roundoff error bound                  *
!                                                                      *
!  c(4)  hstart specification - if not zero, the subroutine  will  use *
!           an  initial  hmag equal to abs(c(4)), except of course for *
!           the restrictions imposed by hmin and hmax  -  otherwise it *
!           uses the default value of hmax*(tol)**(1/6)                *
!                                                                      *
!  c(5)  scale specification - this is intended to be a measure of the *
!           scale of the problem - larger values of scale tend to make *
!           the method more reliable, first  by  possibly  restricting *
!           hmax  (as  described  below) and second, by tightening the *
!           acceptance requirement - if c(5) is zero, a default  value *
!           of  1  is  used.  for  linear  homogeneous  problems  with *
!           constant coefficients, an appropriate value for scale is a *
!           norm  of  the  associated  matrix.  for other problems, an *
!           approximation to  an  average  value  of  a  norm  of  the *
!           jacobian along the trajectory may be appropriate           *
!                                                                      *
!  c(6)  hmax specification - four cases are possible                  *
!           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
!              min(abs(c(6)),2/abs(c(5)))                              *
!           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
!           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
!              2/abs(c(5))                                             *
!           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
!              of 2                                                    *
!                                                                      *
!  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
!           error  return with ind = -1 will be caused when the number *
!           of function evaluations exceeds abs(c(7))                  *
!                                                                      *
!  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
!           interrupt   the  calculations  after  it  has  chosen  its *
!           preliminary value of hmag, and just before choosing htrial *
!           and  xtrial  in  preparation for taking a step (htrial may *
!           differ from hmag in sign, and may  require  adjustment  if *
!           xend  is  near) - the subroutine returns with ind = 4, and *
!           will resume calculation at the point  of  interruption  if *
!           re-entered with ind = 4                                    *
!                                                                      *
!  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
!           interrupt   the  calculations  immediately  after  it  has *
!           decided whether or not to accept the result  of  the  most *
!           recent  trial step, with ind = 5 if it plans to accept, or *
!           ind = 6 if it plans to reject -  y(*)  is  the  previously *
!           accepted  result, while w(*,9) is the newly computed trial *
!           value, and w(*,2) is the unweighted error estimate vector. *
!           the  subroutine  will  resume calculations at the point of *
!           interruption on re-entry with ind = 5 or 6. (the user  may *
!           change ind in this case if he wishes, for example to force *
!           acceptance of a step that would otherwise be rejected,  or *
!           vice versa. he can also restart with ind = 1 or 2.)        *
!                                                                      *
!***********************************************************************
!                                                                      *
!  summary of the components of the communications vector              *
!                                                                      *
!     prescribed at the option       determined by the program         *
!           of the user                                                *
!                                                                      *
!                                    c(10) rreb(rel roundoff err bnd)  *
!     c(1) error control indicator   c(11) dwarf (very small mach no)  *
!     c(2) floor value               c(12) weighted norm y             *
!     c(3) hmin specification        c(13) hmin                        *
!     c(4) hstart specification      c(14) hmag                        *
!     c(5) scale specification       c(15) scale                       *
!     c(6) hmax specification        c(16) hmax                        *
!     c(7) max no of fcn evals       c(17) xtrial                      *
!     c(8) interrupt no 1            c(18) htrial                      *
!     c(9) interrupt no 2            c(19) est                         *
!                                    c(20) previous xend               *
!                                    c(21) flag for xend               *
!                                    c(22) no of successful steps      *
!                                    c(23) no of successive failures   *
!                                    c(24) no of fcn evals             *
!                                                                      *
!  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
!                                                                      *
!***********************************************************************
!                                                                      *
!  an overview of the program                                          *
!                                                                      *
!     begin initialization, parameter checking, interrupt re-entries   *
!  ......abort if ind out of range 1 to 6                              *
!  .     cases - initial entry, normal re-entry, interrupt re-entries  *
!  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
!  v........abort if n.gt.nw or tol.le.0                               *
!  .        if initial entry without options (ind .eq. 1)              *
!  .           set c(1) to c(9) equal to zero                          *
!  .        else initial entry with options (ind .eq. 2)               *
!  .           make c(1) to c(9) non-negative                          *
!  .           make floor values non-negative if they are to be used   *
!  .        end if                                                     *
!  .        initialize rreb, dwarf, prev xend, flag, counts            *
!  .     case 2 - normal re-entry (ind .eq. 3)                         *
!  .........abort if xend reached, and either x changed or xend not    *
!  .        re-initialize flag                                         *
!  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
!  v        transfer control to the appropriate re-entry point.......  *
!  .     end cases                                                  .  *
!  .  end initialization, etc.                                      .  *
!  .                                                                v  *
!  .  loop through the following 4 stages, once for each trial step .  *
!  .     stage 1 - prepare                                          .  *
!***********error return (with ind=-1) if no of fcn evals too great .  *
!  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
!  .        calc hmin, scale, hmax                                  .  *
!***********error return (with ind=-2) if hmin .gt. hmax            .  *
!  .        calc preliminary hmag                                   .  *
!***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
!  .        calc hmag, xtrial and htrial                            .  *
!  .     end stage 1                                                .  *
!  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
!  .     stage 3 - calc the error estimate                          .  *
!  .     stage 4 - make decisions                                   .  *
!  .        set ind=5 if step acceptable, else set ind=6            .  *
!***********interrupt no 2 if requested....................re-entry.v  *
!  .        if step accepted (ind .eq. 5)                              *
!  .           update x, y from xtrial, ytrial                         *
!  .           add 1 to no of successful steps                         *
!  .           set no of successive failures to zero                   *
!**************return(with ind=3, xend saved, flag set) if x .eq. xend *
!  .        else step not accepted (ind .eq. 6)                        *
!  .           add 1 to no of successive failures                      *
!**************error return (with ind=-3) if hmag .le. hmin            *
!  .        end if                                                     *
!  .     end stage 4                                                   *
!  .  end loop                                                         *
!  .                                                                   *
!  begin abort action                                                  *
!     output appropriate  message  about  stopping  the  calculations, *
!        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
!        previous xend,  no of  successful  steps,  no  of  successive *
!        failures, no of fcn evals, and the components of y            *
!     stop                                                             *
!  end abort action                                                    *
!                                                                      *
!***********************************************************************
!
!     ******************************************************************
!     * begin initialization, parameter checking, interrupt re-entries *
!     ******************************************************************
!
!  ......abort if ind out of range 1 to 6
         if (ind.lt.1 .or. ind.gt.6) go to 500
!
!        cases - initial entry, normal re-entry, interrupt re-entries
         go to (5, 5, 45, 1111, 2222, 2222), ind
!        case 1 - initial entry (ind .eq. 1 or 2)
!  .........abort if n.gt.nw or tol.le.0
    5       if (n.gt.nw .or. tol.le.0._kp) go to 500
            if (ind.eq. 2) go to 15
!              initial entry without options (ind .eq. 1)
!              set c(1) to c(9) equal to 0
               do 10 k = 1, 9
                  c(k) = 0._kp
   10          continue
               go to 35
   15       continue
!              initial entry with options (ind .eq. 2)
!              make c(1) to c(9) non-negative
               do 20 k = 1, 9
                  c(k) = abs(c(k))
   20          continue
!              make floor values non-negative if they are to be used
               if (c(1).ne.4._kp .and. c(1).ne.5._kp) go to 30
                  do 25 k = 1, n
                     c(k+30) = abs(c(k+30))
   25             continue
   30          continue
   35       continue
!           initialize rreb, dwarf, prev xend, flag, counts
            c(10) = 2._kp**(-56)
            c(11) = 1.d-35
!           set previous xend initially to initial value of x
            c(20) = x
            do 40 k = 21, 24
               c(k) = 0._kp
   40       continue
            go to 50
!        case 2 - normal re-entry (ind .eq. 3)
!  .........abort if xend reached, and either x changed or xend not
   45       if (c(21).ne.0._kp .and. &
                 (x.ne.c(20) .or. xend.eq.c(20))) go to 500
!           re-initialize flag
            c(21) = 0._kp
            go to 50
!        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
!           transfer control to the appropriate re-entry point..........
!           this has already been handled by the computed go to        .
!        end cases                                                     v
   50    continue
!
!     end initialization, etc.
!
!     ******************************************************************
!     * loop through the following 4 stages, once for each trial  step *
!     * until the occurrence of one of the following                   *
!     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
!     *        stage 4                                                 *
!     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
!     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
!     *        requested, in stage 1 or stage 4                        *
!     ******************************************************************
!
99999 continue
!
!        ***************************************************************
!        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!        * and some parameter  checking,  and  end  up  with  suitable *
!        * values of hmag, xtrial and htrial in preparation for taking *
!        * an integration step.                                        *
!        ***************************************************************
!
!***********error return (with ind=-1) if no of fcn evals too great
            if (c(7).eq.0._kp .or. c(24).lt.c(7)) go to 100
               ind = -1
               return
  100       continue
!
!           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
            if (ind .eq. 6) go to 105
!EXTRADATA call fcn(n, x, y, w(1,1))
               if (present(extradata)) extradata%check=.true.               
               call fcn(n, x, y, w(1,1),extradata)             
               if (present(extradata)) then
                  if (extradata%update) xend = extradata%xend
                  extradata%check=.false.                  
               endif
!END EXTRADATA
               c(24) = c(24) + 1._kp
  105       continue
!
!           calculate hmin - use default unless value prescribed
            c(13) = c(3)
            if (c(3) .ne. 0._kp) go to 165
!              calculate default value of hmin
!              first calculate weighted norm y - c(12) - as specified
!              by the error control indicator c(1)
               temp = 0._kp
               if (c(1) .ne. 1._kp) go to 115
!                 absolute error control - weights are 1
                  do 110 k = 1, n
                     temp = max(temp, abs(y(k)))
  110             continue
                  c(12) = temp
                  go to 160
  115          if (c(1) .ne. 2._kp) go to 120
!                 relative error control - weights are 1/abs(y(k)) so
!                 weighted norm y is 1
                  c(12) = 1._kp
                  go to 160
  120          if (c(1) .ne. 3._kp) go to 130
!                 weights are 1/max(c(2),abs(y(k)))
                  do 125 k = 1, n
                     temp = max(temp, abs(y(k))/c(2))
  125             continue
                  c(12) = min(temp, 1._kp)
                  go to 160
  130          if (c(1) .ne. 4._kp) go to 140
!                 weights are 1/max(c(k+30),abs(y(k)))
                  do 135 k = 1, n
                     temp = max(temp, abs(y(k))/c(k+30))
  135             continue
                  c(12) = min(temp, 1._kp)
                  go to 160
  140          if (c(1) .ne. 5._kp) go to 150
!                 weights are 1/c(k+30)
                  do 145 k = 1, n
                     temp = max(temp, abs(y(k))/c(k+30))
  145             continue
                  c(12) = temp
                  go to 160
  150          continue
!                 default case - weights are 1/max(1,abs(y(k)))
                  do 155 k = 1, n
                     temp = max(temp, abs(y(k)))
  155             continue
                  c(12) = min(temp, 1._kp)
  160          continue
               c(13) = 10._kp*max(c(11),c(10)*max(c(12)/tol,abs(x)))
  165       continue
!
!           calculate scale - use default unless value prescribed
            c(15) = c(5)
            if (c(5) .eq. 0._kp) c(15) = 1._kp
!
!           calculate hmax - consider 4 cases
!           case 1 both hmax and scale prescribed
               if (c(6).ne.0._kp .and. c(5).ne.0._kp) &
                    c(16) = min(c(6), 2._kp/c(5))
!           case 2 - hmax prescribed, but scale not
               if (c(6).ne.0._kp .and. c(5).eq.0._kp) c(16) = c(6)
!           case 3 - hmax not prescribed, but scale is
               if (c(6).eq.0._kp .and. c(5).ne.0._kp) c(16) = 2._kp/c(5)
!           case 4 - neither hmax nor scale is provided
               if (c(6).eq.0._kp .and. c(5).eq.0._kp) c(16) = 2._kp
!
!***********error return (with ind=-2) if hmin .gt. hmax
            if (c(13) .le. c(16)) go to 170
               ind = -2
               return
  170       continue
!
!           calculate preliminary hmag - consider 3 cases
            if (ind .gt. 2) go to 175
!           case 1 - initial entry - use prescribed value of hstart, if
!              any, else default
               c(14) = c(4)
               if (c(4) .eq. 0._kp) c(14) = c(16)*tol**(1./6.)
               go to 185
  175       if (c(23) .gt. 1._kp) go to 180
!           case 2 - after a successful step, or at most  one  failure,
!              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
!              overflow. then avoid reduction by more than half.
               temp = 2._kp*c(14)
               if (tol .lt. (2._kp/.9_kp)**6*c(19)) &
                    temp = .9_kp*(tol/c(19))**(1./6.)*c(14)
               c(14) = max(temp, .5_kp*c(14))
               go to 185
  180       continue
!           case 3 - after two or more successive failures
               c(14) = .5_kp*c(14)
  185       continue
!
!           check against hmax
            c(14) = min(c(14), c(16))
!
!           check against hmin
            c(14) = max(c(14), c(13))
!
!***********interrupt no 1 (with ind=4) if requested
            if (c(8) .eq. 0._kp) go to 1111
               ind = 4
               return
!           resume here on re-entry with ind .eq. 4   ........re-entry..
 1111       continue
!
!           calculate hmag, xtrial - depending on preliminary hmag, xend
            if (c(14) .ge. abs(xend - x)) go to 190
!              do not step more than half way to xend
               c(14) = min(c(14), .5_kp*abs(xend - x))
               c(17) = x + sign(c(14), xend - x)
               go to 195
  190       continue
!              hit xend exactly
               c(14) = abs(xend - x)
               c(17) = xend
  195       continue
!
!           calculate htrial
            c(18) = c(17) - x
!
!        end stage 1
!
!        ***************************************************************
!        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!        * stage 3. w(*,9) is temporary storage until finally it holds *
!        * ytrial.                                                     *
!        ***************************************************************
!
            temp = c(18)/1398169080000._kp
!
            do 200 k = 1, n
               w(k,9) = y(k) + temp*w(k,1)*233028180000._kp
  200       continue
            call fcn(n, x + c(18)/6._kp, w(1,9), w(1,2),extradata)
!
            do 205 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*74569017600._kp &
                    + w(k,2)*298276070400._kp  )
  205       continue
            call fcn(n, x + c(18)*(4._kp/15._kp), w(1,9), w(1,3),extradata)
!
            do 210 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*1165140900000._kp &
                    - w(k,2)*3728450880000._kp &
                    + w(k,3)*3495422700000._kp )
  210       continue
            call fcn(n, x + c(18)*(2._kp/3._kp), w(1,9), w(1,4),extradata)
!
            do 215 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*3604654659375._kp &
                    + w(k,2)*12816549900000._kp &
                    - w(k,3)*9284716546875._kp &
                    + w(k,4)*1237962206250._kp )
  215       continue
            call fcn(n, x + c(18)*(5._kp/6._kp), w(1,9), w(1,5),extradata)
!
            do 220 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*3355605792000._kp &
                    - w(k,2)*11185352640000._kp &
                    + w(k,3)*9172628850000._kp &
                    - w(k,4)*427218330000._kp &
                    + w(k,5)*482505408000._kp  )
  220       continue
            call fcn(n, x + c(18), w(1,9), w(1,6),extradata)
!
            do 225 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*770204740536._kp &
                    + w(k,2)*2311639545600._kp &
                    - w(k,3)*1322092233000._kp &
                    - w(k,4)*453006781920._kp &
                    + w(k,5)*326875481856._kp  )
  225       continue
            call fcn(n, x + c(18)/15._kp, w(1,9), w(1,7),extradata)
!
            do 230 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*2845924389000._kp &
                    - w(k,2)*9754668000000._kp &
                    + w(k,3)*7897110375000._kp &
                    - w(k,4)*192082660000._kp &
                    + w(k,5)*400298976000._kp &
                    + w(k,7)*201586000000._kp  )
  230       continue
            call fcn(n, x + c(18), w(1,9), w(1,8),extradata)
!
!           calculate ytrial, the extrapolated approximation and store
!              in w(*,9)
            do 235 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*104862681000._kp &
                    + w(k,3)*545186250000._kp &
                    + w(k,4)*446637345000._kp &
                    + w(k,5)*188806464000._kp &
                    + w(k,7)*15076875000._kp &
                    + w(k,8)*97599465000._kp   )
  235       continue
!
!           add 7 to the no of fcn evals
            c(24) = c(24) + 7._kp
!
!        end stage 2
!
!        ***************************************************************
!        * stage 3 - calculate the error estimate est. first calculate *
!        * the  unweighted  absolute  error  estimate vector (per unit *
!        * step) for the unextrapolated approximation and store it  in *
!        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!        * specified by the error  control  indicator  c(1).  finally, *
!        * modify  this result to produce est, the error estimate (per *
!        * unit step) for the extrapolated approximation ytrial.       *
!        ***************************************************************
!
!           calculate the unweighted absolute error estimate vector
            do 300 k = 1, n
               w(k,2) = (   w(k,1)*8738556750._kp &
                    + w(k,3)*9735468750._kp &
                    - w(k,4)*9709507500._kp &
                    + w(k,5)*8582112000._kp &
                    + w(k,6)*95329710000._kp &
                    - w(k,7)*15076875000._kp &
                    - w(k,8)*97599465000._kp)/1398169080000._kp
  300       continue
!
!           calculate the weighted max norm of w(*,2) as specified by
!           the error control indicator c(1)
            temp = 0._kp
            if (c(1) .ne. 1._kp) go to 310
!              absolute error control
               do 305 k = 1, n
                  temp = max(temp,abs(w(k,2)))
  305          continue
               go to 360
  310       if (c(1) .ne. 2._kp) go to 320
!              relative error control
               do 315 k = 1, n
                  temp = max(temp, abs(w(k,2)/y(k)))
  315          continue
               go to 360
  320       if (c(1) .ne. 3._kp) go to 330
!              weights are 1/max(c(2),abs(y(k)))
               do 325 k = 1, n
                  temp = max(temp, abs(w(k,2)) &
                       / max(c(2), abs(y(k))) )
  325          continue
               go to 360
  330       if (c(1) .ne. 4._kp) go to 340
!              weights are 1/max(c(k+30),abs(y(k)))
               do 335 k = 1, n
                  temp = max(temp, abs(w(k,2)) &
                       / max(c(k+30), abs(y(k))) )
  335          continue
               go to 360
  340       if (c(1) .ne. 5._kp) go to 350
!              weights are 1/c(k+30)
               do 345 k = 1, n
                  temp = max(temp, abs(w(k,2)/c(k+30)))
  345          continue
               go to 360
  350       continue
!              default case - weights are 1/max(1,abs(y(k)))
               do 355 k = 1, n
                  temp = max(temp, abs(w(k,2)) &
                       / max(1._kp, abs(y(k))) )
  355          continue
  360       continue
!
!           calculate est - (the weighted max norm of w(*,2))*hmag*scale
!              - est is intended to be a measure of the error  per  unit
!              step in ytrial
            c(19) = temp*c(14)*c(15)
!
!        end stage 3
!
!        ***************************************************************
!        * stage 4 - make decisions.                                   *
!        ***************************************************************
!
!           set ind=5 if step acceptable, else set ind=6
            ind = 5
            if (c(19) .gt. tol) ind = 6
!
!***********interrupt no 2 if requested
            if (c(9) .eq. 0._kp) go to 2222
               return
!           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
 2222       continue
!
            if (ind .eq. 6) go to 410
!              step accepted (ind .eq. 5), so update x, y from xtrial,
!                 ytrial, add 1 to the no of successful steps, and set
!                 the no of successive failures to zero
               x = c(17)
               do 400 k = 1, n
                  y(k) = w(k,9)
  400          continue
               c(22) = c(22) + 1._kp
               c(23) = 0._kp
!**************return(with ind=3, xend saved, flag set) if x .eq. xend
               if (x .ne. xend) go to 405
                  ind = 3
                  c(20) = xend
                  c(21) = 1._kp
                  return
  405          continue
               go to 420
  410       continue
!              step not accepted (ind .eq. 6), so add 1 to the no of
!                 successive failures
               c(23) = c(23) + 1._kp
!**************error return (with ind=-3) if hmag .le. hmin
               if (c(14) .gt. c(13)) go to 415
                  ind = -3
                  return
  415          continue
  420       continue
!
!        end stage 4
!
      go to 99999
!     end loop
!
!  begin abort action
  500 continue
      write(6,*)'Computation stopped in dverk with'
      write(6,*)'ind= tol= ',ind,tol
      write(6,*)'x= n= ',x,n
      write(6,*)'c(13)= xend= ',c(13),xend
      write(6,*)'nw= c(16)= c(20)= ',nw, c(16),c(20)
      write(6,*)'c(22)= c(23)= c(24)= ',c(22),c(23),c(24)
      write(6,*)'y(:)= ',y
!      write(6,*) ind, tol, x, n, c(13), xend, nw, c(16), c(20), &
!           c(22), c(23), c(24), (y(k), k = 1, n)
!  505 format( /// 1h0, 58hcomputation stopped in dverk with the following values - / 1h0, 5hind =, i4, 5x, 6htol  =, 1pd13.6, 5x, 11hx         =,&
!           1pd22.15&
!           / 1h , 5hn   =, i4, 5x, 6hhmin =, 1pd13.6, 5x, 11hxend      =,&
!           1pd22.15&
!           / 1h , 5hnw  =, i4, 5x, 6hhmax =, 1pd13.6, 5x, 11hprev xend =,&
!           1pd22.15&
!           / 1h0, 14x, 27hno of successful steps    =, 0pf8.0&
!           / 1h , 14x, 27hno of successive failures =, 0pf8.0&
!           / 1h , 14x, 27hno of function evals      =, 0pf8.0&
!           / 1h0, 23hthe components of y are&
!           // (1h , 1p5d24.15)                                           )
!
      stop
!
!  end abort action
!
    end subroutine dverk
    


    
    function newton(func,dfunc,xguess,tol,extradata,success)
      implicit none
      real(kp) :: newton, xguess, tol
      type(transfert) :: extradata
      logical, optional :: success
      
!successive over-relaxation parameter      
      real(kp), parameter :: omega = 1._kp
      
      real(kp) :: xp, x, df

      integer :: count
      integer, parameter :: maxcount = 1000000
      
      interface
         function func(x,otherdata)
           use infprec, only : kp,transfert         
           implicit none                     
           real(kp) :: func
           real(kp), intent(in) :: x
           type(transfert), optional, intent(inout) :: otherdata
         end function func

         function dfunc(x,otherdata)
           use infprec, only : kp,transfert         
           implicit none                     
           real(kp) :: dfunc
           real(kp), intent(in) :: x
           type(transfert), optional, intent(inout) :: otherdata
         end function dfunc                  
      end interface

      xp = xguess

      count = 0
      
      do

         df = dfunc(xp,extradata)
         if (abs(df).lt.tol) df = sign(tol,df)
         
         x = xp - omega*func(xp,extradata)/df

         if (abs(x-xp).lt.tol) then

            if (present(success)) success = .true.
            exit

         else

            xp = x
            count = count + 1
            
            if (count.gt.maxcount) then
               write(*,*)'newton: iteration count exceeded!'
               write(*,*)'x= f(x)= ',x,func(x,extradata)
               if (present(success)) then
                  success = .false.
               else
                  stop
               endif
               exit
            endif
            
         endif
         
      end do

      newton = x
      
    end function newton



    


    function zbrent(func,x1,x2,tol,extradata)
      INTEGER ITMAX
      real(kp) zbrent,tol,x1,x2,EPS
      type(transfert) :: extradata 
      
      !real(kp) :: func
!      EXTERNAL func
      PARAMETER (ITMAX=1000000,EPS=100._kp*epsilon(1._kp))
      INTEGER iter
      real(kp) a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm      
      
      integer, parameter :: iexmax = 1000
      real(kp), parameter :: tolExpand = 1e-2
      integer :: iex
      logical :: notbracketed

      logical :: display
      
      interface
         function func(x,otherdata)
           use infprec, only : kp,transfert         
           implicit none                     
           real(kp) :: func
           real(kp), intent(in) :: x
           type(transfert), optional, intent(inout) :: otherdata
         end function func
      end interface

      if (x1.gt.x2) then
         stop 'zbrent: min > max on input!'
      endif
      

      notbracketed = .true.
      display = .true.
      iex = 0

      a=x1
      b=x2      

      do while (notbracketed.and.(iex.lt.iexmax))
         fa=func(a,extradata)
         fb=func(b,extradata)
         if((fa.gt.0._kp.and.fb.gt.0._kp).or.(fa.lt.0._kp.and.fb.lt.0._kp)) then
            if (display) then
               write(*,*)
               write(*,*)'------------------------------------------------------'
               write(*,*)'x1=',a,'f(x1)=',fa
               write(*,*)'x2=',b,'f(x2)=',fb            
               write(*,*)'zbrent: expanding interval!'
            endif
            a = a - abs(a)*tolExpand
            b = b + abs(b)*tolExpand
            iex = iex + 1
            notbracketed = .true.
            display = .false.
         else
            notbracketed = .false.
         endif
      enddo

      if (notbracketed) then
         write(*,*)'x1=',a,'f(x1)=',fa
         write(*,*)'x2=',b,'f(x2)=',fb            
         write(*,*)'zbrent: expanding interval!'
         stop 'root must be bracketed for zbrent'
      elseif(iex.ne.0) then
         write(*,*)
         write(*,*)'x1=',a,'f(x1)=',fa
         write(*,*)'x2=',b,'f(x2)=',fb            
         write(*,*)'zbrent: interval expansion number: ',iex
         write(*,*)'------------------------------------------------------'
         write(*,*)
      endif

      c=b
      fc=fb
      do iter=1,ITMAX
        if((fb.gt.0._kp.and.fc.gt.0._kp).or.(fb.lt.0._kp.and.fc.lt.0._kp))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2._kp*EPS*abs(b)+0.5_kp*tol
        xm=.5_kp*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0._kp)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2._kp*xm*s
            q=1._kp-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2._kp*xm*q*(q-r)-(b-a)*(r-1._kp))
            q=(q-1._kp)*(r-1._kp)*(s-1._kp)
          endif
          if(p.gt.0._kp) q=-q
          p=abs(p)
          if(2._kp*p .lt. min(3._kp*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b,extradata)
     enddo
     !stop 'zbrent exceeding maximum iterations'
     write(*,*) 'zbrent exceeding maximum iterations ITMAX=',ITMAX
     zbrent=b
     return
   end function zbrent



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                             !!
!!          SET OF ROUTINES TO SOLVE QUARTIC POLYNOMIAL        !!
!!                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!+
FUNCTION CubeRoot(x) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the Cube Root of a REAL(kp) number. If the argument is
!   negative, then the cube root is also negative.

  REAL(kp),INTENT(IN) :: x
  REAL(kp):: f
  COMPLEX(kp),PARAMETER:: CZERO=(0._kp,0._kp)
  REAL(kp),PARAMETER:: ZERO=0._kp, FOURTH=0.25_kp, HALF=0.5_kp
  REAL(kp),PARAMETER:: ONE=1.0_kp, TWO=2.0_kp, THREE=3.0_kp, FOUR=4.0_kp
  REAL(kp),PARAMETER:: EPS=EPSILON(ONE)
!----------------------------------------------------------------------------
  IF (x < ZERO) THEN
    f=-EXP(LOG(-x)/THREE)
  ELSE IF (x > ZERO) THEN
    f=EXP(LOG(x)/THREE)
  ELSE
    f=ZERO
  END IF
  RETURN
END Function CubeRoot   ! ---------------------------------------------------

!+
SUBROUTINE LinearRoot(a, z)
! ---------------------------------------------------------------------------
! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z 
!     AND STORES THE RESULTS IN Z. It is assumed that a(2) is non-zero.
  REAL(kp),INTENT(IN),DIMENSION(:):: a
  REAL(kp),INTENT(OUT):: z
  COMPLEX(kp),PARAMETER:: CZERO=(0._kp,0._kp)
  REAL(kp),PARAMETER:: ZERO=0._kp, FOURTH=0.25_kp, HALF=0.5_kp
  REAL(kp),PARAMETER:: ONE=1.0_kp, TWO=2.0_kp, THREE=3.0_kp, FOUR=4.0_kp
  REAL(kp),PARAMETER:: EPS=EPSILON(ONE)
!----------------------------------------------------------------------------
  IF (a(2)==0.0) THEN
    z=ZERO
  ELSE
    z=-a(1)/a(2)
  END IF
  RETURN
END Subroutine LinearRoot   ! -----------------------------------------------

!+
SUBROUTINE OneLargeTwoSmall(a1,a2,a4,w, z)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the roots of a cubic when one root, w, is known to be
!   much larger in magnitude than the other two

  REAL(kp),INTENT(IN):: a1,a2,a4
  REAL(kp),INTENT(IN):: w
  COMPLEX(kp),INTENT(OUT),DIMENSION(:):: z

  COMPLEX(kp),PARAMETER:: CZERO=(0._kp,0._kp)
  REAL(kp),PARAMETER:: ZERO=0._kp, FOURTH=0.25_kp, HALF=0.5_kp
  REAL(kp),PARAMETER:: ONE=1.0_kp, TWO=2.0_kp, THREE=3.0_kp, FOUR=4.0_kp


  REAL(kp),DIMENSION(3):: aq
!----------------------------------------------------------------------------
  aq(1)=a1
  aq(2)=a2+a1/w
  aq(3)=-a4*w
  CALL QuadraticRoots(aq, z)
  z(3)=CMPLX(w,ZERO,kp)
  
  IF (AIMAG(z(1)) == ZERO) RETURN
  z(3)=z(2)
  z(2)=z(1)
  z(1)=CMPLX(w,ZERO,kp)
  RETURN
END Subroutine OneLargeTwoSmall   ! -----------------------------------------

!+
SUBROUTINE QuadraticRoots(a, z)
! ---------------------------------------------------------------------------
! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2
!     AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(3) IS NONZERO.

  REAL(kp),INTENT(IN),DIMENSION(:):: a
  COMPLEX(kp),INTENT(OUT),DIMENSION(:):: z
  COMPLEX(kp),PARAMETER:: CZERO=(0._kp,0._kp)
  REAL(kp),PARAMETER:: ZERO=0._kp, FOURTH=0.25_kp, HALF=0.5_kp
  REAL(kp),PARAMETER:: ONE=1.0_kp, TWO=2.0_kp, THREE=3.0_kp, FOUR=4.0_kp
  REAL(kp),PARAMETER:: EPS=EPSILON(ONE)
  INTEGER :: outputCode



  REAL(kp):: d, r, w, x, y
!----------------------------------------------------------------------------
  IF(a(1)==0.0_kp) THEN     ! EPS is a global module constant (private)
    z(1) = CZERO               ! one root is obviously zero
    z(2) = CMPLX(-a(2)/a(3), ZERO,kp)    ! remainder is a linear eq.
    outputCode=21   ! two identical real roots
    RETURN
  END IF

  d = a(2)*a(2) - FOUR*a(1)*a(3)             ! the discriminant
  IF (ABS(d) <= TWO*eps*a(2)*a(2)) THEN
    z(1) = CMPLX(-HALF*a(2)/a(3), ZERO, kp) ! discriminant is tiny
    z(2) = z(1)
    outputCode=22  ! two distinct real roots
    RETURN
  END IF

  r = SQRT(ABS(d))
  IF (d < ZERO) THEN
    x = -HALF*a(2)/a(3)        ! negative discriminant => roots are complex   
    y = ABS(HALF*r/a(3))
    z(1) = CMPLX(x, y, kp)
    z(2) = CMPLX(x,-y, kp)   ! its conjugate
    outputCode=23                        !  COMPLEX ROOTS
    RETURN
  END IF

  IF (a(2) /= ZERO) THEN              ! see Numerical Recipes, sec. 5.5
    w = -(a(2) + SIGN(r,a(2)))
    z(1) = CMPLX(TWO*a(1)/w,  ZERO, kp)
    z(2) = CMPLX(HALF*w/a(3), ZERO, kp)
    outputCode=22           ! two real roots
    RETURN
  END IF

  x = ABS(HALF*r/a(3))   ! a(2)=0 if you get here
  z(1) = CMPLX( x, ZERO, kp)
  z(2) = CMPLX(-x, ZERO, kp)
  outputCode=22
  RETURN
END Subroutine QuadraticRoots   ! -------------------------------------------

!+
SUBROUTINE CubicRoots(a, z)
!----------------------------------------------------------------------------
! PURPOSE - Compute the roots of the real polynomial
!              A(1) + A(2)*Z + A(3)*Z**2 + A(4)*Z**3
  REAL(kp),INTENT(IN),DIMENSION(:):: a
  COMPLEX(kp),INTENT(OUT),DIMENSION(:):: z

  REAL(kp),PARAMETER:: RT3=sqrt(3._kp)    ! (Sqrt(3)
  REAL (kp) :: aq(3), arg, c, cf, d, p, p1, q, q1
  REAL(kp):: r, ra, rb, rq, rt
  REAL(kp):: r1, s, sf, sq, sum, t, tol, t1, w
  REAL(kp):: w1, w2, x, x1, x2, x3, y, y1, y2, y3

  COMPLEX(kp),PARAMETER:: CZERO=(0._kp,0._kp)
  REAL(kp),PARAMETER:: ZERO=0._kp, FOURTH=0.25_kp, HALF=0.5_kp
  REAL(kp),PARAMETER:: ONE=1.0_kp, TWO=2.0_kp, THREE=3.0_kp, FOUR=4.0_kp
  REAL(kp),PARAMETER:: EPS=EPSILON(ONE)


! NOTE -   It is assumed that a(4) is non-zero. No test is made here.
!----------------------------------------------------------------------------
  IF (a(1)==0.0_kp) THEN
    z(1) = CZERO  ! one root is obviously zero
    CALL QuadraticRoots(a(2:4), z(2:3))   ! remaining 2 roots here
    RETURN
  END IF

  p = a(3)/(THREE*a(4))
  q = a(2)/a(4)
  r = a(1)/a(4)
  tol = FOUR*EPS

  c = ZERO
  t = a(2) - p*a(3)
  IF (ABS(t) > tol*ABS(a(2))) c = t/a(4)

  t = TWO*p*p - q
  IF (ABS(t) <= tol*ABS(q)) t = ZERO
  d = r + p*t
  IF (ABS(d) <= tol*ABS(r)) GO TO 110

!           SET  SQ = (A(4)/S)**2 * (C**3/27 + D**2/4)

  s = MAX(ABS(a(1)), ABS(a(2)), ABS(a(3)))
  p1 = a(3)/(THREE*s)
  q1 = a(2)/s
  r1 = a(1)/s

  t1 = q - 2.25_KP*p*p
  IF (ABS(t1) <= tol*ABS(q)) t1 = ZERO
  w = FOURTH*r1*r1
  w1 = HALF*p1*r1*t
  w2 = q1*q1*t1/27.0_KP

  IF (w1 >= ZERO) THEN
    w = w + w1
    sq = w + w2
  ELSE IF (w2 < ZERO) THEN
    sq = w + (w1 + w2)
  ELSE
    w = w + w2
    sq = w + w1
  END IF

  IF (ABS(sq) <= tol*w) sq = ZERO
  rq = ABS(s/a(4))*SQRT(ABS(sq))
  IF (sq >= ZERO) GO TO 40

!                   ALL ROOTS ARE REAL

  arg = ATAN2(rq, -HALF*d)
  cf = COS(arg/THREE)
  sf = SIN(arg/THREE)
  rt = SQRT(-c/THREE)
  y1 = TWO*rt*cf
  y2 = -rt*(cf + rt3*sf)
  y3 = -(d/y1)/y2

  x1 = y1 - p
  x2 = y2 - p
  x3 = y3 - p

  IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)
  IF (ABS(x2) > ABS(x3)) CALL Swap(x2,x3)
  IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)

  w = x3

  IF (ABS(x2) < 0.1_KP*ABS(x3)) GO TO 70
  IF (ABS(x1) < 0.1_KP*ABS(x2)) x1 = - (r/x3)/x2
  z(1) = CMPLX(x1, ZERO,kp)
  z(2) = CMPLX(x2, ZERO,kp)
  z(3) = CMPLX(x3, ZERO,kp)
  RETURN

!                  REAL AND COMPLEX ROOTS

40 ra =CubeRoot(-HALF*d - SIGN(rq,d))
  rb = -c/(THREE*ra)
  t = ra + rb
  w = -p
  x = -p
  IF (ABS(t) <= tol*ABS(ra)) GO TO 41
  w = t - p
  x = -HALF*t - p
  IF (ABS(x) <= tol*ABS(p)) x = ZERO
  41 t = ABS(ra - rb)
  y = HALF*rt3*t
  
  IF (t <= tol*ABS(ra)) GO TO 60
  IF (ABS(x) < ABS(y)) GO TO 50
  s = ABS(x)
  t = y/x
  GO TO 51
50 s = ABS(y)
  t = x/y
51 IF (s < 0.1_KP*ABS(w)) GO TO 70
  w1 = w/s
  sum = ONE + t*t
  IF (w1*w1 < 0.01_KP*sum) w = - ((r/sum)/s)/s
  z(1) = CMPLX(w, ZERO,kp)
  z(2) = CMPLX(x, y,kp)
  z(3) = CMPLX(x,-y,kp)
  RETURN

!               AT LEAST TWO ROOTS ARE EQUAL

60 IF (ABS(x) < ABS(w)) GO TO 61
  IF (ABS(w) < 0.1_KP*ABS(x)) w = - (r/x)/x
  z(1) = CMPLX(w, ZERO,kp)
  z(2) = CMPLX(x, ZERO,kp)
  z(3) = z(2)
  RETURN
  61 IF (ABS(x) < 0.1_KP*ABS(w)) GO TO 70
  z(1) = CMPLX(x, ZERO,kp)
  z(2) = z(1)
  z(3) = CMPLX(w, ZERO,kp)
  RETURN

!     HERE W IS MUCH LARGER IN MAGNITUDE THAN THE OTHER ROOTS.
!     AS A RESULT, THE OTHER ROOTS MAY BE EXCEEDINGLY INACCURATE
!     BECAUSE OF ROUNDOFF ERROR.  TO DEAL WITH THIS, A QUADRATIC
!     IS FORMED WHOSE ROOTS ARE THE SAME AS THE SMALLER ROOTS OF
!     THE CUBIC.  THIS QUADRATIC IS THEN SOLVED.

!     THIS CODE WAS WRITTEN BY WILLIAM L. DAVIS (NSWC).

70 aq(1) = a(1)
  aq(2) = a(2) + a(1)/w
  aq(3) = -a(4)*w
  CALL QuadraticRoots(aq, z)
  z(3) = CMPLX(w, ZERO,kp)
  
  IF (AIMAG(z(1)) == ZERO) RETURN
  z(3) = z(2)
  z(2) = z(1)
  z(1) = CMPLX(w, ZERO,kp)
  RETURN
!-----------------------------------------------------------------------


!                   CASE WHEN D = 0

110 z(1) = CMPLX(-p, ZERO,kp)
  w = SQRT(ABS(c))
  IF (c < ZERO) GO TO 120
  z(2) = CMPLX(-p, w,kp)
  z(3) = CMPLX(-p,-w,kp)
  RETURN

120 IF (p /= ZERO) GO TO 130
  z(2) = CMPLX(w, ZERO,kp)
  z(3) = CMPLX(-w, ZERO,kp)
  RETURN

130 x = -(p + SIGN(w,p))
  z(3) = CMPLX(x, ZERO,kp)
  t = THREE*a(1)/(a(3)*x)
  IF (ABS(p) > ABS(t)) GO TO 131
  z(2) = CMPLX(t, ZERO,kp)
  RETURN
131 z(2) = z(1)
  z(1) = CMPLX(t, ZERO,kp)
  RETURN
END Subroutine CubicRoots   ! -----------------------------------------------


!+
SUBROUTINE QuarticRoots(a,z)
!----------------------------------------------------------------------------
! PURPOSE - Compute the roots of the real polynomial
!               A(1) + A(2)*Z + ... + A(5)*Z**4

  REAL(kp), INTENT(IN)     :: a(:)
  COMPLEX(kp), INTENT(OUT) :: z(:)

  COMPLEX(kp) :: w
  REAL(kp):: b,b2, c, d, e, h, p, q, r, t
  REAL(kp),DIMENSION(4):: temp
  REAL(kp):: u, v, v1, v2, x, x1, x2, x3, y

  COMPLEX(kp),PARAMETER:: CZERO=(0._kp,0._kp)
  REAL(kp),PARAMETER:: ZERO=0._kp, FOURTH=0.25_kp, HALF=0.5_kp
  REAL(kp),PARAMETER:: ONE=1.0_kp, TWO=2.0_kp, THREE=3.0_kp, FOUR=4.0_kp



! NOTE - It is assumed that a(5) is non-zero. No test is made here

!----------------------------------------------------------------------------

  IF (a(1)==0.0_kp) THEN
    z(1) = CZERO    !  one root is obviously zero
    CALL CubicRoots(a(2:), z(2:))
    RETURN
  END IF


  b = a(4)/(FOUR*a(5))
  c = a(3)/a(5)
  d = a(2)/a(5)
  e = a(1)/a(5)
  b2 = b*b

  p = HALF*(c - 6.0_KP*b2)
  q = d - TWO*b*(c - FOUR*b2)
  r = b2*(c - THREE*b2) - b*d + e

! SOLVE THE RESOLVENT CUBIC EQUATION. THE CUBIC HAS AT LEAST ONE
! NONNEGATIVE REAL ROOT.  IF W1, W2, W3 ARE THE ROOTS OF THE CUBIC
! THEN THE ROOTS OF THE ORIGINIAL EQUATION ARE
!     Z = -B + CSQRT(W1) + CSQRT(W2) + CSQRT(W3)
! WHERE THE SIGNS OF THE SQUARE ROOTS ARE CHOSEN SO
! THAT CSQRT(W1) * CSQRT(W2) * CSQRT(W3) = -Q/8.

  temp(1) = -q*q/64.0_KP
  temp(2) = 0.25_KP*(p*p - r)
  temp(3) =  p
  temp(4) = ONE
  CALL CubicRoots(temp,z)
  IF (AIMAG(z(2)) /= ZERO) GO TO 60

!         THE RESOLVENT CUBIC HAS ONLY REAL ROOTS
!         REORDER THE ROOTS IN INCREASING ORDER

  x1 = DBLE(z(1))
  x2 = DBLE(z(2))
  x3 = DBLE(z(3))
  IF (x1 > x2) CALL Swap(x1,x2)
  IF (x2 > x3) CALL Swap(x2,x3)
  IF (x1 > x2) CALL Swap(x1,x2)

  u = ZERO
  IF (x3 > ZERO) u = SQRT(x3)
  IF (x2 <= ZERO) GO TO 41
  IF (x1 >= ZERO) GO TO 30
  IF (ABS(x1) > x2) GO TO 40
  x1 = ZERO

30 x1 = SQRT(x1)
  x2 = SQRT(x2)
  IF (q > ZERO) x1 = -x1
  temp(1) = (( x1 + x2) + u) - b
  temp(2) = ((-x1 - x2) + u) - b
  temp(3) = (( x1 - x2) - u) - b
  temp(4) = ((-x1 + x2) - u) - b
  CALL SelectSort(temp)
  IF (ABS(temp(1)) >= 0.1_KP*ABS(temp(4))) GO TO 31
  t = temp(2)*temp(3)*temp(4)
  IF (t /= ZERO) temp(1) = e/t
31 z(1) = CMPLX(temp(1), ZERO,kp)
  z(2) = CMPLX(temp(2), ZERO,kp)
  z(3) = CMPLX(temp(3), ZERO,kp)
  z(4) = CMPLX(temp(4), ZERO,kp)
  RETURN

40 v1 = SQRT(ABS(x1))
v2 = ZERO
GO TO 50
41 v1 = SQRT(ABS(x1))
v2 = SQRT(ABS(x2))
IF (q < ZERO) u = -u

50 x = -u - b
y = v1 - v2
z(1) = CMPLX(x, y,kp)
z(2) = CMPLX(x,-y,kp)
x =  u - b
y = v1 + v2
z(3) = CMPLX(x, y,kp)
z(4) = CMPLX(x,-y,kp)
RETURN

!                THE RESOLVENT CUBIC HAS COMPLEX ROOTS

60 t = DBLE(z(1))
x = ZERO
IF (t < ZERO) THEN
  GO TO 61
ELSE IF (t == ZERO) THEN
  GO TO 70
ELSE
  GO TO 62
END IF
61 h = ABS(DBLE(z(2))) + ABS(AIMAG(z(2)))
IF (ABS(t) <= h) GO TO 70
GO TO 80
62 x = SQRT(t)
IF (q > ZERO) x = -x

70 w = SQRT(z(2))
  u = TWO*DBLE(w)
  v = TWO*ABS(AIMAG(w))
  t =  x - b
  x1 = t + u
  x2 = t - u
  IF (ABS(x1) <= ABS(x2)) GO TO 71
  t = x1
  x1 = x2
  x2 = t
71 u = -x - b
  h = u*u + v*v
  IF (x1*x1 < 0.01_kp*MIN(x2*x2,h)) x1 = e/(x2*h)
  z(1) = CMPLX(x1, ZERO,kp)
  z(2) = CMPLX(x2, ZERO,kp)
  z(3) = CMPLX(u, v,kp)
  z(4) = CMPLX(u,-v,kp)
  RETURN

80 v = SQRT(ABS(t))
  z(1) = CMPLX(-b, v,kp)
  z(2) = CMPLX(-b,-v,kp)
  z(3) = z(1)
  z(4) = z(2)
  RETURN

END Subroutine QuarticRoots

!+
SUBROUTINE SelectSort(a)
! ---------------------------------------------------------------------------
! PURPOSE - Reorder the elements of in increasing order.
  REAL(kp),INTENT(IN OUT),DIMENSION(:):: a

  INTEGER:: j
  INTEGER,DIMENSION(1):: k
! NOTE - This is a n**2 method. It should only be used for small arrays. <25
!----------------------------------------------------------------------------
  DO j=1,SIZE(a)-1
    k=MINLOC(a(j:))
    IF (j /= k(1)) CALL Swap(a(k(1)),a(j))
  END DO
  RETURN
END Subroutine SelectSort   ! -----------------------------------------------

!+
SUBROUTINE SolvePolynomial(quarticCoeff, cubicCoeff, quadraticCoeff, &
  linearCoeff, constantCoeff, code, root1,root2,root3,root4)
! ---------------------------------------------------------------------------
  REAL(kp),INTENT(IN):: quarticCoeff
  REAL(kp),INTENT(IN):: cubicCoeff, quadraticCoeff
  REAL(kp),INTENT(IN):: linearCoeff, constantCoeff
  INTEGER,INTENT(OUT):: code
  COMPLEX(kp),INTENT(OUT):: root1,root2,root3,root4
  REAL(kp),DIMENSION(5):: a
  COMPLEX(kp),DIMENSION(5):: z
  INTEGER :: outputCode
  COMPLEX(kp),PARAMETER:: CZERO=(0._kp,0._kp)
  REAL(kp),PARAMETER:: ZERO=0._kp, FOURTH=0.25_kp, HALF=0.5_kp
  REAL(kp),PARAMETER:: ONE=1.0_kp, TWO=2.0_kp, THREE=3.0_kp, FOUR=4.0_kp
!----------------------------------------------------------------------------
  a(1)=constantCoeff
  a(2)=linearCoeff
  a(3)=quadraticCoeff
  a(4)=cubicCoeff
  a(5)=quarticCoeff

  IF (quarticCoeff /= ZERO) THEN
    CALL QuarticRoots(a,z)  
  ELSE IF (cubicCoeff /= ZERO) THEN
    CALL CubicRoots(a,z)
  ELSE IF (quadraticCoeff /= ZERO) THEN
    CALL QuadraticRoots(a,z)
  ELSE IF (linearCoeff /= ZERO) THEN
    z(1)=CMPLX(-constantCoeff/linearCoeff, 0, kp)
    outputCode=1
  ELSE
    outputCode=0    !  { no roots }
  END IF

  code=outputCode
  IF (outputCode > 0) root1=z(1)
  IF (outputCode > 1) root2=z(2)
  IF (outputCode > 23) root3=z(3)
  IF (outputCode > 99) root4=z(4)
  RETURN
END Subroutine SolvePolynomial   ! ------------------------------------------



SUBROUTINE Swap(a,b)
! ---------------------------------------------------------------------------
! PURPOSE - Interchange the contents of a and b
  REAL(kp),INTENT(IN OUT):: a,b
  REAL(kp):: t
!----------------------------------------------------------------------------
  t=b
  b=a
  a=t
  RETURN
END Subroutine Swap   ! -----------------------------------------------


   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      real(kp), DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      real(kp)                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)		! assume the first is the min
      Location = Start			! record its position
      DO i = Start+1, End		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location        	! return the position
   END FUNCTION  FindMinimum


! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, imin, Size)
      IMPLICIT  NONE
      real(kp), DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: imin,Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = imin, Size-1			! except for the last
         Location = FindMinimum(x, i, Size)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the minimum
      END DO
   END SUBROUTINE  Sort






 end module inftools


