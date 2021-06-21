!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! timerMod: simple timer module
!
! This module implements a simple CPU timer using Fortran timing 
! intrinsics. 
!
!***********************************************************************

module timerMod

   implicit none
   private

   ! basic timer type
   type Timer
      real :: startTime   ! start time of current interval
      real :: endTime     ! end   time of current interval
      real :: elapsedTime ! total elapsed time for this timer
      logical :: running  ! true if timer has been started
      character (len=120) :: name ! name for this timer
   end type Timer

   public:: &
      timerCreate, &
      timerStart,  &
      timerStop,   &
      timerPrint

!***********************************************************************

contains

!***********************************************************************
! timerCreate
!
! This subroutine creates a timer and initializes the members. The 
! timer is given a name based on user input.
!
!-----------------------------------------------------------------------

   function timerCreate(timerName) result(newTimer)

      type(Timer), intent(out) :: newTimer
      character (*), intent(in) :: timerName ! name to assign this timer

      newTimer%name = ' '
      newTimer%name = trim(timerName)

      newTimer%startTime   = 0.0
      newTimer%endTime     = 0.0
      newTimer%elapsedTime = 0.0
      newTimer%running  = .false.

   end function timerCreate

!***********************************************************************
! timerStart
!
! This subroutine starts a timer.
!
!-----------------------------------------------------------------------

   subroutine timerStart(inTimer)

      type (Timer), intent(inout) :: inTimer ! timer to start

      ! Check if timer has already been started; if so, stop the timer
      if (inTimer%running) call timerStop(inTimer)

      ! Start the timer and call the Fortran intrinsic to get start time
      inTimer%running = .true.
      call cpu_time(inTimer%startTime)

   end subroutine timerStart

!***********************************************************************
! timerStop
!
! This subroutine stops a timer and accumulates total elapsed time.
!
!-----------------------------------------------------------------------

   subroutine timerStop(inTimer)

      type (Timer), intent(inout) :: inTimer ! timer to stop

      ! To minimize overhead, call stop time right away
      call cpu_time(inTimer%stopTime)

      ! Accumulate the elapsed time if the timer was running and reset
      ! If timer wasn't running, ignore the result
      if (inTimer%running) then
         inTimer%running = .false.
         inTimer%elapsedTime = inTimer%elapsedTime &
                             + (inTimer%stopTime - inTimer%startTime)
      endif

   end subroutine timerStop

!***********************************************************************
! timerPrint
!
! This subroutine prints the total elapsed time for a given timer.
!
!-----------------------------------------------------------------------

   subroutine timerPrint(inTimer)

      type (Timer), intent(inout) :: inTimer ! timer to print out

      ! If timer is still running, stop the timer
      if (inTimer%running) call timerStop(inTimer)

      ! Print the elapsed time from this timer
      print *,'Time in timer: ',trim(timerName),' = ', &
              inTimer%elapsedTime, ' seconds.'

   end subroutine timerPrint

!***********************************************************************

end module timerMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
