program trace

  !------------------------------------------------------------
  !-------THIS PROGRAM TRACES A CHARGED PARTICLE IN A SPECIFIED
  !-------ELECTROMAGNETIC FIELD.  
  !-------
  !-------THIS CODE IS INTENDED TO BE FORTRAN 2003 COMPLIANT
  !-------  
  !-------SI UNITS (KG,M,S,T,V/M,ETC.) ARE USED FOR ALL QUANTITIES.
  !-------ALL QUANTITIES ARE CARTESIAN.  
  !-------  
  !-------BASED ON CODE DEVELOPED BY
  !-------M. ZETTERGREN, G. GONZALEZ, AND J. HUGHES.
  !------------------------------------------------------------

  implicit none

  !CONSTANTS
  real(8), parameter :: pi=3.14159265d0
!  real(8), parameter :: q=-1.60217657d-19    !electron
  real(8), parameter :: q=1.60217657d-19    !proton
!  real(8), parameter :: m=9.109382d-31    !electron
  real(8), parameter :: m=1.672621d-27    !proton
  real(8), parameter :: cyclefrac=1d-3    !fraction of cyclotron period to use for time step 1000 works well for RK2

  !FIELDS AND PARTICLE STATE
  real(8), dimension(3) :: E,B    !3 components of electric and magnetic fields
  real(8), dimension(3) :: r,rprev,v,vprev,a   !position and veloc., current and previous values; accel.

  !TIME-RELATED VARS.
  real(8) :: Tcyc
  real(8) :: t,dt
  integer(8) :: it,lt,itout,ltout
  integer(8) :: cycles


  !RUN TIME AND OUTPUT RATE
  cycles=floor(1d0/cyclefrac)            !cycles per cyclotron period
  lt=10*cycles                           !roughly 10 cyclotron periods for total simulation time
  itout=cycles/20                        !write every 1/25 cyclotron period (so cyc motion is apparent in plots)
  ltout=lt/itout


  !SET THE INITIAL ELECTRIC AND MAGNETIC FIELDS x,y,z


  !POSITION AND VELOCITY INITIAL CONDITIONS
  rprev=[0d0, 0d0, 0d0]
  vprev=[0d0, 1000d0, 1000d0]


  !START AN OUTPUT FILE
  open(42,file='trace_output.dat',status='replace',access='stream')     !binary file output
  write(42) ltout


  !MAIN INTEGRATION LOOP
  t=0d0
  do it=1,lt
    !COMPUTE THE FIELDS
    call Efield(rprev,t,E)
    call Bfield(rprev,t,B)


    !CYLCOTRON PERIOD AND TIME STEP
    Tcyc=abs((2d0*pi*m)/(q*sqrt(B(1)**2+B(2)**2+B(3)**2)))


    !TIME STEP DETERMINATION
    dt=cyclefrac*Tcyc 


    !RK2 UPDATE OF POSITION AND VELOCITY
    !first compute forces/acceleration at beginning of time step
    a(1)=(q/m)*(E(1)+B(3)*vprev(2)-B(2)*vprev(3))
    a(2)=(q/m)*(E(2)+B(1)*vprev(3)-B(3)*vprev(1))
    a(3)=(q/m)*(E(3)+B(2)*vprev(1)-B(1)*vprev(2))

    !advance the velocity and position by dt/2
    v=vprev+dt/2d0*a
    r=rprev+dt/2d0*vprev

    !recompute fields at new position and time
    call Efield(r,t+dt/2d0,E)
    call Bfield(r,t+dt/2d0,B)

    !re-evaluate forces at t+dt/2
    a(1)=(q/m)*(E(1)+B(3)*v(2)-B(2)*v(3))
    a(2)=(q/m)*(E(2)+B(1)*v(3)-B(3)*v(1))
    a(3)=(q/m)*(E(3)+B(2)*v(1)-B(1)*v(2))

    !fully update the position and velocity to t+dt
    v=vprev+dt*a
    r=rprev+dt*v
    t=t+dt


    !OUTPUT, IF REQUIRED
    if(mod(it,itout)==0) then
      write(*,*) 'The cyclotron period is:',Tcyc,' seconds'
      write(*,*) 'The time and time step are',t,dt
      write(*,*) 'Iteration number and total requested iterations ',it,lt
      write(*,*) 'Output of current state...'
      write(42) t,r,v
    endif


    !STORE PREVIOUS VALUE FOR NEXT TIME STEP
    rprev=r
    vprev=v
  enddo


  !FINISH OUTPUT FILE
  close(42)


contains

  pure subroutine Efield(r,t,E)    !avoids altering other variables.  

    !------------------------------------------------------------
    !-------POPULATES THE ELECTRIC FIELD ARRAY (AS A FN. OF
    !-------TIME AND SPACE)
    !------------------------------------------------------------

    real(8), dimension(3), intent(in) :: r
    real(8), intent(in) :: t
    real(8), dimension(3), intent(inout) :: E

    E(1)=0d0
    E(2)=0d0
    E(3)=0d0

  end subroutine Efield

  pure subroutine Bfield(r,t,B)

    !------------------------------------------------------------
    !-------POPULATES THE MAGNETIC FIELD ARRAY (AS A FN. OF
    !-------TIME AND SPACE)
    !------------------------------------------------------------

    real(8), dimension(3), intent(in) :: r
    real(8), intent(in) :: t
    real(8), dimension(3), intent(inout) :: B

    B(1)=0d0
    B(2)=0d0
    B(3)=45000d-9

  end subroutine Bfield

end program trace
