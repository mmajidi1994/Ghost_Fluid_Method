program euler_eqs

use m_parameters
use m_startup
use m_bcs
use m_dt
use m_weno
use m_flux
use m_qpqs
use m_riemann_solvers
use leftright
use m_rhs
use m_rk

! main loop

call s_initial_condition
call s_rhs


do k=1,20000!int(nstep)
 
   if (k<50) then 
   cn=0.1d0
   else 
   cn=0.9d0
   end if
  if (mod(k,500)==0) then
  print *, k, dt, time
  end if
  call time_step 
 call s_rk
  p_time = time 
  time = time + dt
  if (time > r_time) then
   dt = r_time -p_time
  call s_rk
   time = p_time + dt
   print *, time    

  write(11,*)" r, rho, u, p, e "
  do i=1,n
  write(11,*) rr(i), rho(i), u(i), p(i), p(i)/((gama-1.0d0)*rho(i))
  end do

  stop
  
  end if 
 
end do


end program

