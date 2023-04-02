MODULE m_parameters

IMPLICIT NONE 

PRIVATE; PUBLIC :: s_parameters

INTEGER :: i,j,k

INTEGER, PARAMETER :: n = 100, dp=3*n/10;

REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: rr, rho, u, p, e;
REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: rol, ror, ul, ur, pl, pr;
REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: al, ar, sl, sr, ss, plr,
sn;
REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:) :: c0,c1,c2,ch0,ch1,ch2,d,dh;

REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: qp0, qp1, qp2;
REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: qp, qpp, qpm, qs, qsp,
qsm;
REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: fl, fr, flx, ds, fls,
frs;

REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: b1, b2, b3, a1, a2, a3;
REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: w1, w2, w3;
REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: qsold, qsnew, qsnew1,
qsnew2;
REAL(KIND(0d0)), ALLOCATABLE, DIMENSION(:,:) :: qsnew3, qsnew4, qsnew5,
qsnew6;

REAL(KIND(0d0)), ALLOCATABLE:: pi, r_max, r_min, interval, delta_r, dr,
gp;
REAL(KIND(0d0)), ALLOCATABLE:: eps, r_time, p_time, time, cn, dt, nstep,
t_r, gama, ep, alpha, snmax;

CONTAINS

  SUBROUTINE s_parameters()
 
 eps=10.0d-15
 pi=3.141592653
 r_max=1.0d0
 r_min=0.0d0
 interval=r_max-r_min
 delta_r= interval/n
 dr=delta_r
 time=0.0d0
 ep=1.0
 gama=1.4
 alpha=dr/dt
 gp = 1.0d0 ! Geometrical Parameter 
 r_time=0.001d0

 c0(0:2) = [1.0/3.0, 5.0/6.0, -1.0/6.0] 
 c1(0:2) = [-1.0/6.0, 5.0/6.0, 1.0/3.0]  
 c2(0:2) = [1.0/3.0, -7.0/6.0, 11.0/6.0]   
 ch0(0:2) = [11.0/6.0, -7.0/6.0, 1.0/3.0]   
 ch1(0:2) = [1.0/3.0, 5.0/6.0, -1.0/6.0]  
 ch2(0:2) = [-1.0/6.0, 5.0/6.0, 1.0/3.0]  
 d(0:2) = [3.0/10.0, 3.0/5.0, 1.0/10.0]   
 dh(0:2) = [1.0/10.0, 3.0/5.0, 3.0/10.0] 

     ALLOCATE(rr(-4:n+4))
     ALLOCATE(rho(-4:n+4))
     ALLOCATE(u(-4:n+4))
     ALLOCATE(p(-4:n+4))
     ALLOCATE(e(-4:n+4))
     ALLOCATE(rol(-4:n+4))
     ALLOCATE(ror(-4:n+4))
     allocate(ul(-4:n+4))
     allocate(ur(-4:n+4))
     ALLOCATE(pl(-4:n+4))
     ALLOCATE(pr(-4:n+4))
     ALLOCATE(al(-4:n+4))
     ALLOCATE(ar(-4:n+4))
     ALLOCATE(sl(-4:n+4))
     ALLOCATE(sr(-4:n+4))
     ALLOCATE(ss(-4:n+4))
     ALLOCATE(plr(-4:n+4))
     ALLOCATE(sn(-4:n+4))

     ALLOCATE(qp0(1:3,-4:n+4))
     ALLOCATE(qp1(1:3,-4:n+4))
     ALLOCATE(qp2(1:3,-4:n+4))
     ALLOCATE(qp(1:3,-4:n+4))
     ALLOCATE(qpp(1:3,-4:n+4))
     ALLOCATE(qpm(1:3,-4:n+4))
     ALLOCATE(qs(1:3,-4:n+4))
     ALLOCATE(qsp(1:3,-4:n+4))
     ALLOCATE(qsm(1:3,-4:n+4))
     ALLOCATE(fl(1:3,-4:n+4))
     ALLOCATE(fr(1:3,-4:n+4))
     ALLOCATE(flx(1:3,-4:n+4))
     ALLOCATE(ds(1:3,-4:n+4))
     ALLOCATE(fls(1:3,-4:n+4))
     ALLOCATE(frs(1:3,-4:n+4))
    
     ALLOCATE(b1(0:2,-4:n+4))
     ALLOCATE(b2(0:2,-4:n+4))
     ALLOCATE(b3(0:2,-4:n+4))
     ALLOCATE(a1(0:2,-4:n+4))
     ALLOCATE(a2(0:2,-4:n+4))
     ALLOCATE(a3(0:2,-4:n+4))
     ALLOCATE(w1(0:2,-4:n+4))
     ALLOCATE(w2(0:2,-4:n+4))
     ALLOCATE(w3(0:2,-4:n+4))

END SUBROUTINE

END MODULE 

