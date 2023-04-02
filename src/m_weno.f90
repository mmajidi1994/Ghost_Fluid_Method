MODULE m_weno! weno reconstruction on primitive varibles 

use m_parameters
use m_bcs

IMPLICIT NONE

PRIVATE; PUBLIC :: s_weno3, s_weno5

CONTAINS
subroutine s_weno3()


do j=1,3
do i=0,n+1

 qp0(j,i)= 0.5*qp(j,i)+0.5*qp(j,i+1)
 qp1(j,i)= -0.5*qp(j,i-1)+1.5*qp(j,i)

end do 
end do 

do i=0,n+1

b1(0,i)=(qp(1,i+1)-qp(1,i))**2.0
b1(1,i)=(qp(1,i)-qp(1,i-1))**2.0
b2(0,i)=(qp(2,i+1)-qp(2,i))**2.0
b2(1,i)=(qp(2,i)-qp(2,i-1))**2.0
b3(0,i)=(qp(3,i+1)-qp(3,i))**2.0
b3(1,i)=(qp(3,i)-qp(3,i-1))**2.0

end do 
do j=0,1
do i= 0,n+1
   a1(j,i)=d3(j)/((eps + b1(j,i))**2.0)
   a2(j,i)=d3(j)/((eps + b2(j,i))**2.0)
   a3(j,i)=d3(j)/((eps + b3(j,i))**2.0)
  end do 
end do 

do j=0,1
do i=0,n+1
   w1(j,i)=a1(j,i)/(a1(0,i)+a1(1,i))
   w2(j,i)=a2(j,i)/(a2(0,i)+a2(1,i))
   w3(j,i)=a3(j,i)/(a3(0,i)+a3(1,i))
  end do
end do 

! (i+1/2)-

do i=0,n+1

qpm(1,i)=w1(0,i)*qp0(1,i)+w1(1,i)*qp1(1,i)
qpm(2,i)=w2(0,i)*qp0(2,i)+w2(1,i)*qp1(2,i)
qpm(3,i)=w3(0,i)*qp0(3,i)+w3(1,i)*qp1(3,i)

end do 

do j=1,3
do i=0,n+1

 qp0(j,i)= 1.5*qp(j,i)-0.5*qp(j,i+1)
 qp1(j,i)= 0.5*qp(j,i-1)+0.5*qp(j,i)

end do 
end do 

do j=0,1
do i= 0,n+1
   a1(j,i)=d3(1-j)/((eps + b1(j,i))**2.0)
   a2(j,i)=d3(1-j)/((eps + b2(j,i))**2.0)
   a3(j,i)=d3(1-j)/((eps + b3(j,i))**2.0)
end do
end do 

do j=0,1
do i= 0,n+1
   w1(j,i)=a1(j,i)/(a1(0,i)+a1(1,i))
   w2(j,i)=a2(j,i)/(a2(0,i)+a2(1,i))
   w3(j,i)=a3(j,i)/(a3(0,i)+a3(1,i))
  end do 
end do 

! (i-1/2)+

do i=0,n+1
  qpp(1,i)=w1(0,i)*qp0(1,i)+w1(1,i)*qp1(1,i)
  qpp(2,i)=w2(0,i)*qp0(2,i)+w2(1,i)*qp1(2,i)
  qpp(3,i)=w3(0,i)*qp0(3,i)+w3(1,i)*qp1(3,i)
end do 

end subroutine


subroutine s_weno5()


do j=1,3
do i=0,n+1

 qs0(j,i)= c0(0)*qp(j,i)+c0(1)*qp(j,i+1)+c0(2)*qp(j,i+2)
 qs1(j,i)= c1(0)*qp(j,i-1)+c1(1)*qp(j,i)+c1(2)*qp(j,i+1)
 qs2(j,i)= c2(0)*qp(j,i-2)+c2(1)*qp(j,i-1)+c2(2)*qp(j,i)

end do 
end do 

do i=0,n+1

b1(0,i)=(13.0/12.0)*(qp(1,i)-2.0*qp(1,i+1)+qp(1,i+2))**2.0&
+(1.0/4.0)*(3.0*qp(1,i)-4.0*qp(1,i+1)+qp(1,i+2))**2.0

b1(1,i)=(13.0/12.0)*(qp(1,i-1)-2.0*qp(1,i)+qp(1,i+1))**2.0&
+(1.0/4.0)*(qp(1,i-1)-qp(1,i+1))**2.0

b1(2,i)=(13.0/12.0)*(qp(1,i-2)-2.0*qp(1,i-1)+qp(1,i))**2.0&
+(1.0/4.0)*(qp(1,i-2)-4.0*qp(1,i-1)+3.0*qp(1,i))**2.0

b2(0,i)=(13.0/12.0)*(qp(2,i)-2.0*qp(2,i+1)+qp(2,i+2))**2.0&
+(1.0/4.0)*(3.0*qp(2,i)-4.0*qp(2,i+1)+qp(2,i+2))**2.0

b2(1,i)=(13.0/12.0)*(qp(2,i-1)-2.0*qp(2,i)+qp(2,i+1))**2.0&
+(1.0/4.0)*(qp(2,i-1)-qp(2,i+1))**2.0

b2(2,i)=(13.0/12.0)*(qp(2,i-2)-2.0*qp(2,i-1)+qp(2,i))**2.0&
+(1.0/4.0)*(qp(2,i-2)-4.0*qp(2,i-1)+3.0*qp(2,i))**2.0

b3(0,i)=(13.0/12.0)*(qp(3,i)-2.0*qp(3,i+1)+qp(3,i+2))**2.0&
+(1.0/4.0)*(3.0*qp(3,i)-4.0*qp(3,i+1)+qp(3,i+2))**2.0

b3(1,i)=(13.0/12.0)*(qp(3,i-1)-2.0*qp(3,i)+qp(3,i+1))**2.0&
+(1.0/4.0)*(qp(3,i-1)-qp(3,i+1))**2.0

b3(2,i)=(13.0/12.0)*(qp(3,i-2)-2.0*qp(3,i-1)+qp(3,i))**2.0&
+(1.0/4.0)*(qp(3,i-2)-4.0*qp(3,i-1)+3.0*qp(3,i))**2.0

end do 

do j=0,2
do i= 0,n+1
   a1(j,i)=d(j)/((eps + b1(j,i))**2.0)
   a2(j,i)=d(j)/((eps + b2(j,i))**2.0)
   a3(j,i)=d(j)/((eps + b3(j,i))**2.0)
  end do 
end do 
do j=0,2
do i=0,n+1
   w1(j,i)=a1(j,i)/(a1(0,i)+a1(1,i)+a1(2,i))
   w2(j,i)=a2(j,i)/(a2(0,i)+a2(1,i)+a2(2,i))
   w3(j,i)=a3(j,i)/(a3(0,i)+a3(1,i)+a3(2,i))
  end do
end do 
! (i+1/2)-
do i=0,n+1
qpm(1,i)=w1(0,i)*qs0(1,i)+w1(1,i)*qs1(1,i)+w1(2,i)*qs2(1,i)
qpm(2,i)=w2(0,i)*qs0(2,i)+w2(1,i)*qs1(2,i)+w2(2,i)*qs2(2,i)
qpm(3,i)=w3(0,i)*qs0(3,i)+w3(1,i)*qs1(3,i)+w3(2,i)*qs2(3,i)
end do 

do j=1,3
do i=0,n+1
 qs0(j,i)= ch0(0)*qp(j,i)+ch0(1)*qp(j,i+1)+ch0(2)*qp(j,i+2)
 qs1(j,i)= ch1(0)*qp(j,i-1)+ch1(1)*qp(j,i)+ch1(2)*qp(j,i+1)
 qs2(j,i)= ch2(0)*qp(j,i-2)+ch2(1)*qp(j,i-1)+ch2(2)*qp(j,i)
end do 
end do 

do j=0,2
do i= 0,n+1
   a1(j,i)=dh(j)/((eps + b1(j,i))**2.0)
   a2(j,i)=dh(j)/((eps + b2(j,i))**2.0)
   a3(j,i)=dh(j)/((eps + b3(j,i))**2.0)
end do
end do 

do j=0,2
do i= 0,n+1
   w1(j,i)=a1(j,i)/(a1(0,i)+a1(1,i)+a1(2,i))
   w2(j,i)=a2(j,i)/(a2(0,i)+a2(1,i)+a2(2,i))
   w3(j,i)=a3(j,i)/(a3(0,i)+a3(1,i)+a3(2,i))
  end do 
end do
! (i-1/2)+
do i=0,n+1
  qpp(1,i)=w1(0,i)*qs0(1,i)+w1(1,i)*qs1(1,i)+w1(2,i)*qs2(1,i)
  qpp(2,i)=w2(0,i)*qs0(2,i)+w2(1,i)*qs1(2,i)+w2(2,i)*qs2(2,i)
  qpp(3,i)=w3(0,i)*qs0(3,i)+w3(1,i)*qs1(3,i)+w3(2,i)*qs2(3,i)
end do 


end subroutine
END MODULE

