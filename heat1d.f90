!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Forward time Central Space (FTCS) scheme to solve heat equation 
!     1st order Euler scheme for time
!     2nd order central scheme for space
! >>> Crank-Nicolson (CN) scheme to solve heat equation 
!     2nd order trapezoidal for time
!     2nd order central scheme for space
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) 
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 8, 2015
!-----------------------------------------------------------------------------!

program heat1d
implicit none
integer::i,k,nx,nt,na,nb,nc,nd
real*8 ::dx,dt,x0,xL,pi,r,t,alpha,Tmax,ta,tb,tc,td
real*8,allocatable ::u(:),x(:)

!Domain
x0 = 0.0d0 !left
xL = 1.0d0 !right

!number of points
nx = 20

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

!maximum time desired
Tmax = 2.0d0

!time step
dt = 0.001d0

!number of points in time
nt = nint(Tmax/dt)

!ploting times:
ta = 0.5d0
tb = 1.0d0
tc = 1.5d0
td = 2.0d0

na = nint(ta/dt)
nb = nint(tb/dt)
nc = nint(tc/dt)
nd = nint(td/dt)

!diffusion coefficient
alpha = 1.0d0

!u: temperature variable 
allocate(u(0:nx))

!initial condition
pi = 4.0d0*datan(1.0d0)
t = 0.0d0
do i=0,nx
u(i) = dsin(pi*x(i))
end do

!boundary conditions: holds for all time
u(0) = 0.0d0
u(nx)= 0.0d0

!Plot initial condition
open(18,file='u.plt')
write(18,*) 'variables ="x","u"'
write(18,100)'zone f=point i=',nx+1,',t="t =',t,'"'
do i=0,nx
write(18,*) x(i),u(i)
end do



r = alpha*dt/(dx*dx)

print*,r
pause

!time integration
do k=1,nt
     
!update with FTCS scheme
    call FTCS(nx,dx,dt,t,alpha,x,u)
    
!update with CN scheme
    !call CN(nx,dx,dt,t,alpha,x,u)
    
    !update t
    t = t+dt 
 


    !plot field
    if (k.eq.na.or.k.eq.nb.or.k.eq.nc.or.k.eq.nd) then
	write(18,100)'zone f=point i=',nx+1,',t="t =',t,'"'
	do i=0,nx
	write(18,*) x(i),u(i)
	end do
    end if

    
	print*,k,t,maxval(u)


end do
  

close(18)

100 format(a16,i8,a10,f8.4,a3)
end


!-----------------------------------------------------------------------------!
!FTCS scheme
!-----------------------------------------------------------------------------!
subroutine FTCS(nx,dx,dt,t,alpha,x,u)
implicit none
integer::nx,i
real*8 ::dx,dt,t,r,alpha,f
real*8 ::u(0:nx),x(0:nx)

r = alpha*dt/(dx*dx)

do i=1,nx-1
u(i) = u(i) + r*(u(i+1)-2.0d0*u(i)+u(i-1)) + dt*f(t,x(i)) 
end do

end 


!-----------------------------------------------------------------------------!
!Crank-Nicolson(CN) scheme
!-----------------------------------------------------------------------------!
subroutine CN(nx,dx,dt,t,alpha,x,u)
implicit none
integer::nx,i
real*8 ::dx,dt,t,alpha,beta,f
real*8 ::u(0:nx),x(0:nx)
real*8,allocatable ::a(:),b(:),c(:),r(:),q(:)

beta = 0.5d0*alpha*dt/(dx*dx)

!Build coefficient matrix:
allocate(a(1:nx-1),b(1:nx-1),c(1:nx-1),r(1:nx-1),q(1:nx-1))

do i=1,nx-1
a(i) = -beta
b(i) = (1.0d0+2.0d0*beta)
c(i) = -beta
r(i) = beta*u(i+1) + (1.0d0-2.0d0*beta)*u(i)+beta*u(i-1) + 0.5d0*dt*(f(t,x(i))+f(t+dt,x(i)))
end do
!apply boundary conditions
r(1)   = r(1) - a(1)*u(0)     !b.c.
r(nx-1) = r(nx-1) - c(nx-1)*u(nx) !b.c.
call tdma(a,b,c,r,q,1,nx-1)

!assign solutions to y
do i=1,nx-1
u(i)=q(i)
end do

end 

!------------------------------------------------------------------!
!Tridiagonal matrix algorithm (TDMA)
!Thomas algorithm
!solution tridiagonal systems
!a: lower diagonal
!b: main diagonal
!c: upper diagonal
!r: source vector
!x: solution vector
!   for indices s(start) to e(end)
!   i: s,s+1,s+2, ....,i,....,e 
!
!Note: a(s) and c(e) are dummy coefficients, not used.
!------------------------------------------------------------------!

subroutine tdma(a,b,c,r,x,s,e)
implicit none
integer s,e,i
real*8, dimension(s:e) ::a,b,c,r,x    

! forward elimination phase
do i=s+1,e
b(i) = b(i) - a(i)/b(i-1)*c(i-1)
r(i) = r(i) - a(i)/b(i-1)*r(i-1)
end do
! backward substitution phase 
x(e) = r(e)/b(e)
do i=e-1,s,-1
x(i) = (r(i)-c(i)*x(i+1))/b(i)
end do

return
end


!-----------------------------------------------------------------------------!
!source term (nonhomogenous forcing term in heat equation)
!-----------------------------------------------------------------------------!
real*8 function f(t,x)
implicit none
real*8::t,x,pi
pi = 4.0d0*datan(1.0d0)
f = (pi*pi - 1.0d0)*dexp(-t)*dsin(pi*x)
end











