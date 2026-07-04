!flux 2d edge flux evaluation module 
!max wood
!version : 0.0.6
!updated : 04-07-26

!module 
module edge_flux
use flux_data_methods
contains 

!VanLeer flux ===============
subroutine vanleer_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen,fx1,fx2,fx3,fx4)
implicit none 

!variables - inout 
real(dp) :: rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen
real(dp) :: fx1,fx2,fx3,fx4

!variables - local
real(dp) ::  ma,c1,c2,vn1,vn2,m1p,m2p,m1,m2
real(dp) :: f11,f21,f31,f41,f12,f22,f32,f42

!get edge normal velocities
vn1 = u1*nx + v1*ny 
vn2 = u2*nx + v2*ny 

!get the speed of sound in each cell 
c1 = speed_of_sound(p1,rho1,gamma)
c2 = speed_of_sound(p2,rho2,gamma)

!get the machs in each cell
m1 = vn1/c1
m2 = vn2/c2
if (abs(m1) .LE. 1.0d0) then 
    m1p = 0.25d0*((m1 + 1.0d0)**2)
else !m1 .LE. -1.0d0
    m1p = 0.5d0*(m1 + abs(m1))
end if 
if (abs(m2) .LE. 1.0d0) then 
    m2p = -0.25d0*((m2 - 1.0d0)**2)
else !m1 .LE. -1.0d0
    m2p = 0.5d0*(m2 - abs(m2))
end if 
ma = m1p + m2p

!evaluate the edge flux
if (ma .GE. 1.0d0) then 
    fx1 = (rho1*vn1)*elen
    fx2 = (rho1*u1*vn1 + nx*p1)*elen
    fx3 = (rho1*v1*vn1 + ny*p1)*elen
    fx4 = (rho1*vn1*e1 + vn1*p1)*elen
elseif (ma .LE. -1.0d0) then 
    fx1 = (rho2*vn2)*elen
    fx2 = (rho2*u2*vn2 + nx*p2)*elen
    fx3 = (rho2*v2*vn2 + ny*p2)*elen
    fx4 = (rho2*vn2*e2 + vn2*p2)*elen
else

    !evaluate flux 1
    f11 = 0.25d0*rho1*c1*(m1 + 1.0d0)**2
    f21 = f11*((nx*(-vn1 + 2.0d0*c1)/gamma) + u1)
    f31 = f11*((ny*(-vn1 + 2.0d0*c1)/gamma) + v1)
    f41 = f11*((((gamma - 1.0d0)*vn1 + 2.0d0*c1)**2)/(2.0d0*(gamma*gamma - 1.0d0)) + 0.5d0*(u1*u1 + v1*v1 - vn1*vn1))

    !evaluate flux 2
    f12 = -0.25d0*rho2*c2*(m2 - 1.0d0)**2
    f22 = f12*((nx*(-vn2 - 2.0d0*c2)/gamma) + u2)
    f32 = f12*((ny*(-vn2 - 2.0d0*c2)/gamma) + v2)
    f42 = f12*((((gamma - 1.0d0)*vn2 - 2.0d0*c2)**2)/(2.0d0*(gamma*gamma - 1.0d0)) + 0.5d0*(u2*u2 + v2*v2 - vn2*vn2))

    !evaluate the edge flux
    fx1 = (f11 + f12)*elen
    fx2 = (f21 + f22)*elen
    fx3 = (f31 + f32)*elen
    fx4 = (f41 + f42)*elen
end if 
return 
end subroutine vanleer_flux

!VanLeer flux (complex) ===============
subroutine vanleer_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen,fx1,fx2,fx3,fx4)
implicit none 

!variables - inout 
complex(dp) :: rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen
complex(dp) :: fx1,fx2,fx3,fx4

!variables - local
complex(dp) ::  ma,c1,c2,vn1,vn2,m1p,m2p,m1,m2
complex(dp) :: f11,f21,f31,f41,f12,f22,f32,f42

!get edge normal velocities
vn1 = u1*nx + v1*ny 
vn2 = u2*nx + v2*ny 

!get the speed of sound in each cell 
c1 = speed_of_sound_cpx(p1,rho1,gamma)
c2 = speed_of_sound_cpx(p2,rho2,gamma)

!get the machs in each cell
m1 = vn1/c1
m2 = vn2/c2
if (real(m1,dp) .LE. 1.0d0) then 
    m1p = 0.25d0*((m1 + 1.0d0)**2)
else !m1 .LE. -1.0d0
    m1p = 0.5d0*(m1 + sqrt(m1*m1))
end if 
if (real(m2,dp) .LE. 1.0d0) then 
    m2p = -0.25d0*((m2 - 1.0d0)**2)
else !m1 .LE. -1.0d0
    m2p = 0.5d0*(m2 - sqrt(m2*m2))
end if 
ma = m1p + m2p

!evaluate the edge flux
if (real(ma,dp) .GE. 1.0d0) then 
    fx1 = (rho1*vn1)*elen
    fx2 = (rho1*u1*vn1 + nx*p1)*elen
    fx3 = (rho1*v1*vn1 + ny*p1)*elen
    fx4 = (rho1*vn1*e1 + vn1*p1)*elen
elseif (real(ma,dp) .LE. -1.0d0) then 
    fx1 = (rho2*vn2)*elen
    fx2 = (rho2*u2*vn2 + nx*p2)*elen
    fx3 = (rho2*v2*vn2 + ny*p2)*elen
    fx4 = (rho2*vn2*e2 + vn2*p2)*elen
else

    !evaluate flux 1
    f11 = 0.25d0*rho1*c1*(m1 + 1.0d0)**2
    f21 = f11*((nx*(-vn1 + 2.0d0*c1)/gamma) + u1)
    f31 = f11*((ny*(-vn1 + 2.0d0*c1)/gamma) + v1)
    f41 = f11*((((gamma - 1.0d0)*vn1 + 2.0d0*c1)**2)/(2.0d0*(gamma*gamma - 1.0d0)) + 0.5d0*(u1*u1 + v1*v1 - vn1*vn1))

    !evaluate flux 2
    f12 = -0.25d0*rho2*c2*(m2 - 1.0d0)**2
    f22 = f12*((nx*(-vn2 - 2.0d0*c2)/gamma) + u2)
    f32 = f12*((ny*(-vn2 - 2.0d0*c2)/gamma) + v2)
    f42 = f12*((((gamma - 1.0d0)*vn2 - 2.0d0*c2)**2)/(2.0d0*(gamma*gamma - 1.0d0)) + 0.5d0*(u2*u2 + v2*v2 - vn2*vn2))

    !evaluate the edge flux
    fx1 = (f11 + f12)*elen
    fx2 = (f21 + f22)*elen
    fx3 = (f31 + f32)*elen
    fx4 = (f41 + f42)*elen
end if 
return 
end subroutine vanleer_flux_cpx

!AUSM flux ===============
subroutine ausm_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen,fx1,fx2,fx3,fx4)
implicit none 

!variables - inout 
real(dp) :: rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen
real(dp) :: fx1,fx2,fx3,fx4

!variables - local
real(dp) ::  ma,c1,c2,p1p,p2p,vn1,vn2,m1p,m2p,m1,m2,pb

!get edge normal velocities
vn1 = u1*nx + v1*ny 
vn2 = u2*nx + v2*ny 

!get the speed of sound in each cell 
c1 = speed_of_sound(p1,rho1,gamma)
c2 = speed_of_sound(p2,rho2,gamma)

!get the machs in each cell
m1 = vn1/c1
m2 = vn2/c2
if (abs(m1) .LE. 1.0d0) then 
    m1p = 0.25d0*((m1 + 1.0d0)**2)
else !m1 .LE. -1.0d0
    m1p = 0.5d0*(m1 + abs(m1))
end if 
if (abs(m2) .LE. 1.0d0) then 
    m2p = -0.25d0*((m2 - 1.0d0)**2)
else !m1 .LE. -1.0d0
    m2p = 0.5d0*(m2 - abs(m2))
end if 
ma = m1p + m2p

!evaluate the boundary pressure
if (abs(m1) .LE. 1.0d0) then 
    p1p = 0.25d0*p1*(m1 + 1.0)*(m1 + 1.0)*(2.0 - m1)
else
    p1p = 0.5d0*p1*(m1 + abs(m1))/m1
end if 
if (abs(m2) .LE. 1.0d0) then 
    p2p = 0.25d0*p2*(m2 - 1.0)*(m2 - 1.0)*(2.0 + m2)
else
    p2p = 0.5d0*p2*(m2 - abs(m2))/m2
end if 
pb = p1p + p2p 

!evaluate the flux
if (ma .GE. 0.0d0) then 
    fx1 = ma*rho1*c1*elen
    fx2 = (ma*rho1*c1*u1 + nx*pb)*elen
    fx3 = (ma*rho1*c1*v1 + ny*pb)*elen
    fx4 = ma*c1*(rho1*e1 + p1)*elen
else
    fx1 = ma*rho2*c2*elen
    fx2 = (ma*rho2*c2*u2 + nx*pb)*elen
    fx3 = (ma*rho2*c2*v2 + ny*pb)*elen
    fx4 = ma*c2*(rho2*e2 + p2)*elen
end if 
return 
end subroutine ausm_flux

!AUSM flux (complex) ===============
subroutine ausm_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen,fx1,fx2,fx3,fx4)
implicit none 

!variables - inout 
complex(dp) :: rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen
complex(dp) :: fx1,fx2,fx3,fx4

!variables - local
complex(dp) ::  ma,c1,c2,p1p,p2p,vn1,vn2,m1p,m2p,m1,m2,pb

!get edge normal velocities
vn1 = u1*nx + v1*ny 
vn2 = u2*nx + v2*ny 

!get the speed of sound in each cell 
c1 = speed_of_sound_cpx(p1,rho1,gamma)
c2 = speed_of_sound_cpx(p2,rho2,gamma)

!get the machs in each cell
m1 = vn1/c1
m2 = vn2/c2
if (abs(real(m1,dp)) .LE. 1.0d0) then 
    m1p = 0.25d0*((m1 + 1.0d0)**2)
else !m1 .LE. -1.0d0
    m1p = 0.5d0*(m1 + sqrt(m1*m1))
end if 
if (abs(real(m2,dp)) .LE. 1.0d0) then 
    m2p = -0.25d0*((m2 - 1.0d0)**2)
else !m1 .LE. -1.0d0
    m2p = 0.5d0*(m2 - sqrt(m2*m2))
end if 
ma = m1p + m2p

!evaluate the boundary pressure
if (abs(real(m1,dp)) .LE. 1.0d0) then 
    p1p = 0.25d0*p1*(m1 + 1.0)*(m1 + 1.0)*(2.0 - m1)
else
    p1p = 0.5d0*p1*(m1 + sqrt(m1*m1))/m1
end if 
if (abs(real(m2,dp)) .LE. 1.0d0) then 
    p2p = 0.25d0*p2*(m2 - 1.0)*(m2 - 1.0)*(2.0 + m2)
else
    p2p = 0.5d0*p2*(m2 - sqrt(m2*m2))/m2
end if 
pb = p1p + p2p 

!evaluate the flux
if (real(ma,dp) .GE. 0.0d0) then 
    fx1 = ma*rho1*c1*elen
    fx2 = (ma*rho1*c1*u1 + nx*pb)*elen
    fx3 = (ma*rho1*c1*v1 + ny*pb)*elen
    fx4 = ma*c1*(rho1*e1 + p1)*elen
else
    fx1 = ma*rho2*c2*elen
    fx2 = (ma*rho2*c2*u2 + nx*pb)*elen
    fx3 = (ma*rho2*c2*v2 + ny*pb)*elen
    fx4 = ma*c2*(rho2*e2 + p2)*elen
end if 
return 
end subroutine ausm_flux_cpx

!AUSM+up flux ===============
subroutine ausmpup_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen,fx1,fx2,fx3,fx4,machinf,R)
implicit none 

!variables - inout 
real(dp) :: rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen,machinf,R
real(dp) :: fx1,fx2,fx3,fx4

!variables - local
real(dp) :: vn1,vn2,m1,m2,mbar2,mo2,mo,fa,rho12,t1,t2,c12,cs1,cs2,ct1,ct2
real(dp) :: m41,m42,m12,mdot12,p51,p52,p12

!get edge normal velocities
vn1 = u1*nx + v1*ny 
vn2 = u2*nx + v2*ny 

!get the speed of sound in each cell and the boundary
t1 = p1/(rho1*R)
cs1 = (((gamma*R*t1)/(gamma - 1.0d0)) + 0.5d0*(u1*u1 + v1*v1))*(2.0*(gamma - 1.0d0))/(gamma + 1.0)
ct1 = cs1/max(sqrt(cs1),vn1)
t2 = p2/(rho2*R)
cs2 = (((gamma*R*t2)/(gamma - 1.0d0)) + 0.5d0*(u2*u2 + v2*v2))*(2.0*(gamma - 1.0d0))/(gamma + 1.0)
ct2 = cs2/max(sqrt(cs2),-vn2)
c12 = min(ct1,ct2)

!get the machs in each cell
m1 = vn1/c12
m2 = vn2/c12

!get the mach blending functions 
mbar2 = (vn1*vn1 + vn2*vn2)/(2.0d0*c12*c12)
mo2 = min(1.0d0,max(mbar2,machinf*machinf))
mo = sqrt(mo2)
fa = mo*(2.0d0 - mo)
m41 = ausmpup_mach_blend(m1,1.0d0)
m42 = ausmpup_mach_blend(m2,-1.0d0)
rho12 = 0.5d0*(rho1 + rho2)

!get the blended mach 
m12 = m41 + m42 - (0.25d0/fa)*max(1.0d0 - mbar2,0.0d0)*((p2 - p1)/(rho12*c12*c12))

!get the mass flux 
if (m12 .GT. 0.0d0) then 
    mdot12 = m12*c12*rho1
else
    mdot12 = m12*c12*rho2
end if 

!get the pressure flux 
p51 = ausmpup_pressure_blend(m1,1.0d0,fa)
p52 = ausmpup_pressure_blend(m2,-1.0d0,fa)
p12 = p51*p1 + p52*p2 - 0.75d0*p51*p52*(rho1 + rho2)*fa*c12*(vn2 - vn1)

!evaluate the flux
if (mdot12 .GT. 0.0d0) then 
    fx1 = mdot12*elen 
    fx2 = (mdot12*u1 + nx*p12)*elen
    fx3 = (mdot12*v1 + ny*p12)*elen
    fx4 = m12*c12*(rho1*e1 + p1)*elen
else
    fx1 = mdot12*elen 
    fx2 = (mdot12*u2 + nx*p12)*elen
    fx3 = (mdot12*v2 + ny*p12)*elen
    fx4 = m12*c12*(rho2*e2 + p2)*elen
end if 
return 
end subroutine ausmpup_flux

!AUSM+up mach function ===============
function ausmpup_mach_blend(m,pm) result(m4)
implicit none 

!variables - inout 
real(dp) :: m,pm,m4

!variables - local
real(dp) :: m1,m2,m2m

!get m1 and m2
m1 = 0.5d0*(m + pm*abs(m))
m2 = pm*0.25d0*((m + pm)**2)

!get m4
if (abs(m) .GE. 1.0d0) then 
    m4 = m1
else
    m2m = -pm*0.25d0*((m - pm)**2)
    m4 = m2*(1.0d0 - pm*2.0d0*m2m)
end if 
return 
end function ausmpup_mach_blend

!AUSM+up pressure function ===============
function ausmpup_pressure_blend(m,pm,fa) result(p5)
implicit none 

!variables - inout 
real(dp) :: m,pm,p5,fa

!variables - local
real(dp) :: m1,m2,m2m,a

!get m1 and m2
m1 = 0.5d0*(m + pm*abs(m))
m2 = pm*0.25d0*((m + pm)**2)

!get p5
if (abs(m) .GE. 1.0d0) then 
    p5 = (1.0d0/m)*m1
else
    a = (3.0d0/16.0d0)*(-4.0d0 + 5.0d0*fa*fa)
    m2m = -pm*0.25d0*((m - pm)**2)
    p5 = m2*((pm*2.0d0 - m) - pm*16.0d0*a*m*m2m)
end if 
return 
end function ausmpup_pressure_blend

!AUSM+up flux (complex) ===============
subroutine ausmpup_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen,fx1,fx2,fx3,fx4,machinf,R)
implicit none 

!variables - inout 
complex(dp) :: rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen,machinf,R
complex(dp) :: fx1,fx2,fx3,fx4

!variables - local
complex(dp) :: vn1,vn2,m1,m2,mbar2,mo2,mo,fa,rho12,t1,t2,c12,cs1,cs2,ct1,ct2
complex(dp) :: m41,m42,m12,mdot12,p51,p52,p12
complex(dp) :: dnm

!get edge normal velocities
vn1 = u1*nx + v1*ny 
vn2 = u2*nx + v2*ny 

!get the speed of sound in each cell and the boundary
t1 = p1/(rho1*R)
cs1 = (((gamma*R*t1)/(gamma - 1.0d0)) + 0.5d0*(u1*u1 + v1*v1))*(2.0*(gamma - 1.0d0))/(gamma + 1.0)
if (real(sqrt(cs1),dp) .GT. real(vn1,dp)) then 
    ct1 = cs1/sqrt(cs1)
else
    ct1 = cs1/vn1
end if 
! ct1 = cs1/max(sqrt(cs1),vn1)
t2 = p2/(rho2*R)
cs2 = (((gamma*R*t2)/(gamma - 1.0d0)) + 0.5d0*(u2*u2 + v2*v2))*(2.0*(gamma - 1.0d0))/(gamma + 1.0)
if (real(sqrt(cs2),dp) .GT. real(-vn2,dp)) then 
    ct2 = cs2/sqrt(cs2)
else
    ct2 = cs2/(-vn2)
end if 
! ct2 = cs2/max(sqrt(cs2),-vn2)
if (real(ct1,dp) .LT. real(ct2,dp)) then 
    c12 = ct1 
else
    c12 = ct2
end if 
! c12 = min(ct1,ct2)

!get the machs in each cell
m1 = vn1/c12
m2 = vn2/c12

!get the mach blending functions 
mbar2 = (vn1*vn1 + vn2*vn2)/(2.0d0*c12*c12)

if (real(mbar2,dp) .GT. real(machinf*machinf,dp)) then 
    mo2 = mbar2
else
    mo2 = machinf*machinf
end if 
if (1.0d0 .LT. real(mo2,dp)) then 
    mo2 = complex(1.0d0,0.0d0) 
end if 
! mo2 = min(1.0d0,max(mbar2,machinf*machinf))
mo = sqrt(mo2)
fa = mo*(2.0d0 - mo)
m41 = ausmpup_mach_blend_cpx(m1,complex(1.0d0,0.0d0))
m42 = ausmpup_mach_blend_cpx(m2,complex(-1.0d0,0.0d0))
rho12 = 0.5d0*(rho1 + rho2)

!get the blended mach 
if (real(1.0d0 - mbar2) .GT. 0.0d0) then
    dnm = 1.0d0 - mbar2
else
    dnm = 0.0
end if 
m12 = m41 + m42 - (0.25d0/fa)*dnm*((p2 - p1)/(rho12*c12*c12))
! m12 = m41 + m42 - (0.25d0/fa)*max(1.0d0 - mbar2,0.0d0)*((p2 - p1)/(rho12*c12*c12))

!get the mass flux 
if (real(m12,dp) .GT. 0.0d0) then 
    mdot12 = m12*c12*rho1
else
    mdot12 = m12*c12*rho2
end if 

!get the pressure flux 
p51 = ausmpup_pressure_blend_cpx(m1,complex(1.0d0,0.0d0),fa)
p52 = ausmpup_pressure_blend_cpx(m2,complex(-1.0d0,0.0d0),fa)
p12 = p51*p1 + p52*p2 - 0.75d0*p51*p52*(rho1 + rho2)*fa*c12*(vn2 - vn1)

!evaluate the flux
if (real(mdot12,dp) .GT. 0.0d0) then 
    fx1 = mdot12*elen 
    fx2 = (mdot12*u1 + nx*p12)*elen
    fx3 = (mdot12*v1 + ny*p12)*elen
    fx4 = m12*c12*(rho1*e1 + p1)*elen
else
    fx1 = mdot12*elen 
    fx2 = (mdot12*u2 + nx*p12)*elen
    fx3 = (mdot12*v2 + ny*p12)*elen
    fx4 = m12*c12*(rho2*e2 + p2)*elen
end if 
return 
end subroutine ausmpup_flux_cpx

!AUSM+up mach function (complex) ===============
function ausmpup_mach_blend_cpx(m,pm) result(m4)
implicit none 

!variables - inout 
complex(dp) :: m,pm,m4

!variables - local
complex(dp) :: m1,m2,m2m

!get m1 and m2
m1 = 0.5d0*(m + pm*abs(m))
m2 = pm*0.25d0*((m + pm)**2)

!get m4
if (abs(m) .GE. 1.0d0) then 
    m4 = m1
else
    m2m = -pm*0.25d0*((m - pm)**2)
    m4 = m2*(1.0d0 - pm*2.0d0*m2m)
end if 
return 
end function ausmpup_mach_blend_cpx

!AUSM+up pressure function (complex) ===============
function ausmpup_pressure_blend_cpx(m,pm,fa) result(p5)
implicit none 

!variables - inout 
complex(dp) :: m,pm,p5,fa

!variables - local
complex(dp) :: m1,m2,m2m,a

!get m1 and m2
m1 = 0.5d0*(m + pm*abs(m))
m2 = pm*0.25d0*((m + pm)**2)

!get p5
if (abs(m) .GE. 1.0d0) then 
    p5 = (1.0d0/m)*m1
else
    a = (3.0d0/16.0d0)*(-4.0d0 + 5.0d0*fa*fa)
    m2m = -pm*0.25d0*((m - pm)**2)
    p5 = m2*((pm*2.0d0 - m) - pm*16.0d0*a*m*m2m)
end if 
return 
end function ausmpup_pressure_blend_cpx

!Roe flux ===============
subroutine roe_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen,fx1,fx2,fx3,fx4)
implicit none 

!variables - inout 
real(dp) :: rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,gamma,nx,ny,elen
real(dp) :: fx1,fx2,fx3,fx4

!variables - local
real(dp) :: srho1,srho2,srhos,h1,h2,drho,dpr,du,dv,dvel,dlam
real(dp) :: rhor,ur,vr,hr,cr,velr,q2r,w1,w2a,w2b,w3,vn
real(dp) :: f1l,f2l,f3l,f4l,f1r,f2r,f3r,f4r
real(dp) :: lamda1,lamda2,lamda3
real(dp) :: r1,r2,r3,r4

!get the enthalpy in each cell 
h1 = e1 + (p1/rho1)
h2 = e2 + (p2/rho2)

!get the cell deltas
du = u2 - u1 
dv = v2 - v1 
dpr = p2 - p1 
drho = rho2 - rho1 
dvel = (u2*nx + v2*ny) - (u1*nx + v1*ny)

!evaluate the Roe average properties
srho1 = sqrt(rho1)
srho2 = sqrt(rho2)
srhos = srho1 + srho2
rhor = sqrt(rho1*rho2)
ur = (srho1*u1 + srho2*u2)/srhos
vr = (srho1*v1 + srho2*v2)/srhos
hr = (srho1*h1 + srho2*h2)/srhos
q2r = ur*ur + vr*vr
cr = sqrt((gamma - 1.0d0)*(hr - 0.5*q2r))
velr = ur*nx + vr*ny 

!evaluate the Roe weightings 
dlam = 0.01*(speed_of_sound(p1,rho1,gamma) + speed_of_sound(p2,rho2,gamma))*0.5
lamda1 = roe_eigen_fix(abs(velr - cr),dlam)
lamda2 = roe_eigen_fix(abs(velr),dlam)
lamda3 = roe_eigen_fix(abs(velr + cr),dlam)
w1 = abs(velr - cr)*((dpr - rhor*cr*dvel)/(2.0d0*cr))
w2a = abs(velr)*(drho - (dpr/(cr*cr)))
w2b = abs(velr)*rhor
w3 = abs(velr + cr)*((dpr + rhor*cr*dvel)/(2.0d0*cr))

!evaluate the Roe offsets 
r1 = w1 + w2a + w3 
r2 = w1*(ur - cr*nx) + w2a*ur + w2b*(du - dvel*nx) + w3*(ur + cr*nx) 
r3 = w1*(vr - cr*ny) + w2a*vr + w2b*(dv - dvel*ny) + w3*(vr + cr*ny) 
r4 = w1*(hr - cr*velr) + w2a*q2r*0.5d0 + w2b*(ur*du + vr*dv - velr*dvel) + w3*(hr + cr*velr)

!evaluate the flux left
vn = u1*nx + v1*ny
f1l = rho1*vn 
f2l = rho1*u1*vn + nx*p1
f3l = rho1*v1*vn + ny*p1
f4l = rho1*vn*e1 + vn*p1

!evaluate the flux right
vn = u2*nx + v2*ny
f1r = rho2*vn 
f2r = rho2*u2*vn + nx*p2
f3r = rho2*v2*vn + ny*p2
f4r = rho2*vn*e2 + vn*p2

!evaluate the Roe flux 
fx1 = 0.5d0*(f1l + f1r - r1)*elen
fx2 = 0.5d0*(f2l + f2r - r2)*elen
fx3 = 0.5d0*(f3l + f3r - r3)*elen
fx4 = 0.5d0*(f4l + f4r - r4)*elen
return 
end subroutine roe_flux

!Roe eigenvalue fix function ===============
function roe_eigen_fix(lamdain,delta) result(lamdaout)
implicit none 

!variables - inout 
real(dp) :: delta
real(dp) :: lamdain,lamdaout

!evaluate
if (lamdain > delta) then 
    lamdaout = lamdain
else
    lamdaout = (lamdain*lamdain + delta*delta)/(2.0d0*delta)
end if
return 
end function roe_eigen_fix

!Jameson flux ===============
subroutine jameson_flux(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,nx,ny,elen,fx1,fx2,fx3,fx4)
implicit none 

!variables - inout 
real(dp) :: rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,nx,ny,elen
real(dp) :: fx1,fx2,fx3,fx4

!variables - local 
real(dp) :: vn
real(dp) :: f1l,f2l,f3l,f4l
real(dp) :: f1r,f2r,f3r,f4r

!flux left
vn = u1*nx + v1*ny 
f1l = rho1*vn 
f2l = rho1*u1*vn + nx*p1
f3l = rho1*v1*vn + ny*p1
f4l = rho1*vn*e1 + vn*p1

!flux right
vn = u2*nx + v2*ny 
f1r = rho2*vn 
f2r = rho2*u2*vn + nx*p2
f3r = rho2*v2*vn + ny*p2
f4r = rho2*vn*e2 + vn*p2

!edge flux
fx1 = 0.5d0*(f1l + f1r)*elen
fx2 = 0.5d0*(f2l + f2r)*elen
fx3 = 0.5d0*(f3l + f3r)*elen
fx4 = 0.5d0*(f4l + f4r)*elen
return 
end subroutine jameson_flux

!Jameson flux (complex) ===============
subroutine jameson_flux_cpx(rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,nx,ny,elen,fx1,fx2,fx3,fx4)
implicit none 

!variables - inout 
complex(dp) :: rho1,u1,v1,e1,p1,rho2,u2,v2,e2,p2,nx,ny,elen
complex(dp) :: fx1,fx2,fx3,fx4

!variables - local 
complex(dp) :: vn
complex(dp) :: f1l,f2l,f3l,f4l
complex(dp) :: f1r,f2r,f3r,f4r

!flux left
vn = u1*nx + v1*ny 
f1l = rho1*vn 
f2l = rho1*u1*vn + nx*p1
f3l = rho1*v1*vn + ny*p1
f4l = rho1*vn*e1 + vn*p1

!flux right
vn = u2*nx + v2*ny 
f1r = rho2*vn 
f2r = rho2*u2*vn + nx*p2
f3r = rho2*v2*vn + ny*p2
f4r = rho2*vn*e2 + vn*p2

!edge flux
fx1 = complex(0.5d0,0.0d0)*(f1l + f1r)*elen
fx2 = complex(0.5d0,0.0d0)*(f2l + f2r)*elen
fx3 = complex(0.5d0,0.0d0)*(f3l + f3r)*elen
fx4 = complex(0.5d0,0.0d0)*(f4l + f4r)*elen
return 
end subroutine jameson_flux_cpx

!farfield prescribed boundary condition (supersonic) ===============
subroutine farfield_supersonic_bc_prescribed(rhob,ub,vb,pb,eb,mesh,c,options,inflow_outflow)
implicit none 

!variables - inout 
integer(in64) :: c,inflow_outflow
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: ub,vb,rhob,pb,eb

!get boundary state
if (inflow_outflow == 1) then !inflow
    ub = options%uinf
    vb = options%vinf
    rhob = options%rhoinf
    pb = options%pinf
elseif (inflow_outflow == -1) then !outflow
    ub = mesh%u(c)
    vb = mesh%v(c)
    rhob = mesh%rho(c)
    pb = mesh%p(c)
end if 
eb = energy(pb,rhob,sqrt(ub*ub + vb*vb),options%gamma)
return 
end subroutine farfield_supersonic_bc_prescribed

!farfield prescribed boundary condition (supersonic) (complex) ===============
subroutine farfield_supersonic_bc_prescribed_cpx(rhob,ub,vb,pb,eb,mesh,c,options,inflow_outflow)
implicit none 

!variables - inout 
integer(in64) :: c,inflow_outflow
type(flux_mesh_cpx) :: mesh 
type(flux_options) :: options 

!variables - local 
complex(dp) :: ub,vb,rhob,pb,eb

!get boundary state
if (inflow_outflow == 1) then !inflow
    ub = complex(options%uinf,0.0d0)
    vb = complex(options%vinf,0.0d0)
    rhob = complex(options%rhoinf,0.0d0)
    pb = complex(options%pinf,0.0d0)
elseif (inflow_outflow == -1) then !outflow
    ub = mesh%u(c)
    vb = mesh%v(c)
    rhob = mesh%rho(c)
    pb = mesh%p(c)
end if 
eb = energy_cpx(pb,rhob,sqrt(ub*ub + vb*vb),complex(options%gamma,0.0d0))
return 
end subroutine farfield_supersonic_bc_prescribed_cpx

!farfield characteristic boundary condition (subsonic) ===============
subroutine farfield_subsonic_bc_characteristic(rhob,ub,vb,pb,eb,mesh,c,e,options,inflow_outflow)
implicit none 

!variables - inout 
integer(in64) :: e,c,inflow_outflow
real(dp) :: rhob,ub,vb,pb,eb
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: a,du,dv,drho,dpr,c1,c2,c3,c4
real(dp) :: vel_in(2),vel_inf(2),vel_in_n(2),vel_inf_n(2),vel_b(2),vel_b_n(2)
real(dp) :: basis_bx(2),basis_by(2)
real(dp) :: Mb2a(2,2),Ma2b(2,2) 

!get the local speed of sound 
a = speed_of_sound(mesh%p(c),mesh%rho(c),options%gamma)

!get the local basis coordinate system wrt to the edge normal 
if (inflow_outflow == 1) then !inflow
    basis_bx(1) = -mesh%edges(e)%nx !nx
    basis_bx(2) = -mesh%edges(e)%ny !ny
elseif (inflow_outflow == -1) then !outflow
    basis_bx(1) = mesh%edges(e)%nx !nx
    basis_bx(2) = mesh%edges(e)%ny !ny
end if 
basis_by(1) = mesh%edges(e)%dx/mesh%edges(e)%length
basis_by(2) = mesh%edges(e)%dy/mesh%edges(e)%length
call get_basis_change_2d(Ma2b,Mb2a,basis_bx,basis_by)

!translate the velocity to this basis 
vel_in(1) = mesh%u(c)
vel_in(2) = mesh%v(c)
vel_inf(1) = options%uinf
vel_inf(2) = options%vinf 
vel_in_n = change_basis(Ma2b,vel_in)
vel_inf_n = change_basis(Ma2b,vel_inf)

!get the characteristic delta variables 
du = vel_in_n(1) - vel_inf_n(1)
dv = vel_in_n(2) - vel_inf_n(2)
drho = mesh%rho(c) - options%rhoinf 
dpr = mesh%p(c) - options%pinf 
c1 = -a*a*drho + dpr 
c2 = mesh%rho(c)*a*dv 
c3 = mesh%rho(c)*a*du + dpr 
c4 = -mesh%rho(c)*a*du + dpr  

!apply the boundary condition 
if (inflow_outflow == 1) then !inflow 
    c1 = 0.0d0 
    c2 = 0.0d0 
    c3 = 0.0d0 
elseif (inflow_outflow == -1) then !outflow 
    c4 = 0.0d0 
end if 

!get the boundary deltas 
drho = (-1.0d0/(a*a))*c1 + (1.0d0/(2.0d0*a*a))*c3 + (1.0d0/(2.0d0*a*a))*c4 
du = (1.0d0/(2.0d0*mesh%rho(c)*a))*(c3 - c4)
dv = (1.0d0/(mesh%rho(c)*a))*c2 
dpr = 0.5d0*(c3 + c4)

!evaluate the boundary velocity in the base coordinate system 
vel_b_n(1) = vel_inf_n(1) + du 
vel_b_n(2) = vel_inf_n(2) + dv
vel_b = change_basis(Mb2a,vel_b_n)

!evaluate the boundary state
ub = vel_b(1)
vb = vel_b(2)
rhob = drho + options%rhoinf
pb = dpr + options%pinf
eb = energy(pb,rhob,sqrt(ub*ub + vb*vb),options%gamma)
return 
end subroutine farfield_subsonic_bc_characteristic

!farfield characteristic boundary condition (subsonic) (complex) ===============
subroutine farfield_subsonic_bc_characteristic_cpx(rhob,ub,vb,pb,eb,mesh,c,e,options,inflow_outflow)
implicit none 

!variables - inout 
integer(in64) :: e,c,inflow_outflow
complex(dp) :: rhob,ub,vb,pb,eb
type(flux_mesh_cpx) :: mesh 
type(flux_options) :: options 

!variables - local 
complex(dp) :: a,du,dv,drho,dpr,c1,c2,c3,c4
complex(dp) :: vel_in(2),vel_inf(2),vel_in_n(2),vel_inf_n(2),vel_b(2),vel_b_n(2)
complex(dp) :: basis_bx(2),basis_by(2)
complex(dp) :: Mb2a(2,2),Ma2b(2,2) 

!get the local speed of sound 
a = speed_of_sound_cpx(mesh%p(c),mesh%rho(c),complex(options%gamma,0.0d0))

!get the local basis coordinate system wrt to the edge normal 
if (inflow_outflow == 1) then !inflow
    basis_bx(1) = -mesh%edges(e)%nx !nx
    basis_bx(2) = -mesh%edges(e)%ny !ny
elseif (inflow_outflow == -1) then !outflow
    basis_bx(1) = mesh%edges(e)%nx !nx
    basis_bx(2) = mesh%edges(e)%ny !ny
end if 
basis_by(1) = mesh%edges(e)%dx/mesh%edges(e)%length
basis_by(2) = mesh%edges(e)%dy/mesh%edges(e)%length
call get_basis_change_2d_cpx(Ma2b,Mb2a,basis_bx,basis_by)

!translate the velocity to this basis 
vel_in(1) = mesh%u(c)
vel_in(2) = mesh%v(c)
vel_inf(1) = options%uinf
vel_inf(2) = options%vinf 
vel_in_n = change_basis_cpx(Ma2b,vel_in)
vel_inf_n = change_basis_cpx(Ma2b,vel_inf)

!get the characteristic delta variables 
du = vel_in_n(1) - vel_inf_n(1)
dv = vel_in_n(2) - vel_inf_n(2)
drho = mesh%rho(c) - options%rhoinf 
dpr = mesh%p(c) - options%pinf 
c1 = -a*a*drho + dpr 
c2 = mesh%rho(c)*a*dv 
c3 = mesh%rho(c)*a*du + dpr 
c4 = -mesh%rho(c)*a*du + dpr  

!apply the boundary condition 
if (inflow_outflow == 1) then !inflow 
    c1 = 0.0d0 
    c2 = 0.0d0 
    c3 = 0.0d0 
elseif (inflow_outflow == -1) then !outflow 
    c4 = 0.0d0 
end if 

!get the boundary deltas 
drho = (-1.0d0/(a*a))*c1 + (1.0d0/(2.0d0*a*a))*c3 + (1.0d0/(2.0d0*a*a))*c4 
du = (1.0d0/(2.0d0*mesh%rho(c)*a))*(c3 - c4)
dv = (1.0d0/(mesh%rho(c)*a))*c2 
dpr = 0.5d0*(c3 + c4)

!evaluate the boundary velocity in the base coordinate system 
vel_b_n(1) = vel_inf_n(1) + du 
vel_b_n(2) = vel_inf_n(2) + dv
vel_b = change_basis_cpx(Mb2a,vel_b_n)

!evaluate the boundary state
ub = vel_b(1)
vb = vel_b(2)
rhob = drho + options%rhoinf
pb = dpr + options%pinf
eb = energy_cpx(pb,rhob,sqrt(ub*ub + vb*vb),complex(options%gamma,0.0d0))
return 
end subroutine farfield_subsonic_bc_characteristic_cpx

!stagnation inflow ===============
subroutine subsonic_stagnation_inflow_bc(rhob,ub,vb,pb,eb,mesh,c,e,options)
implicit none 

!variables - inout 
integer(in64) :: e,c
real(dp) :: rhob,ub,vb,pb,eb
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: sosi,veli2,Hi,jm,alpha,aq,bq,cq
real(dp) :: sosb,sosb1,sosb2,velb,machb,tb

!get the internal enthalpy upstream running riemann invarient 
sosi = speed_of_sound(mesh%p(c),mesh%rho(c),options%gamma)
veli2 = mesh%u(c)**2 + mesh%v(c)**2
Hi = ((sosi*sosi)/(options%gamma - 1.0)) + 0.5d0*veli2 
jm = -sqrt(veli2) + ((2.0d0*sosi)/(options%gamma - 1.0))

!solve for the boundary speed of sound 
alpha = options%gamma - 1.0d0 
aq = (1.0d0/alpha) + (2.0d0/(alpha*alpha))
bq = -2.0d0*jm/alpha
cq = 0.5d0*jm*jm - Hi 
sosb1 = (-bq + sqrt(bq*bq - 4.0*aq*cq))/(2.0d0*aq)
sosb2 = (-bq - sqrt(bq*bq - 4.0*aq*cq))/(2.0d0*aq)
if (sosb1 .GE. sosb2) then 
    sosb = sosb1
else
    sosb = sosb2
end if 

!get the boundary flow state
velb = ((2.0d0*sosb)/(options%gamma - 1.0)) - jm 
machb = velb/sosb
pb = options%p0inf/((1.0d0 + 0.5d0*(options%gamma - 1.0d0)*machb*machb)**(options%gamma/(options%gamma - 1.0d0)))
tb = options%t0inf/(1.0d0 + 0.5d0*(options%gamma - 1.0d0)*machb*machb)
rhob = pb/(options%R*tb)
ub = -mesh%edges(e)%nx*velb
vb = -mesh%edges(e)%ny*velb
eb = energy(pb,rhob,sqrt(ub*ub + vb*vb),options%gamma)

!evaluate boundary flux
if (velb .GE. 0.0d0) then 
    if (machb .GE. 1.0) then 

        print *, 'supersonic stagnation inflow needs implementing '

        stop 
    end if 
else !revert to wall 
    rhob = 0.0
    ub = 0.0
    vb = 0.0
    pb = mesh%p(c)
    eb = 0.0
end if 
return 
end subroutine subsonic_stagnation_inflow_bc

!stagnation inflow (complex) ===============
subroutine subsonic_stagnation_inflow_bc_cpx(rhob,ub,vb,pb,eb,mesh,c,e,options)
implicit none 

!variables - inout 
integer(in64) :: e,c
complex(dp) :: rhob,ub,vb,pb,eb
type(flux_mesh_cpx) :: mesh 
type(flux_options) :: options 

!variables - local 
complex(dp) :: sosi,veli2,Hi,jm,alpha,aq,bq,cq
complex(dp) :: sosb,sosb1,sosb2,velb,machb,tb

!get the internal enthalpy upstream running riemann invarient 
sosi = speed_of_sound_cpx(mesh%p(c),mesh%rho(c),complex(options%gamma,0.0d0))
veli2 = mesh%u(c)**2 + mesh%v(c)**2
Hi = ((sosi*sosi)/(options%gamma - 1.0)) + 0.5d0*veli2 
jm = -sqrt(veli2) + ((2.0d0*sosi)/(options%gamma - 1.0))

!solve for the boundary speed of sound 
alpha = options%gamma - 1.0d0 
aq = (1.0d0/alpha) + (2.0d0/(alpha*alpha))
bq = -2.0d0*jm/alpha
cq = 0.5d0*jm*jm - Hi 
sosb1 = (-bq + sqrt(bq*bq - 4.0*aq*cq))/(2.0d0*aq)
sosb2 = (-bq - sqrt(bq*bq - 4.0*aq*cq))/(2.0d0*aq)
if (real(sosb1,dp) .GE. real(sosb2,dp)) then 
    sosb = sosb1
else
    sosb = sosb2
end if 

!get the boundary flow state
velb = ((2.0d0*sosb)/(options%gamma - 1.0)) - jm 
machb = velb/sosb
pb = options%p0inf/((1.0d0 + 0.5d0*(options%gamma - 1.0d0)*machb*machb)**(options%gamma/(options%gamma - 1.0d0)))
tb = options%t0inf/(1.0d0 + 0.5d0*(options%gamma - 1.0d0)*machb*machb)
rhob = pb/(options%R*tb)
ub = -mesh%edges(e)%nx*velb
vb = -mesh%edges(e)%ny*velb
eb = energy_cpx(pb,rhob,sqrt(ub*ub + vb*vb),complex(options%gamma,0.0d0))

!evaluate boundary flux
if (real(velb,dp) .GE. 0.0d0) then 
    if (real(machb,dp) .GE. 1.0) then 

        print *, 'cpx supersonic stagnation inflow needs implementing '

        stop 
    end if 
else !revert to wall 
    rhob = 0.0
    ub = 0.0
    vb = 0.0
    pb = mesh%p(c)
    eb = 0.0
end if 
return 
end subroutine subsonic_stagnation_inflow_bc_cpx

!pressure outflow boundary condition (subsonic) ===============
subroutine subsonic_pressure_outflow_bc_characteristic(rhob,ub,vb,pb,eb,mesh,c,e,options,pratio)
implicit none 

!variables - inout 
integer(in64) :: e,c
real(dp) :: rhob,ub,vb,pb,eb
real(dp) :: pratio
type(flux_mesh) :: mesh 
type(flux_options) :: options 

!variables - local 
real(dp) :: a,du,dv,drho,dpr,c1,c2,c3,c4
real(dp) :: vel_in(2),vel_in_n(2),vel_b(2),vel_b_n(2)
real(dp) :: basis_bx(2),basis_by(2)
real(dp) :: Mb2a(2,2),Ma2b(2,2) 

!get the local speed of sound 
a = speed_of_sound(mesh%p(c),mesh%rho(c),options%gamma)

!get the local basis coordinate system wrt to the edge normal 
basis_bx(1) = mesh%edges(e)%nx !nx
basis_bx(2) = mesh%edges(e)%ny !ny
basis_by(1) = mesh%edges(e)%dx/mesh%edges(e)%length
basis_by(2) = mesh%edges(e)%dy/mesh%edges(e)%length
call get_basis_change_2d(Ma2b,Mb2a,basis_bx,basis_by)

!translate the velocity to this basis 
vel_in(1) = mesh%u(c)
vel_in(2) = mesh%v(c)
vel_in_n = change_basis(Ma2b,vel_in)

!get the characteristic delta variables 
du = 0.0d0 
dv = 0.0d0 
drho = 0.0d0 
dpr = mesh%p(c) - options%pinf*pratio
c1 = -a*a*drho + dpr 
c2 = mesh%rho(c)*a*dv 
c3 = mesh%rho(c)*a*du + dpr 
c4 = -mesh%rho(c)*a*du + dpr  

!apply the boundary condition 
c4 = 0.0d0 

!get the boundary deltas 
drho = (-1.0d0/(a*a))*c1 + (1.0d0/(2.0d0*a*a))*c3 + (1.0d0/(2.0d0*a*a))*c4 
du = (1.0d0/(2.0d0*mesh%rho(c)*a))*(c3 - c4)
dv = (1.0d0/(mesh%rho(c)*a))*c2 
dpr = 0.5d0*(c3 + c4)

!evaluate the boundary velocity in the local coordinate system 
vel_b_n(1) = vel_in_n(1) + du 
vel_b_n(2) = vel_in_n(2) + dv

!evaluate boundary flux 
if (vel_b_n(1) .LT. 0.0d0) then !if inflow then revert to a wall condition
    rhob = 0.0
    ub = 0.0
    vb = 0.0
    pb = mesh%p(c)
    eb = 0.0
else !normal

    !evaluate the boundary velocity in the base coordinate system 
    vel_b = change_basis(Mb2a,vel_b_n)

    !evaluate the boundary state
    ub = vel_b(1)
    vb = vel_b(2)
    rhob = drho + mesh%rho(c) 
    pb = dpr + options%pinf*pratio
    eb = energy(pb,rhob,sqrt(ub*ub + vb*vb),options%gamma)
end if 
return 
end subroutine subsonic_pressure_outflow_bc_characteristic

!pressure outflow boundary condition (subsonic) (complex) ===============
subroutine subsonic_pressure_outflow_bc_characteristic_cpx(rhob,ub,vb,pb,eb,mesh,c,e,options,pratio)
implicit none 

!variables - inout 
integer(in64) :: e,c
complex(dp) :: rhob,ub,vb,pb,eb
complex(dp) :: pratio
type(flux_mesh_cpx) :: mesh 
type(flux_options) :: options 

!variables - local 
complex(dp) :: a,du,dv,drho,dpr,c1,c2,c3,c4
complex(dp) :: vel_in(2),vel_in_n(2),vel_b(2),vel_b_n(2)
complex(dp) :: basis_bx(2),basis_by(2)
complex(dp) :: Mb2a(2,2),Ma2b(2,2) 

!get the local speed of sound 
a = speed_of_sound_cpx(mesh%p(c),mesh%rho(c),complex(options%gamma,0.0d0))

!get the local basis coordinate system wrt to the edge normal 
basis_bx(1) = mesh%edges(e)%nx !nx
basis_bx(2) = mesh%edges(e)%ny !ny
basis_by(1) = mesh%edges(e)%dx/mesh%edges(e)%length
basis_by(2) = mesh%edges(e)%dy/mesh%edges(e)%length
call get_basis_change_2d_cpx(Ma2b,Mb2a,basis_bx,basis_by)

!translate the velocity to this basis 
vel_in(1) = mesh%u(c)
vel_in(2) = mesh%v(c)
vel_in_n = change_basis_cpx(Ma2b,vel_in)

!get the characteristic delta variables 
du = 0.0d0 
dv = 0.0d0 
drho = 0.0d0 
dpr = mesh%p(c) - options%pinf*pratio
c1 = -a*a*drho + dpr 
c2 = mesh%rho(c)*a*dv 
c3 = mesh%rho(c)*a*du + dpr 
c4 = -mesh%rho(c)*a*du + dpr  

!apply the boundary condition 
c4 = 0.0d0 

!get the boundary deltas 
drho = (-1.0d0/(a*a))*c1 + (1.0d0/(2.0d0*a*a))*c3 + (1.0d0/(2.0d0*a*a))*c4 
du = (1.0d0/(2.0d0*mesh%rho(c)*a))*(c3 - c4)
dv = (1.0d0/(mesh%rho(c)*a))*c2 
dpr = 0.5d0*(c3 + c4)

!evaluate the boundary velocity in the local coordinate system 
vel_b_n(1) = vel_in_n(1) + du 
vel_b_n(2) = vel_in_n(2) + dv

!evaluate boundary flux 
if (real(vel_b_n(1),dp) .LT. 0.0d0) then !if inflow then revert to a wall condition
    rhob = 0.0
    ub = 0.0
    vb = 0.0
    pb = mesh%p(c)
    eb = 0.0
else !normal

    !evaluate the boundary velocity in the base coordinate system 
    vel_b = change_basis_cpx(Mb2a,vel_b_n)

    !evaluate the boundary state
    ub = vel_b(1)
    vb = vel_b(2)
    rhob = drho + mesh%rho(c) 
    pb = dpr + options%pinf*pratio
    eb = energy_cpx(pb,rhob,sqrt(ub*ub + vb*vb),complex(options%gamma,0.0d0))
end if 
return 
end subroutine subsonic_pressure_outflow_bc_characteristic_cpx

end module edge_flux
