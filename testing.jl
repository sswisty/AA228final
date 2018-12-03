include("OrbitMechFcns.jl")

oe₁ = OrbitalElementVec(7500,0.015,1*pi/180,pi/6,0,pi/4)

R,V = OE2ECI(oe₁)

oe₂,ang = ECI2OE(R,V)


ICs = [R;V]
tspan = (0.0,9000.0)
params = [398600.4418;6378.137;1.0826267e-3]
prob = ODEProblem(GravitationalForce,ICs,tspan,params)
sol = solve(prob,saveat = 2)

xdata = [r[1] for r in sol.u]
ydata = [r[2] for r in sol.u]
zdata = [r[3] for r in sol.u]
# xda,yda,zda = [r[1],r[2],r[3] for r in sol.u # if only...

plot1 = plot(xdata,ydata,zdata,xlim=(-8000,8000), ylim=(-8000,8000), zlim=(-8000,8000),
       title = "Orbit via ODEsolver")

# surface(xdata,ydata,zdata)

t = collect(0:5:86400) # one day, every 5 seconds
epoch = UT12MJD(11,21,2018,0,0,0)
Reci,Veci = MeanMotionProp(oe₁,epoch,t/86400)

xdata = [r[1] for r in Reci]
ydata = [r[2] for r in Reci]
zdata = [r[3] for r in Reci]

plot2 = plot(xdata,ydata,zdata,xlim=(-8000,8000), ylim=(-8000,8000), zlim=(-8000,8000),
       title = "Orbit via Mean Motion Prop")

display(plot1)
display(plot2)


# Load Earth topo for ground track
EP = EarthGroundPlot()
GMST = MJD2GMST(epoch.+t/86400)
Recef,Vecef = ECI2ECEF(Reci,Veci,GMST)
ϕ,λ,h = ECEF2GEO(Recef)

plot!(λ,ϕ)
