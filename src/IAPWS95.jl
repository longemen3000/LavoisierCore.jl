struct IAPWS95 <: AbstractHelmholtzModel
    molecularWeight::Array{Float64,1}
    n00::Array{Float64,1}
    gamma00::Array{Float64,1}
    nr1::Array{Float64,1}
    d1::Array{Float64,1}
    t1::Array{Float64,1}
    nr2::Array{Float64,1}
    c2::Array{Float64,1}
    d2::Array{Float64,1}
    t2::Array{Float64,1}
    gamma2::Array{Float64,1}
    nr3::Array{Float64,1}
    d3::Array{Float64,1}
    t3::Array{Float64,1}
    alpha3::Array{Float64,1}
    beta3::Array{Float64,1}
    gamma3::Array{Float64,1}
    epsilon3::Array{Float64,1}
    nr4::Array{Float64,1}
    a4::Array{Float64,1}
    b4::Array{Float64,1}
    B::Array{Float64,1}
    C::Array{Float64,1}
    D::Array{Float64,1}
    A::Array{Float64,1}
    beta4::Array{Float64,1}

    function IAPWS95()
    
    n00= [-8.3204464837497, 6.6832105275932, 3.00632,0.012436, 0.97315, 1.2795, 0.96956, 0.24873]
    gamma00 = [0, 0, 0, 1.28728967, 3.53734222, 7.74073708, 9.24437796,27.5075105]
    
    nr1 = [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1,
    0.31802509345418, -0.26145533859358, -0.78199751687981e-2,
    0.88089493102134e-2]
    d1 = [1, 1, 1, 2, 2, 3, 4]
    t1 = [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1]
    
    nr2 = [-0.66856572307965, 0.20433810950965, -0.66212605039687e-4,
    -0.19232721156002, -0.25709043003438, 0.16074868486251,
    -0.4009282892587e-1, 0.39343422603254e-6, -0.75941377088144e-5,
    0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8,
    .36582165144204e-6, -.13251180074668e-11, -.62639586912454e-9,
    -0.10793600908932, 0.17611491008752e-1, 0.22132295167546,
    -0.40247669763528, 0.58083399985759, 0.49969146990806e-2,
    -0.31358700712549e-1, -0.74315929710341, 0.47807329915480,
    0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1,
    0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
    -0.20393486513704e-1, -0.16554050063734e-2, .19955571979541e-2,
    0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1,
    0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1,
    -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408,
    0.31777497330738, -0.11841182425981]
    c2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6,6, 6]
    d2 = [1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3,
    4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6]
    t2= [4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7, 
    10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23,23, 10, 50, 44, 46, 50]
    gamma2 = fill(1,44)
    
    nr3 = [-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4]
    d3 = fill(3,3)
    t3 = [0, 1, 4]
    alpha3 =  fill(20,3)
    beta3 = [150, 150, 250]
    gamma3 = [1.21, 1.21, 1.25]
    epsilon3 =  fill(1,3)
    
    nr4 =[-0.14874640856724, 0.31806110878444]
    a4 =[3.5, 3.5]
    b4 =[0.85, 0.95]
    B =[0.2, 0.2]
    C =[28, 32]
    D =[700, 800]
    A =[0.32,0.32]
    beta4 =[0.3, 0.3]
    return new([18.015268],n00,gamma00,nr1,d1,t1,nr2,c2,d2,t2,gamma2,nr3,d3,
    t3,alpha3,beta3,gamma3,epsilon3,nr4,a4,b4,B,C,D,A,beta4)
    end
end

@inline function _f0(model::IAPWS95,rho,T)
    
    delta = rho/322
    tau = 647.096/T
    
    res = log(delta)+model.n00[1]+model.n00[2]*tau+model.n00[3]*log(tau)
    
 for i = 4:8
        res=res+model.n00[i]*log(-expm1(-model.gamma00[i]*tau))
    end
    return res
end

function _fr(model::IAPWS95,rho,T)
    delta = rho/322
    tau = 647.096/T
   
    res=zero(promote_type(typeof(rho),typeof(T)))
    for i = 1:7
        res=res+(model.nr1[i]*delta^(model.d1[i])) * (tau^model.t1[i])
    end

    for i = 1:44
        res=res+(model.nr2[i]*delta^(model.d2[i])) * (tau^model.t2[i]) * exp(-delta^model.c2[i])
    end

    for i = 1:3
           res=res+(model.nr3[i]*delta^(model.d3[i])) * (tau^model.t3[i]) * 
           exp(-model.alpha3[i]*abs2(delta-model.epsilon3[i])-model.beta3[i]*abs(tau-model.gamma3[i]))
    end
    delta1m2 = (delta-1)^2# 
    tau1m2 = (tau-1)^2
    
    for i = 1:2
        theta = (1-tau) + model.A[i]*delta1m2^(1/(2*model.beta4[i]))
        del = theta^2 + model.B[i]*delta1m2^model.a4[i]
        psi = exp(- model.C[i]*delta1m2 - model.D[i]*tau1m2)
        res = res+model.nr4[i]*del^model.b4[i]*psi
    end
    return res
end

@inline _f(model::IAPWS95,rho,T) = _fr(model,rho,T)+_f0(model,rho,T)

function core_helmholtz(model::IAPWS95,v,T,x=[1.0]) 
   #R value calculated from molecular weight and specific gas constant
    #return 8.3143713575874*T*_f(model, molar_to_weight(1/v,[model.molecularWeight],[1.0]),T)
    #println(molar_to_weight(1/v,[model.molecularWeight],[1.0]))
    return 8.3143713575874*T*_f(model, 1/molar_to_weight(v,model.molecularWeight,[1.0]),T)
end



#empiric equations

function _p0exp(T) #Vapor–pressure equation, eq 2.5, #SI units

    d=1-T/647.096
    a = [-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502]
    return 22.064e6*exp((647.096/T)*(a[1]*d + a[2]*d^1.5+ a[3]*d^3+ a[4]*d^3.5+ a[5]*d^4+ a[6]*d^7.5))
end
function _dp0dTexp(T) #Vapor–pressure equation derivative, eq 2.5a, #SI units
    d=1-T/647.096
    a = [-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502]
p0 = _p0exp(T)
return-(p0/T)*(log(p0/22.064e6)+ a[1] + 1.5*a[2]*d^0.5+ 3*a[3]*d^2+ 3.5*a[4]*d^2.5+ 4*a[5]*d^3+ 7.5*a[6]*d^6.5)
end

function _rholsatexp(T) #Saturated liquid density equation, eq 2.6, #SI units
    d=1-T/647.096
    b = [1.99274064,1.09965342,-0.510839303,-1.75493479,-45.5170352,-6.74694450e05]
    return 322*(1.0+b[1]*d^(1.0/3.0)+b[2]*d^(2.0/3.0)+b[3]*d^(5.0/3.0)+b[4]*d^(16.0/3.0)+b[5]*d^(43.0/3.0)+b[6]*d^(110.0/3.0))
end


molecular_weight(model::IAPWS95) = model.molecularWeight #MW
compounds_number(model::IAPWS95)= 1
covolumes(model::IAPWS95) =[1.4981e-5] #10000 bar
critical_pressure(model::IAPWS95)  = [2.2064e7]
critical_volume(model::IAPWS95)  = [inv(17873.72799560906)]
critical_temperature(model::IAPWS95)  = [647.096]  
critical_density(model::IAPWS95) = [17873.72799560906]