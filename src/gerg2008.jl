include("utilities.jl")

struct GERG2008
    N::Int64
    x_selected::Array{Float64,1}
    rhoc::Array{Float64,1}
    Tc::Array{Float64,1}
    nrmat::Array{Float64,2}
    zetamat::Array{Float64,2}
    function GERG2008(data)
        Tc = data.criticalTemperature
        N = length(Tc)
        rhoc = data.criticalDensity
        x_selected = data.GERG2008
        nramat = Array{Float64}(undef,(7,21))
        zetamat = Array{Float64}(undef,(7,21))
        zetamat[:,1]=[19.597508817,-83.959667892,3.00088,0.76315,0.00460,8.74432,-4.46921]
        #zetamat[:,2]= "s" #..... to 21
    return new(N,x_selected,rhoc,Tc,nrmat,zetamat)
    end
end


function _f0(model::GERG2008,rho,T,x)
    #idea:x_selected should be in the data as GERG2008number
    # example: x_selected = [4 21 5 3]
    #so the data passed to the function inmediately recognizes the order
 
   
    #delta = rho/rhoc,
    #tau = Tc/T
    

    RR = 8.314472/8.314510
    #nrmat = "ssdasdasd" #matrix of n values (columnar respect of the compunds)
    #zmat = "ssdasdasd" #valid from 4 to 7, fill the rest with NA)
    
    vecres=zeros(model.N)
        
    for ii = 1:model.N
        i = model.x_selected(ii)
        nr1 =  model.nrmat[:,i]
        delta = rho/model.rhoc[i]
        tau = model.Tc[i]/T
        zeta = model.zetamat[:,i]

        ao1 =  nr1[1]+nr1[2]*delta+nr1[3]*log(delta)+
        nr1[4]*log(abs(sinh(zeta(4)*delta)))- nr1[5]*log(cosh(zeta[5]*delta)) +
        nr1[6]*log(abs(sinh(zeta(6)*delta)))- nr1[7]*log(cosh(zeta[7]*delta))
        ao2 = log(delta)
        ao = RR*ao1+ao2
        vecres(i)=x[i]*(ao+log(x[i]))
        end

        return sum(vecres)
end



