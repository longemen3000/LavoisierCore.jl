using JuliaDB, DataFrames, XLSX

#I have an xlsx with a lot of data from DIPPR 801 (2197 compounds and their properties)
#using JuliaDB, i put the data in a julia-friendly format

function main1() #load the table from the xlsx
df = DataFrame(XLSX.readtable("DIPPR 801.xlsx", "DIPPR 801")...)
t = table(df)
println(0)


for i = 1:11
    st = Array{Union{String,Missing}}(undef,2197)
    for j = 1:2197
        x11 = df[i][j]
        st[j]= "$x11"
    end
    t = setcol(t,i,st)
    println(i)
end


for i = 12:length(colnames(t))
    fl  = Array{Union{Float64,Missing}}(undef,2197)
    for j = 1:2197
        fl[j]=df[i][j]
    end
    t = setcol(t,i,fl)
    println(i)
end

return t

end

#function searchchem(t::IndexedTable,query::String)
#    regex_q = Regex(query)
#    for i = 

#0 - from 1 to 100
#1 - 101:200 -> 