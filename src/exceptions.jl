Aij_indices::Vector{Tuple{Int64,Int64,Int64}}



indice1 = [1,1,1,1,1,1,1,2,2,4,4,4,5,5,6]
indice2 = [2,3,4,5,6,7,15,3,4,5,6,7,6,7,7]
indice3 = 1:length(indice1)
unique_indices = intersect([1,2,3,4,5,6,15],xsel)

Aij_indices =[(a1,a2,a3) for (a1,a2,a3) in 
zip(indice1[unique_indices],indice2[unique_indices],indice3[unique_indices])]

for kk in 1:length(model.Aij_indices) # i are CartesianIndices
        
    i1,i2,i0 = model.Aij_indices[kk]
#    i0 = model.Aij_indices[kk]
 #   i1 = model.Aij_indices[kk]
  #  i2 = model.Aij_indices[kk]