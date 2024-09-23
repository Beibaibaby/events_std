

jie_para = 4.0
jei_para = -18.0 * 1.2
jii_para = -16.0
jee_para = 10.0

taui = 10
taue = 15

sqrtK = sqrt(800)

jie = jie_para / (taui * sqrtK)
jei = jei_para / (taue * sqrtK)
jii = jii_para / (taui * sqrtK)
jee = jee_para / (taue * sqrtK)

println("jie: ", jie)
println("jei: ", jei)
println("jii: ", jii)
println("jee: ", jee)

using Statistics
for i in 1:10
   println(rand(5))
end