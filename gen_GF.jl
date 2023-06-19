kG = jED.gen_kGrid(cfg["parameters"]["lattice"], 100)
 
p  = AIMParams(ϵₖ, Vₖ)
Nν = 4*(nBose+nFermi+1)

basis  = jED.Basis(length(ϵₖ) + 1);
νnGrid = jED.OffsetVector([1im * (2*n+1)*π/β for n in 0:Nν-1], 0:Nν-1)
G0W    = GWeiss(νnGrid, μ, p)

model  = AIM(ϵₖ, Vₖ, μ, U)
es     = Eigenspace(model, basis);
GImp, nden = calc_GF_1(basis, es, νnGrid, β, with_density=true)
ΣImp   = Σ_from_GImp(G0W, GImp)
gLoc   = GLoc(ΣImp, μ, νnGrid, kG)
