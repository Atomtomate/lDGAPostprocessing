@testset "SVertexTools" begin
    global gImp
    global sv_up_test_red,sv_do_test_red
    global sv_up_test_full, sv_do_test_full
    F_up_full = permutedims(reshape(sv_up_test_full.data,(11,10,10)),[1,3,2])
    ind_2 = lDGAPostprocessing.indices(sv_up_test_full) 
    ind = indices(sv_up_test_full) 
    @test all(full(sv_up_test_full, ind) .≈ full(sv_up_test_red, ind))
end
