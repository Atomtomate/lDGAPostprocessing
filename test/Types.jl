@testset "basics" begin 
    global sv_read_1
    global sv_read_2
    global sv_01
    global sv_02
    @test typeof(sv_01) === SVertex{Int64}
    @test sv_01.n == UInt32(2)
    @test sv_01.nB == UInt32(5)
    @test sv_01.nB2 == UInt32(25)
    @test sv_01.offset == UInt32(0)
    @test sv_01[1] == 4
    @test sv_01[2] == 5
    @test sv_01[0,0,1] == 4
end
