@testitem "Python bridge reaches the locked stim" tags = [:bridge] begin
    using PythonCall
    stim = FTShadows.pymodule(:stim)
    @test pyconvert(String, stim.__version__) == "1.16.0"
    c = stim.Circuit("H 0\nCNOT 0 1\nM 0 1")
    s = pyconvert(Matrix{Bool}, c.compile_sampler(seed=1).sample(100))
    @test size(s) == (100, 2)
    @test s[:, 1] == s[:, 2]
end
