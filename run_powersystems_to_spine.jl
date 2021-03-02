include("powersystems_to_spine.jl")
psse_path = joinpath(@__DIR__, "WP2019.raw")
#db_url = raw"sqlite:///"*joinpath(@__DIR__, "powersystems_test.sqlite") 
db_url=raw"sqlite:///D:\workspace\Spine\Network_Prune\powersystems_test2.sqlite"
@info db_url
#psse_to_spine(psse_path, db_url)
aggregate_network(db_url)
