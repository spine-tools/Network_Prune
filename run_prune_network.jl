include("prune_network.jl")
db_url="sqlite:///powersystems_test.sqlite"
prune_network(db_url)
