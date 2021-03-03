include("psse_to_spine.jl")
psse_path = "WP2019.raw"
db_url="sqlite:///powersystems_test2.sqlite"
psse_to_spine(psse_path, db_url)