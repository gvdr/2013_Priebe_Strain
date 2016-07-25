@everywhere using DataFrames
@everywhere using Graphs

LinksFile = "Serengeti_Baskerville/Serengeti_Baskerville_Links.csv"
Links = readtable(LinksFile,header=true,separator=',')
Predators = levels(Links[:pred_code])
Prey = levels(Links[:prey_code])
Ids = levels(vcat(Predators,Prey))

N_s = size(Ids)[1]
AdjacencyMatrix = zeros(N_s,N_s)



G_test = simple_inclist(0)
i = 1
for id::String in Ids
	add_vertex!(G_test,ExVertex(i,id))
	i += 1
end

EList = Edge[]
i = 1
for link in Links
	push(EList,Edge(i,
