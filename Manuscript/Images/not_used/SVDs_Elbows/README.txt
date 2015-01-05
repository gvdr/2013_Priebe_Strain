> getElbows(SD_svds_trunc) #de Visser svd truncated so to exclude the zeros
[1]  6 24 35
> getElbows(SB_svds_trunc) #Baskerville svd truncated so to exclude the zeros
[1]  3 12 22
> getElbows(SB_svdd)
[1] 11 29 32
> getElbows(SD_svdd)
[1] 17 37 46

#The best ones are the smaller ones
> getElbows(SD_svds_trunc,n=1)
[1] 6
> getElbows(SB_svds_trunc,n=1)
[1] 3
> getElbows(SB_svdd,n=1)
[1] 11
> getElbows(SD_svdd,n=1)
[1] 17