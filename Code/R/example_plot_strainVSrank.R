library(GGally)
library(ggplot2)

ggplot() +
  geom_point(data=svdd, aes(x=index,y=value)) +
  geom_line(data=svdd, aes(x=index,y=value)) + 
  labs(title="Baskerville's Serengeti Food Web Singular Values")


cols <- c("2"="blue","3"="red","4"="green")
ggplot() +
  geom_point(data=Strain_2, size=3, fill="darkblue", colour="darkblue", aes(x=degree, y=strain), shape=21) +
  stat_smooth(data=Strain_2, fill="blue", colour="blue", alpha=0.1, aes(x=degree, y=strain, colour="2"), method=lm, level = 0.95) +
  geom_point(data=Strain_3, size=3, fill="darkred", colour="darkred", aes(x=degree, y=strain), shape=21) +
  stat_smooth(data=Strain_3, fill="red", colour="red", alpha=0.1, aes(x=degree, y=strain, colour="3"), method=lm, level = 0.95) +
  geom_point(data=Strain_4, size=3, fill="darkgreen", colour="darkgreen", aes(x=degree, y=strain, rank="4"), shape=21) +
  stat_smooth(data=Strain_4, fill="green", colour="green", alpha=0.1, aes(x=degree, y=strain, colour="4"), method=lm, level = 0.95) + 
  labs(title="Strain vs. Node Degree") +
  scale_colour_manual(name="Rank",values=cols)

ggplot() +
  geom_point(data=SD_2, size=3, fill="darkblue", colour="darkblue", aes(x=degree, y=strain), shape=21) +
  stat_smooth(data=SD_2, fill="blue", colour="blue", alpha=0.1, aes(x=degree, y=strain, colour="2"), method=lm, formula=y ~poly(x,3), level = 0.95) +
  geom_point(data=SD_3, size=3, fill="darkred", colour="darkred", aes(x=degree, y=strain), shape=21) +
  stat_smooth(data=SD_3, fill="red", colour="red", alpha=0.1, aes(x=degree, y=strain, colour="3"), method=lm, formula=y ~poly(x,3), level = 0.95) +
  geom_point(data=SD_4, size=3, fill="darkgreen", colour="darkgreen", aes(x=degree, y=strain, rank="4"), shape=21) +
  stat_smooth(data=SD_4, fill="green", colour="green", alpha=0.1, aes(x=degree, y=strain, colour="4"), method=lm, formula=y ~poly(x,3), level = 0.95) + 
  labs(title="Strain vs. Node Degree") +
  scale_colour_manual(name="Rank",values=cols)

ggplot() +
  geom_point(data=SD_2, size=3, fill="darkblue", color="darkblue", aes(x=degree, y=strain, rank="2"), shape=21) +
  stat_smooth(data=SD_2, fill="blue", color="blue", alpha=0.1, aes(x=degree, y=strain), method=glm, formula= y ~ ns(x,3), level = 0.95) +
  geom_point(data=SD_3, size=3, fill="darkred", color="darkred", aes(x=degree, y=strain, rank="3"), shape=21) +
  stat_smooth(data=SD_3, fill="red", color="red", alpha=0.1, aes(x=degree, y=strain), method=glm, formula= y ~ ns(x,3), level = 0.95) +
  geom_point(data=SD_4, size=3, fill="darkgreen", color="darkgreen", aes(x=degree, y=strain, rank="4"), shape=21) +
  stat_smooth(data=SD_4, fill="green", color="green", alpha=0.1, aes(x=degree, y=strain), method=glm, formula= y ~ ns(x,3), level = 0.95)


ggparcoord(Strains, columns=c(2:4), scale="globalminmax",groupColumn = "Node") + geom_line(alpha=0.5)+ xlab("Rank dimension") + ylab("Strain") + theme(legend.position="none")
