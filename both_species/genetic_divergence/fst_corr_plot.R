# making fst correlation plot

library(tidyverse)
library(corrplot)

data = read.csv("both_fst.csv") %>%
  mutate(Fornicata = factor(Fornicata, ordered= TRUE, levels = c("Cape May", "Newport", "Beverly", "Kettle Cove", "Robbinston"))) %>%
  mutate(Plana = factor(Plana, ordered= TRUE, levels = c("Cape May", "Newport", "Beverly", "Kettle Cove")))

ggplot(data, aes(Plana, Fornicata, fill= Fst)) + 
  geom_tile((aes(fill = Fst)), colour="black")+
  geom_text(aes(label = Fst),size=8 ) +
  theme_classic(base_size = 22)+
  scale_fill_gradient2(low = "white", high="orange", midpoint = 0.4819675) +
  guides(fill=guide_legend(title="Fst"))+
  xlab(label = "")+ ylab(label="")

corrplot(data, method = "circle", 
         type = "upper", 
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
         # p.mat = time30dip_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         #title = "T30" #can add title, but adds it really high for some reason
) 