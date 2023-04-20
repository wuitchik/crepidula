# making fst correlation plot

library(tidyverse)
library(corrplot)

data = read.csv("fst_cleaned.csv") %>%
  mutate(Y = factor(Y, levels = c("Cape_May", "Newport", "Beverly", "Kettle_Cove", "Robbinston"))) %>%
  mutate(X = factor(X, levels = c("Cape_May", "Newport", "Beverly", "Kettle_Cove", "Robbinston")))

ggplot(data, aes(X, Y, fill= value)) + 
  geom_tile((aes(fill = value)), colour="black")+
  geom_text(aes(label = value),size=8 ) +
  theme_classic(base_size = 22)+
  scale_fill_gradient2(low = "white", high="orange", midpoint = 0.022) +
  guides(fill=guide_legend(title="Fst"))+
  xlab(label = "")+ ylab(label="")

corrplot(data, method = "circle", 
         type = "upper", 
         order="alphabet", #hierarchical clustering order of factors
         tl.col="black", tl.srt=45, #text label color and rotation
        # p.mat = time30dip_p.mat, sig.level = 0.05, insig = "blank" #combine with significance
         #title = "T30" #can add title, but adds it really high for some reason
) 
