## LT50 following https://lukemiller.org/index.php/2010/02/calculating-lt50-median-lethal-temperature-aka-ld50-quickly-in-r/
library(tidyverse)
library(MASS)

data = read.csv("all_LT50.csv") %>%
  mutate(Temperature = case_when(
    Treatment == "H1" ~ 34,
    Treatment == "H2" ~ 35,
    Treatment == "H3" ~ 36,
    Treatment == "H4" ~ 37
  )) %>%
  mutate(Survivorship = Alive / Initial)

f.data = data %>%
  filter(Species == "Crepidula fornicata", 
         Hour == "24") %>%
  mutate(Population = factor(Population, levels = c("Cape May", "Newport", "Beverly", "Kettle Cove", "Robbinston")))

p.data = data %>%
  filter(Species == "Crepidula plana", 
         Hour == "24")

### Binomial models

# function for printing out LT50 fro glm.dose

print.LT50 = function(x, ...) {
  M = cbind(x, attr(x, "SE"))
  dimnames(M) = list(names(x), c("LT50", "SE"))
  x = M
  NextMethod("print")
}

# Beverly
f.beverly.results = glm(cbind(Alive, Mortalities) ~ Temperature, binomial, 
                        data = f.data %>%
                          filter(Population == "Beverly"))
summary(f.beverly.results)

f.beverly.lt50 = dose.p(f.beverly.results, p = 0.5) 

f.beverly = print.LT50(f.beverly.lt50) %>%
  as.data.frame() %>%
  mutate(Population = "Beverly",
         Species = "Crepidula fornicata")

p.beverly.results = glm(cbind(Alive, Mortalities) ~ Temperature, binomial, 
                        data = p.data %>%
                          filter(Population == "Beverly"))

summary(p.beverly.results)

p.beverly.lt50 = dose.p(p.beverly.results, p = 0.5) 

p.beverly = print.LT50(p.beverly.lt50) %>%
  as.data.frame() %>%
  mutate(Population = "Beverly",
         Species = "Crepidula plana")

# Cape May
f.CapeMay.results = glm(cbind(Alive, Mortalities) ~ Temperature, binomial, 
                        data = f.data %>%
                          filter(Population == "Cape May"))
summary(f.CapeMay.results)

f.CapeMay.lt50 = dose.p(f.CapeMay.results, p = 0.5) 

f.CapeMay = print.LT50(f.CapeMay.lt50) %>%
  as.data.frame() %>%
  mutate(Population = "Cape May",
         Species = "Crepidula fornicata")

p.CapeMay.results = glm(cbind(Alive, Mortalities) ~ Temperature, binomial, 
                        data = p.data %>%
                          filter(Population == "Cape May"))
summary(p.CapeMay.results)

p.CapeMay.lt50 = dose.p(p.CapeMay.results, p = 0.5) 

p.CapeMay = print.LT50(p.CapeMay.lt50) %>%
  as.data.frame() %>%
  mutate(Population = "Cape May",
         Species = "Crepidula plana")

# Kettle Cove
f.KettleCove.results = glm(cbind(Alive, Mortalities) ~ Temperature, binomial, 
                        data = f.data %>%
                          filter(Population == "Kettle Cove"))
summary(f.KettleCove.results)

f.KettleCove.lt50 = dose.p(f.KettleCove.results, p = 0.5) 

f.KettleCove = print.LT50(f.KettleCove.lt50) %>%
  as.data.frame() %>%
  mutate(Population = "Kettle Cove",
         Species = "Crepidula fornicata")

# Robbinston
f.Robbinston.results = glm(cbind(Alive, Mortalities) ~ Temperature, binomial, 
                           data = f.data %>%
                             filter(Population == "Robbinston"))
summary(f.Robbinston.results)

f.Robbinston.lt50 = dose.p(f.Robbinston.results, p = 0.5) 

f.Robbinston = print.LT50(f.Robbinston.lt50) %>%
  as.data.frame() %>%
  mutate(Population = "Robbinston",
         Species = "Crepidula fornicata")

# Newport
f.Newport.results = glm(cbind(Alive, Mortalities) ~ Temperature, binomial, 
                           data = f.data %>%
                             filter(Population == "Newport"))
summary(f.Newport.results)

f.Newport.lt50 = dose.p(f.Newport.results, p = 0.5) 

f.Newport = print.LT50(f.Newport.lt50) %>%
  as.data.frame() %>%
  mutate(Population = "Newport",
         Species = "Crepidula fornicata")




## Plot
f.data = f.data %>%
  mutate(Species = as.factor(Species), 
         Population = factor(Population, levels = c("Robbinston", "Kettle Cove", "Beverly","Newport", "Cape May"))) 

ggplot(f.data, aes(Temperature, Survivorship)) +
  geom_point(outlier.shape = NA) +
  geom_smooth(method = glm,
              method.args = list(family = "binomial")) +
  theme_bw() +
  facet_wrap(. ~ Population, ncol = 1)

ggsave("survivorship.jpeg", last_plot(), height = 5.6, width = 1.5, units = "in")

ggplot(p.data, aes(Temperature, Survivorship)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = glm,
              method.args = list(family = "binomial")) +
  facet_wrap(Population ~ Hour)


# LT50 + SE

lt50 = rbind(p.beverly, p.CapeMay, f.beverly, f.CapeMay, f.Robbinston, f.Newport, f.KettleCove) %>%
  mutate(Species = as.factor(Species), 
         Population = factor(Population, levels = c("Cape May", "Newport", "Beverly", "Kettle Cove", "Robbinston"))) 



f.lt50 = lt50 %>%
  filter(Species == "Crepidula fornicata") %>%
  rename(Site = Population)

# plot
library(forcats)

colours = c("Robbinston" = "#264D59", "Kettle Cove" = "#43978D","Beverly" = "#F9E07F", "Newport" = "#F9AD6A","Cape May" = "#D46C4E")

ggplot(f.lt50, aes(LT50, Site, color = Site)) +
  geom_pointrange(aes(xmin = LT50-SE, xmax = LT50+SE, colour = Site), position = position_dodge(width = 0.5)) + 
  geom_point(size=4, pch=21,aes(fill=Site), colour = "black")  +
  geom_point(size=3.5, pch=21,aes(fill=Site)) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(legend.position = "none") 
 

ggsave("lt50_plot.jpeg", last_plot(), height = 4, width = 4, units = "in")









