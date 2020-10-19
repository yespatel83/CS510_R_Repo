library(ggplot2)

# The bird data, the femur circumference was in cm
# bird slope
(158.4893-0.0912) / (15.85^3 - 0.7079^3)
# 0.03978331

# The dino data, the femur circumference is in mm, above (1cm = 10mm)
# dino slope using the geometric similarity 
(4000 - 55) / (53.4^3 - 10.3^3)
# 0.02609462

dino <- data.frame("Name"=c("Dino1","Dino2","Dino3","Dino4","Dino5","Dino6","Dino7","Dino8","Dino9", "Bird"),
                   "Femur"=c(103, 136, 201, 267, 348, 400, 504, 512, 534, 210),
                   "Weight"=c(55, 115, 311, 640, 1230, 1818, 3300, 3500, 4000, 368.58))

# create the subset
birdy <- subset(dino, Name == "Bird")

# plot the data
p <-ggplot(dino, aes(Femur, Weight)) + geom_point(size=4, colour="slateblue") + 
  geom_point(data=birdy, size=4, color="red") +  # this adds a red point
  geom_text(data=birdy, label="Terror Bird", vjust=3)  # this adds a label for the red point
p
p + labs (title="Plot of Dinosaur Weight vs Femur Circumference", x="Femur circumference(mm)", y="Weight (kg)") + theme (plot.title=element_text(hjust=0.5)) + geom_line(aes(Femur, Weight))

# take the cubic root of the femur length
dino_c <- data.frame("Name"=c("Dino1","Dino2","Dino3","Dino4","Dino5","Dino6","Dino7","Dino8","Dino9"),
                     "Femur"=c((10.3^3), (13.6^3), (20.1^3), (26.7^3), (34.8^3), (40.0^3), (50.4^3), (51.2^3), (53.4^3)),
                     "Weight"=c(55, 115, 311, 640, 1230, 1818, 3300, 3500, 4000))

# plot the data
p_c <-ggplot(dino_c, aes(Femur, Weight)) + geom_point(size=4, colour="slateblue") + 
  geom_smooth(method = "lm", se = FALSE)  #adds the single regression  line
p_c
p_c + labs (title="Plot of Dinosaur Weight vs Femur Circumference cubed", x="Femur circumference(cm) cubed", y="Weight (kg)") + theme (plot.title=element_text(hjust=0.5))

# regression equation
model <- lm(dino_c$Weight ~ dino_c$Femur)
coef(model)

#> coef(model)
#(Intercept) dino_c$Femur 
#101.82656951   0.02548421 
# y = 0.0255(x) + 101.83; x = femur circumference in cm

#femur of terror bird cubed to get femur_cm_cubed
x = 21^3
x
y = 0.02548421*x +  101.82656951
y
#[1] 337.8358


# AN ALTERNATE WAY TO HIGHLIGHT A POINT ON A PLOT
#highlight <- "Bird"
#dino$highlight <- ifelse(dino$Name == highlight, "Dinosaurs", "Terror Bird")
#textdf <- dino[dino$Name == highlight, ]
#mycolors <- c("Dinosaurs" = "red", "Terror Bird" = "blue")
#p <-ggplot(dino, aes(Femur, Weight)) +geom_point(size=4, aes(color=highlight))
#p
#p + labs (title="Plot of Dinosaur Femur vs Weight", x="Femur circumference(mm)", y="Weight (kg)") + theme (plot.title=element_text(hjust=0.5)) + geom_line(aes(Femur, Weight))
