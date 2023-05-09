

# max distance for connecting habitats 
max_dist = 10000 # according to model, 95% of less than this distance

dist = c(0, 2.5, 7.5, 12.5, 17.5)

p1 = c(0, cumsum(c(74, 11, 5, 3))/100)
p2 = c(0, cumsum(c(73, 15, 9, 3))/100)
p3 = c(0, cumsum(c(83, 17, 0, 0))/100)
p4 = c(0, cumsum(c(75, 25, 0, 0))/100)
p5 = c(0, cumsum(c(78, 0, 22, 0))/100)

distWeight = data.frame(
  resp = c(p1, p2, p3, p4, p5),
  dist = rep(dist,5)
)



mod =  nls(resp ~ Vm * dist/(K+dist), data = distWeight, 
           start = list(K = max(distWeight$resp)/2, Vm = max(distWeight$resp)))

weight_pred = data.frame(
  dist = seq(0,20,0.01), 
  pred = predict(mod, newdata = data.frame(dist = seq(0,20,0.01))))

link_plot = data.frame(
  resp = c(p1, p2, p3, p4, p5),
  dist = rep(dist,5),
  col = as.factor(rep(c(1, 0, 0, 0, 0), 5)))



ggplot(data = link_plot, aes(x = dist, y = resp)) +
  
  geom_point(size = 0.8, aes(colour = col)) +
  geom_line(data = weight_pred, aes(dist, y = pred), size = 0.5) +
  scale_colour_manual(values = c("black", "grey"), labels = c( "Measured", "Inferred"), name = "") +
  
  labs(x = "Distance (km)", y = "Cumulative probability") +
  
  theme_bw() +
  
  theme_sets +
  
  theme(legend.position = c(0.75, 0.2))

ggsave("suppFigures/DispersalPlot.jpg", width = 10, height = 10, units = "cm")


