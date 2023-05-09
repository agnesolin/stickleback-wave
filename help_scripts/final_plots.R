#### FINAL PLOTS OF DRIVERS ####


### predation vs connectivity ###


p1 = ggarrange(pred_conn_2, pred_conn_1, common.legend = T, legend = "bottom", labels = "auto", 
               font.label = list(size = 14, family = "sans", face = "bold"), ncol = 2) #, 

ggsave("figures/fig2.pdf", width = 18, height = 8.5, units = "cm", dpi = 300) 
ggsave("figures/fig2.png", width = 18, height = 8.5, units = "cm", dpi = 300) 


### temp vs distance ###

p2 = ggarrange(temp_dist_2, temp_dist_1, common.legend = T, legend = "bottom", labels = "auto", 
               font.label = list(size = 12, family = "sans", face = "bold"), ncol = 2)#, 
              

ggsave("figures/fig3.pdf", width = 18, height = 8.5, units = "cm", dpi = 300) 
ggsave("figures/fig3.png", width = 18, height = 8.5, units = "cm", dpi = 300) 


### fishing vs connectivity ###

ggarrange(fish_conn_1, fish_conn_2, common.legend = T, legend = "bottom", labels = "auto", 
          font.label = list(size = 12, family = "sans", face = "bold"), ncol = 2)#, 

ggsave("suppFigures/fish_conn.jpg", width = 18.5, height = 10, units = "cm", dpi = 300) 


### stickleback plots ###

pred_conn_3
ggsave("suppFigures/pred_conn_stick.jpg", width = 18.5, height = 12, units = "cm", dpi = 300) 

fish_conn_3
ggsave("suppFigures/fish_conn_stick.jpg", width = 18.5, height = 12, units = "cm", dpi = 300) 

temp_dist_3
ggsave("suppFigures/temp_dist_stick.jpg", width = 18.5, height = 12, units = "cm", dpi = 300) 


