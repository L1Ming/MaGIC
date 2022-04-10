library(ggplot2)


## defien the compare means color
p_color = "#082032"

## define my_theme
my_theme = theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1.5),
        axis.ticks = element_line(size=1),
        text = element_text(face = 'bold'),
        plot.title = element_text(hjust=0.5))
my_theme_strip = my_theme +
  theme(strip.background = element_blank())
# save(my_theme, my_theme_strip, file = "Rdata/my_theme.Rdata")