# looking at PCAs
library(cowplot)


# FL
prow1 <- plot_grid( acerv_sub_DEG_PCA_plot + theme(legend.position="none"),
                   mcav_DEG_PCA_plot + theme(legend.position="none"),
                   ofav_DEGPCAplot + theme(legend.position="none"),
                   align = 'vh',
                   labels = c("A", "B", "C"),
                   hjust = -5,
                   nrow = 1
)

legend_b1 <- get_legend(acerv_sub_DEG_PCA_plot + theme(legend.position="right"))

# add the legend to plot
p1 <- plot_grid( prow1, legend_b1, ncol = 2, rel_heights = c(1, .05), rel_widths = c(1, 0.1))
p1
ggsave("Output/Figs/PCA_FL.png", p1, width = 15, height = 5)


# HI
prow2 <- plot_grid( mcap_DEGPCAplot + theme(legend.position="none"),
                    pacuta_DEGPCAplot + theme(legend.position="none"),
                    plob_DEGPCAplot + theme(legend.position="none"),
                    align = 'vh',
                    labels = c("D", "E", "F"),
                    hjust = -5,
                    nrow = 1
)

legend_b2 <- get_legend(mcap_DEGPCAplot + theme(legend.position="right"))

# add the legend to plot
p2 <- plot_grid( prow2, legend_b2, ncol = 2, rel_heights = c(1, .05), rel_widths = c(1, 0.1))
p2
ggsave("Output/Figs/PCA_HI.png", p2, width = 15, height = 5)

# both FL and HI plots 
full <- plot_grid(p1, p2, nrow = 2)
full
ggsave("Output/Figs/PCA_all.png", full, width = 15, height = 10)



