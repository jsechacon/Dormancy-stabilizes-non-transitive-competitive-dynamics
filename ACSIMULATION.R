##!/usr/bin/env Rscript
#-------------------------------------------------------------------------------
# Loading the packages using pacman
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, ggthemes, dplyr, tidyr, patchwork)

#-------------------------------------------------------------------------------

###################################################################################
# EULER-MURAYAMA SCHEME FOR SIMULATING THE MAJORITY VOTING/ALLEN CAHN SCALING LIMIT  
###################################################################################

# Function to simulate the Allen-Cahn SDE with seedbanks via Euler-Maruyama
# The drift and diffusion matrix are arguments to accomodate for arbitrary 
# number of types

AllenCahnSimulation <- function(T = 1, n = 5000, c1 = 1, c2 = 1, s = 1, 
                                x0 = 1/2, y0 = 1/2, sigma = 1,
                                progress = TRUE) {
    dt <- T/n
    X <- matrix(NA, n + 1, 1 + 2)  # time + notypes
    Y <- matrix(NA, n + 1, 2)
    colnames(X) <- c("t", paste0("x", 1:2))
    colnames(Y) <- paste0("y", 1:2)
    
    X[1,] <- c(0, x0, 1 - x0)
    Y[1,] <- c(y0, 1 - y0)
    
    ytox <- function(x, y) c1 * (y - x)
    xtoy <- function(x, y) c2 * (x - y)
    
    mu <- function(x) s * x * (1 - x) * (2 * x - 1)
    Sigma <- function(x) sqrt(sigma * x * (1 - x))
    
    if ( progress ) { pb <- txtProgressBar(max = n, style = 3) }
    
    for (k in 1:n) {
        x <- X[k, 2]
        y <- Y[k, 1]
        
        driftX <- (mu(x) + ytox(x, y)) * dt
        driftY <- xtoy(x, y) * dt
        
        dB <- rnorm(1, 0, sqrt(dt))
        diffusion <- Sigma(x) %*% dB
        
        x_new <- x + driftX + diffusion
        y_new <- y + driftY
        
        x_new[x_new < 0] <- 0
        x_new[x_new > 1] <- 1
        
        y_new[y_new < 0] <- 0
        y_new[y_new > 1] <- 1
        
        X[k + 1, ] <- c(k * dt, x_new, 1 - x_new)
        Y[k + 1, ] <- c(y_new, 1 - y_new)
        
        if ( progress ) { setTxtProgressBar(pb, k) }
    }
    
    if ( progress ) { close(pb) }
    
    return(as.data.frame(X))
}

###################################################################################
################################# PLOTTING  ######################################3
###################################################################################


# Plotting specifications to create a stacked plot
# The palette is colour blind safe (?)

palette_OI <- c(
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#56B4E9",
    "#CC79A7",
    "#F0E442",
    "#E66101",
    "#007373",
    "#808080"
)

plot_stacked <- function(df, title = NULL){
    df_long <- df %>%
        pivot_longer(-t, names_to="type", values_to="freq")
    
    df_long$type <- factor(
        df_long$type, 
        levels = sort(unique(df_long$type)),
        labels = paste0("Type ", seq_along(unique(df_long$type)))
    )
    
    ggplot(df_long, aes(x = t, y = freq, fill = type)) +
        geom_area(alpha = 1) +
        geom_line(
            data = df, inherit.aes = FALSE, aes(x = t, y = 1-x1), 
            alpha = 0.3, linewidth = 0.1, colour = 'black'
        ) + 
        scale_fill_manual(
            values = palette_OI[1:length(unique(df_long$type))],
            labels = levels(df_long$type)
        ) +
        theme_minimal(base_size = 10) +  # axis + labels for poster
        labs(
            x = "Time",
            y = "Frequency",
            title = title
        ) +
        theme(
            plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 11),
            axis.text  = element_text(size = 9),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 9),
            panel.grid.minor = element_blank(),
            legend.key.width  = unit(1.2, "cm"),
            legend.key.height = unit(0.4, "cm"),
            plot.margin = margin(2, 2, 2, 2)
        ) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

#-------------------------------------------------------------------------------
# Examples
set.seed(241125)

print("Simulating Allen-Cahn with s = 1")
ac2.1 <- AllenCahnSimulation(T = 100, n = 1e5, c1 = 1, c2 = 1, s = 1, sigma = 10)
ac2.2 <- AllenCahnSimulation(T = 100, n = 1e5, c1 = 10, c2 = 1, s = 1, sigma = 10)
ac2.3 <- AllenCahnSimulation(T = 100, n = 1e5, c1 = 100, c2 = 1, s = 1, sigma = 10)
ac2.4 <- AllenCahnSimulation(T = 100, n = 1e5, c1 = 1000, c2 = 1, s = 1, sigma = 10)

the_plot <- plot_stacked(ac2.1, "Weak seedbank: c1 = 1") + 
    plot_stacked(ac2.2, "Moderate seedbank: c1 = 10") + 
    plot_stacked(ac2.3, "Strong seedbank: c1 = 100") + 
    plot_stacked(ac2.4, "Very strong seedbank: c1 = 1000") + 
    plot_layout(ncol = 2, guides = 'collect') + plot_annotation(tag_levels = 'A') &
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.width = unit(1.8, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.title = element_blank()
    ) 

ggsave(
    "./diffAllenCahnParams.pdf",
    plot = the_plot,
    width = 10 * 2.5 / 3, height = 10 * 2 / 3, units = "in"
)
