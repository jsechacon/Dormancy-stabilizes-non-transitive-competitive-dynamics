
###################################################################################
########### EULER-MURAYAMA SCHEME FOR SIMULATING THE SCALING LIMIT. ###############
###################################################################################


# Function to simulate the resulting SDE with seedbanks via Euler-Maruyama
# The drift and diffusion matrix are arguments to accomodate for arbitrary 
# number of types

simulate_SDE <- function(T = 1, n = 5000, noTypes = 3, muDrift, SigmaDif,
                         x0 = rep(1/noTypes, noTypes), y0 = x0,
                         c1 = medSeed, c2 = medSeed, sigma = 1, progress = TRUE){
    
    dt <- T/n
    X <- matrix(NA, n + 1, 1 + noTypes)  # time + notypes
    Y <- matrix(NA, n + 1, noTypes)
    colnames(X) <- c("t", paste0("x", 1:noTypes))
    colnames(Y) <- paste0("y", 1:noTypes)
    
    X[1,] <- c(0, x0)
    Y[1,] <- y0
    
    ytox <- function(x, y) c1 * (y - x)
    xtoy <- function(x, y) c2 * (x - y)
    
    if ( progress ) { pb <- txtProgressBar(max = n, style = 3) }
    
    for(k in 1:n){
        x <- X[k, 2:(noTypes + 1)]
        y <- Y[k, ]
        
        driftX <- (muDrift(x) + ytox(x, y)) * dt
        driftY <- xtoy(x, y) * dt
        
        dB <- rnorm(noTypes, 0, sqrt(dt))
        diffusion <- sqrt(sigma) * SigmaDif(x) %*% dB
        
        x_new <- x + driftX + diffusion
        y_new <- y + driftY
        
        x_new[x_new < 0] <- 0
        if(sum(x_new) > 0) x_new <- x_new/sum(x_new)
        
        y_new[y_new < 0] <- 0
        if(sum(y_new) > 0) y_new <- y_new/sum(y_new)
        
        X[k + 1, ] <- c(k * dt, x_new)
        Y[k + 1, ] <- y_new
        
        if ( progress ) { setTxtProgressBar(pb, k) }
    }
    if ( progress ) { close(pb) }
    
    return(as.data.frame(X))
}

sourceCpp("./simulate_wf_seedbank.cpp") ####????###

###################################################################################
############# COMPUTATION OF FIXATION TIMES FOR THE SCALING LIMIT #################
###################################################################################

simulate_SDE_Times_rcpp <- function(
        T = 1, dt = 1e-3, nExper = 1,
        noTypes = 3, muDrift, SigmaDif,
        x0 = rep(1/noTypes, noTypes), y0 = rep(1/noTypes, noTypes), 
        c1 = 1, c2 = 1, sigma = 1,
        progress = TRUE) {
    
    Results <- data.frame(
        FixTime = numeric(nExper),
        FixType = integer(nExper),
        NumberTry = integer(nExper)
    )
    
    if (progress) pb <- txtProgressBar(max = nExper, style = 3)
    
    noSucc <- 0
    noFail <- 0
    Tmax <- T
    
    while (noSucc < nExper) {
        ans <- simulate_wf_seedbank_cpp(
            T = Tmax,
            dt = dt,
            noTypes = noTypes,
            x0 = x0,
            y0 = y0,
            drift_const = 0,       # ignored
            muR = muDrift,
            SigmaR = SigmaDif,
            c1 = c1,
            c2 = c2,
            sigma = sigma
        )
        
        FixTime <- ans$FixTime
        FixType <- ans$FixType
        
        if (FixType != -1) {
            noSucc <- noSucc + 1
            Results[noSucc, ] <- c(FixTime, FixType, noSucc + noFail)
            if (progress) setTxtProgressBar(pb, noSucc)
        } else {
            Tmax <- 1.1 * Tmax
            noFail <- noFail + 1
        }
        
        if (Tmax > 100 * T) {
            dt <- dt * 10
            T <- T * 10
        }
    }
    
    if (progress) close(pb)
    Results
}



run_SDE_Times_vec <- function(
        T = 1, dt = 1e-3, nExper = 1, Tmax = 1e3,
        noTypes = 3, muDrift, SigmaDif,
        x0 = rep(1/noTypes, noTypes),
        y0 = rep(1/noTypes, noTypes),
        c1, c2,
        sigma = 1,
        progress = TRUE
) {
    
    stopifnot(length(c1) == length(c2))
    
    nC <- length(c1)
    nTot <- nC * nExper
    
    FixTime <- rep(NA_real_, nTot)
    FixType <- rep(NA_integer_, nTot)
    c1_out  <- rep(c1, each = nExper)
    
    if (progress)
        pb <- txtProgressBar(min = 0, max = nC, style = 3)
    
    idx <- 1L
    
    for (i in seq_len(nC)) {
        
        aux <- simulate_SDE_Times_rcpp(
            T        = T,
            dt       = dt,
            nExper   = nExper,
            noTypes  = noTypes,
            muDrift  = muDrift,
            SigmaDif = SigmaDif,
            x0       = x0,
            y0       = y0,
            c1       = c1[i],
            c2       = c2[i],
            sigma    = sigma,
            progress = FALSE
        )
        
        meanTime <- mean(aux$FixTime)
        
        rng <- idx:(idx + nExper - 1L)
        FixTime[rng] <- aux$FixTime
        FixType[rng] <- aux$FixType
        idx <- idx + nExper
        
        if (progress) setTxtProgressBar(pb, i)
        
        if (meanTime > Tmax) break
    }
    
    if (progress) close(pb)
    
    data.frame(
        FixTime = FixTime,
        FixType = FixType,
        c1      = c1_out
    )
}




###################################################################################
################## FUNCTIONS TO COMPUTE THE MATRIX ZETA ##########################
###################################################################################

# Create diffusion matrices 
# The one in the notes
createSigmaOrig <- function(noTypes = 3) {
    if (!is.numeric(noTypes) || length(noTypes) != 1 || noTypes < 1) {
        stop("noTypes must be a single positive integer")
    }
    noTypes <- as.integer(noTypes)
    
    Sigma <- function(x) {
        x <- as.numeric(x)
        if (length(x) < noTypes) {
            stop("length(x) must be at least noTypes")
        }
        # Only the first `noTypes` entries of x are used
        x <- x[1:noTypes]
        
        # cumulative sums c_k = sum_{t=1}^k x_t
        csum <- cumsum(x)
        
        # allocate lower-triangular S
        S <- matrix(0, noTypes, noTypes)
        
        # jMax: stop at first index where cumulative sum equals 1,
        # otherwise use all entries
        jMax <- match(TRUE, csum == 1)
        if (is.na(jMax)) jMax <- length(x)
        
        # first diagonal
        S[1, 1] <- sqrt( x[1] * (1 - csum[1]) )
        
        if (jMax >= 2) {
            for (idx in 2:jMax) {
                # first column in row idx (off-diagonal)
                denom1 <- (1 - csum[1])
                if (denom1 > 0) {
                    S[idx, 1] <- - x[idx] * sqrt( x[1] / denom1 )
                } else {
                    S[idx, 1] <- 0
                }
                
                # middle columns j = 2, ..., idx-1
                if (idx > 2) {
                    for (jdx in 2:(idx - 1)) {
                        denom <- (1 - csum[jdx - 1]) * (1 - csum[jdx])
                        if (denom > 0) {
                            S[idx, jdx] <- - x[idx] * sqrt( x[idx] / denom )
                        } else {
                            S[idx, jdx] <- 0
                        }
                    }
                }
                
                # diagonal entry S[idx, idx]
                denom_diag <- (1 - csum[idx - 1])
                val <- 0
                if (denom_diag > 0) {
                    val <- x[idx] * (1 - csum[idx]) / denom_diag
                }
                S[idx, idx] <- ifelse(is.nan(sqrt(val)), 0, sqrt(val))
            }
        }
        
        # return the lower-triangular matrix S
        S
    }
    
    # return the factory result
    Sigma
}


############################ NUMERICALLY STABLE VERSION ####################################

# A modification to avoid division by 0
createSigmaAlt <- function(noTypes = 3) {
    Sigma <- function(x) {
        S <- matrix(0, noTypes, noTypes)
        for (i in 1:noTypes) {
            for (j in 1:noTypes) {
                S[i, j] <- sqrt(x[i]) * (as.numeric(i == j) - sqrt(x[j]))
            }
        }
        return(S)
    }
}

# diffusion matrices for the examples

Sigma3Orig <- createSigmaOrig(noTypes = 3)
Sigma3Alt <- createSigmaAlt(noTypes = 3)
Sigma5Orig <- createSigmaOrig(noTypes = 5)
Sigma5Alt <- createSigmaAlt(noTypes = 5)
Sigma9Orig <- createSigmaOrig(noTypes = 9)
Sigma9Alt <- createSigmaAlt(noTypes = 9)

# Drifts given by selection

mu3 <- function(x){
    c(
        x[1] * (x[3] - x[2]),
        x[2] * (x[1] - x[3]),
        x[3] * (x[2] - x[1])
    )
}

mu5 <- function(x){
    c(
        x[1] * (x[4] + x[5] - x[2] - x[3]),
        x[2] * (x[5] + x[1] - x[3] - x[4]),
        x[3] * (x[1] + x[2] - x[4] - x[5]),
        x[4] * (x[2] + x[3] - x[5] - x[1]),
        x[5] * (x[3] + x[4] - x[1] - x[2])
    )
}

mu9 <- function(x){
    c(
        x[1] * (x[3] + x[7] + x[8] + x[9] - x[2] - x[4] - x[5] - x[6]),
        x[2] * (x[1] + x[7] + x[8] + x[9] - x[3] - x[4] - x[5] - x[6]),
        x[3] * (x[2] + x[7] + x[8] + x[9] - x[1] - x[4] - x[5] - x[6]),
        x[4] * (x[6] + x[1] + x[2] + x[3] - x[5] - x[7] - x[8] - x[9]),
        x[5] * (x[4] + x[1] + x[2] + x[3] - x[6] - x[7] - x[8] - x[9]),
        x[6] * (x[5] + x[1] + x[2] + x[3] - x[4] - x[7] - x[8] - x[9]),
        x[7] * (x[9] + x[4] + x[5] + x[6] - x[8] - x[1] - x[2] - x[3]), 
        x[8] * (x[7] + x[4] + x[5] + x[6] - x[9] - x[1] - x[2] - x[3]),
        x[9] * (x[8] + x[4] + x[5] + x[6] - x[7] - x[1] - x[2] - x[3])
    )
}



###################################################################################
################################# PLOTTING  ######################################3
###################################################################################


#-------------------------------------------------------------------------------
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
        pivot_longer(-t, names_to = "type", values_to = "freq")
    
    df_long$type <- factor(
        df_long$type,
        levels = sort(unique(df_long$type)),
        labels = paste0("Type ", seq_along(unique(df_long$type)))
    )
    
    ggplot(df_long, aes(x = t, y = freq, fill = type)) +
        geom_area(alpha = 1) +
        scale_fill_manual(
            values = palette_OI[seq_along(levels(df_long$type))]
        ) +
        labs(
            x = "Time",
            y = "Frequency",
            title = title
        ) +
        theme_minimal(base_size = 6) +
        theme(
            plot.title = element_text(
                size = 6,
                face = "bold",
                hjust = 0.5,
                lineheight = 0.95,
                margin = margin(b = 6)
            ),
            plot.title.position = "plot",
            axis.title = element_text(size = 5),
            axis.text  = element_text(size = 4),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 7),
            panel.grid.minor = element_blank(),
            legend.key.width  = unit(1.2, "cm"),
            legend.key.height = unit(0.4, "cm"),
            plot.margin = margin(4, 4, 4, 4)
        ) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

plot_times <- function(df, title = NULL){
    
    df$c1[df$c1 == 0] <- 10^(-16.5)
    
    df_summed <- df %>%
        group_by(c1) %>%
        summarise(mean = mean(FixTime, na.rm = TRUE), .groups = "drop")
    
    ggplot(df, aes(x = c1, y = FixTime, group = c1)) +
        geom_violin(
            alpha = 1,
            fill = palette_OI[2],
            colour = palette_OI[2],
            linewidth = 0,
            scale = "width",
            na.rm = TRUE
        ) +
        geom_point(
            data = df_summed,
            aes(x = c1, y = mean),
            colour = palette_OI[1],
            inherit.aes = FALSE
        ) +
        scale_x_continuous(trans = "log10") +
        scale_y_continuous(trans = "log10") +
        labs(
            x = "Seedbank strength (c1, log scale)",
            y = "Fixation time (log scale)",
            title = title
        ) +
        theme_minimal(base_size = 6) +
        theme(
            plot.title = element_text(
                size = 6,
                face = "bold",
                hjust = 0.5,
                lineheight = 0.95,
                margin = margin(b = 6)
            ),
            plot.title.position = "plot",
            axis.title = element_text(size = 5),
            axis.text  = element_text(size = 4),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            plot.margin = margin(4, 4, 4, 4)
        )
}


#-------------------------------------------------------------------------------
# Tournament

nSteps <- 1e5
noSig <- 1e-6
medSig <- 3e-4
bigSig <- 2e-3
noSeed <- 0
medSeed <- 1e-2
bigSeed <- 1e0
startingConfig <- c(2/3, 2/9, 1/9)
startingSeed <- c(2/3, 2/9, 1/3)

# With modified diffusion matrix
set.seed(241125)
print("3 types - no seedbank - no noise")
Asim <- simulate_SDE(
    c1 = noSeed, c2 = noSeed, T = 250, n = nSteps, sigma = noSig, noTypes = 3, 
    muDrift = mu3, SigmaDif = Sigma3Alt, 
    x0 = startingConfig, y0 = startingSeed
)
print("3 types - no seedbank - med noise")
Bsim <- simulate_SDE(
    c1 = noSeed, c2 = noSeed, T = 250, n = nSteps, sigma = medSig, noTypes = 3, 
    muDrift = mu3, SigmaDif = Sigma3Alt, 
    x0 = startingConfig, y0 = startingSeed
)
print("3 types - no seedbank - big noise")
Csim <- simulate_SDE(
    c1 = noSeed, c2 = noSeed, T = 250, n = nSteps, sigma = bigSig, noTypes = 3, 
    muDrift = mu3, SigmaDif = Sigma3Alt, 
    x0 = startingConfig, y0 = startingSeed
)
print("3 types - med seedbank - no noise")
Dsim <- simulate_SDE(
    c1 = medSeed, c2 = medSeed, T = 250, n = nSteps, sigma = noSig, noTypes = 3, 
    muDrift = mu3, SigmaDif = Sigma3Alt, 
    x0 = startingConfig, y0 = startingSeed
)
print("3 types - med seedbank - med noise")
Esim <- simulate_SDE(
    c1 = medSeed, c2 = medSeed, T = 250, n = nSteps, sigma = medSig, noTypes = 3, 
    muDrift = mu3, SigmaDif = Sigma3Alt, 
    x0 = startingConfig, y0 = startingSeed
)
print("3 types - med seedbank - big noise")
Fsim <- simulate_SDE(
    c1 = medSeed, c2 = medSeed, T = 250, n = nSteps, sigma = bigSig, noTypes = 3, 
    muDrift = mu3, SigmaDif = Sigma3Alt, 
    x0 = startingConfig, y0 = startingSeed
)
print("3 types - big seedbank - no noise")
Gsim <- simulate_SDE(
    c1 = bigSeed, c2 = bigSeed, T = 250, n = nSteps, sigma = noSig, noTypes = 3, 
    muDrift = mu3, SigmaDif = Sigma3Alt, 
    x0 = startingConfig, y0 = startingSeed
)
print("3 types - big seedbank - med noise")
Hsim <- simulate_SDE(
    c1 = bigSeed, c2 = bigSeed, T = 250, n = nSteps, sigma = medSig, noTypes = 3, 
    muDrift = mu3, SigmaDif = Sigma3Alt, 
    x0 = startingConfig, y0 = startingSeed
)
print("3 types - big seedbank - big noise")
Isim <- simulate_SDE(
    c1 = bigSeed, c2 = bigSeed, T = 250, n = nSteps, sigma = bigSig, noTypes = 3, 
    muDrift = mu3, SigmaDif = Sigma3Alt, 
    x0 = startingConfig, y0 = startingSeed
)

# the_plot <- plot_stacked(Asim, "No seedbank, small noise") + 
#     plot_stacked(Bsim, "No seedbank, medium noise") + 
#     plot_stacked(Csim, "No seedbank, high noise") + 
#     plot_stacked(Dsim, "Medium seedbank, small noise") + 
#     plot_stacked(Esim, "Medium seedbank, medium noise") + 
#     plot_stacked(Fsim, "Medium seedbank, high noise") + 
#     plot_stacked(Gsim, "Strong seedbank, small noise") + 
#     plot_stacked(Hsim, "Strong seedbank, medium noise") + 
#     plot_stacked(Isim, "Strong seedbank, high noise") + 
#     plot_layout(ncol = 3, guides = 'collect') + plot_annotation(tag_levels = 'A') &
#     theme(
#         legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.box = "horizontal",
#         legend.key.width = unit(1.8, "cm"),
#         legend.key.height = unit(0.5, "cm"),
#         legend.title = element_blank()
#     )

# Images of model

pdf_convert("../3RPS.pdf", pages = 1, dpi = 1200)
pdf_convert("../wfSeedbankScheme.pdf", pages = 1, dpi = 1200)

img1 <- image_read("./3RPS_1.png", density = 1200)
img2 <- image_read("./wfSeedbankScheme_1.png", density = 1200)
img_plot1 <- ggplot() +
    annotation_raster(img1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme_minimal(base_size = 6) +
    theme_void() +
    theme(
        plot.title = element_text(
            size = 6,
            face = "bold",
            hjust = 0.5,
            lineheight = 0.95,
            margin = margin(b = 6)
        ),
        plot.title.position = "plot",
    ) +
    labs(title = "Tournament rule")
img_plot2 <- ggplot() +
    annotation_raster(img2, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme_minimal(base_size = 6) +
    theme_void() +
    theme(
        plot.title = element_text(
            size = 6,
            face = "bold",
            hjust = 0.5,
            lineheight = 0.95,
            margin = margin(b = 6)
        ),
        plot.title.position = "plot",
    ) +
    labs(title = "Realization with tournament rule")

# the_plot <- img_plot1 + img_plot2 +
#     plot_layout(
#         widths = unit(c(2.5, 7.2), c('in', 'in')), 
#         heights = unit(c(2.5), c('in'))
#     ) + plot_annotation(tag_levels = 'A') &
#     theme(
#         legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.box = "horizontal",
#         legend.key.width = unit(1.8, "cm"),
#         legend.key.height = unit(0.5, "cm"),
#         legend.title = element_blank()
#     )

# With modified diffusion matrix
set.seed(241125)
c1vec <- c(0, 10^(seq(from = -16, to = 0, length.out = 100)))
fixTimesSmallNoise <- run_SDE_Times_vec(
    nExper = 1000,
    muDrift = mu3,
    SigmaDif = Sigma3Alt,
    c1 = c1vec,
    sigma = noSig,
    c2 = c1vec,
    progress = TRUE, 
    Tmax = 700
)

fixTimesMedNoise <- run_SDE_Times_vec(
    nExper = 1000,
    muDrift = mu3,
    SigmaDif = Sigma3Alt,
    c1 = c1vec,
    sigma = medSig,
    c2 = c1vec,
    progress = TRUE, 
    Tmax = 700
)

fixTimesBigNoise <- run_SDE_Times_vec(
    nExper = 1000,
    muDrift = mu3,
    SigmaDif = Sigma3Alt,
    c1 = c1vec,
    sigma = bigSeed,
    c2 = c1vec,
    progress = TRUE, 
    Tmax = 700
)

# Plots
pA <- plot_stacked(Asim, "No seedbank, large population size")
pB <- plot_stacked(Bsim, "No seedbank, medium population size")
pC <- plot_stacked(Csim, "No seedbank, small population size")

pD <- plot_stacked(Dsim, "Medium seedbank, large population size")
pE <- plot_stacked(Esim, "Medium seedbank, medium population size")
pF <- plot_stacked(Fsim, "Medium seedbank, small population size")

pG <- plot_stacked(Gsim, "Strong seedbank, large population size")
pH <- plot_stacked(Hsim, "Strong seedbank, medium population size")
pI <- plot_stacked(Isim, "Strong seedbank, small population size")

pT1 <- plot_times(fixTimesSmallNoise, "Large population size")
pT2 <- plot_times(fixTimesMedNoise,   "Medium population size")
pT3 <- plot_times(fixTimesBigNoise,   "Small population size")

layout <- "
AABBBBBBB
CCCDDDEEE
FFFGGGHHH
IIIJJJKKK
LLLMMMNNN
"

the_plot <- img_plot1 + img_plot2 + 
    pA + pB + pC + 
    pD + pE + pF + 
    pG + pH + pI + plot_layout(guides = 'collect') +
    pT1 + pT2 + pT3 + 
    plot_annotation(tag_levels = 'A') +
    plot_layout(design = layout) &
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.width = unit(1.8, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.title = element_blank(), 
        plot.tag = element_text(size = 8)
    )

ggsave(
    "./megaFigure.pdf",
    plot = the_plot,
    width = 17.8, height = 25, units = "cm"
)
