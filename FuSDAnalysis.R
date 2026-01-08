## FuSD Analysis
## Date: 070126
## Author: Emmanuel Ambriz
# Required packages
library(ggplot2)
library(ggnewscale)
library(fdaoutlier)
library(tidyverse)
library(depthTools)
library(ggpubr)
library(fda)
library(ggnewscale)   
library(grid)
library(tibble)
library(RColorBrewer)


rm(list=ls())

load("LSTNDVI18to24.RData")
source("FuSDPackage060126.R")
source("GraphsAux.R")

J=6 # 2024
per=(J*12+1):((J+1)*12)
lst=Dlst[per,]
Fdata=t(lst)
sttdVar=FuSD(Fdata)

### Graphical parameters ####
ylabel="LST"
xlabel="t"
xrangelabs=1:12

## Quantile Plot ####

probs <- seq(0.01, 0.99, length = 150)

Qeval <- function(p) FunctionalSTtDQuantile(p, Fdata, FuSD)[[1]]
Q     <- t(sapply(probs, Qeval))          # (n_alpha x 12)

digits <- 8
keys   <- apply(round(Q, digits), 1, paste, collapse = "\r")
keep   <- ave(probs, keys, FUN = function(p) seq_along(p) == which.max(p))
keep   <- as.logical(keep)
Q      <- Q[keep, , drop = FALSE]
probs  <- probs[keep]

xf      <- seq_len(ncol(Q))
n_curv  <- nrow(Q)
df_all  <- data.frame(
  tt    = rep(xf, times = n_curv),
  st    = as.vector(t(Q)),
  alpha = rep(probs, each = length(xf)),
  id    = rep(seq_len(n_curv), each = length(xf))
)

# Median
qmed <- FunctionalSTtDQuantile(0.5, Fdata, FuSD)[[1]]
df_med <- data.frame(tt = xf, st = qmed)

# --- Plot ---
p_quant <- ggplot(df_all, aes(tt, st, group = id, colour = alpha)) +
  geom_path(linewidth = 1.2) +
  geom_path(data = df_med, aes(tt, st), inherit.aes = FALSE,
            colour = "paleturquoise1", linewidth = 2) +
  # Axis
  scale_x_continuous(breaks = xf, labels = xf, expand = c(0, 0)) +
  scale_y_continuous(limits = c(9.5, 60), breaks = seq(10, 60, 10), expand = c(0, 0)) +
  # Gradient
  scale_color_gradient2(
    name     = expression(atop(alpha, "")),
    low      = "blue",
    mid      = "paleturquoise1",
    high     = "red",
    midpoint = 0.5,
    limits   = c(0, 1),
    breaks   = seq(0, 1, 0.25),
    guide    = guide_colorbar(
      direction      = "vertical",
      title.position = "top",
      title.hjust    = 0.5,
      barheight      = grid::unit(120, "pt"),
      barwidth       = grid::unit(22,  "pt"),
      ticks          = FALSE,
      label.theme    = element_text(size = 10)
    )
  ) +
  labs(x = "t", y = expression(Q[alpha](t))) +
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    legend.position.inside      = c(0.98, 0.5),
    legend.justification = c(1, 0.5),
    legend.box         = "vertical",
    legend.title       = element_text(face = "bold", size = 16),  # α más grande y en negritas
    legend.text        = element_text(size = 14, margin = margin(t = 0, r = 6, b = 0, l = 8),face = "bold"),
    legend.background    = element_rect(fill = scales::alpha("white", 0.6), colour = NA),
    legend.key.width   = grid::unit(8,  "pt"),
    legend.key.height  = grid::unit(18, "pt"),
    legend.spacing.x   = grid::unit(6,  "pt"),
    legend.spacing.y   = grid::unit(4,  "pt"),
    legend.box.spacing = grid::unit(6,  "pt"),
    legend.margin      = margin(0, 6, 0, 6),
    axis.title.x       = element_text(margin = margin(t = 6)),
    axis.title.y       = element_text(margin = margin(r = 6)),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.25, colour = "#F6F7F9"),
    panel.border       = element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
    plot.background    = element_rect(fill = "white", colour = NA),
    plot.margin        = margin(8, 12, 8, 10) # margenes parecidos a Bands
  )


## Outlier Trimming ####

res_dir <- dir_out(Fdata, return_distance = TRUE)

# Detect outliers via directional outlyingness

threshold = sqrt(CerioliOutlierDetection::hr05CriticalValue(
  CerioliOutlierDetection::hr05AdjustedDF(158, 2), 2, 0.01))

idx_shape_outliers <- which(res_dir$distance > threshold)

Fdata2 <- Fdata[-idx_shape_outliers, ]

### Outlier Trimming END

## ALPHA

G=list()

DF=list()

SB=list()

A=rev(seq(0.1,0.9,0.025))

Fdata.df=as.data.frame(Fdata)
names(Fdata.df)=paste0("M",names(Fdata.df))
N=length(sttdVar)

## Main ####

for(i in 1:length(A))
{

alpha=A[i]
  
mn=FunctionalSTtDQuantile(alpha/2,Fdata,FuSD)
md=FunctionalSTtDQuantile(0.5,Fdata,FuSD)
mx=FunctionalSTtDQuantile(1-alpha/2,Fdata,FuSD)

PBind=(sttdVar>=sttdVar[mn[[2]]]) &  (sttdVar <= sttdVar[mx[[2]]])


### Global Skweness
IndNEG=(sttdVar>=sttdVar[mn[[2]]]) &  (sttdVar <= sttdVar[md[[2]]])
IndPOS=(sttdVar>=sttdVar[md[[2]]]) &  (sttdVar <= sttdVar[mx[[2]]])

FDBP=Fdata[IndNEG,]
nFDBP=dim(FDBP)[1]
L2_dist=c()
for(k1 in 1:(nFDBP-1))
{
  x=FDBP[k1,]
  for (k2 in 2:nFDBP) {
    y=FDBP[k2,]
    L2_dist=c(L2_dist,sqrt(sum((x - y)^2)))  
  }
  
}
L2_distNeg=max(L2_dist)

FDBP=Fdata[IndPOS,]
nFDBP=dim(FDBP)[1]
L2_dist=c()
for(k1 in 1:(nFDBP-1))
{
  x=FDBP[k1,]
  for (k2 in 2:nFDBP) {
    y=FDBP[k2,]
    L2_dist=c(L2_dist,sqrt(sum((x - y)^2)))  
  }
  
}
L2_distPos=max(L2_dist)

SKG=(L2_distPos-L2_distNeg)/(L2_distPos+L2_distNeg)

## Quantile-based Skewness 

x=mn[[1]]
y=mx[[1]]
z=md[[1]]

L2_distNeg=sqrt(sum((x - z)^2))
L2_distPos=sqrt(sum((y - z)^2))

SKQ=(L2_distPos-L2_distNeg)/(L2_distPos+L2_distNeg)


PBlab=rep('out.PB',N)
PBlab[PBind]="PB"

FDBP=Fdata[PBind,]

nFDBP=dim(FDBP)[1]

L2_dist=c()

for(k1 in 1:(nFDBP-1))
{
  x=FDBP[k1,]
  for (k2 in 2:nFDBP) {
    y=FDBP[k2,]
    L2_dist=c(L2_dist,sqrt(sum((x - y)^2)))  
  }
  
}
SP2=max(L2_dist)


df_long <- Fdata.df %>%
  mutate(curva = factor(1:N),
         grupo = factor(PBlab)) %>%
  pivot_longer(
    cols = starts_with("M"),   # columns V1, V2, ...
    names_to = "t_index",
    names_prefix = "M",
    values_to = "valor"
  ) %>%
  mutate(t_index = as.numeric(t_index))  

Cuant=data.frame(x=per,minimos=apply(rbind(mn[[1]],mx[[1]]),2,min),mediana=md[[1]],maximos=apply(rbind(mn[[1]],mx[[1]]),2,max))
SP2fpb=sqrt(sum((Cuant$minimos-Cuant$maximos)^2))

IndAdj=PBadjusted(alpha,Fdata,FuSD,SP2,df_long)
Cons=data.frame(x=per,minimos=apply(Fdata[IndAdj[1:2],],2,min),maximos=apply(Fdata[IndAdj[1:2],],2,max))
SP2tpb=sqrt(sum((Cons$minimos-Cons$maximos)^2))

## Tight Skewness 

x=Fdata[IndAdj[1],]
y=Fdata[IndAdj[2],]
z=md[[1]]

L2_distNeg=sqrt(sum((x - z)^2))
L2_distPos=sqrt(sum((y - z)^2))

SKT=(L2_distPos-L2_distNeg)/(L2_distPos+L2_distNeg)



mbd=fbplot(t(Fdata),
       method='MBD',plot=F)$depth


mbdInd=(mbd>quantile(mbd,alpha))
FBP=data.frame(x=per,minimos=apply(Fdata[mbdInd,],2,min),maximos=apply(Fdata[mbdInd,],2,max))
SP2dbb=sqrt(sum((FBP$minimos-FBP$maximos)^2))

mbd=fbplot(t(Fdata2),
           method='MBD',plot=F)$depth

mbdInd=(mbd> quantile(mbd,alpha))
FBP2=data.frame(x=per,minimos=apply(Fdata2[mbdInd,],2,min),maximos=apply(Fdata2[mbdInd,],2,max))
SP2dbb2=sqrt(sum((FBP2$minimos-FBP2$maximos)^2))

## Error Rates

ER=data.frame(df_long,Qmin=rep(Cuant$minimos,N),Qmax=rep(Cuant$maximos,N),
              Cmin=rep(Cons$minimos,N),Cmax=rep(Cons$maximos,N),
              FBPmin=rep(FBP$minimos,N),FBPmax=rep(FBP$maximos,N),
              FBP2min=rep(FBP2$minimos,N),FBP2max=rep(FBP2$maximos,N))



Rates=rbind(data.frame(grupo=ER$grupo,Positive=ER$valor >= ER$Qmin & ER$valor <= ER$Qmax,band="FPB",stringsAsFactors = FALSE),
            data.frame(grupo=ER$grupo,Positive=ER$valor >= ER$Cmin & ER$valor <= ER$Cmax,band="TPB",stringsAsFactors = FALSE),
            data.frame(grupo=ER$grupo,Positive=ER$valor >= ER$FBPmin & ER$valor <= ER$FBPmax,band="DBB",stringsAsFactors = FALSE),
            data.frame(grupo=ER$grupo,Positive=ER$valor >= ER$FBP2min & ER$valor <= ER$FBP2max,band="DBB2s",stringsAsFactors = FALSE))


Rates$band <- factor(Rates$band, levels = c("FPB","TPB","DBB","DBB2s"))


ERates <- Rates %>%
  group_by(band) %>%
  summarise(
    TP = sum(grupo == "PB" & Positive == TRUE),      # inband and predicted inband
    FN = sum(grupo == "PB" & Positive == FALSE),     # inband and predicted out-of-band
    FP = sum(grupo == "out.PB" & Positive == TRUE),  # out-of-band and predicted inband
    TN = sum(grupo == "out.PB" & Positive == FALSE), # out-of-band and predicted out-of-band
    .groups = "drop"
  ) %>%
  mutate(
    Sensitivity = TP / (TP + FN),
    FNR = FN / (FN + TP),   
    Specificity = TN / (TN + FP),
    FPR = FP / (FP + TN)  
  )

SPREADS=c(SP2fpb,SP2tpb,SP2dbb,SP2dbb2)
SK=c(SKQ,SKT,NA,NA)


Precision=ERates$TP/(ERates$TP+ERates$FP)
Recall=ERates$TP/(ERates$TP+ERates$FN)
Beta=(SP2/SPREADS)^2
ERates$BetaScore=((1+Beta)*Precision*Recall)/(Beta*Precision+Recall)
ERates$BetaScore2=(2*Precision*Recall)/(Precision+Recall)


PrecisionN=ERates$TN/(ERates$TN+ERates$FN)
RecallN=ERates$TN/(ERates$TN+ERates$FP)
Beta=(SPREADS/SP2)^2
ERates$BetaScoreN=((1+Beta)*PrecisionN*RecallN)/(Beta*PrecisionN+RecallN)
ERates$BetaScore2N=(2*PrecisionN*RecallN)/(PrecisionN+RecallN)

ERates$BS=(1-alpha)*ERates$BetaScore+alpha*ERates$BetaScoreN

ERates$BS2=(ERates$BetaScore+ERates$BetaScoreN)/2


ERates$spM=SPREADS
ERates$sp=SP2

ERates$skewness=SK
ERates$SKG=SKG


Beta=(SP2/SPREADS)^2
ERates$FscoreAJ=((1+Beta)*Precision*PrecisionN)/(Beta*Precision+PrecisionN)

ERates$Precision=Precision
ERates$PrecisionN=PrecisionN

ERates$Recall=Recall
ERates$RecallN=RecallN

  
DF[[i]]=data.frame(ERates,Alpha=alpha)

  
  BANDS=rbind(Cuant$minimos,Cuant$maximos,Cons$minimos,Cons$maximos,FBP$minimos,FBP$maximos,FBP2$minimos,FBP2$maximos)
  
  SB[[i]]=FuSDfn(BANDS,Fdata)
  
  ## Functional Bands
  if (alpha == 0.5) {
    
    # Functional Bands (aligned with Quantiles)
    # Expected in the environment:
    # df_long: t_index, valor, curva, grupo in {"PB","out.PB"}
    # Cons   : (x, minimos, maximos)   -> TPB
    # FBP    : (x, minimos, maximos)   -> DBB
    # FBP2   : (x, minimos, maximos)   -> DBB2s (two-stage)
    if (!exists("xlabel")) xlabel <- "t"
    if (!exists("ylabel")) ylabel <- "LST"
    
    # Final palette (paper-friendly)
    col_methods <- c(
      "TPB"   = "#E4572E",  # vermillion (hero)
      "DBB"   = "#C51BB6",  # magenta
      "DBB2s" = "#4DAF4A"   # green (two-stage)
    )
    
    # Ensure numeric x
    df_long$t_index <- as.numeric(df_long$t_index)
    Cons$x  <- as.numeric(Cons$x)
    FBP$x   <- as.numeric(FBP$x)
    FBP2$x  <- as.numeric(FBP2$x)
    
    # X range and labels
    xf <- sort(unique(c(df_long$t_index, Cons$x, FBP$x, FBP2$x)))
    if (!exists("xrangelabs")) xrangelabs <- xf
    
    # Y range
    if (!exists("y_limits")) y_limits <- c(9.5, 60)
    if (!exists("y_breaks")) y_breaks <- seq(10, 60, 10)
    
    # Fairness: same visual weight across methods
    method_lw    <- 1.40
    method_alpha <- 0.80
    
    p_bands <- ggplot() +
    # draw spaghetti first to keep data visible
    geom_line(
      data = subset(df_long, grupo == "out.PB"),
      aes(x = t_index, y = valor, group = curva, color = grupo),
      linewidth = 0.45, alpha = 0.25, lineend = "round"
    ) +
      # central curves (white halo + charcoal line)
      geom_line(
        data = subset(df_long, grupo == "PB"),
        aes(x = t_index, y = valor, group = curva),
        color = "white", linewidth = 0.85, alpha = 0.45, lineend = "round"
      ) +
      geom_line(
        data = subset(df_long, grupo == "PB"),
        aes(x = t_index, y = valor, group = curva, color = grupo),
        linewidth = 0.75, alpha = 0.55, lineend = "round"
      ) +
      scale_color_manual(
        name   = "Curves",
        values = c("PB" = "#5A646E", "out.PB" = "#D8DDE7"),
        breaks = c("PB", "out.PB"),
        labels = c("Central", "Non-central")
      ) +
      guides(colour = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1))) +
      
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale("linetype") +
    # DBB (magenta, dotted) — same lw/alpha
    geom_line(
      data = transform(FBP, metodo = "DBB"),
      aes(x = x, y = maximos, color = metodo, linetype = metodo),
      linewidth = method_lw, alpha = method_alpha, lineend = "round"
    ) +
      geom_line(
        data = transform(FBP, metodo = "DBB"),
        aes(x = x, y = minimos, color = metodo, linetype = metodo),
        linewidth = method_lw, alpha = method_alpha, lineend = "round"
      ) +
      
      # DBB2s (green, twodash) — same lw/alpha (no halo, keeps fairness)
      geom_line(
        data = transform(FBP2, metodo = "DBB2s"),
        aes(x = x, y = maximos, color = metodo, linetype = metodo),
        linewidth = method_lw, alpha = method_alpha, lineend = "round"
      ) +
      geom_line(
        data = transform(FBP2, metodo = "DBB2s"),
        aes(x = x, y = minimos, color = metodo, linetype = metodo),
        linewidth = method_lw, alpha = method_alpha, lineend = "round"
      ) +
      
      # TPB (vermillion, solid) — small white halo + line on top
      geom_line(
        data = Cons, aes(x = x, y = maximos),
        inherit.aes = FALSE, color = "white", linewidth = 2.2, alpha = 0.95,
        linetype = "solid", lineend = "round", show.legend = FALSE
      ) +
      geom_line(
        data = Cons, aes(x = x, y = minimos),
        inherit.aes = FALSE, color = "white", linewidth = 2.2, alpha = 0.95,
        linetype = "solid", lineend = "round", show.legend = FALSE
      ) +
      geom_line(
        data = transform(Cons, metodo = "TPB"),
        aes(x = x, y = maximos, color = metodo, linetype = metodo),
        linewidth = method_lw, alpha = 1.0, lineend = "round"
      ) +
      geom_line(
        data = transform(Cons, metodo = "TPB"),
        aes(x = x, y = minimos, color = metodo, linetype = metodo),
        linewidth = method_lw, alpha = 1.0, lineend = "round"
      ) +
      
      # Bands scales
      scale_color_manual(
        name   = "Bands",
        values = col_methods,
        breaks = c("TPB", "DBB", "DBB2s")
      ) +
      scale_linetype_manual(
        name   = "Bands",
        values = c("TPB" = "solid", "DBB" = "dotted", "DBB2s" = "twodash"),
        breaks = c("TPB","DBB","DBB2s")
      ) +
      scale_x_continuous(breaks = xf, labels = xrangelabs, expand = c(0, 0)) +
      scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) +
      labs(x = xlabel, y = ylabel) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position    = "right",
        legend.box         = "vertical",
        legend.title       = element_text(face = "bold", size = 14),
        legend.text        = element_text(size = 11, margin = margin(t = 0, r = 6, b = 0, l = 8)),
        legend.key.width   = grid::unit(14, "pt"),
        legend.key.height  = grid::unit(12, "pt"),
        legend.spacing.x   = grid::unit(6,  "pt"),
        legend.spacing.y   = grid::unit(4,  "pt"),
        legend.box.spacing = grid::unit(6,  "pt"),
        legend.margin      = margin(2, 2, 2, 2),
        axis.title.x       = element_text(margin = margin(t = 6)),
        axis.title.y       = element_text(margin = margin(r = 6)),
        panel.grid.minor   = element_blank(),
        panel.grid.major   = element_line(linewidth = 0.25, colour = "#F6F7F9"),
        panel.border       = element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
        plot.background    = element_rect(fill = "white", colour = NA),
        plot.margin        = margin(8, 4, 8, 8)
      )
    
    # p_bands
    # ggsave("Bands_2024.png", p_bands, width=6, height=5.5, units="in", dpi=400, bg="white")
    # ggsave("Bands_2024.pdf", p_bands, width=6, height=5.5, units="in", device=cairo_pdf)
  }
  
  ## Functional Skewness Plots
  if (alpha == 0.2) {
    thdprop <- 0.5
    # --- Labels / axes like Bands ---
    if (!exists("xlabel"))   xlabel   <- "t"
    if (!exists("ylabel"))   ylabel   <- "LST"
    if (!exists("y_limits")) y_limits <- c(9.5, 60)
    if (!exists("y_breaks")) y_breaks <- seq(10, 60, 10)
    
    # --- Common x (must match everything) ---
    Cons$x <- as.numeric(Cons$x)
    tt  <- Cons$x
    Tn  <- length(tt)
    
    # --- Median aligned to x (12 pts) ---
    med <- as.numeric(md[[1]])
    if (length(med) != Tn) med <- approx(x = seq_along(med), y = med, xout = seq_len(Tn), rule = 2)$y
    
    # --- TPB ↔ median envelope ---
    LeftInf  <- pmin(Cons$minimos, med)
    RightSup <- pmax(Cons$maximos, med)
    df.skewPLOT <- data.frame(tt = tt, LeftInf = LeftInf, RightSup = RightSup, median = med)
    
    # --- Threshold: fraction of domain outside the envelope (default 50%) ---
   
    domTHD <- ceiling(ncol(Fdata) * thdprop)
    
    # --- Indices (exact masks you use) ---
    indneg <- as.numeric(which(IndNEG == TRUE))
    indpos <- as.numeric(which(IndPOS == TRUE))
    
    # --- Find full curves outside ≥ domTHD points (negative / positive) ---
    SkewNeg <- integer(0)
    if (length(indneg)) {
      FEv <- Fdata[indneg, , drop = FALSE]
      INF <- matrix(rep(LeftInf,  length(indneg)), nrow = length(indneg), byrow = TRUE)
      SUP <- matrix(rep(RightSup, length(indneg)), nrow = length(indneg), byrow = TRUE)
      SkewNeg <- indneg[ which(apply(!((INF <= FEv) & (FEv <= SUP)), 1, sum) >= domTHD) ]
    }
    SkewPos <- integer(0)
    if (length(indpos)) {
      FEv <- Fdata[indpos, , drop = FALSE]
      INF <- matrix(rep(LeftInf,  length(indpos)), nrow = length(indpos), byrow = TRUE)
      SUP <- matrix(rep(RightSup, length(indpos)), nrow = length(indpos), byrow = TRUE)
      SkewPos <- indpos[ which(apply(!((INF <= FEv) & (FEv <= SUP)), 1, sum) >= domTHD) ]
    }
    
    col_low_fill <- "#0FA27A"     # lower→median
    col_up_fill  <- "#2C6DF3"     # median→upper
    # Bordes de los ribbons (sólidos, discretos)
    edge_low <-  "#0B7D5E" 
    edge_up  <-   "#1F4ED6" 
    # Median
    col_med  <-   "firebrick" 
    #col_out_low <-  "#79D6C9" #"#79D6C9"    # teal claro
    col_out_low <-"#00C853"
    col_out_up  <-   "navy"     # blue claro
    
    pal <- scales::gradient_n_pal(c("blue","paleturquoise1","red"),
                                  values = c(0, 0.5, 1))
    
    # Ribbons (LEFT=cool, RIGHT=warm)
    col_low_fill <- pal(0.25)  # left region (azul claro)
    edge_low     <- pal(0.35)  # borde left, un poco más intenso
    
    col_up_fill  <- pal(0.80)  # right region (tirando a naranja)
    edge_up      <- pal(0.90)  # borde right, más intenso
    
    # Median
    col_med <- pal(0.50)       # turquesa (el del medio de tu escala)
    
    # Outs 
    col_out_low <- pal(0)
    col_out_up  <- pal(1)

    r_alpha    <- 0.95   # ribbons casi opacos
    edge_lw    <- 0.9   # grosor bordes
    edge_alpha <- 1   # bordes discretos
    cur_lw     <- 1.1   # curvas fuera
    cur_alpha  <- 1   # fuera discretas
    out_lty    <- "43"   #  (opción: "13")
    
    
    
    
    # --- Base plot: ribbons + límites + mediana ---
    p_skew_ribbons <- ggplot() +
      # Ribbons (fills fijos)
      geom_ribbon(data = df.skewPLOT, aes(tt, ymin = LeftInf,  ymax = median),
                  fill = col_low_fill, alpha = r_alpha, show.legend = FALSE) +
      geom_ribbon(data = df.skewPLOT, aes(tt, ymin = median, ymax = RightSup),
                  fill = col_up_fill,  alpha = r_alpha, show.legend = FALSE) +
      # TPB limits (sólidos, discretos)
      geom_line(data = df.skewPLOT, aes(tt, LeftInf),
                colour = edge_low, linewidth = edge_lw, linetype = "solid", alpha = edge_alpha) +
      geom_line(data = df.skewPLOT, aes(tt, RightSup),
                colour = edge_up,  linewidth = edge_lw, linetype = "solid", alpha = edge_alpha) +
      # halo
      geom_line(data = df.skewPLOT, aes(tt, median),
                colour = "#2B3035", linewidth = 2.1, alpha = 0.85,
                lineend = "round", linejoin = "round") +
      # mediana encima
      geom_line(data = df.skewPLOT, aes(tt, median),
                colour = col_med, linewidth = 1.55,
                lineend = "round", linejoin = "round")
    

    
    # --- Flagged Curves ---
    for (i in seq_along(SkewNeg)) {
      df_agr <- data.frame(tt = tt, st = as.numeric(Fdata[SkewNeg[i], ]))
      p_skew_ribbons <- p_skew_ribbons +
        geom_path(data = df_agr, aes(tt, st), inherit.aes = FALSE,
                  colour = col_out_low, linewidth = cur_lw, alpha = cur_alpha,
                  lineend = "round", linetype = out_lty)
    }
    for (i in seq_along(SkewPos)) {
      df_agr <- data.frame(tt = tt, st = as.numeric(Fdata[SkewPos[i], ]))
      p_skew_ribbons <- p_skew_ribbons +
        geom_path(data = df_agr, aes(tt, st), inherit.aes = FALSE,
                  colour = col_out_up, linewidth = cur_lw, alpha = cur_alpha,
                  lineend = "round", linetype = out_lty)
    }
    
    # --- Axes / theme  ---
    p_skew_ribbons <- p_skew_ribbons +
      scale_x_continuous(breaks = tt, labels = if (exists("xrangelabs")) xrangelabs else tt, expand = c(0, 0)) +
      coord_cartesian(clip = "off")+
      scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) +
      labs(x = xlabel, y = ylabel) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position   = "none",
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(linewidth = 0.25, colour = "#F1F3F5"),
        panel.border      = element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
        plot.background   = element_rect(fill = "white", colour = NA),
        plot.margin       = margin(8, 18, 8, 8)
      )
    
    p_skew_ribbons.OUT=p_skew_ribbons
  }
  
 
  if (alpha == 0.2) {
    thdprop <- 1
    # --- Labels / axes like Bands ---
    if (!exists("xlabel"))   xlabel   <- "t"
    if (!exists("ylabel"))   ylabel   <- "LST"
    if (!exists("y_limits")) y_limits <- c(9.5, 60)
    if (!exists("y_breaks")) y_breaks <- seq(10, 60, 10)
    
    # --- Common x (must match everything) ---
    Cons$x <- as.numeric(Cons$x)
    tt  <- Cons$x
    Tn  <- length(tt)
    
    # --- Median aligned to x (12 pts) ---
    med <- as.numeric(md[[1]])
    if (length(med) != Tn) med <- approx(x = seq_along(med), y = med, xout = seq_len(Tn), rule = 2)$y
    
    # --- TPB ↔ median envelope ---
    LeftInf  <- pmin(Cons$minimos, med)
    RightSup <- pmax(Cons$maximos, med)
    df.skewPLOT <- data.frame(tt = tt, LeftInf = LeftInf, RightSup = RightSup, median = med)
    
    # --- Threshold: fraction of domain outside the envelope (default 50%) ---
    
    domTHD <- ceiling(ncol(Fdata) * thdprop)
    
    # --- Indices (exact masks you use) ---
    indneg <- as.numeric(which(IndNEG == TRUE))
    indpos <- as.numeric(which(IndPOS == TRUE))
    
    # --- Find full curves outside ≥ domTHD points (negative / positive) ---
    SkewNeg <- integer(0)
    if (length(indneg)) {
      FEv <- Fdata[indneg, , drop = FALSE]
      INF <- matrix(rep(LeftInf,  length(indneg)), nrow = length(indneg), byrow = TRUE)
      SUP <- matrix(rep(RightSup, length(indneg)), nrow = length(indneg), byrow = TRUE)
      SkewNeg <- indneg[ which(apply(!((INF <= FEv) & (FEv <= SUP)), 1, sum) >= domTHD) ]
    }
    SkewPos <- integer(0)
    if (length(indpos)) {
      FEv <- Fdata[indpos, , drop = FALSE]
      INF <- matrix(rep(LeftInf,  length(indpos)), nrow = length(indpos), byrow = TRUE)
      SUP <- matrix(rep(RightSup, length(indpos)), nrow = length(indpos), byrow = TRUE)
      SkewPos <- indpos[ which(apply(!((INF <= FEv) & (FEv <= SUP)), 1, sum) >= domTHD) ]
    }
    

    col_low_fill <- "#0FA27A"     # lower→median
    col_up_fill  <- "#2C6DF3"     # median→upper
    # Ribbons (sólidos, discretos)
    edge_low <-  "#0B7D5E" 
    edge_up  <-   "#1F4ED6" 
    # Mediana
    col_med  <-   "firebrick" 
    #col_out_low <-  "#79D6C9" #"#79D6C9"    # teal claro
    col_out_low <-  "#00C853" #"#79D6C9"    # teal claro
    col_out_up  <-   "navy"     # blue claro
    

    pal <- scales::gradient_n_pal(c("blue","paleturquoise1","red"),
                                  values = c(0, 0.5, 1))
    
    # Ribbons (LEFT=cool, RIGHT=warm)
    col_low_fill <- pal(0.25)  # left region (azul claro)
    edge_low     <- pal(0.35)  # borde left, un poco más intenso
    
    col_up_fill  <- pal(0.80)  # right region (tirando a naranja)
    edge_up      <- pal(0.90)  # borde right, más intenso
    
    # Median
    col_med <- pal(0.50)       # turquesa (el del medio de tu escala)
    
    # Outs 
    col_out_low <- pal(0.98)
    col_out_up  <- pal(0.98)

    r_alpha    <- 0.95   # ribbons casi opacos
    edge_lw    <- 0.9   # grosor bordes
    edge_alpha <- 1   # bordes discretos
    cur_lw     <- 0.8   # curvas fuera
    cur_alpha  <- 1   # fuera discretas
    out_lty    <- "13"   #  (opción: "13")
    
    # --- Base plot: ribbons + límites + mediana ---
    p_skew_ribbons <- ggplot() +
      # Ribbons (fills fijos)
      geom_ribbon(data = df.skewPLOT, aes(tt, ymin = LeftInf,  ymax = median),
                  fill = col_low_fill, alpha = r_alpha, show.legend = FALSE) +
      geom_ribbon(data = df.skewPLOT, aes(tt, ymin = median, ymax = RightSup),
                  fill = col_up_fill,  alpha = r_alpha, show.legend = FALSE) +
      # TPB limits (sólidos, discretos)
      geom_line(data = df.skewPLOT, aes(tt, LeftInf),
                colour = edge_low, linewidth = edge_lw, linetype = "solid", alpha = edge_alpha) +
      geom_line(data = df.skewPLOT, aes(tt, RightSup),
                colour = edge_up,  linewidth = edge_lw, linetype = "solid", alpha = edge_alpha) +
      # halo
      geom_line(data = df.skewPLOT, aes(tt, median),
                colour = "#2B3035", linewidth = 2.1, alpha = 0.85,
                lineend = "round", linejoin = "round") +
      # mediana encima
      geom_line(data = df.skewPLOT, aes(tt, median),
                colour = col_med, linewidth = 1.55,
                lineend = "round", linejoin = "round")
      

    
    # --- Flagged Curves ---
    for (i in seq_along(SkewNeg)) {
      df_agr <- data.frame(tt = tt, st = as.numeric(Fdata[SkewNeg[i], ]))
      p_skew_ribbons <- p_skew_ribbons +
        geom_path(data = df_agr, aes(tt, st), inherit.aes = FALSE,
                  colour = col_out_low, linewidth = cur_lw, alpha = cur_alpha,
                  lineend = "round", linetype = out_lty)
    }
    for (i in seq_along(SkewPos)) {
      df_agr <- data.frame(tt = tt, st = as.numeric(Fdata[SkewPos[i], ]))
      p_skew_ribbons <- p_skew_ribbons +
        geom_path(data = df_agr, aes(tt, st), inherit.aes = FALSE,
                  colour = col_out_up, linewidth = cur_lw, alpha = cur_alpha,
                  lineend = "round", linetype = out_lty)
    }
    
    # --- Axes / theme  ---
    p_skew_ribbons <- p_skew_ribbons +
      scale_x_continuous(breaks = tt, labels = if (exists("xrangelabs")) xrangelabs else tt, expand =c(0,0) ) +
      coord_cartesian(clip = "off")+
      scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) +
      labs(x = xlabel, y = ylabel) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position   = "none",
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(linewidth = 0.25, colour = "#F1F3F5"),
        panel.border      = element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
        plot.background   = element_rect(fill = "white", colour = NA),
        plot.margin       = margin(8, 18, 8, 8)
      )
    
    
  }
  

}


##### Align both
adj <- match_panel_width_right(p_bands, p_quant)
p_bands_eq <- adj[[1]]
p_quant_eq <- adj[[2]]


##### Axis_tick_boost
# --- 1) Remove ALL legends from Bands (both "Curves" and "Bands") ---
p_bands_eq <- p_bands_eq + theme(legend.position = "none")

# --- 2) Make axis tick labels larger + bold (same tweak for both plots) ---
axis_tick_boost <- theme(
  axis.text.x  = element_text(size = 14, face = "bold"),  # tick labels
  axis.text.y  = element_text(size = 14, face = "bold"),
  axis.title.x = element_text(size = 16,face = "bold"),              # axis titles (optional)
  axis.title.y = element_text(size = 16,face = "bold")
)

p_bands_eq  <- p_bands_eq  + labs(y = NULL)+axis_tick_boost
p_quant_eq  <- p_quant_eq  +  ylab(expression(bold(LST)~(degree*C)))+axis_tick_boost



nameQ=paste0("FQuantiles",max(per/12)+2017)

## Export Quantile Plots and Probability Bands ####
ggsave(paste0(nameQ, ".png"), p_quant_eq, width = 8, height = 5, units = "in", dpi = 400, bg = "white")
ggsave(paste0(nameQ, ".pdf"), p_quant_eq, width = 8, height = 5, units = "in", device = function(filename, ...) grDevices::pdf(filename, useDingbats = FALSE, ...))

nameBand=paste0("B05and",max(per/12)+2017)

ggsave(paste0(nameBand, ".png"), p_bands_eq, width = 6.38, height = 5, units = "in", dpi = 400, bg = "white")
ggsave(paste0(nameBand, ".pdf"), p_bands_eq, width = 6.38, height = 5, units = "in", device = cairo_pdf)



## Export Functional Skewness plots ####

p_skew_ribbons.OUT=p_skew_ribbons.OUT+labs(y = NULL)+ axis_tick_boost

p_skew_ribbons=p_skew_ribbons+ ylab(expression(bold(LST)~(degree*C)))+axis_tick_boost

nameSKEW=paste0("skew",max(per/12)+2017)

ggsave(paste0(nameSKEW,"out", ".png"), p_skew_ribbons.OUT, width = 6.38, height = 5, units = "in", dpi = 400, bg = "white")
ggsave(paste0(nameSKEW,"out", ".pdf"), p_skew_ribbons.OUT, width = 6.38, height = 5, units = "in", device = cairo_pdf)
ggsave(paste0(nameSKEW, ".png"), p_skew_ribbons, width = 6.775, height = 5, units = "in", dpi = 400, bg = "white")
ggsave(paste0(nameSKEW, ".pdf"), p_skew_ribbons, width = 6.775, height = 5, units = "in", device = cairo_pdf)


## Save Functional Skewness plots

## Spread, FBG, Cov
axis_tick_boost <- theme(
  axis.text.x  = element_text(size = 14, face = "bold"),
  axis.text.y  = element_text(size = 14, face = "bold"),
  axis.title.x = element_text(size = 20, face = "bold"),
  axis.title.y = element_text(size = 20, face = "bold")
)

## Spread Comparison ####
# Build input tables 
# DF: list of data.frames con columnas: Alpha, spM, band, sp
dfER <- dplyr::bind_rows(DF) |> dplyr::mutate(gamma = 1 - Alpha)

# Methods and references
spread_tbl <- dplyr::transmute(dfER, gamma, spread = spM, band)
csp_tbl    <- dfER |> dplyr::distinct(gamma, csp = sp) |> dplyr::arrange(gamma)

keep_methods <- c("FPB", "TPB", "DBB")
df_methods <- spread_tbl |>
  dplyr::filter(band %in% keep_methods) |>
  dplyr::mutate(band = factor(band, levels = keep_methods)) |>
  dplyr::arrange(band, gamma)

# Plot 
p_spread <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted",
             color = "#9AA3AD", linewidth = 1.5, alpha = 1, show.legend = FALSE) +
  # referencia: halo + línea twodash
  geom_line(data = csp_tbl, aes(gamma, csp),
            color = "white", linewidth = 2.2, linetype = "twodash",
            alpha = 0.55, lineend = "round", show.legend = FALSE) +
  geom_line(data = csp_tbl, aes(gamma, csp, color = "CSP1"),
            linewidth = 1.6, linetype = "twodash", alpha = 0.95, lineend = "round") +
  # métodos (sólidos)
  geom_line(data = df_methods, aes(gamma, spread, color = band),
            linewidth = 1.6, alpha = 0.95, lineend = "round") +
  geom_vline(xintercept = 0.50, linetype = "dotted",
             color = "#9AA3AD", linewidth = 1.5, alpha = 1, show.legend = FALSE) +
  scale_x_continuous(limits = c(0.05, 0.95),
                     breaks = seq(0.1,0.9,0.1),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 61),
                     breaks = seq(0, 60, 10),
                     labels = as.character(seq(0, 60, 10)),
                     expand = c(0, 0)) +
  scale_color_manual(
    name   = "Method",
    values = c("FPB" = "#1F63FF", "TPB" = "#E4572E", "DBB" = "#C51BB6", "CSP1" = "#2B3035"),
    breaks = c("FPB", "TPB", "DBB", "CSP1"),
    labels = c("FPB", "TPB", "DBB", expression(CSP^{(1)}))
  ) +
  guides(color = guide_legend(
    nrow = 1, title.position = "left",
    override.aes = list(linetype = c("solid","solid","solid","twodash"),
                        linewidth = c(1.8, 1.8, 1.8, 1.6),
                        lineend = "round", alpha = 1)
  )) +
  labs(x = expression(gamma), y = "Spread") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.box.just   = "center",
    legend.title      = element_text(face = "bold"),
    legend.text       = element_text(size = 10),
    legend.spacing.x  = grid::unit(28, "pt"),
    legend.spacing.y  = grid::unit(0,  "pt"),
    legend.key.width  = grid::unit(16, "pt"),
    legend.margin     = margin(t = 2, r = 0, b = 2, l = 0),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.25, colour = "#F6F7F9"),
    panel.border      = element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
    plot.background   = element_rect(fill = "white", colour = NA),
    plot.margin       = margin(8, 10, 8, 8)
  )

## Homologation (single source of truth)
ylims       <- c(-1, 60)
ybreaks     <- seq(0, 60, 10)
top_pad     <- 26
leg_key_h   <- grid::unit(12, "pt")
leg_key_w   <- grid::unit(20, "pt")
leg_box_gap <- grid::unit(4,  "pt")
legend_pad  <- margin(t = 0, r = 6, b = 0, l = 8)

legend_theme <- theme(
  legend.position    = "top",
  legend.box         = "horizontal",
  legend.box.just    = "center",
  legend.title       = element_text(face = "bold", size = 12),
  legend.text        = element_text(size = 11, margin = legend_pad),
  legend.key.width   = leg_key_w,
  legend.key.height  = leg_key_h,
  legend.spacing.x   = grid::unit(28, "pt"),
  legend.spacing.y   = grid::unit(0,  "pt"),
  legend.box.spacing = leg_box_gap,
  legend.margin      = margin(0, 0, 0, 0)
)

panel_theme <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(linewidth = 0.25, colour = "#F2F2F2"),
  panel.border     = element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
  plot.background  = element_rect(fill = "white", colour = NA),
  plot.margin      = margin(t = top_pad, r = 10, b = 8, l = 8)
)


p_spread <- p_spread +
  scale_y_continuous(limits = ylims, breaks = ybreaks, expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  legend_theme + panel_theme

p_spread <- p_spread + theme(legend.position = "none") + axis_tick_boost

nameSP <- paste0("spread_gamma", max(per/12) + 2017)

ggplot2::ggsave(
  filename = paste0(nameSP, ".png"),
  plot     = p_spread,
  width    = 6, height = 5.5, units = "in",
  dpi      = 400, bg = "white"
)
ggplot2::ggsave(
  filename = paste0(nameSP, ".pdf"),
  plot     = p_spread,
  width    = 6, height = 5.5, units = "in",
  device   = function(filename, ...) grDevices::pdf(filename, useDingbats = FALSE, ...),
  bg = "white"
)

## FBetaGamma ####

keep_methods <- c("TPB","DBB")

df_fbg <- dfER |>
  dplyr::filter(band %in% keep_methods) |>
  dplyr::mutate(band = factor(band, levels = keep_methods)) |>
  dplyr::arrange(band, gamma) |>
  dplyr::select(gamma, Fbg = BS, band)

# Shared styling
ylims   <- c(-0.01, 1.00)
ybreaks <- seq(0, 1, 0.1)
top_pad     <- 26
leg_key_h   <- grid::unit(12, "pt")
leg_key_w   <- grid::unit(20, "pt")
leg_box_gap <- grid::unit(4,  "pt")
legend_pad  <- margin(t = 0, r = 6, b = 0, l = 8)

legend_theme <- theme(
  legend.position    = "top",
  legend.box         = "horizontal",
  legend.box.just    = "center",
  legend.title       = element_text(face = "bold", size = 12),
  legend.text        = element_text(size = 11, margin = legend_pad),
  legend.key.width   = leg_key_w,
  legend.key.height  = leg_key_h,
  legend.spacing.x   = grid::unit(28, "pt"),
  legend.spacing.y   = grid::unit(0,  "pt"),
  legend.box.spacing = leg_box_gap,
  legend.margin      = margin(0, 0, 0, 0)
)

panel_theme <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(linewidth = 0.25, colour = "#F2F2F2"),
  panel.border     = element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
  plot.background  = element_rect(fill = "white", colour = NA),
  plot.margin      = margin(t = top_pad, r = 10, b = 8, l = 8)
)

p_fbg <- ggplot(df_fbg, aes(x = gamma, y = Fbg, color = band)) +
  geom_hline(yintercept = 0, linetype = "dotted",
             color = "#9AA3AD", linewidth = 1.5, alpha = 1, show.legend = FALSE) +
  geom_vline(xintercept = 0.50, linetype = "dotted",
             color = "#9AA3AD", linewidth = 1.5, alpha = 1, show.legend = FALSE) +
  geom_line(linewidth = 1.5, alpha = 0.95, lineend = "round") +
  scale_x_continuous(limits = c(0.05, 0.95),
                     breaks = seq(0.1,0.9,0.1),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = ylims,
                     breaks  = ybreaks,
                     labels  = sprintf("%.1f", ybreaks),
                     expand  = c(0, 0)) +
  scale_color_manual(
    name   = "Method",
    values = c("TPB" = "#E4572E", "DBB" = "#C51BB6"),
    breaks = keep_methods
  ) +
  guides(color = guide_legend(
    nrow = 1, title.position = "left",
    override.aes = list(linetype = "solid", linewidth = 1.8, lineend = "round", alpha = 1)
  )) +
  labs(x = expression(gamma), y = expression(bold(F)[bold(beta)]^{bold(gamma)})) +
  theme_minimal(base_size = 12) +
  legend_theme + panel_theme

nameFbG <- paste0("FBG", max(per/12) + 2017)
p_fbg   <- p_fbg + theme(legend.position = "none") + axis_tick_boost

## DBB2s FbetaGamma
keep_methods <- c("DBB2s")

df_fbg <- dfER |>
  dplyr::filter(band %in% keep_methods) |>
  dplyr::mutate(band = factor(band, levels = keep_methods)) |>
  dplyr::arrange(band, gamma) |>
  dplyr::select(gamma, Fbg = BS, band)

v05 <- subset(df_fbg, gamma == 0.5)$Fbg

p_fbg <- p_fbg +
  annotate("point", x = 0.5, y = v05, shape = 23, size = 4.0,
           fill = "#4DAF4A", color = "white", stroke = 1.2)

ggplot2::ggsave(paste0(nameFbG, ".png"), p_fbg, width = 6, height = 5.5, units = "in", dpi = 400, bg = "white")
ggplot2::ggsave(paste0(nameFbG, ".pdf"), p_fbg, width = 6, height = 5.5, units = "in",
                device = function(filename, ...) grDevices::pdf(filename, useDingbats = FALSE, ...), bg = "white")


#### FBetaGamma Comparison ####

keep_methods <- c("TPB","DBB","DBB2s")

df_fbg <- dfER |>
  dplyr::filter(band %in% keep_methods) |>
  dplyr::mutate(band = factor(band, levels = keep_methods)) |>
  dplyr::arrange(band, gamma) |>
  dplyr::select(gamma, Fbg = BS, band)

FBGpaper=subset(df_fbg, gamma == 0.5)



#### F_{beta}^{gamma} components (paper-style facet titles) ####

# - Curves: TPB vs DBB
# - Marker: DBB2s (green diamond) at gamma = 0.5, included in the legend
# - Facets: Rec^C_gamma, Rec^N_gamma, Prec^C_gamma, Prec^N_gamma (no "(B)")
# - Styling consistent with the Fbg figure, with slight y headroom to avoid top clipping

# Methods
# Line methods to compare
keep_methods <- c("TPB", "DBB")

# Legend order (lines + marker method)
keep_methods_all <- c("TPB", "DBB", "DBB2s")

# Facet titles in plotmath (parsed): (Rec^C)_gamma etc.
# Parentheses enforce the correct precedence: superscript first, then subscript gamma.
comp_labs <- c(
  "Recall"     = "Rec^{C}",
  "RecallN"    = "Rec^{N}",
  "Precision"  = "Prec^{C}",
  "PrecisionN" = "Prec^{N}"
)

# Data preparation (long format)
df_comp <- dfER |>
  dplyr::filter(band %in% keep_methods) |>
  dplyr::mutate(band = factor(band, levels = keep_methods_all)) |>
  dplyr::select(gamma, band, Recall, RecallN, Precision, PrecisionN) |>
  tidyr::pivot_longer(
    cols      = c(Recall, RecallN, Precision, PrecisionN),
    names_to  = "component",
    values_to = "value"
  ) |>
  # Keep the raw component keys and let the facet labeller do the math parsing
  dplyr::mutate(component = factor(component, levels = names(comp_labs)))

# DBB2s marker at gamma = 0.5 (one point per facet/component)
df_point <- dfER |>
  dplyr::filter(band == "DBB2s", dplyr::near(gamma, 0.50)) |>
  dplyr::mutate(band = factor("DBB2s", levels = keep_methods_all)) |>
  dplyr::select(gamma, band, Recall, RecallN, Precision, PrecisionN) |>
  tidyr::pivot_longer(
    cols      = c(Recall, RecallN, Precision, PrecisionN),
    names_to  = "component",
    values_to = "value"
  ) |>
  dplyr::mutate(component = factor(component, levels = names(comp_labs)))

# Shared styling (matched to your Fbg plot)
ylims   <- c(-0.01, 1.05)   # small headroom above 1.0 to prevent top clipping
ybreaks <- seq(0, 1, 0.1)   # ticks shown up to 1.0

top_pad     <- 26
leg_key_h   <- grid::unit(12, "pt")
leg_key_w   <- grid::unit(20, "pt")
leg_box_gap <- grid::unit(4,  "pt")
legend_pad  <- margin(t = 0, r = 6, b = 0, l = 8)

legend_theme <- theme(
  legend.position    = "top",
  legend.box         = "horizontal",
  legend.box.just    = "center",
  legend.title       = element_text(face = "bold", size = 12),
  legend.text        = element_text(size = 11, margin = legend_pad),
  legend.key.width   = leg_key_w,
  legend.key.height  = leg_key_h,
  legend.spacing.x   = grid::unit(28, "pt"),
  legend.spacing.y   = grid::unit(0,  "pt"),
  legend.box.spacing = leg_box_gap,
  legend.margin      = margin(0, 0, 0, 0)
)

panel_theme <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(linewidth = 0.25, colour = "#F2F2F2"),
  panel.border     = element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
  plot.background  = element_rect(fill = "white", colour = NA),
  plot.margin      = margin(t = top_pad, r = 10, b = 8, l = 8),
  strip.background = element_blank(),
  strip.text       = element_text(face = "bold", size = 12)
)

axis_tick_boost <- theme(
  axis.ticks        = element_line(linewidth = 0.6, colour = "#2B3035"),
  axis.ticks.length = grid::unit(3, "pt")
)

# Plot
p_comp <- ggplot(df_comp, aes(x = gamma, y = value, color = band)) +
  # Reference lines (as in the Fbg style)
  geom_hline(
    yintercept = 0,
    linetype   = "dotted",
    color      = "#9AA3AD",
    linewidth  = 1.5,
    alpha      = 1,
    show.legend = FALSE
  ) +
  geom_vline(
    xintercept = 0.50,
    linetype   = "dotted",
    color      = "#9AA3AD",
    linewidth  = 1.5,
    alpha      = 1,
    show.legend = FALSE
  ) +
  # Method curves
  geom_line(linewidth = 1.5, alpha = 0.95, lineend = "round") +
  # Facets with parsed math titles
  facet_wrap(
    ~ component,
    ncol = 2,
    labeller = as_labeller(comp_labs, label_parsed)
  ) +
  # Axes
  scale_x_continuous(
    limits = c(0.05, 0.95),
    breaks = seq(0.1, 0.9, 0.1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = ylims,
    breaks = ybreaks,
    labels = sprintf("%.1f", ybreaks),
    expand = c(0, 0)
  ) +
  # Legend colors (DBB2s color is used for the point + legend key)
  scale_color_manual(
    name   = "Method",
    values = c("TPB" = "#E4572E", "DBB" = "#C51BB6", "DBB2s" = "#4DAF4A"),
    breaks = keep_methods_all
  ) +
  # Legend: show lines for TPB/DBB, and a green diamond key for DBB2s
  guides(color = guide_legend(
    nrow = 1,
    title.position = "left",
    override.aes = list(
      linetype = c("solid", "solid", "blank"),
      linewidth = c(1.8, 1.8, 0),
      shape    = c(NA, NA, 23),
      fill     = c(NA, NA, "#4DAF4A"),
      alpha    = c(1, 1, 1)
    )
  )) +
  labs(x = expression(gamma), y = NULL) +
  theme_minimal(base_size = 12) +
  legend_theme + panel_theme + axis_tick_boost

# DBB2s marker (green diamond) at gamma = 0.5 in each facet
# A slightly smaller size keeps the marker visually balanced.
if (nrow(df_point) > 0) {
  p_comp <- p_comp +
    geom_point(
      data = df_point,
      inherit.aes = FALSE,
      aes(x = gamma, y = value, color = band),
      shape  = 23,
      size   = 3.2,
      fill   = "#4DAF4A",
      stroke = 0.8,
      show.legend = TRUE
    )
}

# Save
nameComp <- "FBG_components"

ggsave(
  filename = paste0(nameComp, ".png"),
  plot     = p_comp,
  width    = 12,
  height   = 5.5,
  units    = "in",
  dpi      = 400,
  bg       = "white"
)

ggsave(
  filename = paste0(nameComp, ".pdf"),
  plot     = p_comp,
  width    = 12,
  height   = 5.5,
  units    = "in",
  device   = function(filename, ...) grDevices::pdf(filename, useDingbats = FALSE, ...),
  bg       = "white"
)




## FuSD coverage ####
i <- 1
dfFSD <- data.frame(
  band     = c("FPB","FPB","TPB","TPB","DBB","DBB","DBB2s","DBB2s"),
  quantile = c("min","max","min","max","min","max","min","max"),
  FuSD     = SB[[i]][1:8],
  gamma    = 1 - A[i]
)
for (i in 2:length(A)) {
  dfFSD <- rbind(
    dfFSD,
    data.frame(
      band     = c("FPB","FPB","TPB","TPB","DBB","DBB","DBB2s","DBB2s"),
      quantile = c("min","max","min","max","min","max","min","max"),
      FuSD     = SB[[i]][1:8],
      gamma    = 1 - A[i]
    )
  )
}

# palette
col_methods <- c("TPB" = "#E4572E", "DBB" = "#C51BB6")
col_target  <- "#2B3035"  # referencia unificada

keep_methods <- c("TPB","DBB")

df_bands <- dfFSD |>
  dplyr::filter(band %in% keep_methods) |>
  dplyr::mutate(band = factor(band, levels = keep_methods)) |>
  dplyr::arrange(band, gamma, quantile)

gamma_vals <- sort(unique(df_bands$gamma))
A_vals     <- 1 - gamma_vals
q_lo <- as.numeric(stats::quantile(sttdVar, probs = A_vals/2,     names = FALSE))
q_hi <- as.numeric(stats::quantile(sttdVar, probs = 1 - A_vals/2, names = FALSE))

df_target <- dplyr::bind_rows(
  tibble::tibble(band = "Target", quantile = "min", FuSD = q_lo, gamma = gamma_vals),
  tibble::tibble(band = "Target", quantile = "max", FuSD = q_hi, gamma = gamma_vals)
)

lvls_all <- c(keep_methods, "Target")

p_cov <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted",
             color = "#9AA3AD", linewidth = 1.2, show.legend = FALSE) +
  geom_vline(xintercept = 0.50, linetype = "dotted",
             color = "#9AA3AD", linewidth = 1.2, show.legend = FALSE) +
  # Target (gris twodash con halo fino)
  geom_line(data = df_target,
            aes(gamma, FuSD, group = interaction(band, quantile)),
            color = "white", linewidth = 2.2, alpha = 0.55, show.legend = FALSE) +
  geom_line(data = df_target,
            aes(gamma, FuSD, color = "Target", linetype = "Target",
                group = interaction(band, quantile)),
            linewidth = 1.6, lineend = "round") +
  # Methods (solid)
  geom_line(data = df_bands,
            aes(gamma, FuSD, color = band, linetype = band,
                group = interaction(band, quantile)),
            linewidth = 1.6, alpha = 0.95, lineend = "round") +
  scale_x_continuous(limits = c(0.05, 0.95),
                     breaks = seq(0.1, 0.9, 0.1),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.51, 0.50),
                     breaks = seq(-0.5, 0.5, 0.1),
                     labels = scales::number_format(accuracy = 0.1),
                     expand = c(0, 0)) +
  scale_color_manual(
    name   = "Method",
    values = c(col_methods, "Target" = col_target),
    breaks = lvls_all,
    labels = c("TPB", "DBB", "FuSD coverage")
  ) +
  scale_linetype_manual(
    name   = "Method",
    values = c("TPB" = "solid", "DBB" = "solid", "Target" = "twodash"),
    breaks = lvls_all,
    labels = c("TPB", "DBB", "FuSD coverage")
  ) +
  guides(
    color   = guide_legend(nrow = 1, title.position = "left",
                           override.aes = list(alpha = 1, linewidth = 1.8, lineend = "round")),
    linetype= guide_legend(nrow = 1, title.position = "left")
  ) +
  labs(x = expression(gamma), y = "FuSD") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.box.just   = "center",
    legend.title      = element_text(face = "bold", size = 12),
    legend.text       = element_text(size = 11, margin = margin(t = 0, r = 6, b = 0, l = 8)),
    legend.key.width  = grid::unit(20, "pt"),
    legend.key.height = grid::unit(12, "pt"),
    legend.spacing.x  = grid::unit(28, "pt"),
    legend.spacing.y  = grid::unit(0,  "pt"),
    legend.box.spacing= grid::unit(4,  "pt"),
    legend.margin     = margin(0, 0, 0, 0),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.25, colour = "#F2F2F2"),
    panel.border      = element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
    plot.background   = element_rect(fill = "white", colour = NA),
    plot.margin       = margin(t = 26, r = 10, b = 8, l = 8)
  )

nameFuSDcov <- paste0("FuSDcov", max(per/12) + 2017)
p_cov <- p_cov + theme(legend.position = "none") + axis_tick_boost

# DBB2s @ gamma = 0.5 (min y max)
dbb2s_all  <- subset(dfFSD, band == "DBB2s")
g_half     <- dbb2s_all$gamma[ which.min(abs(dbb2s_all$gamma - 0.5)) ]
v_dbb2s_lo <- dbb2s_all$FuSD[ dbb2s_all$quantile == "min" & dbb2s_all$gamma == g_half ]
v_dbb2s_hi <- dbb2s_all$FuSD[ dbb2s_all$quantile == "max" & dbb2s_all$gamma == g_half ]

x_nudge <- 0.004

p_cov <- p_cov +
  annotate("point", x = g_half + x_nudge, y = v_dbb2s_hi,
           shape = 23, size = 4.0, fill = "#4DAF4A", color = "white", stroke = 1.2) +
  annotate("point", x = g_half + x_nudge, y = v_dbb2s_lo,
           shape = 23, size = 4.0, fill = "#4DAF4A", color = "white", stroke = 1.2)

ggsave(paste0(nameFuSDcov, ".png"), p_cov, width = 6, height = 5.5, dpi = 300, bg = "white")
ggsave(paste0(nameFuSDcov, ".pdf"), p_cov, width = 6, height = 5.5,
       device = function(filename, ...) grDevices::pdf(filename, useDingbats = FALSE, ...), bg = "white")

