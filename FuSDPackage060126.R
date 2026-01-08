## FuSD Package
## Date: 060126
## Author: Emmanuel Ambriz

u_transform <- function(FdataT) {
  t(apply(FdataT, 1, function(row) {
    Fn <- ecdf(row)        # ecdf row
    Fn(row)                # ecdf
  })) |> as.data.frame()
}
FuSD=function(b){
  Fdata=as.matrix(b)
  FdataT=t(Fdata)
  
  fb=fbplot(FdataT,
            method='MBD',
            ylim=c(min(FdataT),
                   max(FdataT)),
            xlab="",ylab="",xaxt = "n",plot=F)
  
  signos=c()
  deepest=as.numeric(fb$medcurve)[1]
  for(i in 1:dim(FdataT)[2])
  {
    sig=sum(FdataT[,i]-FdataT[,deepest])>0
    if(!sig){
      sig=-1
    }
    signos[i]=sig*1    
  }
  
  
  (fb$depth[deepest]-fb$depth)*signos
}
FuSDU=function(b){
  Fdata=as.matrix(b)
  FdataT=t(Fdata)
  
  fb=fbplot(FdataT,
            method='MBD',
            ylim=c(min(FdataT),
                   max(FdataT)),
            xlab="",ylab="",xaxt = "n",plot=F)
  
  signos=c()
  deepest=as.numeric(fb$medcurve)[1]
  
  FdataU=u_transform(FdataT) 
  
  for(i in 1:dim(FdataT)[2])
  {
    sig=sum(FdataU[,i]-FdataU[,deepest])>0
    if(!sig){
      sig=-1
    }
    signos[i]=sig*1    
  }
  
  
  (fb$depth[deepest]-fb$depth)*signos
}
FuSDfn=function(xSample,fndata=xSample){
     o2=as.numeric(MBD(fndata,plotting=F)$MBD)
     deepest=which(o2==max(o2))[1]
     o1=as.numeric(MBD(x=xSample,xRef=fndata, plotting=F)$MBD)
     signos=c()
     for(i in 1:dim(xSample)[1])
     {
          sig=sum(xSample[i,]-fndata[deepest,])>0
          if(!sig){
               sig=-1
          }
          signos[i]=sig*1    
     }
     (max(o2)-o1)*signos
}
FuSDw=function(b){
  Fdata=as.matrix(b)
  FdataT=t(Fdata)
  
  fb=fbplot(FdataT,
            method='MBD',
            ylim=c(min(FdataT),
                   max(FdataT)),
            xlab="",ylab="",xaxt = "n",plot=F)
  
  signos=c()
  deepest=as.numeric(fb$medcurve)[1]
  
  SD=apply(FdataT,1,sd)
  w=SD/sum(SD)
  
  for(i in 1:dim(FdataT)[2])
  {
    sig=sum((FdataT[,i]-FdataT[,deepest])*w)>0
    if(!sig){
      sig=-1
    }
    signos[i]=sig*1    
  }
  
  
  #list((fb$depth[deepest]-fb$depth)*signos,w)
  (fb$depth[deepest]-fb$depth)*signos
}

FuSDw2=function(b){
  Fdata=as.matrix(b)
  FdataT=t(Fdata)
  
  fb=fbplot(FdataT,
            method='MBD',
            ylim=c(min(FdataT),
                   max(FdataT)),
            xlab="",ylab="",xaxt = "n",plot=F)
  
  signos=c()
  deepest=as.numeric(fb$medcurve)[1]
  
  SD=apply(FdataT,1,sd)
  w=(1/SD)/sum(1/SD)
  
  for(i in 1:dim(FdataT)[2])
  {
    sig=sum((FdataT[,i]-FdataT[,deepest])*w)>0
    if(!sig){
      sig=-1
    }
    signos[i]=sig*1    
  }
  
  
  #list((fb$depth[deepest]-fb$depth)*signos,w)
   (fb$depth[deepest]-fb$depth)*signos
}


FunctionalSTtDQuantile=function(alpha,data,STtDfn)
{
     sttdVar=STtDfn(data)
     IndEval=which(ecdf(sttdVar)(sttdVar)>=alpha)
     jEval=which(sttdVar[IndEval]==min(sttdVar[IndEval]))
     Arg=IndEval[jEval]
     #sttdVar[Arg]
     list(as.vector(data[Arg,]),Arg)
}





PBsttd=function(fndata,alpha,xf,STtDfn){
     alpha1=alpha/2
     alpha2=1-alpha/2
     sttdVar=STtDfn(fndata)
     qs=quantile(sttdVar,c(alpha1,alpha2))
     IndVal=which(((sttdVar>=qs[1]) & (sttdVar<=qs[2]))==T)
     minimos=as.numeric(apply(fndata[IndVal,],2,min))
     maximos=as.numeric(apply(fndata[IndVal,],2,max))

     XQ=FunctionalSTtDQuantile(.5,fndata,STtDfn)[[1]]
          
     Cuant=data.frame(x=xf,minimos,mediana=XQ,maximos)
     
     gIQR=ggplot(Cuant,aes(x=x,y=mediana))+
          geom_ribbon(aes(ymin=minimos,ymax=maximos,fill="skyblue"),col="blue")+
          geom_line(size=2)+
          geom_line(aes(x=x,y=maximos),col="blue",size=2)+
          geom_line(aes(x=x,y=minimos),col="blue",size=2)+
          scale_fill_manual(values="skyblue")+
          theme_minimal()+
          theme(legend.position="none",
                legend.title=element_text(size = 15,face="bold"),
                axis.text.x= element_text(size = 12,face = "bold"),
                axis.text.y= element_text(size = 12,face = "bold"),
                legend.key.size = unit(.9, 'cm'), 
                legend.key.height = unit(.9, 'cm'), 
                legend.key.width = unit(.9, 'cm'), legend.text = element_text(size=14,face = "bold"),
                axis.title=element_text(size = 15,face = "bold"),plot.title = element_text(size=20,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.5))+
          labs(y="",x="t")
     
     
     Sp=sqrt(sum((minimos-maximos)^2))
     
     list(Cuant,gIQR,Sp)
     
}

PBadjusted=function(alpha,data,STtDfn,SP,df_longAlpha)
{
  sttdVar=STtDfn(data)
  qeval=quantile(sttdVar,c(alpha/2,1-alpha/2))
  Pos=which(sttdVar>0)
  Neg=which(sttdVar<0)
  sttdPos=sttdVar[Pos]
  sttdNeg=sttdVar[Neg]
  
  Ii=min(sttdVar[Pos])
  S=max(sttdVar[Neg])
  
  # warning if the contention not holds
  
  IndEvalN=as.numeric(Neg[which(sttdNeg>=qeval[1])])
  IndEvalP=as.numeric(Pos[which(sttdPos<=qeval[2])])
  
  BS=c()
  In=c()
  Ip=c()
  for(i in 1:length(IndEvalN))
  {
    for(j in 1:length(IndEvalP))
    {
      PBelements=rbind(data[IndEvalN[i],],data[IndEvalP[j],])
      Cuant=data.frame(minimos=apply(PBelements,2,min),maximos=apply(PBelements,2,max))
      SP2tpb=sqrt(sum((Cuant$minimos-Cuant$maximos)^2))
      ER=data.frame(df_longAlpha,Qmin=rep(Cuant$minimos,N),Qmax=rep(Cuant$maximos,N))
      Rates=data.frame(ER$grupo,Positive=ER$valor >= ER$Qmin & ER$valor <= ER$Qmax)
      
      ERates <- Rates %>%
        summarise(
          TP = sum(ER.grupo == "PB" & Positive == TRUE),      # dentro y predicho dentro
          FN = sum(ER.grupo == "PB" & Positive == FALSE),     # dentro pero predicho fuera
          FP = sum(ER.grupo == "out.PB" & Positive == TRUE),  # fuera pero predicho dentro
          TN = sum(ER.grupo == "out.PB" & Positive == FALSE), # fuera y predicho fuera
          .groups = "drop"
        ) %>%
        mutate(
          Sensitivity = TP / (TP + FN),
          FNR = FN / (FN + TP),   # tasa de falsos negativos # red curves
          Specificity = TN / (TN + FP),
          FPR = FP / (FP + TN)  # tasa de falsos positivos #  blue curves
        )
      
      
      Precision=ERates$TP/(ERates$TP+ERates$FP)
      Recall=ERates$TP/(ERates$TP+ERates$FN)
      Beta=(SP/SP2tpb)^2
      ERates$BetaScore=((1+Beta)*Precision*Recall)/(Beta*Precision+Recall)
      
      PrecisionN=ERates$TN/(ERates$TN+ERates$FN)
      RecallN=ERates$TN/(ERates$TN+ERates$FP)
      Beta=(SP2tpb/SP)^2
      ERates$BetaScoreN=((1+Beta)*PrecisionN*RecallN)/(Beta*PrecisionN+RecallN)
      BS=c(BS,(1-alpha)*ERates$BetaScore+alpha*ERates$BetaScoreN)
      In=c(In,i)
      Ip=c(Ip,j)
    }
  }
  
  c(IndEvalN[In[which(BS==max(BS))]],IndEvalP[Ip[which(BS==max(BS))]],max(BS))
}
plot_functional_quantiles <- function(
    quant_tbl,
    Fdata = NULL,
    xlabel = "t",
    ylabel = "LST",
    xrangelabs = NULL,
    # Colors
    col_inner  = "#1F63FF",   # Q0.25 & Q0.75 (blue)
    col_outer  = "#8F1537",   # Q0.125 & Q0.875 (claret)
    col_median = "#111111",   # Q0.50 (black)
    # Line widths and halos
    lwd_q          = 1.8,     # main line width
    lwd_med        = 1.9,     # median slightly stronger
    halo_w         = 1.2,     # white under-stroke width (all)
    halo_alpha     = 0.70,    # white under-stroke opacity
    # Blue-only gradient/glow (to keep red visible underneath at early x)
    halo_blue_w        = 1.4,
    halo_blue_alpha0   = 0.25,  # start alpha for blue glow
    halo_blue_alpha1   = 0.25,  # end alpha for blue glow
    blue_line_alpha0   = 0.80,  # start alpha for blue line
    blue_line_alpha1   = 0.8,  # end alpha for blue line
    # Background spaghetti
    show_background    = !is.null(Fdata),
    bg_linewidth       = 0.9,
    bg_alpha_outer     = 0.18,
    bg_alpha_central   = 0.34,
    bg_color_light     = "#D8DDE7",
    bg_color_mid       = "#9BA4AE",
    # Axes
    y_headroom         = 1.03,   # 3% headroom
    show_y0line        = TRUE
) {
  # --- Basic checks ---
  req <- c("x","Q0.125","Q0.25","Q0.50","Q0.75","Q0.875")
  stopifnot(all(req %in% names(quant_tbl)))
  n <- nrow(quant_tbl)
  if (is.null(xrangelabs)) xrangelabs <- quant_tbl$x
  stopifnot(length(xrangelabs) == n)
  
  # --- Build long data (use internal index pos = 1..n) ---
  pos <- seq_len(n)
  q_med <- data.frame(pos = pos, y = quant_tbl$`Q0.50`,  qname = "Q0.50",  tier = "Q0.50")
  q_in  <- rbind(
    data.frame(pos = pos, y = quant_tbl$`Q0.25`, qname = "Q0.25"),
    data.frame(pos = pos, y = quant_tbl$`Q0.75`, qname = "Q0.75")
  )
  q_in$tier <- "Q0.25 & Q0.75"
  q_out <- rbind(
    data.frame(pos = pos, y = quant_tbl$`Q0.125`, qname = "Q0.125"),
    data.frame(pos = pos, y = quant_tbl$`Q0.875`, qname = "Q0.875")
  )
  q_out$tier <- "Q0.125 & Q0.875"
  
  # Gradient weights along x: map to [0,1]
  w <- if (n > 1) (pos - min(pos)) / (max(pos) - min(pos)) else rep(1, n)
  blue_glow_alpha <- halo_blue_alpha0 + (halo_blue_alpha1 - halo_blue_alpha0) * w
  blue_line_alpha <- blue_line_alpha0 + (blue_line_alpha1 - blue_line_alpha0) * w
  
  # Attach gradient alphas to each blue series independently
  q_in <- q_in |>
    dplyr::arrange(qname, pos) |>
    dplyr::group_by(qname) |>
    dplyr::mutate(a_glow = blue_glow_alpha, a_line = blue_line_alpha) |>
    dplyr::ungroup()
  
  # Background to long if provided
  df_bg_long <- NULL
  if (show_background) {
    stopifnot(ncol(Fdata) == n)
    df_bg <- as.data.frame(Fdata)
    colnames(df_bg) <- as.character(seq_len(ncol(df_bg)))
    df_bg$curve_id <- seq_len(nrow(df_bg))
    df_bg_long <- tidyr::pivot_longer(df_bg, -curve_id, names_to = "pos", values_to = "y")
    df_bg_long$pos <- as.integer(df_bg_long$pos)
  }
  
  # y upper bound
  y_upper <- max(
    q_med$y, q_in$y, q_out$y,
    if (show_background) max(df_bg_long$y, na.rm = TRUE) else -Inf,
    na.rm = TRUE
  ) * y_headroom
  y_upper <- max(y_upper, 1)
  
  # --- Build plot (paint order: red -> blue -> median) ---
  p <- ggplot2::ggplot()
  
  # Background spaghetti (two passes for soft depth)
  if (show_background) {
    p <- p +
      ggplot2::geom_line(
        data = df_bg_long, ggplot2::aes(pos, y, group = curve_id),
        color = bg_color_light, linewidth = bg_linewidth, alpha = bg_alpha_outer,
        lineend = "round"
      ) +
      ggplot2::geom_line(
        data = df_bg_long, ggplot2::aes(pos, y, group = curve_id),
        color = bg_color_mid, linewidth = bg_linewidth, alpha = bg_alpha_central,
        lineend = "round"
      )
  }
  
  # Baseline at y = 0
  if (show_y0line) {
    p <- p + ggplot2::geom_hline(yintercept = 0,  linetype = "dotted",
                                 color = "#9AA3AD", linewidth = 0.7, alpha = 1, show.legend = FALSE)
  }
  
  # 1) Extremes (red): white halo then colored line
  p <- p +
    ggplot2::geom_line(
      data = q_out, ggplot2::aes(pos, y, group = qname),
      color = "white", linewidth = halo_w, alpha = halo_alpha, lineend = "round"
    ) +
    ggplot2::geom_line(
      data = q_out, ggplot2::aes(pos, y, color = tier, linetype = tier, group = qname),
      linewidth = lwd_q, alpha = 0.98, lineend = "round"
    )
  
  # 2) Quartiles (blue): white halo, blue glow with gradient alpha, then blue line with gradient alpha
  p <- p +
    ggplot2::geom_line(
      data = q_in, ggplot2::aes(pos, y, group = qname),
      color = "white", linewidth = halo_w, alpha = halo_alpha, lineend = "round"
    ) +
    ggplot2::geom_line(
      data = q_in, ggplot2::aes(pos, y, group = qname, alpha = a_glow),
      color = col_inner, linewidth = halo_blue_w, lineend = "round", show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = q_in, ggplot2::aes(pos, y, color = tier, linetype = tier, group = qname, alpha = a_line),
      linewidth = lwd_q, lineend = "round"
    ) +
    ggplot2::scale_alpha_identity(guide = "none")
  
  # 3) Median (black): white halo then dashed black on top
  p <- p +
    ggplot2::geom_line(
      data = q_med, ggplot2::aes(pos, y, group = qname),
      color = "white", linewidth = halo_w + 0.1, alpha = halo_alpha, lineend = "round"
    ) +
    ggplot2::geom_line(
      data = q_med, ggplot2::aes(pos, y, color = tier, linetype = tier, group = qname),
      linewidth = lwd_med, alpha = 1.0, lineend = "round"
    )
  
  # Legend and scales (order: Median → Quartiles → Extremes)
  p <- p +
    ggplot2::scale_color_manual(
      name   = "Quantiles",
      values = c("Q0.50" = col_median,
                 "Q0.25 & Q0.75" = col_inner,
                 "Q0.125 & Q0.875" = col_outer),
      breaks = c("Q0.50","Q0.25 & Q0.75","Q0.125 & Q0.875"),
      labels = c(expression(Q[0.5]),
                 expression(Q[0.25]~"&"~Q[0.75]),
                 expression(Q[0.125]~"&"~Q[0.875]))
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Q0.50" = "longdash", "Q0.25 & Q0.75" = "solid", "Q0.125 & Q0.875" = "solid"),
      breaks = c("Q0.50","Q0.25 & Q0.75","Q0.125 & Q0.875"),
      guide  = "none"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        order = 1, nrow = 1, title.position = "left",
        override.aes = list(
          linewidth = 1.8,
          linetype  = c("longdash","solid","solid"),
          alpha     = 1
        )
      )
    ) +
    ggplot2::scale_x_continuous(breaks = pos, labels = xrangelabs, expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(0, y_upper), expand = c(0,0)) +
    ggplot2::labs(x = xlabel, y = ylabel) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position   = "top",
      legend.box        = "horizontal",
      legend.box.just   = "center",
      legend.title      = ggplot2::element_text(face = "bold"),
      legend.text       = ggplot2::element_text(size = 10),
      legend.spacing.x  = grid::unit(8, "pt"),
      legend.key.width  = grid::unit(18, "pt"),
      legend.margin     = ggplot2::margin(t = 2, r = 0, b = 2, l = 0),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.grid.major  = ggplot2::element_line(linewidth = 0.25, colour = "#F2F2F2"),
      panel.border      = ggplot2::element_rect(colour = "#D5DAE1", fill = NA, linewidth = 0.6),
      plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin       = ggplot2::margin(8, 10, 8, 8)
    )
  
  return(p)
}


