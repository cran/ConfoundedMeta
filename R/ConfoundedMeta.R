
#' Estimates and inference for sensitivity analyses
#'
#' Computes point estimates, standard errors, and confidence interval bounds
#' for (1) \code{prop}, the proportion of studies with true effect sizes above \code{.q} (or below
#' \code{.q} for an apparently preventive \code{.yr}) as a function of the bias parameters;
#' (2) the minimum bias factor on the relative risk scale (\code{Tmin}) required to reduce to
#' less than \code{.r} the proportion of studies with true effect sizes more extreme than
#' \code{.q}; and (3) the counterpart to (2) in which bias is parameterized as the minimum
#' relative risk for both confounding associations (\code{Gmin}).
#' @param .q True effect size that is the threshold for "scientific significance"
#' @param .r For \code{Tmin} and \code{Gmin}, value to which the proportion of large effect sizes is to be reduced
#' @param .muB Mean bias factor on the log scale across studies
#' @param .sigB Standard deviation of log bias factor across studies
#' @param .yr Pooled point estimate (on log scale) from confounded meta-analysis
#' @param .vyr Estimated variance of pooled point estimate from confounded meta-analysis
#' @param .t2 Estimated heterogeneity (tau^2) from confounded meta-analysis
#' @param .vt2 Estimated variance of tau^2 from confounded meta-analysis
#' @param CI.level Confidence level as a proportion
#' @export
#' @details
#' To compute all three point estimates (\code{prop, Tmin, and Gmin}) and inference, all
#' arguments must be non-\code{NULL}. To compute only a point estimate for \code{prop},
#' arguments \code{.r, .vyr}, and \code{.vt2} can be left \code{NULL}. To compute only
#' point estimates for \code{Tmin} and \code{Gmin}, arguments \code{.muB, .vyr}, and \code{.vt2}
#' can be left \code{NULL}. To compute inference for all point estimates, \code{.vyr} and 
#' \code{.vt2} must be supplied. 
#' @import metafor
#' stats 
#' @examples
#' d = metafor::escalc(measure="RR", ai=tpos, bi=tneg,
#' ci=cpos, di=cneg, data=metafor::dat.bcg)
#' 
#' m = metafor::rma.uni(yi= d$yi, vi=d$vi, knha=FALSE,
#'                      measure="RR", method="DL" ) 
#' yr = as.numeric(m$b)  # metafor returns on log scale
#' vyr = as.numeric(m$vb)
#' t2 = m$tau2
#' vt2 = m$se.tau2^2 
#' 
#' # obtaining all three estimators and inference
#' confounded_meta( .q=log(0.90), .r=0.20, .muB=log(1.5), .sigB=0.1,
#'                  .yr=yr, .vyr=vyr, .t2=t2, .vt2=vt2,
#'                  CI.level=0.95 )
#' 
#' # passing only arguments needed for prop point estimate
#' confounded_meta( .q=log(0.90), .muB=log(1.5),
#'                  .yr=yr, .t2=t2, CI.level=0.95 )
#' 
#' # passing only arguments needed for Tmin, Gmin point estimates
#' confounded_meta( .q=log(0.90), .r=0.20,
#'                  .yr=yr, .t2=t2, CI.level=0.95 )


confounded_meta = function( .q, .r=NULL, .muB=NULL, .sigB=0,
                             .yr, .vyr=NULL, .t2, .vt2=NULL,
                            CI.level=0.95 ) {

  # somewhere have option to plot the bias factor distribution, the confounded distribution, and the adjusted distribution
  
  ##### Check for Bad Input #####
  if ( .t2 <= .sigB^2 ) stop("Must have .t2 > .sigB^2")
  if ( is.null(.vyr) | is.null(.vt2) ) warning("Cannot compute inference without .vyr and .vt2. Returning only point estimates.")
  if ( is.null(.r) ) warning("Cannot compute Tmin or Gmin without .r. Returning only prop.")

  ##### Point Estimates: Causative Case #####
  if ( .yr > log(1) ) {
    
    if ( !is.null(.muB) ) {
      # prop above
      Z = ( .q + .muB - .yr ) / sqrt( .t2 - .sigB^2 )
      phat = 1 - pnorm(Z) 
    } else {
      phat = NA
    }
    
    if ( !is.null(.r) ) {
      # min bias factor
      Tmin = exp( qnorm(1-.r) * sqrt(.t2) - .q + .yr )
      
      # min confounding strength
      Gmin = Tmin + sqrt( Tmin^2 - Tmin )
    } else {
      Tmin = Gmin = NA
    }
  }
  
  ##### Point Estimates: Preventive Case #####
  else if ( .yr < log(1) ) {
    
    if ( !is.null(.muB) ) {
      # prop above
      Z = ( .q - .muB - .yr ) / sqrt( .t2 - .sigB^2 )
      phat = pnorm(Z) 
    } else {
      phat = NA
    }
    
    if ( !is.null(.r) ) {
      # min bias factor
      Tmin = exp( .q - .yr - qnorm(.r) * sqrt(.t2) )
      
      # min confounding strength
      Gmin = Tmin + sqrt( Tmin^2 - Tmin )
    } else {
      Tmin = Gmin = NA
    }
  }

  
  ##### Delta Method Inference: P-Hat #####
  # do inference only if given needed SEs
  if ( !is.null(.vyr) & !is.null(.vt2) & !is.null(.muB) ){
    
    # term in numerator depends on whether causative or preventive RR
    num.term = ifelse( .yr > log(1), .q + .muB - .yr, .q - .muB - .yr )
    
    term1.1 = .vyr / (.t2 - .sigB^2 )
    term1.2 = ( .vt2 * (num.term)^2 ) / ( 4 * (.t2 - .sigB^2 )^3 )
    term1 = sqrt( term1.1 + term1.2 )
    
    Z = num.term / sqrt( .t2 - .sigB^2 )
    SE = term1 * dnorm(Z)
    
    # confidence interval
    tail.prob = ( 1 - CI.level ) / 2
    lo.phat = max( 0, phat + qnorm( tail.prob )*SE )
    hi.phat = min( 1, phat - qnorm( tail.prob )*SE )
    
    # TEST ONLY
    # compare to SE from R's DM function - same
    #library(msm)
    #deltamethod( ~ 1 - pnorm( ( q + eval(.muB) - x1 ) / sqrt( x2 - eval(.sigB)^2 ) ), mean = c( .yr, .t2 ), cov = diag( c(.vyr,.vt2) ) )
    #deltamethod( ~ 1 - pnorm( ( log(1.1) + log(1.7) - x1 ) / sqrt( x2 - .1^2 ) ), mean = c( .yr, .t2 ), cov = diag( c(.vyr,.vt2) ) )
    #deltamethod( ~ pnorm( ( log(0.90) - log(1.7) + x1 ) / sqrt( x2 - .1^2 ) ), mean = c( .yr, .t2 ), cov = diag( c(.vyr,.vt2) ) )
    
  } else {
    SE = lo.phat = hi.phat = NA
  }
  
  ##### Delta Method Inference: Tmin and Gmin #####
  # do inference only if given needed SEs and .r
  if ( !is.null(.vyr) & !is.null(.vt2) & !is.null(.r) ){
    
    ##### Tmin #####
    if (.yr > log(1) ) {
      term = ( .vt2 * qnorm(1-.r)^2 ) / ( 4 * .t2 )
      SE.T = exp( sqrt(.t2) * qnorm(1-.r) - .q + .yr ) * sqrt( .vyr + term  )
    } else {
      term = ( .vt2 * qnorm(.r)^2 ) / ( 4 * .t2 )
      SE.T = exp( .q - .yr - sqrt(.t2) * qnorm(.r) ) * sqrt( .vyr + term  )
    }
    
    # TEST ONLY
    # compare to SE from R's DM function
    #   deltamethod( ~ exp( qnorm(1-.r) * sqrt(x2) - .q + x1 ), mean = c( .yr, .t2 ), cov = diag( c(.vyr,.vt2) ) )
    # qnorm(1-.10) = 1.281552
    #   deltamethod( ~ exp( 1.281552 * sqrt(x2) - log(1.1) + x1 ), mean = c( .yr, .t2 ), cov = diag( c(.vyr,.vt2) ) )
    # qnorm(.90) = -1.281552
    #deltamethod( ~ exp( log(0.90) - x1 - (-1.281552) * sqrt(x2) ), mean = c( .yr, .t2 ), cov = diag( c(.vyr,.vt2) ) )  
    
    tail.prob = ( 1 - CI.level ) / 2
    lo.T = max( 1, Tmin + qnorm( tail.prob )*SE.T )  # bias factor can't be < 1
    hi.T = Tmin - qnorm( tail.prob )*SE.T  # but has no upper bound
    
    
    ##### Gmin #####
    SE.G = SE.T * ( 1 + ( 2*Tmin - 1 ) / ( 2 * sqrt( Tmin^2 - Tmin ) ) )
    
    # TEST ONLY
    # compare to SE from R's DM function
    #   deltamethod( ~ x1 + sqrt( x1^2 - x1 ),
    #                mean = c( Tmin ), 
    #                cov = SE.T^2 )
    
    lo.G = max( 1, Gmin + qnorm( tail.prob )*SE.G )  # confounding RR can't be < 1
    hi.G = Gmin - qnorm( tail.prob )*SE.G  # but has no upper bound
    
  } else {  # i.e., user didn't pass parameters needed for inference
    SE.T = SE.G = lo.T = lo.G = hi.T = hi.G = NA
  }
  
  
  # return results
  res = data.frame( Value = c("Prop", "Tmin", "Gmin"), 
                    Est = c( phat, Tmin, Gmin ),
                    SE = c(SE, SE.T, SE.G),
                    CI.lo = c(lo.phat, lo.T, lo.G), 
                    CI.hi = c(hi.phat, hi.T, hi.G) 
                    )
  
  return(res)
}







#' Tables for sensitivity analyses
#'
#' Produces table showing the proportion of true effect sizes more extreme than \code{.q}
#' across a grid of bias parameters \code{.muB} and \code{.sigB} (for \code{.meas == "prop"}).
#' Alternatively, produces a table showing the minimum bias factor (for \code{.meas == "Tmin"})
#' or confounding strength (for \code{.meas == "Gmin"}) required to reduce to less than
#' \code{.r} the proportion of true effects more extreme than \code{.q}.
#' @param .meas \code{prop}, \code{Tmin}, or \code{Gmin}
#' @param .q True effect size that is the threshold for "scientific significance"
#' @param .r For \code{Tmin} and \code{Gmin}, vector of values to which the proportion of large effect sizes is to be reduced
#' @param .muB Mean bias factor on the log scale across studies
#' @param .sigB Standard deviation of log bias factor across studies
#' @param .yr Pooled point estimate (on log scale) from confounded meta-analysis
#' @param .t2 Estimated heterogeneity (tau^2) from confounded meta-analysis
#' @keywords meta-analysis, confounding, sensitivity
#' @export
#' @details
#' For \code{.meas=="Tmin"} or \code{.meas=="Gmin"}, arguments \code{.muB} and
#' \code{.sigB} can be left \code{NULL}; \code{.r} can also be \code{NULL} as
#' it will default to a reasonable range of proportions. Returns a \code{data.frame}
#' whose rows are values of \code{.muB} (for \code{.meas=="prop"}) or of \code{.r} 
#' (for \code{.meas=="Tmin"} or \code{.meas=="Gmin"}). Its columns are values of 
#' \code{.sigB} (for \code{.meas=="prop"}) or of \code{.q} (for \code{.meas=="Tmin"}
#' or \code{.meas=="Gmin"}).
#' Tables for \code{Gmin} will display \code{NaN} for cells corresponding to \code{Tmin}<1,
#' i.e., for which no bias is required to reduce the effects as specified. 
#' 
#' @import ggplot2 
#' @examples
#' sens_table( .meas="prop", .q=log(1.1), .muB=c( log(1.1),
#' log(1.5), log(2.0) ), .sigB=c(0, 0.1, 0.2), 
#' .yr=log(2.5), .t2=0.1 )
#' 
#' sens_table( .meas="Tmin", .q=c( log(1.1), log(1.5) ),
#' .yr=log(1.3), .t2=0.1 ) 
#' 
#' # will have NaNs in cells with Tmin < 1 (no bias needed)
#' sens_table( .meas="Gmin", .r=0.8, .q=c( log(1.1) ),
#' .yr=log(1.3), .t2=0.1 )


sens_table = function( .meas, .q, .r=seq(0.1, 0.9, 0.1), .muB=NULL, .sigB=NULL,
                            .yr, .t2 ) {
  
  ##### Check for Correct Inputs Given Measure ######
  if ( .meas=="prop" & ( is.null(.muB) | is.null(.sigB) ) ) {
    stop( "To compute proportion above .q, provide .muB and .sigB")
  }
  
  if ( .meas=="prop" & length(.q) > 1 ) {
    stop( "To compute proportion above .q, provide only a single value of .q" )
  }
  
  ###### Generate Table #####
  
  # table skeleton
  nrow = ifelse( .meas=="prop", length(.muB), length(.r) )
  ncol = ifelse( .meas=="prop", length(.sigB), length(.q) )
  
  m = matrix( NA, nrow=nrow, ncol=ncol )
  
  # doing this inefficient thing because confounded_meta is not vectorized
  # because returns a dataframe
  for (i in 1:nrow) {
    for (j in 1:ncol) {
      if ( .meas == "prop" ) {
        m[i,j] = confounded_meta( .q=.q, .muB = .muB[i], .sigB = .sigB[j],
                                  .yr=.yr, .t2=.t2 )[1,"Est"]
      } else if ( .meas == "Tmin" ) {
        m[i,j] = confounded_meta( .q=.q[j], .r=.r[i],
                                  .yr=.yr, .t2=.t2 )[2,"Est"]
      } else if ( .meas == "Gmin" ) {
        m[i,j] = confounded_meta( .q=.q[j], .r=.r[i],
                                  .yr=.yr, .t2=.t2 )[3,"Est"]
      }
      
    }
  }
  
  d = data.frame(m)

  if ( .meas=="prop" ) {
    row.names(d) = round( .muB, 3 )
  } else if ( .meas %in% c( "Tmin", "Gmin" ) ) {
    row.names(d) = round( .r, 3 )
  }

  if( .meas=="prop" ) {
    names(d) = round( .sigB, 3 )
  } else if ( .meas %in% c( "Tmin", "Gmin" ) ) {
    names(d) = round( .q, 3 )
  }

  cat("Rows: ", ifelse( .meas=="prop", ".muB", ".r" ) )
  cat("\nColumns: ", ifelse( .meas=="prop", ".sigB", ".q" ) )
  cat("\n\n")
  
  return(d)
}




#' Plots for sensitivity analyses
#'
#' Produces line plots (\code{.type=="line"}) showing the bias factor on the relative risk (RR) scale vs. the proportion
#' of studies with true RRs above \code{.q} (or below it for an apparently preventive relative risk).
#' The plot secondarily includes a X-axis scaled based on the minimum strength of confounding
#' to produce the given bias factor. The shaded region represents a 95\% pointwise confidence band.
#' Alternatively, produces distribution plots (\code{.type=="dist"}) for a specific bias factor showing the observed and 
#' true distributions of RRs with a red line marking exp(\code{.q}).
#' @param .type \code{dist} for distribution plot; \code{line} for line plot (see Details)
#' @param .q True effect size that is the threshold for "scientific significance"
#' @param .muB Single mean bias factor on log scale (only needed for distribution plot)
#' @param .Bmin Lower limit of lower X-axis on the log scale (only needed for line plot)
#' @param .Bmax Upper limit of lower X-axis on the log scale (only needed for line plot)
#' @param .sigB Standard deviation of log bias factor across studies (length 1)
#' @param .yr Pooled point estimate (on log scale) from confounded meta-analysis
#' @param .vyr Estimated variance of pooled point estimate from confounded meta-analysis
#' @param .t2 Estimated heterogeneity (tau^2) from confounded meta-analysis
#' @param .vt2 Estimated variance of tau^2 from confounded meta-analysis
#' @param breaks.x1 Breaks for lower X-axis (bias factor) on RR scale
#' @param breaks.x2 Breaks for upper X-axis (confounding strength) on RR scale
#' @param CI.level Poitnwise confidence level as a proportion
#' @keywords meta-analysis, confounding, sensitivity
#' @details
#' Arguments \code{.vyr} and \code{.vt2} can be left \code{NULL}, in which case no confidence
#' band will appear on the line plot. 
#' @export
#' @import ggplot2 
#' @examples
#' # with variable bias and with confidence band
#' sens_plot( .type="line", .q=log(1.1), .Bmin=log(1), .Bmax=log(4), .sigB=0.1,
#'            .yr=log(1.3), .vyr=0.005, .t2=0.4, .vt2=0.03 )
#' 
#' # with fixed bias and without confidence band
#' sens_plot( .type="line", .q=log(1.1), .Bmin=log(1), .Bmax=log(4),
#'            .yr=log(1.3), .t2=0.4 )
#' 
#' # apparently preventive
#' sens_plot( .type="line", .q=log(0.90), .Bmin=log(1), .Bmax=log(4),
#'            .yr=log(0.6), .vyr=0.005, .t2=0.4, .vt2=0.04 )
#' 
#' # distribution plot: apparently causative
#' # commented out because takes 5-10 seconds to run
#' # sens_plot( .type="dist", .q=log(1.1), .muB=log(2),
#' #           .yr=log(1.3), .t2=0.4 )
#'            
#' # distribution plot: apparently preventive
#' # commented out because takes 5-10 seconds to run
#' # sens_plot( .type="dist", .q=log(0.90), .muB=log(1.5),
#' #           .yr=log(0.7), .t2=0.2 )


sens_plot = function( .type, .q, .muB=NULL, .Bmin=log(1), .Bmax=log(5), .sigB=0,
                      .yr, .vyr=NULL, .t2, .vt2=NULL,
                      breaks.x1=NULL, breaks.x2=NULL,
                      CI.level=0.95 ) {
  
  ##### Check for Bad Input ######
  if ( .type=="dist" ) {

    if( is.null(.muB) ) stop("For type=='dist', must provide .muB")
    
    if ( ( length(.muB) > 1 ) | ( length(.sigB) > 1 ) ) {
      stop( "For type=='dist', .muB and .sigB must be length 1")
    }
  }

  if ( .type=="line" ) {
    
    if ( is.null(.vyr) | is.null(.vt2) ) {
      warning( "No confidence interval because .vyr or .vt2 is NULL")
    }
  }
  
  ##### Distribution Plot ######
  if ( .type=="dist" ) {
    
    # simulate confounded distribution
    reps = 1000000
    RR.c = exp( rnorm( n=reps, mean=.yr, sd=sqrt(.t2) ) )
    
    # simulate unconfounded distribution
    Mt = ifelse( .yr > 0, .yr - .muB, .yr + .muB )
    RR.t = exp( rnorm( n=reps, mean=Mt, sd=sqrt(.t2-.sigB^2) ) )
    
    # get reasonable limits for X-axis
    x.min = min( quantile(RR.c, 0.01), quantile(RR.t, 0.01) )
    x.max = max( quantile(RR.c, 0.99), quantile(RR.t, 0.99) )
    
    temp = data.frame( group = rep( c( "Observed", "True" ), each = reps ), 
                       val = c( RR.c, RR.t ) )
    
    colors=c("black", "orange")
    p = ggplot2::ggplot( data=temp, aes(x=temp$val, group=temp$group ) ) +
      geom_density( aes( fill=temp$group ), alpha=0.4 ) +
      theme_bw() + xlab("Study-specific relative risks") +
      ylab("") + guides(fill=guide_legend(title=" ")) +
      scale_fill_manual(values=colors) +
      geom_vline( xintercept = exp(.q), lty=2, color="red" ) +
      scale_x_continuous( limits=c(x.min, x.max), breaks = seq( round(x.min), round(x.max), 0.5) ) +
      ggtitle("Observed and true relative risk distributions")

    graphics::plot(p)
  }
  
  ##### Line Plot ######
  if ( .type=="line" ) {
    # get mean bias factor values for a bunch of different B's
    t = data.frame( B = seq(.Bmin, .Bmax, .01), phat = NA, lo = NA, hi = NA )
    t$eB = exp(t$B)
    
    for ( i in 1:dim(t)[1] ) {
      # .r is irrelevant here
      cm = confounded_meta(.q, .r=0.10, .muB=t$B[i], .sigB,
                           .yr, .vyr, .t2, .vt2,
                           CI.level=CI.level)
      t$phat[i] = cm$Est[ cm$Value=="Prop" ]
      t$lo[i] = cm$CI.lo[ cm$Value=="Prop" ]
      t$hi[i] = cm$CI.hi[ cm$Value=="Prop" ]
    }

    # compute values of g for the dual X-axis
    if ( is.null(breaks.x1) ) breaks.x1 = seq( exp(.Bmin), exp(.Bmax), .5 )
    if ( is.null(breaks.x2) ) breaks.x2 = round( breaks.x1 + sqrt( breaks.x1^2 + breaks.x1 ), 2)
    
    p = ggplot2::ggplot( t, aes(x=t$eB, y=t$phat ) ) + theme_bw() +
      scale_y_continuous( limits=c(0,1), breaks=seq(0, 1, .1)) +
      scale_x_continuous(  breaks = breaks.x1,
                          sec.axis = sec_axis( ~ . + sqrt(.^2 + .),
                          name = "Minimum strength of both confounding RRs",
                          breaks=breaks.x2 ) ) +
      geom_line(lwd=1.2) +
      xlab("Bias factor (RR scale)") +
      ylab( paste( ifelse( .yr > log(1),
                           paste( "Estimated proportion of studies with true RR >", exp(.q) ),
                          paste( "Estimated proportion of studies with true RR <", exp(.q) ) ) ) )
  
    # can't compute a CI if the bounds aren't there
    no.CI = any( is.na(t$lo) ) | any( is.na(t$hi) )
    
    if ( no.CI ) graphics::plot(p)
    else p + ggplot2::geom_ribbon( aes(ymin=t$lo, ymax=t$hi), alpha=0.15 )   
  
  }
}



#' Convert forest plot or summary table to meta-analytic dataset
#'
#' Given relative risks (RR) and upper bounds of 95\% confidence intervals (CI)
#' from a forest plot or summary table, returns a dataframe ready for meta-analysis
#' (e.g., via the \code{metafor} package) with the log-RRs and their variances.
#' Optionally, the user may indicate studies for which the point estimate is to be
#' interpreted as an odds ratios of a common outcome rather than relative risks;
#' for such studies, the function applies VanderWeele (2017)'s square-root transformation to convert
#' the odds ratio to an approximate risk ratio. 
#' @param .est Vector of study point estimates on RR or OR scale
#' @param .hi Vector of upper bounds of 95\% CIs on RRs
#' @param .sqrt Vector of booleans (TRUE/FALSE) for whether each study measured an odds ratio of a common outcome that should be approximated as a risk ratio via the square-root transformation
#' @export
#' @import stats

scrape_meta = function( .est, .hi, .sqrt=FALSE ){
  
  # take square root for certain elements
  RR = .est
  RR[.sqrt] = sqrt( RR[.sqrt] )
  
  # same for upper CI limit
  hi.RR = .hi
  hi.RR[.sqrt] = sqrt( hi.RR[.sqrt] )
  
  sei = ( log(hi.RR) - log(RR) ) / qnorm(.975)
  
  return( data.frame( yi = log(RR), vyi = sei^2 ) )
}






