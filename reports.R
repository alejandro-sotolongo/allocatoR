tbl_mpt <- function(res, freq = "days") {
  if (inherits(res$b, "xts")) {
    is_bench <- TRUE
  } else {
    is_bench <- FALSE
  }
  a <- freq_to_scaler(freq)
  xcov <- cov(res$xb) * a
  if (nrow(res$xb) >= a) {
    geo_ret <- calc_geo_ret(res$xb, freq)
    rf_geo <- calc_geo_ret(res$rf, freq)
  } else {
    geo_ret <- apply(res$xb + 1, 2, prod) - 1
    rf_geo <- prod(1+res$rf) - 1
  }
  vol <- sqrt(diag(xcov))
  down_vol <- calc_down_vol(res$xb, freq)
  max_dd <- calc_max_drawdown(res$xb)
  sharpe <- calc_sharpe_ratio(res$xb, res$rf, freq)
  sortino <- calc_sortino_ratio(res$xb, res$rf, freq)
  recov <- geo_ret / -max_dd
  per_up <- apply(res$xb, 2, \(x){sum(x >= 0) / length(x)})
  per_down <- 1 - per_up
  avg_up <- apply(res$xb, 2, \(x){mean(x[x >= 0])})
  avg_down <- apply(res$xb, 2, \(x){mean(x[x <0])})
  res_out <- list()
  if (!is_bench) {
    x <- rbind(
      scales::percent(geo_ret, 0.01),
      scales::percent(vol, 0.01), 
      scales::percent(down_vol, 0.01), 
      scales::percent(max_dd, 0.01), 
      scales::number(sharpe, 0.01), 
      scales::number(sortino, 0.01), 
      scales::number(recov, 0.01),
      scales::percent(per_up, 0.01),
      scales::percent(per_down, 0.01),
      scales::percent(avg_up, 0.01),
      scales::percent(avg_down, 0.01),
      scales::number(avg_up / -avg_down, 0.01),
      scales::number((avg_up * per_up) / -(avg_down * per_down), 0.01)
    )
    num <- rbind(geo_ret, vol, down_vol, max_dd, sharpe, sortino, recov, 
      per_up, per_down, avg_up, avg_down, avg_up / -avg_down,
      (avg_up * per_up) / -(avg_down * per_down))
    xdf <- data.frame(
      Metric = c("Geometric Return", "Volatility", "Downside Vol",
                 "Worst Drawdown", "Sharpe Ratio", "Sortino Ratio", "Recovery",
                 "Up Periods", "Down Periods", "Average Up Period",
                 "Average Down Period", "Avg Up / Avg Down",
                 "Wgt Avg Up / Wgt Avg Down"),
      x,
      row.names = NULL
    )
    colnames(xdf) <- c("Metric", colnames(res$xb))
    numdf <- data.frame(
      Metric = c("Geometric Return", "Volatility", "Downside Vol",
                 "Worst Drawdown", "Sharpe Ratio", "Sortino Ratio", "Recovery",
                 "Up Periods", "Down Periods", "Average Up Period",
                 "Average Down Period", "Avg Up / Avg Down",
                 "Wgt Avg Up / Wgt Avg Down"),
      num,
      row.names = NULL
    )
    colnames(numdf) <- c("Metric", colnames(res$xb))
    res_out$xdf <- xdf
    res_out$numdf <- numdf
    return(res_out)
  } else {
    ar <- excess_ret(res$x, res$b)
    acov <- cov(ar) * a
    te <- c(sqrt(diag(acov)), NA)
    xbeta <- calc_uni_beta(res$x, res$b, res$rf)
    up_capt <- c(calc_up_capture(res$x, res$b, freq), NA)
    down_capt <- c(calc_down_capture(res$x, res$b, freq), NA)
    bat_avg <- c(apply(ar, 2, \(x){sum(x >= 0) / length(x)}), NA)
    act_ret <- geo_ret - geo_ret[length(geo_ret)]
    act_ret[length(act_ret)] <- NA
    ir <- act_ret / te
    ir[length(ir)] <- NA
    xcor <- cor(res$xb)[ncol(res$xb), ]
    xalpha <- (geo_ret - rf_geo) - (geo_ret[length(geo_ret)] - rf_geo) * xbeta
    xalpha[length(xalpha)] <- NA
    treynor <- geo_ret / xbeta
    treynor[length(treynor)] <- NA
    max_te_dd <- c(calc_max_drawdown(ar), NA)
    x <- rbind(
      scales::percent(geo_ret, 0.01),
      scales::percent(act_ret, 0.01),
      scales::percent(vol, 0.01), 
      scales::number(xbeta, 0.01),
      scales::percent(te, 0.01),
      scales::percent(down_vol, 0.01), 
      scales::percent(max_dd, 0.01),
      scales::percent(max_te_dd, 0.01),
      scales::number(sharpe, 0.01),
      scales::number(treynor, 0.01),
      scales::number(ir, 0.01),
      scales::number(sortino, 0.01), 
      scales::number(recov, 0.01),
      scales::number(act_ret / -max_te_dd, 0.01),
      scales::percent(up_capt, 0.01),
      scales::percent(down_capt, 0.01),
      scales::number(up_capt / down_capt, 0.01),
      scales::percent(bat_avg, 0.01),
      scales::percent(xalpha, 0.01),
      scales::percent(xcor, 0.01),
      scales::percent(xcor^2, 0.01),
      scales::percent(per_up, 0.01),
      scales::percent(per_down, 0.01),
      scales::percent(avg_up, 0.01),
      scales::percent(avg_down, 0.01),
      scales::number(avg_up / -avg_down, 0.01),
      scales::number((avg_up * per_up) / -(avg_down * per_down), 0.01)
    )
    xdf <- data.frame(
      Metric = c("Geometric Return", "Active Return", "Volatility",
                 "Beta", "Tracking Error", "Downside Vol", "Worst Drawdown",
                 "Worst Active Drawdown", "Sharpe Ratio", "Treynor Ratio",
                 "Info Ratio", "Sortino Ratio", "Recovery Ratio", 
                 "Active Recovery Ratio", "Up Capture", "Down Capture",
                 "Up Capt / Down Capt",
                 "Batting Average", "Alpha", "Correlation", "R-squared",
                 "Up Periods", "Down Periods", "Average Up Period",
                 "Average Down Period", "Avg Up / Avg Down",
                 "Wgt Avg Up / Wgt Avg Down"),
      x,
      row.names = NULL
    )
    num <- rbind(geo_ret, act_ret, vol, xbeta, te, down_vol, max_dd, max_te_dd,
                 sharpe, treynor, ir, sortino, recov, act_ret / -max_te_dd,
                 up_capt, down_capt, up_capt / down_capt, bat_avg, xalpha,
                 xcor, xcor^2, per_up, per_down, avg_up, avg_down, 
                 avg_up / -avg_down, (avg_up * per_up) / -(avg_down * per_down))
    numdf <- data.frame(
      Metric = c("Geometric Return", "Active Return", "Volatility",
                 "Beta", "Tracking Error", "Downside Vol", "Worst Drawdown",
                 "Worst Active Drawdown", "Sharpe Ratio", "Treynor Ratio",
                 "Info Ratio", "Sortino Ratio", "Recovery Ratio", 
                 "Active Recovery Ratio", "Up Capture", "Down Capture",
                 "Up Capt / Down Capt",
                 "Batting Average", "Alpha", "Correlation", "R-squared",
                 "Up Periods", "Down Periods", "Average Up Period",
                 "Average Down Period", "Avg Up / Avg Down",
                 "Wgt Avg Up / Wgt Avg Down"),
      num,
      row.names = NULL
    )
    res_out$xdf <- xdf
    res_out$numdf <- numdf
    return(res_out)
  }
}