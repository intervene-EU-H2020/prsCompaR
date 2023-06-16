#' Title
#'
#' @param plotdata
#' @param metric
#' @param plot_width
#' @param plot_height
#' @param legend_position
#' @param vline
#' @param plot_legend
#' @param dodge_width
#' @param point_size
#' @param nrow
#' @param indicate_best
#' @param show_cv_equal_auto
#' @param to_pdf
#' @param scales
#' @param shape
#' @param x_variable
#' @param y_lab
#' @param strip_text_x_size
#' @param axis_text_y_size
#' @param axis_title_y_size
#' @param adjust_methods_levels_fun
#'
#'
#' @return
#' @export
#'
#' @examples
plot_effect_sizes_full <- function(plotdata,
                                   metric,
                                   plot_width,
                                   plot_height,
                                   legend_position = 'bottom',
                                   vline = NA,
                                   plot_legend = F,
                                   dodge_width = 0.8,
                                   point_size = 0.3,
                                   nrow = 2,
                                   indicate_best = F,
                                   show_cv_equal_auto = T,
                                   scales = 'fixed',
                                   shape = c(1, 16),
                                   x_variable = 'method',
                                   y_lab = NA,
                                   strip_text_x_size = 12,
                                   axis_text_y_size = 12,
                                   axis_title_y_size = 20,
                                   adjust_methods_levels_fun = adjust_methods_levels) {
  if (metric == 'OR') {
    y <- rlang::sym('OR')
    y_min <- rlang::sym('OR_CI_LOW')
    y_max <- rlang::sym('OR_CI_HIGH')
    y_lab <- ifelse(is.na(y_lab), 'OR', y_lab)
    y_intercept <- 1
    data_range <-
      c(min(plotdata[['OR_CI_LOW']]), max(plotdata[['OR_CI_HIGH']]))
  } else if (metric == 'BETA') {
    y <- rlang::sym('BETA')
    y_min <- rlang::sym('BETA_CI_LOW')
    y_max <- rlang::sym('BETA_CI_HIGH')
    y_lab <- ifelse(is.na(y_lab),  latex2exp::TeX('\\beta'), y_lab)
    y_intercept <- 0
    data_range <-
      c(min(plotdata[['BETA_CI_LOW']]), max(plotdata[['BETA_CI_HIGH']]))
  } else if (metric == 'AUC') {
    y <- rlang::sym('AUC_MEDIAN')
    y_min <- rlang::sym('AUC_CI_LOW')
    y_max <- rlang::sym('AUC_CI_HIGH')
    y_lab <- ifelse(is.na(y_lab),  'AUROC', y_lab)
    y_intercept <- 0.5
    data_range <-
      c(min(plotdata[['AUC_CI_LOW']]), max(plotdata[['AUC_CI_HIGH']]))
  }   else if (metric == 'h2l') {
    y <- rlang::sym('h2l')
    y_min <- rlang::sym('h2l_ci_low')
    y_max <- rlang::sym('h2l_ci_high')
    y_lab <- ifelse(is.na(y_lab),  latex2exp::TeX('$r_{liab}^{2}$') , y_lab)
    y_intercept <- 0
    data_range <-
      c(min(plotdata[['h2l_ci_low']]), max(plotdata[['h2l_ci_high']]))
    print(data_range)
  }   else if (metric == 'R2_OBS') {
    y <- rlang::sym('R2_OBS')
    y_min <- rlang::sym('R2_OBS_CI_LOW')
    y_max <- rlang::sym('R2_OBS_CI_HIGH')
    y_lab <- ifelse(is.na(y_lab),  latex2exp::TeX('$r^2$') , y_lab)
    y_intercept <- 0
    data_range <-
      c(min(plotdata[['R2_OBS_CI_LOW']]), max(plotdata[['R2_OBS_CI_HIGH']]))
  }

  y_variable <- as.character(y)
  y_min_variable <- as.character(y_min)
  y_max_variable <- as.character(y_max)


  discrete_cols <- get_discrete_colours()
  if (!('pt.clump' %in% plotdata$method)) {
    plotcolors <- discrete_cols[2:length(discrete_cols)]
  } else {
    plotcolors <- discrete_cols
  }


  if (show_cv_equal_auto) {
    plotdata[method %in% c('prscs', 'lassosum', 'ldpred2', 'megaprs'), auto_equal_cv :=
               !('CV' %in% method_type), by = list(ancestry, phenotype, study, method)]
    tmp <- plotdata[auto_equal_cv == TRUE]
    tmp <- tmp[!(phenotype %in% get_auto_only_phenotypes())]
    tmp$method_type <- 'CV'
    plotdata <- data.table::rbindlist(list(plotdata, tmp))
    rm(tmp)
  }

  plotdata <- rename_phenotypes(plotdata)
  # plotdata$method <- adjust_methods_levels_fun(plotdata$method)
  plotdata$method_type <-
    factor(plotdata$method_type, levels = c('auto', 'CV'))

  library(ggplot2)
  if (x_variable == 'method') {
    plotdata$x_variable <-
      paste(plotdata$method, plotdata$method_type, sep = ':')
    plot_full <-
      ggplot2::ggplot(plotdata,
             aes(
               x = method,
               y = !!y,
               col = method,
               group = x_variable,
               shape = method_type
             ))
  } else if (x_variable == 'biobank') {
    plotdata$x_variable <-
      paste(plotdata$method, plotdata$method_type, sep = ':')
    plot_full <-
      ggplot(plotdata,
             aes(
               x = bbid,
               y = !!y,
               col = method,
               group = x_variable,
               shape = method_type
             ))
  } else {
    stop('x_variable has to be "method" or "biobank"')
  }


  plot_full <- plot_full +
    geom_hline(
      yintercept = y_intercept,
      col = 'black',
      alpha = 0.5,
      size = 0.2
    ) +
    geom_point(position = position_dodge2(width = dodge_width, preserve =
                                            'single'),
               size = point_size) +
    geom_errorbar(
      position = position_dodge2(width = dodge_width, preserve = 'single'),
      width = dodge_width,
      aes(ymin = !!y_min, ymax = !!y_max)
    )

  if (indicate_best) {
    # add a triangle to indicate the best method
    plotdata2 <- data.table::copy(plotdata)
    if ('bbid' %in% colnames(plotdata)) {
      plotdata2[, select := max(get(y_variable)) == get(y_variable), by = list(bbid, ancestry, phenotype)]
    } else {
      plotdata2[, select := max(get(y_variable)) == get(y_variable), by = list(ancestry, phenotype)]
    }
    plotdata2[select == FALSE, c(y_variable, y_min_variable, y_max_variable) :=
                NA, ]

    if (scales == 'fixed') {
      offset <- (data_range[2] - data_range[1]) * 0.05
      plot_full <-
        plot_full + geom_point(
          data = plotdata2,
          aes(y = !!y_max + offset, fill = method),
          position = position_dodge2(width = 0.8, preserve = 'single'),
          shape = 25
        )
    } else {
      offsets <-
        plotdata[, list(offset = (abs(max(
          get(y_max_variable)
        ) - min(
          get(y_min_variable)
        ))) * 0.05), by = list(phenotype)]
      plotdata2 <-
        merge(
          plotdata2,
          offsets,
          all.x = T,
          all.y = F,
          by = 'phenotype'
        )
      plot_full <-
        plot_full + geom_point(
          data = plotdata2,
          aes(y = !!y_max + offset, fill = method),
          position = position_dodge2(width = 0.8, preserve = 'single'),
          shape = 25
        )

    }

  }

  if (!(is.na(vline))) {
    plot_full <-
      plot_full + geom_vline(
        xintercept = vline,
        col = 'black',
        linetype = 1,
        alpha = 0.5,
        size = 0.4
      )
  }

  if (x_variable == 'method') {
    plot_full <-
      plot_full + facet_grid(study + phenotype ~  ancestry + bbid, scales = scales)
  } else if (x_variable == 'biobank') {
    plot_full <-
      plot_full + facet_grid(study + phenotype ~  ancestry + method, scales = scales)
  }

  plot_full <- plot_full +
    scale_color_manual(values = plotcolors) +
    scale_fill_manual(values = plotcolors) +
    scale_shape_manual(values = shape) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      strip.background = element_blank(),
      panel.border = element_rect(color = 'grey80', fill = NA),
      legend.position = legend_position,
      strip.text.x = element_text(size = strip_text_x_size, margin = margin(0, 0, 2, 0, 'pt'))
    ) +
    ylab(y_lab) +
    xlab('')

  if (scales != 'fixed') {
    plot_full <-
      plot_full + theme(axis.text.y = element_text(size = axis_text_y_size))
  }

  if (!is.null(axis_title_y_size)) {
    plot_full <-
      plot_full + theme(axis.title.y = element_text(size = axis_title_y_size))
  }

  if (!plot_legend) {
    plot_full <- plot_full + theme(legend.position = 'none')
  }

  return(plot_full)
}
