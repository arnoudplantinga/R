# Proportion tables ---------------------------------------------------------------------------

#############
# proptable #
#############

# Function to create a frequency table as data frame or pander table.
proptable <- function(
  var1, var2 = NULL, # Variables
  shape = "wide",    # Shape of the data frame, long or wide
  decimals = 1,      # Round
  cumulative = F,    # Cumulative frequencies for single variable
  na.rm = F,         # Delete NAs when shape is "long"
  var = "all",       # Choose variables shown in data frame: all, n, prop, nprop
  output = "df"      # Output: df (data frame), or pander (table, can be used in RMarkdown)
  ) {
  
  library(tidyverse)
  name_var1 <- deparse(substitute(var1)) %>% str_remove(".*\\$")
  name_var2 <- deparse(substitute(var2)) %>% str_remove(".*\\$")

  # 1 variable cumulative
  if (is.null(var2) & cumulative) {
    sum = length(var1)
    df <- data_frame(var1) %>%
      count(var1) %>%
      mutate(n = cumsum(n)) %>%
      mutate(
        prop = n / sum,
        nprop = paste0(n, " (", formatC(prop*100, decimals, format = "f"), "%)")
      )
    
    if (shape == "long") {
      if (na.rm) df <- filter(df, !is.na(var1))
      res <- df
    }
    else {
      res <- df$nprop
      names(res) <- df$var1
    }
  }

  # 1 variable not cumulative
  if (is.null(var2) & !cumulative) {
    df <- data_frame(var1) %>%
      count(var1) %>%
      mutate(prop = prop.table(n)) %>%
      mutate(nprop = paste0(n, " (", formatC(prop*100, decimals, format = "f"), "%)"))
    if (shape == "long") {
      if (na.rm) df <- filter(df, !is.na(var1))
      res <- df
    } else {
      res <- df$nprop
      names(res) <- df$var1
    }
  }   

   # 2 variables
  if (!is.null(var2)) {
    df <- data_frame(var1, var2) %>%
      count(var1, var2) %>%
      group_by(var1) %>%
      mutate(prop = prop.table(n)) %>%
      mutate(nprop = paste0(n, " (", formatC(prop*100, decimals, format = "f"), "%)")) %>% 
      ungroup()
    if (shape == "long" & var == "n") {
      df <- df %>%
        dplyr::select(var1, var2, n)
    }
    if (shape == "long" & var == "prop") {
      df <- df %>%
        dplyr::select(var1, var2, prop)
    }
    if (shape == "long" & var == "nprop") {
      df <- df %>%
      dplyr::select(var1, var2, nprop) %>%
      mutate(nprop = ifelse(is.na(nprop), "0 (0.0%)", nprop))
    }
    if (shape == "wide") {
      df <- df %>% 
        dplyr::select(var1, var2, nprop) %>% 
        spread(var2, nprop) %>% 
        mutate_if(is.character, ~ ifelse(is.na(.), "0 (0.0%)", .))
      names(df)[1] <- name_var1
      names(df)[-1] <- paste0(name_var2, "_", names(df)[-1])
    }
    if (shape == "long") {
      names(df)[1] <- name_var1
      names(df)[2] <- name_var2
    }
    res <- df
  }
  
  if (output == "pander") {
    library(pander)
    return(pander(res, justify = "right", style = "rmarkdown"))
  } else {
    return(res)
  }
}   

############
# proptabr #
############

# Function to create a frequency table using the descr package. Output is text.
proptabr <- function(
  var1, var2 = NULL, # Variables
  output = "text",   # Output: text (table in text format) or pander 
                     #   (table, can be used in RMarkdown)
  format = "SPSS",   # Format used by CrossTable. SPSS (percentages) or SAS (proportions).
  chisq = F,         # Show Chi-squared test?
  ...                # Other arguments passed to CrossTable
  ) {
  
  library(descr)
  
  name_var1 <- deparse(substitute(var1)) %>% str_remove(".*\\$")
  name_var2 <- deparse(substitute(var2)) %>% str_remove(".*\\$")
  
  if (!is.null(var2)) {
    tab <- CrossTable(var1, var2, prop.c = F, prop.t = F, prop.chisq = F, chisq = chisq,
      format = format, cell = F, percent = T,
      dnn = c(name_var1, name_var2), ...)
  } else {
    tab <- CrossTable(var1, format = format, cell.layout = F, chisq = chisq, percent = T)
  }
  
  if (output == "pander") {
    pander(tab, digits = 1)
  } else {
    tab
  }
}

# Robust standard errors ----------------------------------------------------------------------

###########
# robust2 #
###########

# Regression output with robust standard errors
robust2 <- function(mod) {
  robust(mod) %>%
  mutate(p.value = round(p.value, 3))
}   

####################
# stargazer_robust #
####################

# Regression output in a fancy table with robust standard errors
stargazer_robust <- function(reg, type = "htm1", decimals = 2, ...) {
  library(sjstats)
  library(stargazer)
  arguments = list(...)

  if (sum(class(reg) == "Tm") > 0) {
    se_robust <- robust(reg)$std.error %>% list()
    stargazer(reg, type = type, se = se_robust, star.cutoffs = c(.05, .01, .001), digits = decimals,
      ...)
  }
  if (sum(class(reg) == "rlmermod") > 0) {
    s <- summary(reg)
    pvals <- 2 * (1 - pt(abs(as.data.frame(s$coefficients)$'t value'), s$ngrps))
    res <- as.data.frame(s$coefficients)
    res$p_approx <- pvals
    return(round(res, 4))
  }
  if (class(reg) == "Tist") {
    se_robust <- map(reg, ~ robust(.)$std.error)
    stargazer(reg, type = type, se = se_robust, star.cutoffs = c(.05, .01, .001), digits = decimals,
    ...)
  }
  }

# Regression ----------------------------------------------------------------------------------

# Get standardized coefficients
lm_beta <- function (model) {
  b <- summary(model)$coef[-1, 1]
  sx <- sapply(model$model[-1], sd)
  sy <- sapply(model$model[1], sd)
  beta <- b * sx/sy
  return(beta)
}

# Plot fitted values against residuals
plot_residuals <- function(model) {
  require(ggplot2)
  
  fitted <- fitted(model)
  residuals <- residuals(model)
  
  df <- data.frame(fitted, residuals)
  
  ggplot(df, aes(fitted, residuals)) +
    geom_point() +
    geom_smooth(method = 'loess') +
    geom_abline(intercept = 0, slope = 0)
}

# Plot interaction
plot_interaction <- function(X, Y, GROUP) {
  require(effects); require(ggplot2)
  
  name_X <- deparse(substitute(X)) %>% str_remove(".*\\$")
  name_Y <- deparse(substitute(Y)) %>% str_remove(".*\\$")
  name_GROUP <- deparse(substitute(GROUP)) %>% str_remove(".*\\$")
  
  lm1 <- lm(Y ~ X * GROUP)
  d <- data.frame(allEffects(lm1)[[1]])
  d$GROUP <- as.factor(d$GROUP)
  df <- data.frame(X, Y, GROUP)

  ggplot() + 
    geom_line(data = d, aes(X, fit, linetype = GROUP, col = GROUP), size = 1.2) +
    geom_point(data = df, aes(x = X, y = Y, col = as.factor(GROUP)), size = 1.5) +
    xlab(name_X) + ylab(name_Y) +
    guides(linetype = guide_legend(title = name_GROUP),
           col      = guide_legend(title = name_GROUP)) +
    theme_minimal()
}
