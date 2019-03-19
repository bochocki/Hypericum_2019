### ENVIRONMENT SETUP ----------------------------------------------------------
library(cubature)
library(glue)
library(plotly)
library(tidyverse)

### FUNCTION DEFINITIONS -------------------------------------------------------


print_clear <- function(string) {
  # A function to clear then console and then print a string.
  #
  # Arguments:
  #   string: The string to be printed.
  #
  # Returns:
  #   None, but clears the console and prints the string.
  
  cat('\f')
  print(string)}


progress <- function(i, n, start_time) {
  # A completely non-essential function to print a progress bar and time
  # remaining.
  #
  # Arguments:
  #     i: The current iteration of the loop.
  #     n: The total number of iterations in the loop.
  #
  # Returns:
  #     None, but prints a progress bar
  
  # generate 'percent completion' strings
  pct <- i/n * 100
  len <- floor(pct / 2.5)
  bar <- paste0(rep('=', len), collapse='')
  spc <- paste0(rep(' ', 40 - len), collapse='')
  pct_complete <- sprintf('%.1f', pct)
  
  # generate timing strings
  t <- difftime(Sys.time(), start_time, units = 'mins')
  r <- (t / i) * (n - i)
  mt <- sprintf('%02.f', floor(t))
  st <- sprintf('%02.f', (t - floor(t)) * 60)
  mr <- sprintf('%02.f', floor(r))
  sr <- sprintf('%02.f', (r - floor(r)) * 60)

  if (i == n) {
    print_clear(glue('|{bar}{spc}| 100.0% complete!!!\n',
                     'Total time: {mt}:{st}\n'))
  } else {
    print_clear(glue('|{bar}{spc}| {pct_complete}%\n',
                     'Elapsed time: {mt}:{st}\n',
                     'Remaining time: {mr}:{sr}'))
  }
}


revolved_exponential <- function(xy, lambda) {
  # This function returns the probability of dispersing from the origin to a 
  # point in the x-y plane, using the revolved exponential distribution.
  #
  # Arguments:
  #     xy: A vector giving the (x, y) coordinates of the point of interest.
  #     lambda: The exponential distribution's rate parameter.
  #
  # Returns:
  #     Pr(x,y | lambda): the probability of dispersing from (0, 0) to (x, y)
  #     given lambda.
  
  # calculate the scale factor and the distrance from (0, 0) to (x, y)
  scale <- (lambda^2) / (2*pi)
  dist  <- sqrt(xy[1]^2 + xy[2]^2)
  
  # return the probability of dispersing to (x, y)
  return(scale * exp(-lambda * dist))}


exp_distance <- function(L2x2, L2y2) {
  # This is a helper function for `simpson_approx()`, which returns the
  # exponential portion of the revolved exponential kernel. Pre-computing the
  # arguments is more convenient, computationally, than doing everything
  # explicitly in order.
  #
  # NOTE: x and y are squared (and multiplied by lambda^2) **before** being
  #     passed to this function
  #
  # Arguments:
  #     L2x2: A vector of squared x-coordinates,  multiplied by lambda^2.
  #     L2y2: A vector of squared y-coordinates,  multiplied by lambda^2.
  #
  # Returns:
  #     The equivalent of exp(-L * sqrt(x^2 + y^2))
  
  return(exp(-sqrt(L2x2 + L2y2)))
}


simpson_approx <- function(A, B, C, D, lambda) {
  # A function that uses Simpson's Rule to estimate the double-integral of the
  # revolved exponential kernel for a given region of space.
  #
  # This approximation works best over relatively small regions, like the 5x5m
  # grid cells in the current study design.
  #
  # Arguments:
  #     A: The lower bound of the interval over x
  #     B: The upper bound of the interval over x
  #     C: The lower bound of the interval over y
  #     D: The upper bound of the interval over y
  #     lambda: The exponential distribution's rate parameter.
  #
  # Returns:
  #     The total probability of dispersing from (0, 0) to a region in the xy 
  #     plane bounded by {A, B, C, D}:
  #
  #         (A, D) --- (B, D)
  #           |           |
  #           |           |
  #           |           |
  #         (A, C) --- (B, C)
  
  # calculate lamda^2 and two useful scaling factors
  L2 <- lambda^2
  s <- (B-A)*(D-C) / 36
  alpha <- s * L2 / (2*pi)
  
  # pre-calculate lambda^2 times {A^2, B^2, C^2, D^2}, and the midpoints of AB 
  # and CD (squared)
  AB <- L2 * ((A + B) / 2)^2
  CD <- L2 * ((C + D) / 2)^2 
  A  <- L2 * A^2
  B  <- L2 * B^2
  C  <- L2 * C^2
  D  <- L2 * D^2
  
  # calculate the exponential portion of the dispersal kernel for each xy-pair
  # required in the approximation (and scale these values when appropriate)
  AB_CD = exp_distance(AB, CD) * 16
  AB_C  = exp_distance(AB, C) * 4
  AB_D  = exp_distance(AB, D) * 4
  A_CD  = exp_distance(A, CD) * 4
  B_CD  = exp_distance(B, CD) * 4
  A_C   = exp_distance(A, C)
  B_C   = exp_distance(B, C)
  A_D   = exp_distance(A, D)
  B_D   = exp_distance(B, D)
  
  # approximate the integral by summing all pieces and scaling them by alpha
  return(alpha * (AB_CD + AB_C + AB_D + A_CD + B_CD + A_C + B_C + A_D + B_D))
}


dispersal_between_balds <- function(f, id='bald_num', lambda=0.001, width=5,
                                    print_progress=TRUE) {
  # A function to calculate the fraction of seeds that disperse from each bald
  # to every other bald in the landscape, including seeds that remain in the
  # source bald.
  #
  # Arguments:
  #     f: The path to the data file that contains information on lattitude and
  #         longitude for each grid cell, and the bald ID that the cell belongs
  #         to.
  #     id: The column that should be used to identify the bald that each grid
  #         cell belongs to.
  #     lambda: The exponential distribution's rate parameter.
  #     width: The width of each square grid cell (in meters).
  #
  # Returns:
  #     A dataframe where each column describes the fraction of seeds that
  #     disperse from the bald identified by that column to every other bald in
  #     the landscape. For example, the column labeled '10' describes the
  #     fraction of seeds that disperse from Bald 10 to every other bald in the
  #     landscape, including the number of seeds that remain in Bald 10, and the
  #     fraction of seeds that do not arrive in any bald, and are lost to the
  #     surrounding matrix.
  
  # load the data
  df <- read_csv(f, col_types = cols())
  
  # ensure the id column is a character (for indexing)
  df[[id]] <- as.character(df[[id]])
  
  # subset the dataframe to the id column and coordinates
  df <- df[c(id, 'easting', 'northing')]
  
  # gather the unique bald ids and initiate an output dataframe
  names <- unique(df[[id]])
  out <- as.data.frame(matrix(0, nrow=length(names), ncol=length(names)))
  rownames(out) <- names
  colnames(out) <- names
  
  # add columns to df to track the x- and y-distance between patches, and
  # the probability of dispersing that distance
  df$x <- 0
  df$y <- 0
  df$p <- 0
  
  # create a timer to track how long the process takes
  tic <- Sys.time()
  
  # loop over every point in the dataframe
  ni <- nrow(df)
  w <- width / 2
  for (i in 1:ni) {
    
    # get the id and (x, y) coordinates for point [i]
    bald <- df[[id]][i]
    x0 <- df$easting[i]
    y0 <- df$northing[i]
    
    # calculate the x- and y-distance between point [i] and all others
    df$x <- df$easting - x0
    df$y <- df$northing - y0
    
    # calculate the prob. of dispersing to every other patche from point [i]
    df$p <- simpson_approx(df$x-w, df$x+w, df$y-w, df$y+w, lambda = lambda)
    
    # count the total proportion of seeds that each bald receives from [i]
    tmp <- df %>% 
      group_by_(sym(id)) %>%
      summarize(proportion = sum(p))
    
    # add these bald-specific proportions to the relevant row in the output
    out[tmp[[id]], bald] <- out[tmp[[id]], bald] + tmp$proportion
    
    # print a progress update to the console
    if ((i %% 20 == 0) & print_progress) {
      progress(i, ni, tic)
    }
  }
  
  # normalize the values in the output
  out <- out / ni
  
  # calculate the proportion of seeds lost to the matrix, and append to output
  out[nrow(out) + 1, ] <- 1 - colSums(out)
  rownames(out) <- c(names, 'matrix')
  out[id] <- rownames(out)
  
  # print some final processing updates to the user.
  progress(ni, ni, tic)
  
  return(out)
}


theme_black = function(base_size = 12, base_family = '') {
  # A black plotting theme, taken (with slight modifications) from
  # https://gist.github.com/jslefche/eff85ef06b4705e6efbc
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = 'white'),  
      axis.text.y = element_text(size = base_size*0.8, color = 'white'),  
      axis.ticks = element_line(color = 'grey20', size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = 'white',
                                  margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = 'white', angle = 90,
                                  margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, 'lines'),
      
      # Specify legend options
      legend.background = element_rect(color = NA, fill = 'black'),  
      legend.key = element_rect(color = 'white',  fill = 'black'),  
      legend.key.size = unit(1.2, 'lines'),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = 'white'),  
      legend.title = element_text(size = base_size*0.8, face = 'bold',
                                  hjust = 0, color = 'white'),  
      legend.position = 'right',  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = 'vertical',  
      legend.box = NULL, 
      
      # Specify panel options
      panel.background = element_rect(fill = 'black', color  =  NA),  
      panel.border = element_rect(fill = NA, color = 'grey20'),  
      panel.grid.major = element_line(color = 'grey20'),  
      panel.grid.minor = element_blank(),  
      
      # Specify facetting options
      strip.background = element_rect(fill = 'grey20', color = 'grey20'),  
      strip.text.x = element_text(size = base_size*0.8, color = 'white'),  
      strip.text.y = element_text(size = base_size*0.8, color = 'white',
                                  angle = -90),  
      
      # Specify plot options
      plot.background = element_rect(color = 'black', fill = 'black'),  
      plot.title = element_text(size = base_size*1.2, color = 'white',
                                margin = margin(0, 0, 20, 0)),  
      plot.margin = unit(rep(1, 4), 'lines')
    )
}

plot_dispersal <- function(f, disp, bald, id='bald_num', logscale=TRUE,
                           highlight=TRUE) {
  # Draw a plot of the fraction of seeds that disperse to each bald from a focal
  # bald.
  #
  # Arguments:
  #     f: The path to the data file that contains information on lattitude and
  #         longitude for each grid cell, and the bald ID that the cell belongs
  #         to.
  #     disp: The dispersal dataframe returned by `dispersal_between_balds`.
  #     bald: The focal bald for the plot. Should be one of the unique
  #         identifiers given by `id`.
  #     id: The column that should be used to identify the bald that each grid
  #         cell belongs to.
  #     logscale: Should color values be plotted on the log10 scale? Defaults to
  #         TRUE.
  #     highlight: Should the bald designated by `bald` be highlighted? Defaults
  #         to TRUE.
  #
  # Returns:
  #     None, but draws the plot.
  
  # assign a colorscale for the plot
  colorscale = c('dodgerblue', 'white', 'tomato')
  
  # GATHER DATA
  # load the original data
  df <- read_csv(f, col_types = cols())
  
  # ensure the id column is a character (for indexing)
  df[[id]] <- as.character(df[[id]])
  
  # join original data with the calculated dispersal data
  all <- dplyr::left_join(df, disp, by=id)

  # MANIPULATE DATA
  # gather **only the bald data columns** from the dispersal data
  disp_mat <- as.matrix(disp[disp[[id]] != 'matrix', colnames(disp) != id])
  
  # assign plot variables based on whether or not a log-scale plot will be used
  if (logscale) {
    
    # scale transformation
    log_trans <- 'log10'  
    
    # legend label
    prob_scale <- 'Log10(fraction)'
    
    # colorscale limits and values
    li <- c(NA, max(disp_mat))
    va <- c(0, log10(max(disp_mat))/log10(median(disp_mat)), 1)
  } else {
    
    # scale transformation
    log_trans <- 'identity'
    # legend label
    prob_scale <- 'Fraction'
    
    # colorscale limits and values
    li <- c(0, max(disp_mat))
    va <- c(0, median(disp_mat) / max(disp_mat), 1)
  }
  
  # build the plot
  p <- ggplot(all,
              aes_(x = quote(easting),
                   y = quote(northing),
                   colour = as.name(bald))) + 
    theme_black() +
    coord_fixed() +
    ggtitle(glue('Dispersal from Bald {bald}\n')) +
    ylab('Northing') +
    xlab('Easting') +
    geom_point(size = 0.01) +
    scale_colour_gradientn(colours = colorscale,
                           limits = li,
                           values = va,
                           trans = log_trans,
                           name = glue('{prob_scale} of Bald {bald} seeds\n',
                                       'dispersing to patch'))
  
  # if highlighting is to be used, add it here
  if (highlight) {
    p <- p +
      geom_point(data = all[all[[id]] == bald, ],
                 aes(x = easting, y = northing),
                 colour='yellow',
                 size=1) +
      geom_point(data = all[all[[id]] == bald, ],
                 size=0.01)
  }
  
  # draw the plot
  print(p)
}