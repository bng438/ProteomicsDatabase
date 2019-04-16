
# T-test function  ----
pVal <- function(dt, grp1, grp2)
{
  x <- dt[grp1]
  y <- dt[grp2]
  
  # Performs t-test using the mean of x and y
  result <- t.test(x, y, 
                   alternative="two.sided",
                   mu=0, paired=FALSE, var.equal=TRUE)
  
  # Extracts and returns p-values from the results
  return(result$p.value)
}


# Fold-change function  ----
foldChange <- function(dt, grp1, grp2)
{
  x <- dt[grp1] %>% unlist %>% as.numeric() %>% mean()
  y <- dt[grp2] %>% unlist %>% as.numeric() %>% mean()
  
  fold_change <- (x / y)
  return(fold_change)
}



# Applies specified function for all combinations of groups  ----
toAllGroups <- function(dat, num_repl, func)
{
  result <- as.data.frame(matrix(0,nrow(dat)))
  
  for (i in seq(1, (ncol(dat) - (2 * num_repl - 1)), by=num_repl))
  {
    for (j in seq((num_repl + i), (ncol(dat) - num_repl + 1), by=num_repl))
    {
      # Applies the specified function func to all rows in the
      # two groups currently being compared
      result[(length(result) + 1)] <-
        adply(.data=dat, .margins=1, .fun=get(func),
              grp1=c(i : (i + num_repl - 1)),
              grp2=c(j : (j + num_repl - 1))) %>%
        as_tibble() %>% select("V1")
      
      # Determines column heading depending on function specified
      if (strcmp(func, "pVal"))
      {
        func_heading <- "Pval"
      }
      else 
      {
        func_heading <- "Fold-change"
      }
      
      # Renames column heading
      names(result)[length(result)] <-
        paste(func_heading,
              substr(names(dat[i]), 1, nchar(names(dat[i]))-1),
              "vs",
              substr(names(dat[j]), 1, nchar(names(dat[j]))-1))
    }
  }
  return(select(result,-1))
}


# Identifies NAN error values  ----
is.nan.data.frame <- function(x)
{
  do.call(cbind, lapply(x,is.nan))
}


# Replaces NAN error values with 0  ----
removeNan <- function(dat)
{
  dat[is.nan(dat)] <- 0
  return(dat)
}


# Merges datasets such that corresponding columns are next to each other  ----
# I.e. 1st column of dat2 is next to 1st column of dat1, etc...
merge <- function(dat1, dat2)
{
  resorted <- as.data.frame(matrix(0,nrow(dat1)))
  for (i in 1:ncol(dat1))
  {
    resorted[(length(resorted)+1)] <- dat1[i]
    resorted[(length(resorted)+1)] <- dat2[i]
  }
  return(select(resorted,-1))
}


# Rounds values to desired number of decimal places  ----
roundValues <- function(dat,x)
{
  mutate_all(dat, round, x)
}
