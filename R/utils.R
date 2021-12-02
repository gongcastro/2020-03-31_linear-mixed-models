# function for simulating hierarchical data
simulate_rts <- function(
    intercepts, # subject real intercepts
    sd, # real intercept SD
    corr, # correlation between slopes and intercepts
    n_conditions, # number of within participant conditions
    n_trials # number of trials per conditions (cell size)
) {
    
    n_subj <- length(intercepts)
    
    means <- c(rbind({{ intercepts }}, {{ intercepts }} * {{ corr }}))
    
    df <- replicate(
        rtruncnorm(n = 1, a = 0, b = Inf, sd = sd, mean = means),
        n = n_trials
    ) %>%
        as.data.frame() %>% 
        pivot_longer(cols = everything(), names_to = "trial", values_to = "rt") %>%
        mutate(
            participant = rep(paste0("Participant ", 1:n_subj), each = n_trials*n_conditions),
            trial = str_remove(trial, "V"),
            frequency = rep(c("High", "Low"), each = n_trials, times = n_subj)
        )  %>%
        select(participant, frequency, trial, rt)
    
    return(df)
}

# get coefficients
get_coefs <- function(x){
    coefs <- x$coefficients %>%
        as.data.frame() %>%
        rename(coefficient = "Estimate", sem = "Std. Error") %>%
        rownames_to_column("term") %>%
        mutate(term = c("Intercept", "Slope")) %>%
        select(term, coefficient, sem)
    
    return(coefs)
}

# set custom ggplot theme as default
theme_custom <- function(){
    theme(
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey10"),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text (colour = "black", family = "Helvetica"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", colour = NA),
        legend.title = element_blank()
    )
}