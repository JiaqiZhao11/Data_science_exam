#Stan functions


get_diag = function(fit, parameter, prior = FALSE, nsubs = NA){
  
  variable = fit$metadata()$stan_variables
  
  all_parameters = grep(parameter, variable, value = TRUE)
  
  priors = grep("prior", all_parameters, value = TRUE)
  
  posteriors = all_parameters[!all_parameters %in% priors]
  
  search = c("mu","kappa","sd")
  
  
  result <- sapply(search, function(x) grepl(x, priors, fixed = TRUE))
  
  non_hier_prior = priors[which(rowSums(!result) == ncol(result))]
  
  hier_prior = priors[!priors %in% non_hier_prior]
  
  result <- sapply(search, function(x) grepl(x, posteriors, fixed = TRUE))
  
  non_hier_post = posteriors[which(rowSums(!result) == ncol(result))]
  
  hier_post = posteriors[!posteriors %in% non_hier_post]
  
  
  
  if(!prior){
    
    if(is.na(nsubs)){
      nsubs = length(posteriors)
    }
    nplots = ceiling(nsubs/6)
    
    trace_plots = list()
    pairs_plots = list()
    
    for(i in 1:nplots){
      
      if(i != nplots){
        
        indicies = c(((i - 1) * 5 + 1):(i * 5))
        variabless = c(paste0(non_hier_post,"[",indicies,"]"),hier_post)
        
        trace = mcmc_trace(fit$draws(variables = variabless))+
          theme(strip.text = element_text(size = 18), axis.text = element_text(size = 18),axis.title = element_text(size = 18))
        
        
        pairs = mcmc_pairs(fit$draws(variables = variabless), np = nuts_params(fit), pars = variabless,
                   off_diag_args = list(size = 0.75))
        
        trace_plots[[i]] = trace
        pairs_plots[[i]] = pairs
        
      }else{
        
        indicies = c((i * 5):nsubs)
        variabless = c(paste0(non_hier_post,"[",indicies,"]"),hier_post)
        
        
        trace = mcmc_trace(fit$draws(variables = variabless))+
          theme(strip.text = element_text(size = 18), axis.text = element_text(size = 18),axis.title = element_text(size = 18))
        
        
        pairs = mcmc_pairs(fit$draws(variables = variabless), np = nuts_params(fit), pars = variabless,
                           off_diag_args = list(size = 0.75))
        
        trace_plots[[i]] = trace
        pairs_plots[[i]] = pairs
      }
    }
    
  
    return(list(trace_plots = trace_plots, pairs_plots = pairs_plots))
  }
}

ppu_hier = function(fit,parameters,reals){
  
  
  #things to search for 
    
    
    search = c("mu","kappa","sd")
  
    plots = list()
    q = 0
    for(parameter in parameters){
      
      post_parameters = paste0(search,"_",parameter)
      
      post_parameters = post_parameters[post_parameters %in% fit$metadata()$stan_variables]
      
      prior_parameters = paste0("prior_",search,"_",parameter)
      
      prior_parameters = prior_parameters[prior_parameters %in% fit$metadata()$stan_variables]
      
      q = q+1

      
      post = as_draws_df(fit$draws(variables = post_parameters)) %>% 
        mutate(posterior = T)
      
      priors = as_draws_df(fit$draws(variables = prior_parameters)) %>% 
        mutate(posterior = F) %>%
        rename_with(~sub("^prior_", "", .), starts_with("prior_"))
      
      
      df = rbind(post,priors)
      
      if(!is.na(reals)[1]){
        reals = hier %>% select(grep(parameter,names(hier), value = TRUE))
        
        plot1 = df %>% pivot_longer(cols = post_parameters)%>% 
          ggplot(aes(x = value, fill = posterior))+
          geom_histogram(alpha = 0.5, position="identity")+
          facet_wrap(~name, scales = "free")+
          theme_classic()+
          geom_vline(data = reals %>% pivot_longer(cols = post_parameters), aes(xintercept = value))
      
        plots[[q]] = plot1
      }else{
        plot1 = df %>% pivot_longer(cols = post_parameters)%>% 
          ggplot(aes(x = value, fill = posterior))+
          geom_histogram(alpha = 0.5, position="identity")+
          facet_wrap(~name, scales = "free")+
          theme_classic()
        
        plots[[q]] = plot1
        
        }
    }
    
    final_plot <- wrap_plots(plots, ncol = 1)
    
    return(final_plot)
  
    
}


ppu_sub = function(fit,parameter,reals, lims = -10, lime = 10, nsubs = NA){
  
  
  #get rid of these:
  search = c("mu","kappa","sd")
  #get the names of the parameters
  post_parameters = grep(parameter, fit$metadata()$stan_variables, value = TRUE)
  #find where the search term is.
  result <- sapply(search, function(x) grepl(x, post_parameters, fixed = TRUE))
  #find where any of them are
  result <- apply(result, 1, any)
  #get rid of them
  parameters = post_parameters[!result]
  
  
  prior_parameters = grep("prior",parameters, value = TRUE)
  
  post_parameters = parameters[!parameters %in% prior_parameters]
  
  post = as_draws_df(fit$draws(variables = post_parameters)) %>%
    mutate(posterior = T)
  
  priors = as_draws_df(fit$draws(variables = prior_parameters)) %>%
    mutate(posterior = F) %>%
    rename_with(~sub("^prior_", "", .), starts_with("prior_"))
  
  
  df = rbind(post,priors)
  
  if(!is.na(reals)[1]){
    
    nsubs = length(unique(reals$id))
    reals = reals %>% select(parameter,"id") %>% group_by(id) %>% summarize(mean = mean(!!sym(parameter)))
    
    data = df %>% pivot_longer(cols = contains(post_parameters)) %>% drop_na()
    data1 = df %>% pivot_longer(cols = contains(post_parameters)) %>% drop_na() %>% filter(posterior == T)
    
    data %>% 
      ggplot(aes(x = value, fill = posterior))+
      geom_histogram(alpha = 0.5, position="identity")+
      facet_wrap(~name, nrow = round(nsubs/2,0), ncol = 5, scales = "free")+
      geom_vline(data = reals %>% mutate(name = paste0(parameter,"[", row_number(), "]")), aes(xintercept = mean))+
      theme_classic()+
      scale_x_continuous(limits = c(lims, lime), breaks = scales::pretty_breaks(n = 5))
  }else{
    data = df %>% pivot_longer(cols = contains(post_parameters)) %>% drop_na()
    data %>% 
      ggplot(aes(x = value, fill = posterior))+
      geom_histogram(alpha = 0.5, position="identity")+
      facet_wrap(~name, nrow = nsubs/5, ncol = 5, scales = "free")+
      theme_classic()+
      #coord_cartesian(xlim = c(min(data$value)-sd(data1$value),mean(data$value)+2*sd(data1$value)))+
      scale_x_continuous(limits =c(lims, lime), breaks = scales::pretty_breaks(n = 5))
    
    }
}



ppcheck = function(fit, df, n, prior = FALSE){
  
  percept = as.matrix(df %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL))
  
  expectPain = as.matrix(df %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL))
  percept_bin = as.matrix(df %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL))
  
  variables = c("percept","percept_bin","pred")
  
  if(!prior){
    variables = paste0("post_",variables)
    draws = as_draws_df(fit$draws(variables = variables))
    
    
    drawss <- lapply(1:length(unique(df$id)), function(i) {
      draws1 <- draws %>% select(matches(paste0("post_percept\\[\\d+\\,", i, "\\]"))) %>% drop_na()
      colnames(draws1) <- paste0("post_percept[", 1:160, "]")
      return(draws1)
    })
    
  }else{
    variables = paste0("prior_",variables)
    draws = as_draws_df(fit$draws(variables = variables))
    
    
    drawss <- lapply(1:length(unique(df$id)), function(i) {
      draws1 <- draws %>% select(matches(paste0("prior_percept\\[\\d+\\,", i, "\\]"))) %>% drop_na()
      colnames(draws1) <- paste0("prior_percept[", 1:160, "]")
      return(draws1)
    })
  }
    
    
    drawss <- do.call(rbind, drawss)
  
  
  plot = data.frame(t((drawss[rbinom(n,nrow(drawss),0.5),]) )) %>% 
    pivot_longer(everything()) %>% 
    ggplot(aes(x = value, group = name))+
    geom_density(alpha = 0.5, col = "lightblue")+
    theme_classic()+
    geom_density(data = df %>% mutate(name = NA), aes(x = percept),col = "red", size = 1.5)+
    scale_y_continuous(limits = c(0,5), breaks = scales::pretty_breaks(n = 5))
    
  return(plot)
}


