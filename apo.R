
library(doParallel)
library(splines)

#-----------------------------------------------------------------------------------------------------------------
# sample independent normal data by specifying treatment arm means (standard deviation = 1)
# returns a matrix where row denotes treatment arm (arm 1 is control) and column is patient number
sample.data <- function(n_per_arm=100, arm_mus)
{
  num_arms = length(arm_mus)
  means = rep( arm_mus, each=n_per_arm )
  draws = rnorm(n_per_arm*num_arms, mean = means, sd=1)
  patient_data <- matrix(draws, nrow = n_per_arm, ncol = num_arms) #initialize matrix
  return(patient_data)
}

#-----------------------------------------------------------------------------------------------------------------
# calculate the proportion of times that index 'ii' occurs in vector 'vec'
# made to be plugged into the 'apply function'
check.index = function(ii, vec){
	return( mean( vec == ii ) )
}

#-----------------------------------------------------------------------------------------------------------------
# simulate a single 'drop the losers' trial
# returns the result of the test and the chosen winner treatment arm 
simulate.trial <- function(n_per_arm_I, arm_mus, n_per_arm_II, threshold=0.9, post_sample_size=1000, delta=0, epsilon=0)
{
  num_arms = length(arm_mus)
  
  #simulate data via sample.data function
  patient_data <- sample.data(n_per_arm=n_per_arm_I, arm_mus=arm_mus)
  
  SDS = rep( apply(patient_data,MARGIN=2,FUN=sd), each=post_sample_size )
  MEANS = rep( apply(patient_data,MARGIN=2,FUN=mean), each=post_sample_size )
  draws = rt(n=post_sample_size*num_arms, df=n_per_arm_I-1) * SDS/sqrt(n_per_arm_I) + MEANS
  mu_post_draws <- matrix(draws, nrow = post_sample_size, ncol = num_arms)
  
  #find PR(best arm) for each arm
  raw_best_arm <- apply(mu_post_draws,1,which.max) #gives index of "best" arm via largest post mu for each iteration
  pr_best_arm = sapply(1:num_arms,FUN=check.index,vec=raw_best_arm)
  
  #determine winner arm and save its index
  winner_arm_index <- which.max(pr_best_arm) #the urge to label this "windex" is overwhelming
  
  #early stopping rule
  result_string = ""
  stop_for_futility_r1 = (winner_arm_index == 1) # default: if best arm is control
  stop_for_futility_r2 = FALSE
  
  if( stop_for_futility_r1 ){
    result_string = "early stop control is winner"
  } else{
    pr_promising_effect = mean(mu_post_draws[,winner_arm_index] - mu_post_draws[,1] >= delta)
    stop_for_futility_r2 = (pr_promising_effect <= epsilon)
  }
  
  if( stop_for_futility_r2 ){
    result_string = "early stop probability rule"
  }
  
  if( !stop_for_futility_r1 && !stop_for_futility_r2 ){
    
    #simulate data from phase II (winner arm and ctrl only)
    patient_data_phaseII <- sample.data(n_per_arm=n_per_arm_II, arm_mus = c(arm_mus[1], arm_mus[winner_arm_index]))
    
    control_dat <- c(patient_data[,1],patient_data_phaseII[,1]) #ALL control data from phaseI&II
    tx_dat <- c(patient_data[,winner_arm_index],patient_data_phaseII[,2]) #ALL tx data from phaseI&II
    
    #posterior draws from mu from ctrl AND tx
    mu_post_ctrl <- rt(n=post_sample_size, df=(n_per_arm_I+n_per_arm_II-1))*sd(control_dat)/(sqrt(n_per_arm_I+n_per_arm_II)) + mean(control_dat)
    mu_post_tx <-   rt(n=post_sample_size, df=(n_per_arm_I+n_per_arm_II-1))*sd(tx_dat)/(sqrt(n_per_arm_I+n_per_arm_II)) + mean(tx_dat)
    
    #find PR(tx > ctrl + effect size) via mu's
    pr_tx_superior <- mean( mu_post_tx > mu_post_ctrl ) 
    
    if(pr_tx_superior >= threshold)
    {
      result_string = "success"
    } else if(pr_tx_superior < threshold)
    {
      result_string = "failure"
    }
  }
  
  output_object = list(winner=winner_arm_index, result=result_string)
  return( output_object )
}

#-----------------------------------------------------------------------------------------------------------------
# simulate a multiple 'drop the losers' trials
# returns the results of the tests as a vector and the chosen winner treatment arms as a vector
simulate.trial.multiple = function(simulation_runs, n_per_arm_I, arm_mus, n_per_arm_II, threshold=0.9, post_sample_size=1000, delta=0, epsilon=0){
  preout = sapply(
    rep(n_per_arm_I, simulation_runs), FUN=simulate.trial, arm_mus=arm_mus, 
		n_per_arm_II=n_per_arm_II, threshold=threshold, post_sample_size=post_sample_size,
		delta=delta, epsilon=epsilon
  )
  out = list(winners=unlist(preout[1,]), results=unlist(preout[2,]))
  return( out )
}

#-----------------------------------------------------------------------------------------------------------------
# short wrapper function to generate means that increase linearly (use for linear UI setting)
# first mu = 0
# last mu = best_arm_mu
create.linear.means = function(num_arms, best_arm_mu){
  return( ((1:num_arms)-1)*best_arm_mu/(num_arms-1) )
}

#-----------------------------------------------------------------------------------------------------------------
# numerically estimate the rejection threshold required to calibrate the 'target_alpha' type 1 error rate
# generates a search grid of thresholds over which error rates are estimated
# fits a polynomial regression model to estimate the target rejection threshold
# returns said threshold and the search grid used to estimate it
estimate.threshold = function(
  target_alpha = 0.05, num_arms = 5, n_per_arm_I = 50, n_per_arm_II = 50, post_sample_size=1000, dt_start = 0.1, th_start = 0.9, fill_it = 6, 
  shrink = 0.7, shrink_depth = 15, simulation_runs = 1000, seed=-1, delta=0, epsilon=0
){

  if(seed>0){
    set.seed(seed)
  }
  
  # --- error handling --- #
  
  if(shrink_depth < 6){
    shrink_depth = 6
    cat("INSUFFICIENT SHRINK DEPTH! \nDEFAULTING TO A VALUE OF 6 INSTEAD!\n")
  }
  
  # --- adaptive stepwise search --- #
  
  DT = dt_start
  ALPHAS = c()
  THRESHOLDS = c(th_start)
  clev = 1
  LEVELS = c(clev)
  sim_obj = simulate.trial.multiple(
    simulation_runs, n_per_arm_I, arm_mus=rep(0,num_arms), n_per_arm_II, threshold=THRESHOLDS[1], post_sample_size, delta=delta, epsilon=epsilon
  )
  test_outcomes = sim_obj$results
  
  ALPHAS[1] = mean(test_outcomes == "success")
  direction = if(ALPHAS[1] > target_alpha | ALPHAS[1]<0) 1 else -1 #if ALPHAS[1] == 1 then direction is 1
  cat(THRESHOLDS[1],ALPHAS[1],LEVELS[1],direction,"\n")
  
  k=1
  
  while( shrink_depth > clev ){
    k = k+1
    THRESHOLDS[k] = THRESHOLDS[k-1] + DT*direction
	
	if( THRESHOLDS[k] > 1 ){ # set anchor point if necessary
      THRESHOLDS[k] = 1.000001
	  ALPHAS[k] = 0
	  LEVELS[k] = clev
	  direction = -1
	  clev = clev + 1
	  cat(THRESHOLDS[k],ALPHAS[k],LEVELS[k],direction,"\n")
	  
	  k = k+1
	  THRESHOLDS[k] = THRESHOLDS[k-1] + DT*direction
	}
	
	if( THRESHOLDS[k] < 0 ){ # set anchor point if necessary
      THRESHOLDS[k] = -0.000001
	  ALPHAS[k] = 1
	  LEVELS[k] = clev
	  direction = 1
	  clev = clev + 1
	  cat(THRESHOLDS[k],ALPHAS[k],LEVELS[k],direction,"\n")
	  
	  k = k+1
	  THRESHOLDS[k] = THRESHOLDS[k-1] + DT*direction
	}
    
    sim_obj = simulate.trial.multiple(
      simulation_runs, n_per_arm_I, arm_mus=rep(0,num_arms), n_per_arm_II, threshold=THRESHOLDS[k], post_sample_size, delta=delta, epsilon=epsilon
	)
	test_outcomes = sim_obj$results
    
    ALPHAS[k] = mean(test_outcomes == "success")
	LEVELS[k] = clev
    if(ALPHAS[k] > target_alpha){
      if(direction == -1){ # if direction switches, decrease step size
        DT = DT * shrink
		clev = clev + 1
      }
      direction = 1
    } else{
      if(direction == 1){ # if direction switches, decrease step size
        DT = DT * shrink
		clev = clev + 1
      }
      direction = -1
    }
    cat(THRESHOLDS[k],ALPHAS[k],LEVELS[k],direction,"\n")
  }
  
  tab = data.frame(A=ALPHAS,T=THRESHOLDS,L=LEVELS)
  
  # --- snapshot selection --- #
  
  # select the relevant snapshot of the grid
  LG3 = tab[tab$L>2,]
  tab = tab[tab$T >= min(LG3$T),]
  tab = tab[tab$T <= max(LG3$T),]
  tabo = tab[order(tab$T),]
  
  # --- filling iterations --- #
  
  # fill in sparse points of the grid to diminish
  # leverage and improve fit
  tabo$diff = NA
  for(e in 1:fill_it){
    imax = 1
    cmax = -999
    for(i in 2:dim(tabo)[1]){
      tabo$diff[i] = tabo$T[i] - tabo$T[i-1]
      if(tabo$diff[i] > cmax){cmax=tabo$diff[i]; imax=i}
    }
    TH = mean(tabo$T[(imax-1):imax])
    sim_obj = simulate.trial.multiple(
      simulation_runs, n_per_arm_I, arm_mus=rep(0,num_arms), n_per_arm_II, threshold=TH, post_sample_size, delta=delta, epsilon=epsilon
	)
	test_outcomes = sim_obj$results

    ALPH = mean(test_outcomes == "success")
    tabo = rbind(tabo,data.frame(A=ALPH,T=TH,L=NA,diff=NA))
    tabo = tabo[order(tabo$T),]
	cat(TH,ALPH,"NA NA\n")
  }
  
  # --- threshold estimation --- #
  
  A1 = tabo$A
  A2 = A1^2
  A3 = A1^3
  
  m3 = lm(tabo$T ~ A1 + A2 + A3)
  m2 = lm(tabo$T ~ A1 + A2)
  m1 = lm(tabo$T ~ A1)

  pv3 = (anova(m3,m2))$Pr[2]
  pv2 = (anova(m2,m1))$Pr[2]

  tmod = m3
  if(pv3>0.05){
    tmod = m2
    if(pv2>0.05){
      tmod = m1
    }
  }

  cth = predict(tmod,newdata=data.frame(A1=target_alpha,A2=target_alpha^2,A3=target_alpha^3))
  
  #plot(tabo$T,tabo$A)
  #abline(h=target_alpha,col="blue")
  #abline(v=cth,col="blue")
  #anew = seq(from=min(tabo$A),to=max(tabo$A),by=0.001)
  #lines(
  #  predict(tmod,newdata=data.frame(A1=anew,A2=anew^2,A3=anew^3)),
  #  anew, col="red"
  #)
  
  return(list(TH=cth,grid=tabo[,c("A","T")]))
}

#-----------------------------------------------------------------------------------------------------------------
# creates frequency tables of test results and winner indices like the 'table function'
# used instead of 'table' in order to NOT omit labels that do not occur
# returns a breakdown matrix with a column for each winner treatment arm and a row for each test outcome
create.full.frequency.table = function(results,winners,ntx){
  dat = data.frame(r=results,w=winners)
  mat = matrix(ncol=ntx,nrow=4)
  for(i in 1:ntx){
    sdat = dat[dat$w==i,]
	mat[1,i] = sum(sdat$r=="early stop control is winner")
	mat[2,i] = sum(sdat$r=="early stop probability rule")
	mat[3,i] = sum(sdat$r=="failure")
	mat[4,i] = sum(sdat$r=="success")
  }
  colnames(mat) = paste("Tx",1:ntx,sep="")
  row.names(mat) = c("early stop control is winner","early stop probability rule","failure","success")
  return( mat )
}

#-----------------------------------------------------------------------------------------------------------------
# estimate power of a specific alternative hypothesis scenario
# requires specification of treatment means (treatment 1 is control) and a target rejection threshold
# returns the power estimate and a breakdown matrix containing the test result by winner arm frequency table
get.power <- function(
  threshold, arm_mus, n_per_arm_I = 50, n_per_arm_II = 50, post_sample_size = 1000,
  simulation_runs_H1 = 5000, seed=-1, delta=0, epsilon=0
)
{

  if(seed>0){
    set.seed(seed)
  }
  h1_sim_obj <- simulate.trial.multiple(
    simulation_runs_H1, n_per_arm_I=n_per_arm_I, arm_mus=arm_mus, n_per_arm_II=n_per_arm_II, threshold=threshold, post_sample_size=post_sample_size,
    delta=delta, epsilon=epsilon
  )
  h1_test_outcomes = h1_sim_obj$results
  
  breakdown_matrix = create.full.frequency.table(h1_sim_obj$results, h1_sim_obj$winners, ntx=length(arm_mus))
  
  power <- mean(h1_test_outcomes=="success")
  return( list(pow=power, th=threshold, breakdown=breakdown_matrix) )
  
}

#-----------------------------------------------------------------------------------------------------------------
# find the best df value when smoothing the power-curve via cubic splines
# only works when x are fairly evenly spaced and sorted
find.best.df = function(x,y){
	nh = ceiling(length(x)*0.5)
	
	odn = (1:nh)*2-1
	odn = odn[odn <= length(x)]
	evn = (1:nh)*2
	evn = evn[evn <= length(x)]
	x1 = x[odn]
	x2 = x[evn]
	y1 = y[odn]
	y2 = y[evn]
	
	odn_b = odn[odn>min(evn) & odn<max(evn)]
	evn_b = evn[evn>min(odn) & evn<max(odn)]
	x1b = x[odn_b]
	x2b = x[evn_b]
	y1b = y[odn_b]
	y2b = y[evn_b]
	
	dfs = c(3,4,5,6)
	errs = c()
	for(i in 1:length(dfs)){
		dfc = dfs[i]
		
		modstr1 = paste("y1 ~ bs(x1,df=",dfc,")",sep="")
		modstr2 = paste("y2 ~ bs(x2,df=",dfc,")",sep="")
		fit1 = lm(as.formula(modstr1))
		fit2 = lm(as.formula(modstr2))
		
		pred2b = predict(fit1,newdata=data.frame(x1=x2b))
		pred1b = predict(fit2,newdata=data.frame(x2=x1b))
		
		#pred2b[pred2b<0] = 0
		#pred1b[pred1b<0] = 0
		
		err = sqrt(sum((
			c( pred2b, pred1b ) - c(y2b,y1b) 
		)^2))
		
		errs[i] = err
	}
	dat = data.frame(dfs,errs)
	dat = dat[dat$errs==min(errs),]
	out = dat$dfs[1]
	
	return( out )
}

#-----------------------------------------------------------------------------------------------------------------
# provide outcome decompositions for a list of model breakdowns
outcome.decomposition = function(breakdown_list, num_per_arm_I_schedule, arm_mus){
  NP = length(num_per_arm_I_schedule)
  NT = dim(breakdown_list[[1]])[2]
  i = 1
  outm = matrix(NA,ncol=NT,nrow=NP)
  
  best_treatment_columns = ( arm_mus==max(arm_mus) )
  true_treatment_columns = ( arm_mus > arm_mus[1] )
  
  best_treatment = c()
  true_treatment = c()
  futility_control_winner = c()
  futility_probability_rule = c()
  
  for(i in 1:NP){
    kpows = c()
    for(k in 1:NT){
      kpows[k] = breakdown_list[[i]][4,k] / sum(breakdown_list[[i]])
    }
    outm[i,] = kpows
	best_treatment[i] = sum(breakdown_list[[i]][4,best_treatment_columns]) / sum(breakdown_list[[i]])
	true_treatment[i] = sum(breakdown_list[[i]][4,true_treatment_columns]) / sum(breakdown_list[[i]])
    futility_control_winner[i] = sum(breakdown_list[[i]][1,]) / sum(breakdown_list[[i]])
    futility_probability_rule[i] = sum(breakdown_list[[i]][2,]) / sum(breakdown_list[[i]])
  }
  futility_overall = futility_control_winner + futility_probability_rule
  colnames(outm) = colnames(breakdown_list[[1]])
  out = data.frame(
    n.per.arm.I=num_per_arm_I_schedule, outm, Txs.Best.Arms=best_treatment, Txs.True.Effect=true_treatment,
	futility=futility_overall, futility.ctrl.win=futility_control_winner, futility.prob.rule=futility_probability_rule
  )
  return( out )
}


#-----------------------------------------------------------------------------------------------------------------
# two small utility functions used for boostrapping simulation error
success.rate = function(rlist){
  return( mean(rlist=="success") )
}
draw.bootstrap = function(foo,dat,size){
  return( sample(dat,size,replace=TRUE) )
}

#-----------------------------------------------------------------------------------------------------------------
# top level simulation function
# assumes a fixed total sample size of patients can be allocated
# finds the optimal, balanced Stage I/II sample distribution of patients
# helps to answer when it is *best* to do interim analysis. 
optimize.power <- function(
  target_alpha = 0.05, arm_mus, n_patients_TOTAL = 200, post_sample_size=1000,
  simulation_runs_H0 = 1000, simulation_runs_H1 = 5000, threads=2, num_checked_sizes=20,
  seed=-1, delta=0, epsilon=0
){
  num_arms = length(arm_mus)
  
  seed_th_base = -1
  seed_pw_base = -1
  if(seed>0){
    seed_th_base = seed
	seed_pw_base = 9857 + seed
  }
  
  # make a coarse estimate
  posterior_power <- c() #initialize vector
  # mathematically needs to be at >=2 pats per arm in phase I
  if(n_patients_TOTAL<=num_arms*2){return("Need more than 2*(number of arms) patients!")}
  
  # start with 2 pats per arm in phase I 
  num_per_arm_I_schedule <- floor(seq(from=2, to=(n_patients_TOTAL/num_arms)-2, length.out=num_checked_sizes))
  
  # first get thresholds
  cl <- makeCluster(threads, outfile="")
  registerDoParallel(cl)
  threshold_list = foreach(
    i=1:num_checked_sizes, .verbose=FALSE, 
    .export=c(
      "simulate.trial", "sample.data", "check.index", "simulate.trial.multiple", "estimate.threshold", "get.power", "create.full.frequency.table"
    )
  ) %dopar% {
  
    current_seed = -1
    if(seed_th_base>0){
	  current_seed = seed_th_base+i*113
	}
	eout = estimate.threshold(
      target_alpha = target_alpha, num_arms = num_arms, n_per_arm_I = num_per_arm_I_schedule[i], 
      n_per_arm_II = floor((n_patients_TOTAL-(num_per_arm_I_schedule[i]*num_arms))/2), 
      post_sample_size=post_sample_size, simulation_runs = simulation_runs_H0, seed=current_seed,
      delta=delta, epsilon=epsilon
	)
    threshold = eout$TH
    threshold
  }
  stopCluster(cl)
  
  threshold_phase_II = unlist(threshold_list)
  
  th_fit <- lm(threshold_phase_II~bs(num_per_arm_I_schedule,df=4,degree=3))
  threshold_estimated_phase_II = predict(th_fit)
  pred_n1 <- 2:(floor(n_patients_TOTAL/num_arms)-2)
  pred_th <- predict(th_fit,newdata = data.frame(num_per_arm_I_schedule=pred_n1))
  
  cl <- makeCluster(threads, outfile="")
  registerDoParallel(cl)
  posterior_power_objects = foreach(
    i=1:num_checked_sizes, .verbose=FALSE, 
    .export=c(
      "simulate.trial", "sample.data", "check.index", "simulate.trial.multiple", "estimate.threshold", "get.power", "create.full.frequency.table"
    )
  ) %dopar% {
     
	  current_seed = -1
      if(seed_pw_base>0){
	    current_seed = seed_pw_base+i*157
	  }
	  
      outp = get.power(
        threshold = threshold_estimated_phase_II[i],
        arm_mus = arm_mus, n_per_arm_I = num_per_arm_I_schedule[i], 
        n_per_arm_II = floor((n_patients_TOTAL-(num_per_arm_I_schedule[i]*num_arms))/2), 
        post_sample_size=post_sample_size, simulation_runs_H1=simulation_runs_H1,
        seed=current_seed, delta=delta, epsilon=epsilon
    )
    
    outp
  }
  stopCluster(cl)
  
  posterior_power = c()
  breakdown_list = list()
  for(i in 1:num_checked_sizes){
    posterior_power[i] = posterior_power_objects[[i]]$pow
    breakdown_list[[i]] = posterior_power_objects[[i]]$breakdown
  }
  
  decomp_raw = outcome.decomposition(breakdown_list, num_per_arm_I_schedule, arm_mus)
  
  decomp_smooth = matrix(NA,ncol=dim(decomp_raw)[2],nrow=length(pred_n1))
  for(i in 2:dim(decomp_raw)[2]){
    tpostdec = decomp_raw[,i]
	#if(length(unique(tpostdec))>2){
      tbdf = find.best.df(x=num_per_arm_I_schedule, y=tpostdec)  
      tpw_fit = lm(tpostdec ~ bs(num_per_arm_I_schedule,df=tbdf) )
      tpred_pow = predict(tpw_fit,newdata=data.frame(num_per_arm_I_schedule=pred_n1))
      tpred_pow[tpred_pow<0] = 0
      decomp_smooth[,i] = tpred_pow
	#}
  }
  decomp_smooth[,1] = pred_n1
  colnames(decomp_smooth) = colnames(decomp_raw)
  decomp_smooth = data.frame(decomp_smooth)
  
  bdf = find.best.df(x=num_per_arm_I_schedule, y=posterior_power)  
  pw_fit = lm(posterior_power ~ bs(num_per_arm_I_schedule,df=bdf) )
  pred_pow = predict(pw_fit,newdata=data.frame(num_per_arm_I_schedule=pred_n1))
  best_index = which.max(pred_pow)
  best_n_per_arm_I = pred_n1[best_index]
  
  raw_sim = data.frame(
    n.per.arm.I=num_per_arm_I_schedule,threshold=threshold_phase_II,
    threshold.smoothed=threshold_estimated_phase_II,power=posterior_power
  )
  pred_n2 <- floor((n_patients_TOTAL - pred_n1*num_arms)/2)
  pred_nt <- pred_n1*num_arms + pred_n2*2
  predicted_performance = data.frame(
    power=pred_pow, 
    power.low=pred_pow - qnorm(p=0.95) * sqrt(pred_pow*(1-pred_pow)/simulation_runs_H1),
    power.high=pred_pow + qnorm(p=0.95) * sqrt(pred_pow*(1-pred_pow)/simulation_runs_H1),
    n.per.arm.I=pred_n1, n.per.arm.II=pred_n2, 
    n.total=pred_nt, threshold=pred_th,
    pr.futility = decomp_smooth$futility,
    pr.best.tx.success = decomp_smooth$Txs.Best.Arms
  )
  
  best_row = predicted_performance[predicted_performance$power==max(predicted_performance$power),]
  best_row$alpha.low = target_alpha - qnorm(p=0.95) * sqrt(target_alpha*(1-target_alpha)/simulation_runs_H0)
  best_row$alpha.high = target_alpha + qnorm(p=0.95) * sqrt(target_alpha*(1-target_alpha)/simulation_runs_H0)
  
  return( list(
    results.all = predicted_performance,
    results.best = best_row,
    raw.sim = raw_sim,
    raw.decomp = decomp_raw,
    decomposition.smooth = decomp_smooth,
    true.mus = arm_mus
  ))
}

#-----------------------------------------------------------------------------------------------------------------
# generates diagnostic plots of simulation result objects returned by the 'optimize.power' function
plot.sim = function(sim_obj){
  my_red = rgb(red=255, green=48, blue=66, maxColorValue=255)
  my_blue = rgb(red=39, green=103, blue=255, maxColorValue=255)
  my_grey = rgb(red=133, green=137, blue=138, maxColorValue=255)
  par(mfrow=c(1,2))
  
  plot(
    sim_obj$raw.sim$n.per.arm.I, sim_obj$raw.sim$threshold, 
    xlab="Patients per arm (Stage 1)", ylab="Threshold",
    main="Stage 2 Rejection Threshold"
  )
  lines(sim_obj$raw.sim$n.per.arm.I, sim_obj$raw.sim$threshold, col=my_red)
  lines(sim_obj$results.all$n.per.arm.I, sim_obj$results.all$threshold, col=my_blue)
  
  legend(
    "bottomright",
    legend=c("estimates","joined estimates","smoothed fit", "best N Stage 1"),
    lty = c(0,1,1,3),
    lwd = c(0,1,1,2),
    pch = c(1,NA,NA,NA),
    col = c("black",my_red,my_blue,my_grey),
    cex=1
  )
  
  plot(
    sim_obj$raw.sim$n.per.arm.I, sim_obj$raw.sim$power, 
    xlab="Patients per arm (Stage 1)", ylab="Power",
    main="Overall Power"
  )
  lines(sim_obj$raw.sim$n.per.arm.I, sim_obj$raw.sim$power, col=my_red)
  lines(sim_obj$results.all$n.per.arm.I, sim_obj$results.all$power, col=my_blue)
  abline(v=sim_obj$results.best$n.per.arm.I, lty=3, col=my_grey, lwd=2)
}

