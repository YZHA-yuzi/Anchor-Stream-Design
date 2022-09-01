######################################################################
#                                                                    #
#                                                                    #
#  Self-defined functions used for computing 
#  point estimates and interval estimates proposed in the manuscript #
#                                                                    #
#                                                                    #
######################################################################

## A function to compute aggregated cell counts from individual-level data
get_nxbar_obs <- function(dat.obs, Ntot){
  
  ## # selected in both S1 and S2 that are CV negative
  n11.noCV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 1 & dat.obs$case == 0)
  ## # selected in S1 but not in S2 that are CV negative
  n10.noCV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 0 & dat.obs$case == 0)
  ## # selected in S2 but not in S1 that are CV negative
  n01.noCV = sum(dat.obs$y1 == 0 & dat.obs$y2 == 1 & dat.obs$case == 0)
  
  ## # selected in both S1 and S2 that are CV positive
  n11.CV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 1 & dat.obs$case == 1)
  ## # selected in S1 but not in S2 that are CV positive
  n10.CV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 0 & dat.obs$case == 1)
  ## # selected in S2 but not in S1 that are CV positive
  n01.CV = sum(dat.obs$y1 == 0 & dat.obs$y2 == 1 & dat.obs$case == 1)
  
  ## rename variables based on the notation used in Table 2 of anchor stream paper
  n1 = n11.noCV; n2 = n11.CV
  n3 = n10.noCV; n4 = n10.CV
  n5 = n01.noCV; n6 = n01.CV
  n7 = Ntot - n1 - n2 - n3 - n4 - n5 - n6 
  
  xbar11 = mean(subset(dat.obs, y1 == 1 & y2 == 1)$x)
  xbar01 = mean(subset(dat.obs, y1 == 0 & y2 == 1)$x)
  xbar10 = mean(subset(dat.obs, y1 == 1 & y2 == 0)$x)
  xbar2dot = mean(subset(dat.obs, y2 == 1)$x)
  
  return(list(nvec = c(n1, n2, n3, n4, n5, n6, n7),
              xbarvec = c(xbar11, xbar10, xbar01, xbar2dot)))
}


## A function to estimate case count using N_RS, Nhat_psi and Nhat_psihat_star 
## based on individual-level data or aggregated cell counts 
AnchorStream_CaseCount <- function(dat.obs = NULL, Ntot, p2, 
                                   num.post = 10000, 
                                   data.type = "individual",
                                   cellcounts.vec = NULL,
                                   seed){
  
  if(is.null(dat.obs) & data.type == "individual"){
    stop('please provide a data frame which continas individual-level data')
  }
  if(!is.null(dat.obs) & data.type == "individual"){
    inputnames <- colnames(dat.obs)
    if(sum(inputnames %in% c("y1", "y2", "case", "sympt")) != 4){
      stop('please make sure that capture indicators, case status indicator, and symptom idicator are named as y1, y2, case, and sympt')
    }
  }
  if(is.null(cellcounts.vec) & data.type == "aggregated"){
    stop('please provide a vector of cell counts')
  }
  
  ## INPUT: dat.obs: observed data frame
  ## y1: capture indicator of S1
  ## y2: capture indicator of S2
  ## case: 1 tested positive; 0 otherwise
  ## sympt: 1 have symptom; 0 otherwise
  ## INPUT: data.type = {"individual", "aggregated"}
  ## "individual" represents individual-level data
  ## "aggregated" represents cell counts as presented in Table 1 in the manuscript
  ## INPUT: Ntot: total population size
  ## INPUT: seed: random seed
  ## INPUT: p2 (known selection probability into S2)
  ## INPUT: cellcounts.vec (a vector of cell counts if data.type == aggregated)
  
  set.seed(seed)
  
  if(data.type == "individual"){
    ## # cases selected into S2
    ## # cases tested positive in S2 
    num.sel.2 = sum(dat.obs$y2 == 1)
    num.pos.2 = sum(dat.obs$y2 == 1 & dat.obs$case == 1)
    
    ## # selected in both S1 and S2 that are CV negative
    n11.noCV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 1 & dat.obs$case == 0)
    ## # selected in S1 but not in S2 that are CV negative
    n10.noCV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 0 & dat.obs$case == 0)
    ## # selected in S2 but not in S1 that are CV negative
    n01.noCV = sum(dat.obs$y1 == 0 & dat.obs$y2 == 1 & dat.obs$case == 0)
    
    ## # selected in both S1 and S2 that are CV positive
    n11.CV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 1 & dat.obs$case == 1)
    ## # selected in S1 but not in S2 that are CV positive
    n10.CV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 0 & dat.obs$case == 1)
    ## # selected in S2 but not in S1 that are CV positive
    n01.CV = sum(dat.obs$y1 == 0 & dat.obs$y2 == 1 & dat.obs$case == 1)
    n11add.CV = max(n11.CV, 0.5)
    n01add.CV = max(n01.CV, 0.5)
    n10add.CV = max(n10.CV, 0.5)
    
    nc.CV = n11.CV + n10.CV + n01.CV
    n11 = n11.CV + n11.noCV
    n10 = n10.CV + n10.noCV
    n01 = n01.CV + n01.noCV
    nc = n11 + n10 + n01
    
    ## rename variables based on the notation used in Table 2 of anchor stream paper
    npos1 = n11.CV + n10.CV; npos2 = n11.CV + n01.CV 
    n1 = n11.noCV; n2 = n11.CV
    n3 = n10.noCV; n4 = n10.CV
    n5 = n01.noCV; n6 = n01.CV
    n7 = Ntot - n1 - n2 - n3 - n4 - n5 - n6 
  }

  if(data.type == "aggregated"){
    
    n1 = cellcounts.vec[1]; n2 = cellcounts.vec[2]
    n3 = cellcounts.vec[3]; n4 = cellcounts.vec[4]
    n5 = cellcounts.vec[5]; n6 = cellcounts.vec[6]
    n7 = cellcounts.vec[7]
    
    num.sel.2 = n1 + n2 + n5 + n6
    num.pos.2 = n2 + n6
    
    n11.noCV = n1; n11.CV = n2
    n10.noCV = n3; n10.CV = n4
    n01.noCV = n5; n01.CV = n6
    
    n11add.CV = max(n11.CV, 0.5)
    n01add.CV = max(n01.CV, 0.5)
    n10add.CV = max(n10.CV, 0.5)
    
    nc.CV = n11.CV + n10.CV + n01.CV
    n11 = n11.CV + n11.noCV
    n10 = n10.CV + n10.noCV
    n01 = n01.CV + n01.noCV
    nc = n11 + n10 + n01
    npos1 = n11.CV + n10.CV; npos2 = n11.CV + n01.CV 
    
  }
  
  pihat.S1 = (n2 + n4)/(n1 + n2 + n3 + n4)
  pihat.S1bar = (n6 + n7*n6/(n5+n6))/(n5+n6+n7)
  psihat = (n1+n2+n5+n6)/Ntot # pr(selected in S2)
  phihat = (n1+n2+n3+n4)/Ntot # pr(selected in S1)
  
  phat.CV = pihat.S1*phihat + pihat.S1bar*(1-phihat)
  ### Nhat_psihat_star ##
  Nhat.CV = phat.CV*Ntot
  
  # Get Nhat and its SE based on the RS in Stream 2
  phat2 = num.pos.2/num.sel.2
  NhatUse2Only = phat2*Ntot
  phat2 = max(num.pos.2, .5)/num.sel.2
  Vphat2 = phat2*(1-phat2)/num.sel.2
  VNhatUse2Only = Ntot^2*Vphat2
  SEhatNhatUse2Only = sqrt(VNhatUse2Only)
  # Get estimated SE for RS estimator with FPC included
  fpcCochran = min((Ntot - num.sel.2)*num.sel.2/(Ntot*(num.sel.2-1)), 1)
  VNhatUse2OnlyWithFPC = fpcCochran*VNhatUse2Only
  SEhatNhatUse2OnlyWithFPC = sqrt(VNhatUse2OnlyWithFPC)
  
  CI.wald.RS.FPC = NhatUse2Only + c(-1.96, 1.96)*SEhatNhatUse2OnlyWithFPC
  CI.wald.RS.FPC[1] = max(nc.CV, CI.wald.RS.FPC[1])
  
  ## Jeffreys interval based only on Stream 2
  alphpost = num.pos.2 + 0.5
  betpost = num.sel.2 - num.pos.2 + 0.5
  
  LL_JeffreysForP = qbeta(0.025, alphpost, betpost)
  UL_JeffreysForP = qbeta(0.975, alphpost, betpost)
  
  if(num.pos.2 == 0){ LL_JeffreysForP = 0 }
  if(num.pos.2 == num.sel.2){ UL_JeffreysForP = 1 }
  LL_Jeffreys = Ntot*LL_JeffreysForP
  UL_Jeffreys = Ntot*UL_JeffreysForP
  
  ## Adjust the Jeffreys CI to account for FPC
  phat2Only = NhatUse2Only/Ntot
  a = sqrt(fpcCochran)
  b = phat2Only*(1-a)
  LL_JeffreysForPFPC = a*LL_JeffreysForP + b
  UL_JeffreysForPFPC = a*UL_JeffreysForP + b
  if(num.pos.2 == 0){ LL_JeffreysForPFPC = 0 }
  if(num.pos.2 == num.sel.2){ UL_JeffreysForPFPC = 1 }
  LL_JeffreysFPC = Ntot*LL_JeffreysForPFPC
  UL_JeffreysFPC = Ntot*UL_JeffreysForPFPC
  
  ## Allow Stream 2 Only CI to account for observed nc 
  LL_JeffreysFPCKnowNc = max(nc.CV, LL_JeffreysFPC)
  UL_JeffreysFPCKnowNc = max(nc.CV, UL_JeffreysFPC)
  
  CI.RS.Jeffreys = c(LL_JeffreysFPCKnowNc, UL_JeffreysFPCKnowNc)
  
  
  ### MLE of N in eqn. (2) Nhat_psi and 
  ### variance estimator for Nhat_psi in eqn. (3) ##
  NhatKnownPsiML = n11.CV + n10.CV + n01.CV/p2
  SEhatNhatKnownPsiML = sqrt(n01add.CV*(1-p2)/(p2^2))
  CI.wald.Nhatpsi = NhatKnownPsiML + c(-1.96, 1.96)*SEhatNhatKnownPsiML
  CI.wald.Nhatpsi[1] = max(nc.CV, CI.wald.Nhatpsi[1])
  
  
  ## variance estimator of LP estimator ##
  VhatLP12 = (n11add.CV + n10add.CV)*(n11add.CV + n01add.CV)*
    n10add.CV*n01add.CV/(n11add.CV^3)
  VNhatWTDavgF = 1/(1/VNhatUse2OnlyWithFPC + 1/VhatLP12)
  SEhatNhatWTDavgF = sqrt(VNhatWTDavgF)
  LLEstimatedPsiMLnew3 = max(nc.CV, Nhat.CV - 1.96*SEhatNhatWTDavgF)
  ULEstimatedPsiMLnew3 = Nhat.CV + 1.96*SEhatNhatWTDavgF
  
  CI.wald.Nhatpsihatstar = c(LLEstimatedPsiMLnew3, ULEstimatedPsiMLnew3)
  
  ### get Dirichlet-multinomial-based credible intervals (see Appendix A) ###
  get_credint <- function(Ntot, dat.obs, 
                          num.post, data.type,
                          cellcounts.vec){
    ### INPUT: 
    # dat.obs: observed data frame 
    # num.post: number of samples drawn from Dirichlet posterior 
    # data.type: {"individual", "aggregated}
    
    if(data.type == "individual"){
      ## # cases selected into S2
      num.sel.2 = sum(dat.obs$y2 == 1)
      ## # number of cases selected into S2 
      num.pos.2 = sum(dat.obs$y2 == 1 & dat.obs$case == 1)
      
      ## # selected in both S1 and S2 that are CV negative
      n11.noCV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 1 & dat.obs$case == 0)
      ## # selected in S1 but not in S2 that are CV negative
      n10.noCV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 0 & dat.obs$case == 0)
      ## # selected in S2 but not in S1 that are CV negative
      n01.noCV = sum(dat.obs$y1 == 0 & dat.obs$y2 == 1 & dat.obs$case == 0)
      
      ## # selected in both S1 and S2 that are CV positive
      n11.CV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 1 & dat.obs$case == 1)
      ## # selected in S1 but not in S2 that are CV positive
      n10.CV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 0 & dat.obs$case == 1)
      ## # selected in S2 but not in S1 that are CV positive
      n01.CV = sum(dat.obs$y1 == 0 & dat.obs$y2 == 1 & dat.obs$case == 1)
      nc.CV = n11.CV + n10.CV + n01.CV
      n11add.CV = max(n11.CV, 0.5)
      n01add.CV = max(n01.CV, 0.5)
      n10add.CV = max(n10.CV, 0.5)
      
      npos1 = n11.CV + n10.CV; npos2 = n11.CV + n01.CV 
      n1 = n11.noCV; n2 = n11.CV
      n3 = n10.noCV; n4 = n10.CV
      n5 = n01.noCV; n6 = n01.CV
      n7 = Ntot - n1 - n2 - n3 - n4 - n5 - n6 
    }

    if(data.type == "aggregated"){
      
      n1 = cellcounts.vec[1]; n2 = cellcounts.vec[2]
      n3 = cellcounts.vec[3]; n4 = cellcounts.vec[4]
      n5 = cellcounts.vec[5]; n6 = cellcounts.vec[6]
      n7 = cellcounts.vec[7]
      
      num.sel.2 = n1 + n2 + n5 + n6
      num.pos.2 = n2 + n6
      
      n11.noCV = n1; n11.CV = n2
      n10.noCV = n3; n10.CV = n4
      n01.noCV = n5; n01.CV = n6
      
      n11add.CV = max(n11.CV, 0.5)
      n01add.CV = max(n01.CV, 0.5)
      n10add.CV = max(n10.CV, 0.5)
      
      nc.CV = n11.CV + n10.CV + n01.CV
      n11 = n11.CV + n11.noCV
      n10 = n10.CV + n10.noCV
      n01 = n01.CV + n01.noCV
      nc = n11 + n10 + n01
      npos1 = n11.CV + n10.CV; npos2 = n11.CV + n01.CV 
      
    }
    
    pihat.S1 = (n2 + n4)/(n1 + n2 + n3 + n4)
    pihat.S1bar = (n6 + n7*n6/(n5+n6))/(n5+n6+n7)
    psihat = (n1+n2+n5+n6)/Ntot # pr(selected in S2)
    phihat = (n1+n2+n3+n4)/Ntot # pr(selected in S1)
    
    phat.CV = pihat.S1*phihat + pihat.S1bar*(1-phihat)
    Nhat.CV = phat.CV*Ntot
    
    ## estimated psi ##
    psiEstimate = num.sel.2/Ntot
    VhatNhatPsiMLB = n01add.CV*(1-psiEstimate)/(psiEstimate^2)
    
    ## generate a Dirichlet random variable with parameters with 
    # n1 + 0.5, ..., n3 + 0.5
    pstarcond <- rdirichlet(num.post, 
                            alpha = c(n11.CV, n10.CV, n01.CV) + 0.5) #num.post*3
    ptotcheck = rowSums(pstarcond)
    p1starpost = rowSums(pstarcond[,1:2])
    p1post = psiEstimate*p1starpost/(psiEstimate*p1starpost + pstarcond[,3])
    p2giv1starpost = pstarcond[,1]/(pstarcond[,1] + pstarcond[,2])
    
    p11post = p2giv1starpost*p1post
    p10post = (1-p2giv1starpost)*p1post
    p01post = psiEstimate*(1-p1post)
    pcpost = p11post + p10post + p01post 
    
    Nfirst = round(nc.CV/pcpost)
    Nfirst[Nfirst == 0] <- 1
    ncnew = rbinom(num.post, size = Nfirst, prob = pcpost)
    
    n11post = ncnew*pstarcond[,1]
    n10post = ncnew*pstarcond[,2]
    n01post = ncnew*pstarcond[,3]
    
    Npsipost = n11post + n10post + n01post/psiEstimate
    Npsipost = pmax(rep(nc.CV, num.post), Npsipost)
    NpsipostUnAdj = Npsipost
    
    ### For details please refer to eqn. (A.5) in Appendix A ###
    a = sqrt(VNhatWTDavgF/VhatNhatPsiMLB)
    b = Nhat.CV*(1-a)
    Npsipost = a*Npsipost + b
    Npsipost = pmax(rep(nc.CV, num.post), Npsipost)
    
    ### unadjusted CIs ###
    LL.unadj = quantile(NpsipostUnAdj, 0.025, na.rm = T)
    UL.unadj = quantile(NpsipostUnAdj, 0.975, na.rm = T)
    
    ### adjusted CIs ###
    ### compute variance of Chapman estimator ##
    VhatCH12 = (n11.CV + n10.CV + 1)*(n11.CV + n01.CV + 1)*(n10add.CV*n01add.CV)/
      ((n11.CV + 1)^2*(n11.CV + 2))
    ### compute variance of random sample estimator with FPC adjustment ###
    phat2 = num.pos.2/num.sel.2
    NhatUse2Only = phat2*Ntot
    phat2 = max(num.pos.2, .5)/num.sel.2
    Vphat2 = phat2*(1-phat2)/num.sel.2
    VNhatUse2Only = Ntot^2*Vphat2
    # Get estimated SE for RS estimator with FPC included
    fpcCochran = (Ntot - num.sel.2)*num.sel.2/(Ntot*(num.sel.2-1))
    fpcCochran = min(fpcCochran, 1)
    VNhatUse2OnlyWithFPC = fpcCochran*VNhatUse2Only
    SEhatNhatUse2OnlyWithFPC = sqrt(VNhatUse2OnlyWithFPC)
    
    LL.ab = quantile(Npsipost, 0.025, na.rm = T)
    UL.ab = quantile(Npsipost, 0.975, na.rm = T)
    SE.avg = sqrt((VhatCH12 + VNhatUse2OnlyWithFPC)/4)
    LL.avg = max(nc.CV, Nhat.CV - 1.96*SE.avg)
    UL.avg = Nhat.CV + 1.96*SE.avg
    LL.adj = min(LL.ab, (LL.ab+LL.avg)/2)
    UL.adj = max(UL.ab, (UL.ab+UL.avg)/2) 
    
    phatCOVID = Nhat.CV/Ntot
    if(phatCOVID < 0.2 | is.na(phatCOVID)){
      return(list(Nhatpsihatstar = c(LL.unadj, UL.unadj),
                  Nhatpsi = c(LL.unadj, UL.unadj)))
    }else{
      return(list(Nhatpsihatstar = c(LL.adj, UL.adj),
                  Nhatpsi = c(LL.unadj, UL.unadj)))
    }
  }
  re.credint <- get_credint(Ntot = Ntot, 
                            dat.obs = dat.obs, 
                            num.post = num.post,
                            data.type = data.type,
                            cellcounts.vec = cellcounts.vec)
  
  CI.Diri.Nhatpsihatstar <- re.credint$Nhatpsihatstar
  
  re.point = data.frame(Nhat.RS = NhatUse2Only,
                        Nhat.psi = NhatKnownPsiML,
                        Nhat.psihatstar = Nhat.CV)
  re.se = data.frame(SE.RS.FPC = SEhatNhatUse2OnlyWithFPC,
                     SE.Nhat.psi = SEhatNhatKnownPsiML,
                     SE.Nhat.psihatstar = SEhatNhatWTDavgF)
  re.CI = list(CI.RS.wald.FPC = CI.wald.RS.FPC,
               CI.RS.Jeffreys = CI.RS.Jeffreys,
               CI.Nhat.psi.wald = CI.wald.Nhatpsi,
               CI.Nhat.psi.Diri = as.numeric(re.credint$Nhatpsi),
               CI.Nhat.psihatstar.wald = CI.wald.Nhatpsihatstar,
               CI.Nhat.psihatstar.Diri = CI.Diri.Nhatpsihatstar)
  
  re = list(pointest = re.point,
            SE = re.se,
            CI = re.CI)
  ### OUTPUT: a list containing point estimates, standard error, and intervals
  return(re)
}


## A function to estimate continuous means using different estimators 
## for observed data or simulated bootstrap data ##
get_xbar_obs <- function(dat.obs, Ntot, 
                         nc.CV.sim, 
                         nc.noCV.sim,
                         n.vec.ori, 
                         xbar.vec.ori){
  
  # INPUT: dat.obs: observed/simulated bootstrap data; 
  # Ntot: total # of populations
  ## y1: capture indicator of S1
  ## y2: capture indicator of S2
  ## covid: 1 tested positive; 0 otherwise
  ## sympt: 1 have symptom; 0 otherwise
  ## nc.CV.sim: # of cases in observed data set
  ## nc.noCV.sim: # of non-cases in observed data set
  ## n.vec.ori: a vector containing original n1 - n7 values for FPC correction 
  ## xbar.vec.ori: a vector of estimates of x using the observed data 
  
  n1.ori = n.vec.ori[1]; n2.ori = n.vec.ori[2];
  n3.ori = n.vec.ori[3]; n4.ori = n.vec.ori[4];
  n5.ori = n.vec.ori[5]; n6.ori = n.vec.ori[6];
  n7.ori = n.vec.ori[7]
  
  xbar11.ori = xbar.vec.ori[1]; xbar10.ori = xbar.vec.ori[2]
  xbar01.ori = xbar.vec.ori[3]; xbar2dot.ori = xbar.vec.ori[4]
  
  ## # cases selected into S1
  ## # cases tested positive in S1 
  num.sel.1 = sum(dat.obs$y1 == 1)
  num.pos.1 = sum(dat.obs$y1 == 1 & dat.obs$case == 1)
  
  ## # cases selected into S2
  ## # cases tested positive in S2 
  num.sel.2 = sum(dat.obs$y2 == 1)
  num.pos.2 = sum(dat.obs$y2 == 1 & dat.obs$case == 1)
  
  ## # selected in both S1 and S2 that are CV negative
  n11.noCV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 1 & dat.obs$case == 0)
  ## # selected in S1 but not in S2 that are CV negative
  n10.noCV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 0 & dat.obs$case == 0)
  ## # selected in S2 but not in S1 that are CV negative
  n01.noCV = sum(dat.obs$y1 == 0 & dat.obs$y2 == 1 & dat.obs$case == 0)
  
  ## # selected in both S1 and S2 that are CV positive
  n11.CV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 1 & dat.obs$case == 1)
  ## # selected in S1 but not in S2 that are CV positive
  n10.CV = sum(dat.obs$y1 == 1 & dat.obs$y2 == 0 & dat.obs$case == 1)
  ## # selected in S2 but not in S1 that are CV positive
  n01.CV = sum(dat.obs$y1 == 0 & dat.obs$y2 == 1 & dat.obs$case == 1)
  
  n11 = n11.CV + n11.noCV
  n10 = n10.CV + n10.noCV
  n01 = n01.CV + n01.noCV
  nc = n11 + n10 + n01
  
  ## # selected in S1 and tested positive 
  n1.CV = sum(dat.obs$y1 == 1 & dat.obs$case == 1)
  
  ## compute mean of X among cases sampled by S2 
  # estimates of mean.true using random samples
  xbar2dot = mean(subset(dat.obs, y2 == 1)$x)  
  
  ### incorporate FPC adjustment ###
  fpc.2dot = (Ntot - num.sel.2)*(num.sel.2)/
    (Ntot*((num.sel.2)-1))
  fpc.2dot = min(fpc.2dot, 1)
  a.2dot = sqrt(fpc.2dot)
  b.2dot = xbar2dot.ori*(1-a.2dot)
  xbar2dot.adj = a.2dot*xbar2dot + b.2dot
  
  xbar01 = mean(subset(dat.obs, y1 == 0 & y2 == 1)$x)
  xbar1dot = mean(subset(dat.obs, y1 == 1)$x)
  
  xbar11 = mean(subset(dat.obs, y1 == 1 & y2 == 1)$x)
  xbar10 = mean(subset(dat.obs, y1 == 1 & y2 == 0)$x)
  
  dat.obs.CV = subset(dat.obs, case == 1)
  xbar2dot.CV = mean(subset(dat.obs.CV, y2 == 1)$x)
  xbar01.CV = mean(subset(dat.obs.CV, y1 == 0 & y2 == 1)$x)
  xbar1dot.CV = mean(subset(dat.obs.CV, y1 == 1)$x)
  
  dat.obs.noCV = subset(dat.obs, case == 0)
  xbar2dot.noCV = mean(subset(dat.obs.noCV, y2 == 1)$x)
  xbar01.noCV = mean(subset(dat.obs.noCV, y1 == 0 & y2 == 1)$x)
  xbar1dot.noCV = mean(subset(dat.obs.noCV, y1 == 1)$x)
  
  ## rename variables based on the notation used in Table 2 of anchor stream paper
  n1 = n11.noCV; n2 = n11.CV
  n3 = n10.noCV; n4 = n10.CV
  n5 = n01.noCV; n6 = n01.CV
  n7 = Ntot - n1 - n2 - n3 - n4 - n5 - n6 
  
  pihat.S1 = (n2 + n4)/(n1 + n2 + n3 + n4)
  pihat.S1bar = (n6 + n7*n6/(n5+n6))/(n5+n6+n7)
  psihat = (n1+n2+n5+n6)/Ntot # pr(selected in S2)
  phihat = (n1+n2+n3+n4)/Ntot # pr(selected in S1)
  
  phat.CV = pihat.S1*phihat + pihat.S1bar*(1-phihat)
  ## estimated # of COVID-19 cases
  Nhat.CV = phat.CV*Ntot
  ## # of observed COVID-19 cases 
  nc.CV = n2 + n4 + n6
  
  ## for bootstrap samples set Nhat.CV = nc.CV.sim if Nhat.CV < nc.CV.sim
  # nc.CV.sim is the number of cases observed in the original observed data
  # instead of the bootstrap samples
  Nhat.CV = max(Nhat.CV, nc.CV.sim)
  Nhat.noCV = Ntot - Nhat.CV
  Nhat.noCV = max(Nhat.noCV, nc.noCV.sim)
  
  ##### I. compute estimated xbar among populations ####
  # prob of selected into S1 (y1 = 1) 
  phat1dot = num.sel.1/Ntot
  phat0dot = 1 - phat1dot
  phat11 = n11/Ntot
  phat10 = (num.sel.1 - n11)/Ntot
  
  ## estimated mean of x (continuous variable of interest) from
  ## anchor method
  xbar.anchor = xbar1dot*phat1dot + xbar01*phat0dot
  
  #### compute estimated xbar with FPC adjustment for the interval estimation ###
  get_fpcCochran <- function(NtotStar, numsampledStar){
    fpcCochran = (NtotStar - numsampledStar)*numsampledStar/
      (NtotStar*(numsampledStar-1))
    fpcCochran = min(fpcCochran, 1)
    return(fpcCochran)
  }
  
  # I(a). compute estimated xbar with FPC adjustment on xbar01 piece #
  fpcCochran = get_fpcCochran(NtotStar = n5.ori + n6.ori + n7.ori,
                              numsampledStar = n5.ori + n6.ori)
  a = sqrt(fpcCochran); b = xbar01.ori*(1-a)
  xbar01adj = a*xbar01 + b
  
  # I(b) compute estimated xbar with FPC adjustment on xbar11 piece #
  fpcCochranB =  get_fpcCochran(NtotStar = n1.ori + n2.ori + n3.ori + n4.ori,
                                numsampledStar = n1.ori + n2.ori)
  aB = sqrt(fpcCochranB); bB = xbar11.ori*(1-aB)
  xbar11adj = aB*xbar11 + bB
  
  # I(c) compute estimated xbar with FPC adjustment on xbar10 piece #
  fpcCochranC =  get_fpcCochran(NtotStar = n1.ori + n2.ori + n3.ori + n4.ori,
                                numsampledStar = n3.ori + n4.ori)
  aC = sqrt(fpcCochranC); bC = xbar10.ori*(1-aC)
  xbar10adj = aC*xbar10 + bC
  
  xbar.anchor.adj = xbar11adj*phat11 + xbar10adj*phat10 + xbar01adj*phat0dot 
  
  ## compute estimated xbar among cases 
  phat1dot.CV = (n11.CV + n10.CV)/Nhat.CV
  phat0dot.CV = 1 - phat1dot.CV
  xbar.anchor.CV = xbar1dot.CV*phat1dot.CV + xbar01.CV*phat0dot.CV
  
  ## among non-cases
  phat1dot.noCV = (n11.noCV + n10.noCV)/Nhat.noCV
  phat0dot.noCV = 1 - phat1dot.noCV
  xbar.anchor.noCV = xbar1dot.noCV*phat1dot.noCV + xbar01.noCV*phat0dot.noCV
  
  ## difference
  xbar1dot.dif = xbar1dot.CV - xbar1dot.noCV
  xbar2dot.dif = xbar2dot.CV - xbar2dot.noCV
  xbar.anchor.dif = xbar.anchor.CV - xbar.anchor.noCV
  
  re.est = data.frame(xbar2dot = xbar2dot, 
                      xbar2dot.adj = xbar2dot.adj,
                      xbar1dot = xbar1dot, 
                      xbar.anchor = xbar.anchor,
                      xbar.anchor.adj = xbar.anchor.adj,
                      xbar2dot.CV = xbar2dot.CV, xbar1dot.CV = xbar1dot.CV,
                      xbar.anchor.CV = xbar.anchor.CV,
                      xbar2dot.noCV = xbar2dot.noCV, 
                      xbar1dot.noCV = xbar1dot.noCV,
                      xbar.anchor.noCV = xbar.anchor.noCV,
                      xbar2dot.dif = xbar2dot.dif, 
                      xbar1dot.dif = xbar1dot.dif,
                      xbar.anchor.dif = xbar.anchor.dif)
  
  return(re.est)
  
}


## A function to generate non-parametric bootstrap samples
bootsamp <- function(dat.obs){
  ## INPUT: observed data frame
  nobs = nrow(dat.obs)
  sel.index = sample(x = 1:nobs, size = nobs, replace = T)
  dat.boot = dat.obs[sel.index, ]
  ## OUTPUT: a data frame contains the bootstrap samples
  return(dat.boot)
}


## A function to estimate continuous mean
AnchorStream_Continuous <- function(dat.obs, Ntot, seed, nboot){
  ## INPUT:
  ## dat.obs: observed individual-level data
  ## Ntot: total # of populations
  ## seed: random seed
  ## nboot: # of bootstrap replicates
  
  ### get observed # of cases and non-cases
  nc.CV.i = sum(dat.obs$case)
  nc.noCV.i = sum(dat.obs$case == 0)
  re.nxbar.i = get_nxbar_obs(dat.obs = dat.obs, Ntot = Ntot)
  n.vec.i = re.nxbar.i$nvec
  xbar.vec.i = re.nxbar.i$xbarvec
  ### get point estiamtes ###
  re.i = get_xbar_obs(dat.obs = dat.obs, Ntot = Ntot,
                      nc.CV.sim = nc.CV.i, nc.noCV.sim = nc.noCV.i,
                      n.vec.ori = n.vec.i,
                      xbar.vec.ori = xbar.vec.i)
  ### get interval estimations by following procedures presented in Appendix B
  set.seed(seed)
  # (1). generate bootstrap samples based on observed data 
  dat.boot.list <- NULL
  for(iboot in 1:nboot){
    dat.boot.list[[iboot]] <- bootsamp(dat.obs =  dat.obs)
  }
  # (2). get estimates for each of bootstrap replicates 
  re.boot.i <- lapply(dat.boot.list, get_xbar_obs, Ntot = Ntot,
                      nc.CV.sim = nc.CV.i,
                      nc.noCV.sim = nc.noCV.i,
                      n.vec.ori = n.vec.i,
                      xbar.vec.ori = xbar.vec.i)
  re.boot.i <- do.call(rbind.data.frame, re.boot.i)
  # (3). construct percentile intervals 
  sd.boot.i <- apply(re.boot.i, 2, sd, na.rm = T)
  lci.boot.i <- apply(re.boot.i, 2, quantile, 0.025, na.rm = T)
  uci.boot.i <- apply(re.boot.i, 2, quantile, 0.975, na.rm = T)
  
  ### get cleaned results ###
  sel.names.overall <- c("xbar1dot", "xbar2dot.adj", "xbar.anchor.adj")
  sel.names.cases <- c("xbar1dot.CV", "xbar2dot.CV", "xbar.anchor.CV")
  sel.names.noncases <- c("xbar1dot.noCV", "xbar2dot.noCV", "xbar.anchor.noCV")
  sel.names.diff <- c("xbar1dot.dif", "xbar2dot.dif", "xbar.anchor.dif")
  
  est.overall.vec = c(re.i$xbar1dot, re.i$xbar2dot, re.i$xbar.anchor)
  est.cases.vec = c(re.i$xbar1dot.CV, re.i$xbar2dot.CV, re.i$xbar.anchor.CV)
  est.noncases.vec = c(re.i$xbar1dot.noCV, re.i$xbar2dot.noCV, re.i$xbar.anchor.noCV)
  est.diff.vec = c(re.i$xbar1dot.dif, re.i$xbar2dot.dif, re.i$xbar.anchor.dif)
  
  se.overall.vec = sd.boot.i[sel.names.overall]
  se.cases.vec = sd.boot.i[sel.names.cases]
  se.noncases.vec = sd.boot.i[sel.names.noncases]
  se.diff.vec = sd.boot.i[sel.names.diff]
  
  lci.overall.vec = lci.boot.i[sel.names.overall]
  lci.cases.vec = lci.boot.i[sel.names.cases]
  lci.noncases.vec = lci.boot.i[sel.names.noncases]
  lci.diff.vec = lci.boot.i[sel.names.diff]
  
  uci.overall.vec = uci.boot.i[sel.names.overall]
  uci.cases.vec = uci.boot.i[sel.names.cases]
  uci.noncases.vec = uci.boot.i[sel.names.noncases]
  uci.diff.vec = uci.boot.i[sel.names.diff]
  
  tab.overall <- data.frame(est = est.overall.vec,
                            se = se.overall.vec,
                            lci = lci.overall.vec, uci = uci.overall.vec)
  rownames(tab.overall) <- c("xbar1dot", "xbar2dot", "muhatx")
  tab.overall[1, 2:4] <- NA
  tab.cases <- data.frame(est = est.cases.vec,
                          se = se.cases.vec,
                          lci = lci.cases.vec, uci = uci.cases.vec)
  rownames(tab.cases) <- c("xbar1dot.cases", "xbar2dot.cases", "muhatx.cases")
  
  tab.noncases <- data.frame(est = est.noncases.vec,
                             se = se.noncases.vec,
                             lci = lci.noncases.vec, uci = uci.noncases.vec)
  rownames(tab.noncases) <- c("xbar1dot.noncases", "xbar2dot.noncases", "muhatx.noncases")
  
  tab.diff <- data.frame(est = est.diff.vec,
                         se = se.diff.vec,
                         lci = lci.diff.vec, uci = uci.diff.vec)
  rownames(tab.diff) <- c("xbar1dot.diff", "xbar2dot.diff", "muhatx.diff")
  
  return(list(overall = tab.overall,
              cases = tab.cases,
              noncases = tab.noncases,
              difference = tab.diff))
}



# A function to generate individual-level data
# for each record, it contains the capture history, 
# binary indicator of disease status, binary indicator for symptom status
# a continuous variable of interest 
gen_dat <- function(p.CV, p.CV.sympt, p.noCV.sympt, 
                    mu11, mu10, mu01, mu00, 
                    sd11, sd10, sd01, sd00,
                    Ntot, 
                    p2, p1.sympt, p1.nosympt,
                    seed){
  
  # INPUT:
  # p.CV: prevalence
  # p.CV.sympt: proportion of symptom | cases
  # p.noCV.sympt: proportion of symptom | non-cases
  # mu11/sd11: mean/sd of x among cases with symptom
  # mu10/sd10: mean/sd of x among cases without symptom
  # mu01/sd01: mean/sd of x among non-cases with symptom 
  # mu00/sd00: mean/sd of x among non-cases without symptom 
  # Ntot: total # populations
  # p2: prob of being captured by S2 (anchor stream)
  # p1.sympt: prob of being captured by S1 | symptom
  # p1.nosympt: prob of being captured by S1 | no symptom
  # seed: random seed
  
  set.seed(seed) 
  
  prop11 = p.CV*p.CV.sympt
  prop10 = p.CV*(1-p.CV.sympt)
  prop01 = (1-p.CV)*p.noCV.sympt
  prop00 = (1-p.CV)*(1-p.noCV.sympt)
  mu.true.vec = c(mu11, mu10, mu01, mu00)
  prop.true.vec = c(prop11, prop10, prop01, prop00)
  
  mean.true = sum(mu.true.vec*prop.true.vec)
  mean.true.cases = p.CV.sympt*mu11 + (1-p.CV.sympt)*mu10
  
  n.CV = Ntot*p.CV
  
  # Data Generation #
  covid = rep(0, Ntot)
  index.CV = sample(1:Ntot, n.CV)
  covid[index.CV] = 1
  
  sympt = rep(NA, Ntot)
  sympt[covid == 1] <- rbinom(sum(covid == 1), 1, p.CV.sympt)
  sympt[covid == 0] <- rbinom(sum(covid == 0), 1, p.noCV.sympt)
  
  # outcome (continuous variables) 
  x <- rep(NA, Ntot)
  index1 = (covid == 1 & sympt == 1)
  index2 = (covid == 1 & sympt == 0)
  index3 = (covid == 0 & sympt == 1)
  index4 = (covid == 0 & sympt == 0)
  x[index1] <- mu11 + rnorm(sum(index1), sd = sd11)
  x[index2] <- mu10 + rnorm(sum(index2), sd = sd10)
  x[index3] <- mu01 + rnorm(sum(index3), sd = sd01)
  x[index4] <- mu00 + rnorm(sum(index4), sd = sd00)
  
  dat.pop = data.frame(ID = 1:Ntot, covid = covid,
                       sympt = sympt, x = x)
  
  # generate capture histories 
  ## Stream 2 is the anchor stream (random sample)
  ## Stream 1 is the non-random stream 
  y1 <- rep(NA, Ntot)
  y1[dat.pop$sympt == 1] <- rbinom(sum(dat.pop$sympt == 1), 1, p1.sympt)
  y1[dat.pop$sympt == 0] <- rbinom(sum(dat.pop$sympt == 0), 1, p1.nosympt)
  
  y2 <- rep(0, Ntot)
  index.y2 = sample(1:Ntot, Ntot*p2)
  y2[index.y2] <- 1
  
  dat.pop = data.frame(dat.pop, y1 = y1, y2 = y2)
  colnames(dat.pop)[colnames(dat.pop) == "covid"] <- "case"
  return(list(dat.pop = dat.pop, 
              mean.true = mean.true, 
              mean.true.cases = mean.true.cases))
}
