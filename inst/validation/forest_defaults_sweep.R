# Task 4 validation: current(500,gini) vs fast(100,extratrees) forest defaults.
# Production perm.N=1000, reps=10, n in {100,500}, several DGPs (core + realistic).
# Run: Rscript inst/validation/forest_defaults_sweep.R [OUTDIR]
suppressPackageStartupMessages({ library(MASS); library(parallel); library(devtools) })
OUTDIR <- ifelse(length(commandArgs(TRUE))>=1, commandArgs(TRUE)[1], "/tmp/forest_val")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
PLOG <- file.path(OUTDIR,"progress.log"); cat("", file=PLOG)
devtools::load_all(".", quiet=TRUE)            # Phase-1 fastcpt (classifier.args passthrough)
options(fastcpt.num.threads=1L)                # outer-parallel; ranger single-thread
expit <- function(z) 1/(1+exp(-z))

generate_data <- function(scenario,n,p,rho,delta,seed){
  set.seed(seed); Sigma<-matrix(rho,p,p); diag(Sigma)<-1
  X<-pnorm(MASS::mvrnorm(n,rep(0,p),Sigma))
  sig_l<-sd(X[,1]-0.5); sig_m<-sd((X[,1]+X[,2]+X[,3])/3-0.5)
  int_raw<-(X[,1]-0.5)*(X[,2]-0.5); sig_i<-sd(int_raw-mean(int_raw))
  nl_raw<-(X[,1]-0.5)^2; sig_n<-sd(nl_raw-mean(nl_raw))
  prob<-switch(scenario, "null"=rep(0.5,n), "null_unequal"=rep(1/3,n),
    "linear"=expit(delta*(X[,1]-0.5)/sig_l),
    "multi"=expit(delta*((X[,1]+X[,2]+X[,3])/3-0.5)/sig_m),
    "interaction"=expit(delta*(int_raw-mean(int_raw))/sig_i),
    "nonlinear"=expit(delta*(nl_raw-mean(nl_raw))/sig_n))
  W<-rbinom(n,1,prob); colnames(X)<-paste0("X",seq_len(p)); list(X=as.data.frame(X),W=W)
}
generate_realistic <- function(n,delta,seed,scenario){
  set.seed(seed); Sc<-matrix(0.2,4,4); diag(Sc)<-1; cr<-MASS::mvrnorm(n,rep(0,4),Sc)
  age<-18+50*pnorm(cr[,1]); income<-exp(2+1.5*cr[,2]); ideology<-pnorm(cr[,3]); knowledge<-pnorm(cr[,4])
  female<-rbinom(n,1,.52); party<-rbinom(n,1,.45); college<-rbinom(n,1,.35)
  educ<-ordered(sample(1:5,n,TRUE,c(.10,.25,.30,.20,.15))); interest<-ordered(sample(1:4,n,TRUE,c(.20,.30,.30,.20)))
  region<-factor(sample(c("Northeast","South","Midwest","West"),n,TRUE,c(.18,.38,.22,.22)))
  race<-factor(sample(c("White","Black","Hispanic","Asian","Other"),n,TRUE,c(.60,.13,.18,.06,.03)))
  X<-data.frame(age,income,ideology,knowledge,female,party,college,educ,interest,region,race)
  age_c<-(age-mean(age))/sd(age); ideo_c<-(ideology-mean(ideology))/sd(ideology)
  signal<-switch(scenario, "interaction"={s<-ideo_c*(party-mean(party));s/sd(s)},
    "nonlinear"={s<-age_c^2-mean(age_c^2);s/sd(s)},
    "factor"={s<-ifelse(region=="South",1,0)-mean(region=="South");s/sd(s)}, "null"=rep(0,n))
  W<-rbinom(n,1,expit(delta*signal)); list(X=X,W=W)
}

REPS<-10L; PERM<-1000L; ALPHA<-0.05; P<-10L; RHO<-0; DALT<-0.75
CFG<-list(current=list(num.trees=500L),
          fast=list(num.trees=100L, splitrule="extratrees", num.random.splits=1L))
# cells: experiment, scenario, n, delta
core_sc<-c("null","null_unequal","linear","multi","interaction","nonlinear")
cells<-rbind(
  expand.grid(exp="core", scenario=core_sc, n=c(100L,500L), stringsAsFactors=FALSE),
  expand.grid(exp="realistic", scenario=c("interaction","nonlinear","null"), n=500L, stringsAsFactors=FALSE))
cells$delta<-ifelse(cells$scenario %in% c("null","null_unequal"), 0, DALT)
# flatten to (cell, rep, config) tasks
tasks<-do.call(rbind, lapply(seq_len(nrow(cells)), function(ci)
  do.call(rbind, lapply(1:REPS, function(r) data.frame(cells[ci,], rep=r,
    config=names(CFG), stringsAsFactors=FALSE, row.names=NULL)))))
tot<-nrow(tasks); cat(sprintf("START %d tasks (%d cells x %d reps x 2 cfg) perm.N=%d\n",tot,nrow(cells),REPS,PERM),file=PLOG,append=TRUE)
SEEDB<-list(core=1e6, realistic=2e6)
t0all<-proc.time()[3]
run<-function(i){
  tk<-tasks[i,]; dseed<-SEEDB[[tk$exp]] + match(tk$scenario,c(core_sc,"factor"))*1e4 + tk$n*10 + tk$rep
  d<-if(tk$exp=="core") generate_data(tk$scenario,tk$n,P,RHO,tk$delta,dseed) else generate_realistic(tk$n,tk$delta,dseed,tk$scenario)
  t0<-proc.time()[3]
  pv<-tryCatch(unname(fastcpt(Z=d$X,T=d$W,class.methods="forest",classifier.args=CFG[[tk$config]],
        perm.N=PERM,parallel=FALSE,progress=FALSE,R.seed=dseed)$pvals["forest"]),error=function(e)NA_real_)
  el<-proc.time()[3]-t0
  cat(sprintf("%3d/%d | %-9s | %-12s n%-4d r%02d | %-7s p=%.3f | %5.1fs | elapsed %.0fs\n",
      i,tot,tk$exp,tk$scenario,tk$n,tk$rep,tk$config,pv,el,proc.time()[3]-t0all),file=PLOG,append=TRUE)
  data.frame(tk, pval=pv, sec=el)
}
MC<-max(1L,parallel::detectCores()-1L)
res<-do.call(rbind, parallel::mclapply(seq_len(tot), run, mc.cores=MC, mc.preschedule=FALSE))
res$reject<-as.integer(res$pval<ALPHA)
saveRDS(res, file.path(OUTDIR,"res.rds"))
ag<-aggregate(cbind(rate=reject, mean_p=pval, sec=sec)~exp+scenario+n+config, res, mean, na.action=na.pass)
ag<-ag[order(ag$exp,ag$scenario,ag$n,ag$config),]
sink(file.path(OUTDIR,"summary.txt"))
cat(sprintf("VALIDATION current(500,gini) vs fast(100,extratrees)  perm.N=%d reps=%d p=%d delta_alt=%.2f a=%.2f\n",PERM,REPS,P,DALT,ALPHA))
cat("rate = rejection rate (size if null/null_unequal, else power)\n")
key<-unique(ag[,c("exp","scenario","n")])
worst_power_drop<-0; max_size<-0; na_fast<-0
for(k in seq_len(nrow(key))){
  e<-key$exp[k];s<-key$scenario[k];nn<-key$n[k]
  cur<-ag[ag$exp==e&ag$scenario==s&ag$n==nn&ag$config=="current",]
  fst<-ag[ag$exp==e&ag$scenario==s&ag$n==nn&ag$config=="fast",]
  is_null<- s %in% c("null","null_unequal")
  lab<-if(is_null)"size" else "power"
  cat(sprintf("%-9s %-12s n=%-4d  [%s]  current=%.2f  fast=%.2f  (mean_p %.3f/%.3f)  sec %.1f/%.1f\n",
      e,s,nn,lab,cur$rate,fst$rate,cur$mean_p,fst$mean_p,cur$sec,fst$sec))
  if(is_null){ max_size<-max(max_size,cur$rate,fst$rate) } else { worst_power_drop<-max(worst_power_drop, cur$rate-fst$rate) }
}
na_fast<-sum(is.na(res$pval[res$config=="fast"]))
sc_t<-aggregate(sec~config,res,mean)
cat(sprintf("\nspeed: current %.1fs/call, fast %.1fs/call -> %.2fx\n",
    sc_t$sec[sc_t$config=="current"],sc_t$sec[sc_t$config=="fast"],
    sc_t$sec[sc_t$config=="current"]/sc_t$sec[sc_t$config=="fast"]))
cat(sprintf("\nGATE CHECK (reps=10, coarse): max size (either cfg) = %.2f ; worst power drop (current-fast) = %.2f ; fast NA pvals = %d\n",
    max_size, worst_power_drop, na_fast))
verdict <- if(na_fast==0 && max_size<=0.20 && worst_power_drop<=0.15) "LEAN GO" else "LEAN NO-GO / INSPECT"
cat(sprintf("heuristic: %s  (size<=0.20, power_drop<=0.15, no NA)\n", verdict))
cat("\nrun complete\n"); sink()
cat("WROTE summary\n")
