library(Directional)
library(rotations)
library(tidyverse)
library(here)
source(here("R","helper_funs.R"))

dir.create(here("R","figs"))

# sample a random point on a sphere
init_x <- function(n=1,r=1){
  phi <- runif(n,0,pi)
  theta <- runif(n,0,2*pi)
  r*c(sin(phi)*cos(theta),sin(phi)*sin(theta),cos(phi))
}

# base kernel on a sphere
# randomly sample a rotation matrix from SO(3)
base_kernel <- function(position,r=1){
  as.vector(sample_rot(n=1) %*% as.matrix(position,nrow=3))
}

# orbit kernel on a sphere where z defines orbit
orbit_kernel <- function(position,r=1){
  tmp_r <- sphere_radius_z(position[3],r)
  current_theta <- acos(position[1]/tmp_r) 
  # sample a rotation angle
  sample_rotation <- runif(1,0,2*pi)
  updated_theta <- current_theta + sample_rotation
  return(c(tmp_r*cos(updated_theta),tmp_r*sin(updated_theta),position[3]))
}

# naive MCMC algorithm based on sampling uniformly over a sphere
algo_naive <- function(N=10,base_kernel,r=1){
  samples <- matrix(double(N*3),nrow=N)
  # initialization
  samples[1,] <- init_x(n=1,r=r)
  for (i in 2:N){
    # sample a new position on the original space
    samples[i,] <- base_kernel(samples[i-1,],r=r)
  }
  return(samples)
}

# orbit MCMC based on alternation
# between the base kernel and orbit kernel
algo_orbit<- function(N=10,base_kernel,orbit_kernel,r=1){
  samples <- matrix(double(N*3),nrow=N)
  tmp <- double(3)
  # initialization
  samples[1,] <- init_x(n=1,r=r)
  for (i in 2:N){
    # sample a new position on the original space
    tmp <- base_kernel(samples[i-1,],r=r)
    samples[i,] <- orbit_kernel(tmp,r=r)
  }
  return(samples)
}

# orbit MCMC based on mixture of
# the base kernel and orbit kernel
algo_orbit_mixture <- function(N=10,base_kernel,orbit_kernel,
                               r=1,mixture_p=0.5){
  
  samples <- matrix(double(N*3),nrow=N)
  # initialization
  samples[1,] <- init_x(n=1,r=r)
  
  for (i in 2:N){
    if(rbernoulli(1,p=mixture_p)){
      samples[i,] <-  orbit_kernel(samples[i-1,],r=r)
    } else{
      samples[i,] <-  base_kernel(samples[i-1,],r=r)
    }
  }
  return(samples)
}



# sanity check: plot sampled points on a sphere
NN <- 500 # num of samples
naive <- algo_naive(N=NN,base_kernel,r=1)
orbit_alt <- algo_orbit(N=NN,base_kernel,orbit_kernel,r=1)
orbit_mix <- algo_orbit_mixture(N=NN,base_kernel,orbit_kernel,r=1)

# sphereplot(naive,col='red')
# sphereplot(orbit_alt,col='blue')
# sphereplot(orbit_mix,col='brown')
sphereplot3(naive,orbit_alt,orbit_mix)

# eval
NN <- 10000 # num of samples
naive <- algo_naive(N=NN,base_kernel,r=1)
orbit_alt <- algo_orbit(N=NN,base_kernel,orbit_kernel,r=1)
orbit_mix <- algo_orbit_mixture(N=NN,base_kernel,orbit_kernel,r=1)

# compute temperature for all
df <- data.frame(z=naive[,3],temperature=temp(naive[,3]))  %>% 
  mutate(type='naive') %>% 
  rbind(.,data.frame(z=orbit_alt[,3],temperature=temp(orbit_alt[,3]))  %>% 
          mutate(type='orbit_alt')) %>% 
  rbind(.,data.frame(z=orbit_mix[,3],temperature=temp(orbit_mix[,3]))  %>% 
          mutate(type='orbit_mix')) %>% 
  rbind(.,data.frame(data.frame(z=0,
                               temperature=310-80*rbeta(NN,0.5,1),
                               type='truth')))

# temp
num.bins <- 100
truth <- tranformed_density()
ggplot(data=df %>% 
         mutate(type = factor(type,levels=c(lvls,"truth"))) %>% 
         filter(type!="truth"),
       aes(temperature,color=type)) +
  geom_histogram(bins = num.bins,position='identity',
                 alpha=0.5,fill="white") +
  geom_line(data=data.frame(x=truth$x,y=truth$y) %>% 
              filter(x>230 & x < 310) %>% 
              mutate(y=y*10000),aes(x=x,y=y),col='purple') +
  # geom_hline(yintercept=1/num.bins,colour='purple',size=1.5,alpha=0.5)+
  theme_classic(base_size=text_size) +
  scale_color_manual(values=c("naive"="red",
                              "orbit_alt"="blue",
                              "orbit_mix"="brown",
                              "truth"="purple"))+
  theme(legend.position = 'top',
        legend.title=element_blank()) +
  ylab("") +
  scale_x_continuous(breaks=seq(230,310,by=10),
                     label=seq(230,310,by=10))  +
  ylim(0,1000) -> gg.sanity
ggsave(here("R","figs","hist_sanity_check.png"),gg.sanity,"png",dpi = 300)

# formally compare

num.bins <- 100
NN <- 11

result <- data.frame(N=double(NN),
                     naive=double(NN),
                     orbit_alt=double(NN),
                     orbit_mix=double(NN))
result_param <- c()
result_kl <- data.frame(N=double(NN),
                        naive=double(NN),
                        orbit_alt=double(NN),
                        orbit_mix=double(NN))
result_chain <- c()

for(i in 1:NN){
  N <- 2^(i+5)
  naive <- algo_naive(N=N,base_kernel,r=1)
  orbit_alt <- algo_orbit(N=N,base_kernel,orbit_kernel,r=1)
  orbit_mix <- algo_orbit_mixture(N=N,base_kernel,orbit_kernel,r=1)
  df <- data.frame(z=naive[,3],temperature=temp(naive[,3]))  %>% 
    mutate(type='naive') %>% 
    rbind(.,data.frame(z=orbit_alt[,3],temperature=temp(orbit_alt[,3]))  %>% 
            mutate(type='orbit_alt')) %>% 
    rbind(.,data.frame(z=orbit_mix[,3],temperature=temp(orbit_mix[,3]))  %>% 
            mutate(type='orbit_mix')) 
  
  naive.density <- hist(df %>% 
                          filter(type=='naive') %>% 
                          select(temperature) %>% 
                          unlist(),
                        breaks = seq(230,310,length.out = num.bins),plot=F)$density
  orbit_alt.density <- hist(df %>% 
                              filter(type=='orbit_alt') %>% 
                              select(temperature) %>% 
                              unlist(),
                            breaks = seq(230,310,length.out = num.bins),plot=F)$density
  orbit_mix.density <- hist(df %>% 
                              filter(type=='orbit_mix') %>% 
                              select(temperature) %>% 
                              unlist(),
                            breaks = seq(230,310,length.out = num.bins),plot=F)$density
  result[i,] <- c(N,
                  mse(naive.density,true.density),
                  mse(orbit_alt.density,true.density),
                  mse(orbit_mix.density,true.density))
  result_kl[i,] <- c(N,
                  LaplacesDemon::KLD(true.density,naive.density)$sum.KLD.py.px,
                  LaplacesDemon::KLD(true.density,orbit_alt.density)$sum.KLD.py.px,
                  LaplacesDemon::KLD(true.density,orbit_mix.density)$sum.KLD.py.px)
  
  result_param[[i]] <- df %>% 
    group_by(type) %>% 
    group_split() %>% 
    lapply(.,function(tmp){
      tmp2 <- EnvStats::ebeta(tmp$z^2)$parameters
      return(data.frame(alpha=tmp2[1],beta=tmp2[2],type=tmp$type[1]))
    }) %>% 
    do.call(rbind,.) %>% 
    mutate(N=N)
    
  result_chain[[i]] <- df %>%
    mutate(id=rep(1:N,3)) %>% 
    pivot_wider(id_cols=id,
                names_from=type,
                values_from = temperature )
  print(paste0(i," done"))
}
write_rds(list(result=result,
               result_kl=result_kl,
               result_chain=result_chain,
               result_param=result_param),
          here("R","results.rds"))

results <- read_rds(here("R","results.rds"))
result <- results$result

# rmse between true density and estimated density
ggplot(data=pivot_longer(result,cols = 2:4,names_to = 'type',
                         values_to = 'MSE') %>% 
         mutate(RMSE = sqrt(MSE),
                logN = log(N,2),
                type = factor(type,levels=lvls)),
       aes(x=logN,y=RMSE,color=type)) +
  geom_line() + 
  geom_point()+
  scale_color_manual(values = cols)+
  theme_classic(base_size=text_size)+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  ylab("RMSE") +
  xlab("log_2(N)") -> gg.RMSE
ggsave(here("R","figs","RMSE.png"),gg.RMSE,"png",dpi = 300)

ggplot(data=pivot_longer(result_kl,cols = 2:4,names_to = 'type',
                         values_to = 'KL') %>% 
         mutate(logN = log(N,2),
                type = factor(type,levels=lvls)),
       aes(x=logN,y=KL,color=type)) +
  geom_line() + 
  geom_point()+
  scale_color_manual(values = cols)+
  theme_classic(base_size=text_size)+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  ylab("KL") +
  xlab("log_2(N)") +
  scale_x_continuous(breaks=c(11:16),limits = c(11,16)) -> gg.KL
ggsave(here("R","figs","KL.png"),gg.KL,"png",dpi = 300)

# ESS & MCE for the mean temperature
library(mcmcse)
result_chain <- results$result_chain
NN <- 11
chain <- result_chain[[NN]][,-1] %>% as.matrix()
c(ess(chain)[1],ess(chain)[2]/2,ess(chain)[3])

mcerror_bm_mean <- mcse.multi(x = chain, method =  "bm",
                             size = "sqroot", g = NULL)
mcerror_bm_mean$est
sqrt(diag(mcerror_bm_mean$cov))

mcerror_bm_second <- mcse.multi(x = chain, method =  "bm",
                         size = "sqroot", g = function(xx){xx^2})
mcerror_bm_second$est
sqrt(diag(mcerror_bm_second$cov))

# mean temp
burn_in <- 0
df_chain <- chain[-c(1:burn_in),] %>% 
  as.data.frame() %>% 
  mutate(n= row_number(),
         naive=cumsum(naive)/n,
         orbit_alt = cumsum(orbit_alt)/n,
         orbit_mix = cumsum(orbit_mix)/n) %>% 
  pivot_longer(cols=1:3,names_to="type",values_to='MC_Estimate')

ggplot(data=df_chain %>% 
         filter(n>1000) %>% 
         mutate(type = factor(type,levels=lvls),
                n = log(n,2)),aes(x=n,y=MC_Estimate,color=type))+
  geom_line() + 
  ylim(281,285) +
  scale_x_continuous(breaks=c(10:16)) +
  geom_hline(yintercept=310-0.5/1.5*80,
              col='purple',size=2,
             alpha=0.5) +
  ylab("MC estimate of mean temperature") +
  xlab("log_2(N)") +
  theme_classic(base_size=text_size)+
  scale_color_manual(values=cols)+
  theme(legend.position = 'top',
        legend.title=element_blank()) -> gg.mc.sample
ggsave(here("R","figs","MCE.png"),gg.mc.sample,"png",dpi = 300)

library(gridExtra)
library(grid)
tmp1 <- gg.RMSE +
  theme(legend.position='none',
        axis.title.x = element_blank())

tmp2 <- gg.mc.sample +
  theme(legend.position='none',
        axis.title.x = element_blank()) +
  ylab("MC estimate of T")

leg <- lemon::g_legend(gg.RMSE)

combined.gg <- grid.arrange(leg, arrangeGrob(tmp1,tmp2, nrow=1), ncol=1, nrow=2,
             heights=c(1,10),
             bottom=textGrob("log_2(N)",
                             gp = gpar(fontsize = 20)))
ggsave(here("R","figs","combined.png"),combined.gg,"png",dpi=300)
