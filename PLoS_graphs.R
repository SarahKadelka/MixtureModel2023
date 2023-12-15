### Define path to folder where graphs should be saved
folder_graphs = "~/Desktop/Kadelka23_PLoS_Comp_Biol_Figures/" #"~/Documents/AJE_Mixture_model_method/PLoS_submission_figures/"
dir.create(folder_graphs, recursive = TRUE)

## Path to where the provided data sets are saved
folder_data = "~/Desktop/Kadelka23_PLoS_Comp_Biol_Figures/"  #"~/Documents/AJE_Mixture_model_method/PLoS_submission_figures/"


library(sn)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(lubridate)
library(pracma)
library(xtable)
library(SimDesign)


###
# Antibody loads of negative controls
# - skewed normally distributed with parameters
###
xi.pos = 0.9
omega.pos = 0.54
alpha.pos = -8
###
# Antibody loads of negative controls
# - skewed normally distributed with parameters
###
xi.neg = -1.8
omega.neg = 0.75
alpha.neg = 2
###
# Decay rates of antibody loads of cases
# - gamma distributed with parameters
###
shape.decay = 1.93
rate.decay = 600.52

### Define test cases
test.case.1.a = c(1,3,30,10,0.1,0.1,0.5,2,1)/100
test.case.2.a = c(1,0.5,2,0.1,3,30,10,0.1,1)/100
test.case.3.a = c(0.5,1,5.5,5.5,0.1,0.2,8,14,1)/100
test.case.4.a = c(0.5,1,14,8,0.1,0.2,5.5,5.5,1)/100
test.case.5.a = c(0.5,1,22,15,4,0.1,0.5,0.1,0.1,0.4,2,5,12,8,0.5)/100
test.case.6.a = c(0.2,0.1,0.7,1.3,0.2,0.1,0.5,0.8,0.2)/100
test.case.7.a = c(0.5,0.2,0.7,0.02,0.7,1,0.3,0.5,0.4,0.9,1,0.5,0.1,0.8,0.6,0.6,0.07)/100




########################################
## Figure1
########################################
###############
## Figure1A
###############
backgrounds.plot =  rsn(100000, xi = xi.neg, omega = omega.neg, alpha = alpha.neg)
dens.background = density(backgrounds.plot)
i.label.background = which.max(dens.background$y)
peaks.plot =  rsn(100000, xi = xi.pos, omega = omega.pos, alpha = alpha.pos)
dens.peaks = density(peaks.plot)
i.label.peaks = which.max(dens.peaks$y)
fig1A = ggplot() +
  stat_function(data = data.frame(x = c(-3, 2)), aes(x),
                fun = dsn, n = 1001,
                args = list( xi = xi.pos, omega = omega.pos, alpha = alpha.pos),
                size = 2, col = rgb(0,0,0,1)) +
  stat_function(data = data.frame(x = c(-3, 2)), aes(x),
                fun = dsn, n = 1001,
                args = list( xi = xi.neg, omega = omega.neg, alpha = alpha.neg),
                col = rgb(0,0.4,0,1),
                size = 2) +
  geom_histogram(aes(x = peaks.plot, y = ..density..), col = rgb(0,0,0,0.75), fill = rgb(0,0,0,0.2)) +
  geom_histogram(aes(x = backgrounds.plot, y = ..density..), col = rgb(0,0.4,0,0.75), fill = rgb(0,0.4,0,0.2)) +
  xlab("log(AB level)") +
  ylab("Density") +
  scale_y_continuous(breaks = NULL, expand=expansion(mult = c(0,0.05),add = 0)) +
  scale_x_continuous(breaks = NULL, expand=c(0,0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 26))  +
  annotate("text", x = dens.background$x[i.label.background], y = 1.25*dens.background$y[i.label.background] , label = "background\nlevels", size = 8, col = rgb(0,0.4,0)) +
  annotate("text", x = dens.peaks$x[i.label.peaks], y = 1.05*dens.peaks$y[i.label.peaks] , label = "peaks", size = 8, col = rgb(0,0,0))

name.pic = paste0(folder_graphs,"Figure1A.pdf")
pdf(file=name.pic, height=8/2, width = 10/2)
print(fig1A)
dev.off()


###############
## Figure1B
###############
decays.plot = rgamma(100000, shape = shape.decay, rate = rate.decay)
dens.decay = density(decays.plot)
i.label.decay = which.max(dens.decay$y)
fig1B = ggplot() +
  stat_function(data = data.frame(x = c(0, 0.02)), aes(x),
                fun = dgamma, n = 1001,
                args = list( shape = shape.decay, rate = rate.decay),
                size = 2,
                col = rgb(0,0,0,1)) +
  geom_histogram(aes(x = decays.plot, y = ..density..), col = rgb(0,0,0,0.75), fill = rgb(0,0,0,0.2)) +
  xlab("- d/dt(log(AB level(t)))") +
  ylab("Density") +
  scale_y_continuous(breaks = NULL, expand=expansion(mult = c(0,0.05),add = 0)) +
  scale_x_continuous(breaks = NULL, expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 26))+
  annotate("text", x = 3.25*dens.decay$x[i.label.decay], y = 1.05*dens.decay$y[i.label.decay] , label = "decay rates", size = 8, col = rgb(0,0,0))

name.pic = paste0(folder_graphs,"Figure1B.pdf")
pdf(file=name.pic, height=8/2, width = 10/2)
print(fig1B)
dev.off()


###############
## Figure1C
###############
N_total = 100
Apeak.vec = rsn(n = N_total, xi = xi.pos, omega = omega.pos, alpha = alpha.pos)
Aneg.vec = rsn(n = N_total, xi = xi.neg, omega = omega.neg, alpha = alpha.neg)
r.vec = rgamma(n = N_total,  shape = shape.decay, rate = rate.decay)
T1 = 0
tpeaks = rgamma(N_total , shape = 30/2, rate = 0.5)

times = 0:300
A = data.frame(matrix(NA, nrow = length(times), ncol = N_total))%>%mutate(time = times)
for (i in 1:length(times)){
  A[i,1:N_total] =  kadelka23cutoffFree::antibody_dynamics(t = rep(times[i], N_total),
                                                           AB.peak = Apeak.vec,
                                                           AB.background = Aneg.vec,
                                                           decay.rate = r.vec,
                                                           t.growth = T1,
                                                           t.peak = tpeaks,
                                                           growth.dynamics = "exp")
}

fig1C = ggplot() +
  geom_line(data = reshape2::melt(A, id.vars = "time"),
            aes(x = time, y = value, group = variable)) +
  geom_rect(aes(xmin=0,
                xmax=29,
                ymin=min((reshape2::melt(A, id.vars = "time"))$value),
                ymax=max((reshape2::melt(A, id.vars = "time"))$value)),
            fill=rgb(1,0,0,0.2), color=rgb(1,0,0,0.2)) +
  geom_rect(aes(xmin=30,
                xmax=59,
                ymin=min((reshape2::melt(A, id.vars = "time"))$value),
                ymax=max((reshape2::melt(A, id.vars = "time"))$value)),
            fill=rgb(0,0,1,0.2), color=rgb(0,0,1,0.2)) +
  geom_rect(aes(xmin=60,
                xmax=89,
                ymin=min((reshape2::melt(A, id.vars = "time"))$value),
                ymax=max((reshape2::melt(A, id.vars = "time"))$value)),
            fill=rgb(0,1,1,0.2), color=rgb(0,1,1,0.2)) +
  ylab("log(AB level)") +
  xlab("time (in months)") +
  scale_y_continuous(breaks = NULL, expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,300,30), labels = seq(0,10,1), expand = expansion(mult = 0,add = c(0,10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 26))

name.pic = paste0(folder_graphs,"Figure1C.pdf")
pdf(file=name.pic, height=8/2, width = 10/2)
print(fig1C)
dev.off()


###############
## Figure1D
###############
AB.growth.option = "exp0"
time.between.surveys = 30
f.pos.AB.vals.entire.month = function(t1, T1 = 0, growth = AB.growth.option){
  N_total = 100000
  Apeak.vec = rsn(n = N_total, xi = xi.pos, omega = omega.pos, alpha = alpha.pos)
  Aneg.vec = rsn(n = N_total, xi = xi.neg, omega = omega.neg, alpha = alpha.neg)
  r.vec = rgamma(n = N_total,  shape = shape.decay, rate = rate.decay)

  tpeaks = rgamma(N_total , shape = (30)/2, rate = 0.5)
  tpeaks[tpeaks<T1] = T1
  t1.add.days = runif(N_total, min = 0, max = time.between.surveys-1)

  AB_vals = kadelka23cutoffFree::antibody_dynamics(AB.background = Aneg.vec, AB.peak = Apeak.vec, t.growth = T1, t.peak = tpeaks, decay.rate = r.vec, t = t1 + t1.add.days, growth.dynamics = growth)
  return(AB_vals)
}

colors <- c("0-1 month" = "red", "1-2 months" = "blue", "2-3 months" = "cyan",
            "3-4 months" = "green", "4-5 months"  = "orange", "uninfected" = "black")
colors_new <- c("0-1" = "red", "1-2" = "blue", "2-3" = "cyan",
                "3-4" = "green", "4-5"  = "orange", "uninfected" = "black")
fig1D = ggplot() +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = -60, T1 = 0), col = "uninfected"),
               size = 2) +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = 0, T1 = 0), col = "0-1"),
               size = 2) +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = 30, T1 = 0), col = "1-2"),
               size = 2) +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = 60, T1 = 0), col = "2-3"),
               size = 2) +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = 90, T1 = 0), col = "3-4"),
               size = 2) +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = 120, T1 = 0), col = "4-5"),
               size = 2) +
  xlab("log(AB level)") +
  ylab("Density") +
  scale_y_continuous(breaks = NULL, expand = expansion(add = c(0,0.1))) + # c(0,0)) +
  scale_x_continuous(breaks = NULL, limits = c(-2.7, 1.5), expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 26),
        legend.position = "top") +
  scale_color_manual(
    values = (colors_new),
    breaks = c("0-1", "1-2", "2-3", "3-4", "4-5", "uninfected")) +
  guides(color = guide_legend(title = "months since infection", title.position = "top"))

name.pic = paste0(folder_graphs,"Figure1D.pdf")
pdf(file=name.pic, height=8/2, width = 10/2)
print(fig1D)
dev.off()


###############
## Figure1: Combine panels
###############
fig1 = ggarrange(fig1A, fig1B, fig1C, fig1D, nrow =2, ncol = 2,
                 labels = c("A", "B", "C", "D"),
                 font.label = list(size = 26, color = "black")
)
name.pic = paste0(folder_graphs,"Figure1.pdf")
pdf(file=name.pic, height=8, width = 10)
print(fig1)
dev.off()


########################################
## Figure2
########################################
#### Summary graph of epidemic patterns in the seven test scenarios
test.cases.all = data.frame(matrix(NA,nrow =0,ncol = 3), check.names = FALSE)
names(test.cases.all) = c("Test case", "time", "cumulative incidence")
for (i in 1:7){
  test.case.use = eval(parse(text = paste0("test.case.",i,".a")))
  test.cases.all = rbind(test.cases.all, data.frame(`Test case` = i,
                                                    time = 1:length(test.case.use),
                                                    `cumulative incidence` = cumsum(test.case.use),
                                                    check.names = FALSE))
}

fig2 = ggplot(data = test.cases.all, aes(x = time, y = `cumulative incidence`, col = as.character(`Test case`))) +
  geom_point() +
  geom_line() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        legend.position = "top") +
  scale_x_continuous(expand=expansion(mult = 0, add = c(0,0.2))) +
  scale_y_continuous(expand=expansion(mult = 0, add = c(0,0.02))) +
  scale_color_discrete(name ="Test cases", guide = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1))
fig2

name.pic = paste0(folder_graphs,"Figure2.pdf")
pdf(file=name.pic, height=4, width = 6)
print(fig2)
dev.off()




########################################
## Figure3
########################################
AB.both = "exp"
Ttoprod.both = 0
Nstudies = 3000 # 3000 in silico studies
N.brute.force = 0

survey.sizes = c(round(10^seq(1.5,4.5,0.5)),400, 500, 600, 650, 700, 750, 800, 900,1030, 1500, 1600,
                 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100,
                 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200,
                 4500, 5000, 5500, 6000, 6100, 6200, 6300, 6400, 6500, 7000, 7200, 7300, 7400, 7500, 7600, 8000,
                 40000, 60000, 80000, 100000)  # fixed number of participants per survey
test.cases = seq(1,7)

# method in ("cutoff-based", "cutoff-free"), cohort in ("disjoint", "constant")
cohorts = c("disjoint", "constant")
methods  = c("cutoff-based", "cutoff-free")

options(scipen = 7) ### set when scientific notation is supposed to be used
# store powers in data.frame P
P = data.frame(`Test case` = test.cases, matrix(NA, nrow = length(test.cases), ncol = length(survey.sizes)),
               method = rep(methods, each = 2*length(test.cases)),
               cohort = rep(rep(cohorts, each = length(test.cases)),2), check.names = FALSE)
names(P)[2:(length(survey.sizes)+1)] = paste0("N = ",survey.sizes)  # powers will be stored here

times = list()
for(id in test.cases){
  times[[id]] = 1:length(eval(parse(text = paste0("test.case.",id,".a"))))
}
#store quantiles (2.5%, 10%, 25%, 50%, 75%, 90%, 97.5%) in data.frame Q
quantiles = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)

Q1 = data.frame(`Test case` = rep(rep(test.cases, sapply(times, FUN = length)), 4*length(survey.sizes)),
                time = rep(unlist(times), 4*length(survey.sizes)),
                matrix(NA, nrow = length(rep(rep(test.cases, sapply(times, FUN = length)), 4*length(survey.sizes))),
                       ncol = length(quantiles)),
                `survey size` = rep(rep(survey.sizes, each = sum(sapply(times, FUN = length))),4),
                method = rep(methods, each = 2*length(survey.sizes)*sum(sapply(times, FUN = length))),
                cohort = rep(rep(cohorts, each  = length(survey.sizes)*sum(sapply(times, FUN = length))),2),
                check.names = FALSE)
names(Q1)[3:(length(quantiles)+2)] = paste0(100*quantiles,"%")



folder.data.sampleSize.mixture = paste0(folder_data,"datasets/all_results_surveySize_ABdyn_exp0/")
folder.data.constantCohort.mixture = paste0(folder_data,"datasets/all_results_constant_cohort/")
folder.data.cutoff.as.in.mixture = paste0(folder_data,"datasets/all_results_cutofflow/")
folder.data.cutoff.as.in.mixture.constantCohort = paste0(folder_data,"datasets/all_results_cutofflow_constant_cohort/")


for (cohort in cohorts){
  for (method in methods){
    if(cohort == "disjoint" & method == "cutoff-free"){folder.data = folder.data.sampleSize.mixture}
    if(cohort == "disjoint" & method == "cutoff-based"){folder.data = folder.data.cutoff.as.in.mixture}
    if(cohort == "constant" & method == "cutoff-free"){folder.data = folder.data.constantCohort.mixture}
    if(cohort == "constant" & method == "cutoff-based"){folder.data = folder.data.cutoff.as.in.mixture.constantCohort}
    for(id in test.cases){
      for(survey.size in survey.sizes){
        # load cumulative INCIDENCE estimates from fitting 3000 in silico studies
        name.cum.inc.per.survey = paste0(folder.data,
                                         "TC", id,
                                         "_SurveySize", survey.size,
                                         "_Nstudies", Nstudies,
                                         "_NbruteForce", N.brute.force,
                                         "_all_studies.rds")

        if(file.exists(name.cum.inc.per.survey)){
          cum.inc.per.survey = apply(readRDS(name.cum.inc.per.survey), MARGIN  = 2, FUN = cumsum)
          true.new.inf = eval(parse(text = paste0("test.case.",id,".a")))

          power1 = kadelka23cutoffFree::calculate_power(fitted_cumulative_incidences = cum.inc.per.survey,
                                                        test.case = true.new.inf)

          quantiles1 = apply(100*cum.inc.per.survey, MARGIN = 1, FUN = quantile, probs = quantiles)

          P[which(P$`Test case`==id & P$method == method & P$cohort == cohort),
            which(gsub("N = ", "", names(P)) == survey.size)] = round(100*power1[[1]],2)

          Q1[which(Q1$`Test case`==id & Q1$method == method & Q1$cohort == cohort & Q1$`survey size`==survey.size),
             3:(length(quantiles)+2)] = round(t(quantiles1),2)
        }
      }
    }
  }
}


P1 = P %>% filter(method == "cutoff-based" & cohort == "disjoint" & `Test case` %in% 1:7) %>%
  dplyr::select(c("Test case", paste0("N = ",sort(survey.sizes)), "method", "cohort"))
cols = as.numeric(which(sapply(P1, function(x)!all(is.na(x)))))
P1[,cols]



TC = 1:7
surveysize = 1000
colors1 = c("cutoff-free, disjoint" = "red",
            "cutoff-free, constant" = "blue",
            "cutoff-based, disjoint" = "orange",
            "cutoff-based, constant" = "green")
fig3 = ggplot(Q1 %>% filter(!is.na(`50%`) & `Test case` %in% TC &
                              `survey size`  == surveysize &
                              (cohort != "constant" | method != "cutoff-based")) %>%
                mutate(true_cum_inc = rep(rep(unlist(lapply(lapply(parse(text = paste0("test.case.",TC,".a")), FUN = eval),
                                                            FUN = function(x){100*cumsum(x)})),
                                              length(surveysize)),length(unique(paste(method,cohort))) )),
              aes(x = time,
                  y = `50%`,
                  ymin=`2.5%`,
                  ymax = `97.5%`,
                  col = paste(method,cohort, sep = ", "),
                  fill = paste(method,cohort, sep = ", "))) +
  geom_line() +
  geom_ribbon(size = 0, alpha = 0.1) +
  facet_wrap(~paste0("TC = ",`Test case`), scales = "free",
             ncol = 2) +
  scale_color_manual(name = "Method, cohort",
                     values = colors1,
                     breaks = c("cutoff-free, disjoint",
                                "cutoff-free, constant",
                                "cutoff-based, disjoint",
                                "cutoff-based, constant")) +
  scale_fill_manual(name = "Method, cohort",
                    values = colors1,
                    breaks = c("cutoff-free, disjoint",
                               "cutoff-free, constant",
                               "cutoff-based, disjoint",
                               "cutoff-based, constant")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        legend.position = "top")+
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2),
         fill = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2)) +
  xlab("Month") +
  ylab("Cumulative incidence (in %)") +
  geom_point(aes(x = time, y = true_cum_inc), col = "black", pch = 4, size = 2, show.legend = FALSE)
fig3
name.pic = paste0(folder_graphs,"Figure3.pdf")
pdf(file=name.pic, height=10, width = 8)
print(fig3)
dev.off()




########################################
## Figure4
########################################
colors2 = c("cutoff-free, disjoint" = rgb(1,0,0,1),
            "cutoff-free, constant" = rgb(0,0,1,1),
            "cutoff-based, disjoint" = rgb(1,0,0,0.5),
            "cutoff-based, constant" = rgb(0,0,1,0.5))
linetypes = c("cutoff-free" = "solid", "cutoff-based" = "dashed")

D.plot = reshape2::melt(P, id.vars = c("Test case", "method", "cohort"))%>%
  filter(!is.na(value))
D.plot = D.plot %>% filter(`Test case` %in% 1:7)
fig4 = ggplot() +
  scale_x_log10() +
  geom_point(data = D.plot %>% filter(method != "cutoff-based" | cohort != "constant"),
             aes(x = as.numeric(gsub("N = ", "", variable)),
                 y = value,
                 col = paste(method,cohort, sep = ", "),
                 linetype = method), size = 2) +
  geom_line(data = D.plot %>% filter(method != "cutoff-based" | cohort != "constant"),
            aes(x = as.numeric(gsub("N = ", "", variable)),
                y = value,
                col = paste(method,cohort, sep = ", "),
                linetype = method)) +
  facet_wrap(~`Test case`, scales = "free_x", ncol = 2) +
  xlab("Sampled individuals per survey") +
  ylab("Power (in %)") +
  scale_color_manual(name = "Method, cohort",
                     values = colors2,
                     breaks = c("cutoff-free, disjoint",
                                "cutoff-free, constant",
                                "cutoff-based, disjoint",
                                "cutoff-based, constant"))  +
  scale_linetype_manual(
    name = "Method",
    values = linetypes,
    breaks = c("cutoff-free", "cutoff-based")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        legend.position = "top") +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2),
         linetype = FALSE) +
  geom_hline(data = data.frame(x = 1, y = 99), aes(yintercept = 99), col = "black", linetype = "dashed")
fig4

name.pic = paste0(folder_graphs,"Figure4.pdf")
pdf(file=name.pic, height=10, width = 7)
print(fig4)
dev.off()



########################################
## FigureS4
########################################
##############
## FigureS4A
##############
D.plot = reshape2::melt(P, id.vars = c("Test case", "method", "cohort"))%>%
  filter(!is.na(value))
D.plot = D.plot %>% filter(`Test case` %in% 1:5)
figS4A = ggplot() +
  xlim(c(NA,10000)) +
  coord_cartesian(ylim = c(95,100)) +
  geom_point(data = D.plot %>% filter(method != "cutoff-based" | cohort != "constant"),
             aes(x = as.numeric(gsub("N = ", "", variable)),
                 y = value,
                 col = paste(method,cohort, sep = ", "),
                 linetype = method), size = 2) +
  geom_line(data = D.plot %>% filter(method != "cutoff-based" | cohort != "constant"),
            aes(x = as.numeric(gsub("N = ", "", variable)),
                y = value,
                col = paste(method,cohort, sep = ", "),
                linetype = method)) +
  facet_wrap(~`Test case`, scales = "free_x", ncol = 2)+
  xlab("Sampled individuals per survey") +
  ylab("Power (in %)") +
  scale_color_manual(name = "Method, cohort",
                     values = colors2,
                     breaks = c("cutoff-free, disjoint",
                                "cutoff-free, constant",
                                "cutoff-based, disjoint",
                                "cutoff-based, constant"))  +
  scale_linetype_manual(
    name = "Method",
    values = linetypes,
    breaks = c("cutoff-free", "cutoff-based")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 14),
        legend.position = "top") +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2),
         linetype = FALSE) +
  geom_hline(data = data.frame(x = 1, y = 99), aes(yintercept = 99), col = "black", linetype = "dashed")

name.pic = paste0(folder_graphs,"S4AFigure.pdf")
pdf(file=name.pic,  height=7, width=7)
print(figS4A)
dev.off()

##############
## FigureS4B
##############
D.plot = reshape2::melt(P, id.vars = c("Test case", "method", "cohort"))%>%
  filter(!is.na(value))
D.plot = D.plot %>% filter(`Test case` %in% 6:7)
figS4B = ggplot() +
  coord_cartesian(ylim = c(80,100)) +
  geom_point(data = D.plot %>% filter(method != "cutoff-based" | cohort != "constant"),
             aes(x = as.numeric(gsub("N = ", "", variable)),
                 y = value,
                 col = paste(method,cohort, sep = ", "),
                 linetype = method), size = 2) +
  geom_line(data = D.plot %>% filter(method != "cutoff-based" | cohort != "constant"),
            aes(x = as.numeric(gsub("N = ", "", variable)),
                y = value,
                col = paste(method,cohort, sep = ", "),
                linetype = method)) +
  facet_wrap(~`Test case`, scales = "free_x", ncol = 2)+
  xlab("Sampled individuals per survey") +
  ylab("Power (in %)") +
  scale_color_manual(name = "Method, cohort",
                     values = colors2,
                     breaks = c("cutoff-free, disjoint",
                                "cutoff-free, constant",
                                "cutoff-based, disjoint",
                                "cutoff-based, constant"))  +
  scale_linetype_manual(
    name = "Method",
    values = linetypes,
    breaks = c("cutoff-free", "cutoff-based")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 14),
        legend.position = "top") +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2),
         linetype = FALSE) +
  geom_hline(data = data.frame(x = 1, y = 99), aes(yintercept = 99), col = "black", linetype = "dashed")
figS4B

name.pic = paste0(folder_graphs,"S4BFigure.pdf")
pdf(file=name.pic,  height=3.5, width=7)
print(figS4B)
dev.off()



########################################
## Figure5
########################################
####################### SURVEY DATA from PRETE (8 cities), NEW (PRETE) OR OLD (BUSS) VALIDATION DATA
# WEEKLY VS MONTHLY
# SURVEY DATA: UNTIL OCTOBER 2020 vs UNTIIL NOVEMBER 2020
blood_donors = readRDS(paste0(folder_data,"datasets/Bloodbank.rds"))

cities = c("HEMEPAR","FPS","HEMOPE", "HEMOAM", "HEMOBA", "HEMOCE", "HEMORIO", "HEMOMINAS")
colnames1 = c('2.5%','50%','97.5%','growth dynamics','point estimate','time','repeated','oldValidationData','maxmonth','city')
output.summary.all = data.frame(matrix(nrow = 0, ncol = length(colnames1)))
names(output.summary.all) = colnames1
for(city in cities){
  for (maxmonth in c(10,11)){
    survey = (blood_donors%>%filter(blood_center == city, month %in% 3:maxmonth))
    min.week = min(survey$week)
    for (monthly_weekly in c('weekly','monthly')){
      for (oldValidationData in c("FALSE","TRUE")){
        for (ABgrowth1 in c('linear','exp')){
          for (ttoprod1 in c(0,7)){
            filename1 = paste0(folder_data,"datasets/all_results_",monthly_weekly,"_oldValidationData_",
                               oldValidationData,"/",city,"_ABdynTrue",ABgrowth1,ttoprod1,
                               "_Nstudies1000_NbruteForce0_maxmonth",maxmonth,".rds")
            if(file.exists(filename1)){
              output1 = readRDS(file = filename1)
              filename2 = paste0(folder_data,"datasets/all_results_",monthly_weekly,"_oldValidationData_",
                                 oldValidationData,"/",city,"_ABdynTrue",ABgrowth1,ttoprod1,
                                 "_Nstudies1000_NbruteForce0_maxmonth",maxmonth,"_point_estimate.rds")
              output.original = readRDS(file = filename2)
              if(monthly_weekly == "weekly"){
                approx.date.surveys  = as.Date(paste(2020, min.week:(min.week+length(output.original)-1), 1, sep="-"), "%Y-%U-%u")
              }
              if(monthly_weekly == "monthly"){
                approx.date.surveys = seq(as.Date("2020-04-01"),length=maxmonth-2,by="months")-1### last day of month
              }
              output1.summary = data.frame(t(
                apply(apply(output1, MARGIN = 2, FUN = cumsum), MARGIN = 1, FUN = quantile, probs = c(0.025, 0.5, 0.975))),
                check.names = FALSE)%>%mutate(`growth dynamics` = paste0(ABgrowth1,ttoprod1))
              output1.summary$`point estimate` = cumsum(output.original)
              output1.summary$time = approx.date.surveys
              output1.summary$repeated = monthly_weekly
              output1.summary$oldValidationData = oldValidationData
              output1.summary$maxmonth = maxmonth
              output1.summary$city = city
              output.summary.all = rbind(output.summary.all, output1.summary)
            }
          }
        }
      }
    }
  }
}


plot.cities = c("HEMEPAR","FPS","HEMOPE", "HEMOAM", "HEMOBA", "HEMOCE", "HEMORIO", "HEMOMINAS")
growths.assumed = c("linear0","linear7","exp0","exp7")
y.limits = list(c())
y.limits.map <- list("HEMOCE" = c(0,50),
                     "HEMOAM" = c(0,70),
                     "HEMOBA" = c(0,32),
                     "HEMORIO" = c(0,35),
                     "HEMOPE" =c(0,42),
                     "HEMOMINAS" = c(0,15),
                     "HEMEPAR" = c(0,11),
                     "FPS" = c(0,25))
maxmonth.plot = 11
pic.names = c("exp0" = "A", "exp7" = "B", "linear0" = "C", "linear7" = "D")
for(plot.city in c("HEMOAM")){#plot.cities){
  for(growth.assumed in growths.assumed){
    mortality_data = readRDS(file = paste0(folder_data,"datasets/mortality_data_from_mortality_rate_code.rds"))

    citymap <- c("HEMOCE" = "Fortaleza", "HEMOAM" = "Manaus", "HEMOBA" = "Salvador", "HEMORIO" = "Rio de Janeiro",
                 "HEMOPE" = "Recife", "HEMOMINAS" = "Belo Horizonte", "HEMEPAR" = "Curitiba", "FPS" = "Sao Paulo")
    citymap.rev <- c("Fortaleza" = "HEMOCE", "Manaus" = "HEMOAM", "Salvador" = "HEMOBA", "Rio de Janeiro" = "HEMORIO",
                     "Recife" = "HEMOPE", "Belo Horizonte" = "HEMOMINAS", "Curitiba" = "HEMEPAR",  "Sao Paulo" = "FPS")

    cities_plot = c("Fortaleza", "Manaus", "Salvador", "Rio de Janeiro",
                    "Recife", "Belo Horizonte", "Curitiba", "Sao Paulo")
    cities = plot.city
    D.mortality = mortality_data %>%
      filter(DT_SIN_PRI <= as.Date("2020-11-30"))%>%
      filter(city %in% cities_plot) %>%
      mutate(location = citymap.rev[city]) %>%
      filter(location %in% cities)

    D.plot = output.summary.all %>% filter(maxmonth == maxmonth.plot, city == plot.city)
    D.plot = D.plot %>% filter(time <= as.Date("2020-10-31"), `growth dynamics` == growth.assumed) %>% mutate(city = citymap[city])
    ### for monthly estimation use mid (15th) of month
    h1 = D.plot$time[D.plot$repeated == "monthly"]
    D.plot$time[D.plot$repeated == "monthly"] = paste0(year(h1),"-",month(h1),"-",15)
    D.plot$oldValidationData[D.plot$oldValidationData == TRUE] = "(a) convalescents, PCR positive" #"Buss (convalescents, PCR positive)"
    D.plot$oldValidationData[D.plot$oldValidationData == FALSE] = "(b) AB positive repeat blood donors" #"Prete (AB positive repeat blood donors)"

    fig5 = ggplot() +
      labs(x = "Date in 2020", y = "Mortality per million inhabitants,\nCumulative incidence (in %)") +
      geom_ribbon(data = D.mortality, aes(x = DT_SIN_PRI, ymin = 0, ymax= 1E6*MR_mean, group = city), alpha = 0.1, col = rgb(0,0,0,0.3)) +
      geom_point(data = D.plot, aes(x = time, y = 100*`point estimate`, pch = `oldValidationData`)) +
      geom_line(data = D.plot, aes(x = time, y = 100*`50%`, linetype = `oldValidationData`))+
      geom_ribbon(data = D.plot,
                  aes(x = time, ymin = 100*`2.5%`, ymax = 100*`97.5%`,
                      group = `oldValidationData`),
                  alpha = 0.2) +
      geom_line(data = D.plot, aes(x = time, y = 100*`2.5%`,
                                   linetype = `oldValidationData`),
                alpha = 0.3) +
      geom_line(data = D.plot, aes(x = time, y = 100*`97.5%`,
                                   linetype = `oldValidationData`),
                alpha = 0.3) +
      theme_bw() +
      theme(
        legend.position = "top",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 16)) +
      ylim(as.numeric(unlist(y.limits.map[plot.city]))) +
      facet_wrap(~paste0("estimated ",repeated),nrow = 2) +
      scale_linetype_discrete(name = "Validation data", guide = guide_legend(title.position = "left", nrow = 2)) +
      scale_shape_discrete(name = "Validation data", guide = guide_legend(title.position = "left", nrow = 2)) #+

    name.pic = paste0(folder_graphs,"Figure5",as.character(pic.names[growth.assumed]),".pdf")
    pdf(file=name.pic, height=8, width=8)
    print(fig5)
    dev.off()

  }
}



########################################
## FigureS1 - S3
########################################
TCs = c(1,3,6)
figuresS = c(1,2,3)
counter = 1
for (TC in TCs){
  ## read in data from cutoff-free constant
  point_estim = readRDS(paste0(folder_data, "datasets/all_results_boot_constant/all_point_estimates_TC_",TC,"_ABdynTrueexp0_SurveySize1000_Nstudies50.rds"))
  cum_inc_point_esti = apply(point_estim,MARGIN = 2, FUN = cumsum)

  Nstudies = ncol(point_estim)
  length_CIs = matrix(NA, nrow = ncol(point_estim), ncol = nrow(point_estim))

  for(i in 1:Nstudies){
    boot_i = readRDS(paste0(folder_data, "datasets/all_results_boot_constant/bootstrap_study_",i,"_TC_",TC,"_ABdynTrueexp0_SurveySize1000_Nstudies50.rds"))
    CI_i = apply(apply(boot_i, MARGIN = 2, FUN = cumsum), MARGIN = 1, FUN = quantile, probs = c(0.05/2/9, 1-0.05/2/9)) #probs = c(0.025, 0.975))
    length_CI_i = as.numeric(CI_i[2,] - CI_i[1,])
    length_CIs[i,] = length_CI_i
  }

  length_CIs

  L.constant = stack(data.frame(length_CIs)) %>% mutate(survey = gsub("X","",ind))%>% dplyr::select(c(values,survey))
  PE.constant = stack(data.frame(t(cum_inc_point_esti))) %>% mutate(survey = gsub("X","",ind))%>% dplyr::select(c(values,survey))


  ############
  ## read in data from cutoff-based and cutoff-free disjoint
  ############
  point_estim = readRDS(paste0(folder_data, "datasets/all_results_boot_cutoff/all_point_estimates_TC_",TC,"_ABdynTrueexp0_SurveySize1000_Nstudies50.rds"))
  cum_inc_point_esti = apply(point_estim,MARGIN = 2, FUN = cumsum)

  point_estim_COfree = readRDS(paste0(folder_data, "datasets/all_results_boot_cutoff/all_point_estimates_TC_",TC,"_ABdynTrueexp0_SurveySize1000_Nstudies50_cutoff_free_disjoint_fitted_simultaneously_with_cutoff_based.rds"))
  cum_inc_point_esti_COfree = apply(point_estim_COfree,MARGIN = 2, FUN = cumsum)

  Nstudies = ncol(point_estim)
  length_CIs = matrix(NA, nrow = ncol(point_estim), ncol = nrow(point_estim))
  length_CIs_COfree = matrix(NA, nrow = ncol(point_estim_COfree), ncol = nrow(point_estim_COfree))

  for(i in 1:Nstudies){
    boot_i = readRDS(paste0(folder_data,"datasets/all_results_boot_cutoff/bootstrap_study_",i,"_TC_",TC,"_ABdynTrueexp0_SurveySize1000_Nstudies50.rds"))
    CI_i = apply(apply(boot_i, MARGIN = 2, FUN = cumsum), MARGIN = 1, FUN = quantile, probs = c(0.05/2/9, 1-0.05/2/9)) #probs = c(0.025, 0.975))
    length_CI_i = as.numeric(CI_i[2,] - CI_i[1,])
    length_CIs[i,] = length_CI_i

    boot_i = readRDS(paste0(folder_data,"datasets/all_results_boot_cutoff/bootstrap_study_",i,"_TC_",TC,"_ABdynTrueexp0_SurveySize1000_Nstudies50_cutoff_free_disjoint_fitted_simultaneously_with_cutoff_based.rds"))
    CI_i = apply(apply(boot_i, MARGIN = 2, FUN = cumsum), MARGIN = 1, FUN = quantile,  probs = c(0.05/2/9, 1-0.05/2/9)) #probs = c(0.025, 0.975))
    length_CI_i = as.numeric(CI_i[2,] - CI_i[1,])
    length_CIs_COfree[i,] = length_CI_i

  }

  L.cutoff = stack(data.frame(length_CIs)) %>% mutate(survey = gsub("X","",ind))%>% dplyr::select(c(values,survey))
  PE.cutoff = stack(data.frame(t(cum_inc_point_esti))) %>% mutate(survey = gsub("X","",ind))%>% dplyr::select(c(values,survey))

  L.disjoint= stack(data.frame(length_CIs_COfree)) %>% mutate(survey = gsub("X","",ind))%>% dplyr::select(c(values,survey))
  PE.disjoint = stack(data.frame(t(cum_inc_point_esti_COfree))) %>% mutate(survey = gsub("X","",ind))%>% dplyr::select(c(values,survey))



  L = rbind(L.disjoint %>% mutate(method = "disjoint"),
            L.constant %>% mutate(method = "constant"),
            L.cutoff %>% mutate(method = "cutoff")
  )
  PE = rbind(PE.disjoint %>% mutate(method = "disjoint"),
             PE.constant %>% mutate(method = "constant"),
             PE.cutoff %>% mutate(method = "cutoff")
  )

  cum.inc.true = cumsum(eval(parse(text = paste0("test.case.", TC, ".a"))))

  my_comparisons1 <- list( c("CO-free\ndisjoint", "CO-free\nconstant"), c("CO-free\ndisjoint", "CO-based\ndisjoint"), c("CO-free\nconstant", "CO-based\ndisjoint") )

  # ### Test for normality
  # for(method1 in c("disjoint", "constant", "cutoff")){
  #   for(i in 1:9){
  #     normal_data = L%>%filter(survey == i , method == method1)%>% .$values
  #     shapiro.test(normal_data) %>% print()
  #     normal_data = PE%>%filter(survey == i , method == method1)%>% .$values
  #     shapiro.test(normal_data) %>% print()
  #   }
  # }

  p1.a = ggplot(data = PE, aes(x = factor(method, levels = c("disjoint", "constant", "cutoff"),
                                          labels = c("CO-free\ndisjoint", "CO-free\nconstant", "CO-based\ndisjoint")),
                               y = values,
                               fill = factor(method, levels = c("disjoint", "constant", "cutoff"),
                                             labels = c("CO-free\ndisjoint", "CO-free\nconstant", "CO-based\ndisjoint")))) +
    facet_wrap(~survey, nrow = 1) +
    geom_violin(
      alpha = 0.4,
      position = position_dodge(width = 0.7),
      width = 0.7) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_blank(),
          legend.position = "top")  +
    xlab("") +
    ylab("Length of confidence interval") +
    geom_hline(data = data.frame(CI = cum.inc.true) %>% mutate(survey = 1:nrow(.)),
               aes(yintercept = CI), linetype = "dashed") +
    stat_compare_means(comparisons = my_comparisons1,
                       method = "t.test",
                       method.args = list(alternative = "two.sided"),
                       label = "p.signif",
                       vjust = 0.2,
                       size = 3) +
    stat_summary(data = PE,
                 aes(x = factor(method, levels = c("disjoint", "constant", "cutoff"),
                                labels = c("CO-free\ndisjoint", "CO-free\nconstant", "CO-based\ndisjoint")),
                     y = values,
                     group = factor(method, levels = c("disjoint", "constant", "cutoff"),
                                    labels = c("CO-free\ndisjoint", "CO-free\nconstant", "CO-based\ndisjoint"))),
                 fun = "mean",
                 geom = "crossbar",
                 width = 0.5,
                 colour = "red",
                 position = position_dodge(width = 0.7),
                 show.legend = FALSE)+
    guides(fill = guide_legend(title = "Method"))


  p2.a = ggplot(data = L, aes(x = factor(method, levels = c("disjoint", "constant", "cutoff"),labels = c("CO-free\ndisjoint", "CO-free\nconstant", "CO-based\ndisjoint")), y = values, fill = factor(method, levels = c("disjoint", "constant", "cutoff"),labels = c("CO-free\ndisjoint", "CO-free\nconstant", "CO-based\ndisjoint")))) +
    facet_wrap(~survey, nrow = 1) +
    geom_violin(
      alpha = 0.4,
      position = position_dodge(width = 0.7),
      width = 0.7) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text.x = element_blank(),
          legend.position = "top")  +
    xlab("") +
    ylab("Length of confidence interval") +
    stat_compare_means(comparisons = my_comparisons1,
                       method = "t.test",
                       method.args = list(alternative = "two.sided"),
                       label = "p.signif",
                       vjust = 0.2,
                       size = 3) +
    stat_summary(data = L,
                 aes(x = factor(method, levels = c("disjoint", "constant", "cutoff"),
                                labels = c("CO-free\ndisjoint", "CO-free\nconstant", "CO-based\ndisjoint")),
                     y = values,
                     group = factor(method, levels = c("disjoint", "constant", "cutoff"),
                                    labels = c("CO-free\ndisjoint", "CO-free\nconstant", "CO-based\ndisjoint"))),
                 fun = "mean",
                 geom = "crossbar",
                 width = 0.5,
                 colour = "red",
                 position = position_dodge(width = 0.7),
                 show.legend = FALSE)+
    guides(fill = guide_legend(title = "Method"))

  pic1.a = ggarrange(p1.a + ylab("Point estimtes of\ncumulative incidences"), p2.a + ylab("Sidelength of\nconfidence region"), nrow = 2, common.legend = TRUE)
  name.pic = paste0(folder_graphs,"S",figuresS[counter],"Figure.pdf")
  pdf(file=name.pic, height=6, width=10)
  print(pic1.a)
  dev.off()

  counter = counter + 1
}


########################################
# Figure S5
########################################
################
# Figure S5A
################
N.ABs = 1e+4
cutoffs = seq(-3.2,2,0.001)

calc_ROC_AUC = function(mean_peaks, cutoffs = seq(-3.2,2,0.001), N.ABs = 1e+4){
  A.pos1 = rsn(N.ABs, xi = mean_peaks, omega = omega.pos, alpha = alpha.pos)
  A.neg = rsn(N.ABs, xi = xi.neg, omega = omega.neg, alpha = alpha.neg)
  FPR = sapply(cutoffs, FUN = function(x){sum(A.neg > x)/length(A.neg)})
  TPR1 = sapply(cutoffs, FUN = function(x){sum(A.pos1 > x)/length(A.pos1)})
  if(min(FPR) == 0 & max(FPR) == 1){
    AUC_ROC1 = trapz(x = rev(FPR), y = rev(TPR1))
    return(AUC_ROC1)
  }else{
    print(mean_peaks)
    min(A.neg)
    max(A.neg)
    print(min(FPR))
    print(max(FPR))
    print("Error: change cutoffs input.")
  }
}

mean.peaks = c(0.9, 0.6, 0.3, 0, -0.3)
ROC.AUCs = sapply(mean.peaks, FUN = calc_ROC_AUC)

ROC_AUCs = data.frame(mean_peak = mean.peaks,
                      ROC_AUC = ROC.AUCs)


all_cum_inc = data.frame(matrix(NA, nrow = 0, ncol = 11), check.names = FALSE)
names(all_cum_inc) = c("test case", "survey", "mean peak", "ROC-AUC", "survey size", "50%", "2.5%", "97.5%", "true cumulative incidence", "power", "method")

for(TC in c(1,3,6)){
  r.test = eval(parse(text = paste0("test.case.",TC,".a")))
  for(mean.use in mean.peaks){
    for(Nsurveysize in c(1000, 2000, 3000, 5000, 10000)){
      filename1 = paste0(folder_data,"datasets/all_results_sn_dist_mixture_plus_cutoff/all_res_TC_",TC,"_ABdynTrueexp0_SurveySize",Nsurveysize,"_Nstudies1000_sn_location_",mean.use,".rds")
      if(file.exists(filename1)){
        A = readRDS(filename1)
        power = kadelka23cutoffFree::calculate_power(apply(A, MARGIN = 2, FUN = cumsum),test.case =r.test)
        A.add = data.frame(t(apply(apply(A, MARGIN = 2, FUN = cumsum),
                                   MARGIN = 1, FUN = quantile, probs = c(0.025, 0.5, 0.975))),
                           check.names = FALSE) %>%
          dplyr::mutate(`test case` = TC,
                        survey = seq(1,nrow(.)),
                        `mean peak` = mean.use,
                        `ROC-AUC` = ROC_AUCs$ROC_AUC[ROC_AUCs$mean_peak == mean.use],
                        `survey size` = Nsurveysize,
                        `true cumulative incidence` = cumsum(r.test),
                        power = power,
                        method = "cutoff-free, disjoint") %>%
          dplyr::select(names(all_cum_inc))
        all_cum_inc = rbind(all_cum_inc, A.add)
      }
      filename2 = paste0(folder_data,"datasets/all_results_sn_dist_mixture_plus_cutoff/all_res_TC_",TC,"_ABdynTrueexp0_SurveySize",Nsurveysize,"_Nstudies1000_sn_location_",mean.use,"_cutoff.rds")
      if(file.exists(filename2)){
        A = readRDS(filename2)
        power = kadelka23cutoffFree::calculate_power(apply(A, MARGIN = 2, FUN = cumsum),test.case =r.test)
        A.add = data.frame(t(apply(apply(A, MARGIN = 2, FUN = cumsum),
                                   MARGIN = 1, FUN = quantile, probs = c(0.025, 0.5, 0.975))),
                           check.names = FALSE) %>%
          dplyr::mutate(`test case` = TC,
                        survey = seq(1,nrow(.)),
                        `mean peak` = mean.use,
                        `ROC-AUC` = ROC_AUCs$ROC_AUC[ROC_AUCs$mean_peak == mean.use],
                        `survey size` = Nsurveysize,
                        `true cumulative incidence` = cumsum(r.test),
                        power = power,
                        method = "cutoff") %>%
          dplyr::select(names(all_cum_inc))
        all_cum_inc = rbind(all_cum_inc, A.add)
      }
    }
  }
}

M.AUC = all_cum_inc %>%
  filter(survey == 1)
figS5A = ggplot(data = M.AUC,
                aes(x=`ROC-AUC`,y=power, col=as.character(`survey size`), linetype = method)) +
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 16),
        legend.position = "top") +
  facet_wrap(~paste0("Test case ",`test case`), nrow = 1, scales = "free_y") +
  guides(color = guide_legend(title = "Individuals per survey", title.hjust = 0.5,
                              title.position = "top", nrow = 2),
         linetype = guide_legend(title = "Method", title.hjust = 0.5,
                                 title.position = "top", nrow = 2)) +
  scale_linetype_manual(breaks = c("cutoff-free, disjoint", "cutoff"),
                        values=c("solid", "dashed"),
                        labels = c("CO-free, disjoint", "CO-based")) +
  ylab("Power")

name.pic = paste0(folder_graphs,"S5AFigure.pdf")
pdf(file=name.pic, height=4, width=8)
print(figS5A)
dev.off()


################
# Figure S5B
################
a.COfree = all_cum_inc %>%
  filter(survey == 1,
         method == "cutoff-free, disjoint") %>%
  dplyr::select(c("test case", "mean peak", "survey size", "power", "method", "ROC-AUC"))

a.cutoff = all_cum_inc %>%
  filter(survey == 1,
         method == "cutoff") %>%
  dplyr::select(c("test case", "mean peak", "survey size", "power", "method","ROC-AUC"))

A = a.COfree %>% dplyr::mutate(powerCOfree = power)
A = A%>%dplyr::select(-c("power","method"))
A$powerCObased = a.cutoff$power

A$power_rel_improvement_COfree_vs_CObased = (A$powerCOfree - A$powerCObased)/A$powerCObased
A$power_absolute_improvement_COfree_vs_CObased = (A$powerCOfree - A$powerCObased)


### relative improvement is infinite when power using cutoff based method is 0
###### --> do not include in plots
figS5B = ggplot(data = A %>% filter(is.finite(power_rel_improvement_COfree_vs_CObased)),
                aes(x=`ROC-AUC`,y=100*power_rel_improvement_COfree_vs_CObased,
                    col=as.character(`survey size`))) +
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 16),
        legend.position = "top") +
  facet_wrap(~paste0("Test case ",`test case`), nrow = 1, scales = "free_y") +
  guides(color = guide_legend(title = "Individuals per survey", title.hjust = 0.5,
                              title.position = "top", nrow = 2)) +
  ylab("Relative change of power (in %)")

name.pic = paste0(folder_graphs,"S5BFigure.pdf")
pdf(file=name.pic, height=4, width=8)
print(figS5B +  scale_y_log10())
dev.off()



########################################
# Figure S6-S8
########################################
folder.data.disjoint.mixture.diffGrowth = paste0(folder_data,"datasets/all_results_misspec/")

model_misspec_diff_growth = function(id, survey.size, save.graph = FALSE, save.table = FALSE,
                                     figure.name, table.name ){
  ABgrowth = c("exp","linear")
  ttoprod = c(0,7)

  all_powers = data.frame(matrix(NA,nrow=0,ncol=4), check.names = FALSE)
  names(all_powers) = c("Test case", "True dynamics", "Assumed dynamics","power")

  all_powers
  quantiles = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
  all_quantiles = data.frame(matrix(NA,nrow=0,ncol=4+length(quantiles)), check.names = FALSE)
  names(all_quantiles) = c("Test case", "time", "True dynamics", "Assumed dynamics", paste0(100*quantiles,"%"))


  for (dynTrue in ABgrowth){
    for (tTrue in ttoprod){
      for (dynAssumed in ABgrowth){
        for (tAssumed in ttoprod){
          name.df = paste0(folder.data.disjoint.mixture.diffGrowth,"all_res_TC_",id,"_ABdynTrue",dynTrue,tTrue,"_ABdynAssumed",dynAssumed,tAssumed,"_SurveySize",survey.size,"_Nstudies3000.rds")
          if(file.exists(name.df)){
            cum.inc.per.survey = 100*apply(readRDS(name.df), MARGIN  = 2, FUN = cumsum)
            true.cum.inc = cumsum(eval(parse(text = paste0("test.case.",id,".a"))))
            true.test.case = eval(parse(text = paste0("test.case.",id,".a")))

            power1 = kadelka23cutoffFree::calculate_power(fitted_cumulative_incidences = cum.inc.per.survey/100,
                                                          test.case = true.test.case)
            quantiles1 = apply(cum.inc.per.survey, MARGIN = 1, FUN = quantile, probs = quantiles)
            all_powers = rbind(all_powers, data.frame(`Test case` = id, `True dynamics` = paste0(dynTrue,tTrue),
                                                      `Assumed dynamics` = paste0(dynAssumed,tAssumed), power = round(100*power1[[1]],2), check.names = FALSE))
            all_quantiles = rbind(all_quantiles, data.frame(`Test case` = id, time = 1:nrow(t(quantiles1)), `True dynamics` = paste0(dynTrue,tTrue),
                                                            `Assumed dynamics` = paste0(dynAssumed,tAssumed), t(quantiles1), check.names = FALSE))


          }
        }
      }
    }
  }

  all_powers
  powers.show = reshape2::dcast(all_powers, `True dynamics` ~ `Assumed dynamics`)
  all_quantiles %>% dplyr::select(c("Test case", "time", "True dynamics", "Assumed dynamics", "2.5%", "50%", "97.5%"))

  test.case.use = eval(parse(text = paste0("test.case.",id,".a")))
  p1 = ggplot() +
    geom_line(data = all_quantiles, aes(x = time, y = `50%`)) +
    geom_ribbon(data = all_quantiles, aes(x = time, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
    geom_point(data = all_quantiles%>%mutate(true_cum_inc = rep(100*cumsum(test.case.use), length(time)/length(test.case.use))),
               aes(x = time, y =true_cum_inc), pch = 4, size = 2) +
    facet_grid(paste0("True = ",`True dynamics`)~paste0("Assumed = ",`Assumed dynamics`)) +
    xlab("Survey") +
    ylab("Cumulative incidence (in %)") +
    theme_bw() +
    theme(text = element_text(size = 16))
  print(p1)
  if (save.graph){
    name.pic = paste0(folder_graphs,figure.name,".pdf")
    pdf(file=name.pic, height=9, width=10)
    print(p1)
    dev.off()
  }
  if (save.table){
    print(xtable(powers.show),
          file = paste0(folder_graphs,table.name,".tex"), include.rownames = FALSE)
  }
  return(list(all_powers, powers.show, all_quantiles))
}


model_misspec_diff_growth(id = 1, survey.size = 800, save.graph = TRUE, save.table = TRUE, figure.name = "S6Figure", table.name = "TableS1")
model_misspec_diff_growth(id = 3, survey.size = 2500, save.graph = TRUE, save.table = TRUE, figure.name = "S7Figure", table.name = "TableS2")
model_misspec_diff_growth(id = 6, survey.size = 10000, save.graph = TRUE, save.table = TRUE, figure.name = "S8Figure", table.name = "TableS3")


########################################
# Figure S9
########################################
colors_new3 <- c("uninfected" = "grey", "0-1 mpi, exp0" = "red",  "0-1 mpi, exp7" = "blue",
                 "0-1 mpi, linear0" = "cyan", "0-1 mpi, linear7" = "orange", "peaks" = "black")
figS9 = ggplot() +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = -60, T1 = 0), col = "uninfected"),
               size = 2) +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = 0, T1 = 0, growth = "exp"), col = "0-1 mpi, exp0"),
               size = 2) +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = 0, T1 = 7, growth = "exp"), col = "0-1 mpi, exp7"),
               size = 2) +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = 0, T1 = 0, growth = "linear"), col = "0-1 mpi, linear0"),
               size = 2) +
  geom_density(aes(x = f.pos.AB.vals.entire.month(t1 = 0, T1 = 7, growth = "linear"), col = "0-1 mpi, linear7"),
               size = 2) +
  geom_density(aes(x = rsn(n = 100000, xi = xi.pos, omega = omega.pos, alpha = alpha.pos), col = "peaks"), size = 2) +
  xlab("log(AB level)") +
  ylab("Density") +
  scale_y_continuous(breaks = NULL, expand = expansion(mult=c(0,0.05),add = 0)) +
  scale_x_continuous(breaks = NULL, limits = c(-2.7, 1.5), expand = c(0,0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 26),
        legend.position = "top") +
  scale_color_manual(
    values = (colors_new3),
    breaks = c( "uninfected",
                "peaks",
                "0-1 mpi, exp0",
                "0-1 mpi, exp7",
                "0-1 mpi, linear0",
                "0-1 mpi, linear7"
    )) +
  guides(color = guide_legend(title = ""))

name.pic = paste0(folder_graphs,"S9Figure.pdf")
pdf(file=name.pic, height=6, width = 10)
print(figS9)
dev.off()



########################################
# Figure S10 A-C
########################################
folder.data.shifted.peaks = paste0(folder_data,"datasets/all_results_misspec_shifted_peaks/")

model_misspec_shifted_peaks = function(id, survey.size, save.graph = FALSE, save.table = FALSE,
                                       figure.name, table.name){
  shifts = c("exact","10dayslate", "20dayslate", "30dayslate")

  all_powers = data.frame(matrix(NA,nrow=0,ncol=3), check.names = FALSE)
  names(all_powers) = c("Test case", "peaks measured","power")

  all_powers
  quantiles = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
  all_quantiles = data.frame(matrix(NA,nrow=0,ncol=3+length(quantiles)), check.names = FALSE)
  names(all_quantiles) = c("Test case", "time", "peaks measured", paste0(100*quantiles,"%"))

  dynTrue = "exp"
  tTrue = 0
  nstudies1 = 1000

  for (shift in shifts){
    name.df = paste0(folder.data.shifted.peaks,"all_res_TC_",id,"_ABdynTrue",dynTrue,tTrue,"_shifted",shift,"_SurveySize",survey.size,"_Nstudies",nstudies1,".rds")
    if(file.exists(name.df)){
      cum.inc.per.survey = 100*apply(readRDS(name.df), MARGIN  = 2, FUN = cumsum)
      true.cum.inc = cumsum(eval(parse(text = paste0("test.case.",id,".a"))))
      true.test.case = eval(parse(text = paste0("test.case.",id,".a")))

      power1 = kadelka23cutoffFree::calculate_power(fitted_cumulative_incidences = cum.inc.per.survey/100,
                                                    test.case = true.test.case)
      quantiles1 = apply(cum.inc.per.survey, MARGIN = 1, FUN = quantile, probs = quantiles)
      all_powers = rbind(all_powers, data.frame(`Test case` = id, `peaks measured` = shift,
                                                power = round(100*power1[[1]],2), check.names = FALSE))
      all_quantiles = rbind(all_quantiles, data.frame(`Test case` = id, time = 1:nrow(t(quantiles1)), `peaks measured` = shift,
                                                      t(quantiles1), check.names = FALSE))


    }
  }

  all_powers
  all_quantiles %>% dplyr::select(c("Test case", "time", "peaks measured", "2.5%", "50%", "97.5%"))

  test.case.use = eval(parse(text = paste0("test.case.",id,".a")))
  p1 = ggplot() +
    geom_line(data = all_quantiles, aes(x = time, y = `50%`)) +
    geom_ribbon(data = all_quantiles, aes(x = time, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
    geom_point(data = all_quantiles%>%mutate(true_cum_inc = rep(100*cumsum(test.case.use), length(time)/length(test.case.use))),
               aes(x = time, y =true_cum_inc), pch = 4, size = 2) +
    facet_wrap(~factor(`peaks measured`, levels = c("exact","10dayslate","20dayslate","30dayslate"),
                       labels =  c("exact","10 days late","20 days late","30 days late"))) +
    xlab("Survey") +
    ylab("Cumulative incidence (in %)") +
    theme_bw() +
    theme(text = element_text(size = 16))
  print(p1)

  p2 = ggplot() +
    geom_line(data = all_quantiles, aes(x = time, y = `50%`, col = `peaks measured`)) +
    geom_ribbon(data = all_quantiles, aes(x = time, ymin = `2.5%`, ymax = `97.5%`, fill = `peaks measured`), alpha = 0.2) +
    geom_point(data = all_quantiles%>%mutate(true_cum_inc = rep(100*cumsum(test.case.use), length(time)/length(test.case.use))),
               aes(x = time, y =true_cum_inc), pch = 4, size = 2) +
    xlab("Survey") +
    ylab("Cumulative incidence (in %)") +
    theme_bw() +
    scale_color_discrete(breaks = c("exact","10dayslate","20dayslate","30dayslate"),
                         labels = c("exact","10 days late","20 days late","30 days late")) +
    scale_fill_discrete(breaks = c("exact","10dayslate","20dayslate","30dayslate"),
                        labels = c("exact","10 days late","20 days late","30 days late")) +
    theme(text = element_text(size = 16), legend.position = "top") +
    guides(color = guide_legend(title = "Test case",  nrow = 2, label.hjust = 0.5),
           fill = guide_legend(title = "Test case",  nrow = 2, label.hjust = 0.5))

  print(p2)

  if (save.graph){
    name.pic = paste0(folder_graphs,figure.name,".pdf")
    pdf(file=name.pic, height=9, width=10)
    print(p2)
    dev.off()
  }
  if (save.table){
    print(xtable(all_powers),
          file = paste0(folder_graphs,table.name,".tex"), include.rownames = FALSE)
  }
  return(list(all_powers, all_quantiles))
}


model_misspec_shifted_peaks(id = 1, survey.size = 800, save.graph = TRUE, save.table = TRUE, figure.name = "S10AFigure", table.name = "TableS4_column1")
model_misspec_shifted_peaks(id = 3, survey.size = 2500, save.graph = TRUE, save.table = TRUE, figure.name = "S10BFigure", table.name = "TableS4_column2")
model_misspec_shifted_peaks(id = 6, survey.size = 10000, save.graph = TRUE, save.table = TRUE, figure.name = "S10CFigure", table.name = "TableS4_column3")



########################################
# Figure S11
########################################
###################
# Figure S11A
###################
thr = 0.49
df = read.csv(paste0(folder_data,"datasets/Prete_repeat_blood_donors.csv"))
df$donation_date = seq(as.Date("2020-01-05"), as.Date("2021-04-12"), by = "week")[df$week]
df = df %>% filter(covid_result != "", assay != "COV-2IgGII")
df$result <- as.numeric(sapply(stringr::str_split(df$covid_result, pattern = " "), function(x) x[1]))
df = df %>% group_by(donorid) %>% arrange(donation_date) %>%
  mutate(next_result = lead(result), next_date = lead(donation_date),
         halflife = -log(2)*as.numeric(next_date - donation_date)/(log(next_result) - log(result)))

ids_multiple_donations = df %>% group_by(donorid) %>% count() %>% filter(n>1) %>% .$donorid
ids_multiple_positive_donations = df %>% filter(result > thr) %>% group_by(donorid) %>% count() %>% filter(n>1) %>% .$donorid

df1 = df %>%
  summarise(result_max = max(result), donation_date = donation_date, result = result) %>%
  filter(result >= thr, donation_date <= as.Date("2020-05-31")) %>% group_by(donorid) %>%
  filter(result == result_max) %>%
  mutate(group = as.numeric(donation_date - as.Date("2020-03-01"))/7)

figS11A = ggplot(data = df1 %>% filter(donorid %in% ids_multiple_positive_donations)) +
  geom_boxplot(aes(x = donation_date, y = result_max, group = as.character(group))) +
  theme(legend.position = "none") +
  ylab("Maximal AB level") + xlab("Donation time (2020)") +
  theme(text = element_text(size =16))


name.pic = paste0(folder_graphs,"S11AFigure.pdf")
pdf(file=name.pic, height=5, width=5)
print(figS11A)
dev.off()


###################
# Figure S11B
###################
df.comb = rbind(data.frame(value = kadelka23cutoffFree::Buss_validation_data$peaks,
                           data_set = "convalescents (PCR-confirmed)",
                           var = "peaks"),
                data.frame(value = kadelka23cutoffFree::Buss_validation_data$background_levels,
                           data_set = "pre-pandemic controls",
                           var = "backgrounds"),
                data.frame(value = kadelka23cutoffFree::Buss_validation_data$decay_rates,
                           data_set = "convalescents (PCR-confirmed)",
                           var = "decays"),
                data.frame(value = kadelka23cutoffFree::Prete_validation_data$peaks,
                           data_set = "AB positive repeat AB donors",
                           var = "peaks"),
                data.frame(value = kadelka23cutoffFree::Prete_validation_data$background_levels,
                           data_set = "pre-pandemic controls",
                           var = "backgrounds"),
                data.frame(value = kadelka23cutoffFree::Prete_validation_data$decay_rates,
                           data_set = "AB positive repeat AB donors",
                           var = "decays"))

figS11B = ggplot(data = df.comb %>% filter(var == "peaks")) +
  geom_histogram(aes(x = value, y = ..density.., fill = data_set), alpha = 0.2, position = "identity")+
  geom_density(aes(x = value, col = data_set))+
  geom_boxplot(aes(x = value, y = 0, col = data_set), width = 0.2, outlier.shape = NA, show.legend = FALSE) +
  ylab("Density") +
  xlab(expr(log[10]("peak antibody measurement")))+
  theme(text = element_text(size =16),
        legend.position=c(.4,.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) #+

name.pic = paste0(folder_graphs,"S11BFigure.pdf")
pdf(file=name.pic, height=5, width=5)
print(figS11B)
dev.off()


###################
# Figure S11C
###################
figS11C = ggplot(data = df.comb %>% filter(var == "decays" & value > 0 )) +
  geom_histogram(aes(x = value, y = ..density.., fill = data_set), alpha = 0.2, position = "identity")+
  geom_density(aes(x = value, col = data_set))+
  geom_boxplot(aes(x = value, y = 0, col = data_set), width = 0.2, outlier.shape = NA) +
  ylab("Density") +
  xlab("Base-10 decay rate (1/day)")+
  theme(text = element_text(size =16),
        legend.position="none")

name.pic = paste0(folder_graphs,"S11CFigure.pdf")
pdf(file=name.pic, height=5, width=5)
print(figS11C)
dev.off()


###################
# Figure S11D
###################
figS11D = ggplot(data = df.comb %>% filter(var == "decays" & value > 0)) +
  geom_histogram(aes(x = log10(2)/value, y = ..density.., fill = data_set), binwidth = 5, alpha = 0.2, position = "identity")+
  geom_density(aes(x = log10(2)/value, col = data_set), n = 2^15)+
  geom_boxplot(aes(x = log10(2)/value, y = 0, col = data_set), width = 0.001, outlier.shape = NA) +
  coord_cartesian(xlim = c(0, 365)) +
  ylab("Density") +
  xlab("Halflife (in days)")+
  theme(text = element_text(size =16),
        legend.position="none")

name.pic = paste0(folder_graphs,"S11DFigure.pdf")
pdf(file=name.pic, height=5, width=5)
print(figS11D)
dev.off()



###################
# Figure S11E
###################
figS11E = ggplot(data = df.comb %>% filter(var == "backgrounds")) +
  geom_histogram(aes(x = value, y = ..density..), alpha = 0.2, position = "identity")+
  geom_density(aes(x = value))+
  geom_boxplot(aes(x = value, y = 0), width = 0.2, outlier.shape = NA) +
  ylab("Density") +
  xlab(expr(log[10]("peak antibody measurement")))+
  theme(text = element_text(size =16),
        legend.position="none")

name.pic = paste0(folder_graphs,"S11EFigure.pdf")
pdf(file=name.pic, height=5, width=5)
print(figS11E)
dev.off()



########################################
# Figure S12
########################################
density.pos.AB.vals = function(t1, #time after infection (in days)
                               T1, #timetoproduction (in days)
                               case.data,
                               control.data,
                               grads,
                               AB.growth.option){

  grads.fun = approxfun(density(grads))
  pdf.grads.smooth <- function(x){
    pdf <- grads.fun(x)
    pdf[is.na(pdf)]<-min(pdf, na.rm = T)
    return(pdf)
  }

  case.fun <- approxfun(density(case.data))
  pdf.case.smooth <- function(x){
    pdf <- case.fun(x)
    pdf[is.na(pdf)]<-min(pdf, na.rm = T)
    return(pdf)
  }

  control.fun <- approxfun(density(control.data))
  pdf.control.smooth <- function(x){
    pdf <- control.fun(x)
    pdf[is.na(pdf)]<-1e-5
    return(pdf)
  }

  N_total = 100000
  ##----sample N_total peak values
  xmax = max(density(case.data)$x)
  xmin = min(density(case.data)$x)
  dg <- function(x) dunif(x, min = xmin, max = xmax)
  rg <- function(n) runif(n, min = xmin, max = xmax)
  df <- function(x) pdf.case.smooth(x)
  # determine M as max(pdf.case.smooth(seq(xmin,xmax,1e-5)))*(xmax-xmin)
  dat = rejectionSampling(N_total, df=df, dg=dg, rg=rg,M =6)

  ##---sample N_total minimal values
  xmax.control = max(density(control.data)$x)
  xmin.control = min(density(control.data)$x)
  dg.control <- function(x) dunif(x, min = xmin.control, max = xmax.control)
  rg.control <- function(n) runif(n, min = xmin.control, max = xmax.control)
  df.control <- function(x) pdf.control.smooth(x)
  d.contr = rejectionSampling(N_total, df=df.control, dg=dg.control, rg=rg.control,M=4)

  # ##---sample N_total decay rate
  xmax = max(density(grads)$x)
  xmin = min(density(grads)$x)
  dg <- function(x) dunif(x, min = xmin, max = xmax)
  rg <- function(n) runif(n, min = xmin, max = xmax)
  df <- function(x) pdf.grads.smooth(x)
  # determine M as max(pdf.grads.smooth(seq(xmin,xmax,1e-5)))*(xmax-xmin)
  rates = rejectionSampling(N_total, df=df, dg=dg, rg=rg,M =5)
  rates[rates < 0] = min(rates[rates>0])

  ## peak times
  tpeaks = rgamma(N_total , shape = (30)/2, rate = 0.5)
  tpeaks[tpeaks<T1] = T1

  t1.add.days = runif(N_total, min = 0, max = time.between.surveys-1)
  AB_vals = kadelka23cutoffFree::antibody_dynamics(AB.background = d.contr, AB.peak = dat, t.growth = T1, t.peak = tpeaks, decay.rate = rates, t = t1 + t1.add.days, growth.dynamics = AB.growth.option)

  return(AB_vals)
}
density.pos.AB.vals.vec = Vectorize(density.pos.AB.vals, vectorize.args = "t1")

for(validation.data.old in c(FALSE, TRUE)){
  if(validation.data.old == TRUE){
    control.data = kadelka23cutoffFree::Buss_validation_data$background_levels
    case.data = kadelka23cutoffFree::Buss_validation_data$peaks
    grads = kadelka23cutoffFree::Buss_validation_data$decay_rates
    figure.name1 = "S12AFigure.pdf"
    figure.name2 = "S12BFigure.pdf"
  }
  if(validation.data.old == FALSE){
    control.data = kadelka23cutoffFree::Prete_validation_data$background_levels
    case.data = kadelka23cutoffFree::Prete_validation_data$peaks
    grads = kadelka23cutoffFree::Prete_validation_data$decay_rates
    figure.name1 = "S12CFigure.pdf"
    figure.name2 = "S12DFigure.pdf"
  }

  A1.exp0 = density.pos.AB.vals.vec(t1 = (0:(8-1))*time.between.surveys,
                                    T1 = 0,
                                    case.data = case.data,
                                    control.data = control.data,
                                    grads = grads,
                                    AB.growth.option = "exp")
  A1.exp7 = density.pos.AB.vals.vec(t1 = (0:(8-1))*time.between.surveys,
                                    T1 = 7,
                                    case.data = case.data,
                                    control.data = control.data,
                                    grads = grads,
                                    AB.growth.option = "exp")
  A1.linear0 = density.pos.AB.vals.vec(t1 = (0:(8-1))*time.between.surveys,
                                       T1 = 0,
                                       case.data = case.data,
                                       control.data = control.data,
                                       grads = grads,
                                       AB.growth.option = "linear")
  A1.linear7 = density.pos.AB.vals.vec(t1 = (0:(8-1))*time.between.surveys,
                                       T1 = 7,
                                       case.data = case.data,
                                       control.data = control.data,
                                       grads = grads,
                                       AB.growth.option = "linear")
  B = data.frame(exp0 = A1.exp0[,1],
                 exp7 = A1.exp7[,1],
                 linear0 = A1.linear0[,1],
                 linear7 = A1.linear7[,1]) %>% reshape::melt()


  pic.density.first.month.with.background.and.peak = ggplot() +
    geom_density(aes(x = case.data, y = ..density..), fill = rgb(0,0,0,0.2), col = rgb(0,0,0,0), size = 2) +
    geom_density(aes(x = control.data, y = ..density..), fill = rgb(0,0,0,0.4), col = rgb(0,0,0,0), size = 2) +
    geom_density(data = B, aes(x = value, y = ..density.., col = variable), size = 1) +
    xlab("log10(abbott_sc)") +
    scale_color_discrete(name = "growth dynamics") +
    theme(legend.position = "top",
          text = element_text(size = 16)) +
    xlab(expr(log[10]("antibody level"))) +
    xlab("log(AB level)") +
    ylab("Density") +
    scale_y_continuous( expand = expansion(mult=c(0,0.05),add = 0)) +
    scale_x_continuous(limits = c(-2.7, 1.5), expand = c(0,0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 16),
          legend.position = "top") +
    guides(color = guide_legend(title.position = "left",nrow = 2))
  pic.density.first.month.with.background.and.peak

  name.pic = paste0(folder_graphs,figure.name1)
  pdf(file=name.pic, height=5, width=5)
  print(pic.density.first.month.with.background.and.peak)
  dev.off()


  pic.density.months.2to8 = ggplot(data = A1.exp0 %>%
                                     data.frame(.,check.names = FALSE) %>%
                                     setNames(.,paste0((0:(8-1)),"-",(1:8))) %>%
                                     reshape::melt() %>%
                                     dplyr::rename(months_past_infection = variable) %>%
                                     filter(months_past_infection!="0-1")) +
    geom_density(aes(x = value, col = months_past_infection), size = 1) +
    scale_color_discrete(name = "Months past\ninfection")  +
    scale_fill_discrete(name = "Months past\ninfection") +
    xlab(expr(log[10]("antibody level"))) +
    ylab("Density") +
    theme(text = element_text(size = 16),
          legend.position = "top")
  pic.density.months.2to8

  name.pic = paste0(folder_graphs,figure.name2)
  pdf(file=name.pic, height=5, width=5)
  print(pic.density.months.2to8)
  dev.off()
}


########################################
# Figure S13
########################################
D.plot = output.summary.all %>% filter(maxmonth == maxmonth.plot, `growth dynamics` == "exp0")
D.plot = D.plot %>% filter(time <= as.Date("2020-10-31")) %>% mutate(city = citymap[city])
### for monthly estimation use mid (15th) of month
h1 = D.plot$time[D.plot$repeated == "monthly"]
D.plot$time[D.plot$repeated == "monthly"] = paste0(year(h1),"-",month(h1),"-",15)
D.plot$oldValidationData[D.plot$oldValidationData == TRUE] = "(a) convalescents, \nPCR positive"
D.plot$oldValidationData[D.plot$oldValidationData == FALSE] = "(b) AB positive \nrepeat blood donors"

figS13 = ggplot() +
  labs(x = "Date in 2020", y = "Mortality per million inhabitants,\nCumulative incidence (in %)") +
  geom_point(data = D.plot, aes(x = time, y = 100*`point estimate`, col = `repeated`)) +
  geom_line(data = D.plot, aes(x = time, y = 100*`50%`, col = `repeated`))+
  geom_ribbon(data = D.plot, aes(x = time, ymin = 100*`2.5%`, ymax = 100*`97.5%`, fill = `repeated`), alpha = 0.2) +
  theme(
    legend.position = "right",
    axis.line = element_line(colour = "black"),
    text = element_text(size = 16)) +
  facet_grid(city~oldValidationData, scales = "free_y") +
  scale_color_discrete(name = "estimated", guide = guide_legend(title.position = "top", nrow = 2)) +
  scale_fill_discrete(name = "estimated", guide = guide_legend(title.position = "top", nrow = 2)) #+

name.pic = paste0(folder_graphs,"S13Figure.pdf")
pdf(file=name.pic, height=13, width=7)
print(figS13)
dev.off()


########################################
# Figure S14
########################################
plot.city = "HEMOAM"
cities = plot.city
D.mortality = mortality_data %>%
  filter(DT_SIN_PRI <= as.Date("2020-11-30"))%>%
  filter(city %in% cities_plot) %>%
  mutate(location = citymap.rev[city]) %>%
  filter(location %in% cities)

D.plot = output.summary.all %>% filter(maxmonth == maxmonth.plot, city == plot.city)
D.plot = D.plot %>% filter(time <= as.Date("2020-10-31")) %>% mutate(city = citymap[city])
### for monthly estimation use mid (15th) of month
h1 = D.plot$time[D.plot$repeated == "monthly"]
D.plot$time[D.plot$repeated == "monthly"] = paste0(year(h1),"-",month(h1),"-",15)
D.plot$oldValidationData[D.plot$oldValidationData == TRUE] = "(a) convalescents, \nPCR positive"
D.plot$oldValidationData[D.plot$oldValidationData == FALSE] = "(b) AB positive \nrepeat blood donors"

figS14 = ggplot() +
  labs(x = "Date in 2020", y = "Mortality per million inhabitants,\nCumulative incidence (in %)") +
  geom_ribbon(data = D.mortality, aes(x = DT_SIN_PRI, ymin = 0, ymax= 1E6*MR_mean, group = city), alpha = 0.2, col = rgb(0,0,0,0.3)) +
  geom_point(data = D.plot, aes(x = time, y = 100*`point estimate`, col = `growth dynamics`)) +
  geom_line(data = D.plot, aes(x = time, y = 100*`50%`, col = `growth dynamics`))+
  geom_ribbon(data = D.plot, aes(x = time, ymin = 100*`2.5%`, ymax = 100*`97.5%`, fill = `growth dynamics`), alpha = 0.2) +
  theme(
    legend.position = "top",
    axis.line = element_line(colour = "black"),
    text = element_text(size = 16)) +
  ylim(as.numeric(unlist(y.limits.map[plot.city]))) +
  facet_grid(paste0("estimated ",repeated)~oldValidationData) +
  scale_color_discrete(name = "growth dynanmics", guide = guide_legend(title.position = "left", nrow = 2)) +
  scale_fill_discrete(name = "growth dynanmics", guide = guide_legend(title.position = "left", nrow = 2)) #+

name.pic = paste0(folder_graphs,"S14Figure.pdf")
pdf(file=name.pic, height=8, width=8)
print(figS14)
dev.off()


########################################
# Figure S15
########################################
############################### AGE AND SEX STRATIFIED
cities = c("HEMEPAR","FPS","HEMOPE", "HEMOAM", "HEMOBA", "HEMOCE", "HEMORIO", "HEMOMINAS")
agesexgroups.names = c("_weighted", "_M16", "_M25", "_M35", "_M45", "_M55", "_F16", "_F25", "_F35", "_F45", "_F55")
a = c("cum_inc", "q1", "median", "q2")
colnames1= c("month", paste0(rep(a, n = length(agesexgroups.names)), rep(agesexgroups.names, each = length(a))),
             "growth dynamics","time","repeated","oldValidationData","maxmonth","city")
output.summary.all.agesex = data.frame(matrix(nrow = 0, ncol = length(colnames1)))
names(output.summary.all.agesex) = colnames1
for(city in cities){
  for (maxmonth in c(11)){
    survey = (blood_donors%>%filter(blood_center == city, month %in% 3:maxmonth))
    min.week = min(survey$week)
    for (monthly_weekly in c('weekly','monthly')){
      for (oldValidationData in c("FALSE","TRUE")){
        for (ABgrowth1 in c('exp','linear')){
          for (ttoprod1 in c(0,7)){
            filename1 =  paste0(folder_data,"datasets/all_results_withage_",monthly_weekly,"_oldValidationData_",
                                oldValidationData,"/",city,"_ABdynTrue",ABgrowth1,ttoprod1,
                                "_Nstudies1000_NbruteForce0_maxmonth",maxmonth,"_age_sex_weighted.rds")
            if(file.exists(filename1)){
              output1 = readRDS(file = filename1)
              if(monthly_weekly == "weekly"){
                approx.date.surveys  = as.Date(paste(2020, min.week:(min.week+dim(output1)[1]-1), 1, sep="-"), "%Y-%U-%u")
              }
              if(monthly_weekly == "monthly"){
                approx.date.surveys = seq(as.Date("2020-04-01"),length=maxmonth-2,by="months")-1### last day of month
              }
              output1.summary = output1
              output1.summary$`growth dynamics` = paste0(ABgrowth1,ttoprod1)
              output1.summary$time = approx.date.surveys
              output1.summary$repeated = monthly_weekly
              output1.summary$oldValidationData = oldValidationData
              output1.summary$maxmonth = maxmonth
              output1.summary$city = city
              output.summary.all.agesex = rbind(output.summary.all.agesex, output1.summary)
            }
          }
        }
      }
    }
  }
}

plot.city = "HEMOAM"
maxmonth.plot = 11
repeated.plot = "monthly" # or "weekly"
oldValidationData.plot = FALSE
data.for.plot = output.summary.all.agesex%>%filter(city == plot.city, repeated == repeated.plot, oldValidationData == oldValidationData.plot, maxmonth == maxmonth.plot)
D.plot = data.frame(date = data.for.plot$time,
                    `point estimate` = data.for.plot$cum_inc_weighted,
                    `2.5%` = data.for.plot$q1_weighted,
                    `50%` = data.for.plot$median_weighted,
                    `97.5%` = data.for.plot$q2_weighted,
                    `growth dynamics` = data.for.plot$`growth dynamics`,
                    `repeated` = data.for.plot$repeated,
                    oldValidationData = data.for.plot$oldValidationData,
                    maxmonth = data.for.plot$maxmonth,
                    city = data.for.plot$city,
                    age_sex = "weighted",
                    check.names = FALSE)


ages = c(16,25,35,45,55)
groups = c(paste0("M",ages),paste0("F",ages))
for(agesex in groups){
  D.plot.add = data.frame(date = data.for.plot$time,
                          `point estimate` = data.for.plot[,paste0("cum_inc_",agesex)],
                          `2.5%` = data.for.plot[,paste0("q1_",agesex)],
                          `50%` = data.for.plot[,paste0("median_",agesex)],
                          `97.5%` = data.for.plot[,paste0("q2_",agesex)],
                          `growth dynamics` = data.for.plot$`growth dynamics`,
                          `repeated` = data.for.plot$repeated,
                          oldValidationData = data.for.plot$oldValidationData,
                          maxmonth = data.for.plot$maxmonth,
                          city = data.for.plot$city,
                          age_sex = agesex,
                          check.names = FALSE)
  D.plot = rbind(D.plot, D.plot.add)
}
h1 = D.plot$date[D.plot$repeated == "monthly"]
D.plot$date[D.plot$repeated == "monthly"] = paste0(year(h1),"-",month(h1),"-",15)

figS15 = ggplot(data = D.plot %>% filter(month(date) < 11)) +
  geom_point(aes(x = date, y = `point estimate`)) +
  geom_ribbon(aes(x = date, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  geom_line(aes(x = date, y = `50%`)) +
  facet_wrap(~age_sex) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    text = element_text(size = 20)) +
  xlab("Date (in 2020)") +
  ylab("Cumulative incidence")

name.pic = paste0(folder_graphs,"S15Figure.pdf")
pdf(file=name.pic, height=10, width=10)
print(figS15)
dev.off()


########################################
# Figure S16
########################################
#### compare unweighted results with those weighted by age and sex
D1 = output.summary.all %>% mutate(stratified = "no")
D2 = output.summary.all.agesex %>% dplyr::select(c("growth dynamics", "time", "repeated", "oldValidationData",  "maxmonth", "city", "cum_inc_weighted", "q1_weighted", "median_weighted", "q2_weighted"))
names(D2)[names(D2) %in% c("cum_inc_weighted", "q1_weighted", "median_weighted", "q2_weighted")] = c("point estimate", "2.5%", "50%", "97.5%")
D2 = D2 %>%  mutate(stratified = "by age and sex")
D2 = D2 %>% dplyr::select(names(D1))
D.combined = rbind(D1, D2)

plot.city = c("HEMEPAR","FPS","HEMOPE", "HEMOAM", "HEMOBA", "HEMOCE", "HEMORIO", "HEMOMINAS")
maxmonth.plot = 11
repeated.plot = "monthly"
oldValidationData.plot = FALSE
dynamics.plot = "exp0"

D.plot = D.combined %>% filter(maxmonth == maxmonth.plot,
                               city %in% plot.city,
                               repeated == repeated.plot,
                               oldValidationData == oldValidationData.plot,
                               `growth dynamics` == "exp0")
h1 = D.plot$time[D.plot$repeated == "monthly"]
D.plot$time[D.plot$repeated == "monthly"] = paste0(year(h1),"-",month(h1),"-",15)

figS16 = ggplot(data = D.plot %>% filter(month(time) < 11)) +
  geom_point(aes(x = time, y = `point estimate`, col = stratified)) +
  geom_ribbon(aes(x = time, ymin = `2.5%`, ymax = `97.5%`, fill = stratified), alpha = 0.4) +
  geom_line(aes(x = time, y = `50%`, col = stratified)) +
  facet_wrap(~citymap[city], scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.line = element_line(colour = "black"),
    text = element_text(size = 20)) +
  xlab("Date (in 2020)") +
  ylab("Cumulative incidence") +
  guides(color = guide_legend(title = "weighted"),
         fill = guide_legend(title = "weighted"))

name.pic = paste0(folder_graphs,"S16Figure.pdf")
pdf(file=name.pic, height=10, width=10)
print(figS16)
dev.off()




#########################
# Tables S5 - S12
#########################
cities = c("HEMEPAR","FPS","HEMOPE", "HEMOAM", "HEMOBA", "HEMOCE", "HEMORIO", "HEMOMINAS")
table.map = c("HEMEPAR" = "11","FPS" = "12","HEMOPE" = "9", "HEMOAM" = "5",
              "HEMOBA" = "7", "HEMOCE" = "6", "HEMORIO" = "8", "HEMOMINAS" = "10")

for (city.use in cities){
  data.use.week = output.summary.all%>% filter(maxmonth == 11, repeated == "weekly") %>%
    mutate(`2.5%` = round(100*`2.5%`,2), `97.5%` = round(100*`97.5%`,2),
           `50%` = round(100*`50%`,2),
           `point estimate` = round(100*`point estimate`,2)) %>%
    group_by(city, month(time)) %>%
    filter(day(time) == max(day(time))) %>%
    data.frame(., check.names = FALSE) %>%
    mutate(`validation` = ifelse(oldValidationData == FALSE, "(b)", "(a)"),
           `estimate` = paste0("placeholder1",`point estimate`,"placeholder2","(", `2.5%` ,",", `97.5%`,")","placeholder3")) %>%
    filter(city == city.use) %>% mutate(`validation data` = validation, estimated = repeated)


  A.week = data.use.week %>%
    filter(`month(time)` == 3) %>%
    dplyr::select(c("city","validation data", "estimated","growth dynamics")) %>%
    mutate(March = NA, April = NA, May = NA, June = NA, July = NA, August = NA, September = NA, October = NA)

  for(i in 3:10){
    A.week[,i+2] = data.use.week %>% filter(`month(time)` == i) %>% pull(estimate)
  }

  data.use.month = output.summary.all%>%
    filter(maxmonth == 11, repeated == "monthly") %>%
    mutate(`2.5%` = round(100*`2.5%`,2),
           `97.5%` = round(100*`97.5%`,2),
           `50%` = round(100*`50%`,2),
           `point estimate` = round(100*`point estimate`,2)) %>%
    group_by(city, month(time)) %>%
    filter(day(time) == max(day(time))) %>%
    data.frame(., check.names = FALSE) %>%
    mutate(`validation` = ifelse(oldValidationData == FALSE, "(b)", "(a)"),
           `estimate` = paste0("placeholder1",`point estimate`,"placeholder2","(", `2.5%` ,",", `97.5%`,")","placeholder3")) %>%
    filter(city == city.use) %>% mutate(`validation data` = validation, estimated = repeated)

  A.month = data.use.month %>%
    filter(`month(time)` == 3) %>%
    dplyr::select(c("city","validation data", "estimated","growth dynamics")) %>%
    mutate(March = NA, April = NA, May = NA, June = NA, July = NA, August = NA, September = NA, October = NA)
  for(i in 3:10){
    A.month[,i+2] = data.use.month %>% filter(`month(time)` == i) %>% pull(estimate)
  }

  A = rbind(A.week, A.month)
  A$city = citymap[A$city]
  A$city[2:16] = ""
  A$`validation data`[c(2:4, 6:8, 10:12, 14:16)] = ""
  A$estimated[c(2:4, 6:8, 10:12, 14:16)] = ""
  print(xtable(A, caption = paste0("Estimated cumulative incidences in ",gsub(" ","",citymap[city.use]),"."),
               label = paste0("tab:",gsub(" ","",citymap[city.use]),"_cum_inc")),
        include.rownames = FALSE,
        size = "\\footnotesize",
        file = paste0(folder_graphs,"TableS",as.character(table.map[city.use]),".tex"))

}




