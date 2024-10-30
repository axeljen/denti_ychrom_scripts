rm(list = ls())
setwd("/Users/axeljensen/Dropbox/notes/Attachments/denti-lesula/bpp")
library(tidyverse)
library(bppr)
library(ape)
library(coda)
library(ggpubr)
library(psych)
library(MoreTreeTools)
library(treeio)
library(ggrepel)

# specify a basename that will be used for saving plots
NAME <- "output_name"

# will make a directore in current workdir
dir.create(NAME)

# I ran two independent runs for each model, read those two in
mcmc1 <- read_tsv("model_XX_run1_output.mcmc")
mcmc2 <- read_tsv("model_XX_run2_output.mcmc")

# manually specify the order of the tips, from bottom to top
tiporder <- c("C.denti","C.wolfi","C.pogonias","C.mona","C.neglectus","C.mitis","C.nictitans","C.cephus")

### Specify hybrid nodes as they're specified in the control file
hybrid_sources <- c("A","P","Z")
hybrid_dests <- c("B","Q","W")


# fetch topology from the outfile
outfile <- read_file("model_XX_run1_output.txt")
outfile <- str_split(outfile, "\n")

# reduce size as the tree will be towards the end
start <- length(outfile[[1]])-100
end <- length(outfile[[1]])
outfile <- outfile[[1]][start:end]
# get the species tree
next_is_tree <- FALSE
for (line in outfile){
  if (next_is_tree){
    treestring <- line
    top <- read.tree(text = treestring)
    next_is_tree = FALSE
  }
  if (line == "Species tree network:"){
    next_is_tree = TRUE
  }
}

# then we can run the full thing
{
print("Plotting script started.")
# rename the tiplabs of hybrid nodes to source
for (i in 1:length(top$tip.label)) {
  if (top$tip.label[i] %in% c(hybrid_sources, hybrid_dests)) {
    top$tip.label[i] <- paste0(top$tip.label[i], ".source")
  }
}

# calculate the mean values for both and check convergence
## start wotj rim 2
print("Summarizing mcmc samples.")
run1_sum <- mcmc1 %>%
  mutate(run = "run1") %>%
  # pivot all but first column, wich is just the generation number
  pivot_longer(cols = 2:ncol(mcmc1)) %>%
  # group by the parameter, which is in "name" col
  group_by(name) %>%
  # calculate means for each estimate
  summarize(mean_run1 = mean(value))

# summarize second run and join it with first
sum <- mcmc2 %>%
  mutate(run = "run2") %>%
  # pivot all but first column, wich is just the generation number
  pivot_longer(cols = 2:ncol(mcmc2)) %>%
  # group by the parameter, which is in "name" col
  group_by(name) %>%
  # calculate means for each estimate
  summarize(mean_run2 = mean(value)) %>%
  # bind this together with the first run
  left_join(run1_sum,.,by = "name")

## Compare estimates from the two runs

# plot times and Nes separately against each other for more detail
tauplot <- ggplot(sum[grepl("tau_", sum$name),], aes(x = mean_run1, y = mean_run2)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 2.5) +
  theme_bw()

# Nes
thetaplot <- ggplot(sum[grepl("theta_", sum$name),], aes(x = mean_run1, y = mean_run2)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 2.5) +
  theme_bw()


# and phi
phiplot <- ggplot(sum[grepl("phi_", sum$name),], aes(x = mean_run1, y = mean_run2)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 2.5) +
  theme_bw()
phiplot

print("Making the convergence plots.")
# make an arranged plot from these and save
convergence <- ggarrange(tauplot,thetaplot,phiplot, labels = c("Ne estimates","Divergence time estimates","migration proportion"))

ggsave(convergence, filename = paste0(NAME, "/", NAME, "_convergence.jpg"))
ggsave(convergence, filename = paste0(NAME, "/", NAME, "_convergence.jpg"))

############# Scale the mcmc traces to years and Nes

# join the two mcmcs to calculate the stats
mcmc.joint <- rbind(mcmc1, mcmc2)

# mutation rate to use for scaling
mu <- 4.82e-9

# generation time
g <- 10 # here I'm using a generation time of 10 years, which most guenons are spread around in Kuderna et al.

print(paste0("Scaling the mcmc output ugsing mutation rate ", mu, " and generation time of ", g, " years."))

# scale time with the bppr package, setting sds close to 0
times.joint <- msc2time.r(mcmc.joint, u.m = mu, u.sd = 0.000001e-8,
                    g.mean = g, g.sd = 0.000001)
# get the means
means.joint <- as.data.frame(apply(times.joint, 2, mean)) %>%
  mutate(stat = rownames(.))
# calculate 95 % highest posterior probability intervals
conf.joint<- as.data.frame(coda::HPDinterval(coda::as.mcmc(times.joint))) %>%
  mutate(stat = rownames(.))
# join this with the means
joint <- means.joint %>%
  left_join(., conf.joint, by = "stat")

# fetch the scaling factor to get the correct branch lengths in generations/years

# fetch the times for internal nodes
tau.originals <- mcmc.joint %>%
  select(starts_with("tau_")) %>%
  pivot_longer(cols = 1:ncol(.)) %>%
  group_by(name) %>%
  summarise(mean_t = mean(value)) %>%
  mutate(nodename = gsub("tau_","t_",name))

# scale to years
scaled.times <- times.joint %>%
  select(starts_with("t_")) %>%
  pivot_longer(cols = 1:ncol(.)) %>%
  group_by(name) %>%
  summarize(mean_time_years = mean(value),
            confint_upper_years = coda::HPDinterval(coda::as.mcmc(value))[1],
            confint_lower_years = coda::HPDinterval(coda::as.mcmc(value))[2])
colnames(scaled.times) <- c("nodename", colnames(scaled.times)[2:length(colnames(scaled.times))])

# add terminal branches for these
for (node in top$tip.label) {
  dftau <- data.frame(name = paste0("tau_", node),
                      mean_t = 0, nodename = node)
  tau.originals <- rbind(tau.originals, dftau)
  dftimes <- data.frame(nodename = node,
                        mean_time_years = 0, 
                        confint_upper_years = 0,
                        confint_lower_years= 0)
  scaled.times <- rbind(scaled.times,
                        dftimes)
}

# remove the "t_xx" prefix
scaled.times <- scaled.times %>%
  mutate(nodename = gsub("t_\\d{2}", "", nodename))
tau.originals <- tau.originals %>%
  mutate(nodename = gsub("t_\\d{2}", "", nodename))

# bind these together
times.translate <- tau.originals %>%
  left_join(.,scaled.times, by = c("nodename")) %>%
  mutate(scale_factor = (mean_time_years / mean_t))
  # add scaling factor in 1,000,000 units

# add info on type of node
times.translate$nodetype <- NA
for (tip in top$tip.label){
  times.translate$nodetype[times.translate$nodename == tip] <- "tip"
}
for (node in top$node.label){
  times.translate$nodetype[times.translate$nodename == node] <- "node"
}
for (hybsource in hybrid_sources) {
  times.translate$nodetype[times.translate$nodename == hybsource] <- "hybrid_source"
}
for (hybdest in hybrid_dests) {
  times.translate$nodetype[times.translate$nodename == hybdest] <- "hybrid_dest"
}
print("Scaling Nes.")
# scale factor for Nes
theta.originals <- mcmc.joint %>%
  select(starts_with("theta_")) %>%
  pivot_longer(cols = 1:ncol(.)) %>%
  group_by(name) %>%
  summarize(mean_theta = mean(value)) %>%
  mutate(nodename = gsub("theta_\\d*", "",name))

scaled.theta<- times.joint %>%
  select(starts_with("Ne_")) %>%
  pivot_longer(cols = 1:ncol(.)) %>%
  group_by(name) %>%
  summarize(mean_ne = mean(value)) %>%
  mutate(nodename = gsub("Ne_\\d*", "", name))

nes.translate <- theta.originals %>%
  left_join(., scaled.theta, by = "nodename") %>%
  mutate(scale_factor_ne = mean_ne / mean_theta)

print("Creating the master data frame with all branches, ages and Nes.")

# now let's make a "master" table with all necessary info
d <- left_join(times.translate, nes.translate,
               by = "nodename")

########### Break out the Nes for all branches to increase readability

print("Averaging/scaling the Nes of branches affected by hybridization, to only have one Ne per branch.")

# take harmonic mean Ne of branches containing hybrid node
hybrid_branches <- data.frame(from = character(), through_1 = character(),
                              through_2 = character(),
                              to = character())
for (hyb in c(hybrid_sources, hybrid_dests)){
  integrated <- NA
  # get the parent of the node
  parent <- nodelab(top,getParent(top, nodeid(top,hyb)))
  # and the child
  child <- nodelab(top,child(top,hyb))
  # if the parent is a hybrid node, then keep walking till we hit a real node
  realnodefound <- !parent %in% c(hybrid_sources, hybrid_dests)
  while (!realnodefound) {
    print(paste0("Parent of ", hyb, " is ", parent))
    integrated <- hyb
    hyb <- parent
    print(paste0(parent, " is the parent of ", integrated))
    child <- nodelab(top,child(top,integrated))
    parent <- nodelab(top, getParent(top, nodeid(top,parent)))
    # switch their places to have them ordered
    realnodefound <- !parent %in% c(hybrid_sources, hybrid_dests)
  }
  # get the time and ne of parent

  if (length(child) > 1){
    #if there's two children, take the one without source in it
    for (c in child){
      if (!grepl("source",c)){
        child <- c
      }
    }
  }
  # if the child is a hybrid node we'll just skip it alltogether, as it should be included on the parent side
  print(child)
  if (!child %in% c(hybrid_dests,hybrid_sources)){
  # from will be the final parent
  f <- parent
  # through_1 will be the first hybrid node
  t_1 <- hyb
  # through 2 will be the integrated if any
  t_2 <- integrated
  # and last the destination will be the child
  dest <- child
  r <- data.frame(from = f, through_1 = t_1, through_2 = t_2, dest = child)
  print(r)
  hybrid_branches <- rbind(hybrid_branches, r)
  }
}

# now let's adjust all the Nes of the hybrid children to their harmonic mean
d$adjusted_ne <- d$mean_ne
# and we'll also have a column with their "real" parent node, excluding the hybrid nodes
d$parent <- NA
for (tip in top$tip.label){
  d$parent[d$nodename == tip] <- nodelab(top,getParent(top, nodeid(top,tip)))
}
for (node in top$node.label){
  d$parent[d$nodename == node] <- nodelab(top,getParent(top, nodeid(top,node)))
  
}
print("Adjusting the parents and children of hybrid nodes.")
for (i in 1:nrow(hybrid_branches)){
  # so, for each hybrid branch, there's some stuff to check
  destination <- hybrid_branches$dest[i]
  print(destination)
  from <- hybrid_branches$from[i]
  hyb1 <- hybrid_branches$through_1[i]
  print(hyb1)
  if (is.na(hybrid_branches$through_2[i])){
  # grab the first two times and nes, going backwards in time
  t1 <- d$mean_time_years[d$nodename == hyb1] - d$mean_time_years[d$nodename == destination]
  print(t1)
  t2 <- d$mean_time_years[d$nodename == from] - t1
  ne1 <- d$mean_ne[d$nodename == destination]
  ne2 <- d$mean_ne[d$nodename == hyb1]
  ne.adj <- harmonic.mean(c(rep(ne1,t1),rep(ne2,t2)))
  # and then we'll make an harmonic mean ne out of this
  ne.adj <- harmonic.mean(c(rep(ne1,t1),rep(ne2,t2)))
  } else {
  # if we have an internal hybrid node too, we need to incorporate that additional segment
    hyb1 <- hybrid_branches$through_2[i]
    hyb2 <- hybrid_branches$through_1[i]
    t1 <- d$mean_time_years[d$nodename == hyb1] - d$mean_time_years[d$nodename == destination]
    t2 <- d$mean_time_years[d$nodename == hyb2] - t1
    t3 <- d$mean_time_years[d$nodename == from] - t2
    ne1 <- d$mean_ne[d$nodename == destination]
    ne2 <- d$mean_ne[d$nodename == hyb1]
    ne3 <- d$mean_ne[d$nodename == hyb2]
    # harmonic mean of these
    ne.adj <- harmonic.mean(c(rep(ne1,t1),rep(ne2,t2),rep(ne3,t3)))
  }
  d$adjusted_ne[d$nodename == destination] <- ne.adj
  # and set the "corrected" parent too
  d$parent[d$nodename == destination] <- from
}

# and a root where we cut the plot off
root_t <- 8000000
print(paste0("Cutting the tree off at a root of ", root_t, " years ago."))

# drop macaque as we don't want to plot this
print("Removing unwanted branches from the tree.")
top <- drop.tip(top, "M.mulatta")

# and also remove the "true" root from the data table
d <- d %>%
  filter(node != "guenon_rootM.mulatta") %>%
  filter(node != "M.mulatta") %>%
  mutate(parent = if_else(parent == "guenon_rootM.mulatta", "root", parent))

# and anyones containing "source"
falsetips <- top$tip.label[grepl("source",top$tip.label)]
top <- drop.tip(top, falsetips)

# let's prep a df with branches
print("Prepping data frame with branches etc for plotting.")
branches <- data.frame(node = top$tip.label,
                       time = rep(0,length(top$tip.label)),
                       confint.lower = NA,
                       confint.upper = NA,
                       ne = d$adjusted_ne[match(top$tip.label, d$nodename)],
                       parent_node = d$parent[match(top$tip.label, d$nodename)])
# prep internal branches
intbranches <- data.frame(node = top$node.label,
                          time = d$mean_time_years[match(top$node.label, d$nodename)],
                          confint.lower = d$confint_lower_years[match(top$node.label, d$nodename)],
                          confint.upper = d$confint_upper_years[match(top$node.label, d$nodename)],
                          parent_node = d$parent[match(top$node.label, d$nodename)],
                          ne = d$adjusted_ne[match(top$node.label, d$nodename)])
# bind'em together and add an "artificial" root
branches <- rbind(branches, intbranches)

# order them by the order of the tiplabs
print("Ordering tips in the tree.")
#tips <- get_taxa_name(ggtree(top))
nodeorder <- c(tiporder, d$nodename[which(d$nodetype == "node" & d$nodename %in% c(top$node.label, top$tip.label))])
branches <- branches[match(nodeorder, branches$node),]

# add root 
branches <- branches %>%
  bind_rows(data.frame(node = "root",time= root_t, confint.lower = NA,
                       confint.upper = NA, parent_node = NA, ne = NA))

# initate ypos and time_end for ending of branch
branches$ypos <- NA
branches$time_end <- NA

# give all tips a y-pos, dispersed the centers based on total Ne-s to cover and a margin to separate branches with
mar <- 200000
branches$ypos[branches$time == 0] <- cumsum(c(0, rep((sum(branches$ne[which(branches$time == 0)]) / 
                                                       nrow(branches[branches$time == 0,])) + mar,
                                                     nrow(branches[branches$time == 0,]) - 1)))

print("Setting Y-positions of internal branches")
# loop through all the branches and set the ypos of the internal nodes and their endtimes
for (i in c(1,2)){
for (n in 1:nrow(branches)){
  print(branches$node[n])
  max <- nrow(branches)
  if (branches$node[n] != 'root'){
    if (is.na(branches$ypos[which(branches$node == branches$parent_node[n])])){
      for (n2 in (n+1):max){
        if (branches$node[n2] != "root"){
          if (branches$parent_node[n2] == branches$parent_node[n]){
            n_ypos <- branches$ypos[n]
            n2_ypos <- branches$ypos[n2]
            print(n_ypos)
            print(n2_ypos)
            parent_ypos <- min(c(n_ypos,n2_ypos)) + ((max(c(n_ypos,n2_ypos)) - min(c(n_ypos,n2_ypos))) / 2)
            branches$ypos[which(branches$node == branches$parent_node[n])] <- parent_ypos
            # time end of the children branches will be that of the parents
            branches$time_end[n] <- branches$time[which(branches$node == branches$parent_node[n])]
            branches$time_end[n2] <- branches$time[which(branches$node == branches$parent_node[n])]
            break
          }
        }
      }
    }
  }
}
}
# repeating the above step is an ugly workaround for some missing stuff


# add the root branch since that will not be made in the above loop
branches$time_end[which(branches$parent_node=='root')] <- branches$time[which(branches$node == "root")]

print("Making a dataframe with vertical connectors.")
# make a df with the verticle branches (connectors)
connectors <- data.frame(node = branches$node[which(branches$time != 0 & branches$node != "root")])
connectors$ymax <- NA
connectors$ymin <- NA
connectors$time <- NA
for (i in 1:nrow(connectors)){
  nodename <- connectors$node[i]
  ymax <- max(branches$ypos[which(!is.na(match(branches$parent_node, connectors$node[i])))])
  ymin <- min(branches$ypos[which(!is.na(match(branches$parent_node, connectors$node[i])))])
  time <- branches$time[match(connectors$node[i], branches$node)]
  connectors$ymax[i] <- ymax
  connectors$ymin[i] <- ymin
  connectors$time[i] <- time
  }

# add species group column for fill
branches$group <- NA
branches <- branches %>%
  mutate(group = ifelse(node %in% c("C.denti","C.wolfi","C.pogonias","C.mona"), "mona",group)) %>%
  mutate(group = ifelse(node %in% c("C.neglectus"), "neglectus",group)) %>%
  mutate(group = ifelse(node %in% c("C.cephus"), "cephus",group)) %>%
  mutate(group = ifelse(node %in% c("C.mitis","C.nictitans"), "mitis",group))

# make a df for drawing background grid
print("Making the background grid dataframe.")
grid_density <- 1000000
bg <- data.frame(time = seq(from = -grid_density, to = 0- max(branches$time), by = - grid_density),
                 ymin = min(na.omit(branches$ypos)) - max(na.omit(branches$ne)),
                 ymax = max(na.omit(branches$ypos)) + max(na.omit(branches$ne)))

# figure out the gene flow events
print("Prepping a gene flow arrows data frame.")
gf <- data.frame(source = character(), recipient = character(), time = double(),
                 prop = double(), ystart = double(), yend = double(), phy = double())
colnames(means.joint) <- c("value","stat")
means.joint$name <- rownames(means.joint)
means.joint$name <- gsub("t_\\d*", "", means.joint$name)
for (i in 1:length(hybrid_dests)){
  hybdest <- hybrid_dests[i]
  hybsource <- hybrid_sources[i]
  # find the hybdest in the hybrid branches df
  hybdestbranch <- hybrid_branches$dest[which(hybrid_branches$through_1 == hybdest | hybrid_branches$through_2 == hybdest)]
  hybsourcebranch <- hybrid_branches$dest[which(hybrid_branches$through_1 == hybsource | hybrid_branches$through_2 == hybsource)]
  p <- means.joint$value[means.joint$stat == paste0("phi_",hybdest, "<-", hybsource)]
  
  t <- means.joint$value[means.joint$name == hybdest]
  # print(data.frame(source = hybsourcebranch, recipient = hybdestbranch,
  #                  time = t, ystart = NA, yend = NA))
  gf <- rbind(gf, data.frame(source = hybsourcebranch, recipient = hybdestbranch,
                             time = t, ystart = NA, yend = NA, phi = round(p,4),
                             phi_lower = round(conf.joint$lower[conf.joint$stat == paste0("phi_",hybdest, "<-", hybsource)],4),
                             phi_upper = round(conf.joint$upper[conf.joint$stat == paste0("phi_",hybdest, "<-", hybsource)],4)))
}

print("Adding Y-positions to all branches")
# add the ypositions
for (i in 1:nrow(gf)){
  source <- gf$source[i]
  #print(source)
  #print(recipient)
  recipient <- gf$recipient[i]
  ycent_s <- branches$ypos[branches$node == source]
  ycent_r <- branches$ypos[branches$node == recipient]
  if (ycent_s > ycent_r){
    # in this case we should connect the lower boundry of the source branch
    # with the upper of the recipient
    ystart <- ycent_s - branches$ne[branches$node == source] / 2
    yend <- ycent_r + branches$ne[branches$node == recipient] / 2
  }
  else {
    # otherwise switch them around
    ystart <- ycent_r - branches$ne[branches$node == recipient] / 2
    yend <- ycent_s + branches$ne[branches$node == source] / 2
  }
  gf$ystart[i] <- ystart
  gf$yend[i] <- yend
}

print("Making the plot...")
## let's try to draw some branches
p <- ggplot() +
  # draw background grid every million years
  geom_segment(data = bg, aes(x = time, xend = time, y = ymin, yend = ymax),
               linetype = "dashed", linewidth = 0.5, color = "lightgrey") +
  ylim(c(min(branches$ypos) - max(branches$ne - max(branches$ne))),max(branches$ypos) + max(branches$ne) + max(branches$ne)) +
  geom_rect(data = branches, aes(ymin = ypos - (0.5 * ne), ymax = ypos + 0.5 * ne, xmin = 0 - time_end, xmax = 0 -time,
                                 fill = group)) +
  # add connectors
  geom_segment(data = connectors, aes(x = 0-time, xend = 0-time, y = ymin, yend = ymax)) +
  # add gene flow arrows
  geom_segment(data = gf, aes(x = 0-time, xend = 0-time, y = ystart,
                                    yend = yend),
               arrow = arrow(length = unit(0.5, "cm"))) +
  # add labels to gene flow arrows 
  geom_label_repel(data = gf, aes(x = 0-time, y = ystart - ((ystart - yend) /2), 
                                  label = paste0(phi,"\n",phi_lower,"-",phi_upper))) +
  # add error bars/95 confint
  geom_segment(data = branches, aes(x = 0-confint.lower, xend = 0-confint.upper, y = ypos, yend = ypos),
               linewidth = 3, alpha = .5) +
  # and whiskers to those
  geom_segment(data = branches, aes(x = 0-confint.lower, xend = 0-confint.lower, y = ypos - 50000, yend = ypos + 50000)) +
  geom_segment(data = branches, aes(x = 0-confint.upper, xend = 0-confint.upper, y = ypos - 50000, yend = ypos + 50000)) +
  # put a point at internal nodes just for fun
  geom_point(data = branches[branches$time != 0,], aes(x = 0 - time, y = ypos),
             size = 3) +
  # add an ne scale bar
  geom_segment(aes(x = -500000, xend = -500000, y = min(na.omit(branches$ypos) - 300000), 
                   yend = min(na.omit(branches$ypos) - 400000))) +
  geom_segment(aes(x=-500000 - 100000, xend = -500000 + 100000,
                   y = min(na.omit(branches$ypos) - 300000),
                   yend = min(na.omit(branches$ypos) - 300000))) +
  geom_segment(aes(x=-500000 - 100000, xend = -500000 + 100000,
                   y = min(na.omit(branches$ypos) - 400000),
                   yend = min(na.omit(branches$ypos) - 400000))) +
  geom_text(aes(label = "Ne = 100,000", x = -350000, y = -300000),
            hjust = 0, vjust = 1) +
  # add species labels
  geom_text(data = branches[branches$time == 0,], aes(x = 100000, y = ypos, label = node),
            hjust = 0) +
  xlim(c(0-max(branches$time), 2000000)) +
  theme_void() +
  theme(legend.position = "none")
p

ggsave(paste0(NAME, "/",NAME,"tree.pdf"), p)

}