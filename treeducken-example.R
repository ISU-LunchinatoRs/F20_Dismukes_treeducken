# install.packages("treeducken")
library(treeducken)

set.seed(54)
lambda_H <- rexp(n=1)
mu_H <- lambda_H / 2
lambda_C <- rexp(n=1)
time <- 1.0

lambda_S <- rexp(n=1)
mu_S <- lambda_S / 2

H_tips <- treeducken::calculate_expected_leaves_sptree(lambda = lambda_H + lambda_C,
                                                       mu = mu_H,
                                                       t = time)

S_tips <- treeducken::calculate_expected_leaves_sptree(lambda = lambda_S + lambda_C,
                                                       mu = mu_S,
                                                       t = time)
cophy_obj <- treeducken::sim_cophylo_bdp(hbr = lambda_H, hdr = mu_H,
                                         sbr = lambda_S, sdr = mu_S,
                                         cosp_rate = lambda_C, host_exp_rate = 0.0,
                                         time_to_sim = time, numbsim = 1)
plot(cophy_obj[[1]], col = "orange", lty = "dotted")


print(cophy_obj, details = TRUE)
summary(cophy_obj[[1]])

df <- treeducken::cophy_summary_stat(cophy_obj)
D <- treeducken::parafit_stat(cophy_obj[[1]]$host_tree,
                              cophy_obj[[1]]$symb_tree,
                              cophy_obj[[1]]$association_mat)
treeducken::parafit_test(cophy_obj[[1]]$host_tree,
                         cophy_obj[[1]]$symb_tree,
                         cophy_obj[[1]]$association_mat,
                         D = D,
                         reps = 99)


gopher_lice_map <- read.table(system.file("extdata",
                                          "gopher_lice_mapping.txt",
                                          package = "treeducken"),
                              stringsAsFactors = FALSE, header = TRUE)

gopher_lice_assoc_matrix <- convert_assoc_table_to_matrix(gopher_lice_map)
gopher_tree <- ape::read.nexus(system.file("extdata",
                                           "gophers_bd.tre",
                                           package = "treeducken"))
lice_tree <- ape::read.nexus(system.file("extdata",
                                         "lice_bd.tre",
                                         package = "treeducken"))
gopher_lice_cophylo <- convert_to_cophy(hostTree = gopher_tree,
                                        symbTree = lice_tree,
                                        assocMat = gopher_lice_assoc_matrix)
print(gopher_lice_cophylo)
cophy_summary_stat(gopher_lice_cophylo)

plot(gopher_lice_cophylo,
     fsize = 0.5,
     show_tip_label = FALSE,
     gap = 1,
     col = "purple",
     lty = "dashed")



library(paco)
host_tree_pruned <- drop_extinct(cophy_obj[[1]]$host_tree)
symb_tree_pruned <- drop_extinct(cophy_obj[[1]]$symb_tree)
A <- association_mat(cophy_obj[[1]])
host_dist <- cophenetic(host_tree_pruned)
symb_dist <- cophenetic(symb_tree_pruned)
links <- t(A) # paco wants associations with rows as hosts 

# we need to name rows and columns for paco
rownames(links) <- host_tree_pruned$tip.label
colnames(links) <- symb_tree_pruned$tip.label 
D <- paco::prepare_paco_data(H = host_dist, P = symb_dist, HP = links)
D <- paco::add_pcoord(D)
D <- paco::PACo(D, nperm=100, seed = 11, method="r0")
print(D$gof)