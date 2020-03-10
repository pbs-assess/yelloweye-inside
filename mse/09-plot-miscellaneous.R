
# Miscellaneous plots generally associated with conditioning OM
library(MSEtool)

# Plot B/BMSY estimated from state-space surplus production model
SRA_list <- readRDS("mse/scoping/scoping_base.rds")
SRA <- SRA_list[[1]]; ret <- SRA_list[[2]]

SP_mod <- SP_SS(Data = SRA@OM@cpars$Data, AddInd = 1:5, use_r_prior = TRUE, start = list(r_prior = c(0.068, 0.03)))
#plot(SP_mod, dir = getwd(), filename = "mse/scoping/scoping_SP", open_file = FALSE)

# Little to no difference with reconstructed catch (no doubling 1986-2005)
#Data_lowcatch <- SRA@OM@cpars$Data
#Data_lowcatch@Cat[, match(1986:2005, Data_lowcatch@Year)] <- 0.5 * Data_lowcatch@Cat[, match(1986:2005, Data_lowcatch@Year)]
#SP_mod2 <- SP_SS(Data = Data_lowcatch, AddInd = 1:5, use_r_prior = TRUE, start = list(r_prior = c(0.068, 0.03)))

png("mse/figures/SP_fit.png", height = 3.5, width = 4.5, units = "in", res = 400)
par(mar = c(5, 4, 1, 1))
plot(as.numeric(names(SP_mod@B_BMSY)), SP_mod@B_BMSY, typ = "l", lwd = 3, xlab = "Year", ylab = expression(B/B[MSY]),
     ylim = c(0, 2.1))
abline(h = 0, col = "grey")
abline(h = c(0.4, 0.8), lty = 3)
dev.off()
