# gatherCorrelation.R
# Copyright (C) 2016, 2017 flossCoder
#
# This file is part of largeDevCorrelation.
#
# largeDevCorrelation is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# largeDevCorrelation is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

load(paste(directory2, "/", param, "/ana-cor_chi-sqr.RData", sep = ""))

input <- strsplit(param, "_")[[1]] # dirER contains simulation results for the corresponding ER graphs

cinter <- as.numeric(gsub("-", ".", input[length(input) - 1]))
cintra <- as.numeric(gsub("-", ".", input[length(input)]))

c <- (cinter + cintra) / 2.0 # calculate the corresponding er connectivity

# Compare largest components vs. number of components:
cat(sprintf("%s %s %s %s %s %s %s %s %s %s\n",
            toString(cinter),
            toString(cintra),
            toString(c), # average connectivity c
            toString(lc_vs_nc[[2]]$statistic), # X-squared Pearson no MC
            toString(lc_vs_nc[[2]]$parameter), # degrees of fredom Pearson no MC
            toString(lc_vs_nc[[2]]$p.value), # p.value Pearson no MC (might be wrong)
            toString(lc_vs_nc[[3]]$statistic), # X-squared Pearson MC
            toString(lc_vs_nc[[3]]$p.value), # p.value Pearson MC
            toString(lc_vs_nc[[4]]), # correlation coefficient
            toString(lc_vs_nc[[5]])), # qchisq for degree of freedom and 0.95
    file = paste(directory2, "/correlation-lc_vs_nc.dat", sep = ""),
    append = TRUE)

# Compare largest components vs. number of edges:
cat(sprintf("%s %s %s %s %s %s %s %s %s %s\n",
            toString(cinter),
            toString(cintra),
            toString(c), # average connectivity c
            toString(lc_vs_ne[[2]]$statistic), # X-squared Pearson no MC
            toString(lc_vs_ne[[2]]$parameter), # degrees of fredom Pearson no MC
            toString(lc_vs_ne[[2]]$p.value), # p.value Pearson no MC (might be wrong)
            toString(lc_vs_ne[[3]]$statistic), # X-squared Pearson MC
            toString(lc_vs_ne[[3]]$p.value), # p.value Pearson MC
            toString(lc_vs_ne[[4]]), # correlation coefficient
            toString(lc_vs_ne[[5]])), # qchisq for degree of freedom and 0.95
    file = paste(directory2, "/correlation-lc_vs_ne.dat", sep = ""),
    append = TRUE)

# Compare number of components vs number of edges:
cat(sprintf("%s %s %s %s %s %s %s %s %s %s\n",
            toString(cinter),
            toString(cintra),
            toString(c), # average connectivity c
            toString(nc_vs_ne[[2]]$statistic), # X-squared Pearson no MC
            toString(nc_vs_ne[[2]]$parameter), # degrees of fredom Pearson no MC
            toString(nc_vs_ne[[2]]$p.value), # p.value Pearson no MC (might be wrong)
            toString(nc_vs_ne[[3]]$statistic), # X-squared Pearson MC
            toString(nc_vs_ne[[3]]$p.value), # p.value Pearson MC
            toString(nc_vs_ne[[4]]), # correlation coefficient
            toString(nc_vs_ne[[5]])), # qchisq for degree of freedom and 0.95
    file = paste(directory2, "/correlation-nc_vs_ne.dat", sep = ""),
    append = TRUE)
