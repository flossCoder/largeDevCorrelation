# gatherStatisticsKSTest.R
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

load(paste(directory2, "/", param, "/ana-ks.RData", sep = ""))

cat(sprintf("%s %s %s %s %s\n",
            toString(format(cinter), digits=22),
            toString(format(cintra), digits=22),
            toString(format((cintra + cinter) / 2.0, digits=22)), # average connectivity c
            toString(format(ksResult$p.value, digits=22)), # save the p-value of the KS-test
            toString(format(log(ksResult$p.value), digits=22))), # save the logarithm of the p-value of the KS-test
    file = paste(directory2, "/ks-statistics.dat", sep = ""),
    append = TRUE)

cat(sprintf("%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
            toString(cinter),
            toString(cintra),
            toString((cintra + cinter) / 2.0), # average connectivity c
            toString(integralData), # integral of the sbm rate function
            toString(integralER), # integral of the er rate function
            toString(integralDifference), # difference of the integral of the two rate functions
            toString(absIntegralDifference), # difference of the absolut integral of the two rate functions
            toString(integralERAna), # integral of the er rate function
            toString(integralDifferenceAna), # difference of the integral of the two rate functions
            toString(absIntegralDifferenceAna), # difference of the absolut integral of the two rate functions
            toString(integralDataResult), # integral of the sbm rate function
            toString(integralERResult), # integral of the er rate function
            toString(integralDifferenceResult), # difference of the integral of the two rate functions
            toString(absIntegralDifferenceResult)), # difference of the absolut integral of the two rate functions
    file = paste(directory2, "/integral-statistics.dat", sep = ""),
    append = TRUE)

cat(sprintf("%s %s %s %s %s %s %s %s %s\n",
            toString(cinter),
            toString(cintra),
            toString((cintra + cinter) / 2.0), # average connectivity c
            toString(sumDiff), # sum difference of the two large deviation rate functions
            toString(sumAbsDiff), # absolute sum difference of the two large deviation rate functions
            toString(sumDiffAna), # sum difference of the two large deviation rate functions
            toString(sumAbsDiffAna), # absolute sum difference of the two large deviation rate functions
            toString(sumDiffResult), # sum difference of the two large deviation rate functions
            toString(sumAbsDiffResult)), # absolute sum difference of the two large deviation rate functions
    file = paste(directory2, "/sum-statistics.dat", sep = ""),
    append = TRUE)

cat(sprintf("%s %s %s %s %s %s %s %s %s\n",
            toString(cinter),
            toString(cintra),
            toString((cintra + cinter) / 2.0), # average connectivity c
            toString(meanDiff), # mean difference of the two large deviation rate functions
            toString(meanAbsDiff), # mean absolute difference of the two large deviation rate functions
            toString(meanDiffAna), # mean difference of the two large deviation rate functions
            toString(meanAbsDiffAna), # mean absolute difference of the two large deviation rate functions
            toString(meanDiffResult), # mean difference of the two large deviation rate functions
            toString(meanAbsDiffResult)), # mean absolute difference of the two large deviation rate functions
    file = paste(directory2, "/mean-statistics.dat", sep = ""),
    append = TRUE)
