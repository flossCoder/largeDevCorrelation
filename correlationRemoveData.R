# correlationRemoveData.R
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


# for none scientific file processing => needed for file names
options(scipen = 999)

load(paste(directory2, "ana-cor_chi-sqr.RData", sep = ""))

remove(data)

# Save all results as R workspace:
save.image(paste(directory2, "ana-cor_chi-sqr.RData", sep = ""))
