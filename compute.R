## This document is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2, or (at
## your option) any later version.

## This document is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.

## Oscar Perpiñán Lamigueiro 2011

library(nloptr)
library(lattice)
library(latticeExtra)
library(colorspace)
library(solaR)

source('functions.R')

#####################################
#####       WIRING COSTS     ########
#####################################

s=c(4, 6, 10, 16, 25, 35, 50, 70, 95)

##Cable libre de halógenos
RZBuried=c(17.55, 17.6, 19.8, 22.24, 29.51, 38.85, 48.01, 59.34, 73.87)
RZConduit=c(20.4, 20.27, 22.44, 24.96, 31.96, 41.46, 50.89, 62.32, 76.27)

trellis.device(pdf, file='WiringCosts.pdf')
xyplot(RZConduit+RZBuried~s,
       type=c('p', 'r'), par.settings=standard.theme(color=FALSE),
       xlab='Section (mm²)', ylab='Unitary Cost (euro/m)',
       auto.key=list(space='right'))
dev.off()

lmConduit <- lm(RZConduit~s)
lmBuried <- lm(RZBuried~s)

##Example with Lew=50, Lns=30
calcSection(Lew=50, Lns=30, Linv=60)
calcSection(Lew=50, Lns=30, Linv=1)

############################################
#####       Meteorological Data     ########
############################################

##Sevilla, Spain
lat=37.2;
G0dm=c(2766, 3491, 4494, 5912, 6989, 7742, 7919, 7027, 5369, 3562, 2814,
     2179)
Ta=c(10, 14.1, 15.6, 17.2, 19.3, 21.2, 28.4, 29.9, 24.3, 18.2, 17.2, 15.2)
prom=list(G0dm=G0dm, Ta=Ta)
prod0 <- prodGCPV(lat=lat, dataRad=prom, modeTrk='two')
prevG0 <- as(prod0, 'G0')

#################################
#####       Trackers     ########
#################################

##24 trackers, six groups of 4 trackers
struct2x=list(W=23.11, L=9.8, Nrow=6, Ncol=4)
W=struct2x$W


#######################################################################################
## Calculations for a system with Lew=50, Lns=30, and struct2x, with
## the meteorological conditions defined by prevG0 and where the land
## costs are 2.5€/m2 (1000€/Ha during 25 years),
## http://www.solarweb.net/forosolar/aspectos-economicos-legales-administrativos/1404-alquiler-terreno-rustico-huerta-solar-2.html,
## The cost of the equipments is 3.8€/Wp (modules 2.5€/Wp, 0.3€/Wp
## inverters and 1€/Wp trackers)
########################################################################################

calcSystem(50, 30, struct=struct2x, dataRad=prevG0, cLand=2.5, cPV=3.8, cWire=lmConduit)


##############################
####     Optimization     ####
##############################
##Optimize a system with struct2x and prevG0, with costs cLand and
## cPV
cLand=2.5
cPV=4

##NELDER-MEAD

ceNelderMead <- nloptr(x0=c(50, 30), eval_f=calcCostEnergy,
              lb=c(23.11, 9.8),
              struct=struct2x, dataRad=prevG0, cLand=cLand, cPV=cPV, cWire=lmConduit,
               opts = list("algorithm"="NLOPT_LN_NELDERMEAD"))
ceNelderMead

calcSystem(ceNelderMead$solution[1], ceNelderMead$solution[2],
           struct=struct2x, dataRad=prevG0, cLand=cLand, cPV=cPV)

##COBYLA

ceCOBYLA <- nloptr(x0=c(50, 30), eval_f=calcCostEnergy,
              lb=c(23.11, 9.8),
              struct=struct2x, dataRad=prevG0, cLand=cLand, cPV=cPV, cWire=lmConduit,
              opts = list("algorithm"="NLOPT_LN_COBYLA"))
ceCOBYLA

calcSystem(ceCOBYLA$solution[1], ceCOBYLA$solution[2],
           struct=struct2x, dataRad=prevG0, cLand=cLand, cPV=cPV)

###############################################
#### Output of Ce over a grid of distances ####
###############################################

dists <- expand.grid(Lew=seq(40, 50, .5), Lns=seq(20, 30, .5), cLand=seq(1, 4, .5), cPV=seq(2.5, 5, .5))

costs <- apply(dists, 1, function(x, ...)calcCostEnergy(c(x[1], x[2]), cLand=x[3], cPV=x[4],...),
               struct=struct2x, dataRad=prevG0)

costMatrix <- cbind(dists, costs)
costMatrix$GRR <- with(costMatrix, Lew*Lns/with(struct2x, L*W))
costMatrix$rCost <- with(costMatrix, costs/ave(costs, cLand, cPV, FUN=min))

useOuterStrips(
               levelplot(rCost~Lew*Lns|cPV*cLand,
                         data=costMatrix,
                         as.table=TRUE, par.settings=custom.theme.2()))

xx <- loess(rCost~Lew+Lns+cPV+cLand, data=costMatrix)
grid <- expand.grid(Lew=seq(40, 50, .1), Lns=seq(20, 30, .1), cLand=seq(1, 4, .5), cPV=seq(2.5, 5, .5))
rCost <- as.numeric(predict(xx, grid))
grid <- cbind(grid, rCost)
grid$GRR <- with(grid, Lew*Lns/with(struct2x, L*W))


####Graphical output
myTheme <- function(pch=19, cex=0.7, region=sequential_hcl(12, power=3.5), ...){
  theme <- custom.theme(pch=pch, cex=cex, region=region, ...)
  theme$strip.background$col='grey90'
  theme$strip.shingle$col='grey90'
  theme$strip.border$col='transparent'
  theme$grid.pars=list(fontfamily='serif')
  theme
}

##data = grid (loess)
trellis.device(pdf, file='matrix.pdf')
useOuterStrips(
               levelplot(rCost~(Lew/W)*(Lns/W)|cPV*cLand, data=grid, 
                         aspect='iso',
                         contour=TRUE,
                         between=list(x=0.3, y=0.3),
                         xlab=expression(l[ew]),
                         ylab=expression(l[ns]),
                         scales=list(x=list(cex=0.7, rot=30),
                           y=list(cex=0.7)),
                         par.settings=myTheme(),
                         subscripts=TRUE),
               strip=strip.custom(var.name='PV',
                 strip.names=c(TRUE, TRUE),
                 strip.levels=c(TRUE, TRUE),
                 par.strip.text=list(cex=0.7)),
               strip.left=strip.custom(var.name='Land',
                 strip.names=c(TRUE, TRUE),
                 strip.levels=c(TRUE, TRUE),
                 par.strip.text=list(cex=0.7))
               ) +
    layer(panel.contourplot(Lew/W, Lns/W, GRR, at=c(4, 5, 6),
                            labels=list(cex=0.6, col='black'), col='gray40', 
                          region=FALSE, contour=TRUE, subscripts=TRUE), data=grid) +
  layer({idx <- which.min(z[subscripts])
         panel.points(x[subscripts][idx], y[subscripts][idx], pch=4, col='black')
         })
dev.off()

## data= costMatrix (no loess)

trellis.device(pdf, file='minimumLocation.pdf')
xyplot(rCost~Lew/W+Lns/W, outer=TRUE, groups=cut(GRR, 5),
       data=costMatrix, subset=(cLand==2 & cPV==3.5),
       auto.key=list(corner=c(1,0.8), title='GRR', cex.title=1),
       par.settings=list(
         grid.pars=list(fontfamily='serif'),
         strip.background=list(col='transparent'),
         strip.border=list(col='transparent')),
       strip=strip.custom(factor.levels=c(expression(l[ew]), expression(l[ns]))),
       between=list(x=0.3),
       xlab='Normalized distance', ylab='Relative cost',
       scales=list(x=list(rot=45, cex=0.6, relation='free'))
       ) + layer(panel.abline(h=1, lty=3, lwd=0.7))
dev.off()

##########################################################################
####  Optimization for different combinations of cLand and cPV costs  ####
##########################################################################

casos <- expand.grid(cLand=seq(2, 3, .2), cPV=seq(3, 5, .2))
optimResults <- apply(casos, 1, function(x)fOptim(x[1], x[2], struct=struct2x, dataRad=prevG0))
optimResults <- cbind(casos, do.call('rbind', optimResults))

###loess
gridOPTIM <- expand.grid(cLand=seq(2, 3, .01), cPV=seq(3, 5, .01))

lGRR <- loess(GRR~cLand*cPV, data=optimResults)
lenergyCost <- loess(energyCost~cLand*cPV, data=optimResults)
lPVRatio <- loess(PVRatio~cLand*cPV, data=optimResults)
llew <- loess(Lew/W~cLand*cPV, data=optimResults)
llns <- loess(Lns/W~cLand*cPV, data=optimResults)


gridOPTIM$GRR <- as.numeric(predict(lGRR, gridOPTIM))
gridOPTIM$energyCost <- as.numeric(predict(lenergyCost, gridOPTIM))
gridOPTIM$PVRatio <- as.numeric(predict(lPVRatio, gridOPTIM))
gridOPTIM$lew <- as.numeric(predict(llew, gridOPTIM))
gridOPTIM$lns <- as.numeric(predict(llns, gridOPTIM))

##Graphical output
myPalette=custom.theme(
            region = diverge_hcl(5, 
              c = c(100, 0), l = c(50, 90), power = 1)
  )

myPalette$layout.heights =
        list(top.padding = 0,
 	    main.key.padding = 0,
 	    key.axis.padding = 0,
 	    axis.xlab.padding = 0,
 	    xlab.key.padding = 0,
 	    key.sub.padding = 0,
 	    bottom.padding = 0)

myPalette$layout.widths =
        list(left.padding = 0,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 0,
 	    axis.key.padding = 0,
 	    right.padding = 0)

myPalette$grid.pars=list(fontfamily='serif')
  
trellis.device(pdf, file='GRRoptim.pdf')
levelplot(GRR~cLand*cPV, data=gridOPTIM,
          xlab='Land costs', ylab='Equipment costs',
          contour=TRUE, labels=TRUE, pretty=TRUE,
          aspect='iso',
          par.settings=myPalette)
dev.off()

trellis.device(pdf, file='Lnsoptim.pdf')
levelplot(lns~cLand*cPV, data=gridOPTIM,
          xlab='Land costs', ylab='Equipment costs',
          contour=TRUE, labels=TRUE, pretty=TRUE,
          aspect='iso',
          par.settings=myPalette)
dev.off()

trellis.device(pdf, file='Lewoptim.pdf')
levelplot(lew~cLand*cPV, data=gridOPTIM,
          xlab='Land costs', ylab='Equipment costs',
          contour=TRUE, labels=TRUE, pretty=TRUE,
          aspect='iso',
          par.settings=myPalette)
dev.off()
