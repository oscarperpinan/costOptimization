## This document is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2, or (at
## your option) any later version.

## This document is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.

## Oscar Perpiñán Lamigueiro 2011

cu <- function(s, model=lmConduit){##cost for a section
  coefs <- coef(model)
  coefs[1] + coefs[2]*s
  ##predict(model, data.frame(s=s))
}

##section for a group of four trackers feeding an inverter.
calcSection <- function(Lns, Lew, Linv, Ig=52.2, Vnom=506.9, norm=TRUE){
  L1=Lew+Lns+1
  L2=Lns+1
  L3=Lew+1
  L4=1##The junction box is next to the four tracker (1 meter)
  Ls=c(L1, L2, L3, L4)

  It=4*Ig##Current to the inverter, 4 trackers

  DeltaV=1.5/100*Vnom

  F=1+sqrt(Ig*sum(Ls^2)/(It*Linv^2))

  DeltaVinv=DeltaV/F
  sInv=2*Linv*It/(56*DeltaVinv)##section of the inverter wire

  DeltaVs=DeltaV-DeltaVinv

  s=2*Ig*Ls/(56*DeltaVs)##section of each tracker circuit

  result <- c(s, sInv)

  if (norm) {##choose the correspondent normalized wire
    snorm=c(6, 10, 16, 25, 35, 50, 70, 95, 120, 150, 185, 240, 300, 400, 500, 630)
    idx=colSums(outer(snorm, result, '<='))+1
    snorm[idx]
  } else {result}
}


##Calculate sections, energy production and costs for a system defined
##by Lew, Lns and struct. The PV plant is composed of 6 groups of 4
##trackers around an inverter building as described in the paper
##cLand is €/m2 (whole life cycle!)
##cPV is €/Wp
calcSystem <- function(Lew, Lns, struct, dataRad, cLand=2.5, cPV=3.8, cWire=lmConduit){
  ##productividad y energía
  distances=data.frame(Lew=Lew, Lns=Lns, H=0)
  lat=getLat(dataRad, 'deg');
  prodShd <- prodGCPV(lat=lat, modeRad='prev', dataRad=dataRad, modeTrk='two',
                      modeShd='prom', struct=struct, distances=distances)
  Yf <- as.data.frameY(prodShd)$Yf
  Pg=with(struct, Nrow*Ncol)*prodShd@generator$Pg ##Potencia total en Wp
  Eac <- Yf*Pg/1000 ##Energía Anual en kWh!!!
  ##cable
  cable1 <- calcSection(Lew=Lew, Lns=Lns, Linv=2*Lns)
  cable2 <- calcSection(Lew=Lew, Lns=Lns, Linv=1)
  CosteUnitario1 <- cu(cable1, model=cWire)
  CosteUnitario2 <- cu(cable2, model=cWire)
  Ln1=c(Lns+Lew+1, Lns+1, Lew+1, 1, 2*Lns)
  Ln2=c(Lns+Lew+1, Lns+1, Lew+1, 1, 1)
  CosteCable1 <- 4*sum(CosteUnitario1*Ln1)
  CosteCable2 <- 2*sum(CosteUnitario2*Ln2)
  CosteCable <- CosteCable1+CosteCable2
  ##terreno
  GRR=(Lew*Lns)/with(struct, L*W)
  At=(struct$Nrow*Lew)*(struct$Ncol*Lns)
  CosteTerreno <- cLand*At 
  ##equipos
  CosteFV <- cPV*Pg 
  ##total
  Coste <- CosteFV+CosteTerreno+CosteCable
  ##coste energía
  ce <- Coste/(Eac*25)*100      #Coste energía durante 25 años, c€/kWh
  ##entrego resultados
  result <- data.frame(GRR=GRR,
                       Eac=Eac,
                       wireCost=CosteCable,
                       wireRatio=CosteCable/Coste,
                       landCost=CosteTerreno,
                       landRatio=CosteTerreno/Coste,
                       PVCost=CosteFV,
                       PVRatio=CosteFV/Coste,
                       cost=Coste,
                       energyCost=ce)
  result
}

##Simple wrapper to obtain the energy cost from calcSystem
calcCostEnergy <- function(x, struct, dataRad, cLand=2.5, cPV=3.8, cWire=lmConduit){
  output <- calcSystem(Lew=x[1], Lns=x[2], struct=struct, dataRad=dataRad,
             cLand=cLand, cPV=cPV, cWire=cWire)
  res <- output$energyCost
  }

##Function to optimize a system with design parameters (struct, meteoData) and costs (cLand, cPV)
fOptim <- function(cLand, cPV, theta=c(50, 30), ci=c(23.11, 9.8), cWire=lmConduit, struct, dataRad){
  ceOptim <- nloptr(x0=theta, eval_f=calcCostEnergy,
              lb=ci,
              struct=struct, dataRad=dataRad, cLand=cLand, cPV=cPV, cWire=cWire,
              opts = list("algorithm"="NLOPT_LN_COBYLA"))
  res <- calcSystem(Lew=ceOptim$solution[1], Lns=ceOptim$solution[2],
                    struct=struct, dataRad=dataRad,
                    cLand=cLand, cPV=cPV, cWire=lmConduit)
  cbind(res, data.frame(Lew=ceOptim$solution[1], Lns=ceOptim$solution[2]))
}
