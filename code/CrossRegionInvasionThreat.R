#################################################################################
###This function uses annual invasion probabilities from INApest core functions
###To estimate invasion threat to other regions
###Allows national-scale assessment of consequences of inaction 
###and impacts of control
###Outputs are line graphs of N annual invasion events and annual farm-level invasion risk 
###For the sink region and a .csv storing these results
#################################################################################


CrossRegionInvasionThreat = function
(
ThreatSource, ###Source region
ThreatSink = ThreatSink, ###Sink region
LDDrate = LongDistProbPerFarm, ###annual long distance dispersal events per source farm
SourceInvasionProb = SourceInvasionProb, ###Annual farm-level invasion prob from simulations
CrossRegionInvasionDir=CrossRegionInvasionDir , ###Directory where cross-region invasion probs stored
EIDir = InputDir, ##Directory where Ecoclimatic Index data stored
OutputFileNameStem = OutputFileNameStem
)
{

##############################################################
###Assign climate based establishment prob to single farms based on ecoclimatic index
###Default function is logit
##############################################################
SinkEI = data<-read_csv(paste0(EIDir,ThreatSink,"_Simple_points_for_INA.csv"))
prob_est<-as.vector(SinkEI$Probability_Estab)
#if(EI_Prob_CurveType == "SplitLinear")
#	prob_est<-as.vector(data$Probability_Estab_SplitLinear)

##############################################
###Read in farm-level cross-region invasion risk
##############################################
CrossRegionSDDmatrix = readRDS(paste0(CrossRegionInvasionDir,ThreatSource,"_",ThreatSink,"_InvasionProb.rds"))
CrossRegionLDDmatrix = readRDS(paste0(CrossRegionInvasionDir,ThreatSource,"_",ThreatSink,"_InvasionProb_LDDweight.rds"))
CrossRegionLDDmatrix = CrossRegionLDDmatrix*LDDrate
CrossRegionDispersal = CrossRegionSDDmatrix+CrossRegionLDDmatrix
CrossRegionDispersal = sweep(CrossRegionDispersal,1,prob_est,"*")


DispersalEvents = vector(length = ncol(SourceInvasionProb))
for(year in 1:ncol(SourceInvasionProb))
{
AnnualInvasion = sweep(CrossRegionDispersal,2,SourceInvasionProb[,year],"*")
DispersalEvents[year] = sum(AnnualInvasion[,])
}
FarmLevelDispersalRisk = DispersalEvents/nrow(CrossRegionDispersal)
SourceInfestedYear = 1:ncol(SourceInvasionProb)
OutputTable = data.frame(SourceInfestedYear,DispersalEvents,FarmLevelDispersalRisk)
Filename = paste0(OutputFileNameStem,ThreatSource,"_",ThreatSink,".csv")
write.table(OutputTable,Filename,sep = ",",row.names = F)
Title = paste0("Detection Prob. ",DetectionProb," Erradication Prob. ",round(AnnualErradicationProb,digits =3),
		";\nSpread Reduction ",SpreadReduction,"\n",
             ThreatSource, " to ", ThreatSink)
Filename = paste0(OutputFileNameStem,ThreatSource,"_",ThreatSink,".png")
png(Filename)
plot(1:length(DispersalEvents),DispersalEvents,pch = NA,
	xlab = "Years source region infested", ylab = "Mean annual establishment events",
    main = Title)
lines(1:length(DispersalEvents),DispersalEvents, lwd = 3)
dev.off()

Filename = paste0(OutputFileNameStem,ThreatSource,"_",ThreatSink,"FarmRisk.png")
png(Filename)
plot(1:length(DispersalEvents),FarmLevelDispersalRisk,pch = NA,
	xlab = "Years source region infested", ylab = "Mean per-farm annual establishment events",
    main = Title)
lines(1:length(DispersalEvents),FarmLevelDispersalRisk, lwd = 3)
dev.off()
}

########################################################################
###End of function
########################################################################
