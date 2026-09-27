# Example only: apply multi-target PoF to an existing INApest result object
# that contains both pathogen outputs and BiocontrolHistory.
source(file.path('..','src','INApestProofOfFreedom.R'))

# Replace these with your model result and pathogen specification.
# result <- readRDS('my_model_results.rds')
# pathogen <- my_pathogen_spec

# Direct surveillance of adult Q only; all Q stages still count against freedom.
q_surveillance <- INApestBiocontrolSurveillance(
  DetectionProb = 0.20,
  DetectStages = 'adult',
  DetectAgents = 'Q'
)

# Example observed histories. NA means no observation was supplied at that timestep.
# observations <- list(
#   pathogen = c(0, 0, NA, 0),
#   biocontrol = c(0, NA, 0, 0)
# )
#
# fit <- INApestPoF(
#   ModelOutput = result,
#   Pathogen = pathogen,
#   Targets = c('pathogen','biocontrol'),
#   ObservationHistory = observations,
#   BiocontrolSurveillance = q_surveillance,
#   BiocontrolAgents = 'Q'
# )
#
# fit$Summary[, c(
#   'timestep',
#   'PosteriorPoF_Pathogen',
#   'PosteriorPoF_Biocontrol',
#   'PosteriorPoF_Joint'
# )]
