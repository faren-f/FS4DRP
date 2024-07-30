library(progeny)
Progeny_pw_act = function(X){
  
  pathway_activity= progeny(
    t(X),
    scale = TRUE,
    organism = "Human",
    top = 100,
    perm = 1
    
  )
  return(pathway_activity)
  }

