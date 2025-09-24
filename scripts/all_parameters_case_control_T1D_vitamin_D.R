################################################################################
load_all_parameter <- function(){
  all_parameters <- data.frame(
    matrix(ncol = 4, nrow = 0)
  )
  colnames(all_parameters) <- c("variable2categorise",
                                "category_variable",
                                "breaks",
                                "labels")
  
  
  # first_degree_psoriasis_heredity
  variable2categorise = list("first_degree_T1D_heredity")
  category_variable = list("First_degree_T1D_heredity_cat")
  breaks = list(c(0,1,2))
  labels = list(c("No", "Yes"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  )
  
  # week of delivery
  variable2categorise = list("fr_12")
  category_variable = list("Preterm_group")
  breaks = list(c(0,28,32,37,Inf))
  labels = list(c("Extreme preterm", "Very preterm", "Moderate preterm", "Normal"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  )    
  
  # week of delivery
  variable2categorise = list("fr_12")
  category_variable = list("Preterm_dichotomous")
  breaks = list(c(0,37,Inf))
  labels = list(c("Preterm", "Normal"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  )      
  
  # delivery
  variable2categorise = list("fr_13")
  category_variable = list("Delivery_mode")
  breaks = list(c(0,1,2,3))
  labels = list(c("Normal", "caesarean", "other problem"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  ) 
  
  # Did you use any vitamins or mineral supplement during your pregnancy?
  variable2categorise = list("fr_81")
  category_variable = list("Vitamin_mineral_supplement_pregnancy")
  breaks = list(c(0,1,2,3))
  labels = list(c("Yes", "No", "Do not know"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  ) 
  
  # "How often did you eat lake fish?"
  variable2categorise = list("fr_91")
  category_variable = list("Fish_lake_pregnancy")
  breaks = list(c(0,1,2,3,4))
  labels = list(c("Every day", "3\u20135 times/week", "1\u20132 times week", "Seldom"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  ) 
  
  # "How often did you eat Baltic fish ?"
  variable2categorise = list("fr_92")
  category_variable = list("Fish_Baltic_pregnancy")
  breaks = list(c(0,1,2,3,4))
  labels = list(c("Every day", "3\u20135 times/week", "1\u20132 times week", "Seldom"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  ) 
  
  # "How often did you eat other fish?"
  variable2categorise = list("fr_93")
  category_variable = list("Fish_other_pregnancy")
  breaks = list(c(0,1,2,3,4))
  labels = list(c("Every day", "3\u20135 times/week", "1\u20132 times week", "Seldom"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  ) 
  
  # "How often did you any fish?"
  variable2categorise = list("fish_intake_frequency")
  category_variable = list("Fish_any_pregnancy")
  breaks = list(c(0,1,2,3,4))
  labels = list(c("Every day", "3\u20135 times/week", "1\u20132 times week", "Seldom"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  ) 
  
  # smoking during pregnancy
  variable2categorise = list("fr_18")
  category_variable = list("Smoking_during_pregnancy")
  breaks = list(c(0,1,2))
  labels = list(c("No", "Yes"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  ) 
  
  # parity
  variable2categorise = list("parity")
  category_variable = list("Parity")
  breaks = list(c(0,1,2))
  labels = list(c("First parity", "Previous parity"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  ) 
  
  # Gender
  variable2categorise = list("Gender")
  category_variable = list("Sex_cat")
  breaks = list(c(0,1,2))
  labels = list(c("Boy", "Girl"))
  
  all_parameters <- rbind(
    all_parameters, 
    data.frame(
      I(variable2categorise), I(category_variable), I(breaks), I(labels)
    )
  ) 
  
  return(as.data.frame(all_parameters))
}

all_parameters <- load_all_parameter()