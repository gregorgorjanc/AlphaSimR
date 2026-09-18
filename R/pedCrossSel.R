#' @title Pedigree cross with selection
#'
#' @description
#' Creates a \code{\link{Pop-class}} from a generic
#' pedigree and a set of founder individuals. 
#' Selection through the pedigree occurs by matching the 
#' ID individuals with best criteria and the one with 
#' the most offspring in that generation. 
#'
#' @param founderPop a \code{\link{Pop-class}}
#' @param id a vector of unique identifiers for individuals
#' in the pedigree. The values of these IDs are seperate from
#' the IDs in the founderPop if matchID=FALSE.
#' @param mother a vector of identifiers for the mothers
#' of individuals in the pedigree. Must match one of the
#' elements in the id vector or they will be treated as unknown.
#' @param father a vector of identifiers for the fathers
#' of individuals in the pedigree. Must match one of the
#' elements in the id vector or they will be treated as unknown.
#' @param matchID indicates if the IDs in founderPop should be
#' matched to the id argument. See details. 
#' DO NOT WORK
#' @param maxCycle the maximum number of loops to make over the pedigree
#' to sort it.
#' @param simParam an object of 'SimParam' class
#' @param use indicates if the pedigree shall be modified to respect selection. Do not modify pedigree structure.  
#' Match the best individual based on criterion (see sel_trait) to the individual with the most offspring.  
#' Options for use are 'gv', 'ebv' or 'pheno' 
#' ONLY GV WORKS! there is no propagation of ebv and pheno within makeCross2
#' TO ADD 'function'
#' @param sel_trait indicates which trait to use to apply use on the pedigree. 
#' @param sex a vector of 'M' and 'F' indicating the sex of individuals in pedigree. 
#' Necessary to rank sires and dams when param use is used. 
#'
#' @description
#' The way in which the user supplied pedigree is used depends on
#' the value of matchID. If matchID is TRUE, the IDs in the user
#' supplied pedigree are matched against founderNames. If matchID
#' is FALSE, founder individuals in the user supplied pedigree are
#' randomly sampled from founderPop.
#'
#' @examples
 # Create founder haplotypes
 founderPop = quickHaplo(nInd=5, nChr=2, segSites=10)

 # Set simulation parameters
 SP = SimParam$new(founderPop)
# \dontshow{SP$nThreads = 1L}
 SP$addTraitA(nQtlPerChr = 2L, name = 'TraitA')
 SP$addTraitA(nQtlPerChr = 3L, name = 'TraitB')
 # Create population
 pop = newPop(founderPop, simParam=SP)

 # Pedigree 
 #
 #test no recoding 
 #id = 2:11
 #mother = c(0,0,0,0,2,2,3,8,0,8)
 #father = c(0,0,0,0,4,4,5,6,7,7)
 #sex = c("M","M","F","F","F","F","M","F","F","M")
 
 id = 1:10
 mother = c(0,0,0,0,1,1,2,7,0,7)
 father = c(0,0,0,0,3,3,4,5,6,6)
 sex = c("M","M","F","F","F","F","M","F","F","M")
 
 
 pop2 = pedigreeCross2(pop, id, mother, father, sex, use = 'gv', matchID = FALSE, sel_trait = 1, simParam=SP)
#'
#'
#' @export

pedigreeCross2 = function(pop, id, mother, father, matchID=FALSE, use = NULL, sel_trait = 0L, sex = NULL,
                         maxCycle = 100, simParam = NULL){
  
  if(is.null(simParam)){
    simParam = get("SP", envir=.GlobalEnv)
  }
  
  if (!is.null(use)){
    if (is.null(sex)) {
      sex = get("sex", envir=.GlobalEnv)
    }
    if (is.null(sex)){
      stop('Missing sex information')
    }
  }
  
  
  # Coerce input data
  id = as.character(id)
  mother = as.character(mother)
  father = as.character(father)
  
  # Sort pedigree (identifies potential problems)
  # need recoded ped!! ID 2:11 won't give the same results as 1:10 despite same ped struture (issue coming from output dataframe creation and match function within)
  # if not recoded, makeCross2 will give an error as ped and output created will be wrong 
  sortPed = function(id, mother, father, maxCycle=maxCycle){
    nInd = length(id)
    output = data.frame(gen=integer(nInd),
                        id=as.character(id),
                        mother=match(mother, id),
                        father=match(father, id),
                        motherID=as.character(mother),
                        fatherID=as.character(father))
    unsorted = rep(TRUE, nInd)
    for(gen in 1:maxCycle){
      for(i in which(unsorted)){
        if(is.na(output$mother[i])&is.na(output$father[i])){
          # Is a founder
          output$gen[i] = 1
          unsorted[i] = FALSE
        }else if(is.na(output$mother[i])){
          # Mother is a founder
          if(!unsorted[output$father[i]]){
            output$gen[i] = output$gen[output$father[i]] +1 
            unsorted[i] = FALSE
          }
        }else if(is.na(output$father[i])){
          # Father is a founder
          if(!unsorted[output$mother[i]]){
            output$gen[i] = output$gen[output$mother[i]] +1 
            unsorted[i] = FALSE
          }
        }else{
          # Both parents are in the pedigree
          if(!unsorted[output$mother[i]] & !unsorted[output$father[i]]){
            output$gen[i] = pmax(output$gen[output$mother[i]], output$gen[output$father[i]]) + 1
            unsorted[i] = FALSE
          }
        }
      }
    }
    if(any(unsorted)){
      stop("Failed to sort pedigree, may contain loops or require a higher maxGen")
    }
    return(output)
  }

  ped = sortPed(id=id, mother=mother, father=father, 
                maxCycle=100)
  
  # Create list for new population
  output = vector("list", length=length(id))
  selection = rep('NA', length(id))
  
  # Order and assign founders
  unknownMotherAndFather = is.na(ped$mother) & is.na(ped$father)
  unknownMotherOnly =  is.na(ped$mother) & !is.na(ped$father)
  unknownFatherOnly = !is.na(ped$mother) &  is.na(ped$father)
  founderNames = c(unique(id[unknownMotherAndFather]),
                   unique(id[unknownMotherOnly]),
                   unique(id[unknownFatherOnly]))
  
  nFounder = length(founderNames)
  
  if(matchID){
    stop("matchID=TRUE does not work currently!")
    # need to recode founder pop for matchID which may be complex
    # Check that all founders are present
    founderPresent = founderNames %in% pop@id
    if(!all(founderPresent)){
      stop(paste("The following founders are missing:", founderNames[!founderPresent]))
    }
  }
  
  
  # Check that there are enough founders
  if(nFounder>pop@nInd){
    stop(paste("Pedigree requires",nFounder,"founders, but only",pop@nInd,"were supplied"))
  }

  #check that all info is present for selection
  if (!is.null(use)){
    #stop("use does not work currently!")
    if (use == 'gv'){
      if (length(pop@gv) == 0 | any(is.na(pop@gv))){
        stop(paste("There is missing information for 'gv'"))
      }
    } else if (use=='ebv'){
      stop("use 'ebv' does not work currently!")
      
      if (length(pop@ebv) == 0 | any(is.na(pop@ebv))){
        stop(paste("There is missing information for 'ebv'"))
      }
    } else if (use == 'pheno'){
      stop("use 'pheno' does not work currently!")
      
      if (length(pop@pheno) == 0 | any(is.na(pop@pheno))){
        stop(paste("There is missing information for 'pheno'"))
      }
    }
    else {stop(paste("Options for use are 'gv', 'ebv' or 'pheno'"))}
  }
  

  #checking if trait for use to rewrite the pedigree does exist
  if (!is.null(use)){
    trait = as.integer(sel_trait)
    if (pop@nTraits < trait){
      stop(paste('The selected trait was not defined for this population'))
    }
  }
  
  ###necessary for use option 
  #counting the number of offspring for each id (in the next gen and sex specific)
  countsire = function(tmp, x) { 
    tbl= table(tmp$father[tmp$gen == (as.numeric(tmp$gen[x])+1)])
    ret = ifelse(tmp$id[x] %in% names(tbl), tbl[tmp$id[x]], ifelse(tmp$sex[x] == "M", 0, NA))
    return(ret)
  }
  
  countdam = function(tmp, x) { 
    tbl= table(tmp$mother[tmp$gen == (as.numeric(tmp$gen[x])+1)])
    ret = ifelse(tmp$id[x] %in% names(tbl), tbl[tmp$id[x]], ifelse(tmp$sex[x] == "F", 0, NA))
    return(ret)
  }
  
  #switch ID within sire and dam columns 
  create_id = function(tmp, x){
    tbl = tmp[tmp$cat==tmp$cat[x], ]
    tbl = tbl[order(tbl$rank_off), ]
    tbl$new_id = tbl$id[order(tbl$rank_selection)]
    return(tbl$new_id[tbl$id == x])
  }
  
  ###
  
  # Prepare founders for unknownMotherAndFather
  n1 = 1
  n2 = sum(unknownMotherAndFather)
  pop@id[n1:n2] = id[unknownMotherAndFather]
  pop@mother[n1:n2] = rep("0", n2)
  pop@father[n1:n2] = rep("0", n2)
  
  # Prepare founders for unknownMotherOnly
  n = sum(unknownMotherOnly)
  if(n>=1){
    n1 = n2 + 1
    n2 = n2 + n
    pop@id[n1:n2] = paste("-", id[unknownMotherOnly], sep = "")
    pop@mother[n1:n2] = rep("0", n2-n1+1)
    pop@father[n1:n2] = rep("0", n2-n1+1) 
  }
  
  # Prepare founders for unknownFatherOnly
  n = sum(unknownFatherOnly)
  if(n>=1){
    n1 = n2 + 1
    n2 = n2 + n
    pop@id[n1:n2] = paste("-", id[unknownFatherOnly], sep = "")
    pop@mother[n1:n2] = rep("0", n2-n1+1)
    pop@father[n1:n2] = rep("0", n2-n1+1)
  }
  
  #rewrite ped needs to be run after a first run of makeCross (we need 'use' information for the whole ped)
  rewritePed = function(ped, output, sel_trait) {
    #create a tmp df to determine which ID needs to be switched 
    tmp= as.data.frame(cbind(ped$gen, ped$id, ped$mother, ped$father, sex))
    
    if(use == 'gv') {
      for (j in 1:length(output)){
        if (!is.null(output[[j]])){
          selection[j] = output[[j]]@gv[as.integer(sel_trait)]
        } else if (is.null(output[[j]])){
          selection[j] = NA
        }
      }
      tmp$selection = as.numeric(selection)
    } else if (use == "pheno") {
      #for (j in 1:length(output)){
      #  if (!is.null(output[[j]])){
      #    selection[j] = output[[j]]@pheno[as.integer(sel_trait)]
      #  } else if (is.null(output[[j]])){
      #    selection[j] = NA
      #  }
      #}
      #tmp$selection = as.numeric(selection)
    } else if (use == "ebv"){
      
      #for (j in 1:length(output)){
      #  if (!is.null(output[[j]])){
      #    selection[j] = output[[j]]@ebv[as.integer(sel_trait)]
      #  } else if (is.null(output[[j]])){
      #    selection[j] = NA
      #  }
      #}
      #tmp$selection = as.numeric(selection) 
    } 
    ##needs to work on subset of output, won't have full ped in there 
    colnames(tmp) = c('gen', 'id', 'mother', 'father', 'sex', 'selection')
    
    #offspring count column
    tmp$off_father = sapply(1:nrow(tmp), FUN = function(x) countsire(tmp, x))
    tmp$off_mother = sapply(1:nrow(tmp), FUN = function(x) countdam(tmp, x))
    
    #rank based on the number of offspring
    tmp$rank_father =  with(tmp, ave(as.numeric(tmp$off_father), tmp$gen, FUN=function(x) rank(-x, na.last='keep', ties.method = "min")))
    tmp$rank_mother =  with(tmp, ave(as.numeric(tmp$off_mother), tmp$gen, FUN=function(x) rank(-x, na.last='keep', ties.method = "min")))
    
    #create one column for ranking based on number of off (ranking still specific to gen and sex)
    tmp$rank_off = apply(cbind(tmp$rank_father, tmp$rank_mother), 1, function(x) paste(x[!is.na(x)], collapse = ""))
    
    #remove unnecessary columns
    tmp[ , c('rank_father', 'rank_mother', 'off_father', 'off_mother')] = NULL
    
    #create rank based on $selection within gen*sex
    tmp$cat = paste(tmp$gen, tmp$sex)
    tmp = transform(tmp, rank_selection = ave(1:nrow(tmp), cat, FUN=function(x) order(sex[x], -as.numeric(selection[x]), decreasing=F)))
    
    #new_id to be used in pedigree (only on father and mother columns)
    tmp$new_id = sapply(1:nrow(tmp), FUN = function(x) create_id(tmp, x))
    
    #id list 
    dict = as.data.frame(cbind(tmp$id, tmp$new_id))
    colnames(dict) = c('old_id', 'new_id')
    
    #new id for the pedigree
    tmp$father_id = dict$new_id[match(tmp$father, dict$old_id)]
    tmp$mother_id = dict$new_id[match(tmp$mother, dict$old_id)]
    
    #reset original ped formatting
    tmp$mother = tmp$mother_id
    tmp$mother[is.na(tmp$mother)] = 0
    tmp$father = tmp$father_id
    tmp$father[is.na(tmp$father)] = 0
    
    new_ped = tmp[, c('gen', 'id', 'mother', 'father', 'mother_id', 'father_id')]
    
    return(new_ped)
  }
  
  
  # Create individuals
  crossPlan = matrix(c(1,1),ncol=2)
  for (gen in 1:max(ped$gen)){
    for (i in which(ped$gen==gen)){
      if (unknownMotherAndFather[i]){
        # Copy over founder individual
        output[[i]] = pop[id[i]]
      } else if(unknownMotherOnly[i]){
        # Cross founder to newly created individual
        if(!is.null(use)){
          ped = rewritePed(ped, output, sel_trait)
        }
        output[[i]] = makeCross2(pop[paste("-", id[i], sep = "")],
                                 output[[as.numeric(ped$father[i])]],
                                 crossPlan=crossPlan,
                                 simParam=simParam)
        output[[i]]@mother = "0"
      } else if(unknownFatherOnly[i]){
        # Cross newly created individual to founder
        if (!is.null(use)){
          ped = rewritePed(ped, output, sel_trait)
        }
        output[[i]] = makeCross2(output[[as.numeric(ped$mother[i])]],
                                 pop[paste("-", id[i], sep = "")],
                                 crossPlan=crossPlan,
                                 simParam=simParam)
        output[[i]]@father = "0"
      } else {
        # Cross two newly created individuals
        if (!is.null(use)){
          ped = rewritePed(ped, output, sel_trait)
        }
        output[[i]] = makeCross2(output[[as.numeric(ped$mother[i])]],
                                 output[[as.numeric(ped$father[i])]],
                                 crossPlan=crossPlan,
                                 simParam=simParam)
      }
    }
  }
  
  # Collapse list to a population
  output = mergePops(output)
  
  # Copy over names
  output@id = id
  output@mother = mother
  output@father = father
  sex = get("sex",envir=.GlobalEnv)
  output@sex = sex
  
  
  return(output)
  
}




