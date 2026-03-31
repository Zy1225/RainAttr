#' @title Permutation-Based Procedure for Rainfall Enhancement Trial Data
#'
#' @description
#' Implements the permutation-based procedure used in \code{\link{rain_attr}} for inference on attribution and SATE in rainfall enhancement trial data under randomized assignment of ionizer operation.
#' It allows for flexible permutation schemes, including permutations between ionizers, between days, or between gauge-day combinations.
#'
#'
#' @param B_permutation An integer specifying the number of permutation replicates.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permute_between_ionizer Logical. If \code{TRUE}, for each day, a random permutation is performed among the operation statuses of all ionizers that have been deployed on that day.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permute_all_ionizers_between_day Logical. If \code{TRUE}, for each year, a random permutation is performed among the daily operation schedules of all trial days within that year.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permute_between_gaugeday Logical. If \code{TRUE}, for each year, a random permutation is performed among the gauge-day level operation schedule of all gauge-days within that year.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param ionizer_operation_input A data frame containing ionizer operation indicators for each day (row) and each ionizer (column), where 1 indicates that an ionizer is turned on, 0 indicates that it is off and NA indicates that it is not deployed yet.
#' Additionally, this data frame must include two columns with names specified by \code{ionizer_operation_day_column_name} and \code{ionizer_operation_year_column_name}, containing the day and year for each row.
#' Each day must appear only once in this data frame (no duplicated day entries).
#' The ionizer columns must appear in the same order as specified by \code{data_target_column_names} and must be consistent with the column order in \code{gaugeday_downwind_input}.
#'   (User-supplied using \code{\link{permutation_option}})
#' @param gaugeday_downwind_input A binary matrix indicating which gauge-day observations (row) are downwind of which ionizers (column), where 1 indicates that the gauge is downwind of the ionizer on that day, 0 indicates that the gauge is not downwind of the ionizer, and NA indicates that the ionizer has not been deployed yet.
#'   The row order of this matrix must match that of \code{data}. The column order must correspond to the ionizer columns in \code{ionizer_operation_input} (excluding the day and year columns) and be in the same order as specified by \code{data_target_column_names}.
#'   (User-supplied using \code{\link{permutation_option}})
#' @param data_target_column_names A character vector specifying the column names of \code{data} corresponding to the binary target indicators used in the downwind (second stage) LMM fitting.
#'   The order of names in this vector must match the column order of the corresponding ionizers in \code{ionizer_operation_input} (excluding the day and year columns) and in \code{gaugeday_downwind_input}.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param ionizer_operation_year_column_name A character string specifying the column name of \code{ionizer_operation_input} containing the year of each day.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param ionizer_operation_day_column_name A character string specifying the column name of \code{ionizer_operation_input} containing the day of each observation. The same column name should also be found in \code{data}.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permutation_seed An integer specifying the random seed for the permutation-based procedure. Reproducibility is guaranteed only if \code{permutation_parallel} is the same, since parallel execution changes the order of random number generation.
#' (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permutation_parallel Logical. If \code{TRUE}, each permutation run is executed in parallel across multiple workers. If \code{FALSE}, they are run sequentially.
#' (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permutation_parallel_num_worker An integer specifying the number of parallel workers to use when \code{permutation_parallel = TRUE}.
#' (User-configurable permutation option using \code{\link{permutation_option}})
#' @param data A data frame containing the original dataset used in \code{\link{rain_attr}}, along with an additional column containing the fitted values generated from the upwind (first stage) LMM.
#'   Its column names should contain \code{data_target_column_names} (binary target indicators) and \code{ionizer_operation_day_column_name} (day).
#'   The row order of this data frame must match that of \code{gaugeday_downwind_input}.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind_lmm_formula A two sided linear formula object to be used in \link[lme4]{lmer}, describing both the fixed-effects and random intercept part of the downwind (second stage) LMM.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind_propensity_formula A two sided linear formula object to be used in \code{\link{glm}} with \code{family = "binomial"}, for fitting a propensity score model to the treatment indicators of downwind (second stage) observations.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param attr_type A character string specifying the type of attribution estimates. Must be one of \code{"ChambersEtAl"}, \code{"ChambersEtAl_No_Winsorize"}, \code{"ThoEtAl"}, or \code{"No"}. See \code{\link{rain_attr}} for more information.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param x_downwind_name A character vector containing variable names from the right hand side of \code{downwind_lmm_formula}, for those variables that are not related to ionizers (treatment). The intercept is always included and does not need to be specified.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param target_only Logical. If \code{TRUE} the attribution estimates are computed based on only target observations. If \code{FALSE} the attribution estimates are computed based on both treatment and control observations.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param rain_col_name A character string specifying the column name of the raw scale rainfall in \code{ori_data}.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#'
#' @return A list with two components:
#' \describe{
#'   \item{hatattr}{Matrix of permutation replicates of attribution estimates.}
#'   \item{hatsate}{Matrix of permutation replicates of SATE estimates.}
#' }
#'
#' @details
#' This function implements the permutation-based procedure used in \code{\link{rain_attr}}. It is intended for internal use only. Users should not call this function directly.
#' Instead, permutation-based inference should be performed by calling \code{\link{rain_attr}} with \code{permutation = TRUE}.
#'
#  \strong{Additional Information Required for Permutation-Based Procedure} \cr
#' To perform permutation-based procedure, additional information need to be supplied through the following arguments of \code{\link{permutation_option}}:
#' \describe{
#'  \item{\code{ionizer_operation_input}}{Day(group)-level ionizers operation schedule during the rainfall enhancement trial.}
#'  \item{{gaugeday_downwind_input}}{Gauge-day(unit within group)-level information on relative orientation of gauges from ionizers each day.}
#' }
#'
#' \strong{Permutation Steps} \cr
#' The permutation-based procedure considers to randomly permute ionizers' operation statuses via:
#' \enumerate{
#'    \item{If \code{permute_between_ionizer = TRUE}, for each row of the day-level \code{ionizer_operation_input} matrix, a random permutation is performed among the binary indicators in the row that correspond to ionizers that have already been deployed during the year of the row.
#'    This is equivalent to randomly permuting the operation statuses of deployed ionizers for each day. }
#'    \item{If \code{permute_all_ionizers_between_day = TRUE}, for each year, a random permutation is performed among all rows belonging to the year in the day-level \code{ionizer_operation_input} matrix.
#'    This is equivalent to randomly permuting the daily operation schedules of all trial days belonging to each year.  }
#'    \item{The permuted day-level \code{ionizer_operation_input} are then expanded into a gauge-day level binary ionizers' operation indicator matrix, to match the gauge-day level \code{data}. }
#'    \item{If \code{permute_between_gaugeday = TRUE}, for each year, a random permutation is performed among all rows belonging to the year in the gauge-day level binary ionizers' operation indicator matrix.
#'    This is equivalent to randomly permuting the gauge-day level operation schedules of all gauge-days belonging to each year. }
#' }
#'
#' The above three optional permutation steps are performed in sequence, resulting in a final permuted gauge-day level binary ionizers' operation indicator matrix.
#' An elementwise multiplication is carried out between this matrix and \code{gaugeday_downwind_input}, and the results are used to replace the original columns (\code{data_target_column_names}) in \code{data} that contain the binary target indicators used in the downwind (second stage) LMM fitting, where NAs (for cases where the ionizer has not been deployed yet) are replaced by zeros.
#' Therefore, it is important to ensure that
#' \itemize{
#'  \item{The row order of \code{gaugeday_downwind_input} and \code{data} is consistent.}
#'  \item{The column orders of \code{gaugeday_downwind_input} and \code{ionizer_operation_input} is consistent, and matches the order specified by \code{data_target_column_names}.}
#' }
#' Based on these permuted binary target indicators, the original binary indicators \eqn{I_{ij}} for exposure to ionizers (treatment) are also updated to be their permuted counterparts \eqn{I_{ij}^*}, where \eqn{I_{ij}^*} only equals zero if all permuted binary target indicators in its corresponding row equal zero.
#'
#'
#'  \strong{Permutation Distribution of Attribution and SATE} \cr
#' The same procedure described in \code{\link{rain_attr}} is then repeated on the permuted dataset, where the columns containing the binary target indicators and the binary indicators \eqn{I_{ij}} for exposure to ionizers (treatment) are replaced according to the permutation.
#  This includes the fitting of the downwind (second stage) LMM, downwind (second stage) target-only LMM, downwind (second stage) control-only LMM to the subset of observations from \code{data} that are downwind (second stage) and with positive rainfall, where a gauge-day level observation is said to be downwind if the gauge is downwind of at least one deployed ionizer (not necessarily turned on) on that day i.e., at least one of the indicators is equal to one for that row of \code{gaugeday_downwind_input}.
#' This includes the fitting of the downwind (second stage) LMM, downwind (second stage) target-only LMM, downwind (second stage) control-only LMM to the subset of observations from \code{data} that are downwind (second stage) and with positive rainfall, along with the fitting of downwind propensity score model to the subset of observations from \code{data} that are downwind (second stage) with the response being the permuted indicator \eqn{I_{ij}^*} for exposure to ionizers (treatment).
#' Finally, two attribution estimates and the SATE estimates are computed based on the estimation results of these models fitted to the permuted dataset, where \eqn{z_{ij}} (ionizer related covariate vector constructed from the binary target indicators) and \eqn{I_{ij}} are replaced by their permuted counterparts \eqn{z_{ij}^*} and \eqn{I_{ij}^*}, respectively.
#'
#' By repeatedly permuting ionizers' operation schedules, fitting models and computing attribution and SATE estimates for \code{B_permutation} number of times, this function returns the permutation distributions of
#' \itemize{
#'    \item{Two attribution estimates}
#'    \item{Five SATE estimates}
#' }
#'
#' \strong{Permutation-Based P-Value and Permutation-Based Plots} \cr
#' The permutation distributions of attribution and SATE produced by this function are further used in \code{\link{rain_attr}} to:
#' \itemize{
#'    \item{Compute permutation-based p-value as the proportion of permuted estimates that are greater than or equal to the original estimate, i.e.,
#'      \deqn{\frac{1}{B} \sum_{b=1}^{B} 1_{ \{ \hat{\theta}^*_{b} \geq \hat{\theta} \}  } ,}
#'      where \eqn{\hat{\theta}^*_{b}} is the estimate from the \eqn{b}-th permuted dataset and \eqn{\hat{\theta}} is the corresponding estimate from the original dataset.
#'    }
#'    \item{Plot the kernel density estimate with a solid vertical line for the original estimate \eqn{\hat{\theta}}. }
#' }
#'
#' @seealso \code{\link{rain_attr}} for the main function, \code{\link{permutation_option}} for specifying permutation options



permutation_ionizer = function(B_permutation, permute_between_ionizer, permute_all_ionizers_between_day, permute_between_gaugeday,
                               ionizer_operation_input, gaugeday_downwind_input,
                               data_target_column_names, ionizer_operation_year_column_name, ionizer_operation_day_column_name,
                               permutation_seed, permutation_parallel, permutation_parallel_num_worker,
                               data, downwind_lmm_formula, downwind_propensity_formula,
                               attr_type, x_downwind_name, target_only,
                               rain_col_name){
  #data_target_column_names are the column names of 'data', which correspond to the target indicators of all ionizers
  #ionizer_operation_year_column_name is the column name of 'ionizer_operation_input' containing the year of each day
  #ionizer_operation_day_column_name is the column name of 'ionizer_operation_input' containing the day of each observation, which should be the same column name in 'data'
  #Note that data_target_column_names should have the same ordering as colnames(ionizer_operation_input), as well as gaugeday_downwind_input
  #The rows of gaugeday_downwind_input should be ordered as the same as the order of rows of gauge-day observations in 'data'
  if(length(unique(ionizer_operation_input[,ionizer_operation_day_column_name])) != nrow(ionizer_operation_input)){
    stop('ionizer_operation_input has more than one rows associated to the same day')
  }

  year_ionizer_list = lapply(
    unique(ionizer_operation_input[,ionizer_operation_year_column_name]), function(x){
      deployed_ionizers_x = which(colMeans(!is.na(ionizer_operation_input[ ionizer_operation_input[,ionizer_operation_year_column_name] == x  , -which(colnames(ionizer_operation_input) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))])) == 1)
      return(
        setdiff(colnames(ionizer_operation_input),c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))[deployed_ionizers_x]
      )
    }
  )

  names(year_ionizer_list) = unique(ionizer_operation_input[,ionizer_operation_year_column_name])


  perm_attr_matrix = matrix(data = NA, nrow = B_permutation, ncol = 2, dimnames = list(NULL, c('apo','apl')))
  perm_sate_matrix = matrix(data = NA, nrow = B_permutation, ncol = 5, dimnames = list(NULL, c('sate.mb','sate.ipw','sate.ipw.l','sate.ipw.ma','sate.aipw')))

  # ionizer_operation_yearlist = lapply(
  #   names(year_ionizer_list), function(x){
  #     ionizer_operation_input[ionizer_operation_input[,ionizer_operation_year_column_name] == x, -which(colnames(ionizer_operation_input) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name) )  ]
  #   }
  # )
  # names(ionizer_operation_yearlist) = names(year_ionizer_list)

  if(!permutation_parallel){
    if(!is.null(permutation_seed)){
      set.seed(permutation_seed)
    }

    downwind = apply(gaugeday_downwind_input,1,function(x){sum(x, na.rm=TRUE)}) > 0
    positive = ( data[,rain_col_name] > 0)

    test_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = data[downwind & positive,])
    z_downwind_name = setdiff(names(lme4::fixef(test_downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
    perm_downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)
    perm_downwind_propensity_formula = update.formula(downwind_propensity_formula, permuted_target_indicator ~ . )





    for(b in 1:B_permutation){
      tryCatch({
        perm_data = data



        perm_ionizer_operation_day = ionizer_operation_input
        if(permute_between_ionizer){
          for(unique_year in unique(perm_ionizer_operation_day[, ionizer_operation_year_column_name])){
            temp = perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))]
            deployed_ionizers = year_ionizer_list[[as.character(unique_year)]]
            perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))] =
              t(apply(temp,1, FUN = function(x){
                x[colnames(temp) %in% deployed_ionizers] = sample(x[colnames(temp) %in% deployed_ionizers])
                return(x)
              }))
          }
        }

        if(permute_all_ionizers_between_day){
          for(unique_year in unique(perm_ionizer_operation_day[, ionizer_operation_year_column_name])){
            temp = perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))]
            perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))] =
              temp[sample(1:nrow(temp)),]
          }
        }


        # perm_ionizer_operation_yearlist = ionizer_operation_yearlist
        # if(permute_between_ionizer){
        #   for(i in 1:length(perm_ionizer_operation_yearlist)){
        #     deployed_ionizers = year_ionizer_list[[names(perm_ionizer_operation_yearlist)[i]]]
        #     perm_ionizer_operation_yearlist[[i]] = t(apply(perm_ionizer_operation_yearlist[[i]], 1, function(x){
        #       x[colnames(x) %in% deployed_ionizers] = sample(x[colnames(x) %in% deployed_ionizers])
        #       return(x)
        #     }))
        #   }
        # }
        #
        # if(permute_all_ionizers_between_day){
        #   for(i in 1:length(perm_ionizer_operation_yearlist)){
        #     perm_ionizer_operation_yearlist[[i]] = perm_ionizer_operation_yearlist[[i]][sample(1:nrow(perm_ionizer_operation_yearlist[[i]])),]
        #   }
        # }
        #
        # perm_ionizer_operation_day = ionizer_operation_input
        # for(i in 1:length(year_ionizer_list)){
        #   perm_ionizer_operation_day[perm_ionizer_operation_day[,ionizer_operation_year_column_name] == names(year_ionizer_list)[i], -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name) )  ] = perm_ionizer_operation_yearlist[[i]]
        # }

        perm_ionizer_operation_gaugeday = dplyr::left_join(perm_data[, ionizer_operation_day_column_name, drop = FALSE], perm_ionizer_operation_day, by = ionizer_operation_day_column_name)

        if(permute_between_gaugeday){
          for(unique_year in unique(perm_ionizer_operation_gaugeday[, ionizer_operation_year_column_name])){
            temp = perm_ionizer_operation_gaugeday[perm_ionizer_operation_gaugeday[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_gaugeday) %in% c(ionizer_operation_day_column_name, ionizer_operation_year_column_name))  ]
            perm_ionizer_operation_gaugeday[perm_ionizer_operation_gaugeday[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_gaugeday) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name))  ] =
              temp[sample(1: nrow(temp)),]
          }
        }

        perm_data[,data_target_column_names] = (perm_ionizer_operation_gaugeday[, -which(colnames(perm_ionizer_operation_gaugeday) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name))] * gaugeday_downwind_input)

        perm_data[,data_target_column_names][is.na(perm_data[,data_target_column_names])] = 0

        #Refit downwind LMM to permuted data, and then compute hatattr and hatsate, noting we need to becareful with the definition of downwind_positive_target, downwind_positive_control, and downwind_propensity_
        perm_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = perm_data[downwind & positive,])

        perm_downwind_positive_data = perm_data[downwind & positive,]
        target_vec = apply(perm_data[,data_target_column_names],1, sum) > 0
        nontarget_vec = apply(perm_data[,data_target_column_names],1, sum) ==  0
        perm_downwind_positive_target = target_vec[downwind & positive]
        perm_downwind_positive_control = nontarget_vec[downwind & positive]

        perm_hatattr = attr_est(attr_type, perm_downwind_positive_data, rain_col_name, perm_downwind_positive_target, perm_downwind_positive_control,
                                x_downwind_name, target_only = target_only, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)


        perm_downwind_positive_data$permuted_target_indicator = as.logical(perm_downwind_positive_target)


        perm_hatsate = sate_est(perm_downwind_positive_data, perm_downwind_positive_target, perm_downwind_positive_control, perm_downwind_propensity_formula, perm_downwind_separate_formula,
                                x_downwind_name, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

        perm_attr_matrix[b,] = c(perm_hatattr$apo, perm_hatattr$apl)
        perm_sate_matrix[b,] = c(perm_hatsate$estimates$sate.mb, perm_hatsate$estimates$sate.ipw, perm_hatsate$estimates$sate.ipw.l, perm_hatsate$estimates$sate.ipw.ma, perm_hatsate$estimates$sate.aipw)
      },error=function(e){cat(b,"th","Permutation Run Skipped due to ERROR :",conditionMessage(e), "\n")})
    }
  }else{
    downwind = apply(gaugeday_downwind_input,1,function(x){sum(x, na.rm=TRUE)}) > 0
    positive = ( data[,rain_col_name] > 0)

    test_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = data[downwind & positive,])
    z_downwind_name = setdiff(names(lme4::fixef(test_downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
    perm_downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)
    perm_downwind_propensity_formula = update.formula(downwind_propensity_formula, permuted_target_indicator ~ . )

    permutation_cl = parallel::makeCluster(permutation_parallel_num_worker)
    parallel::clusterExport(permutation_cl, varlist = c('attr_est', 'sate_est'))
    doParallel::registerDoParallel(permutation_cl)
    if(!is.null(permutation_seed)){
      doRNG::registerDoRNG(seed = permutation_seed)
    }
    `%dopar%` <- foreach::`%dopar%`


    permutation_results = foreach::foreach(b = 1:B_permutation, .packages = c("dplyr","lme4"), .errorhandling = 'remove') %dopar% {
      perm_data = data



      perm_ionizer_operation_day = ionizer_operation_input
      if(permute_between_ionizer){
        for(unique_year in unique(perm_ionizer_operation_day[, ionizer_operation_year_column_name])){
          temp = perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))]
          deployed_ionizers = year_ionizer_list[[as.character(unique_year)]]
          perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))] =
            t(apply(temp,1, FUN = function(x){
              x[colnames(temp) %in% deployed_ionizers] = sample(x[colnames(temp) %in% deployed_ionizers])
              return(x)
            }))
        }
      }

      if(permute_all_ionizers_between_day){
        for(unique_year in unique(perm_ionizer_operation_day[, ionizer_operation_year_column_name])){
          temp = perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))]
          perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))] =
            temp[sample(1:nrow(temp)),]
        }
      }


      # perm_ionizer_operation_yearlist = ionizer_operation_yearlist
      # if(permute_between_ionizer){
      #   for(i in 1:length(perm_ionizer_operation_yearlist)){
      #     deployed_ionizers = year_ionizer_list[[names(perm_ionizer_operation_yearlist)[i]]]
      #     perm_ionizer_operation_yearlist[[i]] = t(apply(perm_ionizer_operation_yearlist[[i]], 1, function(x){
      #       x[colnames(x) %in% deployed_ionizers] = sample(x[colnames(x) %in% deployed_ionizers])
      #       return(x)
      #     }))
      #   }
      # }
      #
      # if(permute_all_ionizers_between_day){
      #   for(i in 1:length(perm_ionizer_operation_yearlist)){
      #     perm_ionizer_operation_yearlist[[i]] = perm_ionizer_operation_yearlist[[i]][sample(1:nrow(perm_ionizer_operation_yearlist[[i]])),]
      #   }
      # }
      #
      # perm_ionizer_operation_day = ionizer_operation_input
      # for(i in 1:length(year_ionizer_list)){
      #   perm_ionizer_operation_day[perm_ionizer_operation_day[,ionizer_operation_year_column_name] == names(year_ionizer_list)[i], -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name) )  ] = perm_ionizer_operation_yearlist[[i]]
      # }

      perm_ionizer_operation_gaugeday = dplyr::left_join(perm_data[, ionizer_operation_day_column_name, drop = FALSE], perm_ionizer_operation_day, by = ionizer_operation_day_column_name)

      if(permute_between_gaugeday){
        for(unique_year in unique(perm_ionizer_operation_gaugeday[, ionizer_operation_year_column_name])){
          temp = perm_ionizer_operation_gaugeday[perm_ionizer_operation_gaugeday[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_gaugeday) %in% c(ionizer_operation_day_column_name, ionizer_operation_year_column_name))  ]
          perm_ionizer_operation_gaugeday[perm_ionizer_operation_gaugeday[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_gaugeday) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name))  ] =
            temp[sample(1: nrow(temp)),]
        }
      }

      perm_data[,data_target_column_names] = (perm_ionizer_operation_gaugeday[, -which(colnames(perm_ionizer_operation_gaugeday) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name))] * gaugeday_downwind_input)

      perm_data[,data_target_column_names][is.na(perm_data[,data_target_column_names])] = 0

      #Refit downwind LMM to permuted data, and then compute hatattr and hatsate, noting we need to becareful with the definition of downwind_positive_target, downwind_positive_control, and downwind_propensity_
      perm_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = perm_data[downwind & positive,])

      perm_downwind_positive_data = perm_data[downwind & positive,]
      target_vec = apply(perm_data[,data_target_column_names],1, sum) > 0
      nontarget_vec = apply(perm_data[,data_target_column_names],1, sum) ==  0
      perm_downwind_positive_target = target_vec[downwind & positive]
      perm_downwind_positive_control = nontarget_vec[downwind & positive]

      perm_hatattr = attr_est(attr_type, perm_downwind_positive_data, rain_col_name, perm_downwind_positive_target, perm_downwind_positive_control,
                              x_downwind_name, target_only = target_only, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)


      perm_downwind_positive_data$permuted_target_indicator = as.logical(perm_downwind_positive_target)



      perm_hatsate = sate_est(perm_downwind_positive_data, perm_downwind_positive_target, perm_downwind_positive_control, perm_downwind_propensity_formula, perm_downwind_separate_formula,
                              x_downwind_name, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

      perm_attr_b = c(perm_hatattr$apo, perm_hatattr$apl)
      perm_sate_b = c(perm_hatsate$estimates$sate.mb, perm_hatsate$estimates$sate.ipw, perm_hatsate$estimates$sate.ipw.l, perm_hatsate$estimates$sate.ipw.ma, perm_hatsate$estimates$sate.aipw)

      list(
        b = b,
        perm_attr_b = perm_attr_b,
        perm_sate_b = perm_sate_b
      )
    }

    parallel::stopCluster(permutation_cl)

    for(permutation_res in permutation_results){
      perm_attr_matrix[permutation_res$b,] = permutation_res$perm_attr_b
      perm_sate_matrix[permutation_res$b,] = permutation_res$perm_sate_b
    }

  }


  return(list(
    hatattr = perm_attr_matrix,
    hatsate = perm_sate_matrix
  ))
}


#' @title Permutation Options
#'
#' @description
#' This function generates a list of settings controlling how permutation-based procedures are executed within \code{\link{rain_attr}}.
#'
#'
#' @param B_permutation An integer specifying the number of permutation replicates. Default is 10000.
#' @param permute_between_ionizer Logical. If \code{TRUE}, for each day, a random permutation is performed among the operation statuses of all ionizers that have been deployed on that day. Default is \code{TRUE}.
#' @param permute_all_ionizers_between_day Logical. If \code{TRUE}, for each year, a random permutation is performed among the daily operation schedules of all trial days within that year. Default is \code{FALSE}.
#' @param permute_between_gaugeday Logical. If \code{TRUE}, for each year, a random permutation is performed among the gauge-day level operation schedule of all gauge-days within that year. Default is \code{TRUE}.
#' @param ionizer_operation_input A data frame containing ionizer operation indicators for each day (row) and each ionizer (column), where 1 indicates that an ionizer is turned on, 0 indicates that it is off, and NA indicates that it is not deployed yet.
#' Additionally, this data frame must include two columns with names specified by \code{ionizer_operation_day_column_name} and \code{ionizer_operation_year_column_name}, containing the day and year for each row.
#' Each day must appear only once in this data frame (no duplicated day entries).
#' The ionizer columns must appear in the same order as specified by \code{data_target_column_names} and must be consistent with the column order in \code{gaugeday_downwind_input}. Default is \code{ionizer_operation}.
#' @param gaugeday_downwind_input A binary matrix indicating which gauge-day observations (row) are downwind of which ionizers (column), where 1 indicates that the gauge is downwind of the ionizer on that day, 0 indicates that the gauge is not downwind of the ionizer, and NA indicates that the ionizer has not been deployed yet.
#'   The row order of this matrix must match that of the original dataset supplied to \code{\link{rain_attr}}. The column order must correspond to the ionizer columns in \code{ionizer_operation_input} (excluding the day and year columns) and be in the same order as specified by \code{data_target_column_names}. Default is \code{gaugeday_downwind}.
#' @param data_target_column_names A character vector specifying the column names of the original dataset supplied to \code{\link{rain_attr}}, corresponding to the binary target indicators used in the downwind (second stage) LMM fitting.
#'   The order of names in this vector must match the column order of the corresponding ionizers in \code{ionizer_operation_input} (excluding the day and year columns) and in \code{gaugeday_downwind_input}. Default is:
#'   \code{c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10")}.
#' @param ionizer_operation_year_column_name A character string specifying the column name of \code{ionizer_operation_input} containing the year of each day. Default is \code{'Year'}.
#' @param ionizer_operation_day_column_name A character string specifying the column name of \code{ionizer_operation_input} containing the day of each observation. The same column name should also be found in the original dataset supplied to \code{\link{rain_attr}}. Default is \code{'TrialDay'}.
#' @param permutation_seed An integer specifying the random seed for the permutation-based procedure. Reproducibility is guaranteed only if \code{permutation_parallel} is the same, since parallel execution changes the order of random number generation. Default is \code{NULL}, meaning no seed is set internally and users should call \code{set.seed()} beforehand to ensure reproducibility.
#' @param permutation_parallel Logical. If \code{TRUE}, each permutation run is executed in parallel across multiple workers. If \code{FALSE}, they are run sequentially. Default is \code{FALSE}.
#' @param permutation_parallel_num_worker An integer specifying the number of parallel workers to use when \code{permutation_parallel = TRUE}. Default is \code{parallel::detectCores() - 1}.
#'
#' @return A list containing all permutation options, suitable for passing to \code{\link{rain_attr}}.
#'
#' @details
#' This function is used to configure and store all settings needed for performing permutation-based analyses within \code{\link{rain_attr}}.
#' It is important to ensure that
#' \itemize{
#'  \item{The row order of \code{gaugeday_downwind_input} and the original dataset supplied to \code{\link{rain_attr}} is consistent.}
#'  \item{The column orders of \code{gaugeday_downwind_input} and \code{ionizer_operation_input} is consistent, and matches the order specified by \code{data_target_column_names}.}
#' }
#'
#' @seealso \code{\link{rain_attr}} for the main function, \code{\link{permutation_ionizer}} for more details on the permutation-based procedure
#'
#' @examples
#' #Create default permutation options
#' # These are the same permutation settings used in Chambers et al. (2022a) Nudging a Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall Enhancement Trial in Oman. International Statistical Review, 90: 346–373,
#' #  as well as Chambers et al. (2022b) Weighting, Informativeness and Causal Inference, with an Application to Rainfall Enhancement. \emph{Journal of the Royal Statistical Society Series A: Statistics in Society}, 185: 1584–1612.
#' # Specifically: permute_between_ionizer = TRUE, permute_all_ionizers_between_day = FALSE, and permute_between_gaugeday = TRUE
#' perm_options = permutation_opt()
#' str(perm_options)
#'
#' #Permutation option with parallelization over (parallel::detectCores() - 1) number of workers and seed = 1 for reproducibility
#' perm_options_parallel = permutation_opt(
#'   permutation_seed = 1,
#'   permutation_parallel = TRUE,
#'   permutation_parallel_num_worker = parallel::detectCores() - 1
#' )
#' str(perm_options_parallel)
#'
permutation_opt = function(B_permutation = 10000,
                           permute_between_ionizer = T,
                           permute_all_ionizers_between_day = F,
                           permute_between_gaugeday = T,
                           ionizer_operation_input = ionizer_operation,
                           gaugeday_downwind_input = gaugeday_downwind,
                           data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
                           ionizer_operation_year_column_name = 'Year',
                           ionizer_operation_day_column_name = 'TrialDay',
                           permutation_seed = NULL,
                           permutation_parallel = F,
                           permutation_parallel_num_worker = parallel::detectCores() - 1){
  return(list(
    B_permutation = B_permutation,
    permute_between_ionizer = permute_between_ionizer,
    permute_all_ionizers_between_day = permute_all_ionizers_between_day,
    permute_between_gaugeday = permute_between_gaugeday,
    ionizer_operation_input = ionizer_operation_input,
    gaugeday_downwind_input = gaugeday_downwind_input,
    data_target_column_names = data_target_column_names,
    ionizer_operation_year_column_name = ionizer_operation_year_column_name,
    ionizer_operation_day_column_name = ionizer_operation_day_column_name,
    permutation_seed = permutation_seed,
    permutation_parallel = permutation_parallel,
    permutation_parallel_num_worker = permutation_parallel_num_worker
  ))
}

permutation_p_value = function(permutation_result, ori_est){
  temp = sapply(1:ncol(permutation_result), FUN = function(i){
    mean(permutation_result[,i] >= ori_est[i], na.rm = T )
  })
  names(temp) = colnames(permutation_result)
  return(temp)
}

permutation_plot = function(permutation_result, ori_est){
  num_var = ncol(permutation_result)
  output_plot_list = lapply(1:num_var, function(i){
    ggplot2::ggplot() +
      ggplot2::geom_density(ggplot2::aes(x = permutation_result[,i])) +
      ggplot2::geom_vline(xintercept = ori_est[i]) +
      ggplot2::ggtitle(colnames(permutation_result)[i]) +
      ggplot2::xlab('Permuted Values') +
      ggplot2::ylab('Permutation Density') +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  })
  names(output_plot_list) = colnames(permutation_result)
  return(output_plot_list)
}
