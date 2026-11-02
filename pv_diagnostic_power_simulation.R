# PV Diagnostic: Trial-Level Mixed-Model Power Simulation
# -------------------------------------------------------
# Simulation-based power analysis for the PV diagnostic ABX trials.
#
# Primary hypothesis:
#   Does the Bimodal group improve more from Pre to Post than
#   the Unimodal group?
#
# Primary coefficient:
#   phasePost:distributionBimodal
#
# Design represented in the simulation:
#   - 2 distribution groups: Unimodal, Bimodal
#   - 2 phases: Pre, Post
#   - 8 unique PV diagnostic items
#   - each item repeated twice
#   - 16 PV diagnostic trials per participant per phase
#
# Current complete PV participants:
#   PVB01, PVB02, PVB03, PVB05
#   PVU01, PVU02, PVU03, PVU04, PVU05
#
# PVB04 was excluded because Day 2 / Posttest was not completed.


# Packages -----------------------------------------------------------------

library(lme4)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)


# Current PV diagnostic pilot data -----------------------------------------

# Each PV item value is the number of correct responses.
#
# trials_per_item = 1:
#   item was presented once
#
# trials_per_item = 2:
#   item was presented twice

pv_wide <- tribble(

  ~participant, ~distribution, ~phase, ~trials_per_item,
  ~PV_A_1, ~PV_A_2, ~PV_A_3, ~PV_A_4,
  ~PV_R_1, ~PV_R_2, ~PV_R_3, ~PV_R_4,

  # PVB01
  "PVB01", "Bimodal", "Pre",  1,
  0, 1, 0, 1,
  1, 1, 1, 1,

  "PVB01", "Bimodal", "Post", 1,
  1, 1, 1, 1,
  0, 1, 1, 0,

  # PVB02
  "PVB02", "Bimodal", "Pre",  2,
  2, 2, 1, 2,
  1, 0, 2, 0,

  "PVB02", "Bimodal", "Post", 2,
  2, 2, 2, 2,
  2, 1, 0, 2,

  # PVB03
  "PVB03", "Bimodal", "Pre",  2,
  2, 0, 2, 0,
  2, 0, 2, 0,

  "PVB03", "Bimodal", "Post", 2,
  0, 2, 1, 2,
  2, 1, 0, 2,

  # PVB05
  "PVB05", "Bimodal", "Pre",  2,
  2, 2, 1, 2,
  1, 2, 2, 2,

  "PVB05", "Bimodal", "Post", 2,
  1, 2, 2, 2,
  1, 2, 2, 2,

  # PVU01
  "PVU01", "Unimodal", "Pre",  1,
  1, 1, 0, 0,
  0, 1, 0, 1,

  "PVU01", "Unimodal", "Post", 1,
  1, 1, 0, 1,
  1, 0, 0, 0,

  # PVU02
  "PVU02", "Unimodal", "Pre",  2,
  0, 1, 2, 2,
  2, 2, 0, 1,

  "PVU02", "Unimodal", "Post", 2,
  2, 1, 1, 0,
  2, 1, 1, 2,

  # PVU03
  "PVU03", "Unimodal", "Pre",  2,
  2, 0, 2, 0,
  2, 0, 2, 0,

  "PVU03", "Unimodal", "Post", 2,
  1, 0, 1, 0,
  1, 0, 1, 0,

  # PVU04
  "PVU04", "Unimodal", "Pre",  2,
  1, 1, 2, 1,
  0, 1, 1, 2,

  "PVU04", "Unimodal", "Post", 2,
  0, 2, 1, 1,
  1, 2, 2, 1,

  # PVU05
  "PVU05", "Unimodal", "Pre",  2,
  1, 1, 1, 2,
  2, 1, 1, 1,

  "PVU05", "Unimodal", "Post", 2,
  1, 2, 2, 1,
  0, 1, 0, 2
)


# Convert to item-level format ---------------------------------------------

pv <- pv_wide %>%

  pivot_longer(
    cols = starts_with("PV_"),
    names_to = "item",
    values_to = "correct"
  ) %>%

  mutate(
    trials = trials_per_item,
    incorrect = trials - correct,

    phase = factor(
      phase,
      levels = c("Pre", "Post")
    ),

    distribution = factor(
      distribution,
      levels = c("Unimodal", "Bimodal")
    ),

    participant = factor(participant),
    item = factor(item)
  )


# Current participant and group accuracies ---------------------------------

participant_scores <- pv %>%

  group_by(
    participant,
    distribution,
    phase
  ) %>%

  summarise(
    correct = sum(correct),
    trials = sum(trials),
    accuracy = correct / trials,
    .groups = "drop"
  )

print(participant_scores)


# Equal participant weighting
group_scores <- participant_scores %>%

  group_by(
    distribution,
    phase
  ) %>%

  summarise(
    mean_accuracy = mean(accuracy),
    .groups = "drop"
  )

print(group_scores)


# Fit current PV diagnostic model ------------------------------------------

# Items are treated as fixed effects because the earlier item-random-effect
# model had an estimated item variance of approximately zero and was singular.

model_pv_fixeditem <- glmer(

  cbind(correct, incorrect) ~

    phase * distribution +

    item +

    (1 | participant),

  data = pv,

  family = binomial,

  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

summary(model_pv_fixeditem)
isSingular(model_pv_fixeditem)
VarCorr(model_pv_fixeditem)


# Extract Phase x Distribution interaction ---------------------------------

interaction_name <- "phasePost:distributionBimodal"

coef(summary(model_pv_fixeditem))

current_interaction <- coef(summary(model_pv_fixeditem))[
  interaction_name,
  ,
  drop = FALSE
]

print(current_interaction)


# Optional likelihood-ratio test -------------------------------------------

model_pv_no_interaction <- glmer(

  cbind(correct, incorrect) ~

    phase +
    distribution +
    item +

    (1 | participant),

  data = pv,

  family = binomial,

  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

anova(
  model_pv_no_interaction,
  model_pv_fixeditem,
  test = "Chisq"
)


# Extract parameters for simulation ----------------------------------------

beta <- fixef(model_pv_fixeditem)
print(beta)

participant_sd <- as.numeric(
  attr(
    VarCorr(model_pv_fixeditem)$participant,
    "stddev"
  )[1]
)

cat(
  "\nParticipant random-intercept SD =",
  participant_sd,
  "\n\n"
)

item_levels <- levels(pv$item)
print(item_levels)


# Simulate one future PV study ---------------------------------------------

simulate_one_pv_study <- function(n_per_group) {

  # Future design:
  #   2 groups
  #   2 phases
  #   8 PV diagnostic items
  #   each item repeated twice
  #
  # = 16 PV diagnostic trials per participant per phase

  simdat <- expand_grid(

    distribution = c(
      "Unimodal",
      "Bimodal"
    ),

    person_number = seq_len(n_per_group),

    phase = c(
      "Pre",
      "Post"
    ),

    item = item_levels
  ) %>%

    mutate(

      participant = paste0(
        ifelse(
          distribution == "Bimodal",
          "PVB",
          "PVU"
        ),
        "_",
        sprintf("%02d", person_number)
      ),

      trials = 2,

      phase = factor(
        phase,
        levels = c("Pre", "Post")
      ),

      distribution = factor(
        distribution,
        levels = c("Unimodal", "Bimodal")
      ),

      item = factor(
        item,
        levels = item_levels
      ),

      participant = factor(participant)
    )


  # Participant random intercepts
  participant_ids <- levels(simdat$participant)

  participant_effects <- rnorm(
    n = length(participant_ids),
    mean = 0,
    sd = participant_sd
  )

  names(participant_effects) <- participant_ids


  # Fixed-effect design matrix
  X <- model.matrix(
    ~ phase * distribution + item,
    data = simdat
  )

  X <- X[
    ,
    names(beta),
    drop = FALSE
  ]


  # Linear predictor
  eta <-
    as.numeric(X %*% beta) +
    participant_effects[
      as.character(simdat$participant)
    ]

  probability_correct <- plogis(eta)


  # Generate two responses per item
  simdat$correct <- rbinom(
    n = nrow(simdat),
    size = simdat$trials,
    prob = probability_correct
  )

  simdat$incorrect <-
    simdat$trials -
    simdat$correct


  # Fit the same model to the simulated study
  fit <- tryCatch(

    suppressMessages(
      suppressWarnings(

        glmer(

          cbind(
            correct,
            incorrect
          ) ~

            phase * distribution +

            item +

            (1 | participant),

          data = simdat,

          family = binomial,

          control = glmerControl(
            optimizer = "bobyqa",
            optCtrl = list(
              maxfun = 100000
            )
          )
        )
      )
    ),

    error = function(e) NULL
  )


  if (is.null(fit)) {
    return(
      c(
        p = NA_real_,
        singular = NA_real_
      )
    )
  }


  # Check singularity
  singular_fit <- isSingular(
    fit,
    tol = 1e-4
  )


  # Extract interaction p-value
  coef_table <- coef(summary(fit))

  if (
    !(interaction_name %in%
      rownames(coef_table))
  ) {
    return(
      c(
        p = NA_real_,
        singular = as.numeric(singular_fit)
      )
    )
  }

  p_value <- coef_table[
    interaction_name,
    "Pr(>|z|)"
  ]

  return(
    c(
      p = p_value,
      singular = as.numeric(singular_fit)
    )
  )
}


# Power-estimation function -------------------------------------------------

power_pv_fixeditem <- function(
    n_per_group,
    nsim = 500
) {

  simulation_results <- replicate(
    nsim,
    simulate_one_pv_study(
      n_per_group
    )
  )

  simulation_results <- t(simulation_results)

  p_values <- simulation_results[, "p"]

  singular_values <- simulation_results[, "singular"]

  valid_p_values <- p_values[
    !is.na(p_values)
  ]

  valid_singular <- singular_values[
    !is.na(singular_values)
  ]

  tibble(

    PVB_N = n_per_group,

    PVU_N = n_per_group,

    total_PV_N = n_per_group * 2,

    simulations_requested = nsim,

    simulations_successful =
      length(valid_p_values),

    failed_simulations =
      nsim -
      length(valid_p_values),

    singular_simulations =
      sum(valid_singular == 1),

    singular_percent =
      100 *
      mean(valid_singular == 1),

    power =
      mean(valid_p_values < 0.05),

    power_percent =
      100 *
      mean(valid_p_values < 0.05)
  )
}


# Run selected sample sizes -------------------------------------------------

set.seed(123)

result_15 <- power_pv_fixeditem(
  n_per_group = 15,
  nsim = 500
)

print(
  result_15,
  width = Inf
)


set.seed(456)

result_18 <- power_pv_fixeditem(
  n_per_group = 18,
  nsim = 500
)

print(
  result_18,
  width = Inf
)


set.seed(789)

result_20 <- power_pv_fixeditem(
  n_per_group = 20,
  nsim = 500
)

print(
  result_20,
  width = Inf
)


# Combine power results -----------------------------------------------------

final_power_results <- bind_rows(
  result_15,
  result_18,
  result_20
) %>%

  mutate(

    power_percent = round(
      power_percent,
      1
    ),

    singular_percent = round(
      singular_percent,
      1
    )
  )

print(
  final_power_results,
  n = Inf,
  width = Inf
)


# For more stable final estimates, increase nsim from 500 to 1000.
# Running each sample size separately avoids fitting all models in one call.
